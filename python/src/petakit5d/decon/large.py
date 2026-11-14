"""Deconvolution of large files: in memory (``RLdecon_large_in_memory``) and in place
(``RLdecon_large_in_place``, whose block tasks are ``decon.block.decon_zarr_block``,
``RLdecon_for_zarr_block``).

Both cut the volume into batches (``core.blocks.zarr_chunk_coordinates``, MATLAB's
``XR_zarrChunkCoordinatesExtraction``) and deconvolve each batch with its border as RLdecon does
(``decon.file.rl_decon_file`` with the batch as ``rawdata`` and the batch's region as
``decon_bbox``), so a voxel's result depends on the batch it is in (the FFT's periodic boundary
is the batch's): the batches, their borders and their order are MATLAB's.

- **The border** comes from the PSF: the PSF file read, with ``psf_gen`` cleaned and resampled
  (``psf_gen_new``, 'masked', medFactor 1.5) to ``dz`` (``dz sind(skew_angle)`` with ``deskew``,
  the angle's sign kept) and cropped to its support (for every method), and written as
  ``<decon folder>/psfgen/<psf name>.tif`` if that is missing; the raw PSF without ``psf_gen``.
  The border is ``round((size + 10) / 2)`` per axis, at least ``edge_erosion + 1`` with an
  erosion. RLdecon generates the PSF for its decon again, per batch, with ``dz sind(|angle|)``.
- **In memory** (``decon_large_in_memory``): the file (tiff or zarr) is read whole (flipped in z
  with ``flip_zstack``); batches of ``batch_size`` (the block size is the batch size), not of
  the same size (``same_batch_size`` false), in MATLAB's order (y fastest, then x, then z); each
  batch's result (in the result's class) placed into a volume of zeros. With ``edge_erosion``
  each batch's decon erodes its own mask, and the whole result is then multiplied by the eroded
  mask of the whole volume's nonzero voxels: for data without zeros (3D, every size above twice
  the erosion) the volume without ``edge_erosion`` voxels at each face
  (``decon_mask_edge_erosion``), else ``bwdist(~mask) > edge_erosion - 1`` of the mask without
  its outer faces. ``save_maskfile`` keeps it as ``<decon folder>/Masks/<name>_eroded.tif``
  (uint8; an existing one is read instead). An ``erode_maskfile`` that exists sets the erosion to
  0, and the result is multiplied by its values (cast to the result's class). Then the z MIP
  ``MIPs/<name>_MIP_z.tif``, and the result (tiff, or zarr with ``save_zarr`` and ``block_size``)
  under a temporary name, moved into place last.
- **In place** (``decon_large_in_place``): the result is ``<decon_path>/<name>.zarr``; nothing is
  done if it exists. The batches are whole chunks of ``block_size`` with the same size
  (``same_batch_size``), sorted by their size and then their boxes as MATLAB's ``sortrows``
  does (by y start first), and run in tasks of 5 batches (20 for GPU jobs; with the cluster, a
  GPU job or 'omw' and more than 1,000 batches ``max(100, round(n / 5000))``) through the job
  framework into ``<result>_<uuid>``. A task's flag is
  ``<decon_path>/zarr_flag/<name>_<uuid>/blocks_<first>_<last>.mat`` (1-based) with the MIPs of
  its blocks, which ``tools.mip.mip_zarr`` merges into the MIP slabs and MIPs before the result
  is moved into place. The call records its writer (``jobs.writers``), waits for a writer of
  another call of the same result, keeps a result another writer published (also one published
  during its blocks or MIPs), and its blocks skip their work once it has ended. GPU jobs run
  their blocks with ``use_gpu`` (an error without a usable GPU), 4 CPUs each, without the master
  (the MIPs keep the caller's ``master_compute``), and with ``gpu_config_file``. The values the
  blocks get keep MATLAB's trip through its task strings (``%0.20d`` of the pixel sizes, angle,
  background and PSF step, ``%0.20f`` of the OMW parameters, ``%d`` of the damp and scale
  factors, offset and erosion, ``mat2str`` of the apodization bounds).
- **Masks**: with ``mask_fullpaths`` (2D zarr masks of the volume's xy, xz and yz projections:
  MATLAB's ``ny x nx``, ``nx x nz`` and ``ny x nz`` images) a batch whose region has no nonzero
  mask value in one of them is not deconvolved: its result is zeros.

The in-memory batches are only deconvolved, as the in-place blocks: their RLdecon gets
``deskew`` and ``rotate`` false (``deskew`` only sets the z step of the border PSF).

MATLAB quirks kept (each has a test): the in-memory temporary names carry an empty uuid as
``<name>_.tif`` (a new uuid here); the in-place blocks ignore ``flip_zstack`` and always fix the
iterations; the in-place call takes ``rl_method`` as it is (RLdecon's names are case-sensitive:
another name fails in each block's RLdecon); the border PSF of a deskewed file has a negative z
step with the default angle. ``in_place_controller_resources`` gives what the job of an
in-place file asks for (XR_decon_data_wrapper's controller rule).
"""

from __future__ import annotations

import math
import os
import time
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "DeconLargeError",
    "border_psf",
    "border_size",
    "decon_large_in_memory",
    "decon_large_in_place",
    "in_memory_batches",
    "in_place_batches",
    "in_place_controller_resources",
    "in_place_task_size",
    "mask_boxes",
    "masks_empty",
]

_UNDEFINED = "MATLAB:UndefinedFunction"


class DeconLargeError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id ('' for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


def _threads() -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(None)


def _exist(path: str) -> bool:
    from petakit5d.jobs.files import exist

    return exist(path)


def _scalar(v: Any, default: float = 0.0) -> float:
    a = np.asarray(default if v is None else v, dtype=np.float64).reshape(-1)
    return float(a[0]) if a.size else default


def _empty(v: Any) -> bool:
    return v is None or np.size(v) == 0


def _size3(v: Any, what: str) -> tuple[int, int, int]:
    a = np.asarray(v, dtype=np.float64).reshape(-1)
    if a.size != 3 or not np.all(np.isfinite(a)) or np.any(a != np.floor(a)):
        raise ValueError(f"{what} must be 3 integers (z, y, x), got {v!r}")
    return int(a[0]), int(a[1]), int(a[2])


def _toc(t0: float) -> None:
    print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.", flush=True)


def _num(x: float) -> int | float:
    """a MATLAB double as a Python int where it is integral (for the task records)"""
    x = float(x)
    return int(x) if math.isfinite(x) and x == math.floor(x) and abs(x) < 2**53 else x


def _dtype(save_16bit: bool) -> np.dtype:
    return np.dtype(np.uint16) if save_16bit else np.dtype(np.float32)


def _local_box(region: Sequence[int], batch: Sequence[int]) -> tuple[int, ...]:
    """the region of a batch inside the batch (RLdecon's deconBbox)"""
    s = [int(batch[q]) for q in range(3)]
    return tuple(int(region[q]) - s[q % 3] for q in range(6))


# ---------------------------------------------------------------------------------------------
# the border and the batches
# ---------------------------------------------------------------------------------------------


def border_psf(
    psf_fullpath: str,
    decon_path: str,
    *,
    psf_gen: bool = True,
    deskew: bool = False,
    skew_angle: float = -32.45,
    dz_psf: float = 0.1,
    dz: float = 0.1,
    threads: int | None = None,
) -> np.ndarray:
    """The PSF the large-file decon takes its border from: the PSF file, with ``psf_gen``
    cleaned and resampled to ``dz`` (``dz sind(skew_angle)`` with ``deskew``) and cropped to its
    support, written as ``<decon_path>/psfgen/<psf name>.tif`` if that is missing."""
    from petakit5d import io
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.decon import psf as P
    from petakit5d.jobs.publish import remove_path, rename_no_replace
    from petakit5d.matlabcompat.elementary import sind

    pp = io.read_tiff(psf_fullpath, num_threads=threads)
    if not psf_gen:
        return pp
    ratio = float(sind(float(skew_angle))) if deskew else 1.0
    psf = P.generate_psf(pp, float(dz_psf), float(dz) * ratio, 1.5, "masked", threads=threads)
    psf = P.crop_psf_support(psf)
    folder = f"{decon_path}/psfgen/"
    os.makedirs(folder, exist_ok=True)
    _, psf_fsn, _ = fileparts(psf_fullpath)
    final = f"{folder}/{psf_fsn}.tif"
    if not _exist(final):
        # (a temporary name: another call may write the same file)
        tmp = f"{folder}/{psf_fsn}_{get_uuid()}.tif"
        io.write_tiff(psf, tmp, num_threads=threads)
        if not rename_no_replace(tmp, final):
            remove_path(tmp)
    return psf


def border_size(psf_shape: Sequence[int], edge_erosion: float = 0) -> tuple[int, int, int]:
    """``round((size(psf) + 10) / 2)`` per axis (z, y, x), at least ``edge_erosion + 1`` with an
    erosion"""
    from petakit5d.core.dtypes import mround

    if len(psf_shape) != 3:
        raise ValueError(f"the PSF must be a (z, y, x) volume, not of shape {tuple(psf_shape)}")
    b = [float(mround((int(s) + 10) / 2)) for s in psf_shape]
    e = float(edge_erosion)
    if e > 0:
        b = [max(v, e + 1) for v in b]
    return _size3(b, "the border")


def in_memory_batches(
    im_size: Sequence[int], batch_size: Sequence[int], border: Sequence[int]
) -> tuple[np.ndarray, np.ndarray]:
    """The batches of the in-memory decon: (input boxes, output regions), our boxes in MATLAB's
    order (blocks of the batch size, batches not of the same size)."""
    from petakit5d.core.blocks import zarr_chunk_coordinates

    b, r, _ = zarr_chunk_coordinates(
        _size3(im_size, "im_size"),
        batch_size=_size3(batch_size, "batch_size"),
        block_size=_size3(batch_size, "batch_size"),
        same_batch_size=False,
        border_size=_size3(border, "border"),
    )
    return b, r


def in_place_batches(
    im_size: Sequence[int],
    batch_size: Sequence[int],
    block_size: Sequence[int],
    border: Sequence[int],
) -> tuple[np.ndarray, np.ndarray]:
    """The batches of the in-place decon: (input boxes, output regions), our boxes, batches of
    whole blocks with the same size, sorted as MATLAB's ``sortrows([prod(batch sizes),
    batchBBoxes])``: by size, then by the MATLAB box (y, x, z start, then y, x, z end), stable."""
    from petakit5d.core.blocks import zarr_chunk_coordinates

    b, r, _ = zarr_chunk_coordinates(
        _size3(im_size, "im_size"),
        batch_size=_size3(batch_size, "batch_size"),
        block_size=_size3(block_size, "block_size"),
        same_batch_size=True,
        border_size=_size3(border, "border"),
    )
    keys = []
    for k, row in enumerate(b.tolist()):
        z0, y0, x0, z1, y1, x1 = row
        vol = float(z1 - z0) * float(y1 - y0) * float(x1 - x0)
        keys.append((vol, y0, x0, z0, y1, x1, z1, k))
    order = [k[-1] for k in sorted(keys)]
    return b[order], r[order]


def in_place_task_size(num_batch: int, gpu_job: bool, parse_cluster: bool, rl_method: str) -> int:
    """batches per task of the in-place decon: 5 (20 for GPU jobs); with the cluster, a GPU job
    or 'omw' and more than 1,000 batches ``max(100, round(n / 5000))``"""
    from petakit5d.core.dtypes import mround

    size = 20 if gpu_job else 5
    if (parse_cluster or gpu_job or rl_method == "omw") and num_batch > 1000:
        size = max(100, int(mround(num_batch / 5000)))
    return size


def mask_boxes(region: Sequence[int]) -> list[tuple[int, ...]]:
    """the boxes of a region in the 2D masks (xy, xz, yz; MATLAB's ``ny x nx``, ``nx x nz`` and
    ``ny x nz`` images, stored as one-slice volumes): (0, y0, x0, 1, y1, x1), (0, x0, z0, 1, x1,
    z1), (0, y0, z0, 1, y1, z1)"""
    z0, y0, x0, z1, y1, x1 = (int(v) for v in region)
    return [(0, y0, x0, 1, y1, x1), (0, x0, z0, 1, x1, z1), (0, y0, z0, 1, y1, z1)]


def masks_empty(mask_fullpaths: Sequence[str], region: Sequence[int]) -> bool:
    """whether a batch's region has no nonzero value in one of the 2D masks (it is then not
    deconvolved), the masks read in MATLAB's order (xy, xz, yz) until one is empty"""
    from petakit5d import io

    for path, box in zip(mask_fullpaths[:3], mask_boxes(region), strict=False):
        m = io.read_zarr(os.fspath(path), input_bbox=box)
        if not np.any(m != 0):
            return True
    return False


def _use_masks(mask_fullpaths: Any) -> list[str]:
    """the masks, or none (MATLAB: ``isempty(maskFullpaths) || isempty(maskFullpaths{1})``)"""
    if mask_fullpaths is None or isinstance(mask_fullpaths, (str, os.PathLike)):
        masks = [] if not mask_fullpaths else [os.fspath(mask_fullpaths)]
    else:
        masks = [os.fspath(m) for m in mask_fullpaths]
    if not masks or not masks[0]:
        return []
    return masks


# ---------------------------------------------------------------------------------------------
# the in-memory decon
# ---------------------------------------------------------------------------------------------


def _read_data(path: str, threads: int) -> np.ndarray:
    from petakit5d.core.paths import fileparts

    _, _, ext = fileparts(path)
    from petakit5d import io

    if ext in (".tif", ".tiff"):
        return io.read_tiff(path, num_threads=threads)
    if ext == ".zarr":
        return io.read_zarr(path, num_threads=threads)
    raise DeconLargeError("Unrecognized function or variable 'im'.", _UNDEFINED)


def _all_true(v: np.ndarray) -> bool:
    """``all(v, 'all')``, plane by plane (stops at the first plane with a false voxel)"""
    for k in range(v.shape[0]):
        if not v[k].all():
            return False
    return True


def eroded_volume_mask(im: np.ndarray, edge_erosion: float, threads: int | None = None):
    """the in-memory decon's eroded mask of the nonzero voxels of the whole volume: for 3D data
    without zeros and every size above twice the erosion the volume without ``edge_erosion``
    voxels at each face, else ``bwdist(~mask) > edge_erosion - 1`` of the mask with its outer
    faces set false (an image: all false)"""
    from petakit5d.decon.erosion import decon_mask_edge_erosion
    from petakit5d.kernels import ops
    from petakit5d.matlabcompat.distance import bwdist

    e = float(edge_erosion)
    bw = ops.compare(im, "gt", 0.0, num_threads=threads)
    v = bw if bw.ndim == 3 else bw.reshape((1,) * (3 - bw.ndim) + bw.shape)
    if v.shape[0] > 1 and all(s > 2 * e for s in v.shape) and _all_true(v):
        return np.asarray(decon_mask_edge_erosion(v, e)).reshape(bw.shape)
    # (MATLAB's y, x and z faces; a 2D image's only plane is its first and last)
    v[:, [0, -1], :] = False
    v[:, :, [0, -1]] = False
    v[[0, -1], :, :] = False
    np.logical_not(v, out=v)
    d = bwdist(v, threads=threads)
    del v, bw
    return ops.compare(d, "gt", e - 1, num_threads=threads).reshape(im.shape)


def _mask_in_place(im: np.ndarray, mask: np.ndarray, threads: int) -> None:
    """``im = im .* cast(mask, class(im))`` in place (a product with 0 or 1 of the class), in
    z chunks on threads"""
    from concurrent.futures import ThreadPoolExecutor

    a = im.reshape((-1,) + im.shape[-2:]) if im.ndim >= 2 else im.reshape(1, 1, -1)
    m = mask.reshape(a.shape)
    n = a.shape[0]
    step = max(1, -(-n // max(1, threads)))
    chunks = [(k, min(n, k + step)) for k in range(0, n, step)]

    def run(c):
        np.multiply(a[c[0] : c[1]], m[c[0] : c[1]], out=a[c[0] : c[1]])

    if len(chunks) == 1:
        run(chunks[0])
        return
    with ThreadPoolExecutor(len(chunks)) as ex:
        list(ex.map(run, chunks))


def _write_tiff_replace(img: np.ndarray, final: str, threads: int) -> None:
    """writetiff of a file: a temporary name, then moved over the file (MATLAB writes it)"""
    from petakit5d import io
    from petakit5d.core.ids import get_uuid

    tmp = f"{final[:-4]}_{get_uuid()}.tif"
    io.write_tiff(img, tmp, num_threads=threads)
    os.replace(tmp, final)


def _rl_batch(
    rawdata: np.ndarray,
    output_fn: str,
    psf_fn: str,
    xy_pixel_size,
    dz,
    dz_psf,
    local_box,
    kw: dict[str, Any],
) -> np.ndarray:
    """RLdecon of one batch, as the large-file functions call it (no file written)"""
    from petakit5d.decon import file as F

    r = F.rl_decon_file(
        "",
        output_fn,
        psf_fn,
        xy_pixel_size,
        dz,
        dz_psf,
        rawdata=rawdata,
        dsr_combined=True,
        reverse=True,
        decon_bbox=local_box,
        debug=False,
        save_3d_stack=(False, False, False),
        mip_axis=(0, 0, 0),
        **kw,
    )
    return r.decon


@register(
    "RLdecon_large_in_memory",
    Param("frameFullpaths", "@(x) ischar(x) || iscell(x)", "File to deconvolve (tiff or zarr)."),
    Param(
        "deconFullpath",
        "@(x) ischar(x)",
        "PSF file (the second argument; MATLAB's inputParser names it deconFullpath).",
        py="psf_fullpath",
    ),
    Param(
        "psfFullpath",
        "@(x) ischar(x)",
        "Result file (the third argument; MATLAB's inputParser names it psfFullpath).",
        py="decon_fullpath",
    ),
    Param("xyPixelSize", "@isnumeric", "Pixel size (um)."),
    Param("dz", "@isnumeric", "Scan step (um)."),
    Param("save16bit", "@islogical", "Result in uint16 (else single).", convert="same"),
    Param("Rotate", "@islogical", "Accepted; not used (the batches are only deconvolved)."),
    Param("Deskew", "@islogical", "Border PSF at dz sind(SkewAngle); the batches aren't deskewed."),
    Param("SkewAngle", "@isnumeric", "Skew angle (degrees)."),
    Param("flipZstack", "@islogical", "Flip the data in z first."),
    Param("Background", "@isnumeric", "Background ([]: 99).", convert="same"),
    Param("EdgeSoften", "@isnumeric", "Accepted; not used (as in MATLAB)."),
    Param("zEdgeSoften", "@isnumeric", "Accepted; not used (as in MATLAB)."),
    Param("Crop", "@isnumeric", "Accepted; not used (as in MATLAB).", convert="same"),
    Param("zFlip", "@islogical", "Accepted; not used (as in MATLAB)."),
    Param(
        "GenMaxZproj",
        "@isnumeric",
        "Accepted; not used (as in MATLAB; (z, y, x)).",
        convert="flags_yxz",
    ),
    Param("dzPSF", "@isnumeric", "z step of the PSF (um)."),
    Param("DeconIter", "@isnumeric", "Number of iterations."),
    Param("EdgeErosion", "@isnumeric", "Erosion of the edges of the mask of nonzero data."),
    Param("ErodeMaskfile", "@ischar", "Eroded mask tiff multiplied with the result instead."),
    Param("SaveMaskfile", "@islogical", "Save the eroded mask (Masks/<name>_eroded.tif)."),
    Param("dampFactor", "@isnumeric", "Damp factor of the result."),
    Param("scaleFactor", "@isnumeric", "Scale factor of the result ([]: 1)."),
    Param("deconOffset", "@isnumeric", "Offset of the result.", convert="same"),
    Param("maskFullpaths", "@iscell", "2D masks (xy, xz, yz) of the regions to decon."),
    Param("saveZarr", "@islogical", "Save the result as zarr."),
    Param("RLMethod", "@ischar", "'simplified', 'omw', 'original' or 'cudagen'."),
    Param("wienerAlpha", "@isnumeric", "Wiener alpha of the OMW back projector."),
    Param("OTFCumThresh", "@isnumeric", "OTF cumulative sum threshold (OMW)."),
    Param("hannWinBounds", "@isnumeric", "Apodization range of the distance matrix (OMW)."),
    Param("skewed", "@(x) isempty(x) || islogical(x)", "Decon in skewed space (None: automatic)."),
    Param("fixIter", "@islogical", "Fixed number of iterations."),
    Param("errThresh", "@isnumeric", "Accepted; not used (as in MATLAB)."),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x).", convert="size_yxz"),
    Param(
        "blockSize",
        "@isnumeric",
        "Chunk size (z, y, x) of a zarr result (None: writezarr's).",
        convert="size_yxz",
    ),
    Param("Overlap", "@isnumeric", "Accepted; not used (as in MATLAB).", convert="same"),
    Param("CPUMaxMem", "@isnumeric", "Accepted; not used (as in MATLAB)."),
    Param("parseCluster", "@islogical", "Accepted; not used (as in MATLAB)."),
    Param("uuid", "@ischar", "uuid of the temporary names."),
    Param("debug", "@islogical", "Accepted; not used (the batches run without debug)."),
    Param("psfGen", "@islogical", "Generate the PSF (psf_gen) or use it as it is."),
    Param(
        "useGPU",
        "@(x) isempty(x) || islogical(x)",
        "True: on a GPU (an error if none is usable); False: the CPU; None: a GPU if usable.",
    ),
    task="R07",
)
def decon_large_in_memory(
    frame_fullpaths,
    psf_fullpath,
    decon_fullpath,
    xy_pixel_size,
    dz,
    *,
    save_16bit=True,
    rotate=False,
    deskew=False,
    skew_angle=-32.45,
    flip_zstack=False,
    background=None,
    edge_soften=5,
    z_edge_soften=2,
    crop=None,
    z_flip=False,
    gen_max_zproj=(1, 0, 0),
    dz_psf=0.1,
    decon_iter=15,
    edge_erosion=8,
    erode_maskfile="",
    save_maskfile=False,
    damp_factor=1,
    scale_factor=None,
    decon_offset=0,
    mask_fullpaths=(),
    save_zarr=False,
    rl_method="simplified",
    wiener_alpha=0.005,
    otf_cum_thresh=0.9,
    hann_win_bounds=(0.8, 1),
    skewed=None,
    fix_iter=True,
    err_thresh=None,
    batch_size=(1024, 1024, 1024),
    block_size=None,
    overlap=200,
    cpu_max_mem=500,
    parse_cluster=True,
    uuid="",
    debug=False,
    psf_gen=True,
    use_gpu=None,
) -> None:
    """Deconvolve a file that fits in memory batch by batch (RLdecon_large_in_memory): see the
    module. Sizes are (z, y, x)."""
    from petakit5d import io
    from petakit5d.core import dtypes
    from petakit5d.core.arrays import matlab_class
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.decon.file import _flip_z, _move
    from petakit5d.io import ReadBuffer
    from petakit5d.kernels.crop import crop_3d, indexing_3d
    from petakit5d.tools.mip import axis_max

    del edge_soften, z_edge_soften, crop, z_flip, gen_max_zproj, err_thresh, overlap
    del cpu_max_mem, parse_cluster, debug, rotate
    t = _threads()
    frame = os.fspath(frame_fullpaths)
    psf_fn = os.fspath(psf_fullpath)
    decon_fn = os.fspath(decon_fullpath)
    erosion = _scalar(edge_erosion)
    mask_file = os.fspath(erode_maskfile) if erode_maskfile else ""
    use_mask_file = bool(mask_file) and _exist(mask_file)
    if use_mask_file:
        erosion = 0.0
    bg = 99 if _empty(background) else background
    tag = uuid if uuid else get_uuid()

    t0 = time.perf_counter()
    _, fsname, _ = fileparts(frame)
    decon_path, _, _ = fileparts(decon_fn)
    print(f"Start Large-file RL Decon for {frame}...", flush=True)
    psf = border_psf(
        psf_fn,
        decon_path,
        psf_gen=bool(psf_gen),
        deskew=bool(deskew),
        skew_angle=skew_angle,
        dz_psf=dz_psf,
        dz=dz,
        threads=t,
    )

    t1 = time.perf_counter()
    print(f"reading {fsname}...", flush=True)
    im = _read_data(frame, t)
    if flip_zstack:
        im = _flip_z(im, True)
    _toc(t1)

    bw = None
    if erosion > 0:
        mask_fn = mask_tmp = ""
        if save_maskfile:
            mask_path = f"{decon_path}/Masks"
            os.makedirs(mask_path, exist_ok=True)
            mask_fn = f"{mask_path}/{fsname}_eroded.tif"
            mask_tmp = f"{mask_path}/{fsname}_eroded_{tag}.tif"
        if not (save_maskfile and _exist(mask_fn)):
            bw = eroded_volume_mask(im, erosion, t)
        else:
            from petakit5d.kernels import ops

            bw = ops.compare(io.read_tiff(mask_fn, num_threads=t), "gt", 0.0, num_threads=t)
        # a mask of the first time point for the others
        if save_maskfile:
            io.write_tiff(bw.view(np.uint8), mask_tmp, num_threads=t)
            os.replace(mask_tmp, mask_fn)

    if im.ndim != 3:
        raise ValueError(f"the large-file decon takes (z, y, x) volumes, not {im.shape}")
    im_size = im.shape
    dtype = _dtype(bool(save_16bit))
    border = border_size(np.shape(psf), erosion)
    in_boxes, out_boxes = in_memory_batches(im_size, batch_size, border)

    print("Deconvolving chunks...", flush=True)
    masks = _use_masks(mask_fullpaths)
    # the batches are only deconvolved, as the in-place blocks (deskew is for the PSF above)
    kw = {
        "save_16bit": bool(save_16bit),
        "skew_angle": skew_angle,
        "deskew": False,
        "rotate": False,
        "background": bg,
        "decon_iter": decon_iter,
        "rl_method": str(rl_method),
        "skewed": skewed,
        "wiener_alpha": wiener_alpha,
        "otf_cum_thresh": otf_cum_thresh,
        "hann_win_bounds": hann_win_bounds,
        "fix_iter": bool(fix_iter),
        "damp_factor": damp_factor,
        "scale_factor": scale_factor,
        "decon_offset": decon_offset,
        "use_gpu": use_gpu,
        "psf_gen": bool(psf_gen),
        "edge_erosion": _num(erosion),
    }
    n = in_boxes.shape[0]
    done = [False] * n
    imout = np.zeros(im_size, dtype)
    buf = ReadBuffer()
    for i in range(n):
        print(f"Process Batch {i + 1}... ", end="", flush=True)
        tb = time.perf_counter()
        ib = [int(v) for v in in_boxes[i]]
        ob = [int(v) for v in out_boxes[i]]
        # (the result is zeros where a mask is empty)
        if masks and masks_empty(masks, ob):
            done[i] = True
            _toc(tb)
            continue
        shape = tuple(ib[q + 3] - ib[q] for q in range(3))
        in_batch = crop_3d(im, ib, out=buf.view(shape, im.dtype), num_threads=t)
        out = _rl_batch(
            in_batch, decon_fn, psf_fn, xy_pixel_size, dz, dz_psf, _local_box(ob, ib), kw
        )
        del in_batch
        if out is None:
            # (RLdecon skipped the batch: MATLAB can't place its empty result)
            raise DeconLargeError(
                "Unable to perform assignment because the size of the left side is "
                f"{_msize(ob)} and the size of the right side is 0-by-0.",
                "MATLAB:subsassigndimmismatch",
            )
        out = np.asarray(out)
        # the mex takes a batch of the result's class: others are converted first
        if out.dtype != dtype:
            out = np.asarray(dtypes.cast(out, matlab_class(dtype)))
        indexing_3d(imout, ob, out.reshape(tuple(ob[q + 3] - ob[q] for q in range(3))))
        del out
        done[i] = True
        _toc(tb)
    buf.release()
    del im

    if not all(done):
        missing = [k + 1 for k, d in enumerate(done) if not d]
        raise DeconLargeError(f"Deconvolution of chunk {missing} is missing!")

    t3 = time.perf_counter()
    print("Saving combined deconvolved chunks...", flush=True)
    im = imout
    del imout
    if erosion > 0:
        print("Erode edges of deconvolved data w.r.t. raw data...", flush=True)
        _mask_in_place(im, bw, t)
    bw = None
    if use_mask_file:
        print("Erode edges of deconvolved data using a predefined mask...", flush=True)
        m = io.read_tiff(mask_file, num_threads=t)
        im = np.asarray(dtypes.times(im, dtypes.cast(m, matlab_class(im.dtype))))
        del m

    # the MIP comes first: the caller counts the file as done once the result exists
    mip_path = f"{decon_path}/MIPs/"
    os.makedirs(mip_path, exist_ok=True)
    mip = axis_max(im, -3, num_threads=t)
    _write_tiff_replace(mip.reshape(mip.shape[-2:]), f"{mip_path}{fsname}_MIP_z.tif", t)

    if save_zarr:
        tmp = f"{decon_fn[:-5]}_{tag}.zarr"
        if _empty(block_size):
            io.write_zarr(im, tmp, num_threads=t)
        else:
            io.write_zarr(im, tmp, block_size=_size3(block_size, "block_size"), num_threads=t)
    else:
        tmp = f"{decon_fn[:-4]}_{tag}.tif"
        io.write_tiff(im, tmp, num_threads=t)
    _move(tmp, decon_fn)
    _toc(t3)
    _toc(t0)


def _msize(box: Sequence[int]) -> str:
    """MATLAB's 'a-by-b-by-c' of a box's size"""
    z, y, x = (int(box[q + 3]) - int(box[q]) for q in range(3))
    return f"{y}-by-{x}-by-{z}"


# ---------------------------------------------------------------------------------------------
# the in-place decon
# ---------------------------------------------------------------------------------------------


def block_task_values(
    *,
    pixel_size,
    dz,
    skew_angle,
    background,
    dz_psf,
    decon_iter,
    wiener_alpha,
    otf_cum_thresh,
    hann_win_bounds,
    damp_factor,
    scale_factor,
    decon_offset,
    edge_erosion,
) -> dict[str, Any]:
    """The numbers the in-place blocks get, after MATLAB's trip through the task strings:
    ``%0.20d`` (pixel size, dz, angle, background, PSF step), ``%d`` (iterations, damp and
    scale factors, offset, erosion), ``%0.20f`` (OMW alpha and threshold), ``mat2str`` (the
    apodization bounds)."""
    from petakit5d.core.text import through_mat2str, through_text

    def tt(v, fmt):
        return _num(through_text(float(_scalar(v)), fmt))

    hann = np.atleast_1d(through_mat2str(np.asarray(hann_win_bounds, np.float64))).tolist()
    return {
        "xy_pixel_size": tt(pixel_size, "%0.20d"),
        "dz": tt(dz, "%0.20d"),
        "skew_angle": tt(skew_angle, "%0.20d"),
        "background": tt(background, "%0.20d"),
        "dz_psf": tt(dz_psf, "%0.20d"),
        "decon_iter": tt(decon_iter, "%d"),
        "wiener_alpha": tt(wiener_alpha, "%0.20f"),
        "otf_cum_thresh": tt(otf_cum_thresh, "%0.20f"),
        "hann_win_bounds": [_num(v) for v in hann],
        "damp_factor": tt(damp_factor, "%d"),
        "scale_factor": None if _empty(scale_factor) else tt(scale_factor, "%d"),
        "decon_offset": tt(decon_offset, "%d"),
        "edge_erosion": tt(edge_erosion, "%d"),
    }


def in_place_controller_resources(
    data_size: Sequence[int],
    batch_size: Sequence[int],
    block_size: Sequence[int],
    *,
    save_16bit: bool = True,
    gpu_job: bool = False,
    master_compute: bool = True,
    config_file: str = "",
    cpus_per_task: float = 4,
) -> tuple[float, float]:
    """(cpus_per_task, mem_allocate) of the job that runs an in-place large-file decon (the
    controller of its block jobs), XR_decon_data_wrapper's rule: the process and the merge of
    the MIP slabs (``tools.mip.merge_memory`` of the result) with one CPU, and the CPUs of a
    block's job if it deconvolves CPU blocks itself (``master_compute``, or the config file's
    ``masterCompute``; never for GPU jobs). Sizes are (z, y, x)."""
    from petakit5d.jobs.config import cluster_config_value
    from petakit5d.jobs.coordinator import coordinator_resources
    from petakit5d.tools.mip import merge_memory

    computes = (not gpu_job) and bool(
        cluster_config_value(config_file, "masterCompute", master_compute)
    )
    dtype = "uint16" if save_16bit else "single"
    merge = merge_memory(
        _size3(data_size, "data_size"),
        _size3(batch_size, "batch_size"),
        _size3(block_size, "block_size"),
        dtype,
    )
    return coordinator_resources(merge, computes, cpus_per_task, None)


@register(
    "RLdecon_large_in_place",
    Param("frameFullpaths", "@(x) ischar(x) || iscell(x)", "Zarr file to deconvolve."),
    Param("pixelSize", "@isnumeric", "Pixel size (um)."),
    Param("dz", "@isnumeric", "Scan step (um)."),
    Param("deconPath", "@(x) ischar(x) || isempty(x)", "Result folder."),
    Param("PSF", "@(x) ischar(x) || isempty(x)", "PSF file.", convert="same"),
    Param("save16bit", "@islogical", "Result in uint16 (else single).", convert="same"),
    Param("Deskew", "@islogical", "Border PSF at dz sind(SkewAngle)."),
    Param("SkewAngle", "@isnumeric", "Skew angle (degrees)."),
    Param("flipZstack", "@islogical", "Passed to the blocks, which ignore it (as in MATLAB)."),
    Param("Background", "@isnumeric", "Background ([]: 100).", convert="same"),
    Param("dzPSF", "@isnumeric", "z step of the PSF (um)."),
    Param("DeconIter", "@isnumeric", "Number of iterations."),
    Param("RLMethod", "@ischar", "'simplified', 'omw', 'original' or 'cudagen'."),
    Param("wienerAlpha", "@isnumeric", "Wiener alpha of the OMW back projector."),
    Param("OTFCumThresh", "@isnumeric", "OTF cumulative sum threshold (OMW)."),
    Param("hannWinBounds", "@isnumeric", "Apodization range of the distance matrix (OMW)."),
    Param("skewed", "@(x) isempty(x) || islogical(x)", "Decon in skewed space (None: automatic)."),
    Param("fixIter", "@islogical", "Passed to the blocks, which fix the iterations anyway."),
    Param("batchSize", "@isnumeric", "Batch size (z, y, x), in whole chunks.", convert="size_yxz"),
    Param("blockSize", "@isnumeric", "Chunk size (z, y, x) of the result.", convert="size_yxz"),
    Param("dampFactor", "@isnumeric", "Damp factor of the result."),
    Param("scaleFactor", "@isnumeric", "Scale factor of the result ([]: 1)."),
    Param("deconOffset", "@isnumeric", "Offset of the result.", convert="same"),
    Param("EdgeErosion", "@isnumeric", "Erosion of the edges of each batch's mask."),
    Param("maskFullpaths", "@iscell", "2D masks (xy, xz, yz) of the regions to decon."),
    Param("parseCluster", "@islogical", "Run the blocks as Slurm jobs."),
    Param("parseParfor", "@islogical", "Run the blocks on a local process pool."),
    Param("masterCompute", "@islogical", "This process computes blocks (CPU) and MIPs too."),
    Param("jobLogDir", "@ischar", "Folder of the job logs."),
    Param("cpusPerTask", "@isnumeric", "CPUs of a block's job (4 for GPU jobs)."),
    Param("GPUJob", "@islogical", "Deconvolve the blocks on GPUs."),
    Param("uuid", "@ischar", "Id of the call (empty: a new one)."),
    Param("debug", "@islogical", "Passed to the blocks (not used there)."),
    Param("psfGen", "@islogical", "Generate the PSF (psf_gen) or use it as it is."),
    Param(
        "mccMode", "@islogical", "MATLAB's compiled mode (MATLAB-only; batch jobs: at least one)."
    ),
    Param("configFile", "@ischar", "Cluster configuration file (CPU blocks, MIPs)."),
    Param("GPUConfigFile", "@ischar", "Cluster configuration file of the GPU blocks."),
    task="R07",
)
def decon_large_in_place(
    frame_fullpaths,
    pixel_size,
    dz,
    decon_path,
    psf,
    *,
    save_16bit=True,
    deskew=False,
    skew_angle=-32.45,
    flip_zstack=False,
    background=None,
    dz_psf=0.1,
    decon_iter=15,
    rl_method="simplified",
    wiener_alpha=0.005,
    otf_cum_thresh=0.9,
    hann_win_bounds=(0.8, 1),
    skewed=None,
    fix_iter=False,
    batch_size=(1024, 1024, 1024),
    block_size=(256, 256, 256),
    damp_factor=1,
    scale_factor=None,
    decon_offset=0,
    edge_erosion=0,
    mask_fullpaths=(),
    parse_cluster=True,
    parse_parfor=False,
    master_compute=True,
    job_log_dir="../job_logs",
    cpus_per_task=4,
    gpu_job=False,
    uuid="",
    debug=False,
    psf_gen=True,
    mcc_mode=False,
    config_file="",
    gpu_config_file="",
) -> None:
    """Deconvolve a large zarr file in blocks, each written into the result (RLdecon_large_in_
    place): see the module. Sizes are (z, y, x)."""
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts, mkdir_recursive
    from petakit5d.decon.block import decon_zarr_block
    from petakit5d.decon.file import _mkdir_group
    from petakit5d.io import create_zarr, get_image_size, get_zarr_info
    from petakit5d.jobs import framework, writers
    from petakit5d.jobs.publish import keep_other_result, publish_result
    from petakit5d.jobs.slurm import check_slurm_cluster
    from petakit5d.jobs.tasks import Task
    from petakit5d.tools import mip

    t = _threads()
    frame = os.fspath(frame_fullpaths)
    psf_fn = os.fspath(psf) if psf else ""
    decon_path = os.fspath(decon_path) if decon_path else ""
    bg = 100 if _empty(background) else background
    if not uuid:
        uuid = get_uuid()
    gpu_job = bool(gpu_job)
    caller_master = bool(master_compute)
    t0 = time.perf_counter()

    data_path, fsname, _ = fileparts(frame)
    decon_fullpath = f"{decon_path}/{fsname}.zarr"
    if os.path.isdir(decon_fullpath):
        print("The output result exists, skip it!", flush=True)
        return
    print(f"Start Large-file RL Decon for {frame}...", flush=True)
    tmp_path = f"{decon_fullpath}_{uuid}"

    # the PSF and its border
    p = border_psf(
        psf_fn,
        decon_path,
        psf_gen=bool(psf_gen),
        deskew=bool(deskew),
        skew_angle=skew_angle,
        dz_psf=dz_psf,
        dz=dz,
        threads=t,
    )
    print(f"reading {fsname}...", flush=True)
    dtype = "uint16" if save_16bit else "single"
    if parse_cluster:
        parse_cluster = bool(check_slurm_cluster(data_path, job_log_dir).parse_cluster)

    flag_path = f"{decon_path}/zarr_flag/{fsname}_{uuid}/"
    if not os.path.isdir(flag_path):
        mkdir_recursive(flag_path)
    # a writer of another call that may still make the result: wait for it; the record is
    # marked ended when this call ends but by publishing (its blocks are skipped then)
    with writers.recorded_writer(decon_fullpath, uuid, flag_path):
        if writers.wait(decon_fullpath, uuid):
            keep_other_result(decon_fullpath, tmp_path, flag_path)
            return

        im_size = _size3(tuple(get_image_size(frame))[:3], "the image size")
        border = border_size(np.shape(p), _scalar(edge_erosion))
        block = _size3(block_size, "block_size")
        in_boxes, out_boxes = in_place_batches(im_size, batch_size, block, border)

        if not os.path.isdir(tmp_path):
            n_chunks = 1
            for s, b in zip(im_size, block, strict=True):
                n_chunks *= math.ceil(s / b)
            create_zarr(
                tmp_path,
                data_size=im_size,
                block_size=block,
                dtype=dtype,
                dim_separator="/" if n_chunks > 10000 else ".",
            )
        # the tasks save the MIPs of their blocks for the MIP slabs (mip_zarr merges them
        # instead of reading the result again)
        zi = get_zarr_info(tmp_path)
        mbs = [int(v) for v in mip.mip_zarr_batch_size(zi.size, zi.block_size, zi.dtype)]

        num_batch = in_boxes.shape[0]
        task_size = in_place_task_size(num_batch, gpu_job, bool(parse_cluster), str(rl_method))
        cpus = cpus_per_task
        master = caller_master
        if gpu_job:
            cpus = 4
            master = False
        num_tasks = math.ceil(num_batch / task_size)
        values = block_task_values(
            pixel_size=pixel_size,
            dz=dz,
            skew_angle=skew_angle,
            background=bg,
            dz_psf=dz_psf,
            decon_iter=decon_iter,
            wiener_alpha=wiener_alpha,
            otf_cum_thresh=otf_cum_thresh,
            hann_win_bounds=hann_win_bounds,
            damp_factor=damp_factor,
            scale_factor=scale_factor,
            decon_offset=decon_offset,
            edge_erosion=edge_erosion,
        )
        xy_v, dz_v = values.pop("xy_pixel_size"), values.pop("dz")
        masks = _use_masks(mask_fullpaths)
        tasks, outputs = [], []
        for k in range(num_tasks):
            inds = list(range(k * task_size, min((k + 1) * task_size, num_batch)))
            flag = f"{flag_path}/blocks_{inds[0] + 1}_{inds[-1] + 1}.mat"
            outputs.append(flag)
            tasks.append(
                Task.of(
                    decon_zarr_block,
                    inds,
                    frame,
                    psf_fn,
                    tmp_path,
                    flag,
                    in_boxes[inds].tolist(),
                    out_boxes[inds].tolist(),
                    xy_v,
                    dz_v,
                    save_16bit=bool(save_16bit),
                    overwrite=False,
                    flip_zstack=bool(flip_zstack),
                    rl_method=str(rl_method),
                    skewed=skewed,
                    fix_iter=bool(fix_iter),
                    use_gpu=gpu_job,
                    mask_fullpaths=masks,
                    uuid=uuid,
                    debug=bool(debug),
                    psf_gen=bool(psf_gen),
                    mip_batch_size=mbs,
                    **values,
                )
            )
        inputs = [frame] * num_tasks
        done: list[bool] = []
        if not parse_parfor:
            cfg = gpu_config_file if gpu_job else config_file
            # (the blocks stop once the result exists: another writer made it)
            for _ in range(2):
                done = framework.run_tasks(
                    inputs,
                    outputs,
                    tasks,
                    cpus_per_task=cpus,
                    max_job_num=math.inf,
                    task_batch_num=1,
                    master_compute=master,
                    parse_cluster=bool(parse_cluster),
                    gpu_job=gpu_job,
                    job_log_dir=job_log_dir,
                    mcc_mode=mcc_mode,
                    config_file=cfg,
                    final_out_fullpath=decon_fullpath,
                )
                if all(done):
                    break
        else:
            done = framework.run_tasks(
                inputs,
                outputs,
                tasks,
                cluster_type="parfor",
                task_batch_num=1,
                gpu_job=gpu_job,
                uuid=uuid,
            )
        if not all(done):
            raise DeconLargeError("The decon of some chunks are missing!")

        # another writer may have published the result: keep it and its MIPs (the same
        # result), and remove this one's temporary result, which can be incomplete
        if keep_other_result(decon_fullpath, tmp_path, flag_path):
            return

        # the MIPs before the result is moved into place, named after the result and merged
        # from the MIPs of the blocks in the flags (the caller's master_compute)
        if os.path.isdir(tmp_path):
            _mkdir_group(f"{decon_path}/MIPs/")
            mip.mip_zarr(
                tmp_path,
                axis=(1, 1, 1),
                parse_cluster=bool(parse_cluster),
                master_compute=caller_master,
                job_log_dir=job_log_dir,
                mcc_mode=mcc_mode,
                config_file=config_file,
                block_flag_fullpaths=outputs,
            )
        # (one may also have published it during the MIPs)
        publish_result(tmp_path, decon_fullpath, flag_path)
    _toc(t0)
