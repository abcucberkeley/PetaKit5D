"""SIM reconstruction of files (``simReconFrame``), deskewing of phase-interleaved files
(``deskewPhasesFrame``) and the phase split and merge of ``third_parties/sim3D``
(``separatePhases``, ``combinePhases``).

``sim_recon_file``, per tiff or zarr file of raw phases:

1. the volume read (flipped along z with ``flip_zstack``: the page order, phases included);
2. with ``edge_erosion``, a mask of the data (``> 0``) eroded per plane by a disk (saved under
   ``<results>/Masks/`` with ``save_maskfile``) and applied to the raw data (``erode_before``)
   or, resized to the reconstruction (nearest), to the result (``erode_after``);
3. the volume cut into overlapping chunks (``XR_subVolumeCoordinatesExtraction`` on the size
   of one phase, ``chunk_size`` and ``overlap``; their x and y sizes made even), every chunk with
   enough signal (sum above ``int_thresh`` and a nonzero fraction of at least ``occ_thresh``)
   reconstructed (``sim_recon_volume``) and its center written into the result (twice the
   lateral size; a chunk that fails is reported and left 0);
4. ``<results>/<name>_recon.tif`` (single, or uint16 with ``save_16bit``) with negative values
   set to 0, and ``<results>/MIPs/<name>_MIP_z.tif``, the z MIP.

MATLAB quirks kept: an existing ``erode_maskfile`` turns the erosion off altogether (the file is
never read); the MIP is that of the result before negative values are set to 0; the chunks of
x (or y) are all widened by one voxel when any of them is odd; with ``use_gpu`` and not ``ds``
the chunks span the whole width in x (MATLAB's chunk shape for skewed data); the parameters
``edge_soften``, ``z_edge_soften``, ``crop``, ``z_flip``, ``gen_max_zproj``, ``resize_images``,
``cpu_max_mem``, ``large_file``, ``parse_cluster``, ``master_compute``, ``master_cpu``,
``job_log_dir``, ``cpus_per_task``, ``max_trial_num`` and ``unit_wait_time`` are accepted and
not used. The port writes the MIP before it publishes the result (by one rename), where MATLAB
writes it after.
"""

from __future__ import annotations

import os
import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, not_ported, register
from petakit5d.sim.recon import ARRAY_OR_PATH

not_ported(
    "DGXsimReconFrame",
    "A development copy of simReconFrame for one multi-GPU node (its chunks on parallel MATLAB "
    "workers) whose edge erosion is computed but never applied; sim_recon_file does its work.",
    "sim",
)

__all__ = [
    "combine_phases_file",
    "deskew_phases_file",
    "separate_phases_file",
    "sim_recon_file",
]


def _as_list(paths) -> list[str]:
    if isinstance(paths, (str, os.PathLike)):
        return [os.fspath(paths)]
    return [os.fspath(p) for p in paths]


def _read_volume(path: str, ext: str) -> np.ndarray:
    from petakit5d.sim.recon import SimError

    if ext in (".tif", ".tiff"):
        from petakit5d.io import read_tiff

        a = np.asarray(read_tiff(path))
    elif ext == ".zarr":
        from petakit5d.io import read_zarr

        a = np.asarray(read_zarr(path))
    else:
        raise SimError("Unrecognized function or variable 'im_raw'.", "MATLAB:UndefinedFunction")
    return a if a.ndim == 3 else a[None]


def _even_chunks(lo: np.ndarray, hi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """MATLAB's ``if any(rem(hi - lo + 1, 2) == 1), lo = lo - 1; hi = hi + (lo == 0);
    lo = lo + (lo == 0); end`` (every chunk moved when any is odd)"""
    if np.any(np.remainder(hi - lo + 1, 2) == 1):
        lo = lo - 1
        hi = hi + (lo == 0)
        lo = lo + (lo == 0)
    return lo, hi


def _eroded_mask(im_raw: np.ndarray, radius: float) -> np.ndarray:
    """``imerode`` of the data's nonzero mask, padded by one false voxel, by ``strel('disk', r)``"""
    from petakit5d.matlabcompat.morphology import imerode, strel

    bw = np.zeros(tuple(s + 2 for s in im_raw.shape), dtype=bool)
    bw[1:-1, 1:-1, 1:-1] = im_raw > 0
    er = imerode(bw, strel("disk", radius))
    return np.ascontiguousarray(er[1:-1, 1:-1, 1:-1])


def _times_mask(a: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """``a .* cast(mask, class(a))``"""
    return a * mask.astype(a.dtype)


def _recon_to_matlab(v):
    from petakit5d.core.arrays import to_matlab

    return None if v is None else to_matlab(np.asarray(v))


@register(
    "simReconFrame",
    Param("frameFullpaths", "@(x) ischar(x) || iscell(x)", "Tiff or zarr files of raw phases."),
    Param(
        "otf",
        None,
        "The OTF (sim_psf_to_otf: (norientations, norders, z, y, x)) or the PSF's tiff file.",
        convert=ARRAY_OR_PATH,
    ),
    Param("resultsDirName", "@ischar", "Result folder next to each file."),
    Param("islattice", "@islogical", "Lattice light sheet data (else 3-beam SIM)."),
    Param("NA_det", "@isnumeric", "Detection NA."),
    Param("NA_ext", "@isnumeric", "Excitation NA."),
    Param("nimm", "@isnumeric", "Refractive index of the immersion medium."),
    Param("wvl_em", "@isnumeric", "Emission wavelength in um."),
    Param("wvl_ext", "@isnumeric", "Excitation wavelength in um."),
    Param("w", "@isnumeric", "Wiener coefficient for regularization."),
    Param("apodize", "@islogical", "Apodize the final spectrum."),
    Param("DS", "@islogical", "The data are deskewed (else chunks span x on the GPU)."),
    Param("nphases", "@isnumeric", "Phases."),
    Param("norders", "@isnumeric", "Information orders."),
    Param("norientations", "@isnumeric", "Pattern orientations."),
    Param("lattice_period", "@isnumeric", "Lattice period in um (the coarsest)."),
    Param("lattice_angle", "@isnumeric", "Angle of the pattern modulation per orientation (rad)."),
    Param("phase_step", "@isnumeric", "Phase step in um."),
    Param(
        "pxl_dim_data", "@isnumeric", "Voxel size of the data (z, y, x) in um.", convert="size_yxz"
    ),
    Param(
        "pxl_dim_PSF", "@isnumeric", "Voxel size of the PSF (z, y, x) in um.", convert="size_yxz"
    ),
    Param("Background", "@isnumeric", "Background subtracted from data and PSF.", convert="same"),
    Param("normalize_orientations", "@islogical", "Normalize the total intensity per orientation."),
    Param("perdecomp", "@islogical", "Periodic-plus-smooth decomposition against edge effects."),
    Param("edgeTaper", "@islogical", "Taper the edges (Tukey window) against edge effects."),
    Param("edgeTaperVal", "@isnumeric", "Roll-off of the Tukey window."),
    Param("useGPU", "@islogical", "Run on a GPU if one is usable."),
    Param("Overwrite", "@islogical", "Replace existing results."),
    Param("Save16bit", "@islogical", "Save the result as uint16 (else single)."),
    Param(
        "gpuPrecision", "@ischar", "Precision of GPU arrays and phase ramps: 'single', 'double'."
    ),
    Param("flipZstack", "@islogical", "Flip the raw pages along z first."),
    Param("EdgeSoften", "@isnumeric", "Not used."),
    Param("zEdgeSoften", "@isnumeric", "Not used."),
    Param("Crop", "@isnumeric", "Not used."),
    Param("zFlip", "@islogical", "Not used."),
    Param(
        "GenMaxZproj", "@isnumeric", "Not used (the z MIP is always saved).", convert="flags_yxz"
    ),
    Param("ResizeImages", "@isnumeric", "Not used."),
    Param("EdgeErosion", "@isnumeric", "Radius of the disk that erodes the data's mask (0: none)."),
    Param("ErodeBefore", "@islogical", "Apply the eroded mask to the raw data."),
    Param("ErodeAfter", "@islogical", "Apply the eroded mask, resized, to the result."),
    Param("ErodeMaskfile", "@ischar", "An existing mask file turns the erosion off (MATLAB's)."),
    Param("SaveMaskfile", "@islogical", "Save the eroded mask under <results>/Masks/."),
    Param("ChunkSize", "@isvector", "Chunk size (z, y, x) of one phase.", convert="size_yxz"),
    Param("Overlap", "@isnumeric", "Chunk overlap in voxels.", convert="same"),
    Param("maxSubVolume", "@isnumeric", "Largest chunk in voxels."),
    Param("CPUMaxMem", "@isnumeric", "Not used."),
    Param("largeFile", "@islogical", "Not used."),
    Param("parseCluster", "@islogical", "Not used."),
    Param("masterCompute", "@islogical", "Not used."),
    Param("masterCPU", "@islogical", "Not used."),
    Param("jobLogDir", "@ischar", "Not used."),
    Param("cpusPerTask", "@isnumeric", "Not used."),
    Param("uuid", "@ischar", "Id of the temporary files (empty: a new one)."),
    Param("maxTrialNum", "@isnumeric", "Not used."),
    Param("unitWaitTime", "@isnumeric", "Not used."),
    Param("intThresh", "@isnumeric", "A chunk is reconstructed if its sum is above this."),
    Param("occThresh", "@isnumeric", "... and its nonzero fraction at least this."),
    Param("num_threads", None, "Threads (None: all available).", python_only=True),
    Param("fft_backend", None, "CPU FFT backend (None: the default).", python_only=True),
    task="C06",
    matlab_out=_recon_to_matlab,
)
def sim_recon_file(
    frame_fullpaths,
    otf,
    *,
    results_dir_name: str = "sim_recon",
    islattice: bool = True,
    na_det: float = 1.0,
    na_ext: float = 0.55,
    nimm: float = 1.33,
    wvl_em: float = 0.605,
    wvl_ext: float = 0.56,
    w: float = 0.005,
    apodize: bool = True,
    ds: bool = False,
    nphases: float = 5,
    norders: float = 5,
    norientations: float = 1,
    lattice_period: float = 1.2021,
    lattice_angle: Any = 1.5707963267948966,
    phase_step: float = 0.232,
    pxl_dim_data: Sequence[float] = (0.16074803849369898, 0.11, 0.11),
    pxl_dim_psf: Sequence[float] = (0.10716535899579933, 0.11, 0.11),
    background: float = 105,
    normalize_orientations: bool = False,
    perdecomp: bool = False,
    edge_taper: bool = False,
    edge_taper_val: float = 0.1,
    use_gpu: bool = True,
    overwrite: bool = False,
    save_16bit: bool = False,
    gpu_precision: str = "single",
    flip_zstack: bool = False,
    edge_soften: float = 5,
    z_edge_soften: float = 2,
    crop: Any = None,
    z_flip: bool = False,
    gen_max_zproj: Sequence[float] = (1, 0, 0),
    resize_images: Any = None,
    edge_erosion: float = 0,
    erode_before: bool = False,
    erode_after: bool = True,
    erode_maskfile: str = "",
    save_maskfile: bool = False,
    chunk_size: Sequence[float] = (250, 250, 250),
    overlap: float = 10,
    max_sub_volume: float = 300000000,
    cpu_max_mem: float = 500,
    large_file: bool = False,
    parse_cluster: bool = True,
    master_compute: bool = True,
    master_cpu: bool = False,
    job_log_dir: str = "../job_logs",
    cpus_per_task: float = 5,
    uuid: str = "",
    max_trial_num: float = 3,
    unit_wait_time: float = 2,
    int_thresh: float = 1,
    occ_thresh: float = 0.8,
    num_threads: int | None = None,
    fft_backend: str | None = None,
):
    """SIM reconstructions of files of raw phases (``simReconFrame``) into
    ``<folder>/<results_dir_name>/<name>_recon.tif`` and its z MIP; returns the last volume
    reconstructed ``(nz, 2 ny, 2 nx)`` (None when every file was skipped). ``otf`` is an OTF
    array or the PSF's tiff file (``sim_psf_to_otf`` with the same values)."""
    from petakit5d.core import dtypes
    from petakit5d.core.blocks import sub_volume_coordinates_matlab
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.core.text import sprintf
    from petakit5d.io import write_tiff
    from petakit5d.io.tiff import group_writable
    from petakit5d.jobs.publish import publish_result
    from petakit5d.sim.recon import SimError, SimWarning, sim_recon_volume

    del edge_soften, z_edge_soften, crop, z_flip, gen_max_zproj, resize_images, cpu_max_mem
    del large_file, parse_cluster, master_compute, master_cpu, job_log_dir, cpus_per_task
    del max_trial_num, unit_wait_time
    frames = _as_list(frame_fullpaths)
    gpu = False
    if use_gpu:
        from petakit5d.gpu.device import usable_gpu_count

        gpu = usable_gpu_count()[0] > 0
    if isinstance(otf, (str, os.PathLike)):
        empty = os.fspath(otf) == ""
    else:
        empty = otf is None or np.size(otf) == 0
    if empty:
        raise SimError("You should provide an otf file for the frame...\n")
    common = {
        "nphases": nphases,
        "norders": norders,
        "norientations": norientations,
        "lattice_period": lattice_period,
        "phase_step": phase_step,
        "pxl_dim_psf": pxl_dim_psf,
        "background": background,
        "use_gpu": gpu,
        "gpu_precision": gpu_precision,
        "num_threads": num_threads,
        "fft_backend": fft_backend,
    }
    if isinstance(otf, (str, os.PathLike)):
        from petakit5d.sim.otf import sim_psf_to_otf

        otf = sim_psf_to_otf(os.fspath(otf), save_16bit=bool(save_16bit), **common)
    if erode_maskfile and os.path.exists(erode_maskfile):
        edge_erosion = 0
    if not uuid:
        uuid = get_uuid()
    steps = float(nphases) * float(norientations)
    recon_kw = dict(
        islattice=islattice, na_det=na_det, na_ext=na_ext, nimm=nimm, wvl_em=wvl_em,
        wvl_ext=wvl_ext, w=w, apodize=apodize, lattice_angle=lattice_angle,
        pxl_dim_data=pxl_dim_data, normalize_orientations=normalize_orientations,
        perdecomp=perdecomp, edge_taper=edge_taper, edge_taper_val=edge_taper_val, **common,
    )  # fmt: skip
    # MATLAB's [y x z] chunk size
    cs = np.asarray(chunk_size, dtype=np.float64).ravel()
    cs = cs[[1, 2, 0]] if cs.size == 3 else cs
    out_vol = None
    for frame in frames:
        if not os.path.exists(frame):
            warnings.warn(sprintf("%s does not exist! Skip it!", frame), stacklevel=2)
            continue
        pathstr, fsname, ext = fileparts(frame)
        recon_path = f"{pathstr}/{results_dir_name}/"
        os.makedirs(recon_path, exist_ok=True)
        group_writable(recon_path)
        recon_full = f"{recon_path}/{fsname}_recon.tif"
        if os.path.exists(recon_full) and not overwrite:
            print("Reconstruction results already exist, skip it!", flush=True)
            continue
        im_raw = _read_volume(frame, ext)
        if flip_zstack:
            im_raw = np.ascontiguousarray(im_raw[::-1])
        im_bw_erode = None
        if edge_erosion > 0:
            print("Create eroded masks using raw data...", flush=True)
            im_bw_erode = _eroded_mask(im_raw, float(edge_erosion))
            print("Eroding edges for data...", flush=True)
            if erode_before:
                im_raw = _times_mask(im_raw, im_bw_erode)
            if save_maskfile:
                mask_path = f"{recon_path}/Masks"
                os.makedirs(mask_path, exist_ok=True)
                tmp = f"{mask_path}/{fsname}_eroded_{uuid}.tif"
                write_tiff(im_bw_erode.astype(np.uint8), tmp)
                os.replace(tmp, f"{mask_path}/{fsname}_eroded.tif")
            if erode_after:
                from petakit5d.matlabcompat.resize import imresize3

                nzm, nym, nxm = im_bw_erode.shape
                size = (nzm / steps, nym / 0.5, nxm / 0.5)
                # (nearest picks voxels: the logical mask resized as 0/1 uint8 is the same)
                im_bw_erode = imresize3(im_bw_erode.astype(np.uint8), size, "nearest") > 0
        msize = np.array([im_raw.shape[1], im_raw.shape[2], im_raw.shape[0]], dtype=np.float64)
        dtype = np.uint16 if save_16bit else np.float32
        ck_size = cs.copy()
        if not ds and gpu:
            ck_size[0] = float(dtypes.mround(np.float64(ck_size[0] / msize[1] * ck_size[1])))
            ck_size[1] = msize[1]
        xmin, xmax, ymin, ymax, zmin, zmax, nn = sub_volume_coordinates_matlab(
            msize / np.array([1, 1, steps]),
            chunk_size=ck_size,
            overlap_size=overlap,
            max_sub_volume=max_sub_volume,
        )
        xmin, xmax = _even_chunks(np.asarray(xmin, np.float64), np.asarray(xmax, np.float64))
        ymin, ymax = _even_chunks(np.asarray(ymin, np.float64), np.asarray(ymax, np.float64))
        print("Processing image chunks...", flush=True)
        nz_out = msize[2] / steps
        if nz_out != int(nz_out):
            raise SimError("Size inputs must be integers.", "MATLAB:NonIntegerInput")
        im = np.zeros((int(nz_out), int(msize[0] * 2), int(msize[1] * 2)), dtype=dtype)
        isz = np.array([im.shape[1], im.shape[2], im.shape[0]], dtype=np.float64)
        lol, rol = int(np.floor(overlap / 2)), int(np.ceil(overlap / 2))
        for ck in range(int(nn)):
            z0 = int((zmin[ck] - 1) * steps)
            z1 = int(zmax[ck] * steps)
            chunk = im_raw[
                z0:z1, int(ymin[ck]) - 1 : int(ymax[ck]), int(xmin[ck]) - 1 : int(xmax[ck])
            ]
            total = float(np.asarray(_msum(chunk)))
            if total > int_thresh and np.count_nonzero(chunk) / chunk.size >= occ_thresh:
                print(sprintf("processing chunk:%d of %d \n", ck + 1, nn), end="", flush=True)
                try:
                    tim = sim_recon_volume(chunk, otf, **recon_kw)
                except Exception:  # noqa: BLE001 (MATLAB's try/catch around the chunk)
                    print(sprintf("chunk %d failed!!!!!!!!!!!!\n", ck + 1), end="", flush=True)
                    continue
                _place(im, tim, ck, xmin, xmax, ymin, ymax, zmin, zmax, isz, lol, rol)
        if edge_erosion > 0 and erode_after:
            if im_bw_erode.shape != im.shape:
                raise SimError(
                    "Arrays have incompatible sizes for this operation.",
                    "MATLAB:sizeDimensionsMustMatch",
                )
            im = _times_mask(im, im_bw_erode)
        tmp = sprintf("%s_%s_eroded.tif", recon_full[:-4], uuid)
        mip_path = f"{recon_path}/MIPs/"
        os.makedirs(mip_path, exist_ok=True)
        group_writable(mip_path)
        # the z MIP of the result before negative values are set to 0 (first, then the result)
        from petakit5d.matlabcompat import reduce

        mip = np.ascontiguousarray(reduce.max(im, axis=-3)[0])
        write_tiff(mip, f"{mip_path}{fsname}_MIP_z.tif")
        if save_16bit:
            write_tiff(im, tmp)
        else:
            write_tiff(im * (im >= 0), tmp)
        publish_result(tmp, recon_full)
        # MATLAB deletes temporary files that were never written: two warnings
        for never in (
            sprintf("%s_%s_recon.tif", recon_full[:-10], uuid),
            sprintf("%s/%s_%s_MIP_z.tif", recon_path, fsname, uuid),
        ):
            warnings.warn(
                SimWarning(f"File '{never}' not found.", "MATLAB:DELETE:FileNotFound"),
                stacklevel=2,
            )
        out_vol = im
        if os.path.exists(recon_full):
            print(sprintf("Completed sim reconstruction of %s.\n", frame), end="", flush=True)
    return out_vol


def _msum(chunk: np.ndarray):
    """``sum(im_chunk(:))`` (double for integers)"""
    from petakit5d.kernels.ops import sum_double

    if chunk.dtype == np.float32:
        from petakit5d.matlabcompat import reduce

        return np.asarray(reduce.sum(chunk)).ravel()[0]
    return sum_double(chunk)


def _place(im, tim, ck, xmin, xmax, ymin, ymax, zmin, zmax, isz, lol, rol) -> None:
    """``im(yrange, xrange, zrange) = tim(tyrange, txrange, tzrange)``: the chunk's center (its
    overlap halves dropped where it meets other chunks), in ``im``'s class"""
    from petakit5d.core import dtypes

    tsz, tsy, tsx = tim.shape
    ymin_ck = ymin[ck] * 2 - 1
    xmin_ck = xmin[ck] * 2 - 1

    def rng(lo, hi, lo_ck, n_out, half_l, half_r, tsize):
        a = lo_ck + (lo_ck != 1) * half_l
        b = hi - (hi != n_out) * half_r
        r = np.arange(a, b + 1, dtype=np.float64)
        return r[r - lo_ck + 1 <= tsize]

    yr = rng(ymin_ck, ymax[ck] * 2, ymin_ck, isz[0], lol * 2, rol * 2, tsy)
    xr = rng(xmin_ck, xmax[ck] * 2, xmin_ck, isz[1], lol * 2, rol * 2, tsx)
    zr = rng(zmin[ck], zmax[ck], zmin[ck], isz[2], lol, rol, tsz)
    if not (yr.size and xr.size and zr.size):
        return
    ty = (yr - ymin_ck).astype(np.int64)
    tx = (xr - xmin_ck).astype(np.int64)
    tz = (zr - zmin[ck]).astype(np.int64)
    sub = tim[np.ix_(tz, ty, tx)]
    if im.dtype != sub.dtype:
        sub = (
            np.asarray(dtypes.cast(sub, "uint16"))
            if im.dtype == np.uint16
            else sub.astype(im.dtype)
        )
    im[np.ix_(zr.astype(np.int64) - 1, yr.astype(np.int64) - 1, xr.astype(np.int64) - 1)] = sub


# ---------------------------------------------------------------------------------------------
# deskewPhasesFrame, separatePhases, combinePhases
# ---------------------------------------------------------------------------------------------


@register(
    "deskewPhasesFrame",
    Param("dataFile", None, "Tiff file of raw phases (interleaved per z plane)."),
    Param("xyPixelSize", None, "Pixel size in um."),
    Param("dz", None, "Scan step in um."),
    Param("SkewAngle", "@isscalar", "Skew angle in degrees."),
    Param("Reverse", "@islogical", "Reverse scan direction."),
    Param("Rotate", "@islogical", "Rotate after the deskew (DSR)."),
    Param("nphases", "@isnumeric", "Phases."),
    Param("Save16bit", "@islogical", "Save as uint16 (else single)."),
    task="C06",
)
def deskew_phases_file(
    data_file,
    xy_pixel_size,
    dz,
    skew_angle: float = 32.45,
    reverse: bool = True,
    *,
    rotate: bool = False,
    nphases: float = 5,
    save_16bit: bool = False,
) -> None:
    """Deskew a tiff file of raw phases phase by phase (``deskewPhasesFrame``): each phase's
    planes ``deskew_volume``'d, interleaved again (single), rotated with ``rotate``
    (``rotate_volume``, cropped, objective scan) and written to ``<folder>/DS/<name>`` (``DSR``
    with ``rotate``), uint16 with ``save_16bit``."""
    import time

    from petakit5d.core import dtypes
    from petakit5d.core.paths import fileparts
    from petakit5d.core.text import sprintf
    from petakit5d.deskew_rotate.volume import deskew_volume, rotate_volume
    from petakit5d.io import read_tiff, write_tiff
    from petakit5d.io.tiff import group_writable
    from petakit5d.matlabcompat import elementary as E

    if not isinstance(data_file, (str, os.PathLike)):
        data_file = "".join(_as_list(data_file))
    data_file = os.fspath(data_file)
    filepath, name, ext = fileparts(data_file)
    fol = f"{filepath}/"
    fol_str = "DSR" if rotate else "DS"
    if not os.path.isdir(fol + fol_str):
        os.makedirs(fol + fol_str, exist_ok=True)
        group_writable(fol + fol_str)
    t0 = time.perf_counter()
    im = np.asarray(read_tiff(data_file))
    if im.ndim == 2:
        im = im[None]
    nph = int(nphases)
    dsim = None
    for p in range(nph):
        vol = np.ascontiguousarray(im[p::nph])
        out = deskew_volume(vol, skew_angle, dz, xy_pixel_size, reverse)
        if dsim is None:
            dsim = np.zeros((im.shape[0], out.shape[1], out.shape[2]), dtype=np.float32)
        dsim[p::nph] = out
    if rotate:
        theta = skew_angle * np.pi / 180
        dz0 = float(E.sin(np.float64(theta))) * dz
        z_aniso = dz0 / xy_pixel_size
        dsim = rotate_volume(dsim, skew_angle, z_aniso, reverse, crop=True, objective_scan=True)
    dsim = np.asarray(dtypes.to_class(dsim, "uint16")) if save_16bit else dsim.astype(np.float32)
    write_tiff(dsim, f"{fol}{fol_str}/{name}{ext}")
    print(sprintf("Elapsed time is %f seconds.\n", time.perf_counter() - t0), end="", flush=True)


@register(
    "separatePhases",
    Param("dataFile", None, "Tiff file of interleaved phases (path ending in .tif)."),
    Param("nphases", "@isnumeric", "Phases."),
    Param("writeToDisk", "@islogical", "Not used."),
    task="C06",
)
def separate_phases_file(data_file, *, nphases: float = 5, write_to_disk: bool = False) -> None:
    """Write each phase of a tiff file of interleaved phases to ``<name>_phase<p>.tif``
    (``separatePhases``)."""
    from petakit5d.io import read_tiff, write_tiff

    del write_to_disk
    data_file = os.fspath(data_file)
    data = np.asarray(read_tiff(data_file))
    if data.ndim == 2:
        data = data[None]
    nph = int(nphases)
    for p in range(1, nph + 1):
        write_tiff(np.ascontiguousarray(data[p - 1 :: nph]), f"{data_file[:-4]}_phase{p}.tif")


@register(
    "combinePhases",
    Param("dataFile", None, "The name the phase files derive from (<name>_phase<p>.tif)."),
    Param("nphases", "@isnumeric", "Phases."),
    task="C06",
)
def combine_phases_file(data_file, *, nphases: float = 5) -> None:
    """Interleave ``<name>_phase<p>.tif`` again into ``<name>_combined.tif`` (single;
    ``combinePhases``)."""
    from petakit5d.io import read_tiff, write_tiff

    data_file = os.fspath(data_file)
    nph = int(nphases)
    first = np.asarray(read_tiff(f"{data_file[:-4]}_phase1.tif"))
    first = first if first.ndim == 3 else first[None]
    out = np.zeros((first.shape[0] * nph, first.shape[1], first.shape[2]), dtype=np.float64)
    for p in range(1, nph + 1):
        c = first if p == 1 else np.asarray(read_tiff(f"{data_file[:-4]}_phase{p}.tif"))
        c = c if c.ndim == 3 else c[None]
        out[p - 1 :: nph] = c
    write_tiff(out.astype(np.float32), f"{data_file[:-4]}_combined.tif")
