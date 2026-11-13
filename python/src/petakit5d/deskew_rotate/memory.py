"""The peak memory of a deskew/rotate task of the port (``XR_estimateDeskewRotateMemory``,
re-measured for Python): one ``deskew_rotate_file`` call, or one batch of the large-file block
task.

``memory = BASE_GB + MARGIN * arrays``, ``arrays`` the largest set of arrays the task's code
path holds at the same time (computed per path from the sizes and options below), ``BASE_GB``
the Python process itself (the interpreter, NumPy, SciPy, cpp-tiff, cpp-zarr and the port's
modules) with its libraries' buffers, ``MARGIN`` for the readers' and writers' buffers per
thread. The output size is MATLAB's (``deskewRotateFrame3D``'s, with ``XR_deskewRotateFrame``'s
binning of more than 1000 slices before a combined deskew/rotation with resampling), compared
with MATLAB's in oracle cases.

The arrays per path, from the port's code (``deskew_rotate.file``), where ``b`` is the input's
bytes per voxel and ``N`` a voxel count:

- read: ``b N_read`` (partial files are read into one array, without MATLAB's concatenation
  copy; with an input box ``b (N_read + N_in)`` for the crop);
- flip, flat field and input axis order (``preprocess``): ``2 b N_in`` (a new array each);
- combined, deskew kernels (sample scan, linear, no resampling): ``b N_in + o N_out``, the
  input and output turned to MATLAB's layout in place (one page of scratch per thread) and the
  skewed space interpolation done in per-thread line caches (no interpolated volume, which
  MATLAB's mex holds), ``o`` = 2 for uint16 data with ``save16bit``, else 4;
- combined, ``imwarp`` (objective scan, resampling, cubic): ``b N_in + o (N_int + N_out)``
  with the interpolated volume ``N_int`` where the x step is above the threshold (the slices
  with the raw ones), ``imwarp`` keeping its input's class and computing each voxel in place;
- separate deskew and rotation (``imwarp`` twice, in single): ``(b + 4) N_in`` for the
  single copy, ``4 (N_in + N_ds)``, ``6 N_ds`` for the 16-bit DS, ``4 (N_ds + N_out)`` and
  ``6 N_out`` (above 1000 slices the deskewed volume is read back for the rotation: ``6 N_ds``);
- the 16-bit or single copy of a result, the output axis order (``2 o N_out``), the
  rescale (``rescale_rotate``: the result, its double slab copies and the single or double
  result), and zarr writes (the result and a chunk per writer thread, measured).

The constants were measured with the port's ``deskew_rotate_file`` on the 512 x 1800 demo
tile, cropped and repeated along z (one process per task, the peak resident set of the process
tree, 1 and 24 threads): see the module's tests and the task report.

Sizes are ``(nz, ny, nx)`` and ``block_size`` ``(z, y, x)``; ``resample_factor`` ``(z, y, x)``
(or one value, or ``(z, xy)``).
"""

from __future__ import annotations

import math
from collections.abc import Sequence

__all__ = ["BASE_GB", "MARGIN", "estimate_deskew_rotate_memory"]

# the Python process with the port's modules imported and its libraries' buffers (GB)
BASE_GB = 1.0
# the margin on the arrays: thread buffers of the readers, kernels and writers
MARGIN = 1.1
# zarr writes hold about a chunk per thread (measured: 4.4 GB for 16 single chunks of 270 MB at
# 24 threads), counted for up to this many threads (MATLAB's rule)
ZARR_THREADS = 48
_GB = float(2**30)


def _byte_number(data_type) -> int:
    import numpy as np

    from petakit5d.core.arrays import data_type_to_byte_number, matlab_class

    name = data_type if isinstance(data_type, str) else matlab_class(np.dtype(data_type))
    return data_type_to_byte_number(name)


def _resample_xyz(resample_factor) -> list[float] | None:
    """MATLAB's completed ``[x y z]`` factor, None for none or all ones"""
    if resample_factor is None:
        return None
    v = [
        float(x)
        for x in (resample_factor if isinstance(resample_factor, Sequence) else [resample_factor])
    ]
    if not v:
        return None
    if len(v) == 1:
        rs = [v[0], v[0], v[0]]
    elif len(v) == 2:
        rs = [v[1], v[1], v[0]]
    else:
        rs = [v[2], v[1], v[0]]
    return None if all(r == 1 for r in rs) else rs


def estimate_deskew_rotate_memory(
    im_size: Sequence[int],
    xy_pixel_size: float,
    dz: float,
    *,
    skew_angle: float = 32.45,
    objective_scan: bool = False,
    z_stage_scan: bool = False,
    dsr_combined: bool = True,
    rotate: bool = True,
    save_16bit: bool = True,
    interp_method: str = "linear",
    resample_factor: Sequence[float] | None = None,
    read_size: Sequence[int] | None = None,
    partial_file: bool = False,
    preprocess: bool = False,
    out_size: Sequence[int] | None = None,
    output_axis_order: str = "zyx",
    save_zarr: bool = False,
    block_size: Sequence[int] = (500, 500, 500),
    large_file: bool = False,
    data_type: str = "uint16",
    rescale_rotate: bool = False,
    x_step_thresh: float = 2.0,
) -> tuple[float, float, tuple[int, int, int]]:
    """The peak memory (GB, 2^30 bytes) of a deskew/rotate task of the port, the arrays it
    holds at the peak (GB) and the output size ``(nz, ny, nx)``.

    ``im_size``: the volume deskewed/rotated (all partial files together, after the input
    box; one input batch for ``large_file``), ``(nz, ny, nx)``. ``read_size``: what is read
    before the input box (None: ``im_size``). ``partial_file``: read from several files.
    ``preprocess``: flat field, flip or input axis order. ``out_size``: the output size
    (one batch's in large-file mode; None: computed). ``data_type``: the input's class (MATLAB's
    name or a NumPy type). The other options are ``deskew_rotate_file``'s. MATLAB's options and
    defaults, with ``rescale_rotate`` and ``x_step_thresh`` added, except ``block_size``:
    ``deskew_rotate_file``'s chunks ``(500, 500, 500)`` (MATLAB's estimate assumes 256^3; its
    callers pass the chunks they write).
    """
    from petakit5d.deskew_rotate.geometry import _ceil, _cos, _round, _sin, interp_step
    from petakit5d.matlabcompat.elementary import sind

    del partial_file  # partial files are read into one array
    nz, ny, nx = (int(v) for v in im_size)
    n_in = nz * ny * nx
    rd = im_size if read_size is None or len(read_size) == 0 else read_size
    n_read = math.prod(int(v) for v in rd)
    b = _byte_number(data_type)
    linear = interp_method.lower() == "linear"
    angle = 90 - skew_angle if z_stage_scan else skew_angle
    theta = angle * math.pi / 180
    rs = _resample_xyz(resample_factor)

    arrays = [b * n_read]
    if n_read > n_in:
        arrays.append(b * (n_read + n_in))
    if preprocess:
        arrays.append(2 * b * n_in)

    # XR_deskewRotateFrame bins (z) and resizes (xy) inputs of over 1000 slices before a
    # combined deskew/rotation with resampling
    if dsr_combined and rotate and not large_file and rs is not None and nz > 1000:
        import numpy as np

        out_px = [r * xy_pixel_size for r in rs]
        m = min(out_px)
        sd = float(np.asarray(sind(np.array([angle], dtype=np.float64))).reshape(-1)[0])
        pre = [m / xy_pixel_size, m / xy_pixel_size, m / (dz * sd)]
        pre[2] = max(1.0, _round(pre[2]))
        if any(p != 1 for p in pre):
            nz_b = len(range(int(_round(pre[2] / 2)) - 1, nz, int(pre[2])))
            n_bin = ny * nx * nz_b
            ny_rs, nx_rs = int(_round(ny / pre[0])), int(_round(nx / pre[1]))
            arrays.append(b * (n_in + n_bin))
            arrays.append(b * (n_bin + ny_rs * nx_rs * nz_b))
            ny, nx, nz = ny_rs, nx_rs, nz_b
            n_in = ny * nx * nz
            xy_pixel_size = xy_pixel_size * pre[0]
            dz = dz * pre[2]
            rs = [p / m for p in out_px]
            rs = None if all(r == 1 for r in rs) else rs

    dx = _cos(theta) * dz / xy_pixel_size
    if out_size is None or len(out_size) == 0:
        st = _sin(abs(theta))
        if not objective_scan:
            z_aniso = st * dz / xy_pixel_size
            o = [ny, (nx - 1) * _cos(theta) + (nz - 1) * z_aniso / st, (nx - 1) * st - 4]
        else:
            z_aniso = dz / xy_pixel_size
            o = [ny, nx * _cos(theta) + nz * z_aniso * st, nz * z_aniso * _cos(theta) + nx * st]
        o = [_round(v) for v in o]
        if rs is not None:
            o = [_round(o[0] / rs[1]), _round(o[1] / rs[0]), _round(o[2] / rs[2])]
        osz = (int(o[2]), int(o[0]), int(o[1]))
    else:
        osz = tuple(int(v) for v in out_size)
    n_out = math.prod(max(v, 0) for v in osz)

    if large_file and not save_16bit:
        arrays.append((b + 4) * n_in)
        b = 4
    if not save_16bit and b != 2 and b != 4 and not large_file:
        arrays.append((b + 4) * n_in)
        b = 4

    if not dsr_combined or not rotate:
        nx_ds = nx + _ceil((nz - 1) * abs(dx))
        n_ds = int(ny * nx_ds * nz)
        if b != 4:
            arrays.append((b + 4) * n_in)
        arrays.append(4 * (n_in + n_ds))
        arrays.append((4 + 2 * save_16bit) * n_ds)
        if nz > 1000 and rotate:
            arrays.append((2 * save_16bit + 4) * n_ds)
        out_b = 4
        if rotate:
            arrays.append(4 * (n_ds + n_out))
            if save_16bit:
                arrays.append(6 * n_out)
                out_b = 2
        else:
            n_out = 0
    else:
        use_kernels = not objective_scan and rs is None and linear and b in (2, 4)
        out_b = 2 if (save_16bit and b == 2) else 4
        if use_kernels:
            # the in-place page transposes' scratch: a page per thread up to 1/16 of the array
            turn = max(2**26, max(b * n_in, out_b * n_out) / 16)
            arrays.append(b * n_in + out_b * n_out + turn)
        else:
            n_int = 0
            if not objective_scan and abs(dx) > x_step_thresh:
                _, step = interp_step(dz, dx, theta, xy_pixel_size, x_step_thresh)
                n_int = ny * nx * (int(math.floor(_round((nz - 1) / step * 100000) / 100000)) + 1)
            else:
                out_b = b
            arrays.append(b * n_in + out_b * (n_int + n_out))
        if save_16bit and out_b != 2:
            arrays.append((out_b + 2) * n_out)
        elif not save_16bit and out_b != 4:
            arrays.append((out_b + 4) * n_out)
    if output_axis_order.lower() != "zyx" and n_out:
        arrays.append(2 * out_b * n_out)
    final_b = 2 if save_16bit else 4
    if rescale_rotate and n_out:
        # the result and the rescaled array (of the result's class, single or double; double
        # temporaries per slab only)
        arrays.append((out_b + max(out_b, 4)) * n_out)
    if (save_zarr or large_file) and n_out:
        chunk = [min(int(c), max(int(s), 1)) for c, s in zip(block_size, osz, strict=True)]
        nchunks = math.prod(math.ceil(max(s, 1) / c) for s, c in zip(osz, chunk, strict=True))
        arrays.append(final_b * n_out + min(ZARR_THREADS, nchunks) * math.prod(chunk) * final_b)
    array_gb = max(arrays) / _GB
    return BASE_GB + MARGIN * array_gb, array_gb, osz
