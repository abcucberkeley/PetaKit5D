"""The blocked deskew kernels on (z, y, x) volumes (``petakit5d._ext.deskew``'s ``warp`` and
``skewed_warp``): the bits of the mex's loops (the MATLAB-layout kernels, which the oracle cases
of test_deskew.py check against MATLAB) for any rows per block and window of grid slices, with
and without ``fill`` (an output that does not start zeroed: every voxel is written), boxes that
start off a block of rows or end past the input's rows, columns and grid slices running either
way, and any thread count; and the ``out`` of the kernel functions.
"""

from __future__ import annotations

import numpy as np
import pytest

from petakit5d._ext import load
from petakit5d.deskew_rotate.geometry import deskew_rotate_geometry
from petakit5d.kernels import deskew

F32 = np.float32
KINDS = ("single", "u16f", "u16")  # single; uint16 to single; uint16 to uint16 (save16bit)


def _geometry(shape, dz, reverse):
    return deskew_rotate_geometry(
        shape,
        32.45,
        dz,
        0.108,
        reverse=reverse,
        objective_scan=False,
        x_step_thresh=2.0,
        resample_factor=None,
        interp_method="linear",
    )


def _volume(shape, kind, seed):
    rng = np.random.default_rng(seed)
    if kind == "single":
        return (rng.random(shape, dtype=F32) * 1000).astype(F32)
    return rng.integers(0, 65536, shape, dtype=np.uint16)


def _garbage(shape, dtype, seed):
    rng = np.random.default_rng(seed)
    a = np.empty(shape, dtype=dtype)
    a.view(np.uint8).reshape(-1)[:] = rng.integers(0, 256, a.nbytes, dtype=np.uint8)
    return a


def _args(vol, g, box):
    start, out_shape, line_len = deskew._box(box, vol.shape[1])
    return start, out_shape, line_len


def _mex_loops(ext, vol, g, tmat, box, s16, interp):
    """the reference: the mex's loops on MATLAB's layout, the result back in (z, y, x)"""
    start, (oz, oy, ox), line_len = _args(vol, g, box)
    vm = np.ascontiguousarray(vol.transpose(0, 2, 1))
    dt = np.uint16 if vol.dtype == np.uint16 and s16 else F32
    out = np.zeros((oz, ox, oy), dtype=dt)
    if interp:
        ext.skewed_warp_ml(
            vm, float(g.x_step), float(g.int_stepsize), g.reverse, tmat, out, *start, line_len, 2
        )
    else:
        ext.warp_ml(vm, tmat, out, *start, line_len, 2)
    return np.ascontiguousarray(out.transpose(0, 2, 1))


def _blocked(ext, vol, g, tmat, box, s16, interp, out, nt, **kw):
    start, _, line_len = _args(vol, g, box)
    if interp:
        ext.skewed_warp(
            vol,
            float(g.x_step),
            float(g.int_stepsize),
            g.reverse,
            tmat,
            out,
            *start,
            line_len,
            nt,
            **kw,
        )
    else:
        ext.warp(vol, tmat, out, *start, line_len, nt, **kw)
    return out


def _boxes(g, ny, seed):
    """the whole output, a box off the blocks of rows, a box past the input's rows, one inside"""
    oz, oy, ox = g.out_size
    rng = np.random.default_rng(seed)
    boxes = [(0, 0, 0, oz, oy, ox), (1, 3, 2, oz - 1, oy - 1, ox - 2)]
    boxes.append((0, max(0, ny - 5), 0, oz, ny + 7, ox))
    z0, z1 = sorted(int(v) for v in rng.integers(0, oz + 1, 2))
    x0, x1 = sorted(int(v) for v in rng.integers(0, ox + 1, 2))
    y0 = int(rng.integers(0, max(1, oy - 1)))
    boxes.append((z0, y0, x0, z1 + 1, int(rng.integers(y0 + 1, oy + 9)), x1 + 1))
    return boxes


class _G:
    """a geometry with its scan direction (for the calls)"""

    def __init__(self, g, reverse):
        self.__dict__.update({k: getattr(g, k) for k in ("x_step", "int_stepsize", "out_size")})
        self.do_interp = g.do_interp
        self.tmat = g.tmat
        self.reverse = reverse


SHAPES = ((37, 45, 29), (20, 17, 64), (61, 130, 40), (12, 33, 100), (2, 9, 12), (1, 9, 12))
CASES = [
    (shape, dz, rev, kind)
    for shape in SHAPES
    for dz in (0.2, 0.3, 0.5, 1.0)
    for rev in (True, False)
    for kind in KINDS
]


def _cases_ok(shape, dz, rev):
    try:
        g = _geometry(shape, dz, rev)
    except (ValueError, IndexError):
        return None
    if not g.use_fast_method or (g.do_interp and shape[2] < np.ceil(g.x_step) + 2):
        return None
    return _G(g, rev)


@pytest.mark.parametrize(
    "shape, dz, reverse, kind",
    CASES,
    ids=[f"{s[0]}x{s[1]}x{s[2]}-dz{d}-{'rev' if r else 'fwd'}-{k}" for s, d, r, k in CASES],
)
def test_blocked_kernels_equal_the_mex_loops(shape, dz, reverse, kind):
    g = _cases_ok(shape, dz, reverse)
    if g is None:
        pytest.skip("not a kernel geometry")
    ext = load("deskew")
    vol = _volume(shape, kind, 11)
    s16 = kind == "u16"
    tmat = deskew._tmat(g.tmat)
    for bi, box in enumerate(_boxes(g, shape[1], 5)):
        ref = _mex_loops(ext, vol, g, tmat, box, s16, g.do_interp)
        for block_rows, window in ((0, 0), (8, 1), (16, 3), (24, 0), (64, 2), (32, 7)):
            for fill in (False, True):
                out = (
                    _garbage(ref.shape, ref.dtype, bi)
                    if fill
                    else np.zeros(ref.shape, dtype=ref.dtype)
                )
                got = _blocked(
                    ext,
                    vol,
                    g,
                    tmat,
                    box,
                    s16,
                    g.do_interp,
                    out,
                    3 if fill else 1,
                    fill=fill,
                    block_rows=block_rows,
                    window=window,
                )
                assert got.view(np.uint8).tobytes() == ref.view(np.uint8).tobytes(), (
                    box,
                    block_rows,
                    window,
                    fill,
                )


def _mirrored(g, nx, ntotal):
    """the geometry's transform with the columns and the grid slices running the other way"""
    t = np.array(g.tmat, dtype=np.float64)
    t[1, 2], t[1, 3] = -t[1, 2], (nx - 2) - t[1, 3]
    t[2, 1], t[2, 2], t[2, 3] = -t[2, 1], -t[2, 2], (ntotal - 2) - t[2, 3]
    return t


@pytest.mark.parametrize("dz", (0.2, 0.3, 0.8))
@pytest.mark.parametrize("kind", KINDS)
def test_columns_and_grid_slices_running_the_other_way(dz, kind):
    """transforms whose input column falls with the output slice and whose grid slice grows
    along the output lines (the slices are then taken in the other order)"""
    shape = (41, 37, 44)
    g = _cases_ok(shape, dz, True)
    ext = load("deskew")
    vol = _volume(shape, kind, 3)
    s16 = kind == "u16"
    ntotal = int(ext.interp_grid_size(shape[0], float(g.int_stepsize))) if g.do_interp else 41
    for tmat in (_mirrored(g, shape[2], ntotal), deskew._tmat(g.tmat)):
        for box in _boxes(g, shape[1], 9)[:2]:
            ref = _mex_loops(ext, vol, g, tmat, box, s16, g.do_interp)
            for nt in (1, 4):
                out = _garbage(ref.shape, ref.dtype, 2)
                got = _blocked(ext, vol, g, tmat, box, s16, g.do_interp, out, nt, fill=True)
                np.testing.assert_array_equal(got, ref)


@pytest.mark.parametrize("stepsize", (None, 0.9))
@pytest.mark.parametrize("reverse", (True, False))
@pytest.mark.parametrize("kind", ("u16", "single"))
def test_units_starting_at_any_column(reverse, kind, stepsize):
    """windows of one or two grid slices start units at many slices; boxes starting at other
    output slices shift the slices' columns (and the interpolated slices' raw columns, up to xa
    on either side) across the ring's chunks at every offset. A step size of 0.9 puts grid
    slices next to raw ones (z weights 0.1 and 0.9), whose sums read the columns xa away (the
    geometry's own step sizes keep them half way)"""
    shape = (40, 24, 100)
    g = _cases_ok(shape, 0.5, reverse)  # xstep 3.9: xa 4
    if stepsize is not None:
        g.int_stepsize = stepsize
    ext = load("deskew")
    vol = _volume(shape, kind, 12)
    s16 = kind == "u16"
    tmat = deskew._tmat(g.tmat)
    oz, oy, ox = g.out_size
    for z0 in range(34):
        box = (z0, 0, 0, oz, oy, ox)
        ref = _mex_loops(ext, vol, g, tmat, box, s16, True)
        for window in (1, 2):
            out = _garbage(ref.shape, ref.dtype, z0)
            got = _blocked(
                ext, vol, g, tmat, box, s16, True, out, 2, fill=True, block_rows=8, window=window
            )
            assert got.tobytes() == ref.tobytes(), (z0, window)


@pytest.mark.parametrize("kind", KINDS)
def test_threads_and_blocks_do_not_change_results(kind):
    shape = (83, 211, 57)
    g = _cases_ok(shape, 0.3, True)
    ext = load("deskew")
    vol = _volume(shape, kind, 4)
    s16 = kind == "u16"
    tmat = deskew._tmat(g.tmat)
    box = (0, 0, 0, *g.out_size)
    ref = _blocked(ext, vol, g, tmat, box, s16, True, np.zeros(g.out_size, _dt(vol, s16)), 1)
    for nt in (2, 5, 8):
        for block_rows, window in ((0, 0), (16, 9), (48, 40)):
            got = _blocked(
                ext,
                vol,
                g,
                tmat,
                box,
                s16,
                True,
                np.empty(g.out_size, _dt(vol, s16)),
                nt,
                fill=True,
                block_rows=block_rows,
                window=window,
            )
            np.testing.assert_array_equal(got, ref)


def _dt(vol, s16):
    return np.uint16 if vol.dtype == np.uint16 and s16 else F32


@pytest.mark.parametrize("interp", (False, True))
@pytest.mark.parametrize("matlab_layout", (False, True))
def test_kernel_out(interp, matlab_layout):
    """``out``: written (all of it, whatever it held) and returned; a new array without it; the
    wrong class, shape or layout refused"""
    shape = (29, 34, 31)
    g = _cases_ok(shape, 0.3 if interp else 0.2, False)
    vol = _volume(shape, "u16", 6)
    box = (1, 2, 0, g.out_size[0], g.out_size[1] + 4, g.out_size[2] - 1)
    v = np.ascontiguousarray(vol.transpose(0, 2, 1)) if matlab_layout else vol

    def run(**kw):
        if interp:
            f = (
                deskew.skewed_space_interp_volume_deskew_rotate_warp_matlab_layout
                if matlab_layout
                else deskew.skewed_space_interp_volume_deskew_rotate_warp
            )
            return f(v, g.x_step, g.int_stepsize, False, g.tmat, box, True, num_threads=3, **kw)
        f = (
            deskew.volume_deskew_rotate_warp_matlab_layout
            if matlab_layout
            else deskew.volume_deskew_rotate_warp
        )
        return f(v, g.tmat, box, True, num_threads=3, **kw)

    ref = run()
    out = _garbage(ref.shape, ref.dtype, 1)
    got = run(out=out)
    assert got is out
    np.testing.assert_array_equal(got, ref)
    with pytest.raises(ValueError, match="out must be"):
        run(out=np.empty(ref.shape, dtype=F32))
    with pytest.raises(ValueError, match="out must be"):
        run(out=np.empty(ref.shape[::-1], dtype=ref.dtype))
    with pytest.raises(ValueError, match="writeable and C-contiguous"):
        run(out=np.empty(ref.shape[::-1], dtype=ref.dtype).transpose(2, 1, 0))
    with pytest.raises(TypeError, match="NumPy array"):
        run(out=[0])


@pytest.mark.parametrize("layout", ("matlab", "c"))
def test_forced_layouts_write_out(layout):
    """the (z, y, x) functions with a forced kernel layout (the tests' switch) write ``out`` too"""
    shape = (29, 34, 31)
    g = _cases_ok(shape, 0.3, True)
    vol = _volume(shape, "single", 8)
    box = (0, 0, 0, *g.out_size)
    f = deskew.skewed_space_interp_volume_deskew_rotate_warp
    ref = f(vol, g.x_step, g.int_stepsize, True, g.tmat, box, False)
    out = _garbage(ref.shape, ref.dtype, 4)
    got = f(vol, g.x_step, g.int_stepsize, True, g.tmat, box, False, out=out, _kernel_layout=layout)
    assert got is out
    np.testing.assert_array_equal(got, ref)
