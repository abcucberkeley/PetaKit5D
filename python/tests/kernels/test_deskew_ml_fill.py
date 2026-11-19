"""The ``fill`` of the deskew kernels on MATLAB's layout (``petakit5d._ext.deskew``'s ``warp_ml``
and ``skewed_warp_ml``, the mex's loops): into an output that holds any values they give the
bits they give in a zeroed output (the slices and lines outside the input, the rows of a line
past the input's and a box past all of them written as 0), for every class, scan direction,
step size, box and thread count; without ``fill`` they leave those voxels alone, as the mex.
"""

from __future__ import annotations

import numpy as np
import pytest

from petakit5d._ext import load
from petakit5d.deskew_rotate.geometry import deskew_rotate_geometry
from petakit5d.kernels import deskew

F32 = np.float32
KINDS = ("single", "u16f", "u16")  # single; uint16 to single; uint16 to uint16 (save16bit)
SHAPES = ((37, 45, 29), (20, 17, 64), (12, 33, 100), (2, 9, 12), (1, 9, 12))


def _geometry(shape, dz, reverse):
    try:
        g = deskew_rotate_geometry(
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
    except (ValueError, IndexError):
        return None
    if not g.use_fast_method or (g.do_interp and shape[2] < np.ceil(g.x_step) + 2):
        return None
    return g


def _volume(shape, kind, seed):
    rng = np.random.default_rng(seed)
    if kind == "single":
        return (rng.random(shape, dtype=F32) * 1000).astype(F32)
    return rng.integers(0, 65536, shape, dtype=np.uint16)


def _garbage(shape, dtype, seed):
    rng = np.random.default_rng(seed)
    a = np.empty(shape, dtype=dtype)
    a.view(np.uint8).reshape(-1)[:] = rng.integers(1, 256, a.nbytes, dtype=np.uint8)
    return a


def _boxes(g, ny, seed):
    """the whole output; a box off the output's start; boxes past the input's rows (by 7, by
    one); one past the output in z and x (slices and lines outside the input); one past all the
    input's rows (no row computed); a random one"""
    oz, oy, ox = g.out_size
    rng = np.random.default_rng(seed)
    boxes = [
        (0, 0, 0, oz, oy, ox),
        (1, 3, 2, oz - 1, oy - 1, ox - 2),
        (0, max(0, ny - 5), 0, oz, ny + 7, ox),
        (0, 0, 0, oz, ny + 1, ox),  # (one row past the input's)
        (0, 0, 0, oz + 6, oy, ox + 9),
        (2, ny + 1, 1, oz, ny + 4, ox),
    ]
    z0, z1 = sorted(int(v) for v in rng.integers(0, oz + 1, 2))
    x0, x1 = sorted(int(v) for v in rng.integers(0, ox + 1, 2))
    y0 = int(rng.integers(0, max(1, oy - 1)))
    boxes.append((z0, y0, x0, z1 + 1, int(rng.integers(y0 + 1, oy + 9)), x1 + 1))
    return boxes


def _ml(ext, vm, g, reverse, box, s16, out, nt, **kw):
    """the kernel of g (fused with interpolation, else the warp) on MATLAB's layout into out"""
    start, _, line_len = deskew._box(box, vm.shape[2])
    tmat = deskew._tmat(g.tmat)
    if g.do_interp:
        ext.skewed_warp_ml(
            vm, float(g.x_step), float(g.int_stepsize), reverse, tmat, out, *start, line_len, nt,
            **kw,
        )  # fmt: skip
    else:
        ext.warp_ml(vm, tmat, out, *start, line_len, nt, **kw)
    return out


CASES = [
    (shape, dz, rev, kind)
    for shape in SHAPES
    for dz in (0.2, 0.3, 0.5, 1.0)
    for rev in (True, False)
    for kind in KINDS
]


@pytest.mark.parametrize(
    "shape, dz, reverse, kind",
    CASES,
    ids=[f"{s[0]}x{s[1]}x{s[2]}-dz{d}-{'rev' if r else 'fwd'}-{k}" for s, d, r, k in CASES],
)
def test_fill_gives_the_bits_of_a_zeroed_output(shape, dz, reverse, kind):
    g = _geometry(shape, dz, reverse)
    if g is None:
        pytest.skip("not a kernel geometry")
    ext = load("deskew")
    vm = np.ascontiguousarray(_volume(shape, kind, 11).transpose(0, 2, 1))
    s16 = kind == "u16"
    dt = np.uint16 if vm.dtype == np.uint16 and s16 else F32
    for bi, box in enumerate(_boxes(g, shape[1], 5)):
        _, (oz, oy, ox), line_len = deskew._box(box, shape[1])
        ref = _ml(ext, vm, g, reverse, box, s16, np.zeros((oz, ox, oy), dt), 1)
        for nt in (1, 3):
            got = _ml(ext, vm, g, reverse, box, s16, _garbage(ref.shape, dt, bi), nt, fill=True)
            assert got.tobytes() == ref.tobytes(), (box, nt)
        # (in a zeroed output, fill writes the zeros that are there)
        got = _ml(ext, vm, g, reverse, box, s16, np.zeros(ref.shape, dt), 2, fill=True)
        assert got.tobytes() == ref.tobytes(), box
        if line_len == 0:
            assert not ref.any(), box


@pytest.mark.parametrize("dz", (0.2, 0.3))
def test_fill_writes_only_the_voxels_the_loops_leave_alone(dz):
    """without fill the loops leave the slices and lines outside and the rows past the input's
    as they are (two outputs of other values keep their differences there); fill writes those
    as 0 and nothing else differently"""
    shape = (37, 45, 29)
    g = _geometry(shape, dz, True)
    ext = load("deskew")
    vm = np.ascontiguousarray(_volume(shape, "u16", 2).transpose(0, 2, 1))
    oz, _, ox = g.out_size
    box = (0, 0, 0, oz + 6, shape[1] + 7, ox + 9)
    _, (bz, by, bx), _ = deskew._box(box, shape[1])
    g1 = _garbage((bz, bx, by), np.uint16, 3)
    a = _ml(ext, vm, g, True, box, True, g1.copy(), 2)
    # (every voxel of the second output differs from the first's before the call)
    b = _ml(ext, vm, g, True, box, True, g1 ^ np.uint16(1), 2)
    alone = a != b
    assert alone[:, :, shape[1] :].all(), "rows past the input's were written"
    assert 0 < alone.sum() < alone.size
    filled = _ml(ext, vm, g, True, box, True, _garbage((bz, bx, by), np.uint16, 5), 2, fill=True)
    assert not filled[alone].any()
    assert np.array_equal(filled[~alone], a[~alone])


@pytest.mark.parametrize("interp", (False, True))
@pytest.mark.parametrize("kind", KINDS)
def test_kernel_functions_write_a_kept_output_without_zeroing_it(interp, kind, monkeypatch):
    """the *_matlab_layout functions run the loops with fill into ``out`` (no zeroing pass
    first) and give their new-output bits"""
    shape = (29, 34, 31)
    g = _geometry(shape, 0.3 if interp else 0.2, False)
    vm = np.ascontiguousarray(_volume(shape, kind, 6).transpose(0, 2, 1))
    s16 = kind == "u16"
    box = (1, 2, 0, g.out_size[0] + 3, g.out_size[1] + 4, g.out_size[2] - 1)

    def run(**kw):
        if interp:
            return deskew.skewed_space_interp_volume_deskew_rotate_warp_matlab_layout(
                vm, g.x_step, g.int_stepsize, False, g.tmat, box, s16, num_threads=3, **kw
            )
        return deskew.volume_deskew_rotate_warp_matlab_layout(
            vm, g.tmat, box, s16, num_threads=3, **kw
        )

    ref = run()
    zeroed = []
    monkeypatch.setattr(deskew, "_zero", lambda a, n: zeroed.append(a.shape))
    out = _garbage(ref.shape, ref.dtype, 1)
    got = run(out=out)
    assert got is out
    assert got.tobytes() == ref.tobytes()
    assert zeroed == [], "the kept output was zeroed first"
