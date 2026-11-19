"""The range maxima along all three axes in one pass (``kernels.mip.range_max_3d_axes``) and the
block MIPs that use it (``tools.mip.block_mips``).

The specification is ``range_max_3d`` along each axis (its oracle cases test MATLAB's
``range_max_3d_mex``): the same bits for every class, NaN of several payloads, -0 and 0, Inf,
ranges of every length, and every split of the rows into blocks and threads (parts of a range
that several threads share are merged in the order of the rows), at every instruction set; an
image takes z as a size of 1; the input errors are range_max_3d's.
"""

from __future__ import annotations

import os
import subprocess
import sys

import numpy as np
import pytest

from petakit5d.core.dtypes import NAN
from petakit5d.kernels import mip
from petakit5d.testing.synth import CounterRNG
from petakit5d.tools import mip as tmip

DTYPES = (np.uint8, np.uint16, np.float32, np.float64)
AXES = (-3, -2, -1)


def _vol(shape, dtype, name):
    """a deterministic (z, y, x) volume; floats with NaN of three payloads, -0 next to 0, Inf,
    -Inf, lines of NaN only and of negatives only"""
    dt = np.dtype(dtype)
    rng = CounterRNG(11, f"mip_axes/{name}/{dt.name}/{'x'.join(map(str, shape))}")
    if dt.kind != "f":
        return rng.integers(0, int(np.iinfo(dt).max), shape, dtype=dt)
    v = ((rng.uniform(shape) - 0.5) * 1000).astype(dt)
    bits = np.uint32 if dt == np.float32 else np.uint64
    hi = 0xFFC00000 if dt == np.float32 else 0xFFF8000000000000
    lo = 0x7FC00000 if dt == np.float32 else 0x7FF8000000000000
    nans = np.array([hi, lo, hi + 0x123], dtype=bits).view(dt)
    flat = v.reshape(-1)
    pick = rng.integers(0, 10, flat.size, dtype=np.int64).reshape(-1)
    for k in range(3):
        flat[pick == k] = nans[k]
    flat[pick == 3] = -0.0
    flat[pick == 4] = 0.0
    flat[pick == 5] = -np.abs(flat[pick == 5]) - 1
    nz, ny, nx = shape
    v[:, ny // 2, :] = nans[0]  # a plane of NaN along y
    v[nz // 2, :, nx // 2] = nans[1]
    v[-1, -1, :] = -0.0
    v[-1, -1, ::3] = 0.0
    v[0, 0, 0] = np.inf
    v[nz - 1, 0, nx - 1] = -np.inf
    return v


def _same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype and a.shape == b.shape, (a.dtype, a.shape, b.dtype, b.shape)
    return np.array_equal(
        np.ascontiguousarray(a).view(np.uint8), np.ascontiguousarray(b).view(np.uint8)
    )


def _three(v, starts):
    """range_max_3d along z, y and x"""
    return tuple(
        mip.range_max_3d(v, ax, s, num_threads=1) for ax, s in zip(AXES, starts, strict=True)
    )


def _starts(n, kind):
    """range starts along an axis of n values"""
    if kind == "one":
        return [0]
    if kind == "each":
        return list(range(n))
    if kind == "fives":
        return list(range(0, n, 5))
    # irregular: ranges of 1, 2, 3, ... values
    s, k, out = 0, 1, []
    while s < n:
        out.append(s)
        s, k = s + k, k + 1
    return out


@pytest.mark.parametrize("dtype", DTYPES)
@pytest.mark.parametrize(
    "shape", [(7, 23, 19), (1, 5, 300), (40, 3, 1), (2, 1, 1), (1, 1, 1), (5, 64, 33), (3, 9, 70)]
)
def test_the_three_maxima_are_range_max_3d_s(dtype, shape):
    v = _vol(shape, dtype, "three")
    for kinds in (("one",) * 3, ("each",) * 3, ("fives",) * 3, ("odd", "fives", "each"),
                  ("each", "odd", "one"), ("odd",) * 3):  # fmt: skip
        starts = [_starts(n, k) for n, k in zip(shape, kinds, strict=True)]
        got = mip.range_max_3d_axes(v, starts)
        for g, w in zip(got, _three(v, starts), strict=True):
            assert _same(g, w), (kinds, shape)


@pytest.mark.parametrize("dtype", DTYPES)
def test_blocks_and_threads_never_change_a_result(dtype):
    """every block size and thread count: thread boundaries inside ranges of the rows (parts
    merged in their order), a range over several threads, ranges of one row"""
    ext = mip._ext()
    v = _vol((6, 41, 37), dtype, "threads")
    for ys in ([0], [0, 3, 4, 20], list(range(41)), [0, 40], _starts(41, "odd")):
        starts = ([0, 2, 5], ys, [0, 1, 17, 36])
        want = _three(v, starts)
        for block_rows in (1, 2, 3, 7, 40, 41, 100):
            for n in (1, 2, 3, 5, 8, 41, 64):
                out = tuple(np.full(w.shape, 7, dtype=w.dtype) for w in want)
                ext.range_max_axes(v, *out, *starts, n, -1, block_rows)
                for o, w in zip(out, want, strict=True):
                    assert _same(o, w), (ys, block_rows, n)


def test_default_blocks_with_threads(monkeypatch):
    """the kernel's own block size (about 128 KB of a plane: a few rows of wide planes) with
    more threads than blocks and thread boundaries inside the ranges of the rows"""
    monkeypatch.setattr(mip, "_threads", lambda n: 1 if n is None else n)
    v = _vol((4, 57, 9000), np.uint16, "wide")
    starts = ([0, 3], [0, 10, 11, 50], [0, 4096, 8192])
    want = _three(v, starts)
    for n in (1, 2, 4, 7, 24):
        for g, w in zip(mip.range_max_3d_axes(v, starts, num_threads=n), want, strict=True):
            assert _same(g, w), n
    v = _vol((3, 31, 2500), np.float64, "wide")
    starts = ([0], [0, 13], [0, 999, 1000])
    want = _three(v, starts)
    for n in (1, 3, 31):
        for g, w in zip(mip.range_max_3d_axes(v, starts, num_threads=n), want, strict=True):
            assert _same(g, w), n


def test_fewer_threads_than_asked_still_do_every_row():
    """an OpenMP runtime that grants fewer threads than asked (OMP_THREAD_LIMIT): every run of
    rows is still done"""
    code = (
        "import numpy as np\n"
        "from petakit5d.kernels import mip\n"
        "v = (np.arange(4 * 41 * 37, dtype=np.uint32) * np.uint32(2654435761) % 65521)\n"
        "v = v.astype(np.uint16).reshape(4, 41, 37)\n"
        "starts = ([0, 1], [0, 3, 20], [0, 9])\n"
        "want = [mip.range_max_3d(v, a, s, num_threads=1) for a, s in zip((0, 1, 2), starts)]\n"
        "out = [np.full(w.shape, 7, w.dtype) for w in want]\n"
        "mip._ext().range_max_axes(v, *out, *starts, 6, -1, 2)\n"
        "assert all(np.array_equal(o, w) for o, w in zip(out, want))\n"
        "print('ok')\n"
    )
    env = {**os.environ, "OMP_THREAD_LIMIT": "2"}
    r = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True)
    assert r.returncode == 0 and r.stdout.strip() == "ok", r.stdout + r.stderr


@pytest.mark.parametrize("dtype", DTYPES)
def test_every_instruction_set_gives_the_same_bits(dtype):
    ext = mip._ext()
    v = _vol((5, 33, 70), dtype, "isa")
    starts = ([0, 2], [0, 5, 6, 30], [0, 1, 32, 64, 69])
    want = _three(v, starts)
    for isa in range(ext.best_isa() + 1):
        out = tuple(np.empty(w.shape, dtype=w.dtype) for w in want)
        ext.range_max_axes(v, *out, *starts, 3, isa, 4)
        for o, w in zip(out, want, strict=True):
            assert _same(o, w), isa


def test_images_take_z_as_one():
    im = _vol((1, 17, 23), np.float32, "image")[0]
    starts = ([0], [0, 4, 9], [0, 22])
    got = mip.range_max_3d_axes(im, starts)
    assert [g.shape for g in got] == [(17, 23), (3, 23), (17, 2)]
    for g, w in zip(got, _three(im, starts), strict=True):
        assert _same(g, w)
    assert _same(got[0], im)


def test_input_errors_are_range_max_3d_s():
    v = np.zeros((3, 4, 5), np.uint16)
    good = ([0], [0], [0])
    cases = [
        (np.zeros((3, 4, 5), np.int16), good, "uint8, uint16, single or double"),
        (np.zeros((2, 3, 4, 5), np.uint16), good, "at most 3 dimensions"),
        (np.zeros((0, 4, 5), np.uint16), good, "must not be empty"),
        (v, ([0], [0]), "3 axes"),
        (v, ([0], [1], [0]), "increasing integers"),
        (v, ([0, 0], [0], [0]), "increasing integers"),
        (v, ([0], [0, 4], [0]), "increasing integers"),
        (v, ([0], [0], [0, 2.5]), "increasing integers"),
        (v, ([0], [], [0]), "nonempty vector"),
    ]
    for vol, starts, msg in cases:
        with pytest.raises(mip.MatlabError, match=msg) as err:
            mip.range_max_3d_axes(vol, starts)
        assert err.value.identifier == "range_max_3d:inputError"


def test_binding_checks_its_outputs():
    ext = mip._ext()
    v = np.zeros((3, 4, 5), np.uint16)
    good = (np.empty((1, 4, 5), np.uint16), np.empty((3, 1, 5), np.uint16),
            np.empty((3, 4, 1), np.uint16))  # fmt: skip
    ext.range_max_axes(v, *good, [0], [0], [0], 1)
    with pytest.raises(ValueError, match="wrong size"):
        ext.range_max_axes(v, good[0], good[0], good[2], [0], [0], [0], 1)
    with pytest.raises(ValueError, match="begin with 0"):
        ext.range_max_axes(v, *good, [1], [0], [0], 1)
    with pytest.raises(ValueError, match="at least 1"):
        ext.range_max_axes(v, *good, [0], [0], [0], 0)


# ---------------------------------------------------------------------------------------------
# tools.mip.block_mips
# ---------------------------------------------------------------------------------------------


def _block_mips_by_axis(block, box, mbs, dtype):
    """MIP_zarr_block_mips's three range_max_3d calls (the code before the one pass)"""
    out = []
    for p in range(3):
        b0, b1, k = box[p], box[p + 3], mbs[p]
        starts = [s - b0 for s in [b0, *range((b0 // k) * k + k, b1, k)]]
        out.append(tmip._cast(tmip._range_max_parts(block, p, starts), dtype))
    return out


@pytest.mark.parametrize(
    "dtype, zarr_dtype", [(np.uint16, "uint16"), (np.float32, "single"), (np.float32, "uint16")]
)
def test_block_mips_are_those_of_each_axis(dtype, zarr_dtype):
    block = _vol((9, 30, 26), dtype, "block")
    for box, mbs in (
        ((0, 0, 0, 9, 30, 26), (9, 30, 26)),
        ((5, 7, 3, 14, 37, 29), (4, 8, 5)),
        ((2, 0, 1, 11, 30, 27), (2, 16, 3)),
    ):
        got = tmip.block_mips(block, box, mbs, zarr_dtype)
        want = _block_mips_by_axis(block, box, mbs, zarr_dtype)
        for g, w in zip(got, want, strict=True):
            assert _same(g, w), (box, mbs)


def test_block_mips_of_other_classes_take_each_axis():
    block = CounterRNG(12, "mip_axes/int16").integers(-500, 500, (5, 12, 10), dtype=np.int16)
    box, mbs = (1, 2, 3, 6, 14, 13), (2, 5, 4)
    got = tmip.block_mips(block, box, mbs, "int16")
    for g, w in zip(got, _block_mips_by_axis(block, box, mbs, "int16"), strict=True):
        assert _same(g, w)


def test_block_mips_take_one_pass(monkeypatch):
    calls = []
    real = mip.range_max_3d_axes

    def spy(v, starts, **kw):
        calls.append([list(s) for s in starts])
        return real(v, starts, **kw)

    monkeypatch.setattr(mip, "range_max_3d_axes", spy)
    block = _vol((6, 20, 18), np.uint16, "one_pass")
    tmip.block_mips(block, (3, 5, 7, 9, 25, 25), (4, 8, 10), "uint16")
    assert calls == [[[0, 1, 5], [0, 3, 11, 19], [0, 3, 13]]]


def test_nan_is_kept_as_the_first_of_a_range():
    """a range of NaN gives its first NaN's payload in all three, at any thread split"""
    v = np.full((4, 6, 5), NAN, np.float32)
    v[1, 2, 3] = np.array([0x7FC00055], np.uint32).view(np.float32)[0]
    v[0, 0, 0] = np.array([0xFFC00011], np.uint32).view(np.float32)[0]
    starts = ([0, 1], [0, 2, 3], [0, 3])
    want = _three(v, starts)
    for n in (1, 2, 6):
        for g, w in zip(mip.range_max_3d_axes(v, starts, num_threads=n), want, strict=True):
            assert _same(g, w), n
