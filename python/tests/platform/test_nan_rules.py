"""NaN bits that must not depend on the CPU: x86's default NaN (sign bit set, MATLAB's NaN)
for invalid operations, also on ARM, whose default NaN has the sign bit clear."""

import numpy as np
import pytest

from petakit5d.core import dtypes

_NAN64 = 0xFFF8000000000000
_NAN32 = 0xFFC00000

# (op, a, b) whose result is an invalid operation's NaN: Inf - Inf, 0 * Inf, 0 / 0, Inf / Inf
_INVALID = [
    (dtypes.plus, np.inf, -np.inf),
    (dtypes.minus, np.inf, np.inf),
    (dtypes.times, 0.0, np.inf),
    (dtypes.rdivide, 0.0, 0.0),
    (dtypes.rdivide, -np.inf, np.inf),
    (dtypes.ldivide, 0.0, 0.0),
]


def _bits(x):
    x = np.asarray(x)
    return x.view(np.uint64 if x.dtype == np.float64 else np.uint32)


@pytest.mark.parametrize("op, a, b", _INVALID, ids=lambda v: getattr(v, "__name__", None))
def test_invalid_operation_gives_matlabs_nan_double(op, a, b):
    r = op(np.array([a, 1.0, a]), np.array([b, 2.0, b]))
    assert _bits(r)[0] == _NAN64 and _bits(r)[2] == _NAN64
    assert not np.isnan(r[1])
    assert _bits(op(a, b)) == _NAN64  # 0-d


@pytest.mark.parametrize("op, a, b", _INVALID, ids=lambda v: getattr(v, "__name__", None))
@pytest.mark.parametrize("side", ["array_scalar", "scalar_array", "array_array"])
def test_invalid_operation_gives_matlabs_nan_single(op, a, b, side):
    s = np.array([a, 3.0], dtype=np.float32)
    if side == "array_scalar":
        r = op(s, float(b))
    elif side == "scalar_array":
        r = op(float(b), s)
    else:
        r = op(s, np.array([b, 1.0], dtype=np.float32))
    assert r.dtype == np.float32
    assert np.isnan(r[0]) and _bits(r)[0] == _NAN32, hex(int(_bits(r)[0]))


def test_nan_operand_keeps_its_payload_quieted():
    # the first NaN operand, quieted, whatever the other operand or the CPU
    p = np.array([0x7FF0000000000123, 0xFFF8000000000456], dtype=np.uint64).view(np.float64)
    r = dtypes.plus(p, np.array([np.inf, -np.inf]))
    assert list(_bits(r)) == [0x7FF8000000000123, 0xFFF8000000000456]
    r = dtypes.times(np.array([np.inf, 0.0]), p)
    assert list(_bits(r)) == [0x7FF8000000000123, 0xFFF8000000000456]


@pytest.mark.parametrize("n", [2, 50, 100_000], ids=["sequential", "lanes", "blocks"])
@pytest.mark.parametrize("dt", [np.float64, np.float32])
def test_sum_of_both_infinities_is_matlabs_nan(n, dt):
    from petakit5d.kernels import ops

    v = np.ones(n, dtype=dt)
    v[0], v[-1] = np.inf, -np.inf
    for f in (ops.sum_double, ops.sum_double_ref):
        r = np.float64(f(v))
        assert _bits(r) == _NAN64, (f.__name__, hex(int(_bits(r))))


def test_x86_nan_rewrites_only_the_nan_of_invalid_operations(monkeypatch):
    # the ARM branch, also on x86: a positive default NaN where no operand is NaN becomes
    # MATLAB's NaN, an operand's NaN is kept
    monkeypatch.setattr(dtypes, "_cpu_nan_is_x86", lambda: False)
    arm = np.array([0x7FF8000000000000], dtype=np.uint64).view(np.float64)[0]
    payload = np.array([0x7FF8000000000123], dtype=np.uint64).view(np.float64)[0]
    a = np.array([np.inf, payload, 1.0])
    r = dtypes.x86_nan(np.array([arm, payload, 2.0]), a, np.array([-np.inf, 1.0, 1.0]))
    assert list(_bits(r)) == [_NAN64, 0x7FF8000000000123, _bits(np.float64(2.0))]
    r32 = dtypes.x86_nan(np.array([arm, 1.0], dtype=np.float32), np.float32(0.0))
    assert r32.dtype == np.float32 and _bits(r32)[0] == _NAN32
    assert _bits(dtypes.x86_nan(arm, 0.0, np.inf)) == _NAN64  # 0-d


def test_x86_nan_is_the_identity_on_x86():
    r = np.array([1.0, 2.0])
    if dtypes._cpu_nan_is_x86():
        assert dtypes.x86_nan(r, r) is r


# Two places that made the CPU's NaN of an invalid operation, found on ARM. NumPy's error state
# raises at an invalid operation that is not made MATLAB's NaN, so these fail on x86 too without
# the step (with one thread: the error state is per thread).


def test_distance_map_of_a_negative_weight_degree_has_matlabs_nan():
    # weight degree -1 makes the box weights' floor Inf in single: the voxels of distance 0 are
    # Inf .* 0 (tests/stitch/test_distance.py slab.wd-1)
    from petakit5d.stitch import distance as sd
    from petakit5d.testing.synth import CounterRNG

    slab = CounterRNG(0, "e04/slab").integers(1, 1000, (11, 30, 26), dtype=np.uint16)
    slab[:2] = 0
    slab[:, :, :3] = 0
    with np.errstate(invalid="raise"):
        m = sd.tile_distance_map_mip_slab(
            slab, -1, pool_size=(2, 2, 2), dist_bbox=(2, 4, 6, 18, 50, 40), num_threads=1
        )
    nan = np.isnan(m)
    assert nan.any()
    assert set(m[nan].view(np.uint32).tolist()) == {_NAN32}


def test_slab_shift_of_two_zero_weights_has_matlabs_nan():
    # failed slabs weigh 0: 0 / 0 in the registration's shift, and in its peak where all three
    # are 0 (tests/stitch/test_registration.py weights.nan.3, weights.nan.6)
    from petakit5d.stitch import registration as sr

    rsm = np.arange(1.0, 10.0).reshape(3, 3)
    one_seen = np.array([[0.0], [0.0], [0.5]])  # axis 0 is seen by slabs 0 and 1 only
    with np.errstate(invalid="raise"):
        rel, peak = sr._slab_shift(rsm.copy(), one_seen, [2, 1, 0], np.ones(3))
        _, peak0 = sr._slab_shift(rsm.copy(), np.zeros((3, 1)), [2, 1, 0], np.ones(3))
    assert _bits(rel)[0] == _NAN64 and not np.isnan(rel[1:]).any()
    assert not np.isnan(peak)
    assert _bits(np.float64(peak0)) == _NAN64


# The C++ kernels of the integral image, imwarp and resize. MATLAB gives the first NaN operand
# of each operation as written, else x86's NaN for an invalid operation (measured against the
# mex and imwarp). Which operand's NaN a + b gives is the compiler's choice on every system
# (GCC's on x86-64 took the second at some of these sums), and ARM's default NaN is positive:
# both show in these tests, but for imwarp's NaNs that meet, which x86-64 leaves to the
# compiler.

_P32 = np.array([0x7FC00123], np.uint32).view(np.float32)[0]  # a NaN with a payload
_S32 = np.array([0x7F800456], np.uint32).view(np.float32)[0]  # a signalling one
_M32 = np.float32(dtypes.NAN)  # MATLAB's


def _uint(dt):
    return np.uint32 if np.dtype(dt) == np.float32 else np.uint64


def _ii_reference(a, t):
    """integral_image_3d_mex's three passes (cpp/stitch/integral_image.cpp) in MATLAB's
    arithmetic (core.dtypes): running sums along y, then differences along x and z from the
    end"""
    nz, ny, nx = a.shape
    tz, ty, tx = t
    n = np.zeros((nz + tz - 1, ny + ty - 1, nx + tx - 1), np.float32)
    plus, minus = dtypes.plus, dtypes.minus
    n[:nz, 0, :nx] = a[:, 0]
    for y in range(1, n.shape[1]):
        p = n[:nz, y - 1, :nx]
        if y < ty:
            n[:nz, y, :nx] = plus(p, a[:, y])
        elif y < ny:
            n[:nz, y, :nx] = minus(plus(p, a[:, y]), a[:, y - ty])
        else:
            n[:nz, y, :nx] = minus(p, a[:, y - ty])
    for ax, m, w in ((2, nx, tx), (0, nz, tz)):  # x for the input's planes, then z
        v = np.moveaxis(n[:nz] if ax == 2 else n, ax, 0)  # a view
        post = np.zeros(v.shape[1:], np.float32)
        for i in range(v.shape[0] - 1, 0, -1):
            old = v[i].copy()
            if i == v.shape[0] - 1:
                v[i] = v[v.shape[0] - w]
            elif i >= m:
                v[i] = plus(v[i + 1], v[i - w + 1])
            elif i >= w - 1:
                v[i] = plus(minus(v[i + 1], post), v[i - w + 1])
            else:
                v[i] = minus(v[i + 1], post)
            post = old
    return n


def _ii_inputs():
    """(input, window): two NaNs of other bits meeting at each addition of the three passes, in
    both orders and with a signalling NaN, infinities (Inf - Inf) meeting a NaN, and many of
    both"""
    from petakit5d.testing.synth import CounterRNG

    out = []
    for p, q in ((_P32, _M32), (_M32, _P32), (_S32, _P32)):
        for at, t in (((0, 1), (1, 3, 1)), ((1, 2), (1, 2, 1)), ((2, 3), (1, 2, 1))):
            v = np.arange(1, 7, dtype=np.float32)  # pass_y: running sum, slide (op + fy) - ft
            v[at[0]], v[at[1]] = p, q
            out.append((v.reshape(1, 6, 1), t))
        for n, w in ((3, 3), (4, 2)):  # pass_x and pass_z: x >= nx, and nx > x >= tx - 1
            v = np.arange(1, n + 1, dtype=np.float32)
            v[n - 1], v[n - 2] = p, q
            out.append((v.reshape(1, 1, n), (1, 1, w)))
            # (one row: its NaNs are read only as pass_y's first row)
            out.append(
                (np.ascontiguousarray(np.broadcast_to(v[:, None, None], (n, 1, 9))), (w, 1, 1))
            )
        v = np.ones((3, 2, 5), np.float32)
        v[1, 1, 4], v[1, 1, 1], v[2, 0, 3] = np.inf, p, q
        out.append((v, (2, 2, 3)))
    v = CounterRNG(0, "x16/ii").normal((4, 21, 19)).astype(np.float32)
    v.reshape(-1)[::23] = np.inf
    v.reshape(-1)[5::29] = -np.inf
    out.append((v.copy(), (2, 3, 4)))  # infinities only: invalid operations' NaN
    v.reshape(-1)[11::37] = _P32
    v.reshape(-1)[13::41] = _M32
    out.append((v, (2, 3, 4)))
    return out


_II_INPUTS = _ii_inputs()


@pytest.mark.parametrize("k", range(len(_II_INPUTS)))
def test_integral_image_nans_are_matlabs(k):
    from petakit5d.kernels import stitch as S

    a, t = _II_INPUTS[k]
    got = S.integral_image(a, t, num_threads=1)
    want = _ii_reference(a, t)
    assert np.isnan(want).any()
    bad = np.argwhere(got.view(np.uint32) != want.view(np.uint32))
    assert len(bad) == 0, bad[:4].tolist()


def _imwarp(v, T, method):
    from petakit5d.matlabcompat import warp as W

    return W.imwarp(v, T, method, output_size=v.shape, num_threads=1)


def _without(v, nan):
    """v with the NaNs of nan's bits replaced by 1"""
    u = _uint(v.dtype)
    return np.where(v.view(u) == np.asarray(nan, v.dtype).view(u), v.dtype.type(1), v)


@pytest.mark.parametrize("method", ["cubic", "linear"])
@pytest.mark.parametrize("dt", [np.float32, np.float64])
@pytest.mark.parametrize("axis", [0, 1, 2])
@pytest.mark.parametrize("where", ["inside", "start"])
def test_imwarp_meeting_nans_take_the_first_as_written(method, dt, axis, where):
    if dtypes._cpu_nan_is_x86():
        pytest.skip(
            "x86-64's build of imwarp keeps the compiler's NaN where NaNs meet (GCC: the second "
            "operand's at some sums); ARM's computes MATLAB's (docs/platforms.md)"
        )
    # a NaN with a payload, p, and MATLAB's, q, side by side along an axis (or at its start two
    # apart: cubic's end value 3 (f1 - f2) + f3), read at half voxels along it: p comes first in
    # every product and sum that has both, so where both reach a voxel it has p's NaN
    shape = (7, 8, 9)
    p, q = dt(_P32), dt(_M32)
    T = np.eye(4)
    T[3, 2 - axis] = -0.5  # [x y z]'s step
    v = (np.arange(np.prod(shape)).reshape(shape) % 7 + 1).astype(dt)
    i, j = [3, 4, 4], [3, 4, 4]
    if where == "inside":
        j[axis] += 1
    else:
        i[axis], j[axis] = 0, 2
    v[tuple(i)], v[tuple(j)] = p, q
    out = _imwarp(v, T, method)
    p_only, q_only = _imwarp(_without(v, q), T, method), _imwarp(_without(v, p), T, method)
    want = np.where(np.isnan(p_only), p_only, q_only)
    u = _uint(dt)
    assert (out.view(u) == p.view(u)).any() and (out.view(u) == q.view(u)).any()
    assert np.array_equal(out.view(u), want.view(u))


@pytest.mark.parametrize("method", ["cubic", "linear"])
@pytest.mark.parametrize("dt", [np.float32, np.float64])
def test_imwarp_invalid_operations_give_matlabs_nan(method, dt):
    # infinities and no NaN, at whole coordinates (Inf * 0: linear's second sample, cubic's zero
    # weights) and between voxels (Inf - Inf of neighbours): every NaN is x86's
    from petakit5d.testing.synth import CounterRNG

    v = (CounterRNG(0, "x16/w").normal((7, 8, 9)) * 100).astype(dt)
    v.reshape(-1)[::17] = np.inf
    v.reshape(-1)[3::19] = -np.inf
    T = np.eye(4)
    for shift in ((0.0, 0.0, 0.0), (-0.5, -0.5, 0.25)):
        T[3, :3] = shift
        out = _imwarp(v, T, method)
        nan = np.isnan(out)
        assert nan.any()
        assert set(out[nan].view(_uint(dt)).tolist()) == {_NAN32 if dt == np.float32 else _NAN64}


def test_resize_kernels_invalid_operations_give_matlabs_nan():
    # fdlibm's pow: a negative base with a fractional exponent (0 / 0), (+-1)^(+-Inf) (Inf - Inf);
    # the cubic kernel where |x|^3 overflows (Inf * 0) and of +-Inf (Inf - Inf)
    from petakit5d._ext import load
    from petakit5d.matlabcompat import resize as R

    x = np.array([-1.0, -2.5, -0.5])
    for y in (0.5, -2.5, 1e-3):
        o = np.empty_like(x)
        load("resize").fdlibm_pow(x, y, o)
        assert list(_bits(o)) == [_NAN64] * 3
    x = np.array([1.0, -1.0])
    for y in (np.inf, -np.inf):
        o = np.empty_like(x)
        load("resize").fdlibm_pow(x, y, o)
        assert list(_bits(o)) == [_NAN64] * 2
    c = R.cubic(np.array([np.inf, -np.inf, 1e120, -1e200, 1.5]))
    assert list(_bits(c[:4])) == [_NAN64] * 4 and c[4] == R.cubic(np.array([1.5]))[0]
