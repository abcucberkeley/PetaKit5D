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


# The linear filters (C++ filters) and the deskew kernels (C++ deskew). x86 gives MATLAB's NaN for
# Inf * 0 and Inf - Inf, and its fused multiply-adds take the factors' NaN before the addend's
# (the running sum's) in every form; ARM's default NaN is positive, and its fused multiply-add
# takes the addend's NaN first and gives its default NaN for 0 * Inf + NaN. There a value that
# comes out NaN is computed again with x86's NaN for each operation: the filters on ARM (a tap's
# or a chain's NaN before the running sum's, as MATLAB gives them: measured), the deskew kernels
# in every build without the mex's AVX2 code (deskew_scalar on x86-64 too; the mex's formula's
# first NaN operand as written).

_X86_NAN = {np.dtype(np.float32): _NAN32, np.dtype(np.float64): _NAN64}
_P64 = np.array([0x7FF8000000000123], np.uint64).view(np.float64)[0]
_S64 = np.array([0x7FF0000000000456], np.uint64).view(np.float64)[0]


def _infinities(name, shape, dt, frac=0.1):
    """normal values with infinities (no NaN) at a fraction frac of the elements, half of them
    negative: every NaN of an output is an invalid operation's"""
    from petakit5d.testing.synth import CounterRNG

    a = (CounterRNG(0, f"x20/{name}").normal(shape) * 100).astype(dt)
    u = CounterRNG(0, f"x20/{name}/inf").uniform(shape)
    a[u < frac / 2] = np.inf
    a[(u >= frac / 2) & (u < frac)] = -np.inf
    return np.ascontiguousarray(a)


def _all_x86_nans(outs):
    total = 0
    for k, o in enumerate(outs):
        o = np.asarray(o)
        nan = np.isnan(o)
        total += int(nan.sum())
        bits = {hex(b) for b in o[nan].view(_uint(o.dtype)).tolist()}
        assert bits <= {hex(_X86_NAN[o.dtype])}, (k, bits)
    assert total > 0


@pytest.mark.parametrize("dt", [np.float32, np.float64])
def test_filters_invalid_operations_give_matlabs_nan(dt):
    # every path of the filters' kernels: tap correlations fused (convn, conv2, imfilter's general
    # kernels) and plain (a kernel larger than the image, imfilter's sparse path), and the tap
    # programs along x (the padded line, edge taps), y (blocks of rows; padded columns' chains
    # added) and z, with each padding, and imgaussfilt
    from petakit5d.matlabcompat import filters as F
    from petakit5d.testing.synth import CounterRNG

    v = _infinities("fl", (5, 23, 31), dt)
    k = np.abs(CounterRNG(0, "x20/fk").normal((3, 5, 4)))
    k[1, 2, 1] = 0.0  # a zero tap: Inf * 0
    col = np.array([[0.5], [0.0], [2.0], [-1.0], [0.25]])
    outs = [F.convn(v, k.astype(dt), "same"), F.convn(v, k.astype(dt), "full")]
    outs += [F.conv2(v[0], np.array([[1.0, 0.0, -1.0]]).astype(dt), "same")]
    outs += [F.conv2(v[0, :2, :3], k[1].astype(dt))]  # the roles swapped: plain sums
    outs += [F.imfilter(v, k, "symmetric"), F.imfilter(v, k[1], "replicate", "full")]
    for pad in ("replicate", "symmetric", "circular", 3.0):
        outs += [F.imfilter(v, h, pad) for h in (col, col.T, col[:, :, None])]
    outs += [F.imfilter(v, col, "full"), F.imfilter(v, col.T, 0.0)]
    sparse = np.zeros((4, 40, 33), dt)
    sparse.reshape(-1)[::23] = np.inf
    sparse.reshape(-1)[5::29] = -np.inf
    outs += [F.imfilter(sparse, h) for h in (col, col.T, k[1, :3, :3])]
    outs += [F.imgaussfilt3(v, 1.5), F.imgaussfilt(v[1], 1.2), F.imgaussfilt(v[2], 0.6)]
    _all_x86_nans(outs)


def _quiet(x):
    x = np.asarray(x)
    u = _uint(x.dtype)
    q = 0x00400000 if u == np.uint32 else 0x0008000000000000
    return np.array([int(x.view(u)) | q], u).view(x.dtype)[0]


def _tap_ref(w, s, acc, fused):
    """a tap with x86's NaN, a tap's NaN before the running sum's (MATLAB's, measured): fused,
    the factors' NaN, then the sum's (quieted), else MATLAB's NaN of an invalid operation (the
    tests' values add exactly); plain, MATLAB's arithmetic: w * s + acc"""
    if not fused:
        return type(acc)(dtypes.plus(dtypes.times(w, s), acc))
    for x in (w, s, acc):
        if np.isnan(x):
            return _quiet(x)
    r = float(w) * float(s) + float(acc)
    return type(acc)(dtypes.NAN if r != r else r)


def _pad_index(q, n, mode):
    if 0 <= q < n:
        return q
    if mode == 1:
        return 0 if q < 0 else n - 1
    if mode == 2:
        r = q % (2 * n)
        return r if r < n else 2 * n - 1 - r
    if mode == 3:
        return q % n
    return -1


def _program_ref(src, w, prog, rowmap, pre, n_out, mode, value):
    """cpp/filters' line_program with x86's NaN: per output, its program's chains (fused or plain
    taps; -1 adds a chain to the sum, its NaN first, -3 continues from the sum; edge taps fused)"""
    outer, n_in, inner = src.shape
    dt = src.dtype.type
    lo, hi = 1 << 21, 1 << 22
    out = np.empty((outer, n_out, inner), src.dtype)
    for o, i, x in np.ndindex(out.shape):
        acc, total, first, opened = dt(0), dt(0), True, False
        for op in [int(c) for c in prog[rowmap[i]]] + [-2]:
            if op in (-1, -2, -3):
                if opened:
                    total = acc if first else dt(dtypes.plus(acc, total))
                    first, acc, opened = False, dt(0), False
                if op == -2:
                    break
                if op == -3 and not first:
                    acc, first, opened = total, True, True
                continue
            fused = op >= 0
            if op <= -hi:
                t, q, fused = -op - hi, n_in, True
            elif op <= -lo:
                t, q, fused = -op - lo, -1, True
            else:
                t = op if fused else -op - 10
                q = i + t - pre
            r = _pad_index(q, n_in, mode)
            acc = _tap_ref(w[t], dt(value) if r < 0 else src[o, r, x], acc, fused)
            opened = True
        out[o, i, x] = dt(0) if first else total
    return out


def _meeting_values(name, shape, dt):
    """small integers (the sums exact) with NaNs of three kinds (MATLAB's, a payload, a
    signalling one), infinities and zeros, so that they meet in every combination"""
    from petakit5d.testing.synth import CounterRNG

    a = CounterRNG(0, f"x20/{name}").integers(-3, 3, shape).astype(dt)
    u = CounterRNG(0, f"x20/{name}/s").uniform(shape)
    one = dt == np.float32
    special = [_M32 if one else dtypes.NAN, _P32 if one else _P64, np.inf, -np.inf]
    special.append(_S32 if one else _S64)
    for j, sv in enumerate(special):
        a[(u >= 0.06 * j) & (u < 0.06 * (j + 1))] = sv
    return np.ascontiguousarray(a)


_W = [0.5, 0.0, 2.0, -1.0]  # a zero tap (Inf * 0) and a negative one (Inf - Inf of one sign)
_EDGE_LOW, _EDGE_HIGH = 1 << 21, 1 << 22
_FUSED = [[3, 2, 1, 0], [0, 1, 2, 3], [2, 0, 3]]
_EDGE = [[1, 2, -(_EDGE_LOW + 4)], [0, 1, -(_EDGE_HIGH + 5)], [-(_EDGE_LOW + 4), 3]]
_MIXED = [[3, 2, -1, 1, 0], [-13, 2, -1, -11, -3, 0], [-10, -11, -12, -13], [1, -3, -12, 0]]


def _run_program(src, programs, pre, mode, value, extra=()):
    """line_program on src (outer, n, inner), runs of three rows per program, and its reference"""
    from petakit5d.matlabcompat import filters as F

    dt = src.dtype.type
    w = np.array(_W + list(extra), dt)
    outer, n_in, inner = src.shape
    rows = [programs[(i // 3) % len(programs)] for i in range(n_in)]
    prog, rowmap = F._programs(rows)
    out = np.empty_like(src)
    fn = F._ext().line_program_f64 if dt == np.float64 else F._ext().line_program_f32
    fn(src, out, outer, n_in, inner, n_in, w, prog, rowmap, pre, len(_W), mode, value, 2)
    return out, _program_ref(src, w, prog, rowmap, pre, n_in, mode, value)


def _fma_rule_holds_here() -> bool:
    """x86-64 computes the fused chains with its FMA instruction where the filters dispatch to it
    (GCC, Clang); the baseline (MSVC's build) calls the C library's fma"""
    from petakit5d.matlabcompat import filters as F

    return not dtypes._cpu_nan_is_x86() or F._ext().best_isa() >= 1


def _assert_same_bits(got, want):
    u = _uint(got.dtype)
    bad = np.argwhere(got.view(u) != want.view(u))
    assert len(bad) == 0, [
        (tuple(b), hex(int(got[tuple(b)].view(u))), hex(int(want[tuple(b)].view(u))))
        for b in bad[:4]
    ]


@pytest.mark.parametrize("inner", [1, 7])
@pytest.mark.parametrize("dt", [np.float32, np.float64])
def test_filters_fused_chains_take_a_taps_nan_first(dt, inner):
    # fused chains and edge taps: x86's fused multiply-add takes a tap's NaN before the running
    # sum's (a later tap's NaN wins), and keeps the sum's NaN at 0 * Inf
    if not _fma_rule_holds_here():
        pytest.skip(
            "the filters' baseline calls the C library's fma, whose NaN where NaNs meet is the "
            "library's (x86-64 builds without the FMA dispatch: MSVC)"
        )
    src = _meeting_values(f"ff{inner}", (3, 19, inner), dt)
    for programs, pre, mode, value, extra in (
        (_FUSED, 2, 0, 0.0, ()),
        (_FUSED, 1, 2, 0.0, ()),
        (_EDGE, 1, 1, 0.0, (1.5, -0.5)),
        (_EDGE, 2, 0, np.inf, (1.5, -0.5)),
    ):
        got, want = _run_program(src, programs, pre, mode, value, extra)
        assert np.isnan(want).sum() > 10
        _assert_same_bits(got, want)
    # a payload NaN in the sum, then 0 * Inf: the sum's NaN (ARM's instruction gives its default)
    p = _P32 if dt == np.float32 else _P64
    line = np.array([p, np.inf], dt).reshape(1, 2, 1)
    got, want = _run_program(line, [[0, 1]], 0, 0, 0.0)
    assert _bits(got[0, 0, 0]) == _bits(_quiet(p)) == _bits(want[0, 0, 0])


@pytest.mark.parametrize("inner", [1, 7])
@pytest.mark.parametrize("dt", [np.float32, np.float64])
def test_filters_plain_sums_take_a_taps_nan_first(dt, inner):
    # plain taps, chains added (-1) and continued from the sum (-3): a tap's product and a
    # chain's NaN before the running sum's, as MATLAB's plain path and padded columns give them
    # (measured)
    if dtypes._cpu_nan_is_x86():
        pytest.skip(
            "x86-64's filters keep the compiler's operand order in the plain sums; ARM computes "
            "MATLAB's (docs/platforms.md)"
        )
    src = _meeting_values(f"fp{inner}", (3, 19, inner), dt)
    for pre, mode in ((2, 0), (1, 1), (3, 3)):
        got, want = _run_program(src, _MIXED, pre, mode, 0.0)
        assert np.isnan(want).sum() > 10
        _assert_same_bits(got, want)


@pytest.mark.parametrize("fused", [True, False])
@pytest.mark.parametrize("dt", [np.float32, np.float64])
def test_filters_tap_correlation_nans(dt, fused):
    # convn's and imfilter's kernel (one chain per output over the taps): a tap's NaN before
    # the running sum's, fused (x86's instruction) and plain (MATLAB's sparse path, measured)
    from petakit5d.matlabcompat import filters as F

    if fused and not _fma_rule_holds_here():
        pytest.skip("the filters' baseline calls the C library's fma (MSVC builds)")
    if not fused and dtypes._cpu_nan_is_x86():
        pytest.skip("x86-64's filters keep the compiler's operand order in the plain sums")
    P = _meeting_values("tc", (4, 9, 37), dt)
    offs = np.array([[1, 0, 2], [0, 2, 0], [1, 1, 1], [0, 0, 3], [1, 2, 0]], np.int64)
    w = np.array([0.5, 0.0, 2.0, -1.0, 0.25], dt)
    out = np.empty((3, 7, 34), dt)
    fn = F._ext().taps_nd_f64 if dt == np.float64 else F._ext().taps_nd_f32
    fn(P, out, offs, w, fused, 2)
    want = np.empty_like(out)
    for z, y, x in np.ndindex(out.shape):
        acc = dt(0)
        for k in range(len(w)):
            acc = _tap_ref(w[k], P[z + offs[k, 0], y + offs[k, 1], x + offs[k, 2]], acc, fused)
        want[z, y, x] = acc
    assert np.isnan(want).sum() > 50
    _assert_same_bits(out, want)


def _deskew_builds():
    from petakit5d._ext import available

    return ["deskew"] + (["deskew_scalar"] if available("deskew_scalar") else [])


def _deskew_outputs(module, vol, dz, reverse, steps=None):
    """the float kernels of a build of the deskew module on vol (z, y, x), on both layouts: the
    warp, or the fused skewed space interpolation and warp and the interpolation to a step size
    (steps: (xstep, stepsize) instead of the geometry's)"""
    from petakit5d.deskew_rotate.geometry import deskew_rotate_geometry
    from petakit5d.kernels import deskew as DS

    g = deskew_rotate_geometry(
        vol.shape, 32.45, dz, 0.108, reverse=reverse, objective_scan=False, x_step_thresh=2.0,
        resample_factor=None, interp_method="linear",
    )  # fmt: skip
    box = (0, 0, 0, *g.out_size)
    zxy = np.ascontiguousarray(vol.transpose(0, 2, 1))
    kw = {"num_threads": 2}
    real = DS._ext
    DS._ext = lambda m="deskew": real(module if m == "deskew" else m)
    try:
        if not g.do_interp:
            return [
                DS.volume_deskew_rotate_warp(vol, g.tmat, box, False, **kw),
                DS.volume_deskew_rotate_warp_matlab_layout(zxy, g.tmat, box, False, **kw),
            ]
        xs, st = steps or (g.x_step, g.int_stepsize)
        return [
            DS.skewed_space_interp_volume_deskew_rotate_warp(
                vol, xs, st, reverse, g.tmat, box, False, **kw
            ),
            DS.skewed_space_interp_volume_deskew_rotate_warp_matlab_layout(
                zxy, xs, st, reverse, g.tmat, box, False, **kw
            ),
            DS.skewed_space_interp_defined_stepsize(vol, xs, st, reverse, False, **kw),
            DS.skewed_space_interp_defined_stepsize_matlab_layout(
                zxy, xs, st, reverse, False, **kw
            ),
        ]
    finally:
        DS._ext = real


_DESKEW_CASES = [(dz, rev) for dz in (0.5, 0.2) for rev in (False, True)]  # interpolated, or not


@pytest.mark.parametrize("module", _deskew_builds())
@pytest.mark.parametrize(
    "dz, reverse, steps",
    [(dz, rev, None) for dz, rev in _DESKEW_CASES]
    + [(0.5, rev, (2.5, 0.2)) for rev in (False, True)],
)
def test_deskew_invalid_operations_give_matlabs_nan(module, dz, reverse, steps):
    # single data with infinities through the warp (Inf * 0 at whole coordinates, Inf - Inf),
    # the fused kernel's interpolated grid columns and lines, and the interpolation to a step
    # size, on both layouts: every NaN is x86's, the layouts' kernels (blocked, and the mex's
    # loops) give the same bits, and so do the builds (without AVX2, and the AVX2 one). xstep 2.5
    # and step 0.2 shift some interpolated slices by whole pixels: zero weights (Inf * 0).
    vol = _infinities(f"ds{dz}{reverse}", (21, 18, 67), np.float32)
    outs = _deskew_outputs(module, vol, dz, reverse, steps)
    _all_x86_nans(outs)
    for zyx, zxy in zip(outs[0::2], outs[1::2], strict=True):
        _assert_same_bits(zyx, np.ascontiguousarray(zxy.transpose(0, 2, 1)))
    if module != "deskew":
        for a, b in zip(outs, _deskew_outputs("deskew", vol, dz, reverse, steps), strict=True):
            _assert_same_bits(a, b)


@pytest.mark.parametrize("module", _deskew_builds())
@pytest.mark.parametrize("dz, reverse", _DESKEW_CASES)
def test_deskew_meeting_nans_take_the_first_as_written(module, dz, reverse):
    from petakit5d._ext import load

    if load(module).simd() != "scalar":
        pytest.skip(
            "the AVX2 build keeps the mex's code: where NaNs meet, the compiler's choice of "
            "operand (docs/platforms.md)"
        )
    # a payload NaN p and MATLAB's q placed so that p comes first in every sum that has both:
    # q at p's column of the next slice (the warp's 4-term sums take slice z before z + 1 at a
    # column; the interpolation's take slice zorg first: q along p's row there), or at the next
    # column of p's slice (each sum takes column x before x + 1: pairs 3 slices and 5 columns
    # apart, so that no sum has two pairs' NaNs). The fused kernel's lines mix raw and
    # interpolated slices: left to the invalid operations' test.
    shape = (21, 18, 67)
    p, q = np.float32(_P32), np.float32(_M32)
    keep = [0, 1] if dz == 0.2 else [2, 3]  # the warp, or the interpolation to a step size
    for place in ("slices", "columns"):
        v = (np.arange(np.prod(shape)).reshape(shape) % 7 + 1).astype(np.float32)
        if place == "columns":
            for z in range(2, 19, 3):
                for x0 in range(1, 65, 7):  # the rows' ends too: the 2-term sums
                    v[z, 5, x0], v[z, 5, x0 + 1] = p, q
        elif dz == 0.2:
            v[9, 5, 30], v[10, 5, 30] = p, q
        else:
            v[4, 11, 12], v[5, 11, :] = p, q
        outs = [
            _deskew_outputs(module, x, dz, reverse) for x in (v, _without(v, q), _without(v, p))
        ]
        for k in keep:
            out, p_only, q_only = outs[0][k], outs[1][k], outs[2][k]
            want = np.where(np.isnan(p_only), p_only, q_only)
            assert (np.isnan(p_only) & np.isnan(q_only)).any(), (place, k)
            assert np.array_equal(out.view(np.uint32), want.view(np.uint32)), (place, k)
