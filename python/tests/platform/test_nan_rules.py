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
