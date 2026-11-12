import numpy as np
import pytest

from petakit5d.core import dtypes

PRED_HALF_S = np.nextafter(np.float32(0.5), np.float32(0))  # 0.49999997
PRED_HALF_D = np.nextafter(0.5, 0.0)  # 0.49999999999999994


def test_eps_is_double_eps():
    assert dtypes.EPS == 2.0**-52 == dtypes.eps()
    assert dtypes.eps("single") == np.float32(2.0**-23)
    assert dtypes.eps(1.0) == 2.0**-52
    assert dtypes.eps(np.float32(1)) == np.float32(2.0**-23)
    assert dtypes.eps(0.0) == 2.0**-1074
    assert dtypes.eps(np.float32(0)) == np.float32(2.0**-149)
    assert dtypes.eps(np.finfo(np.float64).max) == 2.0**971
    assert np.isnan(dtypes.eps(np.inf)) and np.isnan(dtypes.eps(np.nan))
    np.testing.assert_array_equal(
        dtypes.eps(np.array([1.0, 2.0, 0.75])), [2.0**-52, 2.0**-51, 2.0**-53]
    )


def test_mround_half_away_from_zero():
    x = np.array(
        [
            0.5,
            1.5,
            2.5,
            -0.5,
            -1.5,
            -2.5,
            PRED_HALF_D,
            -PRED_HALF_D,
            np.inf,
            -np.inf,
            np.nan,
            4503599627370497.0,
        ]
    )
    np.testing.assert_array_equal(
        dtypes.mround(x),
        [1, 2, 3, -1, -2, -3, 0, -0.0, np.inf, -np.inf, np.nan, 4503599627370497.0],
    )
    s = dtypes.mround(np.array([0.5, 2.5, PRED_HALF_S, 8388607.5], dtype=np.float32))
    assert s.dtype == np.float32
    np.testing.assert_array_equal(s, [1, 3, 0, 8388608])
    # NumPy's own rounds half to even
    assert np.round(2.5) == 2 and dtypes.mround(2.5) == 3


def test_cast_rounds_saturates_and_maps_nan_to_zero():
    x = np.array([np.nan, np.inf, -np.inf, 65535.5, -0.5, 0.5, 1.5, 70000.0, -3.0])
    np.testing.assert_array_equal(dtypes.cast(x, "uint16"), [0, 65535, 0, 65535, 0, 1, 2, 65535, 0])
    np.testing.assert_array_equal(dtypes.cast(x, "int8"), [0, 127, -128, 127, -1, 1, 2, 127, -3])
    big = np.array([2.0**63, -(2.0**63), 2.0**62, 2.0**64, 2.0**64 - 2048])
    np.testing.assert_array_equal(
        dtypes.cast(big, "int64"), [2**63 - 1, -(2**63), 2**62, 2**63 - 1, 2**63 - 1]
    )
    np.testing.assert_array_equal(
        dtypes.cast(big, "uint64"),
        np.array([2**63, 0, 2**62, 2**64 - 1, 2**64 - 2048], dtype=np.uint64),
    )
    # single values beyond int32 saturate (2^31 is a single, 2^31 - 1 is not)
    np.testing.assert_array_equal(
        dtypes.cast(np.array([2.0**31, -(2.0**31) - 256], np.float32), "int32"),
        [2**31 - 1, -(2**31)],
    )


def test_direct_call_on_single_arrays_rounds_the_value_below_a_half_up():
    # MATLAB: uint16(x) on a single array (2+ elements) is floor(|x| + 0.5) in single;
    # cast, feval and the class function on a scalar round exactly
    x = np.array([PRED_HALF_S, -PRED_HALF_S, 0.5, 1.5, -1.5], dtype=np.float32)
    np.testing.assert_array_equal(dtypes.uint16(x), [1, 0, 1, 2, 0])
    np.testing.assert_array_equal(dtypes.cast(x, "uint16"), [0, 0, 1, 2, 0])
    np.testing.assert_array_equal(dtypes.int16(x), [1, -1, 1, 2, -2])
    np.testing.assert_array_equal(dtypes.int8(x), [1, -1, 1, 2, -2])
    np.testing.assert_array_equal(dtypes.uint8(x), [1, 0, 1, 2, 0])
    # 32 and 64 bits, scalars and double data round exactly
    np.testing.assert_array_equal(dtypes.int32(x), [0, 0, 1, 2, -2])
    np.testing.assert_array_equal(dtypes.uint16(x[:1]), [0])
    assert dtypes.uint16(PRED_HALF_S) == 0
    np.testing.assert_array_equal(dtypes.uint16(np.array([PRED_HALF_D, PRED_HALF_D])), [0, 0])


def test_int_and_logical_conversions():
    i = np.array([-129, -1, 0, 200, 70000], dtype=np.int32)
    np.testing.assert_array_equal(dtypes.cast(i, "int8"), [-128, -1, 0, 127, 127])
    np.testing.assert_array_equal(dtypes.cast(i, "uint8"), [0, 0, 0, 200, 255])
    u = np.array([0, 2**63, 2**64 - 1], dtype=np.uint64)
    np.testing.assert_array_equal(dtypes.cast(u, "int64"), [0, 2**63 - 1, 2**63 - 1])
    np.testing.assert_array_equal(
        dtypes.cast(np.array([-5, 5], np.int64), "uint64"), np.array([0, 5], np.uint64)
    )
    assert dtypes.cast(np.array([2**53 + 1], np.int64), "double")[0] == 2.0**53
    np.testing.assert_array_equal(
        dtypes.logical(np.array([0, -0.0, 0.1, -2, np.inf])), [False, False, True, True, True]
    )
    with pytest.raises(ValueError, match="nologicalnan"):
        dtypes.logical(np.array([np.nan]))
    np.testing.assert_array_equal(dtypes.cast(np.array([True, False]), "uint16"), [1, 0])
    assert dtypes.single(np.array([3.5e38]))[0] == np.inf


def test_single_with_double_is_computed_in_double():
    rng = np.random.default_rng(1)
    a = (rng.random(10000) * 2000 - 1000).astype(np.float32)
    for b in (0.1, 1 / 3, np.pi, 1.1, 123.456):
        r = dtypes.times(a, b)
        assert r.dtype == np.float32
        np.testing.assert_array_equal(r, (a.astype(np.float64) * b).astype(np.float32))
        # NumPy computes in single with b rounded first, which differs
        assert not np.array_equal(r, a * np.float32(b))
    B = rng.random(10000) * 3 + 0.01
    np.testing.assert_array_equal(dtypes.rdivide(a, B), (a / B).astype(np.float32))
    assert dtypes.plus(a, B).dtype == np.float32 and (a + B).dtype == np.float64
    np.testing.assert_array_equal(
        dtypes.ldivide(a, 2.0), (2.0 / a.astype(np.float64)).astype(np.float32)
    )
    assert dtypes.plus(np.float32(1.5), True).dtype == np.float32
    assert dtypes.plus(np.ones(3, bool), np.ones(3, bool)).dtype == np.float64


def test_integer_arithmetic_rounds_and_saturates():
    u = np.array([0, 1, 2, 3, 5, 7, 65535, 65534, 32768], dtype=np.uint16)
    # measured in MATLAB R2024b (probe 1)
    np.testing.assert_array_equal(dtypes.rdivide(u, 2), [0, 1, 1, 2, 3, 4, 32768, 32767, 16384])
    np.testing.assert_array_equal(
        dtypes.rdivide(u, np.uint16(2)), [0, 1, 1, 2, 3, 4, 32768, 32767, 16384]
    )
    np.testing.assert_array_equal(dtypes.plus(u, 0.5), [1, 2, 3, 4, 6, 8, 65535, 65535, 32769])
    np.testing.assert_array_equal(dtypes.minus(u, 0.5), [0, 1, 2, 3, 5, 7, 65535, 65534, 32768])
    np.testing.assert_array_equal(dtypes.minus(0.5, u), [1, 0, 0, 0, 0, 0, 0, 0, 0])
    np.testing.assert_array_equal(dtypes.rdivide(u, 0), [0] + [65535] * 8)
    assert dtypes.rdivide(np.int16(-5), np.int16(0)) == -32768
    assert dtypes.uminus(np.int8(-128)) == 127
    # computed in double: uint16(1) + 0.49999999999999994 is 1.5 in double, so 2
    assert dtypes.plus(np.uint16(1), PRED_HALF_D) == 2
    assert dtypes.plus(u, "a")[1] == 98


def test_int64_arithmetic_is_exact():
    i = np.int64(9007199254740993)
    assert dtypes.plus(i, np.int64(2)) == 9007199254740995
    assert dtypes.times(i, np.int64(3)) == 27021597764222979
    assert dtypes.times(i, 3) == 27021597764222979
    assert dtypes.rdivide(i, 2) == 4503599627370497
    assert dtypes.times(i, 1.5) == 13510798882111490
    m = np.uint64(18446744073709551615)
    assert dtypes.minus(m, np.uint64(1)) == 18446744073709551614
    assert dtypes.rdivide(m, np.uint64(2)) == 9223372036854775808
    assert dtypes.times(np.uint64(18446744073709551614), 0.5) == 9223372036854775807
    assert dtypes.plus(np.int64(2**63 - 1), 1) == 2**63 - 1
    assert dtypes.times(np.int64(5), np.nan) == 0
    assert dtypes.rdivide(np.int64(5), 0.0) == 2**63 - 1


def test_int64_quirk_with_the_double_below_minus_a_half():
    # MATLAB R2024b: a positive int64 plus -0.49999999999999994 gives one too many (kept)
    a = np.array([1, 5, -5, 0], dtype=np.int64)
    np.testing.assert_array_equal(dtypes.minus(a, PRED_HALF_D), [2, 6, -5, 0])
    np.testing.assert_array_equal(dtypes.plus(a, -PRED_HALF_D), [2, 6, -5, 0])
    np.testing.assert_array_equal(dtypes.plus(-PRED_HALF_D, a), [2, 6, -5, 0])
    np.testing.assert_array_equal(dtypes.minus(-PRED_HALF_D, a), [-1, -5, 6, 0])
    np.testing.assert_array_equal(
        dtypes.minus(np.array([1, 5, 0, 0], dtype=np.uint64), PRED_HALF_D),
        np.array([1, 5, 0, 0], dtype=np.uint64),
    )
    assert dtypes.minus(np.int64(5), 0.5) == 5 and dtypes.plus(np.int64(5), PRED_HALF_D) == 5


def test_nan_constant_and_payloads():
    assert np.float64(dtypes.NAN).view(np.uint64) == 0xFFF8000000000000
    payload = np.array([0x7F9C5817], dtype=np.uint32).view(np.float32)  # a signalling NaN
    r = dtypes.plus(payload, np.nan)
    assert r.view(np.uint32)[0] == 0x7FDC5817  # the first NaN operand, quieted (as MATLAB)
    r = dtypes.plus(np.float32(1), np.array([np.nan]))
    assert np.isnan(r[0])


def test_mixed_classes_raise_as_in_matlab():
    u = np.array([1, 2], dtype=np.uint16)
    for other in (np.array([0.5, 0.5]), np.float32(0.3), np.int16(1), True):
        with pytest.raises(TypeError, match="mixedClasses"):
            dtypes.times(u, other)
    assert dtypes.result_class(u, 0.5) == "uint16"
    assert dtypes.result_class(np.float32(1), np.zeros(3)) == "single"
    assert dtypes.result_class(True, 1) == "double"


def test_idivide():
    a = np.array([-7, -5, 5, 7], dtype=np.int16)
    np.testing.assert_array_equal(dtypes.idivide(a, np.int16(2)), [-3, -2, 2, 3])
    np.testing.assert_array_equal(dtypes.idivide(a, np.int16(2), "round"), [-4, -3, 3, 4])
    np.testing.assert_array_equal(dtypes.idivide(a, np.int16(2), "floor"), [-4, -3, 2, 3])
    np.testing.assert_array_equal(dtypes.idivide(a, np.int16(2), "ceil"), [-3, -2, 3, 4])
    assert dtypes.idivide(np.int64(-7), np.int64(2), "floor") == -4
    with pytest.raises(TypeError):
        dtypes.idivide(7.0, 2.0)


def test_implicit_expansion_is_matlabs():
    a = np.ones((3, 4, 5))
    b = np.arange(12.0).reshape(3, 4)  # MATLAB expands it along the third dimension
    r = dtypes.times(a, b)
    assert r.shape == (3, 4, 5)
    np.testing.assert_array_equal(r[:, :, 2], b)
    assert dtypes.plus(np.ones(4), np.ones(4)).shape == (4,)


def test_limits():
    assert dtypes.intmax("uint16") == 65535 and dtypes.intmin("int8") == -128
