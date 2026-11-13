"""Conversions and selects whose bits differed by compiler or processor (docs/platforms.md)."""

import numpy as np
import pytest

from petakit5d.kernels import ops

# 64-bit integers whose float32 differs when rounded through double: 2^53 + 2^29 + 1 is just
# above the halfway point between two floats, but its double (2^53 + 2^29) is on it
_INT64 = np.array(
    [2**53 + 2**29 + 1, -(2**53 + 2**29 + 1), 2**62 + 2**38 + 1, 2**63 - 1, -(2**63), 3, 0],
    dtype=np.int64,
)


def _once(v: int) -> float:
    """v rounded once to the nearest float32 (ties to even), with exact integer arithmetic"""
    if v == 0:
        return 0.0
    s, a = (-1 if v < 0 else 1), abs(v)
    e = a.bit_length() - 24
    if e <= 0:
        return float(v)
    q, r = divmod(a, 1 << e)
    half = 1 << (e - 1)
    if r > half or (r == half and q & 1):
        q += 1
    return s * float(q << e)


@pytest.mark.parametrize("dtype", [np.int64, np.uint64])
def test_cast_of_64_bit_integers_to_single_rounds_once(dtype):
    x = _INT64 if dtype == np.int64 else np.abs(_INT64[_INT64 >= 0]).astype(np.uint64)
    got = ops.convert(x, "single", mode="cast")
    want = np.array([_once(int(v)) for v in x], dtype=np.float32)
    assert got.tobytes() == want.tobytes()
    # 'class' (MATLAB's single()) rounds through double on purpose
    via_double = np.array([np.float32(float(int(v))) for v in x], dtype=np.float32)
    assert ops.convert(x, "single", mode="class").tobytes() == via_double.tobytes()


def test_clamps_keep_negative_zero():
    # min(v, hi) and max(v, lo) give v on a tie, also -0 against +0 (MSVC made it +0)
    x = np.array([-0.0, 0.0, -1.0, 2.0], dtype=np.float32)
    got = ops.scale_offset(x, lo=0.0, hi=0.0)
    assert list(got.view(np.uint32)) == [0x80000000, 0, 0, 0]
