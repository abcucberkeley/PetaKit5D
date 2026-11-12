"""The oracle cases' inputs are the same bytes on every platform (the bytes MATLAB read).

The references in oracle_f02.json are only valid for the inputs they were made from. The first
CI run found inputs built with ``10.0 ** k`` that differ by an ulp between C libraries (glibc
2.28 rounds 10^23 correctly, 2.35 and later do not) and between NumPy's SIMD paths (its AVX-512
power is 1 ulp low at 10^-5); the generators now build powers of ten without ``pow``.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import _f02_gen as gen  # noqa: E402
import f02_oracle_cases as oc  # noqa: E402

# SHA-256 of the column-major bytes MATLAB read (oracle runs of 2026-11-11)
MATLAB_INPUTS = {
    "fmt_values": "0d57495d1ddc0d94935d3030c7be2d20194a0e225a1d7ae5226a089ea385b840",
    "rows_ints": "d36d6d490190143e76efea4eea5416a5fa92289fcaf6abbbee03ab1a61888f9a",
    "rows_floats": "82bcd62e0007998f935ac477201f80a803b19091959a11efa8609ca27eda16ca",
}


def _digest(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a.ravel(order="F")).tobytes()).hexdigest()


def _inputs() -> dict[str, str]:
    return {
        "fmt_values": _digest(oc.fmt_values()),
        "rows_ints": _digest(oc._rows(oc.SEED + 8, 300, True)),
        "rows_floats": _digest(oc._rows(oc.SEED + 8, 300, False)),
    }


def test_inputs_are_the_bytes_matlab_read():
    assert _inputs() == MATLAB_INPUTS


def test_powers_of_ten_do_not_use_pow():
    # correctly rounded; glibc 2.35+ pow(10, 23) gives 0x1.52d02c7e14af7p+76
    assert gen.pow10(23) == 1e23 == float.fromhex("0x1.52d02c7e14af6p+76")
    assert gen.pow10(-5) == 1e-5
    k = np.array([-17, -5, 0, 7, 23])
    np.testing.assert_array_equal(
        gen.pow10_array(k), [1e-17, 9.999999999999999e-06, 1.0, 1e7, 1e23]
    )
    with pytest.raises(ValueError):
        gen.pow10_array(np.array([0.5]))


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="NumPy's x86-64 dispatch groups")
def test_inputs_do_not_depend_on_numpy_simd():
    # NumPy's power on AVX-512 and on AVX2 differ in the last bit for some powers of ten
    code = (
        "import sys; sys.path.insert(0, sys.argv[1]); import test_oracle_inputs as t; "
        "print(t._inputs() == t.MATLAB_INPUTS)"
    )
    env = dict(os.environ, NPY_DISABLE_CPU_FEATURES="AVX512_ICL AVX512_SPR X86_V4")
    r = subprocess.run(
        [sys.executable, "-c", code, str(Path(__file__).parent)],
        env=env,
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        pytest.skip(f"NumPy refused the dispatch setting: {r.stderr.strip().splitlines()[-1:]}")
    assert r.stdout.strip() == "True"
