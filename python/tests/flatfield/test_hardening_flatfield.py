"""Tests the mutation campaign asked for (E06, docs/hardening.md): flat field."""

import numpy as np
import pytest

from petakit5d.core import dtypes
from petakit5d.flatfield import correction, estimation
from petakit5d.io.tiff import write_tiff
from petakit5d.jobs import framework

FSN = (
    "Scan_Iter_0000_0000_{}_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_"
    "000z_0000t.tif"
)


def _bits(a):
    return np.asarray(a, dtype=np.float32).view(np.uint32)


@pytest.mark.parametrize(
    "op, matlab_op", [(np.subtract, dtypes.minus), (np.divide, dtypes.rdivide)]
)
def test_single_op_gives_matlabs_nan_bits(op, matlab_op):
    # an invalid operation (MATLAB's NaN has the sign bit set) and NaN operands (their payloads)
    nan_a = np.array([0x7FC00011], dtype=np.uint32).view(np.float32)[0]
    nan_b = np.array([0xFFC00022], dtype=np.uint32).view(np.float32)[0]
    a = np.array([np.inf, nan_a, 1, nan_a, 4], dtype=np.float32)
    b = np.array([np.inf, 1, nan_b, nan_b, 2], dtype=np.float32)
    if op is np.divide:
        a[0] = b[0] = 0
    r = correction._single_op(op, matlab_op, a, b)
    assert r is not None and r.dtype == np.float32
    np.testing.assert_array_equal(_bits(r), _bits(matlab_op(a, b)))
    assert _bits(r)[0] == 0xFFC00000
    # (no NaN: NumPy's operation)
    clean = correction._single_op(op, matlab_op, a[4:], b[4:])
    np.testing.assert_array_equal(clean, op(a[4:], b[4:]))


def _dataset(folder):
    rng = np.random.default_rng(4)
    for cam in ("CamA", "CamB"):
        dark = np.clip(100 + rng.normal(0, 3, size=(16, 40, 36)), 0, None).astype(np.uint16)
        write_tiff(dark, str(folder / FSN.format(cam)))
    return str(folder)


KW = {"channel_patterns": ("CamA_ch0", "CamB_ch0"), "parse_cluster": False}


def test_a_finished_call_is_not_run_again(tmp_path, monkeypatch):
    real = framework.run_tasks
    calls = []

    def spy(inputs, outputs, tasks, **kw):
        calls.append(kw["mem_allocate"])
        return real(inputs, outputs, tasks, **kw)

    monkeypatch.setattr(framework, "run_tasks", spy)
    estimation.flatfield_background_dataset(_dataset(tmp_path), **KW)
    assert len(calls) == 1


def test_a_failed_call_is_tried_with_twice_the_memory_then_fails(tmp_path, monkeypatch):
    calls = []

    def spy(inputs, outputs, tasks, **kw):
        calls.append(kw["mem_allocate"])
        return [False] * len(outputs)

    monkeypatch.setattr(framework, "run_tasks", spy)
    with pytest.raises(estimation.MatlabError, match="Some output files are missing"):
        estimation.flatfield_background_dataset(_dataset(tmp_path), **KW)
    assert len(calls) == 2 and calls[1] == 2 * calls[0] > 0
