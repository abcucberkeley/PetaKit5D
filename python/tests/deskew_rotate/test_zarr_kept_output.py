"""The large-file deskew/rotation's batches (``deskew_rotate_block``) put the kernels' results
into the memory of the task's first result, kept over its batches, which the kernels on
MATLAB's layout write in full (``fill``): the same result zarr and MIPs as with a new output per
batch, for both scan directions, the warp and the fused kernel, uint16 and single results,
MATLAB's and native inputs, flipped batches; the batches of a task reuse the first one's memory;
the ``imwarp`` paths and the inputs the kernels refuse take new outputs, as before.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.deskew_rotate import volume as dv
from petakit5d.deskew_rotate import zarr as dz
from petakit5d.io import read_tiff, read_zarr, write_zarr
from petakit5d.testing import synth
from petakit5d.testing.synth import CounterRNG

SKEW = 32.45
PX = 0.108
RAW = (29, 37, 61)
FSN = "scan_Iter_0000_CamA_ch0_stack0000_0000000msec_0000t"
KW = {
    "skew_angle": SKEW,
    "batch_size": (1024, 16, 1024),
    "block_size": (64, 16, 64),
    "parse_cluster": False,
}


def _frame(kind):
    rng = CounterRNG(4, f"dsrzarr/kept/{kind}")
    a = rng.integers(0, 3000, RAW, dtype=np.uint16)
    return a.astype(np.float32) if kind == "single" else a


def _input(folder, kind, layout):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    fn = str(folder / f"{FSN}.zarr")
    frame = _frame(kind)
    if layout == "matlab":
        synth.write_zarr(fn, frame, (8, 16, 16), "matlab")
    else:
        write_zarr(frame, fn, block_size=(8, 16, 16))
    return fn


def _results(folder):
    """the result zarr's voxels and the MIP tiffs of a run"""
    d = Path(folder) / "DSR"
    out = {"zarr": read_zarr(str(d / f"{FSN}.zarr"))}
    for name in sorted(os.listdir(d / "MIPs")):
        if name.endswith(".tif"):
            out[name] = read_tiff(str(d / "MIPs" / name))
    return out


@pytest.mark.parametrize(
    "kind, save_16bit, dz_, reverse, layout, flip",
    [
        ("uint16", True, 0.3, True, "matlab", False),
        ("uint16", True, 0.2, False, "matlab", False),
        ("uint16", False, 0.3, True, "native", False),
        ("single", True, 0.3, False, "matlab", True),
        ("single", False, 0.2, True, "native", False),
        ("uint16", True, 0.5, True, "matlab", True),
    ],
)
def test_kept_outputs_give_the_results_of_new_ones(
    tmp_path, monkeypatch, kind, save_16bit, dz_, reverse, layout, flip
):
    # (one task for all batches)
    kw = {**KW, "reverse": reverse, "save_16bit": save_16bit, "flip_zstack": flip, "task_size": 10}
    fn = _input(tmp_path / "kept", kind, layout)
    used = []
    real = dz._kernel_out

    def spy(*a, **k):
        used.append(real(*a, **k))
        return used[-1]

    monkeypatch.setattr(dz, "_kernel_out", spy)
    dz.deskew_rotate_zarr(fn, PX, dz_, **kw)
    assert len(used) > 2 and all(k for k, _ in used), "the kernels did not compute the batches"
    assert used[0][1] is None and any(o is not None for _, o in used), "no result was kept"
    fn2 = _input(tmp_path / "new", kind, layout)
    monkeypatch.setattr(dz, "_kernel_out", lambda *a, **k: (False, None))
    dz.deskew_rotate_zarr(fn2, PX, dz_, **kw)
    a, b = _results(tmp_path / "kept"), _results(tmp_path / "new")
    assert sorted(a) == sorted(b)
    for k in a:
        assert a[k].dtype == b[k].dtype and a[k].shape == b[k].shape, k
        assert a[k].tobytes() == b[k].tobytes(), k


def test_the_batches_of_a_task_reuse_one_output(tmp_path, monkeypatch):
    fn = _input(tmp_path, "uint16", "matlab")
    seen = []
    real = dv.deskew_rotate_volume_matlab_layout

    def spy(*a, out=None, **k):
        res = real(*a, out=out, **k)
        seen.append((None if out is None else out.__array_interface__["data"][0],
                     res.__array_interface__["data"][0]))  # fmt: skip
        return res

    monkeypatch.setattr(dv, "deskew_rotate_volume_matlab_layout", spy)
    # (one task for all batches)
    dz.deskew_rotate_zarr(fn, PX, 0.3, reverse=True, save_16bit=True, task_size=10, **KW)
    assert len(seen) > 2 and seen[0][0] is None
    first = seen[0][1]
    assert all(o == r == first for o, r in seen[1:]), seen


def _vol(kind, shape=RAW):
    v = _frame(kind)[: shape[0], : shape[1], : shape[2]]
    return np.ascontiguousarray(v.transpose(0, 2, 1))


@pytest.mark.parametrize(
    "kw, kind, shape",
    [
        ({"interp_method": "cubic"}, "uint16", RAW),
        ({"resample_factor": (1.5, 1.5, 1.5)}, "uint16", RAW),
        ({"batch_rows": (0, 37)}, "uint16", RAW),
        ({}, "float64", RAW),
        # (a frame narrower than the x step: MATLAB's code, in single)
        ({}, "uint16", (29, 37, 3)),
    ],
)
def test_no_kept_output_off_the_kernels(kw, kind, shape):
    """imwarp would copy its result into a kept output: those paths get none"""
    if kind == "float64":
        v = _vol("single", shape).astype(np.float64)
    else:
        v = _vol(kind, shape)
    args = {
        "reverse": True,
        "resample_factor": None,
        "save_16bit": True,
        "interp_method": "linear",
        "batch_rows": None,
        **kw,
    }
    kept = np.zeros(1 << 20, dtype=np.uint8)
    assert dz._kernel_out(kept, v, (0, 0, 0, 5, 6, 7), SKEW, 0.3, PX, **args) == (False, None)


@pytest.mark.parametrize("kind, save_16bit, dtype", [
    ("uint16", True, np.uint16), ("uint16", False, np.float32), ("single", True, np.float32),
    ("u16_kept", False, np.float32),
])  # fmt: skip
def test_kept_output_has_the_box_in_matlabs_layout(kind, save_16bit, dtype):
    """a view of the kept memory of the box's size in MATLAB's layout and the kernels' class
    (uint16 input without save_16bit: single, as the kernels give it; the batches turn such
    input to single first); none (a new output) without kept memory or with too little"""
    v = _vol("uint16" if kind == "u16_kept" else kind)
    if not save_16bit and kind != "u16_kept":
        v = v.astype(np.float32)
    box = (2, 3, 4, 12, 23, 34)
    kw = {"reverse": False, "resample_factor": None, "save_16bit": save_16bit,
          "interp_method": "linear", "batch_rows": None}  # fmt: skip
    n = 10 * 30 * 20 * np.dtype(dtype).itemsize
    kept = np.zeros(n + 8, dtype=np.uint8)
    k, o = dz._kernel_out(kept, v, box, SKEW, 0.3, PX, **kw)
    assert k and o is not None and o.dtype == dtype and o.shape == (10, 30, 20)
    assert np.shares_memory(o, kept) and o.flags.c_contiguous and o.flags.writeable
    assert dz._kernel_out(kept[: n - 1], v, box, SKEW, 0.3, PX, **kw) == (True, None)
    assert dz._kernel_out(None, v, box, SKEW, 0.3, PX, **kw) == (True, None)
