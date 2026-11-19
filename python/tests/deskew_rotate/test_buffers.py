"""Memory kept from file to file (``deskew_rotate.file.keep_buffers``) and the ``out`` of
``deskew_rotate_volume``: the same files, bit for bit, with and without the kept buffers (tiff
and zarr results, uint16 and single, the warp and the fused kernel, flat field, MATLAB's zarr
layout, files of other sizes, preprocessing copies, the separate path), the buffers reused,
grown, private to the thread that keeps them and freed when the loop ends, and the local loops
of the dataset wrapper and the pipeline running inside them.
"""

from __future__ import annotations

import contextlib
import os
import sys
import threading

import numpy as np
import pytest

from petakit5d.deskew_rotate import file as df
from petakit5d.deskew_rotate import volume as dv
from petakit5d.io import read_tiff, read_zarr
from petakit5d.testing import synth

PX = 0.108
FSN = "Scan_Iter_0000_0000_CamA_ch0_stack0000_488nm_0000000msec_0106060251msecAbs"
SHAPES = ((24, 40, 36), (31, 40, 36), (17, 40, 36), (24, 40, 36))


def _frame(shape, kind, seed):
    rng = np.random.default_rng(seed)
    v = rng.integers(90, 4000, shape).astype(np.uint16)
    return v if kind == "uint16" else (v.astype(np.float32) + 0.25)


def _read(path):
    return read_zarr(path) if str(path).endswith(".zarr") else read_tiff(path)


def _write_files(root, kind, zarr_layout=None):
    """one folder per file (the same name in each), with the flat-field images"""
    paths = []
    for k, shape in enumerate(SHAPES):
        d = root / f"f{k}"
        d.mkdir(parents=True)
        v = _frame(shape, kind, k)
        if zarr_layout is not None:
            synth.write_zarr(d / f"{FSN}.zarr", v, (8, 16, 16), layout=zarr_layout)
            paths.append(str(d / f"{FSN}.zarr"))
        else:
            synth.write_tiff(d / f"{FSN}.tif", v)
            paths.append(str(d / f"{FSN}.tif"))
    bg = np.zeros(SHAPES[0][1:], dtype=np.float32) + 100
    ff = bg + np.linspace(500, 1000, SHAPES[0][2]).astype(np.float32)
    synth.write_tiff(root / "ff.tif", ff)
    synth.write_tiff(root / "bg.tif", bg)
    return paths


def _results(path, save_zarr, separate=False):
    """the result files of a file's call: (name, array)"""
    folder = os.path.dirname(path)
    ext = ".zarr" if save_zarr else ".tif"
    out = []
    for sub in ("DS", "DSR") if separate else ("DSR",):
        out.append((f"{sub}/{FSN}{ext}", _read(f"{folder}/{sub}/{FSN}{ext}")))
        out.append((f"{sub}/MIP", _read(f"{folder}/{sub}/MIPs/{FSN}_MIP_z{ext}")))
    return out


def _run_all(tmp_path, label, kind, keep, zarr_layout=None, **kw):
    """deskew_rotate_file on each file in turn, inside keep_buffers or not; the results"""
    paths = _write_files(tmp_path / label, kind, zarr_layout)
    dz = kw.pop("dz", 0.3)
    if kw.pop("ff", False):
        kw.update(
            ff_correction=True,
            ff_image=str(tmp_path / label / "ff.tif"),
            background_image=str(tmp_path / label / "bg.tif"),
            const_offset=1,
        )
    with df.keep_buffers() if keep else contextlib.nullcontext():
        for p in paths:
            df.deskew_rotate_file(p, PX, dz, reverse=True, num_threads=3, **kw)
    separate = not kw.get("dsr_combined", True)
    return [_results(p, kw.get("save_zarr", False), separate) for p in paths]


def _same(a, b):
    for ra, rb in zip(a, b, strict=True):
        for (na, xa), (nb, xb) in zip(ra, rb, strict=True):
            assert na == nb
            assert xa.dtype == xb.dtype and xa.shape == xb.shape, na
            assert xa.tobytes() == xb.tobytes(), na


@pytest.mark.parametrize("kind", ("uint16", "single"))
@pytest.mark.parametrize("dz", (0.2, 0.3))
@pytest.mark.parametrize("ff", (False, True))
@pytest.mark.parametrize("save_zarr", (False, True))
def test_kept_buffers_give_the_same_files(tmp_path, kind, dz, ff, save_zarr):
    """files of other sizes one after another (the buffers grow, then hold smaller arrays)"""
    kw = {"dz": dz, "ff": ff, "save_zarr": save_zarr, "save_16bit": kind == "uint16"}
    _same(
        _run_all(tmp_path, "fresh", kind, False, **dict(kw)),
        _run_all(tmp_path, "kept", kind, True, **dict(kw)),
    )


@pytest.mark.parametrize(
    "kw",
    [
        {"zarr_layout": "matlab", "save_zarr": True},
        {"zarr_layout": "matlab", "save_zarr": False},
        {"zarr_layout": "native", "save_zarr": True},
        {"input_bbox": (1, 2, 0, 15, 37, 36)},
        {"flip_zstack": True, "ff": True},
        {"input_axis_order": "zxy"},
        {"dsr_combined": False, "rotate": True},
        {"save_16bit": False},
        {"interp_method": "cubic"},
    ],
    ids=lambda kw: "-".join(f"{k}={v}" for k, v in kw.items()),
)
def test_kept_buffers_on_the_other_paths(tmp_path, kw):
    """MATLAB's zarr layout (the mex's loops into a kept buffer set to 0 first), native zarr,
    the preprocessing copies, the separate path, single results of uint16 data and imwarp"""
    kw = dict(kw)
    layout = kw.pop("zarr_layout", None)
    _same(
        _run_all(tmp_path, "fresh", "uint16", False, layout, **dict(kw)),
        _run_all(tmp_path, "kept", "uint16", True, layout, **dict(kw)),
    )


def test_kept_buffers_with_a_frame_narrower_than_the_x_step(tmp_path):
    """frames narrower than the x step (9.4 pixels at dz 1.2: the kernels refuse fewer than 12
    columns and MATLAB's code runs, giving single where the kernels give uint16): no kept output
    for them, the same files"""
    res = []
    for label, keep in (("fresh", False), ("kept", True)):
        d = tmp_path / label
        d.mkdir()
        paths = []
        for k, shape in enumerate(((20, 30, 10), (24, 30, 11), (20, 30, 40))):
            (d / f"f{k}").mkdir()
            synth.write_tiff(d / f"f{k}" / f"{FSN}.tif", _frame(shape, "uint16", k))
            paths.append(str(d / f"f{k}" / f"{FSN}.tif"))
        with df.keep_buffers() if keep else contextlib.nullcontext():
            for p in paths:
                df.deskew_rotate_file(p, PX, 1.2, reverse=True, num_threads=2)
        res.append([_results(p, False) for p in paths])
    _same(*res)


def test_kept_buffers_are_reused_and_freed(tmp_path):
    paths = _write_files(tmp_path, "uint16")
    with df.keep_buffers() as kept:
        assert df._kept() is kept
        df.deskew_rotate_file(paths[0], PX, 0.3, reverse=True, num_threads=2)
        mem = [b._bytes for b in kept.buffers]
        assert all(m is not None for m in mem)
        # the same size again: the same memory, no new arrays
        df.deskew_rotate_file(paths[3], PX, 0.3, reverse=True, num_threads=2)
        assert all(b._bytes is m for b, m in zip(kept.buffers, mem, strict=True))
        # a larger file: the buffers grow; a smaller one keeps them
        df.deskew_rotate_file(paths[1], PX, 0.3, reverse=True, num_threads=2)
        grown = [b._bytes for b in kept.buffers]
        assert sum(m.size for m in grown) > sum(m.size for m in mem)
        df.deskew_rotate_file(paths[2], PX, 0.3, reverse=True, num_threads=2)
        assert all(b._bytes is m for b, m in zip(kept.buffers, grown, strict=True))
        with df.keep_buffers() as inner:
            assert inner is kept
        assert df._kept() is kept
    assert df._kept() is None
    assert all(b._bytes is None for b in kept.buffers)


def test_kept_buffers_belong_to_their_thread():
    seen = []
    with df.keep_buffers():
        t = threading.Thread(target=lambda: seen.append(df._kept()))
        t.start()
        t.join()
        assert df._kept() is not None
    assert seen == [None]


def test_kept_buffers_are_freed_after_an_error(tmp_path):
    with pytest.raises(RuntimeError), df.keep_buffers() as kept:
        kept.read.view((4, 5, 6), np.uint16)
        raise RuntimeError("stop")
    assert df._kept() is None and kept.read._bytes is None


def test_byte_buffer_views():
    """views of any class on the first bytes; the memory grows to the largest view; other()
    is the buffer that does not hold an array"""
    fb = df.FileBuffers()
    a = fb.read.view((3, 4, 5), np.uint16)
    assert a.shape == (3, 4, 5) and a.dtype == np.uint16 and a.flags.c_contiguous
    assert a.flags.writeable
    mem = fb.read._bytes
    b = fb.read.view((2, 4, 5), np.float32)  # 160 bytes, more than the 120 there: it grows
    assert b.dtype == np.float32 and fb.read._bytes is not mem
    c = fb.read.view((1, 2, 3), np.float32)
    assert np.shares_memory(b, c)
    assert fb.other(c) is fb.buffers[1]
    d = fb.other(c).view((2, 2, 2), np.uint16)
    assert fb.other(d) is fb.buffers[0]
    assert fb.other(np.zeros(3)) is fb.buffers[0]
    fb.release()
    assert all(x._bytes is None for x in fb.buffers)


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="the mapping is Linux's")
def test_kept_memory_without_numpy_advice(monkeypatch):
    """without NUMPY_MADVISE_HUGEPAGE, kept memory is an anonymous mapping NumPy does not
    advise; with it set, NumPy's own memory (its advice as the user set it)"""
    import mmap

    def mapped(x):
        # (np.frombuffer's base is a memoryview of the mapping)
        return isinstance(getattr(x.base, "obj", None), mmap.mmap)

    monkeypatch.delenv("NUMPY_MADVISE_HUGEPAGE", raising=False)
    a = df._allocate(1 << 20)
    assert mapped(a) and a.flags.writeable and a.size == 1 << 20
    a[:] = 7
    monkeypatch.setenv("NUMPY_MADVISE_HUGEPAGE", "1")
    b = df._allocate(1 << 20)
    assert not mapped(b) and b.size == 1 << 20


def test_dataset_and_pipeline_keep_buffers_locally(tmp_path, monkeypatch):
    """the files of deskew_rotate_dataset and automatic_processing in local mode run with kept
    buffers (one set for the loop), freed when the loop ends"""
    from petakit5d.deskew_rotate.dataset import deskew_rotate_dataset
    from petakit5d.pipeline import automatic_processing

    seen = []
    real = df._kept

    def spy():
        k = real()
        seen.append(k)
        return k

    monkeypatch.setattr(df, "_kept", spy)
    for run in ("dataset", "pipeline"):
        d = tmp_path / run
        d.mkdir()
        for yi in range(2):
            synth.write_tiff(
                d / f"Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_"
                f"0000000000msecAbs_000x_{yi:03d}y_000z_0000t.tif",
                _frame((12, 20, 24), "uint16", yi),
            )
        seen.clear()
        if run == "dataset":
            deskew_rotate_dataset(
                str(d), channel_patterns=["CamA_ch0"], dz=0.3, reverse=True, parse_cluster=False
            )
        else:
            automatic_processing(
                str(d),
                channel_patterns=["CamA_ch0"],
                dz=0.3,
                reverse=True,
                streaming=False,
                parse_cluster=False,
            )
        assert len(seen) == 2 and seen[0] is not None and seen[0] is seen[1], run
        assert real() is None
        assert all(b._bytes is None for b in seen[0].buffers)


def test_pipeline_gives_the_memory_back_before_a_volume_is_stitched(tmp_path, monkeypatch):
    """locally, the pipeline frees the kept buffers before it stitches a volume whose files are
    all deskewed (two time points of two tiles: twice, each right before that stitching), and
    keeps them across the stitching calls that find a volume incomplete"""
    from petakit5d.pipeline import automatic as A

    data = tmp_path / "d"
    for it in (0, 1):
        vol = _frame((16, 40, 36), "uint16", it)
        tiles = synth.cut_tiles(vol, [(0, 0, 0), (0, 16, 0)], (16, 24, 36), (2 * PX, PX, PX))
        _, csv = synth.write_stitch_tiles(str(data), tiles, cam="A", iter=it)
    events = []
    real_release, real_run = df.release_kept_buffers, A._State._run_here

    def release():
        events.append("free")
        real_release()

    def run_here(self, func_str):
        events.append("stitch" if func_str.startswith("XR_matlab_stitching") else "deskew")
        return real_run(self, func_str)

    monkeypatch.setattr(df, "release_kept_buffers", release)
    monkeypatch.setattr(A._State, "_run_here", run_here)
    A.automatic_processing(
        str(data) + "/",
        xy_pixel_size=PX,
        dz=0.3,
        reverse=True,
        channel_patterns=("CamA_ch0",),
        streaming=False,
        parse_cluster=False,
        stitch=True,
        image_list_fullpaths=[csv],
        axis_order="-x,y,z",
    )
    assert events == [
        "deskew", "stitch", "deskew", "free", "stitch",
        "deskew", "stitch", "deskew", "free", "stitch",
    ], events  # fmt: skip
    assert df._kept() is None


def test_deskew_rotate_volume_out():
    """``out`` and ``num_threads`` of deskew_rotate_volume: the result written into out (every
    voxel) and returned, on the kernel paths and the imwarp paths alike; the wrong class
    refused"""
    vol = _frame((24, 40, 36), "uint16", 3)
    for kw in ({}, {"interp_method": "cubic"}, {"resample_factor": (2, 2, 2)}):
        ref = dv.deskew_rotate_volume(vol, 32.45, 0.3, PX, reverse=True, save_16bit=True, **kw)
        out = np.full(ref.shape, 7, dtype=ref.dtype)
        got = dv.deskew_rotate_volume(
            vol, 32.45, 0.3, PX, reverse=True, save_16bit=True, num_threads=2, out=out, **kw
        )
        assert got is out
        np.testing.assert_array_equal(got, ref)
        with pytest.raises(ValueError, match="out must be"):
            dv.deskew_rotate_volume(
                vol,
                32.45,
                0.3,
                PX,
                reverse=True,
                save_16bit=True,
                out=np.empty(ref.shape, dtype=np.float64),
                **kw,
            )
