"""The tiff conversion's read buffer kept over the files of a loop (``io.convert.keep_read_buffer``)
and the tile function's result in the volume the conversion owns (the edge erosion in place).

The arrays written are those of the conversion without either: for each option (boxes, the flip,
the z blocks, partial files, images, an offset before the erosion, other functions), for files
of growing and shrinking sizes and of other classes in one loop; the buffer is reused, freed when
the loop ends (also after an error), shared by nested blocks and kept per thread; a frame the
caller gives is never written; the dataset wrapper keeps the buffer in local mode only. Also the
memory of kept buffers (``io.imagesize.kept_memory``, ``ByteBuffer``).
"""

from __future__ import annotations

import contextlib
import sys
import threading

import numpy as np
import pytest

from petakit5d.io import convert, read_zarr
from petakit5d.io.imagesize import ByteBuffer, kept_memory
from petakit5d.stitch import tile_functions as TF
from petakit5d.testing.synth import CounterRNG, write_tiff

ERODE = "@(x)erodeVolumeBy2DProjection(x,2)"
SIZES = [(13, 41, 37), (9, 30, 26), (15, 44, 39)]  # (z, y, x): smaller, then larger


def _vol(shape, dtype=np.uint16, name="v"):
    """positive values with zero columns near the edges, so the y projection is not full"""
    rng = CounterRNG(21, f"x11/convert/{name}/{'x'.join(map(str, shape))}")
    v = rng.integers(1, 3000, shape).astype(dtype)
    if len(shape) == 3:
        v[:, :, :2] = 0
        v[:2, :, shape[2] // 3] = 0
        v[-3:, :, -4:] = 0
    if np.dtype(dtype).kind == "f":
        v = v / np.asarray(7, dtype)
    return np.ascontiguousarray(v)


def _files(tmp_path, dtype=np.uint16, sizes=SIZES):
    vols, paths = [], []
    for k, s in enumerate(sizes):
        v = _vol(s, dtype, f"f{k}")
        p = str(tmp_path / f"in_{k}.tif")
        write_tiff(p, v)
        vols.append(v)
        paths.append(p)
    return vols, paths


def _convert_all(tmp_path, paths, tag, keep, **kw):
    """tiff_to_zarr on each file in turn, inside keep_read_buffer or not; the arrays read back"""
    out = []
    with convert.keep_read_buffer() if keep else contextlib.nullcontext():
        for k, p in enumerate(paths):
            z = str(tmp_path / f"{tag}_{k}.zarr")
            convert.tiff_to_zarr(p, z, block_size=(4, 16, 16), **kw)
            out.append(read_zarr(z))
    return out


OPTIONS = [
    {},
    {"usr_fcn": ERODE},
    {"usr_fcn": ERODE, "flip_zstack": True},
    {"usr_fcn": ERODE, "tile_out_bbox": (1, 0, 0, 8, 30, 26)},  # z only: a view of the read
    {"usr_fcn": ERODE, "tile_out_bbox": (1, 2, 3, 8, 28, 24)},
    {"usr_fcn": ERODE, "input_bbox": (1, 2, 3, 9, 29, 25)},
    {"usr_fcn": "@(x)erodeVolumeBy2DProjection(x+3,1)"},
    {"usr_fcn": "@(x)x+100"},
    {"read_whole_tiff": False},
    {"read_whole_tiff": False, "input_bbox": (2, 1, 0, 9, 30, 26)},
]


OPTION_IDS = ["plain", "erode", "flip", "z_crop", "crop", "input_bbox", "offset", "plus",
              "z_blocks", "z_blocks_box"]  # fmt: skip


@pytest.mark.parametrize("dtype", [np.uint16, np.float32])
@pytest.mark.parametrize("kw", OPTIONS, ids=OPTION_IDS)
def test_kept_buffer_gives_the_same_arrays(tmp_path, dtype, kw):
    _, paths = _files(tmp_path, dtype)
    want = _convert_all(tmp_path, paths, "new", False, **kw)
    got = _convert_all(tmp_path, paths, "kept", True, **kw)
    for g, w in zip(got, want, strict=True):
        assert g.dtype == w.dtype and np.array_equal(g, w, equal_nan=True)


def test_the_arrays_are_the_erosion_of_the_files(tmp_path):
    vols, paths = _files(tmp_path)
    got = _convert_all(tmp_path, paths, "kept", True, usr_fcn=ERODE)
    for g, v in zip(got, vols, strict=True):
        assert np.array_equal(g, TF.erode_volume_by_2d_projection(v, 2))


def test_classes_and_images_in_one_loop(tmp_path):
    vols = [_vol((7, 21, 19), np.uint16, "a"), _vol((5, 18, 20), np.float32, "b"),
            _vol((21, 19), np.uint16, "im"), _vol((6, 20, 17), np.uint8, "c")]  # fmt: skip
    paths = []
    for k, v in enumerate(vols):
        paths.append(str(tmp_path / f"m_{k}.tif"))
        write_tiff(paths[-1], v)
    want = _convert_all(tmp_path, paths, "new", False, usr_fcn=ERODE)
    got = _convert_all(tmp_path, paths, "kept", True, usr_fcn=ERODE)
    for g, w in zip(got, want, strict=True):
        assert g.dtype == w.dtype and g.shape == w.shape and np.array_equal(g, w)


def test_partial_files_in_the_loop(tmp_path):
    v = _vol((13, 41, 37), np.uint16, "parts")
    parts = [str(tmp_path / "p1.tif"), str(tmp_path / "p2.tif")]
    write_tiff(parts[0], v[:5])
    write_tiff(parts[1], v[5:])
    with convert.keep_read_buffer():
        convert.tiff_to_zarr(parts, str(tmp_path / "parts.zarr"), usr_fcn=ERODE)
    assert np.array_equal(
        read_zarr(str(tmp_path / "parts.zarr")), TF.erode_volume_by_2d_projection(v, 2)
    )


def test_kept_buffer_is_reused_and_freed(tmp_path):
    _, paths = _files(tmp_path, sizes=[(9, 30, 26), (9, 30, 26), (7, 30, 26), (11, 30, 26)])
    with convert.keep_read_buffer() as kept:
        assert convert._kept_buffer() is kept
        convert.tiff_to_zarr(paths[0], str(tmp_path / "a.zarr"), usr_fcn=ERODE)
        mem = kept._bytes
        assert mem is not None and mem.size == 9 * 30 * 26 * 2
        convert.tiff_to_zarr(paths[1], str(tmp_path / "b.zarr"), usr_fcn=ERODE)
        convert.tiff_to_zarr(paths[2], str(tmp_path / "c.zarr"), usr_fcn=ERODE)
        assert kept._bytes is mem, "a file of the same size or smaller took new memory"
        with convert.keep_read_buffer() as inner:
            assert inner is kept
        assert kept._bytes is mem, "the inner block freed the outer one's buffer"
        convert.tiff_to_zarr(paths[3], str(tmp_path / "d.zarr"), usr_fcn=ERODE)
        assert kept._bytes.size == 11 * 30 * 26 * 2
    assert convert._kept_buffer() is None and kept._bytes is None


def test_the_reads_go_into_the_kept_buffer(tmp_path, monkeypatch):
    """the whole file, its input box and its z blocks"""
    from petakit5d import io

    seen = []
    real = io.read_tiff

    def spy(path, **kw):
        seen.append(kw.get("out"))
        return real(path, **kw)

    monkeypatch.setattr(io, "read_tiff", spy)
    _, paths = _files(tmp_path, sizes=[(9, 30, 26)])
    with convert.keep_read_buffer() as kept:
        convert.tiff_to_zarr(paths[0], str(tmp_path / "a.zarr"))
        convert.tiff_to_zarr(paths[0], str(tmp_path / "b.zarr"), input_bbox=(1, 2, 3, 8, 28, 24))
        convert.tiff_to_zarr(paths[0], str(tmp_path / "c.zarr"), block_size=(4, 16, 16),
                             read_whole_tiff=False)  # fmt: skip
    assert len(seen) == 5 and all(o is kept for o in seen)


def test_kept_buffer_belongs_to_its_thread():
    seen = []
    with convert.keep_read_buffer():
        t = threading.Thread(target=lambda: seen.append(convert._kept_buffer()))
        t.start()
        t.join()
    assert seen == [None]


def test_kept_buffer_is_freed_after_an_error(tmp_path):
    _, paths = _files(tmp_path, sizes=[(5, 20, 18)])
    with pytest.raises(FileNotFoundError), convert.keep_read_buffer() as kept:
        convert.tiff_to_zarr(paths[0], str(tmp_path / "a.zarr"))
        assert kept._bytes is not None
        convert.tiff_to_zarr(str(tmp_path / "missing.tif"), str(tmp_path / "b.zarr"))
    assert convert._kept_buffer() is None and kept._bytes is None


def _erosion_spy(monkeypatch):
    calls = []
    real = TF.erode_volume_by_2d_projection

    def spy(vol, esize, *, out=None, num_threads=None):
        calls.append({"in_place": out is vol, "kept": convert._kept_buffer() is not None and
                      convert._kept_buffer().holds(vol)})  # fmt: skip
        return real(vol, esize, out=out, num_threads=num_threads)

    monkeypatch.setattr(TF, "erode_volume_by_2d_projection", spy)
    return calls


def test_the_erosion_runs_in_place_on_the_conversions_own_volume(tmp_path, monkeypatch):
    calls = _erosion_spy(monkeypatch)
    vols, paths = _files(tmp_path, sizes=[(9, 30, 26)])
    convert.tiff_to_zarr(paths[0], str(tmp_path / "a.zarr"), usr_fcn=ERODE)
    with convert.keep_read_buffer():
        convert.tiff_to_zarr(paths[0], str(tmp_path / "b.zarr"), usr_fcn=ERODE)
        convert.tiff_to_zarr(paths[0], str(tmp_path / "c.zarr"), usr_fcn=ERODE, flip_zstack=True)
    # (the read; the kept buffer; the flipped copy)
    assert calls == [
        {"in_place": True, "kept": False},
        {"in_place": True, "kept": True},
        {"in_place": True, "kept": False},
    ]
    want = TF.erode_volume_by_2d_projection(vols[0], 2)
    assert np.array_equal(read_zarr(str(tmp_path / "a.zarr")), want)
    assert np.array_equal(read_zarr(str(tmp_path / "b.zarr")), want)


FRAME_OPTIONS = [
    {},
    {"flip_zstack": True},
    {"tile_out_bbox": (1, 0, 0, 8, 30, 26)},
    {"tile_out_bbox": (1, 2, 3, 8, 28, 24)},
]


@pytest.mark.parametrize(
    "kw", FRAME_OPTIONS, ids=["as_given", "flipped", "z_crop_view", "crop_copy"]
)
def test_a_frame_the_caller_gives_is_never_written(tmp_path, monkeypatch, kw):
    calls = _erosion_spy(monkeypatch)
    frame = _vol((9, 30, 26), np.uint16, "frame")
    before = frame.copy()
    with convert.keep_read_buffer():
        convert.tiff_to_zarr("", str(tmp_path / "f.zarr"), frame, usr_fcn=ERODE, **kw)
    assert np.array_equal(frame, before)
    # (the erosion writes in place only into the conversion's own copies)
    copied = bool(kw.get("flip_zstack")) or kw.get("tile_out_bbox", (0, 0, 0))[1:3] != (0, 0)
    assert calls == [{"in_place": copied, "kept": False}]
    view = frame[::-1] if kw.get("flip_zstack") else frame
    if "tile_out_bbox" in kw:
        z0, y0, x0, z1, y1, x1 = kw["tile_out_bbox"]
        view = view[z0:z1, y0:y1, x0:x1]
    assert np.array_equal(
        read_zarr(str(tmp_path / "f.zarr")), TF.erode_volume_by_2d_projection(view, 2)
    )


def test_python_functions_are_called_as_before(tmp_path):
    seen = []

    def fn(vol):
        seen.append(vol.shape)
        return vol + np.uint16(1)

    _, paths = _files(tmp_path, sizes=[(5, 20, 18)])
    with convert.keep_read_buffer():
        convert.tiff_to_zarr(paths[0], str(tmp_path / "a.zarr"), usr_fcn=fn)
    assert seen == [(5, 20, 18)]
    assert getattr(TF.parse_tile_function(ERODE), "in_place", False)
    assert not getattr(fn, "in_place", False)


def test_the_wrapper_keeps_the_buffer_in_local_mode_only(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from petakit5d.jobs import framework, slurm

    seen = []

    def fake(inputs, outputs, tasks, **kw):
        seen.append((kw["parse_cluster"], convert._kept_buffer() is not None))
        return [True] * len(list(outputs))

    monkeypatch.setattr(framework, "run_tasks", fake)
    _, paths = _files(tmp_path, sizes=[(5, 20, 18)])
    convert.tiff_to_zarr_dataset("", tiff_fullpaths=paths, parse_cluster=False)
    cluster = SimpleNamespace(parse_cluster=True)
    monkeypatch.setattr(slurm, "check_slurm_cluster", lambda *a, **k: cluster)
    convert.tiff_to_zarr_dataset("", tiff_fullpaths=paths, parse_cluster=True)
    assert seen == [(False, True), (True, False)]
    assert convert._kept_buffer() is None


def test_the_wrapper_converts_the_same_with_the_kept_buffer(tmp_path):
    from petakit5d.io import write_mat

    d = tmp_path / "data"
    d.mkdir()
    vols = []
    for k, s in enumerate(SIZES):
        vols.append(_vol(s, np.uint16, f"w{k}"))
        write_tiff(str(d / f"tile_{k}.tif"), vols[-1])
    fn = str(d / "fun.mat")
    write_mat(fn, {"usrFun": ERODE})
    convert.tiff_to_zarr_dataset(
        str(d), channel_patterns=["tile"], block_size=(4, 16, 16), process_fun_path=[fn],
        parse_cluster=False,
    )  # fmt: skip
    for k, v in enumerate(vols):
        got = read_zarr(str(d / "zarr" / f"tile_{k}.zarr"))
        assert np.array_equal(got, TF.erode_volume_by_2d_projection(v, 2))


# ---------------------------------------------------------------------------------------------
# the memory of kept buffers
# ---------------------------------------------------------------------------------------------


def test_byte_buffer_views_grow_and_tell_their_arrays():
    b = ByteBuffer()
    a = b.view((3, 4, 5), np.uint16)
    assert a.shape == (3, 4, 5) and a.dtype == np.uint16 and a.flags.c_contiguous
    assert a.flags.writeable and b.holds(a)
    mem = b._bytes
    c = b.view((2, 4, 5), np.float32)  # 160 bytes, more than the 120 there: it grows
    assert c.dtype == np.float32 and b._bytes is not mem and not b.holds(a)
    d = b.view((1, 2, 3), np.uint8)
    assert np.shares_memory(c, d) and b.holds(d) and not b.holds(np.zeros(3))
    b.release()
    assert b._bytes is None and not b.holds(d)


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="the mapping is Linux's")
def test_kept_memory_is_not_advised_without_the_numpy_variable(monkeypatch):
    import mmap

    def mapped(x):
        return isinstance(getattr(x.base, "obj", None), mmap.mmap)

    monkeypatch.delenv("NUMPY_MADVISE_HUGEPAGE", raising=False)
    a = kept_memory(1 << 20)
    assert mapped(a) and a.flags.writeable and a.size == 1 << 20 and a.dtype == np.uint8
    monkeypatch.setenv("NUMPY_MADVISE_HUGEPAGE", "0")
    b = kept_memory(1 << 20)
    assert not mapped(b) and b.size == 1 << 20
