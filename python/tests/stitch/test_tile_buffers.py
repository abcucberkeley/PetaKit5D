"""The tile preparation's memory: tile functions written as MATLAB code run in place on a volume
their caller gives up (``owned``), and the zarr tiles of a local loop are read into the kept
buffer (``io.convert.keep_read_buffer``), with the same results.

The tile function with ``owned`` gives the bits it gives without it and writes into the volume
only where it may (the erosion of the volume itself: not after an offset or a flat field, not
into a read-only or strided volume); ``process_zarr_tile`` gives the same zarr arrays with the
kept buffer (both zarr layouts, the options) and never writes a frame its caller gives;
``process_tiles`` keeps the buffer in local mode only.
"""

from __future__ import annotations

import contextlib

import numpy as np
import pytest

from petakit5d.io import convert, read_zarr
from petakit5d.stitch import tile_functions as TF
from petakit5d.stitch import tiles as T
from petakit5d.testing.synth import CounterRNG, write_zarr

ERODE = "@(x)erodeVolumeBy2DProjection(x,2)"


def _vol(shape=(9, 22, 26), dtype=np.uint16, name="v"):
    rng = CounterRNG(22, f"x11/tiles/{name}/{'x'.join(map(str, shape))}")
    v = rng.integers(1, 3000, shape).astype(dtype)
    if len(shape) == 3:
        v[:, :, :2] = 0
        v[:2, :, 9] = 0
        v[6:, :, 20:] = 0
    if np.dtype(dtype).kind == "f":
        v = v / np.asarray(7, dtype)
    return np.ascontiguousarray(v)


@pytest.mark.parametrize(
    "text",
    [ERODE, "@(x)erodeVolumeBy2DProjection(x,1)", "@(x)erodeVolumeBy2DProjection(x+5,2)",
     "@(x)x+100", "@(x)erodeVolumeBy2DProjection(x,0)"],
)  # fmt: skip
@pytest.mark.parametrize("dtype", [np.uint16, np.float32])
def test_owned_gives_the_same_bits(text, dtype):
    f = TF.parse_tile_function(text)
    assert f.in_place is True
    v = _vol(dtype=dtype)
    want = f(v.copy())
    w = v.copy()
    got = f(w, owned=True)
    assert got.dtype == want.dtype and np.array_equal(got, want)
    # only the erosion of the volume itself writes into it (of size 0: the volume as it is)
    in_place = "x+" not in text
    assert (got is w) == in_place
    if not in_place:
        assert np.array_equal(w, v), "the volume was written although the result is new"


def test_owned_writes_only_where_it_may():
    f = TF.parse_tile_function(ERODE)
    v = _vol()
    want = f(v)
    assert np.array_equal(v, _vol()), "without owned the volume was written"
    ro = v.copy()
    ro.flags.writeable = False
    got = f(ro, owned=True)
    assert got is not ro and np.array_equal(got, want)
    big = np.zeros((9, 22, 52), np.uint16)
    big[:, :, ::2] = v
    strided = big[:, :, ::2]
    got = f(strided, owned=True)
    assert got is not strided and np.array_equal(got, want)
    assert np.array_equal(big[:, :, ::2], v)
    im = _vol((22, 26), name="im")
    w = im.copy()
    got = f(w, owned=True)
    assert got is w and np.array_equal(got, f(im))


def _tile(tmp_path, name, v, layout):
    p = str(tmp_path / f"{name}.zarr")
    write_zarr(p, v, (4, 8, 8), layout)
    return p


OPTIONS = [
    {"usr_fcn": ERODE},
    {"usr_fcn": ERODE, "flip_zstack": True},
    {"usr_fcn": ERODE, "input_bbox": (1, 2, 3, 8, 20, 24)},
    {"usr_fcn": ERODE, "tile_out_bbox": (1, 0, 0, 7, 22, 26)},
    {"usr_fcn": "@(x)erodeVolumeBy2DProjection(x+1,1)"},
    {},
]


@pytest.mark.parametrize("layout", ["native", "matlab"])
@pytest.mark.parametrize(
    "kw", OPTIONS, ids=["erode", "flip", "input_bbox", "z_crop", "offset", "plain"]
)
def test_zarr_tiles_with_the_kept_buffer(tmp_path, layout, kw):
    vols = [_vol(name="a"), _vol((10, 23, 27), name="b"), _vol(name="c")]
    paths = [_tile(tmp_path, f"t{k}", v, layout) for k, v in enumerate(vols)]

    def run(tag, keep):
        out = []
        with convert.keep_read_buffer() if keep else contextlib.nullcontext():
            for k, p in enumerate(paths):
                z = str(tmp_path / f"{tag}{k}.zarr")
                T.process_zarr_tile([p], z, block_size=(4, 8, 8), **kw)
                out.append(read_zarr(z))
        return out

    want = run("new", False)
    for g, w in zip(run("kept", True), want, strict=True):
        assert g.dtype == w.dtype and np.array_equal(g, w)
    if kw == {"usr_fcn": ERODE}:
        for g, v in zip(want, vols, strict=True):
            assert np.array_equal(g, TF.erode_volume_by_2d_projection(v, 2))


def test_zarr_tile_reads_into_the_kept_buffer(tmp_path, monkeypatch):
    calls = []
    real = TF.erode_volume_by_2d_projection

    def spy(vol, esize, *, out=None, num_threads=None):
        kept = convert._kept_buffer()
        calls.append((out is vol, kept is not None and kept.holds(vol)))
        return real(vol, esize, out=out, num_threads=num_threads)

    monkeypatch.setattr(TF, "erode_volume_by_2d_projection", spy)
    p = _tile(tmp_path, "t", _vol(), "native")
    T.process_zarr_tile([p], str(tmp_path / "a.zarr"), usr_fcn=ERODE)
    with convert.keep_read_buffer():
        T.process_zarr_tile([p], str(tmp_path / "b.zarr"), usr_fcn=ERODE)
    assert calls == [(True, False), (True, True)]


@pytest.mark.parametrize("kw", [{}, {"tile_out_bbox": (1, 0, 0, 7, 22, 26)}, {"flip_zstack": True}])
def test_a_frame_the_caller_gives_is_never_written(tmp_path, kw):
    frame = _vol(name="frame")
    before = frame.copy()
    with convert.keep_read_buffer():
        T.process_zarr_tile("", str(tmp_path / "f.zarr"), frame, usr_fcn=ERODE, **kw)
    assert np.array_equal(frame, before)
    view = frame[::-1] if kw.get("flip_zstack") else frame
    if "tile_out_bbox" in kw:
        view = view[1:7]
    assert np.array_equal(
        read_zarr(str(tmp_path / "f.zarr")), TF.erode_volume_by_2d_projection(view, 2)
    )


def test_process_tiles_keeps_the_buffer_in_local_mode_only(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from petakit5d.io import write_mat
    from petakit5d.jobs import framework, slurm

    seen = []

    def fake(inputs, outputs, tasks, **kw):
        seen.append((kw["parse_cluster"], convert._kept_buffer() is not None))
        return [True] * len(list(outputs))

    monkeypatch.setattr(framework, "run_tasks", fake)
    d = tmp_path / "d"
    d.mkdir()
    p = _tile(d, "tile_CamA", _vol(), "native")
    fn = str(d / "fun.mat")
    write_mat(fn, {"usrFun": ERODE})
    kw = {"zarr_file": True, "result_dir_name": "proc", "process_fun_path": [[fn]]}
    T.process_tiles([p], parse_cluster=False, **kw)
    cluster = SimpleNamespace(parse_cluster=True)
    monkeypatch.setattr(slurm, "check_slurm_cluster", lambda *a, **k: cluster)
    T.process_tiles([p], parse_cluster=True, **kw)
    assert seen == [(False, True), (True, False)]
    assert convert._kept_buffer() is None
