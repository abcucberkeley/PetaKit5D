"""Tile preparation of the stitching (petakit5d.stitch.tiles, stitch.tile_functions; task S02).

Oracle cases against MATLAB's oracle tree:

- ``fun.*``: tile functions written as MATLAB code, run by MATLAB's ``str2func`` and by the
  port's parser: the flat-field correction (processFFCorrectionFrame) with 4 to 8 arguments,
  its composition with the edge erosion in both orders (the user's form and the demos'
  ``sprintf`` form with ``%d`` and ``%0.10d`` numbers), ``+`` and ``-`` of constants on uint16,
  int8 and single volumes (rounding, saturation, NaN payloads), an image, a non-integer erosion
  size;
- ``tile.*``: stitch_process_zarr_tile on a zarr tile in MATLAB's layout with each option
  (flip, inputBbox, tileOutBbox and all three, resampling, the compressor, the default chunks,
  both standard function forms and the erosion with an offset, a single tile, an image tile and
  its resampling error, an existing result, a frame given instead of the tile);
- ``tiles.*``: stitch_process_tiles on zarr tiles (flipped tiles, boxes per tile location,
  resampling, function files as .mat and .txt, keyword alternatives, a function file of no
  function, none at all) and on tiff tiles (through the tiff wrapper, with the same function
  forms and keywords).

The function files hold absolute paths, so each side writes them in its own work folder before
the call and removes them after it. Zarr results are compared by voxels and metadata after the
layout mapping, the tiff wrapper's parameter files by content.

Then StitchHelperMexTest's edge erosion and the tile-conversion parts of StitchingTest in
intent, the parser's forms and refusals, the tasks' settings and memory (framework spy) and the
publish rule.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from petakit5d.core.dtypes import NAN
from petakit5d.io import read_zarr, write_mat
from petakit5d.stitch import tile_functions as TF
from petakit5d.stitch import tiles as T
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, InPath, Tree, Value, Vol
from petakit5d.testing.synth import CounterRNG, write_tiff, write_zarr

SHAPE = (9, 22, 26)  # (z, y, x): not cubic; MATLAB 22 x 26 x 9
CHUNKS_M = [8, 8, 4]  # MATLAB's [y x z]: remainders on every axis
NAN32_B = np.array([0x7FC00001], dtype=np.uint32).view(np.float32)[0]


def volume(kind: str = "uint16", name: str = "v", shape=SHAPE, holes: bool = True) -> np.ndarray:
    """deterministic positive values with columns of zeros near the edges (the projection is
    not full); single: NaN of two payloads too"""
    rng = CounterRNG(2, f"s02/{name}/{kind}")
    v = rng.integers(1, 3000, shape).astype(np.uint16)
    if holes and len(shape) == 3:
        v[:, :, :2] = 0
        v[:2, :, 9] = 0
        v[6:, :, 20:] = 0
    if kind == "int8":
        return np.ascontiguousarray((v % 200).astype(np.int16) - 100).astype(np.int8)
    if kind == "single":
        v = v.astype(np.float32) / np.float32(7)
        if len(shape) == 3:
            v[4, 5, 6] = NAN
            v[:, 3, 15] = NAN32_B
            v[1, 2, 3] = -0.0
    return np.ascontiguousarray(v)


def _flat_fields(folder: str, plane=SHAPE[1:]) -> None:
    rng = CounterRNG(3, "s02/ff")
    write_tiff(os.path.join(folder, "ff.tif"), rng.integers(500, 1500, plane).astype(np.uint16))
    write_tiff(os.path.join(folder, "bg.tif"), rng.integers(90, 110, plane).astype(np.uint16))


def _m(v: np.ndarray):
    from petakit5d.core.arrays import to_matlab

    return Value(lambda: to_matlab(v))


# ---------------------------------------------------------------------------------------------
# fun.*: tile functions as MATLAB code
# ---------------------------------------------------------------------------------------------

FF = "processFFCorrectionFrame(x,'d/ff.tif','d/bg.tif',1,0.4,false)"
USER = f"@(x)erodeVolumeBy2DProjection({FF},1)"  # the user's flat field and erosion
# the demos' sprintf("...,%d,%d,%s),%0.10d)", ff, bg, 1, 0.4, string(false), 2)
DEMO = (
    "@(x)erodeVolumeBy2DProjection(processFFCorrectionFrame(x,'d/ff.tif','d/bg.tif',1,"
    "4.000000e-01,false),0000000002)"
)


def _py_fun(v, text, _d):
    return TF.parse_tile_function(text)(v)


def fun_case(cid, text, vol, error=False):
    return Case(
        f"fun.{cid}",
        code="f = str2func(in{2}); out = {f(in{1})};",
        args=[Vol(In("v")), text, InPath("d")],
        python=_py_fun,
        out_roles=("vol",),
        inputs={"v": _m(vol), "d": Tree(_flat_fields, "d")},
        expect_error=error,
        py_error=(TF.TileFunctionError, ValueError),
    )


FUN_CASES = [
    fun_case("user", USER, volume()),
    fun_case("demo_sprintf", DEMO, volume()),
    fun_case("ff_only", f"@(x){FF}", volume()),
    fun_case(
        "ff_defaults_strings",
        '@( vol ) processFFCorrectionFrame( vol , "d/ff.tif" , "d/bg.tif" , [] )',
        volume("single"),
    ),
    fun_case(
        "ff_eight_args",
        "@(x)processFFCorrectionFrame(x,'d/ff.tif','d/bg.tif',100,0.2,true,true,false)",
        volume(),
    ),
    fun_case(
        "erode_then_ff",
        "@(x)processFFCorrectionFrame(erodeVolumeBy2DProjection(x,2),'d/ff.tif','d/bg.tif',1)",
        volume("uint16", "e"),
    ),
    fun_case("user_image", USER, np.ascontiguousarray(volume("uint16", "im")[3])),
    fun_case("user_single", USER, volume("single")),
    fun_case("offset_saturates", "@(x)erodeVolumeBy2DProjection(x+64000,2)", volume()),
    fun_case("plus", "@(x)x+100", volume()),
    fun_case("minus_half", "@(x) x - 1.5", volume()),
    fun_case("number_minus_x", "@(x)2.5-x", volume("int8")),
    fun_case("plus_tenth_single", "@(x)x+0.1", volume("single")),
    fun_case("signs_parentheses", "@(x)(-3+x)+(+-2)", volume("single")),
    fun_case("int8_minus_half", "@(x)x-0.5", volume("int8")),
    fun_case("esize_not_integer", "@(x)erodeVolumeBy2DProjection(x,1.5)", volume(), error=True),
]


@oracle.parametrize(FUN_CASES)
def test_functions(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# tile.*: stitch_process_zarr_tile
# ---------------------------------------------------------------------------------------------


def _tile_tree(kind="uint16", shape=SHAPE, out=False):
    def make(folder):
        write_zarr(
            os.path.join(folder, "tile.zarr"), volume(kind, "tile", shape), (4, 8, 8), "matlab"
        )
        _flat_fields(folder, shape[1:])
        if out:  # an existing result
            write_zarr(os.path.join(folder, "out.zarr"), np.ones((2, 3, 4), np.uint16), (2, 3, 4))

    return make


def _py_tile(d, *nv):
    from petakit5d import matlab

    matlab.stitch_process_zarr_tile([d + "tile.zarr"], d + "out.zarr", None, *nv)
    return []


def _py_tile_frame(d, *nv):
    from petakit5d import matlab

    matlab.stitch_process_zarr_tile(
        [d + "tile.zarr"], d + "out.zarr", np.ones((3, 3, 3), np.uint16)
    )
    return []


def tile_case(cid, *nv, kind="uint16", shape=SHAPE, out=False, error=False, resampled=False):
    return Case(
        f"tile.{cid}",
        code=(
            "stitch_process_zarr_tile({[in{1}, 'tile.zarr']}, [in{1}, 'out.zarr'], [], "
            "in{2 : end}); out = {};"
        ),
        nargout=0,
        args=[InPath("d"), *nv],
        python=_needs_resize(_py_tile) if resampled else _py_tile,
        inputs={"d": Tree(_tile_tree(kind, shape, out), "d")},
        expect_error=error,
        py_error=(T.TileError, ValueError),
    )


def _needs_resize(f):
    def g(*a):
        pytest.importorskip("petakit5d._ext.resize", reason="the resize kernels are not built")
        return f(*a)

    return g


BS = ["blockSize", CHUNKS_M]
TILE_CASES = [
    tile_case("plain", *BS),
    tile_case("default_chunks"),
    tile_case("flip", *BS, "flipZstack", True),
    tile_case("input_bbox", *BS, "inputBbox", [3, 4, 2, 20, 24, 8]),
    tile_case("tile_out_bbox", *BS, "tileOutBbox", [2, 3, 1, 15, 20, 6]),
    tile_case(
        "flip_both_boxes_zstd",
        *BS,
        "flipZstack",
        True,
        "inputBbox",
        [3, 4, 2, 20, 24, 8],
        "tileOutBbox",
        [2, 3, 2, 15, 20, 6],
        "compressor",
        "zstd",
    ),
    tile_case("erosion", *BS, "usrFcn", "@(x)erodeVolumeBy2DProjection(x,2)"),
    tile_case(
        "erosion_offset_single",
        *BS,
        "usrFcn",
        "@(x)erodeVolumeBy2DProjection(x+1,1)",
        kind="single",
    ),
    tile_case("user", *BS, "usrFcn", USER),
    tile_case("demo_sprintf_flip", *BS, "flipZstack", True, "usrFcn", DEMO),
    tile_case(
        "everything",
        *BS,
        "flipZstack",
        True,
        "inputBbox",
        [2, 2, 1, 21, 25, 9],
        "tileOutBbox",
        [1, 1, 1, 18, 22, 8],
        "resample",
        [1.5, 1.5, 2],
        "usrFcn",
        "@(x)erodeVolumeBy2DProjection(processFFCorrectionFrame(x,'d/ff.tif','d/bg.tif',1,"
        "0.4,false),1)",
        resampled=True,
    ),
    tile_case("resample_xy_z", *BS, "resample", [2, 0.5], resampled=True),
    tile_case("image", *BS, "usrFcn", USER, shape=(1, 22, 26)),
    tile_case("image_resample", *BS, "resample", [2, 2, 1], shape=(1, 22, 26), error=True),
    tile_case("existing_result", *BS, "usrFcn", USER, out=True),
    # (a frame: MATLAB failed on an undefined sz until 0fe22a5)
    Case(
        "tile.frame",
        code=(
            "stitch_process_zarr_tile({[in{1}, 'tile.zarr']}, [in{1}, 'out.zarr'], "
            "ones(3, 3, 3, 'uint16')); out = {};"
        ),
        nargout=0,
        args=[InPath("d")],
        python=_py_tile_frame,
        inputs={"d": Tree(_tile_tree(), "d")},
        py_error=(T.TileError,),
    ),
]


@oracle.parametrize(TILE_CASES)
def test_zarr_tile(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# tiles.*: stitch_process_tiles
# ---------------------------------------------------------------------------------------------

NAMES = ("tile_CamA_0", "tile_CamA_1", "tile_CamB_0", "tile_CamB_1")
# function files: (name, MATLAB code with <D> for the data folder, keywords or None)
F_USER = (
    "fun.mat",
    "@(x)erodeVolumeBy2DProjection(processFFCorrectionFrame(x,'<D>ff.tif','<D>bg.tif',1,0.4,"
    "false),1)",
    None,
)
F_TXT = ("fun.txt", "@(x)erodeVolumeBy2DProjection(x+1,2)", None)
F_KEY = ("fun_b.mat", F_USER[1], ["CamB"])
F_NONE = ("none.mat", "", None)


def _tiles_tree(fmt):
    def make(folder):
        for k, n in enumerate(NAMES):
            v = volume("uint16", f"tiles/{k}")
            if fmt == "zarr":
                write_zarr(os.path.join(folder, n + ".zarr"), v, (4, 8, 8), "matlab")
            else:
                write_tiff(os.path.join(folder, n + ".tif"), v)
        _flat_fields(folder)

    return make


# MATLAB: write the function files (absolute paths), call, remove them
_TILES_CODE = (
    "d = in{1}; specs = in{3}; pf = cell(size(specs)); "
    "for k = 1 : numel(specs), s = specs{k}; if isempty(s), pf{k} = ''; continue; end; "
    "usrFun = strrep(s{2}, '<D>', d); fn = [d, s{1}]; "
    "if endsWith(fn, '.txt'), fid = fopen(fn, 'w'); fprintf(fid, '%s\\n', usrFun); fclose(fid); "
    "elseif isempty(s{3}), save('-v7.3', fn, 'usrFun'); "
    "else, keywords = s{3}; save('-v7.3', fn, 'usrFun', 'keywords'); end; pf{k} = fn; end; "
    "stitch_process_tiles(cellfun(@(x) [d, x], in{2}, 'unif', 0), 'processFunPath', pf, "
    "in{4 : end}); "
    "for k = 1 : numel(specs), if ~isempty(specs{k}), delete([d, specs{k}{1}]); end; end; "
    "out = {};"
)


def _write_function_files(d: str, specs) -> list[str]:
    paths = []
    for s in specs:
        if not s:
            paths.append("")
            continue
        name, text, kw = s
        text = text.replace("<D>", d)
        fn = d + name
        if name.endswith(".txt"):
            with open(fn, "w", newline="\n") as f:
                f.write(text + "\n")
        else:
            write_mat(fn, {"usrFun": text} if kw is None else {"usrFun": text, "keywords": kw})
        paths.append(fn)
    return paths


def _tiles_kwargs(nv) -> dict:
    """stitch_process_tiles' name-value pairs (MATLAB's values) as process_tiles' keywords"""
    from petakit5d.core.axes import from_yxz
    from petakit5d.core.bbox import box_from_matlab

    kw: dict = {}
    for name, v in zip(nv[::2], nv[1::2], strict=True):
        if name in ("InputBbox", "tileOutBbox"):
            a = np.asarray(v, dtype=np.float64)
            v = (
                tuple(box_from_matlab(r) for r in a)
                if a.ndim == 2 and a.shape[0] > 1
                else box_from_matlab(a.ravel())
            )
            name = "input_bbox" if name == "InputBbox" else "tile_out_bbox"
        elif name == "blockSize":
            v = tuple(int(x) for x in from_yxz(v))
        elif name == "resample":
            v = tuple(from_yxz(v)) if len(v) == 3 else tuple(v[::-1])
        from petakit5d.params import snake_name

        kw[snake_name(name)] = v
    return kw


def _py_tiles(d, names, specs, *nv):
    paths = _write_function_files(d, specs or [])
    kw = _tiles_kwargs(nv)
    try:
        T.process_tiles([d + n for n in names], process_fun_path=[paths] if specs else None, **kw)
    finally:
        for s in specs or []:
            if s:
                os.remove(d + s[0])
    return []


def tiles_case(cid, fmt, specs, *nv, names=NAMES):
    ext = ".zarr" if fmt == "zarr" else ".tif"
    return Case(
        f"tiles.{cid}",
        code=_TILES_CODE,
        nargout=0,
        args=[
            InPath("d"),
            [n + ext for n in names],
            [list(s) if s else None for s in specs] if specs else None,
            *nv,
            "parseCluster",
            False,
        ],
        python=_py_tiles,
        inputs={"d": Tree(_tiles_tree(fmt), "d")},
    )


Z = ["zarrFile", True, "resultDirName", "proc", *BS]
TILES_CASES = [
    tiles_case(
        "zarr_user_flipped_input_bbox",
        "zarr",
        [F_USER],
        *Z,
        "flippedTile",
        [True, False, False, True],
        "InputBbox",
        [2, 3, 1, 20, 24, 8],
    ),
    tiles_case(
        "zarr_txt_boxes_per_location_resample",
        "zarr",
        [F_TXT],
        *Z,
        "locIds",
        [1, 1, 2, 2],
        "tileOutBbox",
        np.array([[1, 1, 1, 22, 26, 9], [3, 2, 2, 18, 25, 7]], dtype=np.float64),
        "resample",
        [1.5, 1.5, 2],
    ),
    tiles_case("zarr_keywords", "zarr", [F_KEY, F_TXT], *Z),
    tiles_case("zarr_function_file_of_no_function", "zarr", [F_NONE], *Z),
    tiles_case("zarr_no_function_files", "zarr", None, *Z),
    tiles_case(
        "tif_user_flipped_boxes",
        "tif",
        [F_USER],
        "resultDirName",
        "zarr",
        *BS,
        "flippedTile",
        [False, True, True, False],
        "InputBbox",
        [2, 3, 1, 20, 24, 8],
        "tileOutBbox",
        [1, 1, 1, 15, 20, 6],
    ),
    tiles_case("tif_keywords_txt", "tif", [F_KEY, F_TXT], "resultDirName", "zarr", *BS),
]
for _c in TILES_CASES:
    if any(isinstance(a, str) and a == "resample" for a in _c.args):
        _c.python = _needs_resize(_c.python)


@oracle.parametrize(TILES_CASES)
def test_tiles(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# StitchHelperMexTest: the edge erosion
# ---------------------------------------------------------------------------------------------


def test_edge_erosion_zeroes_esize_rows_at_both_ends_of_y():
    # edgeErosionWithoutTheMexMatches: esize rows at both ends of y are zero, the next ones
    # are not (MATLAB vol(:, 1 : 3, :) = 0 are the first three x columns)
    rng = CounterRNG(5, "s02/helper")
    vol = rng.integers(1, 1000, (20, 30, 40)).astype(np.uint16)
    vol[:, :, :3] = 0
    for e in (1, 3, 15):
        out = TF.erode_volume_by_2d_projection(vol, e)
        assert np.count_nonzero(out[:, :e]) == 0 and np.count_nonzero(out[:, 30 - e :]) == 0
        if e < 15:
            assert np.count_nonzero(out[:, e]) > 0 and np.count_nonzero(out[:, 29 - e]) > 0


@pytest.mark.parametrize("kind", ["uint16", "single"])
def test_erosion_in_place_is_the_erosion(kind):
    v = volume(kind)
    want = TF.erode_volume_by_2d_projection(v, 2)
    w = v.copy()
    got = TF.erode_volume_by_2d_projection(w, 2, out=w)
    assert got is w and got.tobytes() == want.tobytes()
    assert TF.erode_volume_by_2d_projection(v, 0) is v


def test_erosion_rejects_an_out_of_another_shape():
    with pytest.raises(ValueError, match="out"):
        TF.erode_volume_by_2d_projection(volume(), 1, out=np.empty((2, 2, 2), np.uint16))


# ---------------------------------------------------------------------------------------------
# the parser
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "text",
    [
        "@(x)processFFCorrectionFrame(x,'a.tif','b.tif',1,0.4,false)",
        "@(x)erodeVolumeBy2DProjection(processFFCorrectionFrame(x,'a.tif','b.tif',1,0.4,false),1)",
        "@ ( im ) erodeVolumeBy2DProjection ( im + 2 , 3 )",
        "@(x)x+1000",
        "@(x)((x))",
        "@(x)processFFCorrectionFrame(x,\"it''s.tif\",'b''.tif',[],.4,0,1,true)",
    ],
)
def test_standard_forms_are_read(text):
    assert callable(TF.parse_tile_function(text))


@pytest.mark.parametrize(
    "text",
    [
        "@(x) single(x) / 3",
        "@(x) x(:, end : -1 : 1, :)",
        "@(x) x * 2",
        "@(x) x + x",
        "@(x) -x",
        "@(x) x'",
        "@(x) 3",
        "@(x) imgaussfilt(x, 2)",
        "@(x) erodeVolumeBy2DProjection(2, x)",
        "@(x) processFFCorrectionFrame(x, ff, 'b.tif', 1)",
        "@(x) x + true",
        "x + 1",
        "@(x) x + 1; system('ls')",
        "@(x) erodeVolumeBy2DProjection(x, 1) % comment",
        "",
    ],
)
def test_other_code_is_refused(text):
    with pytest.raises(NotImplementedError) as err:
        TF.parse_tile_function(text)
    assert err.value.identifier == TF.UNSUPPORTED_ID


def test_wrong_argument_counts_are_matlabs_errors():
    with pytest.raises(TF.TileFunctionError) as err:
        TF.parse_tile_function("@(x)processFFCorrectionFrame(x,'a.tif','b.tif')")
    assert err.value.identifier == "MATLAB:minrhs"
    with pytest.raises(TF.TileFunctionError):
        TF.parse_tile_function("@(x)erodeVolumeBy2DProjection(x,1,2)")
    f = TF.parse_tile_function("@(x)erodeVolumeBy2DProjection(x,'a')")
    with pytest.raises(TF.TileFunctionError, match="esize"):
        f(volume())


def test_the_input_is_not_changed():
    v = volume()
    keep = v.copy()
    for text in ("@(x)erodeVolumeBy2DProjection(x,2)", "@(x)erodeVolumeBy2DProjection(x+1,2)"):
        TF.parse_tile_function(text)(v)
        assert v.tobytes() == keep.tobytes()


def test_the_flat_field_result_is_eroded_in_place(tmp_path, monkeypatch):
    # (the erosion of the composition writes into the correction's new array)
    monkeypatch.chdir(tmp_path)
    os.makedirs("d")
    _flat_fields("d")
    v = volume()
    seen = []
    real = TF.erode_volume_by_2d_projection

    def spy(vol, esize, *, out=None, num_threads=None):
        seen.append(out is vol)
        return real(vol, esize, out=out, num_threads=num_threads)

    monkeypatch.setattr(TF, "erode_volume_by_2d_projection", spy)
    TF.parse_tile_function(USER)(v)
    TF.parse_tile_function("@(x)erodeVolumeBy2DProjection(x,1)")(v)
    assert seen == [True, False]


@pytest.mark.parametrize("kind", ["uint16", "int8", "single"])
@pytest.mark.parametrize("text", ["@(x)x+7", "@(x)x-2.5", "@(x)1.25-x", "@(x)x+0.1"])
def test_arithmetic_is_matlabs(kind, text):
    from petakit5d.core import dtypes

    v = volume(kind, shape=(40, 30, 50))  # several threads' ranges
    op, d, first = {
        "@(x)x+7": ("+", 7.0, False),
        "@(x)x-2.5": ("-", 2.5, False),
        "@(x)1.25-x": ("-", 1.25, True),
        "@(x)x+0.1": ("+", 0.1, False),
    }[text]
    f = dtypes.plus if op == "+" else dtypes.minus
    want = np.asarray(f(d, v) if first else f(v, d))
    got = TF.parse_tile_function(text)(v)
    assert got.dtype == v.dtype and got.tobytes() == want.tobytes()


def test_the_matlab_code_is_kept_for_display():
    assert TF.parse_tile_function(" @(x)x+1 ").matlab_code == " @(x)x+1 "


def test_import_paths_still_name_python_functions():
    from petakit5d.io.convert import user_function

    assert user_function("numpy:negative") is np.negative
    assert user_function("@(x)x+1") is not None


# ---------------------------------------------------------------------------------------------
# stitch_process_tiles: tasks, memory, errors (framework spy)
# ---------------------------------------------------------------------------------------------


@pytest.fixture
def framework_spy(monkeypatch):
    """records run_tasks calls (MATLAB's FrameworkSpyFixture); nothing runs"""
    from petakit5d.jobs import framework

    calls = []

    def fake(inputs, outputs, tasks, **kw):
        calls.append({"inputs": list(inputs), "outputs": list(outputs), "tasks": list(tasks), **kw})
        return [False] * len(list(outputs))

    monkeypatch.setattr(framework, "run_tasks", fake)
    return calls


def _zarr_tiles(root, shapes):
    paths = []
    for k, s in enumerate(shapes):
        p = os.path.join(root, f"t{k}_CamA.zarr" if k % 2 == 0 else f"t{k}_CamB.zarr")
        write_zarr(p, np.zeros(s, np.uint16), (4, 8, 8))
        paths.append(p)
    return paths


@pytest.mark.parametrize("big_data", [True, False])
def test_tasks_memory_and_rounds(tmp_path, framework_spy, big_data):
    paths = _zarr_tiles(str(tmp_path), [(9, 22, 26), (5, 10, 12)])
    fn = tmp_path / "f.txt"
    fn.write_text("@(x)x+1\n")
    T.process_tiles(
        paths,
        zarr_file=True,
        result_dir_name="proc",
        process_fun_path=[[str(fn)]],
        resample=(1 / 3, 1, 1),
        parse_cluster=False,
        big_data=big_data,
        cpus_per_task=3,
    )
    one, two = framework_spy
    # the last tile's voxels as single, times 2.5 (and 2 without big data)
    mem = 5 * 10 * 12 * 4 / 1024**3 * 2.5 * (1 if big_data else 2)
    assert one["mem_allocate"] == pytest.approx(mem) and two["mem_allocate"] == pytest.approx(
        2 * mem
    )
    assert (one["cpus_per_task"], two["cpus_per_task"]) == (3, 6)
    assert one["max_trial_num"] == 2 and one["job_log_dir"] == "../job_logs"
    assert one["outputs"] == [f"{tmp_path}/proc/t0_CamA.zarr/", f"{tmp_path}/proc/t1_CamB.zarr/"]
    kw = one["tasks"][0].kwargs
    assert kw["usr_fcn"] == "@(x)x+1" and kw["compressor"] == "zstd"
    assert kw["resample"] == [0.333333333333333, 1.0, 1.0]  # (mat2str's 15 digits)


def test_keywords_pick_the_function(tmp_path, framework_spy):
    paths = _zarr_tiles(str(tmp_path), [(3, 4, 5)] * 2)
    a, b = tmp_path / "a.mat", tmp_path / "b.txt"
    write_mat(a, {"usrFun": "@(x)x+1", "keywords": ["CamB"]})
    b.write_text("@(x)x+2\n")
    T.process_tiles(
        paths, zarr_file=True, result_dir_name="p", process_fun_path=[[str(a), str(b)]],
        parse_cluster=False,
    )  # fmt: skip
    assert [t.kwargs["usr_fcn"] for t in framework_spy[0]["tasks"]] == ["@(x)x+2", "@(x)x+1"]


def test_nothing_to_do_without_function_files_or_boxes(tmp_path, framework_spy):
    paths = _zarr_tiles(str(tmp_path), [(3, 4, 5)])
    T.process_tiles(paths, zarr_file=True, result_dir_name="p", parse_cluster=False)
    T.process_tiles(
        paths, zarr_file=True, result_dir_name="", process_fun_path=[[""]], parse_cluster=False
    )
    assert framework_spy == []


def test_boxes_without_function_files_fail_as_in_matlab(tmp_path, framework_spy):
    # (usrFcn_strs_c{1} of an empty cell)
    paths = _zarr_tiles(str(tmp_path), [(3, 4, 5)])
    with pytest.raises(T.TileError) as err:
        T.process_tiles(
            paths, zarr_file=True, result_dir_name="p", input_bbox=(0, 0, 0, 2, 2, 2),
            parse_cluster=False,
        )  # fmt: skip
    assert err.value.identifier == "MATLAB:badsubscript"


def test_settings_matlab_rejects(tmp_path, framework_spy):
    paths = _zarr_tiles(str(tmp_path), [(3, 4, 5)] * 3)
    common = {"zarr_file": True, "result_dir_name": "p", "parse_cluster": False}
    with pytest.raises(T.TileError) as err:
        T.process_tiles(paths, process_fun_path=str(tmp_path / "f.mat"), **common)
    assert err.value.identifier == "MATLAB:cellRefFromNonCell"
    with pytest.raises(T.TileError, match=r"InputBbox \(3\) does not match .* \(2\)"):
        T.process_tiles(
            paths, loc_ids=(1, 1, 2), input_bbox=((0, 0, 0, 1, 1, 1),) * 3,
            tile_out_bbox=((0, 0, 0, 1, 1, 1),) * 2, process_fun_path=[[""]], **common,
        )  # fmt: skip
    with pytest.raises(T.TileError, match="channel patterns"):
        T.process_tiles(paths, channel_patterns=("CamC",), process_fun_path=[[""]], **common)


def test_tiff_tiles_go_to_the_tiff_wrapper(tmp_path, monkeypatch):
    from petakit5d.io import convert

    seen = {}
    monkeypatch.setattr(convert, "tiff_to_zarr_dataset", lambda *a, **k: seen.update(k))
    fn = str(tmp_path / "f.txt")
    with open(fn, "w") as f:
        f.write("@(x)x+1\n")
    T.process_tiles(
        [str(tmp_path / "a.tif")], result_dir_name="zarr", resample=(2, 1, 1),
        process_fun_path=[[fn]], cpus_per_task=2, parse_cluster=False,
    )  # fmt: skip
    assert seen["resample_factor"] == (2, 1, 1) and seen["channel_patterns"] == (".",)
    assert seen["cpus_per_task"] == 2 and seen["process_fun_path"] == [[fn]]


# ---------------------------------------------------------------------------------------------
# stitch_process_zarr_tile
# ---------------------------------------------------------------------------------------------


def test_zarr_tile_order_of_steps(tmp_path):
    # read the box, flip, crop, function, cast back (single to uint16 rounds)
    v = volume()
    write_zarr(str(tmp_path / "t.zarr"), v, (4, 8, 8), "matlab")
    T.process_zarr_tile(
        str(tmp_path / "t.zarr"), str(tmp_path / "o.zarr"), block_size=(4, 8, 8),
        flip_zstack=True, input_bbox=(1, 2, 3, 8, 20, 24), tile_out_bbox=(0, 1, 1, 6, 15, 18),
        usr_fcn=lambda a: a.astype(np.float32) / np.float32(2),
    )  # fmt: skip
    from petakit5d.core.dtypes import cast

    want = v[1:8, 2:20, 3:24][::-1][0:6, 1:15, 1:18].astype(np.float32) / np.float32(2)
    got = read_zarr(str(tmp_path / "o.zarr"))
    np.testing.assert_array_equal(got, cast(want, "uint16"))


def test_an_image_tile_is_an_image_for_the_function(tmp_path):
    write_zarr(str(tmp_path / "t.zarr"), volume(shape=(1, 6, 7), holes=False), (1, 4, 4), "matlab")
    shapes = []

    def f(a):
        shapes.append(a.shape)
        return a

    T.process_zarr_tile(str(tmp_path / "t.zarr"), str(tmp_path / "o.zarr"), usr_fcn=f)
    assert shapes == [(6, 7)]
    assert read_zarr(str(tmp_path / "o.zarr")).shape == (6, 7)


def test_an_existing_result_is_kept(tmp_path):
    write_zarr(str(tmp_path / "t.zarr"), volume(), (4, 8, 8))
    write_zarr(str(tmp_path / "o.zarr"), np.ones((2, 2, 2), np.uint16), (2, 2, 2))
    T.process_zarr_tile(str(tmp_path / "t.zarr"), str(tmp_path / "o.zarr"))
    assert read_zarr(str(tmp_path / "o.zarr")).shape == (2, 2, 2)
    T.process_zarr_tile(str(tmp_path / "t.zarr"), str(tmp_path / "o.zarr"), overwrite=True)
    assert read_zarr(str(tmp_path / "o.zarr")).shape == SHAPE


def test_stitch_zarr_tile_of_a_frame_or_without_the_frame_input(tmp_path):
    """(TiffZarrConversionTest/stitchZarrTileOfAFrameOrWithoutTheFrameInput) a given frame is
    processed instead of the tile (MATLAB failed on an undefined sz until 0fe22a5), and
    name-value pairs right after the output name are names (MATLAB took the first one as the
    frame); a NaN box end is refused by read_zarr"""
    from petakit5d import matlab
    from petakit5d.core.bbox import NOT_FINITE_ID
    from petakit5d.io import get_zarr_info

    im = CounterRNG(0, "sw8/frame").integers(0, 60000, (12, 40, 36), dtype=np.uint16)
    zr = str(tmp_path / "im.zarr")
    write_zarr(zr, im, (4, 16, 16))
    b = [2, 3, 1, 30, 33, 8]  # MATLAB's box; (0, 1, 2, 8, 30, 33)
    crop = im[0:8, 1:30, 2:33]
    matlab.stitch_process_zarr_tile(zr, str(tmp_path / "given.zarr"), None, "inputBbox", b)
    np.testing.assert_array_equal(read_zarr(str(tmp_path / "given.zarr")), crop)
    matlab.stitch_process_zarr_tile(
        zr, str(tmp_path / "names.zarr"), "inputBbox", b, "blockSize", [8, 8, 4]
    )
    np.testing.assert_array_equal(read_zarr(str(tmp_path / "names.zarr")), crop)
    assert tuple(get_zarr_info(str(tmp_path / "names.zarr")).block_size) == (4, 8, 8)
    T.process_zarr_tile("", str(tmp_path / "frame.zarr"), im, flip_zstack=True)
    np.testing.assert_array_equal(read_zarr(str(tmp_path / "frame.zarr")), im[::-1])
    with pytest.raises(ValueError) as err:
        T.process_zarr_tile(zr, str(tmp_path / "nan.zarr"), input_bbox=(0, 1, 2, 8, NAN, 33))
    assert err.value.identifier == NOT_FINITE_ID


def test_code_the_port_does_not_read_fails_before_any_tile(tmp_path, framework_spy):
    paths = _zarr_tiles(str(tmp_path), [(3, 4, 5)])
    fn = tmp_path / "f.txt"
    fn.write_text("@(x) imgaussfilt(x, 2)\n")
    with pytest.raises(NotImplementedError):
        T.process_tiles(paths, zarr_file=True, process_fun_path=[[str(fn)]], parse_cluster=False)
    assert framework_spy == []
