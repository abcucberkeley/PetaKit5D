"""Stitching of one volume end to end (petakit5d.stitch.volume) against MATLAB's
XR_stitching_frame_zarr_dev_v1 (the oracle tree's).

Oracle cases ``frame.<scene>``: each side writes the same tiles (tiff, or zarr in MATLAB's
layout) and, where the scene has registration, the same registration result
``xcorr/<name>.mat`` (the shifts given, since the registration is S04's), then runs the frame
(MATLAB's function, the port's through its MATLAB-name layer, with the same arguments), a
secondary channel after its primary where the scene has one. Compared: the stitched volumes
and every file the runs wrote (the converted tiles, the pixel-size marker, the stitch and block
information, the per-block JSON files, the distance maps, the MIPs, the results), the zarr
files by voxels and metadata after the axis mapping. The saved inputParser ``ip`` of the
stitch information file is not compared (MATLAB saves the object, the port a struct of its
Results). On Windows the feather scenes skip: their distance maps go into imdist/ itself there
(MATLAB's ispc branch), while the references are MATLAB's on Linux;
test_feather_frame_under_the_windows_rule runs the frame under that rule on every system.

Scenes: DSR tiles in a 2 x 2 grid blended by every method, with and without registration;
skewed space (two axis orders, flipped tiles, registration); 2D tiles; an IO scan; processed
zarr tiles (and zarr tiles cut by a box, the stand-in of S02's processing); a single tile;
an output box; resampled DSR; large files with pooled distance maps; 'separate'; tiles of
another data order; single and uint8 tiles; the tile information file; secondary channels of
the same sizes (feather and not), with a short tile (its own block information and maps), of
another stitched size and with the primary's coordinates; MIP stitching; local workers.

The pytest functions port the volume-level intent of StitchingTest (exact positions give the
ground truth for every method, the canvas convention, the IO scan and skewed-space
placements, processed zarr tiles, registration results applied, secondary channels, large
files, a single tile, the resample factor's axis order, the MIP made before publishing and
merged from the tasks, another run's result kept, local workers) and test the frame's own
logic (task sizes and order, the stitch information and registration files both ways, the
pixel-size marker, resuming a temporary result, the errors).
"""

from __future__ import annotations

import os
import shutil
import warnings
from dataclasses import dataclass, field

import numpy as np
import pytest

from petakit5d.stitch import distance as sd
from petakit5d.stitch import volume as V
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, InPath, Tree, Value, mcell
from petakit5d.testing.compare import tree
from petakit5d.testing.synth import CounterRNG, cut_tiles, stitch_tile_name, write_tiff, write_zarr

PX = 0.108


# ---------------------------------------------------------------------------------------------
# scenes: tiles cut from a ground truth, the frame's arguments (MATLAB's form)
# ---------------------------------------------------------------------------------------------


@dataclass
class Scene:
    gt: tuple  # ground truth shape (z, y, x) or (y, x)
    starts: list  # tile offsets (z0, y0, x0)
    tile: tuple  # tile size (z, y, x)
    args: dict  # MATLAB name-value pairs of the frame (besides the defaults below)
    pixel: tuple = (PX, PX, PX)  # voxel sizes (z, y, x) for the default stage coordinates
    stage: list | None = None  # stage coordinates (x, y, z) per tile
    tidx: list | None = None  # tile indices (x, y, z) per tile
    fmt: str = "tif"
    tile_dir: str = ""
    cls: str = "uint16"
    shift_vox: list | None = None  # registration result: shifts (x, y, z voxels) per tile
    secondary: dict | None = None  # a secondary channel: its arguments (CamB tiles)
    short_tile: int | None = None  # (secondary) this tile is 4 slices shorter
    mips: bool = False  # tiles' z MIPs in MIPs/ (stitchMIP)
    tile_info: bool = False  # the tiles in a tile information file
    seed: int = 1
    extra: dict = field(default_factory=dict)


def _grid(**args):
    return Scene(
        (32, 100, 120),
        [(0, 0, 0), (0, 36, 0), (0, 0, 48), (0, 36, 48)],
        (32, 64, 72),
        args,
    )


_SKEW = {"DSR": False, "dz": 0.3, "SkewAngle": 32.45, "Reverse": True}


def _skewed(axis_order="-x,y,z", **args):
    n = 10
    z1, z2 = (0, n) if axis_order.startswith("-") else (n, 0)
    return Scene(
        (50, 100, 64),
        [(z1, 0, 0), (z1, 36, 0), (z2, 0, 0), (z2, 36, 0)],
        (40, 64, 64),
        {**_SKEW, "axisOrder": axis_order, **args},
        stage=[[0, 0, 0], [0, 36 * PX, 0], [n * 0.3, 0, 0], [n * 0.3, 36 * PX, 0]],
        tidx=[[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]],
    )


_SMALL = {"batchSize": [32, 32, 16], "blockSize": [16, 16, 16]}
_XC = {"xcorrShift": True}
_SHIFT = [[0, 0, 0], [0, 0, 0], [-2, 1, 0], [0, 0, 0]]

SCENES: dict[str, Scene] = {}
for _m in ("none", "feather", "max", "mean", "median"):
    SCENES[f"dsr.{_m}"] = _grid(BlendMethod=_m, **_SMALL)
SCENES["dsr.none.default_batches"] = _grid(BlendMethod="none")
SCENES["dsr.feather.xcorr"] = _grid(BlendMethod="feather", **_XC, **_SMALL)
SCENES["dsr.feather.xcorr"].shift_vox = _SHIFT
SCENES["dsr.mean.xcorr"] = _grid(BlendMethod="mean", **_XC)
SCENES["dsr.mean.xcorr"].shift_vox = _SHIFT
SCENES["dsr.feather.full"] = _grid(BlendMethod="feather", overlapType="full", halfOrder=[1, 2, 3])
SCENES["dsr.minus_x"] = _grid(BlendMethod="max", axisOrder="-x,y,z")
SCENES["dsr.out_bbox"] = _grid(BlendMethod="feather", outBbox=[3, 5, 2, 90, 100, 30])
SCENES["dsr.out_bbox_open"] = _grid(BlendMethod="none", outBbox=[3, 5, 2, 90, np.inf, np.inf])
SCENES["dsr.single"] = _grid(BlendMethod="feather", save16bit=False)
SCENES["dsr.single"].cls = "single"
SCENES["dsr.uint8"] = _grid(BlendMethod="mean")
SCENES["dsr.uint8"].cls = "uint8"
SCENES["dsr.separate"] = _grid(BlendMethod="separate")
SCENES["dsr.resampled"] = Scene(
    (32, 80, 108),
    [(0, 0, 0), (0, 0, 36), (0, 30, 0), (0, 30, 36)],
    (32, 50, 72),
    {"BlendMethod": "feather", "resampleType": "given", "resampleFactor": [2, 1, 1]},
    pixel=(PX, 2 * PX, PX),
)
SCENES["dsr.single_tile"] = Scene(
    (32, 100, 120), [(0, 0, 0)], (32, 64, 72), {"BlendMethod": "none"}
)
SCENES["dsr.single_tile_xcorr"] = Scene(
    (32, 100, 120), [(0, 0, 0)], (32, 64, 72), {"BlendMethod": "feather", **_XC}
)
SCENES["dsr.single_tile_xcorr"].shift_vox = [[0, 0, 0]]
SCENES["dsr.data_order_xyz"] = _grid(BlendMethod="mean", dataOrder="x,y,z")
SCENES["dsr.large_file"] = _grid(
    BlendMethod="feather", largeFile=True, poolSize=[2, 2, 2], **_SMALL
)
SCENES["dsr.large_file_xcorr"] = _grid(
    BlendMethod="feather", largeFile=True, poolSize=[2, 2, 2], **_XC, **_SMALL
)
SCENES["dsr.large_file_xcorr"].shift_vox = _SHIFT
SCENES["dsr.tile_info"] = _grid(BlendMethod="none")
SCENES["dsr.tile_info"].tile_info = True
SCENES["dsr.local_workers"] = _grid(BlendMethod="feather", localWorkerNum=2, **_SMALL)
# two tiles that overlap by 4 voxels in x, of which the registration's shift leaves one plane:
# no overlap any more (the second overlap pass)
SCENES["dsr.touch_after_shift"] = Scene(
    (32, 64, 124), [(0, 0, 0), (0, 0, 60)], (32, 64, 64), {"BlendMethod": "mean", **_XC},
    shift_vox=[[0, 0, 0], [3, 0, 0]],
    # (the second tile's shifted x exactly the first tile's last plane, 63 * 0.108)
    extra={"touch_x": 63 * PX},
)  # fmt: skip
SCENES["skewed.feather"] = _skewed(BlendMethod="feather")
SCENES["skewed.plus_x.none"] = _skewed("x,y,z", BlendMethod="none")
SCENES["skewed.feather.xcorr"] = _skewed(BlendMethod="feather", **_XC)
SCENES["skewed.feather.xcorr"].shift_vox = [[0, 0, 0], [0, -4, 0], [0, 0, 0], [0, -4, 0]]
SCENES["skewed.flipped"] = _skewed(BlendMethod="max", flippedTile=[0, 1, 0, 1])
SCENES["ioscan.none"] = Scene(
    (40, 64, 104),
    [(0, 0, 0), (0, 0, 40), (16, 0, 0), (16, 0, 40)],
    (24, 64, 64),
    {"DSR": False, "IOScan": True, "dz": 0.5, "BlendMethod": "none"},
    pixel=(0.5, PX, PX),
)
SCENES["2d.feather"] = Scene(
    (100, 120),
    [(0, 0), (36, 0), (0, 48), (36, 48)],
    (64, 72),
    {"BlendMethod": "feather"},
)
SCENES["2d.mean.xcorr"] = Scene(
    (100, 120),
    [(0, 0), (36, 0), (0, 48), (36, 48)],
    (64, 72),
    {"BlendMethod": "mean", **_XC, "blockSize": [64, 64, 1], "batchSize": [256, 256, 1]},
    shift_vox=[[0, 0, 0], [1, 0, 0], [0, -1, 0], [0, 0, 0]],
)
SCENES["zarr.dsr"] = _grid(BlendMethod="feather", zarrFile=True, ProcessedDirstr="DSR")
SCENES["zarr.dsr"].fmt, SCENES["zarr.dsr"].tile_dir = "zarr", "DSR"
SCENES["zarr.input_bbox"] = _grid(
    BlendMethod="none", zarrFile=True, InputBbox=[2, 3, 1, 60, 70, 30]
)
SCENES["zarr.input_bbox"].fmt = "zarr"
SCENES["zarr.tile_out_bbox"] = _grid(
    BlendMethod="max", zarrFile=True, tileOutBbox=[1, 2, 3, 62, 70, 31]
)
SCENES["zarr.tile_out_bbox"].fmt = "zarr"
# (skewed space: zarr tiles flipped by the processing, S02's stitch_process_zarr_tile)
SCENES["zarr.skewed_flipped"] = _skewed(BlendMethod="mean", zarrFile=True, flippedTile=[1, 0, 1, 0])
SCENES["zarr.skewed_flipped"].fmt = "zarr"
# a distance map per tile (DSR), partial files (joined in skewed space, ignored for DSR tiles)
SCENES["dsr.feather.per_tile_maps"] = _grid(BlendMethod="feather", singleDistMap=False)
SCENES["dsr.stitch_info_dir"] = _grid(BlendMethod="feather", stitchInfoDir="stitchInfo_t")
SCENES["dsr.part_file_ignored"] = _grid(BlendMethod="none")
SCENES["dsr.part_file_ignored"].extra = {"part_tile": 0}
SCENES["skewed.partial_file"] = _skewed(BlendMethod="max")
SCENES["skewed.partial_file"].extra = {"part_tile": 2}
SCENES["mip.feather"] = _grid(BlendMethod="feather", stitchMIP=[False, False, True])
SCENES["mip.feather"].mips = True
# secondary channels (CamB after CamA)
for _m in ("feather", "none"):
    s = _grid(BlendMethod=_m, **_XC)
    s.shift_vox = _SHIFT
    s.secondary = {"BlendMethod": _m, **_XC}
    SCENES[f"secondary.{_m}"] = s
s = _grid(BlendMethod="feather", **_XC)
s.shift_vox = _SHIFT
s.secondary = {"BlendMethod": "feather", **_XC}
s.short_tile = 1
SCENES["secondary.short_tile"] = s
s = _grid(BlendMethod="mean", **_XC)
s.shift_vox = _SHIFT
s.secondary = {"BlendMethod": "mean", **_XC, "outBbox": [1, 1, 1, 90, 100, 30]}
SCENES["secondary.out_bbox"] = s
s = _grid(BlendMethod="max")
s.secondary = {"BlendMethod": "max", "usePrimaryCoords": True}
SCENES["secondary.primary_coords"] = s
SCENES["secondary.primary_coords"].extra = {"b_stage_offset": 2 * PX}


def _gt(sc: Scene, cam: str = "A") -> np.ndarray:
    rng = CounterRNG(sc.seed, f"s08_{cam}")
    if sc.cls == "single":
        return (rng.uniform(sc.gt) * 1000 + 1).astype(np.float32)
    hi = 250 if sc.cls == "uint8" else 60000
    a = rng.integers(1, hi, sc.gt)
    return a.astype(np.uint8 if sc.cls == "uint8" else np.uint16)


def _tiles(sc: Scene, cam: str = "A") -> list[dict]:
    starts = [s[-2:] for s in sc.starts] if len(sc.gt) == 2 else sc.starts
    tiles = cut_tiles(_gt(sc, cam), starts, sc.tile[-len(sc.gt) :], sc.pixel, sc.stage, sc.tidx)
    if cam == "B" and sc.short_tile is not None:
        t = tiles[sc.short_tile]
        t["im"] = np.ascontiguousarray(t["im"][:-4])
    if cam == "B" and sc.extra.get("b_stage_offset"):
        for t in tiles:
            t["xyz"] = [t["xyz"][0] + sc.extra["b_stage_offset"], t["xyz"][1], t["xyz"][2]]
    return tiles


def _names(sc: Scene, cam: str = "A") -> list[str]:
    return [stitch_tile_name(t["idx"], cam) for t in _tiles(sc, cam)]


def make_root(name: str):
    """the tiles of a scene (both channels), as both sides read them"""

    def make(d):
        sc = SCENES[name]
        for cam in ("A", "B") if sc.secondary is not None else ("A",):
            out = os.path.join(d, sc.tile_dir) if sc.tile_dir else os.fspath(d)
            os.makedirs(out, exist_ok=True)
            for t, n in zip(_tiles(sc, cam), _names(sc, cam), strict=True):
                if sc.fmt == "zarr":
                    write_zarr(os.path.join(out, n + ".zarr"), t["im"], (16, 32, 32), "matlab")
                else:
                    write_tiff(os.path.join(out, n + ".tif"), t["im"])
                if sc.mips:
                    os.makedirs(os.path.join(d, "MIPs"), exist_ok=True)
                    mip = np.ascontiguousarray(t["im"].max(axis=0))
                    write_tiff(os.path.join(d, "MIPs", n + "_MIP_z.tif"), mip)
            if sc.extra.get("part_tile") is not None:
                # a partial file of a tile (its last 6 slices, as an acquisition continued)
                k = sc.extra["part_tile"]
                t = _tiles(sc, cam)[k]
                part = np.ascontiguousarray(t["im"][-6:][::-1])
                write_tiff(os.path.join(out, _names(sc, cam)[k] + "_part001.tif"), part)

    return make


def _layout(sc: Scene, tiles: list[dict]):
    """the frame's coordinates in stitching space (S03's), for the registration result"""
    from petakit5d.core import axes as ax
    from petakit5d.stitch import geometry as G

    a = sc.args
    sizes = np.array([t["im"].shape if t["im"].ndim == 3 else (1, *t["im"].shape) for t in tiles])
    rf = a.get("resampleFactor")
    return G.prepare_layout(
        np.array([t["xyz"] for t in tiles]),
        sizes,
        tile_idx=np.array([[*t["idx"], 0] for t in tiles], dtype=np.float64),
        axis_order=a.get("axisOrder", "x,y,z"),
        data_order=ax.axis_string_from_matlab(a.get("dataOrder", "y,x,z")),
        xy_pixel_size=PX,
        dz=a.get("dz", PX),
        skew_angle=a.get("SkewAngle", 32.45),
        io_scan=a.get("IOScan", False),
        dsr=a.get("DSR", True),
        resample_type=a.get("resampleType", "isotropic"),
        resample_factor=None if rf is None else ax.from_yxz(rf),
    )


def params(name: str) -> dict:
    """the case's struct for MATLAB: per run its tiles (relative to the root), coordinates and
    arguments, and the registration result"""
    sc = SCENES[name]
    runs = []
    xcorr = None
    for cam in ("A", "B") if sc.secondary is not None else ("A",):
        tiles = _tiles(sc, cam)
        names = _names(sc, cam)
        pre = (sc.tile_dir + "/") if sc.tile_dir else ""
        coords = np.array([t["xyz"] for t in tiles], dtype=np.float64)
        tidx = np.array([[*t["idx"], 0] for t in tiles], dtype=np.float64)
        args = {
            "DSR": True,
            "resampleType": "isotropic",
            "xyPixelSize": PX,
            "dz": PX,
            "axisOrder": "x,y,z",
            "xcorrShift": False,
            "parseCluster": False,
            **(sc.args if cam == "A" else {**_without(sc.args, ("outBbox",)), **sc.secondary}),
        }
        if "flippedTile" in args:
            args["flippedTile"] = np.asarray(args["flippedTile"], dtype=np.float64)
        path_args = []
        if cam == "B":
            args["isPrimaryCh"] = False
            path_args = [
                "stitchInfoFullpath",
                f"matlab_stitch/stitchInfo/{_names(sc)[0][:-21]}.mat",
            ]
        flat = []
        for k, v in args.items():
            flat += [k, v]
        if sc.tile_info:
            run = {"tiles": mcell(), "coordinates": np.zeros((0, 0)), "args": flat,
                   "pathArgs": ["tileInfoFullpath", "tile_info.mat"],
                   "tileInfo": {"tiles": [pre + n + ".tif" for n in names], "xyz": coords,
                                "tileIdx": tidx, "flippedTile": np.zeros((0, 0))}}  # fmt: skip
        else:
            flat += ["tileIdx", tidx]
            run = {
                # (processed tiles by their raw tiff names in the data folder)
                "tiles": [
                    n + ".tif" if sc.tile_dir else n + (".zarr" if sc.fmt == "zarr" else ".tif")
                    for n in names
                ],
                "coordinates": coords,
                "args": flat,
                "pathArgs": path_args,
            }
        run["result"] = names[0][:-21]
        runs.append(run)
        if cam == "A":
            if sc.shift_vox is not None:
                lay = _layout(sc, tiles)
                vs = lay.space.pixel_sizes
                # (the rows in the tiles' stitching order)
                d = np.asarray(sc.shift_vox, dtype=np.float64)[lay.sinds] * vs
                xyz_shift = lay.xyz + d
                if "touch_x" in sc.extra:
                    xyz_shift[1, 0] = sc.extra["touch_x"]
                xcorr = {"name": names[0][:-21], "xyz_shift": xyz_shift, "dxyz_shift": d}
    p = {"runs": runs}
    if xcorr is not None:
        p["xcorr"] = xcorr
    return p


def _without(d: dict, keys) -> dict:
    return {k: v for k, v in d.items() if k not in keys}


# (the tiles' paths are absolute, as the wrapper gives them: MATLAB's createzarr fails on a
# relative path, and its job framework counts relative inputs as missing)
FRAME_CODE = r"""
p = in{1}; d = in{2};
r = [d 'matlab_stitch/'];
if ~exist(r, 'dir'), mkdir(r); end
if isfield(p, 'xcorr')
	mkdir([r 'xcorr']);
	xyz_shift = p.xcorr.xyz_shift; dxyz_shift = p.xcorr.dxyz_shift;
	save('-v7.3', [r 'xcorr/' p.xcorr.name '.mat'], 'xyz_shift', 'dxyz_shift');
end
out = cell(1, numel(p.runs));
for k = 1 : numel(p.runs)
	q = p.runs{k};
	args = q.args;
	for j = 1 : 2 : numel(q.pathArgs)
		args = [args, {q.pathArgs{j}, [d q.pathArgs{j + 1}]}];
	end
	fns = cellfun(@(x) [d x], q.tiles, 'unif', 0);
	if isfield(q, 'tileInfo')
		tile_fullpaths = cellfun(@(x) [d x], q.tileInfo.tiles, 'unif', 0);
		xyz = q.tileInfo.xyz; tileIdx = q.tileInfo.tileIdx; flippedTile = q.tileInfo.flippedTile;
		save('-v7.3', [d 'tile_info.mat'], 'tile_fullpaths', 'xyz', 'tileIdx', 'flippedTile');
		fns = {};
	end
	XR_stitching_frame_zarr_dev_v1(fns, q.coordinates, args{:});
	out{k} = readzarr([r q.result '.zarr']);
end
"""


def run_frame(name: str, d: str) -> list:
    """the Python side: the same files and runs through the MATLAB-name layer"""
    from petakit5d import matlab
    from petakit5d.io import read_zarr, write_mat

    p = params(name)
    r = d + "matlab_stitch/"
    os.makedirs(r, exist_ok=True)
    if "xcorr" in p:
        x = p["xcorr"]
        os.makedirs(r + "xcorr", exist_ok=True)
        V.write_xcorr_result(f"{r}xcorr/{x['name']}.mat", x["xyz_shift"], x["dxyz_shift"])
    outs = []
    for q in p["runs"]:
        args = list(q["args"])
        pa = q["pathArgs"]
        for j in range(0, len(pa), 2):
            args += [pa[j], d + pa[j + 1]]
        fns = [] if "tileInfo" in q else [d + t for t in q["tiles"]]
        if "tileInfo" in q:
            ti = q["tileInfo"]
            write_mat(d + "tile_info.mat", {"tile_fullpaths": [d + t for t in ti["tiles"]],
                                            "xyz": ti["xyz"], "tileIdx": ti["tileIdx"],
                                            "flippedTile": None})  # fmt: skip
            fns = []
        matlab.XR_stitching_frame_zarr_dev_v1(fns, q["coordinates"], *args)
        outs.append(read_zarr(f"{r}{q['result']}.zarr"))
    return outs


def _frame_case(name: str) -> Case:
    n = 2 if SCENES[name].secondary is not None else 1
    return Case(
        f"frame.{name}",
        code=FRAME_CODE,
        args=[In("p"), InPath("root")],
        nargout=n,
        python=run_frame,
        py_args=[name, InPath("root")],
        out_roles=("vol",) * n,
        inputs={"p": Value(lambda name=name: params(name)), "root": Tree(make_root(name), "root")},
        tree=tree(ignore_vars=("t", "t1", "ip")),
        marks=("slow",),
    )


CASES = [_frame_case(n) for n in SCENES]

# the golden scenario stitch_feather_xcorr (GoldenOutputTest: two tiles, the second 400 brighter
# and listed 2 voxels off, registered and feather blended by XR_matlab_stitching_wrapper): its
# result, MATLAB's registration result and the tiles, for test_golden_stitch_feather_xcorr
GOLDEN_NAME = "Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs"
GOLDEN = Case(
    "golden.stitch_feather_xcorr.frame",
    code=(
        "root = fileparts(which('setup')); "
        "addpath(fullfile(root, 'tests'), fullfile(root, 'tests', 'golden')); "
        "w = [work, '/']; s = GoldenOutputTest.computeScenario('stitch_feather_xcorr', w); "
        f"x = load([w 'matlab_stitch/xcorr/{GOLDEN_NAME}.mat']); "
        f"t1 = readtiff([w '{GOLDEN_NAME}_000x_000y_000z_0000t.tif']); "
        f"t2 = readtiff([w '{GOLDEN_NAME}_001x_000y_000z_0000t.tif']); "
        "out = {s.stitched, x.xyz_shift, x.dxyz_shift, t1, t2};"
    ),
    nargout=5,
    outputs=("stitched", "xyz_shift", "dxyz_shift", "tile1", "tile2"),
    files=False,
    note="GoldenOutputTest/outputMatchesReference(scenario=stitch_feather_xcorr), the frame",
)
CASES.append(GOLDEN)


# MATLAB on Windows (ispc) writes the feather distance maps into imdist/ itself, not
# imdist/<first tile>/ (compute_tile_distance_transform's ispc branch), and so does the port
# there: the references of the feather scenes, made on Linux, hold the other folder and paths.
# The Windows layout is tested on every system by forcing that branch
# (test_feather_frame_under_the_windows_rule).
_WINDOWS_MAPS = (
    "the feather distance maps go into imdist/ on Windows (MATLAB's ispc branch); these "
    "references are MATLAB's on Linux (imdist/<first tile>/): "
    "test_feather_frame_under_the_windows_rule covers it"
)


def _writes_maps(scene: str) -> bool:
    """whether a scene's frame runs (primary or secondary) make feather distance maps"""
    sc = SCENES[scene]
    return "feather" in (sc.args.get("BlendMethod"), (sc.secondary or {}).get("BlendMethod"))


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    if case.python is None:
        oracle.reference(case)  # (a reference only: test_golden_stitch_feather_xcorr)
        return
    if sd._is_pc() and _writes_maps(case.id.removeprefix("frame.")):
        pytest.skip(_WINDOWS_MAPS)
    oracle.check(case, tmp_path)


def test_golden_stitch_feather_xcorr_given_matlabs_shifts(tmp_path):
    """the golden scenario stitch_feather_xcorr through the port's frame, with MATLAB's
    registration result (S04's registration is not there yet): the wrapper's call of the frame
    (its batch and chunk sizes), MATLAB's stitched volume bit for bit"""
    from petakit5d.core.arrays import from_matlab, to_matlab
    from petakit5d.io import read_zarr

    ref = oracle.reference(GOLDEN)
    data = [getattr(v, "data", None) for v in ref.outputs]
    if any(d is None for d in data):
        pytest.skip(
            "the golden scenario's arrays are in the golden store (the manifest has hashes)"
        )
    stitched, xyz_shift, dxyz_shift, t1, t2 = data
    root = str(tmp_path) + "/"
    fns = []
    for k, t in enumerate((t1, t2)):
        fns.append(f"{root}{GOLDEN_NAME}_{k:03d}x_000y_000z_0000t.tif")
        write_tiff(fns[-1], from_matlab(np.asarray(t)))
    os.makedirs(f"{root}matlab_stitch/xcorr")
    V.write_xcorr_result(f"{root}matlab_stitch/xcorr/{GOLDEN_NAME}.mat", xyz_shift, dxyz_shift)
    V.stitch_volume(
        fns, np.array([[0.0, 0.0, 0.0], [42 * PX, 0.0, 0.0]]),
        tile_idx=np.array([[0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]]), dsr=True,
        resample_type="isotropic", xy_pixel_size=PX, dz=PX, axis_order="x,y,z",
        xcorr_shift=True, xcorr_downsample=(1, 1, 1), xy_max_offset=10, z_max_offset=4,
        blend_method="feather", batch_size=(512, 512, 512), block_size=(256, 256, 256),
        parse_cluster=False,
    )  # fmt: skip
    out = read_zarr(f"{root}matlab_stitch/{GOLDEN_NAME}.zarr")
    got = to_matlab(out)
    assert got.shape == np.asarray(stitched).shape == (38, 114, 38)
    np.testing.assert_array_equal(got, np.asarray(stitched))


# ---------------------------------------------------------------------------------------------
# the frame's own logic
# ---------------------------------------------------------------------------------------------


def test_pixel_info_name():
    assert V.pixel_info_name(0.108, 0.108, 0.108) == "px0.108_py0.108_pz0.108"
    zp = 0.108 * np.sin(np.deg2rad(32.45)) * 0.3 / 0.108
    assert V.pixel_info_name(0.108, 0.108, zp).startswith("px0.108_py0.108_pz0.16")
    assert V.pixel_info_name(0.1234567, 1, 2.5e-7) == "px0.12346_py1_pz2.5e-07"


@pytest.mark.parametrize(
    ("nb", "bs", "cluster", "workers", "expected"),
    [
        (56, (512, 512, 512), False, 1, 56),
        (5000, (512, 512, 512), False, 1, 500),
        (56, (512, 512, 512), True, 1, 20),
        (40, (32, 32, 16), False, 1, 40),
        (100000, (512, 512, 512), True, 1, 20),
        (300000, (256, 256, 256), True, 1, 160),
        (2000000, (1024, 1024, 1024), True, 1, 6),
        (56, (512, 512, 512), False, 2, 7),
        (9, (500, 500, 500), False, 3, 1),
    ],
)
def test_blend_task_size_follows_matlab(nb, bs, cluster, workers, expected):
    # taskRatio = max(0.1, 512^3 / prod(batchSize)); the cluster's 20 per 512^3 at least, and
    # at most 5,000 jobs within round(50 * ratio); locally all, at most round(500 * ratio)
    assert V.blend_task_size(nb, bs, cluster, workers) == expected


def test_task_order_puts_the_largest_first_and_mixes_quarters():
    b = [3, 9, 9, 1, 5]
    assert V._task_order(b) == [1, 2, 4, 0, 3]
    n = 10004
    b = list(range(n))
    order = V._task_order(b)
    desc = list(range(n - 1, -1, -1))
    q = n // 4
    assert order[:4] == [desc[0], desc[q], desc[2 * q], desc[3 * q]]
    assert sorted(order) == list(range(n))


def test_stitch_info_round_trip(tmp_path):
    om = np.array([[False, True, True], [False, False, True], [False, False, False]])
    kw = {
        "params": {"BlendMethod": "feather"},
        "flipped_tile": [True, False, True],
        "overlap_regions": np.arange(18, dtype=float).reshape(3, 6),
        "overlap_matrix": om,
        "dxyz_shift": np.arange(9, dtype=float).reshape(3, 3) / 10,
        "shift_method": "grid",
        "p_im_sz": (34, 102, 122),
        "p_tile_sizes": np.array([[32, 64, 72]] * 3, dtype=float),
        "xyz_voxel_sizes": (0.108, 0.2, 0.3),
        "data_order": "z,y,x",
        "per_block_info_fullpaths": ["a/x_task_size_7/b_1_7.json", "a/x_task_size_7/b_8_9.json"],
        "imdist_fullpaths": ["d/m1.zarr", "d/m2.zarr"],
        "imdist_file_idx": [0, 1, 1],
        "xyz_orig": np.arange(9, dtype=float).reshape(3, 3),
    }
    fn = str(tmp_path / "si.mat")
    V.write_stitch_info(fn, **kw)
    from petakit5d.io import read_mat

    m = read_mat(fn)
    assert set(m) == set(V.STITCH_INFO_VARS)
    np.testing.assert_array_equal(m["pImSz"], [[102, 122, 34]])
    np.testing.assert_array_equal(m["pTileSizes"], [[64, 72, 32]] * 3)
    np.testing.assert_array_equal(m["data_order_mat"], [[2, 1, 3]])
    np.testing.assert_array_equal(m["imdistFileIdx"], [[1], [2], [2]])
    assert m["PerBlockInfoFullpaths"].shape == (2, 1)
    assert m["flippedTile"].dtype == bool
    r = V.read_stitch_info(fn)
    assert r["pImSz"] == (34, 102, 122)
    np.testing.assert_array_equal(r["pTileSizes"], kw["p_tile_sizes"])
    assert r["imdistFileIdx"] == [0, 1, 1]
    assert r["imdistFullpaths"] == kw["imdist_fullpaths"]
    assert r["PerBlockInfoFullpaths"] == kw["per_block_info_fullpaths"]
    np.testing.assert_array_equal(r["overlap_matrix"], om)
    assert r["data_order_mat"] == (2, 1, 0)
    assert (r["xp"], r["yp"], r["zp"]) == (0.108, 0.2, 0.3)
    # one shared map is a scalar; no maps are MATLAB's {} and []
    kw.update(imdist_fullpaths=[], imdist_file_idx=[])
    V.write_stitch_info(fn, **kw)
    m = read_mat(fn)
    assert m["imdistFullpaths"] == [] and m["imdistFileIdx"].shape == (0, 0)
    assert V.read_stitch_info(fn)["imdistFileIdx"] == []
    kw.update(imdist_fullpaths=["d/m.zarr"], imdist_file_idx=[0], flipped_tile=None)
    V.write_stitch_info(fn, **kw)
    m = read_mat(fn)
    assert float(m["imdistFileIdx"]) == 1.0 and m["flippedTile"].shape == (0, 0)
    assert V.read_stitch_info(fn)["flippedTile"] is None


def test_xcorr_result_round_trip(tmp_path):
    fn = str(tmp_path / "x.mat")
    xyz = np.arange(12, dtype=float).reshape(4, 3) * 0.108
    V.write_xcorr_result(fn, xyz + 1, xyz)
    a, b = V.read_xcorr_result(fn)
    np.testing.assert_array_equal(a, xyz + 1)
    np.testing.assert_array_equal(b, xyz)
    V.write_xcorr_result(fn, [1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
    a, _ = V.read_xcorr_result(fn)
    assert a.shape == (1, 3)


def test_the_package_exports_the_frame_and_its_files():
    """petakit5d.stitch exports the frame, its errors and the readers and writers of its files
    lazily, at the frame's planned name, and its table names each registered function it
    exports with its MATLAB name"""
    import importlib
    import re

    import petakit5d.stitch as st
    from petakit5d.params import REGISTRY
    from petakit5d.params.planned import PLANNED

    names = ("stitch_volume", "read_stitch_info", "write_stitch_info", "read_xcorr_result",
             "write_xcorr_result", "StitchVolumeError", "StitchVolumeWarning")  # fmt: skip
    for name in names:
        assert name in st.__all__ and name in dir(st), name
        assert getattr(st, name) is getattr(V, name)
    planned = {m: p for m, p, _, _ in PLANNED}["XR_stitching_frame_zarr_dev_v1"]
    pkg, pname = planned.rsplit(".", 1)
    assert getattr(importlib.import_module(pkg), pname) is V.stitch_volume
    for mname, (where, kind) in REGISTRY.index().items():
        if kind != "register" or not where.startswith("petakit5d.stitch."):
            continue
        pname = REGISTRY.find(mname).python_name.rsplit(".", 1)[-1]
        if pname in st.__all__:
            assert re.search(rf"^``{pname}``\s+{mname}\b", st.__doc__, re.M), (pname, mname)


def test_stitch_mip_flags():
    assert V._stitch_mip_flags(None) is None
    assert V._stitch_mip_flags([False, False, False]) is None
    assert V._stitch_mip_flags(True) == (True, False, False)
    # several: the first in MATLAB's order (y, x, z)
    assert V._stitch_mip_flags((True, True, False)) == (False, True, False)
    assert V._stitch_mip_flags((True, False, True)) == (False, False, True)


def test_process_function_rows():
    assert V._process_fun_rows("") is None
    assert V._process_fun_rows(("",)) is None
    assert V._process_fun_rows(("a.mat", "b.mat")) == [["a.mat", "b.mat"]]
    assert V._process_fun_rows((("a.mat",), ("b.mat",))) == [["a.mat"], ["b.mat"]]
    with pytest.raises(V.StitchVolumeError) as e:
        V._process_fun_rows("a.mat")
    assert e.value.identifier == "MATLAB:cellfun:NotACell"


def test_edge_function_text():
    assert V._edge_function(True, None, 0, 2) == "@(x)erodeVolumeBy2DProjection(x,2)"
    assert V._edge_function(True, (True, False, False), 0, 2) == ""
    assert V._edge_function(False, None, 1, 2) == "@(x)erodeVolumeBy2DProjection(x+1,2)"
    assert V._edge_function(False, None, 0, 1.5) == (
        "@(x)erodeVolumeBy2DProjection(x,1.500000e+00)"
    )


# ---------------------------------------------------------------------------------------------
# runs on small tiles (StitchingTest in intent)
# ---------------------------------------------------------------------------------------------

NAME = "Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs"


def _write(sc: Scene, root: str, cam: str = "A") -> tuple[list[str], np.ndarray, np.ndarray]:
    tiles = _tiles(sc, cam)
    names = _names(sc, cam)
    out = os.path.join(root, sc.tile_dir) if sc.tile_dir else root
    os.makedirs(out, exist_ok=True)
    fns = []
    for t, n in zip(tiles, names, strict=True):
        if sc.fmt == "zarr":
            write_zarr(os.path.join(out, n + ".zarr"), t["im"], (16, 32, 32))
            fns.append(f"{root}/{n}.tif" if sc.tile_dir else f"{out}/{n}.zarr")
        else:
            write_tiff(os.path.join(out, n + ".tif"), t["im"])
            fns.append(f"{root}/{n}.tif")
    coords = np.array([t["xyz"] for t in tiles], dtype=np.float64)
    tidx = np.array([[*t["idx"], 0] for t in tiles], dtype=np.float64)
    return fns, coords, tidx


def _stitch(root: str, fns, coords, tidx, **kw):
    os.makedirs(f"{root}/matlab_stitch", exist_ok=True)
    opts = {"dsr": True, "resample_type": "isotropic", "xy_pixel_size": PX, "dz": PX,
                "axis_order": "x,y,z", "xcorr_shift": False, "parse_cluster": False}  # fmt: skip
    opts.update(kw)
    V.stitch_volume(fns, coords, tile_idx=tidx, **opts)


def _result(root: str, name: str = NAME) -> np.ndarray:
    from petakit5d.io import read_zarr

    return read_zarr(f"{root}/matlab_stitch/{name}.zarr")


@pytest.mark.parametrize("method", ["none", "feather", "max", "mean", "median"])
def test_exact_positions_reproduce_ground_truth(tmp_path, method):
    # the canvas: the first tile corner at voxel 0 and 2 more (empty) voxels along each axis
    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    _stitch(root, fns, coords, tidx, blend_method=method, batch_size=(16, 32, 32),
            block_size=(16, 16, 16))  # fmt: skip
    out = _result(root)
    assert out.shape == (34, 102, 122)
    np.testing.assert_array_equal(out[:32, :100, :120], _gt(sc))
    assert not out[32:].any() and not out[:, 100:].any() and not out[:, :, 120:].any()


def test_io_scan_uses_stage_positions_as_is(tmp_path):
    sc = SCENES["ioscan.none"]
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    _stitch(root, fns, coords, tidx, dsr=False, io_scan=True, dz=0.5, blend_method="none")
    out = _result(root)
    assert out.shape == (42, 66, 106)
    np.testing.assert_array_equal(out[:40, :64, :104], _gt(sc))


@pytest.mark.parametrize("axis_order", ["x,y,z", "-x,y,z"])
def test_skewed_space_placement(tmp_path, axis_order):
    sc = _skewed(axis_order)
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    _stitch(root, fns, coords, tidx, dsr=False, dz=0.3, skew_angle=32.45, axis_order=axis_order,
            blend_method="feather")  # fmt: skip
    out = _result(root)
    assert out.shape == (52, 102, 66)
    np.testing.assert_array_equal(out[:50, :100, :64], _gt(sc))


def test_processed_zarr_tiles(tmp_path):
    sc = _grid()
    sc.fmt, sc.tile_dir = "zarr", "DSR"
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    _stitch(root, fns, coords, tidx, zarr_file=True, processed_dirstr="DSR",
            blend_method="feather")  # fmt: skip
    out = _result(root)
    assert out.shape == (34, 102, 122)
    np.testing.assert_array_equal(out[:32, :100, :120], _gt(sc))
    # (no tile processing: no converted tiles)
    assert sorted(os.listdir(f"{root}/DSR")) == sorted(n + ".zarr" for n in _names(sc))


def test_registration_result_is_applied(tmp_path):
    # the third tile listed 2 voxels off in x and 1 in y, its shift in the registration result
    sc = SCENES["dsr.mean.xcorr"]
    root = str(tmp_path)
    tiles = _tiles(sc)
    for k, t in enumerate(tiles):
        t["xyz"] = [t["xyz"][0] - _SHIFT[k][0] * PX, t["xyz"][1] - _SHIFT[k][1] * PX, t["xyz"][2]]
    fns, coords, tidx = _write(sc, root)
    coords = np.array([t["xyz"] for t in tiles])
    lay = _layout(sc, tiles)
    d = np.asarray(_SHIFT, dtype=np.float64)[lay.sinds] * PX
    os.makedirs(f"{root}/matlab_stitch/xcorr")
    V.write_xcorr_result(f"{root}/matlab_stitch/xcorr/{NAME}.mat", lay.xyz + d, d)
    _stitch(root, fns, coords, tidx, xcorr_shift=True, blend_method="mean")
    out = _result(root)
    # padding of 2 * #tiles voxels per axis on both sides: 4 in y and x, 2 in z
    assert out.shape == (32 + 4 + 2, 100 + 8 + 2, 120 + 8 + 2)
    np.testing.assert_array_equal(out[2:34, 4:104, 4:124], _gt(sc))
    si = V.read_stitch_info(f"{root}/matlab_stitch/stitchInfo/{NAME}.mat")
    np.testing.assert_array_equal(si["dxyz_shift"], d)


def test_registration_needs_its_result_until_s04(tmp_path, monkeypatch):
    monkeypatch.setattr(V, "REGISTRATION_MODULE", "petakit5d.stitch._no_such_registration")
    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    with pytest.raises(NotImplementedError, match="S04"):
        _stitch(root, fns, coords, tidx, xcorr_shift=True)


def test_single_tile_canvas_is_the_tile_plus_two_voxels(tmp_path):
    sc = SCENES["dsr.single_tile"]
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    _stitch(root, fns, coords, tidx, xcorr_shift=True)
    out = _result(root)
    assert out.shape == (34, 66, 74)
    np.testing.assert_array_equal(out[:32, :64, :72], _tiles(sc)[0]["im"])


def test_resample_factor_is_in_the_tiles_axis_order(tmp_path):
    # DSR tiles resampled 2x in y ([2, 1, 1] in [y x z]): the stage y steps are 2 x 0.108 um
    sc = SCENES["dsr.resampled"]
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    _stitch(root, fns, coords, tidx, resample_type="given", resample_factor=(1, 2, 1),
            blend_method="feather")  # fmt: skip
    out = _result(root)
    assert out.shape == (34, 82, 110)
    np.testing.assert_array_equal(out[:32, :80, :108], _gt(sc))


def test_large_file_without_registration(tmp_path):
    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    _stitch(root, fns, coords, tidx, blend_method="feather", large_file=True,
            pool_size=(2, 2, 2), batch_size=(16, 32, 32), block_size=(16, 16, 16))  # fmt: skip
    np.testing.assert_array_equal(_result(root)[:32, :100, :120], _gt(sc))
    with pytest.raises(V.StitchVolumeError, match="poolSize"):
        _stitch(
            root, fns, coords, tidx, large_file=True, result_path=f"{root}/matlab_stitch/b.zarr"
        )


def _primary_and_secondary(root, sc, **kw_b):
    fa, ca, ta = _write(sc, root, "A")
    fb, cb, tb = _write(sc, root, "B")
    _stitch(root, fa, ca, ta, blend_method=sc.args.get("BlendMethod", "none"))
    si = f"{root}/matlab_stitch/stitchInfo/{NAME}.mat"
    _stitch(root, fb, cb, tb, is_primary_ch=False, stitch_info_fullpath=si,
            blend_method=sc.args.get("BlendMethod", "none"), **kw_b)  # fmt: skip
    return _result(root), _result(root, NAME.replace("CamA", "CamB"))


def test_secondary_channel_uses_the_primary_information(tmp_path):
    sc = _grid(BlendMethod="feather")
    root = str(tmp_path)
    a, b = _primary_and_secondary(root, sc)
    np.testing.assert_array_equal(a[:32, :100, :120], _gt(sc, "A"))
    np.testing.assert_array_equal(b[:32, :100, :120], _gt(sc, "B"))
    # the primary's block information and distance maps: none of its own
    info = os.listdir(f"{root}/matlab_stitch/stitchInfo")
    assert not any("CamB" in f and "task_size" in f and not f.endswith(".mat") for f in info)
    assert all("CamB" not in f for f in os.listdir(f"{root}/matlab_stitch/imdist"))


@pytest.mark.parametrize("single_map", [True, False])
def test_feather_frame_under_the_windows_rule(tmp_path, monkeypatch, single_map):
    # compute_tile_distance_transform's ispc branch, forced on any system: the maps go into
    # imdist/ itself ('<stitch>//imdist//<tile>_wr_10.zarr' in the stitch information); the
    # blending of the primary and of its secondary reads them there and gives the volumes and
    # MIPs of the other rule (imdist/<first tile>/, the oracle cases' MATLAB references). The
    # tiles differ in their overlaps (tile k plus 1000 k), so that the weights show.
    sc = _grid()
    sc.cls = "single"
    names_ab = (NAME, NAME.replace("CamA", "CamB"))
    got = {}
    for pc in (False, True):
        monkeypatch.setattr(sd, "_is_pc", lambda pc=pc: pc)
        root = str(tmp_path / ("windows" if pc else "posix"))
        os.makedirs(root)
        for cam in ("A", "B"):
            tiles = _tiles(sc, cam)
            fns = []
            for k, (t, n) in enumerate(zip(tiles, _names(sc, cam), strict=True)):
                fns.append(f"{root}/{n}.tif")
                write_tiff(fns[-1], np.ascontiguousarray(t["im"] + np.float32(1000 * k)))
            coords = np.array([t["xyz"] for t in tiles], dtype=np.float64)
            tidx = np.array([[*t["idx"], 0] for t in tiles], dtype=np.float64)
            kw = {}
            if cam == "B":
                si = f"{root}/matlab_stitch/stitchInfo/{NAME}.mat"
                kw = {"is_primary_ch": False, "stitch_info_fullpath": si}
            _stitch(root, fns, coords, tidx, blend_method="feather", save_16bit=False,
                    single_dist_map=single_map, **kw)  # fmt: skip
        si = V.read_stitch_info(
            f"{root}/matlab_stitch/stitchInfo/{NAME}.mat", ["imdistFullpaths", "imdistFileIdx"]
        )
        maps = si["imdistFullpaths"]
        folder = f"{root}/matlab_stitch//imdist/" + ("" if pc else f"{_names(sc)[0]}/")
        assert len(maps) == (1 if single_map else 4)
        assert all(p.startswith(folder + "/") and os.path.isdir(p) for p in maps)
        assert sorted(os.listdir(folder)) == sorted(
            ["parameters.mat"] + [p[len(folder) + 1 :] for p in maps]
        )
        vols = [_result(root, n) for n in names_ab]
        mips = [open(f"{root}/matlab_stitch/MIPs/{n}_MIP_z.tif", "rb").read() for n in names_ab]
        got[pc] = (si["imdistFileIdx"], vols, mips)
    assert got[True][0] == got[False][0]
    for a, b in zip(got[True][1], got[False][1], strict=True):
        assert a.dtype == np.float32
        np.testing.assert_array_equal(a.view(np.uint32), b.view(np.uint32))
    assert got[True][2] == got[False][2]


def test_secondary_with_a_short_tile_gets_its_own_block_info(tmp_path):
    sc = _grid(BlendMethod="feather")
    sc.short_tile = 1
    root = str(tmp_path)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        a, b = _primary_and_secondary(root, sc)
    ids = [getattr(x.message, "identifier", "") for x in w]
    assert "XR_stitching_frame_zarr_dev_v1:tileSizes" in ids
    # (the third tile in the stitching order, sizes as MATLAB's [ny nx nz])
    text = str(w[ids.index("XR_stitching_frame_zarr_dev_v1:tileSizes")].message)
    name_b = _names(sc, "B")[1]
    assert text == (
        f"Tile 3 ({name_b}) is 64x72x28, the primary's 64x72x32: so this volume gets its own "
        "block information, with the primary's tile positions."
    )
    gt_b = _gt(sc, "B")
    # the short tile's last 4 slices are filled by nothing (tile 2 starts at y 36)
    np.testing.assert_array_equal(b[:28, :100, :120], gt_b[:28])
    info = os.listdir(f"{root}/matlab_stitch/stitchInfo")
    assert any(f.startswith(NAME.replace("CamA", "CamB")) and "block_info" in f for f in info)


def test_secondary_of_another_stitched_size_gets_its_own_block_info(tmp_path):
    sc = _grid(BlendMethod="mean")
    root = str(tmp_path)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        _, b = _primary_and_secondary(root, sc, out_bbox=(0, 0, 0, 30, 90, 100))
    ids = [getattr(x.message, "identifier", "") for x in w]
    assert "XR_stitching_frame_zarr_dev_v1:stitchedSize" in ids
    # (MATLAB's sizes [ny nx nz]: the box's 30 x 90 x 100 against the canvas' 34 x 102 x 122)
    text = str(w[ids.index("XR_stitching_frame_zarr_dev_v1:stitchedSize")].message)
    assert text.startswith("The stitched size 90x100x30 differs from the primary's 102x122x34,")
    assert b.shape == (30, 90, 100)
    np.testing.assert_array_equal(b, _gt(sc, "B")[:30, :90, :100])


def test_secondary_without_the_primary_information_fails(tmp_path):
    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    with pytest.raises(V.StitchVolumeError, match="does not exist"):
        _stitch(root, fns, coords, tidx, is_primary_ch=False, stitch_info_fullpath="none.mat")


def test_an_existing_result_is_kept(tmp_path, capsys):
    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    os.makedirs(f"{root}/matlab_stitch/{NAME}.zarr")
    _stitch(root, fns, coords, tidx)
    assert "already exist, skip it" in capsys.readouterr().out
    assert os.listdir(f"{root}/matlab_stitch/{NAME}.zarr") == []


def test_the_result_folder_must_exist(tmp_path):
    # MATLAB: fclose(fopen(...)) of the pixel-size marker fails without the result folder
    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    with pytest.raises(V.StitchVolumeError) as e:
        V.stitch_volume(fns, coords, tile_idx=tidx, dsr=True, parse_cluster=False,
                        xcorr_shift=False)  # fmt: skip
    assert e.value.identifier == "MATLAB:FileIO:InvalidFid"


def test_other_pixel_size_markers_are_replaced(tmp_path):
    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    os.makedirs(f"{root}/matlab_stitch")
    open(f"{root}/matlab_stitch/px1_py1_pz1", "w").close()
    _stitch(root, fns, coords, tidx)
    marks = [f for f in os.listdir(f"{root}/matlab_stitch") if f.startswith("px")]
    assert marks == ["px0.108_py0.108_pz0.108"]


def test_a_temporary_result_of_the_same_size_is_resumed(tmp_path, monkeypatch):
    from petakit5d.io import create_zarr, read_zarr, write_zarr

    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    tmp = f"{root}/matlab_stitch/{NAME}.zarr_u1"
    create_zarr(tmp, data_size=(34, 102, 122), block_size=(34, 102, 122), dtype="uint16")
    write_zarr(np.full((34, 102, 122), 7, np.uint16), tmp, create=False)
    seen = []
    real = V._blend

    def spy(**k):
        seen.append(read_zarr(tmp)[33, 101, 121])
        real(**k)

    monkeypatch.setattr(V, "_blend", spy)
    _stitch(root, fns, coords, tidx, uuid="u1")
    assert seen == [7]  # kept for the blending tasks
    # another size: made again
    shutil.rmtree(f"{root}/matlab_stitch/{NAME}.zarr")
    create_zarr(tmp, data_size=(3, 4, 5), block_size=(3, 4, 5), dtype="uint16")
    write_zarr(np.full((3, 4, 5), 7, np.uint16), tmp, create=False)
    seen.clear()
    _stitch(root, fns, coords, tidx, uuid="u1")
    assert seen == [0]
    # the same size in other chunks: made again
    shutil.rmtree(f"{root}/matlab_stitch/{NAME}.zarr")
    create_zarr(tmp, data_size=(34, 102, 122), block_size=(34, 51, 122), dtype="uint16")
    write_zarr(np.full((34, 102, 122), 7, np.uint16), tmp, create=False)
    seen.clear()
    _stitch(root, fns, coords, tidx, uuid="u1")
    assert seen == [0]
    # MATLAB's layout (order F) of the same size and chunks: the port writes only its own
    from petakit5d.testing.synth import write_zarr_matlab

    shutil.rmtree(f"{root}/matlab_stitch/{NAME}.zarr")
    write_zarr_matlab(tmp, np.full((102, 122, 34), 7, np.uint16), (102, 122, 34))
    seen.clear()
    _stitch(root, fns, coords, tidx, uuid="u1")
    assert seen == [0]


def test_mip_is_made_before_publishing_and_merged_from_the_tasks(tmp_path, monkeypatch):
    from petakit5d.io import read_tiff
    from petakit5d.jobs import publish

    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    mip = f"{root}/matlab_stitch/MIPs/{NAME}_MIP_z.tif"
    order = []
    real = publish.publish_result

    def spy(tmp, final, flag_path=""):
        if os.fspath(final).endswith(NAME + ".zarr"):
            order.append(os.path.exists(mip))
        return real(tmp, final, flag_path)

    monkeypatch.setattr(publish, "publish_result", spy)
    # (merged from the tasks' flags: the volume is not read again for it)
    from petakit5d.tools import mip as T

    def no_read(*a, **k):
        raise AssertionError("the z MIP was computed from the volume")

    monkeypatch.setattr(T, "save_mip_zarr", no_read)
    _stitch(root, fns, coords, tidx, blend_method="mean", batch_size=(16, 32, 32),
            block_size=(16, 16, 16))  # fmt: skip
    assert order == [True]
    out = _result(root)
    np.testing.assert_array_equal(read_tiff(mip), out.max(axis=0))


def test_large_volume_mips_are_made_by_the_large_file_mips(tmp_path, monkeypatch):
    from petakit5d.io import read_tiff
    from petakit5d.stitch import block as SB

    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    monkeypatch.setattr(SB, "mip_memory", lambda cluster_mode, config_file="": 1e-6)
    _stitch(root, fns, coords, tidx, blend_method="max")
    out = _result(root)
    for ax, lab in ((0, "z"), (1, "y"), (2, "x")):
        m = read_tiff(f"{root}/matlab_stitch/MIPs/{NAME}_MIP_{lab}.tif")
        want = out.max(axis=ax)
        np.testing.assert_array_equal(m, want if ax == 0 else want.T)


def test_result_of_another_run_is_kept(tmp_path, monkeypatch):
    from petakit5d.io import create_zarr

    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    final = f"{root}/matlab_stitch/{NAME}.zarr"
    real = V._blend

    def other_run_first(**k):
        create_zarr(final, data_size=(2, 2, 2), block_size=(2, 2, 2), dtype="uint16")
        real(**k)

    monkeypatch.setattr(V, "_blend", other_run_first)
    _stitch(root, fns, coords, tidx)
    from petakit5d.io import get_zarr_info

    assert get_zarr_info(final).size == (2, 2, 2)
    assert [f for f in os.listdir(f"{root}/matlab_stitch") if f.startswith(NAME + ".zarr_")] == []


def test_result_another_run_publishes_while_the_tasks_run_is_kept(tmp_path, monkeypatch):
    # the other run's result appears after this run's blending tasks finished: it and its
    # MIP are kept, this run makes no MIP and removes its temporary result
    from petakit5d.io import create_zarr, get_zarr_info
    from petakit5d.jobs import framework

    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    final = f"{root}/matlab_stitch/{NAME}.zarr"
    real = framework.run_tasks

    def tasks_then_other_run(*a, **k):
        done = real(*a, **k)
        if k.get("final_out_fullpath") == final:
            create_zarr(final, data_size=(2, 2, 2), block_size=(2, 2, 2), dtype="uint16")
        return done

    monkeypatch.setattr(framework, "run_tasks", tasks_then_other_run)
    _stitch(root, fns, coords, tidx)
    assert get_zarr_info(final).size == (2, 2, 2)
    assert not os.path.exists(f"{root}/matlab_stitch/MIPs/{NAME}_MIP_z.tif")
    assert [f for f in os.listdir(f"{root}/matlab_stitch") if f.startswith(NAME + ".zarr_")] == []


def test_secondary_channel_uses_the_task_size_of_the_primary(tmp_path):
    # the primary's blending tasks of 6 batches (2 local workers), which the secondary's
    # per-block information files are made for; alone it would take all 48 in one task
    sc = _grid()
    root = str(tmp_path)
    fa, ca, ta = _write(sc, root, "A")
    fb, cb, tb = _write(sc, root, "B")
    small = {"batch_size": (16, 32, 32), "block_size": (16, 16, 16)}
    _stitch(root, fa, ca, ta, blend_method="mean", local_worker_num=2, **small)
    si = f"{root}/matlab_stitch/stitchInfo/{NAME}.mat"
    _stitch(root, fb, cb, tb, is_primary_ch=False, stitch_info_fullpath=si,
            blend_method="mean", **small)  # fmt: skip
    info = os.listdir(f"{root}/matlab_stitch/stitchInfo")
    name_b = NAME.replace("CamA", "CamB")
    assert f"{name_b}_block_info_task_size_6.mat" in info
    np.testing.assert_array_equal(_result(root, name_b)[:32, :100, :120], _gt(sc, "B"))


@pytest.mark.xfail(
    reason="MATLAB code in a tile function needs S02's tile-function parser", strict=True,
    raises=NotImplementedError,
)  # fmt: skip
def test_process_function_is_applied_to_tiles(tmp_path):
    # a flat field of 1000, background 100, offset 1: x -> x - 99 (StitchingTest)
    from petakit5d.io import write_mat

    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    os.makedirs(f"{root}/flatfield")
    write_tiff(f"{root}/flatfield/flatfield.tif", np.full((64, 72), 1000, np.uint16))
    write_tiff(f"{root}/flatfield/background.tif", np.full((64, 72), 100, np.uint16))
    fun = (
        f"@(x)erodeVolumeBy2DProjection(processFFCorrectionFrame(x,'{root}/flatfield/"
        f"flatfield.tif','{root}/flatfield/background.tif',1,0.4,false),0)"
    )
    write_mat(f"{root}/processFunction.mat", {"usrFun": fun})
    _stitch(root, fns, coords, tidx, process_fun_path=(f"{root}/processFunction.mat",))
    np.testing.assert_array_equal(_result(root)[:32, :100, :120], _gt(sc) - 99)


def test_local_workers_give_the_same_result(tmp_path):
    sc = _grid()
    outs = []
    for k, workers in enumerate((2, 1)):
        root = str(tmp_path / str(k))
        fns, coords, tidx = _write(sc, root)
        _stitch(root, fns, coords, tidx, blend_method="feather", batch_size=(16, 32, 32),
                block_size=(16, 16, 16), local_worker_num=workers)  # fmt: skip
        outs.append(_result(root))
        si = V.read_stitch_info(f"{root}/matlab_stitch/stitchInfo/{NAME}.mat")
        assert (len(si["PerBlockInfoFullpaths"]) > 1) == (workers > 1)
    np.testing.assert_array_equal(outs[0], outs[1])


def test_without_slurm_the_cluster_setting_runs_here(tmp_path, monkeypatch):
    # checkSlurmCluster: without Slurm the volume is stitched in this process, as locally
    from petakit5d.testing import slurm_fake

    (tmp_path / "no_slurm").mkdir()
    slurm_fake.hide_slurm(monkeypatch, tmp_path / "no_slurm")
    sc = _grid()
    root = str(tmp_path / "data")
    fns, coords, tidx = _write(sc, root)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _stitch(root, fns, coords, tidx, blend_method="max", parse_cluster=True)
    np.testing.assert_array_equal(_result(root)[:32, :100, :120], _gt(sc))


def test_separate_writes_each_tile_on_a_canvas(tmp_path):
    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    _stitch(root, fns, coords, tidx, blend_method="separate")
    from petakit5d.io import read_zarr

    for t, n in zip(_tiles(sc), _names(sc), strict=True):
        a = read_zarr(f"{root}/matlab_stitch/{n}.zarr")
        assert a.shape == (34, 102, 122)
        s = [int(v) for v in np.asarray(t["xyz"]) / PX]
        np.testing.assert_array_equal(a[:32, s[1] : s[1] + 64, s[0] : s[0] + 72], t["im"])
    assert _result(root).shape == (34, 102, 122)


def test_save_multires_is_not_ported(tmp_path):
    with pytest.raises(NotImplementedError, match="saveMultires"):
        V.stitch_volume([str(tmp_path / "a.tif")], np.zeros((1, 3)), save_multires=True)


def test_edge_artifacts_write_the_function_file(tmp_path):
    from petakit5d.io import read_mat

    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    try:
        _stitch(root, fns, coords, tidx, edge_artifacts=2)
    except NotImplementedError:
        pass  # (MATLAB code needs S02's tile-function parser)
    tmp = f"{root}/matlab_stitch/tmp"
    (fn,) = os.listdir(tmp)
    assert fn.startswith("processFunction_") and fn.endswith(".mat")
    assert read_mat(f"{tmp}/{fn}")["usrFun"] == "@(x)erodeVolumeBy2DProjection(x,2)"


def test_mip_stitching(tmp_path):
    sc = _grid()
    sc.mips = True
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    os.makedirs(f"{root}/MIPs")
    for t, n in zip(_tiles(sc), _names(sc), strict=True):
        write_tiff(f"{root}/MIPs/{n}_MIP_z.tif", np.ascontiguousarray(t["im"].max(axis=0)))
    _stitch(root, fns, coords, tidx, blend_method="none", stitch_mip=(True, False, False))
    out = _result(root)
    assert out.shape == (102, 122)
    np.testing.assert_array_equal(out[:100, :120], _gt(sc).max(axis=0))


def test_the_matlab_name_layer_takes_matlabs_values(tmp_path):
    from petakit5d import matlab

    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    os.makedirs(f"{root}/matlab_stitch")
    matlab.XR_stitching_frame_zarr_dev_v1(
        fns, coords, "DSR", True, "resampleType", "isotropic", "axisOrder", "x,y,z",
        "xcorrShift", False, "parseCluster", False, "tileIdx", tidx, "BlendMethod", "max",
        "batchSize", [32, 32, 16], "blockSize", [16, 16, 16], "outBbox", [1, 1, 1, 50, 60, 20],
    )  # fmt: skip
    out = _result(root)
    assert out.shape == (20, 50, 60)
    np.testing.assert_array_equal(out, _gt(sc)[:20, :50, :60])
    si = V.read_stitch_info(f"{root}/matlab_stitch/stitchInfo/{NAME}.mat", ["ip", "pImSz"])
    assert si["pImSz"] == (20, 50, 60)
    r = si["ip"]["Results"]
    np.testing.assert_array_equal(np.asarray(r["batchSize"]).ravel(), [32, 32, 16])
    assert r["dataOrder"] == "y,x,z"


def test_the_stitch_information_keeps_the_classes_of_the_call(tmp_path):
    """ip.Results in the stitch information holds the call's arrays in their classes, as
    MATLAB's inputParser does: a uint16 tile index stays uint16, the coordinates double"""
    sc = _grid()
    root = str(tmp_path)
    fns, coords, tidx = _write(sc, root)
    _stitch(root, fns, coords, tidx.astype(np.uint16), blend_method="max",
            batch_size=(16, 32, 32), block_size=(16, 16, 16))  # fmt: skip
    np.testing.assert_array_equal(_result(root)[:32, :100, :120], _gt(sc))
    si = V.read_stitch_info(f"{root}/matlab_stitch/stitchInfo/{NAME}.mat", ["ip"])
    r = si["ip"]["Results"]
    assert r["tileIdx"].dtype == np.uint16 and r["coordinates"].dtype == np.float64
    np.testing.assert_array_equal(r["tileIdx"], tidx)
    np.testing.assert_array_equal(r["coordinates"], coords)
