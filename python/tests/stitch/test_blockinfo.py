"""Stitching's block information (petakit5d.stitch.blockinfo) against MATLAB's
stitch_process_block_info and stitch_organize_block_info (the oracle tree's).

Oracle cases: 52 tile layouts (grids of 1 to 12 tiles with overlaps and jitter, 2D tiles, tiles
of other sizes, tiles partly or wholly outside the volume, layouts without overlaps, borders,
another data order, small task sizes and folders of few files); each runs
stitch_process_block_info through MATLAB's local job framework, and the files each side writes
are compared: the block information file, every per-block JSON file (MATLAB's text) and the
organizing tasks' flags. A secondary channel uses the primary's folder, and fails without it.

The pytest functions check the conversions of the per-block information (MATLAB's JSON shapes:
one tile, equal tile counts, none), the task sizes and file names (folders of files, the
power-of-two split), the file sizes in task order, the skip and overwrite of organizing tasks,
the publish order (the JSON files before the flag), that MATLAB's overlap test of the tiles
and the batches is an interval overlap, and folders that jobs share.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, field

import numpy as np
import pytest

from petakit5d.stitch import blockinfo as B
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value, Work, mcell
from petakit5d.testing.compare import MCell, mchar
from petakit5d.testing.synth import CounterRNG

# ---------------------------------------------------------------------------------------------
# layouts
# ---------------------------------------------------------------------------------------------


@dataclass
class Layout:
    """tiles at voxel offsets in a stitched volume, with what the frame computes for them:
    the overlaps (second pass, unit voxels), the half and full overlap regions and the pairs"""

    im_sizes: np.ndarray  # (n, 3) array order
    shift: np.ndarray  # (n, 3) x, y, z voxel offsets
    nv_size: tuple[int, int, int]
    batch_size: tuple[int, int, int]
    border: tuple[int, int, int] = (0, 0, 0)
    task_size: int = 10
    max_files: int = 20000
    stitch_2d: bool = False
    data_order: str = "z,y,x"
    overlap_type: str = "half"
    half_order: tuple[int, int, int] = (0, 1, 2)  # MATLAB's default [3, 2, 1]: z first
    om: np.ndarray = field(default=None, repr=False)  # type: ignore[assignment]
    half: list = field(default=None, repr=False)  # type: ignore[assignment]
    full: list = field(default=None, repr=False)  # type: ignore[assignment]
    pairs: np.ndarray = field(default=None, repr=False)  # type: ignore[assignment]

    def __post_init__(self) -> None:
        from petakit5d.stitch import geometry as G

        dop = G.parse_data_order(self.data_order)
        s = np.asarray(self.im_sizes, dtype=np.float64).reshape(-1, 3)
        sh = np.asarray(self.shift, dtype=np.float64).reshape(-1, 3)
        # cuboids in unit voxels (x, y, z), as the frame's with voxel sizes 1
        cub = np.hstack((sh, sh + s[:, list(dop)] - 1))
        om, _ = G.find_overlaps(cub, self.stitch_2d)
        space = G.StitchSpace(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, False)
        self.om, self.half, self.full, self.pairs = G.half_overlap_regions(
            cub,
            om,
            space,
            overlap_type=self.overlap_type,
            half_order=self.half_order,
            stitch_2d=self.stitch_2d,
        )

    @property
    def n(self) -> int:
        return int(np.asarray(self.im_sizes).reshape(-1, 3).shape[0])

    def tile_fns(self, prefix: str = "tiles/") -> list[str]:
        return [f"{prefix}tile_{k + 1:02d}.zarr" for k in range(self.n)]


def _grid_layout(seed: int) -> Layout:
    """a random grid of tiles: 1 to 3 per axis in y and x (1 or 2 in z), overlaps of 1 to 8
    voxels with jitter, sometimes 2D, of other sizes, outside the volume or with a border"""
    rng = CounterRNG(seed, "s07_layout")

    def ri(lo, hi):
        return int(rng.integers(lo, hi, 1)[0])

    two_d = seed % 9 == 4
    gy, gx = ri(1, 3), ri(1, 3)
    gz = 1 if two_d else ri(1, 2)
    t = [1 if two_d else ri(4, 12), ri(9, 24), ri(9, 24)]
    o = [0 if two_d else ri(1, 4), ri(1, 8), ri(1, 8)]
    sizes, shifts = [], []
    for kz in range(gz):
        for ky in range(gy):
            for kx in range(gx):
                ts = list(t)
                if seed % 6 == 1:
                    # tiles of other sizes (secondary channels, a stopped acquisition)
                    ts = [max(1, ts[0] - ri(0, 2)) if not two_d else 1, ts[1] - ri(0, 3), ts[2]]
                pos = [kz * (t[0] - o[0]), ky * (t[1] - o[1]), kx * (t[2] - o[2])]
                jit = [0 if two_d else ri(-1, 1), ri(-2, 2), ri(-2, 2)]
                p = [max(0, pos[q] + jit[q]) for q in range(3)]
                sizes.append(ts)
                shifts.append([p[2], p[1], p[0]])
    sizes = np.asarray(sizes, dtype=np.int64)
    shifts = np.asarray(shifts, dtype=np.int64)
    pad = [0 if two_d else ri(0, 2), ri(0, 3), ri(0, 3)]
    shifts = shifts + np.asarray([pad[2], pad[1], pad[0]])
    ext = [int(np.max(shifts[:, 2 - q] + sizes[:, q])) for q in range(3)]
    nv = [ext[q] + ri(0, 3) for q in range(3)]
    if two_d:
        nv[0] = 1
    if seed % 8 == 5:
        # the volume cuts the tiles at its end
        nv = [max(1, nv[0] - ri(0, 3)), nv[1] - ri(3, 8), nv[2] - ri(3, 8)]
    if seed % 10 == 7:
        # a tile partly before the volume's start (negative offsets)
        shifts[0] = shifts[0] - np.asarray([ri(1, 4), ri(1, 4), 0 if two_d else 1])
    bs = [1 if two_d else ri(3, 14), ri(5, 20), ri(5, 20)]
    border = (0, 0, 0)
    if seed % 11 == 3:
        border = (0 if two_d else 1, ri(1, 2), ri(1, 2))
    task = ri(1, 9)
    max_files = 20000 if seed % 7 != 2 else ri(2, 4)
    otype = ("half", "none", "full", "zero")[seed % 4] if seed % 3 == 0 else "half"
    return Layout(
        im_sizes=sizes,
        shift=shifts,
        nv_size=tuple(nv),
        batch_size=tuple(bs),
        border=border,
        task_size=task,
        max_files=max_files,
        stitch_2d=two_d,
        overlap_type=otype,
    )


def _special_layouts() -> dict[str, Layout]:
    one = np.asarray([[6, 14, 17]])
    return {
        # one tile, one batch, task size larger than the batches
        "single": Layout(one, [[0, 0, 0]], (6, 14, 17), (8, 16, 20), task_size=20),
        # one tile in a volume larger than it, small batches
        "single_padded": Layout(one, [[2, 1, 0]], (7, 19, 22), (3, 5, 7), task_size=4),
        # a tile wholly outside the volume
        "outside": Layout(
            np.asarray([[5, 12, 12], [5, 12, 12]]), [[0, 0, 0], [40, 0, 0]], (5, 12, 12), (5, 6, 6)
        ),
        # tiles that don't overlap
        "no_overlap": Layout(
            np.asarray([[4, 10, 10], [4, 10, 10]]), [[0, 0, 0], [12, 0, 0]], (4, 10, 22), (4, 7, 7)
        ),
        # tiles whose boxes touch in one plane (an overlap of one voxel counts as none)
        "touching": Layout(
            np.asarray([[4, 10, 10], [4, 10, 10]]), [[0, 0, 0], [9, 0, 0]], (4, 10, 19), (4, 5, 5)
        ),
        # tiles of the other data order ('x,y,z': x first)
        "xyz_order": Layout(
            np.asarray([[12, 10, 5], [12, 10, 5]]),
            [[0, 0, 0], [8, 2, 0]],
            (20, 12, 5),
            (7, 6, 5),
            data_order="x,y,z",
        ),
        # a border around each batch
        "border": Layout(
            np.asarray([[6, 16, 16], [6, 16, 16]]),
            [[0, 0, 0], [12, 1, 0]],
            (6, 17, 28),
            (3, 8, 8),
            border=(1, 2, 2),
            task_size=3,
        ),
        # more tasks than files per folder: subfolders, organizing tasks per folder
        "folders": Layout(
            np.asarray([[5, 20, 20], [5, 20, 20], [5, 20, 20]]),
            [[0, 0, 0], [15, 0, 0], [7, 14, 1]],
            (6, 34, 35),
            (2, 4, 5),
            task_size=2,
            max_files=3,
        ),
        # 2D tiles in a row
        "row_2d": Layout(
            np.asarray([[1, 15, 18]] * 3),
            [[0, 0, 0], [14, 1, 0], [28, 0, 0]],
            (1, 16, 46),
            (1, 7, 9),
            stitch_2d=True,
            task_size=5,
        ),
        # tiles stacked in z (the split axis is z)
        "z_stack": Layout(
            np.asarray([[8, 12, 12], [8, 12, 12]]),
            [[0, 0, 0], [0, 0, 5]],
            (13, 12, 12),
            (4, 6, 6),
            task_size=6,
        ),
    }


LAYOUTS: dict[str, Layout] = {**_special_layouts()}
for _seed in range(42):
    LAYOUTS[f"grid{_seed:02d}"] = _grid_layout(_seed)


def _regions_mcell(regions) -> MCell:
    from petakit5d.stitch.geometry import region_to_matlab

    return MCell(
        (len(regions), 1), [mcell(region_to_matlab(a), region_to_matlab(b)) for a, b in regions]
    )


def layout_matlab(L: Layout, tile_prefix: str = "tiles/") -> dict:
    """the frame's arguments of stitch_process_block_info for a layout, in MATLAB's values"""
    from petakit5d.core.axes import to_yxz
    from petakit5d.stitch import geometry as G

    n = L.n
    om = np.asarray(L.om, dtype=bool).reshape(n, n)
    return {
        "int_xyz_shift": np.asarray(L.shift, dtype=np.float64).reshape(n, 3),
        "imSizes": np.asarray(L.im_sizes, dtype=np.float64).reshape(n, 3)[:, [1, 2, 0]],
        "nvSize": np.asarray(to_yxz(L.nv_size), dtype=np.float64),
        "blockSize": np.asarray(to_yxz(L.batch_size), dtype=np.float64),
        "data_order_mat": np.asarray(
            G.data_order_axes_to_matlab(G.parse_data_order(L.data_order)), dtype=np.float64
        ),
        "overlap_matrix": om,
        "ol_region_cell": _regions_mcell(L.full),
        "half_ol_region_cell": _regions_mcell(L.half),
        "overlap_map_mat": G.overlap_map_to_matlab(L.pairs, n).reshape(-1, 2),
        "BorderSize": np.asarray(to_yxz(L.border), dtype=np.float64),
        "tileFns": MCell((n, 1), [mchar(f) for f in L.tile_fns(tile_prefix)]),
        "stitch2D": bool(L.stitch_2d),
        "taskSize": float(L.task_size),
        "maxFileNumPerFolder": float(L.max_files),
    }


def run_block_info(L: Layout, stitch_info: str, *, primary: bool = True, tile_prefix="tiles/"):
    """process_block_info of a layout (as the frame calls it), locally"""
    return B.process_block_info(
        L.shift,
        L.im_sizes,
        L.nv_size,
        L.batch_size,
        L.data_order,
        L.om,
        L.full,
        L.half,
        L.pairs,
        L.border,
        L.tile_fns(tile_prefix),
        stitch_info,
        "vol",
        primary,
        stitch_info_fullpath=f"{stitch_info}/vol.mat",
        stitch_2d=L.stitch_2d,
        task_size=L.task_size,
        max_file_num_per_folder=L.max_files,
        uuid="u0",
        parse_cluster=False,
    )


# ---------------------------------------------------------------------------------------------
# oracle cases
# ---------------------------------------------------------------------------------------------

_INFO_CALL = (
    "[f, p, b] = stitch_process_block_info(L.int_xyz_shift, L.imSizes, L.nvSize, L.blockSize, "
    "L.data_order_mat, L.overlap_matrix, L.ol_region_cell, L.half_ol_region_cell, "
    "L.overlap_map_mat, L.BorderSize, L.tileFns, si, 'vol', primary, 'stitchInfoFullpath', "
    "[si, '/vol.mat'], 'stitch2D', L.stitch2D, 'taskSize', L.taskSize, 'maxFileNumPerFolder', "
    "L.maxFileNumPerFolder, 'uuid', 'u0', 'parseCluster', false); out = {f, p, b};"
)
_INFO_CODE = "L = in{1}; si = [in{2}, 'stitchInfo']; mkdir(si); primary = true; " + _INFO_CALL
# a secondary channel: the primary's per-block folder exists (made by the case) or not
_SECONDARY_CODE = (
    "L = in{1}; si = [in{2}, 'stitchInfo']; mkdir(si); primary = false; "
    "if in{3}, mkdir(sprintf('%s/vol_task_size_%d/', si, L.taskSize)); end; " + _INFO_CALL
)


def p_info(name: str, work: str) -> list:
    si = work + "stitchInfo"
    os.makedirs(si)
    f, p, b = run_block_info(LAYOUTS[name], si)
    return [f, p, b]


def p_secondary(name: str, make_folder: bool, work: str) -> list:
    L = LAYOUTS[name]
    si = work + "stitchInfo"
    os.makedirs(si)
    if make_folder:
        os.makedirs(f"{si}/vol_task_size_{L.task_size}/")
    f, p, b = run_block_info(L, si, primary=False)
    return [f, p, b]


def _layout_value(name: str):
    return lambda: layout_matlab(LAYOUTS[name])


CASES: list[Case] = [
    Case(
        f"info.{name}",
        code=_INFO_CODE,
        args=[In("layout"), Work()],
        nargout=3,
        python=p_info,
        py_args=[name, Work()],
        inputs={"layout": Value(_layout_value(name))},
    )
    for name in LAYOUTS
]
for _name in ("grid03", "folders"):
    CASES.append(
        Case(
            f"info.secondary.{_name}",
            code=_SECONDARY_CODE,
            args=[In("layout"), Work(), True],
            nargout=3,
            python=p_secondary,
            py_args=[_name, True, Work()],
            inputs={"layout": Value(_layout_value(_name))},
        )
    )
CASES.append(
    Case(
        "error.info.secondary_without_primary",
        code=_SECONDARY_CODE,
        args=[In("layout"), Work(), False],
        nargout=3,
        python=p_secondary,
        py_args=["grid03", False, Work()],
        inputs={"layout": Value(_layout_value("grid03"))},
        expect_error=True,
        py_error=B.BlockInfoError,
    )
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# the per-block information in MATLAB's shapes
# ---------------------------------------------------------------------------------------------


def _tile(t, b=(0, 0, 0, 2, 3, 4), regions=()):
    w = tuple(v + 5 for v in b)
    return B.TileBlock(t, b, w, tuple(v + 1 for v in b), tuple(regions), (0, 0, 0))


def test_per_block_info_round_trips_through_matlab_json(tmp_path):
    r = B.RegionBlock((0, 1, 1, 2, 3, 4), (0, 6, 6, 2, 8, 9), (0, 2, 2, 2, 4, 5))
    empty = B.RegionBlock(None, None, None)
    blocks = [
        [],
        [_tile(0)],
        [_tile(0, regions=[r]), _tile(3, regions=[r, empty])],
        [_tile(1), _tile(2)],
    ]
    text = B.block_info_json(blocks)
    fn = tmp_path / "stitch_block_info_blocks_1_4.json"
    fn.write_text(text, encoding="utf-8")
    assert B.read_per_block_info(fn) == blocks
    # MATLAB's PrettyPrint: a batch without tiles is [], one tile an object
    d = json.loads(text)
    assert d[0] == [] and isinstance(d[1], dict) and isinstance(d[2], list)
    assert d[2][0]["mblockInfo"]["bCoords"] == [2, 2, 1, 3, 4, 2]
    assert d[2][1]["mblockInfo"][1]["bCoords"] == []
    assert text.startswith('[\n  [],\n  {\n    "tileInd": 1,\n    "bCoords": [\n')


def test_jsondecode_shapes_give_the_tiles_of_each_batch():
    # MATLAB's jsondecode: one tile per batch is a struct array of batches; equal tile counts
    # a batches x tiles struct array (arrayfun takes its rows)
    from petakit5d.io.matfile import StructArray

    m = B.tile_blocks_to_matlab([[_tile(0)], [_tile(1)]])
    assert B.tile_blocks_from_matlab(m) == [[_tile(0)], [_tile(1)]]
    rows = B.tile_blocks_to_matlab([[_tile(0)], [_tile(2)]])
    sa = StructArray((2, 1), list(rows[0]), rows)
    assert B.tile_blocks_from_matlab(sa) == [[_tile(0)], [_tile(2)]]
    grid = [
        [B._tile_to_matlab(_tile(i)) for i in (0, 1)],
        [B._tile_to_matlab(_tile(i)) for i in (2, 3)],
    ]
    flat = [grid[r][c] for c in range(2) for r in range(2)]
    sa2 = StructArray((2, 2), list(flat[0]), flat)
    assert B.tile_blocks_from_matlab(sa2) == [[_tile(0), _tile(1)], [_tile(2), _tile(3)]]


def test_per_block_info_from_a_mat_file(tmp_path):
    from petakit5d.io import write_mat

    blocks = [[_tile(0)], [], [_tile(1), _tile(2)]]
    cell = np.empty((3, 1), dtype=object)
    for k, v in enumerate(B.tile_blocks_to_matlab(blocks)):
        cell[k, 0] = np.zeros((0, 0)) if v is None else v
    fn = tmp_path / "info.mat"
    write_mat(fn, {"blockInds": np.array([1.0, 3.0]), "stitchBlockInfo_t": cell})
    assert B.read_per_block_info(fn) == blocks


def test_a_malformed_box_is_matlabs_assignment_error():
    with pytest.raises(B.BlockInfoError) as e:
        B.tile_blocks_from_matlab([{"tileInd": 1, "bCoords": [1, 2, 3]}])
    assert e.value.identifier == "MATLAB:subsassigndimmismatch"


# ---------------------------------------------------------------------------------------------
# tasks and files
# ---------------------------------------------------------------------------------------------


def _blk_task_size_matlab(num_blocks, task_size, n_file):
    # stitch_process_block_info.m lines 156-169, literally
    blk = max(1e5, min(1e6, task_size * 500))
    blk = math.ceil(blk / task_size) * task_size
    if math.ceil(num_blocks / task_size) > n_file:
        blk = min(blk, n_file * task_size)
        if blk < n_file * task_size:
            i = 0
            for i in range(1, math.ceil(math.log2(n_file * task_size / blk)) + 1):
                cur = n_file * task_size / 2**i
                if blk > cur or math.floor(cur) != cur:
                    break
            blk = n_file * task_size / 2 ** (i - 1)
    return blk


@pytest.mark.parametrize(
    ("num_blocks", "task_size", "n_file", "expected"),
    [
        (10, 10, 20000, 100000),
        (75, 7, 20000, 100002),  # rounded up to whole tasks
        (10**6, 3000, 20000, 1002000),  # at most 1e6, rounded up to whole tasks
        (100, 2, 3, 6),  # one folder of files per organizing task
        (5 * 10**6, 20, 20000, 200000),  # 400000 per folder; no break: halved once
        (5 * 10**6, 30, 20000, 150000),  # 600000 per folder, halved twice
        (5 * 10**6, 25, 20000, 125000),  # 500000 per folder, halved twice
        (5 * 10**6, 21, 20000, 105000),  # 420000: halved twice
        (5 * 10**7, 3, 70000, 105000),  # 210000: halved once
    ],
)
def test_organizing_task_size_is_matlabs(num_blocks, task_size, n_file, expected):
    got = B.block_task_size(num_blocks, task_size, n_file)
    assert got == expected == _blk_task_size_matlab(num_blocks, task_size, n_file)


def test_per_block_info_paths_split_into_folders():
    folders, paths = B.per_block_info_paths("si/vol_task_size_2/", 2, 13, 3)
    assert folders == [
        "si/vol_task_size_2//task_size_2_tasks_1_3/",
        "si/vol_task_size_2//task_size_2_tasks_4_6/",
        "si/vol_task_size_2//task_size_2_tasks_7_7/",
    ]
    assert (
        paths[0] == "si/vol_task_size_2//task_size_2_tasks_1_3//stitch_block_info_blocks_1_2.json"
    )
    assert (
        paths[-1]
        == "si/vol_task_size_2//task_size_2_tasks_7_7//stitch_block_info_blocks_13_13.json"
    )
    assert len(paths) == 7
    folders, paths = B.per_block_info_paths("p/", 4, 9, 20000)
    assert folders == ["p/"]
    assert paths == [
        "p//stitch_block_info_blocks_1_4.json",
        "p//stitch_block_info_blocks_5_8.json",
        "p//stitch_block_info_blocks_9_9.json",
    ]


def test_block_info_bytes_are_in_task_order(tmp_path):
    folder = f"{tmp_path}/a/"
    os.makedirs(folder)
    for s, n in ((1, 30), (11, 5), (21, 12), (101, 7)):
        with open(f"{folder}stitch_block_info_blocks_{s}_{s + 9}.json", "w") as f:
            f.write("x" * n)
    # (dir's order: 1, 101, 11, 21; by the first batch: 1, 11, 21, 101)
    assert B.block_info_bytes([folder]) == [30, 5, 12, 7]


def test_block_info_file_round_trips(tmp_path):
    from petakit5d.stitch.geometry import tile_canvas_boxes

    L = LAYOUTS["grid05"]
    cb, tb, valid = tile_canvas_boxes(L.shift, L.im_sizes, L.nv_size, L.data_order, L.stitch_2d)
    fn = str(tmp_path / "info.mat")
    B.write_block_info(
        fn,
        tile_fns=L.tile_fns(),
        int_xyz_shift=L.shift,
        im_sizes=L.im_sizes,
        overlap_matrix=L.om,
        half_regions=L.half,
        full_regions=L.full,
        pairs=L.pairs,
        canvas_boxes=cb,
        tile_boxes=tb,
        valid=valid,
        nv_size=L.nv_size,
        block_size=L.batch_size,
        border_size=L.border,
    )
    a = B.read_block_info(fn)
    assert a["tile_fns"] == L.tile_fns()
    assert np.array_equal(a["overlap_matrix"], L.om)
    assert a["half_regions"] == [tuple(map(tuple, p)) for p in L.half]
    assert a["full_regions"] == [tuple(map(tuple, p)) for p in L.full]
    assert np.array_equal(a["pairs"], L.pairs)
    assert np.array_equal(a["canvas_boxes"][valid], cb[valid])
    assert np.array_equal(a["valid"], valid)
    assert a["nv_size"] == tuple(L.nv_size) and a["block_size"] == tuple(L.batch_size)
    assert np.array_equal(a["im_sizes"], L.im_sizes)


def test_matlabs_tile_test_is_an_interval_overlap():
    # ol_4 repeats ol_3 (the tile's start in the batch is missing), but a tile inside the
    # batch has its end in it: the test is that of overlapping intervals
    rng = CounterRNG(3, "s07_overlap")
    a = rng.integers(1, 20, (4000, 4)).reshape(-1, 4)
    lo1, hi1 = np.minimum(a[:, 0], a[:, 1]), np.maximum(a[:, 0], a[:, 1])
    lo2, hi2 = np.minimum(a[:, 2], a[:, 3]), np.maximum(a[:, 2], a[:, 3])
    for k in range(a.shape[0]):
        got = B._matlab_overlap(
            np.array([[lo1[k]] * 3]),
            np.array([[hi1[k]] * 3]),
            np.array([lo2[k]] * 3),
            np.array([hi2[k]] * 3),
        )[0]
        assert got == (lo1[k] <= hi2[k] and lo2[k] <= hi1[k])


def test_organizing_skips_a_done_task_and_redoes_it_with_overwrite(tmp_path, capsys):
    L = LAYOUTS["grid01"]
    si = str(tmp_path / "si")
    os.makedirs(si)
    f, paths, _ = run_block_info(L, si)
    flag = f"{si}/vol_task_size_{L.task_size}//block_flags/x.mat"
    open(flag, "w").close()
    assert not B.organize_block_info([0, 0], L.task_size, f, f"{si}/redo/", flag)
    assert "already exist, skip them!" in capsys.readouterr().out
    os.makedirs(f"{si}/redo/")
    assert B.organize_block_info([0, 0], L.task_size, f, f"{si}/redo/", flag, overwrite=True,
                                 border_size=L.border)  # fmt: skip
    one = B.read_per_block_info(f"{si}/redo//stitch_block_info_blocks_1_1.json")
    assert one == B.read_per_block_info(paths[0])[:1]
    from petakit5d.io import read_mat

    assert list(np.asarray(read_mat(flag)["blockInds"]).ravel()) == [1.0, 1.0]


def test_organizing_writes_its_json_files_before_its_flag(tmp_path):
    from petakit5d.testing.publish_spy import PublishSpy

    L = LAYOUTS["folders"]
    si = str(tmp_path / "si")
    os.makedirs(si)
    with PublishSpy() as spy:
        f, paths, _ = run_block_info(L, si)
    moves = [c.dst for c in spy.calls]
    flags = [k for k, d in enumerate(moves) if "/block_flags/" in d]
    jsons = [k for k, d in enumerate(moves) if d.endswith(".json")]
    assert flags and len(jsons) == len(paths)
    # each organizing task's files are in place before its flag, which comes last
    assert max(jsons) < max(flags)
    assert all(os.path.exists(p) for p in paths)
    assert moves[0].endswith("vol_block_info_task_size_2.mat")


def test_shared_folders_that_exist_are_fine(tmp_path):
    # jobs of other channels or time points make the same folders (MATLAB's mkdir)
    L = LAYOUTS["folders"]
    si = str(tmp_path / "si")
    for sub in ("", "/vol_task_size_2/", "/vol_task_size_2//block_flags/",
                "/vol_task_size_2//task_size_2_tasks_1_3/"):  # fmt: skip
        os.makedirs(si + sub, exist_ok=True)
    f, paths, _ = run_block_info(L, si)
    assert all(os.path.exists(p) for p in paths)


def test_secondary_channel_reports_file_sizes_with_the_cluster(tmp_path, monkeypatch):
    L = LAYOUTS["grid03"]
    si = str(tmp_path / "si")
    os.makedirs(si)
    _, paths, _ = run_block_info(L, si)
    out = B.process_block_info(
        L.shift, L.im_sizes, L.nv_size, L.batch_size, L.data_order, L.om, L.full, L.half,
        L.pairs, L.border, L.tile_fns(), si, "vol2", False, stitch_info_fullpath=f"{si}/vol.mat",
        task_size=L.task_size, uuid="u1", parse_cluster=True,
    )  # fmt: skip
    assert out[1] == paths
    assert out[2] == [os.path.getsize(p) for p in paths]


def test_a_missing_border_is_matlabs_index_error():
    L = LAYOUTS["grid00"]
    from petakit5d.stitch.geometry import tile_canvas_boxes

    cb, tb, valid = tile_canvas_boxes(L.shift, L.im_sizes, L.nv_size, L.data_order, False)
    with pytest.raises(B.BlockInfoError) as e:
        B.organize_blocks(
            (0, 0), cb, tb, valid, L.om, L.half, L.pairs, L.nv_size, L.batch_size, None
        )
    assert e.value.identifier == "MATLAB:badsubscript"
