"""The blending tasks of stitching (petakit5d.stitch.block) against MATLAB's processStitchBlock,
stitch_merge_block_mips and processStitchSeparteTiles (the oracle tree's).

Oracle cases (each side makes the block information, the stitched zarr and runs every blending
task of the volume, through its MATLAB-name layer for processStitchBlock, then merges the MIPs;
the files are compared: the stitched volume, the flags with their MIPs, the JSON files):

- ``blend.<layout>.<method>.<class>``: 6 layouts (a 2 x 2 grid in batches of 1 to 4 tiles, three
  tiles in a row in one batch and in three, two tiles with a narrow overlap in a batch that
  isn't full, 2D tiles, tiles stacked in z, one tile that covers its batch and one that doesn't)
  x 'none', 'max', 'mean', 'median', 'feather' x uint16, uint8 (MATLAB's feather code) and single
  tiles (negative values) with zero regions; distance maps per tile and shared, read whole
  (small) or per batch;
- ``lf.<layout>.<class>``: the large-file feather blending with pooled distance maps: two tiles,
  one major tile that covers two minor ones, the same with a hole in the major tile (minor tiles
  read where it is empty), a minor tile the major one doesn't cover, two major tiles, and a full
  batch;
- ``blend.cast.*``: tiles of another class than the result; ``blend.nomip.*``: flags without
  MIPs; ``blend.task1.*``: one batch per task;
- ``merge.*``: stitch_merge_block_mips of flags as processStitchBlock writes them
  (StitchMergeBlockMipsTest);
- ``separate.*``: each tile on its own canvas (processStitchSeparteTiles, through XR_crop_zarr);
- ``blend.border.*``: batches with a border, cut to their own size; ``blend.Feather.*``: the
  method's name in another case (MATLAB's blocks read past the tiles, and 'Feather' failed at
  the merge, until bef98d1).

The pytest functions port StitchMergeBlockMipsTest and the blending tests of StitchingTest in
intent (the overlap split of 'none' and of 'half' in quarters, exact tiles giving the ground
truth for every method, the blend of every voxel by its definition, uint8 feather tiles as
uint16 ones, the merged MIP of a mean blend that ends in .5, separate tiles), and check the
publish order (the flag after the batches, the MIPs in it), the skips, the batches without
tiles, a region outside its tile, the pooling, the memory decisions, thread independence and
the blending from MATLAB's block information files (with the golden store).
"""

from __future__ import annotations

import math
import os
import shutil
from dataclasses import dataclass, field

import numpy as np
import pytest

from petakit5d.stitch import block as SB
from petakit5d.stitch import blockinfo as B
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, InPath, Tree, Value, Work, mcell
from petakit5d.testing.compare import MCell, mchar
from petakit5d.testing.synth import CounterRNG, write_zarr

# (the layouts and their MATLAB values as in test_blockinfo.py; a test module can't import
# another, whose oracle cases would be registered twice)


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


class _TB:
    layout_matlab = staticmethod(layout_matlab)
    run_block_info = staticmethod(run_block_info)


TB = _TB
METHODS = ("none", "max", "mean", "median", "feather")
CLASSES = ("uint16", "uint8", "single")

# ---------------------------------------------------------------------------------------------
# scenarios: layouts, tiles and distance maps
# ---------------------------------------------------------------------------------------------


@dataclass
class Scene:
    layout: Layout
    chunk: tuple[int, int, int] = (4, 8, 8)
    maps: str = "per_tile"  # 'per_tile', 'shared' or 'pooled'
    holes: bool = True  # zero regions in tiles 2 and 3
    pooled_values: tuple = ()  # per tile (low, high) of the pooled maps
    major_hole: tuple | None = None  # a box of zeros in tile 1 (lf_major_holes)
    pattern: str = "faces"  # pooled map values: by the distance to the faces, or 'cyclic'
    fine: bool = False  # single tiles with fractions that aren't binary (f * w / w != f)


def _L(sizes, shifts, nv, bs, **kw) -> Layout:
    return Layout(np.asarray(sizes), np.asarray(shifts), tuple(nv), tuple(bs), **kw)


SCENES: dict[str, Scene] = {
    "grid2x2": Scene(
        _L(
            [[8, 20, 22]] * 4,
            [[0, 0, 0], [16, 1, 0], [1, 15, 0], [17, 16, 1]],
            (9, 36, 39),
            (5, 12, 13),
            task_size=4,
        ),  # fmt: skip
        maps="shared",
    ),
    "row3": Scene(
        _L(
            [[6, 14, 24]] * 3,
            [[0, 0, 0], [10, 0, 0], [20, 0, 0]],
            (6, 14, 44),
            (6, 14, 44),
            task_size=1,
        ),  # fmt: skip
    ),
    "row3_batches": Scene(
        _L(
            [[6, 14, 24]] * 3,
            [[0, 0, 0], [10, 0, 0], [20, 0, 0]],
            (6, 14, 44),
            (3, 7, 15),
            task_size=5,
        ),  # fmt: skip
    ),
    "partial": Scene(
        _L([[6, 16, 16], [6, 16, 16]], [[0, 0, 0], [13, 0, 0]], (6, 16, 29), (6, 16, 29)),
    ),
    "two_d": Scene(
        _L(
            [[1, 18, 20]] * 4,
            [[0, 0, 0], [15, 1, 0], [0, 14, 0], [16, 15, 0]],
            (1, 33, 36),
            (1, 10, 11),
            stitch_2d=True,
            task_size=3,
        ),  # fmt: skip
        chunk=(1, 8, 8),
    ),
    "z_stack": Scene(
        _L([[10, 12, 12], [10, 12, 12]], [[0, 0, 0], [0, 0, 7]], (17, 12, 12), (8, 12, 12)),
        maps="shared",
    ),
    "single": Scene(_L([[6, 14, 17]], [[0, 0, 0]], (6, 14, 17), (6, 14, 17))),
    "single_batches": Scene(_L([[6, 14, 17]], [[1, 0, 0]], (6, 14, 18), (4, 8, 8), task_size=3)),
}
# the large-file feather blending: maps pooled in z (poolSize [1 1 2]) with these values (the
# major tiles 8 to 9, the minor ones 1 to 3: 8 > 1000^(1/10) * 3)
LF_SCENES: dict[str, Scene] = {
    "lf_two": Scene(
        _L([[8, 14, 16], [8, 14, 16]], [[0, 0, 0], [11, 0, 0]], (8, 14, 27), (8, 14, 14)),
        maps="pooled",
        pooled_values=((1, 3), (1, 3)),
    ),
    "lf_major_valid": Scene(
        _L(
            [[8, 20, 24], [8, 8, 8], [8, 8, 8]],
            [[0, 0, 0], [4, 4, 0], [12, 9, 0]],
            (8, 20, 24),
            (8, 20, 24),
        ),  # fmt: skip
        maps="pooled",
        holes=False,
        pooled_values=((8, 9), (1, 3), (1, 2)),
    ),
    "lf_major_holes": Scene(
        _L(
            [[8, 20, 24], [8, 8, 8], [8, 8, 8]],
            [[0, 0, 0], [4, 4, 0], [12, 9, 0]],
            (8, 20, 24),
            (8, 20, 24),
        ),  # fmt: skip
        maps="pooled",
        holes=False,
        pooled_values=((8, 9), (1, 3), (1, 2)),
        # (the minor tiles are read in the major tile's hole, in the batch's coordinates)
        major_hole=(2, 9, 9, 6, 12, 12),
    ),
    # one minor tile (tile 2) with two major ones: its weight is constant where it is read
    "lf_one_minor": Scene(
        _L(
            [[8, 20, 24], [8, 8, 8], [8, 10, 8]],
            [[0, 0, 0], [4, 4, 0], [16, 0, 0]],
            (8, 20, 24),
            (8, 20, 24),
        ),  # fmt: skip
        maps="pooled",
        holes=False,
        pooled_values=((8, 9), (1, 3), (8, 9)),
        major_hole=(2, 7, 9, 6, 10, 12),
    ),
    "lf_minor_out": Scene(
        _L(
            [[8, 16, 20], [8, 8, 8], [8, 8, 8]],
            [[0, 0, 0], [3, 3, 0], [15, 10, 0]],
            (8, 18, 23),
            (8, 18, 23),
        ),  # fmt: skip
        maps="pooled",
        holes=False,
        pooled_values=((8, 9), (1, 3), (1, 2)),
    ),
    "lf_two_major": Scene(
        _L(
            [[8, 20, 16], [8, 20, 16], [8, 8, 8]],
            [[0, 0, 0], [6, 0, 0], [7, 6, 0]],
            (8, 20, 22),
            (8, 20, 22),
        ),  # fmt: skip
        maps="pooled",
        holes=False,
        pooled_values=((8, 9), (8, 9), (1, 2)),
    ),
    "lf_full": Scene(
        _L(
            [[8, 18, 20], [8, 18, 20], [8, 18, 20]],
            [[0, 0, 0], [1, 1, 0], [2, 0, 0]],
            (8, 19, 22),
            (8, 19, 22),
        ),  # fmt: skip
        maps="pooled",
        holes=False,
        pooled_values=((8, 9), (1, 2), (8, 9)),
    ),
    # the minor tiles are decided by the smallest weight of the full tiles (5) against each
    # tile's largest (cyclic maps): the two partial tiles (largest 3 and 9) stay major, though
    # 9 (the full tile's largest) or their smallest (1) would make them minor
    "lf_ranges": Scene(
        _L(
            [[8, 20, 24], [8, 8, 8], [8, 8, 8]],
            [[0, 0, 0], [4, 4, 0], [12, 9, 0]],
            (8, 20, 24),
            (8, 20, 24),
        ),  # fmt: skip
        maps="pooled",
        holes=False,
        pattern="cyclic",
        pooled_values=((5, 9), (1, 3), (1, 9)),
    ),
}
# an overlap region of exactly 80% of its batch (z 4 of 5): full, by MATLAB's >= 0.8
SCENES_EXTRA: dict[str, Scene] = {
    "exact80": Scene(
        _L([[5, 10, 10], [4, 10, 10]], [[0, 0, 0], [0, 0, 1]], (5, 10, 10), (5, 10, 10)),
        fine=True,
    ),
}
POOL = (2, 1, 1)  # (z, y, x): MATLAB's poolSize [1 1 2]


def _ground_truth(shape) -> np.ndarray:
    return CounterRNG(7, "s07_gt").integers(100, 3000, tuple(shape))


def tile_data(scene: Scene, k: int, cls: str) -> np.ndarray:
    """tile k (0-based) of a scene in class cls: the ground truth at its place, offset by tile
    (so that the blends show), with zero regions in tiles 1 and 2 (no data)"""
    L = scene.layout
    sz = [int(v) for v in L.im_sizes[k]]
    x0, y0, z0 = (int(v) for v in L.shift[k])
    gt = _ground_truth(
        [max(int(L.shift[i, 2]) + int(L.im_sizes[i, 0]) for i in range(L.n)),
         max(int(L.shift[i, 1]) + int(L.im_sizes[i, 1]) for i in range(L.n)),
         max(int(L.shift[i, 0]) + int(L.im_sizes[i, 2]) for i in range(L.n))]
    )  # fmt: skip
    g = gt[z0 : z0 + sz[0], y0 : y0 + sz[1], x0 : x0 + sz[2]]
    if cls == "uint16":
        v = (g + 1000 * k).astype(np.uint16)
    elif cls == "uint8":
        v = (g % 200 + 20 + 10 * (k % 4)).astype(np.uint8)
    else:
        step = 0.1 if scene.fine else 0.25
        v = ((g - 1500) * step + 0.5 * k).astype(np.float32)
    if scene.holes and k == 1:
        v[:, :, :2] = 0
    if scene.holes and k == 2 and v.shape[0] > 2:
        v[-1] = 0
    if scene.major_hole is not None and k == 0:
        h = scene.major_hole
        v[h[0] : h[3], h[1] : h[4], h[2] : h[5]] = 0
    return np.ascontiguousarray(v)


def _face_distance(n: int) -> np.ndarray:
    i = np.arange(n)
    return np.minimum(i + 1, n - i)


def dist_map(shape) -> np.ndarray:
    """a feather weight map (as compute_tile_bwdist's, already raised to its power): integers,
    largest inside the tile"""
    nz, ny, nx = shape
    dy, dx, dz = _face_distance(ny), _face_distance(nx), _face_distance(nz)
    w = (1 + np.minimum(dy[:, None], dx[None, :]))[None] * dz[:, None, None]
    return np.ascontiguousarray(w.astype(np.float32))


def pooled_map(shape, lo_hi, pattern: str = "faces") -> np.ndarray:
    """a pooled map (z pooled by 2) with values lo to hi: growing from the faces, or cyclic in
    y + x (every region of the size of the range has both ends)"""
    nz, ny, nx = shape
    lo, hi = lo_hi
    if pattern == "cyclic":
        w = lo + (np.arange(ny)[:, None] + np.arange(nx)[None, :]) % (hi - lo + 1)
        return np.ascontiguousarray(
            np.broadcast_to(w, (math.ceil(nz / 2), ny, nx)).astype(np.float32)
        )
    dy, dx = _face_distance(ny), _face_distance(nx)
    m = np.minimum(dy[:, None], dx[None, :])
    w = lo + np.minimum(m - 1, hi - lo)
    return np.ascontiguousarray(np.broadcast_to(w, (math.ceil(nz / 2), ny, nx)).astype(np.float32))


def scene_maps(scene: Scene) -> tuple[list[np.ndarray], list[int]]:
    L = scene.layout
    if scene.maps == "shared":
        return [dist_map(L.im_sizes[0])], [0]
    if scene.maps == "pooled":
        return [
            pooled_map(L.im_sizes[k], scene.pooled_values[k], scene.pattern) for k in range(L.n)
        ], list(range(L.n))
    return [dist_map(L.im_sizes[k]) for k in range(L.n)], list(range(L.n))


def make_root(scene: Scene, cls: str):
    """root/tiles/tile_<k>.zarr and root/imdist/d<k>.zarr (MATLAB's layout)"""

    def make(folder):
        folder = os.fspath(folder)
        os.makedirs(os.path.join(folder, "tiles"))
        os.makedirs(os.path.join(folder, "imdist"))
        for k in range(scene.layout.n):
            write_zarr(
                os.path.join(folder, "tiles", f"tile_{k + 1:02d}.zarr"),
                tile_data(scene, k, cls),
                (4, 8, 8),
                "matlab",
            )
        maps, _ = scene_maps(scene)
        for k, m in enumerate(maps):
            write_zarr(os.path.join(folder, "imdist", f"d{k + 1}.zarr"), m, (4, 8, 8), "matlab")

    return make


@dataclass
class Blend:
    scene: Scene
    method: str
    cls: str  # the tiles' class
    dtype: str  # the result's class
    save_mip: bool = True
    task_size: int | None = None
    border: tuple[int, int, int] | None = None

    @property
    def layout(self) -> Layout:
        L = self.scene.layout
        if self.task_size is None and self.border is None:
            return L
        return Layout(
            L.im_sizes,
            L.shift,
            L.nv_size,
            L.batch_size,
            border=L.border if self.border is None else self.border,
            task_size=L.task_size if self.task_size is None else self.task_size,
            max_files=L.max_files,
            stitch_2d=L.stitch_2d,
            overlap_type=L.overlap_type,
        )


def _imdist(b: Blend, root: str = "root/") -> tuple[list[str], list[int]]:
    if b.method.lower() != "feather":
        return [], []
    maps, idx = scene_maps(b.scene)
    return [f"{root}imdist/d{k + 1}.zarr" for k in range(len(maps))], idx


def blend_matlab(b: Blend) -> dict:
    """the case's struct for MATLAB: the layout's block-information arguments and the blending
    options (tiles and maps relative to the root folder)"""
    from petakit5d.core.axes import to_yxz

    L = b.layout
    d = TB.layout_matlab(L, "tiles/")
    paths, idx = _imdist(b, "")
    if not idx:
        file_idx = np.zeros((0, 0))
    elif len(paths) == 1:
        file_idx = np.asarray([[1.0]])
    else:
        file_idx = np.asarray(idx, float).reshape(-1, 1) + 1
    d.update(
        chunkSize=np.asarray(to_yxz(b.scene.chunk), dtype=np.float64),
        dtype=b.dtype,
        method=b.method,
        imdist=MCell((1, len(paths)), [mchar(p) for p in paths]),
        imdistFileIdx=file_idx,
        poolSize=np.asarray([1.0, 1.0, 2.0]) if b.scene.maps == "pooled" else np.zeros((0, 0)),
        wd=10.0,
        saveMIP=bool(b.save_mip),
    )
    return d


def _num_batches(L: Layout) -> int:
    return math.prod(-(-a // c) for a, c in zip(L.nv_size, L.batch_size, strict=True))


def run_blend(b: Blend, *, base: str = "", num_threads=None) -> list:
    """the Python side of a blending case in the folder base (with its root/ of tiles and maps;
    '': the current folder): the block information, the stitched zarr, every task through the
    MATLAB-name layer, the merged MIP"""
    from petakit5d.io import create_zarr

    L = b.layout
    os.makedirs(base + "stitchInfo")
    bf, pf, _ = TB.run_block_info(L, base + "stitchInfo", tile_prefix=base + "root/tiles/")
    create_zarr(base + "vol.zarr", data_size=L.nv_size, block_size=b.scene.chunk, dtype=b.dtype)
    os.makedirs(base + "flags")
    return blend_tasks(b, bf, pf, base=base, num_threads=num_threads)


def blend_tasks(
    b: Blend, bf: str, pf: list[str], *, base: str = "", tile_fns=None, num_threads=None
) -> list:
    """every blending task (processStitchBlock through the MATLAB-name layer), then the merged
    MIP; tile_fns: the tiles' paths given (bf is then '')"""
    from petakit5d import matlab
    from petakit5d.core.axes import to_yxz

    L = b.layout
    paths, idx = _imdist(b, base + "root/")
    nb = _num_batches(L)
    ts = L.task_size
    flags = []
    for t, p in enumerate(pf):
        bi = list(range(t * ts + 1, min((t + 1) * ts, nb) + 1))
        flags.append(f"{base}flags/blocks_{bi[0]}_{bi[-1]}.mat")
        idx1 = [] if not idx else ([1] if len(paths) == 1 else [i + 1 for i in idx])
        opts = [
            "imSize", list(to_yxz(L.nv_size)), "batchSize", list(to_yxz(L.batch_size)),
            "dtype", b.dtype, "BlendMethod", b.method, "BorderSize", list(to_yxz(L.border)),
            "imdistFullpaths", paths, "imdistFileIdx", idx1,
            "poolSize", [1, 1, 2] if b.scene.maps == "pooled" else [], "weightDegree", 10,
            "saveMIP", bool(b.save_mip),
        ]  # fmt: skip
        if num_threads is not None:
            opts += ["numThreads", num_threads]
        fns = [] if tile_fns is None else list(tile_fns)
        matlab.processStitchBlock(bi, bf, p, flags[-1], base + "vol.zarr", [], fns, *opts)
    return [SB.merge_block_mips(flags, L.nv_size, b.dtype)]


# (absolute paths, as the frame gives them: MATLAB's createzarr fails on a relative path)
_BLEND_CODE = (
    "L = in{1}; d = in{2}; w = in{3}; si = [w, 'stitchInfo']; mkdir(si); "
    "tileFns = cellfun(@(x) [d, x], L.tileFns, 'unif', 0); "
    "imdist = cellfun(@(x) [d, x], L.imdist, 'unif', 0); "
    "[bf, pf] = stitch_process_block_info(L.int_xyz_shift, L.imSizes, L.nvSize, L.blockSize, "
    "L.data_order_mat, L.overlap_matrix, L.ol_region_cell, L.half_ol_region_cell, "
    "L.overlap_map_mat, L.BorderSize, tileFns, si, 'vol', true, 'stitchInfoFullpath', "
    "[si, '/vol.mat'], 'stitch2D', L.stitch2D, 'taskSize', L.taskSize, "
    "'maxFileNumPerFolder', L.maxFileNumPerFolder, 'uuid', 'u0', 'parseCluster', false); "
    "nv = [w, 'vol.zarr']; createzarr(nv, dataSize=L.nvSize, blockSize=L.chunkSize, "
    "dtype=L.dtype); mkdir([w, 'flags']); nb = prod(ceil(L.nvSize ./ L.blockSize)); "
    "flags = cell(numel(pf), 1); "
    "for t = 1 : numel(pf), bi = (t - 1) * L.taskSize + 1 : min(t * L.taskSize, nb); "
    "flags{t} = sprintf('%sflags/blocks_%d_%d.mat', w, bi(1), bi(end)); "
    "processStitchBlock(bi, bf, pf{t}, flags{t}, nv, [], [], 'imSize', L.nvSize, "
    "'batchSize', L.blockSize, 'dtype', L.dtype, 'BlendMethod', L.method, "
    "'BorderSize', L.BorderSize, 'imdistFullpaths', imdist, 'imdistFileIdx', "
    "L.imdistFileIdx, 'poolSize', L.poolSize, 'weightDegree', L.wd, 'saveMIP', L.saveMIP); "
    "end; out = {stitch_merge_block_mips(flags, L.nvSize, L.dtype)};"
)

BLENDS: dict[str, Blend] = {}
for _s in SCENES:
    for _m in METHODS:
        for _c in CLASSES:
            BLENDS[f"blend.{_s}.{_m}.{_c}"] = Blend(SCENES[_s], _m, _c, _c)
for _s in LF_SCENES:
    for _c in CLASSES:
        BLENDS[f"lf.{_s}.{_c}"] = Blend(LF_SCENES[_s], "feather", _c, _c)
BLENDS["blend.cast.grid2x2.mean.uint16_to_single"] = Blend(
    SCENES["grid2x2"], "mean", "uint16", "single"
)
BLENDS["blend.cast.grid2x2.feather.single_to_uint16"] = Blend(
    SCENES["grid2x2"], "feather", "single", "uint16"
)
BLENDS["blend.cast.row3.none.uint16_to_uint8"] = Blend(
    SCENES["row3_batches"], "none", "uint16", "uint8"
)
BLENDS["blend.nomip.grid2x2.none.uint16"] = Blend(
    SCENES["grid2x2"], "none", "uint16", "uint16", save_mip=False
)
BLENDS["blend.task1.grid2x2.feather.uint16"] = Blend(
    SCENES["grid2x2"], "feather", "uint16", "uint16", task_size=1
)
for _c in ("uint16", "single"):
    BLENDS[f"blend.exact80.feather.{_c}"] = Blend(SCENES_EXTRA["exact80"], "feather", _c, _c)
BLENDS["blend.border.grid2x2.none"] = Blend(
    SCENES["grid2x2"], "none", "uint16", "uint16", border=(1, 1, 1)
)
BLENDS["blend.Feather.partial"] = Blend(SCENES["partial"], "Feather", "uint16", "uint16")


def p_blend(name: str, work: str):
    return run_blend(BLENDS[name], base=work)[0]


def _blend_case(name: str, b: Blend, **kw) -> Case:
    return Case(
        name,
        code=_BLEND_CODE,
        args=[In("blend"), InPath("root"), Work()],
        nargout=1,
        python=p_blend,
        py_args=[name, Work()],
        inputs={
            "blend": Value(lambda b=b: blend_matlab(b)),
            "root": Tree(make_root(b.scene, b.cls), "root"),
        },
        **kw,
    )


CASES: list[Case] = [_blend_case(n, b) for n, b in BLENDS.items()]

# --- stitch_merge_block_mips: flags as processStitchBlock writes them, merged
_MERGE_CODE = (
    "M = in{1}; flags = {}; mkdir('flags'); for t = 1 : numel(M.tasks), T = M.tasks{t}; "
    "batchInds = T.batchInds; mipBboxes = T.mipBboxes; mips = T.mips; "
    "flags{end + 1} = sprintf('flags/blocks_%d_%d.mat', batchInds(1), batchInds(end)); "
    "if T.empty, fclose(fopen(flags{end}, 'w')); else, "
    "save(flags{end}, 'batchInds', 'mipBboxes', 'mips', '-v6'); end; end; "
    "out = {stitch_merge_block_mips(flags, M.imSize, M.dtype)};"
)


def merge_flags(im: np.ndarray, batch: tuple[int, int, int], task_size: int, empty=()) -> list:
    """the flags of tasks of task_size batches of im ((z, y, x)) in batches of batch (z, y, x),
    as StitchMergeBlockMipsTest's writeFlags: per batch its z MIP and [ymin xmin ymax xmax]"""
    from petakit5d.matlabcompat.reduce import max as mmax

    sub = [-(-a // c) for a, c in zip(im.shape, batch, strict=True)]
    n = math.prod(sub)
    tasks = []
    for s in range(0, n, task_size):
        inds = list(range(s, min(s + task_size, n)))
        mips, boxes = [], []
        for bi in inds:
            sz, sy, sx = SB._batch_sub(bi, sub)
            o0 = [sz * batch[0], sy * batch[1], sx * batch[2]]
            o1 = [min(o0[q] + batch[q], im.shape[q]) for q in range(3)]
            blk = im[o0[0] : o1[0], o0[1] : o1[1], o0[2] : o1[2]]
            mips.append(np.ascontiguousarray(mmax(blk, axis=0)[0]))
            boxes.append([o0[1] + 1, o0[2] + 1, o1[1], o1[2]])
        tasks.append(
            {
                "batchInds": np.asarray([i + 1 for i in inds], dtype=np.float64),
                "mipBboxes": np.asarray(boxes, dtype=np.float64),
                "mips": mips,
                "empty": len(tasks) in empty,
            }
        )
    return tasks


def _merge_im(cls: str) -> np.ndarray:
    rng = CounterRNG(11, "s07_merge")
    if cls == "uint16":
        return rng.integers(1, 60000, (6, 8, 10)).astype(np.uint16)
    a = (rng.normal((6, 8, 10)) - 5).astype(np.float32)
    a[:, 1, 2] = np.float32("nan")
    a[1, 4, 6] = np.float32("nan")
    return a


MERGES = {
    "merge.uint16": ("uint16", (3, 4, 5), 3, ()),
    "merge.single_negative_nan": ("single", (3, 4, 5), 2, ()),
    "merge.flag_without_mips": ("uint16", (3, 4, 5), 3, (1,)),
    "merge.one_task": ("uint16", (6, 8, 10), 1, ()),
}


def _merge_value(name: str) -> dict:
    from petakit5d.testing.compare import to_mvalue

    cls, batch, ts, empty = MERGES[name]
    tasks = merge_flags(_merge_im(cls), batch, ts, empty)
    items = []
    for t in tasks:
        mips = MCell((len(t["mips"]), 1), [to_mvalue(m) for m in t["mips"]])
        items.append(to_mvalue({**t, "mips": mips}))
    return {"tasks": MCell((1, len(items)), items), "imSize": [8.0, 10.0, 6.0], "dtype": cls}


def _mips_cell(mips: list) -> np.ndarray:
    cell = np.empty((len(mips), 1), dtype=object)
    for i, m in enumerate(mips):
        cell[i, 0] = m
    return cell


def p_merge(name: str):
    from petakit5d.io import write_mat

    cls, batch, ts, empty = MERGES[name]
    im = _merge_im(cls)
    os.makedirs("flags")
    flags = []
    for t in merge_flags(im, batch, ts, empty):
        b = [int(v) for v in t["batchInds"]]
        fn = f"flags/blocks_{b[0]}_{b[-1]}.mat"
        flags.append(fn)
        if t["empty"]:
            open(fn, "w").close()
            continue
        data = {
            "batchInds": t["batchInds"],
            "mipBboxes": t["mipBboxes"],
            "mips": _mips_cell(t["mips"]),
        }
        write_mat(fn, data, version="6", compression=False)
    return SB.merge_block_mips(flags, im.shape, cls)


CASES += [
    Case(
        n,
        code=_MERGE_CODE,
        args=[In("m")],
        python=p_merge,
        py_args=[n],
        inputs={"m": Value(lambda n=n: _merge_value(n))},
    )  # fmt: skip
    for n in MERGES
]

# --- processStitchSeparteTiles: each tile padded onto the canvas (XR_crop_zarr tasks)
_SEPARATE_CODE = (
    "S = in{1}; d = in{2}; o = [in{3}, 'out']; mkdir(o); "
    "processStitchSeparteTiles(cellfun(@(x) [d, x], S.tiles, 'unif', 0), o, S.st, S.pImSz, "
    "'batchSize', S.batch, 'blockSize', S.block, 'parseCluster', false, 'uuid', 'u0', "
    "'saveMIP', S.saveMIP); out = {};"
)
SEPARATE = {
    "separate.grid2x2": ("grid2x2", (8, 16, 16), (4, 8, 8), True),
    "separate.two_d": ("two_d", (1, 16, 16), (1, 8, 8), False),
}


def _separate_value(name: str) -> dict:
    from petakit5d.core.axes import to_yxz

    s, batch, block, mip = SEPARATE[name]
    L = SCENES[s].layout
    return {
        "tiles": MCell((L.n, 1), [mchar(f) for f in L.tile_fns("tiles/")]),
        "st": np.asarray(L.shift, dtype=np.float64),
        "pImSz": np.asarray(to_yxz(L.nv_size), dtype=np.float64),
        "batch": np.asarray(to_yxz(batch), dtype=np.float64),
        "block": np.asarray(to_yxz(block), dtype=np.float64),
        "saveMIP": mip,
    }


def p_separate(name: str, work: str) -> list:
    s, batch, block, mip = SEPARATE[name]
    L = SCENES[s].layout
    os.makedirs(work + "out")
    SB.process_separate_tiles(
        L.tile_fns(work + "root/tiles/"), work + "out", L.shift, L.nv_size, batch_size=batch,
        block_size=block, parse_cluster=False, uuid="u0", save_mip=mip,
    )  # fmt: skip
    return []


CASES += [
    Case(
        n,
        code=_SEPARATE_CODE,
        args=[In("s"), InPath("root"), Work()],
        nargout=0,
        python=p_separate,
        py_args=[n, Work()],
        inputs={
            "s": Value(lambda n=n: _separate_value(n)),
            "root": Tree(make_root(SCENES[SEPARATE[n][0]], "uint16"), "root"),
        },
    )  # fmt: skip
    for n in SEPARATE
]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# blending from MATLAB's block information (the golden store's files)
# ---------------------------------------------------------------------------------------------

_FROM_MATLAB = [n for n in BLENDS if n.endswith(".uint16") and ".single_batches." not in n][:40]


@pytest.mark.parametrize("name", _FROM_MATLAB)
def test_blending_reads_matlabs_block_information(name, tmp_path):
    """Python's blending tasks on MATLAB's per-block JSON files (the golden store's) give
    MATLAB's stitched volume and flags (MATLAB's block information file names MATLAB's tile
    paths, so the tiles are given)"""
    from petakit5d.io import create_zarr, read_mat, read_zarr

    store = oracle.Store.from_env()
    if store is None:
        pytest.skip("needs the golden store (PETAKIT5D_GOLDEN_DIR): MATLAB's files")
    case = next(c for c in CASES if c.id == name)
    key = case.key()
    if not store.has(key):
        pytest.skip(f"no reference for {name} in the store")
    files = store.case_dir(key) / "files"
    b = BLENDS[name]
    L = b.layout
    base = f"{tmp_path}/"
    make_root(b.scene, b.cls)(tmp_path / "root")
    shutil.copytree(files / "stitchInfo", tmp_path / "stitchInfo")
    create_zarr(base + "vol.zarr", data_size=L.nv_size, block_size=b.scene.chunk, dtype=b.dtype)
    os.makedirs(base + "flags")
    _, pf = _json_paths(L, base)
    mip = blend_tasks(b, "", pf, base=base, tile_fns=L.tile_fns(base + "root/tiles/"))[0]
    ref = np.asarray(read_zarr(str(files / "vol.zarr")))
    got = np.asarray(read_zarr(base + "vol.zarr"))
    assert got.dtype == ref.dtype and np.array_equal(got.view(np.uint8), ref.view(np.uint8))
    for fn in sorted(os.listdir(base + "flags")):
        a, r = read_mat(f"{base}flags/{fn}"), read_mat(str(files / "flags" / fn))
        assert np.array_equal(np.asarray(a["mipBboxes"]), np.asarray(r["mipBboxes"]))
        ma = list(np.asarray(a["mips"], dtype=object).ravel(order="F"))
        mr = list(np.asarray(r["mips"], dtype=object).ravel(order="F"))
        for x, y in zip(ma, mr, strict=True):
            assert np.array_equal(np.asarray(x), np.asarray(y))
    assert mip is not None


def _json_paths(L, base: str = ""):
    from petakit5d.stitch.blockinfo import per_block_info_paths

    return per_block_info_paths(
        f"{base}stitchInfo/vol_task_size_{L.task_size}/", L.task_size, _num_batches(L), L.max_files
    )


# ---------------------------------------------------------------------------------------------
# StitchMergeBlockMipsTest
# ---------------------------------------------------------------------------------------------


def _write_flags(folder, im, batch, task_size, empty=()) -> list[str]:
    fns = []
    for t in merge_flags(im, batch, task_size, empty):
        b = [int(v) - 1 for v in t["batchInds"]]
        fn = os.path.join(folder, f"blocks_{b[0] + 1}_{b[-1] + 1}.mat")
        fns.append(fn)
        if t["empty"]:
            open(fn, "w").close()
            continue
        boxes = [(int(r[0]) - 1, int(r[1]) - 1, int(r[2]), int(r[3])) for r in t["mipBboxes"]]
        SB.save_block_flag(fn, b, boxes, t["mips"])
    return fns


def test_merged_mip_equals_mip_of_the_image(tmp_path):
    im = _merge_im("uint16")
    fns = _write_flags(tmp_path, im, (3, 4, 5), 3)
    assert np.array_equal(SB.merge_block_mips(fns, im.shape, "uint16"), im.max(axis=0))


def test_single_keeps_negatives_and_nan(tmp_path):
    from petakit5d.matlabcompat.reduce import max as mmax

    im = _merge_im("single")
    fns = _write_flags(tmp_path, im, (3, 4, 5), 2)
    got = SB.merge_block_mips(fns, im.shape, "single")
    exp = np.asarray(mmax(im, axis=0)[0])
    assert got.dtype == np.float32
    assert np.array_equal(np.isnan(got), np.isnan(exp))
    assert np.array_equal(got[~np.isnan(got)], exp[~np.isnan(exp)])
    assert np.isnan(got[1, 2]) and (got < 0).any()


def test_flag_without_mips_gives_no_mip(tmp_path):
    im = _merge_im("uint16")
    fns = _write_flags(tmp_path, im, (3, 4, 5), 3, empty=(1,))
    assert SB.merge_block_mips(fns, im.shape, "uint16") is None


def test_flag_merges_the_mips_of_batches_at_one_place(tmp_path):
    # batches 1 and 3 share their y and x: one MIP for both, the max in batch order
    m1 = np.array([[1, 5], [0, 2]], np.uint16)
    m3 = np.array([[4, 1], [3, 2]], np.uint16)
    m2 = np.array([[7, 7]], np.uint16)
    fn = str(tmp_path / "blocks_1_3.mat")
    SB.save_block_flag(fn, [0, 1, 2], [(0, 0, 2, 2), (2, 0, 3, 2), (0, 0, 2, 2)], [m1, m2, m3])
    from petakit5d.io import read_mat

    a = read_mat(fn)
    assert np.array_equal(np.asarray(a["batchInds"]).ravel(), [1, 2, 3])
    assert np.array_equal(np.asarray(a["mipBboxes"]), [[1, 1, 2, 2], [3, 1, 3, 2]])
    mips = list(np.asarray(a["mips"], dtype=object).ravel(order="F"))
    assert np.array_equal(np.asarray(mips[0]), np.maximum(m1, m3))
    assert np.array_equal(np.asarray(mips[1]).reshape(1, 2), m2)


# ---------------------------------------------------------------------------------------------
# StitchingTest's blending tests, at the level of the blending tasks
# ---------------------------------------------------------------------------------------------


def _stitch(tmp_path, monkeypatch, b: Blend, **kw) -> np.ndarray:
    from petakit5d.io import read_zarr

    make_root(b.scene, b.cls)(tmp_path / "root")
    monkeypatch.chdir(tmp_path)
    run_blend(b, **kw)
    return np.asarray(read_zarr("vol.zarr")).reshape(b.layout.nv_size)


def _two_tiles(overlap_type: str, offset: int) -> Scene:
    # (split along x: MATLAB's halfOrder [1, 2, 3])
    L = Layout(np.asarray([[4, 10, 40]] * 2), np.asarray([[0, 0, 0], [20, 0, 0]]), (4, 10, 60),
               (4, 10, 60), overlap_type=overlap_type, half_order=(2, 1, 0))  # fmt: skip
    return Scene(L, holes=False)


def _constant_root(tmp_path, scene: Scene, values) -> None:
    os.makedirs(tmp_path / "root" / "tiles")
    os.makedirs(tmp_path / "root" / "imdist")
    for k, v in enumerate(values):
        tile = np.full(scene.layout.im_sizes[k], v, np.uint16)
        write_zarr(
            tmp_path / "root" / "tiles" / f"tile_{k + 1:02d}.zarr", tile, (4, 8, 8), "matlab"
        )
    for k, m in enumerate(scene_maps(scene)[0]):
        write_zarr(tmp_path / "root" / "imdist" / f"d{k + 1}.zarr", m, (4, 8, 8), "matlab")


def test_none_blend_takes_half_overlap_from_each_tile(tmp_path, monkeypatch):
    from petakit5d.io import read_zarr

    sc = _two_tiles("none", 0)
    _constant_root(tmp_path, sc, (100, 200))
    monkeypatch.chdir(tmp_path)
    run_blend(Blend(sc, "none", "uint16", "uint16"))
    v = np.asarray(read_zarr("vol.zarr"))
    # overlap x 20..39: the first half from tile 1, the second from tile 2
    assert (v[:, :, :30] == 100).all() and (v[:, :, 30:] == 200).all()


@pytest.mark.parametrize("method", ["max", "mean", "median"])
def test_blending_splits_the_overlap_in_quarters(method, tmp_path, monkeypatch):
    from petakit5d.io import read_zarr

    sc = _two_tiles("half", 0)
    _constant_root(tmp_path, sc, (100, 5100))
    monkeypatch.chdir(tmp_path)
    run_blend(Blend(sc, method, "uint16", "uint16"))
    v = np.asarray(read_zarr("vol.zarr"))
    mid = {"max": 5100, "mean": 2600, "median": 2600}[method]
    # overlap x 20..39 of 20 voxels: 5 from each tile, the middle 10 blended
    assert (v[:, :, :25] == 100).all() and (v[:, :, 35:] == 5100).all()
    assert (v[:, :, 25:35] == mid).all()


@pytest.mark.parametrize("method", METHODS)
def test_exact_positions_reproduce_ground_truth(method, tmp_path, monkeypatch):
    sc = SCENES["grid2x2"]
    L = sc.layout
    gt = _ground_truth(L.nv_size).astype(np.uint16)

    def make(folder):
        os.makedirs(folder / "tiles")
        os.makedirs(folder / "imdist")
        for k in range(L.n):
            x0, y0, z0 = (int(v) for v in L.shift[k])
            nz, ny, nx = (int(v) for v in L.im_sizes[k])
            write_zarr(folder / "tiles" / f"tile_{k + 1:02d}.zarr",
                       np.ascontiguousarray(gt[z0:z0 + nz, y0:y0 + ny, x0:x0 + nx]), (4, 8, 8),
                       "matlab")  # fmt: skip
        write_zarr(folder / "imdist" / "d1.zarr", dist_map(L.im_sizes[0]), (4, 8, 8), "matlab")

    make(tmp_path / "root")
    monkeypatch.chdir(tmp_path)
    run_blend(Blend(sc, method, "uint16", "uint16"))
    from petakit5d.io import read_zarr

    v = np.asarray(read_zarr("vol.zarr"))
    covered = np.zeros(L.nv_size, bool)
    for k in range(L.n):
        x0, y0, z0 = (int(v) for v in L.shift[k])
        nz, ny, nx = (int(v) for v in L.im_sizes[k])
        covered[z0 : z0 + nz, y0 : y0 + ny, x0 : x0 + nx] = True
    assert np.array_equal(v[covered], gt[covered])
    assert (v[~covered] == 0).all()


def _per_voxel(b: Blend, tiles: list[np.ndarray]) -> np.ndarray:
    """the 'none', 'max', 'mean', 'median' blend of every voxel by its definition: the tiles'
    nonzero values outside their half-overlap regions, or where there are none all of them"""
    from petakit5d.stitch.geometry import tile_canvas_boxes

    L = b.layout
    nz, ny, nx = L.nv_size
    cb, tb, valid = tile_canvas_boxes(L.shift, L.im_sizes, L.nv_size, L.data_order, L.stitch_2d)
    vals = np.zeros((L.n, nz, ny, nx), np.float32)
    masked = np.zeros((L.n, nz, ny, nx), bool)
    for k in range(L.n):
        c, t = cb[k], tb[k]
        vals[k, c[0] : c[3], c[1] : c[4], c[2] : c[5]] = tiles[k][
            t[0] : t[3], t[1] : t[4], t[2] : t[5]
        ]
    for p, (i, j) in enumerate(L.pairs):
        for k, r in ((int(i), L.half[p][0]), (int(j), L.half[p][1])):
            off = [int(cb[k][q]) - int(tb[k][q]) for q in range(3)]
            lo = [max(r[q], tb[k][q]) + off[q] for q in range(3)]
            hi = [min(r[q + 3], tb[k][q + 3]) + off[q] for q in range(3)]
            if all(hi[q] > lo[q] for q in range(3)):
                masked[k, lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] = True  # fmt: skip
    out = np.zeros((nz, ny, nx), np.float32)
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                v = vals[:, z, y, x]
                sel = (v != 0) & ~masked[:, z, y, x]
                if not sel.any():
                    sel = v != 0
                if not sel.any():
                    continue
                s = v[sel]
                if b.method == "none":
                    out[z, y, x] = s[0]
                elif b.method == "max":
                    out[z, y, x] = s.max()
                elif b.method == "mean":
                    out[z, y, x] = np.float32(
                        np.add.reduce(s, dtype=np.float32) / np.float32(len(s))
                    )
                else:
                    q = np.sort(s)
                    m = len(q) // 2
                    out[z, y, x] = (
                        q[m] if len(q) % 2 else q[m - 1] + (q[m] - q[m - 1]) / np.float32(2)
                    )
    return out


@pytest.mark.parametrize("method", ["none", "max", "mean", "median"])
def test_blending_gives_each_voxel_its_definition(method, tmp_path, monkeypatch):
    from petakit5d.core import dtypes

    b = Blend(SCENES["grid2x2"], method, "uint16", "uint16")
    v = _stitch(tmp_path, monkeypatch, b)
    tiles = [tile_data(b.scene, k, "uint16") for k in range(b.layout.n)]
    exp = np.asarray(dtypes.cast(_per_voxel(b, tiles), "uint16"))
    assert np.array_equal(v, exp)


def test_feather_blending_of_uint8_tiles_in_a_batch_that_is_not_full(tmp_path, monkeypatch):
    # uint8 tiles (MATLAB's feather code) blend as the same tiles in uint16 (the kernel)
    from petakit5d.io import read_zarr

    sc = SCENES["partial"]
    out = {}
    for cls in ("uint8", "uint16"):
        d = tmp_path / cls
        d.mkdir()
        os.makedirs(d / "root" / "tiles")
        os.makedirs(d / "root" / "imdist")
        for k in range(2):
            # (tile 2 offset by 100, within uint8)
            t = np.minimum(tile_data(sc, k, "uint8"), 150).astype(np.uint16) + 100 * k
            write_zarr(d / "root" / "tiles" / f"tile_{k + 1:02d}.zarr", t.astype(cls), (4, 8, 8),
                       "matlab")  # fmt: skip
        for k, m in enumerate(scene_maps(sc)[0]):
            write_zarr(d / "root" / "imdist" / f"d{k + 1}.zarr", m, (4, 8, 8), "matlab")
        monkeypatch.chdir(d)
        run_blend(Blend(sc, "feather", cls, cls))
        out[cls] = np.asarray(read_zarr("vol.zarr"))
    assert np.array_equal(out["uint8"].astype(np.uint16), out["uint16"])
    assert (out["uint16"] > 150).any()


def test_mip_of_a_mean_blend_is_the_mip_of_the_result(tmp_path, monkeypatch):
    # tiles that differ by 1: their mean ends in .5, which the uint16 result rounds
    from petakit5d.io import read_zarr

    sc = Scene(SCENES["grid2x2"].layout, holes=False, maps="shared")
    _constant_root(tmp_path, sc, (100, 101, 102, 103))
    monkeypatch.chdir(tmp_path)
    mip = run_blend(Blend(sc, "mean", "uint16", "uint16"))[0]
    v = np.asarray(read_zarr("vol.zarr"))
    assert np.array_equal(mip, v.max(axis=0))
    assert (v == 101).any() or (v == 102).any()


def test_separate_blend_writes_each_tile_on_the_canvas(tmp_path, monkeypatch):
    from petakit5d.io import read_zarr

    sc = SCENES["grid2x2"]
    L = sc.layout
    make_root(sc, "uint16")(tmp_path / "root")
    # (absolute paths, as the frame gives: XR_crop_zarr's mkdir_recursive takes a relative path
    # from the root)
    out = f"{tmp_path}/out"
    os.makedirs(out)
    SB.process_separate_tiles(L.tile_fns(f"{tmp_path}/root/tiles/"), out, L.shift, L.nv_size,
                              batch_size=(8, 16, 16), block_size=(4, 8, 8), parse_cluster=False,
                              uuid="u0")  # fmt: skip
    for k in range(L.n):
        v = np.asarray(read_zarr(f"{out}/tile_{k + 1:02d}.zarr"))
        assert v.shape == tuple(L.nv_size)
        x0, y0, z0 = (int(q) for q in L.shift[k])
        nz, ny, nx = (int(q) for q in L.im_sizes[k])
        exp = np.zeros(L.nv_size, np.uint16)
        exp[z0 : z0 + nz, y0 : y0 + ny, x0 : x0 + nx] = tile_data(sc, k, "uint16")
        assert np.array_equal(v, exp)
        assert os.path.exists(f"{out}/MIPs/tile_{k + 1:02d}_MIP_z.tif")


# ---------------------------------------------------------------------------------------------
# the task: publish order, skips, single batches, errors
# ---------------------------------------------------------------------------------------------


def test_flag_comes_after_the_batches_with_the_mips_in_it(tmp_path, monkeypatch):
    from petakit5d import io as Z
    from petakit5d.testing.publish_spy import PublishSpy

    b = Blend(SCENES["grid2x2"], "feather", "uint16", "uint16")
    make_root(b.scene, b.cls)(tmp_path / "root")
    monkeypatch.chdir(tmp_path)
    events = []
    real = Z.write_zarr

    def spy_write(data, path, **kw):
        events.append(("write", os.fspath(path)))
        return real(data, path, **kw)

    monkeypatch.setattr(Z, "write_zarr", spy_write)
    with PublishSpy() as spy:
        real_replace = os.replace

        def replace(src, dst, *a, **kw):
            if os.fspath(dst).startswith("flags/"):
                from petakit5d.io import read_mat

                events.append(("flag", os.fspath(dst), sorted(read_mat(src))))
            return real_replace(src, dst, *a, **kw)

        monkeypatch.setattr(os, "replace", replace)
        run_blend(b)
    flags = [k for k, e in enumerate(events) if e[0] == "flag"]
    writes = [k for k, e in enumerate(events) if e[0] == "write"]
    nb = _num_batches(b.layout)
    assert len(writes) == nb and len(flags) == math.ceil(nb / b.layout.task_size)
    # each flag after its task's batches; with its MIPs when it appears
    ts = b.layout.task_size
    for t, k in enumerate(flags):
        assert writes[min((t + 1) * ts, nb) - 1] < k
        assert sorted(events[k][2]) == ["batchInds", "mipBboxes", "mips"]
    assert spy.calls  # (the JSON files and the block information file)


def test_done_task_is_skipped_and_redone_with_overwrite(tmp_path, monkeypatch, capsys):
    from petakit5d.io import read_mat

    b = Blend(SCENES["single"], "none", "uint16", "uint16")
    _stitch(tmp_path, monkeypatch, b)
    flag = "flags/blocks_1_1.mat"
    t0 = os.path.getmtime(flag)
    from petakit5d.stitch.blockinfo import per_block_info_paths

    pf = per_block_info_paths("stitchInfo/vol_task_size_10/", 10, 1, 20000)[1]
    bf = "stitchInfo/vol_block_info_task_size_10.mat"
    common = {
        "im_size": b.layout.nv_size,
        "batch_size": b.layout.batch_size,
        "blend_method": "none",
        "save_mip": True,
    }
    assert SB.process_stitch_block([0], bf, pf[0], flag, "vol.zarr", **common) is None
    assert "already exist, skip them!" in capsys.readouterr().out
    assert SB.process_stitch_block([0], bf, pf[0], flag, "vol.zarr", overwrite=True, **common) == [
        True
    ]
    assert os.path.getmtime(flag) >= t0 and "mips" in read_mat(flag)
    # a flag folder that's gone: all tasks were done
    assert (
        SB.process_stitch_block([0], bf, pf[0], "gone/blocks_1_1.mat", "vol.zarr", **common) is None
    )


def test_a_batch_without_tiles_is_zeros_and_its_mip_too(tmp_path, monkeypatch):
    from petakit5d.io import read_mat

    b = Blend(SCENES["single_batches"], "feather", "uint16", "uint16")
    v = _stitch(tmp_path, monkeypatch, b)
    assert (v[:, :, 0] == 0).all()  # (x 0 is before the tile)
    t = tile_data(b.scene, 0, "uint16")
    assert np.array_equal(v[:, :, 1:], t)
    a = read_mat("flags/blocks_1_3.mat")
    mips = list(np.asarray(a["mips"], dtype=object).ravel(order="F"))
    assert mips and all(np.asarray(m).dtype == np.uint16 for m in mips)


def test_a_region_outside_its_tile_is_an_error(tmp_path):
    from petakit5d.io import create_zarr
    from petakit5d.stitch.blockinfo import RegionBlock, TileBlock

    write_zarr(tmp_path / "t.zarr", np.ones((4, 6, 6), np.uint16), (4, 8, 8), "matlab")
    os.makedirs(tmp_path / "flags")
    create_zarr(str(tmp_path / "v.zarr"), data_size=(4, 6, 6), block_size=(4, 8, 8),
                dtype="uint16")  # fmt: skip
    bad = RegionBlock((0, 0, 0, 4, 6, 6), (0, 0, 0, 4, 6, 6), (0, 0, 4, 4, 6, 8))
    tile = TileBlock(
        0, (0, 0, 0, 4, 6, 6), (0, 0, 0, 4, 6, 6), (0, 0, 0, 4, 6, 6), (bad,), (0, 0, 0)
    )
    msg = r"The overlap region \[1 5 1 6 8 4\] is outside the region \[6 6 4\] of tile 1\."
    flag, vol = str(tmp_path / "flags/f.mat"), str(tmp_path / "v.zarr")
    with pytest.raises(SB.StitchBlockError, match=msg):
        SB.process_stitch_block(
            [0], "", "", flag, vol, {0: [tile]}, [str(tmp_path / "t.zarr")], im_size=(4, 6, 6),
            batch_size=(4, 6, 6), blend_method="max",
        )  # fmt: skip


def test_an_empty_region_is_skipped_also_outside_its_tile(tmp_path):
    # MATLAB: any(m_bbox(4:6) >= m_bbox(1:3) & (outside)) checks only the axes that aren't
    # empty, and all(...) adds only regions that aren't: an empty one is no error and no region
    from petakit5d.io import create_zarr, read_zarr
    from petakit5d.stitch.blockinfo import RegionBlock, TileBlock

    write_zarr(tmp_path / "t.zarr", np.full((4, 6, 6), 7, np.uint16), (4, 8, 8), "matlab")
    os.makedirs(tmp_path / "flags")
    vol = str(tmp_path / "v.zarr")
    create_zarr(vol, data_size=(4, 6, 6), block_size=(4, 8, 8), dtype="uint16")
    empty = RegionBlock((0, 0, 0, 0, 6, 6), (0, 0, 0, 0, 6, 6), (-1, 0, 0, -1, 6, 6))
    box = (0, 0, 0, 4, 6, 6)
    tile = TileBlock(0, box, box, box, (empty,), (0, 0, 0))
    done = SB.process_stitch_block(
        [0], "", "", str(tmp_path / "flags/f.mat"), vol, {0: [tile]}, [str(tmp_path / "t.zarr")],
        im_size=(4, 6, 6), batch_size=(4, 6, 6), blend_method="max",
    )  # fmt: skip
    assert done == [True]
    assert (np.asarray(read_zarr(vol)) == 7).all()


@pytest.mark.parametrize(
    ("pool", "psz"),
    [((2, 1, 1), (2, 1, 1)), ((3, 2, 2, 5, 4, 4), (3, 4, 4)),
     ((3, 2, 2, 5, 4, 4, 1, 1, 1), (3, 4, 4)), ((1, 2), None)],
)  # fmt: skip
def test_pooled_maps_use_matlabs_pool_axes(pool, psz):
    # MATLAB: psz = [1, 1, poolSize(3)] or poolSize([4, 5, 3]) ([y x z])
    assert SB._pool_psz(pool) == psz


def test_the_dominance_factor_is_matlabs():
    # num2hex(1000 .^ (1/10)) in MATLAB R2024b (probe 1)
    assert SB._DMT.hex() == "0x1.fec982d5bb8b0p+0"


def test_results_do_not_depend_on_the_thread_count(tmp_path, monkeypatch):
    out = []
    for nt in (1, 3):
        d = tmp_path / f"t{nt}"
        d.mkdir()
        b = Blend(LF_SCENES["lf_major_holes"], "feather", "single", "single")
        out.append(_stitch(d, monkeypatch, b, num_threads=nt))
    assert np.array_equal(out[0].view(np.uint32), out[1].view(np.uint32))


def test_in_process_and_mip_memory(monkeypatch, tmp_path):
    from petakit5d.core import memory

    monkeypatch.setattr(memory, "get_system_memory", lambda: 16.0)
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    assert SB.mip_memory(True) == 16.0
    monkeypatch.setenv("SLURM_JOB_ID", "1")
    assert SB.mip_memory(False) == 16.0
    assert SB.mip_memory(True) == 24 * 20.9
    cfg = tmp_path / "c.json"
    cfg.write_text('{"maxCPUNum": 4, "memPerCPU": 2}')
    assert SB.mip_memory(True, str(cfg)) == 16.0
    # 30 batches or more go to jobs; below, without the cluster, in this process
    assert not SB.in_process(30, (512, 512, 512), "uint16", "none", [[10, 10, 10]], False)
    assert SB.in_process(29, (512, 512, 512), "uint16", "none", [[10, 10, 10]], False)
    # in a Slurm job: the batches (20 B per voxel) or the distance maps must fit its memory
    monkeypatch.setattr(memory, "get_system_memory", lambda: 50.0)
    assert SB.in_process(4, (512, 512, 512), "uint16", "none", [[100, 100, 100]], True)
    assert not SB.in_process(4, (2048, 2048, 2048), "uint16", "none", [[10, 10, 10]], True)
    assert not SB.in_process(4, (64, 64, 64), "uint16", "feather", [[600, 1800, 1800]], True)
    assert SB.in_process(4, (64, 64, 64), "uint16", "none", [[600, 1800, 1800]], True)
