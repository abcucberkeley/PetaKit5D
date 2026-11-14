"""Tests of petakit5d.stitch.distance: stitching's distance maps and feather weights (task S06).

The oracle cases run MATLAB's own functions on generated tiles in MATLAB's zarr layout (both
sides read the same files) and compare every map voxel by voxel: distance_weight_single_axis
on its own, the z window (tukeywin and fastPower), compute_tile_bwdist (the single map and the
per-tile map from the block information's overlap regions; tiles with and without empty
regions, 2D and thin tiles, weight degrees 0 to 20, boxes, data orders, the chunk-file
separator), compute_tile_bwdist_mip_slabs (pool sizes of 3 and 6 values, boxes, the separator of
the matrix division) and its map resized for a tile by feather_distance_map_resize_3d, and
compute_tile_distance_transform in local mode (locations and sizes, per-tile maps, boxes, an
earlier folder of other overlap regions, large files with their MIP slabs).

MATLAB tests whose intent is ported here (their class, StitchingTest, is S09's): the distance
maps of featherBlendingWithThinTiles (test_z_window_has_no_zero_weight_for_thin_tiles),
tilesOfDifferentSizesGetTheirOwnDistanceMaps (test_tiles_of_different_sizes_get_their_own_maps)
and largeFileWithoutRegistration (test_large_file_maps_make_their_mip_slabs).
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.axes import from_yxz
from petakit5d.core.bbox import box_from_matlab
from petakit5d.io.matfile import write_mat
from petakit5d.stitch import distance as sd
from petakit5d.stitch.geometry import data_order_axes_from_matlab
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, InPath, Tree, Value
from petakit5d.testing.compare import tree

# --- generated tiles (Python's (z, y, x) form; written in MATLAB's zarr layout) ---------------


def make_tile(shape, seed, kind="full", cls="uint16"):
    """a tile: random values 1..999 ('full'); with empty wedges that move with z ('dsr', as
    deskewed tiles have); with zero blocks and pixels in some slices ('holes'); with one zero
    slice ('zero_slice'); with an empty corner in every slice ('corner'); all zero ('zero');
    single with NaN, -0 and zeros ('special')"""
    nz, ny, nx = shape
    r = synth.CounterRNG(seed, "s06tile")
    if kind == "zero":
        return np.zeros(shape, dtype=cls)
    if cls == "uint16":
        v = r.integers(1, 999, shape, dtype=np.uint16)
    else:
        v = (r.uniform(shape) * 1000 + 1).astype(np.dtype(cls))
    if kind == "dsr":
        for z in range(nz):
            a = (nz - 1 - z) * (nx // 3) // max(nz - 1, 1)
            b = z * (nx // 4) // max(nz - 1, 1)
            v[z, :, :a] = 0
            if b:
                v[z, : ny // 2, nx - b :] = 0
    elif kind == "holes":
        for z in range(0, nz, 3):
            y0, x0 = (7 * z) % max(ny - 4, 1), (5 * z) % max(nx - 4, 1)
            v[z, y0 : y0 + 3, x0 : x0 + 4] = 0
        v[nz // 2 :, ny // 3, nx // 2] = 0
    elif kind == "zero_slice":
        v[nz // 3] = 0
    elif kind == "corner":
        for y in range(ny):
            v[:, y, : max(nx // 2 - y, 0)] = 0
    elif kind == "special":
        v[:, 2, 3] = np.nan
        v[nz // 2, 4:6, 5:9] = -0.0
        v[:, ny - 3, nx // 2] = 0
    return np.ascontiguousarray(v)


MAT_HEADER = b"MATLAB 5.0 MAT-file, Platform: GLNXA64, Created on: Thu Jan  1 00:00:00 2026"


def write_fixed_mat(path, variables):
    """a MAT v6 file whose header has no time or platform (inputs the same everywhere)"""
    write_mat(path, variables, version="6")
    data = bytearray(Path(path).read_bytes())
    data[:116] = MAT_HEADER.ljust(116, b" ")
    Path(path).write_bytes(bytes(data))


def obj_col(items):
    a = np.empty((len(items), 1), dtype=object)
    for k, v in enumerate(items):
        a[k, 0] = v
    return a


def region(x, y, z):
    """MATLAB's 3 x 2 region [s t] (rows x, y, z)"""
    return np.array([x, y, z], dtype=np.float64)


def block_info(names, pairs=(), n=None):
    """block information: tileFns (relative to the work folder), overlap_matrix, ol_region_cell
    and overlap_map_mat of pairs ((i, j), region of i, region of j), i < j, 1-based"""
    n = len(names) if n is None else n
    om = np.zeros((n, n), dtype=bool)
    cells, omap = [], []
    for k, ((i, j), ri, rj) in enumerate(pairs):
        om[i - 1, j - 1] = True
        cells.append([ri, rj])
        omap.append([k + 1, (j - 1) * n + i])
    return {
        "tileFns": obj_col([f"d/tiles/{t}.zarr" for t in names]),
        "overlap_matrix": om,
        "ol_region_cell": obj_col(cells),
        "overlap_map_mat": np.array(omap, dtype=np.float64).reshape(-1, 2),
    }


def grid_pairs(sizes, ox, oy, zr=None, extra=()):
    """the pairs of a 2 x 2 grid of tiles (1 and 2 along x, 3 and 4 below them; ``sizes``: their
    (z, y, x) shapes, or one for all) with ox and oy voxels of overlap; zr: their z range (all
    slices)"""
    sizes = [sizes] * 4 if len(sizes) == 3 and np.isscalar(sizes[0]) else list(sizes)

    def parts(k):
        nz, ny, nx = sizes[k - 1]
        return {
            "xl": [nx - ox + 1, nx], "xr": [1, ox], "yt": [ny - oy + 1, ny], "yb": [1, oy],
            "ax": [1, nx], "ay": [1, ny], "z": zr or [1, nz],
        }  # fmt: skip

    def reg(k, x, y):
        q = parts(k)
        return region(q[x], q[y], q["z"])

    return [
        ((1, 2), reg(1, "xl", "ay"), reg(2, "xr", "ay")),
        ((1, 3), reg(1, "ax", "yt"), reg(3, "ax", "yb")),
        ((1, 4), reg(1, "xl", "yt"), reg(4, "xr", "yb")),
        ((2, 3), reg(2, "xr", "yt"), reg(3, "xl", "yb")),
        ((2, 4), reg(2, "ax", "yt"), reg(4, "ax", "yb")),
        ((3, 4), reg(3, "xl", "ay"), reg(4, "xr", "ay")),
        *extra,
    ]


def tiles_tree(tiles, info, chunks=(4, 16, 16), slabs=None, extra=None):
    """the case's folder: tiles/<name>.zarr (MATLAB's layout), bi.mat, the MIP slabs
    {pool folder: {name: slab}} and extra(folder)"""

    def fill(folder):
        folder = Path(folder)
        (folder / "tiles").mkdir(parents=True)
        for name, make in tiles.items():
            a = make()
            synth.write_zarr(folder / "tiles" / f"{name}.zarr", a, block_size=chunks,
                             layout="matlab")  # fmt: skip
        write_fixed_mat(folder / "bi.mat", info())
        for pool_dir, slab_makers in (slabs or {}).items():
            (folder / "tiles" / pool_dir).mkdir()
            for name, make in slab_makers.items():
                synth.write_zarr(folder / "tiles" / pool_dir / f"{name}_z.zarr", make(),
                                 block_size=chunks, layout="matlab")  # fmt: skip
        if extra is not None:
            extra(folder)

    return fill


# --- oracle cases: distance_weight_single_axis -------------------------------------------------


def py_dwsa(n, ends, b, df, win):
    s, t = (float(v) for v in np.asarray(ends).ravel())
    rng = (None if s == -math.inf else int(s) - 1, None if t == math.inf else int(t))
    return sd.distance_weight_single_axis(n, rng, b, df, win)


def _dwsa_cases():
    specs = [
        (50, [10, 30], 10, 0.99),
        (50, [1, 50], 10, 0.99),  # covers the axis: double ones
        (50, [-5, 80], 10, 0.99),
        (50, [1, 30], 10, 0.99),
        (50, [20, 50], 10, 0.99),
        (50, [3, 47], 10, 0.99),
        (300, [101, 160], 100, 0.99),
        (300, [5, 290], 100, 0.99),
        (300, [150, 151], 100, 0.99),
        (300, [150, 150], 100, 0.99),
        (300, [200, 120], 100, 0.99),
        (1000, [300, 600], 100, 0.99),
        (40, [10, 30], 1, 0.99),
        (40, [10, 30], 3, 0.9801),
        (40, [10, 30], 5, 0.0),
        (40, [10, 30], 0, 0.99),
        (40, [12, math.inf], 6, 0.99),
        (40, [-math.inf, 20], 6, 0.99),
        (40, [45, 60], 6, 0.99),  # starts beyond the axis: MATLAB grows the vector
        (1, [1, 1], 10, 0.99),
        (2, [2, 2], 10, 0.99),
        (7, [3, 5], 50, 0.99),
    ]
    out = []
    for n, ends, b, df in specs:
        cid = f"dwsa.n{n}.{ends[0]}_{ends[1]}.b{b}.d{df}".replace("inf", "Inf")
        out.append(
            Case(
                cid,
                matlab="distance_weight_single_axis",
                args=[n, list(ends), b, df, "hann"],
                python=py_dwsa,
            )
        )
    return out


# --- oracle cases: the z window ----------------------------------------------------------------

_WINDOW_CODE = """
nz = in{1};
nw = max(round(nz * 1.1), nz + 2);
w = tukeywin(nw, 0.5);
ws = max(round(nz * 0.05), 2);
w = w(ws : ws + nz - 1);
out = {w, fastPower(w, in{2})};
"""


def py_window(nz, wd):
    return [sd.z_window(nz).reshape(-1, 1), sd.z_window(nz, wd).reshape(-1, 1)]


def _window_cases():
    nzs = [2, 3, 5, 9, 10, 19, 20, 29, 30, 49, 50, 99, 100, 101, 257, 500, 513, 1001]
    out = []
    for nz in nzs:
        for wd in (10, 1, -1, 3):
            out.append(
                Case(
                    f"window.nz{nz}.wd{wd}",
                    code=_WINDOW_CODE,
                    args=[nz, wd],
                    nargout=2,
                    python=py_window,
                )
            )
    return out


# --- oracle cases: compute_tile_bwdist ---------------------------------------------------------

_BWDIST_CODE = """
p = in{1}; d = in{2};
compute_tile_bwdist([d 'bi.mat'], p.i, [d 'map.zarr'], p.wd, p.single, p.bs, p.shard, ...
    p.cname, p.bbox, p.dom, false);
out = {};
"""


def _bbox_py(b):
    b = np.asarray(b, dtype=np.float64).ravel()
    return None if b.size == 0 else box_from_matlab(b)


def py_bwdist(p, d):
    sd.compute_tile_bwdist(
        f"{d}bi.mat",
        int(p["i"]) - 1,
        f"{d}map.zarr",
        p["wd"],
        bool(p["single"]),
        from_yxz(p["bs"]),
        None if np.size(p["shard"]) == 0 else from_yxz(p["shard"]),
        p["cname"],
        _bbox_py(p["bbox"]),
        data_order_axes_from_matlab(p["dom"]),
        False,
    )
    return []


def bwdist_case(cid, tiles, info, *, i=1, wd=10, single=True, bs=(16, 16, 4), bbox=None,
                dom=(2, 1, 3), shard=None, cname="zstd", chunks=(4, 16, 16), **kw):  # fmt: skip
    p = {
        "i": i, "wd": wd, "single": single, "bs": list(bs), "bbox": [] if bbox is None else bbox,
        "dom": list(dom), "shard": [] if shard is None else list(shard), "cname": cname,
    }  # fmt: skip
    return Case(
        cid,
        code=_BWDIST_CODE,
        args=[In("p"), InPath("d")],
        nargout=0,
        python=py_bwdist,
        inputs={"p": Value(lambda: p), "d": Tree(tiles_tree(tiles, info, chunks), "d")},
        files=True,
        **kw,
    )


def one_tile(shape, seed, kind="full", cls="uint16"):
    tiles = {"t1": lambda: make_tile(shape, seed, kind, cls)}
    return tiles, lambda: block_info(["t1"])


def grid_tiles(shape, seed, kind="full", ox=6, oy=5, zr=None, extra=(), sizes=None,
               regions_of=None):  # fmt: skip
    """four tiles of ``shape`` (or ``sizes``) and their block information; ``regions_of``: the
    shape the regions are made for (the turned tile of another data order)"""
    names = ["t1", "t2", "t3", "t4"]
    sizes = sizes or [shape] * 4
    tiles = {n: (lambda k=k: make_tile(sizes[k], seed + k, kind)) for k, n in enumerate(names)}
    rs = [regions_of] * 4 if regions_of else sizes
    return tiles, lambda: block_info(names, grid_pairs(rs, ox, oy, zr, extra))


def turned_shape(shape, dom):
    """the (z, y, x) shape of a tile turned as compute_tile_bwdist turns it (z last)"""
    if dom[2] == 3:
        return tuple(shape)
    q = sd._mpermute_axes(dom)
    return tuple(shape[q[k]] for k in range(3))


def _bwdist_cases():
    out = []
    s3, s2 = (13, 37, 29), (1, 41, 33)
    # the single map: kinds of tiles, 2D and 3D, weight degrees
    for kind in ("full", "dsr", "holes", "zero_slice", "zero"):
        for wd in (10, 1, 2):
            out.append(bwdist_case(f"single.{kind}.wd{wd}", *one_tile(s3, 11, kind), wd=wd))
    for kind in ("full", "corner", "holes", "zero"):
        out.append(bwdist_case(f"single2d.{kind}", *one_tile(s2, 12, kind)))
    for wd in (0, 3, -1, -2, 7, 20, 31, 10000):
        out.append(bwdist_case(f"single.dsr.wd{wd}", *one_tile(s3, 13, "dsr"), wd=wd))
    out.append(bwdist_case("single.special.single", *one_tile(s3, 14, "special", "float32")))
    out.append(bwdist_case("single.dsr.uint8", *one_tile(s3, 15, "dsr", "uint8")))
    # thin tiles (the z window)
    for nz in (2, 3, 5, 9, 10, 21, 30):
        out.append(bwdist_case(f"single.thin{nz}", *one_tile((nz, 23, 19), 20 + nz, "dsr")))
    # bigger, so that threads share the slices
    out.append(bwdist_case("single.dsr.big", *one_tile((40, 150, 121), 16, "dsr"),
                           bs=(64, 64, 16), chunks=(16, 64, 64)))  # fmt: skip
    # boxes of full weight ([ymin xmin zmin ymax xmax zmax])
    boxes = {
        "inside": [5, 6, 3, 30, 22, 10],
        "edges": [1, 1, 1, 37, 29, 13],
        "beyond": [-50, -50, -10, 200, 200, 100],
        "y_only": [5, 1, 1, 30, 29, 13],
        "z_only": [1, 1, 4, 37, 29, 9],
        "open": [5, 6, 3, math.inf, math.inf, math.inf],
        "thin": [18, 14, 7, 18, 14, 7],
    }
    for name, b in boxes.items():
        out.append(bwdist_case(f"single.bbox.{name}", *one_tile(s3, 17, "dsr"), bbox=b))
    out.append(bwdist_case("single.bbox.wd1", *one_tile(s3, 17, "dsr"), bbox=boxes["inside"],
                           wd=1))  # fmt: skip
    out.append(bwdist_case("single.bbox.wd0", *one_tile(s3, 17, "dsr"), bbox=boxes["inside"],
                           wd=0))  # fmt: skip
    out.append(bwdist_case("single2d.bbox", *one_tile(s2, 18, "dsr"), bbox=[4, 5, 1, 30, 25, 1]))
    # a box far from the tile's ends: the tapers end and the decay (and the 1e-45 floor) shows
    out.append(bwdist_case("single.bbox.far", *one_tile((12, 260, 230), 18, "dsr"),
                           bbox=[5, 5, 2, 40, 40, 10], bs=(64, 64, 8)))  # fmt: skip
    out.append(bwdist_case("single2d.bbox.far", *one_tile((1, 300, 250), 18, "holes"),
                           bbox=[150, 120, 1, 170, 140, 1], bs=(64, 64, 1)))  # fmt: skip
    # data orders (dataOrderMat): the tile turned so that z is last, and back
    for dom in ((1, 2, 3), (3, 2, 1), (3, 1, 2), (1, 3, 2)):
        tag = "".join(str(v) for v in dom)
        out.append(bwdist_case(f"single.dom{tag}", *one_tile(s3, 19, "dsr"), dom=dom))
        out.append(bwdist_case(f"single.dom{tag}.bbox", *one_tile(s3, 19, "dsr"), dom=dom,
                               bbox=[3, 4, 2, 11, 20, 25]))  # fmt: skip
    # compressor and shard size are not used; chunks and the separator
    out.append(bwdist_case("single.lz4_ignored", *one_tile(s3, 21), cname="lz4"))
    out.append(bwdist_case("single.blocks_large", *one_tile(s3, 21), bs=(500, 500, 500)))
    out.append(bwdist_case("single.blocks_odd", *one_tile(s3, 21), bs=(7, 5, 3)))
    out.append(bwdist_case("single.separator_slash", *one_tile((26, 40, 40), 22), bs=(2, 2, 1),
                           marks=("slow",)))  # fmt: skip
    # the per-tile map from the overlap regions
    g3 = (12, 34, 30)
    for kind in ("full", "dsr", "holes"):
        for i in (1, 2, 4):
            out.append(bwdist_case(f"tile.{kind}.i{i}", *grid_tiles(g3, 30, kind), i=i,
                                   single=False))  # fmt: skip
    out.append(bwdist_case("tile.dsr.wd1", *grid_tiles(g3, 31, "dsr"), i=1, single=False, wd=1))
    out.append(bwdist_case("tile.dsr.wd0", *grid_tiles(g3, 31, "dsr"), i=1, single=False, wd=0))
    out.append(bwdist_case("tile.dsr.wd20", *grid_tiles(g3, 31, "dsr"), i=3, single=False,
                           wd=20))  # fmt: skip
    out.append(bwdist_case("tile.dsr.bbox", *grid_tiles(g3, 32, "dsr"), i=2, single=False,
                           bbox=[4, 3, 2, 30, 25, 11]))  # fmt: skip
    out.append(bwdist_case("tile.zrange", *grid_tiles(g3, 33, "dsr", zr=[3, 9]), i=1,
                           single=False))  # fmt: skip
    out.append(bwdist_case("tile.wide_overlap", *grid_tiles(g3, 34, "holes", ox=20, oy=25), i=4,
                           single=False))  # fmt: skip
    # regions that are empty, inside the tile without empty voxels, at its edges, one voxel
    names = ["t1", "t2"]
    tiles2 = {n: (lambda k=k: make_tile(g3, 36 + k, "holes")) for k, n in enumerate(names)}
    regions = [
        ("empty", region([10, 9], [1, 34], [1, 12]), region([1, 5], [1, 34], [1, 12])),
        ("interior", region([10, 20], [10, 20], [4, 8]), region([8, 12], [9, 15], [5, 9])),
        ("edge", region([1, 30], [1, 34], [1, 12]), region([25, 30], [30, 34], [11, 12])),
        ("voxel", region([15, 15], [17, 17], [6, 6]), region([1, 1], [1, 1], [1, 1])),
        ("zero_z", region([1, 8], [1, 34], [1, 1]), region([23, 30], [1, 34], [12, 12])),
    ]
    for name, ri, rj in regions:
        for i in (1, 2):
            out.append(bwdist_case(
                f"tile.region_{name}.i{i}", tiles2,
                lambda ri=ri, rj=rj: block_info(names, [((1, 2), ri, rj)]), i=i, single=False,
            ))  # fmt: skip
    # no overlaps: ones before the window; a pair without its region (MATLAB fails)
    out.append(bwdist_case("tile.no_pairs", tiles2, lambda: block_info(names), i=1, single=False))

    def missing():
        b = block_info(names, [((1, 2), *regions[1][1:])])
        b["overlap_map_mat"] = np.array([[1.0, 7.0]])
        return b

    out.append(bwdist_case("tile.missing_region", tiles2, missing, i=1, single=False,
                           expect_error=True, py_error=(sd.DistanceError,)))  # fmt: skip
    # per-tile maps of 2D tiles, other data orders
    out.append(bwdist_case("tile2d", *grid_tiles((1, 34, 30), 37, "holes", zr=[1, 1]), i=4,
                           single=False))  # fmt: skip
    for dom in ((1, 2, 3), (3, 2, 1), (3, 1, 2)):
        tag = "".join(str(v) for v in dom)
        tt = grid_tiles(g3, 38, "dsr", regions_of=turned_shape(g3, dom))
        out.append(bwdist_case(f"tile.dom{tag}", *tt, i=1, single=False, dom=dom))
    return out


# --- oracle cases: compute_tile_bwdist_mip_slabs -----------------------------------------------

_SLAB_CODE = """
p = in{1}; d = in{2};
compute_tile_bwdist_mip_slabs([d 'bi.mat'], p.i, [d 'map.zarr'], p.wd, true, p.bs, p.shard, ...
    p.cname, p.pool, p.bbox, p.dom, false);
out = {};
"""

_SLAB_RESIZE_CODE = """
p = in{1}; d = in{2};
compute_tile_bwdist_mip_slabs([d 'bi.mat'], p.i, [d 'map.zarr'], p.wd, true, p.bs, p.shard, ...
    p.cname, p.pool, p.bbox, p.dom, false);
m = readzarr([d 'map.zarr']);
r = p.crop;
out = {feather_distance_map_resize_3d(m(r(1) : r(4), r(2) : r(5), r(3) : r(6)), p.tbox, p.wd)};
"""


def py_slab(p, d):
    box = _bbox_py(p["bbox"])
    sd.compute_tile_bwdist_mip_slabs(
        f"{d}bi.mat",
        int(p["i"]) - 1,
        f"{d}map.zarr",
        p["wd"],
        True,
        from_yxz(p["bs"]),
        None if np.size(p["shard"]) == 0 else from_yxz(p["shard"]),
        p["cname"],
        _pool_py(p["pool"]),
        box,
        data_order_axes_from_matlab(p["dom"]),
        False,
    )
    return []


def py_slab_resize(p, d):
    from petakit5d.io import read_zarr

    py_slab(p, d)
    m = read_zarr(f"{d}map.zarr")
    m = sd._as_volume(m)
    z0, y0, x0, z1, y1, x1 = box_from_matlab(p["crop"])
    crop = np.ascontiguousarray(m[z0:z1, y0:y1, x0:x1])
    return sd.feather_distance_map_resize_3d(crop, box_from_matlab(p["tbox"]), p["wd"])


def _pool_py(pool):
    ps = [float(v) for v in np.asarray(pool).ravel()]
    return [int(v) for g in range(0, len(ps), 3) for v in from_yxz(ps[g : g + 3])]


def _pool_dir(pool):
    return "MIP_Slabs_pooling_" + "_".join(str(int(v)) for v in pool)


def slab_case(
    cid, slab_shape, seed, kind="dsr", *, wd=10, pool=(2, 2, 2), bbox=None, dom=(2, 1, 3),
    bs=(16, 16, 4), shard=None, cname="zstd", resize=None, **kw,
):  # fmt: skip
    p = {
        "i": 1, "wd": wd, "bs": list(bs), "bbox": [] if bbox is None else bbox, "dom": list(dom),
        "shard": [] if shard is None else list(shard), "cname": cname, "pool": list(pool),
    }  # fmt: skip
    if resize is not None:
        p["crop"], p["tbox"] = resize
    tiles = {"t1": lambda: make_tile((3, 5, 4), seed)}
    slabs = {_pool_dir(pool): {"t1": lambda: make_tile(slab_shape, seed, kind)}}
    return Case(
        cid,
        code=_SLAB_RESIZE_CODE if resize is not None else _SLAB_CODE,
        args=[In("p"), InPath("d")],
        nargout=1 if resize is not None else 0,
        python=py_slab_resize if resize is not None else py_slab,
        out_roles=("vol",) if resize is not None else None,
        inputs={"p": Value(lambda: p),
                "d": Tree(tiles_tree(tiles, lambda: block_info(["t1"]), slabs=slabs), "d")},
        files=True,
        **kw,
    )  # fmt: skip


def _slab_cases():
    s = (11, 30, 26)
    out = []
    for kind in ("full", "dsr", "holes", "zero_slice"):
        out.append(slab_case(f"slab.{kind}", s, 40, kind))
    out.append(slab_case("slab2d", (1, 30, 26), 41))
    for wd in (1, 3, 0, -1, 20):
        out.append(slab_case(f"slab.wd{wd}", s, 42, wd=wd, bbox=[5, 7, 3, 50, 40, 18]))
    for name, pool in {"222": (2, 2, 2), "335": (3, 3, 5), "6v": (2, 2, 4, 2, 2, 1),
                       "6v_xy": (4, 4, 4, 3, 2, 1)}.items():  # fmt: skip
        out.append(slab_case(f"slab.pool{name}", s, 43, pool=pool))
        out.append(slab_case(f"slab.pool{name}.bbox", s, 43, pool=pool,
                             bbox=[9, 11, 5, 51, 45, 30]))  # fmt: skip
    for name, b in {
        "beyond": [-100, -100, -10, 500, 500, 100],
        "edge": [1, 1, 1, 60, 52, 22],
        "open": [9, 11, 5, math.inf, math.inf, math.inf],
        "odd": [7, 9, 3, 33, 29, 15],
    }.items():
        out.append(slab_case(f"slab.bbox.{name}", s, 44, bbox=b))
    # far from the slab's ends: the decay of 0.99^pool and the floor 1e-50^(1/wd)
    for name, pool in {"222": (2, 2, 2), "6v": (4, 4, 4, 3, 2, 1)}.items():
        out.append(slab_case(f"slab.bbox.far.pool{name}", (6, 130, 120), 44, pool=pool,
                             bbox=[3, 3, 1, 30, 30, 8], bs=(64, 64, 8)))  # fmt: skip
    out.append(slab_case("slab.lz4", s, 45, cname="lz4"))
    # the separator: sz / blockSize (a matrix division) is small, sz ./ blockSize is not
    out.append(slab_case("slab.separator_mrdivide", (26, 40, 40), 46, kind="full",
                         bs=(2, 2, 1), marks=("slow",)))  # fmt: skip
    for dom in ((1, 2, 3), (3, 2, 1), (2, 3, 1)):
        tag = "".join(str(v) for v in dom)
        out.append(slab_case(f"slab.dom{tag}", (10, 10, 10), 47, dom=dom))
        # z not last: MATLAB makes the temporary array of the turned size and fails to write
        # the map into it (the array is left; its .zarray is compared as text, the layouts
        # differ, so it is not compared)
        out.append(slab_case(f"slab.dom{tag}.wide", s, 47, dom=dom, expect_error=dom[2] != 3,
                             py_error=(sd.DistanceError,),
                             tree=tree(ignore=(r"map\.zarr_[^/]*/",))))  # fmt: skip
    # the map resized for a tile's part, raised to the power (the blending's use of it)
    out.append(slab_case("slab.resize", s, 48, bbox=[9, 11, 5, 51, 45, 30],
                         resize=([2, 3, 1, 20, 18, 9], [1, 1, 1, 38, 34, 17])))  # fmt: skip
    out.append(slab_case("slab.resize.wd3", s, 48, wd=3, pool=(3, 3, 5),
                         resize=([1, 1, 1, 30, 26, 11], [1, 1, 1, 90, 78, 55])))  # fmt: skip
    out.append(slab_case("slab2d.resize", (1, 30, 26), 49,
                         resize=([2, 3, 1, 20, 18, 1], [1, 1, 1, 38, 34, 1])))  # fmt: skip
    return out


# --- oracle cases: compute_tile_distance_transform ---------------------------------------------

# (the MATLAB framework changes MATLAB's current folder, so every path is absolute: the block
# information's tile paths are made absolute in each side's work folder first)
_WRAPPER_CODE = """
p = in{1}; d = in{2};
a = load([d 'bi.mat']);
a.tileFns = cellfun(@(x) [d x(3 : end)], a.tileFns, 'unif', 0);
save([d 'bi.mat'], '-struct', 'a');
paths = cellfun(@(x) sprintf('%stiles/%s.zarr', d, x), p.names, 'unif', 0);
[o, idx] = compute_tile_distance_transform([d 'bi.mat'], [d 'stitch/'], paths, p.args{:}, ...
    'parseCluster', false);
out = {o, idx};
"""


def absolute_block_info(d):
    """the block information with its tile paths in the work folder ``d``"""
    from petakit5d.io.matfile import read_mat

    a = read_mat(f"{d}bi.mat")
    a["tileFns"] = obj_col([d + str(x)[2:] for x in a["tileFns"].ravel(order="F")])
    write_mat(f"{d}bi.mat", a, version="7")


def _wrapper_kwargs(args):
    kw = dict(zip(args[0::2], args[1::2], strict=True))
    conv = {
        "singleDistMap": ("single_dist_map", bool),
        "locIds": ("loc_ids", lambda v: list(np.asarray(v, dtype=np.float64).ravel())),
        "distBboxes": ("dist_bboxes", lambda v: [box_from_matlab(r) for r in
                                                 np.asarray(v, dtype=np.float64).reshape(-1, 6)]),
        "blendWeightDegree": ("blend_weight_degree", float),
        "blockSize": ("block_size", lambda v: from_yxz(v)),
        "largeFile": ("large_file", bool),
        "poolSize": ("pool_size", _pool_py),
        "dataOrderMat": ("data_order_axes", data_order_axes_from_matlab),
        "compressor": ("compressor", str),
    }  # fmt: skip
    return {conv[k][0]: conv[k][1](v) for k, v in kw.items()}


def py_wrapper(p, d):
    names = list(p["names"]) if isinstance(p["names"], (list, tuple)) else [p["names"]]
    absolute_block_info(d)
    paths = [f"{d}tiles/{n}.zarr" for n in names]
    o, idx = sd.compute_tile_distance_transform(
        f"{d}bi.mat", f"{d}stitch/", paths, parse_cluster=False, **_wrapper_kwargs(p["args"])
    )
    return [o, np.asarray(idx, dtype=np.float64).reshape(-1, 1)]


def wrapper_case(cid, tiles, info, names, args, *, slabs=None, extra=None, **kw):
    p = {"names": list(names), "args": list(args)}
    return Case(
        cid,
        code=_WRAPPER_CODE,
        args=[In("p"), InPath("d")],
        nargout=2,
        python=py_wrapper,
        out_roles=(None, "index1"),
        inputs={"p": Value(lambda: p), "d": Tree(tiles_tree(tiles, info, slabs=slabs,
                                                            extra=extra), "d")},
        files=True,
        tree=tree(ignore_vars=("ip", "t", "t1"), ignore=(r"(^|/)bi\.mat$",)),
        **kw,
    )  # fmt: skip


def _stale_maps(folder):
    """an earlier imdist folder of the first tile made with other overlap regions"""
    p = Path(folder) / "stitch" / "imdist" / "t1"
    p.mkdir(parents=True)
    write_fixed_mat(p / "parameters.mat", {
        "overlap_matrix": np.zeros((4, 4), dtype=bool),
        "ol_region_cell": obj_col([]),
        "overlap_map_mat": np.zeros((0, 2)),
    })  # fmt: skip
    synth.write_zarr(p / "t1_wr_10.zarr", np.ones((2, 3, 4), np.float32), block_size=(2, 3, 4),
                     layout="matlab")  # fmt: skip


def _wrapper_cases():
    g = (10, 30, 26)
    names = ["t1", "t2", "t3", "t4"]
    out = []
    tiles, info = grid_tiles(g, 60, "dsr")
    out.append(wrapper_case("wrapper.single.one_loc", tiles, info, names,
                            ["locIds", [1, 1, 1, 1]]))  # fmt: skip
    out.append(wrapper_case("wrapper.single.same_sizes", tiles, info, names,
                            ["locIds", [1, 1, 2, 2]]))  # fmt: skip
    sizes = [g, (10, 30, 22), g, (10, 30, 22)]
    tiles_s, info_s = grid_tiles(g, 61, "dsr", sizes=sizes)
    out.append(wrapper_case("wrapper.single.two_sizes", tiles_s, info_s, names,
                            ["locIds", [1, 2, 1, 2]]))  # fmt: skip
    out.append(wrapper_case("wrapper.single.three_locs", tiles_s, info_s, names,
                            ["locIds", [3, 1, 2, 1]]))  # fmt: skip
    out.append(wrapper_case("wrapper.single.three_locs_same", tiles, info, names,
                            ["locIds", [1, 2, 3, 1]]))  # fmt: skip
    out.append(
        wrapper_case(
            "wrapper.single.bbox_row",
            tiles_s,
            info_s,
            names,
            ["locIds", [1, 2, 1, 2], "distBboxes", [4, 5, 2, 25, 20, 9]],
        )
    )
    out.append(
        wrapper_case(
            "wrapper.single.bbox_rows",
            tiles_s,
            info_s,
            names,
            [
                "locIds",
                [1, 2, 1, 2],
                "distBboxes",
                np.array([[4, 5, 2, 25, 20, 9], [1, 1, 1, 30, 22, 10]], float),
            ],
        )
    )
    out.append(wrapper_case("wrapper.single.bbox_rows_wrong", tiles_s, info_s, names,
                            ["locIds", [1, 2, 1, 2], "distBboxes",
                             np.array([[4, 5, 2, 25, 20, 9]] * 3, float)],
                            expect_error=True, py_error=(sd.DistanceError,)))  # fmt: skip
    out.append(wrapper_case("wrapper.tile", tiles, info, names,
                            ["singleDistMap", False, "locIds", [1, 1, 1, 1]]))  # fmt: skip
    out.append(wrapper_case("wrapper.tile.sizes", tiles_s, info_s, names,
                            ["singleDistMap", False, "blendWeightDegree", 3]))  # fmt: skip
    out.append(wrapper_case("wrapper.tile.stale", tiles, info, names,
                            ["singleDistMap", False], extra=_stale_maps))  # fmt: skip
    out.append(wrapper_case("wrapper.single.stale_kept", tiles, info, names,
                            ["locIds", [1, 1, 1, 1]], extra=_stale_maps))  # fmt: skip
    out.append(wrapper_case("wrapper.single.blocks", tiles, info, names,
                            ["locIds", [1, 1, 1, 1], "blockSize", [8, 8, 4], "compressor",
                             "zstd"]))  # fmt: skip
    out.append(wrapper_case("wrapper.single.dom", tiles, info, names,
                            ["locIds", [1, 1, 1, 1], "dataOrderMat", [1, 2, 3]]))  # fmt: skip
    # large files: the MIP slabs are made first
    out.append(wrapper_case("wrapper.large", tiles, info, names,
                            ["locIds", [1, 1, 1, 1], "largeFile", True, "poolSize", [2, 2, 2],
                             "blockSize", [16, 16, 8]]))  # fmt: skip
    out.append(
        wrapper_case(
            "wrapper.large.tile",
            tiles_s,
            info_s,
            names,
            [
                "singleDistMap",
                False,
                "largeFile",
                True,
                "poolSize",
                [2, 2, 2],
                "distBboxes",
                [3, 3, 2, 20, 20, 8],
            ],
            expect_error=True,
            py_error=(IndexError, sd.DistanceError),
        )
    )
    out.append(
        wrapper_case(
            "wrapper.large.sizes",
            tiles_s,
            info_s,
            names,
            [
                "locIds",
                [1, 2, 1, 2],
                "largeFile",
                True,
                "poolSize",
                [2, 2, 3, 2, 2, 1],
                "distBboxes",
                [3, 3, 2, 20, 20, 8],
            ],
        )
    )
    return out


CASES = _dwsa_cases() + _window_cases() + _bwdist_cases() + _slab_cases() + _wrapper_cases()


# MATLAB on Windows (ispc) writes the maps into imdist/ itself, not imdist/<first tile>/, and
# so does the port there: the wrapper's references, made on Linux, hold the other folder. The
# Windows layout is tested on every system by forcing the Windows branch (test_*[windows]).
_WINDOWS_WRAPPER = (
    "compute_tile_distance_transform on Windows writes into imdist/ (MATLAB's ispc branch); "
    "these references are MATLAB's on Linux (imdist/<first tile>/): the [windows] tests cover it"
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    if sd._is_pc() and case.id.startswith("wrapper."):
        pytest.skip(_WINDOWS_WRAPPER)
    oracle.check(case, tmp_path)


# --- unit tests --------------------------------------------------------------------------------


def _ref_map(tile, wd, *, single=True, overlaps=(), box=None):
    """compute_tile_bwdist's map written out as MATLAB's code reads, on whole arrays (a second
    implementation for the tests): standard data order"""
    from petakit5d.core import dtypes
    from petakit5d.matlabcompat.distance import bwdist
    from petakit5d.matlabcompat.elementary import fast_power

    v = sd._as_volume(tile)
    nz, ny, nx = v.shape
    orig = v != 0
    e = v.copy()
    e[:, [0, -1], :] = 0
    e[:, :, [0, -1]] = 0
    if nz > 1:
        e[[0, -1]] = 0
    e = e == 0

    def dist(b):
        return np.asarray(fast_power(dtypes.rdivide(dtypes.plus(bwdist(b), 1.0), 10.0), wd),
                          dtype=np.float32)  # fmt: skip

    if single:
        c = int(dtypes.mround((nz + 1) / 2)) - 1
        dc = dist(e[c])
        m = np.repeat(dc[None], nz, axis=0)
        for z in range(nz):
            if not np.array_equal(e[z], e[c]):
                m[z] = dist(e[z])
    else:
        m = np.ones(v.shape, np.float32)
        for b in overlaps:
            if any(b[k + 3] <= b[k] for k in range(3)):
                continue
            lo = [max(b[k] - 1, 0) for k in range(3)]
            hi = [min(b[k + 3] + 1, v.shape[k]) for k in range(3)]
            crop = e[lo[0] : hi[0], lo[1] : hi[1], lo[2] : hi[2]]
            if not crop.any():
                continue
            dd = np.stack([dist(crop[z]) for z in range(crop.shape[0])])
            m[b[0] : b[3], b[1] : b[4], b[2] : b[5]] = dd[
                b[0] - lo[0] : b[3] - lo[0],
                b[1] - lo[1] : b[4] - lo[1],
                b[2] - lo[2] : b[5] - lo[2],
            ]
    win = sd.z_window(nz, wd)
    m = (m.astype(np.float64) * win[:, None, None]).astype(np.float32)
    m = m * orig.astype(np.float32)
    if box is not None:
        bm = list(sd._box_axis_weights((ny, nx, nz), [float(x) for x in box],
                                       (100, 100, 100), (0.99,) * 3))  # fmt: skip
        wt = (bm[0][:, None] * bm[1][None, :])[None] * bm[2][:, None, None]
        wt = np.fmax(np.asarray(fast_power(wt.astype(np.float32), wd), np.float32), sd._TINY)
        m = m * wt
    return m


def _bits(a):
    return np.ascontiguousarray(a).view(np.uint32 if a.dtype == np.float32 else np.uint64)


@pytest.mark.parametrize("kind", ["full", "dsr", "holes", "zero_slice", "zero"])
@pytest.mark.parametrize("wd", [10, 1, 3])
def test_single_map_equals_the_whole_array_code(kind, wd):
    t = make_tile((13, 37, 29), 3, kind)
    got = sd.tile_distance_map(t, wd, num_threads=3)
    np.testing.assert_array_equal(_bits(got), _bits(_ref_map(t, wd)))


@pytest.mark.parametrize("box", [[5, 6, 3, 30, 22, 10], [-9, -9, -9, 99, 99, 99]])
def test_box_weights_equal_the_whole_array_code(box):
    t = make_tile((13, 37, 29), 4, "dsr")
    got = sd.tile_distance_map(t, 10, dist_bbox=box_from_matlab(box))
    np.testing.assert_array_equal(_bits(got), _bits(_ref_map(t, 10, box=box)))


def test_per_tile_map_equals_the_whole_array_code():
    t = make_tile((12, 34, 30), 5, "holes")
    ov = [
        (0, 0, 24, 12, 34, 30),
        (0, 29, 0, 12, 34, 30),
        (2, 10, 10, 6, 20, 20),
        (3, 3, 3, 2, 9, 9),
    ]
    got = sd.tile_distance_map(t, 10, single_dist_map=False, overlaps=ov, num_threads=2)
    np.testing.assert_array_equal(_bits(got), _bits(_ref_map(t, 10, single=False, overlaps=ov)))


@pytest.mark.parametrize("single", [True, False])
def test_maps_do_not_depend_on_the_thread_count(single):
    t = make_tile((21, 40, 33), 6, "dsr")
    ov = [(0, 0, 25, 21, 40, 33), (0, 33, 0, 21, 40, 33)]
    maps = [
        sd.tile_distance_map(
            t,
            10,
            single_dist_map=single,
            overlaps=ov,
            num_threads=n,
            dist_bbox=(2, 4, 5, 18, 30, 30),
        )  # fmt: skip
        for n in (1, 2, 5, 8)
    ]
    for m in maps[1:]:
        np.testing.assert_array_equal(_bits(m), _bits(maps[0]))


def test_out_is_filled_and_returned():
    t = make_tile((9, 20, 17), 7, "dsr")
    buf = np.full((9, 20, 17), np.float32(7))
    r = sd.tile_distance_map(t, 10, out=buf)
    assert r is buf
    np.testing.assert_array_equal(_bits(buf), _bits(sd.tile_distance_map(t, 10)))
    with pytest.raises(ValueError):
        sd.tile_distance_map(t, 10, out=np.zeros((9, 20, 16), np.float32))


def test_zero_voxels_get_zero_weight_and_data_positive():
    t = make_tile((9, 20, 17), 8, "holes")
    m = sd.tile_distance_map(t, 10)
    assert np.all(m[t == 0] == 0)
    assert np.all(m[t != 0] > 0)


@pytest.mark.parametrize("nz", list(range(1, 41)))
def test_z_window_has_no_zero_weight_for_thin_tiles(nz):
    # StitchingTest/featherBlendingWithThinTiles: the window used to give the first slice weight
    # 0 for 10-29 slices and fail for fewer than 10
    w = sd.z_window(nz, 10)
    assert w.shape == (nz,)
    assert np.all(w > 0) and np.all(w <= 1)


def test_distance_weight_single_axis_shape():
    w = sd.distance_weight_single_axis(60, (9, 40), 5, 0.99)
    assert w.dtype == np.float32 and w.shape == (60,)
    # 1 on the range (its ends included), tapers down to 1e-3 over 5 voxels, then the decay
    assert np.all(w[9:40] == 1)
    assert w[4] == np.float32(1e-3) and w[44] == np.float32(1e-3)
    assert np.all(np.diff(w[:10].astype(np.float64)) > 0)
    assert np.all(np.diff(w[39:].astype(np.float64)) < 0)
    # the box covers the axis: double ones (MATLAB returns before the cast)
    w = sd.distance_weight_single_axis(60, (0, 60))
    assert w.dtype == np.float64 and np.all(w == 1)
    # a box that starts beyond the axis grows the vector, as MATLAB's assignment does
    assert sd.distance_weight_single_axis(40, (44, 60), 6).shape == (45,)


def test_overflowing_powers_give_matlab_nan():
    # Inf distances (fastPower overflows in single) times a window that underflows to 0 in
    # double: MATLAB's NaN (x86's, sign bit set) on every processor
    m = sd.tile_distance_map(make_tile((13, 37, 29), 9), 10000)
    assert np.isinf(m).any()
    nan = np.isnan(m)
    assert nan.any() and np.all(m[nan].view(np.uint32) == 0xFFC00000)


def test_weight_degree_zero_is_the_data_mask():
    t = make_tile((7, 15, 13), 10, "holes")
    for single in (True, False):
        m = sd.tile_distance_map(t, 0, single_dist_map=single, overlaps=[(0, 0, 0, 7, 15, 13)])
        np.testing.assert_array_equal(m, (t != 0).astype(np.float32))


def test_non_integer_weight_degree_fails():
    with pytest.raises(sd.DistanceError):
        sd.tile_distance_map(make_tile((3, 8, 8), 1), 2.5)


def test_data_order_round_trip():
    # a tile in x, y, z order gives the map of its (z, y, x) form turned the same way, when z
    # stays last (no turn) and through MATLAB's permute when it does not
    t = make_tile((9, 21, 17), 11, "dsr")
    m = sd.tile_distance_map(t, 10)
    t_xy = np.ascontiguousarray(np.transpose(t, (0, 2, 1)))
    m_xy = sd.tile_distance_map(t_xy, 10, data_order_axes=data_order_axes_from_matlab([1, 2, 3]))
    np.testing.assert_array_equal(
        _bits(m_xy), _bits(np.ascontiguousarray(np.transpose(m, (0, 2, 1))))
    )
    t_zx = np.ascontiguousarray(np.transpose(t, (1, 0, 2)))  # (y, z, x): MATLAB (z, x, y)
    m_zx = sd.tile_distance_map(t_zx, 10, data_order_axes=data_order_axes_from_matlab([2, 3, 1]))
    np.testing.assert_array_equal(
        _bits(m_zx), _bits(np.ascontiguousarray(np.transpose(m, (1, 0, 2))))
    )


def _write_tile_info(folder: Path, tiles: dict, pairs=(), layout="matlab"):
    (folder / "tiles").mkdir(parents=True, exist_ok=True)
    for name, a in tiles.items():
        synth.write_zarr(folder / "tiles" / f"{name}.zarr", a, block_size=(4, 16, 16),
                         layout=layout)  # fmt: skip
    info = block_info(list(tiles), pairs)
    info["tileFns"] = obj_col([str(folder / "tiles" / f"{n}.zarr") for n in tiles])
    write_fixed_mat(folder / "bi.mat", info)
    return str(folder / "bi.mat")


def _read_map(path):
    from petakit5d.io import read_zarr

    return sd._as_volume(read_zarr(path))


def test_native_and_matlab_tiles_give_the_same_map(tmp_path):
    t = make_tile((9, 30, 26), 12, "dsr")
    maps = []
    for layout in ("matlab", "native"):
        f = tmp_path / layout
        bi = _write_tile_info(f, {"t1": t}, layout=layout)
        sd.compute_tile_bwdist(bi, 0, str(f / "m.zarr"), 10, True, (4, 16, 16), None, "zstd",
                               None, (2, 1, 0))  # fmt: skip
        maps.append(_read_map(f / "m.zarr"))
    np.testing.assert_array_equal(_bits(maps[0]), _bits(maps[1]))
    np.testing.assert_array_equal(_bits(maps[0]), _bits(sd.tile_distance_map(t, 10)))


def test_map_file_ignores_compressor_and_shards(tmp_path):
    bi = _write_tile_info(tmp_path, {"t1": make_tile((5, 20, 18), 13)})
    out = tmp_path / "m.zarr"
    sd.compute_tile_bwdist(bi, 0, str(out), 10, True, (4, 8, 8), (2, 4, 4), "lz4", None, (2, 1, 0))
    meta = json.loads((out / ".zarray").read_text())
    assert meta["compressor"]["cname"] == "zstd" and meta["chunks"] == [4, 8, 8]
    assert meta["dtype"] == "<f4" and meta["shape"] == [5, 20, 18]


def test_existing_map_is_left_alone(tmp_path):
    bi = _write_tile_info(tmp_path, {"t1": make_tile((5, 20, 18), 14)})
    out = tmp_path / "m.zarr"
    out.mkdir()
    (out / "keep").write_text("x")
    sd.compute_tile_bwdist(bi, 0, str(out), 10, True, (4, 8, 8), None, "zstd", None, (2, 1, 0),
                           True)  # fmt: skip
    assert sorted(os.listdir(out)) == ["keep"]


@pytest.mark.parametrize("overwrite", [False, True])
def test_map_that_appears_meanwhile(tmp_path, monkeypatch, capsys, overwrite):
    bi = _write_tile_info(tmp_path, {"t1": make_tile((5, 20, 18), 15)})
    out = tmp_path / "m.zarr"
    real = sd.tile_distance_map

    def racing(*a, **k):
        out.mkdir()
        (out / "other").write_text("x")
        return real(*a, **k)

    monkeypatch.setattr(sd, "tile_distance_map", racing)
    sd.compute_tile_bwdist(bi, 0, str(out), 10, True, (4, 8, 8), None, "zstd", None, (2, 1, 0),
                           overwrite)  # fmt: skip
    left = sorted(os.listdir(tmp_path))
    assert [n for n in left if n.startswith("m.zarr_")] == []
    if overwrite:
        assert (out / ".zarray").is_file() and not (out / "other").exists()
    else:
        assert sorted(os.listdir(out)) == ["other"]
        assert "already exists, skip it!" in capsys.readouterr().out


def test_slab_separator_comes_from_the_matrix_division(tmp_path):
    # MATLAB's prod(ceil(sz / blockSize)): 10,400 chunks, but one quotient of 21, so '.'
    tiles = {"t1": make_tile((3, 5, 4), 16)}
    bi = _write_tile_info(tmp_path, tiles)
    slab_dir = tmp_path / "tiles" / "MIP_Slabs_pooling_2_2_2"
    slab_dir.mkdir()
    synth.write_zarr(slab_dir / "t1_z.zarr", make_tile((26, 40, 40), 16), block_size=(8, 8, 8))
    out = tmp_path / "m.zarr"
    sd.compute_tile_bwdist_mip_slabs(bi, 0, str(out), 10, True, (1, 2, 2), None, "zstd",
                                     (2, 2, 2), None, (2, 1, 0))  # fmt: skip
    assert json.loads((out / ".zarray").read_text()).get("dimension_separator", ".") == "."
    assert sd._mrdivide_rows([40, 40, 26], [2, 2, 1]) == pytest.approx(186 / 9)
    # compute_tile_bwdist divides element by element: '/'
    bi2 = _write_tile_info(tmp_path / "b", {"t1": make_tile((26, 40, 40), 16)})
    out2 = tmp_path / "m2.zarr"
    sd.compute_tile_bwdist(bi2, 0, str(out2), 10, True, (1, 2, 2), None, "zstd", None, (2, 1, 0))
    assert json.loads((out2 / ".zarray").read_text()).get("dimension_separator") == "/"


def test_slab_map_with_shards_is_zarr_v3(tmp_path):
    # the large-file maps take the stitching's shard size: the port writes Zarr v3 with
    # sharding (the user's decision; MATLAB writes a v2 .zarray with a sharding codec)
    bi = _write_tile_info(tmp_path, {"t1": make_tile((3, 5, 4), 18)})
    slab_dir = tmp_path / "tiles" / "MIP_Slabs_pooling_2_2_2"
    slab_dir.mkdir()
    synth.write_zarr(slab_dir / "t1_z.zarr", make_tile((9, 40, 36), 18, "dsr"))
    maps = []
    for name, shard in (("a.zarr", None), ("b.zarr", (4, 8, 8))):
        sd.compute_tile_bwdist_mip_slabs(bi, 0, str(tmp_path / name), 10, True, (8, 16, 16),
                                         shard, "zstd", (2, 2, 2), None, (2, 1, 0))  # fmt: skip
        maps.append(_read_map(tmp_path / name))
    assert (tmp_path / "b.zarr" / "zarr.json").is_file()
    np.testing.assert_array_equal(_bits(maps[0]), _bits(maps[1]))


def test_mip_slab_name():
    assert sd.mip_slab_name((2, 2, 2)) == "MIP_Slabs_pooling_2_2_2"
    assert sd.mip_slab_name((4, 2, 3, 1, 2, 2)) == "MIP_Slabs_pooling_2_3_4_2_2_1"


@pytest.fixture(params=["posix", "windows"])
def rule(request, monkeypatch):
    """the wrapper's folder rule: MATLAB's on Linux and macOS, or on Windows (ispc), forced on
    any system"""
    monkeypatch.setattr(sd, "_is_pc", lambda: request.param == "windows")
    return request.param


def _dist_folder(stitch, first, rule):
    """compute_tile_distance_transform's folder of the maps, with its trailing '/'"""
    return f"{stitch}/imdist/" if rule == "windows" else f"{stitch}/imdist/{first}/"


def _wrapper_tiles(tmp_path, sizes, seed=20):
    names = [f"t{k + 1}" for k in range(len(sizes))]
    tiles = {
        n: make_tile(s, seed + k, "dsr") for k, (n, s) in enumerate(zip(names, sizes, strict=True))
    }
    pairs = grid_pairs(sizes, 5, 5) if len(sizes) == 4 else ()
    bi = _write_tile_info(tmp_path, tiles, pairs)
    return bi, [str(tmp_path / "tiles" / f"{n}.zarr") for n in names]


@pytest.mark.parametrize(
    "loc_ids, sizes, n_maps, idx",
    [
        ([1, 1, 1, 1], [(6, 20, 18)] * 4, 1, [0]),
        ([1, 1, 2, 2], [(6, 20, 18)] * 4, 1, [0]),
        ([1, 2, 3, 1], [(6, 20, 18)] * 4, 1, [0]),
        ([1, 2, 1, 2], [(6, 20, 18), (6, 20, 16)] * 2, 2, [0, 1, 0, 1]),
        ([2, 1, 3, 1], [(6, 20, 18), (6, 20, 16), (6, 20, 18), (6, 20, 16)], 2, [1, 0, 1, 0]),
        ([3, 2, 1, 4], [(6, 20, 18), (6, 20, 16), (6, 18, 18), (6, 20, 16)], 3, [2, 1, 0, 1]),
    ],
)
def test_single_maps_per_location_and_size(tmp_path, rule, loc_ids, sizes, n_maps, idx):
    bi, paths = _wrapper_tiles(tmp_path, sizes)
    stitch = str(tmp_path / "stitch") + "/"
    outs, fidx = sd.compute_tile_distance_transform(
        bi, stitch, paths, loc_ids=loc_ids, parse_cluster=False
    )
    assert len(outs) == n_maps and fidx == idx
    assert all(os.path.isdir(o) for o in outs)
    # MATLAB's paths: '<stitch>/imdist/t1//<tile>_wr_10.zarr' ('<stitch>/imdist//...' on Windows)
    folder = _dist_folder(stitch, "t1", rule)
    names = {f"t{k + 1}_wr_10.zarr" for k in range(4)}
    assert all(o.startswith(folder + "/") and o[len(folder) + 1 :] in names for o in outs)
    assert os.path.isfile(folder + "parameters.mat")


def test_tiles_of_different_sizes_get_their_own_maps(tmp_path):
    # StitchingTest/tilesOfDifferentSizesGetTheirOwnDistanceMaps: without singleDistMap one map
    # per tile, from the tile's overlap regions
    sizes = [(6, 20, 18), (6, 20, 16), (6, 20, 18), (6, 20, 18)]
    bi, paths = _wrapper_tiles(tmp_path, sizes)
    outs, fidx = sd.compute_tile_distance_transform(
        bi, str(tmp_path / "stitch") + "/", paths, single_dist_map=False, parse_cluster=False
    )
    assert fidx == [0, 1, 2, 3] and len(outs) == 4
    for k, o in enumerate(outs):
        assert _read_map(o).shape == sizes[k]


def test_maps_of_other_overlaps_are_removed(tmp_path, capsys, rule):
    sizes = [(6, 20, 18)] * 4
    bi, paths = _wrapper_tiles(tmp_path, sizes)
    stitch = str(tmp_path / "stitch") + "/"
    stale = Path(_dist_folder(stitch, "t1", rule))
    stale.mkdir(parents=True)
    write_fixed_mat(stale / "parameters.mat", {"overlap_matrix": np.zeros((4, 4), bool),
                                               "ol_region_cell": obj_col([]),
                                               "overlap_map_mat": np.zeros((0, 2))})  # fmt: skip
    (stale / "old.zarr").mkdir()
    sd.compute_tile_distance_transform(bi, stitch, paths, single_dist_map=False,
                                       parse_cluster=False)  # fmt: skip
    assert "not in the same setting" in capsys.readouterr().out
    assert not (stale / "old.zarr").exists()
    # the same setting again: kept
    (stale / "old.zarr").mkdir()
    sd.compute_tile_distance_transform(bi, stitch, paths, single_dist_map=False,
                                       parse_cluster=False)  # fmt: skip
    assert (stale / "old.zarr").exists()


@pytest.mark.parametrize("single", [True, False])
def test_windows_folder_holds_the_same_maps(tmp_path, monkeypatch, single):
    # MATLAB's ispc branch: imdist/ itself, the paths '<stitch>/imdist//<tile>_wr_10.zarr', the
    # maps as in imdist/<first tile>/ (whose voxels the oracle cases compare with MATLAB's)
    sizes = [(6, 20, 18), (6, 20, 16), (6, 20, 18), (6, 20, 16)]
    bi, paths = _wrapper_tiles(tmp_path, sizes)
    kw = {"loc_ids": [1, 2, 1, 2]} if single else {"single_dist_map": False}
    got = {}
    for pc in (False, True):
        monkeypatch.setattr(sd, "_is_pc", lambda pc=pc: pc)
        stitch = str(tmp_path / f"stitch_{int(pc)}") + "/"
        outs, fidx = sd.compute_tile_distance_transform(bi, stitch, paths, parse_cluster=False,
                                                        **kw)  # fmt: skip
        folder = _dist_folder(stitch, "t1", "windows" if pc else "posix")
        assert all(o.startswith(folder + "/") for o in outs)
        assert sorted(os.listdir(folder)) == sorted(
            ["parameters.mat"] + [o[len(folder) + 1 :] for o in outs]
        )
        got[pc] = (fidx, [_read_map(o) for o in outs])
    assert got[True][0] == got[False][0]
    for a, b in zip(got[True][1], got[False][1], strict=True):
        np.testing.assert_array_equal(_bits(a), _bits(b))


def test_windows_volumes_share_the_folder(tmp_path, monkeypatch, capsys):
    # on Windows the volumes of one result folder share imdist/: a second volume's maps go next
    # to the first's, and a call with other overlap regions removes them all (MATLAB's rmdir)
    monkeypatch.setattr(sd, "_is_pc", lambda: True)
    stitch = str(tmp_path / "stitch") + "/"
    a = tmp_path / "a"
    bi_a, paths_a = _wrapper_tiles(a, [(6, 20, 18)] * 4)
    outs_a, _ = sd.compute_tile_distance_transform(bi_a, stitch, paths_a, loc_ids=[1, 1, 1, 1],
                                                   parse_cluster=False)  # fmt: skip
    b = tmp_path / "b"
    tiles = {f"u{k}": make_tile((6, 20, 18), 40 + k, "dsr") for k in range(4)}
    bi_b = _write_tile_info(b, tiles, grid_pairs([(6, 20, 18)] * 4, 3, 3))
    paths_b = [str(b / "tiles" / f"{n}.zarr") for n in tiles]
    outs_b, _ = sd.compute_tile_distance_transform(bi_b, stitch, paths_b, loc_ids=[1, 1, 1, 1],
                                                   parse_cluster=False)  # fmt: skip
    assert all(os.path.isdir(o) for o in outs_a + outs_b)
    assert sorted(os.listdir(stitch + "imdist")) == ["parameters.mat", "t1_wr_10.zarr",
                                                     "u0_wr_10.zarr"]  # fmt: skip
    # volume a again, per tile: its overlap regions differ from b's parameters.mat
    sd.compute_tile_distance_transform(bi_a, stitch, paths_a, single_dist_map=False,
                                       parse_cluster=False)  # fmt: skip
    assert "not in the same setting" in capsys.readouterr().out
    assert not any(os.path.isdir(o) for o in outs_b)


def test_large_file_maps_make_their_mip_slabs(tmp_path):
    # StitchingTest/largeFileWithoutRegistration: the distance maps of large files use the MIP
    # slabs, which used to be made only by the registration
    sizes = [(8, 20, 18)] * 4
    bi, paths = _wrapper_tiles(tmp_path, sizes)
    outs, fidx = sd.compute_tile_distance_transform(
        bi, str(tmp_path / "stitch") + "/", paths, loc_ids=[1, 1, 1, 1], large_file=True,
        pool_size=(2, 2, 2), parse_cluster=False,
    )  # fmt: skip
    slab = tmp_path / "tiles" / "MIP_Slabs_pooling_2_2_2" / "t1_z.zarr"
    assert slab.is_dir() and fidx == [0]
    assert outs[0].endswith("/t1_z_pooling_2_2_2_wr_1_for_10.zarr")
    assert _read_map(outs[0]).shape == _read_map(slab).shape


def test_single_map_needs_location_ids(tmp_path):
    bi, paths = _wrapper_tiles(tmp_path, [(6, 20, 18)] * 4)
    with pytest.raises(sd.DistanceError):
        sd.compute_tile_distance_transform(bi, str(tmp_path / "s") + "/", paths,
                                           parse_cluster=False)  # fmt: skip


def test_matlab_name_layer(tmp_path):
    from petakit5d import matlab

    bi = _write_tile_info(tmp_path, {"t1": make_tile((7, 20, 18), 17, "dsr")})
    a, b = tmp_path / "a.zarr", tmp_path / "b.zarr"
    matlab.compute_tile_bwdist(bi, 1, str(a), 10, True, [8, 8, 4], [], "zstd",
                               [3, 4, 2, 15, 14, 6], [2, 1, 3])  # fmt: skip
    sd.compute_tile_bwdist(bi, 0, str(b), 10, True, (4, 8, 8), None, "zstd",
                           (1, 2, 3, 6, 15, 14), (2, 1, 0))  # fmt: skip
    np.testing.assert_array_equal(_bits(_read_map(a)), _bits(_read_map(b)))
