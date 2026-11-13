"""The stitching kernels (petakit5d.kernels.stitch) against the mex functions and MATLAB's code.

Oracle cases (MATLAB R2024b runs the mex of the oracle tree on the same inputs; exact). Inputs
are made in MATLAB's form; the roles give the Python calls the port's form ((t, z, y, x) tile
stacks, (z, y, x) volumes, Python boxes) and map the results back:

- tile.*: tile_blending_3d_mex on TileBlendingTest's random batches (every method, uint8,
  uint16 and single with NaN, Inf and -0, 3D, 2D and singleton axes, 2 to 9 tiles, regions
  that overlap or are empty), a multithreaded batch, the per-voxel test batch, invalid input;
  tile.matlab.*: the same batches against MATLAB's code (tile_blending_3d of double tiles, as
  TileBlendingTest/mexMatchesTheMatlabCode compares them); tile.wrapper.*: other classes.
- feather.*: feather_blending_3d_mex and _with_indexing_mex (partial blocks: regions of a
  larger output, voxels without weight), the wrapper with both, other classes (MATLAB's code),
  weights outside the usual range (the mex's integer conversion), invalid input.
- resize.*: feather_distance_map_resize_3d_mex (sizes not multiples of the map's, 2D, singleton
  axes, non-integer sizes, weight degrees 0, 1, 2, 10, -2, 2.7), the wrapper.
- integral.*: integral_image_3d_mex (uint16, single, 2D, T = 1 and T = size), the wrapper's
  MATLAB code for double, uint8 and logical, invalid input.
- major.*: any_4th_dim_mex and check_major_tile_valid (StitchHelperMexTest's cases, NaN, other
  classes).
- nan.*: replace_nan_with_value (the value's first element as a double, in the input's class).

Cases of MATLAB's own argument class checks (a single or int64 box or size, a column of sizes)
are references only: the port's functions take Python numbers.

test_thread_counts runs every case again with 1 to 4 threads against the same reference.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, size_from_matlab, to_matlab
from petakit5d.core.bbox import box_from_matlab, box_to_matlab
from petakit5d.core.dtypes import NAN, to_class
from petakit5d.kernels import stitch as sk
from petakit5d.testing import oracle
from petakit5d.testing.cases import Box, Case, In, PerAxis, Value, Vol
from petakit5d.testing.synth import CounterRNG

M = "petakit5d.kernels.stitch:"
POW10 = np.array([0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0])
INF = float("inf")


# ---------------------------------------------------------------------------------------------
# inputs, in MATLAB's form (random values in MATLAB's linear order)
# ---------------------------------------------------------------------------------------------


def mfill(draw, shape):
    """a random array of MATLAB size shape, filled in MATLAB's linear order (draw(size) is a
    CounterRNG call on the port's size)"""
    return to_matlab(draw(size_from_matlab(shape)))


def random_batch(sz, n_tiles: int, dtype: str, seed: int):
    """TileBlendingTest.randomBatch: tiles at random parts of the batch (zeros elsewhere), 0 to 3
    regions per tile and an empty one; single values with NaN, Inf, -Inf and -0"""
    rng = CounterRNG(seed, f"batch/{dtype}")
    sz = list(sz) + [1] * (3 - len(sz))
    dt = {"uint8": np.uint8, "uint16": np.uint16, "single": np.float32}[dtype]
    tim = np.zeros(sz + [n_tiles], dtype=dt)
    boxes = []
    for j in range(n_tiles):
        if j == 0:
            s, e = [1, 1, 1], list(sz)
        else:
            s = [int(rng.integers(1, n, 1)[0]) for n in sz]
            e = [int(rng.integers(a, n, 1)[0]) for a, n in zip(s, sz, strict=True)]
        tsz = [b - a + 1 for a, b in zip(s, e, strict=True)]
        if dtype == "uint8":
            v = mfill(lambda q: rng.integers(0, 255, q), tsz).astype(np.float64)
        elif dtype == "uint16":
            v = mfill(lambda q: rng.integers(0, 65535, q), tsz).astype(np.float64)
        else:
            v = mfill(rng.normal, tsz) * POW10[mfill(lambda q: rng.integers(0, 6, q), tsz)]
        v[mfill(rng.uniform, tsz) < 0.3] = 0
        v = v.astype(dt)
        if dtype == "single":
            r = mfill(rng.uniform, tsz)
            v[r < 0.05] = NAN
            v[(r >= 0.05) & (r < 0.07)] = INF
            v[(r >= 0.07) & (r < 0.09)] = -INF
            v[(r >= 0.09) & (r < 0.12)] = -0.0
        tim[s[0] - 1 : e[0], s[1] - 1 : e[1], s[2] - 1 : e[2], j] = v
        for _ in range(int(rng.integers(0, 3, 1)[0])):
            ms = [int(rng.integers(a, b, 1)[0]) for a, b in zip(s, e, strict=True)]
            me = [int(rng.integers(a, b, 1)[0]) for a, b in zip(ms, e, strict=True)]
            boxes.append([j + 1, *ms, *me])
    boxes.append([n_tiles, 2, 1, 1, 1, sz[1], sz[2]])
    return tim, np.array(boxes, dtype=np.float64)


def voxel_test_batch():
    """TileBlendingTest.voxelsBlendTheTilesOutsideTheirRegions: 3 tiles of 10, 21 and 40 with
    zeros; tile 2 counts as 0 in columns 1-2 and tile 3 in column 1, rows 1-2"""
    tim = np.zeros((4, 4, 2, 3), dtype=np.uint16)
    tim[..., 0], tim[..., 1], tim[..., 2] = 10, 21, 40
    tim[0, :, :, 0] = 0
    tim[1:3, 2, 0, 1] = 0
    tim[3, 3, 1, :] = 0
    boxes = np.array([[2, 1, 1, 1, 4, 2, 2], [3, 1, 1, 1, 2, 1, 2]], dtype=np.float64)
    masked = np.zeros(tim.shape, dtype=bool)
    masked[:, 0:2, :, 1] = True
    masked[0:2, 0, :, 2] = True
    return tim, boxes, masked


def feather_tiles(sz, n_tiles: int, dtype: str, seed: int, weights: str = "normal"):
    """tiles with zero regions and single weights (distance-map-like, some 0); 'negative',
    'huge', 'inf' and 'nan' weights reach the mex's conversion outside [0, 65535]"""
    rng = CounterRNG(seed, f"feather/{dtype}/{weights}")
    shape = list(sz) + [1] * (3 - len(sz)) + [n_tiles]
    if dtype == "uint16":
        f = mfill(lambda q: rng.integers(0, 65535, q), shape).astype(np.uint16)
    else:
        f = (mfill(rng.normal, shape) * 1000.0).astype(np.float32)
        r = mfill(rng.uniform, shape)
        f[r < 0.03] = NAN
        f[(r >= 0.03) & (r < 0.05)] = -0.0
    f[mfill(rng.uniform, shape) < 0.25] = 0
    for j in range(1, n_tiles):  # each later tile misses a slab
        f[: shape[0] // (j + 1), :, :, j] = 0
    u = (mfill(rng.uniform, shape) * 2.0).astype(np.float32)
    d = u * u * u  # (IEEE products: the same bytes on every platform)
    d[mfill(rng.uniform, shape) < 0.1] = 0
    if weights == "negative":
        d = d - np.float32(1.0)
    elif weights == "huge":
        d[..., 0] = np.float32(1.0)
        d[..., 1:] = np.float32(-0.999999)
    elif weights == "inf":
        d[mfill(rng.uniform, shape) < 0.2] = INF
    elif weights == "nan":
        d[mfill(rng.uniform, shape) < 0.2] = NAN
    return f, d


def fvol(shape, dtype, seed: int, lo=0.0, hi=1.0, zeros: float = 0.0):
    """a random array of MATLAB size shape"""
    rng = CounterRNG(seed, f"vol/{np.dtype(dtype).name}")
    if np.dtype(dtype).kind == "f":
        v = lo + mfill(rng.uniform, shape) * (hi - lo)
    else:
        v = mfill(lambda q: rng.integers(int(lo), int(hi), q), shape).astype(np.float64)
    if zeros:
        v[mfill(rng.uniform, shape) < zeros] = 0
    return v.astype(dtype)


def V(make):
    return Value(make)


# ---------------------------------------------------------------------------------------------
# the Python side of cases whose MATLAB call takes MATLAB's tile regions or boxes
# ---------------------------------------------------------------------------------------------


def tile_mex(tiles, mboxes, method, num_threads=None):
    """tile_blending_3d_mex: the regions converted to the port's rows inside the call"""
    return sk.tile_blend(tiles, sk.tile_boxes_from_matlab(mboxes), method, num_threads=num_threads)


def tile_wrapper(tiles, mboxes, method, num_threads=None):
    return sk.tile_blending_3d(
        tiles, sk.tile_boxes_from_matlab(mboxes), method, num_threads=num_threads
    )


def tile_of_class(tiles, mboxes, method, cls, num_threads=None):
    """tile_blending_3d of the tiles converted to cls as MATLAB converts them"""
    t = np.ascontiguousarray(tiles != 0) if cls == "logical" else to_class(tiles, cls)
    return tile_wrapper(np.ascontiguousarray(t), mboxes, method, num_threads)


def feather_into(f, d, out, box, num_threads=None):
    sk.feather_blend_into(f, d, out, box, num_threads=num_threads)
    return out


def feather_into_mbox(f, d, out, mbox, num_threads=None):
    """feather_blending_3d_with_indexing_mex with MATLAB's box converted inside the call (a box
    that has no Python form raises there)"""
    return feather_into(f, d, out, box_from_matlab(mbox), num_threads)


def box_or_empty(b):
    """check_major_tile_valid's box as MATLAB's double row ([] when valid)"""
    return None if b is None else np.array([box_to_matlab(b)], dtype=np.float64)


def mex_bounds(b):
    """any_4th_dim_mex's bounds: uint64 [ymin xmin zmin ymax xmax zmax], zeros when valid"""
    return np.zeros((1, 6), np.uint64) if b is None else np.array([box_to_matlab(b)], np.uint64)


def single_of(r):
    return [r.astype(np.float32) if r.dtype != np.float32 else r]


# ---------------------------------------------------------------------------------------------
# tile blending
# ---------------------------------------------------------------------------------------------

METHODS = ["none", "max", "mean", "median"]
TILE_SIZES = [((23, 19, 11), 2), ((17, 13), 3), ((1, 9, 6), 5), ((23, 19, 11), 9)]
ERRORS = (sk.StitchKernelError, sk.StitchKernelTypeError)

CASES: list[Case] = []
for _dt in ["uint8", "uint16", "single"]:
    for _i, (_sz, _nt) in enumerate(TILE_SIZES):
        _seed = 100 * _i + {"uint8": 1, "uint16": 2, "single": 3}[_dt]
        _inp = {
            "tim": V(lambda sz=_sz, nt=_nt, dt=_dt, s=_seed: random_batch(sz, nt, dt, s)[0]),
            "boxes": V(lambda sz=_sz, nt=_nt, dt=_dt, s=_seed: random_batch(sz, nt, dt, s)[1]),
        }
        for _m in METHODS:
            CASES.append(
                Case(
                    f"tile.{_dt}.{_m}.b{_i}",
                    matlab="tile_blending_3d_mex",
                    args=[Vol(In("tim")), In("boxes"), _m],
                    python=tile_mex,
                    inputs=_inp,
                    out_roles=("vol",),
                )
            )
            # the MATLAB code (double tiles are exact copies of uint8, uint16 and single ones)
            CASES.append(
                Case(
                    f"tile.matlab.{_dt}.{_m}.b{_i}",
                    code="out = {tile_blending_3d(double(in{1}), in{2}, in{3})};",
                    args=[Vol(In("tim")), In("boxes"), _m],
                    python=tile_mex,
                    py_out=single_of,
                    inputs=_inp,
                    out_roles=("vol",),
                )
            )

for _dt in ["uint16", "single"]:  # a batch large enough for many threads
    _inp = {
        "tim": V(lambda dt=_dt: random_batch((130, 101, 24), 3, dt, 77)[0]),
        "boxes": V(lambda dt=_dt: random_batch((130, 101, 24), 3, dt, 77)[1]),
    }
    for _m in METHODS:
        CASES.append(
            Case(
                f"tile.{_dt}.{_m}.large",
                matlab="tile_blending_3d_mex",
                args=[Vol(In("tim")), In("boxes"), _m],
                python=tile_mex,
                inputs=_inp,
                out_roles=("vol",),
            )
        )

for _m in METHODS:
    CASES.append(
        Case(
            f"tile.voxels.{_m}",
            matlab="tile_blending_3d",
            args=[Vol(In("tim")), In("boxes"), _m],
            python=tile_wrapper,
            inputs={
                "tim": V(lambda: voxel_test_batch()[0]),
                "boxes": V(lambda: voxel_test_batch()[1]),
            },
            out_roles=("vol",),
        )
    )
    for _cls in ["double", "int16", "logical"]:  # MATLAB's code in single
        CASES.append(
            Case(
                f"tile.wrapper.{_cls}.{_m}",
                code=f"out = {{tile_blending_3d({_cls}(in{{1}}), in{{2}}, in{{3}})}};",
                args=[Vol(In("tim")), In("boxes"), _m],
                python=tile_of_class,
                py_kwargs={"cls": _cls},
                inputs={
                    "tim": V(lambda: random_batch((12, 10, 7), 3, "uint16", 5)[0]),
                    "boxes": V(lambda: random_batch((12, 10, 7), 3, "uint16", 5)[1]),
                },
                out_roles=("vol",),
            )
        )


def _tim4():
    return np.zeros((4, 4, 4, 2), dtype=np.uint16)


_TILE_ERRORS = [  # TileBlendingTest/mexRejectsInvalidInput
    ("tile_region", [[3, 1, 1, 1, 2, 2, 2]], "none"),
    ("zero_start", [[1, 0, 1, 1, 2, 2, 2]], "none"),
    ("past_end", [[1, 1, 1, 1, 5, 2, 2]], "none"),
    ("past_end_x", [[1, 1, 1, 1, 2, 5, 2]], "none"),
    ("past_end_z", [[1, 1, 1, 1, 2, 2, 5]], "none"),
    ("fraction", [[1, 1, 1, 1, 2.5, 2, 2]], "none"),
    ("six_columns", [[1, 1, 1, 1, 2, 2]], "none"),
    ("nan", [[1, 1, NAN, 1, 2, 2, 2]], "none"),
    ("inf", [[1, 1, 1, 1, INF, 2, 2]], "none"),
    ("feather", None, "feather"),
    ("long_name", None, "medianmedianmedian"),
]
for _name, _b, _m in _TILE_ERRORS:
    CASES.append(
        Case(
            f"tile.error.{_name}",
            matlab="tile_blending_3d_mex",
            args=[Vol(In("tim")), np.array(_b, dtype=np.float64) if _b else None, _m],
            python=tile_mex,
            inputs={"tim": V(_tim4)},
            expect_error="tile_blending:inputError",
            py_error=ERRORS,
            # the port's message names its own row form
            error_message=_name != "six_columns",
        )
    )
CASES += [
    Case(
        "tile.error.no_tiles",
        matlab="tile_blending_3d_mex",
        args=[Vol(In("tim")), None, "none"],
        python=tile_mex,
        inputs={"tim": V(lambda: np.zeros((4, 4, 4, 0), dtype=np.uint16))},
        expect_error="tile_blending:inputError",
        py_error=sk.StitchKernelError,
        error_message=True,
    ),
    Case(
        "tile.error.double",
        matlab="tile_blending_3d_mex",
        args=[Vol(In("tim")), None, "none"],
        python=tile_mex,
        inputs={"tim": V(lambda: np.zeros((4, 4, 4, 2)))},
        expect_error="tile_blending:dataTypeError",
        py_error=sk.StitchKernelTypeError,
        error_message=True,
    ),
    Case(
        "tile.error.wrapper_method",
        matlab="tile_blending_3d",
        args=[Vol(In("tim")), None, "None"],
        python=tile_wrapper,
        inputs={"tim": V(_tim4)},
        expect_error=True,
        py_error=ValueError,
        error_message=True,
    ),
    Case(
        "tile.empty_regions_only",
        matlab="tile_blending_3d_mex",
        args=[
            Vol(In("tim")),
            np.array([[1, 3, 1, 1, 2, 4, 4], [2, 1, 2, 1, 4, 1, 4]], float),
            "max",
        ],
        python=tile_mex,
        inputs={"tim": V(lambda: random_batch((4, 4, 4), 2, "uint16", 9)[0])},
        out_roles=("vol",),
    ),
]

# ---------------------------------------------------------------------------------------------
# feather blending
# ---------------------------------------------------------------------------------------------

FEATHER = [  # id, size, tiles
    ("b3d", (21, 17, 9), 3),
    ("b2d", (19, 23), 2),
    ("single_tile", (8, 6, 5), 1),
    ("many", (11, 9, 7), 6),
    ("large", (96, 80, 40), 4),
]
for _dt in ["uint16", "single"]:
    for _k, (_id, _sz, _nt) in enumerate(FEATHER):
        _inp = {
            "f": V(lambda sz=_sz, nt=_nt, dt=_dt, s=_k: feather_tiles(sz, nt, dt, 10 + s)[0]),
            "d": V(lambda sz=_sz, nt=_nt, dt=_dt, s=_k: feather_tiles(sz, nt, dt, 10 + s)[1]),
        }
        CASES.append(
            Case(
                f"feather.{_dt}.{_id}",
                matlab="feather_blending_3d_mex",
                args=[Vol(In("f")), Vol(In("d"))],
                python=M + "feather_blend",
                inputs=_inp,
                out_roles=("vol",),
            )
        )
for _w in ["negative", "huge", "inf", "nan"]:
    for _dt in ["uint16", "single"]:
        CASES.append(
            Case(
                f"feather.{_dt}.weights_{_w}",
                matlab="feather_blending_3d_mex",
                args=[Vol(In("f")), Vol(In("d"))],
                python=M + "feather_blend",
                inputs={
                    "f": V(lambda w=_w, dt=_dt: feather_tiles((13, 11, 5), 3, dt, 31, w)[0]),
                    "d": V(lambda w=_w, dt=_dt: feather_tiles((13, 11, 5), 3, dt, 31, w)[1]),
                },
                out_roles=("vol",),
            )
        )

_WITH_INDEXING = (
    "o = in{3}; feather_blending_3d_with_indexing_mex(in{1}, in{2}, o, in{4}); out = {o};"
)
PARTIAL = [  # id, tile size, tiles, output size, bbox
    ("inner", (9, 7, 5), 3, (20, 18, 12), [4, 6, 3, 12, 12, 7]),
    ("corner_end", (5, 6, 4), 2, (12, 10, 9), [8, 5, 6, 12, 10, 9]),
    ("full", (10, 8, 6), 2, (10, 8, 6), [1, 1, 1, 10, 8, 6]),
    ("plane", (7, 9), 3, (15, 13, 4), [3, 2, 4, 9, 10, 4]),
    ("large", (64, 50, 30), 4, (90, 70, 41), [11, 21, 12, 74, 70, 41]),
]
for _dt in ["uint16", "single"]:
    for _k, (_id, _sz, _nt, _osz, _bb) in enumerate(PARTIAL):
        _inp = {
            "f": V(lambda sz=_sz, nt=_nt, dt=_dt, s=_k: feather_tiles(sz, nt, dt, 40 + s)[0]),
            "d": V(lambda sz=_sz, nt=_nt, dt=_dt, s=_k: feather_tiles(sz, nt, dt, 40 + s)[1]),
            "o": V(
                lambda osz=_osz, dt=_dt, s=_k: fvol(
                    osz, np.uint16 if dt == "uint16" else np.float32, 50 + s, 1, 60000
                )
            ),
        }
        _args = [Vol(In("f")), Vol(In("d")), Vol(In("o")), Box(_bb)]
        CASES.append(
            Case(
                f"feather.indexing.{_dt}.{_id}",
                code=_WITH_INDEXING,
                args=_args,
                python=feather_into,
                inputs=_inp,
                out_roles=("vol",),
            )
        )
        CASES.append(
            Case(
                f"feather.wrapper.indexing.{_dt}.{_id}",
                matlab="feather_blending_3d",
                args=_args,
                nargout=2,
                python=M + "feather_blending_3d",
                inputs=_inp,
                out_roles=("vol", None),
            )
        )
for _dt in ["uint16", "single"]:
    CASES.append(
        Case(
            f"feather.wrapper.{_dt}",
            matlab="feather_blending_3d",
            args=[Vol(In("f")), Vol(In("d"))],
            nargout=2,
            python=M + "feather_blending_3d",
            inputs={
                "f": V(lambda dt=_dt: feather_tiles((14, 12, 6), 3, dt, 61)[0]),
                "d": V(lambda dt=_dt: feather_tiles((14, 12, 6), 3, dt, 61)[1]),
            },
            out_roles=("vol", None),
        )
    )


def _feather_f(fc):
    if fc == "single":
        return feather_tiles((9, 8, 5), 3, "single", 71)[0]
    return feather_tiles((9, 8, 5), 3, "uint16", 71)[0].astype(fc)


# MATLAB's code: other tile classes, double weights; with an output, the blend cast to the
# block's class and written into it (ids feather.wrapper.matlab.<pair>, and with a block of each
# class, .indexing_single and .indexing_uint8, .indexing_double, .indexing_uint16:
# FeatherBlendingTest's cases). Each block is an input of one case only: MATLAB's code writes
# into the block in place, and the cases of a MATLAB pack share their inputs.
for _p, _pair in enumerate(
    ["uint8_single", "double_single", "uint16_double", "single_double", "uint16_single"]
):
    _fc, _dc = _pair.split("_")  # the tiles' and the weights' classes
    _inp = {
        "f": V(lambda fc=_fc: _feather_f(fc)),
        "d": V(lambda dc=_dc: feather_tiles((9, 8, 5), 3, "uint16", 71)[1].astype(dc)),
    }
    CASES.append(
        Case(
            f"feather.wrapper.matlab.{_fc}_{_dc}",
            matlab="feather_blending_3d",
            args=[Vol(In("f")), Vol(In("d"))],
            nargout=2,
            python=M + "feather_blending_3d",
            inputs=_inp,
            out_roles=("vol", None),
        )
    )
    for _oc in sorted({"single", _fc}):
        CASES.append(
            Case(
                f"feather.wrapper.matlab.{_fc}_{_dc}.indexing_{_oc}",
                matlab="feather_blending_3d",
                args=[Vol(In("f")), Vol(In("d")), Vol(In("o")), Box([2, 3, 1, 10, 10, 5])],
                nargout=2,
                python=M + "feather_blending_3d",
                inputs={
                    **_inp,
                    "o": V(lambda oc=_oc, p=_p: fvol((12, 11, 6), oc, 72 + 10 * p, 1, 200)),
                },
                out_roles=("vol", None),
            )
        )

_FB_ID = "feather_blending:inputError"
_FEATHER_ERRORS = [  # MexInputCheckTest/stitchMexCheckTheirInputs
    ("5d", "f5", "d5", _FB_ID),
    ("dmap_more_dims", "f3", "d", _FB_ID),
    ("dmap_fewer_dims", "f", "d3", _FB_ID),
    ("dmap_double", "f", "dd", _FB_ID),
    ("dmap_uint16", "f", "du", _FB_ID),
    ("dmap_size", "f", "dsz", _FB_ID),
    ("double_tiles", "fd", "d", "feather_blending:dataTypeError"),
]
_FEATHER_ERROR_INPUTS = {
    "f": V(lambda: np.ones((4, 3, 2, 2), dtype=np.uint16)),
    "d": V(lambda: np.ones((4, 3, 2, 2), dtype=np.float32)),
    "f5": V(lambda: np.ones((4, 3, 2, 2, 2), dtype=np.uint16)),
    "d5": V(lambda: np.ones((4, 3, 2, 2, 2), dtype=np.float32)),
    "f3": V(lambda: np.ones((4, 3, 2), dtype=np.uint16)),
    "d3": V(lambda: np.ones((4, 3, 2), dtype=np.float32)),
    "dd": V(lambda: np.ones((4, 3, 2, 2))),
    "du": V(lambda: np.ones((4, 3, 2, 2), dtype=np.uint16)),
    "dsz": V(lambda: np.ones((4, 4, 2, 2), dtype=np.float32)),
    "fd": V(lambda: np.ones((4, 3, 2, 2))),
    "o": V(lambda: np.zeros((10, 10, 10), dtype=np.uint16)),
    "o4": V(lambda: np.zeros((10, 10, 10, 2), dtype=np.uint16)),
    "os": V(lambda: np.zeros((10, 10, 10), dtype=np.float32)),
}
for _name, _f, _d, _id in _FEATHER_ERRORS:
    CASES.append(
        Case(
            f"feather.error.{_name}",
            matlab="feather_blending_3d_mex",
            args=[Vol(In(_f)), Vol(In(_d))],
            python=M + "feather_blend",
            inputs={k: _FEATHER_ERROR_INPUTS[k] for k in {_f, _d}},
            expect_error=_id,
            py_error=ERRORS,
            error_message=True,
        )
    )
_INDEXING_ERRORS = [  # name, tiles, weights, output, MATLAB's box; Python: True, msg, refs only
    ("5d", "f5", "d5", "o", [1, 1, 1, 4, 3, 2], True),
    ("dmap_fewer_dims", "f", "d3", "o", [1, 1, 1, 4, 3, 2], True),
    ("4d_output", "f", "d", "o4", [1, 1, 1, 4, 3, 2], True),
    ("single_bbox", "f", "d", "o", np.array([1, 1, 1, 4, 3, 2], np.float32), None),
    ("int64_bbox", "f", "d", "o", np.array([1, 1, 1, 4, 3, 2], np.int64), None),
    ("past_the_rows", "f", "d", "o", [8, 1, 1, 11, 3, 2], True),
    ("nan_bbox", "f", "d", "o", [NAN, 1, 1, 4, 3, 2], False),
    ("five_values", "f", "d", "o", [1, 1, 1, 4, 3], False),
    ("wrong_size", "f", "d", "o", [1, 1, 1, 5, 3, 2], True),
    ("negative_end", "f", "d", "o", [1, 1, 1, -4, 3, 2], True),
    ("output_class", "f", "d", "os", [1, 1, 1, 4, 3, 2], True),
    ("dmap_double", "f", "dd", "o", [1, 1, 1, 4, 3, 2], True),
]
for _name, _f, _d, _o, _bb, _py in _INDEXING_ERRORS:
    CASES.append(
        Case(
            f"feather.error.indexing_{_name}",
            code="feather_blending_3d_with_indexing_mex(in{:}); out = {0};",
            args=[Vol(In(_f)), Vol(In(_d)), Vol(In(_o)), _bb],
            # None: MATLAB's class check of the box, which the port's Python box has no form for
            python=None if _py is None else feather_into_mbox,
            inputs={k: _FEATHER_ERROR_INPUTS[k] for k in {_f, _d, _o}},
            expect_error=_FB_ID,
            # a box without a Python form (NaN, 5 values) raises in its conversion
            py_error=ValueError,
            error_message=bool(_py),
        )
    )
CASES.append(  # the region is written at bbox of the output
    Case(
        "feather.indexing.ones_at_end",
        code=_WITH_INDEXING,
        args=[Vol(In("f")), Vol(In("d")), Vol(In("o")), Box([7, 8, 9, 10, 10, 10])],
        python=feather_into,
        inputs={k: _FEATHER_ERROR_INPUTS[k] for k in ("f", "d", "o")},
        out_roles=("vol",),
    )
)

# ---------------------------------------------------------------------------------------------
# distance map resize
# ---------------------------------------------------------------------------------------------

RESIZE = [  # id, map size, size, wd
    ("up", (5, 7, 3), [17, 20, 9], 2),
    ("wd1", (6, 5, 4), [11, 13, 10], 1),
    ("wd10", (4, 6, 5), [9, 16, 12], 10),
    ("wd0", (4, 4, 3), [7, 9, 5], 0),
    ("wd_negative", (5, 4, 3), [8, 8, 7], -2),
    ("wd_fraction", (5, 4, 3), [8, 8, 7], 2.7),
    ("same", (6, 7, 5), [6, 7, 5], 3),
    ("by_one", (6, 7, 5), [7, 8, 6], 2),
    ("plane", (6, 8), [11, 15, 1], 2),
    ("plane_to_volume", (6, 8), [11, 15, 4], 2),
    ("row", (1, 5, 4), [1, 9, 7], 2),
    ("fractional_size", (4, 3, 2), [8.7, 3.2, 2.5], 2),
    ("large", (40, 33, 20), [256, 201, 77], 10),
]
for _k, (_id, _msz, _sz, _wd) in enumerate(RESIZE):
    _inp = {"dmat": V(lambda msz=_msz, s=_k: fvol(msz, np.float32, 80 + s, 0.05, 1.5, 0.05))}
    CASES.append(
        Case(
            f"resize.{_id}",
            matlab="feather_distance_map_resize_3d_mex",
            args=[Vol(In("dmat")), PerAxis(_sz, "yxz"), _wd],
            python=M + "resize_distance_map",
            inputs=_inp,
            out_roles=("vol",),
        )
    )
CASES += [
    Case(
        "resize.wrapper",
        matlab="feather_distance_map_resize_3d",
        args=[Vol(In("dmat")), Box([3, 5, 2, 19, 21, 14]), 10],
        python=M + "feather_distance_map_resize_3d",
        inputs={"dmat": V(lambda: fvol((8, 9, 6), np.float32, 95, 0.1, 1.3))},
        out_roles=("vol",),
    ),
    Case(
        "resize.wrapper.plane",
        matlab="feather_distance_map_resize_3d",
        args=[Vol(In("dmat")), Box([1, 1, 4, 30, 25, 4]), 2],
        python=M + "feather_distance_map_resize_3d",
        inputs={"dmat": V(lambda: fvol((10, 9), np.float32, 96, 0.1, 1.3))},
        out_roles=("vol",),
    ),
    Case(
        "resize.ones",
        matlab="feather_distance_map_resize_3d_mex",
        args=[Vol(In("dmat")), PerAxis([8, 3, 2.5], "yxz"), 2],
        python=M + "resize_distance_map",
        inputs={"dmat": V(lambda: np.ones((4, 3, 2), dtype=np.float32))},
        out_roles=("vol",),
    ),
    Case(
        "resize.empty",
        matlab="feather_distance_map_resize_3d_mex",
        args=[Vol(In("dmat")), PerAxis([0, 3, 2], "yxz"), 2],
        python=M + "resize_distance_map",
        inputs={"dmat": V(lambda: np.ones((0, 3, 2), dtype=np.float32))},
        out_roles=("vol",),
    ),
]
_RESIZE_ERRORS = [  # name, map, size, identifier, Python: compared (False: reference only)
    ("4d", "d4", [4, 3, 2], _FB_ID, True),
    ("single_size", "d", np.array([4, 3, 2], np.float32), "distmap_resize:inputError", False),
    ("int64_size", "d", np.array([4, 3, 2], np.int64), "distmap_resize:inputError", False),
    ("nan_size", "d", [4, NAN, 2], "distmap_resize:inputError", True),
    ("empty_map", "d0", [4, 3, 2], "distmap_resize:inputError", True),
    ("smaller", "d", [4, 2, 2], "distmap_resize:inputError", True),
    ("column", "d", np.array([[4], [3], [2]], float), "crop:inputError", False),
    ("double_map", "dd", [4, 3, 2], "feather_blending:dataTypeError", True),
]
_RESIZE_ERROR_INPUTS = {
    "d": V(lambda: np.ones((4, 3, 2), dtype=np.float32)),
    "d4": V(lambda: np.ones((4, 3, 2, 2), dtype=np.float32)),
    "d0": V(lambda: np.ones((0, 3, 2), dtype=np.float32)),
    "dd": V(lambda: np.ones((4, 3, 2))),
}
for _name, _in, _sz, _id, _py in _RESIZE_ERRORS:
    CASES.append(
        Case(
            f"resize.error.{_name}",
            matlab="feather_distance_map_resize_3d_mex",
            args=[Vol(In(_in)), PerAxis(_sz, "yxz") if _py else _sz, 2],
            python=M + "resize_distance_map" if _py else None,
            inputs={_in: _RESIZE_ERROR_INPUTS[_in]},
            expect_error=_id,
            py_error=ERRORS,
            error_message=True,
        )
    )

# ---------------------------------------------------------------------------------------------
# integral image
# ---------------------------------------------------------------------------------------------

INTEGRAL = [  # id, size, class, T
    ("u16", (10, 12, 8), "uint16", [3, 4, 2]),
    ("u16_large", (64, 70, 33), "uint16", [17, 9, 12]),
    ("single", (11, 9, 13), "single", [4, 3, 5]),
    ("single_large", (60, 47, 38), "single", [21, 16, 9]),
    ("t_ones", (7, 6, 5), "single", [1, 1, 1]),
    ("t_full", (7, 6, 5), "single", [7, 6, 5]),
    ("plane", (9, 7), "single", [3, 2, 1]),
    ("fraction", (9, 8, 6), "single", [3.9, 2.2, 1.5]),
]
for _k, (_id, _sz, _cls, _t) in enumerate(INTEGRAL):
    _dt = np.uint16 if _cls == "uint16" else np.float32
    _hi = 65535 if _cls == "uint16" else 1000.0
    CASES.append(
        Case(
            f"integral.{_id}",
            matlab="integral_image_3d_mex",
            args=[Vol(In("a")), PerAxis(_t, "yxz")],
            python=M + "integral_image",
            inputs={"a": V(lambda sz=_sz, dt=_dt, hi=_hi, s=_k: fvol(sz, dt, 110 + s, 0, hi))},
            out_roles=("vol",),
        )
    )
CASES.append(
    Case(
        "integral.single_signed",
        matlab="integral_image_3d",
        args=[Vol(In("a")), PerAxis([5, 3, 4], "yxz")],
        python=M + "integral_image_3d",
        inputs={"a": V(lambda: fvol((20, 17, 12), np.float32, 120, -1e4, 1e4))},
        out_roles=("vol",),
    )
)
for _cls, _dt, _hi in [
    ("double", np.float64, 1000.0),
    ("uint8", np.uint8, 255),
    ("logical", np.bool_, 1),
]:
    CASES.append(
        Case(
            f"integral.wrapper.{_cls}",
            matlab="integral_image_3d",
            args=[Vol(In("a")), PerAxis([3, 4, 2], "yxz")],
            python=M + "integral_image_3d",
            inputs={"a": V(lambda dt=_dt, hi=_hi: fvol((10, 12, 8), dt, 121, 0, hi))},
            out_roles=("vol",),
        )
    )
CASES.append(
    Case(
        "integral.wrapper.int16",
        matlab="integral_image_3d",
        args=[Vol(In("a")), PerAxis([3, 4, 2], "yxz")],
        python=M + "integral_image_3d",
        inputs={"a": V(lambda: fvol((10, 12, 8), np.int16, 122, -100, 100))},
        expect_error=True,
        py_error=sk.StitchKernelTypeError,
    )
)
_II_ID = "integral_image_3d:inputError"
_INTEGRAL_ERRORS = [  # StitchHelperMexTest/integralImageMexErrorsHaveIdentifiers, MexInputCheckTest
    ("double", "a", [3, 3, 3], "integral_image_3d:dataTypeError", True),
    ("two_sizes", "a", [3, 3], _II_ID, True),
    ("zero_size", "a", [3, 0, 3], _II_ID, True),
    ("4d", "a4", [3, 3, 3], _II_ID, True),
    ("smaller_than_t", "a1", [5, 1, 1], _II_ID, True),
    ("empty", "a0", [1, 1, 1], _II_ID, True),
    ("single_t", "a1", np.array([2, 2, 2], np.float32), _II_ID, False),
    ("int64_t", "a1", np.array([2, 2, 2], np.int64), _II_ID, False),
    ("nan_t", "a1", [2, NAN, 2], _II_ID, True),
    # y breaks one check and z another: the mex reports y's (it checks y, x, z)
    ("order", "a1", [0, 1, 9], _II_ID, True),
]
_INTEGRAL_ERROR_INPUTS = {
    "a": V(lambda: fvol((10, 12, 8), np.float64, 130)),
    "a4": V(lambda: fvol((10, 12, 8, 2), np.float32, 131)),
    "a1": V(lambda: np.ones((4, 4, 4), dtype=np.float32)),
    "a0": V(lambda: np.ones((4, 0, 4), dtype=np.float32)),
}
for _name, _in, _t, _id, _py in _INTEGRAL_ERRORS:
    CASES.append(
        Case(
            f"integral.error.{_name}",
            matlab="integral_image_3d_mex",
            args=[Vol(In(_in)), PerAxis(_t, "yxz") if _py else _t],
            python=M + "integral_image" if _py else None,
            inputs={_in: _INTEGRAL_ERROR_INPUTS[_in]},
            expect_error=_id,
            py_error=ERRORS,
            error_message=True,
        )
    )
CASES.append(
    Case(
        "integral.ones_full_window",
        matlab="integral_image_3d_mex",
        args=[Vol(In("a1")), PerAxis([4, 4, 4], "yxz")],
        python=M + "integral_image",
        inputs={"a1": _INTEGRAL_ERROR_INPUTS["a1"]},
        out_roles=("vol",),
    )
)

# ---------------------------------------------------------------------------------------------
# the major-tile check
# ---------------------------------------------------------------------------------------------


def _major_full(cls):
    """StitchHelperMexTest/majorTilesWithDataEverywhereAreValid"""
    f = np.zeros((8, 9, 5, 4), dtype=cls)
    f[..., 0] = 7
    f[:, :, 0:2, 3] = 3
    return f


def _major_hole(cls, fill_third: bool = False):
    """StitchHelperMexTest/voxelWithoutDataInTheMajorTilesIsFound"""
    f = np.zeros((8, 9, 5, 3), dtype=cls)
    f[..., 0] = 7
    f[..., 2] = 5
    f[2, 5, 1, 0] = 0
    f[2, 5, 1, 2] = 1 if fill_third else 0
    return f


_BOUNDS_HOLES = [(1, 120, 30), (256, 130, 33), (120, 1, 40), (125, 256, 20), (110, 115, 1)]


def _major_bounds():
    """StitchHelperMexTest/boundsOfVoxelsWithoutDataWithManyThreads"""
    f = np.full((256, 256, 64, 2), 7, dtype=np.uint16)
    f[99:150, 99:150, 9:50, :] = 0
    for y, x, z in [*_BOUNDS_HOLES, (111, 116, 64)]:
        f[y - 1, x - 1, z - 1, :] = 0
    return f


def _major_nan(cls):
    f = _major_hole(cls)
    f[2, 5, 1, 2] = NAN  # NaN is data for the mex, no data for MATLAB's any
    return f


_MAJOR = []
for _cls in ["uint16", "single", "uint8", "double"]:
    _MAJOR += [
        (f"full.{_cls}", lambda c=_cls: _major_full(c), [True, False, True, False]),
        (f"hole.{_cls}", lambda c=_cls: _major_hole(c), [True, False, True]),
        (f"filled.{_cls}", lambda c=_cls: _major_hole(c, True), [True, False, True]),
        (f"no_major.{_cls}", lambda c=_cls: _major_full(c), [False, False, False, False]),
    ]
_MAJOR += [
    ("nan.single", lambda: _major_nan("single"), [True, False, True]),
    ("nan.double", lambda: _major_nan("double"), [True, False, True]),
    ("bounds", _major_bounds, [True, True]),
    ("plane", lambda: np.ascontiguousarray(_major_hole("uint16")[:, :, 1:2, :]), [1, 1, 0]),
]
for _name, _mk, _major in _MAJOR:
    CASES.append(
        Case(
            f"major.{_name}",
            matlab="check_major_tile_valid",
            args=[Vol(In("f")), np.array(_major, dtype=bool)],
            nargout=2,
            python=M + "check_major_tile_valid",
            inputs={"f": V(_mk)},
            out_roles=(None, box_or_empty),
        )
    )
    mex_class = _name.split(".")[-1] in ("uint16", "single") or _name in ("bounds", "plane")
    if mex_class and any(_major):
        CASES.append(
            Case(
                f"major.mex.{_name}",
                matlab="any_4th_dim_mex",
                args=[Vol(In("f")), np.array(_major, dtype=bool)],
                nargout=2,
                python=M + "major_tiles_valid",
                inputs={"f": V(_mk)},
                out_roles=(None, mex_bounds),
            )
        )


def _ones4(shape=(4, 3, 2, 2), cls=np.uint16):
    return np.ones(shape, cls)


_A4 = "any_4th_dim_mex:"
_A4_ERRORS = [
    ("5d", lambda: _ones4((4, 3, 2, 2, 2)), np.ones(2, bool), _A4 + "inputError"),
    ("count", _ones4, np.ones(3, bool), _A4 + "dataTypeError"),
    ("double_flags", _ones4, [1.0, 1.0], _A4 + "dataTypeError"),
    ("no_major", _ones4, np.zeros(2, bool), _A4 + "dataRangeError"),
    ("double", lambda: _ones4(cls=np.float64), np.ones(2, bool), _A4 + "dataTypeError"),
]
for _name, _mk, _flags, _id in _A4_ERRORS:
    CASES.append(
        Case(
            f"major.error.{_name}",
            matlab="any_4th_dim_mex",
            args=[Vol(In("f")), _flags],
            nargout=2,
            python=M + "major_tiles_valid",
            inputs={"f": V(_mk)},
            expect_error=_id,
            py_error=ERRORS,
            error_message=True,
        )
    )

# ---------------------------------------------------------------------------------------------
# NaN replacement
# ---------------------------------------------------------------------------------------------


def _with_nans(cls, seed):
    v = fvol((17, 13, 5), cls, seed, -10, 10)
    r = fvol((17, 13, 5), np.float64, seed + 1)
    v[r < 0.2] = NAN
    v[(r >= 0.2) & (r < 0.25)] = INF
    v[(r >= 0.25) & (r < 0.3)] = -0.0
    return v


_NAN_CASES = [
    ("single_zero", np.float32, 0),
    ("double_five", np.float64, 5),
    ("single_double_five", np.float32, 5),
    ("single_single_five", np.float32, np.float32(5)),
    ("single_double_pi", np.float32, np.pi),
    ("single_uint8_bytes", np.float32, np.array([1, 2, 3, 200], dtype=np.uint8)),
    ("double_single_pair", np.float64, np.array([1.5, -2.25], dtype=np.float32)),
    ("single_minus_zero", np.float32, np.float32(-0.0)),
    # MexInputCheckTest/replaceNanTakesADoubleValueForSingleInput, singleVolumeGetsTheValue...
    ("single_uint8_seven", np.float32, np.uint8(7)),
    ("single_minus_1e8", np.float32, -1e8),
    ("single_logical", np.float32, True),
    ("double_int16", np.float64, np.int16(-3)),
    ("double_row_first", np.float64, [2.5, 3.5]),
]
for _k, (_name, _dt, _val) in enumerate(_NAN_CASES):
    CASES.append(
        Case(
            f"nan.{_name}",
            matlab="replace_nan_with_value",
            args=[Vol(In("x")), _val],
            python=M + "replace_nan_with_value",
            inputs={"x": V(lambda dt=_dt, s=_k: _with_nans(dt, 140 + 2 * s))},
            out_roles=("vol",),
        )
    )
CASES += [  # rows of values: element-wise, no mapping
    Case(
        "nan.single_value_double_input",
        matlab="replace_nan_with_value",
        args=[In("x"), np.float32(5)],
        python=M + "replace_nan_with_value",
        inputs={"x": V(lambda: np.array([1.0, NAN]))},
    ),
    Case(
        "nan.empty_input",
        matlab="replace_nan_with_value",
        args=[In("x"), 1],
        python=M + "replace_nan_with_value",
        inputs={"x": V(lambda: np.zeros((0, 3)))},
    ),
]
_NAN_ERRORS = [  # MexInputCheckTest/pointDetectionMexCheckTheirInputs
    ("empty_value", np.array([1.0, 2.0]), None, "data:inputError"),
    ("char_value", np.array([1.0, NAN]), "a", "data:inputError"),
    ("complex_value", np.array([1.0, NAN], dtype=np.float32), 1j, "data:inputError"),
    ("uint16_input", np.array([1, 2], dtype=np.uint16), 0, "data:TypeError"),
]
for _name, _x, _val, _id in _NAN_ERRORS:
    CASES.append(
        Case(
            f"nan.error.{_name}",
            matlab="replace_nan_with_value",
            args=[In("x"), _val],  # rows of values: element-wise, no mapping
            python=M + "replace_nan_with_value",
            inputs={"x": V(lambda x=_x: x.copy())},
            expect_error=_id,
            py_error=ERRORS,
            error_message=_id == "data:inputError",
        )
    )


# ---------------------------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------------------------


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


_THREADED = [
    c
    for c in CASES
    if not c.expect_error and c.python is not None and not c.id.startswith("tile.matlab.")
]


@oracle.parametrize(_THREADED)
@pytest.mark.parametrize("num_threads", [1, 2, 3, 4])
def test_thread_counts(case, num_threads, tmp_path):
    """the same results with every thread count (the kernels' own split among threads)"""
    kw = dict(case.py_kwargs or {})
    kw["num_threads"] = num_threads
    variant = dataclasses.replace(case, py_kwargs=kw)
    variant._module_globals = case._module_globals
    oracle.check(variant, tmp_path)


# --- TileBlendingTest ---------------------------------------------------------------------------


@pytest.mark.parametrize("method", METHODS)
def test_voxels_blend_the_tiles_outside_their_regions(method):
    """TileBlendingTest/voxelsBlendTheTilesOutsideTheirRegions: each voxel by its definition"""
    mtim, mboxes, mmasked = voxel_test_batch()
    tim, masked = from_matlab(mtim), from_matlab(mmasked)  # (t, z, y, x)
    boxes = sk.tile_boxes_from_matlab(mboxes)
    np.testing.assert_array_equal(boxes, [[1, 0, 0, 0, 2, 4, 2], [2, 0, 0, 0, 2, 2, 1]])
    out = sk.tile_blending_3d(tim, boxes, method)
    expected = np.zeros(tim.shape[1:], dtype=np.float32)
    for idx in np.ndindex(*expected.shape):
        v = tim[(slice(None), *idx)].astype(np.float32)
        sel = (v != 0) & ~masked[(slice(None), *idx)]
        if not sel.any():
            sel = v != 0
        if sel.any():
            vals = np.sort(v[sel]) if method == "median" else v[sel]
            n = vals.size
            if method == "none":
                expected[idx] = vals[0]
            elif method == "max":
                expected[idx] = vals.max()
            elif method == "mean":
                expected[idx] = vals.sum(dtype=np.float32) / np.float32(n)
            elif n % 2:
                expected[idx] = vals[n // 2]
            else:
                expected[idx] = (vals[n // 2 - 1] + vals[n // 2]) / np.float32(2)
    np.testing.assert_array_equal(out.astype(np.float32), expected)
    # the kernel keeps the class of the tiles for 'none' and 'max'
    assert out.dtype == (np.uint16 if method in ("none", "max") else np.float32)


def test_other_classes_use_the_matlab_code():
    """TileBlendingTest/otherClassesUseTheMatlabCode: double tiles are blended in single, as
    single tiles"""
    mtim, mboxes = random_batch((12, 10, 7), 3, "uint16", 5)
    tim, boxes = from_matlab(mtim), sk.tile_boxes_from_matlab(mboxes)
    for method in METHODS:
        out = sk.tile_blending_3d(tim.astype(np.float64), boxes, method)
        ref = sk.tile_blending_3d(tim.astype(np.float32), boxes, method)
        assert out.dtype == np.float32
        np.testing.assert_array_equal(out.view(np.uint32), ref.view(np.uint32))


def test_mex_rejects_invalid_input():
    """TileBlendingTest/mexRejectsInvalidInput (identifiers; the oracle cases compare them with
    MATLAB's)"""
    tim = from_matlab(_tim4())
    for _, b, m in _TILE_ERRORS:
        with pytest.raises(sk.StitchKernelError) as e:
            tile_mex(tim, np.array(b, float) if b else None, m)
        assert e.value.identifier == "tile_blending:inputError"
    with pytest.raises(sk.StitchKernelTypeError) as e:
        sk.tile_blend(tim.astype(np.float64), None, "none")
    assert e.value.identifier == "tile_blending:dataTypeError"
    with pytest.raises(ValueError, match="Unsupported blending method: None."):
        sk.tile_blending_3d(tim, None, "None")


def test_tile_boxes_round_trip():
    _, mboxes = random_batch((23, 19, 11), 9, "uint16", 300)
    rows = sk.tile_boxes_from_matlab(mboxes)
    np.testing.assert_array_equal(sk.tile_boxes_to_matlab(rows), mboxes)
    # a row is the tile's index and its box
    t, *box = rows[0]
    assert int(t) == int(mboxes[0, 0]) - 1
    assert tuple(int(v) for v in box) == box_from_matlab(mboxes[0, 1:])


# --- StitchHelperMexTest (the kernel parts) -----------------------------------------------------


@pytest.mark.parametrize("cls", ["uint16", "single", "uint8", "double"])
def test_major_tiles_with_data_everywhere_are_valid(cls):
    f = from_matlab(_major_full(cls))
    assert sk.check_major_tile_valid(f, np.array([1, 0, 1, 0], bool)) == (True, None)


@pytest.mark.parametrize("cls", ["uint16", "single", "uint8", "double"])
def test_voxel_without_data_in_the_major_tiles_is_found(cls):
    f = from_matlab(_major_hole(cls))
    valid, box = sk.check_major_tile_valid(f, np.array([1, 0, 1], bool))
    assert valid is False
    assert box == (1, 2, 5, 2, 3, 6)  # MATLAB's [3, 6, 2, 3, 6, 2]
    f = from_matlab(_major_hole(cls, True))
    assert sk.check_major_tile_valid(f, np.array([1, 0, 1], bool))[0] is True


@pytest.mark.parametrize("num_threads", [1, 2, 4])
def test_bounds_of_voxels_without_data_with_many_threads(num_threads):
    f = from_matlab(_major_bounds())
    for _ in range(5):
        valid, box = sk.check_major_tile_valid(f, np.ones(2, bool), num_threads=num_threads)
        assert not valid
        assert box == (0, 0, 0, 64, 256, 256)


def test_integral_image_mex_errors_have_identifiers():
    """StitchHelperMexTest/integralImageMexErrorsHaveIdentifiers and RejectsMoreThan3D"""
    a = from_matlab(fvol((10, 12, 8), np.float32, 130))
    for t, ident, a_ in [
        ((3, 3), "integral_image_3d:inputError", a),
        ((3, 3, 0), "integral_image_3d:inputError", a),
        ((3, 3, 3), "integral_image_3d:dataTypeError", a.astype(np.float64)),
        (
            (3, 3, 3),
            "integral_image_3d:inputError",
            from_matlab(fvol((10, 12, 8, 2), "single", 131)),
        ),
    ]:
        with pytest.raises(ERRORS) as e:
            sk.integral_image(a_, t)
        assert e.value.identifier == ident
    # other classes run MATLAB's code (in their class)
    out = sk.integral_image_3d(a.astype(np.float64), (2, 3, 4))
    assert out.dtype == np.float64
    np.testing.assert_allclose(out, sk.integral_image_3d(a, (2, 3, 4)), atol=1e-4)


# --- definitions, independent of MATLAB ---------------------------------------------------------


def test_integral_image_is_the_window_sum():
    """integer-valued inputs: every window sum is exact in single"""
    a = from_matlab(fvol((9, 11, 7), np.uint16, 150, 0, 50))  # (7, 9, 11)
    t = (2, 3, 5)
    out = sk.integral_image(a, t)
    assert out.shape == tuple(a.shape[i] + t[i] - 1 for i in range(3))
    pad = np.pad(a.astype(np.float64), [(n - 1, n - 1) for n in t])
    ref = np.zeros(out.shape)
    for i, j, k in np.ndindex(*out.shape):
        ref[i, j, k] = pad[i : i + t[0], j : j + t[1], k : k + t[2]].sum()
    np.testing.assert_array_equal(out, ref.astype(np.float32))


def test_feather_blending_is_the_weighted_mean():
    mf, md = feather_tiles((7, 6, 5), 3, "uint16", 160)
    f, d = from_matlab(mf), from_matlab(md)
    out = sk.feather_blend(f, d)
    ref = np.zeros(out.shape, dtype=np.uint16)
    for idx in np.ndindex(*out.shape):
        s, n = np.float32(0), np.float32(0)
        for t in range(3):
            if f[(t, *idx)] != 0:
                s = np.float32(s + d[(t, *idx)])
                n = np.float32(n + np.float32(f[(t, *idx)]) * d[(t, *idx)])
        if s != 0:
            ref[idx] = int(float(np.float32(n / s)) + 0.5)
    np.testing.assert_array_equal(out, ref)


def test_distance_map_resize_keeps_the_corners():
    dmat = from_matlab(fvol((5, 4, 3), np.float32, 170, 0.5, 1.5))  # (3, 5, 4)
    out = sk.resize_distance_map(dmat, (7, 9, 10), 1)
    for c in np.ndindex(2, 2, 2):
        src = tuple((dmat.shape[i] - 1) * c[i] for i in range(3))
        dst = tuple((out.shape[i] - 1) * c[i] for i in range(3))
        assert out[dst] == dmat[src]
    same = sk.resize_distance_map(dmat, dmat.shape, 2)
    np.testing.assert_array_equal(same, dmat * dmat)
    # the box form: the box's size
    box = (1, 2, 3, 8, 11, 13)
    np.testing.assert_array_equal(sk.feather_distance_map_resize_3d(dmat, box, 1), out)


def test_replace_nan_in_place_and_layouts():
    x = from_matlab(_with_nans(np.float32, 180))
    ref = np.where(np.isnan(x), np.float32(0), x)
    out = sk.replace_nan_with_value(x, 0)
    assert out is not x and out.flags.c_contiguous
    np.testing.assert_array_equal(out, ref)
    m = to_matlab(x)  # element-wise: any contiguous layout
    np.testing.assert_array_equal(sk.replace_nan_with_value(m, 0), to_matlab(ref))
    y = x.copy()
    assert sk.replace_nan_with_value(y, 0, inplace=True) is y
    np.testing.assert_array_equal(y, ref)


def test_arrays_are_never_copied():
    """MATLAB-layout tiles raise instead of being copied; the indexed blend writes into the
    caller's array"""
    mf, md = feather_tiles((6, 5, 4), 2, "uint16", 190)
    f, d = from_matlab(mf), from_matlab(md)  # (2, 4, 6, 5)
    with pytest.raises(ValueError, match="C-contiguous"):
        sk.feather_blend(np.asfortranarray(f), d)
    with pytest.raises(ValueError, match="C-contiguous"):
        sk.tile_blend(np.asfortranarray(f), None, "max")
    out = np.zeros((6, 8, 7), dtype=np.uint16)
    res, mex = sk.feather_blending_3d(f, d, out, (1, 1, 1, 5, 7, 6))
    assert res is out and mex is True
    np.testing.assert_array_equal(out[1:5, 1:7, 1:6], sk.feather_blend(f, d))
    ro = out.copy()
    ro.flags.writeable = False
    with pytest.raises(ValueError, match="writeable"):
        sk.feather_blend_into(f, d, ro, (1, 1, 1, 5, 7, 6))


# --- FeatherBlendingTest ------------------------------------------------------------------------


@pytest.mark.parametrize(
    "tile_class, weight_class",
    [("uint8", "single"), ("double", "single"), ("uint16", "double"), ("single", "double")],
)
def test_overlap_is_blended_into_the_block(tile_class, weight_class):
    """FeatherBlendingTest/overlapIsBlendedIntoTheBlock: MATLAB's code for these classes writes
    its blend (single), cast to the block's class, into the block at the box"""
    from petakit5d.core.dtypes import cast, single, times

    rng = CounterRNG(400, f"overlap/{tile_class}/{weight_class}")
    shape = (2, 5, 6, 7)
    f = rng.integers(1, 200, shape) * (rng.uniform(shape) > 0.2)
    f[0] = np.maximum(f[0], 1)  # tile 1 has data everywhere
    f = f.astype(tile_class)
    d = (rng.uniform(shape) * 0.9 + 0.1).astype(weight_class)
    block = rng.integers(1, 200, (8, 10, 12)).astype(tile_class)
    box = (1, 2, 3, 6, 8, 10)
    w = d * (f != 0)
    p = times(single(f), w)  # single, as MATLAB's single .* double
    blend = (p[0] + p[1]) / (w[0] + w[1])
    if blend.dtype != np.float32:  # single ./ double: in double, rounded
        blend = blend.astype(np.float32)
    expected = block.copy()
    expected[1:6, 2:8, 3:10] = cast(blend, tile_class)
    out, mex_compute = sk.feather_blending_3d(f, d, block, box)
    assert out is block and not mex_compute and out.dtype == np.dtype(tile_class)
    np.testing.assert_array_equal(out, expected)
