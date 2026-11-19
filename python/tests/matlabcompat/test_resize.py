"""Tests of petakit5d.matlabcompat.resize against MATLAB R2024b: imresize3 and imresize themselves
(uint8, uint16, single, double and the other integer classes; linear, cubic, nearest and the box
and triangle kernels; up and down per axis, non-integer factors, sizes 1 to 600; the scale, size
and NaN forms; antialiasing on and off; Inf, NaN and -0), the internals imresize3_fast calls
(contributions, imresizemex, the cubic and triangle kernels; black box), the lab's
imresize3_fast and imresize3_average, and the port's call sites (tiffToZarr, XR_resampleFrame,
resampleZarrBlock, MIP_block, stitch_process_zarr_tile, XR_deskewRotateFrame, the puncta
removal, feather_distance_map_resize_3d, the multi-resolution registration).

Oracle cases call the built-ins themselves (black box; inputs from our bit-stable generator in
MATLAB's form, roles map them). Exact unless a case names its tolerance (none does).

MATLAB tests: ResampleTest's imresize3FastMatchesImresize3 and
imresize3FastRegionMatchesCroppedImresize3 (36 cases each; the class is T04's) are the
``resample_test.*`` cases: the port's imresize3_fast against MATLAB's own imresize3 on the same
inputs, factors, methods and output regions (whole and four regions, in the class and in
double).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d.core import dtypes
from petakit5d.matlabcompat import resize as R
from petakit5d.testing import oracle
from petakit5d.testing.cases import Box, Case, In, PerAxis, Value, Vol
from petakit5d.testing.synth import CounterRNG

_ext = pytest.importorskip("petakit5d._ext.resize", reason="the resize kernels are not built")

# ---------------------------------------------------------------------------------------------
# inputs in MATLAB's form
# ---------------------------------------------------------------------------------------------

_SCALE = {"uint8": 250.0, "uint16": 60000.0, "single": 1000.0, "double": 1000.0}


def volume(seed: int, msize, cls: str) -> np.ndarray:
    """uniform data of MATLAB size ``msize`` in class ``cls`` (doubles also negative)"""
    u = CounterRNG(seed, "n01").uniform(tuple(int(s) for s in msize))
    if cls == "double":
        a = u * 1000 - 300
    elif cls == "single":
        a = (u * 1000).astype(np.float32)
    elif cls in ("uint8", "uint16"):
        a = dtypes.cast(u * _SCALE[cls], cls)
    elif cls == "int8":
        a = dtypes.cast(u * 250 - 125, cls)
    elif cls == "int16":
        a = dtypes.cast(u * 60000 - 30000, cls)
    elif cls == "int32":
        a = dtypes.cast(u * 4e9 - 2e9, cls)
    elif cls == "uint32":
        a = dtypes.cast(u * 4e9, cls)
    else:
        raise ValueError(cls)
    return np.asfortranarray(np.asarray(a))


def special(seed: int, msize, cls: str) -> np.ndarray:
    """uniform data with Inf, -Inf, NaN (MATLAB's) and a plane of -0 (float classes)"""
    a = volume(seed, msize, cls).copy(order="F")
    r = CounterRNG(seed, "n01-special")
    flat = a.reshape(-1, order="F")
    picks = (r.uniform((6,)) * flat.size).astype(np.int64)
    vals = [np.inf, -np.inf, dtypes.NAN, np.inf, dtypes.NAN, -np.inf]
    for p, v in zip(picks, vals, strict=True):
        flat[p] = v
    if a.ndim == 3:
        a[:, :, -1] = -0.0
    else:
        a[:, -1] = -0.0
    return a


def negzero(msize, cls: str, spot=(1, 1, 1)) -> np.ndarray:
    """-0 everywhere but a 3 at ``spot`` (0-based MATLAB subscripts)"""
    a = np.full(tuple(msize), -0.0, dtype=np.float32 if cls == "single" else np.float64, order="F")
    a[tuple(spot[: a.ndim])] = 3
    return a


# NaN of three payloads (MATLAB's, a positive one, another negative one), Inf and -Inf
PAYLOADS = {
    "single": (0xFFC00000, 0x7FC00001, 0xFFC00005),
    "double": (0xFFF8000000000000, 0x7FF8000000000001, 0xFFF8000000000005),
}


def payloads(seed: int, msize, cls: str, frac: float = 0.06) -> np.ndarray:
    """uniform data with NaN of the three payloads of PAYLOADS, Inf and -Inf (float classes)"""
    a = volume(seed, msize, cls).copy(order="F")
    u = CounterRNG(seed, "n01-payloads").uniform((5, *(int(m) for m in msize)))
    bits = a.view(np.uint32 if cls == "single" else np.uint64)
    for k, b in enumerate(PAYLOADS[cls]):
        bits[u[k] < frac] = b
    a[u[3] < 0.015] = np.inf
    a[u[4] < 0.015] = -np.inf
    return a


def tap_placements(cls: str) -> np.ndarray:
    """every placement of {1.5, two NaN payloads, Inf} over 4 taps: a 4 x 256 matrix"""
    ut = np.uint32 if cls == "single" else np.uint64
    ft = np.float32 if cls == "single" else np.float64
    alpha = np.array(
        [
            np.asarray(1.5, ft).view(ut),
            PAYLOADS[cls][0],
            PAYLOADS[cls][1],
            np.asarray(np.inf, ft).view(ut),
        ],
        dtype=ut,
    )
    idx = np.array(np.meshgrid(*[np.arange(4)] * 4, indexing="ij")).reshape(4, -1, order="F")
    return np.asfortranarray(alpha[idx].view(ft))


def mround(v: float) -> int:
    return int(dtypes.mround(np.float64(v)))


# ---------------------------------------------------------------------------------------------
# the Python sides: MATLAB's name-value pairs in, Python's forms through the roles
# ---------------------------------------------------------------------------------------------

_NV = {
    "antialiasing": "antialiasing",
    "method": "method",
    "scale": "scale",
    "outputsize": "output_size",
}


def _nv(pairs) -> dict:
    kw = {}
    for k in range(0, len(pairs), 2):
        name = str(pairs[k]).lower()
        key = [v for n, v in _NV.items() if n.startswith(name)][0]
        kw[key] = pairs[k + 1]
    return kw


def _split(args) -> tuple[list, dict]:
    """MATLAB's arguments after the image: a scale or size, a method, then name-value pairs"""
    args = list(args)
    pos = []
    if args and not isinstance(args[0], str):
        pos.append(args.pop(0))
    if len(args) % 2 == 1:
        pos.append(args.pop(0))
    return pos, _nv(args)


def py_imresize3(v, *args):
    pos, kw = _split(args)
    return R.imresize3(v, *pos, **kw)


def py_imresize(a, *args):
    pos, kw = _split(args)
    return R.imresize(a, *pos, **kw)


def py_contributions(n, m, scale, kernel, width, aa):
    w, ind = R.contributions(
        int(n),
        int(m),
        scale,
        {"cubic": "cubic", "triangle": "triangle", "box": "box"}[kernel],
        width,
        bool(aa),
    )
    return w, (ind + 1).astype(np.float64)


def py_kernel(kernel, x):
    return R.kernel_values(kernel, x)


def py_imresizemex(a, wt, indt, dim):
    """imresizemex(A, weights', indices', dim) of a matrix: dim 1 is axis -2, 2 axis -1"""
    w = np.ascontiguousarray(np.asarray(wt, dtype=np.float64).T)
    ind = np.ascontiguousarray(np.asarray(indt, dtype=np.float64).T).astype(np.int64) - 1
    return R.resize_axis(a, -2 if int(dim) == 1 else -1, w, ind)


def py_imresizemex_cubic(a, dim, m):
    """imresizemex of a matrix with the cubic weights of contributions (antialiasing)"""
    axis = -2 if int(dim) == 1 else -1
    n = a.shape[axis]
    w, ind = R.contributions(n, int(m), m / n, "cubic", None, True)
    return R.resize_axis(a, axis, w, ind)


def py_fast(v, size, method="cubic", box=None, cls=None):
    return R.imresize3_fast(v, size, method, box, cls)


def py_fast_pair(v, size, method, box):
    return R.imresize3_fast(v, size, method, box), R.imresize3_fast(v, size, method, box, "double")


def py_average(v, rs):
    return R.imresize3_average(v, rs)


# ---------------------------------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------------------------------

CASES: list[Case] = []


def add(case: Case) -> None:
    CASES.append(case)


def resize3_case(cid: str, make, msize_out, method, *extra, seed=None):
    add(
        Case(
            cid,
            matlab="imresize3",
            args=[Vol(In("v")), PerAxis(list(msize_out), "yxz"), method, *extra],
            python=py_imresize3,
            out_roles=("vol",),
            inputs={"v": Value(make)},
        )
    )


def resize2_case(cid: str, make, size_out, method, *extra):
    add(
        Case(
            cid,
            matlab="imresize",
            args=[Vol(In("a")), list(size_out), method, *extra],
            python=py_imresize,
            out_roles=("vol",),
            inputs={"a": Value(make)},
        )
    )


CLASSES = ("uint8", "uint16", "single", "double")

# imresize3: sizes [y x z] (z >= 2: one plane is a MATLAB error), each with three output sizes
SIZES3 = [
    (7, 9, 5),
    (12, 10, 6),
    (16, 16, 8),
    (23, 17, 11),
    (5, 1, 7),
    (1, 6, 4),
    (3, 4, 2),
    (2, 2, 2),
    (40, 36, 12),
    (41, 37, 23),
    (600, 3, 2),
    (4, 600, 2),
    (3, 5, 600),
    (64, 48, 20),
    (1, 1, 9),
    (9, 13, 3),
    (100, 7, 5),
    (2, 50, 30),
    (37, 41, 19),
    (17, 23, 2),
    (6, 6, 6),
    (31, 2, 3),
    (11, 19, 13),
    (8, 3, 47),
]
FACTORS3 = [
    (0.5, 0.5, 0.5),
    (2, 2, 2),
    (1, 1, 2),
    (0.7, 1.6, 3),
    (1.5, 1.3, 2.5),
    (1, 1, 1),
    (3, 0.4, 1),
    (0.33, 1, 0.25),
    (1.1, 0.9, 1.7),
    (0.1, 0.2, 0.5),
    (2.5, 0.75, 0.6),
    (1, 2, 0.5),
]


def out_size(msize, f):
    return [max(1, mround(s * g)) for s, g in zip(msize, f, strict=True)]


for i, ms in enumerate(SIZES3):
    picks = [FACTORS3[(3 * i + k) % len(FACTORS3)] for k in range(3)]
    for j, f in enumerate(picks):
        osz = out_size(ms, f)
        if j == 2 and i % 4 == 0:
            osz[i % 3] = 600 if ms[i % 3] < 300 else 7  # up to 600 and far down
        for cls in CLASSES:
            for method in ("linear", "cubic", "nearest"):
                sid = "x".join(map(str, ms))
                oid = "x".join(map(str, osz))
                resize3_case(
                    f"v3.{cls}.{method}.{sid}.to{oid}.{j}",
                    lambda ms=ms, cls=cls, i=i: volume(100 + i, ms, cls),
                    osz,
                    method,
                )

# imresize3: other forms, options, kernels, classes and values
_A3 = (11, 9, 7)
for cls in CLASSES:
    mk = lambda cls=cls: volume(7, _A3, cls)  # noqa: E731
    for k, (form, extra) in enumerate(
        [
            ([0.5], ()),
            ([1.7], ()),
            ([0.37], ("linear",)),
            ([2], ("nearest",)),
            ([0.3], ("nearest",)),
        ]
    ):
        add(
            Case(
                f"v3form.{cls}.scale{k}",
                matlab="imresize3",
                args=[Vol(In("v")), form[0], *extra],
                python=py_imresize3,
                out_roles=("vol",),
                inputs={"v": Value(mk)},
            )
        )
    for k, args in enumerate(
        [
            ["Scale", PerAxis([0.5, 0.6, 0.7], "yxz"), "Method", "linear"],
            ["Scale", PerAxis([1.5, 0.4, 2.2], "yxz")],
            ["OutputSize", PerAxis([5, 4, 3], "yxz"), "Method", "nearest"],
            [PerAxis([6.2, 4.5, 3.1], "yxz")],
            [PerAxis([6.7, 4.5, 3.5], "yxz"), "linear"],
            [PerAxis([math.nan, math.nan, 3], "yxz")],
            [PerAxis([math.nan, 9, math.nan], "yxz"), "linear"],
            [PerAxis([7, math.nan, math.nan], "yxz"), "nearest"],
            [PerAxis([5, 4, 3], "yxz"), "cubic", "Antialiasing", False],
            [PerAxis([5, 4, 3], "yxz"), "linear", "Antialiasing", False],
            [PerAxis([5, 20, 3], "yxz"), "nearest", "Antialiasing", True],
            [PerAxis([22, 4, 13], "yxz"), "nearest", "Antialiasing", True],
            [PerAxis([5, 4, 3], "yxz"), "box"],
            [PerAxis([22, 4, 13], "yxz"), "box"],
            [PerAxis([5, 4, 3], "yxz"), "triangle"],
            [PerAxis([11, 9, 7], "yxz")],
            [PerAxis([11, 9, 7], "yxz"), "nearest"],
            [PerAxis([5, 4, 3], "yxz"), "Lin"],
            [PerAxis([5, 4, 3], "yxz"), "linear", "Antialiasing", 0],
            [PerAxis([5, 4, 3], "yxz"), "Meth", "linear", "anti", False],
            [0.99],
            [0.99, "linear"],
        ]
    ):
        add(
            Case(
                f"v3form.{cls}.opt{k}",
                matlab="imresize3",
                args=[Vol(In("v")), *args],
                python=py_imresize3,
                out_roles=("vol",),
                inputs={"v": Value(mk)},
            )
        )
for cls in ("int8", "int16", "int32", "uint32"):
    for method in ("linear", "cubic", "nearest"):
        resize3_case(
            f"v3cls.{cls}.{method}", lambda cls=cls: volume(9, (9, 8, 4), cls), (5, 13, 3), method
        )
for cls in ("single", "double"):
    for method in ("linear", "cubic", "nearest"):
        for k, (ms, osz) in enumerate(
            [((9, 8, 6), (5, 13, 4)), ((9, 8, 6), (17, 4, 11)), ((6, 7, 5), (6, 7, 9))]
        ):
            resize3_case(
                f"v3special.{cls}.{method}.{k}",
                lambda cls=cls, ms=ms, k=k: special(20 + k, ms, cls),
                osz,
                method,
            )
    # -0 through passes that copy (scale 1, one-hot weights) and that add
    for k, (ms, osz, method, extra) in enumerate(
        [
            ((6, 5, 4), (6, 5, 4), "linear", ()),
            ((6, 5, 4), (6, 5, 8), "cubic", ()),
            ((6, 5, 4), (12, 10, 8), "nearest", ("Antialiasing", True)),
            ((6, 5, 4), (3, 2, 2), "nearest", ("Antialiasing", True)),
            ((6, 5, 4), (12, 10, 8), "box", ()),
            ((12, 9, 6), (4, 3, 2), "linear", ("Antialiasing", False)),
            ((12, 9, 6), (4, 3, 2), "cubic", ("Antialiasing", False)),
            ((12, 9, 6), (4, 18, 6), "linear", ("Antialiasing", False)),
        ]
    ):
        resize3_case(
            f"v3negzero.{cls}.{k}",
            lambda cls=cls, ms=ms: negzero(ms, cls),
            osz,
            method,
            *extra,
        )
    # Inf at taps of weight 0 (enlarging by 3: rows with a one-hot cubic) and at one-hot passes
    for k, (osz, method, extra) in enumerate(
        [
            ((12, 3, 3), "cubic", ()),
            ((8, 6, 6), "linear", ()),
            ((4, 1, 3), "cubic", ("Antialiasing", False)),
        ]
    ):

        def mk_inf(cls=cls, k=k):
            a = np.ones(
                (12 if k == 2 else 4, 3, 3),
                dtype=np.float32 if cls == "single" else np.float64,
                order="F",
            )
            a[0, 1, 1] = np.inf
            a[2, 0, 2] = dtypes.NAN
            return a

        resize3_case(f"v3inf.{cls}.{k}", mk_inf, osz, method, *extra)

# imresize3: errors
for k, (make, args, ident) in enumerate(
    [
        (lambda: volume(1, (8, 6), "single"), [PerAxis([4, 3, 2], "yxz")], "MATLAB:badsubscript"),
        (
            lambda: volume(1, (8, 6, 1), "single"),
            [PerAxis([4, 3, 2], "yxz"), "linear"],
            "MATLAB:badsubscript",
        ),
        (
            lambda: volume(1, (6, 5, 3, 2), "single"),
            [PerAxis([3, 4, 2], "yxz")],
            "images:imresize3:incorrectDimensions",
        ),
        (
            lambda: np.zeros((0, 3, 4), order="F"),
            [PerAxis([2, 2, 2], "yxz")],
            "MATLAB:imresize3:expectedNonempty",
        ),
        (
            lambda: volume(1, (6, 5, 3), "single"),
            [PerAxis([0, 4, 3], "yxz")],
            "MATLAB:imresize3:expectedPositive",
        ),
        (
            lambda: volume(1, (6, 5, 3), "single"),
            [PerAxis([math.nan, 4, 3], "yxz")],
            "images:imresize3:invalidOutputSize",
        ),
        (
            lambda: volume(1, (6, 5, 3), "single"),
            [PerAxis([5, 4, 3], "yxz"), "bilinear"],
            "images:imresize3:oddNumberArgs",
        ),
    ]
):
    add(
        Case(
            f"v3error.{k}",
            matlab="imresize3",
            args=[Vol(In("v")), *args],
            python=py_imresize3,
            out_roles=("vol",),
            inputs={"v": Value(make)},
            expect_error=ident,
            py_error=(R.MatlabError,),
        )
    )

# imresize: images [y x], two output sizes each
SIZES2 = [
    (7, 9),
    (12, 10),
    (16, 16),
    (23, 17),
    (5, 1),
    (1, 6),
    (3, 4),
    (2, 2),
    (40, 36),
    (41, 37),
    (600, 3),
    (4, 600),
    (64, 48),
    (100, 7),
    (2, 50),
    (37, 41),
    (17, 23),
    (1, 1),
]
FACTORS2 = [
    (0.5, 0.5),
    (2, 2),
    (1, 2),
    (0.7, 1.6),
    (1.5, 1.3),
    (1, 1),
    (3, 0.4),
    (0.33, 1),
    (0.1, 0.25),
]
for i, ms in enumerate(SIZES2):
    for j in range(2):
        f = FACTORS2[(2 * i + j) % len(FACTORS2)]
        osz = out_size(ms, f)
        if j == 1 and i % 5 == 0:
            osz[i % 2] = 600 if ms[i % 2] < 300 else 9
        for cls in CLASSES:
            for method in ("bilinear", "bicubic", "nearest"):
                sid = "x".join(map(str, ms))
                oid = "x".join(map(str, osz))
                resize2_case(
                    f"v2.{cls}.{method}.{sid}.to{oid}.{j}",
                    lambda ms=ms, cls=cls, i=i: volume(300 + i, ms, cls),
                    osz,
                    method,
                )
_A2 = (13, 10)
for cls in CLASSES:
    mk2 = lambda cls=cls: volume(8, _A2, cls)  # noqa: E731
    for k, args in enumerate(
        [
            [0.5],
            [1.7],
            [0.37, "bilinear"],
            [2, "nearest"],
            [[math.nan, 4]],
            [[6, math.nan], "bilinear"],
            [[6, 4], "bilinear", "Antialiasing", False],
            [[6, 4], "nearest", "Antialiasing", True],
            [[26, 4], "nearest", "Antialiasing", True],
            [[6, 4], "Method", "bilinear"],
            ["Scale", 0.5],
            ["Scale", [0.5, 0.3], "Method", "bilinear"],
            ["OutputSize", [6, 4]],
            [[13, 10]],
            [[13, 10], "nearest"],
            [[6.2, 4.5]],
            [[6, 4], "cubic"],
            [[6, 4], "triangle"],
            [[6, 4], "box"],
            [[3, 25], "bil"],
            [[6, 4], "Bicubic"],
            [0.99, "bilinear"],
        ]
    ):
        add(
            Case(
                f"v2form.{cls}.opt{k}",
                matlab="imresize",
                args=[Vol(In("a")), *args],
                python=py_imresize,
                out_roles=("vol",),
                inputs={"a": Value(mk2)},
            )
        )
    # stacks: every plane alike (MATLAB's dimensions 3 and 4)
    for k, (ms, osz, method) in enumerate(
        [((6, 5, 3), (3, 9), "bilinear"), ((6, 5, 3, 2), (3, 9), "bicubic")]
    ):
        resize2_case(
            f"v2stack.{cls}.{k}", lambda ms=ms, cls=cls, k=k: volume(40 + k, ms, cls), osz, method
        )
for cls in ("int8", "int16", "int32", "uint32"):
    for method in ("bilinear", "bicubic", "nearest"):
        resize2_case(
            f"v2cls.{cls}.{method}", lambda cls=cls: volume(10, (9, 8), cls), (5, 13), method
        )
for cls in ("single", "double"):
    for method in ("bilinear", "bicubic", "nearest"):
        for k, osz in enumerate([(5, 13), (17, 4)]):
            resize2_case(
                f"v2special.{cls}.{method}.{k}",
                lambda cls=cls, k=k: special(60 + k, (9, 8), cls),
                osz,
                method,
            )
    for k, (ms, osz, method, extra) in enumerate(
        [
            ((6, 5), (6, 5), "bilinear", ()),
            ((6, 5), (6, 9), "bicubic", ()),
            ((6, 5), (12, 10), "nearest", ("Antialiasing", True)),
            ((12, 9), (4, 3), "bilinear", ("Antialiasing", False)),
        ]
    ):
        resize2_case(
            f"v2negzero.{cls}.{k}", lambda cls=cls, ms=ms: negzero(ms, cls), osz, method, *extra
        )
for k, (make, args, ident) in enumerate(
    [
        (
            lambda: np.zeros((0, 3), order="F"),
            [[2, 2]],
            "MATLAB:resizeParseInputs:expectedNonempty",
        ),
        (
            lambda: volume(1, (6, 5), "single"),
            [[0, 4]],
            "MATLAB:resizeParseInputs:expectedPositive",
        ),
        (
            lambda: volume(1, (6, 5), "single"),
            [[6, 4, 2]],
            "MATLAB:images:imresize:invalidOutputSize",
        ),
        (
            lambda: volume(1, (6, 5), "single"),
            [[math.nan, math.nan]],
            "MATLAB:images:imresize:allNaN",
        ),
        (
            lambda: volume(1, (6, 5), "single"),
            [[5, 3], "linear"],
            "MATLAB:images:imresize:oddNumberArgs",
        ),
    ]
):
    add(
        Case(
            f"v2error.{k}",
            matlab="imresize",
            args=[Vol(In("a")), *args],
            python=py_imresize,
            out_roles=("vol",),
            inputs={"a": Value(make)},
            expect_error=ident,
            py_error=(R.MatlabError,),
        )
    )

# the internals imresize3_fast calls (black box): contributions, the kernels, imresizemex
_CONTRIB = (
    "[w, ind] = matlab.images.internal.resize.contributions(in{1}, in{2}, in{3}, "
    "str2func(['matlab.images.internal.resize.' in{4}]), in{5}, in{6}); out = {w, ind};"
)
for n, m in [
    (1, 1),
    (1, 5),
    (2, 1),
    (2, 3),
    (5, 7),
    (12, 5),
    (23, 47),
    (41, 14),
    (100, 33),
    (600, 13),
    (13, 600),
    (17, 17),
    (9, 4),
    (40, 20),
    (36, 18),
    (512, 256),
    (500, 1000),
]:
    for kernel, width in (("cubic", 4), ("triangle", 2), ("box", 1)):
        for aa in (True, False):
            for k, s in enumerate((m / n, m / n * 1.001)):
                add(
                    Case(
                        f"contrib.{kernel}.{n}to{m}.aa{int(aa)}.s{k}",
                        code=_CONTRIB,
                        args=[n, m, s, kernel, width, aa],
                        nargout=2,
                        python=py_contributions,
                    )
                )
_GRID = (
    np.arange(-3000, 3001) / 1000.0,
    CounterRNG(5, "n01-x").uniform((1, 4000)).reshape(-1) * 6 - 3,
)
for kernel in ("cubic", "triangle", "box"):
    add(
        Case(
            f"kernel.{kernel}",
            code=f"out = {{double(matlab.images.internal.resize.{kernel}(in{{2}}))}};",
            args=[kernel, In("x")],
            python=py_kernel,
            inputs={"x": Value(lambda: np.concatenate(_GRID).reshape(1, -1))},
        )
    )
_MEX = "out = {matlab.images.internal.resize.imresizemex(in{1}, in{2}, in{3}, in{4})};"
_MEX_CUBIC = (
    "n = size(in{1}, in{2}); [w, ind] = matlab.images.internal.resize.contributions(n, in{3}, "
    "in{3} / n, @matlab.images.internal.resize.cubic, 4, true); "
    "out = {matlab.images.internal.resize.imresizemex(in{1}, w', ind', in{2})};"
)
_WQ = np.array(
    [
        [0.1, 0.2, 0.3, 0.4],
        [1 / 3, 1 / 3, 1 / 3, 0],
        [0.5, 0.5, 0, 0],
        [-0.25, 0.75, 0.75, -0.25],
        [1e-8, 1, -1, 1e-8],
        [0.49999999999999994, 0, 0, 0],
        [2.4999999999999996, 0, 0, 0],
        [1.5, 0, 0, 0],
        [-0.49999999999999994, 0, 0, 0],
    ]
)
for cls in ("uint8", "uint16", "int16", "int8", "int32", "uint32", "single", "double"):

    def mk_q(cls=cls):
        a = volume(70, (4, 60), cls)
        if cls in ("single", "double"):
            a[:, -4:] = np.array(
                [[1, 2, 3, 4], [np.inf, 1, 1, 1], [dtypes.NAN, 1, 1, 1], [-0.0, -0.0, -0.0, -0.0]]
            ).T
        else:
            a[:, -2:] = np.array([[1, 2, 3, 4], [1, 1, 1, 1]]).T
        return np.asfortranarray(a)

    add(
        Case(
            f"mex.{cls}.chosen",
            code=_MEX,
            args=[
                Vol(In("a")),
                _WQ.T.copy(),
                np.tile(np.arange(1.0, 5.0)[:, None], (1, _WQ.shape[0])),
                1,
            ],
            python=py_imresizemex,
            out_roles=("vol",),
            inputs={"a": Value(mk_q)},
        )
    )
    for dim in (1, 2):
        add(
            Case(
                f"mex.{cls}.cubic.dim{dim}",
                code=_MEX_CUBIC,
                args=[Vol(In("a")), dim, 10],
                python=py_imresizemex_cubic,
                out_roles=("vol",),
                inputs={
                    "a": Value(
                        lambda cls=cls, dim=dim: volume(71, (23, 9) if dim == 1 else (9, 23), cls)
                    )
                },
            )
        )

# ResampleTest's imresize3FastMatchesImresize3 and imresize3FastRegionMatchesCroppedImresize3:
# the port's imresize3_fast against MATLAB's own imresize3
_FAST_FACTORS = {
    "zOnly": (1, 1, 2),
    "xyOnly": (2, 2, 1),
    "yOnly": (3, 1, 1),
    "zUp": (1, 1, 0.5),
    "all": (2, 2, 2),
    "nonInteger": (1.5, 1.3, 2.5),
    "zBetween": (2, 4, 3),
    "xzDownYUp": (0.7, 1.6, 3),
    "none": (1, 1, 1),
}
_FAST_REGION = (
    "e = imresize3(in{1}, in{2}, in{3}); d = imresize3(double(in{1}), in{2}, in{3}); b = in{4}; "
    "out = {e(b(1) : b(4), b(2) : b(5), b(3) : b(6)), d(b(1) : b(4), b(2) : b(5), b(3) : b(6))};"
)
_RT = (41, 37, 23)
for fname, f in _FAST_FACTORS.items():
    osz = [mround(s / g) for s, g in zip(_RT, f, strict=True)]
    for method in ("linear", "cubic"):
        for cls in ("single", "uint16"):
            mk_rt = lambda cls=cls: np.asfortranarray(  # noqa: E731
                np.asarray(dtypes.cast(CounterRNG(11, "n01-rt").uniform(_RT) * 60000, cls))
            )
            add(
                Case(
                    f"resample_test.fast.{fname}.{method}.{cls}",
                    matlab="imresize3",
                    args=[Vol(In("im")), PerAxis(osz, "yxz"), method],
                    python=py_fast,
                    out_roles=("vol",),
                    inputs={"im": Value(mk_rt)},
                )
            )
            boxes = [
                [1, 1, 1, *osz],
                [1, 2, 1, *[math.ceil(o / 2) for o in osz]],
                [*[math.floor(o / 3) + 1 for o in osz], *osz],
                [*[math.ceil(o / 2) for o in osz], *[math.ceil(o / 2) for o in osz]],
            ]
            for k, b in enumerate(boxes):
                add(
                    Case(
                        f"resample_test.region.{fname}.{method}.{cls}.b{k}",
                        code=_FAST_REGION,
                        args=[Vol(In("im")), PerAxis(osz, "yxz"), method, Box(b)],
                        nargout=2,
                        python=py_fast_pair,
                        out_roles=("vol", "vol"),
                        inputs={"im": Value(mk_rt)},
                    )
                )

# the lab's imresize3_fast itself (regions, classes, the fall-back to imresize3)
for k, (ms, osz, method, b, cls) in enumerate(
    [
        ((30, 26, 14), (45, 20, 30), "linear", None, None),
        ((30, 26, 14), (45, 20, 30), "cubic", [3, 2, 5, 40, 18, 22], None),
        ((30, 26, 14), (45, 20, 30), "cubic", [3, 2, 5, 40, 18, 22], "double"),
        ((30, 26, 14), (15, 26, 7), "linear", [1, 4, 2, 15, 20, 7], "single"),
        ((30, 26, 14), (45, 20, 30), "nearest", [3, 2, 5, 40, 18, 22], "double"),
        ((30, 26, 14), (45, 20, 30), "box", [3, 2, 5, 40, 18, 22], "double"),
        ((30, 26, 14), (12.5, 30.2, 14), "cubic", None, None),
        ((9, 13, 1), (18, 7, 1), "linear", None, None),
    ]
):
    args = [Vol(In("im")), PerAxis(list(osz), "yxz"), method]
    if b is not None or cls is not None:
        args.append(Box(b) if b is not None else None)
    if cls is not None:
        args.append(cls)
    add(
        Case(
            f"fast.{k}",
            matlab="imresize3_fast",
            args=args,
            python=py_fast,
            out_roles=("vol",),
            inputs={
                "im": Value(lambda ms=ms, k=k: volume(80 + k, ms, "uint16" if k % 2 else "single"))
            },
        )
    )

# imresize3_fast leaves the axes that keep their size alone (no mex pass): -0 stays
add(
    Case(
        "fast.negzero_same_size",
        matlab="imresize3_fast",
        args=[Vol(In("im")), PerAxis([6, 5, 4], "yxz"), "linear"],
        python=py_fast,
        out_roles=("vol",),
        inputs={"im": Value(lambda: negzero((6, 5, 4), "single"))},
    )
)

# NaN payloads: where NaN of different payloads meet in one output, MATLAB keeps the first in
# tap order (each pass), and imresize3_average the first of each block in MATLAB's order
_W_TAPS = np.array(
    [[0.25, 0.25, 0.25, 0.25], [-0.1, 0.6, 0.6, -0.1], [0, 1, 0, 0], [0.5, 0, 0, 0.5]]
)
for cls in ("single", "double"):
    for method in ("linear", "cubic", "nearest"):
        for k, (ms, osz) in enumerate(
            [((9, 8, 6), (5, 4, 3)), ((9, 8, 6), (13, 15, 11)), ((12, 7, 5), (9, 3, 12))]
        ):
            resize3_case(
                f"v3payload.{cls}.{method}.{k}",
                lambda cls=cls, ms=ms, k=k: payloads(110 + k, ms, cls),
                osz,
                method,
            )
    for method in ("bilinear", "bicubic"):
        for k, (ms, osz) in enumerate([((9, 8, 3), (5, 3)), ((9, 8), (17, 14))]):
            resize2_case(
                f"v2payload.{cls}.{method}.{k}",
                lambda cls=cls, ms=ms, k=k: payloads(120 + k, ms, cls),
                osz,
                method,
            )
    add(
        Case(
            f"mex.{cls}.payloads",
            code=_MEX,
            args=[Vol(In("a")), _W_TAPS.T.copy(), np.tile(np.arange(1.0, 5.0)[:, None], (1, 4)), 1],
            python=py_imresizemex,
            out_roles=("vol",),
            inputs={"a": Value(lambda cls=cls: tap_placements(cls))},
        )
    )
for method in ("linear", "cubic"):
    add(
        Case(
            f"fast.payloads.{method}",
            matlab="imresize3_fast",
            args=[Vol(In("im")), PerAxis([5, 13, 9], "yxz"), method],
            python=py_fast,
            out_roles=("vol",),
            inputs={"im": Value(lambda: payloads(130, (9, 8, 6), "single"))},
        )
    )
for k, (ms, rs, cls) in enumerate(
    [
        ((12, 10, 6), (2, 2, 2), "single"),
        ((12, 10, 6), (3, 5, 2), "single"),
        ((9, 8, 5), (4, 3, 2), "single"),
        ((13, 11, 7), (1, 1, 2), "double"),
        ((16, 12, 8), (4, 4, 4), "single"),
    ]
):
    add(
        Case(
            f"average.payloads.{k}",
            matlab="imresize3_average",
            args=[Vol(In("im")), PerAxis(list(rs), "yxz")],
            python=py_average,
            out_roles=("vol",),
            inputs={"im": Value(lambda ms=ms, cls=cls, k=k: payloads(140 + k, ms, cls, 0.04))},
        )
    )

# the lab's imresize3_average
for k, (ms, rs, cls) in enumerate(
    [
        ((12, 10, 6), (2, 2, 2), "uint16"),
        ((12, 10, 6), (3, 5, 2), "single"),
        ((40, 36, 12), (3, 5, 5), "uint16"),
        ((41, 37, 23), (2, 2, 2), "single"),
        ((13, 11, 7), (1, 1, 2), "double"),
        ((12, 10, 6), (12, 2, 2), "uint16"),
        ((12, 10, 6), (2, 2, 6), "single"),
        ((9, 8, 5), (4, 3, 2), "uint8"),
        ((100, 64, 24), (4, 4, 4), "uint16"),
        ((7, 9, 4), (7, 9, 4), "single"),
    ]
):
    add(
        Case(
            f"average.{k}",
            matlab="imresize3_average",
            args=[Vol(In("im")), PerAxis(list(rs), "yxz")],
            python=py_average,
            out_roles=("vol",),
            inputs={"im": Value(lambda ms=ms, cls=cls, k=k: volume(90 + k, ms, cls))},
        )
    )

# the port's call sites, on their own classes and shapes
_SITES = [
    # tiffToZarr, XR_resampleFrame, resampleZarrBlock: imresize3_fast(I, round(size ./ rs), method)
    (
        "site.tiffToZarr",
        "imresize3_fast",
        (64, 48, 20),
        "uint16",
        lambda s: [mround(v) for v in np.array(s) / (1, 1, 2)],
        (),
    ),
    (
        "site.resampleFrame.linear",
        "imresize3_fast",
        (40, 36, 12),
        "single",
        lambda s: [mround(v / 2) for v in s],
        ("linear",),
    ),
    (
        "site.resampleFrame.cubic",
        "imresize3_fast",
        (40, 36, 12),
        "single",
        lambda s: [mround(v / 1.5) for v in s],
        ("cubic",),
    ),
    (
        "site.resampleZarrBlock",
        "imresize3_fast",
        (66, 34, 18),
        "uint16",
        lambda s: [mround(v / f) for v, f in zip(s, (2, 2, 1), strict=True)],
        ("linear",),
    ),
    # MIP_block: imresize3(in_batch, round(size(in_batch) ./ dsfactor), 'linear')
    (
        "site.MIP_block",
        "imresize3",
        (64, 64, 32),
        "uint16",
        lambda s: [mround(v / f) for v, f in zip(s, (4, 4, 2), strict=True)],
        ("linear",),
    ),
    # stitch_process_zarr_tile: imresize3(frame, outSize)
    (
        "site.stitch_tile",
        "imresize3",
        (50, 40, 30),
        "uint16",
        lambda s: [mround(v / f) for v, f in zip(s, (1, 1, 2.5), strict=True)],
        (),
    ),
    # the puncta removal and feather_distance_map_resize_3d: enlarging, 'linear'
    ("site.puncta.v3", "imresize3", (4, 5, 3), "single", lambda s: [40, 50, 30], ("linear",)),
    ("site.feather.v3", "imresize3", (10, 8, 6), "double", lambda s: [53, 41, 17], ("linear",)),
    # multires registration: imresize3(region, 1 ./ downSample .* size(region), 'nearest')
    (
        "site.multires",
        "imresize3",
        (50, 47, 21),
        "uint16",
        lambda s: [v / 3 for v in s],
        ("nearest",),
    ),
    # the FSC analysis: imresize3(im, round(size(im) .* [xy xy dz] / psz))
    (
        "site.fsc",
        "imresize3",
        (40, 40, 30),
        "single",
        lambda s: [mround(v * f) for v, f in zip(s, (0.108, 0.108, 0.3), strict=True)],
        (),
    ),
]
for cid, fn, ms, cls, size_of, extra in _SITES:
    add(
        Case(
            cid,
            matlab=fn,
            args=[Vol(In("v")), PerAxis(size_of(ms), "yxz"), *extra],
            python=py_fast if fn == "imresize3_fast" else py_imresize3,
            out_roles=("vol",),
            inputs={"v": Value(lambda ms=ms, cls=cls, cid=cid: volume(len(cid) * 7, ms, cls))},
        )
    )
for cid, ms, cls, osz, extra in [
    # XR_deskewRotateFrame: imresize(frame, round(size(frame, [1, 2]) ./ pre_rs(1 : 2)))
    ("site.deskew.prers", (64, 40, 6), "uint16", (32, 20), ()),
    ("site.deskew.prers.single", (61, 37, 5), "single", (20, 12), ()),
    # XR_resampleFrame on 2-D images: imresize(im, outSize, 'Method', 'bilinear')
    ("site.resampleFrame.2d", (40, 36), "single", (20, 18), ("Method", "bilinear")),
    # the puncta removal, feather_distance_map_resize_3d: imresize(bg, size, 'bilinear')
    ("site.puncta.v2", (6, 5), "single", (60, 41), ("bilinear",)),
    ("site.feather.v2", (12, 9), "double", (50, 31), ("bilinear",)),
]:
    add(
        Case(
            cid,
            matlab="imresize",
            args=[Vol(In("a")), list(osz), *extra],
            python=py_imresize,
            out_roles=("vol",),
            inputs={"a": Value(lambda ms=ms, cls=cls, cid=cid: volume(len(cid) * 11, ms, cls))},
        )
    )


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# tests of our own
# ---------------------------------------------------------------------------------------------


def _ref_cubic(x):
    from petakit5d.matlabcompat import elementary

    a = np.abs(np.asarray(x, dtype=np.float64))
    a2 = a * a
    a3 = elementary.power(a, 3.0)
    return (1.5 * a3 - 2.5 * a2 + 1) * (a <= 1) + (-0.5 * a3 + 2.5 * a2 - 4 * a + 2) * (
        (1 < a) & (a <= 2)
    )


def test_cubic_kernel_is_the_formula_with_fdlibm_pow():
    x = np.concatenate(_GRID)
    assert np.array_equal(R.cubic(x).view(np.uint64), _ref_cubic(x).view(np.uint64))
    # |x|^3 is not |x|*|x|*|x| (MATLAB's power): the kernel differs from that in some values
    a = np.abs(x)
    naive = (1.5 * (a * a * a) - 2.5 * a * a + 1) * (a <= 1)
    naive += (-0.5 * (a * a * a) + 2.5 * a * a - 4 * a + 2) * ((1 < a) & (a <= 2))
    assert np.any(naive != R.cubic(x))


def test_fdlibm_pow_kernel_equals_elementary_power():
    from petakit5d.matlabcompat import elementary

    x = np.concatenate([np.abs(np.concatenate(_GRID)), [0.0, 1.0, 2.0, 1e-310, 7.5, 600.25]])
    for y in (3.0, 2.5, -1.5, 0.75):
        xs = x[x > 0] if y < 0 else x
        out = np.empty_like(xs)
        _ext.fdlibm_pow(xs, y, out)
        assert np.array_equal(
            out.view(np.uint64), np.asarray(elementary.power(xs, y)).view(np.uint64)
        ), y


def _rows(n_in, n_out, scale, kernel, aa):
    w, ind = R.contributions(n_in, n_out, scale, kernel, None, aa)
    return w, ind


def test_contributions_rows_sum_to_about_one_and_indices_mirror():
    for n, m in [(1, 3), (5, 7), (12, 5), (600, 13), (13, 600)]:
        for kernel in ("cubic", "triangle", "box"):
            w, ind = _rows(n, m, m / n, kernel, True)
            assert w.shape == ind.shape and w.shape[0] == m
            assert np.all(np.abs(w.sum(axis=1) - 1) < 1e-12)
            assert ind.min() >= 0 and ind.max() < n
            assert np.all(np.any(w != 0, axis=0))


def test_contributions_widen_the_kernel_only_when_shrinking_with_antialiasing():
    assert (
        _rows(40, 20, 0.5, "cubic", True)[0].shape[1]
        > _rows(40, 20, 0.5, "cubic", False)[0].shape[1]
    )
    assert _rows(20, 40, 2, "cubic", True)[0].shape == _rows(20, 40, 2, "cubic", False)[0].shape


@pytest.mark.parametrize("cls", ["uint8", "uint16", "int16", "single", "double"])
@pytest.mark.parametrize("method", ["linear", "cubic", "nearest"])
def test_results_do_not_depend_on_the_thread_count(cls, method):
    v = np.ascontiguousarray(np.transpose(volume(3, (37, 41, 19), cls), (2, 0, 1)))
    one = R.imresize3(v, (31, 20, 77), method, num_threads=1)
    for nt in (2, 3, 7):
        many = R.imresize3(v, (31, 20, 77), method, num_threads=nt)
        assert many.dtype == one.dtype and np.array_equal(many.view(np.uint8), one.view(np.uint8))


def test_baseline_and_avx2_kernels_agree():
    if _ext.best_isa() < 1:
        pytest.skip("no AVX2 on this CPU")
    a = np.ascontiguousarray(
        np.transpose(volume(4, (33, 29, 9), "single"), (2, 0, 1))
    )  # (9, 33, 29)
    wx, ix = R.contributions(29, 13, 13 / 29, "cubic", None, True)
    wy, iy = R.contributions(33, 50, 50 / 33, "cubic", None, True)
    for outer, n_in, inner, w, ind in ((9 * 33, 29, 1, wx, ix), (9, 33, 29, wy, iy)):
        outs = []
        for isa in (0, 1):
            dst = np.empty(outer * w.shape[0] * inner, dtype=np.float32)
            _ext.resize_single(a.reshape(-1), dst, outer, n_in, inner, w.shape[0], w, ind, 2, isa)
            outs.append(dst)
        assert np.array_equal(outs[0].view(np.uint32), outs[1].view(np.uint32))


def test_out_takes_the_result_and_may_be_reused():
    v = np.ascontiguousarray(np.transpose(volume(5, (20, 18, 10), "uint16"), (2, 0, 1)))
    ref = R.imresize3(v, (5, 9, 10), "cubic")
    out = np.empty_like(ref)
    for _ in range(2):
        r = R.imresize3(v, (5, 9, 10), "cubic", out=out)
        assert r is out and np.array_equal(out, ref)
    with pytest.raises(ValueError):
        R.imresize3(v, (5, 9, 10), "cubic", out=np.empty((5, 9, 9), np.uint16))
    same = R.imresize3(v, v.shape, "linear", out=np.empty_like(v))
    assert np.array_equal(same, v)
    assert R.imresize3(v, v.shape, "linear") is not v


def test_uint16_rounding_at_one_half_below_is_up():
    # MATLAB's resize rounds floor(s + 0.5): 0.49999999999999994 gives 1
    a = np.ones((3, 1), dtype=np.uint16)
    w = np.array([[0.49999999999999994], [2.4999999999999996], [0.5]])
    out = R.resize_axis(a, 0, w, np.zeros((3, 1), np.int64))
    assert out.ravel().tolist() == [1, 2, 1]
    neg = R.resize_axis(-np.ones((3, 1), dtype=np.int16), 0, w, np.zeros((3, 1), np.int64))
    assert neg.ravel().tolist() == [-1, -2, -1]


def test_one_hot_passes_copy_and_others_add():
    z = np.full((4, 6, 5), -0.0, dtype=np.float32)
    assert np.all(np.signbit(R.imresize3(z, (4, 6, 5), "linear")))
    assert np.all(np.signbit(R.imresize3(z, (8, 12, 10), "nearest")))
    assert not np.any(np.signbit(R.imresize3(z, (4, 6, 9), "linear")))
    # imresize3_fast always adds
    assert not np.any(np.signbit(R.imresize3_fast(z, (4, 6, 9), "linear")))


def test_imresize3_needs_three_dimensions():
    with pytest.raises(R.MatlabError) as e:
        R.imresize3(np.zeros((1, 4, 5), np.float32), (2, 2, 2))
    assert e.value.identifier == "MATLAB:badsubscript"
    with pytest.raises(R.MatlabError):
        R.imresize3(np.zeros((4, 5), np.float32), (2, 2, 2))


def test_unported_inputs_are_named():
    with pytest.raises(NotImplementedError):
        R.imresize3(np.zeros((3, 4, 5), bool), (2, 2, 2))
    with pytest.raises(NotImplementedError):
        R.imresize3(np.zeros((3, 4, 5), np.float32), (2, 2, 2), "lanczos3")
    with pytest.raises(TypeError):
        R.imresize(np.zeros((4, 5), np.int64), (2, 2))


def test_imresize3_fast_takes_a_numpy_type_as_its_out_class():
    # outClass has the class_name kind: a NumPy type is its MATLAB class name (str() of the type
    # failed before)
    v = np.ascontiguousarray(np.transpose(volume(7, (21, 17, 13), "uint16"), (2, 0, 1)))
    for cls, numpy_type in (("single", np.float32), ("double", np.dtype("float64"))):
        want = R.imresize3_fast(v, (9, 14, 11), "linear", None, cls)
        got = R.imresize3_fast(v, (9, 14, 11), "linear", None, numpy_type)
        assert got.dtype == want.dtype and np.array_equal(got, want)


def test_imresize3_fast_region_equals_the_crop_of_the_whole():
    v = np.ascontiguousarray(np.transpose(volume(6, (41, 37, 23), "uint16"), (2, 0, 1)))
    whole = R.imresize3_fast(v, (15, 28, 27), "cubic")
    box = (3, 5, 2, 12, 20, 27)
    part = R.imresize3_fast(v, (15, 28, 27), "cubic", box)
    assert np.array_equal(part, whole[3:12, 5:20, 2:27])
    d = R.imresize3_fast(v, (15, 28, 27), "cubic", box, "double")
    assert d.dtype == np.float64


def test_imresize3_average_squeezes_as_matlab():
    v = np.ones((6, 2, 10), np.uint16)
    r = R.imresize3_average(v, (2, 2, 2))
    assert r.shape == (5, 3)  # MATLAB's squeeze drops y = 1: a 5 x 3 image
    assert r.dtype == np.float32 and np.all(r == 1)


def _rule_sum(vals, wts) -> float:
    """MATLAB's sum of one output (the rule of resize.py): the first NaN in tap order stays"""
    acc = 0.0
    for x, w in zip(vals, wts, strict=True):
        x = float(x)
        if x != x:
            b = np.float64(x).view(np.uint64) | np.uint64(0x0008000000000000)
            return float(b.view(np.float64))
        if w == 0 and math.isinf(x):
            return float(np.float64(dtypes.NAN))
        r = acc + w * x
        if r != r:
            return float(np.float64(dtypes.NAN))
        acc = r
    return acc


def test_a_nan_sum_keeps_the_first_nan_in_tap_order():
    b = PAYLOADS["double"]
    na, nb = (np.uint64(v).view(np.float64) for v in b[:2])
    a = np.array([[nb, na], [na, nb], [1.0, np.inf]])  # (3 taps, 2 lines)
    w = np.array([[0.5, 0.5, 0.0], [0.0, 0.5, 0.5]])
    ind = np.array([[0, 1, 2], [0, 1, 2]])
    out = R.resize_axis(a, 0, w, ind).view(np.uint64)
    # line 0: nb first; line 1: na first; Inf * 0 in line 0 after the NaN does not matter
    assert out[:, 0].tolist() == [b[1], b[1]] and out[:, 1].tolist() == [b[0], b[0]]
    inf0 = R.resize_axis(np.array([[np.inf], [1.0]]), 0, np.array([[0.0, 1.0]]), np.array([[0, 1]]))
    assert inf0.view(np.uint64)[0, 0] == b[0]  # Inf * 0: MATLAB's NaN


@pytest.mark.parametrize("cls", ["single", "double"])
def test_nan_payloads_follow_the_rule_on_every_path_and_thread_count(cls):
    ft, ut = (np.float32, np.uint32) if cls == "single" else (np.float64, np.uint64)
    v = np.ascontiguousarray(np.transpose(payloads(7, (31, 29, 9), cls, 0.1), (2, 0, 1)))
    for axis, n_out in ((2, 13), (1, 47), (0, 4)):
        n = v.shape[axis]
        w, ind = R.contributions(n, n_out, n_out / n, "cubic", None, True)
        src = np.moveaxis(v, axis, 0)
        want = np.empty((n_out,) + src.shape[1:], dtype=ft)
        for i in range(n_out):
            for idx in np.ndindex(src.shape[1:]):
                want[(i,) + idx] = ft(_rule_sum([src[(k,) + idx] for k in ind[i]], w[i]))
        want = np.moveaxis(want, 0, axis)
        for nt in (1, 3, 8):
            got = R.resize_axis(v, axis, w, ind, num_threads=nt)
            assert np.array_equal(got.view(ut), want.view(ut)), (axis, nt)
        if _ext.best_isa() >= 1:
            outer = int(np.prod(v.shape[:axis]))
            inner = int(np.prod(v.shape[axis + 1 :]))
            fn = getattr(_ext, f"resize_{cls}")
            for isa in (0, 1):
                dst = np.empty(want.shape, dtype=ft)
                fn(v, dst, outer, n, inner, n_out, w, ind, 2, isa)
                assert np.array_equal(dst.view(ut), want.view(ut)), (axis, isa)
