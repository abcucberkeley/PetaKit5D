"""Tests of petakit5d.matlabcompat.morphology against MATLAB R2024b, bit for bit: strel's
neighborhoods and decompositions, imerode, imdilate, imopen, imclose, bwconvhull and the median
filters.

Oracle cases call MATLAB's built-ins themselves; the MATLAB arguments are in MATLAB's form
((y, x, z) Fortran arrays) and the Python side converts them (``from_matlab``) and its results
back (``to_matlab``) by hand until the oracle's roles (task F14) are in the tree.

MATLAB tests whose intent is covered here (their classes belong to other tasks):
- DeconMaskErosionTest/masksEqualImerode: the imerode code of decon_mask_edge_erosion_before
  (``erosion.*`` cases: 6 shapes x 6 sizes x EdgeErosion 1 to 41, errors included);
- DeconMaskErosionTest/mexEqualsImerodeForOtherNeighborhoods: random, asymmetric and even-sized
  neighborhoods, alone and in sequences with morphop's padding (``nh.*``, ``seq.*``), and the
  line along z (``zline.*``);
- DeconBackProjectorTest/convexHullsMatchBwconvhull, sliceHullsMatchBwconvhull: ``hull.*``;
- DeconBackProjectorTest/medianFilterMatchesMedfilt3: ``medfilt3.*``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.matlabcompat import morphology as M
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value
from petakit5d.testing.compare import tol
from petakit5d.testing.synth import CounterRNG

# ---------------------------------------------------------------------------------------------
# inputs in MATLAB's form (bit-stable generators)
# ---------------------------------------------------------------------------------------------


def uniform(seed: int, shape, stream: str = "n03") -> np.ndarray:
    return CounterRNG(seed, stream).uniform(shape)


def mask(seed: int, shape, p: float, stream: str = "n03-mask") -> np.ndarray:
    """MATLAB's rand(shape) > p, from our generator (Fortran, MATLAB's shape)"""
    return np.asfortranarray(uniform(seed, shape, stream) > p)


def nhood(seed: int, shape) -> np.ndarray:
    """a random 0/1 neighborhood with at least one true element (MATLAB's form)"""
    a = uniform(seed, shape, "n03-nhood") > 0.5
    k = int(CounterRNG(seed, "n03-nhood-k").integers(0, a.size - 1, 1)[0])
    a.reshape(-1, order="F")[k] = True
    return np.asfortranarray(a)


def integers(seed: int, shape, lo: int, hi: int, stream: str = "n03-int") -> np.ndarray:
    return np.asfortranarray(CounterRNG(seed, stream).integers(lo, hi, shape).astype(np.float64))


# ---------------------------------------------------------------------------------------------
# the Python sides (MATLAB's form in and out)
# ---------------------------------------------------------------------------------------------


def strel_args(args) -> list:
    """strel's MATLAB arguments in Python's form (per-axis vectors, neighborhoods)"""
    args = list(args)
    if not args or not isinstance(args[0], str):
        return [from_matlab(np.asarray(args[0]))]
    name = args[0].lower()
    if name == "cuboid":
        m, n, p = (int(v) for v in args[1])
        return [name, (p, m, n)]
    if name == "arbitrary":
        return [name, from_matlab(np.asarray(args[1]))]
    return args


def m_vec(v) -> np.ndarray:
    """a Python per-axis vector ((z,) y, x) as MATLAB's row (y, x(, z)), double"""
    v = list(v)
    return np.array(v[1:] + v[:1] if len(v) == 3 else v, dtype=np.float64)


def py_strel_structure(*args):
    """[getnhood, the decomposition's neighborhoods, isdecompositionorthogonal, getpadsize,
    getneighbors] in MATLAB's form"""
    se = M.strel(*strel_args(args))
    ul, lr = se.getpadsize()
    o = se.getneighbors()
    o = o[:, [1, 2, 0]] if o.shape[1] == 3 else o

    def nhood_m(s):
        # a line shorter than 1 is strel's default object, whose neighborhood is [] (double)
        empty_line = s.type == "line" and s.matlab_nhood.size == 0
        return np.zeros((0, 0)) if empty_line else s.matlab_nhood

    return [
        nhood_m(se),
        [nhood_m(s) for s in se.decompose()],
        se.is_decomposition_orthogonal(),
        m_vec(ul),
        m_vec(lr),
        o.astype(np.float64),
    ]


STREL_CODE = """
se = strel(in{:});
seq = decompose(se);
nh = cell(1, numel(seq));
for k = 1:numel(seq), nh{k} = getnhood(seq(k)); end
[ul, lr] = getpadsize(se);
out = {getnhood(se), nh, isdecompositionorthogonal(se), ul, lr, getneighbors(se)};
"""


def py_morph(func: str, a, *sargs):
    """func(A, strel(sargs{:})) in MATLAB's form"""
    f = getattr(M, func)
    return to_matlab(f(from_matlab(np.asarray(a)), M.strel(*strel_args(sargs))))


MORPH_CODE = "out = {feval(in{1}, in{2}, strel(in{3:end}))};"


def py_morph_nh(func: str, a, nh):
    """func(A, nhood) with the neighborhood array itself (morphop_fast's numeric path)"""
    f = getattr(M, func)
    return to_matlab(f(from_matlab(np.asarray(a)), from_matlab(np.asarray(nh))))


def py_morph_seq(func: str, a, *nhs):
    """func(A, [strel(nh1), strel(nh2), ...]): a strel array, applied as morphop applies it"""
    f = getattr(M, func)
    ses = [M.strel(from_matlab(np.asarray(nh))) for nh in nhs]
    return to_matlab(f(from_matlab(np.asarray(a)), ses))


SEQ_CODE = """
se = strel(in{3});
for k = 4:numel(in), se(k - 2) = strel(in{k}); end
out = {feval(in{1}, in{2}, se)};
"""


def py_zline(a, n):
    """imerode along z by a line of n voxels: MATLAB permutes z first and erodes by true(n, 1)"""
    return to_matlab(M.imerode(from_matlab(np.asarray(a)), np.ones((int(n), 1, 1), dtype=bool)))


ZLINE_CODE = "out = {permute(imerode(permute(in{1}, [3, 1, 2]), true(in{2}, 1)), [2, 3, 1])};"

# ---------------------------------------------------------------------------------------------
# strel cases
# ---------------------------------------------------------------------------------------------


def strel_case(cid: str, args, inputs=None) -> Case:
    return Case(
        cid,
        code=STREL_CODE,
        args=args,
        nargout=6,
        outputs=("nhood", "seq", "orthogonal", "pad_ul", "pad_lr", "offsets"),
        python=py_strel_structure,
        inputs=inputs or {},
    )


STREL_CASES: list[Case] = []
STREL_CASES += [strel_case(f"strel.disk.r{r}", ["disk", r]) for r in range(41)]
# N = 0, 6, 8 (PetaKit5D uses the default 4): the radii where the decomposition changes form
# (r < 3 has none; extra lines from r = 4) and a spread above
STREL_CASES += [
    strel_case(f"strel.disk.r{r}.n{n}", ["disk", r, n])
    for n in (0, 6, 8)
    for r in (0, 1, 2, 3, 4, 5, 7, 10, 15, 20, 27, 33, 40)
]
STREL_CASES += [strel_case(f"strel.disk.r{r}", ["disk", r]) for r in (50, 64, 75, 100)]
STREL_CASES += [strel_case(f"strel.sphere.r{r}", ["sphere", r]) for r in range(11)]
STREL_CASES += [strel_case(f"strel.square.w{w}", ["square", w]) for w in range(10)]
STREL_CASES += [
    strel_case(f"strel.rectangle.{m}x{n}", ["rectangle", [m, n]])
    for m, n in [(1, 5), (5, 1), (3, 4), (0, 3), (6, 6), (2, 9), (1, 1)]
]
STREL_CASES += [strel_case(f"strel.cube.w{w}", ["cube", w]) for w in range(6)]
STREL_CASES += [
    strel_case(f"strel.cuboid.{m}x{n}x{p}", ["cuboid", [m, n, p]])
    for m, n, p in [(2, 3, 4), (1, 1, 5), (5, 1, 1), (3, 3, 1), (1, 4, 1), (4, 4, 4)]
]
STREL_CASES += [
    strel_case(
        f"strel.line.len{str(n).replace('.', 'p')}.deg{str(d).replace('.', 'p')}", ["line", n, d]
    )
    for n in (0, 0.5, 1, 2, 3, 4.5, 7, 10, 21)
    for d in (0, 90, 45, 30, 135, -20, 200, 17.5)
]
STREL_CASES += [
    strel_case(f"strel.periodicline.p{p}.v{v[0]}_{v[1]}", ["periodicline", p, list(v)])
    for p in range(5)
    for v in [(1, 0), (0, 1), (1, 1), (-1, 1), (2, 1), (1, -2), (-2, -1)]
]
STREL_CASES += [strel_case(f"strel.diamond.r{r}", ["diamond", r]) for r in range(13)]
STREL_CASES += [strel_case(f"strel.octagon.r{r}", ["octagon", r]) for r in (0, 3, 6, 9, 12, 15)]
STREL_CASES += [
    strel_case(f"strel.pair.{a}_{b}", ["pair", [a, b]])
    for a, b in [(0, 0), (2, 1), (-3, 2), (0, -4), (5, 0)]
]
STREL_CASES += [
    strel_case(
        f"strel.arbitrary.{k}",
        [In("nh")],
        {"nh": Value(lambda k=k, s=s: nhood(1000 + k, s))},
    )
    for k, s in enumerate(
        [(3, 3), (2, 4), (4, 1), (1, 6), (5, 5), (3, 3, 3), (2, 3, 2), (4, 4, 2), (1, 1, 3), (6, 2)]
    )
]
STREL_CASES += [
    strel_case(
        f"strel.ones.{'x'.join(map(str, s))}",
        [In("nh")],
        {"nh": Value(lambda s=s: np.ones(s, bool, order="F"))},
    )
    for s in [(3, 3), (5, 5), (1, 9), (7, 1), (3, 3, 3), (2, 2), (4, 6, 2), (3, 1, 5)]
]

# ---------------------------------------------------------------------------------------------
# erosion, dilation, opening, closing
# ---------------------------------------------------------------------------------------------

MASK2 = (53, 61)
MASK3 = (41, 46, 6)
SPHERE_MASK = (25, 27, 22)


def morph_case(cid, func, a_make, sargs, marks=()) -> Case:
    return Case(
        cid,
        code=MORPH_CODE,
        args=[func, In("A"), *sargs],
        python=py_morph,
        inputs={"A": Value(a_make)},
        marks=marks,
    )


MORPH_CASES: list[Case] = []
FUNCS = ("imerode", "imdilate", "imopen", "imclose")
for r in list(range(16)) + [20, 25, 30, 40]:
    for f in FUNCS:
        if f in ("imopen", "imclose") and r not in (1, 2, 3, 5, 7, 8, 11, 15, 20):
            continue
        MORPH_CASES.append(
            morph_case(f"{f}.disk{r}.2d", f, lambda r=r: mask(r, MASK2, 0.35), ["disk", r])
        )
        MORPH_CASES.append(
            morph_case(f"{f}.disk{r}.3d", f, lambda r=r: mask(100 + r, MASK3, 0.3), ["disk", r])
        )
for r in range(1, 6):
    for f in FUNCS:
        MORPH_CASES.append(
            morph_case(
                f"{f}.sphere{r}.3d", f, lambda r=r: mask(200 + r, SPHERE_MASK, 0.7), ["sphere", r]
            )
        )
        MORPH_CASES.append(
            morph_case(f"{f}.sphere{r}.2d", f, lambda r=r: mask(210 + r, MASK2, 0.6), ["sphere", r])
        )
OTHER_STRELS = [
    ("square3", ["square", 3]),
    ("square7", ["square", 7]),
    ("rect3x8", ["rectangle", [3, 8]]),
    ("rect1x5", ["rectangle", [1, 5]]),
    ("line9_0", ["line", 9, 0]),
    ("line9_90", ["line", 9, 90]),
    ("line11_45", ["line", 11, 45]),
    ("line10_30", ["line", 10, 30]),
    ("diamond4", ["diamond", 4]),
    ("diamond7", ["diamond", 7]),
    ("octagon6", ["octagon", 6]),
    ("cube3", ["cube", 3]),
    ("cuboid2x3x4", ["cuboid", [2, 3, 4]]),
    ("periodicline2_1_1", ["periodicline", 2, [1, 1]]),
    ("pair2_m1", ["pair", [2, -1]]),
    ("disk3n0", ["disk", 3, 0]),
    ("disk9n8", ["disk", 9, 8]),
    ("disk6n6", ["disk", 6, 6]),
]
for name, sargs in OTHER_STRELS:
    for f in FUNCS:
        MORPH_CASES.append(
            morph_case(f"{f}.{name}.2d", f, lambda n=name: mask(hash_seed(n), MASK2, 0.4), sargs)
        )
        MORPH_CASES.append(
            morph_case(
                f"{f}.{name}.3d", f, lambda n=name: mask(hash_seed(n) + 1, MASK3, 0.3), sargs
            )
        )


def hash_seed(name: str) -> int:
    """a seed from a name, stable across processes"""
    return int.from_bytes(name.encode()[:7].ljust(7, b"\0"), "little") % 1000003


# border cases: masks that are all true, all false, single pixels at borders, thin arrays
EDGE_MASKS = {
    "alltrue2": lambda: np.ones((12, 15), bool, order="F"),
    "allfalse3": lambda: np.zeros((9, 10, 4), bool, order="F"),
    "corner2": lambda: corner_mask((14, 13)),
    "corner3": lambda: corner_mask((9, 8, 5)),
    "row": lambda: mask(31, (1, 40), 0.3),
    "col": lambda: mask(32, (40, 1), 0.3),
    "zcol": lambda: mask(33, (1, 1, 30), 0.3),
    "slab": lambda: mask(34, (1, 30, 7), 0.3),
    "tiny": lambda: mask(35, (2, 2, 2), 0.4),
}


def corner_mask(shape) -> np.ndarray:
    a = np.zeros(shape, bool, order="F")
    a[(0,) * len(shape)] = True
    a[tuple(s - 1 for s in shape)] = True
    a[(0,) + tuple(s - 1 for s in shape[1:])] = True
    return a


for mname, make in EDGE_MASKS.items():
    for sname, sargs in [
        ("disk5", ["disk", 5]),
        ("sphere2", ["sphere", 2]),
        ("line5_90", ["line", 5, 90]),
    ]:
        for f in FUNCS:
            MORPH_CASES.append(morph_case(f"{f}.{sname}.{mname}", f, make, sargs))
# empty images and empty neighborhoods
for f in FUNCS:
    c = morph_case(f"{f}.disk3.empty", f, lambda: np.zeros((0, 5), bool, order="F"), ["disk", 3])
    if f == "imopen":
        # MATLAB's bwpack path needs the row count of the packed image (0 here): it errors
        c.expect_error = "images:morphop:missingPackedM"
        c.py_error = M.MatlabError
    MORPH_CASES.append(c)
    MORPH_CASES.append(
        morph_case(f"{f}.sphere0.3d", f, lambda: mask(36, (8, 9, 4), 0.5), ["sphere", 0])
    )
    MORPH_CASES.append(
        morph_case(f"{f}.line0.2d", f, lambda: mask(37, (8, 9), 0.5), ["line", 0, 0])
    )
    MORPH_CASES.append(
        morph_case(f"{f}.false3x3.2d", f, lambda: mask(38, (8, 9), 0.5), [np.zeros((3, 3), bool)])
    )

# random neighborhoods as arrays (numeric path) and as strel objects, 2D and 3D, even sizes
NH_CASES: list[Case] = []
NH_SHAPES = [
    (2, 2),
    (4, 3),
    (2, 5),
    (6, 1),
    (1, 4),
    (5, 6),
    (3, 3),
    (4, 4),
    (2, 3, 2),
    (3, 2, 3),
    (4, 4, 3),
    (1, 3, 2),
    (5, 5, 1),
    (6, 4, 2),
    (3, 6),
    (2, 2, 2),
]
for k, s in enumerate(NH_SHAPES):
    if k in (4, 6, 9, 11, 12, 14):  # near-duplicates of the shapes kept
        continue
    for dim, ashape in (("2d", (29, 31)), ("3d", (23, 26, 5))):
        for f in FUNCS:
            NH_CASES.append(
                Case(
                    f"nh.{f}.{k}.{dim}",
                    code="out = {feval(in{1}, in{2}, in{3})};",
                    args=[f, In("A"), In("nh")],
                    python=py_morph_nh,
                    inputs={
                        "A": Value(lambda k=k, a=ashape: mask(300 + k, a, 0.3)),
                        "nh": Value(lambda k=k, s=s: nhood(400 + k, s)),
                    },
                )
            )
            NH_CASES.append(
                Case(
                    f"nh.{f}.{k}.{dim}.strel",
                    code=MORPH_CODE,
                    args=[f, In("A"), In("nh")],
                    python=py_morph,
                    inputs={
                        "A": Value(lambda k=k, a=ashape: mask(500 + k, a, 0.3)),
                        "nh": Value(lambda k=k, s=s: nhood(600 + k, s)),
                    },
                )
            )

# strel arrays: sequences of 2 or 3 neighborhoods (padded when not orthogonal)
SEQ_CASES: list[Case] = []
SEQS = [
    [(3, 1), (1, 3)],
    [(1, 5), (3, 3)],
    [(3, 3), (3, 3)],
    [(2, 4), (4, 2), (3, 3)],
    [(5, 1), (2, 2)],
    [(1, 4), (1, 2)],
    [(3, 3, 2), (2, 2)],
    [(4, 1), (1, 4), (2, 2)],
    [(6, 3), (2, 5)],
    [(1, 1, 3), (3, 1)],
]
for k, shapes in enumerate(SEQS):
    for dim, ashape in (("2d", (31, 28)), ("3d", (22, 25, 6))):
        for f in ("imerode", "imdilate"):
            ins = {"A": Value(lambda k=k, a=ashape: mask(700 + k, a, 0.25))}
            for j, s in enumerate(shapes):
                ins[f"nh{j}"] = Value(lambda k=k, j=j, s=s: nhood(800 + 10 * k + j, s))
            SEQ_CASES.append(
                Case(
                    f"seq.{f}.{k}.{dim}",
                    code=SEQ_CODE,
                    args=[f, In("A"), *[In(f"nh{j}") for j in range(len(shapes))]],
                    python=py_morph_seq,
                    inputs=ins,
                )
            )


def one_sided(shape, side: str) -> np.ndarray:
    """a neighborhood (MATLAB's form) true on one side of its center only (and the center),
    so that a sequence of them reaches further on that side: the padding morphop gives the
    sequence then decides the voxels near the borders"""
    a = np.zeros(shape, bool, order="F")
    c = [(s + 1) // 2 - 1 for s in shape]
    idx = np.indices(shape)
    if side == "low":
        sel = np.all([idx[d] <= c[d] for d in range(len(shape))], axis=0)
    else:
        sel = np.all([idx[d] >= c[d] for d in range(len(shape))], axis=0)
    a[sel] = True
    return a


def border_mask(seed: int, shape) -> np.ndarray:
    """random, with the outer frame true (where one-sided sequences differ by their padding)"""
    a = mask(seed, shape, 0.15)
    for d in range(len(shape)):
        sl = [slice(None)] * len(shape)
        for e in (0, 1, shape[d] - 2, shape[d] - 1):
            sl[d] = e
            a[tuple(sl)] = True
    return a


ONE_SIDED = [
    (("low", (1, 5)), ("low", (1, 5))),
    (("high", (5, 1)), ("high", (5, 1)), ("high", (3, 3))),
    (("low", (4, 4)), ("high", (2, 2))),
    (("low", (3, 3)), ("low", (1, 6))),
    (("high", (3, 3, 3)), ("high", (3, 1))),
]
for k, seq in enumerate(ONE_SIDED):
    for dim, ashape in (("2d", (24, 27)), ("3d", (20, 22, 7))):
        for f in ("imerode", "imdilate"):
            ins = {"A": Value(lambda k=k, a=ashape: border_mask(1950 + k, a))}
            for j, (side, s) in enumerate(seq):
                ins[f"nh{j}"] = Value(lambda side=side, s=s: one_sided(s, side))
            SEQ_CASES.append(
                Case(
                    f"seq.onesided.{f}.{k}.{dim}",
                    code=SEQ_CODE,
                    args=[f, In("A"), *[In(f"nh{j}") for j in range(len(seq))]],
                    python=py_morph_seq,
                    inputs=ins,
                )
            )

ZLINE_CASES = [
    Case(
        f"zline.n{n}.{'x'.join(map(str, s))}",
        code=ZLINE_CODE,
        args=[In("A"), n],
        python=py_zline,
        inputs={"A": Value(lambda n=n, s=s: mask(900 + n, s, 0.1))},
    )
    for n, s in [
        (1, (10, 12, 8)),
        (3, (20, 30, 40)),
        (9, (20, 30, 40)),
        (17, (12, 9, 30)),
        (41, (8, 7, 60)),
        (81, (5, 6, 50)),
        (5, (4, 5, 2)),
    ]
]

# ---------------------------------------------------------------------------------------------
# DeconMaskErosionTest: decon_mask_edge_erosion_before's imerode code
# ---------------------------------------------------------------------------------------------

EROSION_CODE = """
A0 = in{1} > 0;
if in{2}, A0 = imdilate(A0, ones(3, 3)); end
res = cell(1, max(in{3}));
for E = in{3}
    try
        mask = A0;
        if sum(mask, 'all') == numel(mask)
            mask([1 : E, end - E + 1 : end], :, :) = false;
            mask(:, [1 : E, end - E + 1 : end], :) = false;
            if ndims(mask) == 3
                mask(:, :, [1 : E, end - E + 1 : end]) = false;
            end
        else
            mask([1, end], :, :) = false;
            mask(:, [1, end], :) = false;
            switch ndims(mask)
                case 2
                    mask = imerode(mask, strel('disk', E - 1));
                case 3
                    mask(:, :, [1, end]) = false;
                    mask = imerode(mask, strel('disk', E - 1));
                    mask = permute(mask, [3, 1, 2]);
                    mask = imerode(mask, true(E * 2 - 1, 1));
                    mask = permute(mask, [2, 3, 1]);
            end
        end
        res{E} = mask;
    catch ME
        res{E} = ME.identifier;
    end
end
out = {res};
"""


class _IndexError(Exception):
    def __init__(self, identifier: str):
        super().__init__(identifier)
        self.identifier = identifier


# MATLAB's error for mask([1 : E, end - E + 1 : end], ...) = false when E > size (the range has
# 0 and below; measured: badsubscript for every size and dimension)
INDEX_ERROR = "MATLAB:badsubscript"


def edge_erosion_before(m: np.ndarray, e: int) -> np.ndarray:
    """decon_mask_edge_erosion_before in Python's form (m: (y, x) or (z, y, x) bool)"""
    m = m.copy()
    nd = len(M._msize(m))
    # MATLAB's dimensions 1, 2, 3 are Python's axes -2, -1, -3
    axes = [-2, -1] + ([-3] if nd == 3 else [])
    if m.all():
        for ax in axes:
            n = m.shape[ax]
            if e > n:
                raise _IndexError(INDEX_ERROR)
            sl = [slice(None)] * m.ndim
            sl[ax] = np.r_[0:e, n - e : n]
            m[tuple(sl)] = False
        return m
    for ax in axes:
        sl = [slice(None)] * m.ndim
        sl[ax] = [0, m.shape[ax] - 1]
        m[tuple(sl)] = False
    m = M.imerode(m, M.strel("disk", e - 1))
    if nd == 3:
        m = M.imerode(m, np.ones((2 * e - 1, 1, 1), dtype=bool))
    return m


def py_erosion_matrix(a, blobs, es):
    """the masks (or error identifiers) for EdgeErosion in es; None (MATLAB's []) for the
    others up to max(es)"""
    m = from_matlab(np.asarray(a)) > 0
    if blobs:
        m = M.imdilate(m, np.ones((3, 3), dtype=bool))
    es = [int(e) for e in np.atleast_1d(es)]
    res: list = [None] * max(es)
    for e in es:
        try:
            res[e - 1] = to_matlab(edge_erosion_before(m, e))
        except _IndexError as err:
            res[e - 1] = err.identifier
    return res


EROSION_SIZES = {
    "odd3D": (37, 41, 23),
    "small": (5, 5, 5),
    "twoSlices": (17, 19, 2),
    "oneRow": (1, 30, 20),
    "oneColumn": (30, 1, 20),
    "image2D": (45, 47),
}


def erosion_mask(shape_name: str, sz) -> np.ndarray:
    """DeconMaskErosionTest.makeMask in MATLAB's form ('blobs': the seeds, dilated by the case)"""
    sz3 = tuple(sz) + (1,) * (3 - len(sz))
    seed = hash_seed(shape_name + "x".join(map(str, sz)))
    if shape_name == "random":
        return mask(seed, sz, 0.3, "n03-erosion")
    if shape_name == "skewed":
        _, x, z = np.meshgrid(
            np.arange(1, sz3[0] + 1),
            np.arange(1, sz3[1] + 1),
            np.arange(1, sz3[2] + 1),
            indexing="ij",
        )
        x = x.astype(np.float64)
        z = z.astype(np.float64)
        a = (x > 0.6 * z + sz3[1] * 0.05) & (x < 0.6 * z + sz3[1] * 0.5)
        return np.asfortranarray(a.reshape(sz, order="F"))
    a = np.zeros(sz, bool, order="F")
    if shape_name == "thin":
        c = max(1, math.floor(sz[1] / 2 + 0.5))
        cols = [c - 1 + d for d in range(min(1, sz[1] - 1) + 1)]
        a[:, cols, ...] = True
        return a
    if shape_name in ("single", "singleFalse"):
        k = math.ceil(a.size / 2) - 1
        a.reshape(-1, order="F")[k] = True
        return a if shape_name == "single" else np.asfortranarray(~a)
    if shape_name == "blobs":
        return mask(seed, sz, 0.99, "n03-erosion")
    raise ValueError(shape_name)


EROSION_CASES = [
    Case(
        f"erosion.{shape}.{name}",
        code=EROSION_CODE,
        args=[In("A"), shape == "blobs", list(range(1, 42))],
        python=py_erosion_matrix,
        inputs={"A": Value(lambda shape=shape, sz=sz: erosion_mask(shape, sz))},
    )
    for shape in ("random", "skewed", "thin", "single", "singleFalse", "blobs")
    for name, sz in EROSION_SIZES.items()
]
EROSION_CASES.append(
    # demoSettingOnASkewedVolume: EdgeErosion 20, zero borders, a few isolated holes
    Case(
        "erosion.demo_skewed_e20",
        code=EROSION_CODE,
        args=[In("A"), False, 20],
        python=py_erosion_matrix,
        inputs={"A": Value(lambda: demo_skewed())},
    )
)


def demo_skewed() -> np.ndarray:
    a = erosion_mask("skewed", (150, 200, 90))
    for y, x, z in [(70, 70, 40), (20, 90, 60), (100, 110, 75)]:
        a[y - 1, x - 1, z - 1] = False
    return a


# ---------------------------------------------------------------------------------------------
# convex hulls
# ---------------------------------------------------------------------------------------------


def py_bwconvhull(a, *args):
    return M.bwconvhull(from_matlab(np.asarray(a)), *args)


def hull_mask(kind: str, k: int) -> np.ndarray:
    s = (48 + 7 * (k % 5), 47 + 5 * (k % 4))  # over 2 KB: hashes in the manifest
    if kind == "sparse":
        return mask(1100 + k, s, 0.97, "n03-hull")
    if kind == "dense":
        return mask(1200 + k, s, 0.6, "n03-hull")
    if kind == "line":
        a = np.zeros(s, bool, order="F")
        n = min(s) - 2
        dy, dx = [(1, 0), (0, 1), (1, 1), (1, 2), (2, 1), (1, -1), (3, 1), (1, 3)][k % 8]
        for i in range(n // max(abs(dy), abs(dx))):
            a[1 + i * dy, (s[1] // 2 + i * dx) % s[1]] = True
        return a
    if kind == "point":
        a = np.zeros(s, bool, order="F")
        a[k % s[0], (3 * k) % s[1]] = True
        return a
    raise ValueError(kind)


HULL_CASES: list[Case] = []
for kind in ("sparse", "dense", "line", "point"):
    for k in range(8):
        HULL_CASES.append(
            Case(
                f"hull.union.{kind}.{k}",
                matlab="bwconvhull",
                args=[In("A")],
                python=py_bwconvhull,
                inputs={"A": Value(lambda kind=kind, k=k: hull_mask(kind, k))},
            )
        )
for k in range(6):
    for conn in (4, 8):
        HULL_CASES.append(
            Case(
                f"hull.objects.c{conn}.{k}",
                matlab="bwconvhull",
                args=[In("A"), "objects", conn],
                python=py_bwconvhull,
                inputs={"A": Value(lambda k=k: mask(1300 + k, (48, 47), 0.93, "n03-hull"))},
            )
        )
HULL_CASES += [
    Case(
        f"hull.union.{name}",
        matlab="bwconvhull",
        args=[In("A")],
        python=py_bwconvhull,
        inputs={"A": Value(make)},
    )
    for name, make in [
        ("empty", lambda: np.zeros((10, 12), bool, order="F")),
        ("row", lambda: mask(1401, (1, 30), 0.8, "n03-hull")),
        ("col", lambda: mask(1402, (30, 1), 0.8, "n03-hull")),
        ("full", lambda: np.ones((7, 9), bool, order="F")),
        ("single", lambda: np.ones((1, 1), bool, order="F")),
    ]
]

SLICE_HULL_CODE = """
B = in{1};
zinds = find(sum(B, [1, 2]))';
xinds = find(sum(B, [1, 3]));
yinds = find(sum(B, [2, 3]))';
for z = zinds, B(:, :, z) = bwconvhull(B(:, :, z), "union"); end
for y = yinds, B(y, :, :) = bwconvhull(squeeze(B(y, :, :)), "union"); end
for x = xinds, B(:, x, :) = bwconvhull(squeeze(B(:, x, :)), "union"); end
out = {B};
"""


def py_slice_hulls(a):
    """the convex hull loops of omw_backprojector_generation (z, then y, then x slices; the
    slice lists from the mask before the loops)"""
    b = from_matlab(np.asarray(a))
    zinds = np.flatnonzero(b.any(axis=(1, 2)))
    yinds = np.flatnonzero(b.any(axis=(0, 2)))
    xinds = np.flatnonzero(b.any(axis=(0, 1)))
    b = M.bwconvhull_slices(b, 0, zinds)
    b = M.bwconvhull_slices(b, 1, yinds)
    b = M.bwconvhull_slices(b, 2, xinds)
    return to_matlab(b)


def box_sum(a: np.ndarray, r: int) -> np.ndarray:
    """the sum of each (2r + 1)^3 box of an integer array (zeros outside), exact"""
    s = a.astype(np.int64)
    for ax in range(3):
        # s at index j of p is s[j - r - 1]; the box of j is p[j + 1 .. j + 2r + 1]
        p = np.pad(s, [(r + 1, r) if d == ax else (0, 0) for d in range(3)])
        c = np.cumsum(p, axis=ax)
        n = s.shape[ax]
        s = np.take(c, np.arange(2 * r + 1, 2 * r + 1 + n), axis=ax) - np.take(
            c, np.arange(n), axis=ax
        )
    return s


def blob_volume(k: int, shape) -> np.ndarray:
    """a few blobs (OTF-mask like): a box-smoothed random integer field above its 80th
    percentile (integers only, so the bytes are the same everywhere)"""
    u = CounterRNG(1500 + k, "n03-blob").integers(0, 99, shape)
    s = box_sum(u, 2)
    t = np.sort(s.reshape(-1))[int(0.8 * s.size)]
    return np.asfortranarray(s > t)


HULL_CASES += [
    Case(
        f"hull.slices.{k}",
        code=SLICE_HULL_CODE,
        args=[In("B")],
        python=py_slice_hulls,
        inputs={"B": Value(lambda k=k, s=s: blob_volume(k, s))},
    )
    for k, s in enumerate([(21, 21, 21), (33, 29, 17), (40, 36, 25), (15, 22, 9)])
]

# ---------------------------------------------------------------------------------------------
# median filters
# ---------------------------------------------------------------------------------------------

CLASSES = {
    "double": np.float64,
    "single": np.float32,
    "uint8": np.uint8,
    "uint16": np.uint16,
    "int16": np.int16,
    "int8": np.int8,
    "uint32": np.uint32,
    "int32": np.int32,
    "logical": np.bool_,
}


def med_data(seed: int, shape, cls: str) -> np.ndarray:
    """random data of a class: integers 0..9 (signed classes -5..9), uniform for floats"""
    if cls == "logical":
        return mask(seed, shape, 0.5, "n03-med")
    if cls in ("double", "single"):
        return np.asfortranarray(uniform(seed, shape, "n03-med").astype(CLASSES[cls]))
    lo = -5 if cls.startswith("int") else 0
    return np.asfortranarray(integers(seed, shape, lo, 9, "n03-med").astype(CLASSES[cls]))


def py_medfilt2(a, size, *opt):
    return M.medfilt2(from_matlab(np.asarray(a)), tuple(int(v) for v in size), *opt)


def py_medfilt3(a, size, *opt):
    m, n, p = (int(v) for v in size)
    return to_matlab(M.medfilt3(from_matlab(np.asarray(a)), (p, m, n), *opt))


MED_CASES: list[Case] = []
for cls in CLASSES:
    # images large enough that the manifests hold hashes only (more than 2 KB)
    for size in ([3, 3], [5, 5], [1, 7], [4, 4]):
        for opt in ((), ("symmetric",), ("indexed",)):
            oname = opt[0] if opt else "zeros"
            MED_CASES.append(
                Case(
                    f"medfilt2.{cls}.{size[0]}x{size[1]}.{oname}",
                    matlab="medfilt2",
                    args=[In("A"), size, *opt],
                    python=py_medfilt2,
                    inputs={
                        "A": Value(
                            lambda c=cls, s=size: med_data(
                                hash_seed(c) + s[0] * 10 + s[1], (47, 53), c
                            )
                        )
                    },
                )
            )
for cls in ("single", "double", "uint16", "uint8", "logical", "int16"):
    for size in ([3, 3, 3], [5, 3, 3], [1, 1, 5], [3, 5, 7]):
        for opt in ("symmetric", "zeros", "replicate"):
            MED_CASES.append(
                Case(
                    f"medfilt3.{cls}.{'x'.join(map(str, size))}.{opt}",
                    matlab="medfilt3",
                    args=[In("A"), size, opt],
                    python=py_medfilt3,
                    inputs={
                        "A": Value(
                            lambda c=cls, s=size: med_data(hash_seed(c) + sum(s), (13, 11, 9), c)
                        )
                    },
                )
            )
MED_CASES += [
    Case(
        f"medfilt3.single.default.{'x'.join(map(str, s))}",
        matlab="medfilt3",
        args=[In("A")],
        python=lambda a: to_matlab(M.medfilt3(from_matlab(np.asarray(a)))),
        inputs={"A": Value(lambda s=s: med_data(1600 + sum(s), s, "single"))},
    )
    for s in [(63, 63, 41), (2, 2, 2), (1, 5, 4), (3, 1, 2), (17, 2, 33), (4, 3, 1)]
]
MED_CASES += [
    Case(
        "medfilt3.double.psf_like",
        matlab="medfilt3",
        args=[In("A"), [3, 3, 3]],
        python=lambda a, s: to_matlab(M.medfilt3(from_matlab(np.asarray(a)), (3, 3, 3))),
        inputs={
            "A": Value(
                lambda: np.asfortranarray(np.floor(uniform(1700, (31, 33, 21), "n03-med") * 4096))
            )
        },
    ),
    Case(
        "medfilt2.tiny.5x5",
        matlab="medfilt2",
        args=[In("A"), [5, 5], "symmetric"],
        python=py_medfilt2,
        inputs={"A": Value(lambda: med_data(1701, (2, 3), "double"))},
    ),
]

# signed zeros and NaN: MATLAB's choice depends on its code path and on the positions in the
# window (characterized with these cases); the values agree, the sign of a zero may not
SIGNED_ZERO = tol(reason="approved 2026-11-12: sign of zero medians (medfilt2/medfilt3)", abs=0)


def zeros_data(seed: int, shape, cls: str, nan: bool) -> np.ndarray:
    vals = np.array([-0.0, 0.0, -1.0, 1.0, 2.0])
    k = CounterRNG(seed, "n03-zeros").integers(0, 4, shape)
    a = vals[k]
    if nan:
        a[uniform(seed, shape, "n03-nan") < 0.12] = np.nan
    return np.asfortranarray(a.astype(CLASSES[cls]))


MED_ZERO_CASES: list[Case] = []
for cls in ("double", "single"):
    for k in range(4):
        MED_ZERO_CASES.append(
            Case(
                f"medfilt2.signedzeros.{cls}.{k}",
                matlab="medfilt2",
                args=[In("A")],
                python=py_medfilt2,
                py_args=[In("A"), (3, 3)],
                inputs={"A": Value(lambda c=cls, k=k: zeros_data(1800 + k, (9, 11), c, False))},
                compare=SIGNED_ZERO,
            )
        )
        MED_ZERO_CASES.append(
            Case(
                f"medfilt3.signedzeros.{cls}.{k}",
                matlab="medfilt3",
                args=[In("A")],
                python=lambda a: to_matlab(M.medfilt3(from_matlab(np.asarray(a)))),
                inputs={"A": Value(lambda c=cls, k=k: zeros_data(1900 + k, (6, 7, 5), c, False))},
                compare=SIGNED_ZERO,
            )
        )

ALL_CASES = (
    STREL_CASES
    + MORPH_CASES
    + NH_CASES
    + SEQ_CASES
    + ZLINE_CASES
    + EROSION_CASES
    + HULL_CASES
    + MED_CASES
    + MED_ZERO_CASES
)


@oracle.parametrize(ALL_CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# tests without MATLAB
# ---------------------------------------------------------------------------------------------


def ref_morph(a: np.ndarray, nh: np.ndarray, erode: bool) -> np.ndarray:
    """erosion (AND, outside true) or dilation (OR of the reflected offsets, outside false) by a
    Python-form neighborhood, center floor((size + 1) / 2) in MATLAB's 1-based terms"""
    nh3 = nh.reshape((1,) * (3 - nh.ndim) + nh.shape)
    a3 = a.reshape((1,) * (3 - a.ndim) + a.shape)
    c = [(s + 1) // 2 - 1 for s in nh3.shape]
    offs = np.argwhere(nh3) - c
    if not erode:
        offs = -offs
    out = np.full(a3.shape, erode)
    if len(offs) == 0:
        return out.reshape(a.shape)
    g = np.abs(offs).max(axis=0)
    p = np.pad(a3, [(int(v), int(v)) for v in g], constant_values=erode)
    for o in offs:
        sl = tuple(slice(int(g[d] + o[d]), int(g[d] + o[d] + a3.shape[d])) for d in range(3))
        out = out & p[sl] if erode else out | p[sl]
    return out.reshape(a.shape)


def test_erosion_and_dilation_equal_the_definition():
    # random neighborhoods (even sizes, 2D and 3D) on random images, against the definition
    rng = np.random.default_rng(1)
    for _ in range(60):
        nd = rng.integers(2, 4)
        a = rng.random(tuple(rng.integers(1, 20, size=nd))) > 0.3
        nh = rng.random(tuple(rng.integers(1, 6, size=rng.integers(2, 4)))) > 0.5
        nh.reshape(-1)[rng.integers(nh.size)] = True
        assert np.array_equal(M.imerode(a, nh), ref_morph(a, nh, True))
        assert np.array_equal(M.imdilate(a, nh), ref_morph(a, nh, False))


@pytest.mark.parametrize("r", [3, 5, 11, 20])
def test_disk_decomposition_gives_the_neighborhood(r):
    # the padded sequence of periodic lines equals the erosion by the whole neighborhood
    se = M.strel("disk", r)
    assert se.minkowski_exact()
    a = np.random.default_rng(r).random((3, 60, 70)) > 0.15
    assert np.array_equal(M.imerode(a, se), ref_morph(a, se.nhood, True))
    assert np.array_equal(M.imdilate(a, se), ref_morph(a, se.nhood, False))


def test_results_do_not_depend_on_the_thread_count():
    rng = np.random.default_rng(2)
    a3 = rng.random((9, 50, 61)) > 0.2
    a1 = rng.random((1, 70, 80)) > 0.2
    for se in (M.strel("disk", 7), M.strel("sphere", 3), M.strel("square", 4)):
        for a in (a3, a1):
            ref = M.imerode(a, se, threads=1)
            for t in (2, 3, 5):
                assert np.array_equal(M.imerode(a, se, threads=t), ref)
                assert np.array_equal(M.imclose(a, se, threads=t), M.imclose(a, se, threads=1))
    v = rng.random((11, 13, 17)).astype(np.float32)
    assert np.array_equal(M.medfilt3(v, threads=1), M.medfilt3(v, threads=4))


def test_strel_forms():
    # Python's form: (z, y, x) neighborhoods, (dz, dy, dx) offsets, (z, y, x) pads
    se = M.strel("cuboid", (2, 3, 4))
    assert se.nhood.shape == (2, 3, 4) and se.matlab_nhood.shape == (3, 4, 2)
    assert se.nhood.flags.c_contiguous
    o = M.strel(np.ones((3, 1, 1), bool)).getneighbors()
    assert o.tolist() == [[-1, 0, 0], [0, 0, 0], [1, 0, 0]]
    ul, lr = M.strel("cuboid", (5, 1, 3)).getpadsize()
    assert ul == (2, 0, 1) and lr == (2, 0, 1)
    assert M.strel("disk", 5).getpadsize() == ((4, 4), (4, 4))


def test_imclose_is_the_closing_of_the_image_extended_with_false():
    rng = np.random.default_rng(3)
    for se in (M.strel("disk", 3), M.strel("sphere", 2), M.strel(rng.random((4, 3)) > 0.4)):
        a = rng.random((6, 30, 33)) > 0.6
        nh = se.nhood
        p = [(n + 1) // 2 for n in nh.reshape((1,) * (3 - nh.ndim) + nh.shape).shape]
        ap = np.pad(a, [(v, v) for v in p])
        c = ref_morph(ref_morph(ap, nh, False), nh, True)
        c = c[p[0] : p[0] + a.shape[0], p[1] : p[1] + a.shape[1], p[2] : p[2] + a.shape[2]]
        assert np.array_equal(M.imclose(a, se), c)


def test_errors_carry_matlab_identifiers():
    with pytest.raises(M.MatlabError) as e:
        M.strel("disk", 3, 5)
    assert e.value.identifier == "images:strel:invalidN"
    with pytest.raises(M.MatlabError) as e:
        M.imerode(np.ones((3, 3), bool), np.array([[0, 2, 1]]))
    assert e.value.identifier == "images:strelcheck:invalidStrelValues"
    with pytest.raises(M.MatlabError) as e:
        M.imopen(np.ones((3, 3), bool), [M.strel("disk", 1), M.strel("disk", 2)])
    assert e.value.identifier == "images:imopen:nonscalarStrel"


def test_median_filters_select_the_median():
    rng = np.random.default_rng(4)
    a = rng.random((7, 8, 9))
    m = M.medfilt3(a)
    p = np.pad(a, 1, mode="symmetric")
    win = np.lib.stride_tricks.sliding_window_view(p, (3, 3, 3)).reshape(a.shape + (27,))
    assert np.array_equal(m, np.sort(win, axis=-1)[..., 13])
    i = rng.integers(0, 50, (20, 30)).astype(np.uint16)
    p2 = np.pad(i, ((1, 1), (2, 2)))
    w2 = np.lib.stride_tricks.sliding_window_view(p2, (3, 5)).reshape(i.shape + (15,))
    assert np.array_equal(M.medfilt2(i, (3, 5)), np.sort(w2, axis=-1)[..., 7])


def test_single_median_network_equals_the_selection():
    rng = np.random.default_rng(5)
    a = (rng.random((12, 21, 37)) * 100).astype(np.float32)
    net = M.medfilt3(a)
    p = np.pad(a, 1, mode="symmetric")
    win = np.lib.stride_tricks.sliding_window_view(p, (3, 3, 3)).reshape(a.shape + (27,))
    assert np.array_equal(net, np.sort(win, axis=-1)[..., 13])
