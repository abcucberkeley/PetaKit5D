"""The edge erosion of the decon masks (petakit5d.decon.erosion) against decon_mask_edge_erosion.m.

Oracle cases (MATLAB runs the oracle tree's decon_mask_edge_erosion on the same masks, in MATLAB's
form; the roles give the Python calls the port's form):

- erosion.full.*: full masks, whose edges are set directly (DeconMaskErosionTest/
  fullMasksAndOtherErosionValues and the thin/oneColumn case of masksEqualImerode), erosions 1 to 41
  on masks thinner than the erosion (all edge there since M07, 834604d; MATLAB's index error
  before), 2D masks, numeric masks, erosions that aren't positive integers, and
  fullMasksThinnerThanTheErosionAreAllEdge's and nonDoubleErosionValuesAreTakenAsDoubles' masks;
- erosion.other.*: masks that aren't full with erosions the mex doesn't take (0, 2.5, NaN: MATLAB's
  imerode code, with strel's errors) and single and integer ones (taken as doubles since M07);
- erosion.mex.*: DeconMaskErosionTest's shapes and sizes, erosions 1 to 41, through the module's
  dispatch to the kernel (the references of K05's kernel cases: the same MATLAB calls).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.decon import erosion as E
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value, Vol
from petakit5d.testing.synth import CounterRNG

# every erosion 1 .. 41 of a mask, an error as its identifier (DeconMaskErosionTest.outcome)
OUTCOME_CODE = (
    "A = in{1}; out = cell(1, 41); for E = 1 : 41, try, out{E} = decon_mask_edge_erosion(A, E); "
    "catch ME, out{E} = ME.identifier; end, end; out = {out};"
)
ONE_OUTCOME_CODE = (
    "try, out = {decon_mask_edge_erosion(in{1}, in{2})}; catch ME, out = {ME.identifier}; end"
)
# K05's kernel cases (decon_mask_edge_erosion's mex path), whose references these share
EROSION_CODE = "A = in{1}; out = {arrayfun(@(E) decon_mask_edge_erosion(A, E), 1 : 41, 'unif', 0)};"
MASK_SIZES = {
    "odd3D": [37, 41, 23],
    "small": [5, 5, 5],
    "twoSlices": [17, 19, 2],
    "oneRow": [1, 30, 20],
    "oneColumn": [30, 1, 20],
    "image2D": [45, 47],
}


def size_name(s) -> str:
    return "x".join(str(v) for v in s)


def make_mask(shape: str, msize, seed: int = 0) -> np.ndarray:
    """DeconMaskErosionTest.makeMask in MATLAB's form (as tests/kernels/test_decon_kernels.py makes
    it, so that the cases share its references)"""
    sz3 = list(msize) + [1] * (3 - len(msize))
    rng = CounterRNG(seed, f"mask/{shape}/{size_name(msize)}")
    py = tuple(int(v) for v in ((sz3[2], sz3[0], sz3[1]) if len(msize) == 3 else msize))
    if shape == "random":
        a = rng.uniform(py) > 0.3
    elif shape == "skewed":
        z = np.arange(1, sz3[2] + 1, dtype=np.float64).reshape(-1, 1, 1)
        x = np.arange(1, sz3[1] + 1, dtype=np.float64).reshape(1, 1, -1)
        a = (x > 0.6 * z + sz3[1] * 0.05) & (x < 0.6 * z + sz3[1] * 0.5)
        a = np.broadcast_to(a, (sz3[2], sz3[0], sz3[1])).reshape(py)
    elif shape == "thin":
        a = np.zeros(py, bool)
        c = max(1, math.floor(sz3[1] / 2 + 0.5)) - 1
        a[..., c : c + min(2, sz3[1])] = True
    elif shape == "single":
        mo = np.zeros(math.prod(sz3), bool)
        mo[math.ceil(mo.size / 2) - 1] = True
        return np.reshape(mo, msize, order="F")
    elif shape == "singleFalse":
        mo = np.ones(math.prod(sz3), bool)
        mo[math.ceil(mo.size / 2) - 1] = False
        return np.reshape(mo, msize, order="F")
    else:  # blobs: imdilate(rand > 0.99, ones(3, 3))
        from petakit5d.matlabcompat.morphology import imdilate

        a = imdilate(rng.uniform(py) > 0.99, np.ones((3, 3), bool))
    return to_matlab(np.ascontiguousarray(a))


def V(make):
    return Value(make)


# ---------------------------------------------------------------------------------------------
# Python sides
# ---------------------------------------------------------------------------------------------


def outcome(mask, e, *, threads=None):
    """the eroded mask in MATLAB's form, or the identifier of the error"""
    try:
        return to_matlab(E.decon_mask_edge_erosion(mask, e, threads=threads))
    except (E.MatlabError, ValueError) as err:
        ident = getattr(err, "identifier", None)
        if ident is None:
            raise
        return ident


def outcomes(mask, *, threads=None):
    return [outcome(mask, e, threads=threads) for e in range(1, 42)]


def eroded(mask, e, *, threads=None):
    return E.decon_mask_edge_erosion(mask, e, threads=threads)


def eroded_all(mask, *, threads=None):
    return [to_matlab(eroded(mask, e, threads=threads)) for e in range(1, 42)]


CASES: list[Case] = []

# full masks (DeconMaskErosionTest/fullMasksAndOtherErosionValues)
_full = V(lambda: to_matlab(np.ones((12, 20, 22), bool)))
for _e in (1, 3, 6):
    CASES.append(
        Case(
            f"erosion.full.20x22x12.e{_e}",
            matlab="decon_mask_edge_erosion",
            args=[Vol(In("a")), _e],
            python=eroded,
            inputs={"a": _full},
            out_roles=("vol",),
        )
    )
CASES.append(
    Case(
        "erosion.full.numeric",
        code="out = {decon_mask_edge_erosion(uint8(in{1}) * 7, 2)};",
        args=[Vol(In("a"))],
        python=lambda a, threads=None: eroded(a.astype(np.uint8) * np.uint8(7), 2, threads=threads),
        inputs={"a": _full},
        out_roles=("vol",),
    )
)
# erosions 1 .. 41 of full masks, also thinner than the erosion (masksEqualImerode's thin/oneColumn
# mask is full; DeconMaskErosionTest's GPU test's full sizes, their CPU outcome)
for _ms in ([30, 1, 20], [3, 20, 20], [20, 3, 20], [20, 20, 3], [1, 30, 20], [3, 20], [30, 1],
            [45, 47], [20, 22, 12], [30, 20, 1]):  # fmt: skip
    CASES.append(
        Case(
            f"erosion.full.outcomes.{size_name(_ms)}",
            code=OUTCOME_CODE,
            args=[Vol(In("a"))],
            python=outcomes,
            inputs={"a": V(lambda m=_ms: np.ones(m, bool))},
        )
    )
# erosions that aren't positive integers, of full masks and others
_skewed = V(lambda: make_mask("skewed", [20, 22, 12]))
for _name, _e in [("e0", 0), ("eneg", -2), ("e2p5", 2.5), ("single3", np.float32(3)),
                  ("int3", np.uint8(3))]:  # fmt: skip
    CASES.append(
        Case(
            f"erosion.full.{_name}",
            code=ONE_OUTCOME_CODE,
            args=[Vol(In("a")), _e],
            python=outcome,
            inputs={"a": _full},
        )
    )
    CASES.append(
        Case(
            f"erosion.other.{_name}",
            code=ONE_OUTCOME_CODE,
            args=[Vol(In("a")), _e],
            python=outcome,
            inputs={"a": _skewed},
        )
    )
for _name, _e in [("single3", np.float32(3)), ("int3", np.uint8(3))]:
    CASES.append(
        Case(
            f"erosion.other.2d.{_name}",
            code=ONE_OUTCOME_CODE,
            args=[Vol(In("a")), _e],
            python=outcome,
            inputs={"a": V(lambda: make_mask("skewed", [45, 47]))},
        )
    )
# NaN (MATLAB's colon gives NaN, which max(., 1) makes index 1 on full masks), int16
for _name, _e in [("enan", math.nan), ("int16_3", np.int16(3))]:
    for _what, _inp in (("full", _full), ("other", _skewed)):
        CASES.append(
            Case(
                f"erosion.{_what}.{_name}",
                code=ONE_OUTCOME_CODE,
                args=[Vol(In("a")), _e],
                python=outcome,
                inputs={"a": _inp},
            )
        )
# DeconMaskErosionTest/nonDoubleErosionValuesAreTakenAsDoubles: a full mask longer than 255 (uint8
# saturated MATLAB's edge indices before M07)
_long = V(lambda: to_matlab(np.ones((12, 300, 22), bool)))
for _name, _e in [("double3", 3), ("single3", np.float32(3)), ("int3", np.uint8(3)),
                  ("int16_3", np.int16(3))]:  # fmt: skip
    CASES.append(
        Case(
            f"erosion.full.300x22x12.{_name}",
            code=ONE_OUTCOME_CODE,
            args=[Vol(In("a")), _e],
            python=outcome,
            inputs={"a": _long},
        )
    )
# DeconMaskErosionTest/fullMasksThinnerThanTheErosionAreAllEdge: all edge, also a mask that isn't
# full (the mex)
for _ms, _e in (([40, 44, 3], 8), ([30, 1, 20], 2), ([2, 20, 10], 3), ([1, 30], 2), ([5, 5, 5], 6)):
    CASES.append(
        Case(
            f"erosion.full.alledge.{size_name(_ms)}.e{_e}",
            code=ONE_OUTCOME_CODE,
            args=[Vol(In("a")), _e],
            python=outcome,
            inputs={"a": V(lambda m=_ms: np.ones(m, bool))},
        )
    )


def _not_quite_full():
    a = np.ones((40, 44, 3), bool)
    a[0, 0, 0] = False
    return a


CASES.append(
    Case(
        "erosion.other.alledge.40x44x3.e8",
        code=ONE_OUTCOME_CODE,
        args=[Vol(In("a")), 8],
        python=outcome,
        inputs={"a": V(_not_quite_full)},
    )
)

# the mex path through the module's dispatch: K05's cases, the same MATLAB calls
for _shape in ["random", "skewed", "thin", "single", "singleFalse", "blobs"]:
    for _sname, _ms in MASK_SIZES.items():
        if make_mask(_shape, _ms).all():
            continue
        CASES.append(
            Case(
                f"erosion.mex.{_shape}.{_sname}",
                code=EROSION_CODE,
                args=[Vol(In("a"))],
                python=eroded_all,
                inputs={"a": V(lambda s=_shape, m=_ms: make_mask(s, m))},
            )
        )


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# without MATLAB
# ---------------------------------------------------------------------------------------------


def test_full_masks_get_their_edges_set():
    a = np.ones((12, 20, 22), bool)
    for e in (1, 3, 6):
        m = E.decon_mask_edge_erosion(a, e)
        inner = np.zeros_like(a)
        inner[e:-e, e:-e, e:-e] = True
        assert np.array_equal(m, inner)
    # a 2D mask (and a volume of one slice) has no z edges
    m = E.decon_mask_edge_erosion(np.ones((1, 20, 22), bool), 3)
    assert m.shape == (1, 20, 22) and m[0, 3:-3, 3:-3].all() and m.sum() == 14 * 16


def test_the_input_is_not_changed():
    a = np.ones((12, 20, 22), bool)
    E.decon_mask_edge_erosion(a, 3)
    assert a.all()
    b = from_matlab(make_mask("skewed", [20, 22, 12]))
    c = b.copy()
    E.decon_mask_edge_erosion(b, 3)
    assert np.array_equal(b, c)


def test_numeric_masks_count_positive_values():
    b = from_matlab(make_mask("skewed", [20, 22, 12]))
    raw = b.astype(np.float32) * np.float32(5)  # zeros outside the mask
    raw[~b & (np.arange(b.shape[-1]) % 2 == 0)] = -1
    raw[0, 0, 0] = np.nan
    assert np.array_equal(E.decon_mask_edge_erosion(raw, 3), E.decon_mask_edge_erosion(b, 3))
    u = b.astype(np.uint16) * np.uint16(7)
    assert np.array_equal(E.decon_mask_edge_erosion(u, 3), E.decon_mask_edge_erosion(b, 3))


def test_erosions_the_mex_does_not_take_run_imerode():
    """MATLAB's imerode code for erosions that aren't integers of at least 1, with strel's errors
    for strel('disk', EdgeErosion - 1) (an integer, then nonnegative)"""
    b = from_matlab(make_mask("skewed", [20, 22, 12]))
    for e, ident in [
        (2.5, "MATLAB:strel:expectedInteger"),
        (0, "MATLAB:strel:expectedNonnegative"),
        (math.nan, "MATLAB:strel:expectedInteger"),
        (-0.5, "MATLAB:strel:expectedInteger"),
    ]:
        with pytest.raises(E.MatlabError) as ei:
            E.decon_mask_edge_erosion(b, e)
        assert ei.value.identifier == ident


def test_non_double_erosion_values_are_taken_as_doubles():
    """DeconMaskErosionTest/nonDoubleErosionValuesAreTakenAsDoubles: single and integer erosions
    give the double's masks (the mex), also on a full mask longer than 255"""
    b = from_matlab(make_mask("skewed", [20, 22, 12]))
    full = np.ones((12, 300, 22), bool)
    for e in (np.float32(3), np.uint8(3), np.int16(3), True):
        d = float(e)
        assert np.array_equal(E.decon_mask_edge_erosion(b, e), E.decon_mask_edge_erosion(b, d))
        assert np.array_equal(
            E.decon_mask_edge_erosion(full, e), E.decon_mask_edge_erosion(full, d)
        )
    m = E.decon_mask_edge_erosion(full, np.uint8(3))
    assert m[3:-3, 3:-3, 3:-3].all() and m.sum() == 6 * 294 * 16


def test_a_full_mask_thinner_than_the_erosion_is_all_edge():
    """DeconMaskErosionTest/fullMasksThinnerThanTheErosionAreAllEdge: the edge indices stay in the
    mask, so a dimension shorter than the erosion is all edge (MATLAB errored before M07); a
    fractional erosion still errors (an index that isn't an integer); NaN takes index 1"""
    for shape, e in (
        ((3, 40, 44), 8),
        ((20, 30, 1), 2),
        ((10, 2, 20), 3),
        ((30,), 2),
        ((5, 5, 5), 6),
    ):
        a = np.ones(shape if len(shape) > 1 else (1, *shape), bool)
        assert not E.decon_mask_edge_erosion(a, e).any()
    a = np.ones((3, 40, 44), bool)
    a[0, 0, 0] = False
    assert not E.decon_mask_edge_erosion(a, 8).any()
    with pytest.raises(E.MatlabError) as ei:
        E.decon_mask_edge_erosion(np.ones((3, 20, 20), bool), 2.5)
    assert ei.value.identifier == "MATLAB:badsubscript" and "position 1" in str(ei.value)
    m = E.decon_mask_edge_erosion(np.ones((3, 20, 20), bool), math.nan)
    inner = np.ones((3, 20, 20), bool)
    inner[0], inner[:, 0], inner[:, :, 0] = False, False, False
    assert np.array_equal(m, inner)


def test_disk_decomposition_is_strels():
    nh, pu, pl = E.disk_decomposition(0)
    assert len(nh) == 1 and nh[0].shape == (1, 1) and pu == pl == (0, 0)
    nh, pu, pl = E.disk_decomposition(2)
    assert len(nh) == 1 and nh[0].shape == (5, 5)
    nh, pu, pl = E.disk_decomposition(10)
    assert len(nh) > 1 and all(v > 0 for v in pu + pl)


def test_results_do_not_depend_on_the_thread_count():
    b = from_matlab(make_mask("skewed", [60, 70, 40]))
    ref = E.decon_mask_edge_erosion(b, 9, threads=1)
    for t in (2, 5, 24):
        assert np.array_equal(E.decon_mask_edge_erosion(b, 9, threads=t), ref)
