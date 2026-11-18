"""Tests of petakit5d.stitch.xcorr: normxcorr2, crop_subregion_by_intensity and the template's
slices of the stitching registration (task S04).

Oracle cases (inputs in MATLAB's form):

- n2.*: MATLAB's normxcorr2 (Image Processing Toolbox) on integer and float images and
  templates of several sizes and classes, a template that is a part of the image (peak 1), equal
  sizes, rows and columns, data with a large offset, images with exactly zero parts. normxcorr2
  is MathWorks' code: the port's is a new implementation (Lewis's method), so the values are
  compared within the proposed tolerance of docs/parity/stitch.md (MATLAB's own result changes
  in the last bits with its thread count); err.n2.*: its errors (identifiers);
- crop.*: crop_subregion_by_intensity (StitchHelperMexTest's long region and its crop around the
  brightest part, regions within the limit, runs of blank slices, single data of 32 x 32 x 1000
  whose threshold MATLAB's threaded sums used to move): the crop and its box, bit for bit.

MATLAB tests ported here (StitchHelperMexTest, S04's): longRegionIsCroppedAroundItsBrightestPart
and cropDoesNotDependOnTheThreadCount (test_long_region_is_cropped_around_its_brightest_part,
test_crop_does_not_depend_on_the_thread_count, crop.*).
"""

from __future__ import annotations

import functools
from decimal import Decimal, getcontext

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.stitch import xcorr as sx
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value, Vol
from petakit5d.testing.compare import Custom, tol
from petakit5d.testing.synth import CounterRNG, texture

N2_REASON = (
    "docs/parity/stitch.md (proposed): normxcorr2, MathWorks' function, a new implementation"
)
# proposed: the coefficients within 1e-7 of MATLAB's. Measured: at most 2.6e-13 where the data
# have contrast; 9.2e-8 on 16-bit data at 2^15 with a contrast of 12, MATLAB's own error there
# (the port's is 1.6e-15 against the exact NCC). MATLAB's results also change by up to 7e-16
# with its thread count.
N2_TOL = tol(reason=N2_REASON, abs=1e-7)


def n2_offset_compare(ref, act) -> list[str]:
    """N2_TOL, but MATLAB's 0 where its rounding takes a coefficient beyond 1 (its true peak of 1
    on the offset data): the port keeps the NCC there"""
    r = np.asarray(ref.data, dtype=np.float64)
    a = np.asarray(act.data, dtype=np.float64)
    if r.shape != a.shape or ref.mclass != act.mclass:
        return [f"class or size: {ref.mclass} {r.shape}, {act.mclass} {a.shape}"]
    lost = (r == 0) & (np.abs(a) > 1 - 1e-7)
    d = np.abs(r - a)[~lost]
    out = [] if d.max(initial=0.0) <= 1e-7 else [f"max |difference| {d.max():.3g} > 1e-7"]
    if int(lost.sum()) != 1:
        out.append(f"{int(lost.sum())} coefficients of 1 that MATLAB sets to 0 (1 expected)")
    return out


# --- inputs (MATLAB's form: images (y, x)) ----------------------------------------------------


def ints(shape, lo, hi, stream, dtype=np.float64):
    return CounterRNG(4, f"n2/{stream}").integers(lo, hi, shape).astype(dtype)


@functools.cache
def n2_pair(name: str):
    """(template, image) of a normxcorr2 case, MATLAB's form"""
    kind, rest = name.split(":", 1)
    if kind == "int":
        my, mx, ny, nx = (int(v) for v in rest.split("x"))
        T = ints((my, mx), 0, 255, f"T{rest}")
        A = ints((ny, nx), 0, 255, f"A{rest}")
        if ny >= my + 3 and nx >= mx + 2:
            A[2 : 2 + my, 1 : 1 + mx] = T  # a peak of 1
        return T, A
    if kind == "tex":
        my, mx, ny, nx = (int(v) for v in rest.split("x"))
        im = texture((ny + 10, nx + 10), 2.0, (100.0, 2000.0), rng=CounterRNG(4, f"tex{rest}"))
        return im[5 : 5 + my, 7 : 7 + mx].copy(), im[:ny, :nx].copy()
    if kind == "offset":
        # 16-bit data at 2^15 with a contrast of a few counts
        im = ints((90, 80), -6, 6, "offset") + 2.0**15
        return im[20:50, 10:40].copy(), im
    if kind == "zeros":
        # an image with zero blocks (deskewed tiles' empty regions)
        T, A = n2_pair("int:12x10x70x90")
        A = A.copy()
        A[:, 40:70] = 0
        A[50:, :] = 0
        return T, A
    if kind == "classes":
        T, A = n2_pair("int:15x12x80x70")
        dt = {
            "uint8": np.uint8,
            "uint16": np.uint16,
            "single": np.float32,
            "int16": np.int16,
            "logical": np.bool_,
        }[rest]
        if dt == np.bool_:
            return T > 127, A > 127
        return T.astype(dt), A.astype(dt)
    raise KeyError(name)


N2_NAMES = [
    "int:15x12x80x70",
    "int:3x4x9x11",
    "int:40x7x60x90",
    "int:1x6x5x20",
    "int:6x1x30x4",
    "int:20x20x20x20",
    "int:2x1x7x5",
    "int:60x50x300x260",
    "tex:31x17x80x77",
    "tex:200x40x260x150",
    "offset:",
    "zeros:",
    "classes:uint8",
    "classes:uint16",
    "classes:single",
    "classes:int16",
    "classes:logical",
]


def n2_cases():
    out = []
    for name in N2_NAMES:
        make = functools.partial(n2_pair, name)
        cid = "n2." + name.replace(":", "_").rstrip("_")
        out.append(
            Case(
                cid,
                matlab="normxcorr2",
                args=[In("T"), In("A")],
                python="petakit5d.stitch.xcorr:normxcorr2",
                inputs={"T": Value(lambda m=make: m()[0]), "A": Value(lambda m=make: m()[1])},
                compare=(
                    Custom(n2_offset_compare, "N2_TOL; MATLAB's lost peak")
                    if name == "offset:"
                    else N2_TOL
                ),
            )
        )
    return out


def err_cases():
    rng = CounterRNG(4, "n2err")
    T = rng.uniform((4, 5))
    A = rng.uniform((20, 20))
    nan_t = T.copy()
    nan_t[2, 3] = np.nan
    inf_a = A.copy()
    inf_a[5, 5] = np.inf
    bad = {
        "const": (np.ones((5, 5)), A, "images:normxcorr2:sameElementsInTemplate"),
        "const0": (np.zeros((5, 5)), A, "images:normxcorr2:sameElementsInTemplate"),
        "bigrows": (rng.uniform((25, 5)), A, "images:normxcorr2:invalidSizeForA"),
        "bigcols": (rng.uniform((5, 25)), A, "images:normxcorr2:invalidSizeForA"),
        "nanT": (nan_t, A, "MATLAB:normxcorr2:expectedFinite"),
        "infA": (T, inf_a, "MATLAB:normxcorr2:expectedFinite"),
        "T3d": (rng.uniform((4, 5, 2)), A, "MATLAB:normxcorr2:expected2D"),
        "A3d": (T, rng.uniform((20, 20, 2)), "MATLAB:normxcorr2:expected2D"),
        "emptyT": (np.zeros((0, 3)), A, "images:normxcorr2:invalidTemplate"),
        "emptyA": (T, np.zeros((0, 20)), "images:normxcorr2:invalidSizeForA"),
        "scalarT": (np.array([[5.0]]), A, "images:normxcorr2:invalidTemplate"),
        "complexT": (T + 1j, A, "MATLAB:normxcorr2:expectedReal"),
    }
    out = []
    for k, (t, a, ident) in bad.items():
        out.append(
            Case(
                f"err.n2.{k}",
                matlab="normxcorr2",
                args=[In("T"), In("A")],
                python="petakit5d.stitch.xcorr:normxcorr2",
                inputs={"T": Value(lambda t=t: t), "A": Value(lambda a=a: a)},
                expect_error=ident,
                py_error=sx.RegistrationError,
            )
        )
    return out


# --- crop_subregion_by_intensity --------------------------------------------------------------


def long_region():
    """StitchHelperMexTest.longRegionIsCroppedAroundItsBrightestPart: 100 + randi(20) of 30 x 40 x
    12 (MATLAB's [y x z]), 5000 in (20:22, 31:33, 5:7)"""
    r = (100 + CounterRNG(4, "long").integers(1, 20, (12, 30, 40))).astype(np.uint16)
    r[4:7, 19:22, 30:33] = 5000
    return to_matlab(r)


def blank_region():
    """a region of 12 x 300 x 20 ([y x z]) whose bright part along x lies between blank runs"""
    r = np.zeros((20, 12, 300), dtype=np.uint16)
    r[:, :, 40:90] = 100 + CounterRNG(4, "blank").integers(0, 20, (20, 12, 50))
    r[:, :, 150:260] = 100 + CounterRNG(4, "blank2").integers(0, 20, (20, 12, 110))
    r[5:9, 3:7, 200:204] = 4000
    r[2:4, 3:5, 60:62] = 3000
    return to_matlab(r)


def thread_region():
    """StitchHelperMexTest.cropDoesNotDependOnTheThreadCount: single 6e7 + 4 randi([-4, 4]) of
    32 x 32 x 1000, slices of bright voxels of one value each, increasing along z"""
    nz = 1000
    rng = CounterRNG(4, "threads")
    r = (6e7 + 4 * rng.integers(-4, 4, (nz, 32, 32))).astype(np.float32)
    for z in range(nz):
        k = 10 + (nz - 1 - z) // 50
        idx = rng.integers(0, 32 * 32 - 1, (k,))
        flat = r[z].T.reshape(-1)  # MATLAB's order within the slice
        flat[np.unique(idx)] = np.float32(6e7 + 4 * (z + 1))
        r[z] = flat.reshape(32, 32).T
    return to_matlab(r)


CROPS = {
    "long10": (long_region, 10, None),
    "long40": (long_region, 40, None),
    "long25": (long_region, 25, None),
    "blank": (blank_region, 60, 20),
    "blank_none": (blank_region, 60, None),
    "blank_wide": (blank_region, 250, 5),
    "threads": (thread_region, 100, None),
    "threads_blank": (thread_region, 100, 3),
}


def crop_cases():
    out = []
    for name, (make, thr, blank) in CROPS.items():
        args = [Vol(In("r")), thr] + ([] if blank is None else [blank])
        out.append(
            Case(
                f"crop.{name}",
                matlab="crop_subregion_by_intensity",
                args=args,
                nargout=2,
                python="petakit5d.stitch.xcorr:crop_subregion_by_intensity",
                inputs={"r": Value(make)},
                out_roles=("vol", "box"),
            )
        )
    return out


CASES = n2_cases() + err_cases() + crop_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# --- normxcorr2 against the exact NCC -----------------------------------------------------------


def exact_ncc(T, A):
    """the exact NCC of integer images (Python integers, the quotient in 40 digits); 0 for flat
    windows"""
    getcontext().prec = 40
    T = T.astype(np.int64)
    A = A.astype(np.int64)
    m, n = T.shape
    M, N = A.shape
    P = np.zeros((M + 2 * m, N + 2 * n), dtype=np.int64)
    P[m : m + M, n : n + N] = A
    nt = m * n
    st, st2 = int(T.sum()), int((T * T).sum())
    dt = nt * st2 - st * st
    out = np.zeros((M + m - 1, N + n - 1))
    for k in range(M + m - 1):
        for j in range(N + n - 1):
            W = P[k + 1 : k + 1 + m, j + 1 : j + 1 + n]
            sa, sa2, sta = int(W.sum()), int((W * W).sum()), int((W * T).sum())
            da = nt * sa2 - sa * sa
            if da:
                num = nt * sta - sa * st
                out[k, j] = float(Decimal(num) / (Decimal(dt) * Decimal(da)).sqrt())
    return out


@pytest.mark.parametrize(
    "name", ["int:15x12x80x70", "int:3x4x9x11", "int:1x6x5x20", "int:6x1x30x4", "zeros:"]
)
def test_normxcorr2_is_the_ncc(name):
    T, A = n2_pair(name)
    C = sx.normxcorr2(T, A)
    E = exact_ncc(T, A)
    assert C.shape == (A.shape[0] + T.shape[0] - 1, A.shape[1] + T.shape[1] - 1)
    assert C.dtype == np.float64
    assert np.abs(C - E).max() < 1e-13
    # the template is a part of the image: its peak is there
    if name.startswith("int:15"):
        k = int(np.argmax(C.T.reshape(-1)))
        assert (k % C.shape[0], k // C.shape[0]) == (2 + T.shape[0] - 1, 1 + T.shape[1] - 1)


def test_normxcorr2_takes_any_class_as_double():
    T, A = n2_pair("int:15x12x80x70")
    C = sx.normxcorr2(T, A)
    for dt in (np.uint8, np.uint16, np.float32, np.int16, np.int32):
        np.testing.assert_array_equal(sx.normxcorr2(T.astype(dt), A.astype(dt)), C)
    assert sx.normxcorr2(T > 100, A > 100).dtype == np.float64


def test_normxcorr2_flat_windows_are_zero():
    """windows of an exactly flat part of the image (any value) give 0; the textured part its NCC
    (MATLAB gives rounding noise up to about 1e-7 in flat windows: docs/parity/stitch.md)"""
    T, _ = n2_pair("int:9x8x50x120")
    rng = CounterRNG(4, "flat")
    for c in (0.0, 1.0, 500.0, 1e8):
        A = np.full((50, 120), c)
        A[:, :40] = rng.uniform((50, 40)) * max(c, 1.0)
        C = sx.normxcorr2(T, A)
        assert np.all(C[8:50, 40 + 7 : 120] == 0)
        assert np.abs(C[8:50, 7:40]).max() > 0.05


def test_normxcorr2_offset_costs_no_precision():
    T, A = n2_pair("offset:")
    C = sx.normxcorr2(T, A)
    E = exact_ncc(T - 2.0**15, A - 2.0**15)
    # (the exact NCC of the data without the offset, which only the padding sees)
    E2 = exact_ncc(T, A)
    assert np.abs(C - E2).max() < 1e-12
    assert E.shape == C.shape


def test_normxcorr2_errors_have_matlabs_identifiers():
    import scipy.sparse as sp

    T = np.arange(20.0).reshape(4, 5)
    with pytest.raises(sx.RegistrationError) as e:
        sx.normxcorr2(sp.csr_matrix(T), np.ones((9, 9)))
    assert e.value.identifier == "MATLAB:normxcorr2:expectedNonsparse"
    with pytest.raises(sx.RegistrationError) as e:
        sx.normxcorr2(np.array([[{}]], dtype=object), np.ones((9, 9)))
    assert e.value.identifier == "MATLAB:normxcorr2:invalidType"


def test_normxcorr2_max_shift_finds_the_offset():
    im = texture((120, 140), 2.0, (100.0, 2000.0), rng=CounterRNG(4, "ms"))
    T = im[30:70, 40:90].copy()
    A = im[10:100, 20:120].copy()
    off, mx, C = sx.normxcorr2_max_shift(T.astype(np.float32), A.astype(np.float32), (0, 30, 30))
    assert off == (0, 20, 20)
    assert mx == pytest.approx(1.0, abs=1e-12)
    assert C.shape == (61, 61)


# --- crop_subregion_by_intensity ------------------------------------------------------------


def test_long_region_is_cropped_around_its_brightest_part():
    region = from_matlab(long_region())
    cropped, box = sx.crop_subregion_by_intensity(region, 10)
    # MATLAB's [16, 27, 1, 25, 36, 10]
    assert box == (0, 15, 26, 10, 25, 36)
    np.testing.assert_array_equal(cropped, region[0:10, 15:25, 26:36])
    kept, box = sx.crop_subregion_by_intensity(region, 40)
    assert kept is region
    assert box == (0, 0, 0, 12, 30, 40)


def test_crop_does_not_depend_on_the_thread_count():
    region = from_matlab(thread_region())
    ref, ref_box = sx.crop_subregion_by_intensity(region, 100, num_threads=1)
    assert ref_box[0] > 0 and ref_box[3] < region.shape[0]
    assert ref_box[1:3] == (0, 0) and ref_box[4:] == (32, 32)
    for t in (2, 4, 24):
        cropped, box = sx.crop_subregion_by_intensity(region, 100, num_threads=t)
        assert box == ref_box
        np.testing.assert_array_equal(cropped, ref)


def test_blank_runs_bound_the_crop():
    """the blank runs are those of the counts of bright voxels (above mean + 3 std): around the
    brightest x (200 .. 203, 0-based) the runs without bright voxels end at 199 and start at 204
    (1-based 200 and 205), which bound the crop; without blankNumTrsh it is 60 wide"""
    region = from_matlab(blank_region())
    _, box = sx.crop_subregion_by_intensity(region, 60, 20)
    assert box == (0, 0, 199, 20, 12, 205)
    _, box = sx.crop_subregion_by_intensity(region, 60)
    assert box[5] - box[2] == 60 and box[2] <= 200 < 204 <= box[5]


def test_one_slice_regions_index_past_the_box():
    """MATLAB's size() of one slice has 2 values: the crop box has 5 and indexing its 6th fails"""
    region = np.ones((1, 30, 20000), dtype=np.uint16)
    with pytest.raises(sx.RegistrationError) as e:
        sx.crop_subregion_by_intensity(region, 10000)
    assert e.value.identifier == "MATLAB:badsubscript"


# --- the template's slices ---------------------------------------------------------------------


def test_cap_region_inds_cuts_the_margins():
    # 6 or more positions: those inside (max(n b1, min(n b2, m / 2)), min(n (1 - b1), max(...)))
    pos = sx.cap_region_inds(np.ones(40, dtype=bool), 40, 10)
    assert pos[0] == 6 and pos[-1] == 34  # (5, 35): 0.25 * 40 > 10 / 2
    pos = sx.cap_region_inds(np.ones(40, dtype=bool), 40, 30)
    assert pos[0] == 11 and pos[-1] == 29  # (10, 30)
    # five or fewer are kept
    m = np.zeros(40, dtype=bool)
    m[[0, 3, 7, 20, 39]] = True
    np.testing.assert_array_equal(sx.cap_region_inds(m, 40, 10), [1, 4, 8, 21, 40])


def test_template_indices_count_within_the_chosen_slices():
    nz = np.zeros((4, 6, 8), dtype=bool)
    nz[1:3, :, 2:7] = True  # slices 2, 3 with data in 30 of 48
    nz[0, 0, 0] = True
    y, x, z = sx.template_indices(nz, 4)
    np.testing.assert_array_equal(z, [False, True, True, False])
    np.testing.assert_array_equal(x, [False, False, True, True, True, True, True, False])
    assert y.all()
    y2, x2, z2 = sx.template_indices(nz, 4, two_d=True)
    assert z2[0] and not z2[1:].any()
    assert not x2.any() and not y2.any()
