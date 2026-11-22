"""Channel offsets in x and z (the lab's XR_correctXZoffsetData3D and XZoffsetEstimation3D):
oracle cases on synthetic deskewed frames, and unit tests of the port's rules (task C04a)."""

from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.arrays import to_matlab
from petakit5d.detection import condition as cond
from petakit5d.detection import offsets
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, In, InPath, Tree, Value, Vol, mcell
from petakit5d.testing.compare import Masks, tol
from petakit5d.testing.compare import tree as tree_spec

DS_SIZE = (14, 10, 36)  # (z, y, x)
# xz_offset.mat holds the channels' folders, which the port makes with MATLAB's filesep
# ('\\' on Windows, as MATLAB there); the references were made on Linux: the files' texts are
# compared with '/' for every separator
SLASHES = ((r"\\", "/"),)
SPOTS = (
    (4.2, 3.1, 9.4), (7.6, 6.2, 14.8), (10.1, 2.7, 21.3), (5.5, 7.4, 26.6), (8.8, 4.9, 18.2),
    (11.4, 6.6, 12.1), (3.7, 5.3, 29.5),
)  # fmt: skip
W1 = "detection.md W1"


def ds_volume(shift=(0, 0, 0), stream="ch1", amp=300.0, size=DS_SIZE) -> np.ndarray:
    """a deskewed-like frame (z, y, x), single: spots shifted by (z, y, x) voxels on a noisy
    background, zero outside the deskewed parallelogram"""
    pos = np.array(SPOTS, dtype=np.float64) + np.array(shift, dtype=np.float64)
    vol = synth.beads(size, pos, (1.4, 1.2), np.full(len(SPOTS), amp), 100.0)
    noise = synth.CounterRNG(2026, f"offsets/{stream}").normal(size)
    v = vol + 3.0 * noise
    nz, _, nx = size
    z = np.arange(nz)[:, None, None]
    x = np.arange(nx)[None, None, :]
    inside = (x >= z) & (x <= nx - 1 - (nz - 1 - z))
    return np.where(inside, v, 0.0).astype(np.float32)


CHANNELS = {"ch1": (0, 0, 0), "ch2": (-1, 0, 2), "ch3": (1, 0, -3)}


def _bump(peaks, size=(15, 11)) -> np.ndarray:
    """a 15 x 11 matrix (x shifts, z shifts) with exact values: 1 - distance^2 / 64 to the
    nearest of the given (row, column) peaks, each peak lifted by its height"""
    i, j = np.meshgrid(np.arange(size[0]), np.arange(size[1]), indexing="ij")
    m = np.full(size, -1.0)
    for (r, c), h in peaks:
        m = np.maximum(m, h - ((i - r) ** 2 + (j - c) ** 2) / 64)
    return m


PICKS = {
    "one": _bump([((9, 3), 1.0)]),
    "tie_in_x": _bump([((6, 5), 1.0), ((8, 5), 1.0)]),
    "tie_in_z": _bump([((7, 4), 1.0), ((7, 6), 1.0)]),
    "four_equal": _bump([((6, 4), 1.0), ((8, 4), 1.0), ((6, 6), 1.0), ((8, 6), 1.0)]),
    "nearer_lower": _bump([((3, 2), 1.0), ((8, 6), 0.5)]),
    "border_peak": _bump([((0, 5), 1.0), ((10, 8), 0.4)]),
    "only_border": _bump([((14, 0), 1.0)]),
    "plateau": np.where(_bump([((7, 5), 1.0)]) > 0.95, 1.0, _bump([((7, 5), 1.0)])),
    "flat": np.full((15, 11), 0.25),
    "far_corner_interior": _bump([((13, 9), 1.0), ((1, 1), 1.0)]),
    # one maximum of two diagonal pixels: equally near MATLAB's middle, not (8, 6)'s
    "diagonal_pair": _bump([((6, 4), 1.0), ((7, 5), 1.0)]),
    # equally near the middle: the first in column-major order, not in row-major order
    "antidiagonal_pair": _bump([((6, 5), 1.0), ((7, 4), 1.0)]),
    # maxima on each of the four sides only: none counts
    "border_sides": _bump([((0, 5), 1.0), ((14, 5), 1.0), ((7, 0), 1.0), ((7, 10), 1.0)]),
}


def write_condition(root: Path, chs, n_frames: int, shifts=None, names=None) -> None:
    shifts = shifts or CHANNELS
    movie = root / "Ex1_syn_z0p5"
    for ch in chs:
        for d in (movie / ch, movie / ch / "DS"):
            d.mkdir(parents=True, exist_ok=True)
        for k in range(n_frames):
            name = (names or {}).get(ch, "syn") + f"_stack{k:04d}_{2000 * k:07d}msec.tif"
            s = shifts[ch]  # one shift, or one per frame
            v = ds_volume(s[k] if isinstance(s, list) else s, f"{ch}/{k}")
            synth.write_tiff(movie / ch / name, v)
            synth.write_tiff(movie / ch / "DS" / name, v)


def make_tree(folder) -> None:
    r = Path(folder)
    write_condition(r / "two", ["ch1", "ch2"], 3)
    write_condition(r / "same", ["ch1", "ch2"], 3, shifts={"ch1": (0, 0, 0), "ch2": (0, 0, 0)})
    write_condition(r / "three", ["ch1", "ch2", "ch3"], 3)
    write_condition(r / "twoframes", ["ch1", "ch2"], 2)
    write_condition(r / "five", ["ch1", "ch2"], 5)
    write_condition(r / "oneframe3", ["ch1", "ch2", "ch3"], 1)
    write_condition(r / "flag", ["ch1", "ch2"], 3)
    synth.write_tiff(r / "flag" / "Ex1_syn_z0p5" / "ch1" / "DS" / "keep.tif", ds_volume())
    open(r / "flag" / "Ex1_syn_z0p5" / "ch1" / "DS" / "xz_offset.mat", "w").close()
    write_condition(r / "lock", ["ch1", "ch2"], 3)
    open(r / "lock" / "Ex1_syn_z0p5" / "ch1" / "DS" / "xz_offset.tmp", "w").close()
    write_condition(r / "tiff", ["ch1", "ch2"], 3)
    for ch in ("ch1", "ch2"):
        d = r / "tiff" / "Ex1_syn_z0p5" / ch
        for p in sorted(d.glob("*.tif")) + sorted((d / "DS").glob("*.tif")):
            p.rename(p.with_suffix(".tiff"))


def make_halves(folder) -> None:
    """two frames whose shifts differ by one voxel in x and z: medians of x.5, of both signs"""
    shifts = {
        "ch1": (0, 0, 0),
        "ch2": [(-1, 0, 2), (-2, 0, 3)],
        "ch3": [(1, 0, -2), (2, 0, -3)],
    }
    write_condition(Path(folder) / "halves", ["ch1", "ch2", "ch3"], 2, shifts=shifts)


CORRECT = (
    "d = XR_loadConditionData3D([in{1} in{2}], in{3}, in{4}); "
    "XR_correctXZoffsetData3D(d, in{5}{:}); out = {};"
)
MARKERS = ("gfp", "rfp", "mcherry")


def _correct(root, rel, chs, mks, opts):
    names = {"Overwrite": "overwrite", "DS": "ds", "sCMOSCameraFlip": "s_cmos_camera_flip"}
    kw = {names[opts[i]]: opts[i + 1] for i in range(0, len(opts), 2)}
    d = cond.load_condition_data(root + rel, list(chs), list(mks))
    offsets.correct_xz_offset_data(d, **kw)


def _cases():
    tree = {"t": Tree(make_tree, "t")}
    cs = []

    def add(cid, rel, chs=("ch1", "ch2"), opts=(), error=None):
        kw = {}
        if error:
            kw = {"expect_error": error, "py_error": (offsets.OffsetError, cond.ConditionError)}
        cs.append(
            Case(
                f"correct.{cid}",
                code=CORRECT,
                args=[InPath("t"), rel, mcell(*chs), mcell(*MARKERS[: len(chs)]), mcell(*opts)],
                py_args=[InPath("t"), rel, list(chs), list(MARKERS[: len(chs)]), list(opts)],
                python=_correct,
                nargout=0,
                inputs=tree,
                tree=tree_spec(masks=SLASHES),
                **kw,
            )
        )

    add("two_channels", "two/Ex1_syn_z0p5/")
    add("no_offset", "same/Ex1_syn_z0p5/")
    add("three_channels", "three/Ex1_syn_z0p5/", ("ch1", "ch2", "ch3"))
    add("two_frames", "twoframes/Ex1_syn_z0p5/")
    add("five_frames", "five/Ex1_syn_z0p5/")
    add("unsorted_channels", "two/Ex1_syn_z0p5/", ("ch2", "ch1"))
    add("flag_skips", "flag/Ex1_syn_z0p5/")
    add("flag_overwrite", "flag/Ex1_syn_z0p5/", opts=("Overwrite", True))
    add("lock_skips", "lock/Ex1_syn_z0p5/")
    add("raw_frames", "two/Ex1_syn_z0p5/", opts=("DS", False))
    # frames that are all .tiff: MATLAB's 3D loader fails (measured)
    add("tiff_names", "tiff/Ex1_syn_z0p5/", error="MATLAB:badsubscript")
    # medians of two frames of x.5 and -x.5 (fix, not floor): its own tree
    cs.append(
        Case(
            "correct.two_frames_half_medians",
            code=CORRECT,
            args=[InPath("h"), "halves/Ex1_syn_z0p5/", mcell("ch1", "ch2", "ch3"),
                  mcell(*MARKERS), mcell()],
            py_args=[InPath("h"), "halves/Ex1_syn_z0p5/", ["ch1", "ch2", "ch3"], list(MARKERS), []],
            python=_correct,
            nargout=0,
            inputs={"h": Tree(make_halves, "h")},
            tree=tree_spec(masks=SLASHES),
        )
    )  # fmt: skip
    add("one_frame_three_channels", "oneframe3/Ex1_syn_z0p5/", ("ch1", "ch2", "ch3"),
        error="MATLAB:badsubscript")  # fmt: skip
    # XZoffsetEstimation3D on volumes in memory (MATLAB's (y, x, z), NaN for no data)
    pairs = {
        "shift": (ds_volume(), ds_volume((-1, 0, 2), "b")),
        "same": (ds_volume(), ds_volume((0, 0, 0), "b")),
        "far": (ds_volume(), ds_volume((3, 0, -5), "b")),
    }
    for name, (a, b) in pairs.items():

        def mk(v):
            m = to_matlab(v.astype(np.float64))
            m[m == 0] = np.nan
            return m

        inputs = {"a": Value(lambda a=a: mk(a)), "b": Value(lambda b=b: mk(b))}
        cs.append(
            Case(
                f"estimate.{name}",
                matlab="XZoffsetEstimation3D",
                args=[Vol(In("a")), Vol(In("b"))],
                python=offsets.xz_offset_estimation_3d,
                py_out=lambda r: [float(r[0]), float(r[1])],
                nargout=2,
                inputs=inputs,
            )
        )
        # the matrix it maximizes, as its loop computes it with corr (W1)
        code = (
            "f1 = in{1}; f2 = in{2}; [~, nx, nz] = size(f1); xs = -7 : 7; zs = -5 : 5; "
            "c = zeros(15, 11); for i = 1 : 15, xo = xs(i); if xo <= 0, a1 = 1 : nx + xo; "
            "a2 = 1 - xo : nx; else, a1 = 1 + xo : nx; a2 = 1 : nx - xo; end; "
            "for j = 1 : 11, zo = zs(j); if zo <= 0, b1 = 1 : nz + zo; b2 = 1 - zo : nz; else, "
            "b1 = 1 + zo : nz; b2 = 1 : nz - zo; end; g1 = f1(:, a1, b1); g2 = f2(:, a2, b2); "
            "k = ~isnan(g1 + g2); c(i, j) = corr(g1(k), g2(k)); end; end; out = {c};"
        )
        cs.append(
            Case(
                f"correlations.{name}",
                code=code,
                args=[Vol(In("a")), Vol(In("b"))],
                python=offsets.xz_correlations,
                inputs=inputs,
                compare=tol(reason=W1, abs=4e-15),
            )
        )
    # the pick from a matrix: XZoffsetEstimation3D's lines after its loop, on matrices made
    # to have several maxima, ties, plateaus and border maxima
    pick = (
        "c = in{1}; BW = imregionalmax(c); BW(1 : end, [1, end]) = false; "
        "BW([1, end], 1 : end) = false; [lmx, lmz] = find(BW); nx = size(c, 1); nz = size(c, 2); "
        "if numel(lmz) > 1, [~, k] = min((lmz - nz / 2) .^ 2 + (lmx - nx / 2) .^ 2); "
        "lmx = lmx(k); lmz = lmz(k); end; xs = -7 : 7; zs = -5 : 5; "
        "if ~isempty(lmz), out = {xs(lmx), zs(lmz)}; else, out = {0, 0}; end"
    )
    for name, mat in PICKS.items():
        cs.append(
            Case(
                f"pick.{name}",
                code=pick,
                args=[In("c")],
                python=offsets.pick_offset,
                py_out=lambda r: [float(r[0]), float(r[1])],
                nargout=2,
                inputs={"c": Value(lambda mat=mat: mat)},
            )
        )
    cs.append(
        Case(
            "estimate.third_argument",
            matlab="XZoffsetEstimation3D",
            args=[np.ones((4, 9, 7)), np.ones((4, 9, 7)), [-1, 0, 1]],
            python=offsets.xz_offset_estimation_3d,
            nargout=2,
            expect_error="MATLAB:UndefinedFunction",
            py_error=(offsets.OffsetError,),
        )
    )
    cs.append(
        Case(
            "estimate.constant",
            matlab="XZoffsetEstimation3D",
            args=[Vol(np.ones((6, 20, 12))), Vol(np.ones((6, 20, 12)))],
            python=offsets.xz_offset_estimation_3d,
            nargout=2,
            expect_error="MATLAB:imregionalmax:expectedNonNaN",
            py_error=(offsets.OffsetError,),
        )
    )
    return cs


CASES = _cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_channel_paths_compare_with_slashes():
    """xz_offset.mat's channels as the port writes them on Windows mask to the Linux texts"""
    m = Masks().plus(*SLASHES)
    assert m.text("<WORK>/t/two/Ex1_syn_z0p5/ch1\\") == "<WORK>/t/two/Ex1_syn_z0p5/ch1/"
    assert m.text("C:\\w\\t/two/ch2\\") == "C:/w/t/two/ch2/"
    assert m.text("<WORK>/t/two/Ex1_syn_z0p5/ch1/") == "<WORK>/t/two/Ex1_syn_z0p5/ch1/"


def test_ranges_are_matlabs():
    """the index ranges of XZoffsetEstimation3D's shifts, 0-based"""
    assert offsets._ranges(-2, 10) == (slice(0, 8), slice(2, 10))
    assert offsets._ranges(0, 10) == (slice(0, 10), slice(0, 10))
    assert offsets._ranges(3, 10) == (slice(3, 10), slice(0, 7))
    assert offsets._ranges(-12, 10) == (slice(0, 0), slice(10, 10))


def test_nearest_of_several_maxima_leans_negative():
    """of two regional maxima at the same distance from (8, 6), MATLAB's 7.5 and 5.5 pick the
    lower index: x offset -1 rather than +1"""
    c = np.zeros((15, 11))
    c[6, 5] = 1.0  # x -1, z 0
    c[8, 5] = 1.0  # x +1, z 0
    k = np.nonzero(c.T)
    d = (k[0] + 1 - 5.5) ** 2 + (k[1] + 1 - 7.5) ** 2
    assert offsets.XZOFF_SHIFTS[int(k[1][np.argmin(d)])] == -1


def test_median_of_one_frame_runs_along_the_row():
    m = np.array([[3.0, -2.0, 5.0]])
    assert offsets._median_rows(m).tolist() == [3.0]
    m = np.array([[3.0, -2.0], [1.0, 4.0], [2.0, 0.0]])
    assert offsets._median_rows(m).tolist() == [2.0, 0.0]


def test_lock_clean_removes_old_tmp_only(tmp_path):
    import os
    import time

    old = tmp_path / "a.tmp"
    new = tmp_path / "sub" / "b.tmp"
    new.parent.mkdir()
    old.write_text("x")
    new.write_text("x")
    t = time.time() - 3600
    os.utime(old, (t, t))
    offsets._lock_clean(str(tmp_path), 30)
    assert not old.exists() and new.exists()


@pytest.mark.parametrize(
    ("flag", "name"),
    [("overwrite", "Overwrite"), ("ds", "DS"), ("s_cmos_camera_flip", "sCMOSCameraFlip")],
)
def test_flags_must_be_logical(flag, name):
    """the registry checks the lab's validators (@islogical), with inputParser's identifier"""
    from petakit5d.params import MatlabArgumentError

    with pytest.raises(MatlabArgumentError, match=f"'{name}'") as e:
        offsets.correct_xz_offset_data(cond.Condition([]), **{flag: 1})
    assert e.value.identifier == "MATLAB:InputParser:ArgumentFailedValidation"
