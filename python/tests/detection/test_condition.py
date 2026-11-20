"""Loading a condition (MATLAB's XR_loadConditionData3D and XR_loadConditionData): oracle cases on
synthetic folder trees, and the rules measured by black-box probes (task C04a).

The ``date`` field is the digits of the folder path, which differ between MATLAB's and Python's
work folders: the oracle cases compare every other field (MATLAB's side drops it in its code,
the Python side in ``_matlab``) and the rule is checked here on MATLAB's probed answers.

The paths: the loader adds MATLAB's ``filesep``, ``\\`` on Windows as MATLAB does there, to the
folders it is given; the references were made on Linux. The cases compare path texts with
``/`` for every separator (``posix_paths``).
"""

import os

import numpy as np
import pytest

from petakit5d.detection import condition as cond
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, Data, InPath, Tree, mcell

DEMO = "$PETAKIT5D_DEMO_DATA/PetaKit5D_point_detection_demo_dataset"
# the dataset's condition folder (its only folder) and its first experiment, found as MATLAB's
# dir lists them
PATH_FIELDS = ("source", "channels", "framePaths", "maskPaths", "framePathsDS", "framePathsDSR")
FIND = (
    "PATHS = {'source', 'channels', 'framePaths', 'maskPaths', 'framePathsDS', 'framePathsDSR'}; "
    "r = in{1}; if r(end) ~= '/', r = [r '/']; end; f = dir(r); "
    "f = f([f.isdir] & ~startsWith({f.name}, '.')); c = [f(1).name '/']; "
    "e = dir([r c 'Ex1*']); rel = {c, [c e(1).name '/']}; rel = rel{in{2}}; "
)
CH = ["ch1", "ch2"]
MK = ["gfp", "rfp"]
SIZE = (3, 6, 8)  # (z, y, x) of a frame


def frame(k: int, c: int, size=SIZE) -> np.ndarray:
    n = int(np.prod(size))
    v = (np.arange(n, dtype=np.int64) * 7 + 11 * k + 101 * c) % 4000
    return v.reshape(size).astype(np.uint16)


def std_names(n: int, dt: int = 2000, t0: int = 0) -> list[str]:
    t = [t0 + k * dt for k in range(n)]
    return [
        f"syn_CamA_ch0_stack{k:04d}_560nm_{t[k]:07d}msec_{1000000 + t[k]:010d}msecAbs.tif"
        for k in range(n)
    ]


def write_movie(movie, chs, names, size=SIZE) -> None:
    per = names if isinstance(names[0], list | tuple) else [names] * len(chs)
    for c, (ch, nm) in enumerate(zip(chs, per, strict=True), start=1):
        d = movie / ch
        d.mkdir(parents=True, exist_ok=True)
        for k, name in enumerate(nm):
            synth.write_tiff(d / name, frame(k, c, size))


def make_tree(folder) -> None:
    """every layout of the cases below, each in its own subfolder"""
    from pathlib import Path

    r = Path(folder)
    write_movie(r / "one" / "Ex1_syn_z0p5", CH, std_names(5))
    for nm in ("Ex2_b_z0p5", "Ex1_a_z0p5", "Ex10_c_z0p5", "xEx3_z0p5", "other_z0p5"):
        write_movie(r / "multi" / "cond" / nm, CH, std_names(3))
    for rel in ("Ex3", "Ex1", "sub/Ex2", "sub/Ex0", "sub/deep/Ex5"):
        write_movie(r / "depths" / "cond" / rel, CH, std_names(2))
    names = [
        "b_CamA_ch0_stack0001_560nm_0004000msec_0000000100msecAbs.tif",
        "a_CamA_ch0_stack0002_560nm_0002000msec_0000000300msecAbs.tif",
        "c_CamA_ch0_stack0000_560nm_0000000msec_0000000200msecAbs.tif",
        "f_stack10_0006000msec.tif",
    ]
    write_movie(r / "order" / "Ex1_z0p5", CH, names)
    t = [0, 2000, 5000, 7000, 12000]
    names = [f"r_stack{k:04d}_{t[k]:07d}msec_{10 + 1000 * k:010d}msecAbs.tif" for k in range(5)]
    write_movie(r / "rate" / "Ex1_z0p5", CH, names)
    for rel in (
        "dz/Ex1_z0p25",
        "dz/Ex1_z1p5",
        "dz/Ex1_Z0P5",
        "dz/cond_z0p3/Ex1",
        "dz/q/Ex1_z0p5_z0p25",
    ):
        write_movie(r / rel, CH, std_names(2))
    write_movie(r / "dz" / "chan" / "Ex1", ["ch1_z0p7", "ch2_z0p7"], std_names(2))
    write_movie(r / "more" / "Ex1_z0p5", CH, [std_names(3), std_names(4)])
    write_movie(r / "fewer" / "Ex1_z0p5", CH, [std_names(4), std_names(3)])
    ext = [
        "e_stack0000_0000000msec.tif",
        "e_stack0001_0002000msec.tiff",
        "e_stack0002_0004000msec.TIF",
    ]
    write_movie(r / "ext" / "Ex1_z0p5", CH, ext)
    write_movie(r / "ext" / "Ex1_z0p5", ["ch1"], [["._e_stack0000_0000000msec.tif"]])
    with open(r / "ext" / "Ex1_z0p5" / "ch1" / "notes.txt", "w", newline="\n") as f:
        f.write("x\n")
    (r / "ext" / "Ex1_z0p5" / "ch1" / "DS").mkdir()
    synth.write_tiff(r / "ext" / "Ex1_z0p5" / "ch1" / "DS" / ext[0], frame(0, 9))
    for nm in ("cell2_z0p5", "cell1_z0p5"):
        write_movie(r / "cells" / nm, CH, std_names(2))
    write_movie(r / "twod" / "Ex1_z0p5", CH, std_names(2), size=(1, 6, 8))
    write_movie(r / "nomsec" / "Ex1_z0p5", CH, ["frame_0001.tif", "frame_0002.tif"])
    write_movie(r / "ms" / "Ex1_40ms_z0p5", CH, std_names(10))
    (r / "empty" / "Ex1_z0p5" / "ch1").mkdir(parents=True)
    (r / "empty" / "Ex1_z0p5" / "ch2").mkdir(parents=True)
    (r / "nothing" / "q").mkdir(parents=True)
    for nm in (
        "Ex12a3_z0p5",
        "Ex2_z0p5",
        "Ex1_b_z0p5",
        "Ex1_a_z0p5",
        "exp3_z0p5",
        "Ex 4_z0p5",
        "Ex05_z0p5",
    ):
        write_movie(r / "nonum" / "cond" / nm, CH, std_names(2))


LOAD = (
    "d = XR_loadConditionData3D([in{1} in{2}], in{3}, in{4}, in{5}{:}); out = {rmfield(d, 'date')};"
)
LOAD2D = (
    "d = XR_loadConditionData([in{1} in{2}], in{3}, in{4}, in{5}{:}); out = {rmfield(d, 'date')};"
)


def posix_paths(v):
    """the value with ``/`` for every ``\\`` in its texts (lists, tuples, dicts, object arrays):
    a path the port made with Windows' filesep as the Linux references hold it"""
    if isinstance(v, str):
        return v.replace("\\", "/")
    if isinstance(v, dict):
        return {k: posix_paths(x) for k, x in v.items()}
    if isinstance(v, list | tuple):
        return type(v)(posix_paths(x) for x in v)
    if isinstance(v, np.ndarray) and v.dtype == object:
        out = np.empty(v.shape, dtype=object)
        for i, x in np.ndenumerate(v):
            out[i] = posix_paths(x)
        return out
    return v


def _no_date(movies) -> list[dict]:
    out = []
    for m in movies:
        d = posix_paths(m.to_matlab())
        del d["date"]
        out.append(d)
    return out


def _options(opts: list) -> dict:
    names = {
        "MovieSelector": "movie_selector",
        "IgnoreExpNumber": "ignore_exp_number",
        "ChannelPatterns": "channel_patterns",
        "Angle": "angle",
        "PixelSize": "pixel_size",
        "dz": "dz",
        "ObjectiveScan": "objective_scan",
        "Reversed": "reversed",
        "FrameRate": "frame_rate",
        "maxlevel": "maxlevel",
        "Parameters": "parameters",
    }
    return {names[opts[i]]: opts[i + 1] for i in range(0, len(opts), 2)}


def _load3d(root, rel, chs, mks, opts):
    return _no_date(cond.load_condition_data(root + rel, list(chs), list(mks), **_options(opts)))


def _load2d(root, rel, chs, mks, opts):
    return _no_date(cond.load_condition_data_2d(root + rel, list(chs), list(mks), **_options(opts)))


def _cases():
    tree = {"t": Tree(make_tree, "t")}
    cs = []

    def add(cid, rel, chs=CH, mks=MK, opts=(), two_d=False, error=None):
        kw = {}
        if error is not None:
            kw = {"expect_error": error, "py_error": (cond.ConditionError,)}
        cs.append(
            Case(
                cid,
                code=LOAD2D if two_d else LOAD,
                args=[InPath("t"), rel, mcell(*chs), mcell(*mks), mcell(*opts)],
                py_args=[InPath("t"), rel, list(chs), list(mks), list(opts)],
                python=_load2d if two_d else _load3d,
                inputs=tree,
                **kw,
            )
        )

    add("load.movie", "one/Ex1_syn_z0p5/")
    add("load.movie_noslash", "one/Ex1_syn_z0p5")
    add("load.one_channel_markers", "one/Ex1_syn_z0p5/", ["ch1"], ["GFP"])
    add("load.marker_case", "one/Ex1_syn_z0p5/", CH, ["GFP", "mCherry"])
    add("load.condition", "multi/cond/")
    add("load.condition_by_name", "multi/cond/", opts=["IgnoreExpNumber", True])
    add("load.depths", "depths/cond/")
    add("load.depths_maxlevel1", "depths/cond/", opts=["maxlevel", 1])
    add("load.order_by_name", "order/Ex1_z0p5/")
    add("load.rate_median", "rate/Ex1_z0p5/")
    add("load.rate_option", "rate/Ex1_z0p5/", opts=["FrameRate", 5])
    for k, rel in enumerate(
        (
            "dz/Ex1_z0p25/",
            "dz/Ex1_z1p5/",
            "dz/Ex1_Z0P5/",
            "dz/cond_z0p3/Ex1/",
            "dz/q/Ex1_z0p5_z0p25/",
        )
    ):
        add(f"load.dz_token{k}", rel)
    add("load.dz_channel_token", "dz/chan/Ex1/", ["ch1_z0p7", "ch2_z0p7"])
    add(
        "load.options",
        "one/Ex1_syn_z0p5/",
        opts=[
            "Angle",
            32.8,
            "PixelSize",
            0.1,
            "Reversed",
            True,
            "dz",
            0.3,
            "ChannelPatterns",
            ["a"],
        ],
    )
    add("load.objective_scan", "one/Ex1_syn_z0p5/", opts=["ObjectiveScan", True, "Angle", 30])
    add("load.selector", "cells/", opts=["MovieSelector", "cell"])
    add("load.more_frames_ch2", "more/Ex1_z0p5/")
    add("load.extensions", "ext/Ex1_z0p5/")
    add("load.err_fewer_frames_ch2", "fewer/Ex1_z0p5/", error="MATLAB:badsubscript")
    add("load.err_2d_frames", "twod/Ex1_z0p5/", error="MATLAB:badsubscript")
    add("load.err_no_msec", "nomsec/Ex1_z0p5/", error="MATLAB:badsubscript")
    add("load.err_empty_channel", "empty/Ex1_z0p5/", error="MATLAB:cellRefFromNonCell")
    add(
        "load.err_missing_channel",
        "one/Ex1_syn_z0p5/",
        ["ch1", "ch3"],
        error="MATLAB:Java:InvalidInput",
    )
    add("load.err_no_movies", "nothing/", error=True)
    add("load.err_unknown_marker", "one/Ex1_syn_z0p5/", CH, ["gfp", "foo"],
        error="MATLAB:services:MissingRequiredCapability")  # fmt: skip
    add("load.err_fewer_markers", "one/Ex1_syn_z0p5/", CH, ["gfp"], error="MATLAB:badsubscript")
    add("load.err_no_number", "nonum/cond/", error="MATLAB:matrix:singleSubscriptNumelMismatch")
    add("load2d.movie", "one/Ex1_syn_z0p5/", two_d=True)
    add("load2d.ms_token", "ms/Ex1_40ms_z0p5/", two_d=True)
    add(
        "load2d.parameters", "one/Ex1_syn_z0p5/", opts=["Parameters", [1.1, 60, 6.5e-6]], two_d=True
    )
    add("load2d.frame_rate", "one/Ex1_syn_z0p5/", opts=["FrameRate", 5], two_d=True)
    add("load2d.no_msec", "nomsec/Ex1_z0p5/", two_d=True)
    # the demo dataset: three experiments of two channels and 60 frames (read in place)
    demo = {"d": Data(DEMO, copy=False)}
    for cid, which, two_d in (("demo.load", 1, False), ("demo.load_experiment", 2, False),
                              ("demo.load2d", 1, True)):  # fmt: skip
        load = (LOAD2D if two_d else LOAD).replace("[in{1} in{2}]", "[r rel]")
        load = load.replace(
            "rmfield(d, 'date')", "rmfield(d, intersect([{'date'}, PATHS], fieldnames(d)))"
        )
        cs.append(
            Case(
                cid,
                code=FIND + load,
                args=[InPath("d"), which, mcell(*CH), mcell(*MK), mcell()],
                py_args=[InPath("d"), which, CH, MK, []],
                python=_demo_load(two_d),
                inputs=demo,
                marks=("demodata", "slow"),
            )
        )
    return cs


def demo_rel(root: str, which: int) -> str:
    """the demo's condition folder (which 1) or its first experiment (2), relative to root"""
    c = sorted(e.name for e in os.scandir(root) if e.is_dir() and not e.name.startswith("."))[0]
    if which == 1:
        return c + "/"
    e = sorted(n for n in os.listdir(os.path.join(root, c)) if n.startswith("Ex1"))[0]
    return f"{c}/{e}/"


def _demo_load(two_d: bool):
    def run(root, which, chs, mks, opts):
        root = root if root.endswith("/") else root + "/"
        rel = demo_rel(root, which)
        out = (_load2d if two_d else _load3d)(root, rel, chs, mks, opts)
        return [{k: v for k, v in d.items() if k not in PATH_FIELDS} for d in out]

    return run


CASES = _cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# -- the rules MATLAB's probes showed (black box, digits-free probe roots) ----------------------


@pytest.mark.parametrize(
    "path,date",
    [
        ("/r/Lda/onefivezero/Ex1_z0p5/", None),
        ("/r/Ldb/150407_cond/Ex1_z0p5/", "150407"),
        ("/r/Ldc/20150407_cond/Ex1_z0p5/", "201504"),
        ("/r/Ldd/x1504_cond/Ex1_z0p5/", None),
        ("/r/Lde/c_123456_x_654321/Ex1_z0p5/", "123456654321"),
        ("/r/Ldf/a1234567/Ex1_z0p5/", "123456"),
        ("/r/Ldg/q/Ex123456_z0p5/", "123456"),
        ("/r/Ldi/x_12345/Ex1_z0p5/", None),
        ("/r/Ldj/2015-04-07/Ex1_z0p5/", None),
        ("/r/Ldk/abc987654321xyz/Ex1_z0p5/", "987654"),
        ("/r/Ldtb/c_111111/Ex2_222222_z0p5/", "111111222222"),
    ],
)
def test_date_digits(path, date):
    """date: the six-digit runs of the movie folder's path (probes L11, D_condition)"""
    digits = cond._SIX.findall(path)
    assert ("".join(digits) if digits else None) == date


def test_date_of_a_loaded_movie(tmp_path):
    root = tmp_path / "x_123456" / "Ex1_z0p5"
    write_movie(root, CH, std_names(2))
    (m,) = cond.load_condition_data(str(root), CH, MK)
    digits = "".join(cond._SIX.findall(str(root) + os.sep))
    assert m.date == digits and "123456" in digits
    assert m.to_matlab()["date"] == digits


@pytest.mark.parametrize(
    "name,dz",
    [
        ("Ex1_z0p25", 0.25), ("Ex1_z1", None), ("Ex1_z1p5", None), ("Ex1_none", None),
        ("Ex1_z0p5_z1", 0.5), ("Ex1_dz0p5", None), ("Ex1_z0.5", None), ("Ex1_z0p5x", 0.5),
        ("Ex1_Z0p5", 0.5), ("Ex1_z05", None), ("Ex1_z0p50", 0.5), ("Ex1_z2p125", None),
        ("z0p5_Ex1", None), ("Ex1_z0p125", 0.125), ("Ex1_z0p", None), ("Ex1_Z0P5", 0.5),
        ("Ex1-z0p5", None), ("Ex1_xz0p5", None), ("Ex1_z0p5.d", 0.5),
    ],
)  # fmt: skip
def test_dz_token(name, dz):
    """dz: '0.' and the digits of the first _z0p<digits> of the path, any case (probes L10, Z_)"""
    got = cond._dz_token(f"/r/q/{name}/ch1/")
    assert (np.isnan(got) and dz is None) or got == dz
    if dz is None:
        assert np.float64(got).view(np.uint64) == 0xFFF8000000000000  # MATLAB's NaN


def test_movie_round_trip():
    m = cond.Movie(
        source="/a/b/Ex1/ch1/",
        channels=["/a/b/Ex1/ch1/", "/a/b/Ex1/ch2/"],
        date=None,
        framerate=2.0,
        image_size=(3, 6, 8),
        movie_length=2,
        markers=["gfp", "rfp"],
        frame_paths=[["/a/b/Ex1/ch1/f1.tif", "/a/b/Ex1/ch1/f2.tif"]] * 2,
        mask_paths=["/a/b/Ex1/ch1/Analysis/Masks/dmask_001.tif"] * 2,
        dz=0.5,
        z_aniso=2.5,
        frame_paths_ds=[["x"], ["y"]],
        frame_paths_dsr=[["x"], ["y"]],
    )
    d = m.to_matlab()
    assert list(d) == [
        "source", "channels", "date", "framerate", "imagesize", "movieLength", "markers",
        "framePaths", "maskPaths", "pixelSize", "dz", "zAniso", "angle", "objectiveScan",
        "reversed", "framePathsDS", "framePathsDSR",
    ]  # fmt: skip
    assert d["imagesize"].tolist() == [6.0, 8.0, 3.0] and d["date"].shape == (0, 0)
    assert d["framePaths"][0].shape == (2, 1)
    back = cond.Movie.from_matlab(d)
    assert back == m


def test_short_path():
    """getShortPath (black box): the last folders of the first channel, without the first part"""
    assert cond.short_path("/r/a/b/c/d/ch1/") == "c/d/ch1"
    assert cond.short_path("/r/a/b/c/d/ch1/", 2) == "d/ch1"
    assert cond.short_path("/r/a/b/c/d/ch1") == "b/c/d"
    assert cond.short_path("rel/b/ch1/") == "b/ch1"
    assert cond.short_path("/ch1/") == "ch1"
    assert cond.short_path("/a/b/") == "a/b"


def test_paths_compare_with_slashes():
    """a condition the port loaded on Windows (MATLAB's filesep added to the folders as given)
    compares as the Linux reference: every separator a slash, nothing else changed"""
    win = {
        "source": "C:\\w\\t/one/Ex1_syn_z0p5/ch1\\",
        "channels": ["C:\\w\\t/one/Ex1_syn_z0p5/ch1\\", "C:\\w\\t/one/Ex1_syn_z0p5/ch2\\"],
        "framePaths": [np.array([["C:\\w\\t/one/Ex1_syn_z0p5/ch1\\a.tif"]], dtype=object)],
        "movieLength": 1.0,
        "markers": ["gfp"],
    }
    got = posix_paths(win)
    assert got["source"] == "C:/w/t/one/Ex1_syn_z0p5/ch1/"
    assert got["channels"] == ["C:/w/t/one/Ex1_syn_z0p5/ch1/", "C:/w/t/one/Ex1_syn_z0p5/ch2/"]
    assert got["framePaths"][0].shape == (1, 1)
    assert got["framePaths"][0][0, 0] == "C:/w/t/one/Ex1_syn_z0p5/ch1/a.tif"
    assert got["movieLength"] == 1.0 and got["markers"] == ["gfp"]
    again = posix_paths(got)  # a no-op on paths with slashes
    assert again["source"] == got["source"] and again["channels"] == got["channels"]


def test_z_aniso_sine():
    """zAniso = dz * sin(angle*pi/180) / pixelSize, bit for bit (probes L34, round 2)"""
    probes = {
        (31.5, 0.104): 2.5120123303651387, (31.5, 0.108): 2.4189748366479114,
        (32.8, 0.1): 2.7085410514136985, (45, 0.104): 3.3995518326276324,
        (30, 0.1): 2.4999999999999996, (59.9, 0.104): 4.1593818296620411,
        # angles where a*pi/180, a/180*pi and pi/180*a round apart (round 2 of the probes)
        (25.2, 0.104): 2.04701582483208, (27.4, 0.104): 2.212498965306979,
        (31.4, 0.104): 2.5048539992335397, (31.6, 0.104): 2.519163009471534,
        (32.9, 0.104): 2.6114156247628397, (33.0, 0.104): 2.6184568991107073,
        (34.6, 0.104): 2.7300180050629868, (36.8, 0.104): 2.879921146709547,
    }  # fmt: skip
    for (a, ps), z in probes.items():
        assert cond._z_aniso(0.5, a, ps, False) == z
    assert cond._z_aniso(0.5, 31.5, 0.104, True) == 0.5 / 0.104


def test_interactive_branches_are_errors(tmp_path):
    with pytest.raises(cond.ConditionError) as e:
        cond.load_condition_data()
    assert e.value.identifier == "MATLAB:services:MissingRequiredCapability"
    with pytest.raises(cond.ConditionError):
        cond.load_condition_data(str(tmp_path), "ch1", MK)
