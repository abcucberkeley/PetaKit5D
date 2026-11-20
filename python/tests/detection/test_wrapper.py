"""The deskew of a condition's frames (deskewData's behaviour) and the detection and tracking
wrapper (the lab's XR_runDetTrack3d): oracle cases on synthetic conditions, the wrapper's stages
replaced by MATLAB's recorded outputs, and its orchestration with stand-ins (task C04a)."""

import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from petakit5d.detection import condition as cond
from petakit5d.detection import wrapper
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, Data, InPath, Tree, mcell

RAW = (9, 24, 40)  # (z, y, x) of a raw frame
LLC = (
    "LLFFCorrection", False, "LowerLimit", 0.4, "LSImageCh1", "", "LSImageCh2", "",
    "LSImageCh3", "", "LSImageCh4", "", "BackgroundCh1", "", "BackgroundCh2", "",
    "BackgroundCh3", "", "BackgroundCh4", "",
)  # fmt: skip


def raw_frame(k: int, c: int, cls=np.uint16) -> np.ndarray:
    """a raw frame (z, y, x): ramps and a bright spot, MATLAB's round of the doubles"""
    z, y, x = np.meshgrid(*(np.arange(1, n + 1, dtype=np.float64) for n in RAW), indexing="ij")
    v = 100 + 3 * y + 2 * x + 5 * z + 7 * (k + 1) + 11 * c
    v = v + 400 * synth.stable_exp(-((y - 13 - k) ** 2 + (x - 20 + c) ** 2 + (z - 5) ** 2 * 4) / 8)
    t = np.trunc(v)
    v = t + np.where(np.abs(v - t) >= 0.5, np.sign(v), 0.0)
    return v.astype(np.uint16) if cls == np.uint16 else (v.astype(np.float32) + np.float32(0.25))


def frame_name(c: int, k: int) -> str:
    cam, wl = ("B", 560) if c == 1 else ("A", 642)
    t = 2000 * k
    return f"synth_Cam{cam}_ch{c - 1}_stack{k:04d}_{wl}nm_{t:07d}msec_{1000000 + t:010d}msecAbs.tif"


def write_raw(movie: Path, n=3, cls=np.uint16) -> None:
    for c in (1, 2):
        d = movie / f"ch{c}"
        d.mkdir(parents=True, exist_ok=True)
        for k in range(n):
            synth.write_tiff(d / frame_name(c, k), raw_frame(k, c, cls))


def marker(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    synth.write_tiff(path, np.zeros((1, 2, 2), dtype=np.float32))


def make_tree(folder) -> None:
    r = Path(folder)
    write_raw(r / "plain" / "Ex1_synth_z0p5")
    write_raw(r / "single" / "Ex1_synth_z0p5", cls=np.float32)
    write_raw(r / "two" / "cond" / "Ex1_synth_z0p5")
    write_raw(r / "two" / "cond" / "Ex2_synth_z0p5", n=2)
    m = r / "partial" / "Ex1_synth_z0p5"
    write_raw(m)
    for k in (0, 2):
        marker(m / "ch1" / "DS" / frame_name(1, k))
    m = r / "ch1done" / "Ex1_synth_z0p5"
    write_raw(m)
    for k in range(3):
        marker(m / "ch1" / "DS" / frame_name(1, k))
    m = r / "others" / "Ex1_synth_z0p5"
    write_raw(m)
    for k in range(3):
        marker(m / "ch1" / "DS" / f"zz_other_{k}.tif")
    with open(m / "ch1" / "DS" / "notes.txt", "w", newline="\n") as f:
        f.write("x\n")


def make_skip_tree(folder) -> None:
    """channel 1's DS folders with as many .tiff files as frames, and with one file more (a
    tree of its own)"""
    r = Path(folder)
    m = r / "tiffs" / "Ex1_synth_z0p5"
    write_raw(m)
    for k in range(3):
        marker(m / "ch1" / "DS" / (frame_name(1, k)[:-4] + ".tiff"))
    m = r / "extra" / "Ex1_synth_z0p5"
    write_raw(m)
    for k in range(3):
        marker(m / "ch1" / "DS" / frame_name(1, k))
    marker(m / "ch1" / "DS" / "zz_extra.tif")


LOAD = "d = XR_loadConditionData3D([in{1} in{2}], {'ch1', 'ch2'}, {'gfp', 'rfp'}, in{3}{:}); "
LOAD_OPTS = {"dz": "dz", "PixelSize": "pixel_size", "Angle": "angle", "Reversed": "reversed",
             "ObjectiveScan": "objective_scan"}  # fmt: skip


def wrapper_opts(overwrite=False, angle=31.5, aname="Analysis", flip=False) -> list:
    """deskewData's options as XR_runDetTrack3d passes them"""
    return [
        "Overwrite", overwrite, "SkewAngle", angle, "aname", aname, "Rotate", False,
        "sCMOSCameraFlip", flip, "Crop", False, *LLC, "zoffsetCorrection", False,
    ]  # fmt: skip


def _deskew(root, rel, load_opts, opts):
    lo = {LOAD_OPTS[load_opts[i]]: load_opts[i + 1] for i in range(0, len(load_opts), 2)}
    d = cond.load_condition_data(root + rel, ["ch1", "ch2"], ["gfp", "rfp"], **lo)
    o = {opts[i]: opts[i + 1] for i in range(0, len(opts), 2)}
    wrapper.deskew_data(
        d,
        overwrite=o.get("Overwrite", False),
        skew_angle=o.get("SkewAngle", 31.5),
        aname=o.get("aname"),
        s_cmos_camera_flip=o.get("sCMOSCameraFlip", False),
        reverse=o.get("Reverse", False),
    )


def _wrapper_deskew(root, rel, load_opts):
    """run_detection_tracking's deskew alone (no detection: an unknown parallel method; no
    tracking), with the scan direction of each movie (the fixed wrapper)"""
    lo = {LOAD_OPTS[load_opts[i]]: load_opts[i + 1] for i in range(0, len(load_opts), 2)}
    d = cond.load_condition_data(root + rel, ["ch1", "ch2"], ["gfp", "rfp"], **lo)
    old = dict(wrapper.STAGES)
    wrapper.STAGES.update(dict.fromkeys(wrapper.STAGES, lambda *_a, **_k: None))
    try:
        wrapper.run_detection_tracking(
            d, sigma=((1.3, 1.5), (1.3, 1.5)), parallel_method="none", track=False
        )
    finally:
        wrapper.STAGES.clear()
        wrapper.STAGES.update(old)


def _without_identifier(*args):
    try:
        _deskew(*args)
    except wrapper.WrapperError as err:
        assert err.identifier == "MATLAB:UndefinedFunction"
        raise RuntimeError(str(err)) from None


def _cases():
    tree = {"t": Tree(make_tree, "t")}
    cs = []

    def add(cid, rel, load_opts=(), opts=None, error=None, inputs=None):
        opts = wrapper_opts() if opts is None else opts
        kw = {"expect_error": error, "py_error": (wrapper.WrapperError, RuntimeError)}
        kw = kw if error else {}
        cs.append(
            Case(
                f"deskew.{cid}",
                code=LOAD + "deskewData(d, in{4}{:}); out = {};",
                args=[InPath("t"), rel, mcell(*load_opts), mcell(*opts)],
                py_args=[InPath("t"), rel, list(load_opts), list(opts)],
                python=_deskew,
                nargout=0,
                inputs=inputs or tree,
                **kw,
            )
        )

    add("wrapper_call", "plain/Ex1_synth_z0p5/")
    add("empty_aname", "plain/Ex1_synth_z0p5/", opts=wrapper_opts(aname=None))
    add(
        "demo_call",
        "plain/Ex1_synth_z0p5/",
        opts=["Overwrite", True, "SkewAngle", 31.5, "aname", "Analysis_update_3", "Rotate", False,
              "sCMOSCameraFlip", False, "Crop", False],
    )  # fmt: skip
    add("angle45", "plain/Ex1_synth_z0p5/", opts=wrapper_opts(angle=45))
    add("movie_fields", "plain/Ex1_synth_z0p5/", ["dz", 0.8, "PixelSize", 0.1, "Angle", 40])
    # a reversed scan: deskewData leaves data.reversed alone (forward, as without the flag),
    # its 'Reverse' option deskews in the scan's direction (detection.md W2)
    add("reversed", "plain/Ex1_synth_z0p5/", ["Reversed", True])
    add("reverse_option", "plain/Ex1_synth_z0p5/", ["Reversed", True],
        opts=[*wrapper_opts(), "Reverse", True])  # fmt: skip
    add("objective_scan", "plain/Ex1_synth_z0p5/", ["ObjectiveScan", True])
    add("single_frames", "single/Ex1_synth_z0p5/")
    add("two_movies", "two/cond/")
    add("partial_ds", "partial/Ex1_synth_z0p5/")
    add("ch1_done_skips", "ch1done/Ex1_synth_z0p5/")
    add("ch1_done_overwrite", "ch1done/Ex1_synth_z0p5/", opts=wrapper_opts(overwrite=True))
    add("other_tifs_skip", "others/Ex1_synth_z0p5/")
    # the count takes '*.tif*' files and must equal the frames' (measured: 3 .tiff skip, 4 not)
    skip = {"t": Tree(make_skip_tree, "t")}
    add("tiff_files_skip", "tiffs/Ex1_synth_z0p5/", inputs=skip)
    add("one_file_more_deskews", "extra/Ex1_synth_z0p5/", inputs=skip)
    # MATLAB:UndefinedFunction (GU_sCMOSCameraFlip) without a parallel pool, the pool's
    # MATLAB:parfor:SourceCodeNotAvailable with one (the oracle session has one): the error's
    # identifier is not compared here
    add("flip_fails", "plain/Ex1_synth_z0p5/", opts=wrapper_opts(flip=True), error=True)
    cs[-1].python = _without_identifier
    # the wrapper's deskew of a reversed scan, fixed (W2): XR_runDetTrack3d's deskewData call
    # (its 'aname' [] and results folder) with 'Reverse' as the movie's direction
    cs.append(
        Case(
            "deskew.wrapper_reversed_scan",
            code=LOAD + "mkdir([d.source 'Analysis']); deskewData(d, 'Overwrite', false, "
            "'SkewAngle', 31.5, 'aname', [], 'Rotate', false, 'sCMOSCameraFlip', false, "
            "'Crop', false, in{4}{:}, 'zoffsetCorrection', false, 'Reverse', logical(d.reversed)); "
            "out = {};",
            args=[InPath("t"), "plain/Ex1_synth_z0p5/", mcell("Reversed", True), mcell(*LLC)],
            py_args=[InPath("t"), "plain/Ex1_synth_z0p5/", ["Reversed", True]],
            python=_wrapper_deskew,
            nargout=0,
            inputs=tree,
        )
    )
    return cs


# -- XR_runDetTrack3d with MATLAB's recorded stage outputs ---------------------------------------

# a synthetic condition whose raw frames deskew into frames with spots, channel 2 shifted by 2
# voxels in x (spots in the raw frames' coordinates: a shift in x is the same after the deskew)
SPOTS = (
    (4.2, 10.3, 9.6),
    (6.4, 24.7, 22.4),
    (8.0, 40.1, 15.2),
    (5.4, 46.5, 28.3),
    (7.1, 16.2, 30.8),
)


def det_raw_frame(k: int, c: int, movie: int = 0) -> np.ndarray:
    """raw frames (z, y, x) = (12, 56, 40) with moving spots on a noisy background (another
    movie: other spots and noise)"""
    pos = []
    for j, p in enumerate(SPOTS):
        if movie and j % (movie + 1) == 0:
            continue
        q = [p[0] + 0.05 * k, p[1] + 0.3 * k + 2.0 * movie, p[2] + 0.2 * k - 1.5 * movie]
        if c == 2:
            q[2] += 2.0
        pos.append(q)
    vol = synth.beads((12, 56, 40), np.array(pos), (1.5, 1.3), np.full(len(pos), 300.0), 100.0)
    stream = f"dettrack/{c}/{k}" if movie == 0 else f"dettrack/m{movie}/{c}/{k}"
    noise = synth.CounterRNG(2026, stream).normal(vol.shape)
    v = vol + 4.0 * noise
    t = np.trunc(v)
    v = t + np.where(np.abs(v - t) >= 0.5, np.sign(v), 0.0)
    return np.clip(v, 0, 65535).astype(np.uint16)


N_DET = 8


def make_det_condition(folder) -> None:
    m = Path(folder) / "Ex1_synth_z0p5"
    for c in (1, 2):
        d = m / f"ch{c}"
        d.mkdir(parents=True, exist_ok=True)
        for k in range(N_DET):
            synth.write_tiff(d / frame_name(c, k), det_raw_frame(k, c))
    with open(m / "synth_Settings.txt", "w", newline="\n") as f:
        f.write('Camera type = "Orca4.0"\n')


DETTRACK = (
    "d = XR_loadConditionData3D(in{1}, {'ch1', 'ch2'}, {'gfp', 'rfp'}); "
    "XR_runDetTrack3d(d, 'Sigma', [1.3 1.5; 1.3 1.5], 'DetectionMethod', 'lowSNR', 'aname', "
    "'Analysis', 'FitMixtures', false, 'BackgroundCorrection', true, 'XZoffsetCorrection', true); "
    "out = {};"
)
# the files the stages of other tasks write (copied from MATLAB's run into the port's)
STAGE_FILES = ("Detection3D.mat", "Detections/", "Masks/", "Amira/", "trackedFeatures.mat",
               "ProcessedTracks.mat", "RotatedTracks.mat")  # fmt: skip


def _ported_tracking() -> dict:
    """the port's tracking stages where they exist (they then track MATLAB's detections, and
    their trackedFeatures.mat is compared with MATLAB's), else none"""
    try:
        from petakit5d.detection import tracking

        return {"track_settings": tracking.track_settings, "run_tracking": tracking.run_tracking}
    except (ImportError, AttributeError):
        return {}


def _recorded(case_id: str):
    def run(root):
        store = oracle.Store.from_env()
        case = next(c for c in CASES if c.id == case_id)
        if store is None or not store.has(case.key()):
            pytest.skip("the stages are replaced by MATLAB's recorded outputs: needs the store")
        src = store.case_dir(case.key()) / "files"
        work = Path.cwd()
        ported = _ported_tracking()
        own = ("trackedFeatures.mat",) if ported else ()

        def copy(*_a, **_k):
            for p in src.rglob("*"):
                rel = p.relative_to(src).as_posix()
                if any(s in rel for s in own):
                    continue
                if any(s in rel for s in STAGE_FILES) and not (work / rel).exists():
                    if p.is_dir():
                        (work / rel).mkdir(parents=True, exist_ok=True)
                    else:
                        (work / rel).parent.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(p, work / rel)

        stand_ins = dict.fromkeys(wrapper.STAGES, copy)
        stand_ins["track_settings"] = lambda **_k: None
        stand_ins.update(ported)
        old = dict(wrapper.STAGES)
        wrapper.STAGES.update(stand_ins)
        try:
            d = cond.load_condition_data(root, ["ch1", "ch2"], ["gfp", "rfp"])
            wrapper.run_detection_tracking(
                d,
                sigma=((1.3, 1.5), (1.3, 1.5)),
                detection_method="lowSNR",
                aname="Analysis",
                fit_mixtures=False,
                background_correction=True,
                x_zoffset_correction=True,
            )
        finally:
            wrapper.STAGES.clear()
            wrapper.STAGES.update(old)

    return run


# a condition whose middle movie has one frame (W3); its spots differ from movie to movie
ONE_FRAME = (("Ex1_synth_z0p5", 6), ("Ex2_synth_z0p5", 1), ("Ex3_synth_z0p5", 6))


def make_one_frame_condition(folder) -> None:
    for i, (name, n) in enumerate(ONE_FRAME):
        m = Path(folder) / name
        for c in (1, 2):
            d = m / f"ch{c}"
            d.mkdir(parents=True, exist_ok=True)
            for k in range(n):
                synth.write_tiff(d / frame_name(c, k), det_raw_frame(k, c, movie=i))
        with open(m / "synth_Settings.txt", "w", newline="\n") as f:
            f.write('Camera type = "Orca4.0"\n')


# the fixed XR_runDetTrack3d on that condition, as MATLAB gives it: the tracking stages on the
# movies of more than one frame with their own folders, the one-frame movie deskewed and
# detected only (each movie loaded alone, as MATLAB's loader asks about a one-frame movie's NaN
# frame rate among others)
DETTRACK_ONE = (
    "r = in{1}; m = in{2}; d = []; for k = 1 : numel(m), "
    "d = [d, XR_loadConditionData3D([r m{k}], {'ch1', 'ch2'}, {'gfp', 'rfp'})]; end; "
    "o = {'Sigma', [1.3 1.5; 1.3 1.5], 'DetectionMethod', 'lowSNR', 'aname', 'Analysis', "
    "'FitMixtures', false, 'BackgroundCorrection', true}; keep = [d.movieLength] ~= 1; "
    "XR_runDetTrack3d(d(keep), o{:}); XR_runDetTrack3d(d(~keep), o{:}, 'Track', false); out = {};"
)


def _recorded_one_frame(root, movies):
    """the port's wrapper on the joined condition: detection, track processing and rotation
    copied from MATLAB's run, tracking the port's where it exists; each stage must get its own
    movie's folder"""
    store = oracle.Store.from_env()
    case = next(c for c in CASES if c.id == "dettrack.one_frame_movie")
    if store is None or not store.has(case.key()):
        pytest.skip("the stages are replaced by MATLAB's recorded outputs: needs the store")
    src = store.case_dir(case.key()) / "files"
    work = Path.cwd()
    ported = _ported_tracking()
    if not ported:
        pytest.skip(
            "the case checks the port's tracking's folders: needs petakit5d.detection.tracking"
        )
    seen = []

    def copy(*a, **k):
        movies = [a[0]] if isinstance(a[0], dict) else list(a[0])
        seen.append(([m["source"] for m in movies], k["results_path"]))
        for p in src.rglob("*"):
            rel = p.relative_to(src).as_posix()
            if "trackedFeatures.mat" in rel:
                continue
            if any(s in rel for s in STAGE_FILES) and not (work / rel).exists():
                if p.is_dir():
                    (work / rel).mkdir(parents=True, exist_ok=True)
                else:
                    (work / rel).parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(p, work / rel)

    old = dict(wrapper.STAGES)
    wrapper.STAGES.update(dict.fromkeys(wrapper.STAGES, copy))
    wrapper.STAGES.update(ported)
    try:
        d = cond.Condition(
            m for rel in movies for m in cond.load_condition_data(root + rel, ["ch1", "ch2"],
                                                                   ["gfp", "rfp"])
        )  # fmt: skip
        wrapper.run_detection_tracking(
            d,
            sigma=((1.3, 1.5), (1.3, 1.5)),
            detection_method="lowSNR",
            aname="Analysis",
            fit_mixtures=False,
            background_correction=True,
        )
    finally:
        wrapper.STAGES.clear()
        wrapper.STAGES.update(old)
    # one stand-in call each (detection of the three movies, processing and rotation of the two
    # of more than one frame), every movie with its own results folder
    for sources, results in seen:
        assert results == [s + "Analysis" + os.sep for s in sources], (sources, results)
    assert [len(s) for s, _ in seen] == [3, 2, 2]


CASES = _cases() + [
    Case(
        "dettrack.recorded_stages",
        code=DETTRACK,
        # the other tasks' stage files are MATLAB's, copied; trackedFeatures.mat is the port's
        # tracking's where it exists (the harness reads v7 files with mat_dtype=True, so
        # MATLAB's integer-valued doubles, stored as small integers, compare as double)
        args=[InPath("cond")],
        python=_recorded("dettrack.recorded_stages"),
        py_args=[InPath("cond")],
        nargout=0,
        inputs={"cond": Tree(make_det_condition, "cond")},
        marks=("slow",),
    ),
    Case(
        "dettrack.one_frame_movie",
        code=DETTRACK_ONE,
        # as above; the tracking the port's (the case needs petakit5d.detection.tracking)
        args=[InPath("cond"), mcell(*[f"{n}/" for n, _ in ONE_FRAME])],
        python=_recorded_one_frame,
        py_args=[InPath("cond"), [f"{n}/" for n, _ in ONE_FRAME]],
        nargout=0,
        inputs={"cond": Tree(make_one_frame_condition, "cond")},
        marks=("slow",),
    ),
]


# the demo dataset: its first experiment's first frames, deskewed and their offsets corrected
# (the raw copies are removed at the end, so that only the outputs are compared and stored)
DEMO = "$PETAKIT5D_DEMO_DATA/PetaKit5D_point_detection_demo_dataset"
N_DEMO = 4
DEMO_CODE = (
    "r = in{1}; if r(end) ~= '/', r = [r '/']; end; f = dir(r); "
    "f = f([f.isdir] & ~startsWith({f.name}, '.')); c = [f(1).name '/']; "
    "e = dir([r c 'Ex1*']); src = [r c e(1).name '/']; m = [pwd '/demo/Ex1_demo_z0p5/']; "
    "for ch = {'ch1', 'ch2'}, mkdir([m ch{1}]); g = dir([src ch{1} '/*.tif']); "
    "for k = 1 : in{2}, t = regexp(g(k).name, '(\\d+)msec', 'tokens', 'once'); "
    "copyfile([src ch{1} '/' g(k).name], sprintf('%s%s/demo_%s_stack%04d_%smsec.tif', m, ch{1}, "
    "ch{1}, k - 1, t{1})); end; end; "
    "d = XR_loadConditionData3D(m, {'ch1', 'ch2'}, {'gfp', 'rfp'}); deskewData(d, in{3}{:}); "
    "XR_correctXZoffsetData3D(d, 'Overwrite', false); "
    "for ch = {'ch1', 'ch2'}, delete([m ch{1} '/*.tif']); end; out = {};"
)


def _demo(root, n, opts):
    import glob
    import os
    import re

    from petakit5d.detection.offsets import correct_xz_offset_data

    root = root if root.endswith("/") else root + "/"
    c = sorted(e.name for e in os.scandir(root) if e.is_dir() and not e.name.startswith("."))[0]
    e = sorted(x for x in os.listdir(root + c) if x.startswith("Ex1"))[0]
    src = f"{root}{c}/{e}/"
    m = Path.cwd() / "demo" / "Ex1_demo_z0p5"
    for ch in ("ch1", "ch2"):
        (m / ch).mkdir(parents=True)
        for k, p in enumerate(sorted(glob.glob(f"{src}{ch}/*.tif"))[:n]):
            t = re.search(r"(\d+)msec", os.path.basename(p)).group(1)
            shutil.copyfile(p, m / ch / f"demo_{ch}_stack{k:04d}_{t}msec.tif")
    d = cond.load_condition_data(str(m), ["ch1", "ch2"], ["gfp", "rfp"])
    o = {opts[i]: opts[i + 1] for i in range(0, len(opts), 2)}
    wrapper.deskew_data(d, overwrite=o["Overwrite"], skew_angle=o["SkewAngle"], aname=o["aname"])
    correct_xz_offset_data(d, overwrite=False)
    for ch in ("ch1", "ch2"):
        for p in (m / ch).glob("*.tif"):
            p.unlink()


CASES.append(
    Case(
        "demo.deskew_offsets",
        code=DEMO_CODE,
        args=[InPath("d"), N_DEMO, mcell(*wrapper_opts())],
        py_args=[InPath("d"), N_DEMO, wrapper_opts()],
        python=_demo,
        nargout=0,
        inputs={"d": Data(DEMO, copy=False)},
        marks=("demodata", "slow"),
    )
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# -- the orchestration, with stand-ins ----------------------------------------------------------


def _joined(root: Path, frames) -> cond.Condition:
    """movies Ex1, Ex2, ... of the given frame counts, each loaded alone and joined (MATLAB's
    loader asks about a one-frame movie's NaN frame rate among others)"""
    out = []
    for i, n in enumerate(frames, start=1):
        m = root / f"Ex{i}_synth_z0p5"
        write_raw(m, n=n)
        out += cond.load_condition_data(str(m) + os.sep, ["ch1", "ch2"], ["gfp", "rfp"])
    return cond.Condition(out)


@pytest.fixture
def condition(tmp_path):
    return _joined(tmp_path, (2, 1))


@pytest.fixture
def calls():
    log = []

    def rec(name):
        def f(*a, **k):
            log.append((name, a, k))

        return f

    old = dict(wrapper.STAGES)
    wrapper.STAGES.update({n: rec(n) for n in wrapper.STAGES})
    yield log
    wrapper.STAGES.clear()
    wrapper.STAGES.update(old)


# XR_runDetTrack3d's detection: 'parfor' on Windows (ispc) whatever ParallelMethod, else the
# method given ('lock_distributed' by default: the frame-parallel detection)
DETECTION = "run_detection" if os.name == "nt" else "run_detection_frame_parallel"


TRACKING = ("run_tracking", "run_track_processing", "rotate_tracks")


def test_stages_in_order(condition, calls):
    wrapper.run_detection_tracking(condition, sigma=((1.3, 1.5), (1.3, 1.5)))
    assert [c[0] for c in calls] == [
        DETECTION, "track_settings", "run_tracking", "run_track_processing", "rotate_tracks",
    ]  # fmt: skip
    det = calls[0][2]
    assert det["overwrite"] is False and det["detection_method"] == "lowSNR"
    paths = [m.source + "Analysis" + os.sep for m in condition]
    assert det["results_path"] == paths and all(Path(p).is_dir() for p in paths)
    # each stage is called once with the condition as MATLAB's struct array (a StructArray of
    # the two movies' dicts) and the movies' folders
    from petakit5d.io.matfile import StructArray

    structs = calls[0][1][0]
    assert isinstance(structs, StructArray) and structs.shape == (1, 2)
    assert [m["movieLength"] for m in structs] == [2.0, 1.0]
    # the one-frame movie is dropped from the tracking; the other keeps its own folder (W3)
    for call in calls[2:]:
        assert isinstance(call[1][0], dict) and call[1][0]["movieLength"] == 2.0
        assert call[2]["results_path"] == [paths[0]]
    assert calls[1][2] == {"radius": (3, 6), "max_gap_length": 2}


def test_tracking_folders_after_a_one_frame_movie(tmp_path, calls):
    """a one-frame movie first: the second movie keeps its own folder (W3); MATLAB's pairing
    (its i-th movie with the i-th folder of all movies) with matlab_quirks"""
    d = _joined(tmp_path, (1, 2))
    paths = [m.source + "Analysis" + os.sep for m in d]
    wrapper.run_detection_tracking(d, sigma=((1.3, 1.5),), parallel_method="none")
    got = [(c[1][0]["source"], c[2]["results_path"]) for c in calls if c[0] in TRACKING]
    assert got == [(d[1].source, [paths[1]])] * 3
    calls.clear()
    # MATLAB passes the folders of all movies: its first (only) movie gets the first folder
    wrapper.run_detection_tracking(d, sigma=((1.3, 1.5),), parallel_method="none",
                                   matlab_quirks=True)  # fmt: skip
    got = [(c[1][0]["source"], c[2]["results_path"]) for c in calls if c[0] in TRACKING]
    assert got == [(d[1].source, paths)] * 3


def test_reversed_scans_deskewed_in_their_direction(tmp_path, calls, monkeypatch):
    """the wrapper deskews each movie in its direction (W2); matlab_quirks: deskewData's
    default, forward; deskew_data(reverse=None) is the deskew of each frame with its movie's
    reverse"""
    from petakit5d.deskew_rotate.volume import deskew_volume
    from petakit5d.io import read_tiff

    m = tmp_path / "Ex1_synth_z0p5"
    write_raw(m, n=1)
    d = cond.load_condition_data(str(m) + os.sep, ["ch1", "ch2"], ["gfp", "rfp"], reversed=True)
    seen = []
    real = wrapper.deskew_data
    monkeypatch.setattr(wrapper, "deskew_data", lambda data, **k: seen.append(k["reverse"]) or data)
    for q in (False, True):
        wrapper.run_detection_tracking(d, sigma=((1.3, 1.5),), parallel_method="none",
                                       track=False, matlab_quirks=q)  # fmt: skip
    assert seen == [None, False]
    for rev, want in ((None, True), (False, False), (True, True)):
        real(d, overwrite=True, reverse=rev)
        raw = np.asarray(read_tiff(d[0].frame_paths[0][0]), dtype=np.float64)
        ds = deskew_volume(raw, 31.5, d[0].dz, d[0].pixel_size, want).astype(np.float32)
        np.testing.assert_array_equal(read_tiff(d[0].frame_paths_ds[0][0]), ds)


def test_parallel_method_and_flags(condition, calls):
    wrapper.run_detection_tracking(
        condition, sigma=((1.3, 1.5),), parallel_method="parfor", track=False,
        overwrite=np.array([True, False, True, False]),
    )  # fmt: skip
    assert [c[0] for c in calls] == ["run_detection"]
    assert calls[0][2]["overwrite"] is False
    calls.clear()
    wrapper.run_detection_tracking(condition, sigma=((1.3, 1.5),), parallel_method="other",
                                   rotate_tracks=False)  # fmt: skip
    pc = ["run_detection"] if os.name == "nt" else []  # Windows: parfor whatever the method
    assert [c[0] for c in calls] == [*pc, "track_settings", "run_tracking", "run_track_processing"]


def test_windows_takes_parfor(condition, calls, monkeypatch):
    """MATLAB's ispc: on Windows the detection is 'parfor' (run_detection) whatever the
    parallel method; elsewhere the method given"""
    monkeypatch.setattr(wrapper, "_is_pc", lambda: True)
    for method in ("lock_distributed", "parfor", "other"):
        calls.clear()
        wrapper.run_detection_tracking(condition, sigma=((1.3, 1.5),), parallel_method=method,
                                       track=False)  # fmt: skip
        assert [c[0] for c in calls] == ["run_detection"]
    monkeypatch.setattr(wrapper, "_is_pc", lambda: False)
    for method, want in (("lock_distributed", ["run_detection_frame_parallel"]),
                         ("parfor", ["run_detection"]), ("other", [])):  # fmt: skip
        calls.clear()
        wrapper.run_detection_tracking(condition, sigma=((1.3, 1.5),), parallel_method=method,
                                       track=False)  # fmt: skip
        assert [c[0] for c in calls] == want


def test_apath_takes_the_last_folders(condition, calls, tmp_path):
    out = str(tmp_path / "res") + "/"
    wrapper.run_detection_tracking(condition, out, "R", sigma=((1.3, 1.5),), track=False)
    p = calls[0][2]["results_path"][0]
    assert p.startswith(out) and p.rstrip("/\\").endswith("R")


def test_stage_homes(monkeypatch):
    """the detection stage is detection.run's run_detection (XR_runDetection3D); the
    frame-parallel one (the default parallel_method) is not ported and says what to use"""
    from petakit5d.detection import run

    monkeypatch.setitem(wrapper.STAGES, "run_detection", None)
    monkeypatch.setitem(wrapper.STAGES, "run_detection_frame_parallel", None)
    assert wrapper._stage("run_detection") is run.run_detection
    if not hasattr(run, "run_detection_frame_parallel"):
        with pytest.raises(NotImplementedError, match="parallel_method='parfor'"):
            wrapper._stage("run_detection_frame_parallel")


def test_missing_stage_raises(condition, monkeypatch):
    # the first stage of the Linux order, the frame-parallel detection, is not ported; Windows
    # takes MATLAB's parfor path (ispc), whose run_detection is, so take the Linux order everywhere
    monkeypatch.setattr(wrapper, "_is_pc", lambda: False)
    old = dict(wrapper.STAGES)
    try:
        wrapper.STAGES.update(dict.fromkeys(wrapper.STAGES))
        with pytest.raises(NotImplementedError, match="C02"):
            wrapper.run_detection_tracking(condition, sigma=((1.3, 1.5),))
    finally:
        wrapper.STAGES.update(old)


def test_interactive_and_bad_flags(condition, calls):
    with pytest.raises(wrapper.WrapperError):
        wrapper.run_detection_tracking(condition, sigma=None)
    with pytest.raises(wrapper.WrapperError):  # no data: MATLAB's interactive loader
        wrapper.run_detection_tracking(sigma=((1.3, 1.5),))
    with pytest.raises(wrapper.WrapperError) as e:
        wrapper.run_detection_tracking(
            condition, sigma=((1.3, 1.5),), overwrite=np.array([True, False])
        )
    assert e.value.identifier == "MATLAB:badsubscript"
