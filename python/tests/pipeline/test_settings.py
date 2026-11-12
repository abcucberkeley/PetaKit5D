"""Setting files: XR_parseSettingFiles and XR_parseSettingFiles_wrapper.

The setting files here are made up: the lines around the stage line only look like the ones
acquisition software writes. Both stage formats (X stage, XZ stage) and the ways the line can
be missed (first or past the 101st line, last line without a newline, CRLF) are cases.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from petakit5d.pipeline.filenames import MatlabError
from petakit5d.pipeline.settings import (
    _fgetl_lines,
    parse_image_setting_files,
    parse_setting_files,
)
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, InPath, Tree, mcell

HEADER = [
    "***** ***** *****   General   ***** ***** *****",
    "Date :\t1/2/2024 3:04:05 PM",
    "Acq Mode :\tZ stack",
    "Version :\tv 0.0.0 test",
    "",
    "***** ***** *****   Waveform   ***** ***** *****",
    "Waveform type :\tLinear",
    "Cycle lasers :\tper Z",
    "Z motion :\tSample piezo",
    "",
    "# of stacks (0) :\t1",
    "Excitation Filter, Laser, Power (%), Exp(ms) (0) :\tN/A\t488\t10\t20",
    "S PZT Offset, Interval (um), # of Pixels for Excitation (0) :\t25\t0\t101",
    "Z PZT Offset, Interval (um), # of Pixels for Excitation (0) :\t0\t0.2\t101",
]
FOOTER = [
    "Z Galvo Offset, Interval (um), # of Pixels for Excitation (0) :\t0\t0\t101",
    "",
    "***** ***** *****   Camera   ***** ***** *****",
    "Camera :\tA",
    "",
]


def stage_line(stage: str = "X", off: str = "0", interval: str = "0.5", n: str = "101") -> str:
    head = f"{stage} Stage Offset, Interval (um), # of Pixels for Excitation (0) :"
    return f"{head}\t{off}\t{interval}\t{n}"


def settings_text(line: str | None = None, *, before=HEADER, after=FOOTER, eol="\n", final=True):
    lines = list(before) + ([line if line is not None else stage_line()]) + list(after)
    return eol.join(lines) + (eol if final else "")


EMPTY = mcell()

# name -> text of a setting file
SETTINGS = {
    "x_stage": settings_text(),
    "x_stage_negative": settings_text(stage_line(interval="-0.8")),
    "xz_stage": settings_text(stage_line("XZ", interval="-1.25")),
    "xz_spaces": settings_text(stage_line("XZ   ", interval="2")),
    "lower_case": settings_text(stage_line("x").lower()),
    "interval_forms_1": settings_text(stage_line(interval=".5")),
    "interval_forms_2": settings_text(stage_line(interval="-12")),
    "interval_forms_3": settings_text(stage_line(interval="1,5")),
    "interval_forms_4": settings_text(stage_line(interval="0p5")),
    "interval_forms_5": settings_text(stage_line(off="7", interval="-0.0", n="0")),
    "interval_tab_inside": settings_text(stage_line(interval="1\t2")),
    "two_digit_offset": settings_text(stage_line(off="12")),
    "no_stage_line": settings_text(""),
    "z_stage_only": settings_text(stage_line("Z")),
    "first_line": settings_text(before=[]),
    "last_line_newline": settings_text(after=[]),
    "last_line_no_newline": settings_text(after=[""], final=False),
    "second_last_no_newline": settings_text(after=["Camera :\tA"], final=False),
    "crlf": settings_text(eol="\r\n"),
    "crlf_stage_line": settings_text(stage_line() + "\r", after=FOOTER),
    "cr_only": settings_text(eol="\r"),
    "cr_last_line": settings_text(after=["x"], eol="\r"),
    "mixed_eol": "Settings\r" + stage_line() + "\r\nEnd\n\n",
    "crlf_no_final": settings_text(after=[""], eol="\r\n", final=False),
    "line_99": settings_text(before=HEADER + ["x"] * (98 - len(HEADER)), after=["", "y"]),
    "line_100": settings_text(before=HEADER + ["x"] * (99 - len(HEADER)), after=["", "y"]),
    "line_101": settings_text(before=HEADER + ["x"] * (100 - len(HEADER)), after=["", "y"]),
    "two_stage_lines": settings_text(after=[stage_line(interval="-3")] + FOOTER),
    "empty": "",
}

SETTING_CASES = [
    "x_stage",
    "x_stage_negative",
    "xz_stage",
    "xz_spaces",
    "lower_case",
    "interval_forms_1",
    "interval_forms_2",
    "interval_forms_3",
    "interval_forms_4",
    "interval_forms_5",
    "interval_tab_inside",
    "line_99",
    "line_100",
    "two_stage_lines",
    "crlf",
    "crlf_stage_line",
    "cr_only",
    "mixed_eol",
    "last_line_newline",
    "last_line_no_newline",
    "second_last_no_newline",
    "crlf_no_final",
    "cr_last_line",
]
SETTING_ERRORS = [
    "two_digit_offset",
    "no_stage_line",
    "z_stage_only",
    "first_line",
    "line_101",
    "empty",
]


def make_settings(root: Path) -> None:
    for name, text in SETTINGS.items():
        (root / f"{name}_Settings.txt").write_bytes(text.encode())
    (root / "folder_Settings.txt").mkdir()


ML_SETTINGS = "d = in{1}; out = {XR_parseSettingFiles(cellfun(@(x) [d x], in{2}, 'unif', 0))};"


def py_settings(d, names):
    return parse_setting_files([d + n for n in names])


def settings_out(info):
    return [info if info else {}]


def _settings_case(cid, names, error=False):
    return Case(
        cid,
        code=ML_SETTINGS,
        args=[InPath("s"), names],
        py_args=[InPath("s"), [] if names is EMPTY else names],
        python=py_settings,
        py_out=settings_out,
        inputs={"s": Tree(make_settings, "s")},
        expect_error=error,
        py_error=MatlabError,
        warnings="ids",
    )


# image folders for the wrapper: the image files can be empty (the wrapper only looks for the
# setting files next to them)


def tile(cam="A", ch=0, x=0, y=0, z=0, t=0, it="0000", sub="_0000", prefix="", suffix=""):
    cam_s = f"_Cam{cam}" if cam else ""
    return (
        f"{prefix}Scan_Iter_{it}{sub}{cam_s}_ch{ch}_CAM1_stack0000_488nm_0000000msec_"
        f"0001234567msecAbs_{x:03d}x_{y:03d}y_{z:03d}z_{t:04d}t{suffix}"
    )


TWICE = (
    tile("A", 0, it="0001", sub="", prefix="a") + "-" + tile("B", 1, it="0002", sub="", prefix="b")
)


def _write(folder: Path, name: str, text: str | None = None) -> None:
    (folder / name).write_bytes(b"" if text is None else text.encode())


def make_images(root: Path) -> None:
    # own: every image with <name>_Settings.txt
    d = root / "own"
    d.mkdir()
    for i, interval in enumerate(["0.5", "-0.5"]):
        n = tile("A", 0, y=i)
        _write(d, n + ".tif")
        _write(d, n + "_Settings.txt", settings_text(stage_line(interval=interval)))
    _write(d, tile("A", 0, y=0, suffix="_part0001") + ".tif")
    # tile: prefix<...>Scan_Iter_<iter>_<x>x_<y>y_<z>z_<t>t_Settings.txt, with a prefix
    d = root / "tile"
    d.mkdir()
    for cam in ("A", "B"):
        _write(d, tile(cam, 0, x=1, prefix="emb_") + ".tif")
    _write(d, "emb_Scan_Iter_0000_0000_001x_000y_000z_0000t_Settings.txt", settings_text())
    # iter: prefix<...>Scan_Iter_<iter>_Settings.txt, with and without a prefix
    d = root / "iter"
    d.mkdir()
    for ch in (0, 1):
        _write(d, tile(None, ch, it="0002", sub="", prefix="cellA_") + ".tif")
    _write(d, "cellA_Scan_Iter_0002_Settings.txt", settings_text(stage_line("XZ", interval="-2")))
    d = root / "noprefix"
    d.mkdir()
    _write(d, tile("A", 0, it="0001", sub="") + ".tif")
    _write(d, tile("A", 0, x=2, it="0001", sub="") + ".tif")
    _write(d, "Scan_Iter_0001_Settings.txt", settings_text())
    _write(d, "Scan_Iter_0001_000x_000y_000z_0000t_Settings.txt", settings_text())
    _write(d, "Scan_Iter_0001_002x_000y_000z_0000t_Settings.txt", settings_text())
    # an image name the grammar matches twice: sprintf takes each field of every match in
    # turn (tmp{x}.prefix, tmp{x}.Iter, ...) and repeats the format for the rest
    d = root / "twice"
    d.mkdir()
    _write(d, TWICE + ".tif")
    (d / "aScan_Iter_b_Settings.txt0001").mkdir()
    _write(d / "aScan_Iter_b_Settings.txt0001", "0002Scan_Iter_", settings_text())
    d = root / "none"
    d.mkdir()
    _write(d, tile("A", 0) + ".tif")
    _write(d, "plain_ch0_x.tif")
    _write(d, "other.tif")
    _write(d, "other_Settings.txt", settings_text())


ML_WRAPPER = (
    "d = in{1}; out = {XR_parseSettingFiles_wrapper(cellfun(@(x) [d x], in{2}, 'unif', 0))};"
)


def py_wrapper(d, names):
    return parse_image_setting_files([d + n for n in names])


def _wrapper_case(cid, names, error=False):
    return Case(
        cid,
        code=ML_WRAPPER,
        args=[InPath("img"), names],
        py_args=[InPath("img"), [] if names is EMPTY else names],
        python=py_wrapper,
        py_out=settings_out,
        inputs={"img": Tree(make_images, "img")},
        expect_error=error,
        py_error=MatlabError,
        warnings="ids",
    )


OWN = [f"own/{tile('A', 0, y=i)}.tif" for i in range(2)]
TILE = [f"tile/{tile(c, 0, x=1, prefix='emb_')}.tif" for c in ("A", "B")]
ITER = [f"iter/{tile(None, c, it='0002', sub='', prefix='cellA_')}.tif" for c in (0, 1)]
NOPREFIX = [
    f"noprefix/{tile('A', 0, it='0001', sub='')}.tif",
    f"noprefix/{tile('A', 0, x=2, it='0001', sub='')}.tif",
]

CASES = (
    [_settings_case(f"settings_{n}", [f"{n}_Settings.txt"]) for n in SETTING_CASES]
    + [_settings_case(f"settings_{n}", [f"{n}_Settings.txt"], True) for n in SETTING_ERRORS]
    + [
        _settings_case("settings_several", [f"{n}_Settings.txt" for n in SETTING_CASES]),
        _settings_case("settings_none", EMPTY),
        _settings_case("settings_first_missing", ["missing_Settings.txt"]),
        _settings_case(
            "settings_middle_missing",
            ["x_stage_Settings.txt", "missing_Settings.txt", "xz_stage_Settings.txt"],
        ),
        _settings_case("settings_last_missing", ["x_stage_Settings.txt", "missing_Settings.txt"]),
        _settings_case("settings_folder", ["folder_Settings.txt"], True),
        _wrapper_case("wrapper_own_settings", OWN),
        _wrapper_case(
            "wrapper_partial_left_out", OWN + [f"own/{tile('A', 0, suffix='_part0001')}.tif"]
        ),
        _wrapper_case("wrapper_tile_settings", TILE),
        _wrapper_case("wrapper_iter_settings", ITER),
        _wrapper_case("wrapper_no_prefix", NOPREFIX),
        _wrapper_case("wrapper_mixed_folders", ITER + OWN + TILE + OWN[:1]),
        _wrapper_case("wrapper_mixed_own_and_tile", OWN + TILE + OWN[:1]),
        _wrapper_case("wrapper_no_setting_file", ["none/" + tile("A", 0) + ".tif"]),
        _wrapper_case("wrapper_no_grammar", ["none/plain_ch0_x.tif", "none/other.tif"]),
        _wrapper_case("wrapper_only_own_without_grammar", ["none/other.tif"]),
        _wrapper_case("wrapper_empty", EMPTY),
        _wrapper_case("wrapper_name_matched_twice", [f"twice/{TWICE}.tif"]),
        _wrapper_case("wrapper_only_partial", [f"own/{tile('A', 0, suffix='_part0001')}.tif"]),
    ]
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# unit tests
# ---------------------------------------------------------------------------------------------


def _file(tmp_path, text, name="a_Settings.txt"):
    p = tmp_path / name
    p.write_bytes(text.encode())
    return str(p)


@pytest.mark.parametrize(
    ("text", "lines"),
    [
        ("a\nb\n", ["a", "b"]),
        ("a\nb", ["a", "b"]),
        ("", []),
        ("a\r\nb\r\n", ["a", "b"]),
        ("a\rb\rc", ["a", "b", "c"]),
        ("a\rb\r", ["a", "b"]),
        ("\n".join(str(i) for i in range(200)), [str(i) for i in range(100)]),
    ],
)
def test_fgetl_loop(tmp_path, text, lines):
    assert _fgetl_lines(_file(tmp_path, text)) == lines


def test_both_stage_formats(tmp_path):
    a = _file(tmp_path, settings_text(stage_line("X", interval="-0.5")), "a_Settings.txt")
    b = _file(tmp_path, settings_text(stage_line("XZ", off="3", interval="2")), "b_Settings.txt")
    info = parse_setting_files([a, b])
    assert info == [
        {"stageOff": 0.0, "StageInterval": -0.5, "numExiciation": 101.0},
        {"stageOff": 3.0, "StageInterval": 2.0, "numExiciation": 101.0},
    ]


def test_missing_files(tmp_path):
    a = _file(tmp_path, settings_text())
    with pytest.warns(UserWarning, match="missing.txt does not exist"):
        assert parse_setting_files([str(tmp_path / "missing.txt")]) == []
    with pytest.warns(UserWarning, match="missing.txt does not exist"):
        info = parse_setting_files([a, str(tmp_path / "missing.txt"), a])
    assert info[1] == {"stageOff": None, "StageInterval": None, "numExiciation": None}
    with pytest.warns(UserWarning):
        assert len(parse_setting_files([a, str(tmp_path / "missing.txt")])) == 1


def test_stage_line_errors(tmp_path):
    with pytest.raises(MatlabError) as e:
        parse_setting_files([_file(tmp_path, settings_text(""))])
    assert e.value.identifier == "MATLAB:badsubscript"


def test_wrapper_reads_each_file_once(tmp_path, monkeypatch):
    make_images(tmp_path)
    from petakit5d.pipeline import settings

    calls = []
    real = settings.parse_setting_files
    monkeypatch.setattr(settings, "parse_setting_files", lambda f: calls.append(f) or real(f))
    root = str(tmp_path) + "/"
    info = parse_image_setting_files([root + n for n in OWN + OWN + TILE])
    assert len(calls) == 1 and len(calls[0]) == 3
    assert [i["StageInterval"] for i in info] == [0.5, -0.5, 0.5, -0.5, 0.5, 0.5]
