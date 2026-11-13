"""The run_mccMaster.sh shim: MATLAB's launcher convention (Runtime folder first, ignored; each
argument quoted and evaluated by the shell), running pk5d-master. And the emulations of the
callers' quoting (tests/cli/callers.py) that the argument corpus is made with."""

from __future__ import annotations

import importlib.util
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from petakit5d.cli import master

HERE = Path(__file__).parent


def _callers():
    if "pk5d_cli_callers" in sys.modules:
        return sys.modules["pk5d_cli_callers"]
    spec = importlib.util.spec_from_file_location("pk5d_cli_callers", HERE / "callers.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pk5d_cli_callers"] = mod
    spec.loader.exec_module(mod)
    return mod


C = _callers()
posix = pytest.mark.skipif(
    os.name == "nt" or shutil.which("sh") is None, reason="the shim is a POSIX sh script"
)

RECORDER = "import json, sys\nprint(json.dumps(sys.argv[1:]))\n"


def _env():
    env = dict(os.environ)
    src = str(Path(master.__file__).resolve().parents[2])
    env["PYTHONPATH"] = os.pathsep.join([src, env.get("PYTHONPATH", "")])
    return env


@pytest.fixture
def shim(tmp_path):
    rec = tmp_path / "recorder.py"
    rec.write_text(RECORDER)
    path = tmp_path / "PetaKit5D" / "mcc" / "linux" / "run_mccMaster.sh"
    path.parent.mkdir(parents=True)
    master.write_shim(path, f"{shlex.quote(sys.executable)} {shlex.quote(str(rec))}")
    return path


def _run(shim, *argv, env=None):
    p = subprocess.run(
        ["sh", str(shim), *argv], capture_output=True, text=True, env=env or _env(), timeout=120
    )
    return p


TOKENS = [
    "XR_crop_dataset",
    "{'/mnt/projects/a','/mnt/projects/b c'}",
    "[1,1,1,100,200,50]",
    "resultDirName",
    "my result",
    "jobLogDir",
    "",
    "channelPatterns",
    "{'it''s'}",
]


@posix
def test_shim_drops_the_runtime_and_passes_the_arguments(shim):
    p = _run(shim, "/opt/MATLAB_Runtime/R2024b", *TOKENS)
    assert p.returncode == 0, p.stderr
    assert json.loads(p.stdout) == TOKENS == C.launcher_eval(TOKENS)


@posix
def test_shim_undoes_pypetakit5ds_escapes(shim):
    tokens = ["XR_MIP_wrapper", "{'/mnt/$HOME/a\\b'}", 'say "hi"', "`date`", "x$y"]
    p = _run(shim, "/opt/rt", *[C.escape_for_launcher(t) for t in tokens])
    assert p.returncode == 0, p.stderr
    assert json.loads(p.stdout) == tokens
    assert C.launcher_eval([C.escape_for_launcher(t) for t in tokens]) == tokens


@posix
def test_shim_makes_the_guis_doubled_quotes_empty(shim):
    a = C.GuiArgs().data_paths(["/mnt/projects/a b"]).opt("jobLogDir", "char", "")
    a.opt("pad", "bool", True)
    cmd = C.gui_fix_empty(f'"{shim}" "/opt/rt" XR_crop_dataset {a.text}')
    words = C.qt_split_command(cmd)
    assert words[:3] == [str(shim), "/opt/rt", "XR_crop_dataset"]
    assert words[4:6] == ["jobLogDir", '""']  # what QProcess hands the launcher
    p = _run(shim, *words[1:])
    assert p.returncode == 0, p.stderr
    want = ["XR_crop_dataset", "{'/mnt/projects/a b'}", "jobLogDir", "", "pad", "true"]
    assert json.loads(p.stdout) == want == C.gui_master_argv("XR_crop_dataset", a)


@posix
def test_shim_usage_without_runtime(shim):
    p = _run(shim)
    assert p.returncode == 0
    assert p.stdout.splitlines()[0] == "Usage:"
    assert "<deployedMCRroot> args" in p.stdout


@posix
def test_shim_runs_pk5d_master_by_default(tmp_path):
    path = master.write_shim(tmp_path / "run_mccMaster.sh")
    text = Path(path).read_text()
    assert shlex.quote(sys.executable) in text and "-m petakit5d.cli.master" in text
    assert os.access(path, os.X_OK)
    p = _run(path, "/opt/rt", "XR_no_such_function", "a b")
    assert p.returncode == 249
    assert p.stdout.splitlines()[1:] == ["Function name: XR_no_such_function", "Variables: a b "]


@posix
def test_shim_command_from_the_environment(shim, tmp_path):
    rec = tmp_path / "other.py"
    rec.write_text("import sys\nprint('other', *sys.argv[1:])\n")
    env = _env()
    env[master.SHIM_ENV] = f"{shlex.quote(sys.executable)} {shlex.quote(str(rec))}"
    p = _run(shim, "/opt/rt", "f", "x")
    assert json.loads(p.stdout) == ["f", "x"]
    p = _run(shim, "/opt/rt", "f", "x", env=env)
    assert p.stdout.strip() == "other f x"


def test_qt_split_command():
    assert C.qt_split_command('a "b c" d') == ["a", "b c", "d"]
    assert C.qt_split_command('x """""" y') == ["x", '""', "y"]
    assert C.qt_split_command('x "" y') == ["x", "y"]  # an empty word is dropped
    assert C.qt_split_command('a"b c"d') == ["ab cd"]


def test_launcher_eval():
    assert C.launcher_eval(["a b", "", '\\"q\\"']) == ["a b", "", '"q"']
    assert C.launcher_eval(['""']) == [""]
    with pytest.raises(ValueError):
        C.launcher_eval(["$HOME"])


def test_gui_arguments():
    a = C.GuiArgs().data_paths(["/p/a", "/p/b"]).array(["1", "2"], False)
    a.opt("channelPatterns", "channels", (None, "CamA, CamB,"))
    a.opt("processFunPath", "cells;", ["/t/a.txt", "/t/b.txt"]).opt("dz", "scalar", C.ftos(0.5))
    assert C.gui_master_argv("f", a) == [
        "f", "{'/p/a','/p/b'}", "[1,2]", "channelPatterns", "{'CamA',' CamB'}",
        "processFunPath", "{'/t/a.txt';'/t/b.txt'}", "dz", "0.500000",
    ]  # fmt: skip


def test_pypetakit5d_arguments():
    argv = C.pypetakit5d_master_argv(
        "XR_crop_dataset", [["/p/a"], [1, 2, 3, 4, 5, 6]], pad=True, channelPatterns=["it's"],
        cpusPerTask=1.5, uuid="",
    )  # fmt: skip
    assert argv[:2] == ["{'/p/a'}", "[1,2,3,4,5,6]"]
    i = argv.index("pad")
    assert argv[i + 1] == "true"
    assert argv[argv.index("channelPatterns") + 1] == "{'it''s'}"
    assert argv[argv.index("cpusPerTask") + 1] == "1.5"
    assert "uuid" not in argv  # empty text is left out
    assert "lastStartCoords" not in argv  # an empty array too
