"""pk5d-master: the compiled master's call convention, output lines and exit statuses.

The compiled master (mccMaster.m; observed on Linux with MATLAB Runtime R2024b) prints
``Hostname: <host> ``, ``Function name: <f>`` and ``Variables: <strings> `` first, the task's
output, then ``Elapsed time is <s> seconds.`` and an empty line; it exits 0, or 249 after an
uncaught error, whose text (``Error using ...``, message, identifier) goes to stderr.
"""

from __future__ import annotations

import importlib.util
import io
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from petakit5d.cli import master

HERE = Path(__file__).parent


def _toy():
    if "pk5d_cli_toy" in sys.modules:
        return sys.modules["pk5d_cli_toy"]
    spec = importlib.util.spec_from_file_location("pk5d_cli_toy", HERE / "toy.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["pk5d_cli_toy"] = mod
    spec.loader.exec_module(mod)
    return mod


TOY = _toy()
PATHS = "{'/mnt/projects/embryo_2024/a','/mnt/projects/embryo_2024/b c'}"


def run(*argv, env=None, monkeypatch=None):
    out, err = io.StringIO(), io.StringIO()
    TOY.CALLS.clear()
    code = master.main(list(argv), registry=TOY.REG, stdout=out, stderr=err)
    return code, out.getvalue(), err.getvalue()


def test_header_lines_are_the_compiled_masters():
    code, out, err = run("XR_toy_dataset", PATHS, "[1,2,3,10,20,30]", "save16bit", "false")
    assert code == 0, err
    lines = out.split("\n")
    assert lines[0] == f"Hostname: {master.hostname()} "
    assert lines[1] == "Function name: XR_toy_dataset"
    assert lines[2] == f"Variables: {PATHS} [1,2,3,10,20,30] save16bit false "
    assert re.fullmatch(r"Elapsed time is [0-9]+\.[0-9]{6} seconds\.", lines[3])  # toc
    assert lines[4:] == ["", ""]


def test_header_without_arguments_has_two_blanks():
    assert master.header("f", []) == "Function name: f\nVariables:  \n"


def test_values_reach_the_function_in_pythons_form():
    code, _, err = run(
        "XR_toy_dataset", PATHS, "[1,2,3,10,20,30]", "BATCH", "[64,128,32]", "save16BIT",
        "false", "channelPatterns", "{'CamA_ch0'}",
    )  # fmt: skip
    assert code == 0, err
    (name, got) = TOY.CALLS[-1]
    assert name == "toy_dataset"
    assert got["data_paths"] == ("/mnt/projects/embryo_2024/a", "/mnt/projects/embryo_2024/b c")
    assert got["input_bbox"] == (2, 0, 1, 30, 10, 20)
    assert got["batch_size"] == (32, 64, 128)
    assert got["save_16bit"] is False
    assert got["channel_patterns"] == ("CamA_ch0",)


def test_text_path_stays_text_and_last_value_wins():
    code, _, err = run(
        "XR_toy_dataset", "/mnt/projects/x", "[1,1,1,2,2,2]", "resultDirName", "a",
        "resultdirname", "b",
    )  # fmt: skip
    assert code == 0, err
    got = TOY.CALLS[-1][1]
    assert got["data_paths"] == "/mnt/projects/x"
    assert got["result_dir_name"] == "b"


def test_no_function_name():
    code, out, err = run()
    assert code == master.EXIT_ERROR == 249
    assert out == f"Hostname: {master.hostname()} \n"
    assert err == "Not enough input arguments.\n\nMATLAB:minrhs\n"


def test_unknown_function():
    code, out, err = run("XR_no_such", "a")
    assert code == 249
    assert out.endswith("Function name: XR_no_such\nVariables: a \n")
    assert err.startswith(
        "Error using mccMaster\nThe parser for function XR_no_such does not exist!"
    )


def test_function_names_are_exact():
    code, _, err = run("xr_toy_dataset", PATHS, "[1,1,1,2,2,2]")
    assert code == 249
    assert "Did you mean XR_toy_dataset?" in err


@pytest.mark.parametrize(
    ("argv", "identifier"),
    [
        ([PATHS], "MATLAB:badsubscript"),
        ([PATHS, "[1,1,1,2,2,2]", "bachSize", "[1,2,3]"], "MATLAB:InputParser:UnmatchedParameter"),
        ([PATHS, "[1,1,1,2,2,2]", "save", "true"], "MATLAB:InputParser:AmbiguousParameter"),
        ([PATHS, "[1,1,1,2,2,2]", "input", "x"], "MATLAB:InputParser:AmbiguousParameter"),
        ([PATHS, "[1,1,1,2,2,2]", "saveZarr"], "MATLAB:InputParser:ParamMissingValue"),
        ([PATHS, "[1,1,1,2,2,2]", "saveZarr", "1"], "MATLAB:InputParser:ArgumentFailedValidation"),
        ([PATHS, "[1,1,1,2,2,2]", "channelPatterns", "CamA"],
         "MATLAB:InputParser:ArgumentFailedValidation"),
    ],
)  # fmt: skip
def test_argument_errors_exit_249_with_matlabs_identifier(argv, identifier):
    code, _, err = run("XR_toy_dataset", *argv)
    assert code == 249
    assert err.startswith("Error using XR_toy_dataset\n")
    assert err.endswith(f"\n\n{identifier}\n")
    assert TOY.CALLS == []


def test_too_few_arguments_message():
    code, _, err = run("XR_toy_dataset", PATHS)
    assert "Index exceeds the number of array elements. Index must not exceed 1." in err


def test_function_error_exits_249_with_traceback_and_identifier():
    code, out, err = run("XR_toy_dataset", PATHS, "[1,1,1,2,2,2]", "fail", "true")
    assert code == 249
    assert "Traceback (most recent call last)" in err
    assert err.endswith("Error using XR_toy_dataset\nToyError: the toy failed\n\npkpy:toy:failed\n")
    assert "Elapsed time" not in out


def test_without_name_value_pairs_extra_strings_are_ignored():
    code, _, err = run("toy_positional_task", "/mnt/in.mat", "0.25", "[3]", "surplus", "x")
    assert code == 0, err
    got = TOY.CALLS[-1][1]
    assert got == {"input_fullpath": "/mnt/in.mat", "threshold": 0.25, "frame": (3,)} or got == {
        "input_fullpath": "/mnt/in.mat",
        "threshold": 0.25,
        "frame": 3,
    }


def test_optional_positional_is_required_on_the_command_line():
    code, _, err = run("toy_positional_task", "/mnt/in.mat", "0.25")
    assert code == 249
    assert "MATLAB:badsubscript" in err


def test_dry_run_binds_without_calling(monkeypatch):
    monkeypatch.setenv(master.DRY_RUN_ENV, "1")
    code, out, err = run("XR_toy_dataset", PATHS, "[1,1,1,2,2,2]", "batchSize", "[4,5,6]")
    assert code == 0, err
    assert TOY.CALLS == []
    assert "Dry run" in out and "batch_size=(6, 4, 5)," in out
    monkeypatch.setenv(master.DRY_RUN_ENV, "0")
    code, _, _ = run("XR_toy_dataset", PATHS, "[1,1,1,2,2,2]")
    assert code == 0 and len(TOY.CALLS) == 1


def test_help():
    code, out, _ = run("--help")
    assert code == 0 and out.startswith("usage: pk5d-master")


def test_hostname_follows_the_environment(monkeypatch):
    if os.name == "nt":
        pytest.skip("POSIX: mccMaster prints $HOSTNAME")
    monkeypatch.setenv("HOSTNAME", "node12")
    assert master.hostname() == "node12"
    monkeypatch.delenv("HOSTNAME")
    assert master.hostname() == os.uname().nodename


def _python_env():
    env = dict(os.environ)
    src = str(Path(master.__file__).resolve().parents[2])
    env["PYTHONPATH"] = os.pathsep.join([src, env.get("PYTHONPATH", "")])
    return env


def test_process_exit_status_and_streams():
    p = subprocess.run(
        [sys.executable, "-m", "petakit5d.cli.master", "XR_no_such_function", "x"],
        capture_output=True,
        text=True,
        env=_python_env(),
        timeout=120,
    )
    assert p.returncode == 249
    assert p.stdout.splitlines()[1:] == ["Function name: XR_no_such_function", "Variables: x "]
    assert "The parser for function XR_no_such_function does not exist!" in p.stderr


def test_no_heavy_imports_before_dispatch():
    """Until the function's module is imported, the master loads neither NumPy nor SciPy."""
    code = (
        "import sys, io\n"
        "from petakit5d.cli import master\n"
        "master.main(['XR_no_such_function'], stdout=io.StringIO(), stderr=io.StringIO())\n"
        "heavy = sorted({m.split('.')[0] for m in sys.modules} & "
        "{'numpy', 'scipy', 'skimage', 'cupy', 'h5py', 'cpptiff', 'cppzarr', 'argparse'})\n"
        "print(heavy)\n"
    )
    p = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, env=_python_env(), timeout=120
    )
    assert p.returncode == 0, p.stderr
    assert p.stdout.strip() == "[]"


def test_console_scripts_are_declared():
    text = (HERE.parents[1] / "pyproject.toml").read_text()
    assert 'pk5d = "petakit5d.cli.main:main"' in text
    assert 'pk5d-master = "petakit5d.cli.master:main"' in text


def test_pairs_naming_positional_arguments_are_ignored():
    """MATLAB quirk: the mcc parser's inputParser accepts a pair for a required argument, but
    the parser passes the function its positional variable."""
    code, _, err = run(
        "XR_toy_dataset", PATHS, "[1,1,1,2,2,2]", "inputBbox", "[1,1,1,5,5,5]", "DATAP", "/x",
        "batchSize", "[1,2,3]",
    )  # fmt: skip
    assert code == 0, err
    got = TOY.CALLS[-1][1]
    assert got["input_bbox"] == (0, 0, 0, 2, 2, 2)
    assert got["data_paths"][0] == "/mnt/projects/embryo_2024/a"
    assert got["batch_size"] == (3, 1, 2)
    code, _, err = run("XR_toy_dataset", PATHS, "[1,1,1,2,2,2]", "inputBbox")
    assert code == 249 and err.endswith("MATLAB:InputParser:ParamMissingValue\n")


def test_function_level_conversion_is_applied():
    code, _, err = run("toy_hooked_task", "25")
    assert code == 0, err
    assert TOY.CALLS[-1] == ("toy_hooked", {"level": 0.25})
