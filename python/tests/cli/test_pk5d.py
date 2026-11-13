"""pk5d: snake_case options with JSON values, list, doc, parse and master."""

from __future__ import annotations

import importlib.util
import io
import json
import sys
from pathlib import Path

import pytest

from petakit5d.cli import main as cli

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


def run(*argv):
    out, err = io.StringIO(), io.StringIO()
    TOY.CALLS.clear()
    code = cli.main(list(argv), registry=TOY.REG, stdout=out, stderr=err)
    return code, out.getvalue(), err.getvalue()


def test_python_values_and_json():
    code, out, err = run(
        "toy_dataset", '["/mnt/projects/a", "/mnt/projects/b"]', "[0, 1, 2, 10, 20, 30]",
        "--batch-size", "[16, 32, 64]", "--save_16bit=false", "--save-zarr",
    )  # fmt: skip
    assert code == 0, err
    got = TOY.CALLS[-1][1]
    assert got["data_paths"] == ("/mnt/projects/a", "/mnt/projects/b")
    assert got["input_bbox"] == (0, 1, 2, 10, 20, 30)  # Python's form: no conversion
    assert got["batch_size"] == (16, 32, 64)
    assert got["save_16bit"] is False and got["save_zarr"] is True
    assert json.loads(out) == {"batch": [16, 32, 64]}


def test_by_matlab_name_and_plain_text_values():
    code, _, err = run("XR_toy_dataset", "/mnt/projects/a", "[0,0,0,1,1,1]", "--result-dir-name",
                       "out dir")  # fmt: skip
    assert code == 0, err
    got = TOY.CALLS[-1][1]
    assert got["data_paths"] == "/mnt/projects/a" and got["result_dir_name"] == "out dir"


def test_double_dash_ends_options():
    code, _, err = run("toy_dataset", "--", "--odd-folder", "[0,0,0,1,1,1]")
    assert code == 0, err
    assert TOY.CALLS[-1][1]["data_paths"] == "--odd-folder"


@pytest.mark.parametrize(
    ("argv", "text"),
    [
        (["toy_dataset", "/a", "[0,0,0,1,1,1]", "--bach-size", "[1,2,3]"],
         "unknown option --bach-size (did you mean --batch_size?)"),
        (["toy_dataset", "/a", "[0,0,0,1,1,1]", "--batchSize", "[1,2,3]"],
         "(did you mean --batch_size?)"),
        (["toy_dataset", "/a"], "missing required arguments: 'input_bbox'"),
        (["toy_dataset", "/a", "[0,0,0,1,1,1]", "--save-zarr", "1"], "save"),
        (["no_such_thing"], "no registered function no_such_thing"),
        (["toy_dataset", "/a", "[0,0,0,1,1,1]", "--fail", "--fail"], "given twice"),
    ],
)  # fmt: skip
def test_usage_errors_exit_2(argv, text):
    code, _, err = run(*argv)
    assert code == 2
    assert err.startswith("pk5d: error: ") and text in err
    assert TOY.CALLS == []


def test_function_error_exits_1_with_traceback():
    code, _, err = run("toy_dataset", "/a", "[0,0,0,1,1,1]", "--fail")
    assert code == 1
    assert "ToyError: the toy failed" in err


def test_list_and_doc():
    code, out, _ = run("list")
    assert code == 0
    assert "toy_dataset     XR_toy_dataset       Record a call" in out and "toy_positional" in out
    code, out, _ = run("doc", "XR_toy_dataset")
    assert code == 0
    assert out.startswith("toy_dataset  (")
    assert (
        "--batch-size  (default [1024, 1024, 1024]; MATLAB batchSize; converted: size_yxz)" in out
    )
    assert "data_paths  (required; MATLAB dataPaths)" in out
    assert "pk5d-master XR_toy_dataset" in out
    assert run("doc")[0] == 2


def test_parse_shows_both_forms():
    code, out, err = run("parse", "XR_toy_dataset", "{'/a'}", "[1,2,3,10,20,30]", "batch",
                         "[1,2,3]")  # fmt: skip
    assert code == 0, err
    d = json.loads(out)
    assert d["python_arguments"]["batch_size"] == [3, 1, 2]
    assert d["matlab_values"]["batchSize"] == [[1.0, 2.0, 3.0]]
    assert d["matlab_values"]["inputBbox"] == [[1.0, 2.0, 3.0, 10.0, 20.0, 30.0]]
    code, _, err = run("parse", "XR_toy_dataset", "{'/a'}")
    assert code == 2 and "MATLAB:badsubscript" in err


def test_master_subcommand():
    code, out, _ = run("master", "XR_toy_dataset", "{'/a'}", "[1,1,1,2,2,2]")
    assert code == 0 and "Function name: XR_toy_dataset" in out
    assert len(TOY.CALLS) == 1


def test_mcc_shim_writes_an_executable(tmp_path):
    p = tmp_path / "run_mccMaster.sh"
    code, out, _ = run("mcc-shim", str(p), "--command", "echo")
    assert code == 0 and p.exists()
    assert "pk5d_master='echo'" in p.read_text()


def test_help_and_version():
    assert run("--help")[0] == 0
    assert run()[0] == 2
    code, out, _ = run("--version")
    assert code == 0 and out.startswith("pk5d (petakit5d ")


def test_parse_value():
    assert cli.parse_value("[1, 2]") == [1, 2]
    assert cli.parse_value("null") is None
    assert cli.parse_value("/mnt/x y") == "/mnt/x y"
    assert cli.parse_value('"123"') == "123"


def test_real_registry_doc():
    """The package registry: a Phase 1 function by its Python and MATLAB names."""
    out = cli.doc_text("XR_crop_dataset")
    assert "crop_dataset" in out and "--batch-size" in out
    assert cli.find_spec("crop_dataset").matlab_name == "XR_crop_dataset"
