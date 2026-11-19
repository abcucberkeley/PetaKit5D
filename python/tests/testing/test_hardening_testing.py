"""Tests the mutation campaign asked for (E06, docs/hardening.md): the test harness itself."""

import argparse
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from petakit5d.testing import cases as T
from petakit5d.testing import compare as C
from petakit5d.testing import fixtures, oracle


def test_an_oracle_run_with_a_nondeterministic_input_stops_before_matlab(monkeypatch, capsys):
    rng = np.random.default_rng()
    case = T.Case("random", matlab="ones", args=[T.In("x")], inputs={"x": T.Value(rng.random)})

    def no_store(k):
        raise RuntimeError("the run went on to the store and MATLAB")

    store = SimpleNamespace(has=no_store)
    monkeypatch.setattr(oracle, "_store_and_oracle", lambda allow: (store, None))
    monkeypatch.setattr(oracle, "collect", lambda args: [("tests/m.py", case)])
    ns = argparse.Namespace(allow_other_oracle=False, pytest_args=[], device="cpu")
    assert oracle.cmd_run(ns) == 2
    out = capsys.readouterr().out
    assert "is not deterministic" in out and "no MATLAB run" in out


def test_gpu_skip_reason_when_cuda_visible_devices_hides_them(monkeypatch):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    assert fixtures.gpu_skip_reason.__wrapped__() == "CUDA_VISIBLE_DEVICES hides all GPUs"
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "-1")
    assert fixtures.gpu_skip_reason.__wrapped__() == "CUDA_VISIBLE_DEVICES hides all GPUs"


def test_a_value_of_no_known_kind_cannot_be_hashed():
    with pytest.raises(TypeError, match="cannot hash object"):
        C.vhash(object())


def test_value_files_read_their_references(tmp_path):
    v = C.to_mvalue(np.arange(6.0).reshape(2, 3))
    inner = tmp_path / "in.pkv"
    C.write_pkv(inner, v)
    outer = tmp_path / "out.pkv"
    C.write_pkv(outer, C.MCell((1, 1), [C.MRef(inner, C.vhash(v))]))
    got = C.read_pkv(outer)
    assert C.vhash(got.items[0]) == C.vhash(v)


def test_a_dict_as_a_map_of_another_length_takes_no_template():
    m = C.MMap(
        "char", "double", C.MCell((1, 1), [C.mchar("x")]), C.MCell((1, 1), [C.to_mvalue(1.0)])
    )
    got = C._dict_to_map({"a": 1.0, "b": 2.0}, m)
    assert [C.vhash(k) for k in got.keys.items] == [C.vhash(C.mchar(k)) for k in ("a", "b")]
    assert [C.vhash(v) for v in got.values.items] == [C.vhash(C.to_mvalue(x)) for x in (1.0, 2.0)]


def test_sparse_matrices_of_v73_mat_files_are_read(tmp_path):
    import scipy.sparse as sp

    from petakit5d.io.matfile import write_mat

    s = sp.csc_matrix(np.array([[0.0, 2.5, 0.0], [1.0, 0.0, 0.0]]))
    path = tmp_path / "s.mat"
    write_mat(path, {"s": s}, version="7.3")
    got = C.read_mat(path).items[0]["s"]
    assert isinstance(got, C.MSparse)
    assert tuple(got.size) == (2, 3)
    assert C.vhash(got) == C.vhash(C.to_mvalue(s))


def test_a_custom_comparison_reports_its_problems():
    spec = C.Custom(lambda r, a: ["they differ"] if C.vhash(r) != C.vhash(a) else [], "mine")
    ref = C.to_mvalue(np.ones(3))
    assert C.compare(ref, np.ones(3), spec).ok
    rep = C.compare(ref, np.zeros(3), spec)
    assert not rep.ok and "they differ" in rep.text("x")


def test_a_missing_json_list_is_empty():
    assert oracle._seq(None) == []
    assert oracle._seq(3) == [3]
    assert oracle._seq([1, 2]) == [1, 2]


@pytest.mark.skipif(os.name == "nt", reason="symbolic links")
def test_a_local_python_copy_points_its_links_at_the_copied_interpreter(tmp_path, monkeypatch):
    from petakit5d.testing import bench

    base, venv, dest = tmp_path / "base", tmp_path / "venv", tmp_path / "dest"
    (base / "bin").mkdir(parents=True)
    (base / "bin" / "python3.12").write_text("interpreter")
    (venv / "bin").mkdir(parents=True)
    (venv / "pyvenv.cfg").write_text(f"home = {base}/bin\nexecutable = {base}/bin/python3.12\n")
    os.symlink(base / "bin" / "python3.12", venv / "bin" / "python")
    monkeypatch.setattr(sys, "base_prefix", str(base))
    monkeypatch.setattr(sys, "prefix", str(venv))
    monkeypatch.setattr(sys, "executable", str(venv / "bin" / "python"))
    monkeypatch.delenv("PYTHONPATH", raising=False)
    monkeypatch.delenv("PETAKIT5D_EXT_DIR", raising=False)
    out = bench.local_python_copy(str(dest))
    assert out["python"] == str(dest / "venv" / "bin" / "python")
    link = os.readlink(dest / "venv" / "bin" / "python")
    assert link == str(dest / "base" / "bin" / "python3.12")
    cfg = (dest / "venv" / "pyvenv.cfg").read_text()
    assert f"home = {dest / 'base' / 'bin'}" in cfg and "executable" not in cfg
