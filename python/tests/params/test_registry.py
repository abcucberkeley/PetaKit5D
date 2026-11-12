"""The registry: discovery by source scan, lazy imports, the package's declarations."""

from __future__ import annotations

import subprocess
import sys
import textwrap

import pytest

from petakit5d.params import REGISTRY, Registry
from petakit5d.params.planned import PLANNED

AREA = '''
from petakit5d.params import Param, not_ported

REG_CALLS = []


def make(reg):
    @reg.register(
        "XR_fake_tool",
        Param("dataPaths", "@(x) ischar(x) || iscell(x)", "Data folders."),
        Param("Overwrite", "@islogical", "Overwrite results."),
    )
    def fake_tool(data_paths, *, overwrite=False):
        """A fake tool."""
        REG_CALLS.append((data_paths, overwrite))
        return data_paths, overwrite

    return fake_tool
'''


@pytest.fixture
def fake_package(tmp_path, monkeypatch):
    """A package 'fakepk' with an area module that registers on import."""
    pkg = tmp_path / "fakepk"
    (pkg / "area").mkdir(parents=True)
    (pkg / "__init__.py").write_text("")
    (pkg / "area" / "__init__.py").write_text("")
    (pkg / "area" / "impl.py").write_text(
        AREA
        + textwrap.dedent(
            """
            from fakepk import REG
            fake_tool = make(REG)
            """
        )
    )
    (pkg / "area" / "gone.py").write_text(
        'from fakepk import REG\nREG.not_ported("XR_fake_old", "legacy", "area")\n'
    )
    (pkg / "__init__.py").write_text(
        "from pathlib import Path\nfrom petakit5d.params import Registry\n"
        "REG = Registry('fakepk', Path(__file__).parent)\n"
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    for m in [m for m in sys.modules if m.startswith("fakepk")]:
        monkeypatch.delitem(sys.modules, m)
    import fakepk

    return fakepk


def test_find_imports_only_the_declaring_module(fake_package):
    reg = fake_package.REG
    # register(" is in the text of impl.py through make(); the scan reads the literal name
    assert reg.index()["XR_fake_tool"] == ("fakepk.area.impl", "register")
    assert reg.index()["XR_fake_old"] == ("fakepk.area.gone", "not_ported")
    assert "fakepk.area.impl" not in sys.modules
    spec = reg.find("XR_fake_tool")
    assert "fakepk.area.impl" in sys.modules and "fakepk.area.gone" not in sys.modules
    assert spec.call_matlab("/mnt/x", "overwrite", True) == ("/mnt/x", True)
    assert reg.names() == ["XR_fake_tool"]
    with pytest.raises(LookupError, match="not ported: legacy"):
        reg.find("XR_fake_old")
    with pytest.raises(LookupError, match="Did you mean XR_fake_tool"):
        reg.find("XR_fake_tol")
    assert [s.matlab_name for s in reg.specs()] == ["XR_fake_tool"]
    assert list(reg.not_ported_all()) == ["XR_fake_old"]
    assert reg.find_python(spec.python_name) is spec
    assert reg.find_python(sys.modules["fakepk.area.impl"].fake_tool) is spec


def test_package_scan_agrees_with_imports():
    """Every declaration the scan finds is registered by its module, and nothing else."""
    REGISTRY.load_all()
    scanned = {k for k, (_, kind) in REGISTRY.index().items() if kind == "register"}
    assert scanned == set(REGISTRY._specs)
    unported = {k for k, (_, kind) in REGISTRY.index().items() if kind == "not_ported"}
    assert unported == set(REGISTRY.not_ported_all())


def test_package_declarations():
    np_ = REGISTRY.not_ported_all()
    assert "XR_stitching_frame_v1" in np_ and "XR_cppDeconFrame3D" in np_
    with pytest.raises(LookupError, match=r"not ported yet \(planned as petakit5d.io.read_tiff"):
        REGISTRY.find("readtiff")
    import petakit5d.matlab as pm

    with pytest.raises(AttributeError, match="not ported"):
        pm.XR_stitching_frame_v1  # noqa: B018
    planned = [m for m, *_ in PLANNED]
    assert len(planned) == len(set(planned))
    assert not set(planned) & set(np_)


def test_planned_names_follow_the_rule():
    for _matlab, py, task, phase in PLANNED:
        assert py.startswith("petakit5d."), py
        name = py.rsplit(".", 1)[-1]
        assert name == name.lower() and not name.startswith("xr_"), py
        assert phase in (1, 2, 3) and task[0].isupper()


def test_import_is_light():
    code = (
        "import sys, petakit5d.params, petakit5d.matlab; "
        "bad = [m for m in ('numpy', 'scipy', 'petakit5d._ext') if m in sys.modules]; "
        "print(bad)"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "[]"


def test_isolated_registries_do_not_mix():
    a, b = Registry(), Registry()
    from petakit5d.params import Param

    @a.register("only_in_a", Param("x", "@isnumeric"))
    def f(x):
        return x

    assert "only_in_a" in a and "only_in_a" not in b and "only_in_a" not in REGISTRY
