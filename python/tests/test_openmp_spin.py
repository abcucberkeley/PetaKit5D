"""The OpenMP spin that importing petakit5d sets (docs/decisions/0007-openmp-spin.md):
GOMP_SPINCOUNT when the user set neither it nor OMP_WAIT_POLICY, in the environment before any
OpenMP runtime is loaded, so that every libgomp of the process and its children read it, with a
mark by which programs that aren't the port (MATLAB) can be started without it."""

import json
import os
import subprocess
import sys

import pytest

import petakit5d
from petakit5d import _ext

SPIN = petakit5d._OPENMP_SPINCOUNT
MARK = petakit5d._OPENMP_MARK
SET_BY_TESTS = ("GOMP_SPINCOUNT", "OMP_WAIT_POLICY", "OMP_DISPLAY_ENV", MARK)
SHOW = f"import os, petakit5d; print(os.environ.get('GOMP_SPINCOUNT'), os.environ.get('{MARK}'))"


def _env(**extra: str) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if k not in SET_BY_TESTS}
    env.update(extra)
    return env


def _run(code: str, tmp_path, env: dict[str, str], first_on_path: str | None = None):
    if first_on_path:
        old = env.get("PYTHONPATH")
        path = [first_on_path, old] if old else [first_on_path]
        env = dict(env, PYTHONPATH=os.pathsep.join(path))
    return subprocess.run(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=True,
        timeout=120,
    )


def _fake_io_packages(tmp_path) -> str:
    # stand-ins named as the I/O wheels, which print the variable while they are imported
    root = tmp_path / "fake"
    for name in ("cppzarr", "cpptiff"):
        (root / name).mkdir(parents=True)
        (root / name / "__init__.py").write_text(
            "import os\nprint(__name__, os.environ.get('GOMP_SPINCOUNT'))\n"
        )
    return str(root)


def test_the_spin_is_set_when_the_user_set_none():
    env = {}
    assert petakit5d._set_openmp_spin(env) is True
    assert env == {"GOMP_SPINCOUNT": SPIN, MARK: SPIN}
    assert int(SPIN) > 0


@pytest.mark.parametrize(
    "user",
    [
        {"GOMP_SPINCOUNT": "300000"},
        {"GOMP_SPINCOUNT": ""},
        {"OMP_WAIT_POLICY": "ACTIVE"},
        {"OMP_WAIT_POLICY": "PASSIVE"},
        {"GOMP_SPINCOUNT": "0", "OMP_WAIT_POLICY": "passive"},
    ],
)
def test_a_users_setting_is_kept(user):
    env = dict(user, OTHER="1")
    assert petakit5d._set_openmp_spin(env) is False
    assert env == dict(user, OTHER="1")


def test_importing_petakit5d_sets_it(tmp_path):
    assert _run(SHOW, tmp_path, _env()).stdout.split() == [SPIN, SPIN]


@pytest.mark.parametrize(
    "user, expected",
    [({"GOMP_SPINCOUNT": "1234"}, "1234"), ({"OMP_WAIT_POLICY": "PASSIVE"}, "None")],
)
def test_importing_petakit5d_keeps_the_users_setting(tmp_path, user, expected):
    assert _run(SHOW, tmp_path, _env(**user)).stdout.split() == [expected, "None"]


@pytest.mark.parametrize(
    "env, expected",
    [
        ({"GOMP_SPINCOUNT": SPIN, MARK: SPIN, "OTHER": "1"}, {"OTHER": "1"}),
        ({"GOMP_SPINCOUNT": SPIN}, {"GOMP_SPINCOUNT": SPIN}),  # the user's own: no mark
        ({"GOMP_SPINCOUNT": "5000", MARK: SPIN}, {"GOMP_SPINCOUNT": "5000"}),  # changed since
        ({"OMP_WAIT_POLICY": "PASSIVE"}, {"OMP_WAIT_POLICY": "PASSIVE"}),
        ({MARK: SPIN}, {}),
        ({}, {}),
    ],
)
def test_a_program_can_be_started_without_the_ports_spin(env, expected):
    before = dict(env)
    assert petakit5d._without_openmp_spin(env) == expected
    assert env == before  # a copy


@pytest.mark.parametrize("user, expected", [({}, "None"), ({"GOMP_SPINCOUNT": "1234"}, "1234")])
def test_a_program_started_without_it_sees_the_users_environment(tmp_path, user, expected):
    # a port process starts a port child, which starts a program that isn't the port
    show, child = tmp_path / "show.py", tmp_path / "child.py"
    show.write_text(
        f"import os\nprint(os.environ.get('GOMP_SPINCOUNT'), os.environ.get('{MARK}'))\n"
    )
    child.write_text(
        "import os, subprocess, sys, petakit5d\n"
        f"subprocess.run([sys.executable, {str(show)!r}],"
        " env=petakit5d._without_openmp_spin(os.environ), check=True)\n"
    )
    code = (
        "import subprocess, sys, petakit5d\n"
        f"subprocess.run([sys.executable, {str(child)!r}], check=True)\n"
    )
    assert _run(code, tmp_path, _env(**user)).stdout.split() == [expected, "None"]


@pytest.mark.parametrize("user, expected", [({}, SPIN), ({"GOMP_SPINCOUNT": "1234"}, "1234")])
def test_the_io_packages_see_it_when_imported_after_petakit5d(tmp_path, user, expected):
    code = "import petakit5d, cppzarr, cpptiff\n"
    out = _run(code, tmp_path, _env(**user), _fake_io_packages(tmp_path)).stdout.splitlines()
    assert out == [f"cppzarr {expected}", f"cpptiff {expected}"]


def test_it_is_in_place_before_the_first_runtime_loads(tmp_path):
    # an import hook notes the variable when each module that brings an OpenMP runtime
    # (an extension module of ours, cppzarr, cpptiff) starts to load
    code = """
import importlib.abc, json, os, sys
seen = {}
class Spy(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path=None, target=None):
        if name.startswith("petakit5d._ext.") or name in ("cppzarr", "cpptiff"):
            seen.setdefault(name, os.environ.get("GOMP_SPINCOUNT"))
        return None
sys.meta_path.insert(0, Spy())
import petakit5d
maps = ""
if os.path.exists("/proc/self/maps"):
    with open("/proc/self/maps") as f:
        maps = f.read()
loaded = sorted({ln.split()[-1] for ln in maps.splitlines()
                 if any(k in ln.rsplit("/", 1)[-1] for k in ("gomp", "libomp", "iomp", "vcomp"))})
from petakit5d import _ext
groups = _ext.groups()
for g in groups:
    _ext.load(g)
for name in ("cppzarr", "cpptiff"):
    try:
        __import__(name)
    except ImportError:
        pass
print(json.dumps({"seen": seen, "before": loaded, "groups": groups}))
"""
    r = json.loads(_run(code, tmp_path, _env()).stdout.splitlines()[-1])
    assert r["before"] == []  # no OpenMP runtime was in the process when the variable was set
    assert r["seen"], "no extension module or I/O package was imported"
    assert set(r["seen"].values()) == {SPIN}, r["seen"]
    assert {f"petakit5d._ext.{g}" for g in r["groups"]} <= set(r["seen"])


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="libgomp's futex waits: Linux")
def test_each_libgomp_reads_it(tmp_path):
    # OMP_DISPLAY_ENV=VERBOSE: every libgomp prints its settings when it is loaded
    if not _ext.available("core"):
        pytest.skip("the core extension module is not built")
    code = (
        "import json, petakit5d\n"
        "from petakit5d import _ext\n"
        "_ext.load('core')\n"
        "try:\n"
        "    import cppzarr, cpptiff\n"
        "except ImportError:\n"
        "    pass\n"
        "with open('/proc/self/maps') as f:\n"
        "    print(json.dumps(sorted({ln.split()[-1] for ln in f"
        " if 'gomp' in ln.rsplit('/', 1)[-1]})))\n"
    )
    r = _run(code, tmp_path, _env(OMP_DISPLAY_ENV="VERBOSE"))
    runtimes = json.loads(r.stdout.splitlines()[-1])
    lines = [ln for ln in r.stderr.splitlines() if "GOMP_SPINCOUNT = " in ln]
    if not runtimes:
        pytest.skip("no libgomp in the process (another OpenMP runtime)")
    assert len(lines) == len(runtimes), (runtimes, r.stderr[-2000:])
    assert all(ln.strip().endswith(f"GOMP_SPINCOUNT = '{SPIN}'") for ln in lines), lines


def test_a_subprocess_sees_it(tmp_path):
    code = (
        "import subprocess, sys, petakit5d\n"
        "c = subprocess.run([sys.executable, '-c', "
        "'import os; print(os.environ.get(\"GOMP_SPINCOUNT\"))'],"
        " capture_output=True, text=True, check=True)\n"
        "print(c.stdout.strip())\n"
    )
    assert _run(code, tmp_path, _env()).stdout.split() == [SPIN]
