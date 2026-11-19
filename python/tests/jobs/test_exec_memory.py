"""The port's task processes start with NumPy's huge-page advice off (``NUMPY_MADVISE_HUGEPAGE=0``,
``jobs.exec.MEMORY_VARIABLES``) unless the user set it: the job entry point, slot task processes
and local pool workers. (E03: on nodes with fragmented free memory the advice made a fresh
process stall for seconds in the kernel's compaction at the first touch of large arrays.)"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import petakit5d
from petakit5d.jobs.tasks import Call, Task, task_record, write_call

SRC_ROOT = str(Path(petakit5d.__file__).resolve().parents[1])

PROBE = """
import os, sys
SEEN = os.environ.get("NUMPY_MADVISE_HUGEPAGE")
NUMPY_BEFORE = "numpy" in sys.modules

def report(out):
    import json
    import numpy._core.multiarray as m
    state = m._set_madvise_hugepage(False) if hasattr(m, "_set_madvise_hugepage") else None
    with open(out, "w") as fh:
        json.dump({"seen": SEEN, "numpy_before": NUMPY_BEFORE, "numpy_madvise": state}, fh)
"""


def _run(tmp_path, env_value):
    (tmp_path / "memprobe.py").write_text(PROBE)
    out = tmp_path / "seen.json"
    rec = task_record(Task("memprobe:report", [str(out)], path=[str(tmp_path)]), "python")
    folder = write_call(Call(["a"], [str(out)], [rec], {}), tmp_path / "c")
    env = {k: v for k, v in os.environ.items() if k != "NUMPY_MADVISE_HUGEPAGE"}
    env["PYTHONPATH"] = os.pathsep.join([SRC_ROOT, env.get("PYTHONPATH", "")])
    if env_value is not None:
        env["NUMPY_MADVISE_HUGEPAGE"] = env_value
    r = subprocess.run(
        [sys.executable, "-m", "petakit5d.jobs.exec", folder, "1"],
        capture_output=True,
        text=True,
        env=env,
        timeout=300,
    )
    assert r.returncode == 0, r.stderr
    return json.loads(out.read_text())


def test_job_entry_point_turns_numpys_huge_page_advice_off(tmp_path):
    got = _run(tmp_path, None)
    assert got["seen"] == "0"
    assert got["numpy_before"] is False
    if got["numpy_madvise"] is not None and sys.platform.startswith("linux"):
        # NumPy read the variable at its import: the advice is off
        assert got["numpy_madvise"] is False


def test_job_entry_point_keeps_the_users_setting(tmp_path):
    got = _run(tmp_path, "1")
    assert got["seen"] == "1"


def test_slot_and_pool_worker_environments(monkeypatch):
    from petakit5d.jobs import local, slots

    monkeypatch.delenv("NUMPY_MADVISE_HUGEPAGE", raising=False)
    assert slots.slot_env(2, None)["NUMPY_MADVISE_HUGEPAGE"] == "0"
    slot = local._Slot(1, 2, None, None)
    assert local._worker_env(slot)["NUMPY_MADVISE_HUGEPAGE"] == "0"
    monkeypatch.setenv("NUMPY_MADVISE_HUGEPAGE", "1")
    assert "NUMPY_MADVISE_HUGEPAGE" not in slots.slot_env(2, None)
    assert local._worker_env(slot)["NUMPY_MADVISE_HUGEPAGE"] == "1"
