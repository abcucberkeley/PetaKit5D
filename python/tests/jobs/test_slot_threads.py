"""proposed for tests/jobs/test_slots.py (J03): each slot's task process gets its share of the
job's CPUs, max(1, CPUs // slots) threads and (Linux) a block of CPUs of its own"""

import json
import os
from pathlib import Path

import pytest

from petakit5d.jobs import slots
from petakit5d.jobs.tasks import Call, Task, task_record, write_call

MODULE = """
import json, os


def seen(out):
    cpus = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
    with open(out, "w") as fh:
        json.dump({"omp": os.environ.get("OMP_NUM_THREADS"), "cpus": cpus,
                   "slot": os.environ.get("PARALLEL_JOBSLOT"),
                   "var": os.environ.get("PETAKIT5D_SLOT_CPUS")}, fh)
"""


def test_slot_shares():
    from petakit5d.core.threads import available_cpus

    n = available_cpus()
    shares = slots.slot_shares(2)
    assert sorted(shares) == [1, 2]
    assert all(t == max(1, n // 2) for t, _ in shares.values())
    if hasattr(os, "sched_getaffinity") and n >= 2:
        a, b = shares[1][1], shares[2][1]
        assert a and b and not set(a) & set(b)
    assert all(t == 1 for t, _ in slots.slot_shares(1000).values())
    assert all(c is None for _, c in slots.slot_shares(1000).values())


@pytest.mark.skipif(not hasattr(os, "sched_getaffinity"), reason="CPU affinity (Linux)")
def test_slot_task_processes_get_their_share(tmp_path, monkeypatch):
    from petakit5d.core.threads import available_cpus

    if available_cpus() < 2:
        pytest.skip("one CPU")
    (tmp_path / "slotseen.py").write_text(MODULE)
    outs = [str(tmp_path / f"o_{k}.json") for k in range(4)]
    recs = [task_record(Task("slotseen:seen", [o], path=[str(tmp_path)]), "python") for o in outs]
    folder = write_call(Call(["a"] * 4, outs, recs, {"task_batch_num": 4}), tmp_path / "c")
    src = str(Path(slots.__file__).resolve().parents[2])
    monkeypatch.setenv("PYTHONPATH", src)
    assert slots.run_slots(folder, range(4), 2, delay=0) == 0
    seen = [json.loads(Path(o).read_text()) for o in outs]
    share = available_cpus() // 2
    assert {x["omp"] for x in seen} == {str(share)}
    assert all(len(x["cpus"]) == share for x in seen)
    assert {tuple(x["cpus"]) for x in seen if x["slot"] == "1"}.isdisjoint(
        {tuple(x["cpus"]) for x in seen if x["slot"] == "2"}
    )
    assert all(x["var"] is None for x in seen), "the variable is the task process's own"
