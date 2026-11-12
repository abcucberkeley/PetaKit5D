"""The parallel rsync wrapper, the memory estimate and the partial-volume grouping (T05).

Ported from tests/tools/ParallelRsyncWrapperTest.m (names in snake_case); the rsync tests use a
fake rsync on PATH and temporary folders only. Oracle cases compare the memory estimate, the
grouping and the shuffle's permutation with MATLAB.
"""

from __future__ import annotations

import json
import math
import os
import random
import re
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core import axes
from petakit5d.jobs import framework
from petakit5d.params.registry import REGISTRY
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, File, InPath, PerAxis, Tree, Work
from petakit5d.testing.compare import MCell, ignore
from petakit5d.testing.slurm_fake import no_slurm  # noqa: F401 (a fixture)
from petakit5d.tools import memory_estimate as me
from petakit5d.tools import partial_volumes as pv
from petakit5d.tools import rsync

needs_bash = pytest.mark.skipif(
    sys.platform.startswith("win"), reason="rsync and bash tasks are not available on Windows"
)

FAKE_RSYNC = '''\
#!{python}
"""fake rsync: copies as rsync -a does (a source with a trailing / copies its contents) and
logs each call (options, source, destination) as a JSON line to $FAKE_RSYNC_LOG"""
import json, os, shutil, sys

args = sys.argv[1:]
opts = [a for a in args if a.startswith("-")]
src, dst = [a for a in args if not a.startswith("-")]
with open(os.environ["FAKE_RSYNC_LOG"], "a") as f:
    f.write(json.dumps([opts, src, dst]) + "\\n")
if src.endswith("/"):
    shutil.copytree(src, dst, dirs_exist_ok=True)
elif os.path.isdir(src):
    shutil.copytree(src, os.path.join(dst, os.path.basename(src)), dirs_exist_ok=True)
else:
    os.makedirs(dst, exist_ok=True)
    shutil.copy2(src, os.path.join(dst, os.path.basename(src)) if dst.endswith("/") else dst)
'''


@pytest.fixture
def fake_rsync(tmp_path, monkeypatch):
    """A fake rsync first on PATH; returns a function that reads its call log."""
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    exe = bin_dir / "rsync"
    exe.write_text(FAKE_RSYNC.format(python=sys.executable))
    exe.chmod(0o755)
    log = tmp_path / "rsync_calls.ndjson"
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_RSYNC_LOG", str(log))
    # the framework's default job log folder is relative to the current folder
    work = tmp_path / "cwd" / "here"
    work.mkdir(parents=True)
    monkeypatch.chdir(work)

    def calls():
        if not log.exists():
            return []
        return [json.loads(x) for x in log.read_text().splitlines()]

    return calls


@pytest.fixture
def spy(monkeypatch):
    """Records run_tasks calls (MATLAB's FrameworkSpyFixture): nothing runs."""
    calls = []

    def fake(inputs, outputs, tasks, **kw):
        calls.append({"inputs": list(inputs), "outputs": list(outputs), "tasks": list(tasks), **kw})
        return [True] * len(tasks)

    monkeypatch.setattr(framework, "run_tasks", fake)
    return calls


def make_files(folder: Path, names) -> list[str]:
    folder.mkdir(parents=True, exist_ok=True)
    for n in names:
        p = folder / n
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(n)
    return [str(folder / n) for n in names]


def batch_args(command: str) -> list[str]:
    """The source and destination arguments of a batch command."""
    words = shlex.split(command)
    assert words[:2] == ["bash", "-c"] and words[2] == rsync.BATCH_SCRIPT
    assert words[3] == "rsync_batch_files"
    return words[4:]


# ---------------------------------------------------------------------------------------------
# parallel_rsync (ParallelRsyncWrapperTest)


@needs_bash
def test_batch_script_is_found_whichever_setup_is_first(tmp_path, fake_rsync, spy):
    # MATLAB took the code folder from the first setup.m on the path; the batch script here
    # is part of the command, so a batch runs from any current folder and install location
    src = tmp_path / "src"
    make_files(src, [f"f{i}.txt" for i in range(1, 5)])
    dest = tmp_path / "dest"
    rsync.parallel_rsync(str(src), str(dest), parse_cluster=False)
    assert spy, "The job framework was not called."
    other = tmp_path / "elsewhere"
    other.mkdir()
    for cmd in spy[0]["tasks"]:
        subprocess.run(["bash", "-c", cmd], cwd=other, check=True, stdout=subprocess.DEVNULL)
    batch_calls = [c for c in fake_rsync() if c[0] == ["-aP"]]
    assert len(batch_calls) == 4
    assert sorted(os.listdir(dest / "src")) == ["f1.txt", "f2.txt", "f3.txt", "f4.txt", "tmp"]
    final = [c for c in fake_rsync() if c[0] == ["-avP"]]
    assert final == [[["-avP"], f"{src}/", f"{dest}/src/"]], "The final rsync is missing."


def test_each_batch_is_checked_by_its_last_file(tmp_path, fake_rsync, spy):
    # 1252 files: batches of 3, the last with 1 file; each is checked by its last file
    src = tmp_path / "src"
    make_files(src, [f"f{i:04d}.txt" for i in range(1, 1253)])
    dest = tmp_path / "dest"
    rsync.parallel_rsync(str(src), str(dest), parse_cluster=False)
    c = spy[0]
    assert len(c["tasks"]) == math.ceil(1252 / 3) == len(c["inputs"]) == len(c["outputs"])
    for cmd, inp, out in zip(c["tasks"], c["inputs"], c["outputs"], strict=True):
        args = batch_args(cmd)
        assert len(args) in (2, 6)
        assert args[-2] == inp, "Batches are not checked by their last file."
        assert out == f"{args[-1]}/{os.path.basename(inp)}"
    assert len(batch_args(c["tasks"][-1])) == 2
    assert len(os.listdir(dest / "src")) == 1252 + 1  # and the tmp folder MATLAB leaves


def test_shuffles_with_rng1_and_leaves_the_global_random_numbers(tmp_path, fake_rsync, spy):
    # batches of one file in the order of rng(1)'s permutation; NumPy's and Python's global
    # random states are left as they were
    src = tmp_path / "src"
    names = make_files(src, [f"f{i:04d}.txt" for i in range(1, 41)])
    np.random.seed(7)  # noqa: NPY002 (the global state the wrapper must leave alone)
    random.seed(7)
    before = (np.random.get_state()[1].copy(), random.getstate())  # noqa: NPY002
    rsync.parallel_rsync(str(src), str(tmp_path / "dest"), parse_cluster=False)
    assert np.array_equal(np.random.get_state()[1], before[0])  # noqa: NPY002
    assert random.getstate() == before[1]
    files = [batch_args(t)[0] for t in spy[0]["tasks"]]
    # rng(1); randperm(40) in MATLAB (oracle case randperm_seed1_40)
    first = [3, 39, 15, 28, 27, 6, 31, 19, 5, 29]
    assert files[:10] == [names[i - 1] for i in first]
    assert sorted(files) == names


def test_randperm_seed1_is_matlabs():
    # MATLAB R2024b: randperm(RandStream('twister', 'Seed', 1), n)
    want = {
        1: [1],
        3: [3, 1, 2],
        5: [3, 5, 4, 1, 2],
        10: [3, 6, 5, 7, 4, 8, 9, 1, 10, 2],
    }
    for n, p in want.items():
        assert [i + 1 for i in rsync.randperm_seed1(n)] == p
    assert rsync.randperm_seed1(0) == []


@pytest.mark.parametrize(
    ("n", "size"),
    [
        (1, 1),
        (249, 1),
        (250, 1),
        (1249, 2),
        (1250, 3),
        (1252, 3),
        (2500, 5),
        (49999, 5),
        (75000, 8),
        (100000, 10),
        (1_000_000, 100),
        (3_000_000, 200),
        (10**8, 200),
    ],
)
def test_batch_size(n, size):
    # max(1, round(min(nF / 500, min(max(5, nF / 10000), 200)))), round half away from zero
    assert rsync.batch_size_for(n) == size


@needs_bash
def test_local_mode_copies_every_file_with_its_batches(tmp_path, fake_rsync):
    # the whole path with the real framework in local mode: the batches copy the files (names
    # with spaces and commas, subfolders), the final rsync runs once, the temporary folder of
    # the call is removed and the source is left as it was
    src = tmp_path / "data set"
    names = ["a.txt", "b c.txt", "d,e.txt", "sub/f.txt", "sub/deeper/g h.txt", "z.tif"]
    make_files(src, names)
    dest = tmp_path / "out"
    rsync.parallel_rsync(str(src), str(dest), parse_cluster=False, master_compute=True)
    calls = fake_rsync()
    batch = [c for c in calls if c[0] == ["-aP"]]
    assert sorted(c[1] for c in batch) == sorted(f"{src}/{n}" for n in names)
    for c in batch:
        rel = os.path.dirname(os.path.relpath(c[1], src))
        assert c[2] == (f"{dest}/data set/{rel}" if rel else f"{dest}/data set") + "/"
    assert [c for c in calls if c[0] == ["-avP"]] == [[["-avP"], f"{src}/", f"{dest}/data set/"]]
    for n in names:
        assert (dest / "data set" / n).read_text() == n
        assert (src / n).read_text() == n
    assert os.listdir(dest / "data set" / "tmp") == []


@needs_bash
def test_existing_files_are_not_copied_again(tmp_path, fake_rsync):
    src = tmp_path / "src"
    make_files(src, ["a.txt", "b.txt", "c.txt"])
    dest = tmp_path / "dest"
    make_files(dest / "src", ["a.txt", "c.txt"])
    rsync.parallel_rsync(str(src), str(dest), parse_cluster=False)
    batch = [c[1] for c in fake_rsync() if c[0] == ["-aP"]]
    assert batch == [f"{src}/b.txt"]


@needs_bash
def test_source_with_trailing_separator_copies_its_contents(tmp_path, fake_rsync, spy):
    src = tmp_path / "src"
    make_files(src, ["a.txt", "s/b.txt"])
    dest = tmp_path / "dest"
    rsync.parallel_rsync(f"{src}/", str(dest), parse_cluster=False)
    assert sorted(spy[0]["outputs"]) == [f"{dest}/a.txt", f"{dest}/s/b.txt"]
    folders = sorted(batch_args(t)[1] for t in spy[0]["tasks"])
    assert folders == [str(dest), f"{dest}/s"]
    assert [c for c in fake_rsync() if c[0] == ["-avP"]] == [[["-avP"], f"{src}/", f"{dest}/"]]


@needs_bash
def test_without_subfolders_only_the_top_files_are_batched(tmp_path, fake_rsync, spy):
    src = tmp_path / "src"
    make_files(src, ["a.txt", "B.txt", ".hidden", "s/b.txt"])
    rsync.parallel_rsync(str(src), str(tmp_path / "d"), include_subdir=False, parse_cluster=False)
    files = sorted(batch_args(t)[0] for t in spy[0]["tasks"])
    assert files == sorted(f"{src}/{n}" for n in ["a.txt", "B.txt", ".hidden"])


def test_framework_gets_the_wrappers_settings(tmp_path, fake_rsync, spy):
    src = tmp_path / "src"
    make_files(src, ["a.txt"])
    cfg = str(tmp_path / "cfg.json")
    rsync.parallel_rsync(
        str(src),
        str(tmp_path / "d"),
        cpus_per_task=2,
        num_stream=7,
        slurm_param="-p compute",
        nice_factor=5,
        parse_cluster=False,
        master_compute=True,
        config_file=cfg,
    )
    c = spy[0]
    assert c["language"] == "bash"
    assert c["slurm_param"] == "-p compute --nice=5"
    assert c["cpus_per_task"] == 2 and c["max_job_num"] == 7
    assert c["master_compute"] is True and c["parse_cluster"] is False
    assert c["config_file"] == cfg
    assert re.fullmatch(rf"{tmp_path}/d/src//tmp/[0-9a-f-]+/", c["tmp_dir"])
    assert not os.path.exists(c["tmp_dir"]), "The call's temporary folder was not removed."
    rsync.parallel_rsync(str(src), str(tmp_path / "e"), parse_cluster=False)
    assert spy[5]["slurm_param"] == " --nice=0" and spy[5]["max_job_num"] == math.inf


def test_source_folder_name_keeps_its_extension(tmp_path, fake_rsync, spy):
    src = tmp_path / "run.v2"
    make_files(src, ["a.txt"])
    dest = tmp_path / "d"
    rsync.parallel_rsync(str(src), f"{dest}//", parse_cluster=False)
    assert spy[0]["outputs"] == [f"{dest}/run.v2/a.txt"]


@needs_bash
def test_batch_runs_in_a_job_from_its_call_folder(tmp_path, fake_rsync, spy):
    # what a Slurm job runs: the job entry point with the call folder, from another folder
    from petakit5d.jobs.tasks import Call, write_call

    src = tmp_path / "src"
    make_files(src, ["a b.txt", "c,d.txt"])
    dest = tmp_path / "dest"
    rsync.parallel_rsync(str(src), str(dest), parse_cluster=False)
    c = spy[0]
    recs = [{"bash": t} for t in c["tasks"]]
    folder = write_call(
        Call(c["inputs"], c["outputs"], recs, {"language": "bash"}), tmp_path / "call"
    )
    root = str(Path(rsync.__file__).resolve().parents[2])
    for b in range(1, len(recs) + 1):
        r = subprocess.run(
            [sys.executable, "-m", "petakit5d.jobs.exec", folder, str(b)],
            cwd=tmp_path / "cwd",
            capture_output=True,
            text=True,
            env=dict(os.environ, PYTHONPATH=root),
            timeout=120,
        )
        assert r.returncode == 0, r.stderr
    assert all(os.path.exists(o) for o in c["outputs"])


def test_up_to_five_rounds_copy_what_is_missing(tmp_path, fake_rsync, spy):
    # the spy copies nothing: five rounds of the same batch, then the final rsync
    src = tmp_path / "src"
    make_files(src, ["a.txt"])
    rsync.parallel_rsync(str(src), str(tmp_path / "d"), parse_cluster=False)
    assert len(spy) == 5
    assert len([c for c in fake_rsync() if c[0] == ["-avP"]]) == 1


def test_large_copies_pass_their_settings_to_the_cluster_check(
    tmp_path, fake_rsync, spy, monkeypatch
):
    # copies of 100000 files or more are checked by jobs (a stub here that reports every file
    # as copied): the wrapper's settings go to the check, and nothing is copied again
    src = tmp_path / "src"
    make_files(src, ["a.txt", "b.txt", "c.txt"])
    monkeypatch.setattr(rsync, "CLUSTER_CHECK_FILES", 3)
    seen = []

    def check(paths, tmp_dir, **kw):
        seen.append((list(paths), tmp_dir, kw))
        return [True] * len(paths)

    monkeypatch.setattr(rsync, "batch_file_exist_cluster", check)
    rsync.parallel_rsync(str(src), str(tmp_path / "d"), parse_cluster=False, config_file="c.json")
    assert len(seen) == 1 and len(seen[0][0]) == 3
    assert re.search(r"/tmp/[0-9a-f-]+//\d{14}/$", seen[0][1])
    assert seen[0][2] == {"parse_cluster": False, "mcc_mode": False, "config_file": "c.json"}
    assert spy == [], "The files were copied again."


@pytest.mark.usefixtures("no_slurm")
def test_cluster_file_check_finds_the_files(tmp_path, monkeypatch):
    # batch_file_exist_cluster: MATLAB's used undefined variables and always failed; without
    # Slurm the framework runs the checks here (default parse_cluster), as with it false
    # (3 per job here; a job of one path writes no result file, as in MATLAB: see the report)
    monkeypatch.setattr(rsync, "CLUSTER_CHECK_BATCH", 3)
    calls = []
    real = framework.run_tasks

    def recording(*a, **kw):
        calls.append(kw)
        return real(*a, **kw)

    monkeypatch.setattr(framework, "run_tasks", recording)
    (tmp_path / "cwd").mkdir()
    monkeypatch.chdir(tmp_path / "cwd")
    names = ["a b.txt", "c,d.txt", "it's.txt", "plain.txt", "missing.txt"]
    make_files(tmp_path, names[:4])
    paths = [str(tmp_path / n) for n in names]
    want = [True, True, True, True, False]
    assert rsync.batch_file_exist_cluster(paths, str(tmp_path / "check1")) == want
    out = rsync.batch_file_exist_cluster(paths, str(tmp_path / "check2"), parse_cluster=False)
    assert out == want
    flags = sorted(os.listdir(tmp_path / "check1" / "file_exist_flags"))
    assert flags == ["file_exist_1_3.mat", "file_exist_4_5.mat", "parameters.mat"]
    assert [c["parse_cluster"] for c in calls] == [True, False]
    assert calls[0]["mem_allocate"] == 20 and calls[0]["master_compute"] is True


@needs_bash
def test_names_with_spaces_commas_and_quotes_are_copied(tmp_path, fake_rsync, spy):
    # MATLAB escaped spaces as '\\ ' (also in the checked files), split its lists at commas
    # and left the final rsync's paths unquoted; 750 files: batches of 2
    import shutil

    from petakit5d.jobs.tasks import Call, write_call

    src = tmp_path / "my src"
    special = ["a b.txt", "c,d.txt", "it's.txt", "sub dir/e,f g.txt"]
    names = special + [f"f{i:04d}.txt" for i in range(1, 747)]
    make_files(src, names)
    dest = tmp_path / "out dir"
    copy = dest / "my src"
    rsync.parallel_rsync(str(src), str(dest), parse_cluster=False)
    c = spy[0]
    assert all(os.path.isfile(p) for p in c["inputs"]), "A batch is checked by a missing file."
    assert sorted(os.path.relpath(p, copy) for p in c["outputs"]) == sorted(
        os.path.relpath(p, src) for p in c["inputs"]
    )
    for n in names:
        assert (copy / n).read_text() == n, "The final rsync"
    final = [x for x in fake_rsync() if x[0] == ["-avP"]]
    assert final == [[["-avP"], f"{src}/", f"{copy}/"]]

    batches = [t for t in c["tasks"] if any(f"{src}/{n}" in batch_args(t) for n in special)]
    files = [os.path.relpath(a, src) for t in batches for a in batch_args(t)[0::2]]
    assert all(len(batch_args(t)) == 4 for t in batches)
    assert set(special) <= set(files), "The batches leave out a file."
    # run as this process runs them, and as a job runs them, in an empty folder
    cwd = tmp_path / "cwd2"
    cwd.mkdir()
    root = str(Path(rsync.__file__).resolve().parents[2])
    recs = [{"bash": t} for t in batches]
    folder = write_call(
        Call(["x"] * len(recs), ["y"] * len(recs), recs, {"language": "bash"}), tmp_path / "call"
    )
    for mode in ("local", "job"):
        shutil.rmtree(copy)
        for b, t in enumerate(batches, start=1):
            if mode == "local":
                cmd = ["bash", "-c", t]
            else:
                cmd = [sys.executable, "-m", "petakit5d.jobs.exec", folder, str(b)]
            r = subprocess.run(
                cmd,
                cwd=cwd,
                capture_output=True,
                text=True,
                timeout=120,
                env=dict(os.environ, PYTHONPATH=root),
            )
            assert r.returncode == 0, r.stderr
            assert not re.search("error|cannot", r.stdout + r.stderr, re.I), r.stdout + r.stderr
        for f in files:
            assert (copy / f).read_text() == f, f"The batches ({mode})"
        assert os.listdir(cwd) == [], "The batches made files elsewhere."


def test_shell_quote_is_matlabs(monkeypatch):
    monkeypatch.setenv("HOME", "/users/someone")
    assert rsync.shell_quote("/a b/it's") == "'/a b/it'\\''s'"
    assert rsync.shell_quote("~/x y") == "'/users/someone/x y'"
    assert rsync.shell_quote("~") == "'/users/someone'"
    assert rsync.shell_quote("~other/x") == "'~other/x'"
    assert shlex.split(rsync.shell_quote("a'b\"c")) == ["a'b\"c"]


@needs_bash
def test_source_under_home_is_listed_and_copied(tmp_path, fake_rsync, spy, monkeypatch):
    # MATLAB's dir and shell_quote expand a leading ~; the checked paths keep it
    home = tmp_path / "home"
    monkeypatch.setenv("HOME", str(home))
    make_files(home / "src", ["a.txt"])
    rsync.parallel_rsync("~/src", "~/out", parse_cluster=False)
    c = spy[0]
    assert c["inputs"] == ["~/src/a.txt"] and c["outputs"] == ["~/out/src/a.txt"]
    assert batch_args(c["tasks"][0]) == [f"{home}/src/a.txt", f"{home}/out/src"]
    assert (home / "out" / "src" / "a.txt").read_text() == "a.txt"
    assert not os.path.exists("~")


def test_rsync_is_registered():
    spec = REGISTRY.find("XR_parallel_rsync_wrapper")
    assert spec.python_name == "petakit5d.tools.rsync.parallel_rsync"
    assert spec.param("SlurmParam").default == ""


# ---------------------------------------------------------------------------------------------
# estimate_computing_memory

GB = 1024.0**3


def test_estimate_deskew_from_the_given_size(tmp_path):
    # imSize (nz, ny, nx) in Python; deskew grows with (nz / 300)^2 above 300 slices
    missing = str(tmp_path / "missing.tif")
    est, gpu, raw, sz = me.estimate_computing_memory(missing, ("deskew",), im_size=(600, 10, 20))
    assert raw == 600 * 10 * 20 * 4 / GB
    assert est.tolist() == [raw * 15 * 4.0]
    assert gpu is None and sz == (600, 10, 20)
    est = me.estimate_computing_memory(missing, ("deskew",), im_size=(299, 10, 20))[0]
    assert est.tolist() == [299 * 10 * 20 * 4 / GB * 15]


# ComputingMemoryEstimateTest (MATLAB sizes [ny nx nz] = [100 200 600] are (600, 100, 200) here)
SZ = (600, 100, 200)
RAW = 100 * 200 * 600 * 4 / GB


def test_each_step_of_several_gets_its_estimate():
    # a step was estimated only if every step matched it (before 4b79f9b)
    est, gpu, raw, sz = me.estimate_computing_memory("", im_size=SZ)
    assert est.tolist() == [RAW * 15 * 4, RAW * 5, RAW * 15]
    assert gpu == RAW * 1.5 and raw == RAW and sz == SZ
    est = me.estimate_computing_memory("", ("rotate", "deskew"), im_size=SZ)[0]
    assert est.tolist() == [RAW * 5, RAW * 15 * 4]


def test_data_size_with_several_steps():
    est = me.estimate_computing_memory(
        "", ("rotate", "deconvolution"), data_size=2**30, cuda_decon=False
    )[0]
    assert est.tolist() == [2 * 5, 2 * 15]
    # with deskew among them, the image size counts
    est = me.estimate_computing_memory("", ("deskew", "rotate"), im_size=SZ, data_size=2**30)[0]
    assert est.tolist() == [RAW * 15 * 4, RAW * 5]


def test_single_steps_as_the_callers_pass_them(tmp_path):
    from petakit5d.testing.synth import write_tiff_matlab as write_tiff

    sz = (301, 64, 48)
    raw = 64 * 48 * 301 * 4 / GB
    assert me.estimate_computing_memory("", ("deskew",), im_size=sz)[0].tolist() == [
        raw * 15 * (301 / 300) * (301 / 300)
    ]
    assert me.estimate_computing_memory("", ("rotate",), im_size=sz)[0].tolist() == [raw * 5]
    est, gpu, _, _ = me.estimate_computing_memory("", ("deconvolution",), im_size=sz)
    assert est.tolist() == [raw * 15] and gpu == raw * 1.5
    gpu = me.estimate_computing_memory("", ("deconvolution",), im_size=sz, cuda_decon=False)[1]
    assert math.isnan(gpu)
    # a tiff file, as XR_cppDeconFrame3D passes it
    fn = tmp_path / "im.tif"
    write_tiff(fn, np.zeros((16, 12, 5), dtype=np.uint16, order="F"))
    est, _, raw, sz = me.estimate_computing_memory(str(fn), ("deconvolution",), cuda_decon=False)
    assert sz == (5, 16, 12) and raw == 16 * 12 * 5 * 4 / GB and est.tolist() == [raw * 15]


def test_estimate_with_text_steps_has_one_value_per_character():
    est = me.estimate_computing_memory("", "deskew", im_size=(5, 6, 7))[0]
    assert est.shape == (6,) and est[0] > 0 and not est[1:].any()


def test_estimate_deconvolution_on_the_gpu_and_cpu():
    est, gpu, raw, sz = me.estimate_computing_memory(
        "", ("deconvolution",), data_size=1e9, mem_factors=(1, 2, 3), gpu_mem_factor=2.5
    )
    assert raw == 1e9 * 2 / GB and est.tolist() == [raw * 3] and gpu == raw * 2.5 and sz is None
    gpu = me.estimate_computing_memory("", ("deconvolution",), data_size=1e9, cuda_decon=False)[1]
    assert math.isnan(gpu) and np.float64(gpu).view(np.uint64) == 0xFFF8000000000000


def test_estimate_errors_as_matlab():
    with pytest.raises(me.MatlabError, match="does not exist") as e:
        me.estimate_computing_memory("", ("deskew",))
    assert e.value.identifier == ""


def test_estimate_reads_the_size_of_the_file(tmp_path):
    from petakit5d.testing.synth import write_tiff_matlab as write_tiff

    f = tmp_path / "a.tif"
    write_tiff(f, np.zeros((4, 5, 3), dtype=np.uint16, order="F"))
    est, gpu, raw, sz = me.estimate_computing_memory(str(f), ("rotate",), im_size=(9, 9, 9))
    assert sz == (3, 4, 5) and raw == 60 * 4 / GB and est.tolist() == [raw * 5]


def test_matlab_layer_conversions():
    # imSize by the registry's size_yxz kind; the returned size by matlab_out
    spec = REGISTRY.find("XR_estimateComputingMemory")
    assert spec.param("imSize").conv.name == "size_yxz"
    assert not hasattr(me, "matlab_in") and spec.matlab_in is None
    assert spec.bind_matlab(["", ("rotate",), "imSize", [10, 20, 30]])["im_size"] == (30, 10, 20)
    assert spec.bind_matlab(["", ("rotate",), "imSize", []]).get("im_size") is None
    assert me.matlab_out((1, 2, 3, (3, 1, 2)))[3] == (1, 2, 3)
    assert me.matlab_out((1, 2, 3, (3, 1, 2, 4)))[3] == (1, 2, 3, 4)
    assert me.matlab_out((1, 2, 3, None))[3] is None


def test_matlab_name_call_gives_matlabs_size(tmp_path):
    from petakit5d import matlab
    from petakit5d.testing.synth import write_tiff_matlab as write_tiff

    est, _, raw, sz = matlab.XR_estimateComputingMemory("", ("deskew",), "imSize", [10, 20, 600])
    assert sz == (10, 20, 600) and raw == 10 * 20 * 600 * 4 / GB
    assert est.tolist() == [raw * 15 * 4.0]
    f = tmp_path / "a.tif"
    write_tiff(f, np.zeros((4, 5, 3), dtype=np.uint16, order="F"))
    assert matlab.XR_estimateComputingMemory(str(f), ("rotate",))[3] == (4, 5, 3)


@pytest.mark.parametrize("layout", ["native", "matlab"])
def test_estimate_reads_z_of_a_non_cubic_zarr(tmp_path, layout):
    # the slice count is the first of get_image_size's (nz, ny, nx), for both zarr layouts
    from petakit5d.testing.synth import write_zarr

    f = tmp_path / "a.zarr"
    write_zarr(f, np.zeros((600, 7, 5), dtype=np.uint16), block_size=(600, 7, 5), layout=layout)
    est, _, raw, sz = me.estimate_computing_memory(str(f), ("deskew",))
    assert sz == (600, 7, 5) and raw == 600 * 7 * 5 * 4 / GB
    assert est.tolist() == [raw * 15 * 4.0]


def test_estimate_of_a_2d_zarr_fails_as_matlab(tmp_path):
    from petakit5d.testing.synth import write_zarr

    f = tmp_path / "a.zarr"
    write_zarr(f, np.zeros((7, 5), dtype=np.uint16), block_size=(7, 5))
    with pytest.raises(me.MatlabError) as e:
        me.estimate_computing_memory(str(f), ("rotate",))
    assert e.value.identifier == "MATLAB:badsubscript"


@pytest.mark.parametrize("layout", ["native", "matlab"])
def test_partial_volume_sizes_of_non_cubic_zarr(tmp_path, layout):
    # sizes are products of get_image_size's values: the same for both zarr layouts
    from petakit5d.testing.synth import write_zarr

    for n, shape in {"v.zarr": (3, 5, 7), "v_part0001.zarr": (2, 5, 7)}.items():
        write_zarr(tmp_path / n, np.zeros(shape, dtype=np.uint16), block_size=shape, layout=layout)
    contain, names, _, sizes = pv.group_partial_volume_files(str(tmp_path), ext=".zarr")
    assert contain.tolist() == [False, True] and names == [["v.zarr", "v_part0001.zarr"]]
    assert sizes[0].tolist() == [3 * 5 * 7 * 2, 2 * 5 * 7 * 2]


# ---------------------------------------------------------------------------------------------
# group_partial_volume_files


def test_dir_datenum():
    # MATLAB's dir on a node in US Pacific time: 2019-02-21 17:11:42 is 737477.71645833331
    import datetime
    import time

    lt = time.localtime(1550797902)
    days = datetime.date(lt.tm_year, lt.tm_mon, lt.tm_mday).toordinal() + 366
    assert (
        pv.dir_datenum(1550797902.9)
        == days + (lt.tm_hour * 3600 + lt.tm_min * 60 + lt.tm_sec) / 86400
    )
    if (lt.tm_year, lt.tm_mon, lt.tm_mday, lt.tm_hour) == (2019, 2, 21, 17):
        assert pv.dir_datenum(1550797902) == 737477.71645833331


def test_grouping_keeps_matlabs_order_and_quirks(tmp_path):
    names = [
        "a.tif",
        "a_part0002.tif",
        "a_part0001.tif",
        "ab.tif",
        "ab_part0001.tif",
        "b.tif",
        "c_part0001.tif",
    ]
    files = [str(tmp_path / n) for n in names]
    for i, f in enumerate(files):
        Path(f).write_bytes(b"x" * (i + 1))
    contain, fnames, dn, ds = pv.group_partial_volume_files("", files)
    assert contain.tolist() == [False, True, True, False, True, False, True]
    # 'b' takes ab's part (the pattern is not anchored); c's part has no main file
    assert fnames == [
        ["a.tif", "a_part0001.tif", "a_part0002.tif"],
        ["ab.tif", "ab_part0001.tif"],
        ["b.tif", "ab_part0001.tif"],
    ]
    # names sorted, sizes in list order
    assert [x.tolist() for x in ds] == [[1, 2, 3], [4, 5], [6, 5]]
    assert len(dn) == 3 and dn[0].shape == (3,)


def test_grouping_without_parts_and_without_files(tmp_path):
    files = [_touch(tmp_path / "x.tif"), _touch(tmp_path / "y.tif")]
    contain, fnames, dn, ds = pv.group_partial_volume_files("", files)
    assert contain.tolist() == [False, False] and fnames == [["x.tif"], ["y.tif"]]
    with pytest.warns(pv.PartialVolumeWarning, match="empty"):
        out = pv.group_partial_volume_files(str(tmp_path / "empty_folder_that_is_missing"))
    assert out[0] is False and out[1] == [] and out[2].size == 0 and out[3].size == 0


def test_grouping_of_a_list_ignores_case_and_takes_the_first_time_point(tmp_path):
    names = ["s_Iter_0000_a.tif", "s_Iter_0000_a_PART0001.tif", "s_Iter_0001_a.tif"]
    files = [_touch(tmp_path / n) for n in names]
    contain, fnames, _, _ = pv.group_partial_volume_files("", files)
    assert contain.tolist() == [False, True, False]
    assert fnames == [["s_Iter_0000_a.tif", "s_Iter_0000_a_PART0001.tif"], ["s_Iter_0001_a.tif"]]
    contain, fnames, _, _ = pv.group_partial_volume_files("", files, only_first_tp="yes")
    assert fnames == [["s_Iter_0000_a.tif", "s_Iter_0000_a_PART0001.tif"]]


def _touch(p: Path) -> str:
    p.write_bytes(b"")
    return str(p)


def test_grouping_is_registered():
    spec = REGISTRY.find("groupPartialVolumeFiles")
    assert spec.python_name == "petakit5d.tools.partial_volumes.group_partial_volume_files"
    spec = REGISTRY.find("XR_estimateComputingMemory")
    assert spec.python_name == "petakit5d.tools.memory_estimate.estimate_computing_memory"


# ---------------------------------------------------------------------------------------------
# oracle cases


def randperm1(n):
    return np.array([i + 1 for i in rsync.randperm_seed1(int(n))], dtype=np.float64)


def _matlab_size(sz):
    # imSize as MATLAB returns it: [] for none, a double row otherwise
    return np.zeros((0, 0)) if sz is None else np.array(sz, dtype=np.float64)


def est_matlab(*args):
    """XR_estimateComputingMemory through the MATLAB-name layer (MATLAB's arguments; imSize
    converted by its kind, the returned size by matlab_out): for the cases about the binding."""
    mem, gpu, raw, sz = REGISTRY.find("XR_estimateComputingMemory").call_matlab(*args)
    return [mem, gpu, raw, _matlab_size(sz)]


_EST_KW = {
    "imSize": "im_size",
    "dataSize": "data_size",
    "memFactors": "mem_factors",
    "cudaDecon": "cuda_decon",
    "GPUMemFactor": "gpu_mem_factor",
}


def est_case(cid, args, nargout=1, **kw):
    """A case of XR_estimateComputingMemory whose Python side is the snake_case function: the
    leading arguments positional, the name-value pairs by their Python names, imSize (nz, ny,
    nx) through its role, the returned size back to MATLAB's [ny nx nz] by its out role."""
    k = next((i for i, a in enumerate(args) if isinstance(a, str) and a in _EST_KW), len(args))
    py_kwargs = {}
    for name, v in zip(args[k::2], args[k + 1 :: 2], strict=True):
        py_kwargs[_EST_KW[name]] = PerAxis(v, "yxz") if name == "imSize" else v
    roles = (None, None, None, lambda sz: _matlab_size(None if sz is None else axes.to_yxz(sz)))
    return Case(
        cid,
        matlab="XR_estimateComputingMemory",
        args=args,
        nargout=nargout,
        python=me.estimate_computing_memory,
        py_args=args[:k],
        py_kwargs=py_kwargs,
        py_out=first(1) if nargout == 1 else None,
        out_roles=roles[:nargout] if nargout == 4 else None,
        **kw,
    )


def first(n):
    def out(r):
        return list(r)[:n]

    return out


def group_matlab(*args):
    from petakit5d import matlab

    return list(matlab.groupPartialVolumeFiles(*args))


def tiff_writer(shape, dtype=np.uint16):
    def make(path):
        from petakit5d.testing.synth import write_tiff_matlab as write_tiff

        write_tiff(path, np.zeros(shape, dtype=dtype, order="F"))

    return make


def zarray(shape, dtype="<u2"):
    return json.dumps(
        {
            "chunks": list(shape),
            "compressor": None,
            "dimension_separator": ".",
            "dtype": dtype,
            "fill_value": 0,
            "filters": None,
            "order": "F",
            "shape": list(shape),
            "zarr_format": 2,
        }
    )


def make_zarr(folder: Path, shape, dtype="<u2"):
    folder.mkdir()
    (folder / ".zarray").write_text(zarray(shape, dtype))


def make_parts_tree(d: Path):
    from petakit5d.testing.synth import write_tiff_matlab as write_tiff

    sizes = {
        "a.tif": (10, 12, 3),
        "a_part0001.tif": (10, 12, 2),
        "a_part0002.tif": (10, 12, 1),
        "ab.tif": (8, 8, 2),
        "ab_part0001.tif": (8, 8, 1),
        "b.tif": (6, 7, 4),
        "c_part0001.tif": (5, 5, 1),
        "z.tif": (4, 4, 2),
    }
    for n, s in sizes.items():
        write_tiff(d / n, np.zeros(s, dtype=np.uint16, order="F"))
    (d / "notes.txt").write_text("not an image")


def make_plain_tree(d: Path):
    from petakit5d.testing.synth import write_tiff_matlab as write_tiff

    for n, s in {"x.tif": (6, 6, 2), "y.tif": (6, 6, 3)}.items():
        write_tiff(d / n, np.ones(s, dtype=np.uint8, order="F"))


def make_channel_tree(d: Path):
    from petakit5d.testing.synth import write_tiff_matlab as write_tiff

    for n in [
        "s_CamA_ch0_Iter_0000_t.tif",
        "s_CamA_ch0_Iter_0001_t.tif",
        "s_CamB_ch1_Iter_0000_t.tif",
        "s_CamB_ch1_Iter_0000_t_part0001.tif",
        "s_CamC_ch2_Iter_0000_t.tif",
    ]:
        write_tiff(d / n, np.zeros((4, 5, 2), dtype=np.float32, order="F"))


def make_zarr_tree(d: Path):
    make_zarr(d / "v.zarr", (10, 11, 12))
    make_zarr(d / "v_part0001.zarr", (10, 11, 3))
    make_zarr(d / "w.zarr", (2, 3, 4))


def make_text_tree(d: Path):
    (d / "a.txt").write_text("a")


SPECIAL_NAMES = [
    "a+b.tif",
    "a+b_part0001.tif",
    "c(1).tif",
    "c(1)_part0002.tif",
    "c(1)_part0001.tif",
    "x_Cam(A_ch0.tif",
    "x_Cam(A_ch0_part0001.tif",
    "zz.tif",
]
SPECIAL_LIST = [Work(f"p/{n}") for n in SPECIAL_NAMES]


def make_special_tree(d: Path):
    from petakit5d.testing.synth import write_tiff_matlab as write_tiff

    for i, n in enumerate(SPECIAL_NAMES):
        write_tiff(d / n, np.zeros((4, 5, i + 1), dtype=np.uint16, order="F"))


def test_grouping_matches_main_names_as_text(tmp_path):
    files = [_touch(tmp_path / n) for n in SPECIAL_NAMES]
    contain, fnames, _, _ = pv.group_partial_volume_files("", files)
    assert contain.tolist() == [False, True, False, True, True, False, True, False]
    assert fnames == [
        ["a+b.tif", "a+b_part0001.tif"],
        ["c(1).tif", "c(1)_part0001.tif", "c(1)_part0002.tif"],
        ["x_Cam(A_ch0.tif", "x_Cam(A_ch0_part0001.tif"],
        ["zz.tif"],
    ]


def test_channel_pattern_that_is_not_a_regexp_is_text(tmp_path):
    make_special_tree(tmp_path)
    _, fnames, _, _ = pv.group_partial_volume_files(str(tmp_path), channel_patterns=("Cam(A",))
    assert fnames == [["x_Cam(A_ch0.tif", "x_Cam(A_ch0_part0001.tif"]]
    with pytest.warns(pv.PartialVolumeWarning):
        pv.group_partial_volume_files(str(tmp_path), channel_patterns=("Cam(B",))


def make_list_tree(d: Path):
    for i, n in enumerate(["m.tif", "m_part0002.tif", "m_part0001.tif", "n.tif"]):
        (d / n).write_bytes(bytes(range(i + 1)) * 3)


# a column cell, as MATLAB's callers pass it (a row cell fails in MATLAB's concatenation)
LIST = [Work(f"p/{n}") for n in ["m.tif", "m_part0002.tif", "m_part0001.tif", "n.tif"]]
LIST_COL = MCell((len(LIST), 1), LIST)

NO_DATES = {"out3": ignore("dates of the work folder's copies differ between runs")}

TIFF_IN = {"t.tif": File(tiff_writer((20, 30, 7)), "t.tif")}

EST_CASES = [
    est_case("est_deskew_tiff", [InPath("t.tif"), ("deskew",)], inputs=TIFF_IN),
    Case(
        # inputParser takes text without a validator as a parameter name, not as the
        # optional steps (the MATLAB-name binding does the same since F13)
        "est_deskew_char_steps",
        matlab="XR_estimateComputingMemory",
        args=[InPath("t.tif"), "deskew"],
        python=est_matlab,
        py_out=first(1),
        inputs=TIFF_IN,
        expect_error="MATLAB:InputParser:ArgumentFailedDefaultValidation",
    ),
    Case(
        # the optional steps given by name, as XR_microscopeAutomaticProcessing calls it
        "est_steps_by_name",
        matlab="XR_estimateComputingMemory",
        args=[InPath("t.tif"), "steps", ("deskew",)],
        python=est_matlab,
        py_out=first(1),
        inputs=TIFF_IN,
    ),
    est_case("est_default_steps", [InPath("t.tif")], inputs=TIFF_IN),
    est_case("est_two_steps", [InPath("t.tif"), ("deskew", "deconvolution")], inputs=TIFF_IN),
    *[
        est_case(f"est_deskew_size_z{z}", ["missing.tif", ("Deskew",), "imSize", [512, 1800, z]])
        for z in (299, 300, 301, 450, 777, 1999, 12345)
    ],
    # one slice of a tiff: the bytes of a file that does not exist
    est_case(
        "est_deskew_size_z1_missing_tiff",
        ["missing.tif", ("deskew",), "imSize", [512, 1800, 1]],
        expect_error="MATLAB:minrhs",
        py_error=me.MatlabError,
    ),
    est_case("est_deskew_size_z1", ["missing.zarr", ("deskew",), "imSize", [512, 1800, 1]]),
    est_case(
        "est_rotate_mem_factors",
        ["", ("rotate",), "imSize", [100, 200, 300], "memFactors", [3, 2.7, 7.5]],
    ),
    est_case(
        "est_decon_cuda_tiff",
        [InPath("t.tif"), ("deconvolution",), "cudaDecon", True, "GPUMemFactor", 2.25],
        nargout=4,
        inputs=TIFF_IN,
    ),
    est_case(
        "est_decon_cpu_tiff",
        [InPath("t.tif"), ("deconvolution",), "cudaDecon", False],
        nargout=4,
        inputs=TIFF_IN,
    ),
    est_case(
        "est_decon_single_slice_tiff_bytes",
        [InPath("t.tif"), ("deconvolution",), "cudaDecon", False],
        nargout=4,
        inputs={"t.tif": File(tiff_writer((33, 17)), "t.tif")},
    ),
    est_case(
        "est_decon_zarr",
        [InPath("z.zarr"), ("deconvolution",)],
        nargout=4,
        inputs={
            "z.zarr": Tree(lambda d: (d / ".zarray").write_text(zarray((10, 11, 1))), "z.zarr")
        },
    ),
    est_case(
        "est_decon_data_size",
        ["", ("deconvolution",), "dataSize", 123456789, "cudaDecon", True],
        nargout=4,
    ),
    est_case(
        "est_missing_file",
        ["missing.tif", ("deskew",)],
        expect_error=True,
        py_error=me.MatlabError,
    ),
    est_case("est_data_size_several_steps", ["", ("rotate", "deconvolution"), "dataSize", 1000]),
    # ComputingMemoryEstimateTest
    Case(
        # the optional steps left out before a name-value pair ('imSize' is a parameter name)
        "est_default_steps_given_size",
        matlab="XR_estimateComputingMemory",
        args=["", "imSize", [100, 200, 600]],
        nargout=4,
        python=est_matlab,
    ),
    est_case("est_rotate_deskew_given_size", ["", ("rotate", "deskew"), "imSize", [100, 200, 600]]),
    est_case(
        "est_data_size_rotate_decon_cpu",
        ["", ("rotate", "deconvolution"), "dataSize", 2**30, "cudaDecon", False],
        nargout=4,
    ),
    est_case(
        "est_data_size_with_deskew",
        ["", ("deskew", "rotate"), "imSize", [100, 200, 600], "dataSize", 2**30],
    ),
]

GROUP_CASES = [
    Case(
        "group_parts_tiff",
        matlab="groupPartialVolumeFiles",
        args=[InPath("p")],
        nargout=4,
        python=group_matlab,
        inputs={"p": Tree(make_parts_tree, "p")},
        compare=NO_DATES,
    ),
    Case(
        "group_plain_tiff",
        matlab="groupPartialVolumeFiles",
        args=[InPath("p")],
        nargout=4,
        python=group_matlab,
        inputs={"p": Tree(make_plain_tree, "p")},
        compare=NO_DATES,
    ),
    Case(
        "group_channel_patterns",
        matlab="groupPartialVolumeFiles",
        args=[InPath("p"), "channelPatterns", ("CamA_ch0", "Cam[B]_ch1")],
        nargout=4,
        python=group_matlab,
        inputs={"p": Tree(make_channel_tree, "p")},
        compare=NO_DATES,
    ),
    Case(
        "group_only_first_tp",
        matlab="groupPartialVolumeFiles",
        args=[InPath("p"), "onlyFirstTP", "yes"],
        nargout=4,
        python=group_matlab,
        inputs={"p": Tree(make_channel_tree, "p")},
        compare=NO_DATES,
    ),
    Case(
        "group_zarr",
        matlab="groupPartialVolumeFiles",
        args=[InPath("p"), "ext", ".zarr"],
        nargout=4,
        python=group_matlab,
        inputs={"p": Tree(make_zarr_tree, "p")},
        compare=NO_DATES,
    ),
    Case(
        "group_no_images",
        matlab="groupPartialVolumeFiles",
        args=[InPath("p")],
        nargout=4,
        python=group_matlab,
        inputs={"p": Tree(make_text_tree, "p")},
        compare=NO_DATES,
    ),
    Case(
        "group_file_list",
        matlab="groupPartialVolumeFiles",
        args=["", LIST_COL],
        py_args=["", LIST],
        nargout=4,
        python=group_matlab,
        inputs={"p": Tree(make_list_tree, "p")},
        compare=NO_DATES,
    ),
    Case(
        "group_file_list_first_tp",
        matlab="groupPartialVolumeFiles",
        args=["", LIST_COL, "onlyFirstTP", ""],
        py_args=["", LIST, "onlyFirstTP", ""],
        nargout=4,
        python=group_matlab,
        inputs={"p": Tree(make_list_tree, "p")},
        compare=NO_DATES,
    ),
    Case(
        # a row list failed in MATLAB's concatenation before b24c975
        "group_file_list_row",
        matlab="groupPartialVolumeFiles",
        args=["", LIST],
        nargout=4,
        python=group_matlab,
        inputs={"p": Tree(make_list_tree, "p")},
        compare=NO_DATES,
    ),
    Case(
        # main names with regular expression characters ('+', '(1)') are matched as text
        "group_names_with_regexp_characters",
        matlab="groupPartialVolumeFiles",
        args=[InPath("p")],
        nargout=4,
        python=group_matlab,
        inputs={"p": Tree(make_special_tree, "p")},
        compare=NO_DATES,
    ),
    Case(
        "group_list_names_with_regexp_characters",
        matlab="groupPartialVolumeFiles",
        args=["", SPECIAL_LIST],
        nargout=4,
        python=group_matlab,
        inputs={"p": Tree(make_special_tree, "p")},
        compare=NO_DATES,
    ),
    Case(
        # a channel pattern that is not a valid regular expression is matched as text
        "group_channel_pattern_not_a_regexp",
        matlab="groupPartialVolumeFiles",
        args=[InPath("p"), "channelPatterns", ("Cam(A", "a\\+b")],
        nargout=4,
        python=group_matlab,
        inputs={"p": Tree(make_special_tree, "p")},
        compare=NO_DATES,
    ),
]

RANDPERM_CASES = [
    Case(
        f"randperm_seed1_{n}",
        code="out = {randperm(RandStream('twister', 'Seed', 1), in{1})};",
        args=[n],
        python=randperm1,
    )
    for n in (1, 7, 40, 1000, 123457)
]


@oracle.parametrize(EST_CASES + GROUP_CASES + RANDPERM_CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_oracle_inputs_are_deterministic(tmp_path):
    for make in (make_parts_tree, make_channel_tree, make_zarr_tree, make_list_tree):
        a, b = tmp_path / f"{make.__name__}_a", tmp_path / f"{make.__name__}_b"
        a.mkdir()
        b.mkdir()
        make(a)
        make(b)
        for p in a.rglob("*"):
            if p.is_file():
                assert p.read_bytes() == (b / p.relative_to(a)).read_bytes()
