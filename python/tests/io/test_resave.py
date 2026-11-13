"""Zarr resaving (petakit5d.io.resave, task T01).

Oracle cases against MATLAB: XR_resave_zarr_wrapper (whole files, batches with a box, every
zarr of the folder, uint8 and single), XR_resaveSingleZarr (chunks and batches that don't divide
the image, a box, one batch per task on a cluster's path), resaveZarrBlock (MATLAB's batch
numbers and boxes, a flag that exists, a missing flag folder). Inputs are zarr arrays in MATLAB's
layout; the results (native arrays) are compared by data and metadata after the layout mapping,
the parameter files by content.

Then the resave methods of TiffZarrConversionTest (tests/tools) in intent: chunk changes,
batches with a box, results of other copies and runs kept, batches grouped without a cluster
(along y, then x and z, up to one chunk per thread, within 1/16 of the process's share of the
memory), flags of grouped batches, the same chunk files with and without grouping, and the
workers of a local pool.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d.io import get_zarr_info, read_zarr, resave, write_zarr
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, InPath, Tree, Work
from petakit5d.testing.compare import tree
from petakit5d.testing.synth import CounterRNG

IM_SIZE = (12, 40, 36)  # MATLAB's ImSize [40 36 12]
DATASET = ("scan_CamA_ch0_t0000", "scan_CamB_ch0_t0000", "scan_CamA_ch1_t0000")
CHANNELS = ("CamA_ch0", "CamB_ch0")
BOX = [5, 7, 3, 30, 33, 10]  # MATLAB's [ymin xmin zmin ymax xmax zmax]


def volume(shape, cls, name):
    rng = CounterRNG(0, f"t01/resave/{name}")
    shape = tuple(int(s) for s in shape)
    if cls == "uint8":
        return rng.integers(0, 255, shape, dtype=np.uint8)
    if cls == "uint16":
        return rng.integers(0, 60000, shape, dtype=np.uint16)
    return np.ascontiguousarray((rng.normal(shape) * 100.0).astype(np.float32))


# ---------------------------------------------------------------------------------------------
# oracle cases
# ---------------------------------------------------------------------------------------------


def call_matlab(name, *args):
    from petakit5d import matlab

    getattr(matlab, name)(*args)
    return []


def dataset_tree(cls="uint16", shape=IM_SIZE, chunks=(4, 16, 16)):
    def make(folder):
        for k, n in enumerate(DATASET):
            a = volume(shape, cls, f"ds/{cls}/{k}")
            synth.write_zarr(Path(folder) / f"{n}.zarr", a, chunks, "matlab")

    return make


LOCAL = ["parseCluster", False, "masterCompute", True]


def wrapper_case(cid, make, *args):
    return Case(
        cid,
        code="XR_resave_zarr_wrapper(in{:}); out = {};",
        nargout=0,
        args=[InPath("d"), *args, *LOCAL],
        python=lambda *a: call_matlab("XR_resave_zarr_wrapper", *a),
        inputs={"d": Tree(make, "d")},
    )


def single_case(cid, make, *args, out="out.zarr"):
    return Case(
        cid,
        code="XR_resaveSingleZarr([in{1}, in{2}], [in{1}, in{3}], in{4 : end}); out = {};",
        nargout=0,
        args=[InPath("d"), "scan_CamA_ch0_t0000.zarr", out, *args],
        python=lambda d, n, o, *a: call_matlab("XR_resaveSingleZarr", d + n, d + o, *a),
        inputs={"d": Tree(make, "d")},
    )


def _py_block(d, sz, bs, *args):
    """resaveZarrBlock into an empty array of MATLAB's size sz and chunks bs ([y x z])"""
    from petakit5d.core.axes import from_yxz
    from petakit5d.io import create_zarr

    create_zarr(d + "out.zarr", data_size=from_yxz(sz), block_size=from_yxz(bs))
    os.makedirs(d + "flags", exist_ok=True)
    return call_matlab("resaveZarrBlock", *args)


def block_case(cid, make, sz, bs, inds, inboxes, outboxes, flag="flags/b.mat", error=False):
    code = (
        "d = in{1}; createzarr([d, 'out.zarr'], dataSize=in{2}, blockSize=in{3});"
        " mkdir([d, 'flags']); resaveZarrBlock(in{4 : end}); out = {};"
    )
    return Case(
        cid,
        code=code,
        nargout=0,
        args=[
            InPath("d"),
            sz,
            bs,
            inds,
            Work("d/scan_CamA_ch0_t0000.zarr"),
            Work("d/out.zarr"),
            Work("d/" + flag),
            np.array(inboxes, dtype=np.float64),
            np.array(outboxes, dtype=np.float64),
        ],
        python=_py_block,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=(resave.ResaveError, ValueError, OSError),
        # (a flag is an empty file named .mat, which the .mat comparison can't read; names and
        # sizes of flags are checked by test_resave_zarr_block_checks_and_skips)
        tree=tree(ignore=[r"flags/b\.mat$"]),
    )


def _cases():
    chans = ["channelPatterns", list(CHANNELS)]
    cases = [
        wrapper_case(
            "wrapper.whole_files", dataset_tree(), "blockSize", [8, 8, 4], "largeFile", False
        ),
        wrapper_case(
            "wrapper.batches_box",
            dataset_tree(),
            "batchSize",
            [16, 16, 4],
            "blockSize",
            [8, 8, 4],
            "inputBbox",
            BOX,
            "channelPatterns",
            ["CamA_ch0"],
            "largeFile",
            True,
        ),
        wrapper_case(
            "wrapper.uint8_channels",
            dataset_tree("uint8"),
            "blockSize",
            [16, 8, 5],
            "batchSize",
            [16, 16, 8],
            *chans,
            "resultDirName",
            "rz",
        ),
        wrapper_case("wrapper.single_default_chunks", dataset_tree("single"), "largeFile", False),
    ]
    big = dataset_tree(shape=(23, 61, 47), chunks=(5, 16, 16))
    cases += [
        single_case(
            "single.uneven",
            big,
            "blockSize",
            [16, 12, 5],
            "batchSize",
            [20, 30, 7],
            "parseCluster",
            False,
            "uuid",
            "u1",
        ),
        single_case(
            "single.box",
            big,
            "blockSize",
            [16, 12, 5],
            "batchSize",
            [16, 12, 5],
            "inputBbox",
            [3, 5, 2, 59, 40, 21],
            "parseCluster",
            False,
            "uuid",
            "u2",
        ),
        single_case("single.default_sizes", big, "parseCluster", False),
    ]
    ims = [61, 47, 23]
    cases += [
        block_case(
            "block.two_batches",
            big,
            [40, 30, 12],
            [16, 12, 5],
            [3, 4],
            [[2, 3, 4, 33, 26, 13], [10, 10, 1, 49, 39, 12]],
            [[1, 1, 1, 32, 24, 10], [1, 1, 1, 40, 30, 12]],
        ),
        block_case(
            "block.whole",
            big,
            ims,
            [16, 12, 5],
            [1],
            [[1, 1, 1, *ims]],
            [[1, 1, 1, *ims]],
        ),
        block_case(
            "block.missing_flag_folder",
            big,
            ims,
            [16, 12, 5],
            [1],
            [[1, 1, 1, *ims]],
            [[1, 1, 1, *ims]],
            flag="noflags/b.mat",
            error=True,
        ),
    ]

    def flagged(folder):
        big(folder)
        os.makedirs(Path(folder) / "flags")
        (Path(folder) / "flags" / "b.mat").write_bytes(b"")

    cases.append(
        block_case(
            "block.flag_exists",
            flagged,
            ims,
            [16, 12, 5],
            [1],
            [[1, 1, 1, *ims]],
            [[1, 1, 1, *ims]],
        )
    )
    return cases


CASES = _cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# TiffZarrConversionTest's resave methods in intent
# ---------------------------------------------------------------------------------------------


def write_dataset(folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    ims = []
    for k, n in enumerate(DATASET):
        a = volume(IM_SIZE, "uint16", f"py/{k}")
        ims.append(a)
        write_zarr(a, str(folder / f"{n}.zarr"), block_size=(4, 16, 16))
    return str(folder) + "/", ims, list(DATASET)


def assert_same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype and a.shape == b.shape, (a.dtype, a.shape, b.dtype, b.shape)
    assert a.tobytes() == b.tobytes()


def test_resave_zarr_changes_chunks(tmp_path):
    path, ims, names = write_dataset(tmp_path / "d")
    resave.resave_zarr_dataset(
        path,
        block_size=(4, 8, 8),
        channel_patterns=(".zarr",),
        large_file=False,
        parse_cluster=False,
    )
    for i in range(3):
        fn = f"{path}resaved_zarr/{names[i]}.zarr"
        assert_same(read_zarr(fn), ims[i])
        assert get_zarr_info(fn).block_size == (4, 8, 8)
    for f in ("parameters.mat", "parameters.json"):
        assert os.path.isfile(f"{path}resaved_zarr/{f}")


def test_matlab_layout_inputs_are_resaved_as_native_arrays(tmp_path):
    im = volume(IM_SIZE, "single", "layout")
    synth.write_zarr(tmp_path / "m.zarr", im, (4, 16, 16), "matlab")
    resave.resave_zarr(str(tmp_path / "m.zarr"), str(tmp_path / "n.zarr"), parse_cluster=False)
    info = get_zarr_info(str(tmp_path / "n.zarr"))
    assert (info.order, info.size, info.block_size) == ("C", IM_SIZE, IM_SIZE)
    assert_same(read_zarr(str(tmp_path / "n.zarr")), im)


@pytest.fixture
def framework_spy(monkeypatch):
    """records run_tasks calls; ``on_call`` runs at each call (MATLAB's FrameworkSpyFixture with
    a function); the tasks don't run, the outputs count as done"""
    from petakit5d.jobs import framework

    spy = {"calls": [], "on_call": None}

    def fake(inputs, outputs, tasks, **kw):
        spy["calls"].append(
            {"inputs": list(inputs), "outputs": list(outputs), "tasks": list(tasks), **kw}
        )
        if spy["on_call"] is not None:
            spy["on_call"]()
        return [True] * len(list(outputs))

    monkeypatch.setattr(framework, "run_tasks", fake)
    return spy


def test_resave_keeps_result_of_another_copy(tmp_path, framework_spy):
    # two copies of a task share the temporary result; the copy that finishes second must not
    # delete the result the first moved into place
    path, _, names = write_dataset(tmp_path / "d")
    out = f"{path}resaved.zarr"
    framework_spy["on_call"] = lambda: os.rename(out + "_test", out)
    resave.resave_zarr(
        f"{path}{names[0]}.zarr", out, block_size=(4, 8, 8), parse_cluster=False, uuid="test"
    )
    assert os.path.isdir(out), "the result of the other copy was deleted"


def test_resave_keeps_result_of_another_run(tmp_path, framework_spy):
    # a copy of another run (its own uuid) that finds the result published keeps it
    path, ims, names = write_dataset(tmp_path / "d")
    out = f"{path}resaved.zarr"
    framework_spy["on_call"] = lambda: write_zarr(ims[0], out, block_size=(4, 8, 8))
    resave.resave_zarr(
        f"{path}{names[0]}.zarr", out, block_size=(4, 8, 8), parse_cluster=False, uuid="test"
    )
    assert_same(read_zarr(out), ims[0])
    assert not os.path.isdir(out + "_test"), "the temporary result was left behind"


def test_resave_large_zarr_in_batches(tmp_path):
    path, ims, names = write_dataset(tmp_path / "d")
    resave.resave_zarr_dataset(
        path,
        batch_size=(4, 16, 16),
        block_size=(4, 8, 8),
        input_bbox=(2, 4, 6, 10, 30, 33),
        channel_patterns=CHANNELS[:1],
        large_file=True,
        parse_cluster=False,
    )
    out = read_zarr(f"{path}resaved_zarr/{names[0]}.zarr")
    assert_same(out, ims[0][2:10, 4:30, 6:33])


def _regions(monkeypatch, src, out, batch, parse_cluster):
    """the output regions of the tasks resave_zarr hands to the framework (chunks (4, 8, 8)),
    one row per read and write (MATLAB's resaveRegions)"""
    from petakit5d.jobs import framework

    calls = []
    monkeypatch.setattr(
        framework, "run_tasks", lambda i, o, t, **kw: calls.append(t) or [True] * len(o)
    )
    resave.resave_zarr(
        src, out, block_size=(4, 8, 8), batch_size=batch, parse_cluster=parse_cluster
    )
    return np.array([r for t in calls[0] for r in t.args[5]], dtype=np.int64)


def _verify_tiles(regions, batch):
    counts = np.zeros(IM_SIZE, np.int64)
    for r in regions:
        counts[r[0] : r[3], r[1] : r[4], r[2] : r[5]] += 1
    assert (counts == 1).all(), "the regions do not cover the image once"
    assert (regions[:, :3] % np.array(batch) == 0).all(), "a region starts inside a batch"
    ends = regions[:, 3:]
    ok = (ends % np.array(batch) == 0) | (ends == np.array(IM_SIZE))
    assert ok.all(), "a region ends inside a batch"


def _unique_sizes(regions):
    return sorted({tuple(int(v) for v in r[3:] - r[:3]) for r in regions})


def test_resave_groups_batches_without_cluster(tmp_path, monkeypatch):
    # neighboring batches together (along y, then x and z, up to one chunk per thread); cluster
    # tasks, and too little or unknown memory, keep one batch per read and write
    path, _, names = write_dataset(tmp_path / "d")
    src = f"{path}{names[0]}.zarr"
    monkeypatch.setattr(resave, "_memory_gb", lambda: 64.0)
    monkeypatch.setattr(resave, "_node_cores", lambda: 24.0)

    # 12 threads, 3 x 5 x 3 batches of two chunks (4, 8, 16): up to 6, 3 along y, 2 along x
    monkeypatch.setattr(resave, "_comp_threads", lambda: 12)
    regions = _regions(monkeypatch, src, f"{path}local2.zarr", (4, 16, 8), False)
    assert regions.shape == (9, 6)
    assert _unique_sizes(regions) == [(4, 40, 4), (4, 40, 16)]
    _verify_tiles(regions, (4, 16, 8))

    # 6 threads, 3 x 5 x 5 batches of (4, 8, 8): 15 groups of the 5 batches along y
    monkeypatch.setattr(resave, "_comp_threads", lambda: 6)
    regions = _regions(monkeypatch, src, f"{path}local.zarr", (4, 8, 8), False)
    assert regions.shape == (15, 6)
    assert _unique_sizes(regions) == [(4, 40, 4), (4, 40, 8)]
    _verify_tiles(regions, (4, 8, 8))

    regions = _regions(monkeypatch, src, f"{path}cluster.zarr", (4, 8, 8), True)
    assert regions.shape == (75, 6)
    _verify_tiles(regions, (4, 8, 8))

    for mem in (math.nan, 1e-6):
        monkeypatch.setattr(resave, "_memory_gb", lambda m=mem: m)
        regions = _regions(monkeypatch, src, f"{path}low_{mem}.zarr", (4, 8, 8), False)
        assert regions.shape == (75, 6), mem


def test_resave_groups_within_the_share_of_the_memory(tmp_path, monkeypatch):
    # at most 1/16 of this process's share of the memory: the memory times its share of the
    # node's CPUs
    path, _, names = write_dataset(tmp_path / "d")
    src = f"{path}{names[0]}.zarr"
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "24")
    # 64 kB: 1/16 is 8 batches of (4, 8, 8) uint16
    monkeypatch.setattr(resave, "_memory_gb", lambda: 2**16 / 2**30)
    monkeypatch.setattr(resave, "_comp_threads", lambda: 24)
    regions = _regions(monkeypatch, src, f"{path}all.zarr", (4, 8, 8), False)
    assert regions.shape == (15, 6)
    # 6 of 24 CPUs: a quarter of the memory, 2 batches
    monkeypatch.setattr(resave, "_comp_threads", lambda: 6)
    regions = _regions(monkeypatch, src, f"{path}quarter.zarr", (4, 8, 8), False)
    assert regions.shape == (45, 6)
    assert _unique_sizes(regions) == [(4, 8, 4), (4, 8, 8), (4, 16, 4), (4, 16, 8)]
    _verify_tiles(regions, (4, 8, 8))


def test_local_batch_size_reads_the_node_and_meminfo(monkeypatch):
    monkeypatch.delenv("SLURM_CPUS_ON_NODE", raising=False)
    assert resave._node_cores() == float(os.cpu_count() or 1)
    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "7")
    assert resave._node_cores() == 7.0
    m = resave._memory_gb()
    assert math.isnan(m) or m > 0


def test_resave_grouped_batches_do_not_take_other_flags(tmp_path, monkeypatch):
    # a run that groups batches, with the uuid of an interrupted run that wrote one batch at a
    # time (and finished its first task, batches 1-5), must not take that task's flag for its
    # own first task (groups 1-5)
    from petakit5d.io import create_zarr

    path, ims, names = write_dataset(tmp_path / "d")
    out = f"{path}resumed.zarr"
    create_zarr(out + "_u", data_size=IM_SIZE, block_size=(4, 8, 8))
    write_zarr(ims[0][:4, :, :8], out + "_u", bbox=(0, 0, 0, 4, 40, 8), create=False)
    flag_path = f"{path}zarr_flag/resumed_u/"
    os.makedirs(flag_path)
    open(flag_path + "batches_1_5.mat", "w").close()
    monkeypatch.setattr(resave, "_memory_gb", lambda: 64.0)
    monkeypatch.setattr(resave, "_comp_threads", lambda: 6)
    resave.resave_zarr(
        f"{path}{names[0]}.zarr",
        out,
        block_size=(4, 8, 8),
        batch_size=(4, 8, 8),
        parse_cluster=False,
        uuid="u",
    )
    assert_same(read_zarr(out), ims[0])
    assert not os.path.isdir(flag_path), "the flag folder was left"


def _chunk_files(fn):
    return sorted(n for n in os.listdir(fn) if os.path.isfile(os.path.join(fn, n)))


def test_resave_grouped_batches_write_the_same_chunks(tmp_path, monkeypatch):
    # batches read and written together give the same chunk files as one batch at a time, also
    # with a crop and batches of two chunks that don't divide the image
    path, ims, names = write_dataset(tmp_path / "d")
    monkeypatch.setattr(resave, "_memory_gb", lambda: 64.0)
    for k, threads in enumerate((1, 6)):
        monkeypatch.setattr(resave, "_comp_threads", lambda t=threads: t)
        resave.resave_zarr_dataset(
            path,
            result_dir_name=f"resaved_{k}",
            batch_size=(4, 16, 8),
            block_size=(4, 8, 8),
            input_bbox=(2, 4, 6, 10, 30, 33),
            channel_patterns=CHANNELS[:1],
            large_file=True,
            parse_cluster=False,
        )
    fns = [f"{path}resaved_{k}/{names[0]}.zarr" for k in range(2)]
    assert_same(read_zarr(fns[1]), ims[0][2:10, 4:30, 6:33])
    files = [_chunk_files(f) for f in fns]
    assert files[0] == files[1]
    assert len(files[0]) == 33  # 2 x 4 x 4 chunks and the metadata
    for f in files[0]:
        a = Path(fns[0], f).read_bytes()
        assert a == Path(fns[1], f).read_bytes(), f


def test_cluster_batches_keep_matlabs_flag_names(tmp_path, framework_spy):
    # one batch per read and write, tasks of five batches, flags named by MATLAB's 1-based
    # batch numbers; 1.5 times a batch per task, then master compute for the second call
    path, _, names = write_dataset(tmp_path / "d")
    resave.resave_zarr(
        f"{path}{names[0]}.zarr",
        f"{path}c.zarr",
        block_size=(4, 8, 8),
        batch_size=(4, 8, 8),
        parse_cluster=True,
        uuid="v",
        cpus_per_task=2,
    )
    c = framework_spy["calls"][0]
    assert [os.path.basename(o) for o in c["outputs"]][:2] == [
        "batches_1_5.mat",
        "batches_6_10.mat",
    ]
    assert len(c["outputs"]) == 15
    assert c["final_out_fullpath"] == f"{path}c.zarr"
    assert c["cpus_per_task"] == 2
    assert c["mem_allocate"] == pytest.approx(4 * 8 * 8 * 2 / 2**30 * 1.5)
    assert "master_compute" not in c


def test_grouped_flags_name_matlabs_group_size(tmp_path, framework_spy, monkeypatch):
    path, _, names = write_dataset(tmp_path / "d")
    monkeypatch.setattr(resave, "_memory_gb", lambda: 64.0)
    monkeypatch.setattr(resave, "_comp_threads", lambda: 6)
    resave.resave_zarr(
        f"{path}{names[0]}.zarr",
        f"{path}g.zarr",
        block_size=(4, 8, 8),
        batch_size=(4, 8, 8),
        parse_cluster=False,
        uuid="w",
    )
    # groups of (4, 40, 8): MATLAB's [y x z] 40 8 4
    flags = [os.path.basename(o) for o in framework_spy["calls"][0]["outputs"]]
    assert flags[0] == "batches_40_8_4_1_5.mat"


def test_resave_wrapper_job_settings(tmp_path, framework_spy):
    # few batches: whole files, 1.9 times a file; many: batches of 1.9 times a batch
    path, _, _ = write_dataset(tmp_path / "d")
    resave.resave_zarr_dataset(path, large_file=False, parse_cluster=False, cpus_per_task=3)
    resave.resave_zarr_dataset(
        path, large_file=True, batch_size=(4, 8, 8), parse_cluster=False, result_dir_name="r2"
    )
    a, b = framework_spy["calls"]
    assert a["mem_allocate"] == pytest.approx(np.prod(IM_SIZE) * 2 / 2**30 * 1.9)
    assert a["cpus_per_task"] == 3
    assert a["tasks"][0].kwargs["batch_size"] == list(IM_SIZE)
    assert a["tasks"][0].kwargs["cpus_per_task"] == 3
    assert b["mem_allocate"] == pytest.approx(4 * 8 * 8 * 2 / 2**30 * 1.9)
    # every file of a call shares its uuid
    assert len({t.kwargs["uuid"] for t in a["tasks"]}) == 1


def test_resave_zarr_block_checks_and_skips(tmp_path, capsys):
    from petakit5d.io import create_zarr

    im = volume(IM_SIZE, "uint16", "block")
    write_zarr(im, str(tmp_path / "in.zarr"), block_size=(4, 16, 16))
    create_zarr(str(tmp_path / "out.zarr"), data_size=IM_SIZE, block_size=(4, 8, 8))
    box = [0, 0, 0, *IM_SIZE]
    with pytest.raises(resave.ResaveError, match="does not exist"):
        resave.resave_zarr_block(
            [0], str(tmp_path / "in.zarr"), str(tmp_path / "out.zarr"), str(tmp_path / "f/b.mat"),
            [box], [box],
        )  # fmt: skip
    os.makedirs(tmp_path / "f")
    flag = str(tmp_path / "f" / "b.mat")
    with pytest.raises(resave.ResaveError, match="doesnot exist"):
        resave.resave_zarr_block(
            [0], str(tmp_path / "no.zarr"), str(tmp_path / "out.zarr"), flag, [box], [box]
        )
    resave.resave_zarr_block(
        [0], str(tmp_path / "in.zarr"), str(tmp_path / "out.zarr"), flag, [box], [box]
    )
    assert os.path.isfile(flag) and os.path.getsize(flag) == 0
    assert_same(read_zarr(str(tmp_path / "out.zarr")), im)
    resave.resave_zarr_block(
        [0, 1], str(tmp_path / "in.zarr"), str(tmp_path / "out.zarr"), flag, [box], [box]
    )
    assert "The batch files (1 - 2) already exist, skip them!" in capsys.readouterr().out


@pytest.mark.slow
@pytest.mark.skipif(not hasattr(os, "sched_getaffinity"), reason="needs CPU affinity (Linux)")
def test_resave_on_local_pool_groups_within_the_workers_threads(tmp_path, capfd):
    # files resaved on a local process pool (localWorkerNum in the configuration file): each
    # worker has its share of the threads, and groups no more batches than that
    import re

    cpus = len(os.sched_getaffinity(0))
    if cpus < 2:
        pytest.skip("needs 2 CPUs")
    path, ims, names = write_dataset(tmp_path / "d")
    cfg = tmp_path / "config.json"
    cfg.write_text('{"localWorkerNum": 2}', newline="\n")
    resave.resave_zarr_dataset(
        path,
        batch_size=(4, 8, 8),
        block_size=(4, 8, 8),
        large_file=True,
        parse_cluster=False,
        config_file=str(cfg),
    )
    for i in range(3):
        assert_same(read_zarr(f"{path}resaved_zarr/{names[i]}.zarr"), ims[i])
    out = capfd.readouterr().out
    threads = max(1, cpus // 2)
    for g in re.findall(r"batches_(\d+)_(\d+)_(\d+)_\d+_\d+", out):
        n = np.prod(np.array([int(v) for v in g]) / np.array([8, 8, 4]))
        assert n <= threads, f"a worker grouped {n} batches"
