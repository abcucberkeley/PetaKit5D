"""Deskew/rotate of large zarr files in blocks (petakit5d.deskew_rotate.zarr, task D04).

Oracle cases against MATLAB: ``XR_deskewRotateZarr`` on a copy of a folder with a zarr file in
MATLAB's layout (batches and chunks that don't divide the result, a task per batch, no MIPs,
single output, an existing result, a box, resample factors of one, two and three values with x
and y different, flipZstack, the error of objective scans, parseParfor on a local pool), its
folder compared with the port's; ``XR_getDeskeRotateBoxesFromMasks`` on MIP masks (with and
without an input box, resample factors of none, one, two and three values, flipZstack).

Then DeskewRotateDataWrapperTest's large-file methods and DeskewRotateMemoryTest's
largeFileBlocksRequestOneBatch in intent: the MIP option of the block tasks, failed blocks not
published, a result another writer published kept (before and during the MIPs), the blocks of
an ended task skipped (a fake Slurm), the wait for a writer of another call, the memory of one
batch; and the port's own rules (the geometry, the MIPs of blocks from MATLAB's layout).
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from petakit5d.deskew_rotate import zarr as dz
from petakit5d.io import create_zarr, read_zarr, write_zarr
from petakit5d.testing import oracle, slurm_fake, synth
from petakit5d.testing.cases import Box, Case, InPath, PerAxis, Tree
from petakit5d.testing.compare import tree
from petakit5d.testing.synth import CounterRNG

SKEW = 32.45
PX = 0.108
DZ = 0.3
RAW = (29, 37, 61)
FSN = "scan_Iter_0000_CamA_ch0_stack0000_0000000msec_0000t"
ERRORS = (dz.DeskewRotateZarrError, RuntimeError, ValueError, OSError, IndexError)


def make_frame(shape=RAW, kind="uint16", seed=0):
    nz, ny, nx = shape
    rng = CounterRNG(seed, f"dsrzarr/{kind}/{'x'.join(map(str, shape))}")
    pos = [
        (2 + (k * 7) % max(1, nz - 4), 3 + (k * 5) % max(1, ny - 6), 4 + (k * 11) % (nx - 8))
        for k in range(9)
    ]
    vol = synth.beads(shape, pos, 1.5, 1000.0, 100.0)
    noise = rng.integers(0, 50, shape, dtype=np.uint16)
    if kind == "single":
        return np.ascontiguousarray(vol + noise.astype(np.float32))
    return np.ascontiguousarray(np.floor(vol).astype(np.uint16) + noise)


# ------------------------------------------------------------------------------------------
# oracle cases
# ------------------------------------------------------------------------------------------


def zarr_tree(kind="uint16", existing=False):
    def make(folder):
        synth.write_zarr(Path(folder) / f"{FSN}.zarr", make_frame(RAW, kind), (8, 16, 16), "matlab")
        if existing:
            synth.write_zarr(
                Path(folder) / "DSR" / f"{FSN}.zarr", np.full((3, 4, 5), 7, np.uint16), (3, 4, 5)
            )

    return make


def mask_tree(folder, shape=RAW, x0=15):
    """the z, y and x MIPs (MATLAB's orientation) of a mask of two boxes"""
    m = np.zeros(shape, dtype=np.uint8)
    m[8:22, 10:30, x0 : x0 + 30] = 1
    m[3, 30:33, x0 + 35 : x0 + 40] = 1
    synth.write_tiff(Path(folder) / "mask_z.tif", m.max(axis=0))
    synth.write_tiff(Path(folder) / "mask_y.tif", np.ascontiguousarray(m.max(axis=1).T))
    synth.write_tiff(Path(folder) / "mask_x.tif", np.ascontiguousarray(m.max(axis=2).T))


def p_zarr(d, *args):
    from petakit5d import matlab

    matlab.XR_deskewRotateZarr(d + f"{FSN}.zarr", *args)
    return []


def zarr_case(cid, make, *args, error=False, code=None):
    return Case(
        f"zarr.{cid}",
        code=code or f"XR_deskewRotateZarr([in{{1}}, '{FSN}.zarr'], in{{2 : end}}); out = {{}};",
        nargout=0,
        args=[
            InPath("d"),
            PX,
            DZ,
            "skewAngle",
            SKEW,
            "reverse",
            True,
            "parseCluster",
            False,
            *args,
        ],  # fmt: skip
        python=p_zarr,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=ERRORS,
        # (a temporary result left by an error is an empty array: its .zarray is the port's
        # native one, which the tree comparison maps only for names ending in .zarr)
        tree=tree(
            masks=[(r"owner=[^\n]*", "owner=<owner>")],
            ignore=[r"\.zarr_[^/]*/\.zarray$"],
        ),
    )


# parseParfor runs the tasks on a pool of 12 (matlab_parfor_generic_computing_wrapper), which
# MATLAB refuses in a job of fewer cores (parallel:cluster:LocalProfileNumWorkersExceeded),
# and the wrapper takes an open pool as it is: the case opens the pool itself, of 12 workers
# at most, so that its reference is the blocks' result in a job of any size, and closes it
_PARFOR = (
    "c = parcluster('local'); p = parpool(c, min(12, c.NumWorkers)); "
    f"try, XR_deskewRotateZarr([in{{1}}, '{FSN}.zarr'], in{{2 : end}}); "
    "catch ME, delete(p); rethrow(ME); end; delete(p); out = {};"
)


def parfor_case(make, *args):
    return zarr_case("parfor", make, *args, "parseParfor", True, code=_PARFOR)


def p_masks(d, border, px, dzv, skew, rev, resample=None, box=None, flip=False):
    names = [d + n for n in ("mask_z.tif", "mask_y.tif", "mask_x.tif")]
    rs = None if resample is None else [float(v) for v in np.asarray(resample).reshape(-1)]
    rs = rs or None
    return dz.deskew_rotate_boxes_from_masks(names, border, px, dzv, skew, rev, rs, box, flip)


def mask_case(cid, border, *extra):
    return Case(
        f"masks.{cid}",
        code=(
            "f = cellfun(@(x) [in{1}, x], {'mask_z.tif', 'mask_y.tif', 'mask_x.tif'}, 'unif', 0);"
            " [o1, o2, o3] = XR_getDeskeRotateBoxesFromMasks(f, in{2 : end}); out = {o1, o2, o3};"
        ),
        nargout=3,
        args=[InPath("d"), PerAxis(border, "yxz"), PX, DZ, SKEW, True, *extra],
        python=p_masks,
        out_roles=("box", "box", None),
        inputs={"d": Tree(mask_tree, "d")},
        files=False,
    )


SMALL = ["batchSize", [16, 1024, 1024], "blockSize", [16, 64, 64]]


def _cases():
    t = zarr_tree()
    return [
        zarr_case("default", t, *SMALL),
        zarr_case("uneven", t, "batchSize", [10, 1024, 1024], "blockSize", [4, 24, 20]),
        zarr_case("task_per_batch", t, *SMALL, "taskSize", 1),
        zarr_case("no_mip", t, *SMALL, "saveMIP", False),
        zarr_case("single_out", t, *SMALL, "save16bit", False),
        zarr_case("single_in", zarr_tree("single"), *SMALL),
        zarr_case("bbox", t, *SMALL, "inputBbox", [3, 5, 2, 30, 58, 27]),
        zarr_case("result_dir", t, *SMALL, "resultDirStr", "out/"),
        zarr_case("existing", zarr_tree(existing=True), *SMALL),
        zarr_case("objective", t, *SMALL, "objectiveScan", True, error=True),
        parfor_case(t, *SMALL),
        # the factor is [x y z], 1 or 2 values completed to 3 (it was read [y x z], and 1 or 2
        # values raised an index error)
        zarr_case("resample_two", t, *SMALL, "resampleFactor", [1.5, 1.5]),
        zarr_case("resample_one", t, *SMALL, "resampleFactor", 2),
        zarr_case("resample_xz", t, *SMALL, "resampleFactor", [1.5, 1, 1.25]),
        zarr_case("resample_yx", t, *SMALL, "resampleFactor", [1, 2, 1]),
        zarr_case("resample_y_frac", t, *SMALL, "resampleFactor", [1, 1.5, 1.25]),
        # each batch flipped along z after its read (the blocks never flipped)
        zarr_case("flip", t, *SMALL, "flipZstack", True),
        zarr_case("flip_bbox", t, *SMALL, "flipZstack", True, "inputBbox", [3, 5, 2, 30, 58, 27]),
        zarr_case("flip_resample", t, *SMALL, "flipZstack", True, "resampleFactor", [1.5, 1, 1.25]),
        mask_case("default", [10, 10, 10], [1, 1, 1]),
        mask_case("box", [10, 10, 10], [1, 1, 1], Box([5, 8, 4, 33, 55, 25])),
        mask_case("resample", [7, 7, 7], [1.5, 1.5, 1.5]),
        mask_case("resample_xy", [10, 12, 14], [1, 1.5, 1.25]),
        mask_case("resample_xz", [10, 12, 14], [1.5, 1, 1.25]),
        mask_case("resample_one", [10, 12, 14], 1.5),
        mask_case("resample_two", [10, 12, 14], [1.5, 1]),
        mask_case("resample_none", [10, 12, 14], []),
        # the y MIP flipped along z after its crop (the input box stays the raw data's)
        mask_case("flip", [10, 10, 10], [1, 1, 1], [], True),
        mask_case("flip_resample", [10, 12, 14], [1, 1.5, 1.25], [], True),
        mask_case("flip_box", [10, 10, 10], [1, 1, 1], Box([5, 8, 4, 33, 55, 25]), True),
    ]


CASES = _cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ------------------------------------------------------------------------------------------
# DeskewRotateDataWrapperTest's large-file methods in intent
# ------------------------------------------------------------------------------------------

KW = {
    "skew_angle": SKEW,
    "reverse": True,
    "batch_size": (1024, 16, 1024),
    "block_size": (64, 16, 64),
    "parse_cluster": False,
}


def write_input(folder, kind="uint16", shape=RAW):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    fn = str(folder / f"{FSN}.zarr")
    frame = make_frame(shape, kind)
    write_zarr(frame, fn, block_size=(8, 16, 16))
    return fn, frame


def same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.dtype == b.dtype and a.shape == b.shape, (a.dtype, b.dtype, a.shape, b.shape)
    np.testing.assert_array_equal(a, b)


@pytest.fixture
def framework_spy(monkeypatch):
    """records run_tasks calls; nothing runs. ``spy.action`` runs at each call, ``spy.done``
    are the done flags returned"""
    from petakit5d.jobs import framework

    calls = []

    def fake(inputs, outputs, tasks, **kw):
        calls.append({"inputs": list(inputs), "outputs": list(outputs), "tasks": list(tasks), **kw})
        if fake.action is not None:
            fake.action()
        return [fake.done] * len(list(outputs))

    fake.action = None
    fake.done = True
    monkeypatch.setattr(framework, "run_tasks", fake)
    fake.calls = calls
    return fake


def test_large_file_blocks_get_the_mip_option_only_for_mips(tmp_path, framework_spy):
    fn, _ = write_input(tmp_path)
    framework_spy.done = False
    for save in (False, True):
        with pytest.raises(dz.DeskewRotateZarrError):
            dz.deskew_rotate_zarr(
                fn, PX, DZ, save_mip=save, result_dir_str=f"DSR_{int(save)}/", **KW
            )
    c0, c1 = framework_spy.calls[0], framework_spy.calls[2]
    assert all("mip_batch_size" not in t.kwargs for t in c0["tasks"])
    assert all(len(t.kwargs["mip_batch_size"]) == 3 for t in c1["tasks"])
    # 3 batches of 16, 16 and 5 rows, 2 per task: MATLAB's flag names, 1-based
    assert [os.path.basename(o) for o in c1["outputs"]] == ["blocks_1_2.mat", "blocks_3_3.mat"]
    assert [t.args[0] for t in c1["tasks"]] == [[0, 1], [2]]


def test_large_file_with_failed_blocks_is_not_published(tmp_path, framework_spy):
    fn, _ = write_input(tmp_path)
    framework_spy.done = False
    with pytest.raises(dz.DeskewRotateZarrError, match="Some blocks are not finished"):
        dz.deskew_rotate_zarr(fn, PX, DZ, save_mip=False, **KW)
    assert not os.path.isdir(tmp_path / "DSR" / f"{FSN}.zarr"), "an incomplete result was published"
    records = list((tmp_path / "DSR" / "zarr_flag").glob(f"{FSN}_*/writer.txt"))
    assert len(records) == 1 and "\nended=1" in records[0].read_text()
    assert len(framework_spy.calls) == 2, "two rounds of tries"


def test_large_file_keeps_a_result_another_writer_published(tmp_path, framework_spy):
    fn, _ = write_input(tmp_path)
    dsr = str(tmp_path / "DSR" / f"{FSN}.zarr")
    other = np.full((5, 5, 5), 7, dtype=np.uint16)
    framework_spy.action = lambda: write_zarr(other, dsr)
    dz.deskew_rotate_zarr(fn, PX, DZ, **KW)
    same(read_zarr(dsr), other)
    assert not list((tmp_path / "DSR").glob(f"{FSN}.zarr_*")), "the temporary result is left"
    assert not list((tmp_path / "DSR" / "zarr_flag").glob(f"{FSN}_*")), "the flags are left"
    assert not list((tmp_path / "DSR").glob("MIPs/*.tif")), "MIPs of the incomplete result"
    assert os.path.normpath(framework_spy.calls[0]["final_out_fullpath"]) == os.path.normpath(dsr)


def test_large_file_keeps_a_result_published_during_its_mips(tmp_path, monkeypatch):
    from petakit5d.tools import mip

    fn, _ = write_input(tmp_path)
    dsr = str(tmp_path / "DSR" / f"{FSN}.zarr")
    other = np.full((5, 5, 5), 7, dtype=np.uint16)
    monkeypatch.setattr(mip, "mip_zarr", lambda *a, **k: write_zarr(other, dsr))
    dz.deskew_rotate_zarr(fn, PX, DZ, **KW)
    same(read_zarr(dsr), other)
    assert not list((tmp_path / "DSR").glob(f"{FSN}.zarr_*"))
    assert not [n for n in os.listdir(dsr) if "zarr_" in n], "moved into the other's"


fake_slurm = slurm_fake.fake_slurm  # (the fixture)


@pytest.mark.skipif(os.name == "nt", reason="The fake Slurm needs a POSIX shell.")
def test_large_file_blocks_of_an_ended_task_are_skipped(tmp_path, framework_spy, fake_slurm):
    from petakit5d.jobs import writers

    writers.reset()
    fn, _ = write_input(tmp_path)
    framework_spy.done = False
    with pytest.raises(dz.DeskewRotateZarrError):
        dz.deskew_rotate_zarr(fn, PX, DZ, save_mip=False, **KW)
    task = framework_spy.calls[0]["tasks"][0]
    flag = framework_spy.calls[0]["outputs"][0]
    record = Path(os.path.dirname(flag)) / "writer.txt"
    assert record.is_file(), "the task did not record itself"
    # a task whose job was cancelled can't say so: its job is gone
    txt = record.read_text()
    lines = [ln for ln in txt.split("\n") if ln and not ln.startswith("ended=")]
    lines = ["owner=777_3:othernode:5" if ln.startswith("owner=") else ln for ln in lines]
    record.write_text("\n".join(lines) + "\n", newline="")
    writers.reset()
    from petakit5d.jobs.tasks import run_record

    assert run_record(task.record()) is None
    assert not os.path.isfile(flag), "the block of an ended task was deskewed"
    writers.reset()


@pytest.mark.skipif(
    not sys.platform.startswith("linux"), reason="Writers outside Slurm are checked on Linux only."
)
def test_large_file_waits_for_a_writer_of_another_call(tmp_path, framework_spy, monkeypatch):
    from petakit5d.jobs import files, writers

    for v in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID"):
        monkeypatch.delenv(v, raising=False)
    writers.reset()
    fn, _ = write_input(tmp_path)
    dsr = f"{tmp_path}/DSR//{FSN}.zarr"
    proc = subprocess.Popen(["sleep", "120"])
    try:
        flag = tmp_path / "DSR" / "zarr_flag" / f"{FSN}_otheruuid"
        flag.mkdir(parents=True)
        owner = f"0_:{writers.node_name()}:{proc.pid}:{writers.process_start(proc.pid)}"
        rec = flag / "writer.txt"
        rec.write_text(
            f"result={writers.normalize_result_path(dsr)}\nuuid=otheruuid\nowner={owner}\n",
            newline="",
        )
        st = os.stat(rec)
        os.utime(rec, ns=(st.st_atime_ns, st.st_mtime_ns - 60 * 10**9))
        other = np.full((5, 5, 5), 7, dtype=np.uint16)

        def spy(groups, wait_time):
            write_zarr(other, dsr)
            return False

        monkeypatch.setattr(files, "wait_for_files", spy)
        dz.deskew_rotate_zarr(fn, PX, DZ, **KW)
    finally:
        proc.kill()
        proc.wait()
    assert framework_spy.calls == [], "the result was made a second time"
    same(read_zarr(dsr), other)


def test_large_file_blocks_request_one_batch(tmp_path, framework_spy):
    # the block tasks hold one batch (input rows with a one-row border, the output region)
    from petakit5d.deskew_rotate.memory import estimate_deskew_rotate_memory

    fn = str(tmp_path / f"{FSN}.zarr")
    write_zarr(np.zeros((80, 64, 96), np.uint16), fn, block_size=(64, 16, 64))
    dz.deskew_rotate_zarr(fn, PX, DZ, save_mip=False, master_compute=True, **KW)
    # 4 batches of 16 output rows; the middle ones read 18 input rows; chunks (z, 16, 64) cut
    # to the 47 output slices
    _, _, out = estimate_deskew_rotate_memory((80, 64, 96), PX, DZ, skew_angle=SKEW)
    expected = estimate_deskew_rotate_memory(
        (80, 18, 96),
        PX,
        DZ,
        skew_angle=SKEW,
        out_size=(out[0], 16, out[2]),
        large_file=True,
        block_size=(out[0], 16, 64),
    )[0]
    assert out[0] == 47
    assert framework_spy.calls[0]["mem_allocate"] == pytest.approx(expected, rel=1e-12)


def test_large_file_skips_an_existing_result(tmp_path, framework_spy, capsys):
    fn, _ = write_input(tmp_path)
    write_zarr(np.ones((2, 2, 2), np.uint16), str(tmp_path / "DSR" / f"{FSN}.zarr"))
    assert dz.deskew_rotate_zarr(fn, PX, DZ, **KW) is None
    assert "The output result exists, skip it!" in capsys.readouterr().out
    assert framework_spy.calls == []


def test_a_temporary_result_of_other_chunks_is_made_again(tmp_path, framework_spy):
    fn, _ = write_input(tmp_path)
    tmp = str(tmp_path / "DSR" / f"{FSN}.zarr_u1")
    create_zarr(tmp, data_size=(5, 5, 5), block_size=(5, 5, 5), dtype="uint16")
    framework_spy.done = False
    with pytest.raises(dz.DeskewRotateZarrError):
        dz.deskew_rotate_zarr(fn, PX, DZ, uuid="u1", **KW)
    from petakit5d.io import get_zarr_info

    g = dz.large_file_geometry((37, 61, 29), PX, DZ, SKEW, [16, 1024, 1024], [16, 64, 64], None,
                               None, None)  # fmt: skip
    info = get_zarr_info(tmp)
    want = g.out_block_size
    assert tuple(info.block_size) == (want[2], want[0], want[1])
    assert tuple(info.size) == (g.out_size[2], g.out_size[0], g.out_size[1])


# ------------------------------------------------------------------------------------------
# the port's own rules
# ------------------------------------------------------------------------------------------


def test_geometry_has_matlabs_batches_and_borders():
    g = dz.large_file_geometry((37, 61, 29), PX, DZ, SKEW, [16, 1024, 1024], [16, 64, 64], None,
                               None, None)  # fmt: skip
    # three batches of rows 1-16, 17-32 and 33-37 with a border row each side within the image
    assert g.out_batch_bboxes[:, [0, 3]].tolist() == [[1, 17], [16, 33], [32, 37]]
    assert g.out_region_bboxes[:, [0, 3]].tolist() == [[1, 16], [17, 32], [33, 37]]
    assert g.in_batch_bboxes[:, [0, 3]].tolist() == [[1, 17], [16, 33], [32, 37]]
    assert g.out_local_bboxes[:, [0, 3]].tolist() == [[1, 16], [2, 17], [2, 6]]
    assert all(g.in_batch_bboxes[:, 4] == 61) and all(g.in_batch_bboxes[:, 5] == 29)
    # a batch is all of x and z, whatever the batch size says
    g2 = dz.large_file_geometry((37, 61, 29), PX, DZ, SKEW, [16, 32, 8], [16, 64, 64], None,
                                None, None)  # fmt: skip
    assert g2.out_batch_bboxes[:, [0, 3]].tolist() == [[1, 17], [16, 33], [32, 37]]
    assert g2.out_batch_size[1:] == g2.out_size[1:]


def test_geometry_reads_the_resample_factor_as_x_y_z():
    # (MATLAB: largeFileReadsTheResampleFactorAsRegularFiles) rs(2), y, sets the rows of a
    # batch's input and divides the result's rows; x and z divide its x and z
    g = dz.large_file_geometry((37, 61, 29), PX, DZ, SKEW, [8, 1024, 1024], [8, 64, 64],
                               [1, 2, 1], None, None)  # fmt: skip
    g1 = dz.large_file_geometry((37, 61, 29), PX, DZ, SKEW, [8, 1024, 1024], [8, 64, 64],
                                None, None, None)  # fmt: skip
    assert g.dsr_out_size == [19, *g1.dsr_out_size[1:]]  # round(37 / 2)
    assert g.in_batch_bboxes[1, 0] == 15  # round((8 - 1) * 2) + 1
    gx = dz.large_file_geometry((37, 61, 29), PX, DZ, SKEW, [8, 1024, 1024], [8, 64, 64],
                                [2, 1, 1], None, None)  # fmt: skip
    assert gx.dsr_out_size[0] == 37 and gx.dsr_out_size[1] == round(g1.dsr_out_size[1] / 2)
    assert gx.in_batch_bboxes[1, 0] == 8


@pytest.mark.parametrize(
    ("given", "xyz"),
    [
        (None, [1, 1, 1]),
        ([], [1, 1, 1]),
        ([1.5], [1.5, 1.5, 1.5]),
        ([2, 1.25], [2, 2, 1.25]),
        ([1, 1.5, 1.25], [1, 1.5, 1.25]),
        ([4, 1, 2, 3], [1, 2, 3]),
    ],
)
def test_the_factor_is_completed_as_the_frame_function_completes_it(given, xyz):
    # [ones(1, 4 - numel(rs)) * rs(1), rs(2 : end)]
    assert dz._rs_xyz(given) == xyz


@pytest.mark.parametrize("factor", [2, (1, 2), (1.25, 1, 1.5), (1, 2, 1)])
def test_large_file_result_equals_the_regular_path_with_its_factor(tmp_path, factor):
    """(MATLAB: largeFileReadsTheResampleFactorAsRegularFiles) one value, (z, xy) and (z, y, x)
    with x and y different: the large file's result is the regular path's (integer factors and
    1 along y give the same values; a fractional y factor only the same size, as in MATLAB)"""
    from petakit5d.deskew_rotate import deskew_rotate_volume
    from petakit5d.deskew_rotate.file import _resample_completed

    # (dz 0.2: no skewed space interpolation, whose uint16 sums depend on a batch's rows; 64
    # rows, MATLAB's test's: each batch's own deskew and rotation centers its rows, which
    # matches the whole result's grid where the batches' and the volume's input rows are even)
    fn, frame = write_input(tmp_path, shape=(29, 64, 61))
    dz.deskew_rotate_zarr(fn, PX, 0.2, save_mip=False, resample_factor=factor, **KW)
    got = read_zarr(str(tmp_path / "DSR" / f"{FSN}.zarr"))
    # (the regular path's factor, completed as XR_deskewRotateFrame completes it)
    want = deskew_rotate_volume(
        frame, SKEW, 0.2, PX, reverse=True, resample_factor=_resample_completed(factor),
        save_16bit=True,
    )  # fmt: skip
    same(got, want)


def test_large_file_is_flipped_as_regular_files(tmp_path):
    """(MATLAB: largeFileIsFlippedAsRegularFiles) each batch is flipped along z after its read:
    the result is the flipped volume's (the blocks never flipped)"""
    from petakit5d.deskew_rotate import deskew_rotate_volume

    fn, frame = write_input(tmp_path)
    dz.deskew_rotate_zarr(fn, PX, 0.2, save_mip=False, flip_zstack=True, **KW)
    want = deskew_rotate_volume(
        np.ascontiguousarray(frame[::-1]), SKEW, 0.2, PX, reverse=True, save_16bit=True
    )
    same(read_zarr(str(tmp_path / "DSR" / f"{FSN}.zarr")), want)


def test_blocks_get_the_flip_and_the_completed_factor(tmp_path, framework_spy):
    fn, _ = write_input(tmp_path)
    dz.deskew_rotate_zarr(fn, PX, DZ, save_mip=False, flip_zstack=True, resample_factor=1.5, **KW)
    tasks = framework_spy.calls[0]["tasks"]
    assert all(t.kwargs["flip_zstack"] is True for t in tasks)
    assert all(list(t.kwargs["resample_factor"]) == [1.5, 1.5, 1.5] for t in tasks)


def test_large_file_parfor_runs_the_blocks_on_a_local_pool(tmp_path, framework_spy):
    """(MATLAB: largeFileParforRunsTheBlocks) parseParfor without a cluster failed on an undefined
    GPUJob; it runs the tasks on a local pool of 12 workers, as XR_crop_zarr"""
    fn, _ = write_input(tmp_path)
    kw = {**KW, "parse_parfor": True}
    dz.deskew_rotate_zarr(fn, PX, DZ, save_mip=False, **kw)
    (call,) = framework_spy.calls
    assert call["cluster_type"] == "parfor" and call["local_worker_num"] == 12
    assert call["gpu_job"] is False and call["task_batch_num"] == 1 and call["uuid"]


def test_masks_follow_the_factor_and_the_flip(tmp_path):
    """(MATLAB: largeFileMIPMaskBoxesFollowTheFactorAndTheFlip) the input box stays the raw
    data's (the data is flipped after its crop); with the flip, the output box and the
    deskewed MIP are those of the flipped masks"""
    raw = tmp_path / "raw"
    flipped = tmp_path / "flipped"
    raw.mkdir()
    flipped.mkdir()
    m = np.zeros(RAW, dtype=np.uint8)
    m[3:20, 10:30, 8:25] = 1
    m[5:27, 12:20, 30:50] = 1
    for folder, mask in ((raw, m), (flipped, m[::-1])):
        synth.write_tiff(folder / "mask_z.tif", mask.max(axis=0))
        synth.write_tiff(folder / "mask_y.tif", np.ascontiguousarray(mask.max(axis=1).T))
        synth.write_tiff(folder / "mask_x.tif", np.ascontiguousarray(mask.max(axis=2).T))
    names = ("mask_z.tif", "mask_y.tif", "mask_x.tif")
    args = ((0, 0, 0), PX, DZ, SKEW, True)
    for rs in ([1, 1, 1], [1.5, 1, 1.25], [1, 1.5, 1.25], [1.5], [1.5, 1]):
        a = dz.deskew_rotate_boxes_from_masks([str(raw / n) for n in names], *args, rs)
        b = dz.deskew_rotate_boxes_from_masks([str(raw / n) for n in names], *args, rs, None, True)
        c = dz.deskew_rotate_boxes_from_masks([str(flipped / n) for n in names], *args, rs)
        assert b[0] == a[0], rs
        assert (c[0][0], c[0][3]) == (RAW[0] - a[0][3], RAW[0] - a[0][0]), rs
        assert b[1] == c[1], rs
        same(b[2], c[2])
        assert b[1] != a[1], rs


@pytest.mark.parametrize("cls", ["uint16", "single"])
def test_block_mips_from_matlabs_layout_equal_the_block_mips(cls):
    from petakit5d.tools.mip import block_mips

    rng = CounterRNG(3, f"dsrzarr/mips/{cls}")
    a = rng.integers(0, 60000, (9, 13, 11), dtype=np.uint16)
    if cls == "single":
        a = a.astype(np.float32)
    box = (3, 17, 5, 12, 30, 16)
    mbs = (4, 8, 6)
    want = block_mips(a, box, mbs, cls)
    got = dz._block_mips_zxy(np.ascontiguousarray(a.transpose(0, 2, 1)), box, mbs, cls)
    for w, g in zip(want, got, strict=True):
        same(g, w)


def test_results_do_not_depend_on_threads(tmp_path, monkeypatch):
    outs = []
    for n in (1, 3):
        monkeypatch.setenv("OMP_NUM_THREADS", str(n))
        d = tmp_path / f"t{n}"
        fn_n, _ = write_input(d)
        dz.deskew_rotate_zarr(fn_n, PX, DZ, **KW)
        outs.append(read_zarr(str(d / "DSR" / f"{FSN}.zarr")))
    same(outs[0], outs[1])


def test_masks_set_the_boxes(tmp_path):
    from petakit5d.deskew_rotate import deskew_rotate_volume

    # (wide in x: the masks' border of 100 voxels leaves part of x out)
    shape = (29, 37, 300)
    mask_tree(tmp_path, shape, 150)
    fn, frame = write_input(tmp_path, shape=shape)
    dz.deskew_rotate_zarr(
        fn,
        PX,
        0.2,
        mask_fullpaths=[str(tmp_path / n) for n in ("mask_z.tif", "mask_y.tif", "mask_x.tif")],
        **KW,
    )
    in_box, out_box, _ = dz.deskew_rotate_boxes_from_masks(
        [str(tmp_path / n) for n in ("mask_z.tif", "mask_y.tif", "mask_x.tif")],
        (100, 100, 100),
        PX,
        0.2,
        SKEW,
        True,
        [1, 1, 1],
    )
    assert in_box[2] == 50 and in_box[5] == 290
    sub = frame[in_box[0] : in_box[3], in_box[1] : in_box[4], in_box[2] : in_box[5]]
    full = deskew_rotate_volume(
        np.ascontiguousarray(sub), SKEW, 0.2, PX, reverse=True, save_16bit=True
    )
    want = full[out_box[0] : out_box[3], :, out_box[2] : out_box[5]]
    same(read_zarr(str(tmp_path / "DSR" / f"{FSN}.zarr")), want)


def test_the_masks_border_follows_the_factor_of_its_axis(tmp_path, framework_spy, monkeypatch):
    """the masks' border is 100 voxels over the factor of its own axis (rs([2, 1, 3]) of the
    [x y z] factor for the [y x z] border), at least 10"""
    fn, _ = write_input(tmp_path)
    seen = []

    def spy(masks, border, *a, **k):
        seen.append(tuple(int(v) for v in border))
        return (0, 0, 0) + RAW, (0, 0, 0, 2, 8, 4), None

    monkeypatch.setattr(dz, "deskew_rotate_boxes_from_masks", spy)
    masks = ["m_z.tif", "m_y.tif", "m_x.tif"]
    dz.deskew_rotate_zarr(
        fn, PX, DZ, save_mip=False, mask_fullpaths=masks, resample_factor=(1.25, 1, 2), **KW
    )
    dz.deskew_rotate_zarr(
        fn, PX, DZ, save_mip=False, mask_fullpaths=masks, resample_factor=(4,),
        result_dir_str="DSR_4/", **KW,
    )  # fmt: skip
    assert seen == [(80, 100, 50), (25, 25, 25)]


def test_the_output_region_of_the_masks_is_written(tmp_path, monkeypatch):
    # the coordinator crops the result to the masks' output region: blocks are computed on
    # that region of the full result (x and z starts moved), the input on the input box
    from petakit5d.deskew_rotate import deskew_rotate_volume

    fn, frame = write_input(tmp_path)
    in_box, out_box = (2, 0, 5, 27, 37, 55), (3, 0, 9, 20, 37, 80)
    monkeypatch.setattr(
        dz, "deskew_rotate_boxes_from_masks", lambda *a, **k: (in_box, out_box, None)
    )
    dz.deskew_rotate_zarr(fn, PX, 0.2, mask_fullpaths=["m_z.tif", "m_y.tif", "m_x.tif"], **KW)
    sub = np.ascontiguousarray(
        frame[in_box[0] : in_box[3], in_box[1] : in_box[4], in_box[2] : in_box[5]]
    )
    full = deskew_rotate_volume(sub, SKEW, 0.2, PX, reverse=True, save_16bit=True)
    want = full[out_box[0] : out_box[3], :, out_box[2] : out_box[5]]
    same(read_zarr(str(tmp_path / "DSR" / f"{FSN}.zarr")), want)
