"""Every task function of the port, run one after another in one process in 5 random orders,
gives what fresh processes give, bit for bit (jobs.isolation's harness; the model is MATLAB's
tests/isolation).

The scenarios build small synthetic datasets and run the port's data wrappers until the
framework call whose tasks they capture, as isolation_prepare.m does: ``run_tasks`` is
replaced by a recorder that lets the earlier calls run and stops the wrapper at the captured
call (the decon frames, which the wrapper runs itself in local mode, are recorded at
``run_record``). Each scenario folder is then as it is just before its tasks run; the folders
are the harness's template, copied into a new work folder before each run.

``test_every_task_function_has_a_scenario`` keeps the list complete (it reads the sources, so it
always runs): a new ``Task.of`` call site in the package needs a scenario, or a listed reason,
here. The harness run (103 tasks of 47 scenarios, about a minute; ``-s`` shows its stages)
runs on Linux, the runners' platform: on macOS the cpp-tiff and cpp-zarr wheels can't be loaded
together in one process yet, which every scenario process does."""

from __future__ import annotations

import ast
import json
import shutil
import sys
import time
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pytest

import petakit5d
from petakit5d.io.tiff import write_tiff
from petakit5d.io.zarr import write_zarr

LINUX = pytest.mark.skipif(
    not sys.platform.startswith("linux"),
    reason="the harness runs where the runners run (and macOS can't load cpp-tiff and cpp-zarr "
    "together yet)",
)
SRC = Path(petakit5d.__file__).resolve().parent
LOCAL = {"parse_cluster": False, "master_compute": True}
FSN = "Scan_Iter_0000_0000_{}_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t"  # noqa: E501
DATASET = ("scan_CamA_ch0_t0000", "scan_CamB_ch0_t0000")

# task functions the scenarios do not run, and why
NOT_RUN = {
    "petakit5d.io.n5:zarr_to_n5_block": "writes blosc-compressed N5 (numcodecs, the optional "
    "n5 extra)",
    "petakit5d.io.n5:n5_to_zarr_block": "reads N5 written with numcodecs (the optional n5 extra)",
    "petakit5d.io.bioformats:bioformats_to_tiff_or_zarr_file": "needs a Bio-Formats reader, "
    "an optional dependency the test environments do not install",
    "petakit5d.tools.psf_analysis:psf_analysis_plot": "draws figures (matplotlib, the "
    "optional figures extra)",
    "petakit5d.tools.psf_analysis:_noop": "does nothing (the task of a file without a PSF)",
    "petakit5d.jobs.files:batch_file_exist": "the rsync wrapper's file check in cluster mode only",
    "petakit5d.stitch.assignment:global_grid_assignment": "made only by the group assignment "
    "in cluster mode (a pure function of its group's .mat file)",
}


class _Stop(BaseException):
    """ends a wrapper at the captured call (a BaseException: the framework's local loop
    retries tasks that raise an Exception)"""


def _frames(folder: Path, names, shape=(12, 40, 36), fmt="tif", seed=0, chunks=(8, 16, 16)):
    folder.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    out = []
    for n in names:
        a = rng.integers(0, 3000, size=shape, dtype=np.uint16)
        if fmt == "zarr":
            write_zarr(a, str(folder / f"{n}.zarr"), block_size=chunks)
        else:
            write_tiff(a, str(folder / f"{n}.tif"))
        out.append(a)
    return str(folder) + "/", out


def capture(build: Callable[[Path], None], folder: Path, stop_at: int = 1) -> list[dict]:
    """the task records of the wrapper's stop_at-th framework call (1-based); earlier calls
    run as they would"""
    from petakit5d.jobs import framework, tasks

    real = framework.run_tasks
    calls = []

    def recorder(inputs, outputs, ts, **kw):
        calls.append(len(ts) if isinstance(ts, (list, tuple)) else 1)
        if len(calls) < stop_at:
            return real(inputs, outputs, ts, **kw)
        lang = kw.get("language", "python")
        raise _Stop([tasks.task_record(t, lang) for t in (ts if isinstance(ts, list) else [ts])])

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(framework, "run_tasks", recorder)
        try:
            build(folder)
        except _Stop as s:
            return s.args[0]
    raise AssertionError(f"{folder.name}: {len(calls)} framework calls, fewer than {stop_at}")


def capture_frames(build: Callable[[Path], None], folder: Path) -> list[dict]:
    """the records the decon and stitching wrappers run themselves in local mode (recorded,
    not run: a frame that waits for another's results would wait for ever here)"""
    from petakit5d.jobs import tasks

    recs = []
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(tasks, "run_record", lambda rec: recs.append(json.loads(json.dumps(rec))))
        build(folder)
    for p in folder.rglob("*.tmp"):
        p.unlink()
    return recs


# --- the scenarios ---------------------------------------------------------------------------
# each builder makes its data in folder and runs a wrapper (local mode)


def b_crop(folder):
    from petakit5d.tools.crop import crop_dataset

    path, _ = _frames(folder, DATASET)
    crop_dataset(path, (2, 4, 6, 10, 30, 33), channel_patterns=("CamA", "CamB"), **LOCAL)


def b_crop_large(folder):
    from petakit5d.tools.crop import crop_dataset

    path, _ = _frames(folder, DATASET[:1], fmt="zarr")
    crop_dataset(
        path, (1, 3, 2, 11, 37, 33), channel_patterns=("CamA",), zarr_file=True,
        large_file=True, batch_size=(8, 16, 16), block_size=(4, 8, 8), **LOCAL,
    )  # fmt: skip


def b_mip(folder):
    from petakit5d.tools.mip import mip_dataset

    path, _ = _frames(folder, DATASET)
    mip_dataset(path, axis=(1, 1, 1), channel_patterns=("CamA", "CamB"), **LOCAL)


def b_mip_zarr(folder):
    from petakit5d.tools.mip import mip_dataset

    path, _ = _frames(folder, DATASET, fmt="zarr")
    mip_dataset(path, axis=(1, 1, 1), channel_patterns=("CamA", "CamB"), zarr_file=True, **LOCAL)


def b_mip_large(folder):
    from petakit5d.tools.mip import mip_dataset

    path, _ = _frames(folder, DATASET[:1], fmt="zarr")
    mip_dataset(
        path, axis=(1, 1, 1), channel_patterns=("CamA",), zarr_file=True, large_file=True,
        batch_size=(4, 16, 16), **LOCAL,
    )  # fmt: skip


def b_resample(folder):
    from petakit5d.tools.resample import resample_dataset

    path, _ = _frames(folder, DATASET)
    resample_dataset(path, (2, 1.5, 1.25), channel_patterns=("CamA", "CamB"), **LOCAL)


def b_resample_large(folder):
    from petakit5d.tools.resample import resample_dataset

    path, _ = _frames(folder, DATASET[:1], fmt="zarr")
    resample_dataset(
        path, (2, 2, 2), channel_patterns=("CamA",), zarr_file=True, large_file=True,
        save_zarr=True, batch_size=(8, 16, 16), block_size=(4, 8, 8), **LOCAL,
    )  # fmt: skip


def _two_channels(folder, fmt="tif"):
    names = [FSN.format(c) for c in ("CamA", "CamB")]
    return _frames(folder, names, shape=(12, 32, 28), fmt=fmt, seed=3)


def b_unmix(folder):
    from petakit5d.tools.unmix import unmix_channels_dataset

    path, _ = _two_channels(folder)
    unmix_channels_dataset(
        path, unmix_factors=(1, -0.1), channel_patterns=("CamA", "CamB"),
        const_background=(100, 100), **LOCAL,
    )  # fmt: skip


def b_unmix_large(folder):
    from petakit5d.tools.unmix import unmix_channels_dataset

    path, _ = _two_channels(folder, "zarr")
    unmix_channels_dataset(
        path, unmix_factors=(1, -0.1), channel_patterns=("CamA", "CamB"), zarr_file=True,
        large_file=True, save_zarr=True, batch_size=(8, 16, 16), block_size=(4, 8, 8), **LOCAL,
    )  # fmt: skip


def b_unmix_gaussian(folder):
    from petakit5d.tools.unmix import unmix_channels_zarr

    path, _ = _two_channels(folder, "zarr")
    fns = [f"{path}{FSN.format(c)}.zarr" for c in ("CamA", "CamB")]
    unmix_channels_zarr(
        fns, (1, -0.1), mode="gaussian", unmix_sigmas=(0, (1, 1, 1)), batch_size=(8, 16, 16),
        block_size=(4, 8, 8), **LOCAL,
    )  # fmt: skip


def b_chromatic(folder):
    from petakit5d.tools.chromatic import chromatic_shift_correction_dataset

    path, _ = _two_channels(folder)
    chromatic_shift_correction_dataset(
        path, chromatic_offset=[[0, 0, 0], [-2, 1, -5]], channel_patterns=("CamA", "CamB"),
        **LOCAL,
    )  # fmt: skip


def b_fft(folder):
    from petakit5d.tools.fft_spectrum import fft_spectrum_dataset

    path, _ = _frames(folder, DATASET)
    fft_spectrum_dataset(path, out_size=(31, 31, 31), channel_patterns=("CamA", "CamB"), **LOCAL)


def b_fsc(folder):
    from petakit5d.testing import synth
    from petakit5d.tools.fsc import fsc_analysis_dataset

    folder.mkdir(parents=True)
    rng = synth.CounterRNG(0, "e06/fsc")
    for k in (2, 4):
        v = 1000.0 + 300.0 * synth.gaussian_filter(rng.spawn(str(k)).normal((40, 40, 40)), 1.0)
        synth.write_tiff(folder / f"Iter_{k:04d}.tif", v.astype(np.float32))
    fsc_analysis_dataset(
        str(folder) + "/", xy_pixel_size=0.1, dz=0.1, dr=2, half_size=(20, 20, 20),
        channel_patterns=("Iter",), channels=(488,), iter_interval=2, **LOCAL,
    )  # fmt: skip


def b_tiff2zarr(folder):
    from petakit5d.io.convert import tiff_to_zarr_dataset

    path, _ = _frames(folder, DATASET)
    tiff_to_zarr_dataset(path, block_size=(4, 16, 16), **LOCAL)


def b_zarr2tiff(folder):
    from petakit5d.io.convert import zarr_to_tiff_dataset

    path, _ = _frames(folder, DATASET, fmt="zarr")
    zarr_to_tiff_dataset(path, channel_patterns=("CamA", "CamB"), **LOCAL)


def b_resave(folder):
    from petakit5d.io.resave import resave_zarr_dataset

    path, _ = _frames(folder, DATASET[:1], fmt="zarr")
    resave_zarr_dataset(path, batch_size=(8, 16, 16), block_size=(4, 8, 8), **LOCAL)


DSR = {"xy_pixel_size": 0.108, "dz": 0.3, "skew_angle": 32.45, "reverse": True}


def b_dsr(folder):
    from petakit5d.deskew_rotate.dataset import deskew_rotate_dataset

    path, _ = _frames(folder, DATASET, shape=(16, 24, 40))
    deskew_rotate_dataset(path, channel_patterns=("CamA", "CamB"), **DSR, **LOCAL)


def b_dsr_large(folder):
    from petakit5d.deskew_rotate.dataset import deskew_rotate_dataset

    path, _ = _frames(folder, DATASET[:1], shape=(16, 24, 40), fmt="zarr")
    deskew_rotate_dataset(
        path, channel_patterns=("CamA",), zarr_file=True, large_file=True,
        batch_size=(1024, 8, 1024), block_size=(64, 8, 64), **DSR, **LOCAL,
    )  # fmt: skip


def b_flatfield(folder):
    from petakit5d.flatfield.estimation import flatfield_background_dataset

    folder.mkdir(parents=True)
    rng = np.random.default_rng(4)
    for cam in ("CamA", "CamB"):
        dark = 100 + rng.normal(0, 3, size=(16, 40, 36))
        write_tiff(np.clip(dark, 0, None).astype(np.uint16), str(folder / f"{FSN.format(cam)}.tif"))
    flatfield_background_dataset(str(folder), channel_patterns=("CamA_ch0", "CamB_ch0"), **LOCAL)


def _decon_data(folder, zarr=False):
    from petakit5d.testing import synth

    data, psf = folder / "data", folder / "psf"
    data.mkdir(parents=True)
    psf.mkdir()
    blur = (2.15, 1.7)
    p = synth.decon_psf((21, 21, 21), (2.0, 1.5))
    write_tiff(np.floor(p + np.float32(0.5)).astype(np.uint16), str(psf / "PSF_CamA.tif"))
    rng = np.random.default_rng(8)
    for t in range(2):
        pos = np.column_stack([rng.uniform(3, 13, 5), rng.uniform(6, 26, 5), rng.uniform(6, 22, 5)])
        im = synth.beads((16, 32, 28), pos, blur, 1000.0, 100.0)
        im = np.floor(im + np.float32(0.5)).astype(np.uint16)
        name = FSN.format("CamA").replace("0000t", f"{t:04d}t")
        if zarr:
            write_zarr(im, str(data / f"{name}.zarr"), block_size=(8, 16, 16))
        else:
            write_tiff(im, str(data / f"{name}.tif"))
    return str(data) + "/", [str(psf / "PSF_CamA.tif")]


DECON = {"xy_pixel_size": 0.108, "dz": 0.2, "dz_psf": 0.2, "background": 100, "decon_iter": 3,
         "channel_patterns": ("CamA_ch0",), "parse_cluster": False, "unit_wait_time": 0,
         "max_trial_num": 1}  # fmt: skip


def _quiet(f, *a, **k):
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f(*a, **k)


def b_decon_rl(folder):
    from petakit5d.decon.dataset import decon_dataset

    data, psfs = _decon_data(folder)
    _quiet(decon_dataset, [data], psf_fullpaths=psfs, rl_method="simplified", **DECON)


def b_decon_omw(folder):
    from petakit5d.decon.dataset import decon_dataset

    data, psfs = _decon_data(folder)
    _quiet(decon_dataset, [data], psf_fullpaths=psfs, rl_method="omw", wiener_alpha=0.005,
           otf_cum_thresh=0.9, skewed=False, **DECON)  # fmt: skip


def b_decon_blocks(folder):
    from petakit5d.decon.dataset import decon_dataset

    data, psfs = _decon_data(folder, zarr=True)
    _quiet(decon_dataset, [data], psf_fullpaths=psfs, rl_method="simplified", zarr_file=True,
           save_zarr=True, large_file=True, large_method="inplace", batch_size=(8, 16, 16),
           block_size=(8, 16, 16), **DECON)  # fmt: skip


STITCH = {"dsr": True, "resample_type": "isotropic", "xy_pixel_size": 0.108, "axis_order": "x,y,z",
          "channel_patterns": ("CamA_ch0",), "primary_ch": "CamA_ch0", "parse_cluster": False,
          "max_trial_num": 1}  # fmt: skip


def _stitch_data(folder):
    from petakit5d.testing.synth import cut_tiles, write_stitch_tiles

    data = folder / "data"
    data.mkdir(parents=True)
    rng = np.random.default_rng(9)
    gt = rng.integers(1000, 60000, size=(16, 60, 72), dtype=np.uint16)
    starts = [(0, 0, 0), (0, 28, 0), (0, 0, 36), (0, 28, 36)]
    csv = None
    for it in range(2):
        csv = write_stitch_tiles(str(data), cut_tiles(gt, starts, (16, 32, 36)), iter=it)[1]
    return str(data) + "/", str(csv)


def b_stitch(folder):
    from petakit5d.stitch.dataset import stitch_dataset

    data, csv = _stitch_data(folder)
    _quiet(stitch_dataset, data, csv, xcorr_shift=True, xcorr_mode="primaryFirst",
           block_size=(16, 16, 16), batch_size=(16, 32, 32), **STITCH)  # fmt: skip


def b_stitch_frames(folder):
    """(each time point registered by itself: no frame waits for the first's results)"""
    from petakit5d.stitch.dataset import stitch_dataset

    data, csv = _stitch_data(folder)
    _quiet(stitch_dataset, data, csv, xcorr_shift=True, xcorr_mode="all",
           block_size=(16, 16, 16), batch_size=(16, 32, 32), **STITCH)  # fmt: skip


def b_stitch_separate(folder):
    """('separate': each tile on a canvas of its own, cropped from the stitched volume's place)"""
    from petakit5d.stitch.dataset import stitch_dataset

    data, csv = _stitch_data(folder)
    _quiet(stitch_dataset, data, csv, xcorr_shift=False, blend_method="separate",
           block_size=(16, 16, 16), batch_size=(16, 32, 32), **STITCH)  # fmt: skip


def b_stitch_large(folder):
    from petakit5d.stitch.dataset import stitch_dataset

    data, csv = _stitch_data(folder)
    _quiet(stitch_dataset, data, csv, xcorr_shift=True, xcorr_mode="primaryFirst",
           block_size=(16, 16, 16), batch_size=(16, 32, 32), large_file=True,
           pool_size=(2, 2, 2), **STITCH)  # fmt: skip


def b_stitch_zarr_box(folder):
    from petakit5d.stitch.dataset import stitch_dataset
    from petakit5d.testing.synth import cut_tiles, write_stitch_tiles

    data = folder / "data"
    data.mkdir(parents=True)
    rng = np.random.default_rng(10)
    gt = rng.integers(1000, 60000, size=(16, 60, 72), dtype=np.uint16)
    starts = [(0, 0, 0), (0, 28, 0), (0, 0, 36), (0, 28, 36)]
    csv = write_stitch_tiles(str(data), cut_tiles(gt, starts, (16, 32, 36)), save_format="zarr")[1]
    _quiet(stitch_dataset, str(data) + "/", str(csv), xcorr_shift=False, zarr_file=True,
           input_bbox=(0, 1, 1, 16, 31, 35), block_size=(16, 16, 16), batch_size=(16, 32, 32),
           **STITCH)  # fmt: skip


def r_stitch_wrapper(folder):
    """the pipeline's stitching task: the whole wrapper as one task"""
    from petakit5d.jobs.tasks import Task, task_record
    from petakit5d.stitch.dataset import stitch_dataset

    data, csv = _stitch_data(folder)
    kw = dict(STITCH, xcorr_shift=False, block_size=[16, 16, 16], batch_size=[16, 32, 32])
    return [task_record(Task.of(stitch_dataset, data, csv, **kw), "python")]


def _puncta(folder):
    from petakit5d.testing import synth

    folder.mkdir(parents=True)
    puncta = np.array([(10, 12, 14), (20, 30, 26)], dtype=np.float64)
    v = synth.beads((32, 40, 36), puncta, 1.5, 20000.0, 0.0).astype(np.float64) + 100
    v[4:28, 6:34, 4:30] += 4900 + 200 * np.random.default_rng(1).normal(size=(24, 28, 26))
    fn = str(folder / "expanded.zarr")
    write_zarr(v.astype(np.float32), fn, block_size=(16, 16, 16))
    return fn


PUNCTA = {"batch_size": (32, 32, 32), "block_size": (16, 16, 16), "parse_cluster": False}


def b_puncta(folder):
    from petakit5d.tools.puncta_removal import puncta_removal_zarr

    puncta_removal_zarr(_puncta(folder), **PUNCTA)


def b_puncta_bg(folder):
    from petakit5d.tools.puncta_removal import puncta_removal_zarr

    puncta_removal_zarr(_puncta(folder), mode="bg_estimation", **PUNCTA)


def b_puncta_interp(folder):
    from petakit5d.tools.puncta_removal import puncta_removal_zarr

    fn = _puncta(folder)
    puncta_removal_zarr(fn, mode="bg_estimation", **PUNCTA)
    bg = folder / "Puncta_background_estimation"
    (bg_fn,) = [str(p) for p in bg.glob("*_background_values.zarr")]
    puncta_removal_zarr(fn, mode="bg_interpolated", bg_fullpath=bg_fn, **PUNCTA)


def b_psf_detection(folder):
    from petakit5d.testing import synth
    from petakit5d.tools.psf_analysis import psf_detection_and_analysis_dataset

    folder.mkdir(parents=True)
    centers = np.array([[20.0, 24.0, 24.0], [44.0, 70.0, 64.0]])
    im = synth.beads((61, 96, 88), centers, (3.0, 1.6), 2000.0, 0.0)
    write_tiff(np.floor(im + np.float32(0.5)).astype(np.uint16), str(folder / "beads.tif"))
    psf_detection_and_analysis_dataset(
        str(folder) + "/", crop_size=(21, 32, 32), dist_thresh=(21, 32, 32),
        channel_patterns=("beads",), channels=(488,), parse_cluster=False,
    )  # fmt: skip


def b_mip_masks(folder):
    from petakit5d.tools.mip import mip_dataset
    from petakit5d.tools.mip_masks import generate_mip_masks

    path, _ = _frames(folder, DATASET, shape=(12, 40, 36), seed=11)
    mip_dataset(path, axis=(1, 1, 1), channel_patterns=("CamA", "CamB"), **LOCAL)
    generate_mip_masks(f"{path}MIPs/", channel_patterns=("CamA", "CamB"), int_thresh=2900,
                       vol_thresh=20, dilate_size=2, **LOCAL)  # fmt: skip


def _sim_frames(folder, psf=True):
    """a synthetic 5-phase lattice SIM file (and the phase-stepped bead for its OTF)"""
    from petakit5d.sim import synthetic

    folder.mkdir(parents=True, exist_ok=True)
    beads = [(2.5, 8.2, 7.1, 900.0), (4.0, 17.3, 21.6, 1400.0), (6.2, 25.4, 9.8, 1100.0)]
    write_tiff(synthetic.phase_stack((10, 36, 32), beads), str(folder / "cell_CamA.tif"))
    if psf:
        write_tiff(synthetic.phase_psf((11, 32, 30)), str(folder / "psf.tif"))
    return str(folder) + "/"


def b_sim_recon(folder):
    from petakit5d.sim.dataset import sim_recon_dataset

    path = _sim_frames(folder)
    sim_recon_dataset(path, [path + "psf.tif"], channel_patterns=("CamA",), use_gpu=False, **LOCAL)


def b_sim_deskew(folder):
    from petakit5d.sim.dataset import deskew_phases_dataset

    path = _sim_frames(folder, psf=False)
    deskew_phases_dataset(path, 0.108, 0.3, channel_patterns=("CamA",), **LOCAL)


def b_sim_rotate(folder):
    from petakit5d.sim.dataset import rotate_3d_dataset

    path = _sim_frames(folder, psf=False)
    rotate_3d_dataset(path, 0.108, 0.3, channel_patterns=("CamA",), **LOCAL)


def b_detection_frames(folder):
    """the frame-parallel detection of a two-channel movie of two DS frames (60 spots on
    noise), 'lowSNR' without the background correction; the settings file gives the
    threshold"""
    import os

    from petakit5d.detection import run
    from petakit5d.testing import synth

    ex = folder / "Ex1_synth_z0p5"
    cols = []
    shape = (24, 72, 68)
    for c in (1, 2):
        ds = ex / f"ch{c}" / "DS"
        ds.mkdir(parents=True)
        paths = []
        for k in range(2):
            rng = np.random.default_rng(10 * c + k)
            centers = rng.uniform(0.0, 1.0, (60, 3)) * (np.array(shape) - 8.0) + 4.0
            amps = rng.uniform(250.0, 600.0, 60)
            im = synth.beads(shape, centers, (1.5, 1.3), amps, 100.0)
            noise = rng.normal(0.0, 4.0, im.shape)
            p = ds / f"synth_ch{c}_stack{k:04d}.tif"
            write_tiff(np.clip(np.round(im + noise), 1, 65535).astype(np.uint16), str(p))
            paths.append(str(p))
        col = np.empty((len(paths), 1), dtype=object)
        col[:, 0] = paths
        cols.append(col)
    with open(ex / "synth_Settings.txt", "w", newline="\n") as f:
        f.write('Camera type = "Orca4.0"\n')
    chans = [str(ex / f"ch{c}") + os.sep for c in (1, 2)]
    movie = {"source": chans[0], "channels": chans, "framePathsDS": cols, "zAniso": 2.5}
    run.run_detection_frame_parallel(
        movie, [[1.3, 1.5], [1.3, 1.5]], [chans[0] + "Analysis" + os.sep],
        detection_method="lowSNR",
    )  # fmt: skip


# name, builder, the captured framework call (1-based; 0: the frames the wrapper runs itself
# in local mode; -1: the builder returns the records), the task functions it runs
_C, _M, _R, _U = "tools.crop:", "tools.mip:", "tools.resample:", "tools.unmix:"
_S, _P = "stitch.", "tools.puncta_removal:puncta_removal_"
SCENARIOS: list[tuple[str, Callable[[Path], None], int, str]] = [
    ("crop", b_crop, 1, _C + "crop_file"),
    ("crop_large", b_crop_large, 1, _C + "crop_file"),
    ("crop_blocks", b_crop_large, 2, _C + "crop_block"),
    ("mip", b_mip, 1, _M + "save_mip_tiff"),
    ("mip_zarr", b_mip_zarr, 1, _M + "save_mip_zarr"),
    ("mip_large", b_mip_large, 1, _M + "mip_zarr"),
    ("mip_blocks", b_mip_large, 2, _M + "mip_block"),
    ("resample", b_resample, 1, _R + "resample_file"),
    ("resample_large", b_resample_large, 1, _R + "resample_zarr"),
    ("resample_blocks", b_resample_large, 2, _R + "resample_zarr_block"),
    ("unmix", b_unmix, 1, _U + "unmix_channels_file"),
    ("unmix_large", b_unmix_large, 1, _U + "unmix_channels_zarr"),
    ("unmix_blocks", b_unmix_large, 2, _U + "unmix_channels_block"),
    ("unmix_gaussian", b_unmix_gaussian, 1, _U + "unmix_channels_gaussian_block"),
    ("chromatic", b_chromatic, 1, _C + "crop_file"),
    ("fft", b_fft, 1, "tools.fft_spectrum:fft_spectrum_file"),
    ("fsc", b_fsc, 1, "tools.fsc:fsc_analysis_file"),
    ("tiff2zarr", b_tiff2zarr, 1, "io.convert:tiff_to_zarr"),
    ("zarr2tiff", b_zarr2tiff, 1, "io.convert:zarr_to_tiff"),
    ("resave", b_resave, 1, "io.resave:resave_zarr"),
    ("resave_blocks", b_resave, 2, "io.resave:resave_zarr_block"),
    ("dsr", b_dsr, 1, "deskew_rotate.file:deskew_rotate_file"),
    ("dsr_large", b_dsr_large, 1, "deskew_rotate.zarr:deskew_rotate_zarr"),
    ("dsr_blocks", b_dsr_large, 2, "deskew_rotate.zarr:deskew_rotate_block"),
    ("flatfield", b_flatfield, 1, "flatfield.estimation:flatfield_background_estimation"),
    ("decon_rl", b_decon_rl, 0, "decon.file:decon_file"),
    ("decon_omw", b_decon_omw, 0, "decon.file:decon_file"),
    ("decon_blocks", b_decon_blocks, 1, "decon.block:decon_zarr_block"),
    ("stitch_frames", b_stitch_frames, 0, _S + "volume:stitch_volume"),
    ("stitch_tiles", b_stitch, 1, "io.convert:tiff_to_zarr"),
    ("stitch_zarr_tiles", b_stitch_zarr_box, 1, _S + "tiles:process_zarr_tile"),
    ("stitch_xcorr", b_stitch, 2, _S + "registration:cross_correlation_registration"),
    ("stitch_blockinfo", b_stitch, 3, _S + "blockinfo:organize_block_info"),
    ("stitch_bwdist", b_stitch, 4, _S + "distance:compute_tile_bwdist"),
    ("stitch_blocks", b_stitch, 5, _S + "block:process_stitch_block"),
    ("stitch_separate", b_stitch_separate, 3, _C + "crop_zarr"),
    ("stitch_large_mips", b_stitch_large, 2, _M + "mip_zarr"),
    ("stitch_large_bwdist", b_stitch_large, 12, _S + "distance:compute_tile_bwdist_mip_slabs"),
    ("stitch_wrapper", r_stitch_wrapper, -1, _S + "dataset:stitch_dataset"),
    ("puncta", b_puncta, 1, _P + "updated_block"),
    ("puncta_bg", b_puncta_bg, 1, _P + "background_estimation_block"),
    ("puncta_interp", b_puncta_interp, 2, _P + "interpolated_background_block"),
    ("psf_detection", b_psf_detection, 1, "tools.psf_analysis:psf_detection_and_cropping"),
    ("mip_masks", b_mip_masks, 2, "tools.mip_masks:generate_single_mip_mask"),
    ("sim_recon", b_sim_recon, 1, "sim.file:sim_recon_file"),
    ("sim_deskew", b_sim_deskew, 1, "sim.file:deskew_phases_file"),
    ("sim_rotate", b_sim_rotate, 1, "deskew_rotate.file:rotate_file"),
    ("detection_frames", b_detection_frames, 1, "detection.run:frame_task"),
]


def build_template(root: Path, names=None) -> dict:
    """builds the scenarios in root (the harness's <work>/out); the spec's tasks and, per
    scenario, its task functions or the error that stopped its builder"""
    tasks: list[dict] = []
    funcs: dict[str, list[str]] = {}
    errors: dict[str, str] = {}
    for name, build, call, _ in SCENARIOS:
        if names and name not in names:
            continue
        t0 = time.perf_counter()
        folder = root / name
        try:
            if call == 0:
                recs = capture_frames(build, folder)
            elif call < 0:
                recs = build(folder)
            else:
                recs = capture(build, folder, call)
        except Exception as e:  # (reported by the test, with the others)
            errors[name] = f"{type(e).__name__}: {e}"
            continue
        if not recs:
            errors[name] = "no tasks"
        tasks += recs
        funcs[name] = sorted({r["func"] for r in recs})
        print(f"{name}: {len(recs)} tasks in {time.perf_counter() - t0:.1f} s", flush=True)
    return {"tasks": tasks, "funcs": funcs, "errors": errors}


def task_call_sites() -> set[str]:
    """the names the package's Task.of(...) and Task("module:function") calls give (a Name or
    an attribute's last part, or the string)"""
    names = set()
    for p in sorted(SRC.rglob("*.py")):
        for node in ast.walk(ast.parse(p.read_text(encoding="utf-8"))):
            if not isinstance(node, ast.Call) or not node.args:
                continue
            f, a = node.func, node.args[0]
            if isinstance(f, ast.Attribute) and f.attr == "of" and _name(f.value) == "Task":
                names.add(a.id if isinstance(a, ast.Name) else getattr(a, "attr", ast.unparse(a)))
            elif isinstance(f, ast.Name) and f.id == "Task":
                if isinstance(a, ast.Constant):
                    names.add(str(a.value).rpartition(":")[2])
                elif isinstance(a, ast.Name):
                    names.add(a.id)
    return names


def _name(node: ast.AST) -> str:
    return node.id if isinstance(node, ast.Name) else ""


# call sites whose first argument is a variable, a constant or a table: what it names
INDIRECT = {
    "func": {"save_mip_tiff", "save_mip_zarr", "n5_to_zarr_block", "zarr_to_n5_block"},
    "FRAME_FUNC": {"decon_file", "stitch_volume"},
    "REGISTRATION_TASK": {"cross_correlation_registration"},
    "TASK_FUNCS[spec.matlab_name]": {"deskew_rotate_file", "stitch_dataset"},
    "task['func']": set(),
    "task": set(),
}


def test_every_task_function_has_a_scenario():
    declared = {f"petakit5d.{f}" for *_, f in SCENARIOS} | set(NOT_RUN)
    short = {f.rpartition(":")[2] for f in declared}
    sites = set()
    for n in task_call_sites():
        sites |= INDIRECT.get(n, {n})
    assert sites <= short, f"task functions without a scenario: {sorted(sites - short)}"
    assert all(f.split(":")[0] for f in declared)


@LINUX
@pytest.mark.slow
def test_every_task_function_equals_fresh_processes_in_5_orders(tmp_path):
    from petakit5d.jobs.isolation import run_harness

    work = tmp_path / "work"
    out = work / "out"
    out.mkdir(parents=True)
    t0 = time.perf_counter()
    built = build_template(out)
    assert not built["errors"], f"scenarios that failed to build: {built['errors']}"
    for name, _, _, func in SCENARIOS:
        assert built["funcs"][name] == [f"petakit5d.{func}"], (name, built["funcs"][name])
    keep = tmp_path / "template"
    shutil.copytree(out, keep)
    shutil.rmtree(work)
    print(f"template: {len(built['tasks'])} tasks in {time.perf_counter() - t0:.0f} s", flush=True)
    prep = {"func": "shutil:copytree", "args": [str(keep), str(work / "out")],
            "kwargs": {"dirs_exist_ok": True}}  # fmt: skip
    spec = {"prepare": prep, "tasks": built["tasks"]}
    res = run_harness(spec, str(work), str(tmp_path / "results"), orders=5, seed=6, jobs=4)
    print(f"harness: {time.perf_counter() - t0:.0f} s", flush=True)
    for r in res.runs:
        print(f"{r['name']}: {r['seconds']:.0f} s, equal {r['equal']}", flush=True)
    assert len(res.runs) == 5
    assert res.ok, res.text()
