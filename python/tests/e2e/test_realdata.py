"""Real-data equivalence of the production workflow, in local and cluster mode, on CPUs and GPUs.

The workflow is the shape of a production run on the demo cell dataset (16 skewed light-sheet
tiles: 4 tiles x 2 cameras x 2 time points), each step one MATLAB wrapper called by its MATLAB
name (``petakit5d.matlab``) with the values a MATLAB script passes (``workflow``):

1. the flat-field function files (``processFFCorrectionFrame`` with the dataset's flat fields
   and camera backgrounds), as the production scripts write them;
2. ``XR_tiffToZarr_wrapper`` with the flat-field correction (``zarr_ff/``);
3. ``XR_matlab_stitching_wrapper`` in skewed space from the raw tiffs, the flat field applied
   while converting the tiles, feather blending, the primary channel's tiles registered by
   cross-correlation (``xcorrShift``; the port's ``stitch.registration``) and the other channel
   placed by its result (``matlab_stitch/``);
4. ``XR_decon_data_wrapper``: large-file OMW decon of the stitched volumes in place, block jobs
   on CPUs (``matlab_decon_omw/``) or GPUs (``matlab_decon_omw_<gpu>/``);
5. ``XR_deskew_rotate_data_wrapper``: large-file deskew/rotation of the deconvolved volumes;
6. ``XR_MIP_wrapper``: the MIPs of the stitched volumes (zarr) and of the deskewed volumes
   (large zarr), along all three axes.

MATLAB runs the same steps from the same description (``write_spec`` writes it as JSON for a
MATLAB driver); ``compare_runs`` compares two runs' trees (names, voxels, metadata, .mat and
JSON values, with each side's root masked) with the listed exceptions (``EXCEPTIONS``,
docs/parity/realdata.md).

Environment (all optional; the tests skip without them):

``PETAKIT5D_DEMO_DATA``
    the demo cell dataset (the folder or its parent); only read, the runs link its files.
``PETAKIT5D_REALDATA_MATLAB``
    a MATLAB run of this workflow: the folder that holds its dataset folder.
``PETAKIT5D_TEST_SLURM_CONFIG``
    a cluster configuration for the cluster-mode test (``cluster``).
"""

from __future__ import annotations

import dataclasses
import glob
import json
import os
import re
import shutil
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pytest

CELL = "PetaKit5D_demo_cell_image_dataset"
CHANNELS = ("CamA_ch0", "CamB_ch0")
IMAGE_LIST = "ImageList_from_encoder.csv"
STITCH = "matlab_stitch"
DECON = "matlab_decon_omw"
ZARR_FF = "zarr_ff"
MIPS = "MIPs_wrapper"
FUNC_DIR = "tmp/processFunc"
PSFS = ("PSF/560_c.tif", "PSF/488_2_c.tif")
STEP_NAMES = (
    "function_files",
    "tiff_to_zarr",
    "stitch",
    "decon",
    "deskew_rotate",
    "mip_stitch",
    "mip_dsr",
)


# ---------------------------------------------------------------------------------------------
# the workflow
# ---------------------------------------------------------------------------------------------


@dataclasses.dataclass
class Step:
    """one MATLAB wrapper call: ``fun(*args, name, value, ...)``; ``columns``: cell arguments
    that MATLAB takes as columns (one row per channel), others are rows"""

    name: str
    fun: str
    args: list[Any]
    kw: dict[str, Any]
    columns: tuple[str, ...] = ()


def _slash(p: str | os.PathLike) -> str:
    p = os.fspath(p).replace("\\", "/")
    return p if p.endswith("/") else p + "/"


def flat_field_files(data: str) -> list[tuple[str, str]]:
    """(flat-field image, camera background) of each channel, as the dataset's ``FF`` folder
    holds them: ``FF/averaged/ff_<cam>_ch0_*.tif`` and ``FF/<folder>/AVG_DF400_<cam>_10ms.tif``"""
    data = _slash(data)
    out = []
    for ch in CHANNELS:
        cam = ch.split("_")[0]
        ff = sorted(glob.glob(f"{data}FF/averaged/ff_{ch}_*.tif"))
        bg = sorted(glob.glob(f"{data}FF/*/AVG_DF400_{cam}_10ms.tif"))
        if not ff or not bg:
            raise FileNotFoundError(f"no flat field or background of {cam} in {data}FF")
        out.append(
            (
                f"{data}FF/averaged/{os.path.basename(ff[0])}",
                f"{data}FF/{os.path.basename(os.path.dirname(bg[0]))}/{os.path.basename(bg[0])}",
            )
        )
    return out


def function_files(data: str) -> list[tuple[str, str]]:
    """(file, usrFun) of the flat-field function of each channel, as a production script
    writes them: ``sprintf("@(x)processFFCorrectionFrame(x,'%s','%s',%d,%d,%s)", ...)`` with a
    constant offset of 1, a lower limit of 0.4 and no flat-field background"""
    from petakit5d.core.text import sprintf

    data = _slash(data)
    out = []
    for i, (ff, bg) in enumerate(flat_field_files(data)):
        usr_fun = sprintf(
            "@(x)processFFCorrectionFrame(x,'%s','%s',%d,%d,%s)", ff, bg, 1, 0.4, "false"
        )
        out.append((f"{data}{FUNC_DIR}/processFunction_c{i + 1:04d}.mat", usr_fun))
    return out


def write_function_files(data: str) -> list[str]:
    """writes the flat-field function files (``usrFun``, char, v7.3) and returns their paths"""
    from petakit5d.io import write_mat

    paths = []
    for fn, usr_fun in function_files(data):
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        write_mat(fn, {"usrFun": usr_fun}, version="7.3")
        paths.append(fn)
    return paths


def workflow(
    data: str,
    *,
    job_log_dir: str,
    parse_cluster: bool = True,
    config_file: str = "",
    gpu_config_file: str = "",
    gpu_name: str = "",
    cpus_per_task: int = 12,
) -> list[Step]:
    """the workflow's wrapper calls on the dataset folder ``data`` with MATLAB's names and
    values; with ``gpu_name`` the decon runs its blocks on GPUs (``gpu_config_file``) into
    ``matlab_decon_omw_<gpu_name>``"""
    data = _slash(data)
    stitch = f"{data}{STITCH}/"
    decon_dir = f"{DECON}_{gpu_name}" if gpu_name else DECON
    funcs = [fn for fn, _ in function_files(data)]
    jobs = {
        "parseCluster": parse_cluster,
        "masterCompute": True,
        "jobLogDir": job_log_dir,
        "configFile": config_file,
    }
    steps = [
        Step(
            "tiff_to_zarr",
            "XR_tiffToZarr_wrapper",
            [data],
            {
                "resultDirName": ZARR_FF,
                "channelPatterns": list(CHANNELS),
                "processFunPath": funcs,
                **jobs,
            },
            columns=("processFunPath",),
        ),
        Step(
            "stitch",
            "XR_matlab_stitching_wrapper",
            [data, data + IMAGE_LIST],
            {
                "DS": False,
                "DSR": False,
                "SkewAngle": 32.45,
                "xyPixelSize": 0.108,
                "dz": 0.3,
                "Reverse": True,
                "axisOrder": "-x,y,z",
                "dataOrder": "y,x,z",
                "xcorrShift": True,
                "ChannelPatterns": list(CHANNELS),
                "PrimaryCh": "CamB_ch0",
                "zarrFile": False,
                "batchSize": [512, 512, 512],
                "blockSize": [256, 256, 256],
                "TileOffset": 0,
                "resampleType": "isotropic",
                "resampleFactor": [],
                "BlendMethod": "feather",
                "resultDir": STITCH,
                "xyMaxOffset": 150,
                "zMaxOffset": 150,
                "xcorrDownsample": [2, 2, 2],
                "InputBbox": [],
                "tileOutBbox": [],
                "parseSettingFile": False,
                "processFunPath": funcs,
                "Save16bit": True,
                **jobs,
                "mccMode": False,
            },
            columns=("processFunPath",),
        ),
        Step(
            "decon",
            "XR_decon_data_wrapper",
            [[stitch]],
            {
                "resultDirName": decon_dir,
                "xyPixelSize": 0.108,
                "dz": 0.3,
                "ChannelPatterns": list(CHANNELS),
                "PSFFullpaths": [data + p for p in PSFS],
                "dzPSF": 0.3,
                "parseSettingFile": False,
                "RLmethod": "omw",
                "wienerAlpha": [0.05, 0.005],
                "OTFCumThresh": [0.88, 0.9],
                "skewed": True,
                "Background": 100,
                "DeconIter": 2,
                "EdgeErosion": 20,
                "Save16bit": True,
                "zarrFile": True,
                "batchSize": [512, 768, 768],
                "blockSize": [256, 256, 256],
                "parseCluster": parse_cluster,
                "largeFile": True,
                "largeMethod": "inplace",
                "GPUJob": bool(gpu_name),
                "debug": False,
                "cpusPerTask": cpus_per_task,
                "jobLogDir": job_log_dir,
                "ConfigFile": config_file,
                "GPUConfigFile": gpu_config_file,
                "mccMode": False,
            },
        ),
        Step(
            "deskew_rotate",
            "XR_deskew_rotate_data_wrapper",
            [[f"{stitch}{decon_dir}/"]],
            {
                "xyPixelSize": 0.108,
                "dz": 0.3,
                "Reverse": True,
                "ChannelPatterns": list(CHANNELS),
                "largeFile": True,
                "zarrFile": True,
                "saveZarr": True,
                "BatchSize": [512, 512, 512],
                "blockSize": [256, 256, 256],
                "Save16bit": True,
                **jobs,
                "mccMode": False,
            },
        ),
        Step(
            "mip_stitch",
            "XR_MIP_wrapper",
            [[stitch]],
            {
                "resultDirName": MIPS,
                "axis": [1, 1, 1],
                "channelPatterns": list(CHANNELS),
                "zarrFile": True,
                "largeFile": False,
                "Save16bit": True,
                **jobs,
                "mccMode": False,
            },
        ),
        Step(
            "mip_dsr",
            "XR_MIP_wrapper",
            [[f"{stitch}{decon_dir}/DSR/"]],
            {
                "resultDirName": MIPS,
                "axis": [1, 1, 1],
                "channelPatterns": list(CHANNELS),
                "zarrFile": True,
                "largeFile": True,
                "batchSize": [2048, 2048, 2048],
                "Save16bit": True,
                **jobs,
                "mccMode": False,
            },
        ),
    ]
    return steps


def select(steps: Sequence[Step], names: Iterable[str] | None) -> list[Step]:
    """the steps named (all if None), in the workflow's order"""
    if names is None:
        return list(steps)
    names = set(names)
    unknown = names - set(STEP_NAMES)
    if unknown:
        raise ValueError(f"unknown steps {sorted(unknown)}; the steps are {STEP_NAMES}")
    return [s for s in steps if s.name in names]


def run_python(data: str, steps: Sequence[Step], *, function_files_too: bool = True) -> dict:
    """runs the steps through the MATLAB-name layer; returns each step's wall time (s)"""
    import time

    from petakit5d import matlab

    times = {}
    if function_files_too:
        t0 = time.perf_counter()
        write_function_files(data)
        times["function_files"] = time.perf_counter() - t0
    for s in steps:
        t0 = time.perf_counter()
        nv = [x for kv in s.kw.items() for x in kv]
        getattr(matlab, s.fun)(*s.args, *nv)
        times[s.name] = time.perf_counter() - t0
    return times


def write_spec(path: str, data: str, steps: Sequence[Step]) -> None:
    """the steps (and the function files) as JSON for a MATLAB driver: ``files`` [{path,
    usrFun}], ``steps`` [{name, fun, args, kw, columns}]; numeric lists are rows, cells are rows
    except the ``columns``"""
    spec = {
        "files": [{"path": p, "usrFun": f} for p, f in function_files(data)],
        "steps": [dataclasses.asdict(s) for s in steps],
    }
    with open(path, "w", newline="\n") as f:
        json.dump(spec, f, indent=1)


def link_dataset(src: str, dst: str) -> str:
    """a run folder for the dataset ``src``: its tiffs linked, its other files copied (the
    wrappers write into the dataset folder); returns the folder"""
    src, dst = _slash(src), _slash(dst)
    for dp, _dns, fns in os.walk(src):
        rel = os.path.relpath(dp, src)
        out = os.path.normpath(os.path.join(dst, rel))
        os.makedirs(out, exist_ok=True)
        for fn in fns:
            p, q = os.path.join(dp, fn), os.path.join(out, fn)
            if os.path.lexists(q):
                continue
            if fn.lower().endswith(".tif") and os.path.getsize(p) > 1 << 24:
                os.symlink(os.path.realpath(p), q)
            else:
                shutil.copyfile(p, q)
    return dst


# ---------------------------------------------------------------------------------------------
# comparison of two runs
# ---------------------------------------------------------------------------------------------

# outputs that are not compared: (regex on the name relative to the dataset folder, why)
NOT_COMPARED = (
    (r"(.*/)?job_logs(/|$)", "job logs (job ids, hosts and times)"),
    (
        r".*\.(png|fig)$",
        "figures: MATLAB draws the back projector's figure when it runs with a display "
        "toolkit; the port draws it with the optional figures extra (demos.md row 5)",
    ),
)


# the GPU decon results (compared by themselves, with ``only``)
GPU_DECON = rf"^{STITCH}/{DECON}_[^/]+/"

# differences between local and cluster mode that MATLAB's modes have too: the stitching's block
# information is made per blending task (all batches in one task locally, 20 per cluster job),
# and the parameter files hold the call's own parseCluster and configFile
MODE_NAMES = (rf"{STITCH}/stitchInfo/.*_task_size_\d+",)
MODE_FIELDS = ("parseCluster", "configFile", "PerBlockInfoFullpaths")


@dataclasses.dataclass(frozen=True)
class Exception_:
    """a listed difference: outputs (regex) allowed to differ within ``spec``"""

    row: int
    pattern: str
    spec: Any
    why: str


def exceptions(kind: str) -> list[Exception_]:
    """the listed differences of docs/parity/realdata.md for a comparison ``kind``:
    'exact' (both CPU FFT backends exact), 'permissive' (the port on its default CPU FFT),
    'gpu' (GPU decon on the same GPU model), 'mode' (the port's cluster against local mode)"""
    from petakit5d.testing.compare import tol

    out: list[Exception_] = []
    if kind == "permissive":
        why = (
            "record 0006: the default CPU FFT (pocketfft) rounds differently from MATLAB's "
            "FFTW; the OMW decon amplifies it"
        )
        spec = tol(reason="realdata row 1 (record 0006)", abs=1, max_frac=1e-3, rel_l2=1e-4)
        # the OTFs and back projectors (single) made by FFTs: record 0006's single-transform rule
        otf = tol(reason="realdata row 1 (record 0006)", rel_max=1e-5, rel_l2=1e-6)
        for pat in (rf"{STITCH}/{DECON}/", rf"{STITCH}/{DECON}_[^/]+/"):
            out.append(Exception_(1, pat + r"[^/]*\.zarr$", spec, why))
            out.append(Exception_(1, pat + r"MIPs/", spec, why))
            out.append(Exception_(1, pat + r"DSR/", spec, why))
            out.append(Exception_(1, pat + r"psfgen/.*_(OTF|back_projector_.*)\.tif$", otf, why))
        # the registration's correlation peaks (the shifts are whole numbers, so they stay exact)
        out.append(
            Exception_(
                3,
                rf"{STITCH}/xcorr/[^/]+/xcorr_tile_[^/]*\.mat$",
                tol(reason="realdata row 3 (P7)", abs=5e-3),
                "record 0006 (P7): the default CPU FFT moves the correlation peaks in their last "
                "bits; the shifts, and so the stitched volumes, are the same",
            )
        )
    return out


def outputs(root: str, inputs: Iterable[str] = ()) -> list[str]:
    """the run's outputs under the dataset folder ``root``: files and zarr arrays that are not
    inputs (empty folders are left out)"""
    from petakit5d.testing.compare import list_tree

    skip = set(inputs)
    rels = [r for r in list_tree(root) if r not in skip]
    return [r for r in rels if r.endswith(".zarr") or not os.path.isdir(os.path.join(root, r))]


def not_compared(rel: str) -> str:
    for pat, why in NOT_COMPARED:
        if re.match(pat, rel):
            return why
    return ""


def compare_runs(
    ref_data: str,
    act_data: str,
    *,
    inputs: Iterable[str] = (),
    kind: str = "exact",
    ref_roots: Sequence[str] = (),
    act_roots: Sequence[str] = (),
    only: str = "",
) -> tuple[Any, dict]:
    """compares the outputs of two runs' dataset folders; returns (report, summary) with the
    differences within listed exceptions as notes. ``*_roots``: each side's own folders
    (masked as <WORK>, <CODE>, ... in names, text and values)"""
    from petakit5d.testing.compare import Masks, Report, compare_trees, scan_tree, tree

    inputs = list(inputs)
    skip = [GPU_DECON] if not only else []
    if kind == "mode":
        skip += MODE_NAMES

    def keep(r: str) -> bool:
        if not_compared(r) or any(re.match(p, r) for p in skip):
            return False
        return not only or re.search(only, r) is not None

    refs = [r for r in outputs(ref_data, inputs) if keep(r)]
    acts = [r for r in outputs(act_data, inputs) if keep(r)]

    def masks(roots: Sequence[str]) -> Masks:
        m = Masks()
        for i, r in enumerate(roots):
            m = m.with_paths(**{("work" if i == 0 else f"root{i}"): r})
        return m

    # outputs under a listed exception are compared within its tolerance (tiff bytes then
    # follow the voxels); all others exactly, tiff bytes included
    exc = exceptions(kind)
    per_file: dict[str, Any] = {}
    for r in refs:
        e = next((e for e in exc if re.match(e.pattern, r)), None)
        if e is not None:
            per_file[r] = per_file[r + "/*"] = e.spec
    covered = {r for r in per_file if not r.endswith("/*")}
    out = Report()
    for part, tiff_bytes in ((False, True), (True, False)):
        rp = [r for r in refs if (r in covered) == part]
        ap = [r for r in acts if (r in covered) == part]
        if not rp and not ap:
            continue
        spec = tree(
            ignore_vars=("t", "t1", "ip"),
            ignore_fields=MODE_FIELDS if kind == "mode" else (),
            text_format=kind != "mode",
            per_file=per_file,
            tiff_bytes=tiff_bytes,
        )
        out.extend(
            compare_trees(
                scan_tree(ref_data, rp),
                scan_tree(act_data, ap),
                spec,
                ref_masks=masks(ref_roots),
                act_masks=masks(act_roots),
            )
        )
    summary = {"ref": len(refs), "act": len(acts), "problems": len(out.problems)}
    return out, summary


def voxel_difference(ref: np.ndarray, act: np.ndarray) -> dict:
    """how two arrays of the same size differ: counts of values off by 1 and by more, the
    largest difference and where it is (index, both values), the relative L2 error"""
    if ref.shape != act.shape:
        raise ValueError(f"sizes differ: {ref.shape} and {act.shape}")
    out = {"n": int(ref.size), "differ": 0, "off_by_1": 0, "off_by_more": 0, "max": 0.0}
    num = den = 0.0
    worst = None
    step = max(1, (1 << 24) // max(1, int(np.prod(ref.shape[1:]))))
    for z0 in range(0, ref.shape[0], step):
        r = np.asarray(ref[z0 : z0 + step], dtype=np.float64)
        a = np.asarray(act[z0 : z0 + step], dtype=np.float64)
        d = np.abs(a - r)
        nz = d > 0
        out["differ"] += int(nz.sum())
        out["off_by_1"] += int((d == 1).sum())
        out["off_by_more"] += int((d > 1).sum())
        num += float(np.sum(d * d))
        den += float(np.sum(r * r))
        m = float(d.max()) if d.size else 0.0
        if m > out["max"]:
            k = np.unravel_index(int(np.argmax(d)), d.shape)
            worst = ((z0 + int(k[0]), *map(int, k[1:])), float(r[k]), float(a[k]))
            out["max"] = m
    out["rel_l2"] = (num**0.5 / den**0.5) if den else (0.0 if num == 0 else float("inf"))
    out["frac"] = out["differ"] / out["n"] if out["n"] else 0.0
    out["worst"] = worst
    return out


# ---------------------------------------------------------------------------------------------
# tests without data
# ---------------------------------------------------------------------------------------------


def fake_dataset(root: Path) -> str:
    """a dataset folder with the FF layout of the demo cell dataset (empty files)"""
    for ch in CHANNELS:
        cam = ch.split("_")[0]
        for rel in (
            f"FF/averaged/ff_{ch}_stack0000_0000t.tif",
            f"FF/cam_folder/AVG_DF400_{cam}_10ms.tif",
        ):
            (root / rel).parent.mkdir(parents=True, exist_ok=True)
            (root / rel).write_bytes(b"")
    return _slash(root)


def test_function_files_are_written_as_a_production_script_writes_them(tmp_path):
    from petakit5d.io import read_mat

    data = fake_dataset(tmp_path)
    files = function_files(data)
    assert [os.path.basename(p) for p, _ in files] == [
        "processFunction_c0001.mat",
        "processFunction_c0002.mat",
    ]
    # MATLAB's %d of 0.4 is '4.000000e-01'
    assert files[0][1] == (
        f"@(x)processFFCorrectionFrame(x,'{data}FF/averaged/ff_CamA_ch0_stack0000_0000t.tif',"
        f"'{data}FF/cam_folder/AVG_DF400_CamA_10ms.tif',1,4.000000e-01,false)"
    )
    assert "FF/cam_folder/AVG_DF400_CamB_10ms.tif',1," in files[1][1]
    assert "ff_CamB_ch0_" in files[1][1]
    paths = write_function_files(data)
    assert read_mat(paths[1])["usrFun"] == files[1][1]


def test_the_flat_field_function_parses_as_the_tile_function(tmp_path):
    """the production form of the function text (MATLAB's %d of 0.4) is one the port's tile
    functions read, with its arguments"""
    from petakit5d.io import write_tiff
    from petakit5d.stitch.tile_functions import parse_tile_function

    data = fake_dataset(tmp_path)
    (ff, bg), _ = flat_field_files(data)
    rng = np.random.default_rng(3)
    write_tiff(rng.uniform(0.5, 1.0, (9, 7)).astype(np.float32), ff)
    write_tiff(np.full((9, 7), 100, np.uint16), bg)
    fn = parse_tile_function(function_files(data)[0][1])
    vol = rng.integers(100, 1000, (4, 9, 7)).astype(np.uint16)
    out = fn(vol)
    from petakit5d.flatfield import flatfield_correct_file

    ref = flatfield_correct_file(vol, ff, bg, 1, 0.4, False)
    assert out.dtype == ref.dtype and np.array_equal(out, ref)


def test_every_step_calls_a_registered_function_with_its_parameters(tmp_path):
    from petakit5d.params.registry import REGISTRY

    data = fake_dataset(tmp_path)
    for gpu in ("", "a100"):
        steps = workflow(data, job_log_dir=str(tmp_path / "logs"), gpu_name=gpu)
        assert [s.name for s in steps] == list(STEP_NAMES[1:])
        for s in steps:
            spec = REGISTRY.find(s.fun)
            names = [p.matlab.lower() for p in spec.params]
            for k in s.kw:
                # inputParser's unique prefixes (the demos pass 'resultDir')
                hits = [n for n in names if n == k.lower()] or [
                    n for n in names if n.startswith(k.lower())
                ]
                assert len(hits) == 1, (s.name, k)
            for c in s.columns:
                assert c in s.kw
    gpu = workflow(data, job_log_dir="logs", gpu_name="a100", gpu_config_file="g.json")
    decon = gpu[2].kw
    assert decon["GPUJob"] and decon["resultDirName"] == "matlab_decon_omw_a100"
    assert gpu[3].args == [[data + "matlab_stitch/matlab_decon_omw_a100/"]]


def test_the_spec_for_matlab_holds_the_same_calls(tmp_path):
    data = fake_dataset(tmp_path)
    steps = workflow(data, job_log_dir="logs", config_file="cpu.json")
    fn = tmp_path / "spec.json"
    write_spec(str(fn), data, steps)
    spec = json.loads(fn.read_text())
    assert [f["usrFun"] for f in spec["files"]] == [u for _, u in function_files(data)]
    assert [s["fun"] for s in spec["steps"]] == [s.fun for s in steps]
    assert spec["steps"][1]["kw"]["processFunPath"] == [p for p, _ in function_files(data)]
    assert spec["steps"][0]["columns"] == spec["steps"][1]["columns"] == ["processFunPath"]
    assert all(not s["columns"] for s in spec["steps"][2:])
    assert list(spec["steps"][2]["kw"]) == list(steps[2].kw)  # the order of the pairs
    assert b"\r\n" not in fn.read_bytes()


def test_steps_are_selected_in_the_workflows_order(tmp_path):
    steps = workflow(fake_dataset(tmp_path), job_log_dir="logs")
    assert [s.name for s in select(steps, ["mip_dsr", "stitch"])] == ["stitch", "mip_dsr"]
    assert len(select(steps, None)) == len(steps)
    with pytest.raises(ValueError, match="unknown steps"):
        select(steps, ["stitching"])


def test_the_dataset_is_linked(tmp_path):
    src = tmp_path / "src"
    (src / "PSF").mkdir(parents=True)
    big = src / "tile.tif"
    with open(big, "wb") as f:
        f.truncate((1 << 24) + 1)
    (src / "PSF" / "psf.tif").write_bytes(b"psf")
    (src / IMAGE_LIST).write_text("a,b\n", newline="\n")
    dst = link_dataset(str(src), str(tmp_path / "run"))
    assert os.path.islink(dst + "tile.tif") and not os.path.islink(dst + "PSF/psf.tif")
    assert Path(dst + IMAGE_LIST).read_text() == "a,b\n"
    assert link_dataset(str(src), dst) == dst  # again: nothing to do


def test_outputs_leave_out_inputs_job_logs_and_empty_folders(tmp_path):
    (tmp_path / "in.tif").write_bytes(b"x")
    (tmp_path / "a" / "b.zarr").mkdir(parents=True)
    (tmp_path / "a" / "b.zarr" / ".zarray").write_text("{}", newline="\n")
    (tmp_path / "empty").mkdir()
    (tmp_path / "job_logs").mkdir()
    (tmp_path / "job_logs" / "job_1.out").write_text("x", newline="\n")
    (tmp_path / STITCH / "xcorr" / "v").mkdir(parents=True)
    (tmp_path / STITCH / "xcorr" / "v" / "xcorr_tile_1.mat").write_bytes(b"x")
    (tmp_path / STITCH / "xcorr" / "v.mat").write_bytes(b"x")
    rels = outputs(str(tmp_path), ["in.tif"])
    assert "a/b.zarr" in rels and "in.tif" not in rels and "empty" not in rels
    # the registration results, per pair and merged, are compared
    assert [r for r in rels if not not_compared(r)] == [
        "a/b.zarr",
        f"{STITCH}/xcorr/v.mat",
        f"{STITCH}/xcorr/v/xcorr_tile_1.mat",
    ]


def test_compare_runs_masks_the_roots_and_applies_the_exceptions(tmp_path):
    from petakit5d.io import write_tiff

    m, p = tmp_path / "m", tmp_path / "p"
    for side in (m, p):
        d = side / STITCH / DECON / "MIPs"
        d.mkdir(parents=True)
        (side / STITCH / "parameters.json").write_text(
            json.dumps({"dataPath": str(side) + "/x"}), newline="\n"
        )
    a = np.arange(2000, dtype=np.uint16).reshape(4, 25, 20)
    b = a.copy()
    b[1, 2, 3] += 1
    write_tiff(a, str(m / STITCH / DECON / "MIPs" / "v_MIP_z.tif"))
    write_tiff(b, str(p / STITCH / DECON / "MIPs" / "v_MIP_z.tif"))
    rep, summary = compare_runs(str(m), str(p), ref_roots=[str(m)], act_roots=[str(p)])
    assert summary["ref"] == summary["act"] == 2
    assert [q.path.split(" ")[0] for q in rep.problems] == [
        f"files/{STITCH}/{DECON}/MIPs/v_MIP_z.tif"
    ]
    for kind in ("mode", "gpu"):  # no exceptions: exact
        rep, _ = compare_runs(str(m), str(p), kind=kind, ref_roots=[str(m)], act_roots=[str(p)])
        assert not rep.ok
    rep, _ = compare_runs(str(m), str(p), kind="permissive", ref_roots=[str(m)], act_roots=[str(p)])
    assert rep.ok
    b[0, 0, 0] += 2  # beyond 1 count
    write_tiff(b, str(p / STITCH / DECON / "MIPs" / "v_MIP_z.tif"))
    rep, _ = compare_runs(str(m), str(p), kind="permissive", ref_roots=[str(m)], act_roots=[str(p)])
    assert not rep.ok


def test_registration_peaks_are_within_p7_on_the_default_backend(tmp_path):
    from petakit5d.io import write_mat

    rel = f"{STITCH}/xcorr/v/xcorr_tile_1.mat"
    for side, peak, shift in (("m", 0.9123456, 3.0), ("p", 0.9123459, 3.0), ("q", 0.9123459, 4.0)):
        f = tmp_path / side / rel
        f.parent.mkdir(parents=True)
        write_mat(
            str(f),
            {
                "relative_shift_mat": np.array([[shift, 0.0, -2.0]]),
                "max_xcorr_mat": np.array([[1.0, 2.0, peak]]),
                "pair_indices": np.array([[1.0, 2.0]]),
            },
            version="7.3",
        )
    m, p, q = (str(tmp_path / s) for s in "mpq")
    assert not compare_runs(m, p, kind="exact")[0].ok
    rep, summary = compare_runs(m, p, kind="permissive")
    assert rep.ok and summary["ref"] == 1, rep.text()
    # a shift is never within it
    assert not compare_runs(m, q, kind="permissive")[0].ok


def test_mode_comparison_leaves_out_what_the_modes_differ_in(tmp_path):
    from petakit5d.io import write_tiff

    trees = {}
    for mode, task in (("local", 56), ("cluster", 20)):
        root = tmp_path / mode
        info = root / STITCH / "stitchInfo"
        info.mkdir(parents=True)
        (info / f"v_block_info_task_size_{task}.mat").write_bytes(b"x")
        params = {"parseCluster": mode == "cluster", "configFile": "", "xyPixelSize": 0.108}
        (root / STITCH / "parameters.json").write_text(json.dumps(params), newline="\n")
        (root / STITCH / f"{DECON}_gpu1").mkdir()
        vol = np.full((2, 3, 4), len(mode), np.uint16)
        write_tiff(vol, str(root / STITCH / f"{DECON}_gpu1" / "a.tif"))
        trees[mode] = str(root)
    rep, summary = compare_runs(trees["local"], trees["cluster"], kind="mode")
    assert rep.ok and summary["ref"] == 1, rep.text()
    rep, _ = compare_runs(trees["local"], trees["cluster"], kind="exact")
    assert len(rep.problems) == 3  # both stitch information files and parseCluster
    rep, summary = compare_runs(trees["local"], trees["cluster"], kind="gpu", only=GPU_DECON)
    assert summary["ref"] == 1 and not rep.ok


def test_voxel_difference_counts_and_locates():
    a = np.zeros((300, 200, 300), np.uint16)
    b = a.copy()
    b[5, 6, 7] = 1
    b[299, 1, 2] = 9
    d = voxel_difference(a, b)
    assert (d["differ"], d["off_by_1"], d["off_by_more"], d["max"]) == (2, 1, 1, 9.0)
    assert d["worst"] == ((299, 1, 2), 0.0, 9.0)
    assert d["rel_l2"] == float("inf")
    assert voxel_difference(a, a)["rel_l2"] == 0.0
    c = np.full((3, 4, 5), 10, np.uint16)
    e = c.copy()
    e[2, 3, 4] = 13
    assert voxel_difference(c, e)["rel_l2"] == 3 / (60 * 100) ** 0.5
    with pytest.raises(ValueError, match="sizes differ"):
        voxel_difference(a, b[:2])


# ---------------------------------------------------------------------------------------------
# tests on the real data
# ---------------------------------------------------------------------------------------------


@pytest.fixture
def dataset() -> str:
    """the demo cell dataset (PETAKIT5D_DEMO_DATA: the folder or its parent), or a skip"""
    value = os.environ.get("PETAKIT5D_DEMO_DATA", "")
    path = Path(value)
    if (path / CELL).is_dir():
        path = path / CELL
    if not value or not (path / IMAGE_LIST).is_file():
        pytest.skip(f"demodata: {value!r} does not contain the demo cell dataset")
    return _slash(path)


def _matlab_run() -> str:
    value = os.environ.get("PETAKIT5D_REALDATA_MATLAB", "")
    if not value or not os.path.isdir(os.path.join(value, CELL)):
        pytest.skip("PETAKIT5D_REALDATA_MATLAB does not hold a MATLAB run of the workflow")
    return _slash(os.path.join(value, CELL))


def _inputs(dataset: str) -> list[str]:
    from petakit5d.testing.compare import list_tree

    return list_tree(dataset)


@pytest.mark.demodata
@pytest.mark.slow
def test_local_workflow_matches_matlab(dataset, tmp_path):
    """the workflow in local mode on the CPU gives MATLAB's run's outputs (bit for bit with
    the exact FFT add-on, else within the listed exceptions)"""
    from petakit5d import fft

    mdata = _matlab_run()
    data = link_dataset(dataset, str(tmp_path / CELL))
    logs = str(tmp_path / "job_logs")
    run_python(data, workflow(data, job_log_dir=logs, parse_cluster=False))
    kind = "exact" if fft.backend() == "exact" else "permissive"
    rep, summary = compare_runs(
        mdata,
        data,
        inputs=_inputs(dataset),
        kind=kind,
        ref_roots=[os.path.dirname(mdata.rstrip("/"))],
        act_roots=[str(tmp_path)],
    )
    assert summary["act"] > 0
    assert rep.ok, rep.text("local workflow against MATLAB")


@pytest.mark.demodata
@pytest.mark.cluster
@pytest.mark.slow
def test_cluster_mode_gives_the_local_tree(dataset, tmp_path):
    """the workflow's trees in cluster mode (jobs from PETAKIT5D_TEST_SLURM_CONFIG) and in local
    mode are the same"""
    config = os.environ["PETAKIT5D_TEST_SLURM_CONFIG"]
    trees = {}
    for mode in ("local", "cluster"):
        root = tmp_path / mode
        data = link_dataset(dataset, str(root / CELL))
        cluster = mode == "cluster"
        steps = workflow(
            data,
            job_log_dir=str(root / "job_logs"),
            parse_cluster=cluster,
            config_file=config if cluster else "",
        )
        run_python(data, steps)
        trees[mode] = (data, str(root))
    (ld, lr), (cd, cr) = trees["local"], trees["cluster"]
    rep, summary = compare_runs(
        ld, cd, inputs=_inputs(dataset), kind="mode", ref_roots=[lr], act_roots=[cr]
    )
    assert summary["act"] > 0
    assert rep.ok, rep.text("cluster against local mode")


@pytest.mark.demodata
@pytest.mark.gpu
@pytest.mark.slow
def test_gpu_decon_matches_matlab(dataset, tmp_path):
    """the large-file OMW decon of the stitched volumes on this GPU, locally, gives MATLAB's GPU
    run's decon (PETAKIT5D_REALDATA_GPU names its result folder, e.g. matlab_decon_omw_a100)"""
    mdata = _matlab_run()
    name = os.environ.get("PETAKIT5D_REALDATA_GPU", "")
    if not name or not os.path.isdir(mdata + STITCH + "/" + name):
        pytest.skip("PETAKIT5D_REALDATA_GPU does not name a GPU decon folder of the MATLAB run")
    data = link_dataset(dataset, str(tmp_path / CELL))
    os.makedirs(data + STITCH, exist_ok=True)
    for p in sorted(glob.glob(mdata + STITCH + "/*.zarr")):
        os.symlink(p, data + STITCH + "/" + os.path.basename(p))
    gpu = name.removeprefix(DECON + "_")
    steps = workflow(
        data, job_log_dir=str(tmp_path / "job_logs"), parse_cluster=False, gpu_name=gpu
    )
    run_python(data, select(steps, ["decon"]), function_files_too=False)
    rep, summary = compare_runs(
        mdata,
        data,
        inputs=_inputs(dataset),
        kind="gpu",
        only=rf"^{STITCH}/{re.escape(name)}/",
        ref_roots=[os.path.dirname(mdata.rstrip("/"))],
        act_roots=[str(tmp_path)],
    )
    assert summary["act"] > 0
    assert rep.ok, rep.text("GPU decon against MATLAB")
