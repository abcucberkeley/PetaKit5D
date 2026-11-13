"""The argument strings the GUI and PyPetaKit5D hand the compiled master, made as they make them.

Both start MATLAB's ``run_mccMaster.sh <Runtime> <function> <arguments>``:

- **PetaKit5D-GUI** (src/matlabthread.cpp, src/matlabhelperfunctions.cpp, the submit functions
  of src/mainwindow.cpp) builds one command line: the helpers ``addCharArrayToArgs`` (text in
  double quotes), ``addBoolToArgs`` (``true``/``false``), ``addScalarToArgs`` (the text of a
  field, or ``std::to_string`` of a number: ``%f``), ``addArrayToArgs`` (``"[a,b]"`` or
  ``"{'a','b'}"``, any separator), ``addChannelPatternsToArgs`` and ``addDataPathsToArgs``.
  It then replaces every blank and two quotes (an empty value) by a blank and six quotes, and
  starts the line with
  ``QProcess::startCommand``, which splits it as ``QProcess::splitCommand`` does.
- **PyPetaKit5D** (_mcc.py and the generated wrappers) passes an argv list: cells
  ``{'a','b'}`` with ``'`` doubled, arrays ``[1,2]`` of ``str()`` of the items, logicals
  ``true``/``false``, scalars by ``str()``, text as it is; it escapes ``\\ " $ ``` for the
  launcher.
- **run_mccMaster.sh** wraps every argument after the Runtime folder in double quotes and
  evaluates the line (``eval "<master>" $args``), which undoes PyPetaKit5D's escapes and makes
  the GUI's ``""`` an empty argument.

``gui_calls()`` and ``pypetakit5d_calls()`` give, per Phase 1 function, the argument lists
mccMaster receives for representative states of the GUI's tabs and calls of the wrappers.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).parent

ROOT = "/mnt/projects/embryo_2024"
PATHS = [f"{ROOT}/cell1", f"{ROOT}/cell 2"]  # a space, as folder names can have
RUNTIME = "/opt/MATLAB_Runtime/R2024b"
LAUNCHER = "/opt/PetaKit5D-GUI/PetaKit5D/mcc/linux/run_mccMaster.sh"


# ----------------------------------------------------------------------------------------------
# the launchers


def qt_split_command(command: str) -> list[str]:
    """QProcess::splitCommand: words split at white space outside double quotes; three
    consecutive quotes are one literal quote; empty words are dropped."""
    args: list[str] = []
    tmp = ""
    quotes = 0
    in_quote = False
    for c in command:
        if c == '"':
            quotes += 1
            if quotes == 3:
                quotes = 0
                tmp += c
            continue
        if quotes:
            if quotes == 1:
                in_quote = not in_quote
            quotes = 0
        if not in_quote and c.isspace():
            if tmp:
                args.append(tmp)
                tmp = ""
        else:
            tmp += c
    if tmp:
        args.append(tmp)
    return args


def gui_fix_empty(command: str) -> str:
    """matlabthread.cpp: a blank and two quotes become a blank and six quotes (std::regex_replace),
    so that QProcess keeps an empty value as a word of two quotes."""
    return command.replace(' ""', ' """"""')


def launcher_eval(tokens: list[str]) -> list[str]:
    """The words of ``eval <master> $args`` in run_mccMaster.sh, where ``args`` is every token
    wrapped in double quotes and joined with spaces (sh double-quote and word rules; tokens
    with ``$`` or a backquote outside an escape are refused: they would be expanded)."""
    line = "".join(f' "{t}"' for t in tokens)
    words: list[str] = []
    cur: str | None = None
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c in " \t\n":
            if cur is not None:
                words.append(cur)
                cur = None
            i += 1
        elif c == '"':
            cur = cur or ""
            i += 1
            while True:
                if i >= n:
                    raise ValueError(f"unterminated quote in {line!r}")
                c = line[i]
                if c == '"':
                    i += 1
                    break
                if c == "\\" and i + 1 < n and line[i + 1] in '$`"\\\n':
                    cur += line[i + 1]
                    i += 2
                    continue
                if c in "$`":
                    raise ValueError(f"the launcher would expand {line[i:]!r}")
                cur += c
                i += 1
        elif c == "\\" and i + 1 < n:
            cur = (cur or "") + line[i + 1]
            i += 2
        elif c in "$`'|&;<>()":
            raise ValueError(f"the launcher would interpret {c!r} in {line!r}")
        else:
            cur = (cur or "") + c
            i += 1
    if cur is not None:
        words.append(cur)
    return words


# ----------------------------------------------------------------------------------------------
# the GUI


class GuiArgs:
    """The GUI's argument builder in mcc mode (prependedString ' ', surroundQuotes '"')."""

    def __init__(self) -> None:
        self.text = ""

    def char(self, s: str, pre: str = " ") -> GuiArgs:  # addCharArrayToArgs, isMcc
        self.text += f'{pre}"{s}"'
        return self

    def bool(self, b: bool) -> GuiArgs:  # addBoolToArgs
        self.text += " true" if b else " false"
        return self

    def scalar(self, s: str) -> GuiArgs:  # addScalarToArgs
        self.text += f" {s}"
        return self

    def array(self, items: list[str], is_char: bool, chars: str = "[]", sep: str = ",") -> GuiArgs:
        body = sep.join(f"'{x}'" if is_char else x for x in items)
        self.text += f' "{chars[0]}{body}{chars[1]}"'
        return self

    def channels(self, patterns: list[str] | None = None, custom: str | None = None) -> GuiArgs:
        """addChannelPatternsToArgs: the checked channel labels, or the custom pattern line
        split at commas (spaces kept)."""
        if custom is not None:
            patterns = custom.split(",")
            if patterns and patterns[-1] == "":
                patterns = patterns[:-1]
        return self.array(list(patterns or []), True, "{}")

    def data_paths(self, paths: list[str]) -> GuiArgs:  # addDataPathsToArgs: the first argument
        self.text += '"{' + ",".join(f"'{p}'" for p in paths) + '}"'
        return self

    def opt(self, name: str, kind: str, value) -> GuiArgs:
        """A name-value pair: the name with addCharArrayToArgs, the value by ``kind``."""
        self.char(name)
        if kind == "char":
            return self.char(value)
        if kind == "bool":
            return self.bool(value)
        if kind == "scalar":
            return self.scalar(value)
        if kind == "num":
            return self.array(value, False)
        if kind == "num;":
            return self.array(value, False, "[]", ";")
        if kind == "cells":
            return self.array(value, True, "{}")
        if kind == "cells;":
            return self.array(value, True, "{}", ";")
        if kind == "channels":
            return self.channels(*value) if isinstance(value, tuple) else self.channels(value)
        raise ValueError(kind)


def gui_master_argv(function: str, args: GuiArgs) -> list[str]:
    """What mccMaster gets from a GUI command line: [function, arguments...]."""
    cmd = f'"{LAUNCHER}" "{RUNTIME}" {function} {args.text}'
    words = qt_split_command(gui_fix_empty(cmd))
    assert words[0] == LAUNCHER and words[1] == RUNTIME
    return launcher_eval(words[2:])


def ftos(x: float) -> str:
    """std::to_string of a double (printf %f)."""
    return f"{x:f}"


CHANNELS = ["CamA_ch0", "CamB_ch0"]


def _gui_main(function: str, s: dict) -> GuiArgs:
    """The main tab's submit (on_submitButton_clicked): deskew/rotate, stitching, decon only
    or the whole pipeline, by ``function``."""
    dsr = function == "XR_deskew_rotate_data_wrapper"
    stitch = function == "XR_matlab_stitching_wrapper"
    decon_only = function == "XR_decon_data_wrapper"
    a = GuiArgs().data_paths(s["paths"])
    if stitch:
        a.array(s["image_lists"], True, "{}")
    a.opt("dz", "scalar", s["dz"])
    a.opt("xyPixelSize", "scalar", ftos(s["xy"]))
    a.opt("channelPatterns", "channels", s["channels"])
    a.opt("skewAngle", "scalar", ftos(s["skew"]))
    if not decon_only:
        a.opt("reverse", "bool", s["reverse"])
        a.opt("objectiveScan", "bool", s["objective_scan"])
    a.opt("parseSettingFile", "bool", s["parse_setting_file"])
    if dsr:
        if s["dsr_combined"]:
            a.opt("DSRCombined", "bool", True).opt("Deskew", "bool", True)
            a.opt("Rotate", "bool", True)
        else:
            a.opt("DSRCombined", "bool", False).opt("Deskew", "bool", s["deskew"])
            a.opt("Rotate", "bool", s["rotate"])
        a.opt("flipZstack", "bool", s["flip"]).opt("save16bit", "bool", s["save16"])
        a.opt("saveMIP", "bool", s["save_mip"]).opt("largeFile", "bool", s["large"])
        a.opt("zarrFile", "bool", s["zarr"]).opt("saveZarr", "bool", s["save_zarr"])
        a.opt("batchSize", "num", s["batch"]).opt("blockSize", "num", s["block"])
        if s.get("bbox"):
            a.opt("InputBbox", "num", s["bbox"])
        if s.get("masks"):
            a.opt("maskFullpaths", "cells", s["masks"])
        a.opt("masterCompute", "bool", s["master_compute"])
    elif stitch:
        a.opt("streaming", "bool", s["streaming"]).opt("multiLoc", "bool", False)
        if s.get("processed_dir"):
            a.opt("processedDirStr", "char", s["processed_dir"])
        a.opt("DS", "bool", s["ds"]).opt("DSR", "bool", s["dsr"])
        a.opt("IOScan", "bool", False).opt("zarrFile", "bool", s["zarr"])
        a.opt("largeFile", "bool", s["large"])
        if s.get("axis_order"):
            a.opt("axisOrder", "char", s["axis_order"])
        a.opt("blockSize", "num", s["block"])
        if s.get("shard"):
            a.opt("shardSize", "num", s["shard"])
        a.opt("resampleType", "char", s["resample_type"])
        if s.get("resample"):
            a.opt("resampleFactor", "num", [ftos(x) for x in s["resample"]])
        if s.get("bbox"):
            a.opt("inputBbox", "num", s["bbox"])
        a.opt("TileOffset", "scalar", "0")
        if s.get("result_dir"):
            a.opt("resultDir", "char", s["result_dir"])
        a.opt("BlendMethod", "char", s["blend"])
        if s.get("overlap_type"):
            a.opt("overlapType", "char", s["overlap_type"])
        a.opt("xcorrShift", "bool", s["xcorr_shift"]).opt("xyMaxOffset", "scalar", "300")
        a.opt("zMaxOffset", "scalar", "50").opt("xcorrDownsample", "num", ["2", "2", "1"])
        a.opt("xcorrThresh", "scalar", "0.25")
        if s.get("out_bbox"):
            a.opt("outBbox", "num", s["out_bbox"])
        a.opt("xcorrMode", "char", "primaryFirst").opt("shiftMethod", "char", "grid")
        a.opt("axisWeight", "num", ["1", ".1", "10"])
        if s.get("primary"):
            a.opt("primaryCh", "char", s["primary"])
        a.opt("usePrimaryCoords", "bool", False).opt("save16bit", "bool", s["save16"])
        a.opt("edgeArtifacts", "scalar", "2")
        a.opt("stitchMIP", "num", ["false", "false", s["stitch_mip_z"]])
        a.opt("onlineStitch", "bool", False)
        if s.get("process_fun"):
            a.opt("processFunPath", "cells;", s["process_fun"])
        a.opt("masterCompute", "bool", s["master_compute"])
    elif decon_only:
        a.opt("save16bit", "bool", s["save16"])
        a.opt("RLMethod", "char", s["rl_method"])
        if s.get("result_dir"):
            a.opt("resultDirName", "char", s["result_dir"])
        if s["rl_method"] == "omw":
            a.opt("wienerAlpha", "num", s["wiener"]).opt("OTFCumThresh", "num", s["otf"])
            a.opt("hannWinBounds", "num", ["0.8", "1.0"])
            if s.get("skewed_manual"):
                a.opt("skewed", "bool", s["skewed"])
        a.opt("background", "scalar", s["background"]).opt("dzPSF", "scalar", s["dz_psf"])
        a.opt("edgeErosion", "scalar", s["erosion"]).opt("ErodeByFTP", "bool", True)
        if s.get("psfs"):
            a.opt("psfFullpaths", "cells", s["psfs"])
        if s["rl_method"] != "original":
            a.opt("debug", "bool", False).opt("saveStep", "scalar", "5")
        a.opt("psfGen", "bool", s["psf_gen"]).opt("GPUJob", "bool", s["gpu"])
        a.opt("deconIter", "scalar", s["iter"]).opt("deconRotate", "bool", False)
        a.opt("dampFactor", "scalar", ftos(1.0))
        if s.get("scale"):
            a.opt("scaleFactor", "num", s["scale"])
        a.opt("deconOffset", "scalar", "0")
        if s.get("masks"):
            a.opt("maskFullpaths", "cells", s["masks"])
        a.opt("batchSize", "num", s["batch"]).opt("blockSize", "num", s["block"])
        a.opt("largeFile", "bool", s["large"]).opt("largeMethod", "char", "inmemory")
        a.opt("zarrFile", "bool", s["zarr"]).opt("saveZarr", "bool", s["save_zarr"])
        if s.get("gpu_config"):
            a.opt("GPUConfigFile", "char", s["gpu_config"])
    else:  # the pipeline (XR_microscopeAutomaticProcessing)
        a.opt("streaming", "bool", s["streaming"]).opt("dzFromEncoder", "bool", False)
        a.opt("resampleType", "char", s["resample_type"])
        if s.get("resample"):
            a.opt("resampleFactor", "num", [ftos(x) for x in s["resample"]])
        a.opt("save16bit", "num", [s["save16_dsr"], s["save16_stitch"]])
        a.opt("zarrFile", "bool", s["zarr"]).opt("saveZarr", "bool", s["save_zarr"])
        a.opt("save3DStack", "bool", True)
        a.opt("DSRCombined", "bool", True).opt("Deskew", "bool", True)
        a.opt("Rotate", "bool", True).opt("Stitch", "bool", s["stitch"])
        a.opt("flipZstack", "bool", s["flip"]).opt("inputAxisOrder", "char", "yxz")
        a.opt("outputAxisOrder", "char", "yxz").opt("FFCorrection", "bool", s["ff"])
        if s.get("ff_images"):
            a.opt("FFImagePaths", "cells", s["ff_images"])
            a.opt("backgroundPaths", "cells", s["backgrounds"])
        a.opt("BKRemoval", "bool", False).opt("LowerLimit", "scalar", ftos(0.4))
        if s.get("result_dir"):
            a.opt("stitchResultDirName", "char", s["result_dir"])
        if s.get("image_lists"):
            a.opt("imageListFullpaths", "cells", s["image_lists"])
        if s.get("axis_order"):
            a.opt("axisOrder", "char", s["axis_order"])
        a.opt("BlendMethod", "char", s["blend"]).opt("xyMaxOffset", "scalar", "300")
        a.opt("zMaxOffset", "scalar", "50").opt("xcorrShift", "bool", True)
        a.opt("xcorrMode", "char", "primaryFirst")
        if s.get("out_bbox"):
            a.opt("boundboxCrop", "num", s["out_bbox"])
        if s.get("primary"):
            a.opt("primaryCh", "char", s["primary"])
        a.opt("stitchMIP", "num", ["false", "false", s["stitch_mip_z"]])
        a.opt("onlineStitch", "bool", False)
        a.opt("minModifyTime", "scalar", "1").opt("maxModifyTime", "scalar", "10")
        a.opt("maxWaitLoopNum", "scalar", "10")
    a.opt("parseCluster", "bool", s["parse_cluster"]).opt("cpusPerTask", "scalar", s["cpus"])
    if s.get("uuid"):
        a.opt("uuid", "char", s["uuid"])
    if not dsr:
        a.opt("maxTrialNum", "scalar", "3").opt("unitWaitTime", "scalar", "1")
    a.opt("mccMode", "bool", True)
    if s.get("config"):
        a.opt("configFile", "char", s["config"])
    if s.get("job_log_dir"):
        a.opt("jobLogDir", "char", s["job_log_dir"])
    return a


MAIN_DEFAULTS = {
    "paths": PATHS[:1],
    "image_lists": [f"{ROOT}/ImageList_from_encoder.csv"],
    "dz": "0.5",
    "xy": 0.108,
    "channels": CHANNELS,
    "skew": 32.45,
    "reverse": True,
    "objective_scan": False,
    "parse_setting_file": False,
    "dsr_combined": True,
    "deskew": True,
    "rotate": True,
    "flip": False,
    "save16": True,
    "save_mip": True,
    "large": False,
    "zarr": False,
    "save_zarr": False,
    "batch": ["1024", "1024", "1024"],
    "block": ["256", "256", "256"],
    "master_compute": True,
    "streaming": False,
    "ds": False,
    "dsr": False,
    "resample_type": "isotropic",
    "blend": "feather",
    "xcorr_shift": True,
    "stitch_mip_z": "false",
    "rl_method": "omw",
    "wiener": ["0.005", "0.005"],
    "otf": ["0.88", "0.88"],
    "background": "100",
    "dz_psf": "0.1",
    "erosion": "0",
    "psf_gen": True,
    "gpu": False,
    "iter": "2",
    "save16_dsr": "true",
    "save16_stitch": "true",
    "stitch": False,
    "ff": False,
    "parse_cluster": False,
    "cpus": "4",
    "job_log_dir": f"{ROOT}/cell1/job_logs/20261113_0930_job",
}

MAIN_ALL = {
    **MAIN_DEFAULTS,
    "paths": PATHS,
    "xy": 0.1,  # printed as 0.100000
    "skew": -31.8,
    "dz": "0.25",
    "channels": (None, "CamA_ch0, CamB_ch1"),  # a custom pattern line: the space is kept
    "dsr_combined": False,
    "deskew": True,
    "rotate": False,
    "flip": True,
    "large": True,
    "zarr": True,
    "save_zarr": True,
    "batch": ["512", "1024", "2048"],
    "block": ["128", "256", "64"],
    "bbox": ["1", "1", "1", "100", "200", "50"],
    "masks": [f"{ROOT}/masks/xy.tif", f"{ROOT}/masks/xz.tif", f"{ROOT}/masks/yz.tif"],
    "master_compute": False,
    "streaming": True,
    "processed_dir": "DSR",
    "ds": True,
    "axis_order": "-x,y,z",
    "shard": ["0", "0", "1"],
    "resample": [2.0, 2.0, 1.0],
    "result_dir": "matlab_stitch_2",
    "overlap_type": "half",
    "out_bbox": ["1", "1", "1", "1000", "2000", "300"],
    "primary": "CamA_ch0",
    "stitch_mip_z": "true",
    "process_fun": [f"{ROOT}/tmp/a_CamA_ch0.txt", f"{ROOT}/tmp/a_CamB_ch1.txt"],
    "rl_method": "omw",
    "skewed_manual": True,
    "skewed": True,
    "psfs": [f"{ROOT}/PSF/CamA.tif", f"{ROOT}/PSF/CamB.tif"],
    "scale": ["1", "2"],
    "gpu": True,
    "gpu_config": f"{ROOT}/configs/gpu.json",
    "save16_stitch": "false",
    "stitch": True,
    "ff": True,
    "ff_images": [f"{ROOT}/FF/CamA.tif", f"{ROOT}/FF/CamB.tif"],
    "backgrounds": [f"{ROOT}/FF/bgA.tif", f"{ROOT}/FF/bgB.tif"],
    "image_lists": [f"{ROOT}/ImageList_from_encoder.csv"],
    "parse_cluster": True,
    "cpus": "24",
    "uuid": "a1b2c3",
    "config": f"{ROOT}/configs/cluster.json",
}


def _gui_tools() -> dict[str, list[tuple[str, GuiArgs]]]:
    out: dict[str, list[tuple[str, GuiArgs]]] = {}

    def crop(paths, bbox, chan, flags, block, batch, cpus, uuid, job_log_dir):
        a = GuiArgs().data_paths(paths).array(bbox, False)
        a.opt("resultDirName", "char", "Cropped").opt("channelPatterns", "channels", chan)
        a.opt("save16bit", "bool", flags[0]).opt("pad", "bool", flags[1])
        a.opt("cropType", "char", "fixed").opt("lastStartCoords", "num", ["0", "0", "0"])
        a.opt("zarrFile", "bool", flags[2]).opt("largeFile", "bool", flags[3])
        a.opt("saveZarr", "bool", flags[4]).opt("blockSize", "num", block)
        a.opt("batchSize", "num", batch).opt("parseCluster", "bool", flags[5])
        a.opt("masterCompute", "bool", True).opt("cpusPerTask", "scalar", cpus)
        if uuid:
            a.opt("uuid", "char", uuid)
        a.opt("mccMode", "bool", True)
        if job_log_dir:
            a.opt("jobLogDir", "char", job_log_dir)
        return a

    out["XR_crop_dataset"] = [
        ("defaults", crop(PATHS[:1], ["1", "1", "1", "512", "512", "100"], CHANNELS,
                          [True, False, False, False, False, False], ["256"] * 3,
                          ["1024"] * 3, "2", "", "")),
        ("options", crop(PATHS, ["10", "20", "3", "500", "600", "70"], (None, "CamA,CamB,"),
                         [False, True, True, True, True, True], ["64", "128", "32"],
                         ["512", "256", "128"], "8", "x9", f"{ROOT}/job_logs")),
    ]  # fmt: skip

    def tiff_zarr(paths, to_zarr, block, flags, chan, cpus, uuid):
        if to_zarr:
            a = GuiArgs().data_paths(paths).opt("blockSize", "num", block)
            a.opt("readWholeTiff", "bool", flags[0]).opt("partialFile", "bool", flags[1])
        else:
            a = GuiArgs().data_paths(paths)
        a.opt("resultDirName", "char", "zarr" if to_zarr else "tiffs")
        a.opt("channelPatterns", "channels", chan).opt("parseCluster", "bool", flags[2])
        a.opt("cpusPerTask", "scalar", cpus).opt("masterCompute", "bool", True)
        if uuid:
            a.opt("uuid", "char", uuid)
        return a.opt("mccMode", "bool", True)

    out["XR_tiffToZarr_wrapper"] = [
        ("defaults", tiff_zarr(PATHS[:1], True, ["256.00"] * 3, [True, False, False], CHANNELS,
                               "1", "")),
        ("options", tiff_zarr(PATHS, True, ["64.00", "128.50", "32.00"], [False, True, True],
                              ["CamA"], "4", "u1")),
    ]  # fmt: skip
    out["XR_zarrToTiff_wrapper"] = [
        ("defaults", tiff_zarr(PATHS[:1], False, [], [False, False, False], CHANNELS, "1", "")),
        ("options", tiff_zarr(PATHS, False, [], [False, False, True], ["CamB_ch1"], "3", "u2")),
    ]

    def mip(paths, chan, axis, flags, batch, cpus):
        a = GuiArgs().data_paths(paths).opt("resultDirName", "char", "MIPs")
        a.opt("channelPatterns", "channels", chan)
        a.opt("axis", "num", [str(int(x)) for x in axis])  # std::to_string(bool)
        a.opt("zarrFile", "bool", flags[0]).opt("largeFile", "bool", flags[1])
        a.opt("batchSize", "num", batch).opt("save16bit", "bool", flags[2])
        a.opt("parseCluster", "bool", flags[3]).opt("masterCompute", "bool", True)
        return a.opt("cpusPerTask", "scalar", cpus).opt("mccMode", "bool", True)

    out["XR_MIP_wrapper"] = [
        ("defaults", mip(PATHS[:1], CHANNELS, [0, 0, 1], [False, False, True, False],
                         ["2048", "2048", "2048"], "1")),
        ("options", mip(PATHS, (None, "ch0"), [1, 1, 1], [True, True, False, True],
                        ["512", "256", "64"], "6")),
    ]  # fmt: skip

    def resample(paths, factor, chan, bbox, interp, flags, block, batch, border, cpus):
        a = GuiArgs().data_paths(paths).array(factor, False)
        a.opt("resultDirName", "char", "resampled").opt("channelPatterns", "channels", chan)
        if bbox:
            a.opt("inputBbox", "num", bbox)
        a.opt("interpMethod", "char", interp).opt("save16bit", "bool", flags[0])
        a.opt("zarrFile", "bool", flags[1]).opt("largeFile", "bool", flags[2])
        a.opt("saveZarr", "bool", flags[3]).opt("blockSize", "num", block)
        a.opt("batchSize", "num", batch).opt("borderSize", "num", border)
        a.opt("parseCluster", "bool", flags[4]).opt("masterCompute", "bool", True)
        return a.opt("cpusPerTask", "scalar", cpus).opt("mccMode", "bool", True)

    out["XR_resample_dataset"] = [
        ("defaults", resample(PATHS[:1], ["2.00", "2.00", "2.00"], CHANNELS, None, "linear",
                              [True, False, False, False, False], ["256"] * 3, ["512"] * 3,
                              ["5", "5", "5"], "1")),
        ("options", resample(PATHS, ["0.50", "1.50", "4.00"], ["CamA_ch0"],
                             ["1", "1", "1", "64", "64", "32"], "cubic",
                             [False, True, True, True, True], ["64", "64", "16"],
                             ["128", "128", "32"], ["0", "0", "0"], "2")),
    ]  # fmt: skip

    def image_list(paths, method, chan, patterns, flags, xy, dz, skew, overlap, extra):
        a = GuiArgs().data_paths(paths).char(method)
        a.opt("channelPatterns", "channels", chan).opt("mapTilename", "bool", flags[0])
        if patterns:
            a.opt("tilePatterns", "cells", patterns)
        a.opt("DS", "bool", flags[1]).opt("DSR", "bool", flags[2])
        a.opt("xyPixelSize", "scalar", xy).opt("dz", "scalar", dz)
        a.opt("skewAngle", "scalar", skew).opt("axisOrder", "char", "x,y,z")
        a.opt("dataOrder", "char", "y,x,z").opt("objectiveScan", "bool", False)
        a.opt("IOScan", "bool", False).opt("zarrFile", "bool", flags[3])
        if overlap:
            a.opt("overlapSize", "num", overlap)
        a.opt("overlapSizeType", "char", "pixel")
        if extra:
            a.opt("tileFilenames", "cells", extra[0]).opt("tileIndices", "num;", extra[1])
            a.opt("tileInterval", "num", extra[2])
        return a

    out["XR_generate_image_list_wrapper"] = [
        ("encoder", image_list(PATHS[:1], "encoder", ["CamA_ch0"], None,
                               [True, False, False, False], "0.108", "0.2", "32.45", None, None)),
        ("tile_list", image_list(PATHS, "tile_list", CHANNELS,
                                 ["0000t", "ch0", "000x", "000y", "000z"],
                                 [False, True, False, True], "0.1", "0.3", "-30", ["10", "12"],
                                 (["tile_a.tif", "tile_b.tif"], ["0 0 1 1 0", "0 0 2 1 0"],
                                  ["1", "2", "3"]))),
    ]  # fmt: skip
    return out


def gui_calls() -> list[tuple[str, str, list[str]]]:
    """(case name, function, mccMaster's arguments) for the GUI's Phase 1 functions."""
    out = []
    for function in (
        "XR_deskew_rotate_data_wrapper",
        "XR_matlab_stitching_wrapper",
        "XR_decon_data_wrapper",
        "XR_microscopeAutomaticProcessing",
    ):
        for name, state in (("defaults", MAIN_DEFAULTS), ("options", MAIN_ALL)):
            argv = gui_master_argv(function, _gui_main(function, state))
            out.append((name, function, argv[1:]))
    for function, calls in _gui_tools().items():
        for name, a in calls:
            out.append((name, function, gui_master_argv(function, a)[1:]))
    return out


# ----------------------------------------------------------------------------------------------
# PyPetaKit5D


def matlab_str(value) -> str:  # _mcc.matlab_str
    return "'" + str(value).replace("'", "''") + "'"


def matlab_cell(items) -> str:
    return "{" + ",".join(matlab_str(i) for i in items) + "}"


def matlab_array(items, separator: str = ",") -> str:
    return "[" + separator.join(str(i) for i in items) + "]"


def matlab_logical_array(items) -> str:
    return "[" + ",".join(str(i).lower() for i in items) + "]"


def escape_for_launcher(token: str) -> str:  # _mcc._escape_for_launcher
    for c in ("\\", '"', "$", "`"):
        token = token.replace(c, "\\" + c)
    return token


def wrappers() -> dict:
    return json.loads((HERE / "data" / "pypetakit5d_wrappers.json").read_text())["wrappers"]


_POS_ENCODERS = {
    "cell": matlab_cell,
    "array": matlab_array,
    "logical_array": matlab_logical_array,
    "text": str,
    "number_text": str,
    "bool_text": lambda b: str(b).lower(),
}


def pypetakit5d_argv(function: str, positional: list, **kwargs) -> list[str]:
    """The tokens a PyPetaKit5D wrapper passes (as run_mcc gets them, before its escaping)."""
    w = wrappers()[function]
    tokens = [
        _POS_ENCODERS[enc](v) for (_, enc), v in zip(w["positional"], positional, strict=True)
    ]
    for key, default, typ in w["options"]:
        value = kwargs.get(key, default)
        if typ == "char":
            if not value:
                continue
            tokens += [key, value]
        elif typ == "cell":
            if not value:
                continue
            tokens += [key, matlab_cell(value)]
        elif typ == "logicalArr":
            tokens += [key, matlab_logical_array(value)]
        elif typ == "logical":
            tokens += [key, str(value).lower()]
        elif typ == "numericArr":
            if not value:
                continue
            if type(value) is not list:
                value = [value]
            sep = (
                ";" if (function, key) == ("XR_generate_image_list_wrapper", "tileIndices") else ","
            )
            tokens += [key, matlab_array(value, sep)]
        elif typ == "numericScalar":
            if type(value) is list:
                if not value:
                    continue
                value = value[0]
            tokens += [key, value]
    return [str(t) for t in tokens]


def pypetakit5d_master_argv(function: str, positional: list, **kwargs) -> list[str]:
    """mccMaster's arguments: the wrapper's tokens escaped by run_mcc, then evaluated by the
    launcher."""
    tokens = pypetakit5d_argv(function, positional, **kwargs)
    words = launcher_eval([function] + [escape_for_launcher(t) for t in tokens])
    assert words == [function, *tokens], (words, tokens)
    return words[1:]


_PYP_POSITIONAL = {
    "XR_crop_dataset": [PATHS[:1], [1, 1, 1, 100, 200, 50]],
    "XR_generate_image_list_wrapper": [PATHS[:1], "encoder"],
    "XR_matlab_stitching_wrapper": [PATHS[:1], [f"{ROOT}/ImageList_from_encoder.csv"]],
    "XR_resample_dataset": [PATHS[:1], [2, 2, 2]],
}

# a second call per wrapper: values of each type as users give them
_PYP_OPTIONS = {
    "XR_crop_dataset": {
        "resultDirName": "crop 1", "pad": True, "lastStartCoords": [10, 20, 5],
        "channelPatterns": ["CamA_ch0", "it's"], "batchSize": [512, 256, 128],
        "cpusPerTask": 4, "uuid": "abc",
    },
    "XR_MIP_wrapper": {
        "axis": [1, 1, 0], "channelPatterns": ["CamA"], "zarrFile": True,
        "batchSize": [100, 200, 300], "cpusPerTask": 1.5, "jobLogDir": f"{ROOT}/logs $HOME",
    },
    "XR_decon_data_wrapper": {
        "xyPixelSize": 0.1, "dz": 1e-1, "Background": 99.5, "EdgeErosion": 3,
        "psfFullpaths": [f"{ROOT}/psf a.tif", f"{ROOT}/psf_b.tif"], "RLMethod": "omw",
        "wienerAlpha": [0.005, 0.01], "OTFCumThresh": [0.9, 0.88], "hannWinBounds": [0.8, 1],
        "GPUJob": True, "deconIter": 15, "scaleFactor": [1e8], "save16bit": False,
    },
    "XR_deskew_rotate_data_wrapper": {
        "xyPixelSize": 0.0975, "dz": 0.25, "skewAngle": -32.45, "flipZstack": True,
        "Rotate": True, "DSRCombined": False, "inputBbox": [1, 2, 3, 40, 50, 60],
        "resampleFactor": [1.5, 1.5, 2], "batchSize": [128, 128, 64], "largeFile": True,
    },
    "XR_generate_image_list_wrapper": {
        "tileFilenames": ["a.tif", "b.tif"], "tileIndices": [[0, 0, 1, 1, 0], [0, 0, 2, 1, 0]],
        "tileInterval": [1, 2, 3],
        "overlapSize": [12, 12, 4], "xyPixelSize": 0.11, "DSR": True,
    },
    "XR_matlab_stitching_wrapper": {
        "streaming": True, "axisOrder": "-x,y,z", "resampleFactor": [2, 2, 1],
        "stitchMIP": [False, False, True], "xcorrDownsample": [2, 2, 1],
        "inputBbox": [1, 1, 1, 10, 20, 30], "primaryCh": "CamA_ch0", "save16bit": False,
    },
    "XR_microscopeAutomaticProcessing": {
        "Stitch": True, "Decon": True, "save16bit": [True, False], "resampleFactor": [2, 2, 2],
        "FFImagePaths": ["/mnt/ff/a.tif", "/mnt/ff/b.tif"], "stitchMIP": [True, False, True],
        "wienerAlpha": [0.005], "GPUJob": True,
    },
    "XR_resample_dataset": {
        "inputBbox": [1, 1, 1, 30, 30, 10], "interpMethod": "cubic", "borderSize": [2, 2, 1],
        "zarrFile": True, "cpusPerTask": 3,
    },
    "XR_resave_zarr_wrapper": {
        "blockSize": [64, 64, 32], "shardSize": [128, 128, 64], "bbox": [1, 1, 1, 20, 20, 5],
        "parseCluster": True,
    },
    "XR_tiffToZarr_wrapper": {
        "blockSize": [32, 32, 16], "flippedTile": [True, False], "partialFile": True,
        "InputBbox": [1, 1, 1, 8, 8, 4], "channelPatterns": ["ch0", "ch1"],
    },
    "XR_zarrToTiff_wrapper": {"channelPatterns": ["ch0"], "parseCluster": True},
}  # fmt: skip


def pypetakit5d_calls() -> list[tuple[str, str, list[str]]]:
    """(case name, function, mccMaster's arguments) for PyPetaKit5D's Phase 1 wrappers."""
    out = []
    for function, w in wrappers().items():
        pos = _PYP_POSITIONAL.get(function, [PATHS[:1]])
        out.append(("defaults", function, pypetakit5d_master_argv(function, pos)))
        names = {o[0] for o in w["options"]}
        opts = {k: v for k, v in _PYP_OPTIONS.get(function, {}).items() if k in names}
        if function == "XR_crop_dataset":
            pos = [PATHS, [5, 6, 7, 50, 60, 70]]
        else:
            pos = [
                PATHS if p[0] == "dataPaths" else v
                for p, v in zip(w["positional"], pos, strict=True)
            ]
        out.append(("options", function, pypetakit5d_master_argv(function, pos, **opts)))
    return out


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", s)
