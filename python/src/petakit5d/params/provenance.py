"""Third-party MATLAB code in the registry: what the repository keeps of it.

The port is BSD-3-Clause. A MATLAB file whose header names an author outside the lab, or
another institution's copyright, is third-party code. Of its functions the repository keeps
only the interface, which anyone can observe from outside: the parameter names, their kinds
(required, optional, name-value), their order and their defaults. Their validator texts are
not kept (decided on 2026-11-20): the stored inputParser extract has ``"validator": null``
and ``"withheld": true`` for them and the function's reason in ``third_party``;
check_registry compares the rest, and the mapping table shows "(third-party: not shown)".

The decisions are made per file from its header (the author lines of its first comment
block, nothing further). A lab-written file that says it is based on outside code is the
lab's (decided on 2026-11-12). A function from a folder that holds outside code needs a
decision here before the extract can hold it.
"""

from __future__ import annotations

import dataclasses

from petakit5d.params.matlab_source import MatlabFunction, MatlabParam

_CME = "its header names an outside author (cmeAnalysis)"
_CME_EDITED = f"{_CME}, with lab edits"

# third-party code: folders (ending in /) and files of the MATLAB tree, with the reason
THIRD_PARTY: dict[str, str] = {
    "imageProcessing/": "an outside-author folder",
    "pointDetection/CMEAnalysis_104_XR/": "an outside-author folder (CMEAnalysis)",
    "pointDetection/CMEAnalysis_104_XR/software/loadTrackSettings.m": _CME,
    "pointDetection/CMEAnalysis_104_XR/software/plotIntensityCohorts.m": _CME,
    "pointDetection/llsmtools/": "an outside-author folder (llsmtools)",
    "pointDetection/deskewData.m": "its header names an outside author",
    "pointDetection/XR_loadConditionData3D.m": _CME_EDITED,
    "pointDetection/XR_runDetection3D.m": _CME_EDITED,
    "pointDetection/XR_runDetection3D_frame_parallel.m": _CME_EDITED,
    "pointDetection/XR_runTrackProcessing3D.m": _CME,
    "pointDetection/XR_runTracking3D.m": _CME_EDITED,
    "pointDetection/rotateTracks3D.m": _CME_EDITED,
    "pointDetection/runLifetimeAnalysis3D.m": _CME,
    "third_parties/BaSiC-master/": "an outside-author folder (BaSiC)",
    "third_parties/u-track3D/": "an outside-author folder (u-track)",
    "third_parties/WBDeconvolution/": "an outside-author folder (WBDeconvolution)",
    "third_parties/sim3D/MaskedTranslationRegistration2D_fit.m": "an outside author's file",
    "third_parties/sim3D/imtranslate_function.m": "a modified MathWorks function",
    "third_parties/sim3D/normxcorr2_masked.m": "an outside author's file",
    "third_parties/sim3D/perdecomp_3D.m": "an outside author's file",
}

# the lab's code in those folders (headers that name the lab's authors only), and sim3D
LAB: dict[str, str] = {
    "pointDetection/GU_estimateSigma3D.m": "a lab author",
    "pointDetection/XR_correctXZoffsetData3D.m": "a lab author",
    "pointDetection/XR_runDetTrack3d.m": "lab authors (based on cmeAnalysis)",
    "pointDetection/src/XR_amiraWriteTracks.m": "a lab author (based on outside code)",
    "pointDetection/src/XR_get_camera_standard_deviation_threshold.m": "a lab author",
    "third_parties/sim3D/": "released under BSD (decided on 2026-11-19), but its outside files",
}

# folders that hold outside code: their files need a decision above
CHECKED_FOLDERS = ("imageProcessing/", "pointDetection/", "third_parties/")

WITHHELD = "(third-party: not shown)"


class ProvenanceError(ValueError):
    """A MATLAB file of a folder with outside code, with no decision on whose code it is."""


def _lookup(table: dict[str, str], file: str) -> str | None:
    """the reason of a file's own entry, else of the first folder that holds it"""
    if file in table:
        return table[file]
    for key, why in table.items():
        if key.endswith("/") and file.startswith(key):
            return why
    return None


def third_party(file: str) -> str | None:
    """Why the MATLAB file ``file`` (relative to the tree, with /) is third-party code; None
    for the lab's. ProvenanceError for a file of a folder with outside code that has no
    decision."""
    f = file.replace("\\", "/")
    why = _lookup(THIRD_PARTY, f)
    if why is not None or _lookup(LAB, f) is not None:
        return why
    if f.startswith(CHECKED_FOLDERS):
        raise ProvenanceError(
            f"{file}: no decision on whose code it is: read the author lines of its header, "
            "nothing further, and add it to petakit5d.params.provenance"
        )
    return None


def interface_only(fn: MatlabFunction) -> MatlabFunction:
    """``fn`` as the repository keeps it: for third-party code without the validator texts
    (its own and its mcc parser's), marked; the lab's code unchanged."""
    why = third_party(fn.file)
    if why is None:
        return fn
    params = [
        MatlabParam(p.kind, p.name, p.default, None, p.withheld or p.validator is not None)
        for p in fn.params
    ]
    parser = fn.parser
    if parser and parser.get("validators"):
        parser = {**parser, "validators": dict.fromkeys(parser["validators"])}
    return dataclasses.replace(fn, params=params, parser=parser, third_party=why)
