"""The OMW back projector (petakit5d.decon.backprojector) against MATLAB's
omw_backprojector_generation and omw_backprojector_files (the oracle tree's, with
omw_backprojector_mex).

Oracle cases (inputs made in MATLAB's form; the roles give the Python calls the port's form):

- bp.<case>.*: DeconBackProjectorTest's nine PSF cases (Gaussian, elongated, 3 and 5 OTF lobes,
  sheared with the skew detected, on grids that aren't cubic) with their alpha, skew and
  parameters: Gaussian tapers, degenerate bounds (MATLAB's code without the mex's window),
  bounds beyond the edge;
- bp.area.*: areaThresholdOfAComponentSizeMatches (OTFAreaThresh the size of the smallest
  lobe, computed on each side as the MATLAB test does);
- bp.var.*: the other branches: the skew detected or not, two labels (one lobe called skewed),
  a 2D PSF (skewed or not), double and uint16 PSFs, bounds the mex doesn't take (l = 0, l > u),
  minIntThrsh above the threshold, sizes odd and even; bp.err.*: MATLAB's errors (an empty or
  zero PSF, an empty mask, no threshold);
- bp.demo.*: the demo datasets' bead PSFs (copies; psf_gen_new first, as RLdecon) with the
  demos' parameters, skewed, not skewed and detected;
- files.*: omw_backprojector_files (the psfgen files, a back projector that exists, a folder
  with a trailing slash, the names of the files of other bounds and of the default ones, its
  second output); decon.*: decon_lucy_omw_function with the back projector file;
- probe.*: MATLAB's arithmetic of the steps on dense values (abs of complex single, the Wiener
  division and its product with the window, the window function, the fold of the thirds).

Cases with FFTs run on both FFT backends: ``.exact`` (MATLAB's bits; skipped without the
add-on or where the C library gives FFTW other twiddle factors) and ``.permissive``
(decision record 0006). test_thread_counts runs cases again at 1, 3 and 24 threads.
"""

from __future__ import annotations

import dataclasses
import math
import os

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.decon import backprojector as B
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, Data, In, InPath, Tree, Value, Vol, Work, mcell
from petakit5d.testing.compare import Custom, tol, tree
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG, decon_psf, stable_exp, stable_sincos_2pi

BACKENDS = ("exact", "permissive")
PERMISSIVE = tol(
    reason="0006: FFT-based outputs with the permissive backend, single precision (max <= 1e-5 "
    "of the largest value, relative L2 <= 1e-6)",
    rel_max=1e-5,
    rel_l2=1e-6,
)
# the Wiener division amplifies the OTF's last bits where |OTF| is small (by up to 1 / alpha):
# measured relative L2 up to 1.9e-6 (max within 1e-5 of the largest value), masks identical
PERMISSIVE_BP = tol(
    reason="accepted (R02, 2026-11-13): the back projector and its masked Wiener filter with the "
    "permissive backend: max <= 1e-5 of the largest value, relative L2 <= 1e-5 (masks identical)",
    rel_max=1e-5,
    rel_l2=1e-5,
)
BP_PERMISSIVE = {0: PERMISSIVE_BP, 1: PERMISSIVE_BP, 2: PERMISSIVE}
# all-NaN back projectors (degenerate bounds make the window NaN): MATLAB's FFTW and the exact
# backend give some of their NaNs other sign bits (the order of NaN operands in the codelets)
NAN_SIGNS = tol(
    reason="accepted (R02, 2026-11-13): all-NaN back projectors (degenerate bounds): FFTW and the "
    "exact backend propagate NaN signs differently, so only the NaN positions are compared",
    abs=0,
)
SIN_32_45 = 0.5365729180004349  # MATLAB's sind(32.45)
COS_30 = 0.8660254037844387  # MATLAB's cosd(30)
OUT_ROLES = ("vol", "vol", "vol", "vol")


def V(make):
    return Value(make)


# ---------------------------------------------------------------------------------------------
# inputs (Python's form; the Values give MATLAB's)
# ---------------------------------------------------------------------------------------------


def _grid(shape):
    """MATLAB's ndgrid((1 : n) - floor(n / 2) - 1) per axis, as (z, y, x) broadcasts"""
    nz, ny, nx = shape
    z = (np.arange(nz) - nz // 2).astype(np.float64).reshape(-1, 1, 1)
    y = (np.arange(ny) - ny // 2).astype(np.float64).reshape(1, -1, 1)
    x = (np.arange(nx) - nx // 2).astype(np.float64).reshape(1, 1, -1)
    return z, y, x


def _cos2pi(u):
    return stable_sincos_2pi(np.mod(u, 1.0))[1]


def bp_psf(kind: str, shape=(48, 44, 52)) -> np.ndarray:
    """DeconBackProjectorTest.makePSF on grids that aren't cubic (single, (z, y, x))"""
    if kind == "gaussian":
        return decon_psf((33, 31, 29), (2.0, 1.5))
    z, y, x = _grid(shape)
    env = stable_exp(-(y**2 / 2 / 1.5**2 + x**2 / 2 / 1.5**2 + z**2 / 2 / 6**2))
    if kind == "elongated":
        p = stable_exp(-(y**2 / 2 / 1.2**2 + x**2 / 2 / 1.6**2 + z**2 / 2 / 4**2))
    elif kind == "lobes3":
        p = env * (1 + 0.9 * _cos2pi(0.25 * z))
    elif kind == "lobes3sheared":
        p = env * (1 + 0.9 * _cos2pi(0.25 * (z * COS_30 + x * 0.5)))
    elif kind == "lobes5":
        p = env * (1 + 0.6 * _cos2pi(0.2 * z) + 0.6 * _cos2pi(0.4 * z))
    else:
        raise ValueError(kind)
    return np.ascontiguousarray(p.astype(np.float32))


def bp_psf_m(kind: str, shape=(48, 44, 52)) -> np.ndarray:
    return to_matlab(bp_psf(kind, shape))


def psf_2d() -> np.ndarray:
    """a 2D PSF (an image, MATLAB's ny x nx): a Gaussian with a sheared modulation"""
    _, y, x = _grid((1, 45, 51))
    p = stable_exp(-(y**2 / 2 / 1.4**2 + x**2 / 2 / 1.8**2)) * (1 + 0.8 * _cos2pi(0.22 * y))
    return np.ascontiguousarray(p[0].astype(np.float32))


# ---------------------------------------------------------------------------------------------
# the Python side
# ---------------------------------------------------------------------------------------------

_NAMES = {
    "otfcumthresh": "otf_cum_thresh",
    "otfareathresh": "otf_area_thresh",
    "hannwinbounds": "hann_win_bounds",
    "minintthrsh": "min_int_thrsh",
}


def _options(opts) -> dict:
    return {_NAMES[str(k).lower()]: v for k, v in zip(opts[::2], opts[1::2], strict=True)}


def py_bp(psf, alpha, skewed, *opts, backend=None, threads=None):
    r = B.generate_backprojector(
        psf, alpha, skewed, backend=backend, threads=threads, **_options(opts)
    )
    return tuple(r)


def _check(case, tmp_path):
    if case.id.endswith(".exact"):
        reason = digests_hold_here()
        if reason:
            pytest.skip(reason)
    oracle.check(case, tmp_path)


def bp_cases(cid, make, alpha, skewed, opts=(), *, error=None, marks=(), compare=None):
    """omw_backprojector_generation(psf, alpha, skewed, opts{:}): one case per backend"""
    out = []
    # degenerate bounds: the window and the back projector are NaN
    nan = any(isinstance(v, list) and len(v) == 2 and v[0] == v[1] for v in opts)
    for be in BACKENDS if not error else (None,):
        spec = None
        if be == "exact" and nan:
            spec = {0: NAN_SIGNS}
        elif be == "permissive":
            spec = compare if compare is not None else BP_PERMISSIVE
        out.append(
            Case(
                f"{cid}.{be}" if be else cid,
                matlab="omw_backprojector_generation",
                args=[Vol(In("psf")), alpha, skewed, *opts],
                nargout=4,
                outputs=["b_omw", "OTF_bp_omw", "abs_OTF_c", "OTF_mask"],
                python=py_bp,
                py_kwargs={"backend": be} if be else {},
                inputs={"psf": V(make)},
                out_roles=OUT_ROLES,
                expect_error=error,
                py_error=B.BackProjectorError if error else Exception,
                marks=marks,
                **({"compare": spec} if spec is not None else {}),
            )
        )
    return out


CASES: list[Case] = []

# DeconBackProjectorTest's nine cases: (PSF, alpha, skewed, parameters)
BP = {
    "gaussian": ("gaussian", 0.005, False, ()),
    "gaussianSkewed": ("gaussian", 0.005, True, ("OTFCumThresh", 0.8)),
    "elongatedWideTaper": (
        "elongated",
        0.02,
        True,
        ("OTFCumThresh", 0.65, "hannWinBounds", [0.8, 1.2]),
    ),
    "elongatedNarrowTaper": ("elongated", 0.01, False, ("hannWinBounds", [0.95, 1.0])),
    "threeLobes": ("lobes3", 0.005, True, ()),
    "threeLobesShearedDetected": ("lobes3sheared", 0.005, None, ("OTFCumThresh", 0.8)),
    "fiveLobes": ("lobes5", 0.005, True, ("OTFCumThresh", 0.85)),
    "degenerateBounds": ("gaussian", 0.005, False, ("hannWinBounds", [0.8, 0.8])),
    "boundsBeyondTheEdge": ("elongated", 0.005, False, ("hannWinBounds", [0.2, 2.0])),
}
for _k, (_kind, _a, _s, _o) in BP.items():
    CASES += bp_cases(f"bp.{_k}", lambda k=_kind: bp_psf_m(k), _a, _s, _o)

# the other branches
VAR = {
    # the skew detected (3 lobes) and not (1 lobe), and skewed false on 3 lobes
    "lobes3.detected": ("lobes3", 0.005, None, ()),
    "lobes3.notSkewed": ("lobes3", 0.005, False, ()),
    "lobes5.detected": ("lobes5", 0.005, None, ("OTFCumThresh", 0.85)),
    "elongated.detected": ("elongated", 0.005, None, ()),
    # bounds the mex's window doesn't take: MATLAB's code for every voxel
    "lobes3.lowerBoundZero": ("lobes3", 0.005, True, ("hannWinBounds", [0.0, 1.0])),
    "elongated.boundsReversed": ("elongated", 0.005, False, ("hannWinBounds", [1.0, 0.8])),
    # the lower bound of the threshold acts
    "elongated.minIntThrsh": ("elongated", 0.005, False, ("minIntThrsh", 0.02)),
    "lobes3.areaThresh": ("lobes3", 0.005, True, ("OTFAreaThresh", 10)),
    "lobes3.odd": ("lobes3.odd", 0.005, True, ("OTFCumThresh", 0.85)),
}


def var_psf(kind):
    if kind == "lobes3.odd":
        return bp_psf_m("lobes3", (47, 41, 45))
    return bp_psf_m(kind)


for _k, (_kind, _a, _s, _o) in VAR.items():
    CASES += bp_cases(f"bp.var.{_k}", lambda k=_kind: var_psf(k), _a, _s, _o)
for _s, _n in ((False, "notSkewed"), (True, "skewed"), (None, "detected")):
    CASES += bp_cases(f"bp.var.2d.{_n}", psf_2d, 0.005, _s, ())
CASES += bp_cases("bp.var.double", lambda: bp_psf_m("lobes3").astype(np.float64), 0.005, True, ())
CASES += bp_cases(
    "bp.var.uint16",
    lambda: np.floor(bp_psf_m("elongated") * 1000 + 0.5).astype(np.uint16),
    0.01,
    False,
    (),
)

# MATLAB's errors
CASES += bp_cases("bp.err.zero", lambda: np.zeros((8, 9, 7), np.float32), 0.005, False, error=True)
CASES += bp_cases(
    "bp.err.emptyMask",
    lambda: bp_psf_m("elongated"),
    0.005,
    False,
    ("minIntThrsh", 2.0),
    error=True,
)
CASES += bp_cases(
    "bp.err.noThreshold",
    lambda: bp_psf_m("elongated"),
    0.005,
    False,
    ("OTFCumThresh", 1.5),
    error="MATLAB:sizeDimensionsMustMatch",
)

# areaThresholdOfAComponentSizeMatches: OTFAreaThresh is the smallest lobe's size after the
# second opening (the threshold acts only in the second area opening)
_AREA_CODE = (
    "psf = in{1}; psfn = single(psf) ./ sum(single(psf), 'all', 'double'); "
    "a = abs(decon_psf2otf(psfn)); v = sort(double(a(:)), 'descend'); "
    "t = max(v(find(cumsum(v) > sum(v) * 0.9, 1, 'first')) / a(1), 1e-3); "
    "M1 = imopen(fftshift(a) > a(1) * t, strel('sphere', 2)); "
    "M2 = imopen(imclose(M1, strel('sphere', 2)), strel('sphere', 2)); "
    "s2 = cellfun(@numel, getfield(bwconncomp(M2, 26), 'PixelIdxList')); "
    "p = min(s2); out = cell(1, 5); "
    "[out{1 : 4}] = omw_backprojector_generation(psf, 0.005, true, 'OTFAreaThresh', p); "
    "out{5} = p;"
)


def py_area(psf, *, backend=None, threads=None):
    """the MATLAB code's steps with the port's functions, then the back projector"""
    from petakit5d.decon.otf import psf2otf
    from petakit5d.kernels import decon as kd
    from petakit5d.matlabcompat.label import bwconncomp
    from petakit5d.matlabcompat.morphology import imclose, imopen, strel

    p = B._single(psf, 1)
    a = B._abs(psf2otf(B._divide_single(p, kd.sum_double(p), 1), backend=backend), 1)
    dc = a.reshape(-1)[0]
    q = B._otf_threshold(a, dc, 0.9, 1e-3, 1, False)
    se = strel("sphere", 2)
    m1 = imopen(np.fft.fftshift(a) > dc * q, se)
    m2 = imopen(imclose(m1, se), se)
    p_area = float(min(c.size for c in bwconncomp(m2, 26).pixel_idx_list))
    r = B.generate_backprojector(
        psf, 0.005, True, otf_area_thresh=p_area, backend=backend, threads=threads
    )
    return (*r, p_area)


for _be in BACKENDS:
    CASES.append(
        Case(
            f"bp.area.lobes3.{_be}",
            code=_AREA_CODE,
            args=[Vol(In("psf"))],
            nargout=5,
            python=py_area,
            py_kwargs={"backend": _be},
            inputs={"psf": V(lambda: bp_psf_m("lobes3"))},
            out_roles=(*OUT_ROLES, None),
            **({"compare": BP_PERMISSIVE} if _be == "permissive" else {}),
        )
    )

# the demo datasets' PSFs: psf_gen_new (no resampling: the demos' dz equals dzPSF), then the
# demos' back projector parameters
DEMO = {
    "rw488": ("PetaKit5D_demo_cell_image_dataset/PSF/488_2_c.tif", 0.3, 0.005, 0.87, (0.8, 1.0)),
    "rw560": ("PetaKit5D_demo_cell_image_dataset/PSF/560_c.tif", 0.3, 0.05, 0.88, (0.8, 1.0)),
    "widefield": (
        "PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/Widefield/WF_PSF.tif",
        0.1,
        0.02,
        0.65,
        (0.8, 1.2),
    ),
    "confocal": (
        "PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/Confocal/Confocal_PSF.tif",
        0.1,
        0.004,
        0.875,
        (0.8, 1.0),
    ),
}
_DEMO_CODE = (
    "psf = psf_gen_new(single(readtiff(in{1})), in{2}, in{2}, 1.5, 'masked'); out = cell(1, 4); "
    "[out{1 : 4}] = omw_backprojector_generation(psf, in{3}, in{4}, 'OTFCumThresh', in{5}, "
    "'hannWinBounds', in{6});"
)


def py_demo(path, dz, alpha, skewed, thresh, bounds, *, backend=None, threads=None):
    from petakit5d.decon.psf import generate_psf
    from petakit5d.io import read_tiff

    psf = generate_psf(np.asarray(read_tiff(path), np.float32), dz, dz, 1.5, "masked")
    r = B.generate_backprojector(
        psf,
        alpha,
        skewed,
        otf_cum_thresh=thresh,
        hann_win_bounds=bounds,
        backend=backend,
        threads=threads,
    )
    return tuple(r)


for _k, (_rel, _dz, _a, _t, _b) in DEMO.items():
    for _s, _n in ((True, "skewed"), (False, "notSkewed"), (None, "detected")):
        for _be in BACKENDS:
            CASES.append(
                Case(
                    f"bp.demo.{_k}.{_n}.{_be}",
                    code=_DEMO_CODE,
                    args=[InPath("psf"), _dz, _a, _s, _t, list(_b)],
                    nargout=4,
                    python=py_demo,
                    py_kwargs={"backend": _be},
                    inputs={"psf": Data("$PETAKIT5D_DEMO_DATA/" + _rel, "psf.tif", content=True)},
                    out_roles=OUT_ROLES,
                    files=False,
                    marks=("slow",),
                    **({"compare": BP_PERMISSIVE} if _be == "permissive" else {}),
                )
            )

# omw_backprojector_files: the psfgen files; decon_lucy_omw_function with the file
FILES_TREE = tree(ignore=[r".*_figure\.png$"])
FILES_PERMISSIVE = tree(
    ignore=[r".*_figure\.png$"],
    per_file={"*_back_projector_*": PERMISSIVE_BP, "*_OTF.tif": PERMISSIVE},
    tiff_bytes=False,
)


def py_files(psf, folder, name, alpha, thresh, skewed, bounds, psf_fn, save_figure=True, **kw):
    import os

    os.makedirs(folder.rstrip("/"), exist_ok=True)
    return B.omw_backprojector_files(
        psf, folder, name, alpha, thresh, skewed, bounds, psf_fn, save_figure, **kw
    ).b_omw


_FILES_CODE = (
    "if ~exist(in{2}, 'dir'), mkdir(in{2}); end; out = {omw_backprojector_files(in{1}, in{2}, "
    "in{3}, in{4}, in{5}, in{6}, in{7}, in{8}, in{9})};"
)
FILES = {
    # (PSF, folder, alpha, threshold, skewed, bounds, figure)
    "lobes3": ("lobes3", "psfgen", 0.005, 0.88, True, [0.8, 1.0], False),
    "elongated.slash": ("elongated", "psfgen/", 0.0123456789, 0.9, False, [0.8, 1.2], False),
    "lobes3.figure": ("lobes3", "psfgen", 0.005, 0.9, None, [0.8, 1.0], True),
}
for _k, (_kind, _sub, _a, _t, _s, _b, _f) in FILES.items():
    for _be in BACKENDS:
        CASES.append(
            Case(
                f"files.{_k}.{_be}",
                code=_FILES_CODE,
                args=[Vol(In("psf")), Work(_sub), "psf_" + _kind, _a, _t, _s, _b, "psf.tif", _f],
                python=py_files,
                py_kwargs={"backend": _be},
                inputs={"psf": V(lambda k=_kind: bp_psf_m(k))},
                out_roles=("vol",),
                tree=FILES_TREE if _be == "exact" else FILES_PERMISSIVE,
                **({"compare": PERMISSIVE_BP} if _be == "permissive" else {}),
            )
        )


def existing_bp(folder) -> None:
    """a psfgen folder with a back projector file (any values: it is read as it is)"""
    from petakit5d.io import write_tiff

    v = CounterRNG(9, "bp_file").uniform((21, 19, 23), np.float32)
    write_tiff(v, str(folder / "psf_x_back_projector_alpha_0.005000_otf_cum_thresh_0.900000.tif"))


# the names of the files (the second output): the bounds in the back projector's name unless
# they are the default in 6 decimals (DeconFunctionTest/backProjectorFilesOfOtherBoundsAreTheirOwn)
_NAMES_CODE = (
    "mkdir(in{2}); [b, f] = omw_backprojector_files(in{1}, in{2}, 'psf', 0.005, 0.9, false, "
    "in{3}, 'psf.tif', false); [b2, f2] = omw_backprojector_files(in{1}, in{2}, 'psf', 0.005, "
    "0.9, false, in{4}, 'psf.tif', false); out = {b, f, b2, f2};"
)


def py_names(psf, folder, bounds, bounds2, *, backend=None):
    import os

    os.makedirs(folder, exist_ok=True)
    kw = {"backend": backend}
    r = B.omw_backprojector_files(
        psf, folder, "psf", 0.005, 0.9, False, bounds, "psf.tif", False, **kw
    )
    r2 = B.omw_backprojector_files(
        psf, folder, "psf", 0.005, 0.9, False, bounds2, "psf.tif", False, **kw
    )
    return [r.b_omw, r.fsns, r2.b_omw, r2.fsns]


for _k, (_b1, _b2) in {
    "other": ([0.8, 1.0], [0.5, 0.9]),
    "near_default": ([0.8, 1.0], [0.8, 1.0 + 2.0**-52]),
    "other_first": ([0.5, 0.9], [0.8, 1.0]),
}.items():
    CASES.append(
        Case(
            f"files.names.{_k}.exact",
            code=_NAMES_CODE,
            args=[Vol(In("psf")), Work("psfgen"), _b1, _b2],
            nargout=4,
            python=py_names,
            py_kwargs={"backend": "exact"},
            inputs={"psf": V(lambda: bp_psf_m("elongated"))},
            out_roles=("vol", None, "vol", None),
            tree=FILES_TREE,
        )
    )

CASES.append(
    Case(
        "files.existing",
        code=(
            "out = {omw_backprojector_files(in{1}, in{2}, 'psf_x', 0.005, 0.9, true, [0.8, 1], "
            "'psf.tif', false)};"
        ),
        args=[Vol(In("psf")), InPath("psfgen")],
        python=lambda psf, folder: (
            B.omw_backprojector_files(
                psf, folder, "psf_x", 0.005, 0.9, True, (0.8, 1.0), "psf.tif", False
            ).b_omw
        ),
        inputs={"psf": V(lambda: bp_psf_m("lobes3")), "psfgen": Tree(existing_bp, "psfgen")},
        out_roles=("vol",),
    )
)

_DECON_CODE = (
    "mkdir(in{3}); b = omw_backprojector_files(in{2}, in{3}, 'psf', 0.005, 0.88, true, "
    "[0.8, 1], 'psf.tif', false); "
    "bp = readtiff([in{3}, '/psf_back_projector_alpha_0.005000_otf_cum_thresh_0.880000.tif']); "
    "out = {decon_lucy_omw_function(in{1}, in{2}, bp, 2, 'Background', 100, 'useGPU', false), "
    "isequal(b, bp)};"
)
U16_PERMISSIVE = tol(
    reason="0006: uint16 deconvolution with the permissive FFT backend (at most 1 count in at "
    "most 0.1% of the voxels)",
    abs=1,
    max_frac=1e-3,
)


def decon_inputs():
    """beads blurred along z, on a background of 100 with noise (uint16, MATLAB's form)"""
    from petakit5d.testing.synth import beads

    shape = (48, 44, 52)
    pos = [[12, 10, 14], [24, 22, 26], [35, 30, 18], [20, 33, 40]]
    im = beads(shape, pos, (2.5, 1.5), 1500.0, 100.0)
    im = im + np.float32(4) * CounterRNG(11, "decon_noise").normal(shape, np.float32)
    return to_matlab(np.floor(np.maximum(im, 0) + np.float32(0.5)).astype(np.uint16))


def py_decon(im, psf, folder, *, backend=None):
    import os

    from petakit5d.decon import rl

    os.makedirs(folder, exist_ok=True)
    b = B.omw_backprojector_files(
        psf, folder, "psf", 0.005, 0.88, True, (0.8, 1.0), "psf.tif", False, backend=backend
    ).b_omw
    r = rl.lucy_omw(im, psf, b, 2, background=100, use_gpu=False, backend=backend)
    return r.decon, True


for _be in BACKENDS:
    CASES.append(
        Case(
            f"decon.lobes3.{_be}",
            code=_DECON_CODE,
            args=[Vol(In("im")), Vol(In("psf")), Work("psfgen")],
            nargout=2,
            python=py_decon,
            py_kwargs={"backend": _be},
            inputs={"im": V(decon_inputs), "psf": V(lambda: bp_psf_m("lobes3"))},
            out_roles=("vol", None),
            tree=FILES_TREE if _be == "exact" else FILES_PERMISSIVE,
            **({"compare": {0: U16_PERMISSIVE}} if _be == "permissive" else {}),
        )
    )


# probes: MATLAB's arithmetic of the steps on dense values
def probe_single(seed: int, shape=(331, 647), lo=-40, hi=8) -> np.ndarray:
    """single values with random signs and exponents in [lo, hi], zeros, -0 and subnormals
    (an image)"""
    rng = CounterRNG(seed, "probe")
    m = rng.uniform(shape, np.float32) + np.float32(0.5)
    e = rng.integers(lo, hi, shape)
    s = np.where(rng.uniform(shape, np.float32) < np.float32(0.5), -1, 1).astype(np.float32)
    v = (s * np.ldexp(m, e.astype(np.int32))).astype(np.float32)
    v[0, :7] = [0.0, -0.0, 1.4e-45, -1.4e-45, 1.2e-38, 3.4e38, -3.4e38]
    v[1, ::5] = 0
    return v


def probe_weights(seed: int, shape=(331, 647)) -> np.ndarray:
    """window values: doubles in [0, 1], 0, 1, 2, 3 and NaN"""
    w = CounterRNG(seed, "weights").uniform(shape)
    w[:, ::7] = 0.0
    w[:, 1::7] = 1.0
    w[:, 2::11] = 2.0
    w[:, 3::13] = 3.0
    w[5, 5] = np.float64("nan")
    return w


def _complex(a, b):
    c = np.empty(a.shape, np.complex64)
    c.real, c.imag = a, b
    return c


def py_abs(a, b):
    return B._abs(_complex(a, b), 4)


def py_wiener(a, b, alpha, m):
    c = _complex(a, b)
    ab = B._abs(c, 4)
    return [
        B._masked_wiener(c, ab, None, float(alpha), 4),
        B._masked_wiener(c, ab, m, float(alpha), 4),
    ]


CASES.append(
    Case(
        "probe.abs",
        code="out = {abs(complex(in{1}, in{2}))};",
        args=[In("a"), In("b")],
        python=py_abs,
        inputs={"a": V(lambda: probe_single(1)), "b": V(lambda: probe_single(2))},
    )
)


def special_pairs():
    """every pair of special single values (zeros, subnormals, the branch bounds of hypotf,
    huge values, Inf, NaN), both signs, as two images"""
    v = np.array(
        [
            0.0,
            1e-45,
            3e-44,
            1e-40,
            1.17e-38,
            2.0**-60,
            2.0**-50,
            2.0**-49,
            1e-20,
            0.3,
            1.0,
            1.5,
            2.0**30,
            2.0**50,
            2.0**51,
            1e30,
            3.4e38,
            np.inf,
            np.float32("nan"),
        ],
        np.float32,
    )
    v = np.concatenate([v, -v])
    a, b = np.meshgrid(v, v, indexing="ij")
    return np.ascontiguousarray(a), np.ascontiguousarray(b)


CASES.append(
    Case(
        "probe.abs.special",
        code="out = {abs(complex(in{1}, in{2}))};",
        args=[In("a"), In("b")],
        python=py_abs,
        inputs={"a": V(lambda: special_pairs()[0]), "b": V(lambda: special_pairs()[1])},
    )
)
for _w, _alpha in (("a005", 0.005), ("a0123", 0.0123456789)):
    CASES.append(
        Case(
            f"probe.wiener.{_w}",
            code=(
                "x = complex(in{1}, in{2}); w = conj(x) ./ (abs(x) .^ 2 + in{3}); "
                "out = {w, in{4} .* w};"
            ),
            args=[In("a"), In("b"), _alpha, In("m")],
            nargout=2,
            python=py_wiener,
            inputs={
                "a": V(lambda: probe_single(3, lo=-30, hi=1)),
                "b": V(lambda: probe_single(4, lo=-30, hi=1)),
                "m": V(lambda: probe_weights(5)),
            },
        )
    )


def special_products():
    """every combination of special real parts, imaginary parts and double factors"""
    s = [0.0, -0.0, 1.0, -1.0, 0.5, 3e-39, np.inf, -np.inf, np.float32("nan")]
    f = [0.0, -0.0, 1.0, -1.0, 0.3, 2.0, 1e-300, np.inf, -np.inf, np.float64("nan")]
    a, b, m = np.meshgrid(
        np.array(s, np.float32), np.array(s, np.float32), np.array(f), indexing="ij"
    )
    return [np.ascontiguousarray(v.reshape(1, -1)) for v in (a, b, m)]


def py_special(a, b, m):
    """complex(a, b) keeps both parts; conj and the window product as the module computes them"""
    c = B.conj_imag(b)
    return [a, b, a, c, _complex(*B.times_complex(m, a, b)), _complex(*B.times_complex(m, a, c))]


def _products_check(ref, act):
    """bit for bit, except where both parts are non-finite, one of them infinite, and the factor
    is +-0 (MATLAB's result depends on the signs there; only a non-finite OTF has such values: a
    known difference)"""
    a, b, m = (v.reshape(-1) for v in special_products())
    skip = ~np.isfinite(a) & ~np.isfinite(b) & (np.isinf(a) | np.isinf(b)) & (m == 0)
    r = np.ascontiguousarray(ref.data).reshape(-1).view(np.uint32).reshape(-1, 2)
    x = (
        np.ascontiguousarray(act.data, dtype=np.complex64)
        .reshape(-1)
        .view(np.uint32)
        .reshape(-1, 2)
    )
    bad = np.flatnonzero((r != x).any(axis=1) & ~skip)
    return [f"{bad.size} products differ, first at {bad[:5].tolist()}"] if bad.size else []


SPECIAL_PRODUCTS = Custom(_products_check, "exact but for both parts non-finite with m = +-0")


CASES.append(
    Case(
        "probe.step.special",
        code=(
            "x = complex(in{1}, in{2}); y = conj(x); "
            "out = {real(x), imag(x), real(y), imag(y), in{3} .* x, in{3} .* y};"
        ),
        args=[In("a"), In("b"), In("m")],
        nargout=6,
        python=py_special,
        compare={4: SPECIAL_PRODUCTS, 5: SPECIAL_PRODUCTS},
        inputs={
            "a": V(lambda: special_products()[0]),
            "b": V(lambda: special_products()[1]),
            "m": V(lambda: special_products()[2]),
        },
    )
)


def probe_x(seed: int) -> np.ndarray:
    """relative distances: [0, 2.5) with the bounds, their neighbors and non-finite values"""
    x = CounterRNG(seed, "x").uniform((257, 389)) * 2.5
    special = [0.0, 0.2, 0.5, 0.8, 0.95, 1.0, 1.2, 2.0, np.inf, -np.inf, np.float64("nan")]
    for i, v in enumerate(special):
        x[0, 3 * i] = v
        x[0, 3 * i + 1] = np.nextafter(v, -np.inf)
        x[0, 3 * i + 2] = np.nextafter(v, np.inf)
    return x


for _lo, _hi in ((0.8, 1.0), (0.8, 1.2), (0.95, 1.0), (0.2, 2.0), (0.8, 0.8), (1.0, 0.8), (0, 1)):
    CASES.append(
        Case(
            f"probe.window.{_lo}_{_hi}",
            code=(
                "x = in{1}; l = in{2}; u = in{3}; "
                "out = {(x < l) + cos(pi * (x - l)/ 2 / (u - l)) .^ 2 .* (x >= l & x < u)};"
            ),
            args=[In("x"), float(_lo), float(_hi)],
            python=B.window_function,
            inputs={"x": V(lambda: probe_x(6))},
        )
    )


def fold_input():
    rng = CounterRNG(7, "fold")
    shape = (3 * 13, 17, 19)
    return to_matlab(rng.uniform(shape) * np.ldexp(1.0, rng.integers(-30, 30, shape)))


CASES.append(
    Case(
        "probe.fold",
        code="m = in{1}; s = size(m); out = {sum(reshape(m, s(1), s(2), s(3) / 3, 3), 4)};",
        args=[Vol(In("m"))],
        python=B._fold,
        inputs={"m": V(fold_input)},
        out_roles=("vol",),
    )
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    _check(case, tmp_path)


_THREADED = [
    c
    for c in CASES
    if c.id.startswith(("bp.threeLobes.", "bp.fiveLobes.", "bp.degenerateBounds.", "bp.var.2d."))
]


@oracle.parametrize(_THREADED)
@pytest.mark.parametrize("threads", [1, 3, 24])
def test_thread_counts(case, threads, tmp_path):
    """backProjectorMatchesWithOneThread: the results don't depend on the thread count
    (MATLAB's references)"""
    variant = dataclasses.replace(case, py_kwargs={**(case.py_kwargs or {}), "threads": threads})
    variant._module_globals = case._module_globals
    _check(variant, tmp_path)


# the golden decon_omw scenario: the oracle tree's back projector and decon of its inputs
GOLDEN = ("backProjector", "decon")
GOLDEN_CASES = {
    f"golden.decon_omw.{what}": Case(
        f"golden.decon_omw.{what}",
        matlab="pkpy_golden",
        args=[what, "decon_omw", mcell(*(("im", "psf") if what == "inputs" else GOLDEN)), Work()],
        nargout=2,
        outputs=list(("im", "psf") if what == "inputs" else GOLDEN),
        helpers=["tests/matlab_map/matlab"],
        files=False,
    )
    for what in ("compute", "inputs")
}


@oracle.parametrize(list(GOLDEN_CASES.values()))
def test_golden_references(case, tmp_path):
    """the references of the golden scenario exist (store or manifest)"""
    oracle.check(case, tmp_path)


@pytest.mark.parametrize("be", BACKENDS)
def test_golden_scenario(be):
    """GoldenOutputTest/outputMatchesReference(scenario=decon_omw), the back projector: the
    port's back projector of the scenario's PSF is the oracle tree's (exact backend: bit for bit;
    permissive: 0006), and the decon with it equals the oracle tree's decon"""
    from petakit5d.decon import rl

    if be == "exact" and digests_hold_here():
        pytest.skip(digests_hold_here())
    if be not in __import__("petakit5d.fft", fromlist=["available"]).available():
        pytest.skip(f"the {be} FFT backend is not built")
    computed = oracle.reference(GOLDEN_CASES["golden.decon_omw.compute"])
    inputs = oracle.reference(GOLDEN_CASES["golden.decon_omw.inputs"])
    if computed.source != "store" or inputs.source != "store":
        pytest.skip("needs the golden store (the manifests hold hashes only)")
    im, psf = (from_matlab(v.data) for v in inputs.outputs)
    bp = B.generate_backprojector(psf, 0.005, False, 0.9, backend=be, verbose=False).b_omw
    ref_bp, ref_decon = (np.asarray(v.data) for v in computed.outputs)
    out = rl.lucy_omw(im, psf, bp, 2, background=100, save_16bit=False, backend=be).decon
    m_bp, m_out = to_matlab(bp), to_matlab(out)
    assert m_bp.dtype == ref_bp.dtype and m_bp.shape == ref_bp.shape
    if be == "exact":
        assert m_bp.tobytes(order="F") == ref_bp.tobytes(order="F")
        assert m_out.tobytes(order="F") == ref_decon.tobytes(order="F")
    else:
        for a, r in ((m_bp, ref_bp), (m_out, ref_decon)):
            d = np.asarray(a, np.float64) - np.asarray(r, np.float64)
            assert np.linalg.norm(d) <= 1e-4 * np.linalg.norm(np.asarray(r, np.float64))


# ---------------------------------------------------------------------------------------------
# DeconBackProjectorTest's other tests (in intent) and the port's own
# ---------------------------------------------------------------------------------------------


def _ellipsoid_mask(shape=(21, 26, 30), bumps=0.1):
    """an ellipsoid with a bumpy surface (DeconBackProjectorTest's window test), (z, y, x)"""
    from petakit5d.matlabcompat.neighbors import cart2sph

    nz, ny, nx = shape
    z = (np.arange(nz) + 1 - (nz + 1) / 2).reshape(-1, 1, 1)
    y = (np.arange(ny) + 1 - (ny + 1) / 2).reshape(1, -1, 1)
    x = (np.arange(nx) + 1 - (nx + 1) / 2).reshape(1, 1, -1)
    z, y, x = np.broadcast_arrays(z, y, x)
    th, ph, rs = cart2sph(x, y, z)
    from petakit5d.matlabcompat import elementary as E

    c, s = E.cos(ph), E.sin(ph)
    ct, st = E.cos(th), E.sin(th)
    lim = (
        7.0
        / np.sqrt(c * c * (ct * ct / 1.2 + st * st / 0.8) + s * s / 2)
        * (1 + bumps * E.cos(5 * th))
    )
    return rs < lim


def _edge(mask):
    from petakit5d.matlabcompat.morphology import imerode, strel

    return mask & ~imerode(mask, strel("sphere", 1))


@pytest.mark.parametrize("bounds", [(0.8, 1.0), (0.8, 1.2), (0.5, 0.6), (0.2, 2.0)])
def test_the_window_without_the_bounds_is_the_same(bounds):
    """withoutTheMexTheMatlabCodeRuns: the window from the mex's bounds and the nearest edge
    points of the other voxels equals MATLAB's computation for every voxel, bit for bit (the
    port runs the latter where the mex's conditions fail)"""
    edge = _edge(_ellipsoid_mask())
    psz = (edge.shape[1], edge.shape[2], edge.shape[0])
    fast = B._window_fast(edge, psz, *bounds, 4)
    assert fast is not None
    full = B._window_full(edge, psz, *bounds, 4)
    assert fast.tobytes() == full.tobytes()


def test_window_bounds_hold_for_sparse_edge_points():
    """windowBoundsHoldForSparseEdgePoints: with few edge points of very different radii, a voxel
    the bounds decide is below l times both radii of the 2 edge points nearest to its direction
    (window 1), or at least u times both (window 0)"""
    from petakit5d.kernels import decon as kd
    from petakit5d.matlabcompat.neighbors import knnsearch

    shape = (36, 41, 38)
    psz = (41, 38, 36)
    t = np.arange(math.prod(shape), dtype=np.int64)
    dx, dy, dz = B._offsets(t, psz)
    yn, xn, zn, rs = B._directions(dx, dy, dz, 4)
    for trial in range(12):
        rng = CounterRNG(100 + trial, "sparse")
        ne = int(rng.integers(3, 300, 1)[0])
        e = rng.normal((ne, 3))
        e = e / np.sqrt((e * e).sum(axis=1, keepdims=True))
        r = 2 + 30 * rng.uniform(ne)
        lo, hi = sorted(0.2 + 1.5 * rng.uniform(2))
        w, u = kd.omw_window(shape, np.ascontiguousarray(e[:, [2, 0, 1]]), r, (lo, hi), threads=3)
        undecided = np.zeros(shape, bool)
        undecided[u] = True
        # MATLAB's linear order of the voxels
        w_t = np.swapaxes(w, 1, 2).reshape(-1)
        dec = ~np.swapaxes(undecided, 1, 2).reshape(-1)
        idx, _ = knnsearch(e, np.stack([yn[dec], xn[dec], zn[dec]], axis=1), 2)
        rlo, rhi = r[idx].min(axis=1), r[idx].max(axis=1)
        one = w_t[dec] == 1
        assert np.all(rs[dec][one] < lo * rlo[one]), trial
        assert np.all(rs[dec][~one] >= hi * rhi[~one]), trial
        assert np.all(w[undecided] == 0)


@pytest.mark.parametrize("threads", [1, 5])
@pytest.mark.parametrize(
    "shape,bounds",
    [((21, 26, 30), (0.8, 1.0)), ((24, 19, 31), (0.5, 0.9)), ((20, 23, 22), (0.2, 2.0))],
)
def test_the_kernels_window_is_the_numpy_window(shape, bounds, threads):
    """the decon kernels' window values (fdlibm in C, the 2 nearest edge points on a grid of
    cells) equal the NumPy computation (matlabcompat's cart2sph and sph2cart, the k-d tree).
    With odd ny and nx and an even nz the voxels on the z axis have the same direction (atan2(0,
    0) is 0), so two edge voxels on it are at distance 0 from the axis voxel of the first plane:
    a zero weight, where MATLAB computes every voxel (both give None)"""
    edge = _edge(_ellipsoid_mask(shape))
    psz = (edge.shape[1], edge.shape[2], edge.shape[0])
    a = B._window_fast(edge, psz, *bounds, threads)
    b = B._window_fast_numpy(edge, psz, *bounds, threads)
    if shape == (24, 19, 31):
        assert a is None and b is None
        return
    assert a is not None and b is not None
    assert a.tobytes() == b.tobytes()


def test_the_kernels_abs_is_hypotf():
    """the decon kernels' abs equals the NumPy hypotf (both MATLAB's, probe.abs*) on dense and
    special values"""
    for x, y in ((probe_single(1), probe_single(2)), special_pairs()):
        c = _complex(x, y)
        assert B._abs(c, 3).tobytes() == B._abs_numpy(c, 3).tobytes()


@pytest.mark.parametrize("psz", [(12, 9, 7), (10, 16, 13), (1, 9, 4)])
def test_directions_from_tables_are_elementwise(psz):
    """the per-plane and per-(hypot, z) tables give cart2sph and sph2cart of every voxel's
    offset, bit for bit"""
    t = np.arange(math.prod(psz), dtype=np.int64)[::-1].copy()
    ref = B._directions(*B._offsets(t, psz), 1)
    out = B._Directions(psz)(t)
    for a, r in zip(out, ref, strict=True):
        assert a.tobytes() == r.tobytes()


def test_knn2_is_knnsearch():
    """the threaded k-d tree search of the window gives matlabcompat's knnsearch (MATLAB's
    order of the squared distances, ties to the lower index), also for duplicate points"""
    from petakit5d.matlabcompat.neighbors import knnsearch

    rng = CounterRNG(21, "knn")
    p = rng.normal((3000, 3))
    p = p / np.sqrt((p * p).sum(axis=1, keepdims=True))
    p[10] = p[20]  # an exact tie
    p[500:510] = p[100]
    q = rng.normal((20000, 3))
    q = q / np.sqrt((q * q).sum(axis=1, keepdims=True))
    q[:50] = p[:50]
    q[50] = p[100]
    for nt in (1, 5):
        i, d = B._knn2(p, q, nt)
        ri, rd = knnsearch(p, q, 2)
        assert np.array_equal(i, ri)
        assert d.tobytes() == rd.tobytes()


def test_knn2_runs_knnsearch_on_the_callers_threads(monkeypatch):
    """where knnsearch searches (few or non-finite points, and the rows with more ties than the
    tree's candidates) it gets the caller's thread count, not every CPU; same result at any"""
    from petakit5d.matlabcompat import neighbors

    real = neighbors.knnsearch
    threads = []

    def spy(*args, **kwargs):
        threads.append(kwargs.get("num_threads"))
        return real(*args, **kwargs)

    monkeypatch.setattr(neighbors, "knnsearch", spy)
    rng = CounterRNG(22, "knn-threads")
    few = rng.normal((40, 3))
    many = rng.normal((3000, 3))
    many[200:400] = many[100]  # 201 equal points: more ties than the tree's 64 candidates
    nonfinite = many.copy()
    nonfinite[7] = np.inf
    q = rng.normal((2000, 3))
    q[:3] = many[100]
    for p in (few, many, nonfinite):
        want = real(p, q, 2)
        for nt in (1, 3):
            threads.clear()
            i, d = B._knn2(p, q, nt)
            assert threads and threads == [nt] * len(threads)
            assert np.array_equal(i, want[0]) and d.tobytes() == want[1].tobytes()


def test_window_function_is_matlabs_expression():
    """(x < l) + cos(pi (x - l) / 2 / (u - l))^2 (l <= x < u), with MATLAB's cos, where x is
    finite (probe.window.* compare it with MATLAB, NaN and Inf too)"""
    from petakit5d.matlabcompat import elementary as E

    x = CounterRNG(3, "w").uniform(5000) * 2.5
    for lo, hi in ((0.8, 1.0), (0.2, 2.0), (1.0, 0.8)):
        c = E.cos(np.pi * (x - lo) / 2 / (hi - lo))
        ref = (x < lo) + c**2 * ((x >= lo) & (x < hi))
        assert B.window_function(x, lo, hi).tobytes() == ref.tobytes()


def _lobes_mask(centers, shape=(30, 21, 23)):
    nz, ny, nx = shape
    z, y, x = np.ogrid[:nz, :ny, :nx]
    m = np.zeros(shape, bool)
    for cz, cy, cx, r in centers:
        m |= (z - cz) ** 2 + (y - cy) ** 2 + (x - cx) ** 2 <= r * r
    return m


def test_skew_detection():
    """skewed PSFs are detected from the y projection: more than one component whose z profile
    covers every slice and peaks away from the middle"""
    # three tilted lobes that together cover every z, the largest at one end
    skew = _lobes_mask([(4, 10, 4, 5), (14, 10, 11, 4), (24, 10, 18, 6)])
    assert B._is_skewed(skew, 2)
    assert not B._is_skewed(_lobes_mask([(15, 10, 11, 6)]), 2)  # one lobe
    # three lobes on one axis leave slices empty between them
    assert not B._is_skewed(_lobes_mask([(4, 10, 11, 3), (15, 10, 11, 3), (26, 10, 11, 3)]), 2)


@pytest.mark.parametrize("peak,skewed", [(5, False), (6, True)])
def test_skew_detection_needs_the_peak_more_than_1_from_the_middle(peak, skewed):
    """abs(pind - (nz + 1) / 2) > 1: a peak exactly one slice from the middle is not skewed"""
    m = np.zeros((9, 5, 20), bool)
    m[:, 2, 0:4] = True  # a bar over every z
    m[peak, 2, 10:16] = True  # a second component, the profile's peak
    assert B._is_skewed(m, 1) is skewed


def test_three_lobes_are_stacked_as_matlab_labels_them():
    """the skewed mask: bwlabeln's labels 3, 2 and 1 along z; one component (labels 0 and 1)
    becomes label 2, the middle third; more than three: the three largest"""
    from petakit5d.matlabcompat.label import bwlabeln

    m = _lobes_mask([(4, 10, 4, 3), (14, 10, 11, 4), (24, 10, 18, 3)])
    lab, n = bwlabeln(m)
    out = B._three_lobes(m, 2)
    nz = m.shape[0]
    assert out.shape == (3 * nz,) + m.shape[1:]
    for k, label in enumerate((3, 2, 1)):
        assert np.array_equal(out[k * nz : (k + 1) * nz], lab == label)
    one = _lobes_mask([(15, 10, 11, 5)])
    out = B._three_lobes(one, 2)
    assert not out[:nz].any() and not out[2 * nz :].any()
    assert np.array_equal(out[nz : 2 * nz], one)
    from petakit5d.kernels import decon as kd

    five = _lobes_mask(
        [(3, 4, 4, 2), (9, 10, 11, 3), (15, 15, 6, 4), (22, 5, 18, 3), (27, 16, 16, 1)]
    )
    assert bwlabeln(five)[1] == 5
    lab, _ = bwlabeln(kd.omw_components(five, "largest", 3))
    out = B._three_lobes(five, 2)
    for k, label in enumerate((3, 2, 1)):
        assert np.array_equal(out[k * nz : (k + 1) * nz], lab == label)


def test_the_running_sum_is_strict_and_carried_over_pieces():
    """find(cumsum(v) > t, 1): strictly above (not at) t, and the running sum continues across the
    pieces of 2^20 values as one sequential cumsum"""
    v = np.ones(8)
    assert B._first_above(v, 2.0) == 2
    assert B._first_above(v, 8.0) == -1
    n = (1 << 20) + 37
    w = np.full(n, 0.5)
    assert B._first_above(w, 0.5 * ((1 << 20) + 9)) == (1 << 20) + 9


def test_the_threshold_is_printed_as_matlab_prints_it(capsys):
    """fprintf('OTF threshold: %0.2d\\n', dc_thresh) of a single value: %e with 2 digits"""
    B.generate_backprojector(bp_psf("gaussian"), 0.005, False)
    assert capsys.readouterr().out.startswith("OTF threshold: ")
    out = capsys.readouterr().out
    assert out == ""
    B.generate_backprojector(bp_psf("gaussian"), 0.005, False, verbose=False)
    assert capsys.readouterr().out == ""


def test_errors_carry_matlabs_identifiers():
    with pytest.raises(B.BackProjectorError) as e:
        B.generate_backprojector(np.zeros((9, 8, 7), np.float32), 0.005, False)
    assert e.value.identifier == "" and "nonzero PSF" in str(e.value)
    with pytest.raises(B.BackProjectorError) as e:
        B.generate_backprojector(bp_psf("elongated"), 0.005, False, 0.9, min_int_thrsh=2.0)
    assert e.value.identifier == "" and "OTF mask is empty" in str(e.value)
    with pytest.raises(B.BackProjectorError) as e:
        B.generate_backprojector(bp_psf("elongated"), 0.005, False, 1.5)
    assert e.value.identifier == "MATLAB:sizeDimensionsMustMatch"


def test_a_volume_of_one_plane_is_an_image():
    p = psf_2d()
    a = B.generate_backprojector(p, 0.005, None, verbose=False)
    b = B.generate_backprojector(p[None], 0.005, None, verbose=False)
    for x, y in zip(a, b, strict=True):
        assert y.shape == (1,) + x.shape
        assert x.tobytes() == y.tobytes()


def test_file_names_are_matlabs():
    n = B.backprojector_file_names("/mnt/projects/psfgen/", "psf_488", 0.0123456789, 0.9)
    stem = "/mnt/projects/psfgen//psf_488"
    assert n["bp"] == f"{stem}_back_projector_alpha_0.012346_otf_cum_thresh_0.900000.tif"
    assert n["figure"] == f"{stem}_back_projector_alpha_0.012346_otf_cum_thresh_0.900000_figure.png"
    assert n["otf"] == f"{stem}_OTF.tif"
    assert n["mask"] == f"{stem}_OTF_mask_otf_cum_thresh_0.900000.tif"
    assert n["fsns"] == [
        "psf_488_back_projector_alpha_0.012346_otf_cum_thresh_0.900000.tif",
        "psf_488_OTF.tif",
        "psf_488_OTF_mask_otf_cum_thresh_0.900000.tif",
    ]


@pytest.mark.parametrize(
    ("bounds", "suffix"),
    [
        ((0.8, 1.0), ""),
        ([0.8, 1], ""),
        ((0.8, 1.0 + 2.0**-52), ""),
        ((0.8000004, 0.9999996), ""),
        ((0.5, 0.9), "_hann_win_bounds_0.500000_0.900000"),
        ((0.8, 1.2), "_hann_win_bounds_0.800000_1.200000"),
        ((0.8000006, 1.0), "_hann_win_bounds_0.800001_1.000000"),
        (np.array([0.8, 1.0], np.float32), ""),
    ],
)
def test_back_projectors_of_other_bounds_have_their_own_name(bounds, suffix):
    """the bounds are in the back projector's name (and its figure's) unless they are the
    default [0.8, 1] in 6 decimals; the OTF and the mask keep their names"""
    n = B.backprojector_file_names("/mnt/projects/psfgen/", "p", 0.005, 0.9, bounds)
    stem = f"/mnt/projects/psfgen//p_back_projector_alpha_0.005000_otf_cum_thresh_0.900000{suffix}"
    assert n["bp"] == stem + ".tif"
    assert n["figure"] == stem + "_figure.png"
    assert n["fsns"][0] == os.path.basename(stem) + ".tif"
    assert n["fsns"][1:] == ["p_OTF.tif", "p_OTF_mask_otf_cum_thresh_0.900000.tif"]


def test_bounds_of_one_value_fail_on_their_second_value():
    with pytest.raises(B.BackProjectorError) as e:
        B.backprojector_file_names("/mnt/projects/psfgen/", "p", 0.005, 0.9, (0.8,))
    assert e.value.identifier == "MATLAB:badsubscript"


def test_a_call_with_other_bounds_makes_its_own_back_projector(tmp_path):
    """the back projector of other bounds isn't the earlier one of the same PSF (it was: its
    name lacked the bounds); each call returns the names of its files"""
    psf = bp_psf("elongated")
    folder = str(tmp_path)
    a = B.omw_backprojector_files(psf, folder, "p", 0.005, 0.9, False, (0.8, 1.0), "p", False)
    b = B.omw_backprojector_files(psf, folder, "p", 0.005, 0.9, False, (0.5, 0.9), "p", False)
    ref = B.generate_backprojector(psf, 0.005, False, 0.9, hann_win_bounds=(0.5, 0.9))
    assert b.b_omw.tobytes() == ref.b_omw.tobytes()
    assert a.b_omw.tobytes() != b.b_omw.tobytes()
    assert a.fsns[0] == "p_back_projector_alpha_0.005000_otf_cum_thresh_0.900000.tif"
    assert b.fsns[0] == (
        "p_back_projector_alpha_0.005000_otf_cum_thresh_0.900000_hann_win_bounds_0.500000_"
        "0.900000.tif"
    )
    assert a.fsns[1:] == b.fsns[1:]
    assert sorted(os.listdir(folder)) == sorted({*a.fsns, *b.fsns})
    # read again by its name
    c = B.omw_backprojector_files(psf, folder, "p", 0.005, 0.9, False, (0.5, 0.9), "p", False)
    assert c.b_omw.tobytes() == b.b_omw.tobytes() and c.fsns == b.fsns


def test_files_are_made_once_and_the_back_projector_last(tmp_path):
    """omw_backprojector_files: the OTF and the mask exist when the back projector is published;
    the next call reads the back projector (cached) and writes nothing"""
    from petakit5d.decon.psf import read_tiff_cached
    from petakit5d.io import read_tiff
    from petakit5d.testing.publish_spy import PublishSpy

    psf = bp_psf("lobes3")
    folder = str(tmp_path)
    names = B.backprojector_file_names(folder, "p", 0.005, 0.88)
    with PublishSpy() as spy:
        b1 = B.omw_backprojector_files(
            psf, folder, "p", 0.005, 0.88, True, (0.8, 1.0), "p.tif", False
        ).b_omw
    files = spy.files_when_moved(names["bp"])
    assert "p_OTF.tif" in files and "p_OTF_mask_otf_cum_thresh_0.880000.tif" in files
    assert sorted(os.listdir(folder)) == sorted(
        os.path.basename(names[k]) for k in ("bp", "otf", "mask")
    )
    r = B.generate_backprojector(psf, 0.005, True, 0.88, verbose=False)
    assert b1.tobytes() == r.b_omw.tobytes()
    assert read_tiff(names["bp"]).tobytes() == r.b_omw.tobytes()
    assert read_tiff(names["otf"]).tobytes() == r.abs_otf_c.tobytes()
    assert read_tiff(names["mask"]).tobytes() == r.otf_mask.astype(np.uint8).tobytes()
    before = {f: os.stat(os.path.join(folder, f)).st_mtime_ns for f in os.listdir(folder)}
    b2 = B.omw_backprojector_files(
        psf, folder, "p", 0.005, 0.88, True, (0.8, 1.0), "p.tif", False
    ).b_omw
    assert b2 is read_tiff_cached(names["bp"])
    assert b2.tobytes() == b1.tobytes()
    assert {f: os.stat(os.path.join(folder, f)).st_mtime_ns for f in os.listdir(folder)} == before


def test_a_back_projector_another_writer_published_is_kept(tmp_path, monkeypatch):
    """a back projector published between the check and the rename is kept, and the temporary
    file removed"""
    from petakit5d.io import read_tiff, write_tiff

    folder = str(tmp_path)
    names = B.backprojector_file_names(folder, "p", 0.005, 0.9)
    other = np.full((5, 6, 7), 3, np.float32)
    real = B.generate_backprojector

    def generate(*a, **k):
        r = real(*a, **k)
        write_tiff(other, names["bp"])
        return r

    monkeypatch.setattr(B, "generate_backprojector", generate)
    B.omw_backprojector_files(
        bp_psf("elongated"), folder, "p", 0.005, 0.9, False, (0.8, 1.0), "p", False
    )
    assert read_tiff(names["bp"]).tobytes() == other.tobytes()
    assert len(os.listdir(folder)) == 3


def test_the_figure_is_drawn_with_matplotlib(tmp_path):
    """the mask's outline on the OTF as a PNG, without effect on the data files"""
    pytest.importorskip("matplotlib")
    r = B.generate_backprojector(bp_psf("lobes3"), 0.005, True, verbose=False)
    p = tmp_path / "f.png"
    assert B._figure(r.abs_otf_c, r.otf_mask, str(p))
    assert p.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_without_matplotlib_no_figure_is_drawn(tmp_path, monkeypatch):
    import builtins

    real = builtins.__import__

    def no_mpl(name, *a, **k):
        if name.startswith("matplotlib"):
            raise ImportError(name)
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_mpl)
    r = B.generate_backprojector(bp_psf("gaussian"), 0.005, False, verbose=False)
    assert not B._figure(r.abs_otf_c, r.otf_mask, str(tmp_path / "f.png"))
    assert not (tmp_path / "f.png").exists()


def test_the_matlab_names_take_and_give_matlabs_arrays():
    """petakit5d.matlab.omw_backprojector_generation: MATLAB's (y, x, z) arrays, [] for the
    detection, name-value pairs in any case"""
    from petakit5d import matlab

    m = bp_psf_m("lobes3")
    out = matlab.omw_backprojector_generation(m, 0.005, [], "otfcumthresh", 0.85)
    ref = B.generate_backprojector(from_matlab(m), 0.005, None, 0.85, verbose=False)
    assert len(out) == 4
    for a, r in zip(out, ref, strict=True):
        assert a.shape == to_matlab(r).shape
        assert np.asarray(a).tobytes(order="F") == to_matlab(r).tobytes(order="F")
