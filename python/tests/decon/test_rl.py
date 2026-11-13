"""Richardson-Lucy deconvolution on the CPU (petakit5d.decon.rl) against MATLAB's
decon_lucy_function and decon_lucy_omw_function (the oracle tree's, on the CPU: the half-spectrum
mex by default, the MATLAB loop for data the mex doesn't take or with useMex false).

Oracle cases (inputs made in MATLAB's form; the roles give the Python calls the port's form):

- rl.* / omw.*: DeconRFFTTest's 11 data sizes (2D and 3D, even and odd, length 64 along y, x or z,
  64 x 230 and 64 x 210, two rows along y); uint16, single, uint8 and double data; 0, 1, 2, 10 and
  30 iterations (and 15 on a noise-free bead whose acceleration factor is clamped to 1); backgrounds
  with fractions; each option of DeconRFFTTest and DeconFunctionTest (dampFactor, EdgeErosion on
  full data and on data with zero borders, deconOffset, scaleFactor, bbox, save16bit, all of them);
  MATLAB's loop (useMex false), data it runs for (double, a background per voxel, one row along y
  or x); debug output (err_mat, the iteration count, the saved estimates) with DeconRFFTTest's
  steps, also where MATLAB's bookkeeping errors; errors.
- Each case runs twice on the Python side: ``.exact`` on the exact FFT backend (MATLAB's bits;
  skipped without the add-on or where the C library gives FFTW other twiddle factors than the
  reference platform) and ``.permissive`` on the default backend (decision record 0006: relative
  L2 <= 1e-4; uint16 off by at most 1 in at most 0.1% of the voxels).
- debug.*.big: err_mat of 89,000 or more voxels is MATLAB's mean at one thread (fixed_order_stat,
  M05 c3efd96, the head of ``test``); the oracle tree d67e1b7 has MATLAB's multithreaded mean
  there, so their references change at switch #4 (marked in SWITCH4).

test_thread_counts runs cases again at 1, 2, 7 and 24 threads against the same references.
"""

from __future__ import annotations

import dataclasses
import math
import warnings

import numpy as np
import pytest

from petakit5d import fft as P
from petakit5d.core.arrays import from_matlab, size_to_matlab, to_matlab
from petakit5d.decon import rl
from petakit5d.testing import oracle
from petakit5d.testing.cases import Box, Case, In, Value, Vol, Work, mcell
from petakit5d.testing.compare import tol, tree
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG, beads, decon_psf, texture

BACKENDS = P.available()
EXACT = "exact" in BACKENDS
PERMISSIVE_TOL = tol(
    reason="0006: deconvolution with the permissive FFT backend (relative L2 <= 1e-4)", rel_l2=1e-4
)
PERMISSIVE_U16 = tol(
    reason="0006: uint16 deconvolution with the permissive FFT backend (at most 1 count in at "
    "most 0.1% of the voxels)",
    abs=1,
    max_frac=1e-3,
)
# debug cases whose d67e1b7 reference changes at oracle switch #4 (M05's fixed-order err_mat):
# until then err_mat is compared within the last bits of MATLAB's multithreaded mean (measured: at
# most 128 ulp, relative 2.5e-14); exact from switch #4 on (drop SWITCH4_TOL there)
SWITCH4: set[str] = set()
SWITCH4_TOL = tol(
    reason="M05 (c3efd96): the oracle tree d67e1b7's debug err_mat is MATLAB's mean at 8 threads, "
    "the port (and the head of test) MATLAB's mean at one thread; exact after oracle switch #4",
    rel=1e-12,
)


def V(make):
    return Value(make)


# ---------------------------------------------------------------------------------------------
# inputs, in MATLAB's form
# ---------------------------------------------------------------------------------------------

BLUR_SIGMA = (math.sqrt(0.8**2 + 2.0**2), math.sqrt(0.8**2 + 1.5**2))  # (z, yx)


def size_name(s) -> str:
    return "x".join(str(v) for v in s)


def _py_size(msize) -> tuple[int, ...]:
    sz = tuple(int(v) for v in msize)
    return (sz[2], sz[0], sz[1]) if len(sz) == 3 else sz


def blurred_beads(msize, seed: int, cls: str = "uint16") -> np.ndarray:
    """DeconRFFTTest.blurredBeads: beads blurred by the test PSF plus a little noise on a background
    of 100 (uint16), or the same in single plus 0.37; MATLAB's form"""
    py = _py_size(msize)
    rng = CounterRNG(seed, "beads")
    n = max(3, round(math.prod(py) / 4000))
    pos = rng.uniform((n, len(py))) * (np.array(py, dtype=np.float64) - 1)
    im = beads(py, pos, BLUR_SIGMA if len(py) == 3 else BLUR_SIGMA[1], 1000.0, 100.0)
    im = im + np.float32(3) * rng.uniform(py, np.float32)
    u = np.floor(im + np.float32(0.5)).astype(np.uint16)  # im >= 0
    if cls == "uint16":
        out = u
    elif cls == "uint8":
        out = np.floor(u.astype(np.float64) / 8 + 0.5).astype(np.uint8)  # MATLAB's uint8(im / 8)
    elif cls == "double":
        out = u.astype(np.float64)
    else:
        out = (u.astype(np.float32) + np.float32(0.37)).astype(np.float32)
    return to_matlab(out)


def bead(msize) -> np.ndarray:
    """DeconRFFTTest/rlAccelerationIsClampedAsInTheMatlabLoop: one noise-free bright bead, whose
    acceleration factor goes above 1; MATLAB's form"""
    py = _py_size(msize)
    im = beads(py, [[17.0, 23.0, 22.0]], (2.0, 1.5), 20000.0, 100.0)
    return to_matlab(np.floor(im + np.float32(0.5)).astype(np.uint16))


def zero_border(msize, seed: int) -> np.ndarray:
    """blurred beads with zero triangles along x-z as skewed data has (the mask of nonzero voxels
    isn't full, so EdgeErosion runs the mex's erosion); MATLAB's form"""
    a = from_matlab(blurred_beads(msize, seed))
    nz, _, nx = a.shape
    z = np.arange(1, nz + 1, dtype=np.float64).reshape(-1, 1, 1)
    x = np.arange(1, nx + 1, dtype=np.float64).reshape(1, 1, -1)
    keep = (x > 0.6 * z + nx * 0.05) & (x < 0.6 * z + nx * 0.95)
    return to_matlab(np.where(keep, a, np.uint16(0)))


def textured(msize) -> np.ndarray:
    """DeconFunctionTest/edgeErosionZeroesTheBorderShell: a texture without zero voxels (uint16)"""
    t = texture(_py_size(msize), 3.0, (200.0, 2000.0))
    return to_matlab(np.floor(t + np.float32(0.5)).astype(np.uint16))


def make_psf(msize, which: str = "f") -> np.ndarray:
    """the PSF (or a narrower one as the back projector) cropped to the data as RLdecon crops, a
    plane for 2D data; MATLAB's form"""
    sigma = (2.0, 1.5) if which == "f" else (1.2, 1.0)
    psf = decon_psf((31, 31, 31), sigma)
    if which == "neg":
        # a back projector with negative side lobes, as OMW back projectors have
        wide = decon_psf((31, 31, 31), (3.0, 2.5)) * np.float32(0.25)
        psf = (psf * np.float32(2.0) - wide).astype(np.float32)
    py = _py_size(msize)
    if len(py) == 2:
        psf = psf[15]
    s = [max(0, (p - n) // 2) for p, n in zip(psf.shape, py, strict=True)]
    t = [a + min(p, n) for a, p, n in zip(s, psf.shape, py, strict=True)]
    psf = psf[tuple(slice(a, b) for a, b in zip(s, t, strict=True))]
    return to_matlab(np.ascontiguousarray(psf))


# ---------------------------------------------------------------------------------------------
# Python sides (MATLAB's arguments mapped by the roles)
# ---------------------------------------------------------------------------------------------


def _iterations_out(r):
    return float(r.iterations) if r.iterations > 0 else np.zeros((0, 0))


def py_rl(im, psf, numit, *, backend, nargout=1, threads=None, **kw):
    r = rl.lucy(im, psf, numit, backend=backend, threads=threads, use_gpu=False, **kw)
    out = (r.decon, r.err_mat, _iterations_out(r))
    return out[0] if nargout == 1 else out[:nargout]


def py_omw(im, psf, psfb, numit, *, backend, nargout=1, threads=None, **kw):
    r = rl.lucy_omw(im, psf, psfb, numit, backend=backend, threads=threads, use_gpu=False, **kw)
    out = (r.decon, r.err_mat)
    return out[0] if nargout == 1 else out[:nargout]


# MATLAB's name-value options and the Python keywords
NAMES = {
    "Background": "background",
    "save16bit": "save_16bit",
    "dampFactor": "damp_factor",
    "scaleFactor": "scale_factor",
    "deconOffset": "decon_offset",
    "EdgeErosion": "edge_erosion",
    "bbox": "bbox",
    "debug": "debug",
    "debug_folder": "debug_folder",
    "saveStep": "save_step",
    "useMex": "use_mex",
}


def _u16(opts: dict) -> bool:
    return bool(opts.get("save16bit", True))


def decon_cases(
    cid: str,
    method: str,
    inputs: dict,
    numit,
    opts: dict,
    *,
    nargout: int = 1,
    files: bool | None = None,
    expect_error: bool | str = False,
    marks=(),
    backends=("exact", "permissive"),
    bg_role: bool = False,
    switch4: bool = False,
) -> list[Case]:
    """the case of a decon call with MATLAB's options (name: value), on each backend; switch4:
    err_mat's reference (d67e1b7) changes at oracle switch #4"""
    vols = [Vol(In("im")), Vol(In("psf"))] + ([Vol(In("psfb"))] if method == "omw" else [])
    args = [*vols, numit, "useGPU", False]
    py_kwargs = {"nargout": nargout}
    for k, v in opts.items():
        if k == "Background" and bg_role:
            args += [k, Vol(In("bg"))]
            py_kwargs[NAMES[k]] = Vol(In("bg"))
            continue
        args += [k, Box(v) if k == "bbox" else v]
        py_kwargs[NAMES[k]] = Box(v) if k == "bbox" else v
    out_roles = ("vol", None, None)[:nargout]
    out = []
    for be in backends:
        trees = {}
        if expect_error:
            spec = None
        elif be == "exact":
            spec = {1: SWITCH4_TOL} if switch4 else None
        else:
            first = PERMISSIVE_U16 if _u16(opts) else PERMISSIVE_TOL
            spec = {0: first, **dict.fromkeys(range(1, nargout), PERMISSIVE_TOL)}
        if files and be == "permissive":
            # the saved estimates (uint16) on the permissive backend
            trees = {"tree": tree(per_file={"*.tif": PERMISSIVE_U16}, tiff_bytes=False)}
        kw = dict(py_kwargs, backend=be)
        c = Case(
            f"{cid}.{be}" if len(backends) > 1 else cid,
            matlab="decon_lucy_function" if method == "simplified" else "decon_lucy_omw_function",
            args=args,
            nargout=nargout,
            python=py_rl if method == "simplified" else py_omw,
            py_args=[*vols, numit],
            py_kwargs=kw,
            inputs=inputs,
            out_roles=out_roles,
            files=files,
            expect_error=expect_error,
            py_error=(rl.DeconError, ValueError, TypeError) if expect_error else Exception,
            marks=marks,
            **({"compare": spec} if spec is not None else {}),
            **trees,
        )
        out.append(c)
    return out


def _inputs(msize, seed, cls="uint16", data=None):
    inp = {
        "im": V(
            data if data is not None else (lambda m=msize, s=seed, c=cls: blurred_beads(m, s, c))
        ),
        "psf": V(lambda m=msize: make_psf(m, "f")),
        "psfb": V(lambda m=msize: make_psf(m, "b")),
    }
    return inp


CASES: list[Case] = []
METHODS = {"rl": ("simplified", 10), "omw": ("omw", 2)}

# DeconRFFTTest's sizes (mexEqualsMatlabLoopForOtherSizes, halfOTF's sizes)
SIZES = [
    [48, 46, 36],
    [47, 45, 33],
    [66, 64, 33],
    [64, 50, 40],
    [10, 23, 64],
    [64, 50],
    [30, 64],
    [64, 230],
    [64, 210],
    [10, 64, 23],
    [2, 64, 33],
]
_seed = 0
for _ms in SIZES:
    _seed += 1
    for _m, (_method, _n) in METHODS.items():
        CASES += decon_cases(
            f"{_m}.size.{size_name(_ms)}",
            _method,
            _inputs(_ms, _seed),
            _n,
            {"Background": 100, "save16bit": False},
        )
        CASES += decon_cases(
            f"{_m}.size.{size_name(_ms)}.uint16",
            _method,
            _inputs(_ms, _seed),
            _n,
            {"Background": 100},
        )

# data classes (mexEqualsMatlabLoop's uint16, single, uint8; DeconFunctionTest's double)
_MS = [48, 46, 36]
for _cls, _bg in [("single", 100.3), ("uint8", 12.5), ("double", 100)]:
    for _m, (_method, _n) in METHODS.items():
        CASES += decon_cases(
            f"{_m}.class.{_cls}",
            _method,
            _inputs(_MS, 40, _cls),
            _n,
            {"Background": _bg, "save16bit": False},
        )

# iterations 0, 1, 2, 30 (and the default background); RL's clamped acceleration factor
for _cls in ("uint16", "single"):
    for _m, (_method, _) in METHODS.items():
        for _n in (0, 1, 2, 30):
            CASES += decon_cases(
                f"{_m}.iter{_n}.{_cls}",
                _method,
                _inputs(_MS, 41, _cls),
                _n,
                {"Background": 100.3, "save16bit": False},
                nargout=3 if _method == "simplified" else 2,
                marks=("slow",) if _n == 30 else (),
            )
for _m, (_method, _n) in METHODS.items():
    CASES += decon_cases(f"{_m}.nobackground", _method, _inputs(_MS, 42), _n, {"save16bit": False})
CASES += decon_cases(
    "rl.clamped",
    "simplified",
    {"im": V(lambda: bead(_MS)), "psf": V(lambda: make_psf(_MS, "f"))},
    15,
    {"Background": 100, "save16bit": False},
)

# options (DeconRFFTTest/optionsAreAppliedAsInTheMatlabLoop)
OPTIONS = {
    "damp": {"dampFactor": 1.2},
    "erosion": {"EdgeErosion": 3},
    "offset": {"deconOffset": 5},
    "scale": {"scaleFactor": 2},
    "bbox": {"bbox": [3, 4, 5, 40, 41, 30]},
    "u16": {"save16bit": True},
    "bg": {"Background": 100.3},
    "all": {
        "EdgeErosion": 2,
        "scaleFactor": 1.5,
        "deconOffset": 2,
        "bbox": [2, 2, 2, 44, 40, 33],
        "save16bit": True,
        "dampFactor": 1.3,
    },
}
for _on, _o in OPTIONS.items():
    for _m, (_method, _n) in METHODS.items():
        _opts = {"Background": 100, "save16bit": False, **_o}
        CASES += decon_cases(f"{_m}.opt.{_on}", _method, _inputs(_MS, 43), _n, _opts)
    # RL with single input (the background subtraction in single)
    CASES += decon_cases(
        f"rl.opt.{_on}.single",
        "simplified",
        _inputs(_MS, 43, "single"),
        10,
        {"Background": 100, "save16bit": False, **_o},
    )
# DeconFunctionTest/loopMatchesThePreviousFormulation's options for each class
DFT_OPTIONS = {
    "uint16Output": {"save16bit": True},
    "singleOutput": {"save16bit": False},
    "dampErosionOffsetScale": {
        "save16bit": False,
        "dampFactor": 1.2,
        "EdgeErosion": 4,
        "deconOffset": 3,
        "scaleFactor": 1.5,
    },
    "bbox": {"save16bit": True, "bbox": [3, 4, 2, 44, 45, 36]},
}
_DFT = [48, 48, 40]
for _cls in ("uint16", "single", "double"):
    for _on, _o in DFT_OPTIONS.items():
        for _m, (_method, _n) in (("rl", ("simplified", 4)), ("omw", ("omw", 2))):
            CASES += decon_cases(
                f"{_m}.dft.{_on}.{_cls}",
                _method,
                _inputs(_DFT, 44, _cls),
                _n,
                {"Background": 100, **_o},
            )
# EdgeErosion on data with zero borders (the mex's erosion), and on a texture without zeros
for _m, (_method, _n) in METHODS.items():
    CASES += decon_cases(
        f"{_m}.opt.erosion.zero_border",
        _method,
        {**_inputs(_MS, 45), "im": V(lambda: zero_border(_MS, 45))},
        _n,
        {"Background": 100, "save16bit": False, "EdgeErosion": 4, "deconOffset": 1},
    )
# a back projector with negative lobes (products with zeros of the data: -0 against +0)
CASES += decon_cases(
    "omw.negbp",
    "omw",
    {**_inputs(_MS, 47), "psfb": V(lambda: make_psf(_MS, "neg"))},
    2,
    {"Background": 100, "save16bit": False},
)

# MATLAB's loop: useMex false, and data the mex doesn't take
for _m, (_method, _n) in METHODS.items():
    for _cls in ("uint16", "single"):
        CASES += decon_cases(
            f"{_m}.loop.{_cls}",
            _method,
            _inputs(_MS, 48, _cls),
            _n,
            {"Background": 100, "save16bit": False, "useMex": False},
        )
    CASES += decon_cases(
        f"{_m}.loop.opts",
        _method,
        _inputs(_MS, 48),
        _n,
        {"Background": 100.3, "useMex": False, **OPTIONS["all"]},
    )
    for _rn, _rs in (("row", [1, 64, 33]), ("column", [64, 1, 33])):
        CASES += decon_cases(
            f"{_m}.loop.{_rn}",
            _method,
            _inputs(_rs, 49),
            _n,
            {"Background": 100, "save16bit": False},
        )
    CASES += decon_cases(
        f"{_m}.loop.voxel_background",
        _method,
        {**_inputs(_MS, 50), "bg": V(lambda: to_matlab(np.full(_py_size(_MS), 100, np.float32)))},
        _n,
        {"Background": None, "save16bit": False},
        bg_role=True,
    )
    CASES += decon_cases(
        f"{_m}.error.mex_double",
        _method,
        _inputs(_MS, 51, "double"),
        _n,
        {"Background": 100, "useMex": True},
        expect_error=True,
        backends=("permissive",),
    )
    CASES += decon_cases(
        f"{_m}.error.numit",
        _method,
        _inputs(_MS, 51),
        2.5,
        {"Background": 100},
        expect_error=True,
        backends=("permissive",),
    )
    CASES += decon_cases(
        f"{_m}.numit_negative",
        _method,
        _inputs(_MS, 51),
        -1,
        {"Background": 100},
        nargout=3 if _method == "simplified" else 2,
    )

# debug output (DeconRFFTTest/debugOutputsEqualTheMatlabLoop): [iterations, saveStep]
DEBUG_CASES = {
    "everyIteration": (6, 1),
    "everySecond": (7, 2),
    "stepPastTheEnd": (2, 5),
    "nonIntegerStep": (9, 2.5),
}
# MATLAB's bookkeeping errors: istp undefined
DEBUG_ERRORS = {("simplified", "stepPastTheEnd"), ("simplified", "nonIntegerStep"),
                ("omw", "nonIntegerStep")}  # fmt: skip
for _m, (_method, _) in METHODS.items():
    for _dn, (_n, _step) in DEBUG_CASES.items():
        for _path, _use in (("kernels", True), ("loop", False)):
            CASES += decon_cases(
                f"{_m}.debug.{_dn}.{_path}",
                _method,
                _inputs(_MS, 52),
                _n,
                {
                    "Background": 100,
                    "save16bit": False,
                    "debug": True,
                    "debug_folder": "dbg/",
                    "saveStep": _step,
                    "useMex": _use,
                },
                nargout=3 if _method == "simplified" else 2,
                files=True,
                expect_error=(_method, _dn) in DEBUG_ERRORS,
            )
    # 89,000 or more voxels: M05's fixed-order mean (the reference changes at switch #4)
    _big = [64, 50, 40]
    _cs = decon_cases(
        f"{_m}.debug.big",
        _method,
        _inputs(_big, 53),
        4,
        {
            "Background": 100,
            "save16bit": False,
            "debug": True,
            "debug_folder": "dbg/",
            "saveStep": 1,
        },
        nargout=3 if _method == "simplified" else 2,
        files=True,
        switch4=True,
    )
    SWITCH4.update(c.id for c in _cs)
    CASES += _cs

# DeconRFFTTest/resultDoesNotDependOnTheThreadCount's data: a 63^3 PSF (89,000 or more elements:
# the normalization's fixed-order sum) and 1 Mi voxels (lambda's blocks)
_BIGS = [128, 128, 64]
for _m, (_method, _n) in METHODS.items():
    CASES += decon_cases(
        f"{_m}.psf63",
        _method,
        {
            "im": V(lambda: blurred_beads(_BIGS, 54)),
            "psf": V(lambda: to_matlab(decon_psf((63, 63, 63), (2.0, 1.5)))),
            "psfb": V(lambda: to_matlab(decon_psf((63, 63, 63), (1.2, 1.0)))),
        },
        _n,
        {"Background": 100, "save16bit": False},
        marks=("slow",),
    )


def _check(case, tmp_path):
    if case.id.endswith(".exact"):
        reason = digests_hold_here()
        if reason:
            pytest.skip(reason)
    oracle.check(case, tmp_path)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    _check(case, tmp_path)


_THREADED = [
    c
    for c in CASES
    if not c.expect_error
    and any(c.id.startswith(p) for p in ("rl.size.", "omw.size.", "rl.opt.all", "omw.opt.all"))
    or c.id.startswith(("rl.iter30.uint16", "rl.debug.everySecond.kernels", "rl.loop.uint16"))
]


@oracle.parametrize(_THREADED)
@pytest.mark.parametrize("threads", [1, 2, 7, 24])
def test_thread_counts(case, threads, tmp_path):
    """the same results with every thread count of ours"""
    kw = dict(case.py_kwargs or {})
    kw["threads"] = threads
    variant = dataclasses.replace(case, py_kwargs=kw)
    variant._module_globals = case._module_globals
    _check(variant, tmp_path)


# ---------------------------------------------------------------------------------------------
# without MATLAB (these run in CI)
# ---------------------------------------------------------------------------------------------


def py_data(msize, seed=1, cls="uint16"):
    return from_matlab(blurred_beads(msize, seed, cls))


def py_psf(msize, which="f"):
    return from_matlab(make_psf(msize, which))


LAYOUTS = [("permissive", "matlab"), ("permissive", "native")] + (
    [("exact", "matlab")] if EXACT else []
)
SHAPES = [(36, 48, 46), (33, 47, 45), (50, 64), (33, 2, 64)]


@pytest.mark.parametrize("be,layout", LAYOUTS)
@pytest.mark.parametrize("shape", SHAPES)
@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_kernels_equal_the_matlab_loop(be, layout, shape, method):
    """DeconRFFTTest/mexEqualsMatlabLoop in the port: in MATLAB's layout the kernels' iterations
    equal MATLAB's loop on full spectra (matlabcompat.fft) on the same backend, values bit for bit
    (the kernels' update makes -0 +0); in the native layout within 0006's tolerance"""
    ms = list(size_to_matlab(shape))
    im, psf, psfb = py_data(ms), py_psf(ms), py_psf(ms, "b")
    kw = {"background": 100.0, "save_16bit": False, "backend": be, "threads": 3}
    if method == "simplified":
        a = rl.lucy(im, psf, 10, layout=layout, **kw).decon
        b = rl.lucy(im, psf, 10, use_mex=False, **kw).decon
    else:
        a = rl.lucy_omw(im, psf, psfb, 2, layout=layout, **kw).decon
        b = rl.lucy_omw(im, psf, psfb, 2, use_mex=False, **kw).decon
    assert a.dtype == b.dtype == np.float32 and a.shape == b.shape == shape
    if layout == "matlab":
        assert np.array_equal(a, b)
    else:
        rel = np.linalg.norm(a.astype(np.float64) - b) / np.linalg.norm(b.astype(np.float64))
        assert rel < 1e-5


@pytest.mark.parametrize("be,layout", LAYOUTS)
@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_results_do_not_depend_on_the_thread_count(be, layout, method):
    """DeconRFFTTest/resultDoesNotDependOnTheThreadCount, DeconAccelTest/
    deconDoesNotDependOnTheThreadCount: 89,000 or more voxels and a 63^3 PSF (fixed-order sums of
    the PSF normalization and of lambda), also with debug output's error metric"""
    ms = [66, 64, 33]
    im = py_data(ms, 3)
    psf = decon_psf((33, 63, 63), (2.0, 1.5))
    psfb = decon_psf((33, 63, 63), (1.2, 1.0))
    outs = []
    for t in (1, 2, 7, 24):
        kw = {"background": 100.0, "save_16bit": False, "backend": be, "threads": t}
        if method == "simplified":
            r = rl.lucy(im, psf, 6, layout=layout, **kw)
        else:
            r = rl.lucy_omw(im, psf, psfb, 2, layout=layout, **kw)
        outs.append(r.decon)
    for o in outs[1:]:
        assert o.tobytes() == outs[0].tobytes()


def test_debug_error_does_not_depend_on_the_thread_count(tmp_path):
    """DeconAccelTest/debugErrorDoesNotDependOnTheThreadCount: err_mat of 89,000 or more voxels"""
    im, psf, psfb = py_data([64, 50, 40], 4), py_psf([64, 50, 40]), py_psf([64, 50, 40], "b")
    ref = None
    for t in (1, 2, 5, 24):
        kw = {"background": 100.0, "debug": True, "threads": t}
        a = rl.lucy(im, psf, 4, debug_folder=str(tmp_path / f"rl{t}"), **kw).err_mat
        b = rl.lucy_omw(im, psf, psfb, 4, debug_folder=str(tmp_path / f"omw{t}"), **kw).err_mat
        c = rl.lucy_omw(
            im, psf, psfb, 4, debug_folder=str(tmp_path / f"loop{t}"), use_mex=False, **kw
        ).err_mat
        if ref is None:
            ref = (a, b, c)
        for x, y in zip((a, b, c), ref, strict=True):
            assert x.tobytes() == y.tobytes()


@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_zero_iterations_give_the_background_subtracted_input(method):
    """DeconRFFTTest/zeroIterationsGiveTheBackgroundSubtractedInput: the subtraction in double,
    rounded to single; MATLAB's max keeps -0 (RL), the kernel's OMW makes it +0"""
    im = py_data([48, 46, 36], 5, "single")
    for x in (im, im.astype(np.uint16)):
        if method == "simplified":
            r = rl.lucy(x, py_psf([48, 46, 36]), 0, background=100.3, save_16bit=False)
            assert r.iterations == 0 and r.err_mat.shape == (0, 4)
        else:
            r = rl.lucy_omw(x, py_psf([48, 46, 36]), py_psf([48, 46, 36], "b"), 0,
                            background=100.3, save_16bit=False)  # fmt: skip
            assert r.err_mat.shape == (0, 4)
        e = np.maximum((x.astype(np.float64) - 100.3).astype(np.float32), np.float32(0))
        assert np.array_equal(r.decon, e)


def test_use_rfft_kernels_takes_what_the_mex_takes():
    """DeconRFFTTest/unsupportedCallsRunTheMatlabLoop: double data, a background per voxel, a
    single row along y or x run the loop; uint16 and debug output the kernels"""
    im = py_data([48, 46, 36])
    can = [
        rl.use_rfft_kernels(None, "omw", im.astype(np.float64), False, 100)[1],
        rl.use_rfft_kernels(None, "omw", im, True, 100)[1],
        rl.use_rfft_kernels(None, "omw", im, False, np.full(im.shape, 100.0))[1],
        rl.use_rfft_kernels(None, "omw", im[:, :1, :], False, 100)[1],
        rl.use_rfft_kernels(None, "omw", im[:, :, :1], False, 100)[1],
        rl.use_rfft_kernels(None, "simplified", im, False, 100)[1],
    ]
    assert can == [False, True, False, False, False, True]
    # two rows or columns are enough; a single slice is an image
    assert rl.use_rfft_kernels(None, "omw", im[:, :2, :2], False, 100) == (True, True)
    assert rl.use_rfft_kernels(None, "omw", im[:1], False, 100) == (True, True)
    assert rl.use_rfft_kernels(True, "omw", im, False, 100) == (True, True)
    assert rl.use_rfft_kernels(False, "omw", im, False, 100) == (False, True)
    with pytest.raises(rl.DeconError, match="does not support the data"):
        rl.use_rfft_kernels(True, "simplified", im.astype(np.float64), False, 100)


class _Counter:
    """counts the OTFs the decon functions make (both paths)"""

    def __init__(self, monkeypatch):
        from petakit5d.kernels import decon as K

        self.n = 0
        orig_otf, orig_full = K.LucyRFFT.otf, rl._psf2otf

        def otf(obj, psf):
            self.n += 1
            return orig_otf(obj, psf)

        def full(*a, **k):
            self.n += 1
            return orig_full(*a, **k)

        monkeypatch.setattr(K.LucyRFFT, "otf", otf)
        monkeypatch.setattr(rl, "_psf2otf", full)


@pytest.mark.parametrize("method", ["simplified", "omw"])
@pytest.mark.parametrize("use_mex", [None, False])
def test_otf_is_cached_for_the_same_psf_bits(monkeypatch, method, use_mex):
    """DeconFunctionTest/otfIsCachedForTheSamePSFBitsAndThreads: the same values (another array)
    use the cached OTF, a PSF that differs in the sign of a zero gets a new one; another thread
    count uses the cached one (the normalization doesn't depend on it in the port)"""
    rl.clear_otf_cache()
    count = _Counter(monkeypatch)
    nf = 1 if method == "simplified" else 2
    im = py_data([48, 46, 36], 6)
    psf = py_psf([48, 46, 36]).copy()
    psf[0, 0, 0] = 0
    psfb = py_psf([48, 46, 36], "b")

    def run(p, t=3):
        kw = {"background": 100.0, "save_16bit": False, "use_mex": use_mex, "threads": t}
        if method == "simplified":
            return rl.lucy(im, p, 3, **kw).decon
        return rl.lucy_omw(im, p, psfb, 2, **kw).decon

    a = run(psf)
    b = run(psf + np.float32(0))
    assert count.n == nf and np.array_equal(a, b)
    neg = psf.copy()
    neg[0, 0, 0] = -0.0
    run(neg)
    assert count.n == 2 * nf
    c = run(neg, t=1)
    assert count.n == 2 * nf and c.tobytes() == run(neg, t=5).tobytes()
    rl.clear_otf_cache()


@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_changed_psf_with_same_energy_is_not_ignored(method):
    """DeconFunctionTest/changedPSFWithSameEnergyIsNotIgnored and
    changedBackProjectorWithSameEnergyIsNotIgnored"""
    rl.clear_otf_cache()
    im = py_data([48, 46, 36], 7)
    psf1 = np.round(py_psf([48, 46, 36]) * np.float32(2**19)).astype(np.float32)
    psf2 = np.roll(psf1, 2, axis=0)
    kw = {"background": 100.0, "save_16bit": False}
    if method == "simplified":
        rl.lucy(im, psf1, 3, **kw)
        out = rl.lucy(im, psf2, 3, **kw).decon
        rl.clear_otf_cache()
        fresh = rl.lucy(im, psf2, 3, **kw).decon
    else:
        rl.lucy_omw(im, psf1, psf1, 2, **kw)
        out = rl.lucy_omw(im, psf1, psf2, 2, **kw).decon
        rl.clear_otf_cache()
        fresh = rl.lucy_omw(im, psf1, psf2, 2, **kw).decon
    assert np.array_equal(out, fresh)
    rl.clear_otf_cache()


def test_psf_normalizations_use_the_fixed_order_sum(monkeypatch):
    """DeconAccelTest/psfNormalizationsUseTheFixedOrderSum: the PSF and the back projector are
    normalized by decon_sum_double (kernels.decon.sum_double)"""
    from petakit5d.kernels import decon as K

    sizes = []
    orig = K.sum_double

    def spy(x, **kw):
        sizes.append(np.size(x))
        return orig(x, **kw)

    monkeypatch.setattr(K, "sum_double", spy)
    rl.clear_otf_cache()
    im = py_data([48, 46, 36], 8)
    psf, psfb = py_psf([48, 46, 36]), py_psf([48, 46, 36], "b")
    rl.lucy(im, psf, 1, background=100.0)
    rl.lucy_omw(im, psf, psfb, 1, background=100.0)
    assert sizes[:3] == [psf.size, psf.size, psfb.size]
    rl.clear_otf_cache()


# GoldenOutputTest's decon scenarios: reference-only cases with the MATLAB calls (and so the keys
# and references) of tests/matlab_map/test_golden.py's cases (F06): the scenario's inputs, the
# oracle tree's result and the stored reference

GOLDEN = {
    "decon_rl_simplified": {"outputs": ["decon"], "inputs": ["im", "psf"]},
    "decon_omw": {"outputs": ["backProjector", "decon"], "inputs": ["im", "psf"]},
}
GOLDEN_CASES = {
    f"golden.{s}.{what}": Case(
        f"golden.{s}.{what}",
        matlab="pkpy_golden",
        args=[what, s, mcell(*spec["inputs" if what == "inputs" else "outputs"]), Work()],
        nargout=len(spec["inputs" if what == "inputs" else "outputs"]),
        outputs=list(spec["inputs" if what == "inputs" else "outputs"]),
        helpers=["tests/matlab_map/matlab"],
        files=False,
    )
    for s, spec in GOLDEN.items()
    for what in ("stored", "compute", "inputs")
}


@oracle.parametrize(list(GOLDEN_CASES.values()))
def test_golden_references(case, tmp_path):
    """the references of the golden scenarios exist (store or manifest)"""
    oracle.check(case, tmp_path)


def _rel_l2(a: np.ndarray, e: np.ndarray) -> float:
    a = np.asarray(a, np.float64).reshape(-1, order="F")
    e = np.asarray(e, np.float64).reshape(-1, order="F")
    return math.sqrt(math.fsum((a - e) ** 2)) / max(math.sqrt(math.fsum(e**2)), 2.0**-52)


@pytest.mark.oracle
@pytest.mark.parametrize("be", ["exact", "permissive"])
@pytest.mark.parametrize("scenario", list(GOLDEN))
def test_golden_scenario(scenario, be):
    """GoldenOutputTest/outputMatchesReference(scenario=decon_rl_simplified, decon_omw): the
    scenario's decon of its inputs equals the oracle tree's result (exact backend: bit for bit;
    permissive: 0006's relative L2 <= 1e-4) and the stored reference within GoldenOutputTest's
    1e-4. OMW takes the oracle tree's back projector (omw_backprojector_generation is R02's)."""
    if be == "exact" and digests_hold_here():
        pytest.skip(digests_hold_here())
    computed = oracle.reference(GOLDEN_CASES[f"golden.{scenario}.compute"])
    inputs = oracle.reference(GOLDEN_CASES[f"golden.{scenario}.inputs"])
    if computed.source != "store" or inputs.source != "store":
        pytest.skip("needs the golden store (the manifests hold hashes only)")
    im, psf = (from_matlab(v.data) for v in inputs.outputs)
    if scenario == "decon_rl_simplified":
        out = rl.lucy(im, psf, 10, background=100, save_16bit=False, backend=be)
        ref = computed.outputs[0].data
    else:
        bp = from_matlab(computed.outputs[0].data)
        out = rl.lucy_omw(im, psf, bp, 2, background=100, save_16bit=False, backend=be)
        ref = computed.outputs[1].data
    m = to_matlab(out.decon)
    assert m.dtype == ref.dtype and m.shape == ref.shape
    if be == "exact":
        assert m.tobytes(order="F") == np.asarray(ref).tobytes(order="F")
    else:
        assert _rel_l2(m, ref) <= 1e-4
    stored = oracle.reference(GOLDEN_CASES[f"golden.{scenario}.stored"])
    if stored.source == "store":
        assert _rel_l2(m, stored.outputs[-1].data) <= 1e-4


# DeconFunctionTest's properties, on beads blurred with a known PSF

IM_SIZE = (40, 48, 48)
PSF_SIGMA = (2.0, 1.5)


def dft_beads(positions=None, amplitude=1000.0, background=100.0):
    if positions is None:
        positions = [[9, 11, 13], [19.2, 23.3, 22.6], [27, 35, 29], [24, 14, 34], [14, 32, 11]]
    sig = (math.sqrt(0.8**2 + 2.0**2), math.sqrt(0.8**2 + 1.5**2))
    im = beads(IM_SIZE, positions, sig, amplitude, background)
    return np.floor(im + np.float32(0.5)).astype(np.uint16)


def dft_psf():
    return decon_psf((31, 31, 31), PSF_SIGMA)


def run_method(method, im, n, **kw):
    kw.setdefault("save_16bit", False)
    if method == "simplified":
        return rl.lucy(im, dft_psf(), n, **kw).decon
    return rl.lucy_omw(im, dft_psf(), decon_psf((31, 31, 31), (1.2, 1.0)), n, **kw).decon


@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_delta_psf_returns_background_subtracted_input(method):
    rng = CounterRNG(9, "delta")
    im = (rng.uniform((16, 24, 24)) * 950 + 50).astype(np.uint16)
    delta = np.zeros((3, 3, 3), np.float32)
    delta[1, 1, 1] = 1
    if method == "simplified":
        out = rl.lucy(im, delta, 5, background=100.0).decon
    else:
        out = rl.lucy_omw(im, delta, delta, 5, background=100.0).decon
    ref = np.floor(np.maximum(im.astype(np.float32) - 100, 0) + np.float32(0.5)).astype(np.uint16)
    assert np.array_equal(out, ref)


def test_flux_is_conserved():
    im = dft_beads()
    out = run_method("simplified", im, 10, background=100.0)
    s = math.fsum(out.astype(np.float64).ravel())
    e = math.fsum(np.maximum(im.astype(np.float64) - 100, 0).ravel())
    assert abs(s - e) <= 1e-4 * e


@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_regions_at_background_stay_zero(method):
    im = dft_beads()
    out = run_method(method, im, 5, background=100.0)
    at = im <= 100
    assert at.sum() > 0.5 * im.size
    assert np.count_nonzero(out[at]) == 0 and out.min() >= 0


@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_output_scales_linearly_with_input(method):
    im = dft_beads(background=10.0)
    a = run_method(method, im, 5, background=0.0)
    b = run_method(method, (im.astype(np.uint32) * 3).astype(np.uint16), 5, background=0.0)
    d = b.astype(np.float64) - 3 * a.astype(np.float64)
    assert np.linalg.norm(d) <= 1e-5 * np.linalg.norm(b.astype(np.float64))


@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_point_source_is_sharpened(method):
    pos = [19.2, 23.3, 22.6]
    im = dft_beads([pos], 1000.0)
    out = run_method(method, im, 10 if method == "simplified" else 2, background=100.0)
    inp = np.maximum(im.astype(np.float32) - 100, 0)
    assert out.max() > 1.5 * inp.max()

    def moments(v):
        lo = [max(round(p) - 8, 0) for p in pos]
        win = v[lo[0] : round(pos[0]) + 9, lo[1] : round(pos[1]) + 9, lo[2] : round(pos[2]) + 9]
        win = win.astype(np.float64)
        g = np.meshgrid(
            *[np.arange(n) + o for n, o in zip(win.shape, lo, strict=True)], indexing="ij"
        )
        s = win.sum()
        c = [float((gi * win).sum() / s) for gi in g]
        var = [float(((gi - ci) ** 2 * win).sum() / s) for gi, ci in zip(g, c, strict=True)]
        return np.array(c), np.array(var)

    ci, vi = moments(inp)
    co, vo = moments(out)
    assert np.all(np.abs(co - ci) <= 0.05)
    assert np.all(vo < 0.8 * vi)


@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_edge_erosion_zeroes_the_border_shell(method):
    im = from_matlab(textured([48, 48, 40]))
    e = 3
    out = run_method(method, im, 3, background=0.0, edge_erosion=e)
    shell = np.ones(im.shape, bool)
    shell[e:-e, e:-e, e:-e] = False
    assert np.count_nonzero(out[shell]) == 0
    inner = out[e:-e, e:-e, e:-e]
    assert np.count_nonzero(inner) / inner.size > 0.99


def test_damp_factor_bounds_the_result():
    d = 1.2
    im = dft_beads()
    out = run_method("simplified", im, 10, background=100.0, damp_factor=d)
    bound = (d * np.maximum(im.astype(np.float32) - 100, 0).astype(np.float64)).astype(np.float32)
    assert (out - bound).max() <= 1e-3
    assert np.count_nonzero((np.abs(out - bound) < 1e-3) & (bound > 0)) > 0


def test_damp_is_matlabs_expression():
    """``J .* (J <= d * I) + min(J, d * I .* (J >= d * I))`` element by element: d * I in double
    rounded to single; at a tie both terms are J (MATLAB's expression doubles it); MATLAB's min
    gives the other value for NaN"""
    d = 1.2
    i = np.array([2.0, 2.0, 2.0, 0.0, 1.0, 3.0, 1.0], np.float32)
    di = (d * i.astype(np.float64)).astype(np.float32)
    j = np.array([di[0], 1.0, 5.0, 0.0, np.nan, -1.0, -0.0], np.float32)
    out = rl._damp(j, i, d, 1)
    # a tie: J + min(J, J); below: J + min(J, 0); above: 0 + min(J, dI); J < 0: J + min(J, 0)
    ref = np.array([2 * di[0], 1.0, di[2], 0.0, 0.0, -2.0], np.float32)
    assert np.array_equal(out[[0, 1, 2, 3, 5]], ref[[0, 1, 2, 3, 5]])
    # NaN J: NaN .* false is NaN, min(NaN, 0) is 0: NaN + 0
    assert np.isnan(out[4])
    # J = -0: -0 .* true + min(-0, +0) (a tie: MATLAB's min gives -0) is -0
    assert out[6] == 0 and np.signbit(out[6])


def test_the_mask_is_the_nonzero_input():
    """I_mask = I ~= 0 (of single(I), before the background): negative and NaN data are in the
    mask, so deconOffset reaches them; zeros are not"""
    im = np.full((4, 6, 5), 7.0, np.float32)
    im[0, 0, :3] = (-5.0, np.nan, 0.0)
    psf = np.zeros((3, 3, 3), np.float32)
    psf[1, 1, 1] = 1
    for method in ("simplified", "omw"):
        kw = {"save_16bit": False, "decon_offset": 3.0, "use_mex": False}
        if method == "simplified":
            out = rl.lucy(im, psf, 0, **kw).decon
        else:
            out = rl.lucy_omw(im, psf, psf, 0, **kw).decon
        assert out[0, 0, :3].tolist() == [3.0, 3.0, 0.0] and out[1, 1, 1] == 10.0


def test_uint16_output_is_matlabs_uint16():
    """save16bit converts as MATLAB's uint16(single): floor(|x| + 0.5) in single, so the single
    just below 0.5 gives 1 (cast would give 0)"""
    x = np.nextafter(np.float32(0.5), np.float32(0))
    im = np.full((4, 6, 5), x, np.float32)
    psf = np.zeros((3, 3, 3), np.float32)
    psf[1, 1, 1] = 1
    for method in ("simplified", "omw"):
        if method == "simplified":
            out = rl.lucy(im, psf, 0, save_16bit=True).decon
        else:
            out = rl.lucy_omw(im, psf, psf, 0, save_16bit=True).decon
        assert out.dtype == np.uint16 and np.all(out == 1)


@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_bbox_crops_the_result(method):
    im = dft_beads()
    full = run_method(method, im, 3, background=100.0, save_16bit=True)
    box = (6, 4, 5, 30, 40, 41)
    out = run_method(method, im, 3, background=100.0, save_16bit=True, bbox=box)
    assert np.array_equal(out, full[6:30, 4:40, 5:41])


def test_overlap_function_runs_once_and_returns_its_output():
    """DeconOverlapTest/overlapFunctionRunsOnceAndReturnsItsOutput (its CPU part)"""
    im = dft_beads()
    calls = []

    def overlap():
        calls.append(1)
        return 42

    a = rl.lucy_omw(im, dft_psf(), dft_psf(), 2, background=100.0)
    b = rl.lucy_omw(im, dft_psf(), dft_psf(), 2, background=100.0, overlap_fcn=overlap)
    assert np.array_equal(a.decon, b.decon) and b.overlap_out == 42 and len(calls) == 1
    a = rl.lucy(im, dft_psf(), 3, background=100.0)
    b = rl.lucy(im, dft_psf(), 3, background=100.0, overlap_fcn=overlap)
    assert np.array_equal(a.decon, b.decon) and b.overlap_out == 42 and len(calls) == 2
    assert rl.lucy(im, dft_psf(), 1).overlap_out is None


@pytest.mark.parametrize("be,layout", LAYOUTS)
def test_a_volume_of_one_slice_is_a_2d_image(be, layout):
    """MATLAB's size of a (1, ny, nx) volume has no third dimension: the 2D transforms (MATLAB's
    layout halves an image along x, a volume along z), the 2D result in the volume's shape"""
    im = py_data([64, 50], 10)
    psf = py_psf([64, 50])
    kw = {"background": 100.0, "backend": be, "layout": layout}
    for fn, args in ((rl.lucy, (psf,)), (rl.lucy_omw, (psf, py_psf([64, 50], "b")))):
        a = fn(im, *args, 5, save_16bit=False, **kw).decon
        b = fn(im[None], *(p[None] for p in args), 5, save_16bit=False, **kw).decon
        assert b.shape == (1, 64, 50) and a.tobytes() == b[0].tobytes()
    c = rl.lucy(im[None], psf[None], 5, bbox=(0, 2, 3, 1, 40, 45), **kw).decon
    assert c.shape == (1, 38, 42)


# the GPU choice (DeconFunctionTest's GPU methods, with gpu.device's selection stubbed)


def _stub_gpus(monkeypatch, found: int, reason: str):
    from petakit5d.gpu import device as G

    def select(index=None, need_bytes=0):
        raise G.GPUUnavailableError(reason)

    monkeypatch.setattr(G, "select_device", select)
    monkeypatch.setattr(G, "gpu_count", lambda refresh=False: G.GPUCount(0, found, reason))


@pytest.mark.parametrize("method", ["simplified", "omw"])
@pytest.mark.parametrize("found", [1, 0])
def test_gpu_request_without_usable_gpu_errors(monkeypatch, method, found):
    """DeconFunctionTest/gpuRequestWithoutUsableGPUErrors and noUsableGPUMessagesGiveMATLABsReason:
    use_gpu True errors (decon:noUsableGPU) with the count and the reason, also after a call that
    cached a CPU OTF"""
    from petakit5d.gpu import device as G

    reason = "Selecting GPU 1 failed: stub reason." if found else "No GPU is found."
    _stub_gpus(monkeypatch, found, reason)
    rl.clear_otf_cache()
    im = dft_beads()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", G.GPUFallbackWarning)
        if method == "simplified":
            rl.decon_lucy(im, dft_psf(), 1, background=100.0)
        else:
            rl.decon_lucy_omw(im, dft_psf(), dft_psf(), 1, background=100.0)
    with pytest.raises(G.NoUsableGPUError) as ei:
        if method == "simplified":
            rl.decon_lucy(im, dft_psf(), 3, background=100.0, use_gpu=True)
        else:
            rl.decon_lucy_omw(im, dft_psf(), dft_psf(), 3, background=100.0, use_gpu=True)
    assert ei.value.identifier == "decon:noUsableGPU"
    assert f"No GPU is available ({found} found, CUDA_VISIBLE_DEVICES=" in str(ei.value)
    assert ("Selecting GPU 1 failed: stub reason." in str(ei.value)) == bool(found)
    rl.clear_otf_cache()


@pytest.mark.parametrize("method", ["simplified", "omw"])
def test_default_falls_back_to_the_cpu_with_a_warning(monkeypatch, method):
    """DeconFunctionTest/defaultFallsBackToTheCPUWithAWarning: without use_gpu (None) a GPU if
    one is usable, else the CPU with a warning (once per new OTF), the CPU's result"""
    from petakit5d.gpu import device as G

    rl.clear_otf_cache()
    im = dft_beads()

    def run(**kw):
        if method == "simplified":
            return rl.decon_lucy(im, dft_psf(), 3, background=100.0, save_16bit=False, **kw)
        return rl.decon_lucy_omw(im, dft_psf(), dft_psf(), 3, background=100.0, save_16bit=False,
                                 **kw)  # fmt: skip

    cpu = run(use_gpu=False).decon
    _stub_gpus(monkeypatch, 1, "Selecting GPU 1 failed: stub reason.")
    rl.clear_otf_cache()
    with pytest.warns(G.GPUFallbackWarning) as wl:
        out = run().decon
    msg = str(wl[0].message)
    assert np.array_equal(out, cpu)
    assert "No GPU is available (1 found" in msg and "so the deconvolution runs on the CPU." in msg
    assert "Selecting GPU 1 failed: stub reason." in msg
    with warnings.catch_warnings():
        warnings.simplefilter("error", G.GPUFallbackWarning)
        assert np.array_equal(run().decon, cpu)  # the cached OTF: no new check, no warning
    rl.clear_otf_cache()


def test_matlab_names_take_and_give_matlabs_values():
    """the MATLAB-name layer: MATLAB's (y, x, z) arrays and box, the outputs in MATLAB's form"""
    from petakit5d.params import find

    im = blurred_beads([48, 46, 36], 11)
    psf = make_psf([48, 46, 36])
    spec = find("decon_lucy_function")
    j, err, k, ov = spec.call_matlab(
        im, psf, 3, "Background", 100, "useGPU", False, "bbox", [2, 3, 4, 40, 41, 30]
    )
    ref = rl.lucy(
        from_matlab(im), from_matlab(psf), 3, background=100.0, bbox=(3, 1, 2, 30, 40, 41)
    )
    assert j.shape == (39, 39, 27) and np.array_equal(j, to_matlab(ref.decon))
    assert err.shape == (3, 4) and k == 3.0 and ov is None
    j, err, ov = find("decon_lucy_omw_function").call_matlab(
        im, psf, psf, 0, "Background", 100, "useGPU", False
    )
    assert j.shape == (48, 46, 36) and err.shape == (0, 4)
