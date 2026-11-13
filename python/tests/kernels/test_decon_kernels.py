"""The deconvolution kernels (petakit5d.kernels.decon) against the mex files and MATLAB's code.

Oracle cases (MATLAB R2024b runs the oracle tree's mex files on the same inputs). Inputs are made
in MATLAB's form; the roles give the Python calls the port's form and map the results back.

- rfft.*: decon_lucy_rfft_mex ('otf', 'omw', 'rlstep', 'fft1'): DeconRFFTTest's sizes (even and odd,
  2D, length 64 along y, x or z, 64 x 230 and 64 x 210, two rows along y) and sizes with 128, uint16
  and single data, backgrounds with fractions, 0 (rfft.omw0.*) to 50 OMW iterations, OMW from an
  estimate, ten RL steps in a row, MATLAB's two FFT modes (grain 32 and 0), MATLAB's own fft for
  'fft1', the mex's argument errors. Each case runs twice on the Python side: ``.exact`` on the
  exact FFT backend in MATLAB's layout (bit for bit; skipped without the add-on, or where the C
  library gives FFTW other twiddle factors than the reference platform), ``.permissive`` on the
  default backend
  with the native half spectra (decision record 0006's tolerance: relative L2 <= 1e-4; the counts
  are in the K05 report).
- accel.*: decon_lucy_accel_mex ('lambda' with its sums, 'predict', 'sum') and the MATLAB functions
  decon_lucy_accel and decon_sum_double around their switches (44, 89000 elements), single and
  double, -0, NaN, Inf, subnormal numbers and NaNs of several payloads, (z, y, x) arrays and
  MATLAB's layout (``.mlayout``); DeconAccelTest's kernel-level cases (the class is R03's).
- erosion.*: decon_mask_edge_erosion (its mex path) for DeconMaskErosionTest's shapes and sizes
  and EdgeErosion 1 to 41, the demo setting, the mex with random neighborhoods and paddings,
  numeric neighborhoods, a line along z, the argument errors.
- omw.*: omw_backprojector_mex's 'convhull', 'window', 'components', 'sortdesc', 'medfilt3'.

test_thread_counts runs the cases again with 1, 2, 7 and 24 threads (as many as there are CPUs)
against the same references.
"""

from __future__ import annotations

import dataclasses
import functools
import math

import numpy as np
import pytest

from petakit5d import fft as P
from petakit5d.core.arrays import from_matlab, matlab_order, size_to_matlab, to_matlab
from petakit5d.core.dtypes import NAN
from petakit5d.kernels import decon as D
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Index1, PerAxis, Value, Vol
from petakit5d.testing.compare import tol
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG, beads, decon_psf

BACKENDS = P.available()
EXACT = "exact" in BACKENDS
INF = float("inf")
PERMISSIVE_TOL = tol(
    reason="0006: deconvolution with the permissive FFT backend (relative L2 <= 1e-4)", rel_l2=1e-4
)
PERMISSIVE_OTF = tol(
    reason="0006: single transforms with the permissive FFT backend", rel_max=1e-5, rel_l2=1e-6
)


def test_the_kernels_use_the_processs_fft_state():
    # petakit5d._ext.decon attaches to petakit5d._ext.fft's backends, planner lock and plan cache at
    # import (tests/fft/test_shared_state.py plans from both modules at once)
    assert D._ext()._fft_state_id() == P._ext()._fft_state_id()


def exact_skip_reason() -> str:
    """'' where the exact backend gives MATLAB's bits here, else why not (the add-on, x86-64 with
    AVX, the reference platform's twiddle factors: petakit5d.testing.fft_reference)"""
    return digests_hold_here()


def V(make):
    return Value(make)


# ---------------------------------------------------------------------------------------------
# inputs, in MATLAB's form
# ---------------------------------------------------------------------------------------------

BLUR_SIGMA = (math.sqrt(0.8**2 + 2.0**2), math.sqrt(0.8**2 + 1.5**2))  # (z, yx)


def blurred_beads(msize, seed: int, cls: str = "uint16") -> np.ndarray:
    """DeconRFFTTest.blurredBeads: beads blurred by the test PSF plus a little noise, on a
    background of 100; MATLAB's form"""
    sz = tuple(int(v) for v in msize)
    py = (sz[2], sz[0], sz[1]) if len(sz) == 3 else sz
    rng = CounterRNG(seed, "beads")
    n = max(3, round(math.prod(sz) / 4000))
    pos = rng.uniform((n, len(py))) * (np.array(py, dtype=np.float64) - 1)
    im = beads(py, pos, BLUR_SIGMA if len(py) == 3 else BLUR_SIGMA[1], 1000.0, 100.0)
    im = im + np.float32(3) * rng.uniform(py, np.float32)
    if cls == "uint16":
        out = np.floor(im + np.float32(0.5)).astype(np.uint16)  # im >= 0
    else:
        out = (im + np.float32(0.37)).astype(np.float32)
    return to_matlab(out)


def make_psf(msize, which: str = "f") -> np.ndarray:
    """the PSF (or a second one for the back projector) cropped to the data as RLdecon crops,
    scaled by a power of two to about unit sum; MATLAB's form"""
    sigma = (2.0, 1.5) if which == "f" else (1.2, 1.0)
    psf = decon_psf((31, 31, 31), sigma)
    psf = np.ascontiguousarray(psf * np.float32(2.0**-19))
    if which == "neg":
        # a back projector with negative side lobes, as OMW back projectors have
        wide = decon_psf((31, 31, 31), (3.0, 2.5)) * np.float32(2.0**-21)
        psf = np.ascontiguousarray(psf * np.float32(2.0) - wide).astype(np.float32)
    if len(msize) == 2:
        psf = psf[15]
        py = tuple(msize)
    else:
        py = (msize[2], msize[0], msize[1])
    s = [max(0, (p - n) // 2) for p, n in zip(psf.shape, py, strict=True)]
    t = [a + min(p, n) for a, p, n in zip(s, psf.shape, py, strict=True)]
    psf = psf[tuple(slice(a, b) for a, b in zip(s, t, strict=True))]
    return to_matlab(np.ascontiguousarray(psf))


# ---------------------------------------------------------------------------------------------
# Python sides of the rfft cases (MATLAB's arguments, mapped by the roles)
# ---------------------------------------------------------------------------------------------


def native_half_to_matlab(h: np.ndarray, shape) -> np.ndarray:
    """MATLAB's half OTF (halved along z; 2D along x) of a native half spectrum (halved along x)"""
    if len(shape) == 2:
        return h
    nz, ny, nx = shape
    nxh = h.shape[-1]
    full = np.empty((nz, ny, nx), dtype=h.dtype)
    full[..., :nxh] = h
    kz = (-np.arange(nz)) % nz
    ky = (-np.arange(ny)) % ny
    kx = nx - np.arange(nxh, nx)
    full[..., nxh:] = np.conj(h[kz][:, ky][:, :, kx])
    return to_matlab(full)[:, :, : nz // 2 + 1]


def kernel(shape, grain, backend, threads) -> D.LucyRFFT:
    return D.LucyRFFT(shape, backend=backend, threads=threads, matlab_parallel=float(grain) > 0)


def rfft_otf(psf, out_size, nt, grain, *, backend, threads=None):
    shape = tuple(int(v) for v in out_size)
    k = kernel(shape, grain, backend, threads)
    h = k.otf(psf)
    return h if k.layout == "matlab" else native_half_to_matlab(h, shape)


def rfft_omw(psf, raw, psfb, n_iter, grain, bg, *, backend, threads=None):
    k = kernel(raw.shape, grain, backend, threads)
    return k.omw(raw, k.otf(psf), k.otf(psfb), n_iter, bg)


def rfft_omw_continued(psf, raw, psfb, n1, n2, grain, bg, *, backend, threads=None):
    k = kernel(raw.shape, grain, backend, threads)
    hf, hb = k.otf(psf), k.otf(psfb)
    return k.omw(raw, hf, hb, n2, bg, k.omw(raw, hf, hb, n1, bg))


def rfft_rlstep(psf, y, im, grain, *, backend, threads=None):
    k = kernel(y.shape, grain, backend, threads)
    return k.rlstep(y, im, k.otf(psf))


def rfft_rl_loop(psf, im, n, grain, *, backend, threads=None):
    """n RL steps in a row, J <- rlstep(J, I, H) from J = I, in the kernel's internal layout"""
    k = kernel(im.shape, grain, backend, threads)
    h = k.otf(psf)
    i = k.to_internal(im)
    j = i
    for _ in range(int(n)):
        j = k.rlstep(j, i, h, internal=True)
    return k.from_internal(j)


def rfft_fft1(x, nt, grain, *, backend, threads=None):
    return D.fft1(x, backend=backend, threads=threads, matlab_parallel=float(grain) > 0)


# MATLAB's code of the composed cases (in{} as listed in args)
OMW_CODE = (
    "sz = size(in{2}); H = decon_lucy_rfft_mex('otf', in{1}, sz, 4, in{5}); "
    "Hb = decon_lucy_rfft_mex('otf', in{3}, sz, 4, in{5}); "
    "out = {decon_lucy_rfft_mex('omw', in{2}, H, Hb, in{4}, in{6}, 4, in{5})};"
)
OMW_CONT_CODE = (
    "sz = size(in{2}); H = decon_lucy_rfft_mex('otf', in{1}, sz, 4, in{6}); "
    "Hb = decon_lucy_rfft_mex('otf', in{3}, sz, 4, in{6}); "
    "J = decon_lucy_rfft_mex('omw', in{2}, H, Hb, in{4}, in{7}, 4, in{6}); "
    "out = {decon_lucy_rfft_mex('omw', in{2}, H, Hb, in{5}, in{7}, J, 4, in{6})};"
)
RLSTEP_CODE = (
    "H = decon_lucy_rfft_mex('otf', in{1}, size(in{2}), 4, in{4}); "
    "out = {decon_lucy_rfft_mex('rlstep', in{2}, in{3}, H, 4, in{4})};"
)
RL_LOOP_CODE = (
    "H = decon_lucy_rfft_mex('otf', in{1}, size(in{2}), 4, in{4}); J = in{2}; "
    "for k = 1 : in{3}, J = decon_lucy_rfft_mex('rlstep', J, in{2}, H, 4, in{4}); end; out = {J};"
)

# MATLAB sizes [ny nx nz]: DeconRFFTTest's, and sizes with 128
RFFT_SIZES = [
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
    [128, 33, 20],
    [31, 128, 18],
    [24, 30, 128],
]
GRAIN0_SIZES = [[66, 64, 33], [64, 230], [10, 23, 64], [31, 128, 18]]


def size_name(s) -> str:
    return "x".join(str(v) for v in s)


def two_backends(case: Case) -> list[Case]:
    """the case's exact and permissive Python sides (one MATLAB reference)"""
    out = []
    for be in ("exact", "permissive"):
        kw = dict(case.py_kwargs or {})
        kw["backend"] = be
        spec = case.compare
        if be == "permissive":
            spec = PERMISSIVE_OTF if case.note == "otf" else PERMISSIVE_TOL
        out.append(dataclasses.replace(case, id=f"{case.id}.{be}", py_kwargs=kw, compare=spec))
    return out


CASES: list[Case] = []
_seed = 0
for _ms in RFFT_SIZES:
    _seed += 1
    _sn = size_name(_ms)
    _inp = {
        "psf": V(lambda m=_ms: make_psf(m, "f")),
        "psfb": V(lambda m=_ms: make_psf(m, "b")),
        "u16": V(lambda m=_ms, s=_seed: blurred_beads(m, s, "uint16")),
        "f32": V(lambda m=_ms, s=_seed: blurred_beads(m, s, "single")),
    }
    for _g in [32] + ([0] if _ms in GRAIN0_SIZES else []):
        CASES += two_backends(
            Case(
                f"rfft.otf.{_sn}.g{_g}",
                matlab="decon_lucy_rfft_mex",
                args=["otf", Vol(In("psf")), [float(v) for v in _ms], 4, _g],
                python=rfft_otf,
                py_args=[Vol(In("psf")), PerAxis([float(v) for v in _ms], "yxz"), 4, _g],
                inputs={"psf": _inp["psf"]},
                note="otf",
            )
        )
        CASES += two_backends(
            Case(
                f"rfft.omw2.uint16.{_sn}.g{_g}",
                code=OMW_CODE,
                args=[Vol(In("psf")), Vol(In("u16")), Vol(In("psfb")), 2, _g, 100],
                python=rfft_omw,
                inputs={k: _inp[k] for k in ("psf", "psfb", "u16")},
                out_roles=("vol",),
            )
        )
        _yinp = {
            "psf": _inp["psf"],
            "y": V(
                lambda m=_ms, s=_seed: (blurred_beads(m, s, "single") - np.float32(99.5)).astype(
                    np.float32
                )
            ),
            "i": V(
                lambda m=_ms, s=_seed: (blurred_beads(m, s, "single") - np.float32(100)).astype(
                    np.float32
                )
            ),
        }
        CASES += two_backends(
            Case(
                f"rfft.rlstep.{_sn}.g{_g}",
                code=RLSTEP_CODE,
                args=[Vol(In("psf")), Vol(In("y")), Vol(In("i")), _g],
                python=rfft_rlstep,
                inputs=_yinp,
                out_roles=("vol",),
            )
        )
    CASES += two_backends(
        Case(
            f"rfft.omw2.single.{_sn}",
            code=OMW_CODE,
            args=[Vol(In("psf")), Vol(In("f32")), Vol(In("psfb")), 2, 32, 100.3],
            python=rfft_omw,
            inputs={k: _inp[k] for k in ("psf", "psfb", "f32")},
            out_roles=("vol",),
        )
    )


def below_background(msize, seed: int) -> np.ndarray:
    """blurred beads with a slab below the background (J is +0 there, and products with it -0
    where the back projection is negative: the mex's update makes them +0)"""
    a = blurred_beads(msize, seed, "uint16")
    if a.ndim == 3:
        a[: a.shape[0] // 3, :, : a.shape[2] // 2] = 60
    else:
        a[: a.shape[0] // 3, :] = 60
    return a


for _ms in ([48, 46, 36], [64, 50]):
    _sn = size_name(_ms)
    CASES += two_backends(
        Case(
            f"rfft.omw2.negbp.{_sn}",
            code=OMW_CODE,
            args=[Vol(In("psf")), Vol(In("u16")), Vol(In("psfb")), 2, 32, 100],
            python=rfft_omw,
            inputs={
                "psf": V(lambda m=_ms: make_psf(m, "f")),
                "psfb": V(lambda m=_ms: make_psf(m, "neg")),
                "u16": V(lambda m=_ms: below_background(m, 90)),
            },
            out_roles=("vol",),
        )
    )

# iterations 0, 1, 5, 50, an estimate to continue from, RL steps in a row
for _ms, _iters in [([48, 46, 36], [0, 1, 5]), ([30, 26, 22], [50]), ([64, 50], [0, 1, 5, 50])]:
    _sn = size_name(_ms)
    _inp = {
        "psf": V(lambda m=_ms: make_psf(m, "f")),
        "psfb": V(lambda m=_ms: make_psf(m, "b")),
        "u16": V(lambda m=_ms: blurred_beads(m, 77, "uint16")),
    }
    for _n in _iters:
        CASES += two_backends(
            Case(
                f"rfft.omw{_n}.uint16.{_sn}",
                code=OMW_CODE,
                args=[Vol(In("psf")), Vol(In("u16")), Vol(In("psfb")), _n, 32, 100.3],
                python=rfft_omw,
                inputs=_inp,
                out_roles=("vol",),
                marks=("slow",) if _n == 50 and len(_ms) == 3 else (),
            )
        )
    CASES += two_backends(
        Case(
            f"rfft.omw.continued.{_sn}",
            code=OMW_CONT_CODE,
            args=[Vol(In("psf")), Vol(In("u16")), Vol(In("psfb")), 2, 3, 32, 100],
            python=rfft_omw_continued,
            inputs=_inp,
            out_roles=("vol",),
        )
    )
    CASES += two_backends(
        Case(
            f"rfft.rl10.{_sn}",
            code=RL_LOOP_CODE,
            args=[Vol(In("psf")), Vol(In("i")), 10, 32],
            python=rfft_rl_loop,
            inputs={
                "psf": _inp["psf"],
                "i": V(
                    lambda m=_ms: np.maximum(
                        blurred_beads(m, 78, "single") - np.float32(100), np.float32(0)
                    )
                ),
            },
            out_roles=("vol",),
        )
    )

# 'fft1': MATLAB's groups of lines in both modes, and this session's fft (its parallel mode)
for _n in [2000, 222, 210]:
    _x = V(lambda n=_n: CounterRNG(n, "fft1").uniform((64, n), np.float32))
    for _g in (32, 0):
        CASES += two_backends(
            Case(
                f"rfft.fft1.{_n}.g{_g}",
                matlab="decon_lucy_rfft_mex",
                args=["fft1", In("x"), 4, _g],
                python=rfft_fft1,
                py_args=[In("x"), 4, _g],
                inputs={"x": _x},
                note="otf",
            )
        )
    CASES += two_backends(
        Case(
            f"rfft.fft1.{_n}.matlab_fft",
            code="Y = fft(in{1}, [], 1); out = {Y(1 : end / 2 + 1, :)};",
            args=[In("x")],
            python=rfft_fft1,
            py_args=[In("x"), 4, 32],
            inputs={"x": _x},
            note="otf",
        )
    )


# the mex's argument errors (DeconRFFTTest/mexRejectsInvalidInputs)
def _err_inputs():
    im = to_matlab(np.ones((10, 16, 12), np.uint16) * np.uint16(150))
    return {"im": V(lambda: im), "y": V(lambda: im.astype(np.float32))}


_RFFT_ID = "decon_lucy_rfft_mex:inputError"
_ERR_H = "H = decon_lucy_rfft_mex('otf', single(ones(3, 3, 3) / 27), [16, 12, 10], 1, 32); "


def rfft_error(im, y, *, which: str, backend="permissive", threads=None):
    """the Python counterparts of the error cases (y: single data, im: uint16)"""
    k = D.LucyRFFT(im.shape, backend=backend, threads=threads)
    h = k.otf(np.full((3, 3, 3), np.float32(1 / 27)))
    if which == "short_otf":
        return k.omw(im, h[..., :-1], h, 2, 100)
    if which == "double":
        return k.omw(im.astype(np.float64), h, h, 2, 100)
    if which == "real_otf":
        return k.omw(im, h.real.copy(), h, 2, 100)
    if which == "rl_uint16":
        return k.rlstep(y, im, h)
    if which == "rl_sizes":
        return k.rlstep(y, y[:-1], h)
    if which == "one_row":
        return D.LucyRFFT((y.shape[0], 1, y.shape[2]))
    if which == "otf_row":
        return D.otf(np.ones((1, 1), np.float32), (10, 1, 12))
    if which == "fft1_volume":
        return D.fft1(y)
    if which == "j0_size":
        return k.omw(im, h, h, 2, 100, y[:-1])
    if which == "j0_double":
        return k.omw(im, h, h, 2, 100, y.astype(np.float64))
    if which == "psf_larger":
        return D.otf(np.ones((5, 5, 5), np.float32), (10, 4, 12))
    raise AssertionError(which)


for _name, _code, _py in [
    (
        "omw.short_otf",
        "decon_lucy_rfft_mex('omw', in{1}, H(:, :, 1 : end - 1), H, 2, 100, 1, 32)",
        "short_otf",
    ),
    ("omw.double", "decon_lucy_rfft_mex('omw', double(in{1}), H, H, 2, 100, 1, 32)", "double"),
    ("omw.real_otf", "decon_lucy_rfft_mex('omw', in{1}, real(H), H, 2, 100, 1, 32)", "real_otf"),
    ("rlstep.uint16", "decon_lucy_rfft_mex('rlstep', in{2}, in{1}, H, 1, 32)", "rl_uint16"),
    (
        "rlstep.sizes",
        "decon_lucy_rfft_mex('rlstep', in{2}, in{2}(:, :, 1 : end - 1), H, 1, 32)",
        "rl_sizes",
    ),
    (
        "rlstep.one_row",
        "decon_lucy_rfft_mex('rlstep', in{2}(1, :, :), in{2}(1, :, :), H, 1, 32)",
        "one_row",
    ),
    ("otf.one_row", "decon_lucy_rfft_mex('otf', single(1), [1, 12, 10], 1, 32)", "otf_row"),
    ("fft1.volume", "decon_lucy_rfft_mex('fft1', in{2}, 1, 32)", "fft1_volume"),
    ("unknown", "decon_lucy_rfft_mex('unknown', in{1}, H, 2, 100, 1, 32)", None),
    (
        "omw.j0_size",
        "decon_lucy_rfft_mex('omw', in{1}, H, H, 2, 100, in{2}(:, :, 1 : end - 1), 1, 32)",
        "j0_size",
    ),
    (
        "omw.j0_double",
        "decon_lucy_rfft_mex('omw', in{1}, H, H, 2, 100, double(in{2}), 1, 32)",
        "j0_double",
    ),
    (
        "otf.psf_larger",
        "decon_lucy_rfft_mex('otf', single(ones(5, 5, 5)), [4, 12, 10], 1, 32)",
        "psf_larger",
    ),
]:
    CASES.append(
        Case(
            f"rfft.error.{_name}",
            code=_ERR_H + _code + ";",
            args=[Vol(In("im")), Vol(In("y"))],
            python=rfft_error if _py else None,
            py_kwargs={"which": _py} if _py else None,
            inputs=_err_inputs(),
            expect_error=_RFFT_ID
            if _py != "psf_larger"
            else "decon_lucy_rfft_mex:outSizeIsSmallerThanPsfSize",
            py_error=D.DeconKernelError,
        )
    )


# ---------------------------------------------------------------------------------------------
# accel: decon_lucy_accel_mex, decon_lucy_accel, decon_sum_double
# ---------------------------------------------------------------------------------------------

ACCEL_SZ = [47, 45, 49]  # 103,635 elements: a last block of 4096 and a last group of four not full
SUM_SIZES = [1, 43, 44, 4097, 88999, 89000, 89001, 4096 * 30, 4096 * 30 + 3, 250047]


def wide_range(n: int, cls: str, seed: int) -> np.ndarray:
    """DeconAccelTest.wideRange: randn .* 2 .^ randi([-40, 40]), a column"""
    rng = CounterRNG(seed, "wide")
    v = np.ldexp(rng.normal(n), rng.integers(-40, 40, n).astype(np.int32))
    return v.astype(np.float32 if cls == "single" else np.float64).reshape(-1, 1)


def accel_inputs(n: int, cls: str, non_finite: bool, seed: int):
    """DeconAccelTest.accelInputs: J_2, Y, J_4, J_3 columns with -0, +0, subnormal numbers (and
    NaN, Inf, -Inf if non_finite)"""
    dt = np.float32 if cls == "single" else np.float64
    rng = CounterRNG(seed, f"accel/{cls}")
    j2 = (rng.uniform(n) * 1000).astype(dt)
    j3 = (rng.uniform(n) * 1000).astype(dt)
    j4 = (rng.normal(n) * 10).astype(dt)
    y = (j2 - dt(0.6) * j4 + rng.normal(n).astype(dt)).astype(dt)
    sp = np.array([0, -0.0, NAN, INF, -INF, 1e-40, -1e-40, 5, -0.0, 0], dtype=np.float64)
    if not non_finite:
        sp[2:5] = [0, -0.0, 3]
    sp = sp.astype(dt)
    j2[:10] = sp
    j3[:10] = sp[[1, 0, 3, 4, 2, 0, 1, 8, 1, 0]]
    y[10:20] = sp
    j4[20:30] = sp
    j2[n - 10 :] = sp
    j3[n - 10 :] = -0.0
    return [a.reshape(-1, 1) for a in (j2, y, j4, j3)]


def in_shape(col: np.ndarray, like: np.ndarray) -> np.ndarray:
    """MATLAB's column (or array) as a (z, y, x) array of like's shape in MATLAB's linear order"""
    msize = size_to_matlab(like.shape) if like.ndim >= 2 else (1, like.size)
    return from_matlab(np.reshape(np.asarray(col).reshape(-1, order="F"), msize, order="F"))


def as_column(a: np.ndarray, matlab_layout: bool) -> np.ndarray:
    if matlab_layout or a.ndim < 2:
        return a.reshape(-1, 1, order="F")
    return np.ascontiguousarray(matlab_order(a)).reshape(-1, 1)


def accel_lambda_case(j2, y, j4, *rest, layout="zyx", threads=None):
    """decon_lucy_accel('lambda', J_2, Y, J_4) or the mex: [lambda, J_4, s]"""
    ml = layout == "matlab"
    if ml:
        a2 = to_matlab(j2) if j2.ndim >= 2 else j2
        ay, a4 = (to_matlab(in_shape(v, j2)) if j2.ndim >= 2 else in_shape(v, j2) for v in (y, j4))
    else:
        a2, ay, a4 = j2, in_shape(y, j2), in_shape(j4, j2)
    r = D.accel_lambda(a2, ay, a4, matlab_layout=ml, threads=threads)
    return [r.lam, as_column(r.j4, ml), np.array([[r.sums[0], r.sums[1]]])]


def accel_lambda_mex(method, j2, y, j4, nt, *, layout="zyx", threads=None):
    return accel_lambda_case(j2, y, j4, layout=layout, threads=threads)


def accel_predict_case(*args, layout="zyx", threads=None):
    """decon_lucy_accel('predict', J_2, J_3, lambda) or the mex (the method first)"""
    if isinstance(args[0], str):
        args = args[1:]
    j2, j3, lam = args[:3]
    if layout == "matlab":
        return from_matlab(D.accel_predict(to_matlab(j2), to_matlab(j3), lam, threads=threads))
    return D.accel_predict(j2, j3, lam, threads=threads)


def sum_case(*args, layout="zyx", threads=None):
    """decon_sum_double(x), or the mex's 'sum' (the method first)"""
    x = args[1] if isinstance(args[0], str) else args[0]
    if layout == "matlab":
        return D.sum_double(to_matlab(x) if x.ndim >= 2 else x, matlab_layout=True, threads=threads)
    return D.sum_double(x, threads=threads)


def with_layouts(case: Case, layouts=("zyx", "matlab")) -> list[Case]:
    out = []
    for lay in layouts:
        kw = dict(case.py_kwargs or {})
        kw["layout"] = lay
        out.append(
            dataclasses.replace(
                case, id=case.id + (".mlayout" if lay == "matlab" else ""), py_kwargs=kw
            )
        )
    return out


def special_vectors(n: int, cls: str, seed: int) -> list[np.ndarray]:
    """DeconAccelTest.specialVectors, plus NaNs of two payloads, and +Inf and -Inf in other blocks
    (their sum is the default NaN) with a NaN"""
    dt = np.float32 if cls == "single" else np.float64
    z = -np.zeros((n, 1), dt)
    p = z.copy()
    p[4096] = 0
    x = wide_range(n, cls, seed)
    xn = x.copy()
    xn[4999] = NAN
    xi = x.copy()
    xi[76] = INF
    xi[8999] = INF
    xm = x.copy()
    xm[2] = -INF
    rng = CounterRNG(seed, "subnormal")
    xd = (rng.normal(n).astype(dt) * np.finfo(dt).tiny / dt(64)).reshape(-1, 1).astype(dt)
    if cls == "single":
        nan_a, nan_b = (
            np.array([v], np.uint32).view(np.float32)[0] for v in (0x7FC00001, 0xFFC12345)
        )
    else:
        nan_a, nan_b = (
            np.array([v], np.uint64).view(np.float64)[0]
            for v in (0x7FF8000000000001, 0xFFF8123456789ABC)
        )
    x2 = x.copy()
    x2[100] = nan_a
    x2[50000] = nan_b
    x3 = x.copy()
    x3[3] = nan_b
    x3[4100] = nan_a
    x3[n - 2] = nan_a
    x4 = x.copy()
    x4[10] = INF
    x4[9000] = -INF
    x4[70000] = nan_a
    x5 = x.copy()
    x5[n - 1] = nan_b  # in the last, short block's tail
    x5[n - 6] = nan_a
    x6 = x.copy()
    x6[200] = nan_a  # two NaNs in one lane of one block, then one in another lane
    x6[204] = nan_b
    x6[201] = nan_b
    return [z, p, np.zeros((n, 1), dt), xn, xi, xm, xd, x, x2, x3, x4, x5, x6]


for _cls in ("single", "double"):
    _n = math.prod(ACCEL_SZ)
    for _c in range(4):
        _nf = _c == 2

        def _make(c=_c, cls=_cls, nf=_nf, n=_n):
            j2, y, j4, _ = accel_inputs(n, cls, nf, 10 + c)
            if c == 1:
                j4 = wide_range(n, cls, 20)
                y = (
                    j2
                    - (0.5 * j4).astype(j4.dtype)
                    - (0.1 * wide_range(n, cls, 21)).astype(j4.dtype)
                ).astype(j4.dtype)
            elif c == 3:
                j4[n - 3 :, 0] = np.array([3e9, -5e9, 7e9]).astype(j4.dtype)
                y[n - 3 :] = j2[n - 3 :] - (0.5 * j4[n - 3 :]).astype(j4.dtype)
            return np.reshape(j2, ACCEL_SZ, order="F"), y, j4

        _inp = {
            "j2": V(lambda m=_make: m()[0]),
            "y": V(lambda m=_make: m()[1]),
            "j4": V(lambda m=_make: m()[2]),
        }
        CASES += with_layouts(
            Case(
                f"accel.lambda.mex.{_cls}.c{_c}",
                matlab="decon_lucy_accel_mex",
                args=["lambda", Vol(In("j2")), In("y"), In("j4"), 3],
                nargout=3,
                python=accel_lambda_mex,
                inputs=_inp,
            )
        )
        CASES += with_layouts(
            Case(
                f"accel.lambda.fn.{_cls}.c{_c}",
                code="[l, J] = decon_lucy_accel('lambda', in{1}, in{2}, in{3}); out = {l, J};",
                args=[Vol(In("j2")), In("y"), In("j4")],
                nargout=2,
                python=lambda *a, **k: accel_lambda_case(*a, **k)[:2],
                inputs=_inp,
            )
        )
    # predict with lambdas of DeconAccelTest
    _pin = {
        "j2": V(
            lambda cls=_cls: np.reshape(
                accel_inputs(math.prod(ACCEL_SZ), cls, True, 30)[0], ACCEL_SZ, order="F"
            )
        ),
        "j3": V(
            lambda cls=_cls: np.reshape(
                accel_inputs(math.prod(ACCEL_SZ), cls, True, 30)[3], ACCEL_SZ, order="F"
            )
        ),
    }
    for _k, _lam in enumerate([0.0, -0.0, 1.0, 0.5, 0.123456789012345678, 1 - 2.0**-52, 2.0**-60]):
        CASES += with_layouts(
            Case(
                f"accel.predict.mex.{_cls}.l{_k}",
                matlab="decon_lucy_accel_mex",
                args=["predict", Vol(In("j2")), Vol(In("j3")), _lam, 3],
                python=accel_predict_case,
                inputs=_pin,
                out_roles=("vol",),
            )
        )
    # below 89000 elements: MATLAB's code (its own sums)
    _sm = {
        "j2": V(lambda cls=_cls: accel_inputs(88999, cls, False, 40)[0]),
        "j3": V(lambda cls=_cls: accel_inputs(88999, cls, False, 40)[3]),
        "j4": V(lambda cls=_cls: wide_range(88999, cls, 41)),
        "y": V(
            lambda cls=_cls: (
                accel_inputs(88999, cls, False, 40)[0]
                - (0.5 * wide_range(88999, cls, 41)).astype(
                    np.float32 if cls == "single" else np.float64
                )
                - (0.1 * wide_range(88999, cls, 42)).astype(
                    np.float32 if cls == "single" else np.float64
                )
            )
        ),
    }
    CASES += with_layouts(
        Case(
            f"accel.lambda.fn.small.{_cls}",
            code="[l, J] = decon_lucy_accel('lambda', in{1}, in{2}, in{3}); out = {l, J};",
            args=[Vol(In("j2")), In("y"), In("j4")],
            nargout=2,
            python=lambda *a, **k: accel_lambda_case(*a, **k)[:2],
            inputs={k: _sm[k] for k in ("j2", "y", "j4")},
        )
    )
    CASES += with_layouts(
        Case(
            f"accel.predict.fn.small.{_cls}",
            code="out = {decon_lucy_accel('predict', in{1}, in{2}, 0.3)};",
            args=[Vol(In("j2")), Vol(In("j3"))],
            python=lambda j2, j3, **k: accel_predict_case(j2, j3, 0.3, **k),
            inputs={k: _sm[k] for k in ("j2", "j3")},
            out_roles=("vol",),
        )
    )
    # decon_sum_double around the switches; 3D arrays summed as x(:)
    for _ns in SUM_SIZES:
        CASES += with_layouts(
            Case(
                f"accel.sum.{_cls}.n{_ns}",
                matlab="decon_sum_double",
                args=[Vol(In("x"))],
                python=sum_case,
                inputs={"x": V(lambda n=_ns, cls=_cls: wide_range(n, cls, 50 + n % 97))},
            )
        )
    CASES += with_layouts(
        Case(
            f"accel.sum.{_cls}.volume",
            matlab="decon_sum_double",
            args=[Vol(In("x"))],
            python=sum_case,
            inputs={
                "x": V(
                    lambda cls=_cls: np.reshape(
                        wide_range(250047, cls, 60), [63, 63, 63], order="F"
                    )
                )
            },
        )
    )
    for _k in range(13):
        CASES += with_layouts(
            Case(
                f"accel.sum.mex.{_cls}.special{_k}",
                matlab="decon_lucy_accel_mex",
                args=["sum", Vol(In("x")), 2],
                python=sum_case,
                inputs={
                    "x": V(lambda k=_k, cls=_cls: special_vectors(math.prod(ACCEL_SZ), cls, 70)[k])
                },
            )
        )
    CASES += with_layouts(
        Case(
            f"accel.sum.mex.{_cls}.volume_nan",
            matlab="decon_lucy_accel_mex",
            args=["sum", Vol(In("x")), 2],
            python=sum_case,
            inputs={
                "x": V(
                    lambda cls=_cls: np.reshape(
                        special_vectors(math.prod(ACCEL_SZ), cls, 71)[9], ACCEL_SZ, order="F"
                    )
                )
            },
        )
    )

# lambda's clamp (DeconAccelTest/lambdaIsClampedAsMatlabDoes), single columns


def clamp_inputs(case: str):
    n = math.prod(ACCEL_SZ)
    rng = CounterRNG(80, case)
    j2 = (rng.uniform(n, np.float32) * np.float32(100)).reshape(-1, 1)
    j4 = (rng.uniform(n, np.float32) - np.float32(0.5)).reshape(-1, 1)
    one = np.ones((n, 1), np.float32)
    if case == "aboveOne":
        return j2, j2 - np.float32(2), one
    if case == "negative":
        return j2, j2 + np.float32(1), one
    if case == "nan":
        j4[6] = INF
        return j2, j2.copy(), j4
    if case == "negativeZero":
        return j2, j2.copy(), -one
    if case == "zeroDenominator":
        return j2, j2 - np.float32(1), np.zeros((n, 1), np.float32)
    if case == "inf":
        j2[4] = INF
        y = j2 - np.float32(1)
        y[4] = 0
        return j2, y, one
    return j2, (j2 - np.float32(0.5) * j4).astype(np.float32), j4


for _lc in ["aboveOne", "negative", "nan", "negativeZero", "zeroDenominator", "inf", "inside"]:
    CASES += with_layouts(
        Case(
            f"accel.clamp.{_lc}",
            code="[l, J] = decon_lucy_accel('lambda', in{1}, in{2}, in{3}); out = {l, J};",
            args=[Vol(In("j2")), In("y"), In("j4")],
            nargout=2,
            python=lambda *a, **k: accel_lambda_case(*a, **k)[:2],
            inputs={
                k: V(lambda c=_lc, i=i: clamp_inputs(c)[i]) for i, k in enumerate(["j2", "y", "j4"])
            },
        )
    )


# the mex's argument errors (DeconAccelTest/mexRejectsInvalidInputs)
def accel_error(which: str, threads=None):
    a = CounterRNG(90).uniform((10, 1), np.float32)
    if which == "class":
        return D.accel_lambda(a, a.astype(np.float64), a)
    if which == "size":
        return D.accel_lambda(a, a[:9], a)
    if which == "complex":
        return D.accel_lambda(a.astype(np.complex64), a, a)
    if which == "uint16":
        return D.accel_lambda(*(a.astype(np.uint16),) * 3)
    if which == "lambda_single":
        return D.accel_predict(a, a, np.float32(0.5))
    if which == "lambda_vector":
        return D.accel_predict(a, a, np.array([0.5, 0.5]))
    if which == "predict_class":
        return D.accel_predict(a, a.astype(np.float64), 0.5)
    raise AssertionError(which)


_ACCEL_ID = "decon_lucy_accel_mex:inputError"
for _w, _code in [
    ("class", "decon_lucy_accel_mex('lambda', a, double(a), a, 1)"),
    ("size", "decon_lucy_accel_mex('lambda', a, a(1 : 9), a, 1)"),
    ("complex", "decon_lucy_accel_mex('lambda', complex(a), a, a, 1)"),
    ("uint16", "decon_lucy_accel_mex('lambda', uint16(a), uint16(a), uint16(a), 1)"),
    ("lambda_single", "decon_lucy_accel_mex('predict', a, a, single(0.5), 1)"),
    ("lambda_vector", "decon_lucy_accel_mex('predict', a, a, [0.5, 0.5], 1)"),
    ("predict_class", "decon_lucy_accel_mex('predict', a, double(a), 0.5, 1)"),
]:
    CASES.append(
        Case(
            f"accel.error.{_w}",
            code=f"a = in{{1}}; {_code};",
            args=[In("a")],
            python=accel_error,
            py_args=[],
            py_kwargs={"which": _w},
            inputs={"a": V(lambda: CounterRNG(90).uniform((10, 1), np.float32))},
            expect_error=_ACCEL_ID,
            py_error=D.DeconKernelError,
        )
    )
CASES.append(
    Case(
        "accel.sum.mex.empty",
        matlab="decon_lucy_accel_mex",
        args=["sum", np.zeros((0, 1), np.float32), 2],
        python=lambda *a, **k: D.sum_double(np.zeros(0, np.float32)),
    )
)


# ---------------------------------------------------------------------------------------------
# erosion: decon_mask_edge_erosion (its mex path) and decon_mask_edge_erosion_mex
# ---------------------------------------------------------------------------------------------

EROSION_CODE = "A = in{1}; out = {arrayfun(@(E) decon_mask_edge_erosion(A, E), 1 : 41, 'unif', 0)};"
MASK_SIZES = {
    "odd3D": [37, 41, 23],
    "small": [5, 5, 5],
    "twoSlices": [17, 19, 2],
    "oneRow": [1, 30, 20],
    "oneColumn": [30, 1, 20],
    "image2D": [45, 47],
}


def make_mask(shape: str, msize, seed: int = 0) -> np.ndarray:
    """DeconMaskErosionTest.makeMask in MATLAB's form: random, skewed data's zero triangles along
    x-z, a thin plane, one true or one false voxel, blobs"""
    sz3 = list(msize) + [1] * (3 - len(msize))
    rng = CounterRNG(seed, f"mask/{shape}/{size_name(msize)}")
    py = tuple(int(v) for v in ((sz3[2], sz3[0], sz3[1]) if len(msize) == 3 else msize))
    if shape == "random":
        a = rng.uniform(py) > 0.3
    elif shape == "skewed":
        z = np.arange(1, sz3[2] + 1, dtype=np.float64).reshape(-1, 1, 1)
        x = np.arange(1, sz3[1] + 1, dtype=np.float64).reshape(1, 1, -1)
        a = (x > 0.6 * z + sz3[1] * 0.05) & (x < 0.6 * z + sz3[1] * 0.5)
        a = np.broadcast_to(a, (sz3[2], sz3[0], sz3[1])).reshape(py)
    elif shape == "thin":
        a = np.zeros(py, bool)
        c = max(1, math.floor(sz3[1] / 2 + 0.5)) - 1
        a[..., c : c + min(2, sz3[1])] = True
    elif shape == "single":
        a = np.zeros(py, bool)
        mo = np.zeros(math.prod(sz3), bool)
        mo[math.ceil(mo.size / 2) - 1] = True
        return np.reshape(mo, msize, order="F")
    elif shape == "singleFalse":
        mo = np.ones(math.prod(sz3), bool)
        mo[math.ceil(mo.size / 2) - 1] = False
        return np.reshape(mo, msize, order="F")
    else:  # blobs: imdilate(rand > 0.99, ones(3, 3))
        from petakit5d.matlabcompat.morphology import imdilate

        a = imdilate(rng.uniform(py) > 0.99, np.ones((3, 3), bool))
    return to_matlab(np.ascontiguousarray(a))


@functools.cache
def disk_decomposition(r: int):
    """decon_mask_edge_erosion.m's disk_decomposition: the neighborhoods of strel('disk', r)'s
    decomposition and morphop's padding for them (before and after swapped for erosion)"""
    from petakit5d.matlabcompat.morphology import strel

    se = strel("disk", r)
    seq = se.decompose()
    nh = [np.asarray(s.nhood) for s in seq]
    pu = pl = (0, 0)
    if len(seq) > 1 and not se.is_decomposition_orthogonal():
        before, after = se.getpadsize()
        pu, pl = tuple(after[:2]), tuple(before[:2])
    return nh, pu, pl


def edge_erosion(mask: np.ndarray, e: int, threads=None) -> np.ndarray:
    """decon_mask_edge_erosion(mask, EdgeErosion)'s mex path for masks that aren't full (R03 ports
    the function; this is its kernel call)"""
    m = mask if mask.dtype == np.bool_ else mask > 0
    nh, pu, pl = disk_decomposition(e - 1)
    rz = e - 1 if m.ndim == 3 else 0
    return D.mask_edge_erosion(m, nh, pu, pl, rz, True, threads=threads)


def erosion_all(mask, *, threads=None):
    return [edge_erosion(mask, e, threads) for e in range(1, 42)]


for _shape in ["random", "skewed", "thin", "single", "singleFalse", "blobs"]:
    for _sname, _ms in MASK_SIZES.items():
        if make_mask(_shape, _ms).all():
            # a full mask (thin with one column) has its edges set without the mex (R03's)
            continue
        CASES.append(
            Case(
                f"erosion.{_shape}.{_sname}",
                code=EROSION_CODE,
                args=[Vol(In("a"))],
                python=erosion_all,
                inputs={"a": V(lambda s=_shape, m=_ms: make_mask(s, m))},
                out_roles=(("each", "vol"),),
            )
        )


def demo_mask():
    """DeconMaskErosionTest/demoSettingOnASkewedVolume: EdgeErosion 20, zero borders as in skewed
    data, more slices than threads and a few isolated holes"""
    a = make_mask("skewed", [150, 200, 90])
    for y, x, z in zip([70, 20, 100], [70, 90, 110], [40, 60, 75], strict=True):
        a[y - 1, x - 1, z - 1] = False
    return a


CASES.append(
    Case(
        "erosion.demo",
        matlab="decon_mask_edge_erosion",
        args=[Vol(In("a")), 20],
        python=lambda a, e, threads=None: edge_erosion(a, e, threads),
        inputs={"a": V(demo_mask)},
        out_roles=("vol",),
    )
)
CASES.append(
    Case(
        "erosion.demo.numeric",
        matlab="decon_mask_edge_erosion",
        args=[Vol(In("a")), 20],
        python=lambda a, e, threads=None: edge_erosion(a, e, threads),
        inputs={
            "a": V(
                lambda: (
                    demo_mask().astype(np.uint16)
                    * to_matlab(CounterRNG(3).integers(1, 1000, (90, 150, 200), np.uint16))
                ).astype(np.uint16)
            )
        },
        out_roles=("vol",),
    )
)


def random_nhoods(trial: int):
    """DeconMaskErosionTest/mexEqualsImerodeForOtherNeighborhoods: 1 to 3 random (asymmetric,
    even-sized) neighborhoods, a random mask, and morphop's padding for them"""
    from petakit5d.matlabcompat.morphology import Strel

    rng = CounterRNG(trial, "nhoods")
    ns = int(rng.integers(1, 3, 1)[0])
    nh = []
    for _ in range(ns):
        sh = tuple(int(v) for v in rng.integers(1, 6, 2))
        e = rng.uniform(sh) > 0.5
        e.reshape(-1)[int(rng.integers(0, e.size - 1, 1)[0])] = True
        nh.append(e)
    msz = [int(v) for v in rng.integers(1, 30, 2)] + [int(rng.integers(2, 4, 1)[0])]
    a = rng.uniform((msz[2], msz[0], msz[1])) > 0.2
    pu = pl = (0, 0)
    if ns > 1:
        se = Strel(np.ones((1, 1), bool), [Strel(x) for x in nh])
        if not se.is_decomposition_orthogonal():
            before, after = se.getpadsize()
            pu, pl = tuple(after[:2]), tuple(before[:2])
    return to_matlab(np.ascontiguousarray(a)), nh, [float(v) for v in pu], [float(v) for v in pl]


def erosion_mex(a, nh, pu, pl, rz, clear, nt, *, threads=None):
    return D.mask_edge_erosion(a, nh, pu, pl, rz, clear, threads=threads)


for _t in range(24):
    _a, _nh, _pu, _pl = random_nhoods(_t)
    for _nt in (1, 3):
        CASES.append(
            Case(
                f"erosion.mex.random.t{_t}.n{_nt}",
                matlab="decon_mask_edge_erosion_mex",
                args=[Vol(In("a")), _nh, _pu, _pl, 0, False, _nt],
                python=erosion_mex,
                inputs={"a": V(lambda t=_t: random_nhoods(t)[0])},
                out_roles=("vol",),
            )
        )
CASES.append(
    Case(
        "erosion.mex.zline",
        matlab="decon_mask_edge_erosion_mex",
        args=[Vol(In("a")), [np.array([[True]])], [0.0, 0.0], [0.0, 0.0], 4, False, 2],
        python=erosion_mex,
        inputs={"a": V(lambda: to_matlab(CounterRNG(5).uniform((40, 20, 30)) > 0.1))},
        out_roles=("vol",),
    )
)
CASES.append(
    Case(
        "erosion.mex.numeric_nhood",
        matlab="decon_mask_edge_erosion_mex",
        args=[Vol(In("a")), [np.array([[0.0, 2.0, 1.0]])], [0.0, 0.0], [0.0, 0.0], 0, False, 1],
        python=erosion_mex,
        inputs={"a": V(lambda: to_matlab(CounterRNG(6).uniform((5, 10, 12)) > 0.5))},
        out_roles=("vol",),
    )
)
CASES.append(
    Case(
        "erosion.mex.empty",
        matlab="decon_mask_edge_erosion_mex",
        args=[Vol(In("a")), [np.array([[True]])], [0.0, 0.0], [0.0, 0.0], 1, True, 1],
        python=erosion_mex,
        inputs={"a": V(lambda: np.zeros((0, 3, 2), bool))},
        out_roles=("vol",),
    )
)
CASES.append(
    Case(
        "erosion.mex.clear_border_pads",
        matlab="decon_mask_edge_erosion_mex",
        args=[
            Vol(In("a")),
            [np.ones((3, 5), bool), np.ones((4, 2), bool)],
            [2.0, 1.0],
            [0.0, 3.0],
            2,
            True,
            3,
        ],
        python=erosion_mex,
        inputs={"a": V(lambda: to_matlab(CounterRNG(7).uniform((9, 21, 18)) > 0.15))},
        out_roles=("vol",),
    )
)


def erosion_error(which: str, threads=None):
    a = CounterRNG(8).uniform((5, 10, 12)) > 0.5
    t = [np.array([[True]])]
    z = (0, 0)
    if which == "double":
        return D.mask_edge_erosion(a.astype(np.float64), t, z, z, 1, True)
    if which == "4d":
        return D.mask_edge_erosion(np.ones((2, 3, 3, 3), bool), t, z, z, 1, True)
    if which == "not_cell":
        return D.mask_edge_erosion(a, np.array([[True]]), z, z, 1, True)
    if which == "false_nhood":
        return D.mask_edge_erosion(a, [np.zeros((3, 3), bool)], z, z, 1, True)
    if which == "3d_nhood":
        return D.mask_edge_erosion(a, [np.ones((2, 2, 2), bool)], z, z, 1, True)
    if which == "neg_pad":
        return D.mask_edge_erosion(a, t, (0, -1), z, 1, True)
    if which == "scalar_pad":
        return D.mask_edge_erosion(a, t, z, 0, 1, True)
    if which == "neg_rz":
        return D.mask_edge_erosion(a, t, z, z, -1, True)
    if which == "frac_rz":
        return D.mask_edge_erosion(a, t, z, z, 1.5, True)
    if which == "nan_nhood":
        return D.mask_edge_erosion(a, [np.array([[1.0, NAN]])], z, z, 0, True)
    raise AssertionError(which)


_EROSION_ID = "decon_mask_edge_erosion_mex:inputError"
for _w, _code in [
    ("double", "decon_mask_edge_erosion_mex(double(A), {true}, [0, 0], [0, 0], 1, true, 1)"),
    ("4d", "decon_mask_edge_erosion_mex(true(3, 3, 3, 2), {true}, [0, 0], [0, 0], 1, true, 1)"),
    ("not_cell", "decon_mask_edge_erosion_mex(A, true, [0, 0], [0, 0], 1, true, 1)"),
    ("false_nhood", "decon_mask_edge_erosion_mex(A, {false(3)}, [0, 0], [0, 0], 1, true, 1)"),
    ("3d_nhood", "decon_mask_edge_erosion_mex(A, {true(2, 2, 2)}, [0, 0], [0, 0], 1, true, 1)"),
    ("neg_pad", "decon_mask_edge_erosion_mex(A, {true}, [0, -1], [0, 0], 1, true, 1)"),
    ("scalar_pad", "decon_mask_edge_erosion_mex(A, {true}, [0, 0], 0, 1, true, 1)"),
    ("neg_rz", "decon_mask_edge_erosion_mex(A, {true}, [0, 0], [0, 0], -1, true, 1)"),
    ("frac_rz", "decon_mask_edge_erosion_mex(A, {true}, [0, 0], [0, 0], 1.5, true, 1)"),
    ("nan_nhood", "decon_mask_edge_erosion_mex(A, {[1, NaN]}, [0, 0], [0, 0], 0, true, 1)"),
]:
    CASES.append(
        Case(
            f"erosion.error.{_w}",
            code=f"A = in{{1}}; {_code};",
            args=[Vol(In("a"))],
            python=erosion_error,
            py_args=[],
            py_kwargs={"which": _w},
            inputs={"a": V(lambda: to_matlab(CounterRNG(8).uniform((5, 10, 12)) > 0.5))},
            # (the mex's call of logical() raises MATLAB's error, not the mex's)
            expect_error="MATLAB:nologicalnan" if _w == "nan_nhood" else _EROSION_ID,
            py_error=D.DeconKernelError,
        )
    )


# ---------------------------------------------------------------------------------------------
# the OMW back projector's steps: omw_backprojector_mex
# ---------------------------------------------------------------------------------------------


def otf_like_mask(py_shape, seed: int, holes: bool = True) -> np.ndarray:
    """a blob like an OTF support (an ellipsoid with a bumpy surface, a few holes and specks);
    (z, y, x)"""
    rng = CounterRNG(seed, "otfmask")
    nz, ny, nx = py_shape
    z = (np.arange(nz) - (nz - 1) / 2).reshape(-1, 1, 1) / (0.35 * nz)
    y = (np.arange(ny) - (ny - 1) / 2).reshape(1, -1, 1) / (0.42 * ny)
    x = (np.arange(nx) - (nx - 1) / 2).reshape(1, 1, -1) / (0.38 * nx)
    r = z * z + y * y + x * x
    m = r + 0.25 * rng.uniform(py_shape) < 1.0
    if holes:
        m &= rng.uniform(py_shape) > 0.02
        m |= rng.uniform(py_shape) > 0.995
    return m


def convhull_case(a, zi, yi, xi, nt, *, threads=None):
    return D.omw_convhull(a, zi, yi, xi, threads=threads)


def _hull_inputs(shape, seed):
    def make():
        return to_matlab(otf_like_mask(shape, seed))

    return make


for _k, (_shape, _sel) in enumerate(
    [((21, 30, 26), "all"), ((17, 33, 29), "some"), ((2, 9, 7), "all"), ((40, 24, 31), "all")]
):
    _m = otf_like_mask(_shape, 100 + _k)

    def _inds(m=_m, sel=_sel):
        z = np.nonzero(m.any(axis=(1, 2)))[0]
        y = np.nonzero(m.any(axis=(0, 2)))[0]
        x = np.nonzero(m.any(axis=(0, 1)))[0]
        if sel == "some":
            z, y, x = z[::2], y[1::3], x[::4]
        return [(v + 1).astype(np.float64).reshape(1, -1) for v in (z, y, x)]

    _zi, _yi, _xi = _inds()
    CASES.append(
        Case(
            f"omw.convhull.{_k}",
            matlab="omw_backprojector_mex",
            args=["convhull", Vol(In("a")), Index1(_zi), Index1(_yi), Index1(_xi), 4],
            python=lambda mode, *a, **k: convhull_case(*a, **k),
            inputs={"a": V(_hull_inputs(_shape, 100 + _k))},
            out_roles=("vol",),
        )
    )


def window_inputs(shape, seed: int):
    """the edge of an OTF-like mask (voxels of the mask with a neighbor outside it along an axis)
    as unit vectors E ([y x z] rows, exact arithmetic) and radii R, as omw_window_with_mex makes
    them (MATLAB's form)"""
    m = otf_like_mask(shape, seed, holes=False)
    edge = m.copy()
    edge[1:-1, 1:-1, 1:-1] &= ~(
        m[:-2, 1:-1, 1:-1]
        & m[2:, 1:-1, 1:-1]
        & m[1:-1, :-2, 1:-1]
        & m[1:-1, 2:, 1:-1]
        & m[1:-1, 1:-1, :-2]
        & m[1:-1, 1:-1, 2:]
    )
    edge_m = to_matlab(edge)
    ey, ex, ez = np.nonzero(edge_m)  # MATLAB order is find's: sort by linear index
    order = np.lexsort((ey, ex, ez))
    ey, ex, ez = ey[order], ex[order], ez[order]
    ny, nx, nz = edge_m.shape
    v = np.stack([ey + 1 - (ny + 1) / 2, ex + 1 - (nx + 1) / 2, ez + 1 - (nz + 1) / 2], axis=1)
    rad = np.sqrt(v[:, 0] * v[:, 0] + v[:, 1] * v[:, 1] + v[:, 2] * v[:, 2])
    keep = rad > 0
    v, rad = v[keep], rad[keep]
    return v / rad[:, None], rad.reshape(-1, 1)


def window_case(mode, psz, e, r, bounds, nt, *, threads=None):
    mask, inds = D.omw_window(
        psz, np.ascontiguousarray(e[:, [2, 0, 1]]), r, bounds, threads=threads
    )
    return [mask, inds]


def _lin_col(v, msize):
    from petakit5d.core.arrays import zyx_to_matlab_index

    k = zyx_to_matlab_index(tuple(v), (msize[2], msize[0], msize[1]))
    return (np.asarray(k, np.float64) + 1).reshape(-1, 1)


for _k, (_shape, _b) in enumerate(
    [
        ((21, 30, 26), (0.8, 1.0)),
        ((33, 29, 41), (0.8, 1.0)),
        ((15, 15, 15), (0.5, 0.9)),
        ((25, 40, 18), (0.9, 1.2)),
    ]
):
    _msz = [float(_shape[1]), float(_shape[2]), float(_shape[0])]
    CASES.append(
        Case(
            f"omw.window.{_k}",
            matlab="omw_backprojector_mex",
            args=["window", PerAxis(_msz, "yxz"), In("e"), In("r"), list(_b), 4],
            nargout=2,
            python=window_case,
            inputs={
                "e": V(lambda s=_shape, k=_k: window_inputs(s, 200 + k)[0]),
                "r": V(lambda s=_shape, k=_k: window_inputs(s, 200 + k)[1]),
            },
            out_roles=("vol", lambda v, m=_msz: _lin_col(v, m)),
        )
    )
CASES.append(
    Case(
        "omw.window.axes",
        matlab="omw_backprojector_mex",
        args=[
            "window",
            PerAxis([8.0, 8.0, 8.0], "yxz"),
            np.eye(3),
            np.array([[1.0], [2.0], [3.0]]),
            [0.8, 1.0],
            1,
        ],
        nargout=2,
        python=window_case,
        out_roles=("vol", lambda v: _lin_col(v, [8, 8, 8])),
    )
)


def components_mask(seed: int, ties: bool):
    rng = CounterRNG(seed, "components")
    a = rng.uniform((23, 31, 27)) > 0.92
    from petakit5d.matlabcompat.morphology import imdilate

    a = imdilate(a, np.ones((1, 2, 2), bool))
    if ties:
        a[:] = False
        a[2:5, 2:5, 2:5] = True  # 27 voxels, first
        a[10:13, 10:13, 10:13] = True  # 27 voxels, second
        a[18:20, 20:25, 3:6] = True  # 30 voxels
        a[0, 30, 26] = True
    return to_matlab(a)


def components_case(mode, a, cmode, n, nt, *, threads=None):
    return D.omw_components(a, cmode, n, threads=threads)


for _k, _ties in enumerate([False, True]):
    for _mode, _n in [
        ("minsize", 0),
        ("minsize", 5),
        ("minsize", 27),
        ("largest", 1),
        ("largest", 2),
        ("largest", 3),
        ("largest", 100),
    ]:
        CASES.append(
            Case(
                f"omw.components.m{_k}.{_mode}{_n}",
                matlab="omw_backprojector_mex",
                args=["components", Vol(In("a")), _mode, float(_n), 4],
                python=components_case,
                inputs={"a": V(lambda s=_k, t=_ties: components_mask(300 + s, t))},
                out_roles=("vol",),
            )
        )


def sort_inputs(which: str):
    rng = CounterRNG(400, which)
    a = (rng.uniform((20, 31, 17), np.float32) * np.float32(3)).astype(np.float32)
    a[a < np.float32(0.5)] = 0
    if which == "nan":
        a[3, 4, 5] = NAN
    elif which == "negzero":
        a[0, 0, 1] = -0.0
    elif which == "negative":
        a[5, 5, 5] = -1
    elif which == "inf":
        a[1, 2, 3] = INF
    return to_matlab(a)


for _w in ["plain", "inf", "nan", "negzero", "negative"]:
    CASES.append(
        Case(
            f"omw.sortdesc.{_w}",
            matlab="omw_backprojector_mex",
            args=["sortdesc", Vol(In("a")), 4],
            python=lambda mode, a, nt, threads=None: D.omw_sortdesc(a, threads=threads),
            py_out=lambda v: [v.reshape(-1, 1) if v is not None else None],
            inputs={"a": V(lambda w=_w: sort_inputs(w))},
        )
    )
    CASES.append(
        Case(
            f"omw.medfilt3.{_w}",
            matlab="omw_backprojector_mex",
            args=["medfilt3", Vol(In("a")), 4],
            python=lambda mode, a, nt, threads=None: D.omw_medfilt3(a, threads=threads),
            py_out=lambda v: [to_matlab(v) if v is not None else None],
            inputs={"a": V(lambda w=_w: sort_inputs(w))},
        )
    )


def omw_error(which: str, threads=None):
    m = np.ones((4, 4, 4), bool)
    e = np.eye(3)
    r = np.array([1.0, 2.0, 3.0])
    if which == "hull_2d":
        return D.omw_convhull(np.ones((4, 4), bool), [], [], [])
    if which == "hull_repeat":
        return D.omw_convhull(m, [0, 0], [], [])
    if which == "hull_range":
        return D.omw_convhull(m, [4], [], [])
    if which == "window_bounds":
        return D.omw_window((8, 8, 8), e, r, (1.0, 0.8))
    if which == "window_unit":
        return D.omw_window((8, 8, 8), e * 2, r, (0.8, 1.0))
    if which == "window_one_point":
        return D.omw_window((8, 8, 8), e[:1], r[:1], (0.8, 1.0))
    if which == "components_mode":
        return D.omw_components(m, "smallest", 1)
    if which == "components_n":
        return D.omw_components(m, "minsize", 1.5)
    raise AssertionError(which)


_OMW_ID = "omw_backprojector_mex:invalidInput"
for _w, _code in [
    ("hull_2d", "omw_backprojector_mex('convhull', true(4, 4), [], [], [], 1)"),
    ("hull_repeat", "omw_backprojector_mex('convhull', true(4, 4, 4), [1, 1], [], [], 1)"),
    ("hull_range", "omw_backprojector_mex('convhull', true(4, 4, 4), 5, [], [], 1)"),
    (
        "window_bounds",
        "[a, b] = omw_backprojector_mex('window', [8, 8, 8], eye(3), [1; 2; 3], [1, 0.8], 1)",
    ),
    (
        "window_unit",
        "[a, b] = omw_backprojector_mex('window', [8, 8, 8], 2 * eye(3), [1; 2; 3], [0.8, 1], 1)",
    ),
    (
        "window_one_point",
        "[a, b] = omw_backprojector_mex('window', [8, 8, 8], [1, 0, 0], 1, [0.8, 1], 1)",
    ),
    ("components_mode", "omw_backprojector_mex('components', true(4, 4, 4), 'smallest', 1, 1)"),
    ("components_n", "omw_backprojector_mex('components', true(4, 4, 4), 'minsize', 1.5, 1)"),
]:
    CASES.append(
        Case(
            f"omw.error.{_w}",
            code=f"{_code};",
            python=omw_error,
            py_kwargs={"which": _w},
            expect_error=_OMW_ID,
            py_error=D.DeconKernelError,
        )
    )
# MexInputCheckTest/omwBackProjectorWindowNeedsBothOutputs: MATLAB's nargout and sparse rules have
# no Python counterpart (both outputs are always returned); references only
for _w, _code in [
    (
        "window_one_output",
        "a = omw_backprojector_mex('window', [8, 8, 8], eye(3), [1; 2; 3], [0.8, 1], 1)",
    ),
    (
        "window_sparse",
        "[a, b] = omw_backprojector_mex('window', [8, 8, 8], sparse(eye(3)), [1; 2; 3], "
        "[0.8, 1], 1)",
    ),
    ("hull_sparse", "omw_backprojector_mex('convhull', true(4, 4, 4), sparse(1), [], [], 1)"),
]:
    CASES.append(Case(f"omw.error.{_w}", code=f"{_code};", expect_error=_OMW_ID))


# ---------------------------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------------------------


def _check(case, tmp_path):
    if case.id.endswith(".exact"):
        reason = exact_skip_reason()
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
    and c.python is not None
    and not c.id.startswith("erosion.mex.random")
    and not c.id.startswith("accel.error")
    and "fft1" not in c.id
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


# --- the steps against compositions of the port's own FFTs (no MATLAB; these run in CI) ----------

EPS32 = np.float32(2.0**-52)
needs_exact = pytest.mark.skipif(not EXACT, reason="the exact FFT backend is not installed")


def mprod(s, h, conj=False):
    """MATLAB's complex product (no fused multiply-add), h .* s or conj(h) .* s"""
    sr, si, hr, hi = s.real, s.imag, h.real, h.imag
    out = np.empty(s.shape, dtype=s.dtype)
    if conj:
        out.real = sr * hr + si * hi
        out.imag = si * hr - sr * hi
    else:
        out.real = sr * hr - si * hi
        out.imag = sr * hi + si * hr
    return out


def positive(w):
    return np.where(w > 0, w, np.float32(0)).astype(np.float32)


def padded_psf(psf, shape):
    p = np.zeros(shape, np.float32)
    p[tuple(slice(0, n) for n in psf.shape)] = psf
    return np.roll(p, tuple(-(n // 2) for n in psf.shape), axis=tuple(range(p.ndim)))


class Composition:
    """decon_lucy_omw_function's and decon_lucy_function's MATLAB loops, on matlabcompat's fftn and
    ifftn (MATLAB's driver) or on HalfSpectrum's native transforms"""

    def __init__(self, shape, backend, native):
        from petakit5d.matlabcompat import fft as F

        self.F, self.be, self.native = F, backend, native
        self.hs = P.HalfSpectrum(shape, backend=backend, threads=3) if native else None
        self.shape = shape

    def otf(self, psf):
        if self.native:
            return self.hs.otf(psf)
        return self.F.fftn(padded_psf(psf, self.shape), backend=self.be, threads=3)

    def conv(self, x, h, conj=False):
        if self.native:
            return self.hs.convolve(x, h, conj)
        s = self.F.fftn(x, backend=self.be, threads=3)
        return self.F.ifftn(mprod(s, h, conj), backend=self.be, threads=3).real.astype(np.float32)

    def omw(self, raw, psf, psfb, n, bg):
        hf, hb = self.otf(psf), self.otf(psfb)
        i = positive((raw.astype(np.float64) - bg).astype(np.float32))
        j = i.copy()
        for _ in range(n):
            v = self.conv(j, hf)
            v = self.conv(i / np.where(v > EPS32, v, EPS32), hb)
            j = positive(v * j)
        return j

    def rlstep(self, y, i, psf):
        h = self.otf(psf)
        v = self.conv(y, h)
        v = self.conv(positive(i) / np.where(v > EPS32, v, EPS32), h, True)
        return positive(v * y)


STRUCT_SIZES = [(33, 66, 64), (20, 47, 45), (64, 23, 10), (64, 50), (47, 30), (33, 2, 64)]
LAYOUTS = [("permissive", "matlab"), ("permissive", "native")] + (
    [("exact", "matlab")] if EXACT else []
)


def py_inputs(shape, seed=1):
    ms = list(size_to_matlab(shape))
    raw = from_matlab(blurred_beads(ms, seed, "uint16"))
    psf = from_matlab(make_psf(ms, "f"))
    psfb = from_matlab(make_psf(ms, "b"))
    return raw, psf, psfb


@pytest.mark.parametrize("be,layout", LAYOUTS)
@pytest.mark.parametrize("shape", STRUCT_SIZES)
def test_steps_are_the_matlab_loops_on_the_backends_ffts(be, layout, shape):
    """DeconRFFTTest/mexEqualsMatlabLoop in the port: in MATLAB's layout the mex's loops equal the
    MATLAB loops on matlabcompat's fftn and ifftn (MATLAB's driver) on the same backend; in the
    native layout they equal the loops on HalfSpectrum's transforms; bit for bit"""
    raw, psf, psfb = py_inputs(shape)
    k = D.LucyRFFT(shape, backend=be, layout=layout, threads=4)
    ref = Composition(shape, be, layout == "native")
    hf, hb = k.otf(psf), k.otf(psfb)
    out = k.omw(raw, hf, hb, 2, 100.3)
    assert out.dtype == np.float32 and out.shape == shape
    assert out.tobytes() == ref.omw(raw, psf, psfb, 2, 100.3).tobytes()
    y = (raw.astype(np.float32) - np.float32(99.5)).astype(np.float32)
    i = (raw.astype(np.float32) - np.float32(100)).astype(np.float32)
    assert k.rlstep(y, i, hf).tobytes() == ref.rlstep(y, i, psf).tobytes()


@pytest.mark.parametrize("be,layout", LAYOUTS)
def test_half_otf_is_the_first_half_of_the_otf(be, layout):
    """DeconRFFTTest/halfOTFIsTheFirstHalfOfTheOTF: MATLAB's half (along z; 2D along x) of fftn of
    the padded, shifted PSF (equal values; the mex's fill writes -im where fftn's writes 0 - im, so
    zeros may differ in sign)"""
    from petakit5d.matlabcompat import fft as F

    for shape in [(36, 48, 46), (33, 47, 45), (64, 50), (23, 10, 64)]:
        _, psf, _ = py_inputs(shape)
        h = D.LucyRFFT(shape, backend=be, layout=layout).otf(psf)
        full = F.fftn(padded_psf(psf, shape), backend=be)
        if layout == "matlab":
            m = to_matlab(full)
            ref = m[:, :, : shape[0] // 2 + 1] if len(shape) == 3 else m[:, : shape[1] // 2 + 1]
            assert h.flags.f_contiguous and h.shape == ref.shape
            np.testing.assert_array_equal(h, ref)
        else:
            ref = full[..., : shape[-1] // 2 + 1]
            assert np.max(np.abs(h - ref)) <= 1e-5 * np.max(np.abs(ref))
    with pytest.raises(D.DeconKernelError) as e:
        D.otf(np.ones((5, 5, 5), np.float32), (10, 4, 12), backend=be, layout=layout)
    assert e.value.identifier == "decon_lucy_rfft_mex:outSizeIsSmallerThanPsfSize"


@pytest.mark.parametrize("be,layout", LAYOUTS)
@pytest.mark.parametrize("shape", [(33, 66, 64), (64, 50)])
def test_results_do_not_depend_on_the_thread_count(be, layout, shape):
    """DeconRFFTTest/resultsDoNotDependOnTheThreadCount, in both of MATLAB's FFT modes"""
    raw, psf, psfb = py_inputs(shape, 2)
    y = (raw.astype(np.float32) - np.float32(99)).astype(np.float32)
    for parallel in (True, False):
        outs = []
        for t in (1, 2, 7, 24):
            k = D.LucyRFFT(shape, backend=be, layout=layout, threads=t, matlab_parallel=parallel)
            hf, hb = k.otf(psf), k.otf(psfb)
            outs.append(
                (hf.tobytes(), k.omw(raw, hf, hb, 2, 100).tobytes(), k.rlstep(y, y, hf).tobytes())
            )
        assert all(o == outs[0] for o in outs)


@pytest.mark.parametrize("be,layout", LAYOUTS)
def test_the_update_makes_negative_zeros_positive(be, layout):
    """the mex's max(v .* J, 0) is +0 where v .* J is -0 (J = +0 below the background, a negative
    back projection), where MATLAB's max keeps -0"""
    shape = (36, 48, 46)
    raw = from_matlab(below_background([48, 46, 36], 90))
    psf, psfb = from_matlab(make_psf([48, 46, 36], "f")), from_matlab(make_psf([48, 46, 36], "neg"))
    k = D.LucyRFFT(shape, backend=be, layout=layout, threads=3)
    out = k.omw(raw, k.otf(psf), k.otf(psfb), 1, 100)
    ref = Composition(shape, be, layout == "native")
    hb = ref.otf(psfb)
    i = positive((raw.astype(np.float64) - 100).astype(np.float32))
    v = ref.conv(i, ref.otf(psf))
    v = ref.conv(i / np.where(v > EPS32, v, EPS32), hb)
    w = v * i
    assert np.any((w == 0) & np.signbit(w))  # the products of -0 are there
    assert not np.any((out == 0) & np.signbit(out))
    assert out.tobytes() == ref.omw(raw, psf, psfb, 1, 100).tobytes()


@pytest.mark.parametrize("be,layout", LAYOUTS)
def test_omw_continues_from_an_estimate(be, layout):
    """DeconRFFTTest/omwContinuesFromAnEstimate: J0 continues the iterations as one call runs
    them"""
    shape = (36, 48, 46)
    raw, psf, psfb = py_inputs(shape, 3)
    k = D.LucyRFFT(shape, backend=be, layout=layout, threads=3)
    h, hb = k.otf(psf), k.otf(psfb)
    j5 = k.omw(raw, h, hb, 5, 100)
    j2 = k.omw(raw, h, hb, 2, 100)
    assert k.omw(raw, h, hb, 3, 100, j2).tobytes() == j5.tobytes()
    j1 = k.omw(raw, h, hb, 1, 100)
    j11 = k.omw(raw, h, hb, 1, 100, j1)
    assert k.omw(raw, h, hb, 3, 100, j11).tobytes() == j5.tobytes()
    assert k.omw(raw, h, hb, 0, 100, j2).tobytes() == j2.tobytes()


@pytest.mark.parametrize("be,layout", LAYOUTS)
def test_zero_iterations_give_the_background_subtracted_input(be, layout):
    """DeconRFFTTest/zeroIterationsGiveTheBackgroundSubtractedInput: single - double in double,
    rounded to single, then max(., 0); a count below 1 or NaN runs no iteration, 2.7 runs 2"""
    shape = (36, 48, 46)
    raw, psf, psfb = py_inputs(shape, 4)
    im = (raw.astype(np.float32) + np.float32(0.37)).astype(np.float32)
    k = D.LucyRFFT(shape, backend=be, layout=layout)
    h, hb = k.otf(psf), k.otf(psfb)
    ref = positive((im.astype(np.float64) - 100.3).astype(np.float32))
    for n in (0, -1, 0.5, float("nan")):
        assert k.omw(im, h, hb, n, 100.3).tobytes() == ref.tobytes()
    u = positive((raw.astype(np.float64) - 100.3).astype(np.float32))
    assert k.omw(raw, h, hb, 0, 100.3).tobytes() == u.tobytes()
    assert k.omw(raw, h, hb, 2.7, 100.3).tobytes() == k.omw(raw, h, hb, 2, 100.3).tobytes()


@pytest.mark.parametrize("be,layout", LAYOUTS)
def test_internal_arrays_keep_the_layout_across_steps(be, layout):
    """RL's loop on internal arrays (no transposes per step) gives the steps' results"""
    shape = (33, 66, 64)
    raw, psf, _ = py_inputs(shape, 5)
    k = D.LucyRFFT(shape, backend=be, layout=layout, threads=3)
    h = k.otf(psf)
    i = (raw.astype(np.float32) - np.float32(100)).astype(np.float32)
    j = i
    ji = k.to_internal(i)
    ii = ji
    for _ in range(3):
        j = k.rlstep(j, i, h)
        ji = k.rlstep(ji, ii, h, internal=True)
    assert ji.shape == k.internal_shape
    assert k.from_internal(ji).tobytes() == j.tobytes()
    if layout == "matlab":
        assert ji.flags.f_contiguous and np.array_equal(from_matlab(ji), j)


def test_layouts_round_trip():
    rng = CounterRNG(9)
    for shape in [(5, 7, 9), (1, 3, 4), (6, 1, 5), (7, 9)]:
        for dt in (np.float32, np.uint16, np.complex64, np.uint8):
            a = (rng.uniform(shape) * 100).astype(dt)
            m = D.to_matlab_layout(a, threads=3)
            assert m.flags.f_contiguous and np.array_equal(m, to_matlab(a))
            assert np.array_equal(D.from_matlab_layout(m), a)


@needs_exact
def test_lines_are_grouped_as_matlab_does():
    """DeconRFFTTest/linesAreGroupedAsMatlabDoes: 'fft1' runs MATLAB's groups of lines: its own fft
    of a group of four or of one line gives the reference"""
    from petakit5d.matlabcompat import fft as F

    for n_cols in (2000, 222, 210):
        x = CounterRNG(n_cols, "fft1").uniform((64, n_cols), np.float32)

        def by_pieces(pieces, x=x):
            out = []
            for a, b in pieces:
                j = a
                while j < b:
                    m = 4 if b - j >= 4 else 1
                    out.append(
                        F.fft(np.ascontiguousarray(x[:, j : j + m]), axis=-2, backend="exact")[:33]
                    )
                    j += m
            return np.concatenate(out, axis=1)

        def halve(a, b):
            if b - a <= 32:
                return [(a, b)]
            mid = a + (b - a) // 2
            return halve(a, mid) + halve(mid, b)

        split = by_pieces(halve(0, n_cols) if 64 * n_cols * 7 > 98304 else [(0, n_cols)])
        serial = by_pieces([(0, n_cols)])
        assert D.fft1(x, backend="exact").tobytes() == split.tobytes()
        assert D.fft1(x, backend="exact", matlab_parallel=False).tobytes() == serial.tobytes()


def test_the_native_layout_needs_the_permissive_backend():
    if EXACT:
        with pytest.raises(ValueError):
            D.LucyRFFT((4, 4, 4), backend="exact", layout="native")
    assert D.LucyRFFT((4, 4, 4), backend="permissive").layout == "native"
    with pytest.raises(D.DeconKernelError):
        D.LucyRFFT((4, 1, 4))
    with pytest.raises(D.DeconKernelError):
        D.LucyRFFT((4,))


# --- acceleration and sums ------------------------------------------------------------------------


def accel_arrays(shape, dt, seed=11):
    g = CounterRNG(seed)
    j2 = (g.uniform(shape) * 1000).astype(dt)
    j4 = (g.normal(shape) * 10).astype(dt)
    y = (j2 - dt(0.6) * j4 + g.normal(shape).astype(dt)).astype(dt)
    return j2, y, j4


@pytest.mark.parametrize("dt", [np.float32, np.float64])
@pytest.mark.parametrize("shape", [(49, 47, 45), (3, 200, 160), (88999,), (40, 41, 47)])
def test_lambda_sums_in_matlabs_order_in_both_layouts(dt, shape):
    """the products in the class, summed in MATLAB's linear order (blocks of 4096 from 89,000
    elements, MATLAB's sum below): (z, y, x) arrays and MATLAB's arrays give the same bits, at any
    thread count"""
    from petakit5d.kernels import ops

    j2, y, j4 = accel_arrays(shape, dt)
    ref = None
    for t in (1, 2, 7, 24):
        r = D.accel_lambda(j2, y, j4, threads=t)
        if len(shape) > 1:
            rm = D.accel_lambda(
                to_matlab(j2), to_matlab(y), to_matlab(j4), matlab_layout=True, threads=t
            )
            assert rm.sums == r.sums and rm.lam == r.lam
            assert rm.j4.flags.f_contiguous and np.array_equal(from_matlab(rm.j4), r.j4)
        ref = ref or r
        assert r.sums == ref.sums and r.j4.tobytes() == ref.j4.tobytes()
    with np.errstate(all="ignore"):
        d = (j2 - y).astype(dt)
    assert ref.j4.tobytes() == d.tobytes()
    if j2.size < D.SUM_BLOCKS_FROM:
        assert ref.sums[0] == ops.sum_double(d * j4)
    assert 0 < ref.lam < 1


def test_sum_double_switches_and_blocks():
    """DeconAccelTest/sumIsMatlabsAtOneThread in the port: below 89,000 elements MATLAB's sum; from
    there blocks of 4096 in four lanes (a reference written out), at any thread count and layout"""
    from petakit5d.kernels import ops

    for n in (1, 43, 44, 4097, 88999, 89000, 89001, 4096 * 30 + 3, 250047):
        for cls in ("single", "double"):
            x = wide_range(n, cls, n % 89).reshape(-1)
            got = {D.sum_double(x, threads=t) for t in (1, 3, 24)}
            assert len(got) == 1
            s = got.pop()
            if n < D.SUM_BLOCKS_FROM:
                assert s == ops.sum_double(x)
                continue
            v = x.astype(np.float64)
            bs = []
            for b0 in range(0, n, 4096):
                lanes = [-0.0] * 4
                for t, val in enumerate(v[b0 : b0 + 4096].tolist()):
                    lanes[t & 3] = lanes[t & 3] + val
                bs.append(((lanes[0] + lanes[1]) + lanes[2]) + lanes[3])
            ref = -0.0
            for b in bs:
                ref = ref + b
            assert s == ref
            if n % 7 == 0:
                vol = np.reshape(x, (n // 7, 7))
                flat = np.ascontiguousarray(matlab_order(vol)).reshape(-1)
                assert D.sum_double(vol) == D.sum_double(flat)
    assert math.copysign(1, D.sum_double(-np.zeros(100000, np.float32))) == -1
    assert D.sum_double(np.zeros(0, np.float32)) == 0


def _bits(v: float) -> int:
    return int(np.array([v]).view(np.uint64)[0])


def test_a_nan_sum_has_the_payload_of_the_mexs_order():
    """a sum that is NaN is summed again with the running sum as the first operand of every
    addition (x86 returns the first NaN operand): the first NaN in MATLAB's order wins, across
    blocks and, within a block, the lanes added as ((l1 + l2) + l3) + l4"""
    n = 100000
    a = np.array([0x7FF8000000000001], np.uint64).view(np.float64)[0]
    b = np.array([0xFFF8123456789ABC], np.uint64).view(np.float64)[0]
    x = np.ones(n, np.float64)
    x[10], x[50000] = a, b
    assert _bits(D.sum_double(x)) == _bits(a)
    x[10], x[50000] = b, a
    assert _bits(D.sum_double(x, threads=5)) == _bits(b)
    y = np.ones(n, np.float64)
    y[6], y[5] = a, b  # lanes 2 and 1 of the first block: lane 1's NaN comes first
    assert _bits(D.sum_double(y)) == _bits(b)
    w = np.ones(n, np.float64)
    w[8], w[12] = a, b  # one lane: the running sum (a) is the first operand when b comes
    assert _bits(D.sum_double(w)) == _bits(a)


@pytest.mark.parametrize("cls", ["single", "double"])
def test_predict_is_the_matlab_expression(cls):
    """DeconAccelTest/predictIsTheMatlabExpression: lambda * x in double rounded to the class, max
    is 0 for NaN and keeps -0; any layout"""
    dt = np.float32 if cls == "single" else np.float64
    j2, _, _, j3 = (a.reshape(-1) for a in accel_inputs(103635, cls, True, 30))
    for lam in (0.0, -0.0, 1.0, 0.5, 0.123456789012345678, 1 - 2.0**-52, 2.0**-60):
        out = D.accel_predict(j2, j3, lam, threads=3)
        with np.errstate(all="ignore"):
            u = (j2 + (lam * (j2 - j3).astype(np.float64)).astype(dt)).astype(dt)
        ref = np.where(u >= 0, u, dt(0)).astype(dt)
        assert out.tobytes() == ref.tobytes()
        v, w = np.reshape(j2, (49, 47, 45)), np.reshape(j3, (49, 47, 45))
        got = from_matlab(D.accel_predict(to_matlab(v), to_matlab(w), lam))
        assert got.tobytes() == D.accel_predict(v, w, lam).tobytes()


@pytest.mark.parametrize(
    "case,expected",
    [
        ("aboveOne", 1.0),
        ("negative", 0.0),
        ("nan", 1.0),
        ("negativeZero", -0.0),
        ("zeroDenominator", 0.0),
        ("inf", 1.0),
    ],
)
def test_lambda_is_clamped_as_matlab_does(case, expected):
    """DeconAccelTest/lambdaIsClampedAsMatlabDoes: max(min(lambda, 1), 0), min(NaN, 1) is 1, -0
    kept"""
    j2, y, j4 = clamp_inputs(case)
    assert _bits(D.accel_lambda(j2, y, j4).lam) == _bits(expected)


@pytest.mark.parametrize("shape", [(49, 47, 45), (40, 41, 47)])
@pytest.mark.parametrize("ml", [False, True])
def test_accel_writes_into_out_in_place(shape, ml):
    """out= writes the new J_4 into J_4 itself, and the prediction into J_3, with the results of
    new arrays (RL's loop reuses its buffers)"""
    j2, y, j4 = accel_arrays(shape, np.float32)
    j3 = (j2 * np.float32(0.9)).astype(np.float32)
    if ml:
        j2, y, j4, j3 = (to_matlab(a) for a in (j2, y, j4, j3))
    ref = D.accel_lambda(j2, y, j4, matlab_layout=ml)
    pref = D.accel_predict(j2, j3, ref.lam)
    r = D.accel_lambda(j2, y, j4, matlab_layout=ml, out=j4)
    assert r.j4 is j4 and r.sums == ref.sums and j4.tobytes(order="A") == ref.j4.tobytes(order="A")
    p = D.accel_predict(j2, j3, ref.lam, out=j3)
    assert p is j3 and j3.tobytes(order="A") == pref.tobytes(order="A")
    # a NaN sum is summed again from the same J_4, also in place
    j2, y, j4 = accel_arrays(shape, np.float32, 12)
    a, b = (np.array([v], np.uint32).view(np.float32)[0] for v in (0x7FC00001, 0xFFC12345))
    j4[1, 2, 3], j4[30, 3, 3] = a, b
    if ml:
        j2, y, j4 = (to_matlab(v) for v in (j2, y, j4))
    ref = D.accel_lambda(j2, y, j4, matlab_layout=ml)
    r = D.accel_lambda(j2, y, j4, matlab_layout=ml, out=j4)
    assert np.isnan(ref.sums[0]) and _bits(r.sums[0]) == _bits(ref.sums[0])
    assert _bits(r.sums[1]) == _bits(ref.sums[1]) and j4.tobytes(order="A") == ref.j4.tobytes(
        order="A"
    )
    with pytest.raises(ValueError):
        D.accel_predict(j2, j3, 0.5, out=np.zeros(3, np.float32))


def test_accel_rejects_what_the_mex_rejects():
    for w in (
        "class",
        "size",
        "complex",
        "uint16",
        "lambda_single",
        "lambda_vector",
        "predict_class",
    ):
        with pytest.raises(D.DeconKernelError) as e:
            accel_error(w)
        assert e.value.identifier == "decon_lucy_accel_mex:inputError", w


# --- mask erosion --------------------------------------------------------------------------------


def imerode_reference(mask, e):
    """the previous code of decon_mask_edge_erosion (imerode), on N03's imerode"""
    from petakit5d.matlabcompat.morphology import imerode, strel

    b = mask.copy()
    if b.ndim == 3:
        b[:, [0, -1], :] = False
        b[:, :, [0, -1]] = False
        b[[0, -1]] = False
        b = imerode(b, strel("disk", e - 1))
        return imerode(b, np.ones((2 * e - 1, 1, 1), bool))
    b[[0, -1], :] = False
    b[:, [0, -1]] = False
    return imerode(b, strel("disk", e - 1))


@pytest.mark.parametrize("shape", ["random", "skewed", "thin", "blobs"])
@pytest.mark.parametrize("sname", ["odd3D", "twoSlices", "oneRow", "oneColumn", "image2D"])
def test_masks_equal_imerode(shape, sname):
    """DeconMaskErosionTest/masksEqualImerode against N03's imerode"""
    a = from_matlab(make_mask(shape, MASK_SIZES[sname]))
    for e in (1, 2, 3, 4, 5, 9, 13, 21):
        assert np.array_equal(edge_erosion(a, e), imerode_reference(a, e)), e


def test_erosion_does_not_depend_on_the_thread_count():
    """DeconMaskErosionTest/resultDoesNotDependOnTheThreadCount: whole slices per thread, or the
    rows of each slice split among the threads when there are few slices"""
    for msize in ([60, 70, 40], [60, 70, 3], [80, 90]):
        a = from_matlab(make_mask("skewed", msize))
        a = a & (CounterRNG(1).uniform(a.shape) > 0.001)
        ref = edge_erosion(a, 9, 1)
        for t in (2, 5, 7, 24):
            assert np.array_equal(edge_erosion(a, 9, t), ref)
        assert np.array_equal(ref, imerode_reference(a, 9))


def test_numeric_neighborhoods_are_converted_as_logical_does():
    from petakit5d.matlabcompat.morphology import imerode

    a = CounterRNG(6).uniform((5, 10, 12)) > 0.5
    out = D.mask_edge_erosion(a, [np.array([0.0, 2.0, 1.0])], (0, 0), (0, 0), 0, False)
    assert np.array_equal(out, imerode(a, np.array([[False, True, True]])))
    with pytest.raises(D.DeconKernelError):
        D.mask_edge_erosion(a, [np.array([1.0, NAN])], (0, 0), (0, 0), 0, False)


def test_erosion_rejects_what_the_mex_rejects():
    for w in (
        "double",
        "4d",
        "not_cell",
        "false_nhood",
        "3d_nhood",
        "neg_pad",
        "scalar_pad",
        "neg_rz",
        "frac_rz",
        "nan_nhood",
    ):
        with pytest.raises(D.DeconKernelError) as e:
            erosion_error(w)
        expected = (
            "MATLAB:nologicalnan" if w == "nan_nhood" else "decon_mask_edge_erosion_mex:inputError"
        )
        assert e.value.identifier == expected, w


# --- back projector steps ------------------------------------------------------------------------


def test_components_match_bwareaopen_and_the_largest_ties_go_to_the_first():
    from petakit5d.matlabcompat.label import bwareaopen, bwconncomp

    for seed, ties in ((300, False), (301, True)):
        a = from_matlab(components_mask(seed, ties))
        for n in (0, 5, 27, 31):
            assert np.array_equal(D.omw_components(a, "minsize", n, threads=3), bwareaopen(a, n))
        cc = bwconncomp(a, 26)
        sizes = np.array([p.size for p in cc.pixel_idx_list])
        for k in (1, 2, 3):
            ref = np.zeros(a.size, bool)
            for i in np.argsort(-sizes, kind="stable")[:k]:
                ref[cc.pixel_idx_list[i]] = True
            assert np.array_equal(D.omw_components(a, "largest", k), ref.reshape(a.shape))
        for t in (1, 2, 7):
            assert np.array_equal(
                D.omw_components(a, "largest", 2, threads=t), D.omw_components(a, "largest", 2)
            )


def test_sort_and_median_return_none_where_the_mex_returns_empty():
    from petakit5d.matlabcompat.morphology import medfilt3

    a = from_matlab(sort_inputs("plain"))
    s = D.omw_sortdesc(a, threads=4)
    assert s.tobytes() == np.sort(a.astype(np.float64).reshape(-1))[::-1].tobytes()
    for w in ("nan", "negzero", "negative"):
        assert D.omw_sortdesc(from_matlab(sort_inputs(w))) is None
    assert np.array_equal(D.omw_medfilt3(a), medfilt3(a))
    for w in ("nan", "negzero"):
        assert D.omw_medfilt3(from_matlab(sort_inputs(w))) is None


def test_hulls_are_bwconvhull_of_the_slices():
    from petakit5d.matlabcompat.morphology import bwconvhull

    a = otf_like_mask((17, 33, 29), 101)
    zi, yi, xi = [0, 5, 8], [3, 10], [2, 28]
    out = D.omw_convhull(a, zi, yi, xi, threads=3)
    ref = a.copy()
    for z in zi:
        ref[z] = bwconvhull(ref[z])
    for y in yi:
        ref[:, y, :] = bwconvhull(np.ascontiguousarray(ref[:, y, :].T)).T
    for x in xi:
        ref[:, :, x] = bwconvhull(np.ascontiguousarray(ref[:, :, x].T)).T
    assert np.array_equal(out, ref)


def test_window_is_decided_only_where_the_bounds_hold():
    """DeconBackProjectorTest/windowIsDecidedOnlyWhereTheBoundsHold in the kernel: where the mask is
    1, r / rn < l, and voxels left out of the indices have r / rn >= u, for rn between the radii
    of the two nearest edge points (brute force)"""
    shape = (15, 17, 13)
    em, r = window_inputs(shape, 202)
    e = np.ascontiguousarray(em[:, [2, 0, 1]])
    lo, hi = 0.8, 1.0
    mask, inds = D.omw_window(shape, e, r, (lo, hi), threads=3)
    undecided = np.zeros(shape, bool)
    undecided[inds] = True
    rr = r.reshape(-1)
    for z, y, x in np.ndindex(*shape):
        p = np.array([z - (shape[0] - 1) / 2, y - (shape[1] - 1) / 2, x - (shape[2] - 1) / 2])
        rad = math.sqrt(float(p @ p))
        if rad == 0:
            assert undecided[z, y, x]
            continue
        d = np.sqrt(((e - p / rad) ** 2).sum(axis=1))
        k = np.argsort(d, kind="stable")[:2]
        if mask[z, y, x] == 1:
            assert rad / rr[k].min() < lo and not undecided[z, y, x]
        elif not undecided[z, y, x]:
            assert rad / rr[k].max() >= hi
    assert set(np.unique(mask)) <= {0.0, 1.0}
    k = (inds[0] * shape[2] + inds[2]) * shape[1] + inds[1]  # MATLAB's linear order
    assert np.all(np.diff(k) > 0)


def test_omw_rejects_what_the_mex_rejects():
    for w in (
        "hull_2d",
        "hull_repeat",
        "hull_range",
        "window_bounds",
        "window_unit",
        "window_one_point",
        "components_mode",
        "components_n",
    ):
        with pytest.raises(D.DeconKernelError) as e:
            omw_error(w)
        assert e.value.identifier == "omw_backprojector_mex:invalidInput", w
