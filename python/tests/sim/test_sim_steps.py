"""The steps of the SIM reconstruction (petakit5d.sim, task C06) against MATLAB R2024b's
``third_parties/sim3D`` functions.

Oracle cases (MATLAB has no tests of sim3D; the inputs are ours, from the bit-stable generators,
in MATLAB's form): make_forward_separation_matrix and its pseudo-inverse, cylMask, tukwin,
fourierShift3D, apodizeEllipse, resampleOTF, the black boxes written from published methods
and measurements (perdecomp_3D, normxcorr2_masked, MaskedTranslationRegistration2D_fit,
imtranslate_function), MATLAB's fminsearch and Lsq_GaussianFit_3D. Exact where the arithmetic
is reproduced; otherwise the tolerances of docs/parity/sim.md (proposed rows R1 to R4).

Then unit tests of the methods' properties.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.core.dtypes import mround
from petakit5d.sim import masked_ncc as M
from petakit5d.sim import neldermead as FM
from petakit5d.sim import perdecomp as PD
from petakit5d.sim import recon as R
from petakit5d.sim import translate as T
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value, Vol
from petakit5d.testing.compare import tol
from petakit5d.testing.synth import CounterRNG, stable_exp

# tolerances (docs/parity/sim.md, proposed)
R1_DOUBLE = tol(
    reason="proposed (C06, R1): double-precision steps, rounding of FFTs and sums", rel_max=1e-13
)
R1_SINGLE = tol(
    reason="proposed (C06, R1): single-precision steps (MATLAB's single sin, cos, FFTs)",
    rel_max=2e-6,
    rel_l2=1e-6,
)
R2_FIT = tol(
    reason="proposed (C06, R2): fminsearch on sums MATLAB adds on its threads",
    rel=1e-6,
)
R3_DOUBLE = tol(
    reason="proposed (C06, R3): masked NCC by FFT in double (sums over small overlaps)",
    rel_max=1e-9,
)
R3_SINGLE = tol(
    reason="proposed (C06, R3): masked NCC by FFT in single precision", rel_max=5e-5, abs=1e-5
)


def uniform(seed: int, shape, scale: float = 1.0, offset: float = 0.0) -> np.ndarray:
    """MATLAB's form (Fortran order), double"""
    u = CounterRNG(seed, "c06").uniform(shape)
    return np.asfortranarray(offset + scale * u)


def complex_volume(seed: int, shape) -> np.ndarray:
    a = uniform(seed, shape, 2.0, -1.0) + 1j * uniform(seed + 1, shape, 2.0, -1.0)
    return np.asfortranarray(a)


def mshape(v: np.ndarray) -> tuple[int, ...]:
    """MATLAB's [ny nx nz] of a (z, y, x) shape"""
    return (v.shape[1], v.shape[2], v.shape[0])


# ---------------------------------------------------------------------------------------------
# separation matrix, cylMask, tukwin, fourierShift3D, apodizeEllipse, resampleOTF
# ---------------------------------------------------------------------------------------------

_SEP = (
    "m = make_forward_separation_matrix(in{1}, in{2}, in{3}, in{4}); "
    "out = {m, cast(pinv(m), in{5})};"
)


def py_sep(nph, nor, period, step, prec):
    m = R.separation_matrix(nph, nor, period, step)
    return [m, R.inverse_separation_matrix(m, prec)]


def _dk(msize, pxl):
    return list(1.0 / (np.array(msize, dtype=np.float64) * np.array(pxl, dtype=np.float64)))


_PXL_DATA = [0.11, 0.11, 0.3 * math.sin(32.4 * math.pi / 180)]
_PXL_PSF = [0.11, 0.11, 0.2 * math.sin(32.4 * math.pi / 180)]


def py_cylmask(na_det, wvl_em, na_ext, wvl_ext, nimm, norders, ny, nx, nz, dk, islattice, use_gpu):
    del use_gpu
    m = R.cyl_mask(na_det, wvl_em, na_ext, wvl_ext, nimm, norders, (nz, ny, nx), dk, islattice)
    return m


def py_tukwin(data, r, use_gpu):
    del use_gpu
    return np.ascontiguousarray(np.broadcast_to(R.tukey_window(data.shape, r), data.shape))


def py_fshift(img, shifts, use_gpu, prec):
    del use_gpu
    return R.fourier_shift(img, shifts, prec)


def py_apod(k, dk, pv, angle, na_det, nimm, na_ext, wvl_em, wvl_ext, islattice, use_gpu):
    del use_gpu
    out, apo = R.apodize_ellipse(k, dk, np.asarray(pv), angle, na_det, nimm, na_ext, wvl_em,
                                 wvl_ext, islattice)  # fmt: skip
    return [out, apo]


def py_resample(otf, pxl_psf, data, pxl_data, nph, nor, nori, use_gpu, prec):
    del nor, use_gpu, prec
    nimgs, ny, nx = data.shape
    out = R.resample_otf(
        otf, [pxl_psf[2], pxl_psf[0], pxl_psf[1]], (nimgs // (nph * nori), ny, nx),
        [pxl_data[2], pxl_data[0], pxl_data[1]],
    )  # fmt: skip
    return out[0]


def _step_cases() -> list[Case]:
    cs = []
    for cid, args in [
        ("sep.5", [5, 5, 1.2021, 0.232, "single"]),
        ("sep.3", [3, 3, 0.52932, 0.170, "single"]),
        ("sep.5.double", [5, 5, 1.4, 0.28, "double"]),
        ("sep.7", [7, 7, 1.1, 0.15, "single"]),
    ]:
        cs.append(Case(cid, code=_SEP, args=args, nargout=2, python=py_sep,
                       compare={1: R1_DOUBLE if args[-1] == "double" else R1_SINGLE}))  # fmt: skip
    for cid, nor, lattice, shape in [
        ("cylmask.lattice", 5, True, (30, 26, 11)),
        ("cylmask.3beam", 5, False, (30, 26, 11)),
        ("cylmask.3orders", 3, True, (24, 31, 8)),
    ]:
        dk = _dk(shape, _PXL_DATA)
        cs.append(
            Case(
                cid,
                matlab="cylMask",
                args=[1.0, 0.605, 0.55, 0.56, 1.33, nor, *shape, dk, lattice, False],
                python=py_cylmask,
                out_roles=("vol",),
            )
        )
    for cid, shape, r in [("tukwin.0.1", (30, 26, 4), 0.1), ("tukwin.0.25", (17, 22, 3), 0.25)]:
        cs.append(
            Case(
                cid,
                matlab="tukwin",
                args=[Vol(In("z")), r, False],
                python=py_tukwin,
                inputs={"z": Value(lambda shape=shape: np.zeros(shape, order="F"))},
                out_roles=("vol",),
            )
        )
    for cid, seed, kind, prec, sh in [
        ("fshift.complex.single", 11, "c", "single", [2.3, -1.7, 0.6]),
        ("fshift.complex.double", 12, "c", "double", [-3.6, 4.25, -1.1]),
        ("fshift.real.single", 13, "r", "single", [0.5, 1.5, 0.0]),
    ]:
        shape = (14, 12, 7)
        make = (
            (lambda seed=seed, shape=shape: complex_volume(seed, shape))
            if kind == "c"
            else (lambda seed=seed, shape=shape: uniform(seed, shape, 10.0))
        )
        cs.append(
            Case(
                cid,
                matlab="fourierShift3D",
                args=[Vol(In("v")), sh, False, prec],
                python=py_fshift,
                inputs={"v": Value(make)},
                out_roles=("vol",),
                compare=R1_SINGLE if prec == "single" else R1_DOUBLE,
            )
        )
    pv = [[-7.2, 0.4, 0.0], [-3.6, 0.2, 0.0], [0, 0, 0], [3.6, -0.2, 0.0], [7.2, -0.4, 0.0]]
    shape = (20, 18, 6)
    dk = _dk((10, 9, 6), _PXL_DATA)
    cs.append(
        Case(
            "apod.lattice",
            matlab="apodizeEllipse",
            args=[
                Vol(In("k")),
                dk,
                np.array(pv),
                math.pi / 2,
                1.0,
                1.33,
                0.55,
                0.605,
                0.56,
                True,
                False,
            ],
            nargout=2,
            python=py_apod,
            inputs={"k": Value(lambda: complex_volume(21, shape).astype(np.complex64))},
            out_roles=("vol", "vol"),
            compare={0: R1_SINGLE, 1: R1_DOUBLE},
        )
    )
    cs.append(
        Case(
            "apod.3beam",
            matlab="apodizeEllipse",
            args=[Vol(In("k")), dk, np.array(pv), 0.3, 1.0, 1.33, 0.55, 0.605, 0.56, False, False],
            nargout=2,
            python=py_apod,
            inputs={"k": Value(lambda: complex_volume(22, shape).astype(np.complex64))},
            out_roles=("vol", "vol"),
            compare={0: R1_SINGLE, 1: R1_DOUBLE},
        )
    )
    cs.append(
        Case(
            "resample.linear",
            matlab="resampleOTF",
            args=[Vol(In("otf")), _PXL_PSF, Vol(In("data")), _PXL_DATA, 5, 5, 1, False, "single"],
            python=py_resample,
            inputs={
                "otf": Value(lambda: complex_volume(31, (16, 14, 9, 5))),
                "data": Value(lambda: np.zeros((20, 18, 35), dtype=np.uint16, order="F")),
            },
            out_roles=("vol",),
            compare=R1_DOUBLE,
        )
    )
    return cs


# ---------------------------------------------------------------------------------------------
# the black boxes: perdecomp_3D, the masked NCC and registration, imtranslate_function
# ---------------------------------------------------------------------------------------------


def py_perdecomp(u, use_gpu, prec):
    return PD.perdecomp_3d(u, use_gpu=use_gpu, gpu_precision=prec)


def reg_images(case: str):
    """fixed, moving, masks (MATLAB's (y, x) images): a smooth texture, the moving one shifted"""
    from petakit5d.matlabcompat import elementary as E

    wave = np.outer(E.cos(np.arange(1, 41) / 3), E.sin(np.arange(1, 37) / 4))
    f = uniform(41, (40, 36)) + 0.2 * wave
    fm = np.ones(f.shape, dtype=bool)
    gm = np.ones(f.shape, dtype=bool)
    if case == "int":
        g = np.roll(f, (3, -2), axis=(0, 1))
    else:
        ky = np.fft.ifftshift(np.arange(-20, 20)) / 40
        kx = np.fft.ifftshift(np.arange(-18, 18)) / 36
        # (exact frequencies; the ramp from the stable sine and cosine)
        u = -(ky[:, None] * 1.3 + kx[None, :] * -0.7)
        from petakit5d.testing.synth import stable_sincos_2pi

        s, c = stable_sincos_2pi(u - np.floor(u))
        spec = _dft2(f) * (c + 1j * s)
        g = _idft2(spec).real
        fm[4:12, 7:20] = False
        gm[24:33, 1:9] = False
    if case == "small":
        g = g[3:33, 2:31].astype(np.float32)
        gm = gm[3:33, 2:31]
        f = f.astype(np.float32)
    return [np.asfortranarray(a) for a in (f, g, fm, gm)]


def _dft2(a):
    """a 2-D DFT by exact-order matrix products of stable twiddles (bit-stable inputs)"""
    from petakit5d.testing.synth import stable_sincos_2pi

    def mat(n):
        k = np.outer(np.arange(n), np.arange(n)) % n / n
        s, c = stable_sincos_2pi(k)
        return c - 1j * s

    my, mx = mat(a.shape[0]), mat(a.shape[1])
    return _matmul(_matmul(my, a.astype(np.complex128)), mx.T)


def _idft2(a):
    return np.conj(_dft2(np.conj(a))) / a.size


def _matmul(a, b):
    """a product in a fixed order (no BLAS: the same bits everywhere)"""
    out = np.zeros((a.shape[0], b.shape[1]), dtype=np.complex128)
    for k in range(a.shape[1]):
        out = out + a[:, k : k + 1] * b[k : k + 1, :]
    return out


def py_reg(f, g, fm, gm, ratio, use_gpu, prec):
    t, mc, c, n = M.masked_translation_registration(f, g, fm, gm, ratio, use_gpu, prec)
    return [t, mc, c, n]


def py_ncc(f, g, fm, gm, use_gpu, prec):
    c, n = M.normxcorr2_masked(f, g, fm, gm, use_gpu, prec)
    return [c, n]


def py_translate(v, t):
    return T.imtranslate(v, t)


def _blackbox_cases() -> list[Case]:
    cs = []
    for cid, seed, shape, prec, single in [
        ("perdecomp.odd.double", 51, (9, 8, 7), "double", False),
        ("perdecomp.odd.single", 51, (9, 8, 7), "single", False),
        ("perdecomp.even.double", 52, (10, 6, 6), "double", False),
        ("perdecomp.single_input", 53, (6, 7, 4), "single", True),
    ]:

        def make(seed=seed, shape=shape, single=single):
            ramp = np.arange(np.prod(shape)).reshape(shape, order="F")
            u = uniform(seed, shape, 100.0) + ramp
            return np.asfortranarray(u.astype(np.float32 if single else np.float64))

        cs.append(
            Case(
                cid,
                matlab="perdecomp_3D",
                args=[Vol(In("u")), False, prec],
                python=py_perdecomp,
                inputs={"u": Value(make)},
                out_roles=("vol",),
                compare=R1_SINGLE if prec == "single" else R1_DOUBLE,
            )
        )
    for cid, case, ratio, prec in [
        ("reg.int.double", "int", 0.5, "double"),
        ("reg.int.single", "int", 0.5, "single"),
        ("reg.frac.double", "frac", 0.5, "double"),
        ("reg.small.single", "small", 0.5, "single"),
    ]:
        names = ("f", "g", "fm", "gm")
        cs.append(
            Case(
                cid,
                matlab="MaskedTranslationRegistration2D_fit",
                args=[In("f"), In("g"), In("fm"), In("gm"), ratio, False, prec],
                nargout=4,
                python=py_reg,
                inputs={
                    n: Value(lambda case=case, i=i: reg_images(case)[i])
                    for i, n in enumerate(names)
                },
                compare=dict.fromkeys((0, 1, 2), R3_DOUBLE if prec == "double" else R3_SINGLE),
            )
        )
    cs.append(
        Case(
            "reg.window_outside",
            matlab="MaskedTranslationRegistration2D_fit",
            args=[In("f"), In("g"), In("fm"), In("gm"), 0, False, "double"],
            nargout=4,
            python=py_reg,
            inputs={
                n: Value(lambda i=i: reg_images("frac")[i])
                for i, n in enumerate(("f", "g", "fm", "gm"))
            },
            expect_error="MATLAB:badsubscript",
            py_error=(R.SimError,),
        )
    )
    cs.append(
        Case(
            "ncc.frac.double",
            matlab="normxcorr2_masked",
            args=[In("f"), In("g"), In("fm"), In("gm"), False, "double"],
            nargout=2,
            python=py_ncc,
            inputs={
                n: Value(lambda i=i: reg_images("frac")[i])
                for i, n in enumerate(("f", "g", "fm", "gm"))
            },
            compare={0: R3_DOUBLE},
        )
    )
    for cid, t in [("translate.frac", [1.3, -2.6, 0.4]), ("translate.int", [2, -1, 1]),
                   ("translate.neg", [-0.25, -3.75, -0.5])]:  # fmt: skip
        cs.append(
            Case(
                cid,
                matlab="imtranslate_function",
                args=[Vol(In("v")), t],
                python=py_translate,
                inputs={"v": Value(lambda: uniform(61, (12, 10, 8), 10.0))},
                out_roles=("vol",),
                compare=R1_DOUBLE,
            )
        )
    return cs


# ---------------------------------------------------------------------------------------------
# fminsearch and Lsq_GaussianFit_3D
# ---------------------------------------------------------------------------------------------

_FMIN = (
    "f = {code}; opts = optimset('Display', 'off'{opts}); [x, fv, ef, o] = fminsearch(f, in{{1}}, "
    "opts); out = {{x, fv, ef, o.iterations, o.funcCount}};"
)


def _rosen(x):
    return 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2


def _quad(x):
    d = x - np.array([1.0, -2, 0, 3])
    v = d * d * np.array([1.0, 2, 3, 4])
    s = 0.0
    for t in v:
        s = s + t
    return s


def _plateau(x):
    # piecewise constant: the expanded point is often as good as the reflected one (a tie that
    # decides the path)
    return float(np.floor(x[0] * x[0] + x[1] * x[1]))


def py_fmin(fun, kw):
    def run(x0):
        r = FM.fminsearch(fun, x0, **kw)
        return [r.x, r.fval, float(r.exitflag), float(r.iterations), float(r.func_count)]

    return run


def py_gfit(bead, display):
    from petakit5d.sim.otf import gaussian_fit_3d

    return gaussian_fit_3d(bead, display)


def bead(seed: int, shape=(13, 11, 9), center=(6.4, 5.7, 4.8), sigma=(1.6, 1.3, 2.1)) -> np.ndarray:
    """a Gaussian bead plus noise, MATLAB's form (y, x, z), 1-based center"""
    yy, xx, zz = np.meshgrid(*(np.arange(1, n + 1, dtype=np.float64) for n in shape), indexing="ij")
    e = (yy - center[0]) ** 2 / (2 * sigma[0] ** 2) + (xx - center[1]) ** 2 / (2 * sigma[1] ** 2)
    e = e + (zz - center[2]) ** 2 / (2 * sigma[2] ** 2)
    return np.asfortranarray(1000 * stable_exp(-e) + 30 * uniform(seed, shape))


def _fit_cases() -> list[Case]:
    rosen = "@(x) 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2"
    quad = "@(x) sum((x - [1, -2, 0, 3]).^2 .* [1, 2, 3, 4])"
    tight = ", 'MaxIter', 10000, 'MaxFunEvals', 10000, 'TolX', 1e-6, 'TolFun', 1e-6"
    cs = [
        Case("fminsearch.rosen", code=_FMIN.format(code=rosen, opts=""), args=[[-1.2, 1.0]],
             nargout=5, python=py_fmin(_rosen, {})),
        Case("fminsearch.quad", code=_FMIN.format(code=quad, opts=tight),
             args=[[0.0, 0.0, 0.0, 0.0]], nargout=5,
             python=py_fmin(_quad, {"tol_x": 1e-6, "tol_fun": 1e-6, "max_iter": 10000,
                                    "max_fun_evals": 10000})),
        Case("fminsearch.maxiter", code=_FMIN.format(code=rosen, opts=", 'MaxIter', 30"),
             args=[[-1.2, 1.0]], nargout=5, python=py_fmin(_rosen, {"max_iter": 30})),
        Case("fminsearch.zero_start", code=_FMIN.format(code=rosen, opts=""), args=[[0.0, 0.0]],
             nargout=5, python=py_fmin(_rosen, {})),
        Case("fminsearch.plateau", code=_FMIN.format(code="@(x) floor(x(1)^2 + x(2)^2)",
             opts=""), args=[[-1.2, -2.6]], nargout=5, python=py_fmin(_plateau, {})),
    ]  # fmt: skip
    for cid, seed in [("gaussfit.bead", 71), ("gaussfit.bead2", 72)]:
        cs.append(
            Case(
                cid,
                matlab="Lsq_GaussianFit_3D",
                args=[Vol(In("b")), False],
                python=py_gfit,
                inputs={"b": Value(lambda seed=seed: bead(seed))},
                compare=R2_FIT,
            )
        )
    return cs


CASES = _step_cases() + _blackbox_cases() + _fit_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# properties
# ---------------------------------------------------------------------------------------------


def test_separation_matrix_inverse_separates_the_orders():
    m = R.separation_matrix(5, 5, 1.2021, 0.232)
    inv = R.inverse_separation_matrix(m, "double")
    assert np.allclose(inv @ m, np.eye(5), atol=1e-12)
    # orders [-2 -1 0 1 2]: the zero order is all ones, the others conjugate pairs
    assert np.array_equal(m[:, 2], np.ones(5))
    assert np.allclose(m[:, 0], np.conj(m[:, 4])) and np.allclose(m[:, 1], np.conj(m[:, 3]))


def test_perdecomp_periodic_component_has_no_border_jumps():
    rng = CounterRNG(5, "c06-pd")
    u = np.cumsum(rng.uniform((6, 9, 11)), axis=-1) * 10
    p = PD.perdecomp_3d(u)
    s = u - p
    # the smooth part takes the jumps: the periodic part's jump across x is the interior's size
    jump_u = np.abs(u[..., -1] - u[..., 0]).mean()
    jump_p = np.abs(p[..., -1] - p[..., 0]).mean()
    assert jump_p < 0.2 * jump_u
    assert abs(float(s.mean())) < 1e-9 * float(np.abs(u).max())


def test_imtranslate_integer_shift_is_a_move_with_zeros():
    v = np.arange(4 * 5 * 6, dtype=np.float64).reshape(4, 5, 6)
    out = T.imtranslate(v, [2, -1, 1])
    want = np.zeros_like(v)
    want[1:, :-1, 2:] = v[:-1, 1:, :-2]
    assert np.array_equal(out, want)


def test_imtranslate_fraction_blends_neighbours_inside_and_is_zero_at_the_edge():
    v = np.zeros((1, 1, 5))
    v[0, 0] = [0, 10, 20, 30, 40]
    out = T.imtranslate(v, [0.25, 0, 0])
    # out[i] = 0.75 v[i] + 0.25 v[i - 1], 0 where v[i - 1] is outside
    assert np.allclose(out[0, 0], [0, 7.5, 17.5, 27.5, 37.5])
    assert out[0, 0, 0] == 0


def test_masked_registration_finds_the_shift_and_its_sign():
    f, g, fm, gm = (np.ascontiguousarray(a) for a in reg_images("frac"))
    t, mc, c, n = M.masked_translation_registration(f, g, fm, gm, 0.5, False, "double")
    # moved by (1.3, -0.7) in (y, x): the translation is [x y], refined by the window's centroid
    # (biased toward the window's center on these textures, as MATLAB's)
    assert abs(t[0] - -0.7) < 0.4 and abs(t[1] - 1.3) < 0.4
    assert 0.5 < mc <= 1
    assert c.shape == n.shape == (79, 71)
    # whole counts, at least eps where nothing overlaps
    assert n.dtype == np.float64 and float(n.min()) >= np.finfo(np.float64).eps
    assert np.array_equal(n[n >= 1], np.floor(n[n >= 1]))


def test_masked_registration_single_is_single():
    f, g, fm, gm = (np.ascontiguousarray(a) for a in reg_images("int"))
    t, mc, c, n = M.masked_translation_registration(f, g, fm, gm, 0.5, False, "single")
    assert t.dtype == np.float32 and c.dtype == np.float32 and n.dtype == np.float32
    assert int(mround(np.float64(t[0]))) == -2


def test_fminsearch_stops_at_max_iter_with_exitflag_0():
    r = FM.fminsearch(_rosen, [-1.2, 1.0], max_iter=5)
    assert r.exitflag == 0 and r.iterations == 5


def test_fminsearch_output_function_sees_every_iteration():
    states = []
    FM.fminsearch(_rosen, [-1.2, 1.0], max_iter=4, output_fcn=lambda *a: states.append(a[-1]))
    assert states[0] == "init" and states[-1] == "done" and states.count("iter") == 5


def test_tukey_window_tapers_x_and_y_only():
    w = R.tukey_window((3, 20, 30), 0.2)
    assert w.shape == (1, 20, 30)
    assert w[0, 10, 0] == 0 and w[0, 0, 15] == 0 and w[0, 10, 15] == 1


def test_cyl_mask_zero_order_is_the_detection_cylinder():
    dk = _dk((30, 26, 11), _PXL_DATA)
    m = R.cyl_mask(1.0, 0.605, 0.55, 0.56, 1.33, 5, (11, 30, 26), dk, True)
    assert m.shape == (5, 11, 30, 26) and m.dtype == bool
    # the zero order reaches furthest in z (its axial extent is the largest)
    assert m[2].sum() >= m[1].sum() >= m[0].sum()
    assert np.array_equal(m[1], m[3]) and np.array_equal(m[0], m[4])


@pytest.mark.parametrize("prec", ["single", "double"])
def test_fourier_shift_by_whole_voxels_moves_circularly(prec):
    v = np.zeros((4, 6, 8))
    v[1, 2, 3] = 1
    out = R.fourier_shift(v, [1, 2, 1], prec)
    assert np.unravel_index(np.argmax(np.abs(out)), out.shape) == (2, 3, 5)
    assert out.dtype == (np.complex64 if prec == "single" else np.complex128)


def test_resample_otf_onto_its_own_grid_is_the_otf():
    o = complex_volume(3, (9, 8, 5, 3))
    po = from_matlab(o)
    out = R.resample_otf(po, (0.2, 0.11, 0.11), (5, 9, 8), (0.2, 0.11, 0.11))
    assert np.allclose(to_matlab(out[0]), o, atol=1e-15)
