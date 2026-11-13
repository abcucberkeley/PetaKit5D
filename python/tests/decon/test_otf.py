"""The OTF of a PSF (petakit5d.decon.otf) against MATLAB's decon_psf2otf (the oracle tree's).

Oracle cases (inputs in MATLAB's form; the roles give the Python calls the port's form): PSFs of
odd and even sizes padded to larger sizes and to their own, images (with an image and a volume as
the output size), single, double, uint16 and logical classes, a delta, all-zero and empty PSFs,
and the errors of output sizes smaller than the PSF. Each case runs on both FFT backends:
``.exact`` (MATLAB's bits; skipped without the add-on or where the C library gives FFTW other
twiddle factors) and ``.permissive`` (decision record 0006's tolerances).

Tests without MATLAB: DeconFunctionTest's decon_psf2otf tests in intent, the placement against
padarray and circshift, and the OTF of the decon loop (``rl._psf2otf``) being this one.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d import fft as F
from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.decon import otf as O
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, PerAxis, Value, Vol
from petakit5d.testing.compare import tol
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG, decon_psf

BACKENDS = ("exact", "permissive")
PERMISSIVE = {
    "single": tol(
        reason="0006: FFT-based outputs with the permissive backend, single precision (max <= "
        "1e-5 of the largest value, relative L2 <= 1e-6)",
        rel_max=1e-5,
        rel_l2=1e-6,
    ),
    "double": tol(
        reason="0006: FFT-based outputs with the permissive backend, double precision (max <= "
        "1e-13 of the largest value, relative L2 <= 1e-14)",
        rel_max=1e-13,
        rel_l2=1e-14,
    ),
}


def rand_psf(shape, cls="single", seed=0) -> np.ndarray:
    """MATLAB's form of a random PSF (z, y, x) or (y, x) of a class"""
    u = CounterRNG(seed, "otf").uniform(shape)
    if cls == "uint16":
        a = np.floor(u * 1000).astype(np.uint16)
    elif cls == "logical":
        a = u > 0.7
    else:
        a = u.astype(np.float32 if cls == "single" else np.float64)
    return to_matlab(a)


def delta(shape) -> np.ndarray:
    a = np.zeros(shape, np.float32)
    a[tuple(s // 2 for s in shape)] = 1
    return to_matlab(a)


def py_otf(psf, out_size=None, *, backend=None, threads=None):
    return O.psf2otf(psf, out_size, backend=backend, threads=threads)


def otf_cases(cid, make, msize=None, *, cls="single", error=False):
    """decon_psf2otf(psf[, outSize]) on each backend; msize: MATLAB's outSize ([y x z])"""
    args = [Vol(In("psf"))] + ([] if msize is None else [PerAxis(msize, "yxz")])
    out = []
    for be in (None,) if error else BACKENDS:
        out.append(
            Case(
                f"{cid}.{be}" if be else cid,
                matlab="decon_psf2otf",
                args=args,
                python=py_otf,
                py_kwargs={"backend": be} if be else {},
                inputs={"psf": Value(make)},
                out_roles=("vol",),
                expect_error=error,
                py_error=ValueError if error else Exception,
                **({"compare": PERMISSIVE[cls]} if be == "permissive" else {}),
            )
        )
    return out


CASES: list[Case] = []
# psf2otfOfCenteredDeltaIsFlat, psf2otfPadsToTheOutputSize and DeconRFFTTest's shapes
CASES += otf_cases("delta", lambda: delta((7, 8, 9)))
CASES += otf_cases("delta.padded", lambda: delta((7, 8, 9)), [16, 18, 13])
CASES += otf_cases("rand5.padded", lambda: rand_psf((5, 5, 5)), [16, 18, 13])
CASES += otf_cases("rand.odd_even", lambda: rand_psf((9, 10, 11), seed=1), [20, 23, 30])
CASES += otf_cases("rand.own_size", lambda: rand_psf((9, 10, 11), seed=2))
CASES += otf_cases("rand.same_size", lambda: rand_psf((9, 10, 11), seed=2), [10, 11, 9])
CASES += otf_cases("bead.data_size", lambda: to_matlab(decon_psf((31, 31, 31))), [48, 46, 36])
CASES += otf_cases("bead.thin", lambda: to_matlab(decon_psf((31, 31, 31))[13:18]), [64, 50, 5])
CASES += otf_cases("image", lambda: rand_psf((9, 12), seed=3), [20, 31])
CASES += otf_cases("image.own_size", lambda: rand_psf((9, 12), seed=3))
CASES += otf_cases("image.into_volume", lambda: rand_psf((9, 12), seed=4), [20, 31, 6])
CASES += otf_cases("row", lambda: rand_psf((1, 7), seed=5), [1, 20])
CASES += otf_cases("double", lambda: rand_psf((9, 10, 11), "double", 6), [20, 23, 30], cls="double")
CASES += otf_cases("uint16", lambda: rand_psf((9, 10, 11), "uint16", 7), [20, 23, 30], cls="double")
CASES += otf_cases(
    "logical", lambda: rand_psf((9, 10, 11), "logical", 8), [20, 23, 30], cls="double"
)
CASES += otf_cases(
    "negative", lambda: to_matlab(-from_matlab(rand_psf((5, 6, 7), seed=9))), [16, 18, 13]
)
CASES += otf_cases("zero", lambda: to_matlab(np.zeros((5, 6, 7), np.float32)), [16, 18, 13])
CASES += otf_cases("zero.double", lambda: to_matlab(np.zeros((5, 6), np.float64)), [16, 18])
CASES += otf_cases("zero.uint16", lambda: to_matlab(np.zeros((5, 6, 7), np.uint16)), [16, 18, 13])
CASES += otf_cases("empty", lambda: np.zeros((0, 0), np.float32), [16, 18, 13])
CASES += otf_cases("neg_zero", lambda: to_matlab(-np.zeros((3, 4, 5), np.float32)), [6, 7, 8])
CASES += otf_cases("err.smaller", lambda: rand_psf((9, 9, 9)), [8, 8, 8], error=True)
CASES += otf_cases("err.smaller_z", lambda: rand_psf((9, 9, 9)), [16, 16, 8], error=True)
CASES += otf_cases("err.two_dims", lambda: rand_psf((9, 9, 9)), [16, 16], error=True)
CASES += otf_cases("err.image_smaller", lambda: rand_psf((9, 12)), [20, 11], error=True)


def _check(case, tmp_path):
    if case.id.endswith(".exact"):
        reason = digests_hold_here()
        if reason:
            pytest.skip(reason)
    oracle.check(case, tmp_path)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    _check(case, tmp_path)


_THREADED = [c for c in CASES if c.id.startswith(("rand.odd_even.", "bead.data_size."))]


@oracle.parametrize(_THREADED)
@pytest.mark.parametrize("threads", [1, 3, 24])
def test_thread_counts(case, threads, tmp_path):
    """the results don't depend on the thread count (MATLAB's references)"""
    import dataclasses

    variant = dataclasses.replace(case, py_kwargs={**(case.py_kwargs or {}), "threads": threads})
    variant._module_globals = case._module_globals
    _check(variant, tmp_path)


def _backends():
    return [b for b in BACKENDS if b in F.available()]


@pytest.mark.parametrize("be", _backends())
def test_psf2otf_of_centered_delta_is_flat(be):
    """DeconFunctionTest/psf2otfOfCenteredDeltaIsFlat"""
    sz = (7, 8, 9)
    d = np.zeros(sz, np.float32)
    d[tuple(s // 2 for s in sz)] = 1
    otf = O.psf2otf(d, backend=be)
    assert otf.shape == sz
    assert np.abs(otf.astype(np.complex128) - 1).max() <= 1e-6


@pytest.mark.parametrize("be", _backends())
def test_psf2otf_pads_to_the_output_size(be):
    """DeconFunctionTest/psf2otfPadsToTheOutputSize: the size, the sum at the origin, the PSF back
    at its place after the inverse transform, and the error for a smaller size"""
    from petakit5d.matlabcompat import fft as MF

    psf = CounterRNG(1, "pad").uniform((5, 5, 5), np.float32)
    otf = O.psf2otf(psf, (13, 16, 18), backend=be)
    assert otf.shape == (13, 16, 18)
    assert math.isclose(
        float(otf[0, 0, 0].real), math.fsum(psf.astype(np.float64).ravel()), rel_tol=1e-6
    )
    back = np.roll(MF.ifftn(otf, backend=be).real, (2, 2, 2), axis=(0, 1, 2))
    assert np.allclose(back[:5, :5, :5], psf, rtol=1e-5, atol=1e-6)
    with pytest.raises(ValueError, match="outSizeIsSmallerThanPsfSize") as e:
        O.psf2otf(CounterRNG(2).uniform((9, 9, 9), np.float32), (8, 8, 8))
    assert e.value.identifier == ""


@pytest.mark.parametrize(
    "psf_shape,out_shape",
    [
        ((5, 5, 5), (13, 16, 18)),
        ((4, 6, 8), (4, 6, 8)),
        ((9, 12), (20, 31)),
        ((9, 12), (6, 20, 31)),
        ((1, 9, 12), (20, 31)),
        ((3, 1, 1), (5, 2, 2)),
        ((2, 2, 2), (2, 3, 4)),
    ],
)
def test_pad_shift_is_padarray_and_circshift(psf_shape, out_shape):
    """the placement equals MATLAB's padarray(psf, outSize - psfSize, 'post') and
    circshift(psf, -floor(psfSize / 2)) (NumPy's pad and roll on the three dimensions)"""
    psf = CounterRNG(3, "place").uniform(psf_shape, np.float32) + np.float32(1)
    ps = (1,) * (3 - len(psf_shape)) + psf_shape
    os_ = (1,) * (3 - len(out_shape)) + out_shape
    p = np.pad(psf.reshape(ps), [(0, o - s) for o, s in zip(os_, ps, strict=True)])
    p = np.roll(p, tuple(-(s // 2) for s in ps), axis=(0, 1, 2)).reshape(out_shape)
    got = O.pad_shift_psf(psf, out_shape)
    assert got.shape == out_shape and got.dtype == psf.dtype and np.array_equal(got, p)


def test_all_zero_psf_gives_zeros_of_its_class():
    for dt in (np.float32, np.float64, np.uint16):
        z = O.psf2otf(np.zeros((3, 4, 5), dt), (6, 7, 8))
        assert z.dtype == dt and z.shape == (6, 7, 8) and not z.any()
    # -0 is zero
    assert not np.iscomplexobj(O.psf2otf(-np.zeros((3, 4, 5), np.float32)))


def test_the_decon_loop_uses_this_otf(monkeypatch):
    """decon's MATLAB loop (rl._psf2otf) gets its OTF from decon.otf.psf2otf"""
    from petakit5d.decon import rl

    calls = []
    orig = O.psf2otf

    def counted(*a, **k):
        calls.append(a[0].shape)
        return orig(*a, **k)

    monkeypatch.setattr(O, "psf2otf", counted)
    p = decon_psf((31, 31, 31))
    h = rl._psf2otf(p, (40, 48, 50), None, 2)
    assert calls == [(31, 31, 31)]
    assert np.array_equal(h, orig(p, (40, 48, 50)))


def test_an_empty_psf_gives_zeros_of_the_output_size():
    """MATLAB's [] (0 x 0) has no size to check: zeros of the output size"""
    z = O.psf2otf(np.zeros((0, 0), np.float32), (16, 18, 13))
    assert z.shape == (16, 18, 13) and z.dtype == np.float32 and not z.any()
