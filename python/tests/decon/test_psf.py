"""PSF preparation (petakit5d.decon.psf) against MATLAB's psf_gen_new, the PSF crops of RLdecon and
XR_decon_data_wrapper, decon_readtiff_cached and XR_rotate_PSF (the oracle tree's).

Oracle cases (inputs made in MATLAB's form; the roles give the Python calls the port's form):

- gen.bp.*: the PSFs of DeconBackProjectorTest's nine cases (Gaussian, elongated, 3 and 5 OTF
  lobes, sheared, on grids that aren't cubic) as DeconBackProjectorTest makes raw PSFs from them
  (``psf * 10000 + 100 + randi([0, 20])``), each with a z step ratio and method of its own;
- gen.bead.*: bead PSFs with a halo and noise (the psf_gen_masked scenario's) at z step ratios
  below, at and above 1, odd and even slice counts, down to one and to no slice left;
- gen.2d.*: images (DeconFunctionTest's hot pixel and line), negative values (MATLAB's complex
  square root), resampled to one and to no plane;
- gen.case.*: the other branches: no positive voxel in the edge slices (also with sums that
  cancel), z padding before and after, the peak outside every large component, unknown and
  upper-case method names, input classes, medFactor; gen.err.*: MATLAB's errors;
- gen.demo.*: the demo datasets' bead PSFs (copies; demodata);
- crop.*: RLdecon's crops (to the support, to the data size with its warning), also after
  psf_gen_new on the demo PSFs; read.*: decon_readtiff_cached; rotate.*: XR_rotate_PSF's files.

Cases whose PSF is resampled (z step ratio above 1) run on both FFT backends: ``.exact`` (MATLAB's
bits; skipped without the add-on or where the C library gives FFTW other twiddle factors) and
``.permissive`` (decision record 0006's single-precision tolerance: the spectrum is single).
Everything else has no FFT and is compared bit for bit.

Tests without MATLAB port DeconFunctionTest's psf_gen_new tests, DeconBackProjectorTest's
psf generation tests and DeconReadCacheTest in intent.
"""

from __future__ import annotations

import math
import os
import warnings

import numpy as np
import pytest

from petakit5d import fft as F
from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.decon import psf as P
from petakit5d.testing import oracle
from petakit5d.testing.cases import (
    Case,
    Data,
    File,
    In,
    InPath,
    PerAxis,
    Tree,
    Value,
    Vol,
    Work,
    mcell,
)
from petakit5d.testing.compare import tol
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import (
    CounterRNG,
    beads,
    decon_psf,
    stable_exp,
    stable_sincos_2pi,
    write_tiff,
)

BACKENDS = ("exact", "permissive")
# the resampling's spectrum is single precision: 0006's single tolerance
PERMISSIVE = tol(
    reason="0006: FFT-based outputs with the permissive backend, single precision (max <= 1e-5 "
    "of the largest value, relative L2 <= 1e-6)",
    rel_max=1e-5,
    rel_l2=1e-6,
)
# one plane left (new_nz = 1): the plane's inverse transform is negative but for rounding noise, so
# MATLAB's result is that noise above 0 (at most 1.6e-4); the permissive backend's noise differs at
# about 7e-8 of the raw PSF's largest value, as everywhere: proposed, 1e-6 of it (about 2,200)
PERMISSIVE_ONE_PLANE = tol(
    reason="proposed (R01): a PSF resampled to one plane is rounding noise; with the permissive "
    "backend |difference| <= 1e-6 of the raw PSF's largest value",
    abs=2.2e-3,
)
SIN_32_45 = 0.5365729180004349  # MATLAB's sind(32.45)
COS_30 = 0.8660254037844387  # MATLAB's cosd(30)

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


def bp_psf(kind: str) -> np.ndarray:
    """DeconBackProjectorTest.makePSF on grids that aren't cubic (single, (z, y, x))"""
    if kind == "gaussian":
        return decon_psf((33, 31, 29), (2.0, 1.5))
    z, y, x = _grid((48, 44, 52))
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


def bp_raw(kind: str, seed: int) -> np.ndarray:
    """psfGenerationMatchesMedfilt3's raw PSF: psf * 10000 + 100 + single(randi([0, 20]))"""
    p = bp_psf(kind)
    noise = CounterRNG(seed, "bp_noise").integers(0, 20, p.shape).astype(np.float32)
    return to_matlab((p * np.float32(10000) + np.float32(100) + noise).astype(np.float32))


def bead_psf(shape, center=None, sigma=(2.2, 1.6), seed=1, cls="uint16", halo=True) -> np.ndarray:
    """the psf_gen_masked scenario's PSF: a bead on a background of 100 with a wide dim halo and
    noise, rounded (MATLAB's form)"""
    if center is None:
        center = [s // 2 for s in shape]
    p = beads(shape, [center], sigma, 2000.0, 100.0).astype(np.float64)
    if halo:
        p = p + beads(shape, [center], (1.2, 12.0), 50.0, 0.0)
    p = p + 3.0 * CounterRNG(seed, "bead_noise").normal(shape)
    r = np.floor(p + 0.5)
    out = (
        r.astype(np.uint16)
        if cls == "uint16"
        else r.astype(np.float32 if cls == "single" else np.float64)
    )
    return to_matlab(out)


def image_2d(neg: bool = False) -> np.ndarray:
    """psfGen2DMedianMethodRemovesIsolatedPoints's image (41 x 45 here), a hot pixel and a line
    from the edge; ``neg``: single with noise below zero in the border rows"""
    im = beads((41, 45), [[17.0, 22.0]], 2.0, 5000.0, 100.0)
    im[4, 34] = 3000
    im[17, 29:] = 3000
    if neg:
        im = (
            im - np.float32(100) + np.float32(6) * CounterRNG(3, "neg").normal(im.shape, np.float32)
        )
        return to_matlab(im.astype(np.float32))
    return to_matlab(im)


def centered(shape=(31, 29, 33), center=(13, 12, 17), zero_edges=True) -> np.ndarray:
    """psfGenCentersOddSizePSF: a bead without background, the edge slices zero"""
    p = beads(shape, [center], (2.0, 1.5), 1e4, 0.0)
    if zero_edges:
        p[:5] = 0
        p[-5:] = 0
    return to_matlab(p)


def cancelling_edges() -> np.ndarray:
    """single data without a positive value in the edge slices (zero, and -0.5 in the last):
    plane 5 holds only +1 and -1, whose sum is 0, so the first nonzero plane is the bead's
    (15), and the slices are not padded (a test for any nonzero value would pad them)"""
    p = beads((31, 29, 33), [[22, 14, 16]], (2.0, 1.5), 1e4, 0.0)
    p[:15] = 0
    p[26:] = 0
    p[30] = -0.5
    p[5, 3, 4], p[5, 20, 7] = 1, -1
    return to_matlab(p.astype(np.float32))


def elongated_bead(peak_z: int) -> np.ndarray:
    """a bead elongated along z with its peak near one end (padding of z slices)"""
    shape = (41, 40, 44)
    p = beads(shape, [[peak_z, 20, 22]], (9.0, 1.6), 2000.0, 100.0).astype(np.float64)
    p = p + 3.0 * CounterRNG(7, "elong").normal(shape)
    return to_matlab(np.floor(p + 0.5).astype(np.uint16))


def unpadded(peak_z: int) -> np.ndarray:
    """a bead on a background of 1 with its peak near one end (padding without a background
    step: method 'none')"""
    p = beads((31, 30, 34), [[peak_z, 15, 17]], (2.0, 1.5), 1000.0, 1.0)
    return to_matlab(p)


def spot_2d(amp: float) -> np.ndarray:
    """a small spot in a 41 x 45 image on a background of 100 with noise: its masked component
    has 37 pixels for amp 400, between the 'median' method's 30 and the 'masked' method's 50"""
    im = beads((41, 45), [[20.0, 22.0]], 1.5, amp, 100.0)
    im = im + np.float32(3) * CounterRNG(4, "spot").normal(im.shape, np.float32)
    return to_matlab(np.floor(im + np.float32(0.5)).astype(np.float32))


def negative_background() -> np.ndarray:
    """single data whose first and last ten rows are 0 (the background the 'masked' method
    subtracts is 0) with negative values in between: the masked voxels become -0, which
    psf_gen_new keeps"""
    p = beads((31, 30, 34), [[15, 15, 17]], (2.0, 1.5), 2000.0, 0.0)
    p[:, 10:20, :] -= 5
    p[:, :10, :] = 0
    p[:, 20:, :] = 0
    p[:5, 12, 3] = 1  # positive values in the edge slices, outside the border rows
    return to_matlab(p.astype(np.float32))


def clipped_bead() -> np.ndarray:
    """a bead PSF with its background subtracted and clipped at 0: zeros in the edge slices"""
    b = from_matlab(bead_psf((41, 40, 44))).astype(np.int32) - 100
    return to_matlab(np.maximum(b, 0).astype(np.uint16))


def nan_bead() -> np.ndarray:
    """a bead PSF (single) with a NaN in its first slice, far from the bead: the plane counts as
    nonzero (its sum is NaN), so the slices are padded, and max(psf, 0) makes the NaN 0"""
    p = from_matlab(bead_psf((41, 40, 44), (24, 20, 22))).astype(np.float32)
    p[0, 5, 5] = np.float32("nan")
    return to_matlab(p)


def tied_maxima() -> np.ndarray:
    """two equal maxima in one plane, at (y, x) (2, 10) and (8, 3): MATLAB's first is (8, 3)
    (x before y)"""
    p = from_matlab(unpadded(15)).copy()
    p[15, 2, 10] = p[15, 8, 3] = 5000
    return to_matlab(p)


def small(shape, value=0.0) -> np.ndarray:
    p = np.full(shape, value, dtype=np.float32)
    if value:
        p[tuple(s // 2 for s in shape)] = 1000
    return to_matlab(p)


# ---------------------------------------------------------------------------------------------
# Python sides
# ---------------------------------------------------------------------------------------------


def py_gen(psf, dz_psf, dz_data, med_factor=1.5, method="masked", *, backend=None, threads=None):
    return P.generate_psf(
        psf, dz_psf, dz_data, med_factor, method, backend=backend, threads=threads
    )


def _check(case, tmp_path):
    if case.id.endswith(".exact"):
        reason = digests_hold_here()
        if reason:
            pytest.skip(reason)
    oracle.check(case, tmp_path)


def gen_cases(
    cid,
    make,
    dz_psf,
    dz_data,
    med=1.5,
    method="masked",
    *,
    nargs=5,
    error=False,
    marks=(),
    permissive=PERMISSIVE,
):
    """psf_gen_new(psf, dz_psf, dz_data[, medFactor[, method]]): one case, or one per backend
    when the PSF is resampled"""
    args = [Vol(In("psf")), dz_psf, dz_data, med, method][:nargs]
    fft = dz_data / dz_psf > 1 and not error
    out = []
    for be in BACKENDS if fft else (None,):
        out.append(
            Case(
                f"{cid}.{be}" if be else cid,
                matlab="psf_gen_new",
                args=args,
                python=py_gen,
                py_kwargs={"backend": be} if be else {},
                inputs={"psf": Value(make)},
                out_roles=("vol",),
                expect_error=error,
                py_error=(P.PSFError, ValueError, IndexError) if error else Exception,
                marks=marks,
                **({"compare": permissive} if be == "permissive" else {}),
            )
        )
    return out


CASES: list[Case] = []

# DeconBackProjectorTest's nine PSF cases: (kind, dz_psf, dz_data, method)
BP = {
    "gaussian": ("gaussian", 0.3, 0.3 * SIN_32_45, "masked"),
    "gaussianSkewed": ("gaussian", 0.1, 0.1, "masked"),
    "elongatedWideTaper": ("elongated", 0.1, 0.2, "masked"),
    "elongatedNarrowTaper": ("elongated", 0.2, 0.2, "median"),
    "threeLobes": ("lobes3", 0.1, 0.15, "masked"),
    "threeLobesShearedDetected": ("lobes3sheared", 0.3, 0.3 * SIN_32_45, "masked"),
    "fiveLobes": ("lobes5", 0.1, 0.25, "masked"),
    "degenerateBounds": ("gaussian", 0.1, 0.2, "median"),
    "boundsBeyondTheEdge": ("elongated", 0.1, 0.3, "masked"),
}
for _seed, (_k, (_kind, _dp, _dd, _m)) in enumerate(BP.items()):
    CASES += gen_cases(f"gen.bp.{_k}", lambda k=_kind, s=_seed: bp_raw(k, s), _dp, _dd, 1.5, _m)
# psfGenerationMatchesMedfilt3's two PSFs with both methods
for _kind in ("lobes3sheared", "elongated"):
    CASES += gen_cases(
        f"gen.bp.{_kind}.median", lambda k=_kind: bp_raw(k, 2), 0.3, 0.3 * SIN_32_45, 1.5, "median"
    )

# bead PSFs: z step ratios below, at and above 1, odd and even slice counts
for _nz in (41, 40):
    _shape = (_nz, 63, 61)
    for _r in (0.5, 1, 1.25, 2, 3.3, 30, 100):
        for _m in ("masked", "median"):
            if _m == "median" and _r not in (1, 2):
                continue
            CASES += gen_cases(
                f"gen.bead.nz{_nz}.r{_r}.{_m}",
                lambda s=_shape: bead_psf(s),
                0.1,
                0.1 * _r,
                1.5,
                _m,
                permissive=PERMISSIVE_ONE_PLANE if _r == 30 else PERMISSIVE,
            )
# the golden scenario's PSF
CASES += gen_cases("gen.bead.golden", lambda: bead_psf((41, 63, 63), (20, 31, 31)), 0.2, 0.2)
CASES += gen_cases("gen.bead.default_args", lambda: bead_psf((41, 63, 61)), 0.1, 0.1, nargs=3)
CASES += gen_cases(
    "gen.bead.median_factor2", lambda: bead_psf((41, 63, 61)), 0.1, 0.1, 2.0, "median"
)
CASES += gen_cases("gen.bead.single", lambda: bead_psf((40, 50, 47), cls="single"), 0.1, 0.2)
CASES += gen_cases("gen.bead.double", lambda: bead_psf((40, 50, 47), cls="double"), 0.1, 0.1)
CASES += gen_cases("gen.bead.upper", lambda: bead_psf((40, 50, 47)), 0.1, 0.1, 1.5, "MASKED")
CASES += gen_cases("gen.bead.none", lambda: bead_psf((40, 50, 47)), 0.1, 0.1, 1.5, "none")
CASES += gen_cases("gen.bead.rows12", lambda: bead_psf((24, 12, 47)), 0.1, 0.1)
CASES += gen_cases("gen.bead.slices5", lambda: bead_psf((5, 40, 47)), 0.1, 0.2)
# the peak's component is smaller than 300 voxels: L(peakInd) is 0
CASES += gen_cases(
    "gen.case.small_peak", lambda: bead_psf((31, 40, 44), sigma=(0.6, 0.5), halo=False), 0.1, 0.1
)
CASES += gen_cases(
    "gen.case.small_peak.median",
    lambda: bead_psf((31, 40, 44), sigma=(0.6, 0.5), halo=False),
    0.1,
    0.1,
    1.5,
    "median",
)
# no positive voxel in the edge slices: only the centering (also with sums that cancel)
CASES += gen_cases("gen.case.centered", centered, 0.2, 0.2)
CASES += gen_cases(
    "gen.case.centered_even", lambda: to_matlab(decon_psf((32, 30, 34), (2.0, 1.5))), 0.2, 0.2
)
CASES += gen_cases("gen.case.cancelling", cancelling_edges, 0.2, 0.2)
CASES += gen_cases("gen.case.cancelling.r2", cancelling_edges, 0.1, 0.2)
# -0 kept (a background of 0), zeros in the edge slices for the 'median' method's median
CASES += gen_cases("gen.case.negative_zero", negative_background, 0.2, 0.2)
CASES += gen_cases("gen.case.median_clipped", clipped_bead, 0.1, 0.1, 1.5, "median")
# NaN (a plane whose sum is NaN; max(psf, 0)), tied maxima
CASES += gen_cases("gen.case.nan_masked", nan_bead, 0.1, 0.1)
CASES += gen_cases("gen.case.nan_none", nan_bead, 0.1, 0.1, 1.5, "none")
CASES += gen_cases("gen.case.tied_maxima", tied_maxima, 0.1, 0.1, 1.5, "none")
# z padding after and before the slices
CASES += gen_cases("gen.case.pad_post", lambda: elongated_bead(31), 0.1, 0.1)
CASES += gen_cases("gen.case.pad_pre", lambda: elongated_bead(9), 0.1, 0.1)
CASES += gen_cases("gen.case.pad_post.r2", lambda: elongated_bead(31), 0.1, 0.2)
CASES += gen_cases("gen.case.pad_post.none", lambda: unpadded(24), 0.1, 0.1, 1.5, "none")
CASES += gen_cases("gen.case.pad_pre.none", lambda: unpadded(5), 0.1, 0.1, 1.5, "none")
CASES += gen_cases(
    "gen.case.pad_pre.none.even",
    lambda: to_matlab(from_matlab(unpadded(5))[:30]),
    0.1,
    0.1,
    1.5,
    "none",
)
# images
CASES += gen_cases("gen.2d.masked", image_2d, 0.2, 0.2)
CASES += gen_cases("gen.2d.median", image_2d, 0.2, 0.2, 1.5, "median")
CASES += gen_cases("gen.2d.masked.negative", lambda: image_2d(True), 0.2, 0.2)
CASES += gen_cases("gen.2d.median.negative", lambda: image_2d(True), 0.2, 0.2, 1.5, "median")
CASES += gen_cases("gen.2d.masked.small_spot", lambda: spot_2d(400.0), 0.2, 0.2)
CASES += gen_cases("gen.2d.median.small_spot", lambda: spot_2d(400.0), 0.2, 0.2, 1.5, "median")
CASES += gen_cases("gen.2d.masked.r2", image_2d, 0.1, 0.2)
CASES += gen_cases("gen.2d.masked.r3", image_2d, 0.1, 0.3)
CASES += gen_cases("gen.2d.none", image_2d, 0.2, 0.2, 1.5, "none")
# MATLAB's errors
CASES += gen_cases("gen.err.rows8", lambda: bead_psf((24, 8, 47)), 0.1, 0.1, error=True)
CASES += gen_cases("gen.err.slices3", lambda: bead_psf((3, 40, 47)), 0.1, 0.1, error=True)
CASES += gen_cases("gen.err.slices4", lambda: bead_psf((4, 40, 47)), 0.1, 0.1, error=True)
CASES += gen_cases("gen.err.zero", lambda: small((20, 30, 31)), 0.1, 0.1, error=True)
CASES += gen_cases(
    "gen.err.2d.rows8", lambda: to_matlab(from_matlab(image_2d())[:8]), 0.1, 0.1, error=True
)
CASES += gen_cases(
    "gen.err.2d.median.zero", lambda: small((30, 31)), 0.1, 0.1, 1.5, "median", error=True
)

# the demo datasets' bead PSFs (copies), as RLdecon reads them
DEMO = {
    "rw488": ("PetaKit5D_demo_cell_image_dataset/PSF/488_2_c.tif", 0.1, 0.3 * SIN_32_45),
    "rw560": ("PetaKit5D_demo_cell_image_dataset/PSF/560_c.tif", 0.1, 0.3 * SIN_32_45),
    "confocal": (
        "PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/Confocal/Confocal_PSF.tif",
        0.2,
        0.2,
    ),
    "widefield": (
        "PetaKit5D_2P_Confocal_Phase_Widefield_demo_datasets/Widefield/WF_PSF.tif",
        0.1,
        0.2,
    ),
}


def py_gen_file(path, dz_psf, dz_data, method, *, backend=None):
    from petakit5d.io import read_tiff

    return P.generate_psf(read_tiff(path), dz_psf, dz_data, 1.5, method, backend=backend)


def py_gen_crop_file(path, dz_psf, dz_data, method, msize, *, backend=None):
    from petakit5d.io import read_tiff

    p = P.generate_psf(
        np.asarray(read_tiff(path), np.float32), dz_psf, dz_data, 1.5, method, backend=backend
    )
    return P.crop_psf_to_data(P.crop_psf_support(p), msize)


_RLDECON_CROP = (
    "psf = single(readtiff(in{1})); imSize = in{5}; "
    "psf = psf_gen_new(psf, in{2}, in{3}, 1.5, in{4}); pSize = size(psf, 1 : 3); "
    "py = find(squeeze(sum(psf, [2, 3]))); px = find(squeeze(sum(psf, [1, 3]))); "
    "pz = find(squeeze(sum(psf, [1, 2]))); "
    "cropSz = [min(py(1) - 1, pSize(1) - py(end)), min(px(1) - 1, pSize(2) - px(end)), "
    "min(pz(1) - 1, pSize(3) - pz(end))] - 1; cropSz = max(0, cropSz); "
    "bbox = [cropSz + 1, pSize - cropSz]; "
    "psf = psf(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6)); pSize = size(psf, 1 : 3); "
    "if any(pSize > imSize), "
    "warning('The psf size %s is larger than the data size %s, crop it!', mat2str(pSize), "
    "mat2str(imSize)); s = max(0, floor((pSize - imSize) / 2)) + 1; "
    "t = s + min(pSize, imSize) - 1; psf = psf(s(1) : t(1), s(2) : t(2), s(3) : t(3)); end; "
    "out = {psf};"
)

for _k, (_rel, _dp, _dd) in DEMO.items():
    for _m in ("masked", "median"):
        fft = _dd / _dp > 1
        for be in BACKENDS if fft else (None,):
            CASES.append(
                Case(
                    f"gen.demo.{_k}.{_m}" + (f".{be}" if be else ""),
                    code="out = {psf_gen_new(single(readtiff(in{1})), in{2}, in{3}, 1.5, in{4})};",
                    args=[InPath("psf"), _dp, _dd, _m],
                    python=py_gen_file,
                    py_kwargs={"backend": be} if be else {},
                    inputs={"psf": Data("$PETAKIT5D_DEMO_DATA/" + _rel, "psf.tif", content=True)},
                    out_roles=("vol",),
                    files=False,
                    **({"compare": PERMISSIVE} if be == "permissive" else {}),
                )
            )
    # RLdecon's PSF preparation for the simplified method, data of 101 slices
    for be in BACKENDS if _dd / _dp > 1 else (None,):
        CASES.append(
            Case(
                f"crop.demo.{_k}" + (f".{be}" if be else ""),
                code=_RLDECON_CROP,
                args=[InPath("psf"), _dp, _dd, "masked", PerAxis([256, 300, 101], "yxz")],
                python=py_gen_crop_file,
                py_kwargs={"backend": be} if be else {},
                inputs={"psf": Data("$PETAKIT5D_DEMO_DATA/" + _rel, "psf.tif", content=True)},
                out_roles=("vol",),
                files=False,
                warnings="ids",
                **({"compare": PERMISSIVE} if be == "permissive" else {}),
            )
        )

# ---------------------------------------------------------------------------------------------
# RLdecon's and XR_decon_data_wrapper's crops
# ---------------------------------------------------------------------------------------------

_SUPPORT = (
    "psf = in{1}; pSize = size(psf, 1 : 3); "
    "py = find(squeeze(sum(psf, [2, 3]))); px = find(squeeze(sum(psf, [1, 3]))); "
    "pz = find(squeeze(sum(psf, [1, 2]))); "
    "cropSz = [min(py(1) - 1, pSize(1) - py(end)), min(px(1) - 1, pSize(2) - px(end)), "
    "min(pz(1) - 1, pSize(3) - pz(end))] - 1; cropSz = max(0, cropSz); "
    "bbox = [cropSz + 1, pSize - cropSz]; "
    "psf = psf(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6)); out = {psf};"
)
_WRAPPER_SUPPORT = (
    "psf = in{1}; py = find(squeeze(sum(psf, [2, 3]))); px = find(squeeze(sum(psf, [1, 3]))); "
    "pz = find(squeeze(sum(psf, [1, 2]))); "
    "cropSz = [min(py(1) - 1, size(psf, 1) - py(end)), min(px(1) - 1, size(psf, 2) - px(end)), "
    "min(pz(1) - 1, size(psf, 3) - pz(end))] - 1; cropSz = max(0, cropSz); "
    "bbox = [cropSz + 1, size(psf, 1 : 3) - cropSz]; "
    "psf = psf(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6)); out = {psf};"
)
_TO_DATA = (
    "psf = in{1}; imSize = in{2}; pSize = size(psf, 1 : 3); "
    "if any(pSize > imSize), "
    "warning('The psf size %s is larger than the data size %s, crop it!', mat2str(pSize), "
    "mat2str(imSize)); s = max(0, floor((pSize - imSize) / 2)) + 1; "
    "t = s + min(pSize, imSize) - 1; psf = psf(s(1) : t(1), s(2) : t(2), s(3) : t(3)); end; "
    "out = {psf};"
)


def offcenter(shape, center, neg=False) -> np.ndarray:
    """a PSF with zero planes around a bead off the center (and negative values whose sums along
    some planes cancel)"""
    p = np.zeros(shape, np.float32)
    cz, cy, cx = center
    p[cz - 3 : cz + 4, cy - 2 : cy + 3, cx - 4 : cx + 2] = beads(
        (7, 5, 6), [[3, 2, 4]], (1.5, 1.0), 100.0
    )
    if neg:
        # row y 0 sums to 0 (so the support along y starts at the bead); x column 0 and z
        # plane 0 don't
        p[0, 0, 0], p[3, 0, 5] = 7, -7
    return to_matlab(p)


SUPPORT_INPUTS = {
    "bead": lambda: to_matlab(decon_psf((33, 31, 35), (2.0, 1.5))),
    "offcenter": lambda: offcenter((21, 25, 23), (8, 15, 9)),
    "offcenter_neg": lambda: offcenter((21, 25, 23), (8, 15, 9), True),
    "full": lambda: bead_psf((11, 12, 13)),
    "image": image_2d,
    "double": lambda: to_matlab(decon_psf((33, 31, 35), (2.0, 1.5)).astype(np.float64)),
}
for _k, _make in SUPPORT_INPUTS.items():
    for _name, _code in (("support", _SUPPORT), ("wrapper", _WRAPPER_SUPPORT)):
        CASES.append(
            Case(
                f"crop.{_name}.{_k}",
                code=_code,
                args=[Vol(In("psf"))],
                python=P.crop_psf_support,
                inputs={"psf": Value(_make)},
                out_roles=("vol",),
            )
        )
CASES.append(
    Case(
        "crop.support.zero",
        code=_SUPPORT,
        args=[Vol(In("psf"))],
        python=P.crop_psf_support,
        inputs={"psf": Value(lambda: small((9, 10, 11)))},
        out_roles=("vol",),
        expect_error=True,
        py_error=P.PSFError,
    )
)
for _k, _msize in {
    "smaller_z": [64, 64, 21],
    "smaller_all": [20, 25, 11],
    "odd_even": [30, 30, 30],
    "larger": [128, 128, 128],
    "image": [16, 20, 1],
}.items():
    CASES.append(
        Case(
            f"crop.to_data.{_k}",
            code=_TO_DATA,
            args=[Vol(In("psf")), PerAxis(_msize, "yxz")],
            python=P.crop_psf_to_data,
            inputs={
                "psf": Value(
                    image_2d if _k == "image" else (lambda: to_matlab(decon_psf((33, 31, 35))))
                )
            },
            out_roles=("vol",),
            warnings="ids",
        )
    )

# ---------------------------------------------------------------------------------------------
# decon_readtiff_cached, XR_rotate_PSF
# ---------------------------------------------------------------------------------------------


def _tiff(make):
    def write(path):
        write_tiff(path, from_matlab(make()))

    return write


CASES.append(
    Case(
        "read.cached",
        code="a = decon_readtiff_cached(in{1}); b = decon_readtiff_cached(in{1}); "
        "out = {a, isequal(a, b)};",
        args=[InPath("psf")],
        nargout=2,
        python=lambda p: [
            P.read_tiff_cached(p),
            bool(np.array_equal(P.read_tiff_cached(p), P.read_tiff_cached(p))),
        ],
        inputs={"psf": File(_tiff(lambda: bead_psf((21, 30, 31))), "psf.tif")},
        out_roles=("vol", None),
        files=False,
    )
)


def rotate_case(cid, make, *opts, existing=False):
    """XR_rotate_PSF(<work>/psfs/psf.tif, ...) (the Rotated folder's files)"""

    def tree(folder):
        write_tiff(os.path.join(folder, "psf.tif"), from_matlab(make()))
        if existing:
            os.makedirs(os.path.join(folder, "Rotated"))
            write_tiff(os.path.join(folder, "Rotated", "psf.tif"), np.full((3, 4, 5), 7, np.uint16))

    return Case(
        cid,
        code="XR_rotate_PSF([in{1}, 'psf.tif'], in{2 : end}); out = {};",
        nargout=0,
        args=[InPath("d"), *opts],
        python=py_rotate,
        inputs={"d": Tree(tree, "psfs")},
    )


def py_rotate(d, *opts):
    from petakit5d import matlab

    matlab.XR_rotate_PSF(d + "psf.tif", *opts)
    return []


_ROT = ["xyPixelSize", 0.108, "dz", 0.3, "SkewAngle", 32.45, "uuid", "u0"]
_ROT_PSF = lambda: bead_psf((31, 36, 34), sigma=(2.0, 1.5), halo=False)  # noqa: E731
CASES += [
    rotate_case("rotate.default", _ROT_PSF, *_ROT),
    rotate_case("rotate.reverse", _ROT_PSF, *_ROT, "Reverse", True),
    rotate_case("rotate.objective_scan", _ROT_PSF, *_ROT, "objectiveScan", True),
    rotate_case("rotate.single", _ROT_PSF, *_ROT, "Save16bit", False, "HighPrctile", 90),
    rotate_case(
        "rotate.ties",
        lambda: to_matlab(
            np.full((31, 36, 34), 100, np.uint16)
            + from_matlab(bead_psf((31, 36, 34), halo=False)) // 1000
        ),
        *_ROT,
        "HighPrctile",
        50,
    ),
    rotate_case("rotate.existing", _ROT_PSF, *_ROT, existing=True),
    rotate_case("rotate.overwrite", _ROT_PSF, *_ROT, "Overwrite", True, existing=True),
]

# GoldenOutputTest's psf_gen_masked scenario: reference-only cases with the MATLAB calls (and so the
# keys and references) of tests/matlab_map/test_golden.py's cases (F06)
GOLDEN_CASES = {
    f"golden.psf_gen_masked.{what}": Case(
        f"golden.psf_gen_masked.{what}",
        matlab="pkpy_golden",
        args=[what, "psf_gen_masked", mcell("psf"), Work()],
        nargout=1,
        outputs=["psf"],
        helpers=["tests/matlab_map/matlab"],
        files=False,
    )
    for what in ("stored", "compute", "inputs")
}


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    _check(case, tmp_path)


@oracle.parametrize(list(GOLDEN_CASES.values()))
def test_golden_references(case, tmp_path):
    """the references of the psf_gen_masked scenario exist (store or manifest)"""
    oracle.check(case, tmp_path)


def _rel_l2(a, e) -> float:
    a = np.asarray(a, np.float64).ravel()
    e = np.asarray(e, np.float64).ravel()
    return math.sqrt(math.fsum((a - e) ** 2)) / max(math.sqrt(math.fsum(e**2)), 2.0**-52)


@pytest.mark.oracle
def test_golden_scenario():
    """GoldenOutputTest/outputMatchesReference(scenario=psf_gen_masked): psf_gen of the scenario's
    input equals the oracle tree's result bit for bit (no resampling, so no FFT) and the stored
    reference within GoldenOutputTest's 1e-4"""
    computed = oracle.reference(GOLDEN_CASES["golden.psf_gen_masked.compute"])
    inputs = oracle.reference(GOLDEN_CASES["golden.psf_gen_masked.inputs"])
    if computed.source != "store" or inputs.source != "store":
        pytest.skip("needs the golden store (the manifests hold hashes only)")
    out = P.generate_psf(from_matlab(inputs.outputs[0].data), 0.2, 0.2, 1.5, "masked")
    m = to_matlab(np.asarray(out, np.float32))
    ref = np.asarray(computed.outputs[0].data)
    assert m.dtype == ref.dtype and m.shape == ref.shape
    assert m.tobytes(order="F") == ref.tobytes(order="F")
    stored = oracle.reference(GOLDEN_CASES["golden.psf_gen_masked.stored"])
    if stored.source == "store":
        assert _rel_l2(m, stored.outputs[0].data) <= 1e-4


_THREADED = [
    c for c in CASES if c.id.startswith(("gen.bp.", "gen.bead.nz41.r2.masked", "gen.2d.masked."))
]


@oracle.parametrize(_THREADED)
@pytest.mark.parametrize("threads", [1, 3, 24])
def test_thread_counts(case, threads, tmp_path):
    """the results don't depend on the thread count (MATLAB's references)"""
    import dataclasses

    variant = dataclasses.replace(case, py_kwargs={**(case.py_kwargs or {}), "threads": threads})
    variant._module_globals = case._module_globals
    _check(variant, tmp_path)


# ---------------------------------------------------------------------------------------------
# DeconFunctionTest's psf_gen_new tests (in intent)
# ---------------------------------------------------------------------------------------------


def test_psf_gen_centers_odd_size_psf():
    """psfGenCentersOddSizePSF: without background (zero edge slices) the PSF is only shifted
    circularly so that its peak is at the center"""
    psf = from_matlab(centered())
    out = P.psf_gen(psf, 0.2, 0.2)
    assert out.dtype == np.float32
    assert np.array_equal(out, np.roll(psf, (15 - 13, 14 - 12, 16 - 17), axis=(0, 1, 2)))


def test_psf_gen_removes_background():
    """psfGenRemovesBackground: 'masked' keeps the component around the peak and subtracts the
    background of the border rows"""
    shape = (41, 43, 39)
    psf = beads(shape, [[20, 21, 19]], (2.0, 1.5), 5000.0, 100.0)
    psf = psf + np.float32(3) * CounterRNG(5, "bg").normal(shape, np.float32)
    ideal = beads(shape, [[20, 21, 19]], (2.0, 1.5), 5000.0, 0.0)
    out = P.psf_gen(psf, 0.2, 0.2)
    assert out.min() >= 0
    assert not out[:, :5, :].any(), "background far from the peak remains"
    assert math.isclose(float(out.max()), 5000, rel_tol=0.03)
    assert math.isclose(
        math.fsum(out.astype(np.float64).ravel()),
        math.fsum(ideal.astype(np.float64).ravel()),
        rel_tol=0.05,
    )


def test_psf_gen_2d_median_method_removes_isolated_points():
    """psfGen2DMedianMethodRemovesIsolatedPoints: the 2D 'median' method subtracts medFactor times
    the median of the positive pixels and keeps the component around the peak, without the hot
    pixel or the line from the edge, and moves the peak to the center"""
    psf = from_matlab(image_2d())
    out = P.psf_gen(psf, 0.2, 0.2, 1.5, "median")
    from petakit5d.matlabcompat import stats

    core = np.maximum(psf - np.float32(1.5 * float(stats.median(psf[psf > 0]))), 0).astype(
        np.float32
    )
    expected = np.roll(core, (20 - 17, 22 - 22), axis=(0, 1))
    assert out.shape == (41, 45)
    assert out[20, 22] == core.max()
    assert np.all((out == expected) | (out == 0))
    assert out[4 + 3, 34] == 0, "the isolated hot pixel remains"
    assert not out[20, 30:39].any(), "the line from the PSF remains"
    core[4, 34] = 0
    core[17, 29:] = 0
    assert math.isclose(
        math.fsum(out.astype(np.float64).ravel()),
        math.fsum(core.astype(np.float64).ravel()),
        rel_tol=0.01,
    )


@pytest.mark.parametrize("be", BACKENDS)
def test_psf_gen_resamples_z_to_the_data_step(be):
    """psfGenResamplesZToTheDataStep: a PSF sampled at half the data's z step is resampled by
    Fourier truncation to half the slices and half the width; no resampling for a smaller data
    step"""
    if be not in F.available():
        pytest.skip(f"the {be} FFT backend is not installed")
    psf = decon_psf((31, 31, 33), (2.0, 1.5))
    out = P.generate_psf(psf, 0.1, 0.2, backend=be)
    assert out.shape == (16, 31, 33) and out.dtype == np.float64

    def sd(v):
        z = np.arange(1, v.size + 1, dtype=np.float64)
        c = math.fsum(z * v) / math.fsum(v)
        return math.sqrt(math.fsum((z - c) ** 2 * v) / math.fsum(v))

    s_in = sd(psf.astype(np.float64).sum(axis=(1, 2)))
    s_out = sd(out.sum(axis=(1, 2)))
    assert math.isclose(s_out, s_in / 2, rel_tol=0.1)
    assert np.array_equal(P.psf_gen(psf, 0.2, 0.1), psf)


def test_psf_gen_keeps_even_size_psf_at_the_fft_center():
    """psfGenKeepsEvenSizePSFAtTheFFTCenter: the peak of an even-size PSF ends at size // 2"""
    psf = decon_psf((32, 30, 34), (2.0, 1.5))
    out = P.psf_gen(psf, 0.2, 0.2)
    from petakit5d.matlabcompat import reduce

    assert reduce.max(out)[1] == (16, 15, 17)


def test_psf_generation_matches_medfilt3_order_filter():
    """DeconBackProjectorTest/psfGenerationMatchesMedfilt3 and psfGenerationDoesNotCallMedfilt3 in
    intent: the 'masked' median of single volumes (the mex's selection network) equals the
    general order filter (medfilt3's code), so the PSFs are the same either way"""
    from petakit5d.matlabcompat import morphology as M

    for kind in ("lobes3sheared", "elongated"):
        raw = from_matlab(bp_raw(kind, 2))
        fast = M.medfilt3(raw, (3, 3, 3))
        slow = M._order_filter(raw, (3, 3, 3), 13, "symmetric", None)
        assert np.array_equal(fast, slow)


def test_psf_gen_does_not_depend_on_the_thread_count():
    psf = from_matlab(bead_psf((41, 63, 61)))
    ref = P.generate_psf(psf, 0.1, 0.1, threads=1)
    for t in (2, 5, 24):
        assert np.array_equal(P.generate_psf(psf, 0.1, 0.1, threads=t), ref)


def test_psf_gen_keeps_the_input():
    psf = from_matlab(bead_psf((41, 63, 61), cls="single"))
    before = psf.copy()
    P.psf_gen(psf, 0.1, 0.2)
    assert np.array_equal(psf, before)


def test_component_of_exactly_the_area_is_kept():
    """bwareaopen keeps a component of exactly P voxels: L(peakInd) is its label at 300 voxels,
    and 0 (every voxel outside the large components) at 299"""
    for n in (300, 299):
        first = np.zeros((20, 30, 40), bool)
        first.reshape(-1)[:n] = True  # one component along the first rows
        mask = first.copy()
        mask[15:18, 5:8, 5:8] = True  # another one, of 27 voxels
        bw = P._peak_component(mask, (0, 0, 0), 300, 3)
        assert np.array_equal(bw, first if n == 300 else np.ones_like(mask))


def test_closing_on_the_box_equals_the_whole_image():
    """the closing computed on the mask's box (widened by the neighborhood) equals MATLAB's
    imclose of the whole mask, also where the box meets the border"""
    from petakit5d.matlabcompat import morphology as M

    rng = CounterRNG(11, "close")
    for shape in ((30, 40, 44), (1, 40, 44)):
        for lo, hi in (
            ((5, 8, 9), (20, 30, 35)),
            ((0, 0, 0), (12, 9, 44)),
            ((22, 30, 30), (30, 40, 44)),
        ):
            bw = np.zeros(shape, bool)
            sl = tuple(
                slice(min(a, s - 1), min(b, s)) for a, b, s in zip(lo, hi, shape, strict=True)
            )
            bw[sl] = rng.uniform(bw[sl].shape) > 0.6
            se = M.strel("sphere", 3) if shape[0] > 1 else M.strel("disk", 3)
            whole = M.imclose(bw if shape[0] > 1 else bw[0], se).reshape(shape)
            assert np.array_equal(P._close(bw, 4), whole)


def _sparse(shape, lo, hi, frac, seed):
    """a volume of zeros with positive values (and some -0) in a box"""
    rng = CounterRNG(seed, "sparse")
    p = np.zeros(shape, np.float32)
    sl = tuple(slice(a, b) for a, b in zip(lo, hi, strict=True))
    sub = p[sl]
    u = rng.uniform(sub.shape)
    sub[u < frac] = (rng.uniform(sub.shape, np.float32) * 100 + 1)[u < frac]
    sub[(u > 0.95)] = -0.0
    return p


@pytest.mark.parametrize(
    "shape,lo,hi,frac",
    [
        ((30, 40, 44), (10, 12, 15), (18, 25, 22), 0.7),  # inside
        ((30, 40, 44), (0, 0, 0), (6, 7, 9), 0.8),  # the first corner
        ((30, 40, 44), (25, 33, 38), (30, 40, 44), 0.8),  # the last corner
        ((30, 40, 44), (0, 10, 40), (30, 13, 44), 0.75),  # through z, at the x border
        ((30, 40, 44), (10, 12, 15), (18, 25, 22), 0.2),  # no window with 14 positives
        ((30, 40, 44), (0, 0, 0), (30, 40, 44), 0.6),  # dense: medfilt3 itself
        ((3, 5, 4), (0, 0, 0), (3, 5, 4), 0.9),
        ((30, 40, 44), (0, 0, 0), (0, 0, 0), 0.5),  # no positive value
    ],
)
def test_median_mask_and_peak_equals_medfilt3(shape, lo, hi, frac):
    """the 'median' method's mask and peak from the windows with 14 or more positive values equal
    medfilt3's (symmetric padding) > 0 and its first maximum"""
    from petakit5d.matlabcompat import morphology as M
    from petakit5d.matlabcompat import reduce

    for seed in range(3):
        p = _sparse(shape, lo, hi, frac, seed)
        med = M.medfilt3(p.astype(np.float64), (3, 3, 3))
        mask, peak = P._median_mask_and_peak(p, 4)
        assert np.array_equal(mask, med > 0)
        assert peak == reduce.max(med)[1]


@pytest.mark.parametrize("nz", [9, 10, 41, 40, 1, 2])
@pytest.mark.parametrize("ratio", [1.01, 1.5, 2, 2.5, 3.3, 7, 30, 100])
def test_truncation_index_arithmetic(nz, ratio, monkeypatch):
    """the planes of the truncated spectrum as MATLAB's uint16 index arithmetic gives them: the
    first ceil(new_nz / 2) planes, then the last ones, the second over the first"""
    psf = np.zeros((nz, 4, 5), np.float32)
    psf[:, 1, 2] = np.arange(1, nz + 1)
    from petakit5d.matlabcompat import fft as MF

    spec = MF.fftn(psf)
    seen = {}

    def ifftn(x, **kw):
        seen["trunc"] = x
        return orig(x, **kw)

    orig = MF.ifftn
    monkeypatch.setattr(MF, "ifftn", ifftn)
    P._truncate_z(psf, ratio, None, 1)
    t = seen["trunc"]
    new_nz = int(math.floor(nz / ratio + 0.5))
    h = (new_nz + 1) // 2
    exp = np.zeros((new_nz, 4, 5), np.complex128)
    exp[:h] = spec[:h]
    if new_nz:
        exp[new_nz - h :] = spec[nz - h :]
    assert t.shape == exp.shape and np.array_equal(t, exp)


# ---------------------------------------------------------------------------------------------
# the crops
# ---------------------------------------------------------------------------------------------


def test_crop_psf_support_keeps_the_center():
    psf = from_matlab(offcenter((21, 25, 23), (8, 15, 9)))
    out = P.crop_psf_support(psf)
    # nonzero planes z 5..11, y 13..17, x 5..10: min(first, n - 1 - last) - 1 cut on both sides
    assert out.shape == (21 - 2 * 4, 25 - 2 * 6, 23 - 2 * 4)
    assert np.array_equal(out, psf[4:17, 6:19, 4:19])
    # a row whose values cancel has no support, also next to the bead
    neg = from_matlab(offcenter((21, 25, 23), (8, 15, 9), True))
    assert P.crop_psf_support(neg).shape == (21, 13, 23)


def test_crop_psf_to_data_warns_and_cuts_about_the_center():
    psf = decon_psf((33, 31, 35))
    with pytest.warns(
        UserWarning, match=r"The psf size \[31 35 33\] is larger than the data size \[64 64 21\]"
    ):
        out = P.crop_psf_to_data(psf, (21, 64, 64))
    assert np.array_equal(out, psf[6:27])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert P.crop_psf_to_data(psf, (40, 40, 40)) is not None


# ---------------------------------------------------------------------------------------------
# DeconReadCacheTest (decon_readtiff_cached)
# ---------------------------------------------------------------------------------------------


@pytest.fixture
def fresh_cache():
    P.clear_read_cache()
    yield
    P.clear_read_cache()


def _stub_reads(monkeypatch):
    """readtiff stubbed: zeros(2, 'single') for every read, as pktest.FunctionStubFixture"""
    import petakit5d.io

    calls = []

    def stub(fn, *a, **k):
        calls.append(fn)
        return np.zeros((2, 2), np.float32)

    monkeypatch.setattr(petakit5d.io, "read_tiff", stub)
    return calls


def test_back_projector_is_read_once(tmp_path, fresh_cache, monkeypatch):
    """backProjectorIsReadOnce: the second read of a file comes from the cache"""
    fn = tmp_path / "PSF_back_projector_alpha_0.005000_otf_cum_thresh_0.900000.tif"
    bp = CounterRNG(1, "bp").uniform((7, 9, 8), np.float32)
    write_tiff(fn, bp)
    assert np.array_equal(P.read_tiff_cached(fn), bp)
    calls = _stub_reads(monkeypatch)
    assert np.array_equal(P.read_tiff_cached(fn), bp)
    assert calls == []


def test_changed_files_are_read_again(tmp_path, fresh_cache, monkeypatch):
    """changedFilesAreReadAgain: another size, or the same size with another modification time"""
    fn = tmp_path / "psf.tif"
    a = CounterRNG(2, "a").uniform((7, 9, 8), np.float32)
    write_tiff(fn, a)
    assert np.array_equal(P.read_tiff_cached(fn), a)
    b = CounterRNG(2, "b").uniform((8, 9, 8), np.float32)
    write_tiff(fn, b)
    assert np.array_equal(P.read_tiff_cached(fn), b)
    assert np.array_equal(P.read_tiff_cached(fn), b)
    st = os.stat(fn)
    os.utime(fn, ns=(st.st_atime_ns, st.st_mtime_ns + 5_000_000_000))
    _stub_reads(monkeypatch)
    assert np.array_equal(P.read_tiff_cached(fn), np.zeros((2, 2), np.float32))


def test_several_files_are_kept(tmp_path, fresh_cache, monkeypatch):
    """severalFilesAreKept: the last four files read stay in memory, the first is read again"""
    fns = [tmp_path / f"f{i}.tif" for i in range(1, 6)]
    ims = []
    for i, fn in enumerate(fns):
        ims.append(
            (CounterRNG(3, i).uniform((3, 4, 5), np.float32) + np.float32(i + 1)).astype(np.float32)
        )
        write_tiff(fn, ims[-1])
        assert np.array_equal(P.read_tiff_cached(fn), ims[-1])
    calls = _stub_reads(monkeypatch)
    for fn, im in zip(fns[1:], ims[1:], strict=True):
        assert np.array_equal(P.read_tiff_cached(fn), im)
    assert calls == []
    assert np.array_equal(P.read_tiff_cached(fns[0]), np.zeros((2, 2), np.float32))


def test_a_hit_does_not_refresh_the_order(tmp_path, fresh_cache, monkeypatch):
    """MATLAB keeps the files in the order they were read: reading a cached file again doesn't
    keep it longer"""
    fns = [tmp_path / f"f{i}.tif" for i in range(5)]
    for i, fn in enumerate(fns[:4]):
        write_tiff(fn, np.full((2, 3, 4), i, np.uint16))
        P.read_tiff_cached(fn)
    P.read_tiff_cached(fns[0])  # a hit
    write_tiff(fns[4], np.full((2, 3, 4), 4, np.uint16))
    P.read_tiff_cached(fns[4])
    calls = _stub_reads(monkeypatch)
    P.read_tiff_cached(fns[0])
    assert calls == [os.fspath(fns[0])]


def test_missing_files_error(tmp_path, fresh_cache):
    """missingFilesError"""
    with pytest.raises(OSError):
        P.read_tiff_cached(tmp_path / "missing.tif")
    with pytest.raises(OSError):
        P.read_tiff_cached(tmp_path)


def test_cached_arrays_are_read_only(tmp_path, fresh_cache):
    fn = tmp_path / "psf.tif"
    write_tiff(fn, np.ones((3, 4, 5), np.float32))
    a = P.read_tiff_cached(fn)
    with pytest.raises(ValueError):
        a[0, 0, 0] = 2
    assert P.read_tiff_cached(fn)[0, 0, 0] == 1


def test_cache_is_safe_across_threads(tmp_path, fresh_cache):
    from concurrent.futures import ThreadPoolExecutor

    fns = [tmp_path / f"f{i}.tif" for i in range(6)]
    for i, fn in enumerate(fns):
        write_tiff(fn, np.full((2, 3, 4), i, np.uint16))
    with ThreadPoolExecutor(8) as ex:
        res = list(ex.map(lambda k: int(P.read_tiff_cached(fns[k % 6])[0, 0, 0]), range(200)))
    assert res == [k % 6 for k in range(200)]


# ---------------------------------------------------------------------------------------------
# the MATLAB-name layer
# ---------------------------------------------------------------------------------------------


def test_matlab_names_take_and_give_matlabs_values():
    """petakit5d.matlab.psf_gen_new takes and gives MATLAB's arrays"""
    from petakit5d import matlab

    m = bead_psf((21, 30, 31))
    out = matlab.psf_gen_new(m, 0.1, 0.1)
    assert np.array_equal(out, to_matlab(P.psf_gen(from_matlab(m), 0.1, 0.1)))
    out = matlab.psf_gen_new(m, 0.1, 0.2, 2.0)
    assert np.array_equal(out, to_matlab(P.psf_gen(from_matlab(m), 0.1, 0.2, 2.0)))


@pytest.mark.xfail(
    strict=True,
    reason="the MATLAB-name layer binds the positional arguments of a function without an "
    "inputParser (psf_gen_new) as inputParser's optionals, which never take text (F03)",
)
def test_matlab_names_take_the_method_positionally():
    from petakit5d import matlab

    m = bead_psf((21, 30, 31))
    assert np.array_equal(
        matlab.psf_gen_new(m, 0.1, 0.1, 1.5, "median"),
        to_matlab(P.psf_gen(from_matlab(m), 0.1, 0.1, 1.5, "median")),
    )
