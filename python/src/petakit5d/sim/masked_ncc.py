"""Masked normalized cross-correlation and the translation it finds (``normxcorr2_masked`` and
``MaskedTranslationRegistration2D_fit`` in ``third_parties/sim3D``).

A new implementation from the published method (D. Padfield, "Masked object registration in the
Fourier domain", IEEE Trans. Image Process. 21(5), 2012). The MATLAB files (another author's)
are not read; the details marked "measured" were found from their outputs on chosen images. For
a fixed image ``f`` with mask ``mf`` and a moving image ``g`` with mask ``mg``, every relative
offset of the two has the normalized cross-correlation of the pixels where both masks are set,

    C = (S_fg - S_f S_g / N) / sqrt((S_ff - S_f^2 / N) (S_gg - S_g^2 / N)),

where ``N`` is the number of those pixels and the ``S`` are sums over them, all computed for
every offset at once as correlations by FFT (the moving image and its mask turned by 180
degrees). ``C`` and ``N`` have the size of the full correlation, ``(Mf + Mg - 1, Nf + Ng - 1)``.

The registration (measured):

- everything in ``gpu_precision``, on the CPU too; ``N`` (rounded, at least double's ``eps``)
  is of that precision;
- offsets with an overlap below ``overlap_ratio`` times the largest overlap get ``C = 0``;
- the peak is the first maximum of ``C`` in MATLAB's order, refined by the centroid of ``C`` over
  the 11 x 11 window around it (a window that leaves ``C`` fails, as MATLAB's indexing does);
- the translation is ``[x y]``, the moving image's offset from the fixed one: ``moving``'s size
  minus the peak's 1-based subscripts minus the centroid's offsets.
"""

from __future__ import annotations

import numpy as np

__all__ = ["CENTROID_HALF_WIDTH", "masked_translation_registration", "normxcorr2_masked"]

# half width of the window whose centroid refines the peak (measured: 11 x 11)
CENTROID_HALF_WIDTH = 5


def _ncc(f, g, mf, mg, xp, dt):
    """C and the overlap N (rounded, in ``dt``) of images and 0/1 masks already in ``dt``"""
    f = f * mf
    g = g * mg
    gr = g[::-1, ::-1]
    mgr = mg[::-1, ::-1]
    shape = (f.shape[0] + g.shape[0] - 1, f.shape[1] + g.shape[1] - 1)

    def fft(a):
        return xp.fft.fft2(a, s=shape)

    def corr(a, b):
        return xp.fft.ifft2(a * b).real.astype(dt)

    ff, mff = fft(f), fft(mf)
    gg, mgg = fft(gr), fft(mgr)
    # the overlap counts, rounded, at least eps (double's, in the precision: measured)
    n = xp.maximum(xp.round(corr(mff, mgg)), dt(np.finfo(np.float64).eps))
    sf = corr(ff, mgg)
    sg = corr(mff, gg)
    num = corr(ff, gg) - sf * sg / n
    vf = xp.maximum(corr(fft(f * f), mgg) - sf * sf / n, dt(0))
    vg = xp.maximum(corr(mff, fft(gr * gr)) - sg * sg / n, dt(0))
    den = xp.sqrt(vf * vg)
    big = float(xp.max(xp.abs(den))) if den.size else 0.0
    tol = dt(1000) * dt(np.spacing(dt(big)))
    c = xp.zeros(shape, dtype=dt)
    ok = den > tol
    c[ok] = num[ok] / den[ok]
    return xp.clip(c, dt(-1), dt(1)), n


def _inputs(fixed, moving, fixed_mask, moving_mask, use_gpu: bool, gpu_precision: str):
    from petakit5d.sim.recon import _precision

    dt = _precision(gpu_precision)
    if use_gpu:
        from petakit5d.gpu._cupy import import_cupy

        xp = import_cupy()
    else:
        xp = np
    f = xp.asarray(fixed).astype(dt)
    g = xp.asarray(moving).astype(dt)
    mf = (xp.asarray(fixed_mask) != 0).astype(dt)
    mg = (xp.asarray(moving_mask) != 0).astype(dt)
    return f, g, mf, mg, xp, dt


def normxcorr2_masked(
    fixed,
    moving,
    fixed_mask,
    moving_mask,
    use_gpu: bool = False,
    gpu_precision: str = "double",
):
    """``[C, numberOfOverlapMaskedPixels] = normxcorr2_masked(fixedImage, movingImage,
    fixedMask, movingMask, useGPU, gpuPrecision)``: the masked NCC of every offset and the
    number of overlapping masked pixels, ``(Mf + Mg - 1, Nf + Ng - 1)`` (``(y, x)`` images)."""
    f, g, mf, mg, xp, dt = _inputs(fixed, moving, fixed_mask, moving_mask, use_gpu, gpu_precision)
    return _ncc(f, g, mf, mg, xp, dt)


def _window_error(index: int, size: int, position: int):
    from petakit5d.sim.recon import SimError

    if index < 1:
        return SimError(
            f"Index in position {position} is invalid. Array indices must be positive integers "
            "or logical values.",
            "MATLAB:badsubscript",
        )
    return SimError(
        f"Index in position {position} exceeds array bounds. Index must not exceed {size}.",
        "MATLAB:badsubscript",
    )


def masked_translation_registration(
    fixed,
    moving,
    fixed_mask,
    moving_mask,
    overlap_ratio: float = 3 / 10,
    use_gpu: bool = False,
    gpu_precision: str = "single",
):
    """``[transform, maxC, C, numberOfOverlapMaskedPixels] =
    MaskedTranslationRegistration2D_fit(fixedImage, movingImage, fixedMask, movingMask,
    overlapRatio, useGPU, gpuPrecision)``: the translation ``[x, y]`` (MATLAB's order) of
    ``moving`` against ``fixed`` under their masks, refined to a fraction of a pixel, the
    correlation at the peak, the correlation image (0 where the overlap is below
    ``overlap_ratio`` of the largest) and the overlap counts. Arrays are ``(y, x)`` images."""
    f, g, mf, mg, xp, dt = _inputs(fixed, moving, fixed_mask, moving_mask, use_gpu, gpu_precision)
    c, n = _ncc(f, g, mf, mg, xp, dt)
    nmax = float(xp.max(n)) if n.size else 0.0
    c[n < float(overlap_ratio) * nmax] = 0
    ch = np.asarray(xp.asnumpy(c) if xp is not np else c)
    # the first maximum in MATLAB's order (y fastest)
    k = int(np.argmax(np.swapaxes(ch, 0, 1).ravel()))
    xpk, ypk = divmod(k, ch.shape[0])
    w = CENTROID_HALF_WIDTH
    for pos, (p, size) in enumerate(((ypk, ch.shape[0]), (xpk, ch.shape[1])), start=1):
        if p - w < 0:
            raise _window_error(p - w + 1, size, pos)
        if p + w >= size:
            raise _window_error(p + w + 1, size, pos)
    s = ch[ypk - w : ypk + w + 1, xpk - w : xpk + w + 1]
    # the centroid's offsets from the peak, summed in MATLAB's order
    r = np.arange(-w, w + 1, dtype=dt)
    sv = np.swapaxes(s, 0, 1).ravel()
    yr = np.tile(r, 2 * w + 1)
    xr = np.repeat(r, 2 * w + 1)
    tot = _msum(sv)
    cy = dt(_msum(sv * yr) / tot)
    cx = dt(_msum(sv * xr) / tot)
    ty = dt(dt(moving.shape[0] - (ypk + 1)) - cy)
    tx = dt(dt(moving.shape[1] - (xpk + 1)) - cx)
    return np.array([tx, ty], dtype=dt), dt(ch[ypk, xpk]), c, n


def _msum(v: np.ndarray):
    """MATLAB's sum of a vector in its order, in the vector's class"""
    from petakit5d.matlabcompat import reduce

    return np.asarray(reduce.sum(np.ascontiguousarray(v))).ravel()[0]
