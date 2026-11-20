"""Structured illumination reconstruction of one volume (``simRecon`` of ``third_parties/sim3D``)
and its steps: the separation matrix, the OTF resampling, the cylindrical OTF masks, the edge
taper, the Fourier shift, the elliptical apodization.

The data are the raw phases of a deskewed volume, ``(nz * nphases * norientations, ny, nx)``,
the phases (then the orientations) interleaved per z plane as the camera writes them; the OTF is
``(norientations, norders, z, y, x)`` complex (MATLAB's ``(y, x, z, orders, orientations)``; one
orientation may be ``(norders, z, y, x)``). The steps, as MATLAB's:

1. the background subtracted in the data's class (negatives set to 0);
2. the OTF orders interpolated (``interp3`` 'linear', 0 outside) onto the data's frequency grid,
   each orientation's zero order normalized to unit energy;
3. per orientation, each phase tapered (Tukey window in x and y, ``edge_taper``) and made
   periodic (``perdecomp``: the periodic-plus-smooth decomposition), transformed, and the orders
   separated with the pseudo-inverse of the separation matrix;
4. per negative order the pattern's shift vector guessed from ``lattice_period`` and
   ``lattice_angle``, refined once by a masked normalized cross-correlation of the order's overlap
   with the zero order, and its complex weight ``B`` fitted by least squares on the overlap; the
   positive orders are the conjugates;
5. each order Wiener-filtered (``w_ver`` 1: by its own OTF; 0: by the sum of all shifted OTFs,
   translated back; 2: by that sum after the assembly), padded to twice the lateral size, moved
   to its place in frequency space and summed;
6. apodized by an ellipsoid that covers the extended support, transformed back: the real part.

Precision is MATLAB's: on the CPU the volumes are double, but the phase ramps of the Fourier
shifts and the inverse separation matrix are of ``gpu_precision`` (single by default), and the
steps that use them give single results (MATLAB's single with double: computed in double,
rounded), so the result is single; on the GPU every array is of ``gpu_precision``. The result is
MATLAB's within the tolerances of ``docs/parity/sim.md``: single precision sines and cosines,
MATLAB's BLAS products and its linear interpolation are not reproduced bit for bit.

MATLAB quirks kept: every order's overlap is taken with the mask of the third order (the zero
order only for five orders: fewer than three orders fail); the apodization uses the first
orientation's shift vectors and fails for several orientations; ``B_TLS`` (total least squares)
is not used, so only its SVD's failure on NaN or Inf is kept.

Shift vectors (``p_vec``) and frequency steps (``dk``) stay in MATLAB's ``[y x z]`` order inside
(internal vectors compared with MATLAB's, like the deskew transforms); pixel sizes and shapes at
the boundary are ``(z, y, x)``.
"""

from __future__ import annotations

import math
import threading
from collections.abc import Sequence
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "SimError",
    "apodize_ellipse",
    "centered_frequencies",
    "cyl_mask",
    "fourier_shift",
    "inverse_separation_matrix",
    "resample_otf",
    "separation_matrix",
    "sim_recon_volume",
    "tukey_window",
]

_PI = math.pi


class SimError(RuntimeError):
    """An error MATLAB raises here; ``identifier`` is its error id (empty for ``error(msg)``)."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


class SimWarning(UserWarning):
    """A warning MATLAB gives here, with its identifier."""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


# ---------------------------------------------------------------------------------------------
# arrays, precision, FFTs
# ---------------------------------------------------------------------------------------------


def _precision(gpu_precision: str) -> type:
    if gpu_precision == "single":
        return np.float32
    if gpu_precision == "double":
        return np.float64
    raise SimError(
        f"Unsupported class name '{gpu_precision}': gpuPrecision must be 'single' or 'double'.",
        "MATLAB:cast:UnsupportedClass",
    )


def _complex_of(real: Any) -> type:
    return np.complex64 if np.dtype(real) == np.float32 else np.complex128


class _Ops:
    """The array module and FFTs of a path: NumPy with the port's FFTs (MATLAB's driver,
    ``petakit5d.fft``) on the CPU, CuPy with cuFFT on the GPU."""

    def __init__(self, gpu: bool, gpu_precision: str, threads: int, backend: str | None):
        self.gpu = gpu
        self.precision = gpu_precision
        self.p = _precision(gpu_precision)
        self.pc = _complex_of(self.p)
        self.nt = threads
        self.backend = backend
        if gpu:
            from petakit5d.gpu._cupy import import_cupy

            self.xp = import_cupy()
        else:
            self.xp = np

    def fftn(self, x):
        if self.gpu:
            return self.xp.fft.fftn(x)
        from petakit5d.matlabcompat import fft

        return fft.fftn(x, threads=self.nt, backend=self.backend)

    def ifftn(self, x):
        if self.gpu:
            return self.xp.fft.ifftn(x)
        from petakit5d.matlabcompat import fft

        return fft.ifftn(x, threads=self.nt, backend=self.backend)

    def fftshift(self, x):
        return self.xp.fft.fftshift(x)

    def ifftshift(self, x):
        return self.xp.fft.ifftshift(x)

    def host(self, x):
        return self.xp.asnumpy(x) if self.gpu else x


def _round_to(a, dtype):
    """``a`` (computed in double) as ``dtype``: MATLAB's single with double rule"""
    return a if a.dtype == np.dtype(dtype) else a.astype(dtype)


def _wide(dtype) -> type:
    return np.complex128 if np.dtype(dtype).kind == "c" else np.float64


def _is_single(dtype) -> bool:
    return np.dtype(dtype) in (np.dtype(np.float32), np.dtype(np.complex64))


_POOLS: dict[int, Any] = {}
_POOL_LOCK = threading.Lock()


def _pool(nt: int):
    """a thread pool of ``nt`` workers, kept for later calls (threads idle between them)"""
    with _POOL_LOCK:
        ex = _POOLS.get(nt)
        if ex is None:
            from concurrent.futures import ThreadPoolExecutor

            ex = _POOLS[nt] = ThreadPoolExecutor(max_workers=nt, thread_name_prefix="pk5d-sim")
    return ex


def _par(n: int, fn, nt: int) -> None:
    """``fn(lo, hi)`` over pieces of ``range(n)`` on up to ``nt`` threads; each element's
    arithmetic is the same whatever the pieces, so results never depend on ``nt``"""
    if nt <= 1 or n < 2:
        fn(0, n)
        return
    pieces = min(n, 2 * nt)
    edges = [n * i // pieces for i in range(pieces + 1)]
    ex = _pool(min(nt, pieces))
    for f in [ex.submit(fn, edges[i], edges[i + 1]) for i in range(pieces)]:
        f.result()


def _flat(a) -> np.ndarray:
    return a.reshape(-1)


def _binary(a, b, op, nt: int = 1):
    """MATLAB's element-wise ``a op b`` of arrays of one shape (or with a scalar ``b``): single
    with double computed in double and rounded to single; NumPy's arrays in pieces on ``nt``
    threads, CuPy's at once"""
    scalar = not hasattr(b, "shape") or b.shape == ()
    if hasattr(b, "dtype"):
        b_dt = np.dtype(b.dtype)
    else:
        b_dt = np.dtype(np.complex128 if isinstance(b, complex) else np.float64)
    # (a Python number is a double, as in MATLAB)
    mixed = _is_single(a.dtype) != _is_single(b_dt)
    kind = "c" if a.dtype.kind == "c" or b_dt.kind == "c" else "f"
    out_dt = (
        (np.complex64 if kind == "c" else np.float32) if mixed else np.result_type(a.dtype, b_dt)
    )
    if not isinstance(a, np.ndarray):
        if mixed:
            y = b if scalar else b.astype(_wide(b_dt))
            return op(a.astype(_wide(a.dtype)), y).astype(out_dt)
        return op(a, b)
    out = np.empty(a.shape, dtype=out_dt)
    fa = _flat(np.ascontiguousarray(a))
    fb = b if scalar else _flat(np.ascontiguousarray(b))
    fo = _flat(out)

    def run(lo, hi):
        x = fa[lo:hi]
        y = fb if scalar else fb[lo:hi]
        if mixed:
            fo[lo:hi] = op(x.astype(_wide(x.dtype)), y if scalar else y.astype(_wide(b_dt)))
        else:
            fo[lo:hi] = op(x, y)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        _par(fo.size, run, nt)
    return out


def _times(a, b, nt: int = 1):
    """MATLAB's element-wise product (complex or real)"""
    import operator

    return _binary(a, b, operator.mul, nt)


def _plus(a, b, nt: int = 1):
    """MATLAB's sum of two arrays"""
    import operator

    return _binary(a, b, operator.add, nt)


def _div_parts(z, d, nt: int = 1):
    """a complex (or real) array divided by a real scalar or real array as MATLAB divides: each
    part (NumPy's complex division by a real multiplies by its reciprocal)"""
    if z.dtype.kind != "c":
        import operator

        return _binary(z, d, operator.truediv, nt)
    out = z.copy() if isinstance(z, np.ndarray) else z.copy()
    v = out.view(out.real.dtype).reshape(*out.shape, 2)
    scalar = not hasattr(d, "ndim") or d.ndim == 0
    if not isinstance(out, np.ndarray):
        v /= d if scalar else d[..., None]
        return out
    fv = v.reshape(-1, 2)
    fd = d if scalar else _flat(np.ascontiguousarray(d))
    # (single data with a double divisor: divided in double, rounded)
    widen = fv.dtype == np.float32 and not _is_single(np.asarray(d).dtype)

    def run(lo, hi):
        y = fd if scalar else fd[lo:hi, None]
        if widen:
            fv[lo:hi] = fv[lo:hi].astype(np.float64) / y
        else:
            fv[lo:hi] /= y

    with np.errstate(divide="ignore", invalid="ignore"):
        _par(fv.shape[0], run, nt)
    return out


def centered_frequencies(n: float) -> np.ndarray:
    """MATLAB's ``[-ceil((n-1)/2):floor((n-1)/2)]`` (double)"""
    return np.arange(-math.ceil((n - 1) / 2), math.floor((n - 1) / 2) + 1, dtype=np.float64)


def _yxz_of(v: Sequence[float]) -> np.ndarray:
    """a (z, y, x) vector as MATLAB's [y x z]"""
    a = np.asarray(v, dtype=np.float64).ravel()
    if a.size != 3:
        raise SimError("A pixel size needs three values (z, y, x).", "MATLAB:dimagree")
    return a[[1, 2, 0]]


def _colon(a: float, d: float, b: float) -> np.ndarray:
    """MATLAB's ``a : d : b`` for doubles (d > 0): MATLAB's count, the values from both ends,
    the middle one the mean of the ends"""
    from petakit5d.core.dtypes import mround

    if math.isinf(d):
        return np.array([a]) if a <= b else np.zeros(0)
    if d <= 0 or b < a:
        return np.zeros(0)
    tol = 2.0 * np.finfo(np.float64).eps * max(abs(a), abs(b))
    if a == math.floor(a) and d == math.floor(d):
        n = int(math.floor((b - a) / d))
    else:
        n = int(mround(np.float64((b - a) / d)))
        if a + n * d - b > tol:
            n -= 1
    last = a + n * d
    if last - b >= -tol:
        last = b
    out = np.empty(n + 1)
    for k in range(n // 2 + 1):
        out[k] = a + k * d
        out[n - k] = last - k * d
    if n % 2 == 0:
        out[n // 2] = (a + last) / 2
    return out


def _cos(x: float) -> float:
    from petakit5d.matlabcompat import elementary as E

    return float(E.cos(np.float64(x)))


def _sin(x: float) -> float:
    from petakit5d.matlabcompat import elementary as E

    return float(E.sin(np.float64(x)))


def _asin(x: float) -> float:
    if not -1.0 <= x <= 1.0:
        raise SimError(
            "The numerical aperture must not exceed the immersion medium's refractive index.",
            "PetaKit5D:sim:naAboveIndex",
        )
    return math.asin(x)


# ---------------------------------------------------------------------------------------------
# make_forward_separation_matrix, pinv
# ---------------------------------------------------------------------------------------------


def separation_matrix(nphases: int, norders: int, lattice_period: float, phase_step: float):
    """``make_forward_separation_matrix``: the ``(nphases, norders)`` complex matrix of phase
    factors, orders ``[..., -2, -1, 0, +1, +2, ...]`` (a matrix of parameters: ``M[i, j]`` is
    MATLAB's ``M(i + 1, j + 1)``)."""
    nph, nor = int(nphases), int(norders)
    delta_phi = phase_step / lattice_period * 2 * _PI
    m = np.zeros((nph, nor), dtype=np.complex128)
    for ii in range(1, nph + 1):
        counter = 0.0
        phase_sign = 1.0
        for jj in range(1, nor + 1):
            # exp(1i*phase_sign*counter*(ii-1)*delta_phi): the products part by part, as MATLAB
            # multiplies a complex value by reals; exp of (+-0, t) is (cos t, sin t)
            re, im = 0.0 * phase_sign, 1.0 * phase_sign
            for f in (counter, float(ii - 1), delta_phi):
                re, im = re * f, im * f
            e = math.exp(re)
            # (a zero imaginary part makes MATLAB's value real: exp of it has a +0 one)
            m[ii - 1, jj - 1] = complex(e * _cos(im), e * _sin(im)) if im != 0 else complex(e, 0.0)
            counter += jj % 2
            phase_sign = -1.0 + 2.0 * (jj % 2)
    # [sep_matrix(:, end:-2:1), sep_matrix(:, 2:2:end)]
    return np.concatenate([m[:, ::-1][:, ::2], m[:, 1::2]], axis=1)


def inverse_separation_matrix(sep: np.ndarray, gpu_precision: str = "single") -> np.ndarray:
    """``cast(pinv(sep_matrix), gpuPrecision)``: the pseudo-inverse from the SVD (singular values
    up to ``max(size) * eps(norm)`` dropped), in ``gpu_precision``."""
    from petakit5d.core.dtypes import eps

    a = np.asarray(sep, dtype=np.complex128)
    u, s, vh = np.linalg.svd(a, full_matrices=False)
    tol = max(a.shape) * float(eps(float(s[0]))) if s.size else 0.0
    r = int(np.count_nonzero(s > tol))
    pinv = (vh[:r].conj().T / s[:r]) @ u[:, :r].conj().T
    return pinv.astype(_complex_of(_precision(gpu_precision)))


# ---------------------------------------------------------------------------------------------
# resampleOTF, cylMask, tukwin, fourierShift3D, apodizeEllipse
# ---------------------------------------------------------------------------------------------


def _linear_weights(grid: np.ndarray, q: np.ndarray):
    """per query point the left node, the fraction to the right node, and whether it is inside
    the grid (an increasing grid; the last node belongs to the last interval)"""
    n = grid.size
    inside = (q >= grid[0]) & (q <= grid[-1])
    i0 = np.clip(np.searchsorted(grid, q, side="right") - 1, 0, max(n - 2, 0))
    if n == 1:
        return i0, np.zeros(q.shape), inside
    t = (q - grid[i0]) / (grid[i0 + 1] - grid[i0])
    return i0, t, inside


def _interp_axis(v, axis: int, grid: np.ndarray, q: np.ndarray, xp, nt: int = 1):
    """linear interpolation of ``v`` along ``axis`` onto ``q`` (0 outside the grid): ``a w0 +
    b w1`` of the two nodes (NumPy's in pieces of another axis on ``nt`` threads)"""
    i0, t, inside = _linear_weights(grid, q)
    shape = [1] * v.ndim
    shape[axis] = q.size
    w0 = np.where(inside, 1.0 - t, 0.0).reshape(shape)
    w1 = np.where(inside, t, 0.0).reshape(shape)
    if xp is not np:
        a = xp.take(v, xp.asarray(i0), axis=axis) * xp.asarray(w0)
        if grid.size == 1:
            return a
        return a + xp.take(v, xp.asarray(i0 + 1), axis=axis) * xp.asarray(w1)
    out_shape = list(v.shape)
    out_shape[axis] = q.size
    out = np.empty(out_shape, dtype=np.result_type(v.dtype, np.float64))
    split = 0 if axis != 0 else 1

    def run(lo, hi):
        sl = [slice(None)] * v.ndim
        sl[split] = slice(lo, hi)
        part = v[tuple(sl)]
        a = np.take(part, i0, axis=axis) * w0
        if grid.size > 1:
            a = a + np.take(part, i0 + 1, axis=axis) * w1
        out[tuple(sl)] = a

    _par(v.shape[split], run, nt)
    return out


def resample_otf(
    otf,
    pxl_dim_psf: Sequence[float],
    data_shape: Sequence[int],
    pxl_dim_data: Sequence[float],
    *,
    _ops: _Ops | None = None,
):
    """``resampleOTF``: the OTF orders ``(norientations, norders, z, y, x)`` interpolated
    (``interp3`` 'linear', 0 outside) from the OTF's frequency grid onto that of data of
    ``data_shape`` ``(nz, ny, nx)``; pixel sizes ``(z, y, x)`` in um. Complex double on the CPU,
    of the path's precision on the GPU."""
    ops = _ops or _Ops(False, "single", 1, None)
    xp = ops.xp
    o = xp.asarray(otf)
    if o.ndim == 4:
        o = o[None]
    nzp, nyp, nxp = o.shape[2:]
    psf, dat = _yxz_of(pxl_dim_psf), _yxz_of(pxl_dim_data)
    nzd, nyd, nxd = data_shape
    dk_psf = 1.0 / (np.array([nyp, nxp, nzp], dtype=np.float64) * psf)
    dk_dat = 1.0 / (np.array([nyd, nxd, nzd], dtype=np.float64) * dat)
    gy, gx, gz = (centered_frequencies(n) * d for n, d in zip((nyp, nxp, nzp), dk_psf, strict=True))
    qy, qx, qz = (centered_frequencies(n) * d for n, d in zip((nyd, nxd, nzd), dk_dat, strict=True))
    out_dtype = ops.pc if ops.gpu else np.complex128
    out = xp.zeros((o.shape[0], o.shape[1], qz.size, qy.size, qx.size), dtype=out_dtype)
    for jj in range(o.shape[0]):
        for ii in range(o.shape[1]):
            v = o[jj, ii].astype(np.complex128)
            v = _interp_axis(v, 2, gx, qx, xp, ops.nt)
            v = _interp_axis(v, 1, gy, qy, xp, ops.nt)
            out[jj, ii] = _interp_axis(v, 0, gz, qz, xp, ops.nt)
    return out


def cyl_mask(
    na_det: float,
    wvl_em: float,
    na_ext: float,
    wvl_ext: float,
    nimm: float,
    norders: int,
    data_shape: Sequence[int],
    dk_data: Sequence[float],
    islattice: bool,
    *,
    single: bool = False,
) -> np.ndarray:
    """``cylMask``: per order ``(norders, nz, ny, nx)`` the cylinder of the theoretical lateral
    and axial OTF extent; ``dk_data`` in MATLAB's [y x z] order. ``single``: MATLAB's GPU path,
    whose frequency step is single (frequencies and their radius in single)."""
    nz, ny, nx = (int(v) for v in data_shape)
    no = int(norders)
    alpha_det = _asin(na_det / nimm)
    kmag_det = nimm / wvl_em
    alpha_ext = _asin(na_ext / nimm)
    kmag_ext = nimm / wvl_ext
    radial_extent = 2 * na_det / wvl_em
    dt = np.float32 if single else np.float64
    dk = [float(v) for v in dk_data]
    kx = (centered_frequencies(nx) * dk[1]).astype(dt)
    ky = (centered_frequencies(ny) * dk[0]).astype(dt)
    kz = (centered_frequencies(nz) * dk[2]).astype(dt)
    radial = np.sqrt(kx[None, :] * kx[None, :] + ky[:, None] * ky[:, None]) <= radial_extent
    axial_k = np.sqrt(kz * kz)
    mask = np.zeros((no, nz, ny, nx), dtype=bool)
    for ii in range(1, no + 1):
        if islattice:
            ext = kmag_det * (1 - _cos(alpha_det)) + (
                (math.floor(no / 2) - abs((math.ceil(no / 2) - ii) * 0.5)) * kmag_ext
            ) * _sin(alpha_ext)
        else:
            ext = kmag_det * (1 - _cos(alpha_det)) + ((1 - ii % 2) * kmag_ext) * (
                1 - _cos(alpha_ext)
            )
        mask[ii - 1] = radial[None] & (axial_k <= ext)[:, None, None]
    return mask


def tukey_window(shape: Sequence[int], r: float) -> np.ndarray:
    """``tukwin``: Tukey windows in x and y (``r`` the tapered fraction), their minimum, equal
    along z: a ``(1, ny, nx)`` double array that broadcasts to ``shape`` ``(nz, ny, nx)``."""
    from petakit5d.matlabcompat import elementary as E

    _, ny, nx = (int(v) for v in shape)
    h = r / 2
    c = 2 * _PI / r

    def axis_window(n: int) -> np.ndarray:
        u = _colon(0.0, 1.0 / (n - 1) if n > 1 else math.inf, 1.0)
        w = np.ones(u.size)
        i1 = (0 <= u) & (u <= h)
        i2 = (h <= u) & (u <= (1 - h))
        i3 = ((1 - h) <= u) & (u <= 1)
        w[i1] = np.minimum(w[i1], 0.5 * (1 + E.cos(c * (u[i1] - h))))
        w[i2] = np.minimum(w[i2], 1.0)
        w[i3] = np.minimum(w[i3], 0.5 * (1 + E.cos(c * (u[i3] - 1 + h))))
        return w

    wx, wy = axis_window(nx), axis_window(ny)
    if wx.size != nx or wy.size != ny:
        # (the colon's count differs from the size: MATLAB's logical indexing fails)
        raise SimError(
            "The logical indices contain a true value outside of the array bounds.",
            "MATLAB:matrix:indexExceedsDims",
        )
    return np.minimum(wx[None, :], wy[:, None])[None]


def _shift_terms(n: int, shift: float, p, xp):
    """``cast(k, p) / n`` times the shift along one axis (each product in double, rounded to p)"""
    k = centered_frequencies(n).astype(p).astype(np.float64)
    kn = (k / n).astype(p).astype(np.float64)
    return xp.asarray((kn * shift).astype(p))


def fourier_shift(img, shifts: Sequence[float], gpu_precision: str = "single", *, _ops=None):
    """``fourierShift3D``: ``img`` ``(nz, ny, nx)`` shifted by ``shifts`` pixels (MATLAB's [y x z]
    order) by a phase ramp on its centred spectrum; the frequencies and the ramp are of
    ``gpu_precision``, and on the CPU a double image gives a single result for single precision
    (MATLAB's single with double)."""
    ops = _ops or _Ops(False, gpu_precision, 1, None)
    xp = ops.xp
    p = _precision(gpu_precision)
    nz, ny, nx = img.shape
    s = [float(v) for v in shifts]
    # the ramp is applied to the spectrum in the FFT's own order: ifftshift(fftshift(F) .* R) is
    # F .* ifftshift(R), the same products without the two copies
    ty = xp.asarray(np.fft.ifftshift(_shift_terms(ny, s[0], p, np)))
    tx = xp.asarray(np.fft.ifftshift(_shift_terms(nx, s[1], p, np)))
    tz = np.fft.ifftshift(_shift_terms(nz, s[2], p, np))
    tyx = ty[:, None] + tx[None, :]
    f = ops.fftn(img)
    if f.dtype.kind != "c":
        f = f.astype(_complex_of(f.dtype))
    pc = _complex_of(p)
    out_dt = np.complex64 if (_is_single(f.dtype) or pc == np.complex64) else np.complex128
    if ops.gpu:
        # (kyy*s1 + kxx*s2) + kzz*s3 in p; exp(-1i*2*pi*arg) is (cos t, sin t) of
        # t = single(-2 pi arg)
        arg = tyx[None] + xp.asarray(tz)[:, None, None]
        theta = (arg.astype(np.float64) * (-2.0 * _PI)).astype(p)
        ramp = xp.cos(theta) + 1j * xp.sin(theta)
        g = _times(f, ramp)
        del arg, theta, ramp
    else:
        g = np.empty(f.shape, dtype=out_dt)
        # single with double: the product in double, rounded
        wide = f.dtype != pc

        def run(lo, hi):
            for z in range(lo, hi):
                arg = tyx + tz[z]
                th = (arg.astype(np.float64) * (-2.0 * _PI)).astype(p).astype(np.float64)
                ramp = np.empty(arg.shape, pc)
                ramp.real = np.cos(th)
                ramp.imag = np.sin(th)
                if wide:
                    g[z] = f[z].astype(np.complex128) * ramp.astype(np.complex128)
                else:
                    g[z] = f[z] * ramp

        _par(nz, run, ops.nt)
    f = None
    return ops.ifftn(g)


def apodize_ellipse(
    data_k_space,
    dk_data: Sequence[float],
    p_vec_guess: np.ndarray,
    angle,
    na_det: float,
    nimm: float,
    na_ext_max: float,
    wvl_em: float,
    wvl_ext: float,
    islattice: bool,
    *,
    return_apo: bool = True,
    _ops=None,
):
    """``apodizeEllipse``: the spectrum ``(nz, ny, nx)`` times the ellipsoidal apodization (1 at
    the center, 0 at and beyond the extended support), and the apodization (double; None
    without ``return_apo``).
    ``dk_data`` and the shift vectors ``p_vec_guess`` ``(norders, 3, norientations)`` are in
    MATLAB's [y x z] order; only the first orientation's vectors count (MATLAB's ``p(:, 1)``)."""
    ops = _ops or _Ops(False, "single", 1, None)
    xp = ops.xp
    ang = np.asarray(angle, dtype=np.float64).ravel()
    if ang.size != 1:
        # ([cos(angle), -sin(angle), 0; ...] with several angles)
        raise SimError(
            "Dimensions of arrays being concatenated are not consistent.",
            "MATLAB:catenate:dimensionMismatch",
        )
    alpha_det = _asin(na_det / nimm)
    alpha_ext = _asin(na_ext_max / nimm)
    kmag_det = nimm / wvl_em
    kmag_ext = nimm / wvl_ext
    nz, ny, nx = data_k_space.shape
    dk = [float(v) for v in dk_data]
    ky = centered_frequencies(ny) * dk[0]
    kx = centered_frequencies(nx) * dk[1]
    kz = centered_frequencies(nz) * dk[2]
    pv = np.asarray(p_vec_guess, dtype=np.float64)
    pv = pv.reshape(pv.shape[0], -1, order="F")
    a = float(np.max(np.sqrt((pv[:, 0] * dk[0]) ** 2 + (pv[:, 1] * dk[1]) ** 2)))
    a = a + kmag_det * _sin(alpha_det) * 2
    b = kmag_det * _sin(alpha_det) * 2
    if islattice:
        c = kmag_ext * _sin(alpha_ext) * 2 + kmag_det * (1 - _cos(alpha_det))
    else:
        c = kmag_ext * (1 - _cos(alpha_ext)) * 2 + kmag_det * (1 - _cos(alpha_det))
    ca, sa = _cos(ang[0]), _sin(ang[0])
    # the rows of [cos -sin 0; sin cos 0; 0 0 1] * [kx; ky; kz] in index order
    r1 = ca * kx[None, :] + (-sa) * ky[:, None]
    r2 = sa * kx[None, :] + ca * ky[:, None]
    q = r1 * r1 / (a * a) + r2 * r2 / (b * b)
    t3 = kz * kz / (c * c)
    if ops.gpu:
        apo = -(xp.asarray(q)[None] + xp.asarray(t3)[:, None, None]) + 1
        apo[apo < 0] = 0
        return _times(data_k_space, apo), apo
    # per z plane: apo = -(q + t3(z)) + 1, 0 where negative, times the spectrum
    dks = data_k_space
    single = _is_single(dks.dtype)
    out = np.empty(dks.shape, dtype=dks.dtype)
    apo = np.empty(dks.shape) if return_apo else None

    def run(lo, hi):
        for z in range(lo, hi):
            az = -(q + t3[z]) + 1
            az[az < 0] = 0
            if apo is not None:
                apo[z] = az
            # (single data with the double apodization: the product in double, rounded)
            out[z] = (dks[z].astype(_wide(dks.dtype)) if single else dks[z]) * az

    _par(nz, run, ops.nt)
    return out, apo


# ---------------------------------------------------------------------------------------------
# simRecon
# ---------------------------------------------------------------------------------------------


def _matlab_flat(a, mask, xp):
    """``a(mask(:))``: the masked values in MATLAB's linear order (y fastest, then x, then z)"""
    return xp.swapaxes(a, -1, -2)[xp.swapaxes(mask, -1, -2)]


def _sum_vector(v, gpu: bool):
    """MATLAB's sum of a vector (complex: part by part), in MATLAB's order at one thread on
    the CPU; on the GPU CuPy's"""
    if gpu:
        return complex(v.sum()) if v.dtype.kind == "c" else float(v.sum())
    from petakit5d.matlabcompat import reduce

    def one(a) -> float:
        return float(np.asarray(reduce.sum(np.ascontiguousarray(a).ravel())).ravel()[0])

    if v.dtype.kind == "c":
        return complex(one(v.real), one(v.imag))
    return one(v)


def _sum_z(a, gpu: bool):
    """``sum(a, 3)``: each (y, x) summed along z in index order"""
    if gpu:
        return a.sum(axis=0)
    acc = a[0].copy()
    for z in range(1, a.shape[0]):
        acc += a[z]
    return acc


def _abs(a, gpu: bool):
    """MATLAB's abs: of complex single msun's hypotf (the CPU), else the magnitude"""
    if not gpu and a.dtype == np.complex64:
        from petakit5d.kernels import ops

        return ops.complex_abs(a)
    return abs(a)


def _subtract_background(data: np.ndarray, background: float) -> np.ndarray:
    """``data = data - Background; data(data < 0) = 0`` in the data's class"""
    from petakit5d.core import dtypes

    if data.dtype.kind in "iu" and float(background) == math.floor(float(background)):
        # (a whole background: the difference is exact, then saturated at the class's limits)
        info = np.iinfo(data.dtype)
        d = data.astype(np.int64) - int(background)
        return np.clip(d, info.min, info.max).astype(data.dtype)
    if data.dtype.kind in "iub":
        return np.asarray(dtypes.minus(data, float(background)))
    d = _round_to(data.astype(np.float64) - float(background), data.dtype)
    d[d < 0] = 0
    return d


def _mround(v: float) -> int:
    return int(math.floor(abs(v) + 0.5)) * (1 if v >= 0 else -1)


def _shift_slices(shift: Sequence[int], msz: Sequence[int]):
    """source and destination of ``cylmask_shift(sn:tn) = cylmask(so:to)`` (``shift`` MATLAB's
    [y x z], ``msz`` (nz, ny, nx))"""
    sy, sx, sz = shift
    src, dst = [], []
    for n, s in ((msz[0], sz), (msz[1], sy), (msz[2], sx)):
        so, to = max(1, 1 - s), min(n, n - s)
        sn, tn = max(1, s + 1), min(n, n + s)
        src.append(slice(so - 1, max(so - 1, to)))
        dst.append(slice(sn - 1, max(sn - 1, tn)))
    return tuple(src), tuple(dst)


def _pad_lateral(x, pad: Sequence[int], xp, nt: int = 1):
    """``padarray(x, padrange, 'both')`` with zeros (``pad`` MATLAB's [y x z]); NumPy's copied in
    pieces of z on ``nt`` threads"""
    py, px, pz = (int(v) for v in pad)
    nz, ny, nx = x.shape
    out = xp.zeros((nz + 2 * pz, ny + 2 * py, nx + 2 * px), dtype=x.dtype)
    if xp is not np:
        out[pz : pz + nz, py : py + ny, px : px + nx] = x
        return out

    def run(lo, hi):
        out[pz + lo : pz + hi, py : py + ny, px : px + nx] = x[lo:hi]

    _par(nz, run, nt)
    return out


def _volume_to_matlab(v):
    from petakit5d.core.arrays import to_matlab

    return to_matlab(np.asarray(v))


def _array_or_path_in(v):
    """the MATLAB-name layer's conversion of an image array, or a file name (kept)"""
    import os

    if isinstance(v, (str, os.PathLike)) or v is None:
        return v
    from petakit5d.core.arrays import from_matlab

    a = np.asarray(v)
    return a if a.size == 0 else from_matlab(a)


def _array_or_path_out(v):
    import os

    if isinstance(v, (str, os.PathLike)) or v is None:
        return v
    a = np.asarray(v)
    return a if a.size == 0 else _volume_to_matlab(a)


ARRAY_OR_PATH = (_array_or_path_in, _array_or_path_out)


@register(
    "simRecon",
    Param(
        "data",
        None,
        "Raw phases of a deskewed volume (nz * phases * orientations, y, x), phases interleaved.",
        convert="volume",
    ),
    Param(
        "otf",
        None,
        "OTF orders (norientations, norders, z, y, x), complex (sim_psf_to_otf).",
        convert="volume",
    ),
    Param("islattice", "@islogical", "Lattice light sheet data (else 3-beam SIM)."),
    Param("NA_det", "@isnumeric", "Detection NA."),
    Param("NA_ext", "@isnumeric", "Excitation NA."),
    Param("nimm", "@isnumeric", "Refractive index of the immersion medium."),
    Param("wvl_em", "@isnumeric", "Emission wavelength in um."),
    Param("wvl_ext", "@isnumeric", "Excitation wavelength in um."),
    Param("w", "@isnumeric", "Wiener coefficient for regularization."),
    Param(
        "wVer",
        "@isnumeric",
        "Wiener filter: 1 each order's OTF, 0 the sum of the shifted OTFs, 2 that sum after.",
    ),
    Param("apodize", "@islogical", "Apodize the final spectrum."),
    Param("normalize_orientations", "@islogical", "Normalize the total intensity per orientation."),
    Param("perdecomp", "@islogical", "Periodic-plus-smooth decomposition against edge effects."),
    Param("edgeTaper", "@islogical", "Taper the edges (Tukey window) against edge effects."),
    Param("edgeTaperVal", "@isnumeric", "Roll-off of the Tukey window."),
    Param("nphases", "@isnumeric", "Phases."),
    Param("norders", "@isnumeric", "Information orders."),
    Param("norientations", "@isnumeric", "Pattern orientations."),
    Param("lattice_period", "@isnumeric", "Lattice period in um (the coarsest)."),
    Param("lattice_angle", "@isnumeric", "Angle of the pattern modulation per orientation (rad)."),
    Param("phase_step", "@isnumeric", "Phase step in um."),
    Param(
        "pxl_dim_data", "@isnumeric", "Voxel size of the data (z, y, x) in um.", convert="size_yxz"
    ),
    Param(
        "pxl_dim_PSF", "@isnumeric", "Voxel size of the PSF (z, y, x) in um.", convert="size_yxz"
    ),
    Param("Background", "@isnumeric", "Background subtracted from the data.", convert="same"),
    Param("Save16bit", "@islogical", "Not used here (sim_recon_file saves)."),
    Param(
        "gpuPrecision", "@ischar", "Precision of GPU arrays and phase ramps: 'single', 'double'."
    ),
    Param("useGPU", "@islogical", "Run on a GPU if one is usable."),
    Param("num_threads", None, "Threads (None: all available).", python_only=True),
    Param("fft_backend", None, "CPU FFT backend (None: the default).", python_only=True),
    task="C06",
    matlab_out=_volume_to_matlab,
)
def sim_recon_volume(
    data,
    otf,
    *,
    islattice: bool = True,
    na_det: float = 1.0,
    na_ext: float = 0.55,
    nimm: float = 1.33,
    wvl_em: float = 0.605,
    wvl_ext: float = 0.56,
    w: float = 0.005,
    w_ver: float = 1,
    apodize: bool = True,
    normalize_orientations: bool = False,
    perdecomp: bool = True,
    edge_taper: bool = True,
    edge_taper_val: float = 0.1,
    nphases: float = 5,
    norders: float = 5,
    norientations: float = 1,
    lattice_period: float = 1.2021,
    lattice_angle: Any = 1.5707963267948966,
    phase_step: float = 0.232,
    pxl_dim_data: Sequence[float] = (0.16074803849369898, 0.11, 0.11),
    pxl_dim_psf: Sequence[float] = (0.10716535899579933, 0.11, 0.11),
    background: float = 105,
    save_16bit: bool = False,
    gpu_precision: str = "single",
    use_gpu: bool = True,
    num_threads: int | None = None,
    fft_backend: str | None = None,
):
    """A super-resolved SIM volume from the raw phases of a deskewed volume and the OTF
    (``simRecon``): ``(nz, 2 (ny // 2) + ny, 2 (nx // 2) + nx)`` (twice the lateral size), single
    (double with ``gpu_precision='double'``). With ``use_gpu`` it runs on a usable GPU (CuPy),
    else on the CPU; MATLAB's ``useGPU`` is the same (true with a GPU)."""
    from petakit5d.core.threads import resolve_threads

    del save_16bit
    gpu = False
    if use_gpu:
        from petakit5d.gpu.device import usable_gpu_count

        gpu = usable_gpu_count()[0] > 0
    ops = _Ops(gpu, gpu_precision, resolve_threads(num_threads), fft_backend)
    if gpu:
        from petakit5d.gpu.device import select_device
        from petakit5d.gpu.memory import task_scope

        select_device()
        with task_scope():
            return ops.host(_sim_recon(ops, data, otf, locals()))
    return _sim_recon(ops, data, otf, locals())


def _sim_recon(ops: _Ops, data, otf, k: dict[str, Any]) -> Any:
    from petakit5d.sim.masked_ncc import masked_translation_registration
    from petakit5d.sim.perdecomp import perdecomp_3d

    xp, gpu, p, pc = ops.xp, ops.gpu, ops.p, ops.pc
    nph, nor, nori = (float(k[n]) for n in ("nphases", "norders", "norientations"))
    nph_i, nor_i, nori_i = int(nph), int(nor), int(nori)
    pxd, pxp = _yxz_of(k["pxl_dim_data"]), _yxz_of(k["pxl_dim_psf"])
    angles = np.asarray(k["lattice_angle"], dtype=np.float64).ravel()
    o = np.asarray(otf) if not gpu else otf
    if o.ndim == 4:
        o = o[None]
    if o.ndim != 5:
        raise SimError(
            "The OTF needs (norientations, norders, z, y, x) axes.", "MATLAB:badsubscript"
        )
    nzp, nyp, nxp = o.shape[2:]
    dk_psf = 1.0 / (np.array([nyp, nxp, nzp], dtype=np.float64) * pxp)
    if gpu:
        dk_psf = dk_psf.astype(p)

    # background, sizes
    d = np.asarray(data) if not hasattr(data, "__cuda_array_interface__") else ops.host(data)
    if d.ndim == 2:
        d = d[None]
    d = _subtract_background(d, k["background"])
    nimgs, ny, nx = d.shape
    nz_f = nimgs / (nph * nori)
    dk_data = 1.0 / (np.array([ny, nx, nz_f], dtype=np.float64) * pxd)
    if gpu:
        dk_data = dk_data.astype(p)
    if nz_f != math.floor(nz_f):
        raise SimError("Size inputs must be integers.", "MATLAB:NonIntegerInput")
    nz = int(nz_f)
    step = nph_i * nori_i
    if k["normalize_orientations"]:
        d = _normalize_orientations(d, nph_i, nori_i, nz)
    if gpu:
        d = xp.asarray(d)

    # resampled OTFs, each orientation's zero order to unit energy
    o_scaled = resample_otf(o, k["pxl_dim_psf"], (nz, ny, nx), k["pxl_dim_data"], _ops=ops)
    zero = math.ceil(nor / 2) - 1
    pdk = np.prod(dk_psf)
    for jj in range(nori_i):
        o0 = o_scaled[jj, zero]
        energy = _sum_vector(((o0.real * o0.real + o0.imag * o0.imag) * pdk).ravel(), gpu)
        o_scaled[jj] = _div_parts(o_scaled[jj], math.sqrt(energy) if energy >= 0 else math.nan)
        del o0
    if nor_i < 3 and math.floor(nor / 2) >= 1:
        raise SimError(
            f"Index in position 4 exceeds array bounds. Index must not exceed {nor_i}.",
            "MATLAB:badsubscript",
        )
    cylmask = xp.asarray(
        cyl_mask(
            k["na_det"], k["wvl_em"], k["na_ext"], k["wvl_ext"], k["nimm"], nor_i, (nz, ny, nx),
            dk_data, k["islattice"], single=gpu and p == np.float32,
        )
    )  # fmt: skip

    # separated orders
    inv_sep = inverse_separation_matrix(
        separation_matrix(nph_i, nor_i, k["lattice_period"], k["phase_step"]), ops.precision
    )
    store = pc if gpu else np.complex128
    dk_sep = xp.zeros((nori_i, nor_i, nz, ny, nx), dtype=store)
    pdd = np.prod(dk_data)
    window = (
        xp.asarray(tukey_window((nz, ny, nx), k["edge_taper_val"])) if k["edge_taper"] else None
    )
    m = xp.asarray(inv_sep if gpu else inv_sep.astype(np.complex128))
    for jj in range(nori_i):
        dks = []
        for ii in range(nph_i):
            dr = d[ii + jj * nph_i :: step].astype(p if gpu else np.float64)
            if window is not None:
                dr = _round_to(dr * window, dr.dtype)
            if k["perdecomp"]:
                dr = perdecomp_3d(dr, use_gpu=gpu, gpu_precision=ops.precision, _ops=ops)
                dr = dr.astype(p if gpu else np.float64)
            dk = ops.fftshift(ops.ifftn(ops.ifftshift(dr)))
            dks.append(_div_parts(dk if dk.dtype == store else dk.astype(store), pdd, ops.nt))
            del dr, dk
        # reshape(Dk, [], nphases) * inv_sep_matrix.': order j is sum_k inv(j, k) Dk_k; on the
        # CPU double data with the single matrix give single values (computed in double)
        for j in range(nor_i):
            if gpu:
                acc = dks[0] * m[j, 0]
                for kk in range(1, nph_i):
                    acc = acc + dks[kk] * m[j, kk]
                dk_sep[jj, j] = acc
                del acc
                continue
            out = _flat(dk_sep[jj, j])
            flats = [_flat(a) for a in dks]
            row = [m[j, kk] for kk in range(nph_i)]
            single = p == np.float32

            def run(lo, hi, out=out, flats=flats, row=row, single=single):
                acc = flats[0][lo:hi] * row[0]
                for kk in range(1, nph_i):
                    acc = acc + flats[kk][lo:hi] * row[kk]
                out[lo:hi] = acc.astype(np.complex64) if single else acc

            _par(out.size, run, ops.nt)
        del dks
    del d

    # shift vectors and complex weights
    nneg = math.floor(nor / 2)
    p_vec = np.zeros((nor_i, 3, nori_i))
    b_pls = np.zeros((nor_i, nori_i), dtype=np.complex128)
    msz = (nz, ny, nx)
    ks = (1 / k["lattice_period"]) / np.asarray(dk_data, dtype=np.float64)
    if gpu:
        ks = ks.astype(p).astype(np.float64)
    for jj in range(nori_i):
        a_j = float(angles[jj]) if jj < angles.size else math.nan
        unit = np.array([_sin(a_j), _cos(a_j), 0.0])
        for kk in range(1, nneg + 1):
            order = kk - math.ceil(nor / 2)
            pv = (order * ks) * unit
            if gpu:
                pv = pv.astype(p).astype(np.float64)
            p_vec[kk - 1, :, jj] = pv
            overlap = d0om = dmo0 = None
            for qq in (1, 2):
                pv = p_vec[kk - 1, :, jj]
                shift_dk = fourier_shift(dk_sep[jj, kk - 1], pv, ops.precision, _ops=ops)
                shift_om = fourier_shift(o_scaled[jj, kk - 1], pv, ops.precision, _ops=ops)
                # the order's mask moved by the rounded shift (0 where it enters)
                src, dst = _shift_slices([_mround(v) for v in pv], msz)
                cm_shift = xp.zeros(msz, dtype=bool)
                cm_shift[dst] = cylmask[kk - 1][src]
                overlap = cylmask[2] & cm_shift
                del cm_shift
                dmo0 = _times(shift_dk, o_scaled[jj, zero])
                d0om = _times(dk_sep[jj, zero], shift_om)
                del shift_dk, shift_om
                if qq == 1:
                    m2 = overlap.max(axis=0)
                    t = masked_translation_registration(
                        _abs(_sum_z(d0om, gpu), gpu),
                        _abs(_sum_z(dmo0, gpu), gpu),
                        m2,
                        m2,
                        0.5,
                        use_gpu=gpu,
                        gpu_precision=ops.precision,
                    )[0]
                    p_vec[kk - 1, :, jj] = p_vec[kk - 1, :, jj] - np.array(
                        [float(t[1]), float(t[0]), 0.0]
                    )
            ai = _matlab_flat(d0om, overlap, xp)
            bi = _matlab_flat(dmo0, overlap, xp)
            if not (bool(xp.all(xp.isfinite(ai))) and bool(xp.all(xp.isfinite(bi)))):
                # (tls's SVD of [ai bi]; its result is not used)
                raise SimError(
                    "Input to SVD must not contain NaN or Inf.", "MATLAB:svd:matrixWithNaNInf"
                )
            num = _sum_vector(xp.conj(ai) * bi, gpu)
            den = _sum_vector(_abs(ai, gpu) ** 2, gpu)
            with np.errstate(divide="ignore", invalid="ignore"):
                b = complex(np.complex128(num).real / np.float64(den),
                            np.complex128(num).imag / np.float64(den))  # fmt: skip
            b_pls[kk - 1, jj] = complex(np.complex64(b)) if ai.dtype == np.complex64 else b
            del ai, bi, d0om, dmo0, overlap
        c = math.ceil(nor / 2) - 1
        p_vec[c, :, jj] = 0
        b_pls[c, jj] = 1
        if nneg:
            p_vec[c + 1 :, :, jj] = -p_vec[nneg - 1 :: -1, :, jj]
            b_pls[c + 1 :, jj] = np.conj(b_pls[nneg - 1 :: -1, jj])
    del cylmask
    return _assemble(ops, dk_sep, o_scaled, b_pls, p_vec, msz, dk_data, k)


def _normalize_orientations(d: np.ndarray, nph: int, nori: int, nz: int) -> np.ndarray:
    """each orientation's planes divided by its total intensity relative to the largest
    (``normalize_orientations``), in the data's class"""
    from petakit5d.core import dtypes
    from petakit5d.matlabcompat import reduce

    idx = [
        [jj * nph + kk * nph * nori + q for kk in range(nz) for q in range(nph)]
        for jj in range(nori)
    ]
    tot = [reduce.sum(d[i]) for i in idx]
    big = max(tot)
    out = d.copy()
    for jj in range(nori):
        if d.dtype.kind in "iub":
            out[idx[jj]] = dtypes.rdivide(d[idx[jj]], float(tot[jj]) / float(big))
        else:
            ratio = np.asarray(tot[jj] / big, dtype=d.dtype)
            out[idx[jj]] = _round_to(d[idx[jj]].astype(np.float64) / float(ratio), d.dtype)
    return out


def _assemble(ops: _Ops, dk_sep, o_scaled, b_pls, p_vec, shape, dk_data, k):
    """the Wiener filter, the shifts into place, the sum, the apodization, the inverse
    transform"""
    xp, gpu, p = ops.xp, ops.gpu, ops.p
    nz, ny, nx = shape
    nori, nor = dk_sep.shape[:2]
    w2 = k["w"] ** 2
    size = np.array([ny, nx, nz])
    padrange = np.floor((np.array([2, 2, 1]) * size - size) / 2).astype(np.int64)
    w_ver = float(k["w_ver"])

    def b_of(ii, jj):
        b = b_pls[ii, jj]
        return xp.asarray(np.complex64(b) if gpu and p == np.float32 else np.complex128(b))

    if w_ver == 2:
        # (MATLAB clears denom after the filtering and fails at its use in the assembly: the
        # error, without the work before it)
        raise SimError("Reference to a cleared variable denom.", "MATLAB:refClearedVar")
    denom = None
    if w_ver in (0, 2):
        for qq in range(nori):
            for kk in range(nor):
                big = _pad_lateral(o_scaled[qq, kk] * b_of(kk, qq), padrange, xp)
                sh = fourier_shift(big, p_vec[kk, :, qq], ops.precision, _ops=ops)
                a2 = _abs(sh, gpu)
                a2 = a2 * a2
                denom = a2.astype(np.float64) if denom is None and not gpu else (
                    a2 if denom is None else _plus(denom, a2)
                )  # fmt: skip
                del big, sh, a2
    elif w_ver != 1:
        # (MATLAB leaves Dk_sep_scaled 3-D for other values and fails at its 4-D index)
        raise SimError("Index in position 4 exceeds array bounds.", "MATLAB:badsubscript")

    acc = None
    for jj in range(nori):
        for ii in range(nor):
            if w_ver == 1 and not gpu:
                # Dk .* conj(O .* B) ./ (abs(O .* B).^2 + w^2), zero in the padding
                scaled = _pad_lateral(
                    _wiener_own_otf(dk_sep[jj, ii], o_scaled[jj, ii], b_of(ii, jj), w2, ops.nt),
                    padrange,
                    xp,
                    ops.nt,
                )
                sh = fourier_shift(scaled, p_vec[ii, :, jj], ops.precision, _ops=ops)
                del scaled
                if sh.dtype.kind != "c":
                    sh = sh.astype(_complex_of(sh.dtype))
                # (the first sum adds MATLAB's double zeros: a -0 part becomes +0)
                acc = _plus(sh, 0.0, ops.nt) if acc is None else _plus(acc, sh, ops.nt)
                del sh
                continue
            num = dk_sep[jj, ii] * xp.conj(o_scaled[jj, ii] * b_of(ii, jj))
            if w_ver == 1:
                # numerator ./ (abs(O .* B).^2 + w^2), zero in the padding
                ob = _abs(o_scaled[jj, ii] * b_of(ii, jj), gpu)
                scaled = _pad_lateral(_div_parts(num, ob * ob + w2), padrange, xp)
                del ob
            elif w_ver == 0:
                from petakit5d.sim.translate import imtranslate

                tr = -np.array([p_vec[ii, 1, jj], p_vec[ii, 0, jj], p_vec[ii, 2, jj]])
                shifted = imtranslate(denom, tr, _xp=xp)
                scaled = _div_parts(_pad_lateral(num, padrange, xp), shifted + w2)
                del shifted
            else:
                scaled = _pad_lateral(num, padrange, xp)
            del num
            sh = fourier_shift(scaled, p_vec[ii, :, jj], ops.precision, _ops=ops)
            del scaled
            if sh.dtype.kind != "c":
                sh = sh.astype(_complex_of(sh.dtype))
            if acc is None:
                acc = sh if gpu else _plus(sh, 0.0, ops.nt)
            else:
                acc = _plus(acc, sh, ops.nt)
            del sh
    if k["apodize"]:
        acc, _ = apodize_ellipse(
            acc, dk_data, p_vec, k["lattice_angle"], k["na_det"], k["nimm"], k["na_ext"],
            k["wvl_em"], k["wvl_ext"], k["islattice"], _ops=ops, return_apo=False,
        )  # fmt: skip
    # real(fftshift(fftn(ifftshift(.)))) * prod(dk_data): the real part and the scaling before
    # the (element-moving) fftshift
    r = ops.fftn(ops.ifftshift(acc))
    del acc
    re = r.real if r.dtype.kind == "c" else r
    pdd = float(np.prod(dk_data))
    if gpu:
        return ops.fftshift(re * re.dtype.type(pdd))
    return np.fft.fftshift(_times(np.ascontiguousarray(re), pdd, ops.nt))


def _wiener_own_otf(dk, o, b, w2: float, nt: int) -> np.ndarray:
    """``Dk .* conj(O .* B) ./ (abs(O .* B).^2 + w^2)`` of complex double arrays (the parts
    divided), in pieces on ``nt`` threads"""
    out = np.empty(dk.shape, np.complex128)
    fd, fo, fout = _flat(dk), _flat(o), _flat(out)

    def run(lo, hi):
        ob = fo[lo:hi] * b
        num = fd[lo:hi] * np.conj(ob)
        a = np.abs(ob)
        den = a * a + w2
        v = num.view(np.float64).reshape(-1, 2)
        v /= den[:, None]
        fout[lo:hi] = num

    with np.errstate(divide="ignore", invalid="ignore"):
        _par(fout.size, run, nt)
    return out
