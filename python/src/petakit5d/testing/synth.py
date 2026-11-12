"""Synthetic test data, bit-stable across NumPy versions, platforms and CPUs.

The generators of the MATLAB suite (``tests/+pktest``: beads, texture, deconPSF,
writeStitchTiles, centroid(s), the tests' randomImage) with random numbers of our own:

- ``CounterRNG`` is Philox4x32-10 (Salmon et al., SC'11), a counter-based generator. Value i
  of a stream is a pure function of (seed, stream, i), so arrays have the same content however
  they are generated, in one call or in pieces.
- Floating-point results use only operations IEEE 754 rounds correctly (+, -, *, /, sqrt,
  floor, frexp, ldexp, conversions), in a fixed order, with exp, log, sin and cos of our own
  built from them. NumPy's and the C library's exp/log/sin/cos differ between platforms and
  SIMD paths in the last bit; NumPy's sums change order with its version.

Arrays are the port's (docs/conventions.md): volumes ``(z, y, x)`` C-ordered, images
``(y, x)``; sizes, positions, sigmas and offsets in ``(z, y, x)``, positions and offsets
0-based. Random arrays are filled in MATLAB's linear order (y fastest, then x, then z), so
``core.arrays.to_matlab`` of any result is, bit for bit, the array these generators gave in
MATLAB's ``(y, x, z)`` form before the convention change: stored input hashes stay valid.
MATLAB never draws these numbers itself: oracle cases hand it the arrays (``to_matlab``) or the
files written here, so both sides read the same data.

Files: tiff as cpp-tiff writes a ``(z, y, x)`` array (unchanged on disk); zarr native
(``(z, y, x)``, C order) by default, or in MATLAB's layout (``(y, x, z)``, F order) on request,
for MATLAB and for the tests of MATLAB-layout reading. ``write_paired_inputs`` writes one
volume for both sides of a benchmark (each reads its own layout).
"""

from __future__ import annotations

import hashlib
import math
import os
from collections.abc import Sequence
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np

__all__ = [
    "CounterRNG",
    "array_hash",
    "beads",
    "centroid",
    "centroids",
    "cut_tiles",
    "decon_psf",
    "gaussian_filter",
    "random_image",
    "stable_exp",
    "stable_log",
    "stable_sincos_2pi",
    "PairedInput",
    "stitch_tile_name",
    "texture",
    "write_paired_inputs",
    "write_stitch_tiles",
    "write_tiff",
    "write_tiff_matlab",
    "write_zarr",
    "write_zarr_matlab",
]

# ------------------------------------------------------------------------------------------
# Philox4x32-10
# ------------------------------------------------------------------------------------------

_M32 = 0xFFFFFFFF
_PHILOX_M0 = 0xD2511F53
_PHILOX_M1 = 0xCD9E8D57
_PHILOX_W0 = 0x9E3779B9
_PHILOX_W1 = 0xBB67AE85


def philox4x32(counter: np.ndarray, key: tuple[int, int], rounds: int = 10) -> np.ndarray:
    """Philox4x32 of counters (n x 4 uint32 words) under a key (2 uint32 words): n x 4 uint32"""
    c = np.asarray(counter, dtype=np.uint64)
    if c.ndim != 2 or c.shape[1] != 4:
        raise ValueError("counter must be n x 4")
    c0, c1, c2, c3 = (c[:, i].copy() for i in range(4))
    k0, k1 = int(key[0]) & _M32, int(key[1]) & _M32
    m0, m1, mask, s32 = np.uint64(_PHILOX_M0), np.uint64(_PHILOX_M1), np.uint64(_M32), np.uint64(32)
    for r in range(rounds):
        if r:
            k0 = (k0 + _PHILOX_W0) & _M32
            k1 = (k1 + _PHILOX_W1) & _M32
        # 32 x 32 -> 64 bit products are exact in uint64
        p0 = m0 * c0
        p1 = m1 * c2
        hi0, lo0 = p0 >> s32, p0 & mask
        hi1, lo1 = p1 >> s32, p1 & mask
        c0, c1, c2, c3 = hi1 ^ c1 ^ np.uint64(k0), lo1, hi0 ^ c3 ^ np.uint64(k1), lo0
    return np.stack([c0, c1, c2, c3], axis=1).astype(np.uint32)


def _stream_words(stream: int | str) -> tuple[int, int]:
    if isinstance(stream, str):
        v = int.from_bytes(hashlib.sha256(stream.encode()).digest()[:8], "little")
    else:
        v = int(stream)
        if not 0 <= v < 2**64:
            raise ValueError("stream must be in [0, 2^64)")
    return v & _M32, v >> 32


class CounterRNG:
    """Random numbers from Philox4x32-10, bit-stable everywhere.

    seed: the key (0 <= seed < 2^64). stream: an int or a name (hashed), the upper half of the
    counter, so streams of one seed never overlap. Each block of the stream (one counter) gives
    4 uint32 words: 4 ``uint32``/``uniform(float32)``/``integers`` values, or 2
    ``uniform(float64)``/``normal`` values. A call uses whole blocks from the current position
    (``position``, in blocks) and moves it on; ``at(block)`` returns a copy at another
    position. Shapes are the port's (``(z, y, x)``, ``(y, x)``, leading extra axes); arrays are
    C-ordered and filled in MATLAB's linear order (y fastest, then x, then z, then the leading
    axes): value k of the call sits at MATLAB linear index k, so ``to_matlab`` of the result is
    the column-major fill of MATLAB's size.
    """

    def __init__(self, seed: int = 0, stream: int | str = 0, position: int = 0):
        seed = int(seed)
        if not 0 <= seed < 2**64:
            raise ValueError("seed must be in [0, 2^64)")
        self.seed = seed
        self.stream = stream
        self._key = (seed & _M32, seed >> 32)
        self._s = _stream_words(stream)
        self.position = int(position)

    def at(self, position: int) -> CounterRNG:
        return CounterRNG(self.seed, self.stream, position)

    def spawn(self, name: str | int) -> CounterRNG:
        """an independent generator of the same seed (a stream named after this one's and name)"""
        return CounterRNG(self.seed, f"{self.stream}/{name}")

    def blocks(self, n: int) -> np.ndarray:
        """the next n blocks: n x 4 uint32"""
        n = int(n)
        idx = np.arange(self.position, self.position + n, dtype=np.uint64)
        ctr = np.empty((n, 4), dtype=np.uint64)
        ctr[:, 0] = idx & np.uint64(_M32)
        ctr[:, 1] = idx >> np.uint64(32)
        ctr[:, 2] = self._s[0]
        ctr[:, 3] = self._s[1]
        self.position += n
        return philox4x32(ctr, self._key)

    def _take(self, shape: int | Sequence[int], per_block: int) -> tuple[tuple[int, ...], int, int]:
        shape = (int(shape),) if np.isscalar(shape) else tuple(int(s) for s in shape)
        if any(s < 0 for s in shape):
            raise ValueError("negative size")
        n = math.prod(shape)
        return shape, n, -(-n // per_block)

    @staticmethod
    def _shape(flat: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
        if len(shape) < 2:
            return np.ascontiguousarray(flat.reshape(shape))
        # MATLAB's linear order is the C order of the (..., x, y) view
        swapped = flat.reshape((*shape[:-2], shape[-1], shape[-2]))
        return np.ascontiguousarray(np.swapaxes(swapped, -1, -2))

    def uint32(self, shape: int | Sequence[int]) -> np.ndarray:
        shape, n, nb = self._take(shape, 4)
        return self._shape(self.blocks(nb).reshape(-1)[:n], shape)

    def uniform(self, shape: int | Sequence[int], dtype: Any = np.float64) -> np.ndarray:
        """uniform in [0, 1): float64 with 53 random bits (2 words), float32 with 24 (1 word)"""
        dtype = np.dtype(dtype)
        if dtype == np.float32:
            shape, n, nb = self._take(shape, 4)
            w = self.blocks(nb).reshape(-1)[:n]
            u = (w >> np.uint32(8)).astype(np.float32) * np.float32(2.0**-24)
            return self._shape(u, shape)
        if dtype != np.float64:
            raise TypeError("uniform: float32 or float64")
        shape, n, nb = self._take(shape, 2)
        return self._shape(_u53(self.blocks(nb).reshape(-1, 2))[:n], shape)

    def normal(self, shape: int | Sequence[int], dtype: Any = np.float64) -> np.ndarray:
        """standard normal (Box-Muller on two 53-bit uniforms; float32: the float64 rounded)"""
        dtype = np.dtype(dtype)
        if dtype not in (np.float32, np.float64):
            raise TypeError("normal: float32 or float64")
        shape, n, nb = self._take(shape, 2)
        b = self.blocks(nb)
        u1 = 1.0 - _u53(b[:, :2])  # (0, 1]
        u2 = _u53(b[:, 2:])
        r = np.sqrt(-2.0 * stable_log(u1))
        s, c = stable_sincos_2pi(u2)
        z = np.stack([r * c, r * s], axis=1).reshape(-1)[:n]
        return self._shape(z.astype(dtype), shape)

    def integers(
        self, low: int, high: int, shape: int | Sequence[int], dtype: Any = np.int64
    ) -> np.ndarray:
        """integers in [low, high] (both included, as MATLAB's randi), high - low < 2^32"""
        low, high = int(low), int(high)
        span = high - low + 1
        if not 0 < span <= 2**32:
            raise ValueError("integers: need 0 <= high - low < 2^32")
        shape, n, nb = self._take(shape, 4)
        w = self.blocks(nb).reshape(-1)[:n].astype(np.uint64)
        # the high word of word * span: exact in uint64 (span <= 2^32)
        if span == 2**32:
            v = w
        else:
            v = (w * np.uint64(span)) >> np.uint64(32)
        if not -(2**62) <= low < 2**62:
            raise ValueError("integers: |low| < 2^62")
        out = v.astype(np.int64) + low
        return self._shape(out.astype(np.dtype(dtype), copy=False), shape)


def _u53(words: np.ndarray) -> np.ndarray:
    """53-bit uniforms in [0, 1) from pairs of words (n x 2 uint32), as MT19937's genrand_res53"""
    a = (words[:, 0] >> np.uint32(5)).astype(np.float64)
    b = (words[:, 1] >> np.uint32(6)).astype(np.float64)
    return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0)


# ------------------------------------------------------------------------------------------
# elementary functions from correctly rounded operations (float64, elementwise)
# ------------------------------------------------------------------------------------------

_LN2_HI = 6.93147180369123816490e-01  # ln 2 with 21 trailing zero bits: k * _LN2_HI is exact
_LN2_LO = 1.90821492927058770002e-10
_INV_LN2 = 1.44269504088896338700e00
# 1 / k! for k = 13 .. 0 (Horner order)
_EXP_C = tuple(1.0 / math.factorial(k) for k in range(13, -1, -1))
# 1 / (2k + 1) for k = 11 .. 0
_LOG_C = tuple(1.0 / (2 * k + 1) for k in range(11, -1, -1))
_SQRT_HALF = 0.70710678118654752440
_HALF_PI = 1.57079632679489661923


def _horner(coefs: Sequence[float], x: np.ndarray) -> np.ndarray:
    p = np.full_like(x, coefs[0])
    for c in coefs[1:]:
        p = p * x + c
    return p


def stable_exp(x: Any) -> np.ndarray:
    """exp(x) in float64, the same bits everywhere (within about 1 ulp of exp).

    Results below 2^-1000 are 0, so subnormal handling (flush-to-zero modes) never matters.
    """
    x = np.asarray(x, dtype=np.float64)
    big = x > 709.0
    tiny = x < -690.0
    xc = np.clip(x, -690.0, 709.0)
    k = np.floor(xc * _INV_LN2 + 0.5)
    r = (xc - k * _LN2_HI) - k * _LN2_LO
    p = np.ldexp(_horner(_EXP_C, r), k.astype(np.int32))
    p = np.where(tiny, 0.0, p)
    p = np.where(big, np.inf, p)
    return np.where(np.isnan(x), np.nan, p)


def stable_log(x: Any) -> np.ndarray:
    """natural log in float64 for finite x > 0, the same bits everywhere (about 1 ulp)"""
    x = np.asarray(x, dtype=np.float64)
    if np.any(~(x > 0)) or np.any(~np.isfinite(x)):
        raise ValueError("stable_log: finite positive input only")
    m, e = np.frexp(x)  # x = m 2^e, m in [0.5, 1)
    low = m < _SQRT_HALF
    m = np.where(low, m * 2.0, m)
    e = (e - low).astype(np.float64)
    s = (m - 1.0) / (m + 1.0)
    s2 = s * s
    lm = 2.0 * s * _horner(_LOG_C, s2)
    return e * _LN2_HI + (e * _LN2_LO + lm)


def _sin_cos_small(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """sin and cos for 0 <= a <= pi/4 (Taylor series to degree 17/18)"""
    a2 = a * a
    sin_c = tuple((-1.0) ** k / math.factorial(2 * k + 1) for k in range(8, -1, -1))
    cos_c = tuple((-1.0) ** k / math.factorial(2 * k) for k in range(9, -1, -1))
    return a * _horner(sin_c, a2), _horner(cos_c, a2)


def stable_sincos_2pi(u: Any) -> tuple[np.ndarray, np.ndarray]:
    """sin(2 pi u) and cos(2 pi u) for u in [0, 1) with at most 53 significant bits after the
    binary point (as the 53-bit uniforms): the quadrant split is exact"""
    u = np.asarray(u, dtype=np.float64)
    t = u * 4.0
    q = np.floor(t)
    f = t - q  # in [0, 1), exact
    flip = f > 0.5
    g = np.where(flip, 1.0 - f, f)  # in [0, 0.5], exact
    s0, c0 = _sin_cos_small(g * _HALF_PI)
    # angle in the quadrant: f * pi / 2; flip swaps sin and cos
    s1 = np.where(flip, c0, s0)
    c1 = np.where(flip, s0, c0)
    q = q.astype(np.int64) & 3
    sin = np.choose(q, [s1, c1, -s1, -c1])
    cos = np.choose(q, [c1, -s1, -c1, s1])
    return sin, cos


def _round_half_away(x: np.ndarray) -> np.ndarray:
    """MATLAB's round (half away from zero), exact"""
    t = np.trunc(x)
    f = x - t
    return t + np.where(np.abs(f) >= 0.5, np.sign(x), 0.0)


# ------------------------------------------------------------------------------------------
# generators (pktest)
# ------------------------------------------------------------------------------------------


def _size(sz: Sequence[int]) -> tuple[int, ...]:
    sz = tuple(int(s) for s in sz)
    if len(sz) not in (2, 3) or any(s < 0 for s in sz):
        raise ValueError("size: (nz, ny, nx) or (ny, nx), non-negative")
    return sz


def _sigma_z_yx(sigma: float | Sequence[float]) -> tuple[float, float]:
    sg = np.asarray(sigma, dtype=np.float64).reshape(-1)
    if sg.size == 1:
        return float(sg[0]), float(sg[0])
    if sg.size == 2:
        return float(sg[0]), float(sg[1])
    if sg.size == 3 and sg[1] == sg[2]:
        return float(sg[0]), float(sg[1])
    raise ValueError(
        "sigma: a number, (sigma_z, sigma_yx), or (sigma_z, sigma_y, sigma_x) with equal y and x"
    )


def beads(
    sz: Sequence[int],
    positions: Any,
    sigma: float | Sequence[float],
    amplitude: float | Sequence[float] = 1000.0,
    background: float = 0.0,
) -> np.ndarray:
    """single precision volume of Gaussian spots (pktest.beads)

    sz: (nz, ny, nx), or (ny, nx) for an image. positions: #spots x 3, centers (z, y, x), or
    #spots x 2 (y, x) for an image; 0-based (voxel i is at coordinate i), can be fractional.
    sigma: a number or (sigma_z, sigma_yx) in voxels. amplitude: a number or one per spot. Each
    spot is computed in double and added in single, in the order of the rows, as MATLAB does;
    the value at a voxel is MATLAB's for the 1-based center position + 1.
    """
    sz = _size(sz)
    image = len(sz) == 2
    shape = (1, *sz) if image else sz
    pos = np.atleast_2d(np.asarray(positions, dtype=np.float64))
    if pos.size == 0:
        pos = pos.reshape(0, 2 if image else 3)
    if image and pos.shape[1] == 2:
        pos = np.concatenate([np.zeros((pos.shape[0], 1)), pos], axis=1)
    if pos.shape[1] != 3:
        raise ValueError("positions: #spots x 3 (z, y, x), or #spots x 2 (y, x) for an image")
    s_z, s_yx = _sigma_z_yx(sigma)
    amp = np.broadcast_to(np.asarray(amplitude, dtype=np.float64).reshape(-1), (pos.shape[0],))
    z = np.arange(shape[0], dtype=np.float64)[:, None, None]
    y = np.arange(shape[1], dtype=np.float64)[None, :, None]
    x = np.arange(shape[2], dtype=np.float64)[None, None, :]
    vol = np.zeros(shape, dtype=np.float32) + np.float32(background)
    d1 = 2.0 * (s_yx * s_yx)
    d2 = 2.0 * (s_z * s_z)
    for i in range(pos.shape[0]):
        p = pos[i]
        # i - p rounds as MATLAB's (i + 1) - (p + 1): the same real number
        dz = z - p[0]
        dy = y - p[1]
        dx = x - p[2]
        arg = (-(dy * dy + dx * dx)) / d1 - (dz * dz) / d2
        vol = vol + (amp[i] * stable_exp(arg)).astype(np.float32)
    return np.ascontiguousarray(vol[0] if image else vol)


def _gaussian_kernel(sigma: float) -> np.ndarray:
    half = int(math.ceil(2.0 * sigma))
    t = np.arange(-half, half + 1, dtype=np.float64)
    w = stable_exp(-(t * t) / (2.0 * sigma * sigma))
    total = 0.0
    for v in w.tolist():
        total += v
    return w / total


def _matlab_dim_axes(ndim: int) -> list[int]:
    # MATLAB's dimensions 1, 2, 3, ... (y, x, z, ...) as axes of our array
    if ndim < 2:
        return list(range(ndim))
    return [ndim - 2, ndim - 1, *range(ndim - 3, -1, -1)]


def gaussian_filter(vol: np.ndarray, sigma: float) -> np.ndarray:
    """separable Gaussian smoothing in float64 with replicated borders and the kernel size of
    imgaussfilt/imgaussfilt3 (2 ceil(2 sigma) + 1), taps summed in a fixed order, passes along
    y, then x, then z (MATLAB's dimension order). The same bits everywhere, not MATLAB's bits.
    """
    a = np.asarray(vol, dtype=np.float64)
    w = _gaussian_kernel(float(sigma))
    half = (w.size - 1) // 2
    for ax in _matlab_dim_axes(a.ndim):
        n = a.shape[ax]
        if n == 1:
            continue
        idx = np.clip(np.arange(-half, n + half), 0, n - 1)
        p = np.take(a, idx, axis=ax)
        out = np.zeros_like(a)
        for k in range(w.size):
            sl = [slice(None)] * a.ndim
            sl[ax] = slice(k, k + n)
            out = out + w[k] * p[tuple(sl)]
        a = out
    return np.ascontiguousarray(a)


def texture(
    sz: Sequence[int],
    sigma: float = 2.0,
    int_range: Sequence[float] = (100.0, 2000.0),
    rng: CounterRNG | None = None,
) -> np.ndarray:
    """smooth random texture in single (pktest.texture): normal noise, Gaussian smoothing (2D
    for an image size (ny, nx) or nz == 1), scaled to int_range. sz: (nz, ny, nx) or (ny, nx).
    rng: default CounterRNG(0) (the MATLAB tests seed rng(0) for every test)."""
    sz = _size(sz)
    rng = CounterRNG(0) if rng is None else rng
    vol = gaussian_filter(rng.normal(sz), sigma)
    lo, hi = vol.min(), vol.max()
    vol = (vol - lo) / (hi - lo)
    r0, r1 = float(int_range[0]), float(int_range[1])
    return np.ascontiguousarray((r0 + vol * (r1 - r0)).astype(np.float32))


def decon_psf(
    sz: Sequence[int] = (31, 31, 31),
    sigma: Sequence[float] = (2.0, 1.5),
    amplitude: float = 10000.0,
) -> np.ndarray:
    """single Gaussian PSF (z, y, x) with its peak at the FFT center sz // 2 (0-based), the
    first and last 5 z slices zero (pktest.deconPSF). sigma: (sigma_z, sigma_yx)."""
    sz = _size(sz)
    if len(sz) != 3:
        raise ValueError("decon_psf: size (nz, ny, nx)")
    psf = beads(sz, [[s // 2 for s in sz]], sigma, amplitude, 0.0)
    psf[:5] = 0
    psf[max(psf.shape[0] - 5, 0) :] = 0
    return psf


_INT_MAX = {
    "uint8": 2**8 - 1,
    "int8": 2**7 - 1,
    "uint16": 2**16 - 1,
    "int16": 2**15 - 1,
    "uint32": 2**32 - 1,
    "int32": 2**31 - 1,
}


def random_image(sz: Sequence[int], dtype: Any, rng: CounterRNG | None = None) -> np.ndarray:
    """the MATLAB tests' randomImage: logical rand > 0.5; single/double randn * 1000 (cast);
    integer classes randi(intmax) (1 .. intmax). MATLAB class names or NumPy dtypes. sz in
    the port's order ((nz, ny, nx), (ny, nx), leading extra axes)."""
    rng = CounterRNG(0) if rng is None else rng
    name = {"logical": "bool", "single": "float32", "double": "float64"}.get(str(dtype), None)
    dt = np.dtype(name if name else dtype)
    sz = tuple(int(s) for s in sz)
    if dt == np.bool_:
        return rng.uniform(sz) > 0.5
    if dt.kind == "f":
        return np.ascontiguousarray((rng.normal(sz) * 1000.0).astype(dt))
    key = dt.name
    if key not in _INT_MAX:
        raise TypeError(f"random_image: no {dt}")
    return rng.integers(1, _INT_MAX[key], sz, dtype=dt)


def centroid(vol: np.ndarray, center: Sequence[float], radius: int) -> np.ndarray:
    """intensity weighted centroid (z, y, x), 0-based, of the spot near center (z, y, x), from
    the voxels of a cube of half width radius around round(center), minus the window's minimum
    (pktest.centroid). An image (y, x) is a volume with one z plane. Sums in float64 with
    math.fsum (exact, order-free)."""
    v = np.asarray(vol)
    if v.ndim == 2:
        v = v[None]
    c = np.asarray(center, dtype=np.float64).reshape(-1)
    if c.size == 2:
        c = np.concatenate([[0.0], c])
    # MATLAB rounds the 1-based center
    c = (_round_half_away(c + 1.0) - 1.0).astype(np.int64)
    lo = np.maximum(c - radius, 0)
    hi = np.minimum(c + radius + 1, v.shape)
    win = v[lo[0] : hi[0], lo[1] : hi[1], lo[2] : hi[2]].astype(np.float64)
    win = win - win.min()
    g = np.meshgrid(*(np.arange(lo[i], hi[i], dtype=np.float64) for i in range(3)), indexing="ij")
    w = win.reshape(-1)
    total = math.fsum(w.tolist())
    return np.array([math.fsum((gi.reshape(-1) * w).tolist()) / total for gi in g])


def centroids(vol: np.ndarray, positions: Any, radius: int) -> np.ndarray:
    """centroids (#spots x 3, (z, y, x), 0-based) of the spots near positions (pktest.centroids)"""
    pos = np.atleast_2d(np.asarray(positions, dtype=np.float64))
    return np.array([centroid(vol, p, radius) for p in pos]).reshape(-1, 3)


# ------------------------------------------------------------------------------------------
# files
# ------------------------------------------------------------------------------------------


def write_tiff(path: str | os.PathLike, vol: np.ndarray) -> None:
    """write a (z, y, x) volume (one page per z plane) or a (y, x) image as a tiff; the file is
    the one MATLAB's writetiff makes of the same data"""
    import cpptiff

    cpptiff.write_tiff(os.fspath(path), np.ascontiguousarray(vol), num_threads=1)


def write_tiff_matlab(path: str | os.PathLike, vol_yxz: np.ndarray) -> None:
    """write_tiff of an array in MATLAB's form ((y, x, z) or (y, x))"""
    from petakit5d.core.arrays import from_matlab

    write_tiff(path, from_matlab(vol_yxz))


ZARR_LAYOUTS = ("native", "matlab")


def write_zarr(
    path: str | os.PathLike,
    vol: np.ndarray,
    block_size: Sequence[int] = (32, 32, 32),
    layout: str = "native",
) -> None:
    """write a (z, y, x) volume (or a (y, x) image) as a Zarr v2 array.

    layout 'native': the port's layout (shape (nz, ny, nx), order C, chunks block_size).
    'matlab': MATLAB's writezarr layout of the same data (shape (ny, nx, nz), order F, chunks
    block_size in (y, x, z)). block_size is (z, y, x) ((y, x) for an image) either way.
    """
    import cppzarr

    if layout not in ZARR_LAYOUTS:
        raise ValueError(f"layout: one of {ZARR_LAYOUTS}")
    a = np.ascontiguousarray(vol)
    block = [int(b) for b in block_size][-a.ndim :] if a.ndim else []
    if layout == "matlab":
        from petakit5d.core.arrays import to_matlab

        chunks = block[-2:] + block[-3::-1] if a.ndim >= 2 else block
        cppzarr.write_zarr(os.fspath(path), to_matlab(a), chunks=chunks, order="F", num_threads=1)
    else:
        cppzarr.write_zarr(os.fspath(path), a, chunks=block, order="C", num_threads=1)


def write_zarr_matlab(
    path: str | os.PathLike, vol_yxz: np.ndarray, block_size_yxz: Sequence[int] = (32, 32, 32)
) -> None:
    """write an array in MATLAB's form ((y, x, z), chunks (y, x, z)) as MATLAB's writezarr does"""
    import cppzarr

    a = np.asfortranarray(vol_yxz)
    chunks = [int(b) for b in block_size_yxz][: a.ndim]
    cppzarr.write_zarr(os.fspath(path), a, chunks=chunks, order="F", num_threads=1)


class PairedInput(NamedTuple):
    """one benchmark or oracle input for both sides: the path MATLAB reads and the path the
    port reads (the same file for tiff)"""

    matlab: str
    python: str


def write_paired_inputs(
    folder: str | os.PathLike,
    name: str,
    vol: np.ndarray,
    fmt: str = "zarr",
    block_size: Sequence[int] = (32, 32, 32),
) -> PairedInput:
    """write a (z, y, x) volume once per side, each in its own layout: a tiff is shared
    (<folder>/<name>.tif); a zarr is MATLAB's layout in <folder>/matlab/<name>.zarr and native
    in <folder>/python/<name>.zarr (same name, same data, chunks block_size in (z, y, x)).
    Paths use '/' separators."""
    base = os.fspath(folder).replace("\\", "/").rstrip("/")
    if fmt == "tif":
        p = f"{base}/{name}.tif"
        write_tiff(p, vol)
        return PairedInput(p, p)
    if fmt != "zarr":
        raise ValueError("fmt: 'tif' or 'zarr'")
    out = {}
    for side, layout in (("matlab", "matlab"), ("python", "native")):
        Path(base, side).mkdir(parents=True, exist_ok=True)
        out[side] = f"{base}/{side}/{name}.zarr"
        write_zarr(out[side], vol, block_size, layout)
    return PairedInput(out["matlab"], out["python"])


def stitch_tile_name(
    idx: Sequence[int], cam: str = "A", iter: int = 0, ch: int = 0, laser: int = 488
) -> str:
    """the LLSM style tile name of pktest.writeStitchTiles (without extension); idx is the tile
    index (x, y, z), by name as in the file names"""
    cam_str = f"_Cam{cam}" if cam else ""
    x, y, z = (int(v) for v in idx)
    return (
        f"Scan_Iter_{int(iter):04d}_0000{cam_str}_ch{int(ch)}_CAM1_stack0000_{int(laser)}nm_"
        f"0000000msec_0000000000msecAbs_{x:03d}x_{y:03d}y_{z:03d}z_0000t"
    )


def _csv_number(v: float) -> str:
    # writetable writes doubles with up to 15 significant digits
    return f"{float(v):.15g}"


def write_stitch_tiles(
    data_path: str | os.PathLike,
    tiles: Sequence[dict],
    cam: str = "A",
    iter: int = 0,
    ch: int = 0,
    laser: int = 488,
    image_list: str = "ImageList.csv",
    tile_dir: str = "",
    save_format: str = "tif",
    zarr_layout: str = "native",
) -> tuple[list[str], str]:
    """write tiles and add them to the stitching image list (pktest.writeStitchTiles)

    tiles: dicts with 'im' (the tile, (z, y, x)), 'idx' ((x, y, z) tile indices for the name)
    and 'xyz' (stage coordinates in um, x, y, z). Rows are appended to an existing image list.
    The list always names .tif files. zarr tiles are native, or MATLAB's layout with
    zarr_layout='matlab' (chunks 32 x 32 x 32 either way). Returns the tile names (no
    extension) and the image list path (with the data path as given, '/' separated).
    """
    fmt = save_format.lower()
    if fmt not in ("tif", "zarr"):
        raise ValueError("save_format: 'tif' or 'zarr'")
    if zarr_layout not in ZARR_LAYOUTS:
        raise ValueError(f"zarr_layout: one of {ZARR_LAYOUTS}")
    base = os.fspath(data_path)
    out = Path(base, tile_dir) if tile_dir else Path(base)
    out.mkdir(parents=True, exist_ok=True)
    names, rows = [], []
    for t in tiles:
        fsn = stitch_tile_name(t["idx"], cam, iter, ch, laser)
        if fmt == "tif":
            write_tiff(out / f"{fsn}.tif", t["im"])
        else:
            write_zarr(out / f"{fsn}.zarr", t["im"], (32, 32, 32), zarr_layout)
        names.append(fsn)
        rows.append(",".join([f"{fsn}.tif", *(_csv_number(v) for v in t["xyz"])]))
    list_path = base + image_list if base.endswith(("/", "\\")) else os.path.join(base, image_list)
    lp = Path(list_path)
    if lp.exists():
        old = lp.read_text().splitlines()
        lines = old + rows
    else:
        lines = ["Filename,StageX_um_,StageY_um_,StageZ_um_", *rows]
    with open(lp, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return names, list_path


def cut_tiles(
    gt: np.ndarray,
    starts: Any,
    tile_size: Sequence[int],
    pixel_sizes: Sequence[float] = (0.108, 0.108, 0.108),
    stage_xyz: Any = None,
    tile_idx: Any = None,
) -> list[dict]:
    """tiles of tile_size (z, y, x) at 0-based voxel offsets starts ((z0, y0, x0) per row) of a
    ground truth volume (z, y, x) (StitchingTest.cutTiles); for an image (y, x), offsets and
    sizes (y, x). Stage coordinates default to the offsets times the voxel sizes (pixel_sizes
    in (z, y, x)), in x, y, z order as in the image list; tile indices (x, y, z) to the rank of
    each offset among the distinct offsets of its axis."""
    st = np.atleast_2d(np.asarray(starts, dtype=np.int64))
    if st.shape[1] == 2:
        st = np.concatenate([np.zeros((st.shape[0], 1), np.int64), st], axis=1)
    ps = np.asarray(pixel_sizes, dtype=np.float64)
    if ps.size == 2:
        ps = np.concatenate([[ps[0]], ps])
    xyz_starts = st[:, [2, 1, 0]]
    if stage_xyz is None:
        stage = xyz_starts.astype(np.float64) * ps[[2, 1, 0]]
    else:
        stage = np.atleast_2d(np.asarray(stage_xyz, dtype=np.float64))
    if tile_idx is None:
        tidx = np.stack(
            [np.unique(xyz_starts[:, i], return_inverse=True)[1].reshape(-1) for i in range(3)], 1
        )
    else:
        tidx = np.atleast_2d(np.asarray(tile_idx, dtype=np.int64))
    ts = [int(v) for v in tile_size]
    if len(ts) == 2:
        ts = [1, *ts]
    tiles = []
    for i in range(st.shape[0]):
        s = st[i]
        if gt.ndim == 2:
            im = gt[s[1] : s[1] + ts[1], s[2] : s[2] + ts[2]]
        else:
            im = gt[s[0] : s[0] + ts[0], s[1] : s[1] + ts[1], s[2] : s[2] + ts[2]]
        tiles.append(
            {
                "im": np.ascontiguousarray(im),
                "idx": [int(v) for v in tidx[i]],
                "xyz": [float(v) for v in stage[i]],
            }
        )
    return tiles


# ------------------------------------------------------------------------------------------
# hashes
# ------------------------------------------------------------------------------------------


def array_hash(a: np.ndarray) -> str:
    """SHA-256 of an array's dtype, shape and values in column-major order, little-endian"""
    a = np.asarray(a)
    h = hashlib.sha256(f"{a.dtype.str.lstrip('<>|=')}:{a.shape}:".encode())
    h.update(np.asarray(a, dtype=a.dtype.newbyteorder("<")).tobytes(order="F"))
    return h.hexdigest()
