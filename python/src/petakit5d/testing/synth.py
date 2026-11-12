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

Arrays are MATLAB's: volumes ``(y, x, z)``, Fortran order, positions 1-based ``[y, x, z]``;
random arrays are filled in column-major (MATLAB linear index) order. MATLAB never draws these
numbers itself: oracle cases hand it the arrays (or the files written here), so both sides
read the same bytes.
"""

from __future__ import annotations

import hashlib
import math
import os
from collections.abc import Sequence
from pathlib import Path
from typing import Any

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
    "stitch_tile_name",
    "texture",
    "write_stitch_tiles",
    "write_tiff",
    "write_zarr",
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
    position. Arrays are filled in column-major order (MATLAB's linear index), Fortran-ordered.
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
        return np.asfortranarray(flat.reshape(shape, order="F"))

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


def _size3(sz: Sequence[int]) -> tuple[int, ...]:
    sz = tuple(int(s) for s in sz)
    if len(sz) not in (2, 3) or any(s < 0 for s in sz):
        raise ValueError("size: 2 or 3 non-negative values")
    return sz


def beads(
    sz: Sequence[int],
    positions: Any,
    sigma: float | Sequence[float],
    amplitude: float | Sequence[float] = 1000.0,
    background: float = 0.0,
) -> np.ndarray:
    """single precision volume of Gaussian spots (pktest.beads)

    sz: [ny, nx, nz]. positions: #spots x 3, centers [y, x, z], 1-based, can be fractional.
    sigma: scalar or [sigmaYX, sigmaZ] in voxels. amplitude: scalar or one per spot. Each spot
    is computed in double and added in single, in the order of the rows, as MATLAB does.
    """
    sz = _size3(sz)
    if len(sz) == 2:
        sz = (*sz, 1)
    pos = np.atleast_2d(np.asarray(positions, dtype=np.float64))
    if pos.size == 0:
        pos = pos.reshape(0, 3)
    if pos.shape[1] != 3:
        raise ValueError("positions: #spots x 3")
    sg = np.asarray(sigma, dtype=np.float64).reshape(-1)
    if sg.size == 1:
        sg = np.repeat(sg, 2)
    amp = np.broadcast_to(np.asarray(amplitude, dtype=np.float64).reshape(-1), (pos.shape[0],))
    y = np.arange(1, sz[0] + 1, dtype=np.float64)[:, None, None]
    x = np.arange(1, sz[1] + 1, dtype=np.float64)[None, :, None]
    z = np.arange(1, sz[2] + 1, dtype=np.float64)[None, None, :]
    vol = np.asfortranarray(np.zeros(sz, dtype=np.float32) + np.float32(background))
    d1 = 2.0 * (sg[0] * sg[0])
    d2 = 2.0 * (sg[1] * sg[1])
    for i in range(pos.shape[0]):
        p = pos[i]
        dy = y - p[0]
        dx = x - p[1]
        dz = z - p[2]
        arg = (-(dy * dy + dx * dx)) / d1 - (dz * dz) / d2
        vol = vol + (amp[i] * stable_exp(arg)).astype(np.float32)
    return np.asfortranarray(vol)


def _gaussian_kernel(sigma: float) -> np.ndarray:
    half = int(math.ceil(2.0 * sigma))
    t = np.arange(-half, half + 1, dtype=np.float64)
    w = stable_exp(-(t * t) / (2.0 * sigma * sigma))
    total = 0.0
    for v in w.tolist():
        total += v
    return w / total


def gaussian_filter(vol: np.ndarray, sigma: float) -> np.ndarray:
    """separable Gaussian smoothing in float64 with replicated borders and the kernel size of
    imgaussfilt/imgaussfilt3 (2 ceil(2 sigma) + 1), taps summed in a fixed order. The same
    bits everywhere, not MATLAB's bits."""
    a = np.asarray(vol, dtype=np.float64)
    w = _gaussian_kernel(float(sigma))
    half = (w.size - 1) // 2
    for ax in range(a.ndim):
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
    return np.asfortranarray(a)


def texture(
    sz: Sequence[int],
    sigma: float = 2.0,
    int_range: Sequence[float] = (100.0, 2000.0),
    rng: CounterRNG | None = None,
) -> np.ndarray:
    """smooth random texture in single (pktest.texture): normal noise, Gaussian smoothing (2D
    for a 2-element size or nz == 1), scaled to int_range. rng: default CounterRNG(0) (the
    MATLAB tests seed rng(0) for every test)."""
    sz = _size3(sz)
    rng = CounterRNG(0) if rng is None else rng
    vol = rng.normal(sz)
    if len(sz) == 2 or sz[2] == 1:
        vol = gaussian_filter(vol.reshape(sz[:2], order="F"), sigma).reshape(sz, order="F")
    else:
        vol = gaussian_filter(vol, sigma)
    lo, hi = vol.min(), vol.max()
    vol = (vol - lo) / (hi - lo)
    r0, r1 = float(int_range[0]), float(int_range[1])
    return np.asfortranarray((r0 + vol * (r1 - r0)).astype(np.float32))


def decon_psf(
    sz: Sequence[int] = (31, 31, 31),
    sigma: Sequence[float] = (1.5, 2.0),
    amplitude: float = 10000.0,
) -> np.ndarray:
    """single Gaussian PSF with its peak at the FFT center floor(sz / 2) + 1 (1-based), the
    first and last 5 z slices zero (pktest.deconPSF)"""
    sz = _size3(sz)
    center = [s // 2 + 1 for s in sz]
    psf = beads(sz, [center], sigma, amplitude, 0.0)
    psf[:, :, :5] = 0
    psf[:, :, max(psf.shape[2] - 5, 0) :] = 0
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
    integer classes randi(intmax) (1 .. intmax). MATLAB class names or NumPy dtypes."""
    rng = CounterRNG(0) if rng is None else rng
    name = {"logical": "bool", "single": "float32", "double": "float64"}.get(str(dtype), None)
    dt = np.dtype(name if name else dtype)
    sz = tuple(int(s) for s in sz)
    if dt == np.bool_:
        return rng.uniform(sz) > 0.5
    if dt.kind == "f":
        return np.asfortranarray((rng.normal(sz) * 1000.0).astype(dt))
    key = dt.name
    if key not in _INT_MAX:
        raise TypeError(f"random_image: no {dt}")
    return rng.integers(1, _INT_MAX[key], sz, dtype=dt)


def centroid(vol: np.ndarray, center: Sequence[float], radius: int) -> np.ndarray:
    """intensity weighted centroid [y, x, z] (1-based) of the spot near center, from the
    voxels of a cube of half width radius around round(center), minus the window's minimum
    (pktest.centroid). Sums in float64 with math.fsum (exact, order-free)."""
    c = _round_half_away(np.asarray(center, dtype=np.float64)).astype(np.int64)
    shp = np.asarray(vol.shape[:3] + (1,) * (3 - min(vol.ndim, 3)))[:3]
    lo = np.maximum(c - radius, 1)
    hi = np.minimum(c + radius, shp)
    v = np.asarray(vol)
    if v.ndim == 2:
        v = v[:, :, None]
    win = v[lo[0] - 1 : hi[0], lo[1] - 1 : hi[1], lo[2] - 1 : hi[2]].astype(np.float64)
    win = win - win.min()
    g = np.meshgrid(
        *(np.arange(lo[i], hi[i] + 1, dtype=np.float64) for i in range(3)), indexing="ij"
    )
    w = win.reshape(-1, order="F")
    total = math.fsum(w.tolist())
    return np.array([math.fsum((gi.reshape(-1, order="F") * w).tolist()) / total for gi in g])


def centroids(vol: np.ndarray, positions: Any, radius: int) -> np.ndarray:
    """centroids (#spots x 3) of the spots near positions (pktest.centroids)"""
    pos = np.atleast_2d(np.asarray(positions, dtype=np.float64))
    return np.array([centroid(vol, p, radius) for p in pos]).reshape(-1, 3)


# ------------------------------------------------------------------------------------------
# files
# ------------------------------------------------------------------------------------------


def write_tiff(path: str | os.PathLike, vol: np.ndarray) -> None:
    """write a (y, x[, z]) array as a tiff, one page per z slice (MATLAB's writetiff layout)"""
    import cpptiff

    a = np.asarray(vol)
    data = a if a.ndim == 2 else np.transpose(a, (2, 0, 1))
    cpptiff.write_tiff(os.fspath(path), np.ascontiguousarray(data), num_threads=1)


def write_zarr(
    path: str | os.PathLike, vol: np.ndarray, block_size: Sequence[int] = (32, 32, 32)
) -> None:
    """write a (y, x, z) array as a Zarr v2 array in F order (MATLAB's writezarr layout)"""
    import cppzarr

    a = np.asfortranarray(vol)
    cppzarr.write_zarr(
        os.fspath(path), a, chunks=[int(b) for b in block_size][: a.ndim], order="F", num_threads=1
    )


def stitch_tile_name(
    idx: Sequence[int], cam: str = "A", iter: int = 0, ch: int = 0, laser: int = 488
) -> str:
    """the LLSM style tile name of pktest.writeStitchTiles (without extension); idx is [x, y, z]"""
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
) -> tuple[list[str], str]:
    """write tiles and add them to the stitching image list (pktest.writeStitchTiles)

    tiles: dicts with 'im' (the tile), 'idx' ([x, y, z] tile indices for the name) and 'xyz'
    (stage coordinates in um). Rows are appended to an existing image list. The list always
    names .tif files. Returns the tile names (no extension) and the image list path (with the
    data path as given, '/' separated).
    """
    fmt = save_format.lower()
    if fmt not in ("tif", "zarr"):
        raise ValueError("save_format: 'tif' or 'zarr'")
    base = os.fspath(data_path)
    out = Path(base, tile_dir) if tile_dir else Path(base)
    out.mkdir(parents=True, exist_ok=True)
    names, rows = [], []
    for t in tiles:
        fsn = stitch_tile_name(t["idx"], cam, iter, ch, laser)
        if fmt == "tif":
            write_tiff(out / f"{fsn}.tif", t["im"])
        else:
            write_zarr(out / f"{fsn}.zarr", t["im"])
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
    """tiles of tile_size ([y, x, z]) at 0-based voxel offsets starts ([y0, x0, z0] per row) of
    a ground truth volume (StitchingTest.cutTiles). Stage coordinates default to the offsets
    times the voxel sizes ([y, x, z] sizes), in [x, y, z] order; tile indices [x, y, z] to the
    rank of each offset among the distinct offsets of its axis."""
    st = np.atleast_2d(np.asarray(starts, dtype=np.int64))
    ps = np.asarray(pixel_sizes, dtype=np.float64)
    xyz_starts = st[:, [1, 0, 2]]
    if stage_xyz is None:
        stage = xyz_starts.astype(np.float64) * ps[[1, 0, 2]]
    else:
        stage = np.atleast_2d(np.asarray(stage_xyz, dtype=np.float64))
    if tile_idx is None:
        tidx = np.stack(
            [np.unique(xyz_starts[:, i], return_inverse=True)[1].reshape(-1) for i in range(3)], 1
        )
    else:
        tidx = np.atleast_2d(np.asarray(tile_idx, dtype=np.int64))
    ts = [int(v) for v in tile_size]
    tiles = []
    for i in range(st.shape[0]):
        s = st[i]
        if gt.ndim == 2:
            im = gt[s[0] : s[0] + ts[0], s[1] : s[1] + ts[1]]
        else:
            im = gt[s[0] : s[0] + ts[0], s[1] : s[1] + ts[1], s[2] : s[2] + ts[2]]
        tiles.append(
            {
                "im": np.asfortranarray(im),
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
