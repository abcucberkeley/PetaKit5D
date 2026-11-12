"""Bit-stable inputs for the core oracle cases (the same bytes on every NumPy version).

Only PCG64's raw 64-bit stream is used (``random_raw``, fixed by its algorithm); every value is
derived from it with integer and exactly rounded float operations.
"""

from __future__ import annotations

import numpy as np

INT_CLASSES = ("int8", "uint8", "int16", "uint16", "int32", "uint32", "int64", "uint64")


def raw(seed: int, n: int) -> np.ndarray:
    return np.random.Generator(np.random.PCG64(seed)).bit_generator.random_raw(n).astype(np.uint64)


def uniform(seed: int, n: int) -> np.ndarray:
    """Doubles in [0, 1) with 53 random bits."""
    return (raw(seed, n) >> np.uint64(11)).astype(np.float64) * 2.0**-53


def _pred(x, dt):
    return np.nextafter(dt(x), dt(-np.inf))


def _succ(x, dt):
    return np.nextafter(dt(x), dt(np.inf))


def float_edges(dtype) -> np.ndarray:
    """Edge values for conversions to every integer class: halves, their neighbours, limits."""
    dt = np.float32 if dtype == np.float32 else np.float64
    v = [
        0.0,
        -0.0,
        np.nan,
        np.inf,
        -np.inf,
        2.0**63,
        -(2.0**63),
        2.0**64,
        2.0**53 + 2,
        2.0**24 + 2,
        1e300 if dt is np.float64 else 3e38,
    ]
    for h in (
        0.5,
        1.5,
        2.5,
        3.5,
        126.5,
        127.5,
        254.5,
        255.5,
        32766.5,
        32767.5,
        65534.5,
        65535.5,
        8388607.5,
        2147483647.5,
        4294967295.5,
    ):
        for s in (1.0, -1.0):
            x = dt(s * h)
            v += [x, _pred(x, dt), _succ(x, dt)]
    for cls in INT_CLASSES:
        info = np.iinfo(cls)
        for b in (float(info.min), float(info.max)):
            x = dt(b)
            v += [x, _pred(x, dt), _succ(x, dt), dt(b + 1), dt(b - 1)]
    out = np.array(v, dtype=np.float64).astype(dt)
    return out


def float_values(seed: int, n: int, dtype) -> np.ndarray:
    """Mixed floats: fractions in the uint16 range, exact halves, wide integers, random bit
    patterns (NaN, Inf and subnormals included), then the edge values."""
    dt = np.float32 if dtype == np.float32 else np.float64
    q = n // 4
    u = uniform(seed, 3 * q)
    r = raw(seed + 1, n - 3 * q)
    frac = u[:q] * 140000.0 - 70000.0
    halves = np.floor(u[q : 2 * q] * 140000.0) - 70000.0 + 0.5
    wide = np.floor(u[2 * q :] * 2.0**34) - 2.0**33 + np.floor(u[2 * q :] * 4) / 4
    if dt is np.float32:
        bits = (r & np.uint64(0xFFFFFFFF)).astype(np.uint32).view(np.float32)
        body = np.concatenate(
            [frac.astype(np.float32), halves.astype(np.float32), wide.astype(np.float32), bits]
        )
    else:
        bits = r.view(np.float64)
        body = np.concatenate([frac, halves, wide, bits])
    return np.concatenate([body, float_edges(dt)]).astype(dt)


def int_values(seed: int, n: int, cls: str) -> np.ndarray:
    """Integers over a class's whole range, then its limits and small values."""
    dt = np.dtype(cls)
    r = raw(seed, n)
    body = r.view(np.uint64).astype(np.uint64)
    body = (body >> np.uint64(64 - 8 * dt.itemsize)).astype(np.uint64)
    if dt.kind == "i":
        vals = (
            body.astype(np.int64) - (1 << (8 * dt.itemsize - 1))
            if dt.itemsize < 8
            else body.view(np.int64)
        )
    else:
        vals = body
    info = np.iinfo(dt)
    edges = [
        info.min,
        info.max,
        0,
        1,
        2,
        3,
        5,
        7,
        100,
        127,
        128,
        255,
        256,
        32767,
        32768,
        65535,
        65536,
    ]
    edges += [-1, -2, -3, -5, -127, -128, -129, -32768] if dt.kind == "i" else []
    edges = [e for e in edges if info.min <= e <= info.max]
    return np.concatenate([np.asarray(vals).astype(dt), np.array(edges, dtype=object).astype(dt)])


def int_edges(cls: str) -> np.ndarray:
    """The edge set of a class for operation tables (about 24 values)."""
    info = np.iinfo(cls)
    v = [
        0,
        1,
        2,
        3,
        5,
        7,
        10,
        100,
        127,
        128,
        255,
        256,
        1000,
        32767,
        32768,
        65535,
        65536,
        2**31 - 1,
        2**31,
        2**32 - 1,
        2**53 + 1,
        2**62 + 3,
    ]
    v += [-1, -2, -3, -7, -100, -128, -129, -32768, -(2**31), -(2**53) - 1]
    v += [info.min, info.max, info.max - 1, info.min + 1, info.max // 2, info.max // 2 + 1]
    v = sorted({x for x in v if info.min <= x <= info.max})
    return np.array(v, dtype=object).astype(cls)


DOUBLE_SCALARS = np.array(
    [
        0.0,
        0.5,
        -0.5,
        1.5,
        2.5,
        -2.5,
        1.0 / 3,
        -2.0 / 3,
        0.1,
        0.7,
        3.0,
        255.5,
        1e10,
        -1e10,
        1e-20,
        2.0**53,
        np.nan,
        np.inf,
        -np.inf,
        0.49999999999999994,
    ],
    dtype=np.float64,
)


def bbox_params(seed: int, n: int) -> np.ndarray:
    """Parameter rows for XR_zarrChunkCoordinatesExtraction:
    [vol(3), batch(3), block(3), border(3), sameBatchSize, useBbox, bbox(6)]."""
    r = (raw(seed, n * 20) % np.uint64(1 << 30)).astype(np.int64).reshape(n, 20)
    vol = 1 + r[:, 0:3] % 1500
    batch = 1 + r[:, 3:6] % 1100
    block = 1 + r[:, 6:9] % 300
    border = r[:, 9:12] % 70
    same = r[:, 12] % 2
    use = r[:, 13] % 3 == 0
    s = 1 + r[:, 14:17] % vol
    t = s + r[:, 17:20] % (vol - s + 1)
    bbox = np.concatenate([s, t], axis=1)
    rows = np.concatenate([vol, batch, block, border, same[:, None], use[:, None], bbox], axis=1)
    return rows.astype(np.float64)


def subvol_params(seed: int, n: int) -> np.ndarray:
    """[vol(3), chunk(3), overlap, maxSubVolume (in 1e6), useGoodNum] rows."""
    r = (raw(seed, n * 9) % np.uint64(1 << 30)).astype(np.int64).reshape(n, 9)
    vol = 2 + r[:, 0:3] % 4000
    chunk = 50 + r[:, 3:6] % 1500
    overlap = r[:, 6] % 80
    maxsub = (1 + r[:, 7] % 3000) * 1e6
    good = r[:, 8] % 2
    return np.concatenate(
        [vol, chunk, overlap[:, None], maxsub[:, None], good[:, None]], axis=1
    ).astype(np.float64)


def blocks_params(seed: int, n: int) -> np.ndarray:
    """bboxToBlocks rows: [bbox(6), bboxStart(3), blockSize(3), imageSize(3), borderSize(3)]."""
    r = (raw(seed, n * 18) % np.uint64(1 << 30)).astype(np.int64).reshape(n, 18)
    img = 1 + r[:, 0:3] % 2000
    blk = 1 + r[:, 3:6] % 400
    s = 1 + r[:, 6:9] % img
    t = s + r[:, 9:12] % (img - s + 1)
    start = 1 + r[:, 12:15] % 50
    border = np.where(r[:, 15:18] % 3 == 0, 0, r[:, 15:18] % 40)
    return np.concatenate([s, t, start, blk, img, border], axis=1).astype(np.float64)
