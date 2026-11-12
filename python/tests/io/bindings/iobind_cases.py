"""Cases and helpers for the cpp-tiff / cpp-zarr binding checks (task IO1).

Volumes are (y, x, z) Fortran-ordered, as in MATLAB. cpp-tiff's Python package
takes and returns (z, y, x) C-ordered arrays; the helpers here convert. Inputs
for the MATLAB side are raw little-endian files in Fortran order, so both sides
read the same bytes.
"""

from __future__ import annotations

import json
import os
import zlib
from concurrent.futures import ThreadPoolExecutor
from importlib import metadata
from pathlib import Path

import numpy as np

# numpy dtype name -> MATLAB class
MATLAB_CLASS = {
    "uint8": "uint8",
    "int8": "int8",
    "uint16": "uint16",
    "int16": "int16",
    "uint32": "uint32",
    "int32": "int32",
    "uint64": "uint64",
    "int64": "int64",
    "float32": "single",
    "float64": "double",
    "bool": "logical",
}
TIFF_DTYPES = [
    "uint8",
    "int8",
    "uint16",
    "int16",
    "uint32",
    "int32",
    "uint64",
    "int64",
    "float32",
    "float64",
]
# createzarr/writezarr accept these four; the others go through parallelWriteZarr directly
ZARR_WRAPPER_DTYPES = ["uint8", "uint16", "float32", "float64"]
ZARR_OTHER_DTYPES = ["int8", "int16", "uint32", "int32", "uint64", "int64", "bool"]
TIFF_COMPRESSIONS = ["lzw", "zstd", "none"]
ZARR_COMPRESSORS = ["zstd", "lz4", "lz4hc", "blosclz", "zlib", "gzip"]


def versions() -> dict:
    out = {}
    for pkg in ("cpp-tiff", "cpp-zarr", "numpy"):
        try:
            out[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            out[pkg] = None
    return out


def seed_of(name: str) -> int:
    return zlib.crc32(name.encode())


def make_data(dtype: str, shape, pattern: str = "random", seed: int = 0) -> np.ndarray:
    """A (y, x[, z]) Fortran-ordered array; 'random' covers the type's full range
    (floats with NaN, infinities, -0 and subnormals), 'smooth' compresses like
    microscope data, 'blocks' is zero outside its first octant."""
    shape = tuple(int(s) for s in shape)
    n = int(np.prod(shape))
    rng = np.random.default_rng(seed)
    dt = np.dtype(dtype)
    if pattern == "random":
        if dt.kind == "b":
            v = rng.integers(0, 2, n).astype(bool)
        elif dt.kind in "ui":
            info = np.iinfo(dt)
            v = rng.integers(info.min, info.max, size=n, dtype=dt, endpoint=True)
            if n >= 2:
                v[:2] = [info.min, info.max]
        else:
            v = (rng.standard_normal(n) * 1000).astype(dt)
            specials = np.array(
                [
                    np.nan,
                    np.inf,
                    -np.inf,
                    -0.0,
                    np.finfo(dt).smallest_subnormal,
                    np.finfo(dt).max,
                    np.finfo(dt).min,
                    0.5,
                ],
                dtype=dt,
            )
            k = min(n, specials.size)
            v[:k] = specials[:k]
    elif pattern in ("smooth", "blocks"):
        s3 = shape + (1,) * (3 - len(shape))
        out = np.empty(s3, dtype=dt, order="F")
        y = np.arange(s3[0], dtype=np.float64)[:, None]
        x = np.arange(s3[1], dtype=np.float64)[None, :]
        base = 1000 + 500 * np.sin(y / 7) + 300 * np.cos(x / 11)
        for z in range(s3[2]):
            # (plane by plane: the volumes over 4 GB fit in memory)
            p = base + 20 * z + rng.integers(0, 8, s3[:2])
            if dt.kind == "f":
                p = p.astype(dt)
            elif dt.kind == "b":
                p = p > 1000
            else:
                p = (np.floor(p).astype(np.int64) % (int(np.iinfo(dt).max) + 1)).astype(dt)
            if pattern == "blocks":
                p[(y >= s3[0] / 2) | (x >= s3[1] / 2) | (z >= s3[2] / 2)] = 0
            out[:, :, z] = p
        return out.reshape(shape, order="F")
    elif pattern == "zeros":
        v = np.zeros(n, dt)
    else:
        raise ValueError(pattern)
    return v.reshape(shape, order="F")


def bits_equal(a: np.ndarray, b: np.ndarray) -> bool:
    """Same class, shape and bits (NaN payloads and -0 included), any memory order."""
    if a.dtype != b.dtype or a.shape != b.shape:
        return False
    if a.dtype.kind == "f":
        u = np.dtype(f"u{a.dtype.itemsize}")
        return np.array_equal(
            a.view(u)
            if a.flags.c_contiguous or a.flags.f_contiguous
            else np.ascontiguousarray(a).view(u),
            b.view(u)
            if b.flags.c_contiguous or b.flags.f_contiguous
            else np.ascontiguousarray(b).view(u),
        )
    return np.array_equal(a, b)


def write_raw(path, a: np.ndarray) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    np.asfortranarray(a).astype(a.dtype.newbyteorder("<"), copy=False).ravel(order="F").tofile(path)


def read_raw(path, dtype, shape) -> np.ndarray:
    return (
        np.fromfile(path, dtype=np.dtype(dtype).newbyteorder("<"))
        .astype(dtype)
        .reshape(shape, order="F")
    )


# --- tiff: (y, x, z) Fortran <-> cpp-tiff's (z, y, x) C order --------------------------------


def yxz_to_zyx(a: np.ndarray) -> np.ndarray:
    """The view cpp-tiff's write_tiff takes for a (y, x[, z]) array (it copies it to C order)."""
    return a.transpose(2, 0, 1) if a.ndim == 3 else a


def zyx_to_yxz(a: np.ndarray, num_threads: int = 1) -> np.ndarray:
    """(z, y, x) C-ordered (or (y, x)) to a (y, x, z) Fortran-ordered copy, split by pages."""
    if a.ndim == 2:
        return np.asfortranarray(a)
    nz, ny, nx = a.shape
    out = np.empty((ny, nx, nz), dtype=a.dtype, order="F")
    if num_threads <= 1 or nz < 2:
        out[...] = a.transpose(1, 2, 0)
        return out
    edges = np.linspace(0, nz, min(num_threads, nz) + 1).astype(int)

    def part(i):
        out[:, :, edges[i] : edges[i + 1]] = a[edges[i] : edges[i + 1]].transpose(1, 2, 0)

    with ThreadPoolExecutor(len(edges) - 1) as ex:
        list(ex.map(part, range(len(edges) - 1)))
    return out


def tiff_write_yxz(fn, a: np.ndarray, compression: str = "lzw", num_threads=None) -> None:
    """Through the public write_tiff (which copies the transposed view to C order)."""
    import cpptiff

    cpptiff.write_tiff(str(fn), yxz_to_zyx(a), compression=compression, num_threads=num_threads)


def tiff_write_yxz_direct(fn, a: np.ndarray, compression: str = "lzw", num_threads=None) -> None:
    """Through the extension's own entry with transpose=True, as the mex writes:
    the Fortran (y, x, z) memory is passed as an (x, y, z) array, without a copy.
    Not part of cpp-tiff's public API (see the gap list)."""
    from cpptiff.cpptiff import pybind11_write_tiff

    a = np.asfortranarray(a)
    view = a.transpose(1, 0, 2) if a.ndim == 3 else a.T
    pybind11_write_tiff(
        str(fn), view, True, compression, 0, 0 if num_threads is None else int(num_threads)
    )


def tiff_read_yxz(fn, z_range=None, num_threads=None, convert_threads: int = 1) -> np.ndarray:
    import cpptiff

    return zyx_to_yxz(cpptiff.read_tiff(str(fn), z_range, num_threads=num_threads), convert_threads)


# --- the case matrix -------------------------------------------------------------------------


def tiff_cases(big: bool = False) -> list[dict]:
    cases = []
    shapes = {"s3d": (37, 53, 11), "s2d": (37, 53), "row": (1, 64, 5), "col": (64, 1, 5)}
    for dt in TIFF_DTYPES:
        for sname, shape in shapes.items():
            for comp in TIFF_COMPRESSIONS:
                cases.append(
                    {
                        "name": f"t_{dt}_{sname}_{comp}",
                        "dtype": dt,
                        "shape": shape,
                        "compression": comp,
                        "pattern": "random",
                    }
                )
    extra = [
        ("strips_u16r", "uint16", (1031, 517, 3), "random"),
        ("strips_u16s", "uint16", (1031, 517, 3), "smooth"),
        ("strips_u8s", "uint8", (1031, 517, 3), "smooth"),
        ("strips_f32s", "float32", (1031, 517, 3), "smooth"),
        ("tall_u16s", "uint16", (2050, 3, 2), "smooth"),
        ("pix_u16", "uint16", (1, 1, 1), "random"),
        ("pix2d_u16", "uint16", (1, 1), "random"),
        ("pages70k_u8", "uint8", (4, 4, 70000), "random"),
    ]
    for name, dt, shape, pat in extra:
        for comp in TIFF_COMPRESSIONS:
            cases.append(
                {
                    "name": f"t_{name}_{comp}",
                    "dtype": dt,
                    "shape": shape,
                    "compression": comp,
                    "pattern": pat,
                }
            )
    if big:
        # over 4 GB of pixels: uncompressed (offsets past 2^32) and LZW of random data
        cases.append(
            {
                "name": "t_big_u16_none",
                "dtype": "uint16",
                "shape": (2048, 2048, 520),
                "compression": "none",
                "pattern": "random",
                "big": True,
            }
        )
        cases.append(
            {
                "name": "t_big_u16_lzw",
                "dtype": "uint16",
                "shape": (2048, 2048, 520),
                "compression": "lzw",
                "pattern": "random",
                "big": True,
            }
        )
    for c in cases:
        c["kind"] = "tiff"
        nz = c["shape"][2] if len(c["shape"]) == 3 else 1
        # MATLAB range [first, last] (1-based, inclusive); Python z_range [first-1, last)
        c["range"] = [2, nz - 1] if nz >= 3 else None
    return cases


def zarr_cases(big: bool = False) -> list[dict]:
    im = (37, 41, 13)
    cases = []

    def add(name, **kw):
        c = {
            "kind": "zarr",
            "name": name,
            "dtype": "uint16",
            "shape": im,
            "chunks": (16, 16, 5),
            "cname": "zstd",
            "order": "F",
            "sep": ".",
            "pattern": "random",
            "writer": "writezarr",
            "sparse": False,
        }
        c.update(kw)
        cases.append(c)

    for dt in ZARR_WRAPPER_DTYPES:
        for bname, bs in {
            "nondiv": (16, 16, 5),
            "single": (37, 41, 13),
            "larger": (64, 64, 64),
        }.items():
            add(f"z_{dt}_{bname}", dtype=dt, chunks=bs)
    for dt in ZARR_OTHER_DTYPES:
        add(f"z_{dt}_direct", dtype=dt, writer="direct")
    for cn in ZARR_COMPRESSORS:
        add(f"z_cname_{cn}", cname=cn, writer="createzarr")
    for order in ("F", "C"):
        for csname, cs in {"zsmaller": (16, 8, 5), "zlarger": (8, 16, 13)}.items():
            add(f"z_order{order}_{csname}", order=order, chunks=cs, writer="createzarr")
    add("z_sep_nested", sep="/")
    add("z_2d", shape=(21, 17), chunks=(500, 500, 500))
    add("z_pix", shape=(1, 1, 1), chunks=(500, 500, 500))
    add("z_sparse_false", pattern="blocks", shape=(32, 32, 10), chunks=(16, 16, 5), sparse=False)
    add("z_sparse_true", pattern="blocks", shape=(32, 32, 10), chunks=(16, 16, 5), sparse=True)
    add("z_zeros", pattern="zeros", shape=(32, 32, 10), chunks=(16, 16, 5))
    add("z_smooth_big", pattern="smooth", shape=(1100, 700, 30), chunks=(256, 256, 16))
    add("z_shard", chunks=(32, 32, 10), shard=(16, 16, 5))
    add("z_meta_only", writer="createonly")
    # region writes into an existing array: MATLAB bboxes (1-based, inclusive)
    add(
        "z_regions",
        writer="regions",
        regions=[[5, 7, 2, 30, 20, 9], [10, 1, 1, 37, 15, 13], [17, 17, 6, 17, 17, 6]],
    )
    add(
        "z_regions_zero",
        writer="regions",
        regions=[[1, 1, 1, 37, 41, 13], [1, 1, 1, 16, 16, 5]],
        region_patterns=["random", "zeros"],
    )
    # a double region written into a uint16 array (MATLAB rounds and saturates)
    add(
        "z_convert",
        writer="convert",
        shape=(8, 8, 2),
        chunks=(8, 8, 2),
        regions=[[1, 1, 1, 8, 8, 2]],
    )
    if big:
        add(
            "z_big_u16", pattern="smooth", shape=(2048, 2048, 520), chunks=(256, 256, 256), big=True
        )
    for c in cases:
        c["bboxes"] = (
            [[3, 5, 2, 30, 33, 11], [17, 17, 6, 17, 17, 6], [1, 1, 7, 37, 41, 7]]
            if tuple(c["shape"]) == im
            else []
        )
    return cases


def convert_values(shape) -> np.ndarray:
    """Doubles for a uint16 array: halves, negatives, values over the range, NaN, infinities."""
    special = np.array(
        [
            2.5,
            3.5,
            -0.5,
            -1.0,
            70000.0,
            np.nan,
            np.inf,
            -np.inf,
            65535.4,
            65535.5,
            0.49999999999999994,
            1.5,
            100.7,
            1e10,
            -1e10,
            0.5,
        ],
        dtype=np.float64,
    )
    n = int(np.prod(shape))
    v = np.resize(special, n) + np.repeat(np.arange((n + 15) // 16), 16)[:n] * 1000.25
    v[:16] = special
    return v.reshape(shape, order="F")


def matlab_uint16(v: np.ndarray) -> np.ndarray:
    """MATLAB's uint16(double): round half away from zero, saturate, NaN to 0.
    (floor(|v| + 0.5) would round 0.49999999999999994 up; v - trunc(v) is exact.)"""
    with np.errstate(invalid="ignore"):
        t = np.trunc(v)
        r = t + np.sign(v) * (np.abs(v - t) >= 0.5)
    r = np.where(np.isnan(v), 0.0, r)
    return np.clip(r, 0, 65535).astype(np.uint16)


def matlab_shape(shape) -> tuple:
    """A shape as MATLAB's size reports it: trailing singleton axes past the second dropped."""
    shape = list(shape) + [1] * max(0, 2 - len(shape))
    while len(shape) > 2 and shape[-1] == 1:
        shape.pop()
    return tuple(shape)


def bits_equal_m(a: np.ndarray, b: np.ndarray) -> bool:
    """bits_equal up to trailing singleton axes (how MATLAB compares sizes)."""
    if matlab_shape(a.shape) != matlab_shape(b.shape):
        return False
    return bits_equal(
        a.reshape(matlab_shape(a.shape), order="A"), b.reshape(matlab_shape(b.shape), order="A")
    )


def manifest_write(folder, cases) -> None:
    Path(folder).mkdir(parents=True, exist_ok=True)
    out = []
    for c in cases:
        d = dict(c)
        d["mclass"] = MATLAB_CLASS[c["dtype"]]
        d["shape"] = list(c["shape"])
        for k in ("chunks", "shard"):
            if k in d:
                d[k] = list(d[k])
        out.append(d)
    with open(os.path.join(folder, "manifest.json"), "w") as f:
        json.dump({"versions": versions(), "cases": out}, f, indent=1)
