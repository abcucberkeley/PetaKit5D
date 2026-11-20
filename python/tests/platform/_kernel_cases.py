"""Fixed-input cases of the port's kernels for the cross-platform hash test (test_kernel_hashes.py).

Each case runs kernels on inputs made from exact arithmetic (the counter-based generator of
petakit5d.testing.synth, literals, and special values: MATLAB's NaN, a NaN with a payload, both
infinities, -0, subnormals) and returns their outputs. Its digest pins those outputs
(kernel_hashes.json, made on Linux x86-64), so that a kernel that gives other bits on macOS or
Windows fails by name. The cases in RECORDED hold outputs whose last bits follow the system:
FFTs (pocketfft's and FFTW's twiddle factors come from the C library's cos and sin), dense
LAPACK (OpenBLAS's kernels for the CPU, or Accelerate) and the point detection's
Anderson-Darling statistic (the C library's erfc and log); their digests are recorded, not
pinned, and docs/platforms.md has the measured differences. The cases in ARM_NAN (none now)
would hold the outputs of kernels that still give ARM's NaN for an invalid operation (MATLAB's
is x86's); the test expects them to differ there. No output sees two NaNs of other bits
(``specials``), except where the kernel fixes their order itself: their sum's NaN is the
compiler's choice of operand otherwise.

The cases cover every binding of the extension modules (ALL_BINDINGS less EXEMPT_BINDINGS) and
every function of petakit5d.matlabcompat (less EXEMPT_FUNCTIONS); test_kernel_hashes.py checks
that by recording the calls.

usage (from python/):
    python tests/platform/_kernel_cases.py                         the digests (JSON)
    python tests/platform/_kernel_cases.py --write                 write kernel_hashes.json
    python tests/platform/_kernel_cases.py --save F.npz            every case's outputs
    python tests/platform/_kernel_cases.py --compare A.npz B.npz   the differences, per case
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import math
import struct
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np

from petakit5d.core.dtypes import NAN
from petakit5d.testing.synth import CounterRNG

HASHES = Path(__file__).with_name("kernel_hashes.json")

CASES: dict[str, Callable[[], Any]] = {}
# case -> why its digest is recorded, not pinned: FFTs (twiddle factors from the C library),
# LAPACK (OpenBLAS's kernels or Accelerate), the C library's erfc and log (detection.ad_stat)
RECORDED: dict[str, str] = {}
# case -> why its bits differ where the CPU's default NaN is not x86's (ARM): kernels that still
# give the CPU's NaN of an invalid operation (docs/platforms.md); the test expects the failure
ARM_NAN: dict[str, str] = {}


def case(name: str, recorded: str | None = None, arm_nan: str | None = None):
    def add(f):
        if name in CASES:
            raise ValueError(f"case {name} twice")
        CASES[name] = f
        if recorded is not None:
            RECORDED[name] = recorded
        if arm_nan is not None:
            ARM_NAN[name] = arm_nan
        return f

    return add


# -- inputs ------------------------------------------------------------------------------------


def gen(name: str) -> CounterRNG:
    return CounterRNG(0, f"e04/{name}")


def u16(name: str, shape, hi: int = 65535) -> np.ndarray:
    return gen(name).integers(0, hi, shape, dtype=np.uint16)


def ints(name: str, shape, lo: int, hi: int, dtype) -> np.ndarray:
    return gen(name).integers(lo, hi, shape, dtype=dtype)


def double(name: str, shape, scale: float = 1000.0) -> np.ndarray:
    """normal values times scale (Box-Muller on the counter-based generator: the same bits
    everywhere)"""
    return np.ascontiguousarray(gen(name).normal(shape) * scale)


def single(name: str, shape, scale: float = 1000.0) -> np.ndarray:
    return double(name, shape, scale).astype(np.float32)


def uniform(name: str, shape, lo: float = 0.0, hi: float = 1.0, dtype=np.float64) -> np.ndarray:
    u = gen(name).uniform(shape, dtype)
    return np.ascontiguousarray(u * (dtype(hi) - dtype(lo)) + dtype(lo)).astype(dtype)


def logical(name: str, shape, p: float = 0.5) -> np.ndarray:
    return np.ascontiguousarray(gen(name).uniform(shape) < p)


def _payload_nan(dtype) -> Any:
    bits = np.array([0x7FF8000000000123], dtype=np.uint64) if dtype == np.float64 else None
    if bits is None:
        return np.array([0x7FC00123], dtype=np.uint32).view(np.float32)[0]
    return bits.view(np.float64)[0]


def specials(a: np.ndarray, every: int = 97, payload: bool = False) -> np.ndarray:
    """a copy of float data with MATLAB's NaN, both infinities, -0 and +0, subnormals and the
    largest finite values placed every ``every`` elements, and with ``payload`` a NaN with a
    payload in the second place.

    No output may see two NaNs of other bits: where they meet in a sum, the result is the
    compiler's choice of operand, which another build can change (the stable-ABI build
    changed imwarp's). MATLAB's NaN and an invalid operation's are the same bits (x86's), so
    ``payload`` is for kernels whose outputs each see one input value, that fix the order
    themselves (the ordered sums, reduce), or whose footprint can't reach another NaN from the
    payload's place (said where it is used); other kernels with a footprint take the payload NaN
    on an input of its own (``payload_only``)."""
    a = np.array(a, copy=True)
    dt = a.dtype.type
    fi = np.finfo(a.dtype)
    tiny = dt(fi.smallest_subnormal) * dt(3)
    vals = [NAN, _payload_nan(a.dtype), np.inf, -np.inf, -0.0, 0.0, tiny, -tiny, fi.max, -fi.max]
    f = a.reshape(-1)
    for k, v in enumerate(vals):
        i = 3 + k * every
        if i < f.size and (payload or k != 1):
            f[i] = dt(v)
    return a


def payload_only(a: np.ndarray, every: int = 97) -> np.ndarray:
    """a copy of float data with only the NaN with a payload, where ``specials`` places it"""
    a = np.array(a, copy=True)
    if 3 + every < a.size:
        a.reshape(-1)[3 + every] = _payload_nan(a.dtype)
    return a


def cplx(re: np.ndarray, im: np.ndarray) -> np.ndarray:
    """complex data with these parts, set without arithmetic (re + 1j * im makes an invalid
    operation's NaN of 0 * Inf, whose sign is the CPU's)"""
    c = np.empty(re.shape, np.complex64 if re.dtype == np.float32 else np.complex128)
    c.real, c.imag = re, im
    return c


def no_nan(a: np.ndarray) -> np.ndarray:
    """``specials`` without the NaNs"""
    a = specials(a, payload=True)
    a[np.isnan(a)] = 0
    return a


def threads() -> int:
    from petakit5d.core.threads import resolve_threads

    return resolve_threads(None)


# -- digests -----------------------------------------------------------------------------------


def _feed(h, o) -> None:
    if isinstance(o, np.ndarray | np.generic):
        a = np.asarray(o)
        if a.dtype.hasobject:
            raise TypeError("no digest of object arrays")
        h.update(b"A" + a.dtype.str.encode() + repr(a.shape).encode())
        h.update(np.ascontiguousarray(a).tobytes())
    elif isinstance(o, bool):
        h.update(b"T" if o else b"F")
    elif isinstance(o, int):
        h.update(b"I" + str(o).encode())
    elif isinstance(o, float):
        h.update(b"D" + struct.pack("<d", o))
    elif o is None:
        h.update(b"N")
    elif isinstance(o, str):
        h.update(b"S" + o.encode())
    elif isinstance(o, list | tuple):
        h.update(b"(" + str(len(o)).encode())
        for v in o:
            _feed(h, v)
        h.update(b")")
    elif isinstance(o, dict):
        h.update(b"{" + str(len(o)).encode())
        for k in sorted(o):
            _feed(h, str(k))
            _feed(h, o[k])
        h.update(b"}")
    else:
        raise TypeError(f"no digest of {type(o).__name__}")


def digest(o) -> str:
    """the first 20 hex digits of the SHA-256 of the outputs: arrays by dtype, shape and C-order
    bytes (NaN payloads and the sign of zero included), containers element by element"""
    h = hashlib.sha256()
    _feed(h, o)
    return h.hexdigest()[:20]


def leaves(o, path: str = ""):
    """(path, array) of every value of the outputs"""
    if isinstance(o, list | tuple):
        for i, v in enumerate(o):
            yield from leaves(v, f"{path}/{i}")
    elif isinstance(o, dict):
        for k in sorted(o):
            yield from leaves(o[k], f"{path}/{k}")
    elif o is None:
        yield path, np.array("None")
    else:
        yield path, np.asarray(o)


def run(name: str):
    # the special values make NumPy warn of invalid operations and overflows, on purpose
    with np.errstate(all="ignore"):
        return CASES[name]()


# -- core: kernels.ops -------------------------------------------------------------------------


def _ops():
    from petakit5d.kernels import ops

    return ops


VOL = (5, 37, 70)  # (z, y, x), never cubic; 70 columns: vector loops and their tails
BIG = (9, 103, 107)  # 99,189 elements: the fixed-order sums' blocks of 4,096


@case("ops.convert.single_to_classes")
def _():
    K = _ops()
    x = specials(single("conv", VOL, 40000.0), payload=True)
    classes = ("uint8", "int8", "uint16", "int16", "int32", "uint32", "int64", "uint64", "double")
    return [K.convert(x, c) for c in classes] + [K.convert(x, c, mode="cast") for c in classes]


@case("ops.convert.double_and_64_bit_integers")
def _():
    K = _ops()
    d = specials(double("convd", VOL, 1e5), payload=True)
    hi = ints("convhi", (3000,), -(2**30), 2**30, np.int64)
    lo = ints("convlo", (3000,), 0, 2**32 - 1, np.int64)
    i64 = hi * np.int64(2**32) + lo  # exact, up to 2^62
    i64[:4] = [2**53 + 2**29 + 1, -(2**53 + 2**29 + 1), 2**63 - 1, -(2**63)]
    u64 = np.abs(i64[i64 > 0]).astype(np.uint64)
    out = [K.convert(d, c) for c in ("single", "uint16", "int32", "int16")]
    for v in (i64, u64):
        out += [K.convert(v, "single", mode="class"), K.convert(v, "single", mode="cast")]
        out += [K.convert(v, "double"), K.convert(v, "int16")]
    out.append(K.convert(u16("convu", VOL), "single"))
    return out


@case("ops.scale_offset")
def _():
    K = _ops()
    x = specials(single("so", VOL), payload=True)
    u = u16("so16", VOL)
    d = specials(double("sod", VOL), payload=True)
    m = logical("som", VOL)
    return [
        K.scale_offset(x, 1.25, -3.5),
        K.scale_offset(x, 0.1, 7.3, hi=400.0, lo=-0.0, mask=m),
        K.scale_offset(x, 2.0, 0.0, hi=math.nan, lo=0.0, out_class="uint16"),
        K.scale_offset(x, -0.7, 1.0e-3, out_class="int8", mode="cast"),
        K.scale_offset(u, 0.3, 100.5, hi=60000.0, lo=12.0),
        K.scale_offset(u, 1.7, -50.0, out_class="single"),
        K.scale_offset(d, 1e-3, 0.5, hi=1.0, lo=-1.0, mask=m),
    ]


@case("ops.clamp_compare")
def _():
    K = _ops()
    x = specials(single("cl", VOL), payload=True)
    d = specials(double("cld", VOL), payload=True)
    u = u16("clu", VOL)
    out = [K.clamp(x, -0.0, 500.0), K.clamp(x, None, 0.0), K.clamp(d, -100.25, None)]
    out += [K.clamp(u, 100.4, 60000.6), K.clamp(x, math.nan, 3.0)]
    for op in ("eq", "ne", "lt", "le", "gt", "ge"):
        for v in (float(x.flat[11]), -0.0, math.nan, math.inf):
            out.append(K.compare(x, op, v))
        out.append(K.compare(u, op, float(u.flat[5])))
        out.append(K.compare(d, op, 1.0e-300))
    return out


@case("ops.complex_abs")
def _():
    K = _ops()
    re = specials(single("cre", VOL), payload=True)
    im = specials(single("cim", VOL), every=89, payload=True)
    c = cplx(re, im)
    return [K.complex_abs(c), K.complex_abs(c.astype(np.complex128)), K.complex_abs(re)]


@case("ops.sum_double")
def _():
    from petakit5d.core.arrays import to_matlab

    K = _ops()
    out = []
    for n in (7, 43, 1000, 88_999, 99_189):
        out.append(K.sum_double(single(f"sd{n}", (n,))))
    v = single("sdv", BIG)
    out += [K.sum_double(v), K.sum_double(to_matlab(v), matlab_layout=True)]
    out.append(K.sum_double(double("sdd", BIG)))
    out.append(K.sum_double(no_nan(single("sdi", BIG))))  # both infinities: MATLAB's NaN
    out.append(K.sum_double(specials(single("sdn", BIG), payload=True)))  # NaN payloads
    out.append(K.sum_double(-np.zeros(BIG, np.float32)))
    out += [K.sum_double(u16("sdu", BIG)), K.sum_double(logical("sdl", BIG))]
    return out


@case("ops.sum_ordered_and_slices")
def _():
    from petakit5d.core.arrays import to_matlab

    K = _ops()
    v = single("so3", BIG)
    d = double("so3d", BIG)
    out = [K.sum_ordered(v, single=True), K.sum_ordered(v, single=False)]
    out += [
        K.sum_ordered(v, single=True, center=12.5),
        K.sum_ordered(d, single=False, center=-3.25),
    ]
    out += [K.sum_ordered(to_matlab(v), single=True, matlab_layout=True)]
    out += [K.sum_ordered(to_matlab(d), single=False, center=0.5, matlab_layout=True)]
    x3 = single("ss", (7, 1500, 5))
    x3d = double("ssd", (3, 2100, 4))
    c3 = single("ssc", (7, 5))
    for s in (True, False):
        for ch in (True, False):
            out.append(K.sum_slices(x3, single=s, chunked=ch))
    out += [K.sum_slices(x3, single=True, chunked=True, center=c3)]
    out += [K.sum_slices(x3d, single=False, chunked=True, center=double("ssdc", (3, 4)))]
    return out


# -- crop --------------------------------------------------------------------------------------


@case("crop.copy_box")
def _():
    from petakit5d.kernels import crop as C

    v = u16("cr", (7, 40, 33))
    f = single("crf", (7, 40, 33))
    t4 = single("cr4", (3, 5, 20, 17))
    out = [C.crop_3d(v, (1, 3, 2, 6, 38, 30)), C.crop_3d(f, (0, 0, 5, 7, 40, 6))]
    out.append(C.crop_3d(v[0], (2, 3, 30, 31)))
    out.append(C.crop_4d(t4, (0, 1, 2, 3, 3, 4, 19, 15)))
    o = np.zeros((9, 44, 35), np.uint16)
    out.append(C.indexing_3d(o, (1, 2, 1, 8, 40, 34), v[:7, :38, :33]).copy())
    o4 = np.zeros((4, 6, 22, 20), np.float32)
    out.append(C.indexing_4d(o4, (1, 0, 1, 2, 3, 5, 18, 17), t4, (0, 0, 1, 1, 2, 5, 18, 16)))
    return out


# -- mip ---------------------------------------------------------------------------------------


def _mip_inputs():
    return {
        "bool": logical("mipb", VOL),
        "uint8": ints("mip8", VOL, 0, 255, np.uint8),
        "uint16": u16("mip16", VOL),
        "single": specials(single("mipf", VOL), payload=True),
        "double": specials(double("mipd", VOL), payload=True),
        "int16": ints("mipi", VOL, -(2**15), 2**15 - 1, np.int16),
    }


@case("mip.max_pooling_min_box")
def _():
    from petakit5d.kernels import mip as P

    out = []
    for v in _mip_inputs().values():
        out.append(P.max_pooling_3d(v, (2, 3, 4)))
        out.append(P.max_pooling_3d(v, (1, 5, 9)))
        out.append(P.min_bbox_3d(v, (1, 2, 3, 5, 30, 61)))
        out.append(P.min_bbox_3d(v, (0, 0, 0, None, None, None)))
    out.append(P.min_bbox_3d(-np.zeros(VOL, np.float32), (0, 0, 0, 5, 37, 70)))
    last = u16("miplast", VOL, 60000) + np.uint16(1)
    last[-1, -1, -1] = 0  # the only minimum, at the open ends
    out.append(P.min_bbox_3d(last, (0, 0, 0, None, None, None)))
    out.append(
        P.min_bbox_3d(ints("mipi32", VOL, -(2**31), 2**31 - 1, np.int32), (0, 1, 1, 4, 9, 9))
    )
    return out


@case("mip.range_max")
def _():
    from petakit5d.kernels import mip as P

    out = []
    for name, v in _mip_inputs().items():
        if name in ("bool", "int16"):
            continue
        for axis, starts in ((-3, [0]), (-3, [0, 2, 3]), (-2, [0, 9, 10, 30]), (-1, [0, 33])):
            out.append(P.range_max_3d(v, axis, starts))
        out.append(P.range_max_3d_axes(v, ([0, 3], [0, 10, 21], [0, 7, 64])))
    return out


# -- flat field --------------------------------------------------------------------------------


@case("flatfield.correct")
def _():
    from petakit5d.kernels import flatfield as FF

    vol = u16("ffv", VOL)
    flat = uniform("fff", VOL[1:], 0.5, 1.5, np.float32)
    bg = uniform("ffb", VOL[1:], 80.0, 120.0, np.float32)
    flat2, bg2 = flat.copy(), bg.copy()
    flat2.flat[[3, 50, 700]] = [0.0, np.inf, -0.0]
    bg2.flat[[4, 51]] = [np.inf, -np.inf]
    s = specials(single("ffs", VOL), payload=True)
    out = [
        FF.flat_field_correction_3d(vol, flat, bg, cast_data_type=True),
        FF.flat_field_correction_3d(vol, flat, bg, cast_data_type=False),
        FF.flat_field_correction_3d(vol, flat, bg, cast_data_type=True, const_offset=7.25),
        FF.flat_field_correction_3d(vol, flat2, bg2, cast_data_type=False),
        FF.flat_field_correction_3d(s, flat, bg),
        FF.flat_field_correction_3d(s, flat2, bg2, const_offset=-0.1),
        FF.flat_field_correction_3d(vol[0], flat, bg, cast_data_type=True),
    ]
    return out


# -- stitch ------------------------------------------------------------------------------------


def _tiles(name: str, dtype, shape=(3, 4, 21, 19)):
    if dtype == np.uint16:
        t = u16(name, shape)
    elif dtype == np.uint8:
        t = ints(name, shape, 0, 255, np.uint8)
    else:
        t = np.abs(single(name, shape))
    t[logical(name + "z", shape, 0.3)] = 0  # no data
    return t


@case("stitch.feather_blending")
def _():
    from petakit5d.kernels import stitch as S

    w = uniform("fbw", (3, 4, 21, 19), 0.0, 2.0, np.float32)
    w[logical("fbw0", w.shape, 0.1)] = 0
    out = []
    for dt in (np.uint16, np.float32):
        t = _tiles(f"fb{np.dtype(dt).name}", dt)
        out.append(S.feather_blend(t, w))
        o = np.full((6, 30, 25), 7, dtype=dt)
        S.feather_blend_into(t, w, o, (1, 2, 3, 5, 23, 22))
        out.append(o)
    dm = uniform("dm", (4, 9, 7), 0.0, 1.0, np.float32)
    out += [S.resize_distance_map(dm, (7, 20, 15), 2.0), S.resize_distance_map(dm, (4, 9, 8), 1.0)]
    return out


@case("stitch.tile_blending")
def _():
    from petakit5d.kernels import stitch as S

    boxes = [(0, 0, 0, 0, 2, 10, 9), (1, 1, 5, 4, 4, 21, 19), (2, 0, 0, 0, 0, 0, 0)]
    out = []
    for dt in (np.uint8, np.uint16, np.float32):
        t = _tiles(f"tb{np.dtype(dt).name}", dt)
        for method in ("none", "max", "mean", "median"):
            out.append(S.tile_blend(t, boxes, method))
            out.append(S.tile_blend(t, None, method))
    return out


@case("stitch.integral_image_valid_nan")
def _():
    from petakit5d.kernels import stitch as S

    out = [S.integral_image(u16("ii", (6, 20, 17)), (3, 4, 5))]
    out.append(S.integral_image(np.abs(single("iif", (6, 20, 17))), (6, 1, 17)))
    t = _tiles("mt", np.uint16)
    out.append(S.major_tiles_valid(t, np.array([True, False, True])))
    t2 = t.copy()
    t2[0] = 1
    out.append(S.major_tiles_valid(t2, np.array([True, False, False])))
    out.append(S.replace_nan_with_value(specials(single("rn", VOL), payload=True), 3.75))
    out.append(S.replace_nan_with_value(specials(double("rnd", VOL), payload=True), 1e-8))
    rn2 = specials(single("rn2", VOL), payload=True)
    out.append(S.replace_nan_with_value(rn2, 0.1))  # rounded to single
    return out


# -- deskew ------------------------------------------------------------------------------------


@contextlib.contextmanager
def _deskew_module(name: str):
    """the kernels of kernels.deskew from another build of the module (deskew_scalar)"""
    from petakit5d.kernels import deskew as DS

    real = DS._ext

    def other(module: str = "deskew"):
        return real(name if module == "deskew" else module)

    DS._ext = other
    try:
        yield
    finally:
        DS._ext = real


def _deskew_builds() -> list[str]:
    from petakit5d._ext import available

    return ["deskew"] + (["deskew_scalar"] if available("deskew_scalar") else [])


def _geometry(shape, dz, reverse):
    from petakit5d.deskew_rotate.geometry import deskew_rotate_geometry

    return deskew_rotate_geometry(
        shape,
        32.45,
        dz,
        0.108,
        reverse=reverse,
        objective_scan=False,
        x_step_thresh=2.0,
        resample_factor=None,
        interp_method="linear",
    )


def _deskew_outputs():
    from petakit5d.kernels import deskew as DS

    out = []
    for shape, dz, reverse in (((20, 17, 64), 0.5, False), ((37, 45, 29), 0.2, True)):
        g = _geometry(shape, dz, reverse)
        for kind in ("single", "uint16"):
            vol = np.abs(single("dsr", shape)) if kind == "single" else u16("dsr16", shape)
            oz, oy, ox = g.out_size
            for box in ((0, 0, 0, oz, oy, ox), (1, 3, 2, oz - 1, oy - 1, ox - 2)):
                for s16 in (False, True) if kind == "uint16" else (False,):
                    zxy = np.ascontiguousarray(vol.transpose(0, 2, 1))
                    if g.do_interp:
                        xs, st = g.x_step, g.int_stepsize
                        out.append(
                            DS.skewed_space_interp_volume_deskew_rotate_warp(
                                vol, xs, st, reverse, g.tmat, box, s16
                            )
                        )
                        out.append(
                            DS.skewed_space_interp_volume_deskew_rotate_warp_matlab_layout(
                                zxy, xs, st, reverse, g.tmat, box, s16
                            )
                        )
                        out.append(
                            DS.skewed_space_interp_defined_stepsize(vol, xs, st, reverse, s16)
                        )
                        out.append(
                            DS.skewed_space_interp_defined_stepsize_matlab_layout(
                                zxy, xs, st, reverse, s16
                            )
                        )
                    else:
                        out.append(DS.volume_deskew_rotate_warp(vol, g.tmat, box, s16))
                        out.append(
                            DS.volume_deskew_rotate_warp_matlab_layout(zxy, g.tmat, box, s16)
                        )
        out.append(DS.interp_grid_size(shape[0], g.int_stepsize if g.do_interp else 1.0))
    return out


@case("deskew.warp_and_interpolation")
def _():
    """the same bits from every build of the module (deskew_scalar: the code without AVX2)"""
    ref = None
    for name in _deskew_builds():
        with _deskew_module(name):
            got = _deskew_outputs()
        if ref is None:
            ref = got
        elif digest(got) != digest(ref):
            raise AssertionError(f"{name} differs from deskew")
    return ref


@case("deskew.swap_axes_and_page_turn")
def _():
    from petakit5d.deskew_rotate.file import _turn_pages
    from petakit5d.kernels import deskew as DS

    swapped = []
    for name in _deskew_builds():
        with _deskew_module(name):
            swapped.append(DS._swap(single("sw", (6, 23, 41)), threads()))
    if any(s.tobytes() != swapped[0].tobytes() for s in swapped):
        raise AssertionError("the deskew builds swap differently")
    out = swapped[:1]
    vols = (
        u16("tp", (5, 31, 47)),
        single("tpf", (3, 64, 17)),
        ints("tp8", (4, 9, 13), 0, 255, np.uint8),
    )
    for v in vols:
        out.append(_turn_pages(v.copy(), threads()).copy())
    return out


# -- unmix -------------------------------------------------------------------------------------


@case("unmix.channels_sum")
def _():
    from petakit5d.tools.unmix import _subtract, unmix_channels_sum

    chans = [u16("um1", VOL), u16("um2", VOL), u16("um3", VOL)]
    sub = [
        _subtract(specials(single("ums", VOL), payload=True), b, None) for b in (100.25, 0.1, -1e-3)
    ]
    f = [1.0, -0.37, 0.1234567]
    out = [
        unmix_channels_sum(chans, f, [100.0, 98.5, 0.0], 3.0, "uint16"),
        unmix_channels_sum(chans, f, None, None, "single"),
        unmix_channels_sum(chans[:2], [0.5, 0.75], [10.0, 20.0], None, ""),
        unmix_channels_sum(
            [single("um4", VOL), single("um5", VOL)], [1.0, -1.0], None, 0.0, "double"
        ),
        unmix_channels_sum(
            [ints("um6", VOL, 0, 255, np.uint8)] * 2, [0.3, 0.3], [1.0, 2.0], 5.0, "uint8"
        ),
    ]
    return out + sub


@case("xcorr.prepare_finish")
def _():
    """the NCC's steps around its FFTs (the zero-mean inputs and their sums; the normalization of
    a correlation), on a correlation made up here"""
    from petakit5d._ext import load
    from petakit5d.kernels import xcorr as X

    ext = load("xcorr")
    out = []
    for f in (single, double):
        a = np.abs(f("xpa", (12, 30, 26)))
        t = np.ascontiguousarray(a[2:10, 5:25, 4:22] + a.dtype.type(3.0))
        T0, A0, stats = ext.prepare(t, a, threads())
        full = X._full_shape(t, a)
        win = X._win((1, 2, 2, 8, 20, 18), full, "", "")
        corr = f("xpc", (win[5] - win[2] + 1, win[3] - win[0] + 1, win[4] - win[1] + 1), 50.0)
        out += [T0, A0, stats, ext.finish(corr, T0, A0, win, stats, threads())]
    return out


# -- decon (the parts without FFTs) ------------------------------------------------------------


@case("decon.acceleration")
def _():
    from petakit5d.core.arrays import to_matlab
    from petakit5d.kernels import decon as DK

    out = []
    for f in (single, double):
        for shape in (VOL, BIG):  # MATLAB's code below 89,000 elements, the kernel above
            j2, y = f("j2", shape), f("y", shape)
            j4 = j2 - y
            out.append(DK.accel_lambda(j2, y, j4.copy()))
            mj = [to_matlab(a) for a in (j2, y, j4)]
            out.append(DK.accel_lambda(*mj, matlab_layout=True))
        j2, j3 = f("j2", VOL), f("j3", VOL)
        out.append(DK.accel_predict(j2, j3, 0.37))
        out.append(DK.accel_predict(specials(j2, payload=True), j3, 1.0))
        out.append(DK.sum_double(f("ds", BIG)))
        out.append(DK.sum_double(to_matlab(f("ds2", BIG)), matlab_layout=True))
    return out


@case("decon.erosion_and_layouts")
def _():
    from petakit5d.kernels import decon as DK

    m = logical("me", (6, 29, 31), 0.85)
    nh = [np.ones((3, 1), bool), np.ones((1, 3), bool)]
    out = [
        DK.mask_edge_erosion(m, nh, (0, 0), (0, 0), 1, True),
        DK.mask_edge_erosion(m, nh, (0, 0), (0, 0), 2, False),
        DK.mask_edge_erosion(m[0], [np.ones((5, 5), bool)], (0, 0), (0, 0), 0, True),
    ]
    v = single("tl", (6, 29, 31))
    a = DK.to_matlab_layout(v)
    out += [a, DK.from_matlab_layout(a), DK.to_matlab_layout(u16("tl16", (29, 31)))]
    return out


@case("decon.back_projector_parts")
def _():
    from petakit5d.kernels import decon as DK

    m = logical("omw", (7, 19, 23), 0.2)
    m[3, 5:12, 4:15] = True
    out = [DK.omw_convhull(m, [0, 3, 6], [5, 10], [2, 22])]
    out += [DK.omw_components(m, "minsize", 5), DK.omw_components(m, "largest", 2)]
    out += [DK.omw_sortdesc(np.abs(single("sd", (9, 11, 13))))]
    out += [DK.omw_medfilt3(single("mf", (7, 19, 23))), DK.omw_medfilt3(single("mf2", (19, 23)))]
    r3 = 1.0 / math.sqrt(3.0)
    e = np.array(
        [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
        + [[a * r3, b * r3, c * r3] for a in (1, -1) for b in (1, -1) for c in (1, -1)],
        dtype=np.float64,
    )
    r = np.array([3.0, 3.5, 6.0, 5.5, 8.0, 7.25] + [5.0, 4.5, 6.5, 6.0, 5.25, 5.0, 4.0, 7.0])
    out.append(DK.omw_window((11, 17, 21), e, r, (0.6, 1.4)))
    out.append(_window_values())
    return out


def _window_values():
    """the back projector's window values on an ellipsoidal shell of edge points"""
    from petakit5d.decon import backprojector as BP

    z, y, x = np.meshgrid(np.arange(11), np.arange(17), np.arange(21), indexing="ij")
    r2 = 4 * (z - 5) ** 2 + (y - 8) ** 2 + (x - 10) ** 2
    edge = np.ascontiguousarray((r2 >= 30) & (r2 <= 40))
    return BP._window_fast(edge, edge.shape, 0.6, 1.4, threads())


# -- point detection: the Gaussian fits and their tests (C++ detection) ------------------------
# Windows as XR_fitGaussians3D cuts them around candidates (detection.fit): (2 wz + 1, 2 w + 1,
# 2 w + 1) voxels, w = ceil(2 s) and wz = ceil(2 sz), clipped at the volume's faces, NaN for the
# voxels left out (another mask component's, a deskewed frame's empty wedge: 26,814 of MATLAB's
# 42,412 calls had such voxels). The data are Gaussian spots on a background with noise, in a
# camera's counts or in half counts (the background-corrected volume); the starts are the
# callers': the candidate voxel with the prefit's A and c (the master channel), or a fitted
# position with an Ac fit's A, which can be negative (the other channels).

DET_SIGMA = ((1.3, 1.6), (1.5, 1.9))  # [s sz] of the demo's two channels (tests/detection)


def _spots(name: str, shape, spots, c: float, sigma, noise: float, counts: bool = True):
    """a (z, y, x) window of spots A exp(-((x - x0)^2 + (y - y0)^2) / (2 s^2) - (z - z0)^2 /
    (2 sz^2)) on the background c with normal noise (synth's exp: the same bits everywhere), in
    counts (rounded) or half counts; spots are ((x0, y0, z0), A) in the window's 0-based
    coordinates"""
    from petakit5d.testing.synth import stable_exp

    s, sz = sigma
    Z, Y, X = np.meshgrid(*(np.arange(n, dtype=np.float64) for n in shape), indexing="ij")
    w = np.full(shape, c) + noise * gen(name).normal(shape)
    for (x0, y0, z0), A in spots:
        dx, dy, dz = X - x0, Y - y0, Z - z0
        w = w + A * stable_exp(-(dx * dx + dy * dy) / (2 * s * s) - dz * dz / (2 * sz * sz))
    return np.ascontiguousarray(np.floor(w + 0.5) if counts else np.floor(w) + 0.5)


def _fit_cd(windows, p0, max_iter: int = 500, ad_crit: float | None = None):
    """fit_cd's binding on windows packed as detection.fit packs them (rows in MATLAB's order,
    NaN past a window's voxels) at threads() threads, checked against 1 and 3 threads (each fit
    runs on one thread): ([the voxel counts, prm, pstd, C, the stats but A^2], A^2). A^2 (stats
    column 4) is in detection.ad_stat: the C library's erfc and log compute it."""
    from petakit5d._ext import load
    from petakit5d.detection import fit

    vals, dims, n = fit._pack(windows)
    p = np.ascontiguousarray(np.asarray(p0, dtype=np.float64).reshape(-1, 7))
    crit = fit.AD_CRITICAL[(2, 0.05)] if ad_crit is None else float(ad_crit)
    k = p.shape[0]
    runs = []
    for nt in (threads(), 1, 3):
        o = [np.empty((k, 7)), np.empty((k, 7)), np.empty((k, 25)), np.empty((k, 6))]
        load("detection").fit_cd(vals, dims, p, *o, int(max_iter), crit, nt)
        runs.append(o)
    if any(digest(r) != digest(runs[0]) for r in runs[1:]):
        raise AssertionError("fit_cd: other bits at 1 or 3 threads")
    prm, pstd, C, st = runs[0]
    return [n, prm, pstd, C, np.delete(st, 4, axis=1)], st[:, 4].copy()


def _fit_spots():
    """full windows of one spot, amplitudes 15 to 300 over backgrounds of 100 to 1000, in counts
    and half counts, with each channel's sigma, started from the candidate voxel with A and c
    near the spot's; a spot on a zero background with signed zeros"""
    amps = (15.0, 25.0, 40.0, 60.0, 80.0, 100.0, 130.0, 160.0, 200.0, 240.0, 280.0, 300.0)
    bgs = (100.0, 100.0, 150.0, 100.0, 300.0, 100.0, 1000.0, 100.0, 120.0, 500.0, 100.0, 200.0)
    wins, p0 = [], []
    for j, sigma in enumerate(DET_SIGMA):
        w, wz = (math.ceil(2 * v) for v in sigma)
        shape = (2 * wz + 1, 2 * w + 1, 2 * w + 1)
        off = uniform(f"dso{j}", (len(amps), 3), -0.5, 0.5)
        for i, (A, c) in enumerate(zip(amps, bgs, strict=True)):
            xyz = (w + off[i, 0], w + off[i, 1], wz + off[i, 2])
            noise = 2.0 + 2.0 * (i % 4)
            wins.append(_spots(f"ds{j}.{i}", shape, [(xyz, A)], c, sigma, noise, i % 3 != 2))
            p0.append([w, w, wz, 0.9 * A, *sigma, c + 2.0])
    z0 = _spots("dsz", (9, 7, 7), [((3.1, 2.9, 4.2), 30.0)], 0.0, DET_SIGMA[0], 2.0)
    z0.reshape(-1)[np.flatnonzero(z0.reshape(-1) == 0)[::2]] = -0.0
    wins.append(z0)
    p0.append([3, 3, 4, 25.0, *DET_SIGMA[0], 0.0])
    return wins, p0


def _wedge(shape) -> np.ndarray:
    """a deskewed frame's empty wedge in a window: x < z - 2"""
    Z, _, X = np.meshgrid(*(np.arange(n) for n in shape), indexing="ij")
    return X < Z - 2


def _fit_nan():
    """windows with voxels left out (NaN, MATLAB's and with a payload: any NaN is left out):
    another mask component's corner, a deskewed frame's empty wedge, scattered voxels (10 to
    70%); 48 to 42 voxels left (the fit's sums: in index order up to 43, in four lanes above,
    with their tails), 8 and 7 (the fewest the fit takes), 6 and none (NaN results); windows
    clipped at the volume's faces, a volume of 2 planes"""
    sigma = DET_SIGMA[0]
    wins, p0 = [], []
    full = (9, 7, 7)
    Z, Y, X = np.meshgrid(*(np.arange(n) for n in full), indexing="ij")
    regions = [
        (X - 6) ** 2 + (Y - 6) ** 2 + (Z - 7) ** 2 // 2 < 8,
        _wedge(full),
        logical("dnr1", full, 0.1),
        logical("dnr2", full, 0.35),
        logical("dnr3", full, 0.7),
    ]
    for i, out in enumerate(regions):
        w = _spots(f"dn{i}", full, [((3.3, 2.8, 4.2), 120.0)], 100.0, sigma, 5.0)
        w[out] = NAN
        w.reshape(-1)[np.flatnonzero(out.reshape(-1))[::3]] = _payload_nan(np.float64)
        wins.append(w)
        p0.append([3, 3, 4, 100.0, *sigma, 100.0])
    small = _spots("dnc", (3, 4, 4), [((1.6, 1.4, 1.2), 70.0)], 100.0, sigma, 3.0)
    order = np.arange(48).reshape(3, 4, 4).transpose(0, 2, 1).reshape(-1)  # MATLAB's voxel order
    for left in (48, 47, 46, 45, 44, 43, 42, 8, 7, 6, 0):
        w = small.copy()
        w.reshape(-1)[order[left:]] = NAN
        wins.append(w)
        p0.append([1, 1, 1, 50.0, *sigma, 100.0])
    # clipped at the faces: (shape, the spot, the candidate voxel)
    for i, (shape, xyz, start) in enumerate(
        (
            ((5, 7, 7), (3.2, 2.7, 0.6), (3, 3, 0)),
            ((9, 4, 7), (3.4, 0.3, 4.1), (3, 0, 4)),
            ((9, 7, 4), (3.4, 2.9, 3.8), (3, 3, 4)),
            ((5, 4, 4), (0.4, 0.2, 0.3), (0, 0, 0)),
            ((2, 7, 7), (3.1, 3.2, 0.7), (3, 3, 1)),
        )
    ):
        w = _spots(f"dnf{i}", shape, [(xyz, 90.0)], 100.0, sigma, 4.0, i % 2 == 0)
        wins.append(w)
        p0.append([*start, 80.0, *sigma, 100.0])
    return wins, p0


def _pairs() -> np.ndarray:
    """bright components of 2 voxels on a lattice (x and y even, z in pairs 3 apart) over weak
    noise: more starts than the fit takes (16)"""
    full = (9, 7, 7)
    Z, Y, X = np.meshgrid(*(np.arange(n) for n in full), indexing="ij")
    w = _spots("drp", full, [], 100.0, DET_SIGMA[0], 1.0)
    dots = (X % 2 == 0) & (Y % 2 == 0) & (Z % 3 != 2)
    w[dots] += 50.0 + np.floor(uniform("drpd", full, 0.0, 20.0))[dots]
    return w


def _fit_restarts():
    """starts the fit leaves for the window's LoG maxima and bright components: a negative A (an
    Ac fit's, on another channel), also at a fitted position, an A above the data's range, a
    start on the window's face with the step outward (the spot's centre outside the clipped
    window), two spots, noise only, a larger window of strong noise, the candidate on the LoG
    maximum (left out of the starts), voxels left out, a dark spot (31 iterations without a
    better fit), more than 16 starts (cut to 16)"""
    sigma = DET_SIGMA[0]
    full, mid = (9, 7, 7), (3, 3, 4)
    wins, p0 = [], []

    def add(name, shape, spots, start, A, noise=5.0, counts=True, out=None):
        w = _spots(name, shape, spots, 100.0, sigma, noise, counts)
        if out is not None:
            w[out] = NAN
        wins.append(w)
        p0.append([*start, A, *sigma, 100.0])

    for i, A0 in enumerate((-20.0, -3.5, -0.25)):
        add(f"drn{i}", full, [((3.2, 2.9, 4.1), 80.0 + 40 * i)], mid, A0)
    add("drf", full, [((3.27, 2.91, 4.06), 120.0)], (3.27, 2.91, 4.06), -12.5)
    add("dra", full, [((2.8, 3.3, 3.7), 90.0)], mid, 5000.0, 4.0)
    add("dro", (9, 7, 4), [((-0.7, 3.1, 4.2), 150.0)], (0, 3, 4), 100.0, 4.0)
    add("dry", (9, 4, 7), [((3.0, 3.8, 4.0), 150.0)], mid, 100.0, 4.0)
    for A0 in (-10.0, 50.0):
        add("dr2", full, [((1.5, 2.0, 3.0), 90.0), ((5.2, 4.5, 5.5), 60.0)], mid, A0, 4.0, False)
    for A0 in (-3.0, 2.0):
        add("drz", full, [], mid, A0, 6.0)
    add("drm", (9, 9, 9), [((4.0, 4.0, 4.0), 10.0)], (4, 4, 4), -1.0, 25.0)
    add("drl", full, [((3.0, 3.0, 4.0), 200.0)], mid, -5.0, 2.0)
    add("drw", full, [((3.4, 3.1, 4.3), 110.0)], mid, -8.0, out=_wedge(full))
    add("drs", (3, 4, 4), [((1.6, 1.4, 1.2), 70.0)], (1, 1, 1), -2.0, 3.0)
    for A0 in (30.0, -30.0):
        add("drd", full, [((3.2, 2.9, 4.1), -60.0)], mid, A0, 3.0)
    for A0 in (-5.0, 20.0):
        wins.append(_pairs())
        p0.append([*mid, A0, *sigma, 100.0])
    return wins, p0


@case("detection.fit_cd.spots")
def _():
    return _fit_cd(*_fit_spots())[0]


@case("detection.fit_cd.nan_voxels")
def _():
    return _fit_cd(*_fit_nan())[0]


@case("detection.fit_cd.restarts")
def _():
    return _fit_cd(*_fit_restarts())[0]


@case("detection.fit_cd.iterations")
def _():
    """the iteration limit (none: the start's own results; 1, 2, 5) and the critical value of the
    Anderson-Darling test"""
    wins, p0 = _fit_spots()
    rw, rp = _fit_restarts()
    wins, p0 = wins[::4] + rw[::3], p0[::4] + rp[::3]
    out = [_fit_cd(wins, p0, m)[0] for m in (0, 1, 2, 5)]
    return out + [_fit_cd(wins, p0, 500, 0.5)[0]]


@case(
    "detection.ad_stat",
    recorded="the C library's erfc and log, which the Anderson-Darling statistic (adtest_mex's, "
    "and the fits' A^2) takes",
)
def _():
    """the Anderson-Darling statistic of samples in their given order (tests/detection's sorted
    normal samples and one out of order, against three sigmas), of 1 to 441 values (a window's),
    of counts against their mean, with ties and signed zeros, of a 2-D array (its C order); the
    fits' A^2 of the sorted residuals"""
    from petakit5d._ext import load
    from petakit5d.detection import fit

    out = []
    g = gen("ads").normal((6, 294))
    for i, s in enumerate((0.7, 0.9, 1.0, 1.1, 1.3, 1.0)):
        x = np.sort(g[i] * s) if i < 5 else g[i] * s
        for sig in (0.8, 1.0, 1.25):
            out += [fit.ad_statistic(x, 0.0, sig), fit.adtest(x, 0.0, sig)]
    for n in (1, 2, 7, 43, 44, 441):
        out.append(fit.ad_statistic(np.sort(double(f"adn{n}", (n,), 1.0)), 0.0, 1.0))
    c = np.sort(np.floor(double("adc", (441,), 5.0) + 100.0))
    out += [fit.ad_statistic(c, 100.0, 5.0), fit.ad_statistic(c, 99.5, 4.75)]
    t = np.sort(np.floor(double("adt", (200,), 2.0)))
    t[t == 0] = np.where(np.arange(int(np.sum(t == 0))) % 2 == 0, -0.0, 0.0)
    out += [fit.ad_statistic(t, 0.0, 2.0), fit.adtest(t, 0.0, 2.0)]
    out.append(
        load("detection").ad_stat(np.sort(double("ad2", (294,), 1.0)).reshape(7, 42), 0.0, 1.0)
    )
    return out + [_fit_cd(*f())[1] for f in (_fit_spots, _fit_nan, _fit_restarts)]


def _bits(*words: int) -> np.ndarray:
    return np.array(words, dtype=np.uint64).view(np.float64)


@case("detection.exp")
def _():
    """the fits' exp (fdlibm's, e at 1, as MATLAB's) on each of its paths, at both sides of their
    bounds: |x| below 2^-28, up to 0.5 ln 2, below 1.5 ln 2, above; the overflow and underflow
    thresholds, subnormal results; the special values; the fits' arguments (-20 to 0)"""
    from petakit5d._ext import load

    D = load("detection")
    # the high words that bound the paths (|x| < 2^-28, <= 0.5 ln 2, < 1.5 ln 2, >= 709.78)
    words = []
    for h in (0x3E300000, 0x3FD62E43, 0x3FF0A2B2, 0x40862E42):
        words += [(h << 32) - 1, h << 32, (h << 32) + 1, ((h - 1) << 32) + 0x80000000]
    # the overflow and underflow thresholds and their neighbours, 1 (e) and its neighbours
    for t in np.array([7.09782712893383973096e02, -7.45133219101941108420e02, 1.0]).view(np.uint64):
        words += [int(t) - 1, int(t), int(t) + 1]
    b = _bits(*words)
    x = np.concatenate(
        [
            b,
            -b,
            -uniform("dexs", (300,), 708.0, 746.0),  # subnormal results
            uniform("dexa", (3000,), -746.0, 710.0),
            -uniform("dexf", (2000,), 0.0, 20.0),
            specials(uniform("dexx", (1000,), -3.0, 3.0), payload=True),
        ]
    )
    return np.array([D.exp(float(v)) for v in x])


# -- matlabcompat: filters, resize, warp, morphology, labels, distances -----------------------


@case("filters.convolutions")
def _():
    from petakit5d.matlabcompat import filters as F

    v = single("fv", (6, 23, 31))
    d = double("fvd", (6, 23, 31))
    im = specials(single("fim", (37, 70)), payload=True)  # no other NaN in the filters' reach
    k3 = np.abs(single("fk", (3, 5, 4)))
    out = [F.gaussian_kernel_1d(1.7, 11), F.fspecial("gaussian", (7, 5), 1.3)]
    out += [F.imgaussfilt3(v, 2.5), F.imgaussfilt3(d, (1.0, 2.0, 0.5)), F.imgaussfilt(im, 1.3)]
    out += [F.imgaussfilt(single("fim2", (37, 70)), 100)]
    k3 = k3.astype(np.float64)
    out += [F.imfilter(v, k3, "conv", "same", "symmetric"), F.imfilter(d, k3, "replicate")]
    out += [F.imfilter(u16("fu", (6, 23, 31)), np.ones((3, 3, 3)) / 27.0)]
    out += [F.convn(v, k3, "same"), F.convn(d, k3, "full"), F.convn(v, k3, "valid")]
    out += [F.conv2(im, np.ones((3, 3)), "same"), F.conv2(d[1], double("fk2", (4, 5)), "valid")]
    out += [F.conv2(np.array([1.0, 2.0, 1.0]), np.array([1.0, 0.0, -1.0]), d[0], "same")]
    return out


@case("resize.imresize")
def _():
    from petakit5d.matlabcompat import resize as R

    v = single("rv", (9, 21, 17))
    u = u16("ru", (9, 21, 17))
    out = [R.cubic(np.arange(-25, 26) / 10.0), R.kernel_values("triangle", np.arange(-8, 9) / 4.0)]
    out += [R.kernel_values("box", np.arange(-8, 9) / 4.0)]
    out += [R.contributions(17, 40, 40 / 17, "cubic", 4.0, True)]
    for method in ("nearest", "linear", "cubic"):
        out += [R.imresize3(v, (5, 30, 9), method), R.imresize3(u, 0.6, method)]
        m2 = {"linear": "bilinear", "cubic": "bicubic"}.get(method, method)
        out += [R.imresize(v[0], (40, 11), m2), R.imresize(u[0], 1.7, m2)]
    out += [R.imresize(ints("ri8", (21, 17), -128, 127, np.int8), 0.5, "bicubic")]
    out += [R.imresize(ints("ri32", (21, 17), -(2**31), 2**31 - 1, np.int32), 1.5, "bilinear")]
    out += [R.imresize(ints("ru32", (21, 17), 0, 2**32 - 1, np.uint32), 0.7, "bicubic")]
    out += [R.imresize(ints("ri16", (21, 17), -(2**15), 2**15 - 1, np.int16), 2.0, "bilinear")]
    out += [R.imresize(ints("ru8", (21, 17), 0, 255, np.uint8), 0.4, "nearest")]
    out += [R.imresize(double("rd", (21, 17)), 0.45, "bicubic")]
    out += [R.imresize(ints("ru8b", (21, 17), 0, 255, np.uint8), 1.3, "bilinear")]
    out += [R.imresize(double("rdn", (21, 17)), 1.3, "nearest")]
    out += [R.imresize3(double("r3d", (9, 21, 17)), 0.5, "nearest")]
    out += [R.imresize3_fast(v, (5, 30, 9), "linear", None, "single")]
    out += [R.imresize3_fast(u, (4, 11, 30), "cubic", None, "uint16")]
    out += [R.imresize3_average(v, (2, 2, 2)), R.imresize3_average(u, (3, 4, 2))]
    w, idx = R.contributions(21, 9, 9 / 21, "triangle", 2.0, True)
    out += [R.resize_axis(v, -2, w, idx)]
    # the kernel's own pow (fdlibm's, inside the cubic weights)
    from petakit5d._ext import load

    x = np.concatenate([uniform("rpw", (999,), 0.0, 30.0), [0.0, -0.0, np.inf, NAN]])
    for y in (3.0, 2.5, -1.5, 0.5, 1.0e-3):
        o = np.empty_like(x)
        load("resize").fdlibm_pow(x, y, o)
        out.append(o)
    xn = -x  # negative bases: whole exponents only (an invalid operation otherwise)
    for y in (3.0, -2.0, 0.0):
        o = np.empty_like(xn)
        load("resize").fdlibm_pow(xn, y, o)
        out.append(o)
    return out


@case("specials.resize_pow")
def _():
    """the resize kernels' fdlibm pow of negative bases with fractional exponents (0 / 0); the
    cubic weights take |x|^3 only, so no output has it"""
    from petakit5d._ext import load

    x = -uniform("rpn", (99,), 0.5, 30.0)
    out = []
    for y in (2.5, -1.5, 0.5):
        o = np.empty_like(x)
        load("resize").fdlibm_pow(x, y, o)
        out.append(o)
    return out


@case("warp.imwarp")
def _():
    from petakit5d.matlabcompat import warp as W

    T = np.array(
        [[0.875, 0.0625, 0.0, 0.0], [-0.125, 1.0625, 0.25, 0.0], [0.0, 0.03125, 0.9375, 0.0],
         [1.5, -2.25, 0.75, 1.0]]
    )  # fmt: skip
    out = [W.mrdivide(np.arange(16.0).reshape(4, 4) + np.eye(4), T)]
    out.append(W.composite_matrix(T, (8, 22, 18), None, None))
    for vol in (
        single("wv", (9, 21, 17)),
        u16("wu", (9, 21, 17)),
        ints("w8", (9, 21, 17), 0, 255, np.uint8),
    ):
        for method in ("nearest", "linear", "cubic"):
            out.append(W.imwarp(vol, T, method, output_size=(8, 22, 18), fill_value=0.0))
    out.append(
        W.imwarp(
            double("wd", (9, 21, 17)), T, "linear", output_size=(9, 21, 17), fill_value=math.nan
        )
    )
    wi = ints("wi", (9, 21, 17), -(2**15), 2**15 - 1, np.int16)
    out.append(W.imwarp(wi, T, "cubic", output_size=(7, 20, 19)))
    return out


@case("morphology.binary_and_median")
def _():
    from petakit5d.matlabcompat import morphology as M

    b3 = logical("mb", (7, 23, 29), 0.6)
    b2 = logical("mb2", (37, 70), 0.6)
    se3 = M.strel(np.ones((3, 3, 3), bool))
    se2 = M.strel(np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], bool))
    out = [M.imerode(b3, se3), M.imdilate(b3, se3), M.imopen(b2, se2), M.imclose(b2, se2)]
    out += [M.bwconvhull(logical("mh", (21, 26), 0.05), "union")]
    out += [M.bwconvhull(logical("mh2", (21, 26), 0.03), "objects", 8)]
    out += [M.bwconvhull_slices(b3, 0, [0, 3]), M.bwconvhull_slices(b3, 2, [1, 28])]
    for v in (
        specials(single("md", (21, 26)), payload=True),  # no other NaN in the windows' reach
        u16("md16", (21, 26)),
        ints("md8", (21, 26), -128, 127, np.int8),
        double("mdd", (21, 26)),
        logical("mdb", (21, 26)),
        ints("mdu8", (21, 26), 0, 255, np.uint8),
        ints("mdi16", (21, 26), -(2**15), 2**15 - 1, np.int16),
        ints("mdi32", (21, 26), -(2**31), 2**31 - 1, np.int32),
        ints("mdu32", (21, 26), 0, 2**32 - 1, np.uint32),
    ):
        out += [M.medfilt2(v, (3, 3)), M.medfilt2(v, (4, 5), "symmetric")]
    out += [M.medfilt3(single("m3", (7, 23, 29))), M.medfilt3(u16("m316", (7, 23, 29)), (3, 5, 3))]
    out += [M.medfilt3(single("m3s", (7, 23, 29)), (3, 3, 3), "zeros")]
    out += [M.medfilt3(specials(single("m3n", (7, 23, 29)), payload=True), (3, 3, 3), "zeros")]
    return out


@case("labels_and_distances")
def _():
    from petakit5d.matlabcompat import distance as D
    from petakit5d.matlabcompat import label as L

    b3 = logical("lb", (7, 23, 29), 0.3)
    b2 = logical("lb2", (37, 70), 0.4)
    out = [L.conndef(3, "minimal"), L.conndef(2, "maximal"), L.connectivity_matrix(18)]
    out += [L.bwlabeln(b3), L.bwlabeln(b3, 6), L.bwlabel(b2, 4), L.bwlabel(b2, 8)]
    cc = L.bwconncomp(b3)
    out += [cc.num_objects, [np.asarray(p) for p in cc.pixel_idx_list], L.labelmatrix(cc)]
    out += [L.bwareaopen(b3, 6), L.bwareaopen(b2, 4, 4), L.imfill(b2, "holes")]
    out += [D.bwdist(b3), D.bwdist(b2), D.bwdist(u16("dist", (5, 9, 11), 3))]
    return out


@case("labels.int64_kernels")
def _():
    """the label kernels' int64 builds (for 2^31 voxels or more) on a small volume: the int32
    builds' labels"""
    from petakit5d._ext import load
    from petakit5d.core.threads import resolve_threads
    from petakit5d.matlabcompat import label as L

    m = load("morph")
    nt = resolve_threads(None)
    v = logical("l64", (7, 23, 29), 0.3)
    off = L._offsets(L.connectivity_matrix(26))
    moff = L._matlab_offsets(L.connectivity_matrix(18))
    out = []
    for kind in ("any", "first"):
        labs = []
        for dt in (np.int32, np.int64):
            lab = np.empty(v.shape, dt)
            n = getattr(m, f"label_{kind}_{np.dtype(dt).name}")(v, lab, off, nt)
            groups = getattr(m, f"group_by_label_{np.dtype(dt).name}")(lab, n)
            labs.append((lab.astype(np.int64), n, [np.asarray(g, np.int64) for g in groups]))
        if digest(labs[0]) != digest(labs[1]):
            raise AssertionError(f"label_{kind}: int64 and int32 builds differ")
        out.append(labs[1])
    ls = []
    for f in (m.bwlabeln_int32, m.bwlabeln_int64):
        o = np.empty(v.shape, np.float64)
        ls.append((o, f(v, o, moff, nt)))
    if digest(ls[0]) != digest(ls[1]):
        raise AssertionError("bwlabeln: int64 and int32 builds differ")
    return out + [ls[1]]


# -- matlabcompat: elementary functions, sums and statistics, solvers --------------------------


def _angles():
    """angles of every octant, large arguments (argument reduction) and special values"""
    a = uniform("ang", (4000,), -20.0, 20.0)
    big = np.array([1e5, -3.0e5 + 0.125, 2.0**19 * math.pi, 1e-300, -1e-310, 2**-1074, 0.0, -0.0])
    return np.concatenate([a, big, [NAN, np.inf, -np.inf]])


@case("elementary.trigonometry")
def _():
    from petakit5d.matlabcompat import elementary as E

    x = _angles()
    deg = np.concatenate([uniform("deg", (2000,), -720.0, 720.0), np.arange(-360.0, 361.0, 15.0)])
    out = [E.sin(x), E.cos(x), E.tan(x), E.cot(x), E.atan(x), E.fdlibm_sin(x), E.fdlibm_cos(x)]
    out += [E.sind(deg), E.cosd(deg), E.atan2(x[:-3], x[::-1][3:]), E.atan2(-0.0, -0.0)]
    small = uniform("ksm", (500,), -0.78, 0.78)
    out += [E.kernel_sin(small, 0.0, 0), E.kernel_cos(small, 1e-18), E.kernel_tan(small, 0.0, 1)]
    out += [E.kernel_tan(small, 0.0, -1), [E.rem_pio2(float(v)) for v in x[:50]]]
    out += [E.sin(x.astype(np.float32)), E.cos(x.astype(np.float32))]
    return out


@case("elementary.exp_power_sqrt")
def _():
    from petakit5d.matlabcompat import elementary as E

    x = np.concatenate(
        [uniform("ex", (3000,), -745.0, 709.0), [1.0, 0.0, -0.0, NAN, np.inf, -np.inf]]
    )
    p = uniform("pw", (3000,), 0.0, 50.0)
    out = [E.exp(x), E.exp(x.astype(np.float32)), E.power(p, 2.5), E.power(p, -1.0 / 3.0)]
    out += [E.power(p, 3.0), E.power(-p, 2.0), E.power(p.astype(np.float32), 0.5)]
    out += [E.fast_power(uniform("fp", (500,), 0.0, 1.0), 10), E.fast_power(np.float32(0.97), 7)]
    out += [E.sqrt(p), E.sqrt(-p[:5]), E.hypot(x[:3000], p), E.hypot(np.float32(3), np.float32(4))]
    out += [E.tukeywin(110, 0.5), E.tukeywin(7, 1.0), E.tukeywin(5, 0.0)]
    a = u16("pad", (3, 4, 5))
    out += [E.padarray(a, (1, 0, 2), method="replicate", direction="post")]
    out += [E.padarray(a, (1, 2, 1), 7), E.padarray(a[0], (1, 2), method="symmetric")]
    out += [E.circshift(a, (1, -2), (-3, -1))]
    return out


@case("reduce.sums_and_orders")
def _():
    from petakit5d.matlabcompat import reduce as Rd

    v = specials(single("rdv", BIG), payload=True)
    d = double("rdd", (6, 37, 70))
    u = u16("rdu", (6, 37, 70))
    out = [Rd.default_axis(d), Rd.sum(v), Rd.sum(d, -3), Rd.sum(d, -2), Rd.sum(u, -1), Rd.sum(u)]
    out += [Rd.sum(v, omitnan=True), Rd.sum(d, None, outtype="native"), Rd.mean(d), Rd.mean(u, -3)]
    out += [Rd.mean(v, omitnan=True), Rd.cumsum(d, -2), Rd.cumsum(v.reshape(-1)[:5000], 0)]
    out += [Rd.max(v), Rd.min(v), Rd.max(d, -1), Rd.min(u, -3), Rd.max(v, omitnan=False)]
    out += [Rd.maximum(d, np.float32(1.5)), Rd.minimum(u, 300.0), Rd.maximum(v, -0.0)]
    out += [Rd.sort(d, -1), Rd.sort(v.reshape(-1), descend=True), Rd.unique(u[0])]
    out += [Rd.unique(ints("un", (500,), 0, 50, np.int32), stable=True)]
    return out


@case("stats")
def _():
    from petakit5d.matlabcompat import stats as St

    v = single("stv", (6, 37, 70))
    d = double("std", (6, 37, 70))
    u = u16("stu", (6, 37, 70), 1000)
    out = [
        St.var(v),
        St.var(d, 1, -3),
        St.std(d),
        St.std(v, 0, -2),
        St.std(specials(v, payload=True), omitnan=True),
    ]
    out += [
        St.median(v),
        St.median(d, -1),
        St.median(u, -3),
        St.median(specials(v, payload=True), omitnan=True),
    ]
    out += [St.prctile(d.reshape(-1), [0, 1, 50, 99.9, 100]), St.prctile(v, 37.5, -2)]
    out += [St.prctile(np.array([-np.inf, np.inf, 1.0]), [10, 50, 90])]
    out += [St.mode(u), St.mode(u, -1, with_count=True)]
    return out


def _shift_problem():
    """a shift assignment of 6 tiles: 9 pair rows, bounds on the same differences"""
    n = 6
    pairs = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 2), (1, 3), (2, 4), (3, 5)]
    C = np.zeros((len(pairs), n))
    for i, (a, b) in enumerate(pairs):
        C[i, a], C[i, b] = 1.0, -1.0
    d = np.array([1.5, -2.25, 0.75, 3.125, -1.0, -0.5, -1.25, 4.0, 2.0])
    return C, d, np.vstack([C, -C]), np.full(2 * len(pairs), 10.0)


@case("lsq.assignment")
def _():
    """the assignment without equality rows whose bounds hold: sparse LU and sparse products"""
    from petakit5d.matlabcompat import lsq

    C, d, A, b = _shift_problem()
    out = [lsq.assignment(C, d, A, b, np.zeros((0, C.shape[1])), np.zeros(0))]
    out += [lsq.near_round_boundary(np.array([0.5, 1.49999, 2.0]), 1e-4)]
    out += [lsq.near_half_integer(np.array([0.5, 1.2, -2.5000001]), 1e-6)]
    return out


@case(
    "lsq.lsqlin",
    recorded="LAPACK: the exact optimum (equality rows, or a bound the least squares break) takes "
    "null spaces by a dense SVD, whose last bits follow the BLAS's kernels for the CPU; the shifts "
    "rounded from it differ only near a rounding boundary, which is flagged",
)
def _():
    from petakit5d.matlabcompat import lsq

    C, d, A, b = _shift_problem()
    Aeq = np.zeros((1, C.shape[1]))
    Aeq[0, 0], Aeq[0, 5] = 1.0, -1.0
    beq = np.array([0.25])
    out = [lsq.lsqlin(C, d, A, b, Aeq, beq), lsq.assignment(C, d, A, b, Aeq, beq)]
    out.append(np.asarray(lsq.solve(C, d, A, b, Aeq, beq).x))
    tight = np.full(A.shape[0], 1.0)
    return out + [lsq.assignment(C, d, A, tight, np.zeros((0, C.shape[1])), np.zeros(0))]


@case("neighbors_and_graph")
def _():
    from petakit5d.matlabcompat import graph as G
    from petakit5d.matlabcompat import neighbors as N

    X = double("kx", (300, 3), 10.0)
    Y = double("ky", (40, 3), 10.0)
    out = [N.knnsearch(X, Y, 3), N.knnsearch(X, X[:5], 1)]
    az, el, r = N.cart2sph(X[:, 0], X[:, 1], X[:, 2])
    out += [az, el, r, N.sph2cart(az, el, r), N.sph2cart(np.inf, 0.0, np.inf)]
    g = G.graph([0, 1, 2, 3, 0, 4], [1, 2, 3, 0, 4, 5], [1.5, 2.0, 0.5, 3.0, 1.0, 2.5])
    t, pred = G.minspantree(g)
    return out + [t.edges, t.weights, np.asarray(pred), G.dfsearch(g, 0), G.neighbors(g, 0)]


@case("neighbors.rangesearch")
def _():
    """rangesearch as the detection counts neighbours (points in a volume) and as tests/detection
    calls it (a close pair); a lattice, whose distances are exact, with r on them (d <= r, ties
    in row order) and duplicates at r = 0; queries apart from the points (some with none in
    range), 2-D points, no points, one query as a vector, coordinates near 1e6, signed zeros"""
    from petakit5d.matlabcompat import neighbors as N

    P = double("rsp", (60, 3), 2.0)
    P[10] = P[11] + 0.25
    out = [N.rangesearch(P, P, 1.5)]
    V = uniform("rsv", (300, 3)) * np.array([64.0, 64.0, 24.0])
    out.append(N.rangesearch(V, V, 4.25))
    z, y, x = np.meshgrid(np.arange(3.0), np.arange(4.0), np.arange(5.0), indexing="ij")
    L = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
    L = np.concatenate([L, L[[7, 30]]])
    out += [N.rangesearch(L, L, r) for r in (0.0, 1.0, math.sqrt(2.0), math.sqrt(3.0), 2.0)]
    out.append(N.rangesearch(L, uniform("rsq", (25, 3), -2.0, 6.0), 1.25))
    out.append(N.rangesearch(double("rs2", (80, 2), 1.0), double("rs2q", (9, 2), 1.0), 0.75))
    out += [N.rangesearch(np.zeros((0, 3)), L[:3], 1.0), N.rangesearch(L, L[17], 1.5)]
    B = 1e6 + uniform("rsb", (50, 3), 0.0, 2.0)
    out.append(N.rangesearch(B, B, 0.5))
    S = np.array([[0.0, -0.0, 0.0], [-0.0, 0.0, -0.0], [1.0, -0.0, 0.0], [0.0, 0.0, -1.0]])
    out.append(N.rangesearch(S, S, 1.0))
    return out


def _sigma_r(name: str, n: int) -> np.ndarray:
    """residual sigmas as the bad-point filter sees them (tests/detection's: |normal| 4 + 3)"""
    return np.abs(gen(name).normal(n)) * 4.0 + 3.0


@case("signal.histcounts")
def _():
    """histcounts as the bad-point filter calls it (residual sigmas, 60 bins and its coarsening
    by ceil(0.9 nb); the oracle cases' 341 and 17); integers on the edges (10, 7, 2 and 1 bins:
    widths of 10 and 100, powers of ten), negative data, signed zeros at the minimum, a narrow
    range near 1e6, values near 1e-300, two values"""
    from petakit5d.matlabcompat import signal as Sg

    sr = _sigma_r("hcs", 800)
    out = []
    nb = 60
    for _ in range(6):
        out.append(Sg.histcounts(sr, nb))
        nb = math.ceil(0.9 * nb)
    out += [Sg.histcounts(sr, n) for n in (341, 17)] + [Sg.histcounts_edges(sr, 45)]
    iv = np.concatenate([np.arange(101.0), np.arange(0.0, 101.0, 10.0)])
    out += [Sg.histcounts(iv, n) for n in (10, 7, 2, 1)]
    out.append(Sg.histcounts(double("hcn", (500,), 50.0) - 20.0, 13))
    z = np.array([-0.0, 0.0, 1.0, 2.5, 3.0, 0.5, -0.0])
    out += [Sg.histcounts(z, 4), Sg.histcounts_edges(z, 3)]
    out.append(Sg.histcounts(1e6 + uniform("hcb", (300,), 0.0, 1e-3), 25))
    out.append(Sg.histcounts(uniform("hct", (200,), 1.0, 2.0) * 1e-300, 9))
    return out + [Sg.histcounts(np.array([3.0, 7.5]), 3)]


@case("signal.findpeaks")
def _():
    """findpeaks as the bad-point filter calls it (histogram counts, prominence 5) and without a
    prominence: flat tops (one peak, at the first sample), plateaus that rise or reach the end
    (none), prominences at the limit and below it, negative values and -0, vectors of 0 to 3
    samples, normal values"""
    from petakit5d.matlabcompat import signal as Sg

    out = []
    for n in (60, 37):
        cnt, _ = Sg.histcounts(_sigma_r("fps", 800), n)
        out += [Sg.findpeaks(cnt, 5.0), Sg.findpeaks(cnt)]
    y = np.array([0, 1, 1, 1, 0, 2, 3, 3, 2, 5, 5, 6, 4, 4, 4, 7, 7, 3, 9, 9], dtype=np.float64)
    out += [Sg.findpeaks(y), Sg.findpeaks(y, 2.0), Sg.findpeaks(y, 3.0)]
    out += [Sg.findpeaks(np.array(v), 5.0) for v in ([0.0, 5.0, 0.0], [1.0, 6.0, 2.0])]
    out += [Sg.findpeaks(np.array([0.0, 4.0, 1.0, 4.0, 0.0]), 3.0)]
    out += [Sg.findpeaks(np.array([-3.0, -0.0, -2.0, -1.0, -5.0, -0.5, -0.5, -4.0]))]
    out += [Sg.findpeaks(np.arange(n, dtype=np.float64)) for n in (0, 1, 2)]
    out += [Sg.findpeaks(np.array([1.0, 3.0, 2.0]))]
    v = double("fpn", (500,), 1.0)
    return out + [Sg.findpeaks(v), Sg.findpeaks(v, 0.5)]


@case("signal.medfilt1")
def _():
    """medfilt1 as the bad-point filter smooths histogram counts (zeros past the ends), ties, 1
    to 3 samples, normal values with the special values apart (no window holds two NaNs or two
    zeros: the median's pick among equal values would be the sort's)"""
    from petakit5d.matlabcompat import signal as Sg

    cnt, _ = Sg.histcounts(_sigma_r("mfs", 800), 60)
    out = [Sg.medfilt1(cnt), Sg.medfilt1(np.array([4.0, 4.0, 1.0, 4.0, 4.0, 7.0, 7.0]))]
    out += [Sg.medfilt1(np.array(v)) for v in ([5.0], [5.0, 2.0], [1.0, 9.0, 3.0])]
    v = double("mfn", (61,), 10.0)
    tiny = np.finfo(np.float64).smallest_subnormal * 3
    put = (NAN, _payload_nan(np.float64), np.inf, -np.inf, -0.0, tiny, -tiny)
    v[[5, 12, 20, 21, 30, 40, 50]] = put
    v[29], v[31] = -1.0, 1.0  # -0 the middle of its window
    return out + [Sg.medfilt1(v), Sg.medfilt1(double("mfl", (1000,), 1.0))]


# -- special values through weighted sums ------------------------------------------------------
# Inf and NaN in float data meet zero weights and opposite infinities in the interpolations and
# filters: invalid operations, whose NaN has the sign bit set on x86 and clear on ARM (the
# kernels give x86's everywhere).


@case("specials.deskew")
def _():
    from petakit5d.kernels import deskew as DS

    out = []
    for shape, dz, reverse in (((20, 17, 64), 0.5, False), ((37, 45, 29), 0.2, True)):
        g = _geometry(shape, dz, reverse)
        raw = np.abs(single("dsx", shape))
        box = (0, 0, 0, *g.out_size)
        for vol in (specials(raw, every=41), payload_only(raw, every=41)):
            if g.do_interp:
                xs, st = g.x_step, g.int_stepsize
                out.append(
                    DS.skewed_space_interp_volume_deskew_rotate_warp(
                        vol, xs, st, reverse, g.tmat, box, False
                    )
                )
                out.append(DS.skewed_space_interp_defined_stepsize(vol, xs, st, reverse, False))
            else:
                out.append(DS.volume_deskew_rotate_warp(vol, g.tmat, box, False))
    return out


@case("specials.warp_and_resize")
def _():
    from petakit5d.matlabcompat import resize as R
    from petakit5d.matlabcompat import warp as W

    T = np.array(
        [[0.875, 0.0625, 0.0, 0.0], [-0.125, 1.0625, 0.25, 0.0], [0.0, 0.03125, 0.9375, 0.0],
         [1.5, -2.25, 0.75, 1.0]]
    )  # fmt: skip
    out = []
    for put in (specials, payload_only):
        v = put(single("sx", (9, 21, 17)), every=37)
        d = put(double("sxd", (9, 21, 17)), every=37)
        out += [W.imwarp(x, T, "linear", output_size=(8, 22, 18)) for x in (v, d)]
        out += [R.imresize3(x, (5, 30, 9), m) for x in (v, d) for m in ("linear", "cubic")]
        out += [R.imresize(v[0], 1.7, "bicubic"), R.imresize(d[0], 0.6, "bilinear")]
    return out


@case("specials.warp_cubic")
def _():
    """imwarp's cubic interpolation where Inf meets zero and negative weights"""
    from petakit5d.matlabcompat import warp as W

    T = np.array(
        [[0.875, 0.0625, 0.0, 0.0], [-0.125, 1.0625, 0.25, 0.0], [0.0, 0.03125, 0.9375, 0.0],
         [1.5, -2.25, 0.75, 1.0]]
    )  # fmt: skip
    out = []
    for put in (specials, payload_only):
        v = put(single("sx", (9, 21, 17)), every=37)
        d = put(double("sxd", (9, 21, 17)), every=37)
        out += [W.imwarp(x, T, "cubic", output_size=(8, 22, 18)) for x in (v, d)]
    return out


@case("specials.filters")
def _():
    from petakit5d.matlabcompat import filters as F

    k3 = np.abs(double("sfk", (3, 5, 4)))
    out = []
    for put in (specials, payload_only):
        v = put(single("sf", (6, 23, 31)), every=53)
        d = put(double("sfd", (6, 23, 31)), every=53)
        out += [
            F.imgaussfilt3(v, 1.5),
            F.imgaussfilt(d[0], 1.2),
            F.convn(v, k3.astype(np.float32), "same"),
            F.convn(d, k3, "full"),
            F.imfilter(v, k3, "symmetric"),
        ]
    return out


@case("specials.blending")
def _():
    from petakit5d.kernels import stitch as S

    w = uniform("sbw", (3, 4, 21, 19), 0.0, 2.0, np.float32)
    w.reshape(-1)[[5, 100, 400]] = [np.inf, 0.0, np.inf]  # not at the payload NaN's voxel
    out = []
    for put in (specials, payload_only):
        t = put(np.abs(single("sb", (3, 4, 21, 19))), every=61)
        out += [S.feather_blend(t, w)]
        out += [S.tile_blend(t, None, m) for m in ("none", "max", "mean", "median")]
    return out


@case("specials.integral_image")
def _():
    """integral_image_3d's differences of running sums: Inf - Inf"""
    from petakit5d.kernels import stitch as S

    raw = np.abs(single("sb", (3, 4, 21, 19)))
    return [S.integral_image(put(raw, every=61)[0], (2, 3, 4)) for put in (specials, payload_only)]


def infinities(name: str, shape, dtype, frac: float = 0.1) -> np.ndarray:
    """float data with infinities (no NaN) at a fraction frac of the elements, half of them
    negative: every NaN of an output is then an invalid operation's"""
    a = double(name, shape).astype(dtype)
    u = gen(name + "/inf").uniform(shape)
    a[u < frac / 2] = np.inf
    a[(u >= frac / 2) & (u < frac)] = -np.inf
    return np.ascontiguousarray(a)


def meeting_nans(dtype) -> np.ndarray:
    """a (7, 8, 9) volume whose NaNs of other bits (a payload, MATLAB's) sit side by side
    along each axis and at the ends, and two infinities"""
    v = (np.arange(7 * 8 * 9).reshape(7, 8, 9) % 7 + 1).astype(dtype)
    p, q = _payload_nan(dtype), dtype(NAN)
    v[3, 4, 4], v[3, 4, 5] = p, q  # along x
    v[1, 2, 6], v[2, 2, 6] = q, p  # along z
    v[5, 0, 1], v[5, 2, 1] = p, q  # at the start of y
    v[6, 7, 8], v[4, 7, 8] = q, p  # at the end of z
    v[3, 6, 2], v[3, 6, 3] = np.inf, -np.inf
    return v


def _shift(dx: float, dy: float, dz: float) -> np.ndarray:
    T = np.eye(4)
    T[3, :3] = [dx, dy, dz]
    return T


@case("specials.meeting_nans")
def _():
    """NaNs of other bits that meet in the integral image's sums, on purpose: the kernel computes
    an input with NaN in MATLAB's order itself (the first NaN operand as written), so these bits
    don't depend on the compiler"""
    from petakit5d.kernels import stitch as S

    v = meeting_nans(np.float32)
    return [S.integral_image(v, t) for t in ((2, 3, 4), (3, 1, 1), (1, 2, 1), (7, 8, 9))]


@case("specials.infinities")
def _():
    """float data with many infinities through the integral image, imwarp at whole and half
    coordinates (Inf * 0, Inf - Inf) and resize's cubic kernel and pow where |x|^3 overflows"""
    from petakit5d._ext import load
    from petakit5d.kernels import stitch as S
    from petakit5d.matlabcompat import resize as R
    from petakit5d.matlabcompat import warp as W

    a = infinities("sii", (6, 20, 17), np.float32)
    out = [S.integral_image(a, t) for t in ((2, 3, 4), (1, 1, 1), (6, 20, 17))]
    for dt in (np.float32, np.float64):
        v = infinities("swi", (9, 21, 17), dt)
        for T in (np.eye(4), _shift(-0.5, 0, 0.25)):
            for m in ("cubic", "linear"):
                out.append(W.imwarp(v, T, m, output_size=(9, 21, 17)))
    x = np.array([np.inf, -np.inf, 1e120, -1e200, 1e300, 0.5, 1.5, NAN, _payload_nan(np.float64)])
    out.append(R.cubic(x))
    x = np.array([1.0, -1.0, -1.0, -2.5, -0.0, -np.inf])
    for y in (np.inf, -np.inf, 0.5, -2.5):
        o = np.empty_like(x)
        load("resize").fdlibm_pow(x, y, o)
        out.append(o)
    return out


@case("specials.filters_infinities")
def _():
    """the linear filters' multiply-add chains where Inf meets a zero tap and the other
    infinity"""
    from petakit5d.matlabcompat import filters as F

    out = []
    for dt in (np.float32, np.float64):
        v = infinities("sfi", (6, 23, 31), dt)
        k = np.abs(double("sfik", (3, 5, 4)))
        k[1, 2, 1] = 0.0  # a zero tap
        out += [F.convn(v, k.astype(dt), "same"), F.imfilter(v, k, "symmetric")]
        out += [F.imgaussfilt3(v, 1.5), F.conv2(v[0], np.array([[1.0, 0.0, -1.0]]), "same")]
    return out


@case("specials.conv_overlap")
def _():
    """convn and conv2 of kernels with Inf and NaN taps: sums over the overlap of the two arrays
    only (no Inf * 0 with the zero padding), with MATLAB's NaN rules written out, so NaNs of
    other bits meet here on purpose: the swapped roles of a kernel larger than the image (plain
    sums) and kernels with such taps (fused)"""
    from petakit5d.matlabcompat import filters as F

    out = []
    for dt in (np.float32, np.float64):
        small = np.floor(double("sco", (3, 4), 3.0)).astype(dt)
        small.reshape(-1)[[0, 5, 6, 9, 11]] = [dt(NAN), _payload_nan(dt), dt(np.inf), -np.inf, -0.0]
        big = double("scb", (5, 6)).astype(dt)
        big.reshape(-1)[[7, 20]] = [dt(NAN), _payload_nan(dt)]
        out += [F.conv2(small, big), F.conv2(small, big, "same"), F.convn(small, big)]
        img = double("sci", (6, 23, 31)).astype(dt)
        img.reshape(-1)[[100, 2000]] = [_payload_nan(dt), dt(NAN)]
        k = np.abs(double("sck", (3, 5, 4))).astype(dt)
        k[0, 0, 0], k[1, 2, 1], k[2, 4, 3] = np.inf, _payload_nan(dt), -np.inf
        out += [F.convn(img, k, "same"), F.convn(img, k, "full"), F.conv2(img[0], k[1], "full")]
    return out


@case("specials.deskew_infinities")
def _():
    """the deskew/rotate warp and the skewed space interpolation of single data with
    infinities: Inf * 0 at whole coordinates, Inf - Inf"""
    from petakit5d.kernels import deskew as DS

    out = []
    for shape, dz, reverse in (((20, 17, 64), 0.5, False), ((37, 45, 29), 0.2, True)):
        g = _geometry(shape, dz, reverse)
        vol = infinities("sdi", shape, np.float32)
        box = (0, 0, 0, *g.out_size)
        if g.do_interp:
            xs, st = g.x_step, g.int_stepsize
            out.append(
                DS.skewed_space_interp_volume_deskew_rotate_warp(
                    vol, xs, st, reverse, g.tmat, box, False
                )
            )
            out.append(DS.skewed_space_interp_defined_stepsize(vol, xs, st, reverse, False))
        else:
            out.append(DS.volume_deskew_rotate_warp(vol, g.tmat, box, False))
    return out


# -- FFT-based outputs (recorded, not pinned) ---------------------------------------------------


@case("fft.matlab_fft", recorded="FFT: fft, ifft, fftn, ifftn (single and double, real, complex)")
def _():
    from petakit5d.matlabcompat import fft as FT

    s = single("ft", (6, 23, 30))
    d = double("ftd", (6, 23, 30))
    c = cplx(s, single("fti", (6, 23, 30)))
    out = [FT.fft(s), FT.fft(d, None, -1), FT.fft(c, 40, -3), FT.ifft(FT.fft(s))]
    out += [FT.ifft(FT.fft(d), None, None, True), FT.fftn(s), FT.fftn(d, (8, 24, 31))]
    out += [FT.ifftn(FT.fftn(s)), FT.ifftn(FT.fftn(d), None, True), FT.fftn(c)]
    return out


@case("fft.half_spectra", recorded="FFT: the half spectra of FFT convolutions, OTFs")
def _():
    from petakit5d import fft as PF

    out = []
    # MATLAB's half spectra too where the FFTW add-on is installed (its digest is then another)
    backends = ["permissive"] + (["exact"] if "exact" in PF.available() else [])
    for be in backends:
        for dt, f in ((np.float32, single), (np.float64, double)):
            hs = PF.HalfSpectrum((6, 23, 30), dt, backend=be)
            s = hs.forward(f("hs", (6, 23, 30)))
            h = hs.otf(np.abs(f("hk", (3, 5, 4))))
            out += [s, h, hs.inverse(s.copy()), hs.convolve(f("hx", (6, 23, 30)), h, True)]
    return out


@case("filters.edgetaper", recorded="FFT: edgetaper (the PSF's and the image's FFTs)")
def _():
    from petakit5d.matlabcompat import filters as F

    psf = np.abs(double("etp", (5, 7)))
    return [F.edgetaper(single("eti", (37, 70)), psf), F.edgetaper(double("etd", (37, 70)), psf)]


@case("fft.half_multiply")
def _():
    """the product of two half spectra (MATLAB's complex product: no FFT) and their sizes"""
    from petakit5d import fft as PF
    from petakit5d._ext import load

    out = [load("fft").half_shape([6, 23, 30]), load("fft").half_shape([23, 30])]
    for dt, f in ((np.float32, single), (np.float64, double)):
        hs = PF.HalfSpectrum((6, 23, 30), dt, backend="permissive")
        shape = hs.spectrum_shape
        s = cplx(specials(f("hm1", shape), payload=True), f("hm2", shape))
        h = cplx(f("hm3", shape), specials(f("hm4", shape), every=53, payload=True))
        out += [hs.multiply(s.copy(), h), hs.multiply(s.copy(), h, conjugate=True)]
    return out


@case("decon.otf_and_rl", recorded="FFT: decon's OTF, RL and OMW steps (single)")
def _():
    from petakit5d.kernels import decon as DK
    from petakit5d.testing.synth import decon_psf

    psf = decon_psf((15, 15, 15), (2.0, 1.5))
    raw = u16("rl", (16, 33, 30), 3000)
    H = DK.otf(psf, raw.shape)
    y = raw.astype(np.float32)
    out = [H, DK.rlstep(y, y, H)]
    out.append(DK.omw(raw, H, H, 3, 100.0, None))
    out.append(DK.omw(y, H, H, 2, 0.0, y))
    kw = {"layout": "matlab", "backend": "permissive"}  # the mex's loops on pocketfft
    Hm = DK.otf(psf, raw.shape, **kw)
    out += [Hm, DK.rlstep(y, y, Hm, **kw), DK.omw(raw, Hm, Hm, 3, 100.0, None, **kw)]
    out.append(DK.omw(y, Hm, Hm, 2, 0.0, y, **kw))
    out.append(DK.fft1(single("f1", (31, 6)), backend="permissive"))
    return out


@case("specials.fft_and_decon", recorded="FFT: FFTs and an RL step of data with Inf and NaN")
def _():
    from petakit5d.kernels import decon as DK
    from petakit5d.matlabcompat import fft as FT
    from petakit5d.testing.synth import decon_psf

    H = DK.otf(decon_psf((15, 15, 15), (2.0, 1.5)), (16, 33, 30))
    out = []
    for put in (specials, payload_only):
        s = put(single("sft", (6, 23, 30)), every=211)
        y = np.abs(put(single("srl", (16, 33, 30), 100.0), every=401))
        out += [FT.fft(s), FT.fftn(s.astype(np.float64)), DK.rlstep(y, y, H)]
    return out


@case("xcorr.ncc", recorded="FFT: the normalized cross-correlations of the registration")
def _():
    from petakit5d.kernels import xcorr as X

    a = np.abs(single("xa", (12, 30, 26)))
    t = np.ascontiguousarray(a[2:10, 5:25, 4:22] + np.float32(3.0))
    out = list(X.normxcorr3_max_shift(t, a, (2, 4, 3), True))
    out += list(X.normxcorr3_max_shift(t, a, (2, 4, 3), False))
    out += list(X.normxcorr3_window(t, a, (1, 2, 2, 8, 20, 18), True))
    out += list(X.normxcorr3_max_shift(t.astype(np.float64), a.astype(np.float64), (2, 4, 3), True))
    return out


# -- what the cases cover ------------------------------------------------------------------------

# bindings that are not data kernels, with the reason
EXEMPT_BINDINGS = {
    "core.affinity_cpu_count": "machine information",
    "core.build_info": "build information",
    "core.fma_contracted": "a build check (tests/test_smoke.py)",
    "core.omp_max_threads": "machine information",
    "core.omp_num_procs": "machine information",
    "core.ops_team_cpus": "thread placement",
    "core.parallel_threads": "thread counts",
    "core.subnormals_flushed": "a build check (tests/test_smoke.py)",
    "deskew.simd": "build information",
    "deskew_scalar.simd": "build information",
    "filters.best_isa": "the CPU's instruction set",
    "flatfield.best_isa": "the CPU's instruction set",
    "io.simd": "build information",
    "mip.best_isa": "the CPU's instruction set",
    "resize.best_isa": "the CPU's instruction set",
    "unmix.best_isa": "the CPU's instruction set",
    "fft.backend_info": "build information",
    "fft.exact_loaded": "the FFTW add-on's state",
    "fft.load_exact": "loads the FFTW add-on",
    "fft.unload_exact": "unloads the FFTW add-on",
    "fft.plan_string": "the FFTW add-on's plans",
    "fft.line_alone": "MATLAB's grouping of lines (tests/fft)",
}
for _g in ("decon", "fft", "xcorr"):
    for _n in ("_fft_cached_plan", "_fft_plan_lengths", "_fft_state_id"):
        EXEMPT_BINDINGS[f"{_g}.{_n}"] = "the FFT plan cache"
EXEMPT_BINDINGS["fft._shared_state"] = "the FFT state shared by the modules"

# matlabcompat functions that compute no data or need none here, with the reason
EXEMPT_FUNCTIONS = {
    "lsq.free_group_offsets": "integer bookkeeping of the stitching groups (tests/stitch)",
}


def exact_only_bindings() -> list[str]:
    """the bindings that only the exact FFT backend (the FFTW add-on) calls, where it is not
    installed: MATLAB's half spectra"""
    from petakit5d import fft

    if "exact" in fft.available():
        return []
    return ["fft.half_forward", "fft.half_inverse"]


def binding_names() -> list[str]:
    """every binding (``group.name``) of the extension modules found here"""
    from petakit5d._ext import groups, load

    out = []
    for g in groups():
        m = load(g)
        for n in dir(m):
            if n.startswith("__") or n.isupper():
                continue
            out.append(f"{g}.{n}")
    return sorted(out)


# -- command line --------------------------------------------------------------------------------


def _leaf_difference(a: np.ndarray, b: np.ndarray) -> str | None:
    """None when a and b have the same bits, else how they differ: the values that differ,
    those that are NaN or Inf on one side, and for the others the largest difference in ulps of
    the largest |value| of a (the measure of an FFT's rounding) and relative to it"""
    if a.dtype != b.dtype or a.shape != b.shape:
        return f"{a.dtype}{a.shape} against {b.dtype}{b.shape}"
    if a.tobytes() == b.tobytes():
        return None
    if a.dtype.kind not in "fc":
        return f"{int(np.count_nonzero(a != b))} of {a.size} values differ"
    ra = np.ascontiguousarray(a).reshape(-1).view(a.real.dtype)
    rb = np.ascontiguousarray(b).reshape(-1).view(b.real.dtype)
    w = ra.dtype.itemsize
    bad = (ra.view(np.uint8).reshape(-1, w) != rb.view(np.uint8).reshape(-1, w)).any(axis=1)
    special = bad & ~(np.isfinite(ra) & np.isfinite(rb))
    fin = bad & ~special
    text = f"{int(bad.sum())} of {ra.size} {ra.dtype} values differ"
    if special.any():
        both_nan = special & np.isnan(ra) & np.isnan(rb)
        text += f", {int(special.sum())} NaN or Inf ({int(both_nan.sum())} NaN on both sides)"
    if fin.any():
        top = float(np.max(np.abs(ra[np.isfinite(ra)]))) or 1.0
        d = float(np.max(np.abs(ra[fin].astype(np.float64) - rb[fin].astype(np.float64))))
        ulp = float(np.spacing(ra.dtype.type(top)))
        text += f", by at most {d / ulp:.3g} ulp of the largest |value| ({d / top:.3g} of it)"
    return text


def compare(fa: str, fb: str) -> None:
    """the differences of two --save files, per case and per differing output"""
    A, B = np.load(fa), np.load(fb)
    names = sorted({k.split("/", 1)[0] for k in A.files} | {k.split("/", 1)[0] for k in B.files})
    for case_name in names:
        keys = sorted(k for k in set(A.files) | set(B.files) if k.split("/", 1)[0] == case_name)
        lines = []
        for k in keys:
            if k not in A.files or k not in B.files:
                lines.append(f"    {k}: only in {'the first' if k in A.files else 'the second'}")
                continue
            d = _leaf_difference(A[k], B[k])
            if d is not None:
                lines.append(f"    {k}: {d}")
        tag = "REC" if case_name in RECORDED else "   "
        print(
            f"{tag} {case_name}: "
            + ("the same bits" if not lines else f"{len(lines)} outputs differ")
        )
        for line in lines:
            print(line)


def main(argv: list[str]) -> int:
    if argv[:1] == ["--compare"]:
        compare(argv[1], argv[2])
        return 0
    got = {n: digest(run(n)) for n in CASES}
    if argv[:1] == ["--write"]:
        import platform

        if platform.machine().lower() not in ("x86_64", "amd64") or sys.platform != "linux":
            print("the pins are made on Linux x86-64")
            return 1
        doc = {
            "made_with": "python tests/platform/_kernel_cases.py --write (Linux x86-64)",
            "pinned": {n: h for n, h in sorted(got.items()) if n not in RECORDED},
            "recorded": {n: h for n, h in sorted(got.items()) if n in RECORDED},
        }
        target = Path(argv[1]) if len(argv) > 1 else HASHES
        target.write_text(json.dumps(doc, indent=1) + "\n", newline="\n")
        return 0
    if argv[:1] == ["--save"]:
        arrays = {}
        for n in CASES:
            for p, a in leaves(run(n)):
                arrays[n + p] = a
        np.savez_compressed(argv[1], **arrays)
        return 0
    print(json.dumps(got, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
