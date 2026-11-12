"""Tests of petakit5d.matlabcompat.fft: MATLAB R2024b's fft, ifft, fftn and ifftn, bit for bit.

MATLAB tests whose intent is ported: tests/decon/DeconRFFTTest.m linesAreGroupedAsMatlabDoes
(MATLAB's groups of lines in both FFT modes) and resultsDoNotDependOnTheThreadCount (for the
transforms themselves); its other methods test the decon mex (the decon task's). The oracle packs
(MATLAB's own fft, ifft, fftn and ifftn: 1D lengths 2 to 1024, MATLAB's rules, the 3D sizes of the
demos, the tests and real data) run with MATLAB on the cluster; here the 1D and rule packs are
checked against digests of MATLAB's outputs, where our plans are MATLAB's (x86-64, AVX) and the C
library's cos and sin give FFTW the reference platform's twiddle factors.
"""

import functools
import hashlib
import math
import os
import platform

import numpy as np
import pytest

from petakit5d.matlabcompat import fft as F

# --- bit-stable inputs (the oracle packs' generator) -------------------------------------------


def splitmix(seed: int, n: int) -> np.ndarray:
    """n uint64 values of splitmix64 started at seed."""
    with np.errstate(over="ignore"):
        z = np.uint64(seed) + np.arange(1, n + 1, dtype=np.uint64) * np.uint64(0x9E3779B97F4A7C15)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return z ^ (z >> np.uint64(31))


def values(seed: int, n: int, kind: str) -> np.ndarray:
    """float64 values: 'u' uniform [-1, 1), 'img' integers 0..4095 (camera data), 'wide' about
    2^-20 .. 2^19."""
    u = (splitmix(seed, n) >> np.uint64(11)).astype(np.float64) * (1.0 / 9007199254740992.0)
    if kind == "u":
        return 2.0 * u - 1.0
    if kind == "img":
        return np.floor(u * 4096.0)
    e = (splitmix(seed + 77, n) % np.uint64(40)).astype(np.float64) - 20.0
    return (2.0 * u - 1.0) * np.exp2(e)


def make_array(seed, shape, cls, cplx, kind="u", hermitian=""):
    """A Fortran-ordered array; hermitian 'lines0' (along dimension 0) or 'nd' makes it exactly
    conjugate symmetric ((z + conj(mirror z)) / 2)."""
    dt = np.float32 if cls == "single" else np.float64
    n = int(np.prod(shape))
    re = values(seed, n, kind).astype(dt).reshape(shape, order="F")
    if not cplx:
        return np.asfortranarray(re)
    z = np.empty(shape, dtype=np.complex64 if cls == "single" else np.complex128, order="F")
    z.real = re
    z.imag = values(seed + 1, n, kind).astype(dt).reshape(shape, order="F")
    if hermitian:
        axes = (0,) if hermitian == "lines0" else tuple(range(z.ndim))
        m = z
        for ax in axes:
            m = np.roll(np.flip(m, axis=ax), 1, axis=ax)
        z = np.asfortranarray(((z + np.conj(m)) / dt(2)).astype(z.dtype))
    return z


def case_input(c):
    kind = c.get("kind", "u")
    cls = c["cls"]
    if cls == "uint16":
        v = values(c["seed"], int(np.prod(c["shape"])), "img").astype(np.uint16)
        return np.asfortranarray(v.reshape(c["shape"], order="F"))
    dt = np.float32 if cls == "single" else np.float64
    if kind == "zeros":
        return np.zeros(c["shape"], dtype=dt, order="F")
    if kind.startswith("zeros") and len(kind) == 7:
        # complex zeros with the given signs of the real and imaginary parts
        z = np.zeros(
            c["shape"], dtype=np.complex64 if cls == "single" else np.complex128, order="F"
        )
        z.real = -0.0 if kind[5] == "-" else 0.0
        z.imag = -0.0 if kind[6] == "-" else 0.0
        return z
    if kind == "sepsym":
        # a product of integer vectors, even-symmetric along dimensions 0 and 1: exact zeros in the
        # imaginary parts of fftn's self-conjugate planes, which its Hermitian fill copies
        shape = c["shape"]
        x = np.ones(shape, dtype=dt, order="F")
        for k, n in enumerate(shape):
            v = np.floor(values(c["seed"] + k, n, "u") * 4.0 + 5.0)
            if k < 2 and n > 1:
                v[1:] = v[1:] + v[1:][::-1]
            sh = [1] * len(shape)
            sh[k] = n
            x = x * v.astype(dt).reshape(sh)
        return np.asfortranarray(x)
    if kind in ("evensym", "evensym0"):
        x = make_array(c["seed"], c["shape"], cls, False, "u")
        y = x.copy(order="F")
        if x.shape[0] > 1:
            y[1:] = (x[1:] + x[1:][::-1]) / x.dtype.type(2)
        return np.asfortranarray(y + 0j) if kind == "evensym0" else y
    return make_array(c["seed"], c["shape"], cls, c["cplx"], kind, c.get("herm", ""))


# --- the packs ---------------------------------------------------------------------------------


def case(func, cls, cplx, shape, seed, **kw):
    return {"func": func, "cls": cls, "cplx": cplx, "shape": shape, "seed": seed, **kw}


def pack_1d():
    """1D lengths 2..1024 as n x 7 matrices along dimension 1 (lines 0-3 in a group, 4-6 alone),
    every kind, both precisions; rows along dimension 2 for every 7th length."""
    cases = []
    for n in range(2, 1025):
        for cls in ("single", "double"):
            k = ("u", "img", "wide")[n % 3]
            s = 1000 * n + (0 if cls == "single" else 500)
            cases.append(case("fft", cls, False, [n, 7], s, kind=k, dim=1))
            cases.append(case("fft", cls, True, [n, 7], s + 10, kind=k, dim=1))
            cases.append(case("ifft", cls, True, [n, 7], s + 20, kind=k, dim=1))
            cases.append(case("ifft", cls, True, [n, 7], s + 30, kind=k, dim=1, herm="lines0"))
            if n % 7 == 0:
                cases.append(case("fft", cls, False, [7, n], s + 40, kind=k, dim=2))
                cases.append(case("ifft", cls, True, [7, n], s + 50, kind=k, dim=2))
    return cases


def pack_edge():
    """MATLAB's rules: real results, the symmetry test, padding and truncation, length-1
    dimensions, n and siz, 'symmetric', the split threshold, integer input."""
    c = []
    for cls in ("single", "double"):
        b = 0 if cls == "single" else 100000
        for n in (1, 2, 3, 4, 6, 8, 64, 97, 100):
            s = b + n
            c.append(case("fft", cls, False, [n, 1], s, kind="u"))
            c.append(case("fft", cls, False, [1, n], s + 1, kind="u"))
            c.append(case("fft", cls, True, [n, 1], s + 2, kind="u"))
            c.append(case("ifft", cls, True, [n, 1], s + 3, kind="u"))
            c.append(case("ifft", cls, True, [n, 1], s + 4, kind="u", herm="lines0"))
            c.append(case("ifft", cls, False, [n, 1], s + 5, kind="u"))
            c.append(case("ifft", cls, True, [n, 1], s + 6, kind="u", sym=True))
            c.append(case("fft", cls, False, [n, 1], s + 7, kind="zeros"))
            c.append(case("fft", cls, False, [n, 1], s + 8, kind="evensym"))
            c.append(case("fft", cls, True, [n, 1], s + 9, kind="evensym0"))
            c.append(case("ifft", cls, False, [n, 1], s + 10, kind="evensym"))
        sh = [12, 10, 8]
        for d in (1, 2, 3):
            for m in (1, 5, 13):
                s = b + d * 10 + m
                c.append(case("fft", cls, False, sh, s + 200, n=m, dim=d))
                c.append(case("fft", cls, True, sh, s + 300, n=m, dim=d))
                c.append(case("ifft", cls, True, sh, s + 400, n=m, dim=d))
                c.append(case("ifft", cls, True, sh, s + 500, n=m, dim=d, sym=True))
            c.append(case("fft", cls, False, sh, b + 600 + d, dim=d))
        c.append(case("fft", cls, False, sh, b + 610, dim=4))
        c.append(case("fft", cls, False, [12, 10], b + 611, n=4, dim=3))
        c.append(case("ifft", cls, True, [3, 1], b + 620, herm="lines0", n=4))
        c.append(case("ifft", cls, True, [8, 3], b + 621, herm="lines0", n=6))
        c.append(case("ifft", cls, True, [1, 1], b + 622, herm="lines0", n=2))
        for siz in (
            [9, 7, 5],
            [16, 7, 8],
            [12, 10, 8, 2],
            [1, 10, 8],
            [12, 1, 8],
            [13, 11, 9],
            [12, 10],
        ):
            s = b + sum(siz) * 7 + len(siz)
            c.append(case("fftn", cls, False, sh, s + 700, siz=siz))
            c.append(case("fftn", cls, True, sh, s + 800, siz=siz))
            c.append(case("ifftn", cls, True, sh, s + 900, siz=siz))
            c.append(case("ifftn", cls, True, sh, s + 1000, siz=siz, herm="nd"))
            c.append(case("ifftn", cls, True, sh, s + 1100, siz=siz, sym=True))
        c.append(case("ifftn", cls, True, sh, b + 1200, herm="lines0"))
        c.append(case("ifftn", cls, True, sh, b + 1201, herm="nd"))
        c.append(case("ifftn", cls, False, sh, b + 1202))
        c.append(case("ifft", cls, True, sh, b + 1203, herm="nd", dim=3))
        for shape in (
            [2, 2, 2],
            [1, 7, 5],
            [7, 1, 5],
            [1, 1, 5],
            [6, 5, 4, 3],
            [15, 12],
            [2, 64, 33],
            [1, 1],
        ):
            s = b + int(np.prod(shape)) * 3 + len(shape)
            c.append(case("fftn", cls, False, shape, s + 1300))
            c.append(case("fftn", cls, True, shape, s + 1400))
            c.append(case("ifftn", cls, True, shape, s + 1500, herm="nd"))
            c.append(case("ifftn", cls, True, shape, s + 1600))
            c.append(case("fft", cls, False, shape, s + 1700))
        for N in (219, 220, 221, 222, 300, 1000, 5001):
            c.append(case("fft", cls, False, [64, N], b + 2000 + N))
            c.append(case("fft", cls, True, [64, N], b + 3000 + N))
            c.append(case("ifft", cls, True, [64, N], b + 4000 + N, herm="lines0"))
        c.append(case("fft", cls, False, [60, 222], b + 5000, n=64))
        c.append(case("fft", cls, False, [70, 222], b + 5001, n=64))
        for shape in ([64, 6, 37], [6, 64, 37], [6, 37, 64], [64, 64, 64], [6, 64, 1000]):
            for d in (1, 2, 3):
                s = b + int(np.prod(shape)) + d
                c.append(case("fft", cls, False, shape, s + 6000, dim=d))
                c.append(case("fft", cls, True, shape, s + 7000, dim=d))
                c.append(case("ifft", cls, True, shape, s + 8000, dim=d))
            c.append(case("fftn", cls, False, shape, s + 9000))
            c.append(case("ifftn", cls, True, shape, s + 9500, herm="nd"))
    for cls in ("single", "double"):
        b = 0 if cls == "single" else 100000
        for shape in ([12, 10, 8], [9, 7, 5], [6, 5, 4, 3], [15, 12], [1, 7, 5], [8, 1], [7, 6]):
            s = b + int(np.prod(shape)) + len(shape)
            c.append(case("ifftn", cls, True, shape, s + 12000, sym=True))
            c.append(case("ifftn", cls, True, shape, s + 13000, herm="nd", sym=True))
        for shape in ([4, 4, 4], [4, 6], [8, 1], [5, 3, 2]):
            for k in ("zeros--", "zeros+-", "zeros-+"):
                s = b + int(np.prod(shape)) + len(shape)
                c.append(case("ifftn", cls, True, shape, s + 14000, kind=k))
                c.append(case("ifftn", cls, True, shape, s + 15000, kind=k, sym=True))
                c.append(case("ifft", cls, True, shape, s + 16000, kind=k))
                c.append(case("fftn", cls, True, shape, s + 17000, kind=k))
        # vectors: MATLAB's last dimension of a column is the second (length 1)
        for shape in ([8, 1], [64, 1], [97, 1], [1, 8], [1, 97]):
            s = b + 19000 + 10 * shape[0] + shape[1]
            c.append(case("ifftn", cls, True, shape, s, herm="nd"))
            c.append(case("ifftn", cls, True, shape, s + 1, sym=True))
            c.append(case("fftn", cls, False, shape, s + 2))
            c.append(case("ifft", cls, True, shape, s + 3, herm="nd"))
        # zeros that fftn's Hermitian fill copies (MATLAB's conjugate there is 0 - im)
        for shape in ([16, 16, 8], [12, 10, 6], [9, 8, 4]):
            c.append(case("fftn", cls, False, shape, b + 20000 + shape[0], kind="sepsym"))
        c.append(case("fft", cls, False, [1, 1], b + 18000, n=3))
        c.append(case("ifft", cls, True, [1, 1], b + 18001, n=4))
        c.append(case("fftn", cls, False, [1, 1], b + 18002, siz=[2, 3]))
    c.append(case("fft", "uint16", False, [37, 5], 99001, kind="img"))
    c.append(case("fftn", "uint16", False, [12, 10, 8], 99002, kind="img"))
    c.append(case("ifftn", "uint16", False, [12, 10, 8], 99003, kind="img"))
    return c


def pack_3d():
    """fftn and ifftn at the sizes of the demos (images, PSFs), the decon tests and a real decon
    batch."""
    sizes = [
        [64, 50, 40],
        [47, 45, 33],
        [2, 64, 33],
        [64, 230],
        [64, 210],
        [128, 512, 81],
        [160, 512, 101],
        [201, 201, 214],
        [201, 201, 315],
        [551, 551, 231],
        [571, 501, 501],
        [1024, 1024, 50],
        [2304, 2304],
        [1800, 512, 500],
        [512, 512, 500],
        [3082, 11530, 4],
        [3072, 11520, 4],
    ]
    c = []
    for shape in sizes:
        k = int(np.prod(shape)) % 100000
        big = np.prod(shape) > 2e8
        for cls in ("single",) if big else ("single", "double"):
            b = k * 10 + (0 if cls == "single" else 1)
            c.append(case("fftn", cls, False, shape, b, kind="img"))
            c.append(case("ifftn", cls, True, shape, b + 2, herm="nd"))
            if not big:
                c.append(case("fftn", cls, True, shape, b + 4))
                c.append(case("ifftn", cls, True, shape, b + 6))
    return c


PACKS = {"1d": pack_1d, "edge": pack_edge, "3d": pack_3d}


def run_case(c, x, parallel=True, threads=None):
    """The case's call (the manifest's dim is MATLAB's, 1-based)."""
    kw = {"matlab_parallel": parallel, "threads": threads}
    dim = c.get("dim")
    d0 = None if dim is None else dim - 1
    f = c["func"]
    if f == "fft":
        return F.fft(x, c.get("n"), d0, **kw)
    if f == "ifft":
        return F.ifft(x, c.get("n"), d0, c.get("sym", False), **kw)
    if f == "fftn":
        return F.fftn(x, c.get("siz"), **kw)
    return F.ifftn(x, c.get("siz"), c.get("sym", False), **kw)


def matlab_shape(shape):
    s = list(shape)
    while len(s) > 2 and s[-1] == 1:
        s.pop()
    while len(s) < 2:
        s.append(1)
    return s


def case_record(cid, y):
    """What the digests cover: class, realness, MATLAB's size and the bytes of a result."""
    cls = {
        np.float32: "single",
        np.complex64: "single",
        np.float64: "double",
        np.complex128: "double",
    }
    raw = np.asfortranarray(y).tobytes(order="F")
    dims = " ".join(str(v) for v in matlab_shape(y.shape))
    real = int(not np.iscomplexobj(y))
    return f"{cid} {cls[y.dtype.type]} {real} {dims} {hashlib.sha256(raw).hexdigest()}"


def pack_digest(pack, parallel, threads=None, step=1):
    """sha256 over the records of a pack's cases ('<id> ERROR' where MATLAB raises); step > 1 takes
    every step-th case."""
    lines = []
    for i, c in enumerate(PACKS[pack]()):
        if i % step:
            continue
        cid = f"{pack}{i:05d}"
        x = case_input(c)
        try:
            y = run_case(c, x, parallel, threads)
        except (ValueError, TypeError):
            lines.append(f"{cid} ERROR")
            continue
        lines.append(case_record(cid, y))
    return hashlib.sha256("\n".join(lines).encode()).hexdigest(), len(lines)


# digests of MATLAB R2024b Update 5's outputs (oracle jobs on the cluster), per pack and FFT mode
MATLAB_DIGESTS = {
    ("1d", "multi"): ("064aaea2bd05091adc79162040964ccb0c72d8ff9d6b4fa290c4169fef1de50b", 8768),
    ("1d", "single"): ("064aaea2bd05091adc79162040964ccb0c72d8ff9d6b4fa290c4169fef1de50b", 8768),
    ("edge", "multi"): ("83eb081d0ceba789dcc5832bdabcdfcde6f6a00e698e13daaef209165cf6ed4b", 779),
    ("edge", "single"): ("924c3366505a159ec9d14670c80348ee9fecd1952db91db09d6247ef2ee08301", 779),
}

# sha256 of the plans of MATLAB's libmwfftw3f / libmwfftw3 (whitespace collapsed, one line per
# plan) for every n in 1..1024, both precisions, kinds r2c c2c-1 c2c+1 c2r, 1 and 4 lines per
# execute
MATLAB_PLANS_DIGEST = "238f71138c3fc44f1e4ee2f59f906394ef30577c29d25fe103989a91ccb48bd9"

# FFTW makes its twiddle factors with the C library's cos and sin (kernel/trig.c, in double for
# both precisions), so their last bits, and with them the FFT's, follow the C library: its version
# and the code path it picks for the CPU. MATLAB's own FFTW calls the same functions. The digests
# above were recorded where cos and sin give these values for FFTW's arguments (twiddle_fingerprint;
# glibc 2.28 on Intel x86-64); elsewhere MATLAB's results would differ from them too.
REFERENCE_TWIDDLE_FINGERPRINT = "5442ee1a3bc9ca08a71ef7c693dab89979f4788915be2b27e3525416adce53bb"
_K2PI = 6.2831853071795864769252867665590057683943388  # FFTW's K2PI, rounded to double as C does


def twiddle_arguments(n):
    """The arguments FFTW's real_cexp passes to cos and sin for 2 pi m / n, m = 0 .. n - 1: 4m and
    4n reduced to the first octant, then (K2PI * m) / n in double."""
    m = 4 * np.arange(n, dtype=np.int64)
    n4 = 4 * n
    m = np.where(m > n4 - m, n4 - m, m)
    m = np.where(m - n > 0, m - n, m)
    m = np.where(m > n - m, n - m, m)
    return (_K2PI * np.unique(m).astype(np.float64)) / float(n4)


def twiddle_fingerprint(nmax=2048, cos=math.cos, sin=math.sin):
    """sha256 of cos and sin at every argument of twiddle_arguments(n), n = 1 .. nmax (the lengths
    of the packs and of FFTW's sub-plans for them)."""
    h = hashlib.sha256()
    for n in range(1, nmax + 1):
        t = twiddle_arguments(n).tolist()
        h.update(np.array([cos(v) for v in t] + [sin(v) for v in t], dtype=np.float64).tobytes())
    return h.hexdigest()


@functools.cache
def digests_hold_here():
    """'' where the MATLAB digests apply, else why not (FFTW's plans or twiddle factors differ)."""
    if not matlab_platform():
        return "FFTW's plans differ from MATLAB's here (not x86-64 with AVX)"
    fp = twiddle_fingerprint()
    if fp != REFERENCE_TWIDDLE_FINGERPRINT:
        return (
            "the C library's cos and sin give FFTW other twiddle factors than where the digests "
            f"were recorded (fingerprint {fp[:16]}, not {REFERENCE_TWIDDLE_FINGERPRINT[:16]}): "
            "MATLAB's own FFTW would differ from them here too"
        )
    return ""


def plans_digest():
    ext = F._ext()
    lines = []
    for prec in (False, True):
        for kind in ("r2c", "c2c-1", "c2c+1", "c2r"):
            for hm in (1, 4):
                for n in range(1, 1025):
                    s = " ".join(ext.plan_string(prec, kind, n, hm).split())
                    lines.append(f"{'d' if prec else 'f'} {kind} hm={hm} n={n} | {s}")
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def matlab_platform():
    """Whether our plans can be MATLAB's here: x86-64 with AVX (FFTW picks its codelets from the
    CPU)."""
    if platform.machine().lower() not in ("x86_64", "amd64"):
        return False
    return ext_plan("r2c", 64, 4) == '(rdft2-r2hc-direct-64-x4 "r2cf_64")' and "avx" in ext_plan(
        "c2c-1", 64, 4
    )


def ext_plan(kind, n, hm, prec=False):
    return " ".join(F._ext().plan_string(prec, kind, n, hm).split())


needs_matlab_platform = pytest.mark.skipif(
    "not matlab_platform()", reason="plans differ from MATLAB here"
)


# --- MATLAB's driver rules ---------------------------------------------------------------------


def mark_alone(N, n, seg, grain):
    """decon_lucy_rfft_mex's mark_alone and stage_grain: the lines MATLAB transforms alone."""
    bits = int(n).bit_length()
    g = grain if grain > 0 and n * N * bits > 98304 else 0
    alone = np.zeros(N, dtype=bool)

    def mark(a, b):
        if g > 0 and b - a > g:
            mid = a + (b - a) // 2
            mark(a, mid)
            mark(mid, b)
            return
        s = a
        while s < b:
            e = b
            if seg > 0 and (s // seg + 1) * seg < b:
                e = (s // seg + 1) * seg
            alone[e - (e - s) % 4 : e] = True
            s = e

    mark(0, N)
    return alone


@pytest.mark.parametrize(
    "N,n,seg,grain",
    [
        (300, 64, 0, 32),
        (222, 64, 0, 32),
        (219, 64, 0, 32),
        (210, 64, 0, 32),
        (5001, 64, 0, 32),
        (222, 64, 6, 32),
        (2368, 6, 64, 32),
        (256000, 512, 512, 32),
        (262144, 500, 0, 32),
        (300, 64, 0, 0),
        (37 * 6, 64, 6, 0),
        (7, 1024, 0, 32),
        (1, 8, 0, 32),
        (4, 5, 2, 32),
    ],
)
def test_line_rule_matches_the_decon_mex(N, n, seg, grain):
    ext = F._ext()
    ours = np.array([ext.line_alone(j, N, n, seg, grain) for j in range(N)])
    np.testing.assert_array_equal(ours, mark_alone(N, n, seg, grain))


def halve(a, b, g):
    if b - a <= g:
        return [(a, b)]
    m = a + (b - a) // 2
    return halve(a, m, g) + halve(m, b, g)


def fft_by_pieces(X, pieces):
    """MATLAB's lines of a piece: groups of four from its start, the rest alone; each computed by a
    call that MATLAB (and we) run as one group or alone."""
    out = []
    for a, b in pieces:
        j = a
        while j + 4 <= b:
            out.append(F.fft(np.asfortranarray(X[:, j : j + 4]), dim=0))
            j += 4
        while j < b:
            out.append(F.fft(np.asfortranarray(X[:, j : j + 1]), dim=0))
            j += 1
    return np.concatenate(out, axis=1)


@pytest.mark.parametrize("N", [2000, 222, 210])
def test_lines_are_grouped_as_matlab_does(N):
    # DeconRFFTTest/linesAreGroupedAsMatlabDoes: MATLAB's parallel mode halves the lines into pieces
    # of at most 32 (if n * N * (bits of n) > 98304), its single-thread mode does not
    X = make_array(N, [64, N], "single", False)
    pieces = halve(0, N, 32) if 64 * N * 7 > 98304 else [(0, N)]
    ref_split = fft_by_pieces(X, pieces)
    ref_serial = fft_by_pieces(X, [(0, N)])
    if N == 2000 and np.array_equal(ref_split, ref_serial):
        pytest.skip("FFTW gives the same bits for groups of four and single lines here")
    assert F.fft(X, dim=0).tobytes() == ref_split.tobytes()
    assert F.fft(X, dim=0, matlab_parallel=False).tobytes() == ref_serial.tobytes()


@pytest.mark.parametrize(
    "case",
    [
        ("fft", [64, 300], False, None),
        ("fft", [6, 64, 37], True, 1),
        ("fftn", [64, 50, 40], False, None),
        ("fftn", [47, 45, 33], True, None),
        ("ifftn", [64, 50, 40], True, None),
        ("ifft", [37, 222], True, 1),
    ],
)
@pytest.mark.parametrize("cls", ["single", "double"])
def test_results_do_not_depend_on_the_thread_count(case, cls):
    func, shape, cplx, dim = case
    x = make_array(7, shape, cls, cplx)
    if func == "ifftn":
        x = make_array(7, shape, cls, True, hermitian="nd")
    f = getattr(F, func)
    kw = {} if dim is None else {"dim": dim}
    ref = f(x, threads=1, **kw)
    for t in (2, 3, 8):
        assert f(x, threads=t, **kw).tobytes() == ref.tobytes()


@pytest.mark.parametrize("cls,tol", [("single", 2e-6), ("double", 1e-14)])
def test_values_agree_with_numpy(cls, tol):
    x = make_array(3, [24, 18, 10], cls, True)
    r = make_array(4, [24, 18, 10], cls, False)
    ref = np.fft.fftn(x.astype(np.complex128))
    assert np.max(np.abs(F.fftn(x) - ref)) <= tol * np.max(np.abs(ref))
    ref = np.fft.fft(r.astype(np.float64), axis=1)
    assert np.max(np.abs(F.fft(r, dim=1) - ref)) <= tol * np.max(np.abs(ref))
    y = F.ifftn(F.fftn(r))
    assert y.dtype == r.dtype and np.max(np.abs(y - r)) <= 10 * tol
    z = F.ifft(x, 30, 0)
    assert z.shape == (30, 18, 10)
    assert np.max(np.abs(z - np.fft.ifft(x.astype(np.complex128), 30, axis=0))) <= tol


def test_classes_follow_matlab():
    assert F.fft(np.arange(8, dtype=np.uint16)).dtype == np.complex128
    # imaginary parts all zero
    assert F.fft(np.array([True, False, True, False])).dtype == np.float64
    assert F.fftn(np.zeros((3, 4), np.float32)).dtype == np.float32
    assert F.fft(np.ones((2, 1), np.float32)).dtype == np.float32
    assert F.fft(np.arange(3.0)).dtype == np.complex128
    assert F.fftn(make_array(1, [2, 2, 2], "double", False)).dtype == np.float64
    with pytest.raises(TypeError):
        F.fft(np.array(["a", "b"]))
    with pytest.raises(TypeError):
        F.fft(np.arange(4, dtype=np.int64))


def test_conjugate_symmetric_input_gives_a_real_result():
    X = F.fft(make_array(5, [16, 3], "single", False), dim=0)
    assert X.dtype == np.complex64
    assert F.ifft(X, dim=0).dtype == np.float32
    Y = X.copy(order="F")
    Y[1, 0] = complex(Y[1, 0].real, Y[1, 0].imag * (1 + np.finfo(np.float32).eps))
    assert F.ifft(Y, dim=0).dtype == np.complex64  # one line fails: every line is complex
    Y = X.copy(order="F")
    Y[0, 0] = complex(Y[0, 0].real, -0.0)
    assert F.ifft(Y, dim=0).dtype == np.float32  # -0 counts as 0
    Y[0, 0] = complex(Y[0, 0].real, np.finfo(np.float32).eps)
    assert F.ifft(Y, dim=0).dtype == np.complex64
    Y = X.copy(order="F")
    Y[1, 1] = Y[15, 1] = np.nan
    assert F.ifft(Y, dim=0).dtype == np.complex64  # NaN is not equal to itself
    # N-D: data symmetric along dimension 0 only is not
    A = make_array(6, [8, 6, 4], "double", False)
    assert F.ifftn(F.fft(A, dim=0)).dtype == np.complex128
    assert F.ifftn(F.fftn(A)).dtype == np.float64


def test_fftn_of_real_data_is_exactly_hermitian():
    for shape in ([9, 7, 5], [12, 10, 8], [1, 7, 6], [16, 1]):
        X = F.fftn(make_array(8, shape, "single", False, "img"))
        m = X
        for ax in range(X.ndim):
            m = np.roll(np.flip(m, axis=ax), 1, axis=ax)
        assert np.array_equal(X.real, m.real) and np.array_equal(X.imag, -m.imag)


def test_fftn_fill_conjugate_is_zero_minus_im():
    # MATLAB's fill writes 0 - im: a zero imaginary part copied from a mirror is +0, never -0
    for cls in ("single", "double"):
        x = case_input(case("fftn", cls, False, [16, 16, 8], 7, kind="sepsym"))
        X = F.fftn(x)
        assert np.iscomplexobj(X)
        flat = X.reshape(-1, order="F")
        idx = np.arange(X.size).reshape(X.shape, order="F")
        m = idx
        for ax in range(X.ndim):
            m = np.roll(np.flip(m, axis=ax), 1, axis=ax)
        mi = m.reshape(-1, order="F")
        filled = (mi < np.arange(X.size)) & (flat.imag[mi] == 0)
        assert filled.sum() > 0
        assert not np.signbit(flat.imag[filled]).any()


def test_symmetric_flag_is_the_fill_of_fftn():
    # ifftn(X, 'symmetric') makes the (padded) data conjugate symmetric as fftn's fill does, then
    # inverts it
    for siz in (None, [13, 11, 9], [9, 7, 5]):
        X = make_array(9, [12, 10, 8], "double", True)
        P = X if siz is None else np.asfortranarray(_pad(X, siz))
        H = _fill(P)
        ref = F.ifftn(H)
        assert ref.dtype == np.float64
        assert F.ifftn(X, siz, True).tobytes() == ref.tobytes()


def _pad(a, siz):
    out = np.zeros(siz, dtype=a.dtype, order="F")
    sl = tuple(slice(0, min(s, t)) for s, t in zip(a.shape, siz, strict=True))
    out[sl] = a[sl]
    return out


def _fill(a):
    """fftn's Hermitian fill: an element whose mirror has the lower linear index takes the mirror's
    conjugate, a self-conjugate element its real part."""
    idx = np.arange(a.size).reshape(a.shape, order="F")
    m = idx
    for ax in range(a.ndim):
        m = np.roll(np.flip(m, axis=ax), 1, axis=ax)
    flat = a.reshape(-1, order="F").copy()
    mi = m.reshape(-1, order="F")
    lower = mi < np.arange(a.size)
    src = flat[mi[lower]]
    # MATLAB's conjugate here is 0 - im (a zero imaginary part becomes +0)
    flat.real[lower] = src.real
    flat.imag[lower] = flat.real.dtype.type(0) - src.imag
    self_ = mi == np.arange(a.size)
    flat.imag[self_] = 0
    return np.asfortranarray(flat.reshape(a.shape, order="F"))


def test_padding_and_truncation_come_first():
    x = make_array(10, [12, 10, 8], "single", False)
    assert F.fft(x, 15, 1).tobytes() == F.fft(_pad(x, [12, 15, 8]), None, 1).tobytes()
    assert F.fft(x, 7, 2).tobytes() == F.fft(np.asfortranarray(x[:, :, :7]), None, 2).tobytes()
    assert F.fftn(x, [14, 9, 11]).tobytes() == F.fftn(_pad(x, [14, 9, 11])).tobytes()


def test_length_one_dimensions():
    x = make_array(11, [12, 10], "single", False)
    assert F.fft(x, dim=2).tobytes() == x.tobytes()  # no transform: the input
    assert F.fft(np.array([[5.0]]), 3).shape == (1, 3)  # a scalar is a row
    assert F.fftn(np.array([[2.0]])).tobytes() == np.array([[2.0]]).tobytes()
    with pytest.raises(ValueError):
        F.fftn(make_array(1, [4, 4, 4], "single", False), [4, 4])


def test_dimensions_beyond_the_input_are_length_1():
    # the copy for a transform along a dimension beyond the input's read its size past the end
    for k in range(40):
        shape = [3 + k % 5, 4 + k % 3, 2 + k % 4][: 1 + k % 3]
        x = make_array(30 + k, shape, "single" if k % 2 else "double", bool(k % 4 == 1))
        nd = len(shape)
        assert F.fft(x, dim=nd + k % 2).tobytes() == x.tobytes()
        assert F.ifft(x, dim=nd + 1).tobytes() == np.asfortranarray(x).tobytes()
        y = F.fftn(x, shape + [1, 1])
        assert y.shape == tuple(shape + [1, 1])
        assert y.tobytes() == F.fftn(x).tobytes()


# MATLAB R2024b's outcomes of the same calls (ones, zeros, the scalar 5, 'symmetric', dims from 1):
# an error, or class, realness, MATLAB's size and sha256 of the values' hex digits (num2hex, column
# order). siz needs ndims(X) elements, 2 at least; an empty siz is the input's size; n and siz must
# be non-negative integers.
MATLAB_ARGUMENT_RULES = """
fftn(ones(4, 1), [8])                          | error
fftn(ones(1, 4), [8])                          | error
fftn(five, [3])                                | error
fftn(five, [])                                 | double real 1x1 5c399a544777c03b
fftn(ones(4, 1), [])                           | double real 4x1 3326a5df512f4fc0
fftn(ones(4, 3), [])                           | double real 4x3 101fa2aa3bf016f9
fftn(ones(4, 3), [8])                          | error
fftn(ones(4, 3, 2), [4, 3])                    | error
fftn(ones(4, 1), [8, 2])                       | double complex 8x2 30b38e59e40bdceb
fftn(ones(4, 1), [0, 1])                       | double real 0x1 -
fftn(ones(4, 3), [4, -1])                      | error
fftn(ones(4, 3), [2.5, 3])                     | error
fft(ones(4, 3), 0)                             | double real 0x3 -
fft(ones(4, 3), 2.5)                           | error
fft(ones(4, 3), -1)                            | error
ifftn(ones(4, 1), [8])                         | error
ifftn(five, [])                                | double real 1x1 5c399a544777c03b
ifftn(ones(4, 3), [])                          | double real 4x3 a4028499a61557f7
fftn(ones(4, 3), [4, 3, 2])                    | double real 4x3x2 3b12c2a817e2a230
ifftn(ones(4, 3), [4, 3, 2], True)             | double real 4x3x2 bd2977f6480bbc3e
fftn(ones(4, 1, 3), [8])                       | error
fft(zeros(0, 3))                               | double real 0x3 -
fftn(zeros(0, 3))                              | double real 0x3 -
fft(zeros(3, 0))                               | double real 3x0 -
ifft(zeros(0, 3))                              | double real 0x3 -
ifftn(zeros(4, 0))                             | double real 4x0 -
fft(ones(4, 3), 8, 2)                          | double real 4x3x8 329c4fbca541496f
fft(ones(4, 3), 0, 2)                          | double real 4x3x0 -
fftn(five, [1])                                | error
fftn(five, [1, 1])                             | double real 1x1 5c399a544777c03b
ifftn(five, [3], True)                         | error
fftn(ones(1, 4), [1, 8])                       | double complex 1x8 5021c2463f072f34
fftn(ones(1, 4), [8, 1])                       | double real 8x1 4409c369140abdeb
fftn(ones(4, 3), [4, 3, 1, 1])                 | double real 4x3 101fa2aa3bf016f9
fft(ones(4, 3), None, 1)                       | double real 4x3 e3466cf4650ee6e1
fftn(zeros(0, 3), [2, 3])                      | double real 2x3 cb0216e7ae909ac5
ifft(five, 1)                                  | double real 1x1 5c399a544777c03b
fft(np.ones((4, 1), dtype=bool))               | double real 4x1 3326a5df512f4fc0
fftn(np.array([[1], [2], [3]], dtype=np.int8)) | double complex 3x1 adde030bc0191df2
fft(np.full((1, 1), 5, dtype=np.float32), 2)   | single real 1x2 d5c178b5d582d28e
"""


def _outcome(y):
    v = np.asfortranarray(y).ravel(order="F")
    cls = "double" if v.dtype in (np.float64, np.complex128) else "single"
    size = "x".join(str(s) for s in matlab_shape(y.shape))
    real = "complex" if np.iscomplexobj(v) else "real"
    if np.iscomplexobj(v):
        v = v.view(np.float64 if cls == "double" else np.float32)
    w = 16 if cls == "double" else 8
    digits = "".join(f"{int(u):0{w}x}" for u in v.view(np.uint64 if w == 16 else np.uint32))
    digest = hashlib.sha256(digits.encode()).hexdigest()[:16] if v.size else "-"
    return f"{cls} {real} {size} {digest}"


def test_arguments_follow_matlab():
    names = {
        "fft": F.fft,
        "ifft": F.ifft,
        "fftn": F.fftn,
        "ifftn": F.ifftn,
        "np": np,
        "ones": lambda *s: np.ones(s, order="F"),
        "zeros": lambda *s: np.zeros(s, order="F"),
        "five": np.full((1, 1), 5.0),
    }
    # the values' digests hold where the packs' do (FFTW's plans and twiddle factors); class,
    # realness and size everywhere
    values = digests_hold_here() == ""
    for line in MATLAB_ARGUMENT_RULES.strip().splitlines():
        call, want = (p.strip() for p in line.split("|"))
        if want == "error":
            with pytest.raises(ValueError):
                eval(call, names)
        elif values:
            assert _outcome(eval(call, names)) == want, call
        else:
            assert _outcome(eval(call, names)).rsplit(" ", 1)[0] == want.rsplit(" ", 1)[0], call
    # a 1-D array is a column: siz needs two elements; integer values of any type are counts
    with pytest.raises(ValueError):
        F.fftn(np.ones(4), [8])
    x = make_array(40, [6, 5], "double", True)
    assert F.fftn(np.ones(4), [8, 1]).tobytes() == F.fftn(np.ones((4, 1)), [8, 1]).tobytes()
    assert F.fft(x, 8.0, np.int64(1)).tobytes() == F.fft(x, 8, 1).tobytes()
    assert F.fftn(x, np.array([7.0, 6.0])).tobytes() == F.fftn(x, [7, 6]).tobytes()
    with pytest.raises(ValueError):
        F.fft(x, None, 0.5)


def test_the_two_matlab_modes_differ_only_in_the_lines_run_alone():
    X = make_array(12, [64, 2000], "single", False)
    a = F.fft(X, dim=0)
    b = F.fft(X, dim=0, matlab_parallel=False)
    split = mark_alone(2000, 64, 0, 32)
    serial = mark_alone(2000, 64, 0, 0)
    same = split == serial
    assert a[:, same].tobytes() == b[:, same].tobytes()
    if a[:, ~same].tobytes() == b[:, ~same].tobytes():
        pytest.skip("FFTW gives the same bits for groups of four and single lines here")


@needs_matlab_platform
def test_plans_are_matlabs():
    assert ext_plan("r2c", 64, 4) == '(rdft2-r2hc-direct-64-x4 "r2cf_64")'
    assert ext_plan("c2c-1", 64, 4) == '(dft-direct-64-x4 "n1fv_64_avx")'
    assert ext_plan("c2c+1", 64, 1) == (
        '(dft-ct-dit/8 (dftw-direct-8/24 "t3bv_8_avx") (dft-direct-8-x8 "n2bv_8_avx"))'
    )
    assert ext_plan("c2r", 97, 1) == "(rdft2-rdft-hc2r-97/1-0 (rdft-generic-hc2r-97) (rdft2-nop))"
    if MATLAB_PLANS_DIGEST:
        assert plans_digest() == MATLAB_PLANS_DIGEST


@pytest.mark.oracle
@pytest.mark.parametrize("pack", ["1d", "edge"])
@pytest.mark.parametrize("mode", ["multi", "single"])
def test_oracle_pack_equals_matlab(pack, mode):
    if (pack, mode) not in MATLAB_DIGESTS:
        pytest.skip("no MATLAB digest recorded for this pack")
    reason = digests_hold_here()
    if reason:
        pytest.skip(reason)
    digest, count = pack_digest(pack, mode == "multi")
    assert (digest, count) == MATLAB_DIGESTS[(pack, mode)]


@pytest.mark.parametrize("pack", ["1d", "edge"])
def test_pack_results_do_not_depend_on_our_thread_count(pack):
    # MATLAB's split of a stage follows its own FFT mode (the grain), never our thread count, so the
    # packs' digests hold on a runner with any number of CPUs
    for mode in (True, False):
        d = {t: pack_digest(pack, mode, threads=t, step=7) for t in (1, 3, 8)}
        assert d[1] == d[3] == d[8], d


def test_twiddle_fingerprint_sees_another_c_library():
    def bumped_sin(v):
        return math.nextafter(math.sin(v), math.inf) if v > 0.5 else math.sin(v)

    a = twiddle_fingerprint(64)
    assert a == twiddle_fingerprint(64)
    assert twiddle_fingerprint(64, sin=bumped_sin) != a
    # the arguments are FFTW's: in the first octant, 0 and pi/4 included where n is a multiple of 8
    t = twiddle_arguments(8)
    assert t[0] == 0 and t[-1] == _K2PI / 8 and np.all(t <= math.pi / 4 + 1e-15)


@pytest.mark.skipif(not os.environ.get("PETAKIT5D_TEST_FFTW"), reason="PETAKIT5D_TEST_FFTW not set")
def test_switch_to_another_fftw():
    # the driver on MATLAB's own FFTW gives the same bits as on ours
    info = F.fftw_info()
    assert info["name"] != "builtin" and info["single"].startswith("fftw-3.3")
    x = make_array(13, [64, 300], "single", False)
    y = make_array(14, [47, 45, 33], "double", True)
    a = (F.fft(x), F.fftn(y), F.ifftn(F.fftn(x)))
    try:
        F.use_fftw()
        assert F.fftw_info()["name"] == "builtin"
        b = (F.fft(x), F.fftn(y), F.ifftn(F.fftn(x)))
    finally:
        F.use_fftw(*F._library_paths(os.environ["PETAKIT5D_TEST_FFTW"]))
    assert all(u.tobytes() == v.tobytes() for u, v in zip(a, b, strict=True))
