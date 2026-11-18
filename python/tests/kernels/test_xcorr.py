"""The registration's NCC in the window of the allowed shifts (petakit5d.kernels.xcorr) against
normxcorr3_window_mex, normxcorr3_fast and normxcorr3_max_shift of the oracle tree.

Oracle cases (inputs in MATLAB's form; the roles give the Python calls (z, y, x) arrays and Python
boxes and map the results back). Each case runs on both CPU FFT backends:

- exact (the FFTW add-on): bit for bit (test_matlab_exact; skipped without the add-on, and where
  its plans can't be MATLAB's), and again with 1 and 24 threads (test_thread_counts);
- permissive (pocketfft, the default): within the FFT tolerances of decision record 0006
  (test_matlab_permissive): the correlation to 1e-5 of its Cauchy-Schwarz bound, the window sums
  exactly (no FFT), the NCC and the peak values to 5e-3, the integer peaks exactly.

Cases:

- mex.<shape>.<window>.<grain>: normxcorr3_window_mex('ncc', T, A, bbox, 4, grain) on
  RegistrationXcorrTest's pairs (even and odd sizes, a length of 64 along y, x or z, lines MATLAB
  splits among its threads, thin templates, one slice, 2D, equal sizes) in its windows (the whole
  output, the first and the last element, an edge, a centre part, random parts), MATLAB's two FFT
  modes; mexfast.*: the same with fastFFT;
- fast.<shape>.<window>: normxcorr3_fast(T, A, 'full', bbox); whole.<shape>.<shape-arg>: without
  a box ('full', 'same', 'valid'); fastfft.*: with fastFFT; double.*, thin.*: MATLAB's code
  without the mex (double data, an output of one row or one column);
- shift.*: normxcorr3_max_shift on textures with a known shift, symmetric and asymmetric bounds,
  a large offset with low contrast, flat and faint regions, 2D, fastFFT;
- stats.*, stats_wide.*: the template's and the image's means, sums, standard deviation and
  variance as normxcorr3_fast takes them (fixed_order_stat: MATLAB's order at one thread), across
  the thresholds of MATLAB's summation rule, on offset data and on data over 80 binary orders;
- threads.*: RegistrationXcorrTest/nccDoesNotDependOnTheThreadCount's pair (double and single,
  both past MATLAB's threaded sums);
- err.*: the errors of the mex and of normxcorr3_fast.
"""

from __future__ import annotations

import dataclasses
import functools

import numpy as np
import pytest

from petakit5d import fft as pfft
from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.core.axes import from_yxz
from petakit5d.core.bbox import box_to_matlab
from petakit5d.core.dtypes import mround
from petakit5d.kernels import xcorr
from petakit5d.testing import oracle
from petakit5d.testing.cases import Box, Case, In, Value, Vol
from petakit5d.testing.compare import Custom, MArray, tol
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG, stable_exp, texture

M = "petakit5d.kernels.xcorr:"
REASON = "decision record 0006: the permissive FFT backend against MATLAB's FFTW"

# RegistrationXcorrTest's pairs: template and image sizes, MATLAB's [y x z]
SHAPES = {
    "evenOdd": ([20, 17, 9], [31, 33, 15]),
    "y64": ([32, 20, 5], [33, 25, 12]),
    "x64": ([10, 30, 7], [20, 35, 9]),
    "z64": ([7, 9, 31], [9, 13, 34]),
    "split": ([40, 64, 6], [60, 101, 20]),
    "thinTemplate": ([1, 5, 3], [2, 7, 4]),
    "singleSlice": ([30, 20, 1], [50, 60, 3]),
    "twoD": ([30, 40, 1], [50, 60, 1]),
    "equal": ([16, 16, 16], [16, 16, 16]),
}


# ---------------------------------------------------------------------------------------------
# inputs (Python's form; the cases give MATLAB to_matlab of them)
# ---------------------------------------------------------------------------------------------


def zyx(m):
    """the port's shape of MATLAB's size [y x z] (an image for nz = 1, as MATLAB sees it)"""
    y, x, z = m
    return (y, x) if z == 1 else (z, y, x)


@functools.cache
def pair(name: str):
    """RegistrationXcorrTest.pair: an image and a template with structure and zero borders (zero
    lines in the padding and in the data), single"""
    st, sa = SHAPES[name]
    rng = CounterRNG(6, f"pair/{name}")
    A = (rng.uniform(zyx(sa)) * 100.0).astype(np.float32)
    A[..., : min(2, sa[0]), :] = 0  # A(1 : min(2, end), :, :) = 0
    T = (rng.uniform(zyx(st)) * 50.0).astype(np.float32)
    T[..., -1] = 0  # T(:, end, :) = 0
    if st[2] > 2:
        T[0] = 0  # T(:, :, 1) = 0
    return T, A


def full_size(name: str) -> list[int]:
    st, sa = SHAPES[name]
    return [p + q - 1 for p, q in zip(st, sa, strict=True)]


def windows(name: str):
    """RegistrationXcorrTest.windows: MATLAB bboxes of the whole output, the first and the last
    element, an edge, a centre part and three random parts"""
    sz = full_size(name)
    out = {
        "whole": [1, 1, 1, *sz],
        "first": [1, 1, 1, 1, 1, 1],
        "last": [*sz, *sz],
        "edge": [1, sz[1], 1, sz[0], sz[1], sz[2]],
    }
    s = [max(1, n // 3) for n in sz]
    out["centre"] = [*s, *(min(n, a + n // 4) for a, n in zip(s, sz, strict=True))]
    rng = CounterRNG(7, f"windows/{name}")
    for k in range(3):
        a = [int(rng.integers(1, n, 1)[0]) for n in sz]
        b = [int(rng.integers(lo, n, 1)[0]) for lo, n in zip(a, sz, strict=True)]
        out[f"rand{k}"] = [*a, *b]
    return out


def tex(sz, sigma=2.0, int_range=(100.0, 2000.0), stream="tex"):
    return texture(sz, sigma, int_range, rng=CounterRNG(0, stream))


def moved_texture():
    """RegistrationXcorrTest's moved texture: A = im(11:70, 6:65, 6:35), T = im(25:54, 18:47,
    13:27) of a 90 x 80 x 40 texture; T starts at (15, 13, 8) in A"""
    im = tex((40, 90, 80))
    return im[12:27, 24:54, 17:47].copy(), im[5:35, 10:70, 5:65].copy()


def low_contrast_texture():
    """lowContrastOnALargeOffset: 16-bit data shifted by 2^15 with a contrast of 12 counts"""
    im = mround(tex((40, 90, 80), 2.0, (2.0**15 - 6, 2.0**15 + 6)))
    return im[12:27, 24:54, 17:47].copy(), im[5:35, 10:70, 5:65].copy()


def flat_texture():
    """flatWindowsAreZero: a wide range of intensities, exactly flat from x = 31 on"""
    im = tex((40, 90, 80)).astype(np.float64)
    im = stable_exp(12.0 * (im - 100.0) / 1900.0).astype(np.float32)
    im[..., 30:] = np.float32(1000.1)
    return im[12:27, 24:54, 17:47].copy(), im[5:35, 10:70, 5:65].copy()


def faint_texture():
    """faintRegionsHaveTheirNCC: a contrast of a few counts next to hundreds"""
    im = tex((40, 90, 80))
    f = tex((40, 90, 80), 2.0, (0.0, 6.0), stream="faint")
    im[..., 30:] = np.float32(1000) + f[..., 30:]
    return im[12:27, 24:54, 17:47].copy(), im[5:35, 10:70, 5:65].copy()


def texture_2d():
    """maxShiftOf2DImages: A = im(11:90, 6:85), T = im(31:70, 21:60) of a 120 x 100 texture"""
    im = tex((120, 100))
    return im[30:70, 20:60].copy(), im[10:90, 5:85].copy()


def large_texture():
    """a template of 108,000 voxels (MATLAB's sums of it depend on its thread count) in an image of
    256,000"""
    im = tex((50, 100, 100), stream="large")
    return im[10:40, 20:80, 20:80].copy(), im[5:45, 10:90, 10:90].copy()


def offset_pair(sz_t, sz_a, off, seed, sigma=1.5):
    """a texture image and a template cut from it at the 0-based offset off ((z, y, x))"""
    im = tex(sz_a, sigma, stream=f"offset/{seed}")
    sl = tuple(slice(o, o + n) for o, n in zip(off, sz_t, strict=True))
    return im[sl].copy(), im


def stats_pair(n_t: int, n_a: int, seed: int):
    """vectors of n_t and n_a values with an offset (the sums' rounding shows)"""
    rng = CounterRNG(seed, "stats")
    T = (rng.uniform((n_t, 1)) * 1000.0 + 30000.0).astype(np.float32)
    A = (rng.uniform((n_a, 1)) * 1000.0 + 30000.0).astype(np.float32)
    return T, A


def stats_wide(n_t: int, n_a: int, seed: int):
    """vectors of n_t and n_a values of both signs over 80 binary orders of magnitude (no sum of
    them is exact, so the order of the additions shows in every sum)"""
    rng = CounterRNG(seed, "stats_wide")

    def draw(n):
        u = rng.uniform((n, 1)) - 0.5
        return np.ldexp(u, rng.integers(-40, 40, (n, 1)).astype(np.int32)).astype(np.float32)

    return draw(n_t), draw(n_a)


def mval(make, k):
    """a Value input of MATLAB's form of the k-th array that make() gives"""
    return Value(lambda: to_matlab(make()[k]))


# ---------------------------------------------------------------------------------------------
# the Python sides
# ---------------------------------------------------------------------------------------------


def py_window(T, A, box, fast=False, **kw):
    return xcorr.normxcorr3_window(T, A, box, fast, **kw)


def py_shift(T, A, max_shifts, fast=False, **kw):
    off, mx, C = xcorr.normxcorr3_max_shift(T, A, max_shifts, fast, **kw)
    return [np.asarray(off, dtype=np.float64), mx, C]


def py_stats(T, A, **kw):
    """normxcorr3_fast's statistics, as the kernel computes them (MATLAB's form: column vectors
    are (n, 1) arrays, one slice)"""
    _, _, st = xcorr._ext().prepare(T.reshape(1, *T.shape), A.reshape(1, *A.shape), 1)
    return list(st)


def shifts_zyx(m):
    """MATLAB's maxShifts ([y x z] row, or [lower; upper] rows) in the port's (z, y, x)"""
    m = np.asarray(m, dtype=np.float64)
    if m.ndim == 1:
        return np.asarray(from_yxz(m))
    return np.stack([np.asarray(from_yxz(r)) for r in m])


# ---------------------------------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------------------------------

EXACT_KW = {"backend": "exact"}


def corr_bound(make):
    """the permissive backend's tolerance of a correlation: 1e-5 of sum(|T|) max(|A|) (the bound
    of every value of the full correlation), as decision record 0006 bounds a transform by its
    largest value"""

    def fn(ref: MArray, act: MArray) -> list[str]:
        T, A = make()
        bound = float(np.abs(T.astype(np.float64)).sum() * np.abs(A).max())
        d = np.abs(ref.data.astype(np.float64) - act.data.astype(np.float64))
        if ref.data.shape != act.data.shape or ref.mclass != act.mclass:
            return [f"class or size: {ref.mclass} {ref.data.shape}, {act.mclass} {act.data.shape}"]
        m = float(d.max(initial=0.0))
        return [] if m <= 1e-5 * bound else [f"max |difference| {m:.3g} > 1e-5 * {bound:.3g}"]

    return Custom(fn, "within 1e-5 of the correlation's bound (permissive FFT)")


NCC_TOL = tol(reason=REASON + ": the NCC (well-conditioned windows)", abs=5e-3)


def mex_cases():
    out = []
    for name in SHAPES:
        make = functools.partial(pair, name)
        for wname, bbox in windows(name).items():
            for grain in (32, 0):
                out.append(
                    Case(
                        f"mex.{name}.{wname}.{grain}",
                        matlab="normxcorr3_window_mex",
                        args=["ncc", Vol(In("T")), Vol(In("A")), Box(bbox), 4, grain],
                        nargout=4,
                        python=py_window,
                        py_args=[Vol(In("T")), Vol(In("A")), Box(bbox)],
                        py_kwargs={**EXACT_KW, "matlab_parallel": grain == 32},
                        inputs={"T": mval(make, 0), "A": mval(make, 1)},
                        out_roles=("vol", "vol", "vol", "vol"),
                        note="permissive:mex",
                    )
                )
        # fastFFT (MATLAB's mex takes the grain 0 then)
        for wname in ("whole", "centre", "rand0"):
            bbox = windows(name)[wname]
            out.append(
                Case(
                    f"mexfast.{name}.{wname}",
                    matlab="normxcorr3_window_mex",
                    args=["ncc", Vol(In("T")), Vol(In("A")), Box(bbox), 4, 32, True],
                    nargout=4,
                    python=py_window,
                    py_args=[Vol(In("T")), Vol(In("A")), Box(bbox), True],
                    py_kwargs=dict(EXACT_KW),
                    inputs={"T": mval(make, 0), "A": mval(make, 1)},
                    out_roles=("vol", "vol", "vol", "vol"),
                    note="permissive:mex",
                )
            )
    return out


def fast_cases():
    out = []
    for name in SHAPES:
        make = functools.partial(pair, name)
        inputs = {"T": mval(make, 0), "A": mval(make, 1)}
        for wname, bbox in windows(name).items():
            out.append(
                Case(
                    f"fast.{name}.{wname}",
                    matlab="normxcorr3_fast",
                    args=[Vol(In("T")), Vol(In("A")), "full", Box(bbox)],
                    python=M + "normxcorr3_fast",
                    py_kwargs=dict(EXACT_KW),
                    inputs=inputs,
                    out_roles=("vol",),
                    note="permissive:ncc",
                )
            )
        for shape in ("full", "same", "valid"):
            out.append(
                Case(
                    f"whole.{name}.{shape}",
                    matlab="normxcorr3_fast",
                    args=[Vol(In("T")), Vol(In("A")), shape],
                    python=M + "normxcorr3_fast",
                    py_kwargs=dict(EXACT_KW),
                    inputs=inputs,
                    out_roles=("vol",),
                    note="permissive:ncc",
                )
            )
        bbox = windows(name)["centre"]
        out.append(
            Case(
                f"fastfft.{name}",
                matlab="normxcorr3_fast",
                args=[Vol(In("T")), Vol(In("A")), "full", Box(bbox), True],
                python=M + "normxcorr3_fast",
                py_kwargs=dict(EXACT_KW),
                inputs=inputs,
                out_roles=("vol",),
                note="permissive:ncc",
            )
        )

    # MATLAB's code without the mex: double data, outputs of one row or one column
    def dbl(name):
        T, A = pair(name)
        return T.astype(np.float64), A.astype(np.float64)

    for name in ("evenOdd", "twoD", "thinTemplate"):
        make = functools.partial(dbl, name)
        bbox = windows(name)["centre"]
        for tag, args in (("box", ["full", Box(bbox)]), ("full", ["full"]), ("valid", ["valid"])):
            out.append(
                Case(
                    f"double.{name}.{tag}",
                    matlab="normxcorr3_fast",
                    args=[Vol(In("T")), Vol(In("A")), *args],
                    python=M + "normxcorr3_fast",
                    py_kwargs=dict(EXACT_KW),
                    inputs={"T": mval(make, 0), "A": mval(make, 1)},
                    out_roles=("vol",),
                    note="permissive:ncc",
                )
            )

    def thin(kind):
        rng = CounterRNG(8, f"thin/{kind}")
        if kind == "row":  # one row: the full output is 1 long along y
            return (
                (rng.uniform((3, 1, 5)) * 9.0).astype(np.float32),
                (rng.uniform((7, 1, 9)) * 9.0).astype(np.float32),
            )
        return (  # one column
            (rng.uniform((4, 6, 1)) * 9.0).astype(np.float32),
            (rng.uniform((8, 11, 1)) * 9.0).astype(np.float32),
        )

    for kind in ("row", "column"):
        make = functools.partial(thin, kind)
        t, a = make()
        full = [p + q - 1 for p, q in zip(t.shape, a.shape, strict=True)]
        box = (1, 0, 0, full[0] - 1, full[1], full[2])
        out.append(
            Case(
                f"thin.{kind}",
                matlab="normxcorr3_fast",
                args=[Vol(In("T")), Vol(In("A")), "full", Box(list(box_to_matlab(box)))],
                python=M + "normxcorr3_fast",
                py_kwargs=dict(EXACT_KW),
                inputs={"T": mval(make, 0), "A": mval(make, 1)},
                out_roles=("vol",),
                note="permissive:ncc",
            )
        )
    return out


def shift_cases():
    out = []

    def add(cid, make, max_shifts, fast=False, peak=True):
        args = [Vol(In("T")), Vol(In("A")), np.asarray(max_shifts, dtype=np.float64)]
        py_args = [Vol(In("T")), Vol(In("A")), shifts_zyx(max_shifts)]
        if fast:
            args.append(True)
            py_args.append(True)
        out.append(
            Case(
                cid,
                matlab="normxcorr3_max_shift",
                args=args,
                nargout=3,
                python=py_shift,
                py_args=py_args,
                py_kwargs=dict(EXACT_KW),
                inputs={"T": mval(make, 0), "A": mval(make, 1)},
                out_roles=("yxz", None, "vol"),
                note="permissive:shift" if peak else "permissive:shift_nopeak",
            )
        )

    bounds = [[-20, -15, -8], [25, 20, 9]]
    add("shift.moved", moved_texture, bounds)
    add("shift.moved.symmetric", moved_texture, [20, 15, 8])
    add("shift.moved.fastfft", moved_texture, bounds, fast=True)
    add("shift.moved.fractional", moved_texture, [[-19.5, -14.2, -7.7], [24.3, 19.9, 8.5]])
    add("shift.low_contrast", low_contrast_texture, bounds)
    add("shift.flat", flat_texture, [[-20, -15, -8], [25, 30, 9]])
    add("shift.faint", faint_texture, [[-20, -15, -8], [25, 30, 9]])
    add("shift.twoD", texture_2d, [[-30, -30, -1], [30, 30, 1]])
    add("shift.large_template", large_texture, [[-15, -15, -8], [25, 25, 13]])
    # synthetic overlaps of odd sizes with known offsets (MATLAB's [y x z] sizes)
    for k, (st, sa, off) in enumerate(
        [
            ((9, 21, 17), (13, 33, 29), (2, 7, 5)),
            ((5, 40, 33), (11, 61, 47), (3, 11, 9)),
            ((1, 25, 31), (1, 44, 52), (0, 6, 13)),
            ((12, 12, 12), (12, 12, 12), (0, 0, 0)),
        ]
    ):
        make = functools.partial(offset_pair, st, sa, off, k)
        add(f"shift.offset{k}", make, [n // 2 for n in (sa[1], sa[2], sa[0])])
    return out


def stats_cases():
    out = []
    sizes = [(5, 43), (43, 44), (44, 100), (1000, 88999), (88999, 89000), (100000, 250001)]
    wide = [(43, 44), (44, 45), (88999, 89000), (200000, 300001)]
    for k, (gen, nt, na) in enumerate(
        [(stats_pair, *s) for s in sizes] + [(stats_wide, *s) for s in wide]
    ):
        make = functools.partial(gen, nt, na, k)
        tag = "stats" if gen is stats_pair else "stats_wide"
        out.append(
            Case(
                f"{tag}.{nt}.{na}",
                # as normxcorr3_fast takes them: MATLAB's order at one thread (fixed_order_stat)
                code=(
                    "T = in{1}; A = in{2}; cT = cast(fixed_order_stat('mean', T, 'double'), "
                    "'like', T); cA = cast(fixed_order_stat('mean', A, 'double'), 'like', A); "
                    "T0 = T - cT; out = {double(cT), double(cA), "
                    "fixed_order_stat('sum', T0, 'double'), fixed_order_stat('std', double(T0)), "
                    "double(fixed_order_stat('var', A, 1))};"
                ),
                args=[In("T"), In("A")],
                nargout=5,
                python=py_stats,
                inputs={"T": Value(lambda m=make: m()[0]), "A": Value(lambda m=make: m()[1])},
                note="stats",
            )
        )
    return out


def error_cases():
    out = []

    def make():
        T, A = pair("evenOdd")
        return T[:4, :5, :6].copy(), A[:7, :9, :11].copy()

    inputs = {"T": mval(make, 0), "A": mval(make, 1)}
    bbox = [1, 1, 1, 13, 16, 10]
    mid = "normxcorr3_window_mex:inputError"

    def mex_err(cid, targ, aarg, b, pybox):
        out.append(
            Case(
                cid,
                matlab="normxcorr3_window_mex",
                args=["ncc", targ, aarg, b, 1, 32],
                python=py_window,
                py_args=[targ, aarg, pybox],
                expect_error=mid,
                py_error=xcorr.XcorrError,
                inputs=inputs,
                note="error",
            )
        )

    # the Python boxes: MATLAB's, converted where they convert (the full output is (10, 13, 16))
    T_, A_ = Vol(In("T")), Vol(In("A"))
    mex_err("err.mex.swapped", A_, T_, bbox, (0, 0, 0, 10, 13, 16))
    mex_err("err.mex.zero_start", T_, A_, [0, 1, 1, 13, 16, 10], (0, -1, 0, 10, 13, 16))
    mex_err("err.mex.past_end", T_, A_, [1, 1, 1, 14, 16, 10], (0, 0, 0, 10, 14, 16))
    mex_err("err.mex.s_after_t", T_, A_, [3, 1, 1, 2, 16, 10], (0, 2, 0, 10, 2, 16))
    mex_err("err.mex.fraction", T_, A_, [1.5, 1, 1, 13, 16, 10], (0, 0.5, 0, 10, 13, 16))
    mex_err("err.mex.five", T_, A_, bbox[:5], (0, 0, 0, 10, 13))

    def fast_err(cid, args, py_args):
        out.append(
            Case(
                cid,
                matlab="normxcorr3_fast",
                args=args,
                python=M + "normxcorr3_fast",
                py_args=py_args,
                expect_error=True,
                py_error=xcorr.XcorrError,
                inputs=inputs,
                note="error",
            )
        )

    fast_err(
        "err.fast.bbox_needs_full",
        [Vol(In("T")), Vol(In("A")), "same", Box([1, 1, 1, 2, 2, 2])],
        [Vol(In("T")), Vol(In("A")), "same", Box([1, 1, 1, 2, 2, 2])],
    )
    fast_err(
        "err.fast.bad_bbox",
        [Vol(In("T")), Vol(In("A")), "full", [1, 1, 1, 13, 16, 11]],
        [Vol(In("T")), Vol(In("A")), "full", (0, 0, 0, 11, 13, 16)],
    )
    fast_err("err.fast.larger_template", [Vol(In("A")), Vol(In("T"))], [Vol(In("A")), Vol(In("T"))])
    return out


def thread_pair(cls):
    """nccDoesNotDependOnTheThreadCount: a template of 92,250 and an image of 231,000 elements (both
    past MATLAB's threaded sums), uniform in [0, 100), double or single"""
    rng = CounterRNG(9, "threads")
    T = rng.uniform((41, 50, 45)) * 100.0
    A = rng.uniform((50, 70, 66)) * 100.0
    dt = np.float64 if cls == "double" else np.float32
    return T.astype(dt), A.astype(dt)


def thread_cases():
    out = []
    for cls in ("double", "single"):
        make = functools.partial(thread_pair, cls)
        out.append(
            Case(
                f"threads.{cls}",
                matlab="normxcorr3_fast",
                args=[Vol(In("T")), Vol(In("A")), "full", Box([30, 25, 20, 90, 85, 70])],
                python=M + "normxcorr3_fast",
                py_kwargs=dict(EXACT_KW),
                inputs={"T": mval(make, 0), "A": mval(make, 1)},
                out_roles=("vol",),
                note="permissive:ncc",
            )
        )
    return out


CASES = mex_cases() + fast_cases() + shift_cases() + stats_cases() + thread_cases() + error_cases()


# ---------------------------------------------------------------------------------------------
# running them on each backend
# ---------------------------------------------------------------------------------------------


def exact_skip_reason() -> str:
    """'' where the exact backend gives MATLAB's bits, else why not (no add-on, or other plans or
    twiddle factors than where the references were made)"""
    return digests_hold_here()


def variant(case: Case, **changes) -> Case:
    v = dataclasses.replace(case, **changes)
    v._module_globals = case._module_globals
    return v


def permissive(case: Case) -> Case:
    """the case on the permissive backend, with its tolerances"""
    kind = case.note.split(":")[-1]
    kw = {**(case.py_kwargs or {}), "backend": "permissive"}
    if kind == "mex":
        mk = case.inputs

        def make():
            return from_matlab(mk["T"].make()), from_matlab(mk["A"].make())

        compare = {0: corr_bound(make)}
    elif kind == "ncc":
        compare = NCC_TOL
    elif kind == "shift":
        compare = {1: NCC_TOL, 2: NCC_TOL}
    elif kind == "shift_nopeak":
        compare = {0: Custom(lambda r, a: [], "not compared"), 1: NCC_TOL, 2: NCC_TOL}
    else:
        compare = case.compare
    return variant(case, py_kwargs=kw, compare=compare)


def uses_fft(case: Case) -> bool:
    return case.note.startswith("permissive")


@oracle.parametrize(CASES)
def test_matlab_exact(case, tmp_path):
    """bit for bit with the exact FFT backend (and the cases without FFTs on any backend)"""
    if uses_fft(case) and exact_skip_reason():
        pytest.skip(exact_skip_reason())
    if not uses_fft(case) and case.py_kwargs and case.py_kwargs.get("backend") == "exact":
        case = variant(case, py_kwargs={**case.py_kwargs, "backend": "permissive"})
    oracle.check(case, tmp_path)


@oracle.parametrize([c for c in CASES if uses_fft(c)])
def test_matlab_permissive(case, tmp_path):
    """within decision record 0006's tolerances on the permissive FFT backend"""
    oracle.check(permissive(case), tmp_path)


THREAD_CASES = [
    c for c in CASES if c.id.startswith(("mex.split", "fast.split", "shift.moved", "threads."))
]


@pytest.mark.parametrize("threads", [1, 24])
@oracle.parametrize(THREAD_CASES)
def test_thread_counts(case, threads, tmp_path):
    """MATLAB's bits with 1 and 24 threads (exact backend)"""
    if exact_skip_reason():
        pytest.skip(exact_skip_reason())
    oracle.check(variant(case, py_kwargs={**case.py_kwargs, "threads": threads}), tmp_path)


# ---------------------------------------------------------------------------------------------
# RegistrationXcorrTest, in intent (both backends)
# ---------------------------------------------------------------------------------------------


@pytest.fixture(params=["permissive", "exact"])
def be(request):
    if request.param not in pfft.available():
        pytest.skip("the exact FFT backend (the FFTW add-on) is not installed")
    return request.param


def py_box(mbox):
    from petakit5d.core.bbox import box_from_matlab

    return box_from_matlab(mbox)


def crop_box(C, box):
    if C.ndim == 2:
        return C[box[1] : box[4], box[2] : box[5]]
    return C[box[0] : box[3], box[1] : box[4], box[2] : box[5]]


def bits(a):
    a = np.ascontiguousarray(a)
    return a.view(np.uint32 if a.dtype == np.float32 else np.uint64)


@pytest.mark.parametrize("name", list(SHAPES))
def test_window_equals_the_whole_output(name, be):
    """RegistrationXcorrTest/windowEqualsTheWholeOutput: every window equals that part of the
    whole output, bit for bit"""
    T, A = pair(name)
    C = xcorr.normxcorr3_fast(T, A, backend=be)
    for mbox in windows(name).values():
        box = py_box(mbox)
        Cw = xcorr.normxcorr3_fast(T, A, "full", box, backend=be)
        np.testing.assert_array_equal(bits(Cw), bits(crop_box(C, box)), err_msg=str(mbox))


def box_sums_reference(X, lo, hi):
    """box_sums of normxcorr3_fast's MATLAB code: cumulative sums along y, then x, then z (MATLAB
    layout), P(hi) - P(lo - 1)"""
    S = X
    for d in range(3):
        P = np.cumsum(S, axis=d)
        W = np.take(P, np.asarray(hi[d]) - 1, axis=d)
        m = np.asarray(lo[d]) > 1
        if m.any():
            L = np.take(P, np.asarray(lo[d])[m] - 2, axis=d)
            idx = [slice(None)] * 3
            idx[d] = np.flatnonzero(m)
            W[tuple(idx)] = W[tuple(idx)] - L
        S = W
    return S


@pytest.mark.parametrize("name", ["evenOdd", "split", "thinTemplate", "twoD"])
def test_window_sums_are_the_cumulative_sums(name, be):
    """RegistrationXcorrTest/mexPartsEqualMatlabs (the window sums): MATLAB's box_sums in double,
    bit for bit, on any backend"""
    T, A = pair(name)
    st, sa = SHAPES[name]
    for mbox in windows(name).values():
        _, sA, sA2, sT = xcorr.normxcorr3_window(T, A, py_box(mbox), backend=be)
        loT, hiT, loA, hiA = [], [], [], []
        for d in range(3):
            k = np.arange(mbox[d], mbox[d + 3] + 1)
            loT.append(np.maximum(1, st[d] - k + 1))
            hiT.append(np.minimum(st[d], st[d] - k + sa[d]))
            loA.append(np.maximum(1, k - st[d] + 1))
            hiA.append(np.minimum(k, sa[d]))
        Xm = to_matlab(A).astype(np.float64).reshape(sa, order="F")
        Tm = to_matlab(T).astype(np.float64).reshape(st, order="F")
        for got, ref in (
            (sA, box_sums_reference(Xm, loA, hiA)),
            (sA2, box_sums_reference(Xm * Xm, loA, hiA)),
            (sT, box_sums_reference(Tm, loT, hiT)),
        ):
            np.testing.assert_array_equal(
                bits(to_matlab(got).reshape(ref.shape, order="F")), bits(ref)
            )


def corr_reference(T0, A0):
    """the full correlation of T0 and A0 (port's arrays) in double, directly from its definition"""
    from numpy.fft import irfftn, rfftn

    t3 = T0.reshape((1,) * (3 - T0.ndim) + T0.shape).astype(np.float64)
    a3 = A0.reshape((1,) * (3 - A0.ndim) + A0.shape).astype(np.float64)
    full = [p + q - 1 for p, q in zip(t3.shape, a3.shape, strict=True)]
    ax = (0, 1, 2)
    c = irfftn(rfftn(a3, full, ax) * rfftn(t3[::-1, ::-1, ::-1], full, ax), full, ax)
    return c


@pytest.mark.parametrize("name", ["evenOdd", "z64", "split", "twoD"])
def test_correlation_is_the_correlation(name, be):
    """the mex's correlation equals the correlation of double precision to the FFTs' rounding
    (the exact backend's bits are MATLAB's: the oracle cases)"""
    T, A = pair(name)
    ref = corr_reference(T, A)
    bound = float(np.abs(T.astype(np.float64)).sum() * np.abs(A).max())
    for mbox in windows(name).values():
        box = py_box(mbox)
        c, *_ = xcorr.normxcorr3_window(T, A, box, backend=be)
        r = ref[box[0] : box[3], box[1] : box[4], box[2] : box[5]].reshape(c.shape)
        assert np.abs(c - r).max() <= 1e-5 * bound


@pytest.mark.parametrize("name", list(SHAPES))
def test_mex_path_equals_the_matlab_code(name, be):
    """RegistrationXcorrTest/withoutTheMexTheMatlabCodeRuns: MATLAB's code (fftn and ifftn of the
    whole output, then cropped) gives the mex's bits"""
    T, A = pair(name)
    ext = xcorr._ext()
    t3, a3 = T.reshape(xcorr._shape3(T)), A.reshape(xcorr._shape3(A))
    full = xcorr._full_shape(T, A)
    for mbox in list(windows(name).values())[:5]:
        win = [mbox[0], mbox[1], mbox[2], mbox[3], mbox[4], mbox[5]]
        C, imax = xcorr._ncc_by_fftn(ext, t3, a3, full, win, 4, xcorr._backend(be), True, be)
        Cm = xcorr.normxcorr3_fast(T, A, "full", py_box(mbox), backend=be)
        if be == "exact":
            np.testing.assert_array_equal(bits(C.reshape(Cm.shape)), bits(Cm))
        else:
            assert np.abs(C.reshape(Cm.shape) - Cm).max() <= 1e-4


def ncc_reference(T, A, box):
    """the NCC of each window of the box from its elements in double (the image zero-padded): 0
    for flat windows; and which windows are flat"""
    t = T.astype(np.float64).reshape(xcorr._shape3(T))
    a = A.astype(np.float64).reshape(xcorr._shape3(A))
    nt = t.size
    pad = np.zeros([q + 2 * p for p, q in zip(t.shape, a.shape, strict=True)])
    pad[tuple(slice(p, p + q) for p, q in zip(t.shape, a.shape, strict=True))] = a
    t0 = t - t.mean()
    st = t.std(ddof=1)
    shape = [box[3 + i] - box[i] for i in range(3)]
    C = np.zeros(shape)
    flat = np.zeros(shape, dtype=bool)
    from numpy.lib.stride_tricks import sliding_window_view

    win = sliding_window_view(pad, t.shape)
    for idx in np.ndindex(*shape):
        w = win[tuple(1 + b + i for b, i in zip(box[:3], idx, strict=True))]
        sw = w.std(ddof=1)
        flat[idx] = sw == 0
        if sw > 0:
            C[idx] = (t0 * w).sum() / ((nt - 1) * st * sw)
    return C, flat


@pytest.mark.parametrize(
    "make, expected",
    [
        (moved_texture, (7, 14, 12)),
        (low_contrast_texture, (7, 14, 12)),
        (flat_texture, (7, 14, 12)),
        (faint_texture, (7, 14, 12)),
    ],
    ids=["moved", "lowContrastOnALargeOffset", "flatWindowsAreZero", "faintRegionsHaveTheirNCC"],
)
def test_max_shift_finds_the_shift(make, expected, be):
    """RegistrationXcorrTest/maxShiftEqualsTheWholeOutputsMax, lowContrastOnALargeOffset,
    flatWindowsAreZero, faintRegionsHaveTheirNCC: the shift of the template, a peak of 1, the NCC
    of double precision in the window (0 for flat windows)"""
    T, A = make()
    bounds = shifts_zyx([[-20, -15, -8], [25, 30, 9]])
    off, mx, C = xcorr.normxcorr3_max_shift(T, A, bounds, backend=be)
    assert off == expected
    assert abs(float(mx) - 1) <= 1e-5
    assert np.abs(C).max() <= 1 + 1e-5
    # the same as the whole output's maximum in the window
    box = xcorr.max_shift_box(T.shape, A.shape, bounds)
    Cw = crop_box(xcorr.normxcorr3_fast(T, A, backend=be), box)
    np.testing.assert_array_equal(bits(Cw), bits(C))
    # the NCC of double precision, on a part of the window around the peak
    sub = (2, 10, 10, 9, 18, 18)
    ref, flat = ncc_reference(T, A, tuple(b + s for b, s in zip(box[:3] * 2, sub, strict=True)))
    got = C[sub[0] : sub[3], sub[1] : sub[4], sub[2] : sub[5]].astype(np.float64)
    assert np.abs(got - ref).max() <= 1e-4
    assert not got[flat].any()


def test_flat_windows_are_zero(be):
    """flatWindowsAreZero: windows of exactly flat data give 0 (not their sums' rounding)"""
    T, A = flat_texture()
    bounds = shifts_zyx([[-20, -15, -8], [25, 30, 9]])
    _, _, C = xcorr.normxcorr3_max_shift(T, A, bounds, backend=be)
    box = xcorr.max_shift_box(T.shape, A.shape, bounds)
    # windows inside the image, in its flat part (A's x >= 25)
    sub = (8, 20, 40, 12, 26, 46)
    ref, flat = ncc_reference(T, A, tuple(b + s for b, s in zip(box[:3] * 2, sub, strict=True)))
    assert flat.sum() > 100
    got = C[sub[0] : sub[3], sub[1] : sub[4], sub[2] : sub[5]]
    assert not got[flat].any()
    assert np.abs(got.astype(np.float64) - ref).max() <= 1e-4


def test_max_shift_of_2d_images(be):
    """RegistrationXcorrTest/maxShiftOf2DImages: 2D images (MATLAB's maxShifts of three values)"""
    T, A = texture_2d()
    off, mx, C = xcorr.normxcorr3_max_shift(
        T, A, shifts_zyx([[-30, -30, -1], [30, 30, 1]]), backend=be
    )
    assert off == (0, 20, 15)
    assert abs(float(mx) - 1) <= 1e-3
    assert C.ndim == 2


def test_symmetric_bounds_are_minus_to_plus(be):
    T, A = moved_texture()
    a = xcorr.normxcorr3_max_shift(T, A, (8, 20, 15), backend=be)
    b = xcorr.normxcorr3_max_shift(T, A, [[-8, -20, -15], [8, 20, 15]], backend=be)
    assert a[0] == b[0]
    np.testing.assert_array_equal(bits(a[2]), bits(b[2]))


def test_fast_fft_is_close_in_every_window(be):
    """RegistrationXcorrTest/fastFFTIsCloseInEveryWindow: smaller FFTs change the correlation only
    by rounding; the window sums stay the same bits"""
    for name in SHAPES:
        T, A = pair(name)
        ref = corr_reference(T, A)
        for mbox in windows(name).values():
            box = py_box(mbox)
            c, s1, s2, s3 = xcorr.normxcorr3_window(T, A, box, True, backend=be)
            _, r1, r2, r3 = xcorr.normxcorr3_window(T, A, box, backend=be)
            r = ref[box[0] : box[3], box[1] : box[4], box[2] : box[5]].reshape(c.shape)
            assert np.abs(c - r).max() <= 1e-5 * np.abs(ref).max()
            for got, want in ((s1, r1), (s2, r2), (s3, r3)):
                np.testing.assert_array_equal(bits(got), bits(want))


def test_fast_fft_finds_the_same_shift(be):
    T, A = moved_texture()
    bounds = shifts_zyx([[-20, -15, -8], [25, 20, 9]])
    off, mx, _ = xcorr.normxcorr3_max_shift(T, A, bounds, backend=be)
    off_f, mx_f, _ = xcorr.normxcorr3_max_shift(T, A, bounds, True, backend=be)
    assert off_f == off
    assert abs(float(mx_f) - float(mx)) <= 1e-5


def test_valid_and_same_shapes(be):
    """RegistrationXcorrTest/validAndSameShapes: 'valid' starts along x at the template's size
    along x (MATLAB d67e1b7's fix), 'same' is the image's size"""
    T, A = pair("evenOdd")
    T, A = T[:4, :6, :9].copy(), A[:11, :20, :22].copy()
    C = xcorr.normxcorr3_fast(T, A, backend=be)
    Cv = xcorr.normxcorr3_fast(T, A, "valid", backend=be)
    assert Cv.shape == tuple(a - t + 1 for a, t in zip(A.shape, T.shape, strict=True))
    np.testing.assert_array_equal(bits(Cv), bits(C[3:-3, 5:-5, 8:-8]))
    Cs = xcorr.normxcorr3_fast(T, A, "same", backend=be)
    np.testing.assert_array_equal(bits(Cs), bits(C[1 : 1 + 11, 2 : 2 + 20, 4 : 4 + 22]))


@pytest.mark.parametrize("name", ["split", "evenOdd", "twoD"])
def test_results_do_not_depend_on_the_thread_count(name, be):
    """RegistrationXcorrTest/resultsDoNotDependOnTheThreadCount (both of MATLAB's FFT modes)"""
    T, A = pair(name)
    box = py_box(windows(name)["centre"])
    ext = xcorr._ext()
    win = xcorr._win(box, xcorr._full_shape(T, A), "", "")
    t3, a3 = T.reshape(xcorr._shape3(T)), A.reshape(xcorr._shape3(A))
    for grain in (32, 0):
        runs = [ext.ncc(t3, a3, win, n, grain, False, xcorr._backend(be)) for n in (1, 3, 24)]
        for r in runs[1:]:
            np.testing.assert_array_equal(bits(r[0]), bits(runs[0][0]))
            assert r[1] == runs[0][1] and r[2] == runs[0][2]


def test_double_data_runs_the_matlab_code(be):
    """RegistrationXcorrTest/unsupportedCallsRunTheMatlabCode: double data give a double NCC in
    the window, the part of the whole output; a box needs SHAPE 'full'"""
    T, A = pair("evenOdd")
    T, A = T[:5, :12, :10].astype(np.float64), A[:9, :20, :22].astype(np.float64)
    C = xcorr.normxcorr3_fast(T, A, backend=be)
    box = (0, 1, 2, 13, 30, 31)
    Cw = xcorr.normxcorr3_fast(T, A, "full", box, backend=be)
    assert Cw.dtype == np.float64
    np.testing.assert_array_equal(bits(Cw), bits(crop_box(C, box)))
    with pytest.raises(xcorr.XcorrError):
        xcorr.normxcorr3_fast(T, A, "same", (0, 0, 0, 2, 2, 2), backend=be)


def test_window_rejects_invalid_inputs():
    """RegistrationXcorrTest/mexRejectsInvalidInputs, with the mex's identifier"""
    T, A = pair("evenOdd")
    T, A = T[:4, :5, :6].copy(), A[:7, :9, :11].copy()
    box = (0, 0, 0, 10, 13, 16)
    bad = [
        (T.astype(np.float64), A, box),
        (A, T, box),
        (T, A, (0, 0, 0, 10, 13, 17)),
        (T, A, (3, 0, 0, 2, 13, 16)),
        (T, A, (0, 0.5, 0, 10, 13, 16)),
        (T, A, box[:5]),
        (T, A[None, None], box),
        (T[:, :1, :1], A[:, :1, :1], (0, 0, 0, 10, 1, 1)),
    ]
    for t, a, b in bad:
        with pytest.raises(xcorr.XcorrError) as err:
            xcorr.normxcorr3_window(t, a, b)
        assert err.value.identifier == "normxcorr3_window_mex:inputError"


def test_max_shift_box_is_matlabs_window():
    """the window of normxcorr3_max_shift: max(1, ceil(lower + szT)) .. min(floor(upper + szT),
    szT + szA - 1), as a Python box"""
    box = xcorr.max_shift_box((15, 30, 30), (30, 60, 60), [[-8, -20, -15.5], [9.7, 25, 100]])
    assert box == (6, 9, 14, 24, 55, 89)
    assert xcorr.max_shift_box((30, 40), (50, 60), (1, 30, 30)) == (0, 0, 9, 1, 60, 70)


def test_first_maximum_is_first_in_matlabs_order(be):
    """at a tie the maximum is the first in MATLAB's linear order (y fastest, then x, then z), not
    the first in C order"""
    T = np.zeros((3, 4, 4), dtype=np.float32)
    T[1, 1:3, 1:3] = 1
    A = np.zeros((6, 9, 9), dtype=np.float32)
    A[1, 2:4, 5:7] = 1  # at offset (0, 1, 4): first in C order
    A[1, 6:8, 1:3] = 1  # at offset (0, 5, 0): first in MATLAB's order (smaller x)
    off, mx, C = xcorr.normxcorr3_max_shift(T, A, (5, 8, 8), backend=be)
    box = xcorr.max_shift_box(T.shape, A.shape, (5, 8, 8))
    ties = np.argwhere(C == C.max())
    first = min(ties.tolist(), key=lambda i: (i[0], i[2], i[1]))
    assert off == tuple(int(i + s + 1 - n) for i, s, n in zip(first, box[:3], T.shape, strict=True))
    assert mx == C.max()


def test_normxcorr2_max_shift_window_and_first_maximum():
    """normxcorr2_max_shift around a given normxcorr2: the window of the bounds along y and x, the
    first maximum in MATLAB's order (NaN ignored), the offset (0, y, x)"""
    T = np.zeros((4, 5))
    A = np.zeros((10, 12))
    C = np.zeros((13, 16))
    C[0, 0] = np.nan
    C[5, 7] = 0.75
    C[6, 3] = 0.75  # smaller x: first in MATLAB's order
    off, mx, Cw = xcorr.normxcorr2_max_shift(T, A, (3, 4), normxcorr2=lambda t, a: C)
    # rows max(1, ceil(-3 + 4)) .. min(floor(3 + 4), 13), columns 1 .. min(4 + 5, 16)
    assert Cw.shape == (7, 9)
    assert mx == 0.75
    assert off == (0, 6 + 1 - 4, 3 + 1 - 5)
