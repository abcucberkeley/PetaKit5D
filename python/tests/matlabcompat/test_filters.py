"""Tests of petakit5d.matlabcompat.filters against MATLAB R2024b: fspecial('gaussian'), convn,
conv2, imfilter, imgaussfilt, imgaussfilt3 and edgetaper, on the options the port's callers use
(unmix_channels_gaussian_block / XR_unmix_channels_frame: imgaussfilt3 of single blocks with a
scalar or [y x z] sigma; the puncta removal blocks: imgaussfilt3 of uint16 blocks, sigma 2.5,
and imgaussfilt of their MIP with sigma 50; the PSF analysis: imgaussfilt3(F, 2) and its
'spatial' form; XR_phase_image_flat_field_correction: imgaussfilt of single images, sigma 100;
RLdecon: edgetaper with fspecial('gaussian', nTapering + 30, nTapering); point detection:
imfilter(vol, k, 'conv', 'same', 'symmetric')).

Oracle cases call the built-ins themselves (black box; inputs from our bit-stable generator in
MATLAB's form, roles map them). Exact unless the case names its tolerance (proposed in
docs/parity/matlabcompat.md).

MATLAB tests whose intent is covered here (their classes belong to other tasks): none is N04's;
UnmixChannelsTest (B02), PunctaRemovalTest (B08), PSFAnalysisTest (B06), FSCAnalysisTest (B05)
and PhaseFlatFieldCorrectionTest (B07) use these filters as their oracles and cite the cases
``ig3.*``, ``ig2.*`` and ``site.*``.
"""

from __future__ import annotations

import numpy as np
import pytest

from petakit5d.core import dtypes
from petakit5d.matlabcompat import filters as F
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, PerAxis, Value, Vol
from petakit5d.testing.compare import tol
from petakit5d.testing.synth import CounterRNG

# ---------------------------------------------------------------------------------------------
# inputs in MATLAB's form
# ---------------------------------------------------------------------------------------------


def uniform(seed: int, shape, cls: str = "double") -> np.ndarray:
    u = CounterRNG(seed, "n04").uniform(shape)
    if cls == "double":
        a = u
    elif cls == "single":
        a = u.astype(np.float32)
    elif cls == "uint16":
        a = dtypes.cast(u * 60000, "uint16")
    elif cls == "uint8":
        a = dtypes.cast(u * 250, "uint8")
    elif cls == "int16":
        a = dtypes.cast(u * 60000 - 30000, "int16")
    elif cls == "uint32":
        a = dtypes.cast(u * 4e9, "uint32")
    elif cls == "int32":
        a = dtypes.cast(u * 4e9 - 2e9, "int32")
    elif cls == "logical":
        a = u > 0.5
    else:
        raise ValueError(cls)
    return np.asfortranarray(np.asarray(a))


def sparse(seed: int, shape, frac: float, cls: str = "double") -> np.ndarray:
    """nonzeros in [0.5, 1.5) at a share frac of the places"""
    r = CounterRNG(seed, "n04-sparse")
    keep = r.uniform(shape) < frac
    v = 0.5 + r.spawn("v").uniform(shape)
    a = np.where(keep, v, 0.0)
    return np.asfortranarray(a.astype(np.float32 if cls == "single" else np.float64))


def kernel(seed: int, shape) -> np.ndarray:
    return np.asfortranarray(0.5 + CounterRNG(seed, "n04-kernel").uniform(shape))


def beads_volume(seed: int, shape, cls: str) -> np.ndarray:
    """a smooth background with sparse bright spots (puncta-like), in MATLAB's form"""
    r = CounterRNG(seed, "n04-beads")
    a = 100 + 50 * r.uniform(shape)
    spots = r.spawn("spots").uniform(shape) > 0.995
    a = a + spots * (1000 + 3000 * r.spawn("amp").uniform(shape))
    return np.asfortranarray(
        np.asarray(dtypes.cast(a, cls))
        if cls not in ("single", "double")
        else a.astype(np.float32 if cls == "single" else np.float64)
    )


# ---------------------------------------------------------------------------------------------
# the Python sides: MATLAB's name-value pairs in, Python's forms through the roles
# ---------------------------------------------------------------------------------------------


def _nv(pairs) -> dict:
    kw = {}
    for k in range(0, len(pairs), 2):
        name = str(pairs[k]).lower()
        v = pairs[k + 1]
        key = {"filtersize": "filter_size", "padding": "padding", "filterdomain": "filter_domain"}[
            name
        ]
        kw[key] = v
    return kw


def py_imgaussfilt3(v, sigma=0.5, *pairs):
    return F.imgaussfilt3(v, sigma, **_nv(pairs))


def py_imgaussfilt(a, sigma=0.5, *pairs):
    return F.imgaussfilt(a, sigma, **_nv(pairs))


def py_imfilter(a, h, *options):
    return F.imfilter(a, h, *options)


def py_convn(a, b, shape="full"):
    return F.convn(a, b, shape)


def py_conv2(*args):
    return F.conv2(*args)


def py_fspecial(kind, *args):
    return F.fspecial(kind, *args)


def py_edgetaper(i, psf):
    return F.edgetaper(i, psf)


def py_rldecon_taper(im, ntaper):
    """RLdecon's softening of the x-y borders, plane by plane (Python's (z, y, x))"""
    k = F.fspecial("gaussian", ntaper + 30, ntaper)
    out = np.empty_like(im)
    for z in range(im.shape[0]):
        out[z] = F.edgetaper(im[z], k)
    return out


# ---------------------------------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------------------------------

CASES: list[Case] = []
SIZE3 = {"a": (37, 41, 19), "b": (67, 45, 31), "c": (100, 64, 24)}
SIZE2 = {"a": (37, 41), "b": (67, 45), "c": (100, 130)}


def add(case: Case) -> None:
    CASES.append(case)


# fspecial('gaussian')
for hs in (3, 5, 7, 11, 15, 31, 40, [5, 7], [7, 5], [1, 9]):
    for s in (0.5, 1.0, 2.0, 3.3, 10.0):
        name = hs if np.isscalar(hs) else f"{hs[0]}x{hs[1]}"
        add(
            Case(
                f"fspecial.{name}.s{s}",
                matlab="fspecial",
                args=["gaussian", hs, s],
                python=py_fspecial,
            )
        )
add(Case("fspecial.default", matlab="fspecial", args=["gaussian"], python=py_fspecial))

# convn and conv2 (2-D and 3-D, every shape, double and single)
_KER = {
    "k33": (3, 3),
    "k54": (5, 4),
    "k17": (1, 7),
    "k61": (6, 1),
    "k333": (3, 3, 3),
    "k452": (4, 5, 2),
    "k115": (1, 1, 5),
}
for cls in ("double", "single"):
    for kn, ks in _KER.items():
        for shp in ("full", "same", "valid"):
            src = "v3" if len(ks) == 3 else "v2"
            add(
                Case(
                    f"convn.{cls}.{kn}.{shp}",
                    matlab="convn",
                    args=[Vol(In(src)), Vol(In("k")), shp],
                    python=py_convn,
                    out_roles=("vol",),
                    inputs={
                        "v3": Value(lambda cls=cls: uniform(31, (23, 29, 19), cls)),
                        "v2": Value(lambda cls=cls: uniform(32, (37, 41), cls)),
                        "k": Value(lambda ks=ks: kernel(33, ks)),
                    },
                )
            )
            if len(ks) == 2:
                add(
                    Case(
                        f"conv2.{cls}.{kn}.{shp}",
                        matlab="conv2",
                        args=[In("v2"), In("k"), shp],
                        python=py_conv2,
                        inputs={
                            "v2": Value(lambda cls=cls: uniform(32, (37, 41), cls)),
                            "k": Value(lambda ks=ks: kernel(33, ks)),
                        },
                    )
                )
    for shp in ("full", "same", "valid"):
        add(
            Case(
                f"conv2uv.{cls}.{shp}",
                matlab="conv2",
                args=[In("u"), In("v"), In("a"), shp],
                python=py_conv2,
                inputs={
                    "u": Value(lambda: kernel(34, (5, 1))),
                    "v": Value(lambda: kernel(35, (1, 7))),
                    "a": Value(lambda cls=cls: uniform(36, (37, 41), cls)),
                },
            )
        )
    add(
        Case(
            f"convn.{cls}.bigkernel.same",
            matlab="convn",
            args=[In("a"), In("k"), "same"],
            python=py_convn,
            inputs={
                "a": Value(lambda cls=cls: uniform(37, (4, 5), cls)),
                "k": Value(lambda: kernel(38, (7, 9))),
            },
        )
    )

# imgaussfilt3: classes, sigmas (MATLAB's [y x z]), paddings, filter sizes, sizes
_SIG3 = {
    "s0.5": 0.5,
    "s1": 1.0,
    "s2": 2.0,
    "s2.5": 2.5,
    "aniso": [1.2, 2.7, 0.9],
    "unmix": [1.0, 1.0, 1.0],
}
for cls in ("single", "double", "uint16"):
    for sn, s in _SIG3.items():
        for size in ("a", "b"):
            add(
                Case(
                    f"ig3.{cls}.{sn}.{size}",
                    matlab="imgaussfilt3",
                    args=[Vol(In("v")), PerAxis(s, "yxz")],
                    python=py_imgaussfilt3,
                    out_roles=("vol",),
                    inputs={"v": Value(lambda cls=cls, size=size: uniform(41, SIZE3[size], cls))},
                )
            )
for cls in ("single", "double", "uint8", "int16", "int32"):
    for pad in ("replicate", "symmetric", "circular", 0, 3.5):
        pn = pad if isinstance(pad, str) else f"v{pad}"
        add(
            Case(
                f"ig3.{cls}.pad.{pn}",
                matlab="imgaussfilt3",
                args=[Vol(In("v")), PerAxis([1.0, 1.5, 2.0], "yxz"), "Padding", pad],
                python=py_imgaussfilt3,
                out_roles=("vol",),
                inputs={"v": Value(lambda cls=cls: uniform(42, SIZE3["c"], cls))},
            )
        )
for fs in ([3, 3, 3], [5, 1, 1], [1, 7, 1], [1, 1, 9], [3, 5, 7], 1):
    fsn = fs if np.isscalar(fs) else "x".join(map(str, fs))
    add(
        Case(
            f"ig3.single.fs{fsn}",
            matlab="imgaussfilt3",
            args=[Vol(In("v")), 1.5, "FilterSize", PerAxis(fs, "yxz"), "FilterDomain", "spatial"],
            python=py_imgaussfilt3,
            out_roles=("vol",),
            inputs={"v": Value(lambda: uniform(43, SIZE3["b"], "single"))},
        )
    )
add(
    Case(
        "ig3.single.frequency",
        matlab="imgaussfilt3",
        args=[Vol(In("v")), 1.5, "FilterDomain", "frequency"],
        python=py_imgaussfilt3,
        out_roles=("vol",),
        inputs={"v": Value(lambda: uniform(44, SIZE3["a"], "single"))},
        compare=tol(reason="proposed: imgaussfilt frequency domain (FFT rounding)", rel_max=1e-6),
    )
)
add(
    Case(
        "ig3.single.sparse",
        matlab="imgaussfilt3",
        args=[Vol(In("v")), 1.0],
        python=py_imgaussfilt3,
        out_roles=("vol",),
        inputs={"v": Value(lambda: sparse(45, SIZE3["b"], 0.02, "single"))},
    )
)
# y sizes of whole vectors of the column pass (16 values single, 8 double) with kernels of 3 mod
# 8 taps (sigma 1/3 as the unmixing task strings write it, 2.5 as in the puncta removal blocks,
# 4.5), tails with a one-tap remainder after two chunks (sigma 2 double, 4 single), and y sizes
# of two 8 KB blocks
_IG3Y = [
    ("single", "s1_3.y48", (48, 40, 24), [0.3333333333, 1.0, 1.0]),
    ("single", "s1_3.y64", (64, 40, 12), [0.3333333333, 1.0, 1.0]),
    ("single", "s2.5.y32", (32, 36, 20), 2.5),
    ("single", "s2.5.y48", (48, 48, 48), 2.5),
    ("uint16", "s2.5.y48", (48, 48, 40), 2.5),
    ("single", "s4.5.y64", (64, 30, 22), 4.5),
    ("single", "s4.y40", (40, 30, 22), 4.0),
    ("double", "s2.y52", (52, 40, 24), 2.0),
    ("single", "s2.5.y2064", (2064, 6, 5), [2.5, 1.0, 1.0]),
    ("double", "s2.y1040", (1040, 6, 5), [2.0, 1.0, 1.0]),
]
for cls, name, shape, s in _IG3Y:
    add(
        Case(
            f"ig3.{cls}.{name}",
            matlab="imgaussfilt3",
            args=[Vol(In("v")), PerAxis(s, "yxz")],
            python=py_imgaussfilt3,
            out_roles=("vol",),
            inputs={"v": Value(lambda cls=cls, shape=shape: uniform(46, shape, cls))},
        )
    )

# imgaussfilt (2-D): classes, sigmas, paddings, filter sizes, the frequency domain
_SIG2 = {
    "s0.5": 0.5,
    "s1": 1.0,
    "s2": 2.0,
    "s3": 3.0,
    "s1.3": 1.3,
    "aniso": [0.8, 1.7],
    "aniso2": [1.2, 2.7],
}
for cls in ("single", "double", "uint16"):
    for sn, s in _SIG2.items():
        add(
            Case(
                f"ig2.{cls}.{sn}",
                matlab="imgaussfilt",
                args=[In("a"), s],
                python=py_imgaussfilt,
                inputs={"a": Value(lambda cls=cls: uniform(51, SIZE2["c"], cls))},
            )
        )
for cls in ("single", "double"):
    for pad in ("replicate", "symmetric", "circular", 0, 3.5):
        pn = pad if isinstance(pad, str) else f"v{pad}"
        add(
            Case(
                f"ig2.{cls}.pad.{pn}",
                matlab="imgaussfilt",
                args=[In("a"), 2.2, "Padding", pad],
                python=py_imgaussfilt,
                inputs={"a": Value(lambda cls=cls: uniform(52, SIZE2["b"], cls))},
            )
        )
    add(
        Case(
            f"ig2.{cls}.frequency",
            matlab="imgaussfilt",
            args=[In("a"), 2.2, "FilterDomain", "frequency"],
            python=py_imgaussfilt,
            inputs={"a": Value(lambda cls=cls: uniform(53, SIZE2["b"], cls))},
            compare=tol(
                reason="proposed: imgaussfilt frequency domain (FFT rounding)", rel_max=1e-6
            ),
        )
    )
for fs in ([5, 1], [1, 5], [3, 7], 9):
    fsn = fs if np.isscalar(fs) else "x".join(map(str, fs))
    add(
        Case(
            f"ig2.double.fs{fsn}",
            matlab="imgaussfilt",
            args=[In("a"), 1.5, "FilterSize", fs],
            python=py_imgaussfilt,
            inputs={"a": Value(lambda: uniform(54, SIZE2["b"], "double"))},
        )
    )
add(
    Case(
        "ig2.single.stack",
        matlab="imgaussfilt",
        args=[Vol(In("v")), 1.8],
        python=py_imgaussfilt,
        out_roles=("vol",),
        inputs={"v": Value(lambda: uniform(55, (40, 30, 3), "single"))},
    )
)

# imfilter: vector kernels along each dimension with every option, general kernels, classes,
# sparse arrays (the plain path)
_OPTS = {
    "default": [],
    "conv": ["conv"],
    "rep": ["replicate"],
    "symconv": ["symmetric", "conv"],
    "circ": ["circular"],
    "full": ["full"],
    "v2.5": [2.5, "conv"],
}
for cls in ("double", "single", "uint16"):
    for kn, ks in (("col9", (9, 1)), ("col5", (5, 1)), ("row7", (1, 7)), ("z5", (1, 1, 5))):
        for on, opts in _OPTS.items():
            add(
                Case(
                    f"imf.{cls}.{kn}.{on}",
                    matlab="imfilter",
                    args=[Vol(In("v")), Vol(In("k")), *opts],
                    python=py_imfilter,
                    out_roles=("vol",),
                    inputs={
                        "v": Value(lambda cls=cls: uniform(61, (64, 40, 6), cls)),
                        "k": Value(lambda ks=ks: kernel(62, ks)),
                    },
                )
            )
for cls in ("double", "single"):
    for frac in (0.02, 0.5):
        add(
            Case(
                f"imf.{cls}.sparse{frac}.col9",
                matlab="imfilter",
                args=[In("a"), In("k")],
                python=py_imfilter,
                inputs={
                    "a": Value(lambda cls=cls, frac=frac: sparse(63, (200, 50), frac, cls)),
                    "k": Value(lambda: kernel(64, (9, 1))),
                },
            )
        )
    for kn, ks in (("k33", (3, 3)), ("k55", (5, 5)), ("k333", (3, 3, 3))):
        add(
            Case(
                f"imf.{cls}.{kn}.symconv",
                matlab="imfilter",
                args=[Vol(In("v")), Vol(In("k")), "conv", "same", "symmetric"],
                python=py_imfilter,
                out_roles=("vol",),
                inputs={
                    "v": Value(lambda cls=cls: uniform(65, (40, 30, 12), cls)),
                    "k": Value(lambda ks=ks: kernel(66, ks)),
                },
                compare=tol(reason="proposed: imfilter general kernels (tap order)", ulp=4),
            )
        )
# the kernel lengths around the switch from the chunked vector code to one chain
for cls in ("double", "single"):
    for m in (17, 19, 21):
        add(
            Case(
                f"imf.{cls}.col{m}.default",
                matlab="imfilter",
                args=[Vol(In("v")), Vol(In("k"))],
                python=py_imfilter,
                out_roles=("vol",),
                inputs={
                    "v": Value(lambda cls=cls: uniform(68, (64, 40, 6), cls)),
                    "k": Value(lambda m=m: kernel(69, (m, 1))),
                },
            )
        )


def exact_count(seed: int, shape, k: int) -> np.ndarray:
    """an array with exactly k nonzeros in [0.5, 1.5)"""
    n = int(np.prod(shape))
    r = CounterRNG(seed, "n04-count")
    order = np.argsort(r.uniform(n), kind="stable")
    a = np.zeros(n)
    a[order[:k]] = 0.5 + r.spawn("v").uniform(k)
    return np.asfortranarray(a.reshape(shape, order="F"))


# the switch to the plain path: fewer than a tenth nonzero values (2,560 values)
for k in (250, 262):
    add(
        Case(
            f"imf.double.nnz{k}.col9",
            matlab="imfilter",
            args=[In("a"), In("k")],
            python=py_imfilter,
            inputs={
                "a": Value(lambda k=k: exact_count(70, (64, 40), k)),
                "k": Value(lambda: kernel(64, (9, 1))),
            },
        )
    )
add(
    Case(
        "imf.logical.k33",
        matlab="imfilter",
        args=[In("a"), In("k")],
        python=py_imfilter,
        inputs={
            "a": Value(lambda: uniform(67, (37, 41), "logical")),
            "k": Value(lambda: np.full((3, 3), 1.0 / 9)),
        },
    )
)
_PADS = {"rep": ["replicate"], "sym": ["symmetric"], "circ": ["circular"], "v2.5": [2.5]}
# explicitly padded columns in 8 KB blocks of padded rows (1,024 double, 2,048 single), the last
# vector of columns of whole vectors (single, 3 mod 8 taps), tails with a one-tap remainder after
# two or more chunks
for cls, n, m in (
    ("double", 1024, 11),
    ("double", 3000, 19),
    ("single", 2048, 19),
    ("single", 2070, 11),
    ("single", 48, 3),
    ("double", 100, 13),
    ("single", 100, 25),
):
    for on, opts in _PADS.items():
        add(
            Case(
                f"imf.{cls}.n{n}.col{m}.{on}",
                matlab="imfilter",
                args=[In("a"), In("k"), *opts],
                python=py_imfilter,
                inputs={
                    "a": Value(lambda cls=cls, n=n: uniform(71, (n, 6), cls)),
                    "k": Value(lambda m=m: kernel(72, (m, 1))),
                },
            )
        )

# edgetaper
for cls in ("double", "single", "uint16"):
    for pn, (hs, s) in (("g15", (15, 3.0)), ("g40", (40, 10.0))):
        add(
            Case(
                f"edgetaper.{cls}.{pn}",
                matlab="edgetaper",
                args=[In("i"), In("p")],
                python=py_edgetaper,
                inputs={
                    "i": Value(lambda cls=cls: uniform(71, (128, 150), cls)),
                    "p": Value(lambda hs=hs, s=s: F.fspecial("gaussian", hs, s)),
                },
                compare=tol(reason="proposed: edgetaper (FFT rounding)", rel_max=1e-6, abs=1),
            )
        )

# errors
_ERR = [
    (
        "err.ig3.sigma0",
        "imgaussfilt3",
        [Vol(In("v")), PerAxis([1, 1, 0], "yxz")],
        "MATLAB:imgaussfilt3:expectedPositive",
    ),
    ("err.ig2.sigmaneg", "imgaussfilt", [In("a"), -1], "MATLAB:imgaussfilt:expectedPositive"),
    ("err.ig3.sigma2", "imgaussfilt3", [Vol(In("v")), [1, 2]], "images:imgaussfilt:invalidLength3"),
    ("err.ig2.sigma3", "imgaussfilt", [In("a"), [1, 2, 3]], "images:imgaussfilt:invalidLength"),
    (
        "err.ig3.even",
        "imgaussfilt3",
        [Vol(In("v")), 1, "FilterSize", 4],
        "MATLAB:validateThreeDFilterSize:expectedOdd",
    ),
    (
        "err.ig2.even",
        "imgaussfilt",
        [In("a"), 1, "FilterSize", 4],
        "MATLAB:validateTwoDFilterSize:expectedOdd",
    ),
    ("err.ig2.logical", "imgaussfilt", [In("l"), 1], "MATLAB:imgaussfilt:invalidType"),
    ("err.imf.singlekernel", "imfilter", [In("a"), In("ks")], "MATLAB:conv2:NonDoubleFilter"),
    (
        "err.edgetaper.big",
        "edgetaper",
        [In("a"), In("p40")],
        "images:edgetaper:psfMustBeSmallerThanHalfImage",
    ),
]
for cid, fn, args, ident in _ERR:
    py = {
        "imgaussfilt3": py_imgaussfilt3,
        "imgaussfilt": py_imgaussfilt,
        "imfilter": py_imfilter,
        "edgetaper": py_edgetaper,
    }[fn]
    add(
        Case(
            cid,
            matlab=fn,
            args=args,
            python=py,
            expect_error=ident if ident else True,
            py_error=(F.MatlabError,),
            inputs={
                "v": Value(lambda: uniform(81, (10, 12, 8), "single")),
                "a": Value(lambda: uniform(82, (37, 41), "double")),
                "l": Value(lambda: uniform(83, (12, 10), "logical")),
                "ks": Value(lambda: np.ones((3, 3), dtype=np.float32)),
                "p40": Value(lambda: F.fspecial("gaussian", 40, 10.0)),
            },
        )
    )

# the call sites' options on typical data
add(
    Case(
        "site.unmix.single.s1",
        matlab="imgaussfilt3",
        args=[Vol(In("v")), PerAxis([1.0, 1.0, 1.0], "yxz")],
        python=py_imgaussfilt3,
        out_roles=("vol",),
        inputs={"v": Value(lambda: beads_volume(91, (160, 130, 40), "single"))},
    )
)
add(
    Case(
        "site.unmix.single.aniso",
        matlab="imgaussfilt3",
        args=[Vol(In("v")), PerAxis([1.5, 1.5, 0.8], "yxz")],
        python=py_imgaussfilt3,
        out_roles=("vol",),
        inputs={"v": Value(lambda: beads_volume(92, (130, 160, 40), "single"))},
    )
)
add(
    Case(
        "site.puncta.uint16.s2.5",
        matlab="imgaussfilt3",
        args=[Vol(In("v")), 2.5],
        python=py_imgaussfilt3,
        out_roles=("vol",),
        inputs={"v": Value(lambda: beads_volume(93, (140, 150, 36), "uint16"))},
    )
)
add(
    Case(
        "site.puncta.mip.s50",
        matlab="imgaussfilt",
        args=[In("m"), 50.0],
        python=py_imgaussfilt,
        inputs={"m": Value(lambda: beads_volume(94, (300, 340), "uint16"))},
    )
)
add(
    Case(
        "site.psf.double.s2",
        matlab="imgaussfilt3",
        args=[Vol(In("v")), 2.0],
        python=py_imgaussfilt3,
        out_roles=("vol",),
        inputs={"v": Value(lambda: beads_volume(95, (65, 65, 81), "double"))},
    )
)
add(
    Case(
        "site.psf.double.s2.spatial",
        matlab="imgaussfilt3",
        args=[Vol(In("v")), 2.0, "FilterDomain", "spatial"],
        python=py_imgaussfilt3,
        out_roles=("vol",),
        inputs={"v": Value(lambda: beads_volume(96, (17, 65, 21), "double"))},
    )
)
add(
    Case(
        "site.phaseff.single.s100.spatial",
        matlab="imgaussfilt",
        args=[In("a"), 100.0],
        python=py_imgaussfilt,
        inputs={"a": Value(lambda: uniform(97, (512, 600), "single"))},
    )
)
add(
    Case(
        "site.phaseff.single.s100.frequency",
        matlab="imgaussfilt",
        args=[In("a"), 100.0],
        python=py_imgaussfilt,
        inputs={"a": Value(lambda: uniform(98, (1500, 1500), "single"))},
        compare=tol(reason="proposed: imgaussfilt frequency domain (FFT rounding)", rel_max=1e-6),
        marks=("slow",),
    )
)
for cls in ("uint16", "single"):
    add(
        Case(
            f"site.rldecon.taper10.{cls}",
            code=(
                "im = in{1}; k = fspecial('gaussian', in{2} + 30, in{2});\n"
                "for iz = 1 : size(im, 3), im(:, :, iz) = edgetaper(im(:, :, iz), k); end\n"
                "out = {im};"
            ),
            args=[Vol(In("v")), 10],
            python=py_rldecon_taper,
            out_roles=("vol",),
            inputs={"v": Value(lambda cls=cls: beads_volume(99, (100, 120, 4), cls))},
            compare=tol(reason="proposed: edgetaper (FFT rounding)", rel_max=1e-6, abs=1),
        )
    )
add(
    Case(
        "site.pointdet.imfilter.symconv",
        matlab="imfilter",
        args=[Vol(In("v")), Vol(In("k")), "conv", "same", "symmetric"],
        python=py_imfilter,
        out_roles=("vol",),
        inputs={
            "v": Value(lambda: beads_volume(100, (48, 40, 20), "double")),
            "k": Value(lambda: kernel(101, (5, 5, 3))),
        },
        compare=tol(reason="proposed: imfilter general kernels (tap order)", ulp=4),
    )
)


# Measured differences (no exact model yet), proposed tolerances (docs/parity/matlabcompat.md,
# for the orchestrator): the largest distance in ulp and the share of values that differ, each
# with a margin over the measurement on these cases.
_COLUMN_TAIL = (
    "proposed: imfilter's column pass with zero padding or 'full', rows of the last vectors "
    "(double 7 taps or more, single 11 or more; 'full' from 5)"
)
_SMALL_2D = (
    "proposed: imgaussfilt's small (< 49 values) 2-D kernels, imfilter's 2-D vector code "
    "(tap order not measured)"
)
_GENERAL = "proposed: imfilter with general 2-D and 3-D kernels (tap order not measured)"
_LONG = (
    "proposed: conv2 with long kernels (imgaussfilt's separable path, kernels of 21 or more taps)"
)
_SPARSE = "proposed: imfilter's plain path for sparse arrays, borders of padded arrays"
_INT = (
    "proposed: integer results one count off where the double result is within rounding of a half"
)
TOLERANCES = {
    "imf.double.col5.full": tol(reason=_COLUMN_TAIL, ulp=2, max_frac=0.01),
    "imf.double.col9.full": tol(reason=_COLUMN_TAIL, ulp=2, max_frac=0.01),
    "imf.single.col5.full": tol(reason=_COLUMN_TAIL, ulp=2, max_frac=0.02),
    "imf.single.col9.full": tol(reason=_COLUMN_TAIL, ulp=3, max_frac=0.04),
    "ig3.single.sparse": tol(reason=_SPARSE, ulp=4, max_frac=0.08),
    "imf.double.col17.default": tol(reason=_COLUMN_TAIL, ulp=3, max_frac=0.05),
    "imf.double.col19.default": tol(reason=_COLUMN_TAIL, ulp=3, max_frac=0.07),
    "imf.double.col21.default": tol(reason=_COLUMN_TAIL, ulp=3, max_frac=0.07),
    "imf.single.col19.default": tol(reason=_COLUMN_TAIL, ulp=5, max_frac=0.15),
    "imf.single.col21.default": tol(reason=_COLUMN_TAIL, ulp=5, max_frac=0.15),
    "site.phaseff.single.s100.spatial": tol(reason=_LONG, ulp=32, rel_max=2e-6),
    "site.puncta.mip.s50": tol(reason=_LONG, abs=1, max_frac=1e-4),
    "site.pointdet.imfilter.symconv": tol(reason=_GENERAL, ulp=16),
}
for _cid in (
    "ig2.double.aniso",
    "ig2.double.fs3x7",
    "ig2.double.s0.5",
    "ig2.double.s1",
    "ig2.single.aniso",
    "ig2.single.s0.5",
    "ig2.single.s1",
):
    TOLERANCES[_cid] = tol(reason=_SMALL_2D, ulp=12)
for _cid in ("ig2.uint16.aniso", "ig2.uint16.s0.5", "ig2.uint16.s1"):
    TOLERANCES[_cid] = tol(reason=_SMALL_2D, abs=1, max_frac=0.006)
for _cid in (
    "imf.double.k55.symconv",
    "imf.double.k333.symconv",
    "imf.single.k55.symconv",
    "imf.single.k333.symconv",
):
    TOLERANCES[_cid] = tol(reason=_GENERAL, ulp=10)
for _cid in (
    "imf.uint16.row7.full",
    "imf.uint16.z5.circ",
    "imf.uint16.z5.conv",
    "imf.uint16.z5.default",
    "imf.uint16.z5.full",
    "imf.uint16.z5.rep",
    "imf.uint16.z5.v2.5",
):
    TOLERANCES[_cid] = tol(reason=_INT, abs=1, max_frac=0.001)
for _c in CASES:
    if _c.id in TOLERANCES:
        _c.compare = TOLERANCES[_c.id]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# unit tests
# ---------------------------------------------------------------------------------------------


def _rand(shape, dt, seed=1):
    return CounterRNG(seed, "n04-unit").uniform(shape).astype(dt)


@pytest.mark.parametrize("dt", [np.float32, np.float64])
def test_results_do_not_depend_on_the_thread_count(dt):
    v = _rand((23, 70, 41), dt)
    a = F.imgaussfilt3(v, (1.0, 2.5, 1.5), num_threads=1)
    for nt in (2, 3, 7):
        assert np.array_equal(F.imgaussfilt3(v, (1.0, 2.5, 1.5), num_threads=nt), a)
    img = _rand((120, 90), dt)
    b = F.imgaussfilt(img, 3.0, num_threads=1)
    assert np.array_equal(F.imgaussfilt(img, 3.0, num_threads=5), b)
    k = _rand((3, 4, 5), np.float64)
    c = F.convn(v, k, "same", num_threads=1)
    assert np.array_equal(F.convn(v, k, "same", num_threads=4), c)


def _fma1(wj, xj, acc, dt):
    """fma(wj, xj, acc) by exact arithmetic"""
    from fractions import Fraction

    flat = [
        float(Fraction(float(wj)) * Fraction(float(v)) + Fraction(float(a)))
        for v, a in zip(xj.ravel(), acc.ravel(), strict=True)
    ]
    return np.array(flat, dtype=np.float64).astype(dt).reshape(acc.shape)


def _chain(x, w, order, dt, fused):
    acc = np.zeros(x.shape[1:], dtype=dt)
    for j in order:
        if fused:
            acc = _fma1(w[j], x[j], acc, dt)
        else:
            acc = (acc + (x[j] * dt(w[j])).astype(dt)).astype(dt)
    return acc


@pytest.mark.parametrize("dt", [np.float32, np.float64])
def test_column_pass_follows_its_row_programs(dt):
    """the kernel's chains equal a reference chain row by row (body, tail and border rows)"""
    n, m = 37, 9
    x = _rand((n, 5), dt, 3)
    w = _rand((m,), dt, 4) + dt(0.5)
    rows = F._column_rows(n, m, (m - 1) // 2, dt == np.float64, True)
    y = F._line_pass(x, w, 0, rows, (m - 1) // 2, n, 0, 0.0, 2)
    xp = np.pad(x, (((m - 1) // 2, m - 1 - (m - 1) // 2), (0, 0)))
    for i in (0, 5, 20, 31, 33, 36):
        ref = _chain(xp[i : i + m], w, rows[i], dt, True)
        assert np.array_equal(y[i], ref), i
    # the programs: body rows in chunks of 2H interleaved, the vector tail with the top chunk
    # ascending, the border rows descending
    vec = 8 if dt == np.float64 else 16
    half = 2 if dt == np.float64 else 4
    body = F._body_order(m, half)
    c = min(m, 2 * half)
    assert rows[0] == body
    assert rows[n - 1] == sorted(body[:c], reverse=True) + body[c:]
    assert rows[vec * ((n - 4) // vec)] != rows[0]
    padded = F._column_rows(n, m, (m - 1) // 2, dt == np.float64, False)
    assert padded[n - 1] == F._padded_tail(m, 2 * half) and padded[0] == body


def _program_ref(x, w, prog, dt):
    """a tap program: fused chains split by -1 and added in order; the taps after a -3 continue
    from the sum so far"""
    total = acc = None
    for op in [*prog, -1]:
        if op in (-1, -3):
            if acc is not None:
                total = acc if total is None else (total + acc).astype(dt)
            acc = None
            if op == -3:
                acc, total = total, None
            continue
        start = np.zeros(x.shape[1:], dtype=dt) if acc is None else acc
        acc = _fma1(w[op], x[op], start, dt)
    return total


@pytest.mark.parametrize("dt", [np.float32, np.float64])
def test_padded_columns_in_blocks(dt):
    """an explicitly padded column in blocks of 8 KB of padded rows: each block's whole vectors
    in body order, then its tail; single's last vector with a three-tap remainder; a one-tap
    remainder fused last into the sum of two or more chunks"""
    double = dt == np.float64
    vec, half, block = (8, 2, 1024) if double else (16, 4, 2048)
    body, tail = F._body_order(11, half), F._padded_tail(11, 2 * half)
    rows = F._column_rows(block + 30, 11, 5, double, False)
    first = block - 10
    nb = vec * (first // vec)
    assert rows[:nb] == [body] * nb and rows[nb:first] == [tail] * (first - nb)
    second = len(rows) - first
    nb2 = vec * (second // vec)
    assert rows[first : first + nb2] == [body] * nb2
    assert rows[first + nb2 :] == [tail] * (second - nb2)
    whole = F._column_rows(4 * vec, 11, 5, double, False)
    if double:
        assert whole == [body] * (4 * vec)
    else:
        assert whole[: 3 * vec] == [body] * (3 * vec)
        assert whole[-1] == [10, 6, 9, 5, 8, 4, 7, 3, -1, 0, 1, 2]
    m = 4 * half + 1
    top, low = list(range(1 + 2 * half, m)), list(range(1, 1 + 2 * half))
    assert F._padded_tail(m, 2 * half) == [*top, -1, *low, -3, 0]
    # the kernel runs those programs as the reference does
    n = 2 * vec + 3
    x = _rand((n, 4), dt, 7)
    w = _rand((m,), dt, 8) + dt(0.5)
    rows = F._column_rows(n, m, (m - 1) // 2, double, False)
    y = F._line_pass(x, w, 0, rows, (m - 1) // 2, n, 1, 0.0, 2)
    xp = np.pad(x, (((m - 1) // 2, m - 1 - (m - 1) // 2), (0, 0)), mode="edge")
    for i in (0, n - 3, n - 1):
        assert np.array_equal(y[i], _program_ref(xp[i : i + m], w, rows[i], dt)), i


def test_body_order_interleaves_the_halves_of_each_chunk():
    assert F._body_order(9, 2) == [8, 6, 7, 5, 4, 2, 3, 1, 0]
    assert F._body_order(11, 2) == [10, 8, 9, 7, 6, 4, 5, 3, 2, 0, 1]
    assert F._body_order(9, 4) == [8, 4, 7, 3, 6, 2, 5, 1, 0]
    assert F._body_order(15, 4) == [14, 10, 13, 9, 12, 8, 11, 7, 6, 2, 5, 1, 4, 0, 3]
    assert F._body_order(5, 4) == [4, 0, 3, 2, 1]


def test_plain_path_for_sparse_arrays():
    """fewer than a tenth nonzero values: products rounded, added in the correlation's order"""
    a = np.zeros((64, 40))
    rng = CounterRNG(5, "n04-unit-sparse")
    idx = rng.integers(0, a.size - 1, 200)
    a.reshape(-1)[idx] = 0.5 + rng.spawn("v").uniform(200)
    w = 0.5 + _rand((9, 1), np.float64, 6)
    y = F.imfilter(a, w)
    xp = np.pad(a, ((4, 4), (0, 0)))
    ref = np.zeros_like(a)
    for j in range(9):
        ref = ref + xp[j : j + 64] * w[j, 0]
    assert np.array_equal(y, ref)
    assert F._sparse(a, 9) and not F._sparse(a, 17)
    b = np.zeros(2560)
    b[:255] = 1.0
    assert F._sparse(b, 9)
    b[:270] = 1.0
    assert not F._sparse(b, 9)
    dense = a.copy()
    dense.reshape(-1)[: a.size // 10 + 1] = 1.0
    assert not F._sparse(dense, 9)


def test_gaussian_kernels():
    g = F.gaussian_kernel_1d(2.0, 9)
    assert g.shape == (9,) and np.array_equal(g, g[::-1])
    assert abs(g.sum() - 1) < 1e-15
    h = F.fspecial("gaussian", 7, 1.5)
    assert h.shape == (7, 7) and np.array_equal(h, h.T)
    assert np.array_equal(F.fspecial("gaussian", [5, 7], 2.0).shape, (5, 7))
    with pytest.raises(F.MatlabError):
        F.fspecial("disk", 3)


def test_sizes_classes_and_errors():
    v = _rand((5, 6, 7), np.float32)
    assert F.imgaussfilt3(v, 1).dtype == np.float32
    u = dtypes.cast(v * 1000, "uint16")
    assert F.imgaussfilt3(u, 1).dtype == np.uint16
    assert F.imgaussfilt(np.zeros((0, 4), np.float32), 1).shape == (0, 4)
    for bad, ident in (
        (lambda: F.imgaussfilt3(v, (1, 0, 1)), "MATLAB:imgaussfilt3:expectedPositive"),
        (lambda: F.imgaussfilt(v[0] > 0.5, 1), "MATLAB:imgaussfilt:invalidType"),
        (lambda: F.imgaussfilt(v[0].astype(np.complex64), 1), "MATLAB:imgaussfilt:expectedReal"),
        (
            lambda: F.imgaussfilt3(v, 1, filter_size=4),
            "MATLAB:validateThreeDFilterSize:expectedOdd",
        ),
        (lambda: F.imfilter(v, np.ones((3, 1), np.float32)), "MATLAB:conv2:NonDoubleFilter"),
    ):
        with pytest.raises(F.MatlabError) as e:
            bad()
        assert e.value.identifier == ident
    assert F._match("Symm", ("replicate", "symmetric", "circular")) == "symmetric"


def test_auto_domain_rule():
    assert not F._auto_frequency_2d((1448, 1448), (401, 401), "single")
    assert F._auto_frequency_2d((1500, 1500), (121, 121), "single")
    assert not F._auto_frequency_2d((1500, 1500), (101, 101), "single")
    assert not F._auto_frequency_2d((1500, 1500), (181, 181), "double")
    assert F._auto_frequency_2d((1500, 1500), (201, 201), "uint16")


@pytest.mark.parametrize("isa", [0, 1])
def test_instruction_sets_give_the_same_bits(isa):
    ext = F._ext()
    if isa > ext.best_isa():
        pytest.skip("this CPU lacks AVX2 with FMA")
    x = _rand((3, 50, 40), np.float32, 7)
    w = (0.5 + _rand((9,), np.float32, 8)).astype(np.float32)
    rows = F._column_rows(50, 9, 4, False, True)
    prog, rowmap = F._programs(rows)
    out = np.empty_like(x)
    ext.line_program_f32(x, out, 3, 50, 40, 50, w, prog, rowmap, 4, 9, 0, 0.0, 2, isa)
    ref = F._line_pass(x, w, 1, rows, 4, 50, 0, 0.0, 2)
    assert np.array_equal(out, ref)
