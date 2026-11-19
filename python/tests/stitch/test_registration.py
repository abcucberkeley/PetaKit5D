"""Tests of petakit5d.stitch.registration: the registration of stitching tile pairs (task S04).

Oracle cases on synthetic tile pairs, written in MATLAB's layouts (zarr in MATLAB's order, tiff),
so both sides read the same files:

- r3.<k>: cross_correlation_registration_3d on pairs of textured, deskewed-like (empty
  wedges), sparse, offset low-contrast and noise tiles, uint16, uint8 and single, tiff and zarr,
  in all six data orders, overlapping along x, y or z (some also offset in another axis), the
  second tile listed off its position (within the max offsets) and its stage position off the
  voxel grid, with and without max pooling (downSample), several max offsets and voxel sizes,
  either tile leading; flat ones: a flat template (no NCC: the first allowed shift), a flat image,
  an empty template (no registration: [0 0 0] and 1); crop.<k>: long templates cropped by
  intensity (dimNumThrsh, blankNumTrsh); axes.<order>: RegistrationXcorrTest's
  pairRegistrationTakesTheAxesOfTheDataOrder in each of its data orders;
- r2.<k>: cross_correlation_registration_2d on 2D tiff pairs ('y,x,z' and 'x,y,z'), textured,
  sparse, with empty rows and columns (the template's masks), an empty template;
- slab.<k>: cross_correlation_registration_3d_mip_slabs on the three MIP slabs of tile pairs
  (pool sizes of 3 and 6 values; its per-axis registrations run as tasks in local mode on both
  sides);
- wrap.<k>: cross_correlation_registration_wrapper's result file for a leading tile with two and
  three pairs (3D, 2D, large files);
- err.*: a tile of one slice (MATLAB's two-value size), an unknown extension, a constant 2D
  template (normxcorr2's error), a 2D tile whose z is not its third axis.

3D cases run on both CPU FFT backends: exact (the FFTW add-on): bit for bit (shifts and peaks);
permissive (the default): the same shifts, the peaks within decision record 0006's 5e-3. 2D cases
compare the shifts exactly and the peaks within the proposed normxcorr2 tolerance
(docs/parity/stitch.md).

MATLAB tests ported here: RegistrationXcorrTest/pairRegistrationTakesTheAxesOfTheDataOrder
(axes.*, test_pair_registration_takes_the_axes_of_the_data_order); StitchingTest's registration
methods (S09's class) at the level of the registration: xcorrRecoversPerturbedTilePosition,
localShiftMethod, registration2DTiles, skewedSpaceRegistration, largeFileStitching and
largeFileRegistrationRecoversPerturbedTilePosition (test_registration_recovers_*, the
test_stitching_* functions: the stitching wrapper with the port's registration).
"""

from __future__ import annotations

import dataclasses
import functools
import os
from pathlib import Path

import numpy as np
import pytest

from petakit5d import fft as pfft
from petakit5d.core.arrays import from_matlab, to_matlab
from petakit5d.core.axes import from_xyz, from_yxz
from petakit5d.io.matfile import read_mat
from petakit5d.stitch import registration as sr
from petakit5d.stitch.geometry import check_cuboids_overlaps, data_order_axes_from_matlab
from petakit5d.stitch.xcorr import RegistrationError
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, Index1, InPath, PerAxis, Tree, Work
from petakit5d.testing.compare import EXACT, tol, tree
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG

REASON = "decision record 0006: the permissive FFT backend against MATLAB's FFTW"
NCC_TOL = tol(reason=REASON + ": the NCC's peak", abs=5e-3)
N2_REASON = (
    "docs/parity/stitch.md (proposed): normxcorr2, MathWorks' function, a new implementation"
)
# proposed: 2D peaks within 1e-7 of MATLAB's (measured: 5.3e-8 on offset data, where it is
# MATLAB's own error; 1e-13 elsewhere)
N2_PEAK = tol(reason=N2_REASON, abs=1e-7)

# MATLAB's data orders (the names of the tile arrays' dimensions)
ORDERS = ("y,x,z", "x,y,z", "z,y,x", "z,x,y", "y,z,x", "x,z,y")


def dom_of(order: str) -> list[int]:
    """MATLAB's data_order_mat: [~, dom] = sort(strrep(order, ',', ''))"""
    s = order.replace(",", "")
    return [s.index(c) + 1 for c in "xyz"]


def perm_of(order: str) -> list[int]:
    """the permutation of a [y x z] tile into the order's dimensions:
    [~, perm] = sort(dom([2, 1, 3]))"""
    d = dom_of(order)
    v = [d[1], d[0], d[2]]
    return [int(k) + 1 for k in np.argsort(v, kind="stable")]


# ---------------------------------------------------------------------------------------------
# synthetic tiles (MATLAB's form: [y x z] arrays of the ground truth)
# ---------------------------------------------------------------------------------------------


def content(kind: str, shape_yxz, seed: int) -> np.ndarray:
    """a ground truth volume [y x z] (float64): smooth texture, with deskewed-like empty wedges,
    sparse blobs, low contrast on an offset, noise, or one constant"""
    ny, nx, nz = shape_yxz
    zyx = (nz, ny, nx)
    rng = CounterRNG(seed, f"s04/{kind}")
    if kind in ("texture", "dsr", "holes"):
        v = synth.texture(zyx, 2.0, (100.0, 3000.0), rng=rng).astype(np.float64)
        if kind == "dsr":
            for z in range(nz):
                a = (nz - 1 - z) * (nx // 3) // max(nz - 1, 1)
                v[z, :, :a] = 0
                b = z * (nx // 4) // max(nz - 1, 1)
                if b:
                    v[z, : ny // 2, nx - b :] = 0
        elif kind == "holes":
            v[:, ny // 3 : ny // 3 + 3, :] = 0
            v[:, :, nx // 2 : nx // 2 + 2] = 0
    elif kind == "sparse":
        k = max(4, (ny * nx * nz) // 4000)
        pos = np.stack(
            [
                rng.uniform((k,)) * (nz - 1),
                rng.uniform((k,)) * (ny - 1),
                rng.uniform((k,)) * (nx - 1),
            ],
            axis=1,
        )
        v = synth.beads(zyx, pos, (1.5, 2.0), 2000.0, 0.0).astype(np.float64)
    elif kind == "offset":
        v = mround_(synth.texture(zyx, 2.0, (2.0**15 - 6, 2.0**15 + 6), rng=rng))
    elif kind == "noise":
        v = rng.uniform(zyx) * 1000.0
    elif kind == "flat":
        v = np.full(zyx, 500.0)
    elif kind == "zero":
        v = np.zeros(zyx)
    else:
        raise KeyError(kind)
    return to_matlab(np.ascontiguousarray(v))


def mround_(v):
    from petakit5d.core.dtypes import mround

    return mround(np.asarray(v, dtype=np.float64))


def as_class(v: np.ndarray, cls: str) -> np.ndarray:
    from petakit5d.core.dtypes import cast

    if cls == "uint8":
        return np.asarray(cast(v / 12.0, "uint8"))
    return np.asarray(cast(v, cls)) if cls == "uint16" else v.astype(np.float32)


@dataclasses.dataclass(frozen=True)
class Pair:
    """a tile pair: sizes [y x z] (MATLAB's), the overlap axis (0 x, 1 y, 2 z) and length, the
    true offset of tile 2 in its other axes, its listed position's error p (voxels, x y z), the
    stage jitter (voxels), the content and class, the data order, the file type, which tile
    leads, voxel sizes (x y z), max offsets (x y z), downSample ([d1 d2 d3]), the crop limits"""

    size: tuple[int, int, int]
    axis: int
    ov: int
    lateral: tuple[int, int, int] = (0, 0, 0)
    p: tuple[int, int, int] = (0, 0, 0)
    jitter: tuple[float, float, float] = (0.0, 0.0, 0.0)
    kind: str = "texture"
    cls: str = "uint16"
    order: str = "y,x,z"
    ext: str = ".zarr"
    swap: bool = False
    vs: tuple[float, float, float] = (1.0, 1.0, 1.0)
    mo: tuple[float, float, float] = (10.0, 10.0, 4.0)
    ds: tuple[float, float, float] = (1.0, 1.0, 1.0)
    thr: float = 10000
    blank: float | None = None
    kind2: str | None = None  # another content for tile 1 (a flat image)
    seed: int = 0
    two_d: bool = False

    @property
    def offset2(self) -> tuple[int, int, int]:
        """tile 2's true origin (x, y, z voxels) in the ground truth"""
        ny, nx, nz = self.size
        n_xyz = (nx, ny, nz)
        o = list(self.lateral)
        o[self.axis] = n_xyz[self.axis] - self.ov
        return tuple(o)


def tiles(pr: Pair, gsize=None):
    """the two tiles ([y x z] arrays of their class, before the data order) of a ground truth of
    gsize ([y x z]; default: just large enough) and their listed origins (x, y, z voxels)"""
    ny, nx, nz = pr.size
    o2 = pr.offset2
    if gsize is None:
        gsize = (ny + o2[1], nx + o2[0], nz + o2[2])
    g = content(pr.kind, gsize, pr.seed)
    t1 = g[:ny, :nx, :nz]
    t2 = g[o2[1] : o2[1] + ny, o2[0] : o2[0] + nx, o2[2] : o2[2] + nz]
    if pr.kind2 is not None:
        t1 = content(pr.kind2, (ny, nx, nz), pr.seed + 1)
    t1 = as_class(np.ascontiguousarray(t1), pr.cls)
    t2 = as_class(np.ascontiguousarray(t2), pr.cls)
    listed2 = tuple(o2[k] + pr.p[k] for k in range(3))
    return t1, t2, (0, 0, 0), listed2


def cuboids(pr: Pair):
    """the tiles' cuboids [x0 y0 z0 x1 y1 z1] (um) at their listed (jittered) positions, and the
    overlap; tile 1 is the leading one unless swap"""
    ny, nx, nz = pr.size
    _, _, l1, l2 = tiles(pr)
    vs = np.asarray(pr.vs)
    n = np.asarray([nx, ny, nz], dtype=np.float64)
    c = []
    for lst, jit in ((l1, (0.0, 0.0, 0.0)), (l2, pr.jitter)):
        x0 = (np.asarray(lst, dtype=np.float64) + np.asarray(jit)) * vs
        c.append(np.concatenate([x0, x0 + (n - 1) * vs]))
    if pr.swap:
        c = c[::-1]
    is_ov, ov = check_cuboids_overlaps(c[0], c[1], pr.two_d)
    if not np.asarray(is_ov).ravel()[0]:
        raise ValueError(f"no overlap: {pr}")
    return [float(v) for v in c[0]], [float(v) for v in c[1]], [float(v) for v in ov.ravel()]


def write_tile(path: Path, a_yxz: np.ndarray, order: str, ext: str) -> None:
    f = np.transpose(a_yxz, [k - 1 for k in perm_of(order)])
    if ext == ".zarr":
        synth.write_zarr_matlab(path, f, (16, 16, 16))
    else:
        synth.write_tiff_matlab(path, np.ascontiguousarray(f))


def pair_tree(pr: Pair):
    def fill(folder):
        folder = Path(folder)
        t1, t2, _, _ = tiles(pr)
        if pr.swap:
            t1, t2 = t2, t1
        write_tile(folder / f"t1{pr.ext}", t1, pr.order, pr.ext)
        write_tile(folder / f"t2{pr.ext}", t2, pr.order, pr.ext)

    return fill


# ---------------------------------------------------------------------------------------------
# the Python sides
# ---------------------------------------------------------------------------------------------


def py_reg3d(d, ext, c1, c2, co, vs, dom, ds, mo, thr, blank, *, backend=None, out=""):
    r = sr.cross_correlation_registration_3d(
        f"{d}t1{ext}",
        f"{d}t2{ext}",
        out,
        c1,
        c2,
        co,
        vs,
        data_order_axes_from_matlab(dom),
        down_sample=ds,
        max_offset=mo,
        dim_num_thrsh=thr,
        blank_num_trsh=blank,
        backend=backend,
    )
    return [] if r is None else list(r)


def py_reg2d(d, ext, c1, c2, co, vs, dom, ds, mo):
    return list(
        sr.cross_correlation_registration_2d(
            f"{d}t1{ext}",
            f"{d}t2{ext}",
            "",
            c1,
            c2,
            co,
            vs,
            data_order_axes_from_matlab(dom),
            down_sample=ds,
            max_offset=mo,
        )
    )


def py_slabs(d, out, c1, c2, co, vs, dom, ds, mo, mip, pool, *, backend=None):
    with pfft.using(backend or pfft.backend()):
        rel, mx, rsm, mxm = sr.cross_correlation_registration_3d_mip_slabs(
            f"{d}t1.zarr",
            f"{d}t2.zarr",
            out,
            c1,
            c2,
            co,
            vs,
            data_order_axes_from_matlab(dom),
            down_sample=ds,
            max_offset=mo,
            mip_dir_str=mip,
            pool_size=pool,
            parse_cluster=False,
        )
    return [rel, mx, rsm, mxm]


def py_wrapper(d, ext, n, out, pairs, c1, c2, co, vs, dom, *, backend=None, **kw):
    with pfft.using(backend or pfft.backend()):
        sr.cross_correlation_registration(
            f"{d}t1{ext}",
            [f"{d}t{k + 2}{ext}" for k in range(n)],
            out,
            pairs,
            c1,
            c2,
            co,
            vs,
            data_order_axes_from_matlab(dom),
            parse_cluster=False,
            **kw,
        )


_REG3D = """
d = in{1}; e = in{2};
[rs, mx] = cross_correlation_registration_3d([d 't1' e], [d 't2' e], '', in{3:end});
out = {rs, mx};
"""
_REG2D = """
d = in{1}; e = in{2};
[rs, mx] = cross_correlation_registration_2d([d 't1' e], [d 't2' e], '', in{3:end});
out = {rs, mx};
"""
_SLABS = """
d = in{1};
[rs, mx, rsm, mxm] = cross_correlation_registration_3d_mip_slabs([d 't1.zarr'], [d 't2.zarr'], ...
    in{2:end});
out = {rs, mx, rsm, mxm};
"""
_WRAP = """
d = in{1}; e = in{2}; n = in{3};
f2 = arrayfun(@(k) sprintf('%st%d%s', d, k + 1, e), 1 : n, 'unif', 0);
cross_correlation_registration_wrapper([d 't1' e], f2, in{4:end});
out = {};
"""


def reg3d_case(cid: str, pr: Pair, note: str = "permissive:shift", **kw) -> Case:
    c1, c2, co = cuboids(pr)
    dom = dom_of(pr.order)
    args = [
        InPath("d"),
        pr.ext,
        c1,
        c2,
        co,
        list(pr.vs),
        dom,
        "downSample",
        list(pr.ds),
        "MaxOffset",
        list(pr.mo),
        "dimNumThrsh",
        pr.thr,
    ]
    if pr.blank is not None:
        args += ["blankNumTrsh", pr.blank]
    return Case(
        cid,
        code=_REG3D,
        args=args,
        nargout=2,
        python=py_reg3d,
        py_args=[
            InPath("d"),
            pr.ext,
            c1,
            c2,
            co,
            list(pr.vs),
            dom,
            PerAxis(list(pr.ds), "yxz"),
            PerAxis(list(pr.mo), "xyz"),
            pr.thr,
            pr.blank,
        ],
        py_kwargs={"backend": "exact"},
        inputs={"d": Tree(pair_tree(pr), "d")},
        files=False,
        note=note,
        **kw,
    )


def reg2d_case(cid: str, pr: Pair, **kw) -> Case:
    c1, c2, co = cuboids(pr)
    dom = dom_of(pr.order)
    return Case(
        cid,
        code=_REG2D,
        args=[
            InPath("d"),
            pr.ext,
            c1,
            c2,
            co,
            list(pr.vs),
            dom,
            "downSample",
            list(pr.ds),
            "MaxOffset",
            list(pr.mo),
        ],
        nargout=2,
        python=py_reg2d,
        py_args=[
            InPath("d"),
            pr.ext,
            c1,
            c2,
            co,
            list(pr.vs),
            dom,
            PerAxis(list(pr.ds), "yxz"),
            PerAxis(list(pr.mo), "xyz"),
        ],
        inputs={"d": Tree(pair_tree(pr), "d")},
        files=False,
        compare={1: N2_PEAK},
        **kw,
    )


# ---------------------------------------------------------------------------------------------
# the cases
# ---------------------------------------------------------------------------------------------

_DS = ((1, 1, 1), (1, 1, 1), (1, 1, 1), (2, 2, 1), (2, 2, 2), (1, 2, 2), (3, 2, 1), (2, 1, 2))
_MO = ((10, 10, 4), (8, 6, 3), (20, 20, 10), (6, 6, 6), (5, 12, 2), (12, 5, 8))
_VS = ((1.0, 1.0, 1.0), (0.108, 0.108, 0.3), (0.5, 0.25, 0.7), (1.0, 1.0, 1.0))
_KINDS = ("texture", "texture", "texture", "dsr", "holes", "sparse", "offset", "noise")
_CLASSES = ("uint16", "uint16", "single", "uint8")


class Picker:
    """deterministic choices (CounterRNG)"""

    def __init__(self, stream: str):
        self.rng = CounterRNG(5, stream)

    def int(self, lo: int, hi: int) -> int:
        return int(self.rng.integers(lo, hi, (1,))[0])

    def of(self, seq):
        return seq[self.int(0, len(seq) - 1)]

    def frac(self, scale: float) -> float:
        return float(np.round(self.rng.uniform((1,))[0] * 2 * scale - scale, 3))


def random_pair(k: int, two_d: bool = False) -> Pair:
    pk = Picker(f"pair{'2' if two_d else '3'}/{k}")
    axis = pk.of((0, 1) if two_d else (0, 1, 2))
    mo = pk.of(_MO)
    ds = pk.of(((1, 1, 1), (1, 1, 1), (2, 2, 1), (1, 2, 1))) if two_d else pk.of(_DS)
    ny, nx = pk.int(28, 52), pk.int(28, 52)
    nz = 1 if two_d else pk.int(16, 34)
    size = (ny, nx, nz)
    n_xyz = (nx, ny, nz)
    ov = pk.int(10, min(18, n_xyz[axis] - 6))
    lateral = [0, 0, 0]
    if pk.int(0, 2) == 0:
        other = pk.of([a for a in ((0, 1) if two_d else (0, 1, 2)) if a != axis])
        lateral[other] = pk.int(1, 3)
    p = [0, 0, 0]
    for a in (0, 1) if two_d else (0, 1, 2):
        m = int(mo[a]) - 1
        p[a] = pk.int(-min(m, 4), min(m, 4))
    jit = (0.0, 0.0, 0.0)
    if pk.int(0, 3) == 0:
        jit = (pk.frac(0.4), pk.frac(0.4), 0.0 if two_d else pk.frac(0.4))
    kind = pk.of(_KINDS[:-1] if two_d else _KINDS)
    cls = pk.of(_CLASSES)
    if kind == "offset":
        cls = "uint16" if two_d else pk.of(("uint16", "single"))
    order = pk.of(("y,x,z", "x,y,z")) if two_d else pk.of(ORDERS)
    ext = ".tif" if two_d else pk.of((".zarr", ".zarr", ".tif"))
    return Pair(
        size=size,
        axis=axis,
        ov=ov,
        lateral=tuple(lateral),
        p=tuple(p),
        jitter=jit,
        kind=kind,
        cls=cls,
        order=order,
        ext=ext,
        swap=pk.int(0, 1) == 1,
        vs=pk.of(_VS),
        mo=mo,
        ds=ds,
        seed=k,
        two_d=two_d,
    )


N_RANDOM_3D = 160
N_RANDOM_2D = 80


def r3_cases():
    out = [reg3d_case(f"r3.{k:03d}", random_pair(k)) for k in range(N_RANDOM_3D)]
    base = Pair(size=(36, 40, 20), axis=0, ov=14, p=(3, -2, 1))
    flats = {
        # a flat template: no NCC anywhere (0), the first allowed shift
        "flat_template": dataclasses.replace(base, kind="flat"),
        # a flat image under a textured template
        "flat_image": dataclasses.replace(base, kind2="flat"),
        "flat_image_single": dataclasses.replace(base, kind2="flat", cls="single"),
        # an empty template: [0 0 0] and 1
        "empty_template": dataclasses.replace(base, kind="zero"),
        "empty_template_tif": dataclasses.replace(base, kind="zero", ext=".tif"),
        # a template with data in less than a quarter of its area
        "sparse_eighth": dataclasses.replace(base, kind="sparse", seed=3),
        "flat_template_pooled": dataclasses.replace(base, kind="flat", ds=(2, 2, 2)),
        "zero_image": dataclasses.replace(base, kind2="zero"),
    }
    out += [reg3d_case(f"r3.{k}", pr) for k, pr in flats.items()]
    # long templates cropped by intensity
    crops = {
        "x": Pair(size=(30, 90, 20), axis=1, ov=12, p=(2, 1, 0), thr=40),
        "y": Pair(size=(90, 30, 20), axis=0, ov=12, p=(1, 2, 1), thr=40),
        "z": Pair(size=(30, 30, 70), axis=0, ov=12, p=(2, 0, 1), thr=30, mo=(8, 8, 6)),
        "blank": Pair(size=(30, 90, 20), axis=1, ov=12, p=(2, 1, 0), thr=40, blank=5, kind="holes"),
        "pooled": Pair(size=(30, 90, 20), axis=1, ov=12, p=(2, 2, 0), thr=40, ds=(2, 2, 1)),
        "two_axes": Pair(size=(80, 90, 20), axis=2, ov=10, p=(2, 1, 1), thr=40, mo=(6, 6, 4)),
        "dsr": Pair(
            size=(30, 90, 24),
            axis=1,
            ov=12,
            p=(-2, 1, 0),
            thr=45,
            blank=3,
            kind="dsr",
            order="x,y,z",
        ),
    }
    out += [reg3d_case(f"crop.{k}", pr) for k, pr in crops.items()]
    return out


def axes_cases():
    """RegistrationXcorrTest.pairRegistrationTakesTheAxesOfTheDataOrder: tiles gt(1:40, 1:40, :)
    and gt(5:44, 25:66, :) of a 44 x 66 x 16 uint16 volume; tile 2 listed at x 30, y 0 instead of
    24, 4; max offsets 10, 10, 2; x downsampled by 2"""
    out = []
    for order in ("y,x,z", "x,y,z", "z,y,x", "z,x,y"):
        perm = perm_of(order)
        ds = [[1, 2, 1][k - 1] for k in perm]
        dom = dom_of(order)
        c1 = [0.0, 0.0, 0.0, 39.0, 39.0, 15.0]
        c2 = [30.0, 0.0, 0.0, 71.0, 39.0, 15.0]
        co = [30.0, 0.0, 0.0, 39.0, 39.0, 15.0]
        args = [
            InPath("d"),
            ".zarr",
            c1,
            c2,
            co,
            [1, 1, 1],
            dom,
            "downSample",
            ds,
            "MaxOffset",
            [10, 10, 2],
        ]
        out.append(
            Case(
                f"axes.{order.replace(',', '')}",
                code=_REG3D,
                args=args,
                nargout=2,
                python=py_reg3d,
                py_args=[
                    InPath("d"),
                    ".zarr",
                    c1,
                    c2,
                    co,
                    [1, 1, 1],
                    dom,
                    PerAxis(ds, "yxz"),
                    PerAxis([10, 10, 2], "xyz"),
                    10000,
                    None,
                ],
                py_kwargs={"backend": "exact"},
                inputs={"d": Tree(functools.partial(axes_tree, order), "d")},
                files=False,
                note="permissive:shift",
            )
        )
    return out


@functools.cache
def axes_gt() -> np.ndarray:
    return CounterRNG(5, "axes").integers(1000, 60000, (44, 66, 16)).astype(np.uint16)


def axes_tree(order: str, folder) -> None:
    gt = axes_gt()  # [y x z]
    t = (gt[0:40, 0:40, :], gt[4:44, 24:66, :])
    for j in range(2):
        write_tile(Path(folder) / f"t{j + 1}.zarr", np.ascontiguousarray(t[j]), order, ".zarr")


def r2_cases():
    out = [reg2d_case(f"r2.{k:03d}", random_pair(k, two_d=True)) for k in range(N_RANDOM_2D)]
    base = Pair(size=(48, 52, 1), axis=0, ov=14, p=(3, -2, 0), ext=".tif", two_d=True)
    special = {
        "empty_template": dataclasses.replace(base, kind="zero"),
        "holes": dataclasses.replace(base, kind="holes"),
        "holes_xyz": dataclasses.replace(base, kind="holes", order="x,y,z", swap=True),
        "dsr": dataclasses.replace(base, kind="dsr"),
        "sparse": dataclasses.replace(base, kind="sparse", seed=5),
        "pooled": dataclasses.replace(base, ds=(2, 2, 1), p=(4, -2, 0)),
        "single": dataclasses.replace(base, cls="single", vs=(0.108, 0.108, 0.3)),
    }
    out += [reg2d_case(f"r2.{k}", pr) for k, pr in special.items()]
    return out


# --- MIP slabs ------------------------------------------------------------------------------


# the slab of the files' dimension i is named after it, as XR_MIP_zarr names its MIPs, whatever
# the data order (MATLAB took the names by axis until 0fe22a5, right only for 'y,x,z')
SLAB_NAMES = ("y", "x", "z")


def slab_pool(pool_m) -> tuple[list[float], list[float]]:
    """poolSize ([y x z], 3 or 6 values) as the slabs' pool and the secondary pool"""
    p = list(pool_m)
    return (p[:3], p[3:6]) if len(p) >= 6 else (p[:3], [1, 1, 1])


def mip_dir_of(pool_m) -> str:
    return "MIP_Slabs_pooling_" + "_".join(str(int(v)) for v in pool_m)


def slab_of(a_file: np.ndarray, i: int, pool, pool1) -> np.ndarray:
    """the slab of dimension i (0-based) of a file array: max pooled by pool[i] along i, by pool1
    along the others (as large-file registration takes the slabs' voxel sizes)"""
    f = [int(v) for v in pool1]
    f[i] = int(pool[i])
    ny, nx, nz = a_file.shape
    pad = [(0, (-n) % q) for n, q in zip(a_file.shape, f, strict=True)]
    v = np.pad(a_file, pad)
    v = v.reshape(v.shape[0] // f[0], f[0], v.shape[1] // f[1], f[1], v.shape[2] // f[2], f[2])
    return np.ascontiguousarray(v.max(axis=(1, 3, 5)))


def slab_tree(pr: Pair, pool_m):
    def fill(folder):
        folder = Path(folder)
        t1, t2, _, _ = tiles(pr)
        if pr.swap:
            t1, t2 = t2, t1
        mip = folder / mip_dir_of(pool_m)
        mip.mkdir()
        pool, pool1 = slab_pool(pool_m)
        for j, t in ((1, t1), (2, t2)):
            f = np.transpose(t, [k - 1 for k in perm_of(pr.order)])
            for i in range(3):
                synth.write_zarr_matlab(
                    mip / f"t{j}_{SLAB_NAMES[i]}.zarr", slab_of(f, i, pool, pool1), (16, 16, 16)
                )

    return fill


def slab_case(cid: str, pr: Pair, pool_m, **kw) -> Case:
    c1, c2, co = cuboids(pr)
    dom = dom_of(pr.order)
    mip = mip_dir_of(pool_m)
    args = [
        InPath("d"),
        Work("x.mat"),
        c1,
        c2,
        co,
        list(pr.vs),
        dom,
        "downSample",
        list(pr.ds),
        "MaxOffset",
        list(pr.mo),
        "mipDirStr",
        mip,
        "poolSize",
        list(pool_m),
        "parseCluster",
        False,
    ]
    pool_py = [v for g in range(0, len(pool_m), 3) for v in from_yxz(list(pool_m[g : g + 3]))]
    return Case(
        cid,
        code=_SLABS,
        args=args,
        nargout=4,
        python=py_slabs,
        py_args=[
            InPath("d"),
            Work("x.mat"),
            c1,
            c2,
            co,
            list(pr.vs),
            dom,
            PerAxis(list(pr.ds), "yxz"),
            PerAxis(list(pr.mo), "xyz"),
            mip,
            pool_py,
        ],
        py_kwargs={"backend": "exact"},
        inputs={"d": Tree(slab_tree(pr, pool_m), "d")},
        files=False,
        note="permissive:slab",
        **kw,
    )


def slab_cases():
    out = []
    pools = ((2, 2, 2), (2, 2, 2), (3, 3, 2), (2, 2, 2, 1, 1, 1), (2, 2, 2, 2, 2, 1))
    for k in range(24):
        pk = Picker(f"slab/{k}")
        axis = pk.of((0, 1))
        mo = pk.of(((20, 20, 10), (12, 12, 6), (16, 10, 8)))
        p = (pk.int(-4, 4), pk.int(-4, 4), pk.int(-2, 2))
        pr = Pair(
            size=(pk.int(56, 72), pk.int(56, 72), pk.int(32, 40)),
            axis=axis,
            ov=pk.int(14, 20),
            p=p,
            kind=pk.of(("texture", "texture", "dsr", "sparse")),
            order=pk.of(ORDERS[:2]),
            mo=mo,
            swap=pk.int(0, 1) == 1,
            vs=pk.of(_VS),
            seed=100 + k,
        )
        out.append(slab_case(f"slab.{k:02d}", pr, pk.of(pools)))
    # StitchingTest.largeFileStitching's tiles: 64 x 72 x 32, pool 2, max offsets 20, 20, 10
    pr = Pair(size=(64, 72, 32), axis=0, ov=28, p=(-2, 0, 0), mo=(20, 20, 10), seed=150)
    out.append(slab_case("slab.stitching_test", pr, (2, 2, 2)))
    return out


# --- the wrapper's result files ---------------------------------------------------------------


def common_size(prs) -> tuple[int, int, int]:
    """a ground truth ([y x z]) that holds the tiles of all the pairs"""
    ny, nx, nz = prs[0].size
    o = np.max([pr.offset2 for pr in prs], axis=0)
    return (ny + int(o[1]), nx + int(o[0]), nz + int(o[2]))


def wrap_tree(prs, ext, slabs=None):
    def fill(folder):
        folder = Path(folder)
        gs = common_size(prs)
        t1, _, _, _ = tiles(prs[0], gs)
        write_tile(folder / f"t1{ext}", t1, prs[0].order, ext)
        for k, pr in enumerate(prs):
            _, t2, _, _ = tiles(pr, gs)
            write_tile(folder / f"t{k + 2}{ext}", t2, pr.order, ext)
        if slabs is not None:
            mip = folder / mip_dir_of(slabs)
            mip.mkdir()
            pool, pool1 = slab_pool(slabs)
            arrays = [t1] + [tiles(pr, gs)[1] for pr in prs]
            for j, t in enumerate(arrays):
                f = np.transpose(t, [k - 1 for k in perm_of(prs[0].order)])
                for i in range(3):
                    synth.write_zarr_matlab(
                        mip / f"t{j + 1}_{SLAB_NAMES[i]}.zarr",
                        slab_of(f, i, pool, pool1),
                        (16, 16, 16),
                    )

    return fill


def wrap_case(cid, prs, ext=".zarr", stitch_2d=False, slabs=None):
    """the wrapper for tile 1 and tiles 2.. (each pair's tile 1 has the same content: the pairs
    share the ground truth's first tile)"""
    n = len(prs)
    cub = [cuboids(pr) for pr in prs]
    c1 = cub[0][0]
    # (matrices: MATLAB's n x 6 and n x 2 double arrays)
    c2 = np.array([c[1] for c in cub])
    co = np.array([c[2] for c in cub])
    pairs = np.array([[1, k + 2] for k in range(n)], dtype=np.float64)
    pr = prs[0]
    dom = dom_of(pr.order)
    nv = [
        "Stitch2D",
        stitch_2d,
        "downSample",
        list(pr.ds),
        "MaxOffset",
        list(pr.mo),
        "parseCluster",
        False,
    ]
    kw = {
        "stitch_2d": stitch_2d,
        "down_sample": PerAxis(list(pr.ds), "yxz"),
        "max_offset": PerAxis(list(pr.mo), "xyz"),
        "backend": "exact",
    }
    if slabs is not None:
        nv += ["largeFile", True, "mipDirStr", mip_dir_of(slabs), "poolSize", list(slabs)]
        kw.update(
            large_file=True,
            mip_dir_str=mip_dir_of(slabs),
            pool_size=list(from_yxz(list(slabs[:3]))),
        )
    args = [InPath("d"), ext, n, Work("xcorr_tile_1.mat"), pairs, c1, c2, co, list(pr.vs), dom, *nv]
    return Case(
        cid,
        code=_WRAP,
        args=args,
        nargout=0,
        python=py_wrapper,
        py_args=[
            InPath("d"),
            ext,
            n,
            Work("xcorr_tile_1.mat"),
            Index1(pairs),
            c1,
            c2,
            co,
            list(pr.vs),
            dom,
        ],
        py_kwargs=kw,
        inputs={"d": Tree(wrap_tree(prs, ext, slabs), "d")},
        files=True,
        # (2D peaks: normxcorr2's proposed tolerance; shifts are integers, so it keeps them exact)
        tree=tree(
            values=N2_PEAK if stitch_2d else EXACT,
            ignore=[r"(^|/)job_logs?(/|$)", r"\.tmp$", r"(^|/)tmp(/|$)"],
        ),
        note="permissive:files",
    )


def wrap_cases():
    a = Pair(size=(40, 44, 20), axis=0, ov=14, p=(2, -1, 1), seed=200)
    b = dataclasses.replace(a, axis=1, p=(1, 2, -1))
    c = dataclasses.replace(a, axis=2, ov=10, p=(-1, 1, 1))
    a2 = Pair(size=(48, 52, 1), axis=0, ov=14, p=(3, -2, 0), ext=".tif", two_d=True, seed=201)
    b2 = dataclasses.replace(a2, axis=1, p=(-2, 3, 0))
    big = Pair(size=(64, 64, 32), axis=0, ov=20, p=(-2, 1, 0), mo=(20, 20, 10), seed=202)
    big_b = dataclasses.replace(big, axis=1, p=(1, -2, 0))
    return [
        wrap_case("wrap.two", [a, b]),
        wrap_case("wrap.three", [a, b, c]),
        wrap_case("wrap.tif", [a, c], ext=".tif"),
        wrap_case("wrap.2d", [a2, b2], ext=".tif", stitch_2d=True),
        wrap_case("wrap.large_one", [big], slabs=(2, 2, 2)),
        wrap_case("wrap.large_two", [big, big_b], slabs=(2, 2, 2)),
    ]


# --- errors -----------------------------------------------------------------------------------


def err_cases():
    one = Pair(size=(40, 44, 1), axis=0, ov=14, p=(2, 1, 0), ext=".tif", two_d=True)
    out = [
        dataclasses.replace(
            reg3d_case("err.one_slice", one),
            expect_error="MATLAB:badsubscript",
            py_error=RegistrationError,
        ),
    ]
    base = Pair(size=(36, 40, 20), axis=0, ov=14, p=(3, -2, 1), ext=".tiff")
    out.append(
        dataclasses.replace(
            reg3d_case("err.extension", base),
            py_error=RegistrationError,
            expect_error="MATLAB:UndefinedFunction",
        )
    )
    flat2 = Pair(size=(48, 52, 1), axis=0, ov=14, p=(3, -2, 0), ext=".tif", two_d=True, kind="flat")
    out.append(
        dataclasses.replace(
            reg2d_case("err.flat_template_2d", flat2),
            py_error=RegistrationError,
            expect_error="images:normxcorr2:sameElementsInTemplate",
        )
    )
    return out


MAT_HEADER = b"MATLAB 5.0 MAT-file, Platform: GLNXA64, Created on: Thu Jan  1 00:00:00 2026"


def write_slab_results(folder: Path, rows, peaks, names=("y", "x", "z")):
    """the three slab registrations' results (x_<axis>.mat), as their tasks write them (MAT v6
    files whose header has no time or platform: the same inputs everywhere)"""
    from petakit5d.io.matfile import write_mat

    for k, n in enumerate(names):
        path = folder / f"x_{n}.mat"
        write_mat(
            path,
            {
                "relative_shift_mat": np.asarray(rows[k], dtype=np.float64).reshape(1, 3),
                "max_xcorr_mat": float(peaks[k]),
                "pair_indices": np.asarray([[k + 1, k + 2]], dtype=np.float64),
            },
            version="6",
        )
        data = bytearray(path.read_bytes())
        data[:116] = MAT_HEADER.ljust(116, b" ")
        path.write_bytes(bytes(data))


WEIGHTS = {
    # a slab without data (peak 1 at no shift, 0.3 over the median) weighs 0
    "failed": ([[0, 0, 0], [2, 4, 1], [3, 5, 1]], [1.0, 0.6, 0.7]),
    # 0.2 over the median or less: it weighs
    "not_failed": ([[0, 0, 0], [2, 4, 1], [3, 5, 1]], [1.0, 0.85, 0.82]),
    "weighted": ([[2, 3, 1], [3, 4, 2], [5, 4, 1]], [0.9, 0.3, 0.6]),
    # no weight for an axis: NaN
    "nan": ([[2, 3, 1], [3, 4, 2], [5, 4, 1]], [0.0, 0.0, 0.6]),
    "halves": ([[1, 2, 3], [2, 3, 4], [3, 4, 5]], [0.5, 0.5, 0.5]),
}


def weights_tree(rows, peaks, pool_m):
    def fill(folder):
        folder = Path(folder)
        mip = folder / mip_dir_of(pool_m)
        mip.mkdir()
        for j in (1, 2):
            for i in range(3):
                synth.write_zarr_matlab(
                    mip / f"t{j}_{SLAB_NAMES[i]}.zarr",
                    np.ones((10, 12, 8), np.uint16),
                    (16, 16, 16),
                )
        write_slab_results(folder, rows, peaks, SLAB_NAMES)

    return fill


_SLABS_DONE = """
d = in{1};
[rs, mx, rsm, mxm] = cross_correlation_registration_3d_mip_slabs([d 't1.zarr'], [d 't2.zarr'], ...
    [d 'x.mat'], in{2:end});
out = {rs, mx, rsm, mxm};
"""


def py_slabs_done(d, c1, c2, co, vs, dom, mip, pool):
    r = sr.cross_correlation_registration_3d_mip_slabs(
        f"{d}t1.zarr",
        f"{d}t2.zarr",
        f"{d}x.mat",
        c1,
        c2,
        co,
        vs,
        data_order_axes_from_matlab(dom),
        mip_dir_str=mip,
        pool_size=pool,
        parse_cluster=False,
    )
    return list(r)


def weights_cases():
    out = []
    c1 = [0.0, 0.0, 0.0, 20.0, 22.0, 15.0]
    c2 = [14.0, 0.0, 0.0, 34.0, 22.0, 15.0]
    co = [14.0, 0.0, 0.0, 20.0, 22.0, 15.0]
    for name, (rows, peaks) in WEIGHTS.items():
        for pool_m in ((2, 2, 2), (2, 2, 2, 2, 1, 3)):
            pool_py = [
                v for g in range(0, len(pool_m), 3) for v in from_yxz(list(pool_m[g : g + 3]))
            ]
            mip = mip_dir_of(pool_m)
            args = [
                InPath("d"),
                c1,
                c2,
                co,
                [1, 1, 1],
                [2, 1, 3],
                "mipDirStr",
                mip,
                "poolSize",
                list(pool_m),
                "parseCluster",
                False,
            ]
            out.append(
                Case(
                    f"weights.{name}.{len(pool_m)}",
                    code=_SLABS_DONE,
                    args=args,
                    nargout=4,
                    python=py_slabs_done,
                    py_args=[InPath("d"), c1, c2, co, [1, 1, 1], [2, 1, 3], mip, pool_py],
                    inputs={"d": Tree(weights_tree(rows, peaks, pool_m), "d")},
                    files=True,
                    tree=tree(ignore=[r"(^|/)job_logs?(/|$)", r"\.tmp$"]),
                )
            )
    return out


CASES_3D = r3_cases() + axes_cases()
CASES_2D = r2_cases()
CASES_SLAB = slab_cases()
CASES_WRAP = wrap_cases()
CASES_ERR = err_cases()
CASES_WEIGHTS = weights_cases()
CASES = CASES_3D + CASES_2D + CASES_SLAB + CASES_WRAP + CASES_ERR + CASES_WEIGHTS


# ---------------------------------------------------------------------------------------------
# running them on each backend
# ---------------------------------------------------------------------------------------------


def variant(case: Case, **changes) -> Case:
    v = dataclasses.replace(case, **changes)
    v._module_globals = case._module_globals
    return v


def uses_fft(case: Case) -> bool:
    return case.note.startswith("permissive")


def permissive(case: Case) -> Case:
    """the case on the permissive backend: the shifts exact, the peaks within record 0006"""
    kind = case.note.split(":")[-1]
    kw = {**(case.py_kwargs or {}), "backend": "permissive"}
    if kind == "shift":
        compare = {1: NCC_TOL}
        return variant(case, py_kwargs=kw, compare=compare)
    if kind == "slab":
        return variant(case, py_kwargs=kw, compare={1: NCC_TOL, 3: NCC_TOL})
    # files: shifts are integers, so 5e-3 keeps them exact
    return variant(case, py_kwargs=kw, tree=dataclasses.replace(case.tree, values=NCC_TOL))


@oracle.parametrize(CASES)
def test_matlab_exact(case, tmp_path):
    """bit for bit with the exact FFT backend (2D and error cases: as declared)"""
    if uses_fft(case) and digests_hold_here():
        pytest.skip(digests_hold_here())
    oracle.check(case, tmp_path)


@oracle.parametrize([c for c in CASES if uses_fft(c)])
def test_matlab_permissive(case, tmp_path):
    """the same shifts and the peaks within decision record 0006 on the permissive backend"""
    oracle.check(permissive(case), tmp_path)


# ---------------------------------------------------------------------------------------------
# the intent of the MATLAB tests (both backends)
# ---------------------------------------------------------------------------------------------


@pytest.fixture(params=["permissive", "exact"])
def be(request):
    if request.param not in pfft.available():
        pytest.skip("the exact FFT backend (the FFTW add-on) is not installed")
    return request.param


def run_pair(tmp_path, pr: Pair, backend=None, out=""):
    folder = tmp_path / "d"
    folder.mkdir(parents=True)
    pair_tree(pr)(folder)
    c1, c2, co = cuboids(pr)
    fn = sr.cross_correlation_registration_2d if pr.two_d else sr.cross_correlation_registration_3d
    kw = (
        {}
        if pr.two_d
        else {"backend": backend, "dim_num_thrsh": pr.thr, "blank_num_trsh": pr.blank}
    )
    return fn(
        str(folder / f"t1{pr.ext}"),
        str(folder / f"t2{pr.ext}"),
        out,
        c1,
        c2,
        co,
        pr.vs,
        data_order_axes_from_matlab(dom_of(pr.order)),
        down_sample=from_yxz(list(pr.ds)),
        max_offset=from_xyz(list(pr.mo)),
        **kw,
    )


def test_pair_registration_takes_the_axes_of_the_data_order(tmp_path, be):
    """the same pair in each data order gives the shift of 'y,x,z' data: tile 2 is 6 voxels
    before its stage position in x and 4 in y"""
    shifts = []
    for order in ("y,x,z", "x,y,z", "z,y,x", "z,x,y"):
        d = tmp_path / order.replace(",", "")
        d.mkdir()
        axes_tree(order, d)
        perm = perm_of(order)
        ds = [[1, 2, 1][k - 1] for k in perm]
        rel, mx = sr.cross_correlation_registration_3d(
            str(d / "t1.zarr"),
            str(d / "t2.zarr"),
            "",
            [0, 0, 0, 39, 39, 15],
            [30, 0, 0, 71, 39, 15],
            [30, 0, 0, 39, 39, 15],
            [1, 1, 1],
            data_order_axes_from_matlab(dom_of(order)),
            down_sample=from_yxz(ds),
            max_offset=from_xyz([10, 10, 2]),
            backend=be,
        )
        shifts.append(rel)
        assert mx > 0.9
    np.testing.assert_array_equal(np.abs(shifts[0]), [6, 4, 0])
    for s in shifts[1:]:
        np.testing.assert_array_equal(s, shifts[0])


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("ext", [".zarr", ".tif"])
def test_registration_recovers_the_listed_offset(tmp_path, be, order, ext):
    """StitchingTest.xcorrRecoversPerturbedTilePosition at the pair level: tile 2 listed off its
    position (x, y, z) is found back, in every data order and file type"""
    pr = Pair(size=(40, 44, 24), axis=0, ov=16, p=(4, -2, 1), order=order, ext=ext, seed=7)
    rel, mx = run_pair(tmp_path, pr, be)
    np.testing.assert_array_equal(np.abs(rel), np.abs(pr.p))
    assert mx > 0.99


def test_registration_signs_follow_the_tiles_order(tmp_path, be):
    """either tile may lead: the shift's size is the same"""
    pr = Pair(size=(40, 44, 24), axis=0, ov=16, p=(4, -2, 1), seed=8)
    a, _ = run_pair(tmp_path / "a", pr, be)
    b, _ = run_pair(tmp_path / "b", dataclasses.replace(pr, swap=True), be)
    np.testing.assert_array_equal(np.abs(a), np.abs(pr.p))
    np.testing.assert_array_equal(np.abs(b), np.abs(pr.p))


def test_downsampled_registration_moves_in_steps_of_the_pooling(tmp_path, be):
    pr = Pair(size=(40, 44, 24), axis=0, ov=16, p=(4, -2, 2), ds=(2, 2, 2), seed=9)
    rel, _ = run_pair(tmp_path, pr, be)
    np.testing.assert_array_equal(np.abs(rel), [4, 2, 2])


def test_a_template_without_data_is_not_registered(tmp_path, be):
    pr = Pair(size=(36, 40, 20), axis=0, ov=14, p=(3, -2, 1), kind="zero")
    rel, mx = run_pair(tmp_path, pr, be)
    np.testing.assert_array_equal(rel, [0, 0, 0])
    assert mx == 1 and isinstance(mx, float)


def test_registration_2d_recovers_the_listed_offset(tmp_path):
    """StitchingTest.registration2DTiles at the pair level"""
    for k, order in enumerate(("y,x,z", "x,y,z")):
        pr = Pair(
            size=(48, 52, 1),
            axis=k,
            ov=14,
            p=(3, -2, 0),
            ext=".tif",
            two_d=True,
            order=order,
            seed=10,
        )
        d = tmp_path / str(k)
        d.mkdir()
        rel, mx = run_pair(d, pr)
        np.testing.assert_array_equal(np.abs(rel), [3, 2, 0])
        # (listed 2 voxels before its position along y, tile 2's template overhangs tile 1)
        assert mx > (0.99 if k == 0 else 0.7)


def test_native_zarr_tiles_give_the_same_shift(tmp_path, be):
    """the port's own (z, y, x) zarr tiles register as MATLAB's layout does"""
    pr = Pair(size=(40, 44, 24), axis=1, ov=16, p=(1, 3, -1), seed=11, kind="dsr")
    ref = run_pair(tmp_path, pr, be)
    t1, t2, _, _ = tiles(pr)
    d = tmp_path / "native"
    d.mkdir()
    for j, t in ((1, t1), (2, t2)):
        synth.write_zarr(d / f"t{j}.zarr", from_matlab(t), (16, 16, 16))
    c1, c2, co = cuboids(pr)
    rel, mx = sr.cross_correlation_registration_3d(
        str(d / "t1.zarr"),
        str(d / "t2.zarr"),
        "",
        c1,
        c2,
        co,
        pr.vs,
        "z,y,x",
        backend=be,
        max_offset=from_xyz(list(pr.mo)),
    )
    np.testing.assert_array_equal(rel, ref[0])
    assert mx == ref[1]


def test_result_file_and_an_existing_one(tmp_path, be):
    pr = Pair(size=(36, 40, 20), axis=0, ov=14, p=(3, -2, 1))
    out = str(tmp_path / "r.mat")
    rel, mx = run_pair(tmp_path, pr, be, out=out)
    a = read_mat(out)
    np.testing.assert_array_equal(np.asarray(a["relative_shift"]).ravel(), rel)
    assert a["max_xcorr"] == mx and a["max_xcorr"].dtype == np.float32
    assert sorted(os.listdir(tmp_path)) == ["d", "r.mat"]
    # an existing result is kept
    folder = tmp_path / "d"
    c1, c2, co = cuboids(pr)
    assert (
        sr.cross_correlation_registration_3d(
            str(folder / "t1.zarr"), str(folder / "t2.zarr"), out, c1, c2, co, pr.vs, "z,y,x"
        )
        is None
    )


def test_wrapper_writes_matlabs_result_file(tmp_path, be):
    a = Pair(size=(40, 44, 20), axis=0, ov=14, p=(2, -1, 1), seed=200)
    b = dataclasses.replace(a, axis=1, p=(1, 2, -1))
    d = tmp_path / "d"
    d.mkdir()
    wrap_tree([a, b], ".zarr")(d)
    cub = [cuboids(pr) for pr in (a, b)]
    out = str(tmp_path / "xcorr_tile_1.mat")
    with pfft.using(be):
        sr.cross_correlation_registration(
            str(d / "t1.zarr"),
            [str(d / "t2.zarr"), str(d / "t3.zarr")],
            out,
            [[0, 1], [0, 2]],
            cub[0][0],
            [c[1] for c in cub],
            [c[2] for c in cub],
            a.vs,
            (2, 1, 0),
            max_offset=from_xyz(list(a.mo)),
            parse_cluster=False,
        )
    m = read_mat(out)
    np.testing.assert_array_equal(m["pair_indices"], [[1, 2], [1, 3]])
    np.testing.assert_array_equal(np.abs(m["relative_shift_mat"]), [[2, 1, 1], [1, 2, 1]])
    assert m["max_xcorr_mat"].shape == (2, 1)


def test_wrapper_is_the_registered_task_of_the_shift_assignment():
    from petakit5d.params import find
    from petakit5d.stitch import assignment

    assert assignment.REGISTRATION_TASK == (
        "petakit5d.stitch.registration:cross_correlation_registration"
    )
    spec = find("cross_correlation_registration_wrapper")
    assert spec.python_name == "petakit5d.stitch.registration.cross_correlation_registration"


def test_mip_slab_registration_recovers_the_offset(tmp_path, be):
    """StitchingTest.largeFileRegistrationRecoversPerturbedTilePosition at the pair level: the
    weighted shifts of the three slabs give the listed offset"""
    pr = Pair(size=(64, 72, 32), axis=0, ov=28, p=(-2, 0, 0), mo=(20, 20, 10), seed=150)
    d = tmp_path / "d"
    d.mkdir()
    slab_tree(pr, (2, 2, 2))(d)
    c1, c2, co = cuboids(pr)
    out = str(tmp_path / "x.mat")
    with pfft.using(be):
        rel, mx, rsm, mxm = sr.cross_correlation_registration_3d_mip_slabs(
            str(d / "t1.zarr"),
            str(d / "t2.zarr"),
            out,
            c1,
            c2,
            co,
            pr.vs,
            "z,y,x",
            max_offset=from_xyz(list(pr.mo)),
            mip_dir_str=mip_dir_of((2, 2, 2)),
            pool_size=(2, 2, 2),
            parse_cluster=False,
        )
    np.testing.assert_array_equal(np.abs(rel), [2, 0, 0])
    assert rsm.shape == (3, 3) and mxm.shape == (3, 1)
    assert {p.name for p in tmp_path.iterdir()} >= {"x.mat", "x_y.mat", "x_x.mat", "x_z.mat"}
    # an existing result is read back
    again = sr.cross_correlation_registration_3d_mip_slabs(
        "a/t1.zarr", "a/t2.zarr", out, c1, c2, co, pr.vs, "z,y,x"
    )
    np.testing.assert_array_equal(np.asarray(again[0]).ravel(), rel)


def test_mip_slab_registration_takes_the_axes_of_the_data_order(tmp_path, be):
    """(RegistrationXcorrTest/mipSlabRegistrationTakesTheAxesOfTheDataOrder) the same pair in
    each data order gives the shift of 'y,x,z' data, 9 voxels in x and 8 in y within max offsets
    of 10, 10 and 2: the slabs are named after the files' dimensions and take each axis's pooling
    (MATLAB named them by axis and halved y's bound by x's pooling and the reverse)"""
    gt = CounterRNG(0, "sw8/slab_axes").integers(1000, 60000, (32, 48, 80), dtype=np.uint16)
    gt = to_matlab(gt)  # [y x z]
    t = (gt[0:40, 0:44, :], gt[8:48, 30:80, :])
    shifts = []
    for order in ("y,x,z", "x,y,z", "z,y,x", "z,x,y"):
        d = tmp_path / order.replace(",", "")
        d.mkdir()
        fns = [str(d / f"t{j + 1}.zarr") for j in range(2)]
        for j in range(2):
            write_tile(Path(fns[j]), np.ascontiguousarray(t[j]), order, ".zarr")
        mip = sr.compute_tile_mip_slabs(fns, (2, 2, 2), parse_cluster=False)
        with pfft.using(be):
            rel, mx, _, _ = sr.cross_correlation_registration_3d_mip_slabs(
                fns[0],
                fns[1],
                str(d / "xcorr.mat"),
                [0, 0, 0, 43, 39, 31],
                [21, 0, 0, 70, 39, 31],
                [21, 0, 0, 43, 39, 31],
                [1, 1, 1],
                data_order_axes_from_matlab(dom_of(order)),
                max_offset=from_xyz([10, 10, 2]),
                mip_dir_str=mip,
                pool_size=(2, 2, 2),
                major_downsample=1,
                parse_cluster=False,
            )
        shifts.append(np.asarray(rel, dtype=np.float64).ravel())
    np.testing.assert_array_equal(np.abs(shifts[0]), [9, 8, 0])
    for s in shifts[1:]:
        np.testing.assert_array_equal(s, shifts[0])


def test_failed_slab_weighs_nothing(tmp_path):
    """the slabs' results: a failed slab (peak 1 at no shift while 0.3 over the median) weighs
    0; each axis is the average of its two slabs' shifts, weighted by their squared peaks"""
    rows, peaks = WEIGHTS["failed"]
    weights_tree(rows, peaks, (2, 2, 2))(tmp_path)
    rel, mx, rsm, mxm = py_slabs_done(
        f"{tmp_path}/",
        [0, 0, 0, 20, 22, 15],
        [14, 0, 0, 34, 22, 15],
        [14, 0, 0, 20, 22, 15],
        [1, 1, 1],
        [2, 1, 3],
        mip_dir_of((2, 2, 2)),
        (2, 2, 2),
    )
    # x from the y and z slabs, y from the x (weight 0) and z slabs, z from y and x
    np.testing.assert_array_equal(mxm.ravel(), [0, 0.6, 0.7])
    np.testing.assert_array_equal(rel, [3, 5, 1])
    assert mx == pytest.approx((0.6**3 + 0.7**3) / (0.6**2 + 0.7**2))


# ---------------------------------------------------------------------------------------------
# StitchingTest's registration methods, with the port's registration in the stitching
# ---------------------------------------------------------------------------------------------


PX = 0.108


@pytest.fixture
def s05_axes(monkeypatch):
    """stand-in for the request to S05: the registration tasks get the data order's axes (the
    Python form of data_order_mat) instead of its string; a no-op once
    assignment._registration_task passes them"""
    from petakit5d.jobs.tasks import Task
    from petakit5d.stitch import assignment as sa
    from petakit5d.stitch.geometry import parse_data_order

    orig = sa._registration_task

    def task(*args):
        t = orig(*args)
        a = list(t.args)
        if isinstance(a[8], str):
            a[8] = list(parse_data_order(a[8]))
        return Task(t.func, a, t.kwargs)

    monkeypatch.setattr(sa, "_registration_task", task)


def random_volume(shape, seed: int = 1):
    """StitchingTest.randomVolume: nonzero random uint16 values"""
    return CounterRNG(seed, "s04/gt").integers(1000, 60000, shape).astype(np.uint16)


def stitch_tiles(data: Path, gt, starts, tile, stage=None, tidx=None, pixel=(PX, PX, PX), **kw):
    """the tiles of gt at starts ((z0, y0, x0)) and the stitching wrapper with StitchingTest's
    defaults (DSR tiles, isotropic, 'x,y,z') and the registration"""
    from petakit5d.stitch.dataset import stitch_dataset

    tiles = synth.cut_tiles(gt, starts, tile, pixel, stage, tidx)
    csv = synth.write_stitch_tiles(str(data), tiles)[1]
    args = {
        "dsr": True,
        "resample_type": "isotropic",
        "xy_pixel_size": PX,
        "axis_order": "x,y,z",
        "channel_patterns": ("CamA_ch0",),
        "primary_ch": "CamA_ch0",
        "parse_cluster": False,
        "xcorr_shift": True,
        "xcorr_downsample": (1, 1, 1),
        "xy_max_offset": 20,
        "z_max_offset": 10,
        **kw,
    }
    stitch_dataset(str(data) + "/", csv, **args)
    name = synth.stitch_tile_name((0, 0, 0)).rsplit("_000x", 1)[0]
    return name


def stitched(data: Path, name: str) -> np.ndarray:
    from petakit5d.io import read_zarr

    return read_zarr(f"{data}/matlab_stitch/{name}.zarr")


def xcorr_shifts(data: Path, name: str) -> np.ndarray:
    a = read_mat(f"{data}/matlab_stitch/xcorr/{name}.mat", ["dxyz_shift"])
    from petakit5d.core.dtypes import mround

    return mround(np.asarray(a["dxyz_shift"], dtype=np.float64) / PX)


def crop_to_data(v: np.ndarray) -> np.ndarray:
    nz = np.nonzero(v)
    return v[tuple(slice(int(a.min()), int(a.max()) + 1) for a in nz)]


@pytest.mark.parametrize("shift_method", ["grid", "local"])
def test_stitching_recovers_the_perturbed_tile(tmp_path, s05_axes, shift_method):
    """StitchingTest.xcorrRecoversPerturbedTilePosition and localShiftMethod: the second tile is
    listed 4 voxels off in x; the registration finds it, and the stitched volume is the ground
    truth (padded by 2 voxels per tile on each side)"""
    gt = random_volume((32, 64, 104))
    name = stitch_tiles(
        tmp_path,
        gt,
        [(0, 0, 0), (0, 0, 40)],
        (32, 64, 64),
        stage=[[0, 0, 0], [44 * PX, 0, 0]],
        tidx=[[0, 0, 0], [1, 0, 0]],
        blend_method="none",
        shift_method=shift_method,
    )
    np.testing.assert_array_equal(xcorr_shifts(tmp_path, name), [[0, 0, 0], [-4, 0, 0]])
    out = stitched(tmp_path, name)
    assert out.shape == (32 + 2 * 2 + 2, 64 + 2 * 2 + 2, 104 + 2 * 4 + 2)
    np.testing.assert_array_equal(out[2:34, 2:66, 4:108], gt)


def test_stitching_two_value_max_offset_is_y_and_x(tmp_path, s05_axes):
    """StitchingTest.twoValueMaxOffsetIsXAndY: the second tile, listed 6 voxels off in x, is
    found within max offsets of 10 in x and 2 in y (MATLAB's [10, 2]; two values garbled the
    stitching string until 0fe22a5)"""
    gt = random_volume((32, 64, 104))
    name = stitch_tiles(
        tmp_path,
        gt,
        [(0, 0, 0), (0, 0, 40)],
        (32, 64, 64),
        stage=[[0, 0, 0], [46 * PX, 0, 0]],
        tidx=[[0, 0, 0], [1, 0, 0]],
        blend_method="none",
        xy_max_offset=(2, 10),
        z_max_offset=2,
    )
    np.testing.assert_array_equal(xcorr_shifts(tmp_path, name), [[0, 0, 0], [-6, 0, 0]])


def test_stitching_registers_2d_tiles(tmp_path, s05_axes):
    """StitchingTest.registration2DTiles: 2D registration with xcorrMode 'primary'"""
    gt = random_volume((64, 104))
    name = stitch_tiles(
        tmp_path,
        gt,
        [(0, 0), (0, 40)],
        (64, 64),
        stage=[[0, 0, 0], [44 * PX, 0, 0]],
        tidx=[[0, 0, 0], [1, 0, 0]],
        dsr=False,
        io_scan=True,
        xcorr_mode="primary",
        blend_method="feather",
        block_size=(1, 64, 64),
        batch_size=(1, 256, 256),
    )
    np.testing.assert_array_equal(crop_to_data(stitched(tmp_path, name)).reshape(gt.shape), gt)


@pytest.mark.slow
def test_stitching_registers_in_skewed_space(tmp_path, s05_axes):
    """StitchingTest.skewedSpaceRegistration: the y = 1 tiles are listed 4 voxels off in y
    ('-x,y,z', feather)"""
    n, dz = 10, 0.3
    gt = random_volume((50, 100, 64))
    starts = [(0, 0, 0), (0, 36, 0), (n, 0, 0), (n, 36, 0)]
    stage = [[0, 0, 0], [0, 40 * PX, 0], [n * dz, 0, 0], [n * dz, 40 * PX, 0]]
    tidx = [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]]
    name = stitch_tiles(
        tmp_path,
        gt,
        starts,
        (40, 64, 64),
        stage=stage,
        tidx=tidx,
        dsr=False,
        dz=dz,
        skew_angle=32.45,
        reverse=True,
        axis_order="-x,y,z",
        blend_method="feather",
    )
    np.testing.assert_array_equal(crop_to_data(stitched(tmp_path, name)), gt)


@pytest.mark.slow
def test_stitching_large_files(tmp_path, s05_axes):
    """StitchingTest.largeFileStitching: MIP slabs made and registered (2 x 2 tiles)"""
    gt = random_volume((32, 100, 120))
    name = stitch_tiles(
        tmp_path,
        gt,
        [(0, 0, 0), (0, 36, 0), (0, 0, 48), (0, 36, 48)],
        (32, 64, 72),
        blend_method="feather",
        large_file=True,
        pool_size=(2, 2, 2),
        batch_size=(16, 32, 32),
        block_size=(16, 16, 16),
    )
    out = stitched(tmp_path, name)
    assert out.shape == (32 + 2 * 2 + 2, 100 + 2 * 4 + 2, 120 + 2 * 4 + 2)
    np.testing.assert_array_equal(out[2:34, 4:104, 4:124], gt)


@pytest.mark.slow
def test_stitching_large_files_recovers_the_perturbed_tile(tmp_path, s05_axes):
    """StitchingTest.largeFileRegistrationRecoversPerturbedTilePosition: the third tile is listed
    2 voxels off in x; each pair registers its own tiles"""
    gt = random_volume((32, 100, 120))
    stage = [[0, 0, 0], [0, 36 * PX, 0], [50 * PX, 0, 0], [48 * PX, 36 * PX, 0]]
    tidx = [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]]
    name = stitch_tiles(
        tmp_path,
        gt,
        [(0, 0, 0), (0, 36, 0), (0, 0, 48), (0, 36, 48)],
        (32, 64, 72),
        stage=stage,
        tidx=tidx,
        blend_method="feather",
        large_file=True,
        pool_size=(2, 2, 2),
        batch_size=(16, 32, 32),
        block_size=(16, 16, 16),
    )
    rows = sorted(map(tuple, xcorr_shifts(tmp_path, name).tolist()))
    assert rows == [(-2, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)]
    np.testing.assert_array_equal(stitched(tmp_path, name)[2:34, 4:104, 4:124], gt)
