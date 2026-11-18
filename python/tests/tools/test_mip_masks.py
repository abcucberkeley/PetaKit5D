"""MIP masks and the figure tools (petakit5d.tools.mip_masks, petakit5d.tools.figures, task B10).

No MATLAB tests cover XR_generate_MIP_masks, gererate_single_MIP_mask,
XR_visualize_tile_positions or visualize_OTF_and_mask_outline; DeconFrameTest's
otfMaskVisualizationSavesTheFigure covers XR_visualize_OTF_mask_segmentation (ported below).

Oracle cases against MATLAB:

- gererate_single_MIP_mask: tiff and zarr MIPs, tiff and zarr masks, uint16 and single masks,
  equal and unequal dilate sizes (one of them 0), the defaults, an empty mask (box ``[]``), an
  existing result, single input with a threshold that is not a single value, and an all-zero
  MIP (MATLAB's ``clim`` fails in the figure). The masks are compared byte for byte (tiff) or
  by data and metadata (zarr), the box files by content; the figures exist on both sides.
- XR_generate_MIP_masks: tiff and zarr datasets of MIPs along y, x and z of two channels and one
  that matches no pattern, per-axis thresholds and dilate sizes, scalars, zarr masks, an axis
  left out (MATLAB's misaligned axes), missing MIPs of an axis (warnings).
- the figures: visualize_OTF_and_mask_outline on volumes (even and odd sizes, a double mask) and
  an image, XR_visualize_OTF_mask_segmentation on a PSF file (both FFT backends) and
  XR_visualize_tile_positions on image lists with and without camera: the data the figures
  show (images, points, limits, color limits, labels, titles, texts and their positions, the
  view), read from MATLAB's figure handles and from matplotlib's figures, and that the files
  are written. Pixels are not compared. The figure cases need matplotlib (skipped without it).
"""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core.arrays import to_matlab
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, File, In, InPath, Tree, Value, Vol, Work
from petakit5d.testing.compare import Custom, MArray, MCell, tol, tree
from petakit5d.testing.fft_reference import digests_hold_here
from petakit5d.testing.synth import CounterRNG
from petakit5d.tools import figures, mip_masks

HAS_MPL = figures.available()
needs_mpl = pytest.mark.skipif(not HAS_MPL, reason="matplotlib (the figures extra) is missing")

# ---------------------------------------------------------------------------------------------
# data (exact arithmetic: bit-stable everywhere)
# ---------------------------------------------------------------------------------------------


def mip_image(shape, name, n_blobs=4, cls="uint16", spikes=6):
    """a (y, x) MIP: noise in [40, 70), blobs ``900 / (1 + d^4)`` (d: distance / radius) and a
    few single bright pixels (removed by the opening)"""
    rng = CounterRNG(0, f"b10/{name}")
    ny, nx = shape
    y = np.arange(ny, dtype=np.float64)[:, None]
    x = np.arange(nx, dtype=np.float64)[None, :]
    im = 40.0 + 30.0 * rng.uniform(shape)
    u = rng.uniform((n_blobs, 3))
    for k in range(n_blobs):
        cy, cx = 8 + (ny - 16) * u[k, 0], 8 + (nx - 16) * u[k, 1]
        r = 4 + 10 * u[k, 2]
        d2 = ((y - cy) ** 2 + (x - cx) ** 2) / (r * r)
        im = im + 900.0 / (1.0 + d2 * d2)
    p = rng.integers(0, ny * nx - 1, spikes)
    im.reshape(-1)[p] = 2000.0
    if cls == "single":
        return im.astype(np.float32)
    return np.floor(im).astype(np.uint16)


def tif_file(shape, name, **kw):
    def make(path):
        synth.write_tiff(path, mip_image(shape, name, **kw))

    return make


def zarr_tree(shape, name, fname="im.zarr", **kw):
    def make(folder):
        a = mip_image(shape, name, **kw)
        synth.write_zarr(Path(folder) / fname, a[None], block_size=(1, 32, 32), layout="matlab")

    return make


# MIP names: along z (ny x nx), y (nx x nz) and x (ny x nz) of a (nz, ny, nx) = (30, 70, 90) stack
NZ, NY, NX = 30, 70, 90
MIP_SHAPES = {"x": (NY, NZ), "y": (NX, NZ), "z": (NY, NX)}
STEMS = ("scan_CamA_ch0_t0000", "scan_CamB_ch0_t0000", "scan_CamA_ch2_t0000")


def mip_dataset(fmt="tif", axes=("x", "y", "z")):
    """a folder of MIPs <stem>_MIP_<axis>.tif (.zarr); the third stem matches no default channel
    pattern"""

    def make(folder):
        for s, stem in enumerate(STEMS):
            for a in axes:
                im = mip_image(MIP_SHAPES[a], f"ds/{fmt}/{s}/{a}")
                p = Path(folder) / f"{stem}_MIP_{a}.{fmt}"
                if fmt == "zarr":
                    synth.write_zarr(p, im[None], block_size=(1, 32, 32), layout="matlab")
                else:
                    synth.write_tiff(p, im)

    return make


# ---------------------------------------------------------------------------------------------
# oracle cases: gererate_single_MIP_mask and XR_generate_MIP_masks
# ---------------------------------------------------------------------------------------------

M = "petakit5d.matlab"
PLOT = r".*_plot\.png$"
FILES = tree(ignore=[PLOT])


def call_matlab(name, *args):
    from petakit5d import matlab

    getattr(matlab, name)(*args)
    return []


def py_single(fn, fnout, *args):
    """gererate_single_MIP_mask's MATLAB arguments: dilateSize in the image's (rows, columns),
    blockSize [y x z]"""
    from petakit5d.core.axes import from_yxz

    a = list(args)
    if len(a) >= 4:
        a[3] = tuple(from_yxz(tuple(a[3])))
    mip_masks.generate_single_mip_mask(fn, fnout, *a)
    return []


def py_single_mpl(fn, fnout, *args):
    pytest.importorskip("matplotlib")
    return py_single(fn, fnout, *args)


def py_single_plot(fn, fnout, *args):
    pytest.importorskip("matplotlib")
    py_single(fn, fnout, *args)
    return [os.path.isfile(f"{os.path.splitext(fnout)[0]}_plot.png")]


def single_case(cid, make, *args, out="m.tif", python=py_single, **kw):
    return Case(
        cid,
        matlab="gererate_single_MIP_mask",
        nargout=0,
        args=[InPath("im"), Work(out), *args],
        python=python,
        inputs={"im": File(make, "im.tif")},
        tree=FILES,
        **kw,
    )


def _py_zarr_in(d, fnout, *args):
    return py_single(d + "im.zarr", fnout, *args)


def _single_cases():
    s = (75, 110)
    img = tif_file(s, "s1")
    u = "u0b10"
    cases = [
        single_case("single.tif.equal", img, 100, 120, [4, 4], [256, 256, 1], True, u),
        single_case("single.tif.unequal", img, 100, 120, [3, 8], [256, 256, 1], True, u),
        single_case("single.tif.rows_only", img, 100, 120, [6, 0], [256, 256, 1], True, u),
        single_case("single.tif.cols_only", img, 100, 120, [0, 5], [256, 256, 1], True, u),
        single_case("single.tif.single_out", img, 100, 120, [2, 2], [256, 256, 1], False, u),
        single_case(
            "single.tif.small_objects_removed",
            tif_file(s, "s2", n_blobs=6),
            150,
            400,
            [1, 2],
            [256, 256, 1],
            True,
            u,
        ),
        single_case("single.zarr_out", img, 100, 120, [3, 1], [32, 48, 1], True, u, out="m.zarr"),
        single_case(
            "single.zarr_out_single", img, 100, 120, [3, 3], [32, 48, 1], False, u, out="m.zarr"
        ),
        Case(
            "single.zarr_in",
            code="gererate_single_MIP_mask([in{1}, 'im.zarr'], in{2 : end}); out = {};",
            nargout=0,
            args=[InPath("d"), Work("m.tif"), 100, 120, [2, 5], [256, 256, 1], True, u],
            python=_py_zarr_in,
            inputs={"d": Tree(zarr_tree((64, 81), "s3"), "d")},
            tree=FILES,
        ),
        single_case(
            "single.single_input_thresh",
            tif_file(s, "s4", cls="single"),
            100.3,
            120,
            [2, 3],
            [256, 256, 1],
            True,
            u,
        ),
        # the defaults: volThresh 1e5 removes every object of a small MIP (box [])
        single_case("single.defaults", img, python=py_single_mpl),
        single_case(
            "single.empty_mask",
            tif_file(s, "s5", n_blobs=0, spikes=0),
            100,
            10,
            [1, 1],
            [256, 256, 1],
            True,
            u,
        ),
        Case(
            "single.figure_written",
            code=(
                "gererate_single_MIP_mask(in{:}); "
                "out = {exist([in{2}(1 : end - 4), '_plot.png'], 'file') == 2};"
            ),
            args=[InPath("im"), Work("m.tif"), 100, 120, [4, 4], [256, 256, 1], True, u],
            python=py_single_plot,
            inputs={"im": File(img, "im.tif")},
            tree=FILES,
        ),
        # an all-zero MIP: prctile(im(:), 99) is 0 and clim([0, 0]) fails
        single_case(
            "single.all_zero_figure_fails",
            lambda p: synth.write_tiff(p, np.zeros(s, np.uint16)),
            100,
            10,
            [1, 1],
            [256, 256, 1],
            True,
            u,
            python=py_single_mpl,
            expect_error=True,
            py_error=(figures.FigureError,),
        ),
    ]
    # an existing result: nothing is written
    cases.append(
        Case(
            "single.exists_skips",
            matlab="gererate_single_MIP_mask",
            nargout=0,
            args=[InPath("im"), InPath("old"), 100, 120, [4, 4], [256, 256, 1], True, u],
            python=py_single,
            inputs={"im": File(img, "im.tif"), "old": File(tif_file((5, 6), "old"), "old.tif")},
            tree=FILES,
        )
    )
    return cases


def wrapper_case(cid, make, *args, expect_error=False):
    return Case(
        cid,
        expect_error=expect_error,
        py_error=(mip_masks.MipMaskError,),
        code="XR_generate_MIP_masks(in{:}); out = {};",
        nargout=0,
        args=[InPath("d"), *args, "parseCluster", False, "masterCompute", True],
        python=lambda *a: call_matlab("XR_generate_MIP_masks", *a),
        inputs={"d": Tree(make, "d")},
        tree=FILES,
    )


def _wrapper_cases():
    per_axis = ["intThresh", [100, 110, 120], "volThresh", 150, "dilateSize", [2, 3, 4]]
    return [
        wrapper_case("ds.tif.per_axis", mip_dataset(), *per_axis),
        wrapper_case(
            "ds.tif.scalars", mip_dataset(), "intThresh", 105, "volThresh", 100, "dilateSize", 3
        ),
        # x left out: the axes of the files are not filtered with them (MATLAB's quirk)
        wrapper_case("ds.tif.axis_left_out", mip_dataset(), *per_axis, "axis", [1, 0, 1]),
        wrapper_case(
            "ds.tif.zarr_masks",
            mip_dataset(),
            *per_axis,
            "saveZarr",
            True,
            "blockSize",
            [32, 32, 1],
            "save16bit",
            False,
        ),
        wrapper_case("ds.zarr", mip_dataset("zarr"), *per_axis, "zarrFile", True),
        wrapper_case("ds.tif.missing_axes_warn", mip_dataset(axes=("z",)), *per_axis),
        # one file in all: MATLAB's fileparts of a 1x1 cell gives a cell, which it wraps in
        # another cell, and sprintf fails
        wrapper_case(
            "ds.tif.one_file_fails",
            mip_dataset(axes=("z",)),
            *per_axis,
            "channelPatterns",
            ["CamA_ch0"],
            "ResultDirStr",
            "m2",
            expect_error=True,
        ),
    ]


def _py_greater(x, t):
    return mip_masks._greater(x, t)


PROBES = [
    # single data against a double threshold
    Case(
        "probe.single_gt_double",
        code="out = {in{1} > in{2}};",
        args=[
            np.array(
                [np.float32(100.3), np.nextafter(np.float32(100.3), np.float32(200)), 100, 101],
                dtype=np.float32,
            ),
            100.3,
        ],
        python=_py_greater,
    ),
]

# ---------------------------------------------------------------------------------------------
# oracle cases: the figures (their data from MATLAB's handles and matplotlib's figures)
# ---------------------------------------------------------------------------------------------

# MATLAB's single log10 is not correctly rounded (it is 1-2 ulp off the rounded double log10 at
# 5-10% of the values, and no single-precision formula reproduces it); the port's is the double
# log10 rounded to single (proposed exception, figure data only: B10's report)
LOG10_ULP = tol(
    reason="proposed (B10): figure images, log10 of single within 2 ulp of MATLAB's", ulp=2
)

# MATLAB: the data of the panels (axes with an image, in creation order), the colorbar and
# the texts of the figure (the title; labels are hidden handles)
FIG_DATA = textwrap.dedent(
    """
    ax = flipud(findobj(fig, 'type', 'axes'));
    ax = ax(arrayfun(@(a) ~isempty(findobj(a, 'type', 'image')), ax));
    n = numel(ax); c = cell(1, n); x = c; y = c; cl = c; xl = c; yl = c; xs = c; ys = c;
    for k = 1 : n
        im = findobj(ax(k), 'type', 'image'); sc = findobj(ax(k), 'type', 'scatter');
        c{k} = im.CData; x{k} = sc.XData; y{k} = sc.YData; cl{k} = ax(k).CLim;
        xl{k} = ax(k).XLim; yl{k} = ax(k).YLim; xs{k} = ax(k).XLabel.String;
        ys{k} = ax(k).YLabel.String;
    end
    cb = findobj(fig, 'type', 'colorbar');
    tx = findobj(fig, 'type', 'text'); ts = cell(1, 0);
    if ~isempty(tx), ts = {tx.String}; end
    fd = {c, x, y, cl, xl, yl, xs, ys, cb.Label.String, numel(cb), ts};
    """
).strip()
FIG_OUTPUTS = [
    "cdata",
    "x",
    "y",
    "clim",
    "xlim",
    "ylim",
    "xlabel",
    "ylabel",
    "cblabel",
    "ncolorbar",
    "texts",
]
PANELS = ("pk5d:otf_xy", "pk5d:otf_yz", "pk5d:otf_xz")


def _tex(s: str) -> str:
    """matplotlib's mathtext label as MATLAB's TeX text"""
    return s.replace("$", "")


def fig_data(fig) -> list:
    """FIG_DATA of a matplotlib figure"""
    axes = [a for a in fig.axes if a.get_gid() in PANELS]
    c, x, y, cl, xl, yl, xs, ys = ([] for _ in range(8))
    for a in axes:
        im = a.images[0]
        off = np.asarray(a.collections[0].get_offsets(), dtype=np.float64).reshape(-1, 2)
        c.append(np.ma.getdata(im.get_array()))
        x.append(off[:, 0].copy())
        y.append(off[:, 1].copy())
        cl.append([float(v) for v in im.get_clim()])
        xl.append([float(v) for v in a.get_xlim()])
        yl.append([float(v) for v in a.get_ylim()])
        xs.append(_tex(a.get_xlabel()))
        ys.append(_tex(a.get_ylabel()))
    cbs = [a for a in fig.axes if a.get_gid() == "pk5d:colorbar"]
    texts = [t.get_text() for a in fig.axes if a.get_gid() == "pk5d:title" for t in a.texts]
    return [c, x, y, cl, xl, yl, xs, ys, _tex(cbs[0].get_ylabel()), float(len(cbs)), texts]


def otf_volume(shape, name, double_mask=False):
    """(|OTF|, mask) in Python's form: ``1 / (1 + q)^2`` of the quadratic form q of the
    distance to the FFT center (plus a little noise), the mask ``q < 2.2`` (as double: 1 inside,
    2 in its inner part)"""
    rng = CounterRNG(0, f"b10/otf/{name}")
    grids = np.meshgrid(*[np.arange(n, dtype=np.float64) - n // 2 for n in shape], indexing="ij")
    scale = [0.12, 0.2, 0.3][-len(shape) :]
    q = sum((g * s) ** 2 for g, s in zip(grids, scale, strict=True))
    a = (1.0 / ((1.0 + q) * (1.0 + q)) * (1.0 + 0.01 * rng.uniform(shape))).astype(np.float32)
    m = q < 2.2
    if double_mask:
        m = m.astype(np.float64) + (q < 0.7)
    return a, m


def _otf_input(shape, name, which, double_mask=False):
    return Value(lambda: to_matlab(otf_volume(shape, name, double_mask)[which]))


def py_otf(abs_otf, mask, visible):
    pytest.importorskip("matplotlib")
    fig = figures.visualize_otf_and_mask_outline(abs_otf, mask, visible)
    return fig_data(fig)


def otf_case(cid, shape, double_mask=False):
    return Case(
        cid,
        code="fig = visualize_OTF_and_mask_outline(in{1}, in{2}, in{3}); "
        + FIG_DATA
        + " out = fd; close(fig);",
        args=[Vol(In("a")), Vol(In("m")), False],
        nargout=len(FIG_OUTPUTS),
        outputs=FIG_OUTPUTS,
        python=py_otf,
        inputs={
            "a": _otf_input(shape, cid, 0, double_mask),
            "m": _otf_input(shape, cid, 1, double_mask),
        },
        compare={"cdata": LOG10_ULP},
    )


OTF_CASES = [
    otf_case("otf.even", (22, 40, 36)),
    otf_case("otf.odd", (21, 41, 35)),
    otf_case("otf.double_mask", (20, 38, 34), double_mask=True),
    otf_case("otf.image", (40, 54)),
]

# XR_visualize_OTF_mask_segmentation on a PSF file
SEG_PNG = r".*figures/psf\.png$"


def psf_noisy() -> np.ndarray:
    """a Gaussian PSF with 5% deterministic noise (no tied |OTF| values at the mask's
    threshold)"""
    p = synth.decon_psf((31, 31, 31), (2.0, 1.5))
    u = CounterRNG(21, "b10/psf_noise").uniform(p.shape, np.float32)
    p = p * (np.float32(1) + np.float32(0.05) * u)
    return np.floor(p + np.float32(0.5)).astype(np.uint16)


def py_seg(psf_fn, thresh, skewed, *opts, backend):
    pytest.importorskip("matplotlib")
    from petakit5d import fft

    kw = dict(zip(opts[::2], opts[1::2], strict=True))
    with fft.using(backend):
        fig = figures.visualize_otf_mask_segmentation(
            psf_fn,
            thresh,
            None if isinstance(skewed, list) else skewed,
            min_int_thrsh=kw.get("minIntThrsh", 1e-3),
            visible=kw.get("visible", True),
            save_fig=kw.get("saveFig", False),
        )
    png = os.path.join(os.path.dirname(psf_fn), "figures", "psf.png")
    return [*fig_data(fig), os.path.isfile(png)]


def _clipped_close(limit):
    """the images as shown: values clipped to the color limits [-3, 0] within ``limit``"""

    def fn(ref, act) -> list[str]:
        if not isinstance(ref, MCell) or not isinstance(act, MCell):
            return ["not a cell"]
        if len(ref.items) != len(act.items):
            return [f"{len(ref.items)} images in MATLAB, {len(act.items)} in Python"]
        out = []
        for k, (r, a) in enumerate(zip(ref.items, act.items, strict=True)):
            if not isinstance(r, MArray) or not isinstance(a, MArray):
                out.append(f"image {k}: not an array")
                continue
            if r.mclass != a.mclass or tuple(r.size) != tuple(a.size):
                out.append(f"image {k}: {r.mclass} {r.size} against {a.mclass} {a.size}")
                continue
            rv = np.clip(np.asarray(r.data, np.float64), -3, 0)
            av = np.clip(np.asarray(a.data, np.float64), -3, 0)
            d = float(np.max(np.abs(rv - av))) if rv.size else 0.0
            if d > limit:
                out.append(f"image {k}: largest difference in [-3, 0] {d:g} > {limit:g}")
        return out

    return Custom(fn, f"images within {limit:g} in [-3, 0] (proposed: permissive FFT backend)")


def seg_cases(cid, *args):
    out = []
    for be in ("exact", "permissive"):
        out.append(
            Case(
                f"{cid}.{be}",
                code="XR_visualize_OTF_mask_segmentation(in{:}); fig = gcf; "
                + FIG_DATA
                + " out = [fd, {exist([fileparts(in{1}), '/figures/psf.png'], 'file') == 2}];"
                + " close(fig);",
                args=[InPath("psf"), *args],
                nargout=len(FIG_OUTPUTS) + 1,
                outputs=[*FIG_OUTPUTS, "png"],
                python=py_seg,
                py_kwargs={"backend": be},
                inputs={"psf": File(lambda p: synth.write_tiff(p, psf_noisy()), "psf.tif")},
                tree=tree(ignore=[SEG_PNG]),
                compare={"cdata": _clipped_close(1e-4) if be == "permissive" else LOG10_ULP},
            )
        )
    return out


SEG_CASES = [
    *seg_cases("seg.unskewed_saved", 0.88, False, "visible", True, "saveFig", True),
    *seg_cases("seg.auto_skewed", 0.65, [], "minIntThrsh", 1e-3, "visible", True),
]

# XR_visualize_tile_positions on an image list
TILE = (
    "Scan_Iter_{it:04d}{cam}_ch{ch}_CAM1_stack0000_{nm}nm_0000000msec_0001234567msecAbs_"
    "{x:03d}x_{y:03d}y_{z:03d}z_0000t.tif"
)


def image_list(cam=True, sub_iter=False, mixed=False):
    """an image list csv; sub_iter: names with a sub-iteration (Scan_Iter_0000_0000_, as the
    demo's), mixed: names with and without a camera"""

    def make(path):
        lines = ["Filepath,Filename,StageX (um),StageY (um),StageZ (um)"]
        channels = (("_CamA", 0), ("_CamB", 0), ("_CamA", 1)) if cam else (("", 0), ("", 1))
        if mixed:
            channels = (("_CamA", 0), ("", 1))
        for it in (0, 1):
            for c, ch in channels:
                for x, y, z in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 1)):
                    name = TILE.format(it=it, cam=c, ch=ch, nm=488 + 72 * ch, x=x, y=y, z=z)
                    if sub_iter:
                        name = name.replace(f"Iter_{it:04d}_", f"Iter_{it:04d}_0000_")
                    sx, sy, sz = 100.5 * x + 0.25 * it, 80.25 * y + 0.5 * ch, 20.125 * z
                    lines.append(f"/data/{name},{name},{sx},{sy},{sz}")
        # LF on every platform: the case keys hash the inputs' bytes
        Path(path).write_text("\n".join(lines) + "\n", newline="\n")

    return make


TILE_OUTPUTS = ["x", "y", "z", "texts", "positions", "view", "title", "xl", "yl", "zl", "eps"]


def py_tiles(fn):
    pytest.importorskip("matplotlib")
    fig = figures.visualize_tile_positions(fn)
    ax = next(a for a in fig.axes if a.get_gid() == "pk5d:tiles")
    xs, ys, zs = (np.asarray(v, dtype=np.float64) for v in ax.collections[0]._offsets3d)
    texts = [t.get_text() for t in ax.texts]
    pos = np.array([t.get_position_3d() for t in ax.texts], dtype=np.float64).reshape(-1, 3)
    eps = os.path.join(os.path.dirname(fn), "visualize_tile_positions.eps")
    return [
        xs,
        ys,
        zs,
        texts,
        pos,
        [float(ax.azim) + 90, float(ax.elev)],
        ax.get_title(),
        ax.get_xlabel(),
        ax.get_ylabel(),
        ax.get_zlabel(),
        os.path.isfile(eps),
    ]


def tiles_case(cid, cam, **kw):
    make = image_list(cam, kw.pop("sub_iter", False), kw.pop("mixed", False))
    return Case(
        cid,
        code=(
            "XR_visualize_tile_positions(in{1}); fig = gcf; ax = fig.CurrentAxes; "
            "s = findobj(ax, 'type', 'scatter'); tx = flipud(findobj(ax, 'type', 'text')); "
            "out = {s.XData, s.YData, s.ZData, {tx.String}, vertcat(tx.Position), ax.View, "
            "ax.Title.String, ax.XLabel.String, ax.YLabel.String, ax.ZLabel.String, "
            "exist([fileparts(in{1}), '/visualize_tile_positions.eps'], 'file') == 2}; close(fig);"
        ),
        args=[InPath("list")],
        nargout=len(TILE_OUTPUTS),
        outputs=TILE_OUTPUTS,
        python=py_tiles,
        inputs={"list": File(make, "ImageList_from_encoder.csv")},
        tree=tree(ignore=[r".*\.eps$"]),
        **kw,
    )


TILE_CASES = [
    tiles_case("tiles.cam", True),
    tiles_case("tiles.no_cam", False),
    # MATLAB's expression has no sub-iteration: the names don't match and it fails
    tiles_case(
        "tiles.sub_iter_fails",
        True,
        sub_iter=True,
        expect_error=True,
        py_error=(figures.FigureError,),
    ),
    # neither expression applies to every name: MATLAB's expression is undefined
    tiles_case(
        "tiles.mixed_fails", True, mixed=True, expect_error=True, py_error=(figures.FigureError,)
    ),
]

CASES = _single_cases() + _wrapper_cases() + PROBES + OTF_CASES + SEG_CASES + TILE_CASES


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    if case.id.endswith(".exact"):
        reason = digests_hold_here()
        if reason:
            pytest.skip(reason)
    if os.name == "nt" and case.id.startswith("ds.") and not case.expect_error:
        # simplifyPath gives the data and result folders with '\' on Windows (from dir, as
        # MATLAB's there), so the fn and fnout the _bbox.mat files hold are '<WORK>\d/...' and
        # '<WORK>\d\masks/...'; these references are MATLAB's on Linux, with '/'
        pytest.skip("the fn and fnout in _bbox.mat hold Windows folders; references are Linux's")
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# DeconFrameTest.otfMaskVisualizationSavesTheFigure, in intent
# ---------------------------------------------------------------------------------------------


@needs_mpl
def test_otf_mask_visualization_saves_the_figure(tmp_path, capsys):
    psf_fn = str(tmp_path / "psf.tif")
    synth.write_tiff(psf_fn, psf_noisy())
    r = figures.visualize_otf_mask_segmentation(psf_fn, 0.88, False, visible=False, save_fig=True)
    assert r is None
    png = tmp_path / "figures" / "psf.png"
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    # automatic skewed space detection and a minimum intensity threshold
    fig = figures.visualize_otf_mask_segmentation(psf_fn, 0.65, None, min_int_thrsh=1e-3)
    assert [a.get_gid() for a in fig.axes].count("pk5d:title") == 1
    assert "OTFCumThresh: 0.650000, skewed: " in capsys.readouterr().out


# ---------------------------------------------------------------------------------------------
# the masks
# ---------------------------------------------------------------------------------------------


def test_mask_does_not_depend_on_threads():
    im = mip_image((90, 130), "threads")
    ref = mip_masks.mip_mask(im, 100, 120, (3, 7), num_threads=1)
    for t in (2, 5):
        np.testing.assert_array_equal(mip_masks.mip_mask(im, 100, 120, (3, 7), num_threads=t), ref)


def test_mask_follows_the_image_axes():
    """dilate_size is (rows, columns): a non-square image transposed gives the transposed mask
    with the sizes swapped"""
    im = mip_image((80, 120), "axes", n_blobs=3)
    a = mip_masks.mip_mask(im, 100, 50, (2, 9))
    b = mip_masks.mip_mask(np.ascontiguousarray(im.T), 100, 50, (9, 2))
    np.testing.assert_array_equal(a.T, b)
    assert not np.array_equal(a, mip_masks.mip_mask(im, 100, 50, (9, 2)))


def test_mask_with_one_dilate_size_fails_as_matlab():
    with pytest.raises(mip_masks.MipMaskError) as e:
        mip_masks.mip_mask(mip_image((40, 50), "one"), 100, 10, (3,))
    assert e.value.identifier == "MATLAB:badsubscript"


def test_publish_order_box_and_figure_before_the_mask(tmp_path, monkeypatch):
    """the box file (and the figure) exist when the mask is renamed into place"""
    fn = str(tmp_path / "im.tif")
    synth.write_tiff(fn, mip_image((60, 70), "order"))
    seen = []
    real = os.replace

    def spy(src, dst):
        if str(dst).endswith("m.tif"):
            seen.append(sorted(p.name for p in tmp_path.iterdir()))
        return real(src, dst)

    monkeypatch.setattr(mip_masks.os, "replace", spy)
    mip_masks.generate_single_mip_mask(fn, str(tmp_path / "m.tif"), 100, 50, (2, 2), uuid="u1")
    assert seen and "m_bbox.mat" in seen[0]
    if HAS_MPL:
        assert "m_plot.png" in seen[0]


def test_no_figure_without_matplotlib(tmp_path, monkeypatch):
    monkeypatch.setattr(figures, "available", lambda: False)
    fn = str(tmp_path / "im.tif")
    synth.write_tiff(fn, np.zeros((30, 40), np.uint16))
    # (an all-zero MIP: no figure, so MATLAB's clim error does not happen)
    mip_masks.generate_single_mip_mask(fn, str(tmp_path / "m.tif"), 100, 10, (1, 1))
    names = sorted(p.name for p in tmp_path.iterdir())
    assert names == ["im.tif", "m.tif", "m_bbox.mat"]


def _mask_files(folder):
    """the files of a mask call and their contents (the mask's voxels, the box file's values
    but its paths)"""
    from petakit5d.io import read_mat, read_tiff

    out = {}
    for p in sorted(Path(folder).rglob("*")):
        if p.suffix == ".tif" and "_MIP_" not in p.name and p.name != "im.tif":
            out[p.relative_to(folder).as_posix()] = read_tiff(str(p)).tobytes()
        elif p.name.endswith("_bbox.mat"):
            m = read_mat(str(p))
            out[p.relative_to(folder).as_posix()] = repr(
                {k: np.asarray(v).tolist() for k, v in m.items() if k not in ("fn", "fnout")}
            )
        elif p.is_file():
            out[p.relative_to(folder).as_posix()] = None
    return out


@needs_mpl
def test_job_tasks_draw_no_figure_direct_calls_do(tmp_path):
    """MATLAB's job tasks run -nojvm: a task of the job framework saves no _plot.png, a direct
    call does; the other files are the same"""
    from petakit5d.jobs.tasks import Task, run_record

    a, b = tmp_path / "direct", tmp_path / "task"
    for d in (a, b):
        d.mkdir()
        synth.write_tiff(str(d / "im.tif"), mip_image((60, 70), "jobs"))
    args = (100, 50, (2, 2), (1, 256, 256), True, "u1")
    mip_masks.generate_single_mip_mask(str(a / "im.tif"), str(a / "m.tif"), *args)
    t = Task.of(mip_masks.generate_single_mip_mask, str(b / "im.tif"), str(b / "m.tif"), *args)
    run_record(t.record())
    assert (a / "m_plot.png").is_file() and not (b / "m_plot.png").exists()
    fa, fb = _mask_files(a), _mask_files(b)
    assert fa.pop("m_plot.png") is None
    assert fa == fb and set(fb) == {"im.tif", "m.tif", "m_bbox.mat"}
    # an all-zero MIP: the figure's clim fails in a direct call, not in a task (no figure)
    z = tmp_path / "zero"
    z.mkdir()
    synth.write_tiff(str(z / "im.tif"), np.zeros((30, 40), np.uint16))
    with pytest.raises(figures.FigureError):
        mip_masks.generate_single_mip_mask(str(z / "im.tif"), str(z / "d.tif"), 100, 10, (1, 1))
    t = Task.of(mip_masks.generate_single_mip_mask, str(z / "im.tif"), str(z / "t.tif"), 100, 10)
    run_record(t.record())
    assert (z / "t.tif").is_file() and not (z / "t_plot.png").exists()


@needs_mpl
def test_dataset_tasks_draw_no_figures(tmp_path):
    """generate_mip_masks runs its masks as framework tasks: no figures (in local mode here),
    and the masks are those of direct calls"""
    d = tmp_path / "d"
    d.mkdir()
    mip_dataset()(str(d))
    mip_masks.generate_mip_masks(
        str(d), int_thresh=(100, 110, 120), vol_thresh=150, dilate_size=(2, 3, 4), uuid="u2",
        parse_cluster=False,
    )  # fmt: skip
    masks = d / "masks"
    assert sorted(masks.glob("*_plot.png")) == []
    made = sorted(p.name for p in masks.glob("*_MIP_*.tif"))
    assert len(made) == 6
    # the same mask by a direct call: with its figure, the same voxels
    from petakit5d.io import read_tiff

    e = tmp_path / "e"
    e.mkdir()
    name = "scan_CamA_ch0_t0000_MIP_z.tif"
    # the z MIP: the threshold of z, the dilate sizes of (y, x) (see the module)
    mip_masks.generate_single_mip_mask(str(d / name), str(e / name), 100, 150, (3, 4))
    assert (e / name.replace(".tif", "_plot.png")).is_file()
    np.testing.assert_array_equal(read_tiff(str(e / name)), read_tiff(str(masks / name)))


def test_modules_import_without_matplotlib():
    code = (
        "import sys; sys.modules['matplotlib'] = None\n"
        "from petakit5d.tools import figures, mip_masks\n"
        "import petakit5d.matlab as m\n"
        "m.XR_generate_MIP_masks, m.XR_visualize_OTF_mask_segmentation\n"
    )
    subprocess.run([sys.executable, "-c", code], check=True, env=dict(os.environ))


# ---------------------------------------------------------------------------------------------
# the figures
# ---------------------------------------------------------------------------------------------


@needs_mpl
def test_mip_mask_figure_shows_the_mip_and_the_mask():
    from petakit5d.matlabcompat.stats import prctile

    im = mip_image((50, 70), "fig")
    bw = mip_masks.mip_mask(im, 100, 20, (1, 1))
    fig = figures.mip_mask_figure(im, bw)
    a, b = (next(x for x in fig.axes if x.get_gid() == g) for g in ("pk5d:mip", "pk5d:mask"))
    np.testing.assert_array_equal(np.ma.getdata(a.images[0].get_array()), im)
    np.testing.assert_array_equal(np.ma.getdata(b.images[0].get_array()), bw)
    assert a.images[0].get_clim() == (0, float(prctile(im, 99)[0]))
    assert a.get_xlim() == (0, 70) and a.get_ylim() == (50, 0)


@needs_mpl
def test_figures_are_saved_without_a_display(tmp_path):
    """no pyplot for figures that are not shown, so no display or GUI backend is needed"""
    code = textwrap.dedent(
        f"""
        import os, sys
        os.environ.pop("DISPLAY", None); os.environ.pop("MPLBACKEND", None)
        import numpy as np
        from petakit5d.tools import figures
        a = np.ones((6, 8, 10), np.float32); m = np.zeros(a.shape, bool); m[2:4, 3:6, 4:7] = True
        fig = figures.visualize_otf_and_mask_outline(a, m, False)
        fig.savefig({str(tmp_path / "f.png")!r})
        assert "matplotlib.pyplot" not in sys.modules
        """
    )
    subprocess.run([sys.executable, "-c", code], check=True)
    assert (tmp_path / "f.png").read_bytes()[:4] == b"\x89PNG"


@needs_mpl
def test_tile_positions_figure_is_eps(tmp_path):
    fn = tmp_path / "ImageList_from_encoder.csv"
    image_list(True)(fn)
    figures.visualize_tile_positions(str(fn))
    assert (tmp_path / "visualize_tile_positions.eps").read_bytes()[:10] == b"%!PS-Adobe"


def test_tile_positions_without_cam(tmp_path):
    fn = tmp_path / "list.csv"
    image_list(False)(fn)
    tp = figures.tile_positions(str(fn))
    assert tp["xyz"].shape == (4, 3)
    assert tp["names"][3] == "001x_001y_001z"
    np.testing.assert_array_equal(tp["xyz"][:, 0], [0, 100.5, 0, 100.5])


def test_segmentation_title_as_matlab_prints_it():
    assert figures._segmentation_title(0.9, None) == "OTFCumThresh: 0.900000, skewed: "
    assert figures._segmentation_title(0.65, True) == "OTFCumThresh: 0.650000, skewed: true"


def test_missing_psf_file_fails_as_matlab(tmp_path):
    with pytest.raises(figures.FigureError, match="does not exist"):
        figures.visualize_otf_mask_segmentation(str(tmp_path / "none.tif"))
    (tmp_path / "p.mat").write_bytes(b"")
    with pytest.raises(figures.FigureError, match="Unknown format"):
        figures.visualize_otf_mask_segmentation(str(tmp_path / "p.mat"))
