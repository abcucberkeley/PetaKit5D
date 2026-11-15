"""Figure tools (MATLAB's ``XR_visualize_tile_positions``, ``visualize_OTF_and_mask_outline``,
``XR_visualize_OTF_mask_segmentation`` and the figure of ``gererate_single_MIP_mask``), drawn with
matplotlib, an optional dependency (``pip install matplotlib``): without it ``available()`` is
False, the MIP masks are saved without their figure (as MATLAB without its JVM) and the other
functions raise ``ImportError``. Figures that are not shown are made without pyplot, so they
are saved on machines without a display.

The figures show MATLAB's data: the same images (``imagesc``'s CData, rows and columns as MATLAB
draws them), points, limits, color limits, labels and titles; fonts, colors and layout follow
MATLAB's where matplotlib has them, pixels are not the same. MATLAB's TeX labels are mathtext
here (``k_x`` is ``$k_x$``, ``log_{10}`` is ``log$_{10}$``). Each data axes carries a gid
(``pk5d:otf_xy``, ``pk5d:otf_yz``, ``pk5d:otf_xz``, ``pk5d:mip``, ``pk5d:mask``,
``pk5d:tiles``), the title axes ``pk5d:title``.

- ``visualize_otf_and_mask_outline(abs_otf, otf_mask)``: the log10 of the fftshifted |OTF|
  (``(z, y, x)``, or an image ``(y, x)``) normalized by each slice's maximum, through its center
  (``size // 2``) in xy, yz and xz, in [-3, 0] (``hot``), with the outline of the OTF mask
  (``mask & ~imerode(mask, strel('disk', 1))``, 1-based points in MATLAB's ``find`` order).
  The xy panel shows x along its rows: MATLAB's ``log10(...)'`` of the ``(y, x)`` slice.
- ``visualize_otf_mask_segmentation(psf_fn, ...)``: the PSF file cleaned (``psf_gen``: dz 0.5,
  median factor 1.5, 'masked') and normalized, its OMW back projector's |OTF| and mask
  (alpha 0.01, Hann bounds (0.8, 1.0)) in that figure, titled with the threshold;
  ``<psf folder>/figures/<name>.png`` with ``save_fig``.
- ``visualize_tile_positions(image_list)``: the stage positions of the tiles of the first
  channel and camera at iteration 0 of an image list csv file, the points at x + 7, y + 7 (as
  MATLAB plots them) and the names ``<x>x_<y>y_<z>z`` at the positions, saved as
  ``<folder>/visualize_tile_positions.eps``.
- ``mip_mask_figure(im, bw)``: a MIP (color limits [0, its 99th percentile]) and its mask.

Kept from MATLAB: a MIP whose 99th percentile is not above 0 (an all-zero MIP) makes MATLAB's
``clim`` fail, and so the mask task, when the figure is drawn.
"""

from __future__ import annotations

import os
import re
from typing import Any

import numpy as np

from petakit5d.params import Param, register

__all__ = [
    "FigureError",
    "available",
    "mip_mask_figure",
    "save_mip_mask_figure",
    "visualize_otf_and_mask_outline",
    "visualize_otf_mask_segmentation",
    "visualize_tile_positions",
]

_FONT = "Helvetica"
_FONT_SIZE = 14


class FigureError(ValueError):
    """an error MATLAB raises, with its identifier"""

    def __init__(self, message: str, identifier: str = ""):
        super().__init__(message)
        self.identifier = identifier


def available() -> bool:
    """matplotlib can be imported"""
    import importlib.util

    return importlib.util.find_spec("matplotlib") is not None


def _figure(width_px: int, height_px: int, visible: bool = False):
    """a white figure of MATLAB's size in pixels (100 per inch); through pyplot when it is to be
    shown"""
    size = (width_px / 100, height_px / 100)
    if visible:
        import matplotlib.pyplot as plt

        return plt.figure(figsize=size, dpi=100, facecolor="w")
    from matplotlib.figure import Figure

    return Figure(figsize=size, dpi=100, facecolor="w")


def _save(fig, path: str, **kw: Any) -> None:
    """save under a temporary name in the same folder, then rename"""
    from petakit5d.core.ids import get_uuid

    root, ext = os.path.splitext(path)
    tmp = f"{root}_{get_uuid()}{ext}"
    fig.savefig(tmp, format=ext[1:], facecolor=fig.get_facecolor(), **kw)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------------------------
# the figure of a MIP mask
# ---------------------------------------------------------------------------------------------


def _clim_hi(im: np.ndarray) -> float:
    from petakit5d.matlabcompat.stats import prctile

    return float(np.asarray(prctile(np.asarray(im), 99)).reshape(-1)[0])


def mip_mask_figure(im: np.ndarray, bw: np.ndarray):
    """``gererate_single_MIP_mask``'s figure: the MIP (``imagesc``, color limits [0, its 99th
    percentile]) and its mask, side by side, the axes limits [0, nx] and [0, ny]"""
    a = np.asarray(im)
    if a.ndim != 2:
        raise NotImplementedError("the figure of a MIP mask needs a 2D image")
    hi = _clim_hi(a)
    if not hi > 0:
        raise FigureError(
            "Limits must be a 2-element vector of increasing numeric values.",
            "MATLAB:rulerFunctions:InvalidNumericLimits",
        )
    ny, nx = a.shape
    fig = _figure(1200, 600)
    for k, (data, gid) in enumerate(((a, "pk5d:mip"), (np.asarray(bw), "pk5d:mask"))):
        ax = fig.add_subplot(1, 2, k + 1)
        ax.set_gid(gid)
        img = ax.imshow(
            data,
            cmap="viridis",
            origin="upper",
            extent=(0.5, nx + 0.5, ny + 0.5, 0.5),
            interpolation="nearest",
        )
        if k == 0:
            img.set_clim(0, hi)
        else:
            img.set_clim(0, 1)
        ax.set_aspect("equal")
        ax.set_xlim(0, nx)
        ax.set_ylim(ny, 0)
    return fig


def save_mip_mask_figure(im: np.ndarray, bw: np.ndarray, fig_fn: str) -> None:
    """``mip_mask_figure`` saved as png (screen resolution, MATLAB's ``-r0``)"""
    fig = mip_mask_figure(im, bw)
    _save(fig, fig_fn, dpi=100)


# ---------------------------------------------------------------------------------------------
# OTF and OTF mask
# ---------------------------------------------------------------------------------------------


def _log_normalized(s: np.ndarray) -> np.ndarray:
    """``log10(s ./ max(s(:)))`` of a single slice (the quotient in single, its log10 in double
    rounded to single)"""
    from petakit5d.core.dtypes import NAN

    a = np.asarray(s)
    with np.errstate(divide="ignore", invalid="ignore"):
        finite = a[~np.isnan(a)]
        m = finite.max() if finite.size else a.dtype.type(NAN)
        q = a / m
        out = np.log10(q.astype(np.float64)).astype(a.dtype)
    return out


def _outline(mask2d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``[y, x] = find(m & ~imerode(m, strel('disk', 1)))`` of a mask in MATLAB's rows and
    columns: the 1-based rows and columns (double) in column-major order"""
    from petakit5d.matlabcompat.morphology import imerode, strel

    m = np.asarray(mask2d) != 0
    supp = m & ~imerode(m, strel("disk", 1))
    cols, rows = np.nonzero(supp.T)
    return rows.astype(np.float64) + 1, cols.astype(np.float64) + 1


def _otf_panels(abs_otf: np.ndarray, otf_mask: np.ndarray) -> list[dict[str, Any]]:
    """per panel: the image (MATLAB's CData), the outline points and the labels"""
    a = np.asarray(abs_otf)
    m = np.asarray(otf_mask)
    if a.ndim == 2:
        a, m = a[None], m[None]
    if a.ndim != 3 or m.shape != a.shape:
        raise FigureError("the |OTF| and its mask: (z, y, x) volumes or (y, x) images of one size")
    nz, ny, nx = a.shape
    cz, cy, cx = nz // 2, ny // 2, nx // 2
    # (MATLAB's slices in its rows and columns; the images are their transposes)
    slices = [
        ("pk5d:otf_xy", a[cz].T, m[cz], "$k_x$", "$k_y$"),
        ("pk5d:otf_yz", a[:, :, cx], m[:, :, cx].T, "$k_x$", "$k_z$"),
        ("pk5d:otf_xz", a[:, cy, :], m[:, cy, :].T, "$k_y$", "$k_z$"),
    ]
    out = []
    for gid, img, mask2d, xl, yl in slices[: 1 if nz == 1 else 3]:
        x, y = _outline(mask2d)
        out.append(
            {
                "gid": gid,
                "cdata": _log_normalized(np.ascontiguousarray(img)),
                "x": x,
                "y": y,
                "xlabel": xl,
                "ylabel": yl,
            }
        )
    return out


def visualize_otf_and_mask_outline(abs_otf: np.ndarray, otf_mask: np.ndarray, visible: bool = True):
    """``visualize_OTF_and_mask_outline``: the figure (see the module) of an fftshifted |OTF| and
    its OTF mask (``(z, y, x)``, or images ``(y, x)``); returns the matplotlib figure. With
    ``visible`` it is a pyplot figure (shown by an interactive backend)."""
    panels = _otf_panels(abs_otf, otf_mask)
    a = np.asarray(abs_otf)
    nz, ny, nx = (1, *a.shape) if a.ndim == 2 else a.shape
    ws = np.array([ny, ny, nx], dtype=np.float64)
    ws = ws / ws.sum() * 0.85
    hs = np.array([nx, nz, nz], dtype=np.float64)
    hs = hs / hs.max() * 0.7
    aps = [
        [0.04, 0.2, ws[0], hs[0]],
        [ws[0] + 0.075, 0.2, ws[1], hs[1]],
        [ws[0] + ws[1] + 0.11, 0.2, ws[2], hs[2]],
    ]
    if nz == 1:
        aps = [[0.04, 0.2, 0.8, 0.7]]
    fig = _figure(1920, 1080, visible) if nz != 1 else _figure(800, 800, visible)
    img = None
    ax = None
    for p, pos in zip(panels, aps, strict=False):
        ax = fig.add_axes(pos)
        ax.set_gid(p["gid"])
        rows, cols = p["cdata"].shape
        img = ax.imshow(
            p["cdata"],
            cmap="hot",
            vmin=-3,
            vmax=0,
            origin="lower",
            extent=(0.5, cols + 0.5, 0.5, rows + 0.5),
            interpolation="nearest",
        )
        ax.scatter(p["x"], p["y"], marker=".")
        ax.set_aspect("equal")
        ax.set_xlim(1, cols)
        ax.set_ylim(1, rows)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(p["xlabel"], fontname=_FONT, fontsize=_FONT_SIZE)
        ax.set_ylabel(p["ylabel"], fontname=_FONT, fontsize=_FONT_SIZE)
        for s in ax.spines.values():
            s.set_linewidth(1)
        ax.tick_params(direction="out")
    cb = fig.colorbar(img, ax=ax)
    cb.ax.set_gid("pk5d:colorbar")
    cb.set_label("log$_{10}$(OTF Intensity)", fontname=_FONT, fontsize=_FONT_SIZE)
    cb.ax.tick_params(labelsize=_FONT_SIZE, direction="out")
    cb.outline.set_linewidth(1)
    return fig


def _segmentation_title(otf_cum_thresh: float, skewed: bool | None) -> str:
    """``sprintf('OTFCumThresh: %f, skewed: %s', OTFCumThresh, string(skewed))``"""
    from petakit5d.core.text import sprintf

    sk: Any = "" if skewed is None else ("true" if skewed else "false")
    return sprintf("OTFCumThresh: %f, skewed: %s", float(otf_cum_thresh), sk)


def _read_psf(psf_fn: str) -> np.ndarray:
    from petakit5d.core.paths import fileparts
    from petakit5d.jobs.files import exist

    if not exist(psf_fn):
        raise FigureError(
            f"The PSF file {psf_fn} does not exist, please double check to make sure the path "
            "is correct!\n"
        )
    ext = fileparts(psf_fn)[2]
    if ext in (".tif", ".tiff"):
        from petakit5d.io import read_tiff

        p = read_tiff(psf_fn)
    elif ext == ".zarr":
        from petakit5d.io import read_zarr

        p = read_zarr(psf_fn)
    else:
        raise FigureError(f"Unknown format for PSF file {psf_fn}\n!")
    from petakit5d.core.dtypes import single

    return np.ascontiguousarray(single(p))


@register(
    "XR_visualize_OTF_mask_segmentation",
    Param("psfFn", "@ischar", "PSF file (tiff or zarr)."),
    Param(
        "OTFCumThresh",
        "@(x) isscalar(x)",
        "Accumulated energy quantile to segment the OTF mask (smaller: tighter masks).",
    ),
    Param(
        "skewed",
        "@(x) isempty(x) || islogical(x)",
        "PSF in the skewed space; None (MATLAB's []): detected.",
    ),
    Param("minIntThrsh", "@(x) isnumeric(x)", "Minimum intensity threshold of the segmentation."),
    Param("visible", "@(x) islogical(x)", "Make the figure visible (else it is closed)."),
    Param("saveFig", "@(x) islogical(x)", "Save the figure to <psf folder>/figures/<name>.png."),
    task="B10",
)
def visualize_otf_mask_segmentation(
    psf_fn: str,
    otf_cum_thresh: float = 0.9,
    skewed: bool | None = None,
    *,
    min_int_thrsh: float = 0.001,
    visible: bool = True,
    save_fig: bool = False,
):
    """The OTF mask outline segmented with ``otf_cum_thresh`` over the OTF of a PSF file
    (``XR_visualize_OTF_mask_segmentation``; see the module). Returns the matplotlib figure, or
    None when it is not ``visible`` (closed, as MATLAB closes it)."""
    from petakit5d.core.paths import fileparts
    from petakit5d.decon.backprojector import generate_backprojector
    from petakit5d.decon.psf import generate_psf
    from petakit5d.matlabcompat.reduce import sum as msum

    psf = _read_psf(psf_fn)
    psf = generate_psf(psf, 0.5, 0.5, 1.5, "masked")
    psf = np.asarray(psf, dtype=np.float32)
    with np.errstate(divide="ignore", invalid="ignore"):
        psf = psf / np.float32(msum(psf))

    title = _segmentation_title(otf_cum_thresh, skewed)
    print(title, flush=True)
    bp = generate_backprojector(
        psf,
        0.01,
        skewed,
        otf_cum_thresh=otf_cum_thresh,
        hann_win_bounds=(0.8, 1.0),
        min_int_thrsh=min_int_thrsh,
    )
    fig = visualize_otf_and_mask_outline(bp.abs_otf_c, bp.otf_mask, visible)

    tax = fig.add_axes([0, 0.95, 1, 0.05])
    tax.set_gid("pk5d:title")
    tax.set_facecolor("none")
    tax.set_xticks([])
    tax.set_yticks([])
    for s in tax.spines.values():
        s.set_color("white")
    tax.text(
        0.5,
        0,
        title,
        fontsize=14,
        fontweight="bold",
        horizontalalignment="center",
        verticalalignment="bottom",
        transform=tax.transData,
    )

    if save_fig:
        pstr, fsn, _ = fileparts(psf_fn)
        fig_path = f"{pstr}/figures/"
        os.makedirs(fig_path, exist_ok=True)
        _save(fig, f"{fig_path}/{fsn}.png", dpi=150)

    # (a figure that is not shown is not pyplot's: nothing else holds it)
    return fig if visible else None


# ---------------------------------------------------------------------------------------------
# tile positions
# ---------------------------------------------------------------------------------------------

_CAM_RX = (
    r"Scan_Iter_(?P<Iter>\d+)_Cam(?P<Cam>\w+)_ch(?P<ch>\d+)_CAM1_stack(?P<stack>\d+)_"
    r"(?P<laser>\d+)nm_(?P<abstime>\d+)msec_(?P<fpgatime>\d+)msecAbs_(?P<x>\d+)x_(?P<y>\d+)y_"
    r"(?P<z>\d+)z_(?P<t>\d+)t.tif"
)
_NOCAM_RX = (
    r"Scan_Iter_(?P<Iter>\d+)_ch(?P<ch>\d+)_CAM1_stack(?P<stack>\d+)_(?P<laser>\d+)nm_"
    r"(?P<abstime>\d+)msec_(?P<fpgatime>\d+)msecAbs_(?P<x>\d+)x_(?P<y>\d+)y_(?P<z>\d+)z_"
    r"(?P<t>\d+)t.tif"
)


def tile_positions(image_list_file_name: str) -> dict[str, Any]:
    """What ``XR_visualize_tile_positions`` plots: ``xyz`` (n x 3, the stage x, y, z of the
    tiles of iteration 0, the first row's channel and camera) and ``names``
    (``<x>x_<y>y_<z>z``, 3 digits each)."""
    from petakit5d.core.text import sprintf, str2double
    from petakit5d.stitch.imagelist import read_table

    tab = read_table(image_list_file_name, ",")
    fn = [str(v) for v in tab["Filename"]]
    flags = re.ASCII
    if all(re.search(r"_Cam\w_ch", f, flags) for f in fn):
        rx, specify_cam = _CAM_RX, True
    elif all(re.search(r"_ch[0-9]_", f, flags) for f in fn):
        rx, specify_cam = _NOCAM_RX, False
    else:
        raise FigureError(
            "Unrecognized function or variable 'expression'.", "MATLAB:UndefinedFunction"
        )
    rows = []
    for f in fn:
        m = re.search(rx, f, flags | re.IGNORECASE)
        if m is None:
            raise FigureError("Not enough input arguments.", "MATLAB:minrhs")
        g = m.groupdict()
        cam = g["Cam"] if specify_cam else "N"
        if len(cam) != 1:
            raise FigureError(
                "Unable to perform assignment because the left and right sides have a different "
                "number of elements.",
                "MATLAB:index:assignmentDimensionMismatch",
            )
        r = {k: str2double(v) for k, v in g.items() if k != "Cam"}
        r["Cam"] = cam
        rows.append(r)
    it = np.array([r["Iter"] for r in rows])
    ch = np.array([r["ch"] for r in rows])
    cam = [r["Cam"] for r in rows]
    sel = [k for k in range(len(rows)) if it[k] == 0 and ch[k] == ch[0] and cam[k] == cam[0]]
    cols = [
        np.asarray(tab[n], dtype=np.float64) for n in ("StageX_um_", "StageY_um_", "StageZ_um_")
    ]
    xyz = np.stack([c[sel] for c in cols], axis=1) if sel else np.zeros((0, 3))
    names = [sprintf("%03dx_%03dy_%03dz", rows[k]["x"], rows[k]["y"], rows[k]["z"]) for k in sel]
    return {"xyz": xyz, "names": names}


def visualize_tile_positions(image_list_file_name: str):
    """``XR_visualize_tile_positions``: a 3D scatter plot of the tile positions of an image list
    csv file (see the module), saved as ``<its folder>/visualize_tile_positions.eps``; returns
    the matplotlib figure (MATLAB leaves it open)."""
    from petakit5d.core.paths import fileparts

    tp = tile_positions(image_list_file_name)
    xyz = tp["xyz"]
    fig = _figure(1200, 1000)
    ax = fig.add_axes([0.15, 0.15, 0.7, 0.7], projection="3d")
    ax.set_gid("pk5d:tiles")
    ax.scatter(xyz[:, 0] + 7, xyz[:, 1] + 7, xyz[:, 2], s=50, marker="*")
    for (x, y, z), txt in zip(xyz, tp["names"], strict=True):
        ax.text(x, y, z, txt, fontsize=18)
    ax.set_xlabel("x", fontsize=18)
    ax.set_ylabel("y", fontsize=18)
    ax.set_zlabel("z", fontsize=18)
    # MATLAB's view([-45, 45]): azimuth from -y, matplotlib's from +x
    ax.view_init(elev=45, azim=-45 - 90)
    ax.set_title("Tile positions", fontsize=18)
    ax.tick_params(labelsize=18)

    figure_location = fileparts(image_list_file_name)[0]
    _save(fig, f"{figure_location}/visualize_tile_positions.eps")
    return fig
