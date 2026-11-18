"""NaN box ends: taken only where MATLAB takes them (the user's decision of 2026-11-18).

MATLAB treats a NaN box end unlike Inf: ``min(bbox(4 : 6), imSize)`` makes both the image's
end, ``NaN > imSize`` is false where ``Inf > imSize`` is true, and a size ``NaN`` fails where it
is used. The port keeps a NaN end as ``math.nan`` (MATLAB's NaN on the way back) and:

- the crop functions take it where MATLAB does (XR_crop_frame of tiff and zarr files, padded or
  not; XR_crop_dataset 'fixed' and 'moving'), with the open (Inf) end beside it for contrast;
- every other function raises MATLAB's error: ``MATLAB:colon:nonFiniteEndpoint`` where the end
  reaches MATLAB's indexing, ``MATLAB:nonaninf`` where the size of the box makes the batches of
  a zarr file, ``MATLAB:imref3d:expectedInteger`` for deskewRotateFrame3D's output reference;
  for the large-file crop, XR_crop_dataset 'center' and padding by ``Inf``, MATLAB's errors of
  those paths.

Oracle cases against MATLAB (the error cases compare the error's identifier and message, not the
files left behind),
then the conversions of ``petakit5d.core.bbox`` and the helpers' refusal.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from petakit5d.core import bbox as B
from petakit5d.core.dtypes import NAN
from petakit5d.testing import oracle, synth
from petakit5d.testing.cases import Case, File, InPath, Tree, Work
from petakit5d.testing.synth import CounterRNG

NANB = [2, 3, 1, math.nan, 30, 8]  # MATLAB's [ymin xmin zmin ymax xmax zmax], ymax NaN
NANZ = [2, 3, 1, 35, 30, math.nan]  # zmax NaN
INFB = [2, 3, 1, math.inf, 30, 8]
IM = (12, 40, 50)  # (z, y, x); MATLAB's [40 50 12]
LOCAL = ["parseCluster", False, "masterCompute", True]
SMALL = ["batchSize", [16, 16, 4], "blockSize", [8, 8, 4]]
ONE = ["channelPatterns", ["CamA_ch0"]]
COLON = "MATLAB:colon:nonFiniteEndpoint"
NONANINF = "MATLAB:nonaninf"
PADSIZE = "MATLAB:padarray:expectedInteger"


def volume(name, shape=IM):
    return CounterRNG(0, f"x07/{name}").integers(0, 60000, shape, dtype=np.uint16)


def tif_file(name):
    def make(path):
        synth.write_tiff(path, volume(name))

    return make


def zarr_tree(name):
    def make(folder):
        synth.write_zarr(
            Path(folder) / "im.zarr", volume(name), block_size=(4, 16, 16), layout="matlab"
        )

    return make


def tif_dataset(name, nt=1, fmt="tif"):
    """nt time points of one channel (tiff or zarr files in MATLAB's layout)"""

    def make(folder):
        for t in range(nt):
            p = Path(folder) / f"scan_CamA_ch0_t{t:04d}.{fmt}"
            if fmt == "zarr":
                synth.write_zarr(p, volume(f"{name}/{t}"), block_size=(4, 16, 16), layout="matlab")
            else:
                synth.write_tiff(p, volume(f"{name}/{t}"))

    return make


def call_matlab(name, *args):
    from petakit5d import matlab

    getattr(matlab, name)(*args)
    return []


def _in_zarr(name):
    def f(d, *args):
        return call_matlab(name, d + "im.zarr", *args)

    f.__name__ = f"py_{name}_zarr"
    return f


def _named(name):
    def f(*args):
        return call_matlab(name, *args)

    f.__name__ = f"py_{name}"
    return f


def file_case(cid, name, out, *args, error=False):
    """``name(<im.tif>, <out>, args...)`` (out None: no output path)"""
    a = [InPath("im")] + ([Work(out)] if out else []) + list(args)
    return Case(
        cid,
        matlab=name,
        nargout=0,
        args=a,
        python=_named(name),
        inputs={"im": File(tif_file(cid), "im.tif")},
        expect_error=error,
        files=None if not error else False,
        error_message=bool(error),
    )


def zarr_case(cid, name, out, *args, error=False):
    """``name(<d>/im.zarr, <out>, args...)``"""
    a = [InPath("d")] + ([Work(out)] if out else []) + list(args)
    return Case(
        cid,
        code=f"{name}([in{{1}}, 'im.zarr'], in{{2 : end}}); out = {{}};",
        nargout=0,
        args=a,
        python=_in_zarr(name),
        inputs={"d": Tree(zarr_tree(cid), "d")},
        expect_error=error,
        files=None if not error else False,
        error_message=bool(error),
    )


def ds_case(cid, name, make, *args, error=False):
    """``name(<d>, args...)``: a dataset function on a folder"""
    return Case(
        cid,
        code=f"{name}(in{{:}}); out = {{}};",
        nargout=0,
        args=[InPath("d"), *args],
        python=_named(name),
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        files=None if not error else False,
        error_message=bool(error),
    )


def _crop_cases():
    c = "XR_crop_frame"
    d = "XR_crop_dataset"
    large = ["zarrFile", True, "largeFile", True, *SMALL]
    return [
        # XR_crop_frame: min(bbox(4 : 6), imSize) cuts both ends at the image's end
        file_case("crop_frame.tif.nan_end", c, "c.tif", NANB, *LOCAL),
        file_case("crop_frame.tif.nan_z_end", c, "c.tif", NANZ, *LOCAL),
        file_case("crop_frame.tif.inf_end", c, "c.tif", INFB, *LOCAL),
        # with pad: NaN > imSize is false (no padding), Inf pads by Inf
        file_case("crop_frame.tif.nan_end_pad", c, "c.tif", NANB, "pad", True, *LOCAL),
        file_case(
            "crop_frame.tif.inf_end_pad", c, "c.tif", INFB, "pad", True, *LOCAL, error=PADSIZE
        ),
        zarr_case("crop_frame.zarr.nan_end", c, "c.tif", NANB, "zarrFile", True, *LOCAL),
        zarr_case(
            "crop_frame.zarr.nan_end_zarr_out", c, "c.zarr", NANZ, "zarrFile", True,
            "saveZarr", True, *LOCAL,
        ),
        # large files (XR_crop_zarr): the batches of a size NaN; Inf is out of bound
        zarr_case("crop_frame.large.nan_end", c, "c.zarr", NANB, *large, *LOCAL, error=NONANINF),
        zarr_case(
            "crop_frame.large.nan_end_pad", c, "c.zarr", NANB, "pad", True, *large, *LOCAL,
            error=NONANINF,
        ),
        zarr_case("crop_frame.large.inf_end", c, "c.zarr", INFB, *large, *LOCAL, error=True),
        zarr_case(
            "crop_frame.large.inf_end_pad", c, "c.zarr", INFB, "pad", True, *large, *LOCAL,
            error=True,
        ),
        # XR_crop_dataset
        ds_case("crop_dataset.fixed.nan_end", d, tif_dataset("f"), NANB, *ONE, *LOCAL),
        ds_case(
            "crop_dataset.fixed.nan_end_pad", d, tif_dataset("fp"), NANZ, "pad", True, *ONE,
            *LOCAL,
        ),
        ds_case(
            "crop_dataset.zarr.nan_end", d, tif_dataset("z", fmt="zarr"), NANB, "zarrFile",
            True, *ONE, *LOCAL,
        ),
        # moving: the check of the last time point passes with NaN, fails with Inf
        ds_case(
            "crop_dataset.moving.nan_end", d, tif_dataset("m", 3), NANB, "cropType", "moving",
            "lastStartCoords", [6, 5, 3], *ONE, *LOCAL,
        ),
        ds_case(
            "crop_dataset.moving.inf_end", d, tif_dataset("mi", 3), INFB, "cropType",
            "moving", "lastStartCoords", [6, 5, 3], *ONE, *LOCAL, error=True,
        ),
        # center: the start is NaN (or -Inf)
        ds_case(
            "crop_dataset.center.nan_end", d, tif_dataset("c"), NANB, "cropType", "center",
            *ONE, *LOCAL, error=COLON,
        ),
        ds_case(
            "crop_dataset.center.nan_end_pad", d, tif_dataset("cp"), NANB, "cropType",
            "center", "pad", True, *ONE, *LOCAL, error=COLON,
        ),
        ds_case(
            "crop_dataset.center.inf_end", d, tif_dataset("ci"), INFB, "cropType", "center",
            *ONE, *LOCAL, error=COLON,
        ),
        ds_case(
            "crop_dataset.center.inf_end_pad", d, tif_dataset("cip"), INFB, "cropType",
            "center", "pad", True, *ONE, *LOCAL, error=PADSIZE,
        ),
        ds_case(
            "crop_dataset.large.nan_end", d, tif_dataset("l", fmt="zarr"), NANB, *large, *ONE,
            *LOCAL, error=NONANINF,
        ),
    ]  # fmt: skip


def _error_cases():
    return [
        file_case(
            "deskew_rotate_frame.nan_end", "XR_deskewRotateFrame", None, 0.108, 0.5,
            "inputBbox", NANB, "Save16bit", True, error=COLON,
        ),
        file_case(
            "rotate_frame_3d.nan_end", "XR_RotateFrame3D", None, 0.108, 0.5, "bbox", NANB,
            error=COLON,
        ),
        file_case(
            "resample_frame.nan_end", "XR_resampleFrame", "r.tif", [2, 2, 2], "inputBbox",
            NANB, error=COLON,
        ),
        zarr_case(
            "resample_single_zarr.nan_end", "XR_resampleSingleZarr", "r.zarr", [2, 2, 2],
            "inputBbox", NANB, "parseCluster", False, error=NONANINF,
        ),
        file_case("save_mip_tiff.nan_end", "saveMIP_tiff", "m.tif", "inputBbox", NANB, error=COLON),
        zarr_case(
            "mip_zarr.nan_end", "XR_MIP_zarr", None, "inputBbox", NANB, "parseCluster", False,
            "batchSize", [16, 16, 4], error=NONANINF,
        ),
        file_case(
            "tiff_to_zarr.nan_end", "tiffToZarr", "t.zarr", [], "inputBbox", NANB, error=COLON
        ),
        file_case(
            "tiff_to_zarr.nan_tile_out_end", "tiffToZarr", "t.zarr", [], "tileOutBbox", NANB,
            error=COLON,
        ),
        zarr_case(
            "resave_single_zarr.nan_end", "XR_resaveSingleZarr", "r.zarr", "inputBbox", NANB,
            "parseCluster", False, error=NONANINF,
        ),
        ds_case(
            "deskew_rotate_wrapper.nan_end", "XR_deskew_rotate_data_wrapper", tif_dataset("w1"),
            "inputBbox", NANB, *ONE, "Deskew", True, "Rotate", True, *LOCAL, error=COLON,
        ),
        ds_case(
            "mip_wrapper.nan_end", "XR_MIP_wrapper", tif_dataset("w2"), "inputBbox", NANB,
            *ONE, *LOCAL, error=COLON,
        ),
        ds_case(
            "tiff_to_zarr_wrapper.nan_end", "XR_tiffToZarr_wrapper", tif_dataset("w3"),
            "inputBbox", NANB, *ONE, *LOCAL, error=COLON,
        ),
        ds_case(
            "resave_zarr_wrapper.nan_end", "XR_resave_zarr_wrapper",
            tif_dataset("w4", fmt="zarr"), "inputBbox", NANB, *ONE, *LOCAL, error=NONANINF,
        ),
        ds_case(
            "resample_dataset.nan_end", "XR_resample_dataset", tif_dataset("w5"), [2, 2, 2],
            "inputBbox", NANB, *ONE, *LOCAL, error=COLON,
        ),
    ]  # fmt: skip


def py_imresize3_fast(v, out_size, method, out_bbox):
    """(the MATLAB-name layer reads 'linear' as a parameter name: imresize3_fast's positional
    arguments are given in Python's form)"""
    from petakit5d.core.arrays import from_matlab
    from petakit5d.core.axes import from_yxz
    from petakit5d.matlabcompat.resize import imresize3_fast

    box = B.box_from_matlab(out_bbox)
    return [imresize3_fast(from_matlab(v), tuple(from_yxz(out_size)), method, box)]


def vol_case(cid, name, *args, error, python=None):
    """``out = name(<a (y, x, z) single volume>, args...)``"""
    v = np.ascontiguousarray(np.transpose(volume(cid).astype(np.float32), (1, 2, 0)))
    return Case(
        cid,
        matlab=name,
        nargout=1,
        args=[v, *args],
        python=python or _named(name),
        expect_error=error,
        error_message=True,
    )


def _more_error_cases():
    dsr = ["Deskew", True, "Rotate", True]
    return [
        zarr_case(
            "deskew_rotate_frame.zarr.nan_end", "XR_deskewRotateFrame", None, 0.108, 0.5,
            "inputBbox", NANB, "Save16bit", True, error=COLON,
        ),
        zarr_case(
            "deskew_rotate_zarr.nan_end", "XR_deskewRotateZarr", None, 0.108, 0.5,
            "inputBbox", NANB, "parseCluster", False, error=NONANINF,
        ),
        ds_case(
            "deskew_rotate_wrapper.zarr.nan_end", "XR_deskew_rotate_data_wrapper",
            tif_dataset("w6", fmt="zarr"), "inputBbox", NANB, *ONE, *dsr, "zarrFile", True,
            *LOCAL, error=COLON,
        ),
        ds_case(
            "deskew_rotate_wrapper.large.nan_end", "XR_deskew_rotate_data_wrapper",
            tif_dataset("w7", fmt="zarr"), "inputBbox", NANB, *ONE, *dsr, "zarrFile", True,
            "largeFile", True, *LOCAL, error=NONANINF,
        ),
        vol_case(
            "imresize3_fast.nan_end", "imresize3_fast", [80, 100, 24], "linear", NANB,
            error=COLON, python=py_imresize3_fast,
        ),
        vol_case(
            "deskew_rotate_frame_3d.nan_end", "deskewRotateFrame3D", 32.45, 0.5, 0.108, "bbox",
            NANB, error="MATLAB:imref3d:expectedInteger",
        ),
    ]  # fmt: skip


CASES = _crop_cases() + _error_cases() + _more_error_cases()


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    rep = oracle.check(case, tmp_path)
    # Python's error carries MATLAB's identifier (the harness only notes a missing one)
    assert not [n for n in rep.notes if "no identifier attribute" in n], rep.notes


# ---------------------------------------------------------------------------------------------
# the conversions and the helpers
# ---------------------------------------------------------------------------------------------


def test_nan_end_survives_the_conversions():
    b = B.box_from_matlab([2, 3, 1, math.nan, 30, 8])
    assert b[:4] == (0, 1, 2, 8) and math.isnan(b[4]) and b[5] == 30
    m = B.box_to_matlab(b)
    assert m[:3] == (2, 3, 1) and m[4:] == (30, 8)
    # the NaN's bits are kept both ways (MATLAB's NaN and Python's)
    for nan in (NAN, math.nan):
        m = B.box_to_matlab(B.box_from_matlab([2, 3, 1, nan, 30, 8]))
        assert np.float64(m[3]).view(np.uint64) == np.float64(nan).view(np.uint64)
    assert B.has_nan_end(b) and not B.has_nan_end(B.box_from_matlab([1, 1, 1, 5, 5, math.inf]))
    assert B.open_nan_ends(b) == (0, 1, 2, 8, None, 30)
    # 2D
    assert math.isnan(B.box_from_matlab([1, 2, math.nan, 4])[2])


def test_nan_start_is_refused():
    with pytest.raises(ValueError, match="box starts"):
        B.box_from_matlab([math.nan, 1, 1, 5, 5, 5])


@pytest.mark.parametrize(
    "f",
    [
        lambda b: B.check_box(b),
        lambda b: B.box_shape(b, (4, 5, 6)),
        lambda b: B.box_slices(b),
        lambda b: B.clamp_box(b, (4, 5, 6)),
        lambda b: B.pad_widths(b, (4, 5, 6)),
        lambda b: B.crop(np.zeros((4, 5, 6)), b),
    ],
)
def test_helpers_refuse_a_nan_end_as_matlab_indexing(f):
    with pytest.raises(B.BoxEndError) as e:
        f((0, 0, 0, 2, math.nan, 3))
    assert e.value.identifier == "MATLAB:colon:nonFiniteEndpoint"
    assert isinstance(e.value, ValueError)


def test_registry_kind_keeps_a_nan_end():
    from petakit5d.params.kinds import KINDS

    k = KINDS["bbox"]
    b = k.from_matlab([1, 1, 1, math.inf, math.inf, math.nan])  # zmax NaN
    assert math.isnan(b[3]) and b[4] is None and b[5] is None
    back = k.to_matlab(b)
    assert back[3] == math.inf and back[4] == math.inf and math.isnan(back[5])
