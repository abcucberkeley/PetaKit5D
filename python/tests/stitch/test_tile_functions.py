"""The stitching's edge erosion as MATLAB code (petakit5d.stitch.tile_functions, a stand-in of
task S02 made by P02): erodeVolumeBy2DProjection against MATLAB's, and the two forms of MATLAB
code the frame writes for it read by ``parse_tile_function``."""

from __future__ import annotations

import numpy as np
import pytest

from petakit5d.core.arrays import to_matlab
from petakit5d.core.dtypes import NAN
from petakit5d.stitch import tile_functions as TF
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value, Vol
from petakit5d.testing.synth import CounterRNG

SHAPE = (11, 13, 17)  # (z, y, x): not cubic


def volume(kind: str = "uint16", holes: bool = True, seed: int = 0) -> np.ndarray:
    """random positive values; with holes, columns along y of zeros near the edges (so the
    projection is not full) and a zero corner"""
    v = CounterRNG(seed, f"tf/{kind}").integers(1, 3000, SHAPE).astype(np.uint16)
    if holes:
        v[:, :, :2] = 0
        v[:3, :, 9] = 0
        v[8:, :, 12:] = 0
    if kind == "single":
        v = v.astype(np.float32)
        if holes:
            v[5, 4, 6] = NAN
            v[:, :, 15] = NAN
    return np.ascontiguousarray(v)


def _matlab(v):
    return Value(lambda: to_matlab(v))


def erode_case(cid, vol, esize):
    return Case(
        f"erode.{cid}",
        matlab="erodeVolumeBy2DProjection",
        args=[Vol(In("v")), esize],
        python=TF.erode_volume_by_2d_projection,
        out_roles=("vol",),
        inputs={"v": _matlab(vol)},
    )


def code_case(cid, code, vol):
    return Case(
        f"code.{cid}",
        code=f"f = str2func('{code}'); out = {{f(in{{1}})}};",
        args=[Vol(In("v"))],
        python=lambda v: TF.parse_tile_function(code)(v),
        out_roles=("vol",),
        inputs={"v": _matlab(vol)},
    )


CASES = [
    erode_case("holes_2", volume(), 2),
    erode_case("holes_3", volume(), 3),
    erode_case("full_2", volume(holes=False), 2),
    erode_case("zero", volume(), 0),
    erode_case("single_nan", volume("single"), 2),
    erode_case("image", np.ascontiguousarray(volume()[0]), 2),
    code_case("plain", "@(x)erodeVolumeBy2DProjection(x,2)", volume()),
    code_case("offset", "@(x)erodeVolumeBy2DProjection(x+1,2)", volume()),
    code_case("offset_single", "@(x)erodeVolumeBy2DProjection(x+1,3)", volume("single")),
]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_other_matlab_code_is_s02s():
    with pytest.raises(NotImplementedError, match="S02"):
        TF.parse_tile_function("@(x)processFFCorrectionFrame(x,'a.tif','b.tif',1,0.4,false)")


def test_the_offset_saturates_uint16():
    v = np.full(SHAPE, 65535, np.uint16)
    out = TF.parse_tile_function("@(x)erodeVolumeBy2DProjection(x+1,1)")(v)
    assert out.dtype == np.uint16 and out[5, 6, 7] == 65535 and out[0, 6, 7] == 0


@pytest.mark.parametrize("kind", ["uint16", "single"])
def test_threads_do_not_change_the_result(kind):
    v = np.ascontiguousarray(np.tile(volume(kind), (40, 20, 4)))  # big enough to be split
    one = TF.erode_volume_by_2d_projection(v, 2, num_threads=1)
    many = TF.erode_volume_by_2d_projection(v, 2, num_threads=7)
    assert one.tobytes() == many.tobytes()


def test_a_full_projection_erodes_its_border_and_the_y_rows():
    v = np.ones(SHAPE, np.uint16)
    out = TF.erode_volume_by_2d_projection(v, 2)
    want = np.zeros(SHAPE, np.uint16)
    want[2:-2, 2:-2, 2:-2] = 1
    np.testing.assert_array_equal(out, want)


def test_an_image_keeps_its_values():
    """(MATLAB's vol .* (vol ~= 0)) negative values stay"""
    im = np.array([[-1.5, 0, 2], [3, -0.0, -4]], np.float32)
    out = TF.erode_volume_by_2d_projection(im, 2)
    assert out.tobytes() == im.tobytes()
