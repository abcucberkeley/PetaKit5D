"""ND2 and CZI conversion without Bio-Formats (petakit5d.io.bioformats, task B09).

MATLAB has no tests of XR_bioformats_to_tiff_or_zarr_wrapper. The oracle cases run MATLAB's
wrapper (bfopen, as a black box) and the port's on public sample files (Bio-Formats' sample
images, CC BY 4.0), named by ``PETAKIT5D_FORMATS_SAMPLES`` (a folder with ``positions/`` (an ND2
of 5 positions and 2 channels), ``timeseries/`` (an ND2 of 50 time points and 2 channels),
``scenes/`` (a CZI of 3 channels and 21 z slices) and ``single/`` (a legacy ND2, BF007.nd2,
where the readers differ: bfopen gives 3 planes, nd2 one plane of other values; not a case,
see test_sample_files)); they skip without it, and the Python side needs the formats extra.
Then the plane order, the channel dealing, the names, tiff and zarr results, existing results,
and the errors, on stand-in readers.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from petakit5d.io import bioformats as bf
from petakit5d.io import read_tiff, read_zarr
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, Data, InPath

SAMPLES = "PETAKIT5D_FORMATS_SAMPLES"
try:
    import czifile  # noqa: F401
    import nd2  # noqa: F401

    HAS_READERS = True
except ImportError:
    HAS_READERS = False


def py_wrapper(*args):
    from petakit5d import matlab

    matlab.XR_bioformats_to_tiff_or_zarr_wrapper(*args)
    return []


def bf_case(cid, sub, *nv, error=False):
    return Case(
        cid,
        code="XR_bioformats_to_tiff_or_zarr_wrapper(in{:}); out = {};",
        nargout=0,
        args=[InPath("d"), *nv],
        python=py_wrapper,
        inputs={"d": Data(f"${SAMPLES}/{sub}", "d")},
        expect_error=error,
        py_error=(bf.BioformatsError,),
    )


CASES = [
    bf_case("bf_nd2_positions_tiff", "positions", "nChannels", 2),
    bf_case(
        "bf_nd2_timeseries_zarr",
        "timeseries",
        "nChannels",
        2,
        "saveZarr",
        True,
        "blockSize",
        [256, 256, 16],
    ),
    bf_case("bf_czi_scene_tiff", "scenes", "nChannels", 3),
    bf_case(
        "bf_nd2_three_channels_error",
        "positions",
        "nChannels",
        3,
        error="MATLAB:getReshapeDims:notDivisible",
    ),
]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    if not os.environ.get(SAMPLES):
        pytest.skip(f"{SAMPLES} (the sample files) is not set")
    if not HAS_READERS:
        pytest.skip("nd2 and czifile (the formats extra) are missing")
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# the plane order and the channels
# ---------------------------------------------------------------------------------------------


def test_planes_are_in_xyczt_order():
    t, z, c, y, x = 2, 3, 4, 2, 5
    a = np.arange(t * z * c * y * x).reshape(t, z, c, y, x)
    p = bf._planes(a, "TZCYX", "f")
    assert p.shape == (t * z * c, y, x)
    np.testing.assert_array_equal(p[1], a[0, 0, 1])
    np.testing.assert_array_equal(p[c], a[0, 1, 0])
    np.testing.assert_array_equal(p[z * c], a[1, 0, 0])
    # another axis order and singleton axes of other kinds
    b = np.ascontiguousarray(a.transpose(2, 0, 1, 3, 4))[np.newaxis, ..., np.newaxis]
    np.testing.assert_array_equal(bf._planes(b, "BCTZYXS", "f"), p)
    # missing axes count as one
    np.testing.assert_array_equal(bf._planes(a[0, 0], "CYX", "f"), a[0, 0])


def test_other_axes_raise():
    with pytest.raises(bf.BioformatsError, match="axes S") as e:
        bf._planes(np.zeros((2, 4, 5, 3)), "CYXS", "f")
    assert e.value.identifier == "petakit5d:bioformats:unsupportedAxes"
    with pytest.raises(bf.BioformatsError, match="unexpected axes"):
        bf._planes(np.zeros((2, 4, 5)), "CZY", "f")


def test_planes_are_dealt_to_channels_as_matlabs_reshape():
    p = np.arange(6 * 2 * 3).reshape(6, 2, 3)
    a, b = bf._channels(p, 2)
    np.testing.assert_array_equal(a, p[0::2])
    np.testing.assert_array_equal(b, p[1::2])
    (one,) = bf._channels(p[:1], 1)
    assert one.shape == (2, 3)  # squeeze: one plane is an image
    with pytest.raises(bf.BioformatsError) as e:
        bf._channels(p[:5], 2)
    assert e.value.identifier == "MATLAB:getReshapeDims:notDivisible"


def test_unsupported_formats_name_the_supported_ones(tmp_path):
    with pytest.raises(bf.BioformatsError, match=r"only \.nd2, \.czi files") as e:
        bf.bioformats_to_tiff_or_zarr_dataset(str(tmp_path), data_formats=(".nd2", ".lif"))
    assert "Bio-Formats is not part of the port" in str(e.value)
    with pytest.raises(bf.BioformatsError, match=r"only \.nd2, \.czi"):
        bf.read_series(str(tmp_path / "a.oib"))


# ---------------------------------------------------------------------------------------------
# the wrapper, on stand-in readers
# ---------------------------------------------------------------------------------------------


def _vol(shape, k):
    return (np.arange(int(np.prod(shape)), dtype=np.uint16).reshape(shape) * (k + 1)) % 60000


@pytest.fixture
def fake_reader(monkeypatch):
    series = {"one": [_vol((8, 6, 7), 0)], "two": [_vol((4, 6, 7), 1), _vol((4, 6, 7), 2)]}

    def read(path):
        name = os.path.basename(path).split(".")[0]
        return iter(series[name])

    monkeypatch.setattr(bf, "read_series", read)
    return series


def test_wrapper_names_and_tiffs(tmp_path, fake_reader):
    d = tmp_path / "data"
    d.mkdir()
    for n in ("one.nd2", "two.czi", "other.tif"):
        (d / n).write_bytes(b"")
    bf.bioformats_to_tiff_or_zarr_dataset(str(d))
    out = d / "tiffs"
    names = sorted(n for n in os.listdir(out))
    assert names == [
        "one_ch0.tif",
        "one_ch1.tif",
        "two_fov0001_ch0.tif",
        "two_fov0001_ch1.tif",
        "two_fov0002_ch0.tif",
        "two_fov0002_ch1.tif",
    ]
    np.testing.assert_array_equal(read_tiff(str(out / "one_ch1.tif")), fake_reader["one"][0][1::2])
    np.testing.assert_array_equal(
        read_tiff(str(out / "two_fov0002_ch0.tif")), fake_reader["two"][1][0::2]
    )


def test_wrapper_zarr_existing_results_and_overwrite(tmp_path, fake_reader):
    d = tmp_path / "data"
    d.mkdir()
    (d / "one.nd2").write_bytes(b"")
    (d / "z").mkdir()
    bf.bioformats_to_tiff_or_zarr_dataset(
        str(d), result_dir_name="z", save_zarr=True, block_size=(2, 4, 4)
    )
    from petakit5d.io import get_zarr_info

    fn = str(d / "z" / "one_ch0.zarr")
    assert get_zarr_info(fn).block_size == (2, 4, 4)
    np.testing.assert_array_equal(read_zarr(fn), fake_reader["one"][0][0::2])
    # an existing result is kept, overWrite replaces it
    fake_reader["one"][0] = fake_reader["one"][0] + 1
    bf.bioformats_to_tiff_or_zarr_dataset(str(d), result_dir_name="z", save_zarr=True)
    np.testing.assert_array_equal(read_zarr(fn), (fake_reader["one"][0] - 1)[0::2])
    bf.bioformats_to_tiff_or_zarr_dataset(
        str(d), result_dir_name="z", save_zarr=True, over_write=True
    )
    np.testing.assert_array_equal(read_zarr(fn), fake_reader["one"][0][0::2])
    assert sorted(os.listdir(d / "z")) == ["one_ch0.zarr", "one_ch1.zarr"]


def test_matlab_name_converts_the_block_size(tmp_path, fake_reader):
    from petakit5d import matlab
    from petakit5d.io import get_zarr_info

    d = tmp_path / "data"
    d.mkdir()
    (d / "one.nd2").write_bytes(b"")
    matlab.XR_bioformats_to_tiff_or_zarr_wrapper(
        str(d), "saveZarr", True, "blockSize", [4, 5, 2], "resultDirName", "r"
    )
    assert get_zarr_info(str(d / "r" / "one_ch0.zarr")).block_size == (2, 4, 5)


@pytest.mark.skipif(not HAS_READERS, reason="nd2 and czifile (the formats extra) are missing")
def test_sample_files():
    root = os.environ.get(SAMPLES)
    if not root:
        pytest.skip(f"{SAMPLES} (the sample files) is not set")
    (pos,) = os.listdir(f"{root}/positions")
    series = list(bf.read_series(f"{root}/positions/{pos}"))
    assert [s.shape for s in series] == [(2, 1040, 1392)] * 5
    (czi,) = os.listdir(f"{root}/scenes")
    (scene,) = bf.read_series(f"{root}/scenes/{czi}")
    assert scene.shape == (63, 512, 672)
    # a legacy ND2 where nd2 and Bio-Formats differ: nd2 reads one plane (bfopen three)
    (single,) = bf.read_series(f"{root}/single/BF007.nd2")
    assert single.shape == (1, 156, 164)
