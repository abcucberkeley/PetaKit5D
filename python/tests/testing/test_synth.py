"""petakit5d.testing.synth: bit-stable generators of synthetic data.

test_generator_hashes pins the bits of every generator: the same hashes on Linux, macOS and
Windows, with every NumPy version (CI runs it on all three systems). If a generator changes
on purpose, every oracle case with its data changes too: regenerate the hashes with
`python tests/testing/test_synth.py` and say so in the commit.
"""

from __future__ import annotations

import csv
import math
import sys

import numpy as np
import pytest

from petakit5d.testing import synth
from petakit5d.testing.fixtures import limit_threads, pk_temp_folder, rng, thread_count  # noqa: F401
from petakit5d.testing.synth import CounterRNG


def _generators() -> dict:
    """name -> function making the array whose hash is pinned"""
    g = {
        "uint32": lambda: CounterRNG(0).uint32((5, 7, 3)),
        "uint32_seed_stream": lambda: CounterRNG(2**63 + 12345, "decon").uint32(1001),
        "uniform64": lambda: CounterRNG(1).uniform((9, 8, 7)),
        "uniform32": lambda: CounterRNG(1).uniform((9, 8, 7), np.float32),
        "normal64": lambda: CounterRNG(3).normal((33, 17, 5)),
        "normal32": lambda: CounterRNG(3).normal(1000, np.float32),
        "integers_uint16": lambda: CounterRNG(4).integers(1, 65535, (31, 9), np.uint16),
        "integers_signed": lambda: CounterRNG(4).integers(-1000, 1000, 257, np.int32),
        # inputs made with exact operations only (linspace and geomspace aren't)
        "exp": lambda: synth.stable_exp(-700.0 + np.arange(4097) * (1409.0 / 4096.0)),
        "log": lambda: synth.stable_log(
            np.ldexp(
                1.0 + (np.arange(4097) % 97) / 97.0, (np.arange(4097) // 4 - 512).astype(np.int32)
            )
        ),
        "sincos": lambda: np.stack(synth.stable_sincos_2pi(np.arange(4096) / 4096.0)),
        "beads": lambda: synth.beads(
            (32, 48, 40), [[4.0, 12, 8], [10.6, 24, 20], [30, 36, 32]], [1.5, 2], 1000, 100
        ),
        "beads_2d_amplitudes": lambda: synth.beads(
            (40, 30), [[10.3, 7.2, 1], [30, 20, 1]], 2.0, [500, 2000]
        ),
        "texture3d": lambda: synth.texture((32, 48, 40), 3, (0, 200)),
        "texture2d": lambda: synth.texture((64, 50), 2, (100, 2000), CounterRNG(7)),
        "texture_nz1": lambda: synth.texture((20, 30, 1), 1.5),
        "decon_psf": lambda: synth.decon_psf(),
        "decon_psf_even": lambda: synth.decon_psf((32, 32, 24), (1.2, 1.8), 1),
        "gaussian_filter": lambda: synth.gaussian_filter(CounterRNG(9).normal((17, 13, 11)), 1.3),
    }
    for dt in ("logical", "uint8", "uint16", "int16", "uint32", "single", "double"):
        g[f"random_image_{dt}"] = lambda dt=dt: synth.random_image((23, 19, 11), dt)
    return g


GENERATORS = _generators()

# made on Linux x86-64 (NumPy 2.5.4); the same on every system
EXPECTED = {
    "uint32": "3084a6b8a6ae3830b745ceff11953d720094a5be8a5096ebc1bb89248c89388c",
    "uint32_seed_stream": "5fc12cb9582c3fa8040f19b67d5902463e43b1e20d88606b8bb1d7101d7915ee",
    "uniform64": "b0bf0b5390991e5c694429fd4da0628ececa7b9f93323e3a34c51467e4d31949",
    "uniform32": "21f90f2134f332694d02f8eafe3389c03e564092c5aef88af25bbb67770ecb57",
    "normal64": "1eb4169ea8c757f9389fb53b8af00f22791b61b662f9bf981b37766bbb64637b",
    "normal32": "e850d97daed2e820cf300ed3601d4f507ee6c41999c5be036e251504bb431442",
    "integers_uint16": "1067dc88609da25934eca8312c58ac54f62285a3d8d6e5e9fad47122b4c28681",
    "integers_signed": "3d039c22efdc5cf78d36e8cfbcb91ad491a1a0d7c0a9d510b62f954fff96da25",
    "exp": "fef02ca0c21168f73a55046024db5e8c73667aeab9cc7c75ec3a780e02ad2603",
    "log": "5d6305ecb208c007a1404631747f012bdb31646ecafe27cf9de4eefd4052ba4c",
    "sincos": "c3258c8f0ffe3e0925514f99ff11d86bf54896962ffa7698ebc4b7ca2a5d60de",
    "beads": "224f38cdb4d1aa0fe0a4012bc42b8254261baea2340ce7384204409ad007de12",
    "beads_2d_amplitudes": "f7c2a619627a7f1f59c2e7a627266427df1cb6a4d158bf5ba156eab1ed8c3412",
    "texture3d": "b2b3344a0a20282132d7ee0c40a59d2a307204030c04fc6dc737f5a8201c9c2c",
    "texture2d": "accb2bf474f9f579727638f3f254116372ec453832033cb477da3bf4bdd33fa9",
    "texture_nz1": "c29e3fa918534aa897562a82198973356eccc8f77776ea9b8e4adc11a7ab12cf",
    "decon_psf": "567790005802651c75bc0d41291ba5732ffc92f73c22f3f1a8e33ff80108786f",
    "decon_psf_even": "fa4bdc44e0ce49ac1424f9508d0a32bf43b6b0f8761ab6c43bd4b5ebc9f52258",
    "gaussian_filter": "7e874a94edff45c27c27d1e0a1f35e5c995f170f7abb09100905c550c1020e1c",
    "random_image_logical": "68d1e2a19c08f5a655f69410f0d4e08f88c0d4711addcaa60dee0c7c2404788d",
    "random_image_uint8": "825318c0dbba854096008c3dc7b44c02f0c3c856f2a5be73a742d25c1b9bd275",
    "random_image_uint16": "d3b3669e7d72cb8ed7436803689513bfff07e309f0989d38a4cf4cd3a8eaef9f",
    "random_image_int16": "890f542b6d01c224d803b6efcf65986af2bd2969d03e46b51d1ce5299dd7993e",
    "random_image_uint32": "1e8d17a6728f4e718dd5f6222e298e1b9ac54328dcb1eccde68dacdf6eb561fa",
    "random_image_single": "1287e595f1c3c35a19cda6a6da897754c02e2ed9e1ea8b52a683ce50c8e4e859",
    "random_image_double": "b70f90ccf707c875418226a18f808c8e0451486b778015dd00e0b0b9e7e71865",
}


@pytest.mark.parametrize("name", list(GENERATORS))
def test_generator_hashes(name):
    a = GENERATORS[name]()
    assert synth.array_hash(a) == EXPECTED[name], (
        f"{name} changed (dtype {a.dtype}, shape {a.shape})"
    )


def test_philox_known_answers():
    # Random123's known-answer vectors of philox4x32-10
    cases = [
        ([0, 0, 0, 0], (0, 0), [0x6627E8D5, 0xE169C58D, 0xBC57AC4C, 0x9B00DBD8]),
        (
            [0xFFFFFFFF] * 4,
            (0xFFFFFFFF, 0xFFFFFFFF),
            [0x408F276D, 0x41C83B0E, 0xA20BC7C6, 0x6D5451FD],
        ),
        (
            [0x243F6A88, 0x85A308D3, 0x13198A2E, 0x03707344],
            (0xA4093822, 0x299F31D0),
            [0xD16CFE09, 0x94FDCCEB, 0x5001E420, 0x24126EA1],
        ),
    ]
    for ctr, key, out in cases:
        got = synth.philox4x32(np.array([ctr], dtype=np.uint64), key)[0]
        assert [int(v) for v in got] == out


def test_values_depend_only_on_their_index():
    # one call, or the same blocks in pieces, or from a later position
    whole = CounterRNG(5, "s").uint32(4 * 100)
    g = CounterRNG(5, "s")
    parts = np.concatenate([g.uint32(4 * 30), g.uint32(4 * 70)])
    assert np.array_equal(whole, parts)
    assert np.array_equal(CounterRNG(5, "s").at(30).uint32(4 * 70), whole[120:])
    # column-major fill: element (i, j, k) is value i + ny (j + nx k)
    v = CounterRNG(5, "s").uint32((4, 5, 20))
    assert np.array_equal(v.reshape(-1, order="F"), whole)
    assert v.flags.f_contiguous


def test_streams_and_seeds_differ():
    a = CounterRNG(0).uint32(64)
    assert not np.array_equal(a, CounterRNG(1).uint32(64))
    assert not np.array_equal(a, CounterRNG(0, 1).uint32(64))
    assert not np.array_equal(a, CounterRNG(0, "a").uint32(64))
    assert not np.array_equal(
        CounterRNG(0).spawn("x").uint32(64), CounterRNG(0).spawn("y").uint32(64)
    )


def test_rng_fixture_is_seed_zero(rng):  # noqa: F811
    assert np.array_equal(rng.uint32(8), CounterRNG(0).uint32(8))


def test_distributions():
    u = CounterRNG(11).uniform(200_000)
    assert u.min() >= 0 and u.max() < 1
    assert abs(u.mean() - 0.5) < 0.005
    u32 = CounterRNG(11).uniform(200_000, np.float32)
    assert u32.dtype == np.float32 and u32.min() >= 0 and u32.max() < 1
    z = CounterRNG(12).normal(200_000)
    assert abs(z.mean()) < 0.01 and abs(z.std() - 1) < 0.01
    assert abs(np.mean(z**4) - 3) < 0.1
    i = CounterRNG(13).integers(1, 6, 60_000)
    assert i.min() == 1 and i.max() == 6
    assert np.all(np.abs(np.bincount(i)[1:] / 60_000 - 1 / 6) < 0.01)
    assert CounterRNG(0).integers(0, 2**32 - 1, 10, np.uint32).dtype == np.uint32
    with pytest.raises(ValueError):
        CounterRNG(0).integers(0, 2**32, 1)
    with pytest.raises(ValueError):
        CounterRNG(-1)


def test_elementary_functions_are_accurate():
    x = np.linspace(-690.0, 709.0, 100_001)
    assert np.max(np.abs(synth.stable_exp(x) - np.exp(x)) / np.spacing(np.exp(x))) <= 2
    assert synth.stable_exp(np.array([-1000.0]))[0] == 0  # no subnormals
    assert synth.stable_exp(np.array([1000.0]))[0] == np.inf
    assert synth.stable_exp(np.array([0.0]))[0] == 1
    y = np.geomspace(1e-300, 1e300, 100_001)
    assert (
        np.max(np.abs(synth.stable_log(y) - np.log(y)) / np.maximum(1, np.abs(np.log(y)))) < 4e-16
    )
    assert synth.stable_log(np.array([1.0]))[0] == 0
    u = np.arange(1 << 14) / (1 << 14)
    s, c = synth.stable_sincos_2pi(u)
    assert np.max(np.abs(s - np.sin(2 * np.pi * u))) < 2e-15
    assert np.max(np.abs(c - np.cos(2 * np.pi * u))) < 2e-15
    with pytest.raises(ValueError):
        synth.stable_log(np.array([0.0]))


def test_beads_match_the_formula():
    pos = [[10.3, 12.0, 7.5], [20, 5, 3]]
    vol = synth.beads((24, 18, 12), pos, [1.5, 2.0], [1000, 300], 100)
    y, x, z = np.meshgrid(np.arange(1, 25), np.arange(1, 19), np.arange(1, 13), indexing="ij")
    ref = np.full(vol.shape, 100, np.float32)
    for (py, px, pz), a in zip(pos, [1000, 300], strict=True):
        ref = ref + (
            a * np.exp(-((y - py) ** 2 + (x - px) ** 2) / 4.5 - (z - pz) ** 2 / 8.0)
        ).astype(np.float32)
    assert vol.dtype == np.float32 and vol.flags.f_contiguous and vol.shape == (24, 18, 12)
    assert np.max(np.abs(vol - ref)) <= 4 * np.spacing(np.float32(1100))
    assert (
        synth.beads((5, 6, 7), np.zeros((0, 3)), 1.0, 1, 3).tolist()
        == np.full((5, 6, 7), 3, np.float32).tolist()
    )


def test_centroids_find_the_beads():
    pos = np.array([[10.3, 12.6, 9.2], [22, 30, 14]])
    vol = synth.beads((32, 40, 24), pos, [1.2, 1.6], 1000, 50)
    c = synth.centroids(vol, pos + 0.4, 5)
    assert np.max(np.abs(c - pos)) < 0.02


def test_texture_spans_its_range():
    t = synth.texture((16, 20, 12), 2, (100, 2000))
    assert t.dtype == np.float32 and t.flags.f_contiguous
    assert t.min() == 100 and t.max() == 2000
    assert np.array_equal(t, synth.texture((16, 20, 12), 2, (100, 2000), CounterRNG(0)))
    assert not np.array_equal(t, synth.texture((16, 20, 12), 2, (100, 2000), CounterRNG(1)))
    # smooth: neighbors are closer than independent values
    assert np.std(np.diff(t, axis=0)) < 0.5 * np.std(t)
    t2 = synth.texture((30, 25), 2)
    assert t2.shape == (30, 25) and t2.min() == 100 and t2.max() == 2000


def test_gaussian_filter_keeps_constants_and_sums():
    c = np.full((9, 7, 5), 3.0)
    assert np.allclose(synth.gaussian_filter(c, 1.7), 3.0, rtol=1e-15)
    d = np.zeros((41, 1, 1))
    d[20] = 1
    f = synth.gaussian_filter(d, 2.0)
    assert math.isclose(math.fsum(f.ravel().tolist()), 1.0, rel_tol=1e-14)
    assert np.count_nonzero(f) == 2 * math.ceil(2 * 2.0) + 1  # the kernel size of imgaussfilt


def test_decon_psf_peak_and_zero_slices():
    p = synth.decon_psf()
    assert p.shape == (31, 31, 31) and p.dtype == np.float32
    assert np.unravel_index(np.argmax(p), p.shape) == (15, 15, 15)
    assert p.max() == 10000
    assert not p[:, :, :5].any() and not p[:, :, -5:].any() and p[:, :, 5].any()
    e = synth.decon_psf((32, 32, 24))
    assert np.unravel_index(np.argmax(e), e.shape) == (16, 16, 12)


def test_random_image_classes():
    for name, dt in (
        ("logical", np.bool_),
        ("uint8", np.uint8),
        ("uint16", np.uint16),
        ("int16", np.int16),
        ("single", np.float32),
        ("double", np.float64),
    ):
        im = synth.random_image((6, 5, 4), name)
        assert im.dtype == dt and im.shape == (6, 5, 4) and im.flags.f_contiguous
        if dt not in (np.bool_, np.float32, np.float64):
            assert im.min() >= 1
    assert np.array_equal(
        synth.random_image((6, 5), np.uint16), synth.random_image((6, 5), "uint16")
    )
    with pytest.raises(TypeError):
        synth.random_image((2, 2), "int64")


def test_results_do_not_depend_on_threads(thread_count):  # noqa: F811
    # NumPy's element-wise loops are single threaded; this pins that nothing else is used
    assert synth.array_hash(GENERATORS["texture3d"]()) == EXPECTED["texture3d"]


def test_temp_folder_style(pk_temp_folder):  # noqa: F811
    assert pk_temp_folder.endswith("/") and "\\" not in pk_temp_folder


def test_stitch_tiles_names_and_image_list(pk_temp_folder):  # noqa: F811
    pytest.importorskip("cpptiff")
    gt = synth.texture((24, 60, 10), 2, (100, 3000)).astype(np.uint16)
    # offsets [y0, x0, z0] and voxel sizes [y, x, z]; stage coordinates and indices are [x, y, z]
    tiles = synth.cut_tiles(gt, [[0, 0, 0], [4, 30, 2]], (20, 30, 8), (0.1, 0.108, 0.3))
    assert [t["idx"] for t in tiles] == [[0, 0, 0], [1, 1, 1]]
    assert tiles[1]["xyz"] == [30 * 0.108, 4 * 0.1, 2 * 0.3]
    assert np.array_equal(tiles[1]["im"], gt[4:24, 30:60, 2:10])
    names, list_path = synth.write_stitch_tiles(pk_temp_folder, tiles)
    assert names[1] == (
        "Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs_"
        "001x_001y_001z_0000t"
    )
    assert list_path == pk_temp_folder + "ImageList.csv"
    rows = list(csv.reader(open(list_path)))
    assert rows[0] == ["Filename", "StageX_um_", "StageY_um_", "StageZ_um_"]
    assert rows[2] == [names[1] + ".tif", "3.24", "0.4", "0.6"]
    # rows are appended; no camera, another iteration, a tile folder
    names2, _ = synth.write_stitch_tiles(pk_temp_folder, tiles[:1], cam="", iter=1, tile_dir="DSR")
    assert names2[0].startswith("Scan_Iter_0001_0000_ch0_CAM1_")
    assert len(list(csv.reader(open(list_path)))) == 4
    import cpptiff

    back = cpptiff.read_tiff(f"{pk_temp_folder}{names[1]}.tif")
    assert np.array_equal(np.transpose(back, (1, 2, 0)), tiles[1]["im"])


def test_zarr_tiles(pk_temp_folder):  # noqa: F811
    cppzarr = pytest.importorskip("cppzarr")
    vol = synth.random_image((40, 33, 9), "uint16")
    names, _ = synth.write_stitch_tiles(
        pk_temp_folder, [{"im": vol, "idx": [0, 0, 0], "xyz": [0, 0, 0]}], save_format="zarr"
    )
    back = cppzarr.read_zarr(f"{pk_temp_folder}{names[0]}.zarr")
    assert np.array_equal(back, vol)


if __name__ == "__main__":
    # print the hash table (run on Linux, paste into EXPECTED)
    for k, f in GENERATORS.items():
        print(f'    "{k}": "{synth.array_hash(f())}",')
    sys.exit(0)
