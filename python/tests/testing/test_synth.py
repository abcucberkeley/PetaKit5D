"""petakit5d.testing.synth: bit-stable generators of synthetic data.

The generators return the port's arrays ((z, y, x), C order) with sizes, positions and sigmas
in (z, y, x), positions 0-based. Their MATLAB form (to_matlab) is, bit for bit, what they gave
before the convention change, so the hashes of stored inputs stay valid:

- test_generator_hashes pins the 26 hashes made before the change (of MATLAB's form);
- test_matlab_form_is_unchanged runs every generator on 20 parameter sets each, given in the
  old (MATLAB) form and converted, and compares with the hashes the old code gave
  (`python tests/testing/test_synth.py --old <old synth.py>` prints them).

The same hashes on Linux, macOS and Windows, with every NumPy version (CI runs these on all
three systems). If a generator changes on purpose, every oracle case with its data changes
too: say so in the commit.
"""

from __future__ import annotations

import csv
import math
import sys

import numpy as np
import pytest

from petakit5d.core.arrays import to_matlab
from petakit5d.testing import synth
from petakit5d.testing.fixtures import (  # noqa: F401
    limit_threads,
    paired_input,
    pk_temp_folder,
    rng,
    thread_count,
    write_volume,
    zarr_layout,
)
from petakit5d.testing.synth import CounterRNG


def _pos(m) -> np.ndarray:
    """MATLAB's 1-based [y, x, z] positions as 0-based (z, y, x) (subtracting 1 is exact)"""
    m = np.atleast_2d(np.asarray(m, dtype=np.float64))
    return m[:, [2, 0, 1]] - 1.0


def _generators() -> dict:
    """name -> function making the array whose hash (of its MATLAB form) is pinned"""
    g = {
        "uint32": lambda: CounterRNG(0).uint32((3, 5, 7)),
        "uint32_seed_stream": lambda: CounterRNG(2**63 + 12345, "decon").uint32(1001),
        "uniform64": lambda: CounterRNG(1).uniform((7, 9, 8)),
        "uniform32": lambda: CounterRNG(1).uniform((7, 9, 8), np.float32),
        "normal64": lambda: CounterRNG(3).normal((5, 33, 17)),
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
            (40, 32, 48), _pos([[4.0, 12, 8], [10.6, 24, 20], [30, 36, 32]]), [2, 1.5], 1000, 100
        ),
        "beads_2d_amplitudes": lambda: synth.beads(
            (1, 40, 30), _pos([[10.3, 7.2, 1], [30, 20, 1]]), 2.0, [500, 2000]
        ),
        "texture3d": lambda: synth.texture((40, 32, 48), 3, (0, 200)),
        "texture2d": lambda: synth.texture((64, 50), 2, (100, 2000), CounterRNG(7)),
        "texture_nz1": lambda: synth.texture((1, 20, 30), 1.5),
        "decon_psf": lambda: synth.decon_psf(),
        "decon_psf_even": lambda: synth.decon_psf((24, 32, 32), (1.8, 1.2), 1),
        "gaussian_filter": lambda: synth.gaussian_filter(CounterRNG(9).normal((11, 17, 13)), 1.3),
    }
    for dt in ("logical", "uint8", "uint16", "int16", "uint32", "single", "double"):
        g[f"random_image_{dt}"] = lambda dt=dt: synth.random_image((11, 23, 19), dt)
    return g


GENERATORS = _generators()

# made on Linux x86-64 (NumPy 2.5.4) before the convention change, of MATLAB's form
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
    assert synth.array_hash(to_matlab(a)) == EXPECTED[name], (
        f"{name} changed (dtype {a.dtype}, shape {a.shape})"
    )


# ------------------------------------------------------------------------------------------
# 20 parameter sets per generator, in the old (MATLAB) form; the new call converts them
# ------------------------------------------------------------------------------------------

# MATLAB sizes [y, x, z, ...]: images, single voxels, lines along each axis, 4D, odd sizes
_SIZES = [
    (5, 7, 3), (1, 1, 1), (7, 1, 1), (1, 7, 1), (1, 1, 7), (4, 6), (6, 4), (13,), (3, 5, 7, 2),
    (17, 13, 11), (2, 2, 2), (9, 8, 7), (31, 9), (1, 40), (40, 1), (10, 10, 10), (3, 4, 5),
    (16, 1, 4), (1, 16, 4), (23, 19, 11),
]  # fmt: skip
_VOL_SIZES = [  # 3D (or 2D for texture) sizes for the image generators
    (12, 10, 8), (9, 9, 9), (1, 20, 6), (20, 1, 6), (15, 11, 1), (8, 13, 5), (16, 16, 16),
    (21, 17, 3), (5, 30, 4), (30, 5, 7), (7, 7, 20), (11, 12, 13), (2, 3, 4), (24, 18, 12),
    (13, 9, 10), (6, 6, 6), (19, 14, 9), (10, 25, 2), (4, 4, 30), (17, 11, 15),
]  # fmt: skip
_DTYPES = ["logical", "uint8", "uint16", "int16", "uint32", "int8", "int32", "single", "double"]
_INT_RANGES = [
    (1, 65535, "uint16"), (-1000, 1000, "int32"), (0, 1, "uint8"), (1, 255, "uint8"),
    (0, 2**32 - 1, "uint32"), (-128, 127, "int8"), (-(2**31), 2**31 - 1, "int64"),
]  # fmt: skip


def _frac(k: int, den: int) -> float:
    return (k * 37 % den) / 10.0


def param_sets() -> list[tuple[str, dict]]:
    """(id, spec) of the 20 sets of every generator, each spec in MATLAB's form"""
    out = []
    for i in range(20):
        size = _SIZES[i]
        seed, stream = i * 7919 % 1000, f"s{i % 3}" if i % 2 else i
        for gen in ("uint32", "uniform64", "uniform32", "normal64", "normal32"):
            out.append(
                (f"{gen}_{i:02d}", {"gen": gen, "size": size, "seed": seed, "stream": stream})
            )
        lo, hi, dt = _INT_RANGES[i % len(_INT_RANGES)]
        out.append(
            (f"integers_{i:02d}", {"gen": "integers", "size": size, "seed": seed, "stream": stream,
                                   "low": lo, "high": hi, "dtype": dt})
        )  # fmt: skip
        vs = _VOL_SIZES[i]
        out.append(
            (f"random_image_{i:02d}", {"gen": "random_image", "size": vs, "seed": seed,
                                       "dtype": _DTYPES[i % len(_DTYPES)]})
        )  # fmt: skip
        n_spots = i % 5
        pos = [
            [1 + _frac(i + 3 * k + 1, vs[0] * 10 + 1), 1 + _frac(i + 3 * k + 2, vs[1] * 10 + 1),
             1 + _frac(i + 3 * k + 3, vs[2] * 10 + 1)]
            for k in range(n_spots)
        ]  # fmt: skip
        sigma = [1.5 + 0.1 * (i % 4), 2.0 + 0.3 * (i % 3)] if i % 3 else 1.2 + 0.1 * i
        amp = [500.0 + 100 * k for k in range(n_spots)] if i % 2 else 1000.0
        bead_size = vs if i % 4 else vs[:2]
        out.append(
            (f"beads_{i:02d}", {"gen": "beads", "size": bead_size, "pos": pos, "sigma": sigma,
                                "amp": amp, "bg": float(i % 3) * 50})
        )  # fmt: skip
        tex_size = vs if i % 5 else vs[:2]
        out.append(
            (f"texture_{i:02d}", {"gen": "texture", "size": tex_size, "sigma": 1.0 + 0.25 * (i % 6),
                                  "range": (float(i), 100.0 + 50 * i), "seed": seed})
        )  # fmt: skip
        psf_size = (vs[0] + 10, vs[1] + 10, vs[2] + 10)
        out.append(
            (f"decon_psf_{i:02d}", {"gen": "decon_psf", "size": psf_size,
                                    "sigma": [1.0 + 0.1 * i, 1.5 + 0.05 * i], "amp": 10.0**(i % 5)})
        )  # fmt: skip
        out.append(
            (f"gaussian_filter_{i:02d}", {"gen": "gaussian_filter", "size": size, "seed": seed,
                                          "sigma": 0.6 + 0.2 * (i % 7)})
        )  # fmt: skip
    return out


def _from_msize(size) -> tuple[int, ...]:
    """MATLAB's size [y, x, z, d4, ...] as the port's (..., d4, z, y, x); 1-D stays"""
    size = tuple(size)
    return size if len(size) < 2 else (*size[:1:-1], size[0], size[1])


def run_new(spec: dict) -> np.ndarray:
    """the generator of the new API on a spec of MATLAB's form"""
    g = spec["gen"]
    size = _from_msize(spec["size"])
    if g in ("uint32", "uniform64", "uniform32", "normal64", "normal32", "integers"):
        r = CounterRNG(spec["seed"], spec["stream"])
        if g == "uint32":
            return r.uint32(size)
        if g == "integers":
            return r.integers(spec["low"], spec["high"], size, spec["dtype"])
        dt = np.float32 if g.endswith("32") else np.float64
        return (r.uniform if g.startswith("uniform") else r.normal)(size, dt)
    if g == "random_image":
        return synth.random_image(size, spec["dtype"], CounterRNG(spec["seed"]))
    if g == "beads":
        sz = (1, *size) if len(size) == 2 else size  # MATLAB's 2-element size made a 1-plane volume
        sg = spec["sigma"]
        sg = sg if np.isscalar(sg) else [sg[1], sg[0]]
        pos = _pos(spec["pos"]) if spec["pos"] else np.zeros((0, 3))
        return synth.beads(sz, pos, sg, spec["amp"], spec["bg"])
    if g == "texture":
        return synth.texture(size, spec["sigma"], spec["range"], CounterRNG(spec["seed"]))
    if g == "decon_psf":
        return synth.decon_psf(size, [spec["sigma"][1], spec["sigma"][0]], spec["amp"])
    if g == "gaussian_filter":
        return synth.gaussian_filter(CounterRNG(spec["seed"]).normal(size), spec["sigma"])
    raise ValueError(g)


def run_old(old, spec: dict) -> np.ndarray:
    """the generator of the old API (module old: synth.py before the change, MATLAB's form)"""
    g = spec["gen"]
    size = spec["size"]
    if g in ("uint32", "uniform64", "uniform32", "normal64", "normal32", "integers"):
        r = old.CounterRNG(spec["seed"], spec["stream"])
        if g == "uint32":
            return r.uint32(size)
        if g == "integers":
            return r.integers(spec["low"], spec["high"], size, spec["dtype"])
        dt = np.float32 if g.endswith("32") else np.float64
        return (r.uniform if g.startswith("uniform") else r.normal)(size, dt)
    if g == "random_image":
        return old.random_image(size, spec["dtype"], old.CounterRNG(spec["seed"]))
    if g == "beads":
        pos = spec["pos"] if spec["pos"] else np.zeros((0, 3))
        return old.beads(size, pos, spec["sigma"], spec["amp"], spec["bg"])
    if g == "texture":
        return old.texture(size, spec["sigma"], spec["range"], old.CounterRNG(spec["seed"]))
    if g == "decon_psf":
        return old.decon_psf(size, spec["sigma"], spec["amp"])
    if g == "gaussian_filter":
        return old.gaussian_filter(old.CounterRNG(spec["seed"]).normal(size), spec["sigma"])
    raise ValueError(g)


PARAM_SETS = param_sets()

# hashes of the old code's outputs (`python tests/testing/test_synth.py --old <old synth.py>`)
EXPECTED_SETS: dict[str, str] = {
    "uint32_00": "3084a6b8a6ae3830b745ceff11953d720094a5be8a5096ebc1bb89248c89388c",
    "uniform64_00": "dec4994596a43186efbfc070e6809046d8ed4f08ac72a8604dcff4f4f9805684",
    "uniform32_00": "be4aabd08ef659dcb3b1ccb7668b1278f098636a207b450ac4cd0bb5824abf55",
    "normal64_00": "ec18dca567635d6dfbcefe05f98e5837a6513f1c41c4e88071d31fc87336a881",
    "normal32_00": "911e081c6865ec574aa9fe3a4088cfac3cb118f68c46bcf86caec3c9ab4ebcb6",
    "integers_00": "89d704a78e6cfb6f2229bbfa746b80f9ff635091de1b440eb9a04abe5796f102",
    "random_image_00": "e985745531b0ed25f1d04c265f062f9b668dac82978d5cae240fbc4d1a21732e",
    "beads_00": "d580080d6214fbbc7280d198905d1351ada683394ed7193c89b16aac59e2e8c5",
    "texture_00": "047f35428b58818fb3f7dd5b8f45be2f0757852c7a9b86af1aefe38170fdd750",
    "decon_psf_00": "9e5eadbdd2895a28f22183722b23b6afbfef834f1406893d97c7725c0fe46ff7",
    "gaussian_filter_00": "dd9ec804bcd769dfce4ce8ad5438f339330eb22b8495643967c97fcebb358f72",
    "uint32_01": "c258531bc9fb7c3033a5357bf6d9812197006891f06018149b03f8e3f0b9f8c4",
    "uniform64_01": "66d2556f4a6ab1f28a068c785d0705fbe7d7ad0bc79db7a8cc61b8591ca831e1",
    "uniform32_01": "959ed04f27b1a2846bf8d0a9965a111f5e5cf41b354aa61ab402be988a27950e",
    "normal64_01": "380d4878c1483e10db59d43128bece2da01e4017736b48a534c85d086544cfed",
    "normal32_01": "1e3c5693bd8f75c646fb1914454ef22b69a056519ab217e9b4ac63d7cd8770ec",
    "integers_01": "c3dce2eb6e07dc13450f7f4e06132c5874073b7efbcac16548253f2dcf99c0d9",
    "random_image_01": "6a6d4767f03382bf71f82baec82a6a0577664c2c71dbe998ae9d7009716b2347",
    "beads_01": "ebde59fea2a9f333c9af4f5bd242b1b840bcdb09b0f65c6d4b8c36fb672c0d64",
    "texture_01": "d72d2020132d01121d4f9e8583b23770060b2aa75356e6a1b2f1244eaa37fa63",
    "decon_psf_01": "ad6bab41bbff6c1e6556752637cc6e35b435c60897625b1f501431e8c8c2a8a7",
    "gaussian_filter_01": "f715cf047f229db52889cccc618303a5d9b2f9238f5459c61a90258d3be8da4f",
    "uint32_02": "fa1fd1cc4e48097d9c562e63ec4983371f5d8b800e14c998f778d4d181db6341",
    "uniform64_02": "9fb36e02374f0c4c9bbfcc057c6102617f2efba230d9c0eca5cd5a3655a368a9",
    "uniform32_02": "ca937ffe0020899ecfe80cf4bb50908a61c977a716b4349d2c1fe641b12c2563",
    "normal64_02": "08645273e2cc8bb66b411ec823edb4a4c63772baca6ee7869c3345861d5c784a",
    "normal32_02": "1166c5482e4a03dfad8b71d23aca2bdd7a9084162a8ff3a5e60f596998fcddac",
    "integers_02": "613e81b1148b1887d959d04398c4f69a7c68b1976f382cac18bf5ca87a8fe591",
    "random_image_02": "dfb832fb60cb12d7f989e8c85485a8b27b2c11796b082e7bbbe6baf8391dae38",
    "beads_02": "45f6d035b3b4ebb447667ca9ce963804fef3c8beb7c62874cbac7f0525b379dd",
    "texture_02": "93b02f19772ec5c93ffc3bf888907c0289728e0723d4d9ada91677021861e61e",
    "decon_psf_02": "73b64b676dc509f602d4155c1309d0183d8268ffa50fea9066998fb80a548a67",
    "gaussian_filter_02": "5c4e5cb972136cfccc5527be0da44df9eab358ac4c1410ec44f94c6c45063a40",
    "uint32_03": "89ce56c41948c42352e25ab28680515ee485bbe9600e3d9931ded1d9c492055e",
    "uniform64_03": "c669057a6f8addfd94c9d431473852eb62c48da2b473d45d177abf1d32ce90a4",
    "uniform32_03": "429699249815bc671576d616a5144869d6873c48bc9b998ebaec2722071d7477",
    "normal64_03": "8407c3974f38623542ce746ab9f712d16363609d76d20d1da2901b86e0b5afbe",
    "normal32_03": "1e850da33ec0f3ede925fbc6050237340a25739303fac32afb82ca9d63ff4c95",
    "integers_03": "10afb640929fbb741bc29962fc60ab43776db45ab407fcde60c4721f914ce022",
    "random_image_03": "34cb3f405b7b8f4f3ab70a03a0fcf459eb7130162438b2739368e5a5117338fb",
    "beads_03": "b15dbec73f032078036fe9921fbc3aa7762d03f729c44fc3227c016e9e9b937e",
    "texture_03": "3e567ea1bf2231958ca54ad34f2f37022fdf8b358258ad6f3a7fb832b8c0e0f4",
    "decon_psf_03": "5f5f87228c55d397694efe85bdaeaabc627aff60202404e3404ceaf83adf7d9b",
    "gaussian_filter_03": "2e78e6bd6944e350f700438ab3c9d9f4df72e2ad39db77b7785e202d237e8273",
    "uint32_04": "2cf946356dbd0b8ca019e75f2b5f75dbd22a11d34c7687bc6d6d1c70b3394610",
    "uniform64_04": "abf8297629bbe30714a59c2aba813cb93472060e3febabf893f5a979d2c877ce",
    "uniform32_04": "3533326b006251a2051ac118e18a1b2dde7b40eb7444cb5b21b23fca3ef380f3",
    "normal64_04": "ccb6c3b66db8bbd8f09c2bf9a1ba27e20da925015b97f31896d01ab5177fc0f2",
    "normal32_04": "3db742dc6e5668fcd4e010be8682f1585d6e7f4bc6ad8200c3013bff5bcb1ed2",
    "integers_04": "2cf946356dbd0b8ca019e75f2b5f75dbd22a11d34c7687bc6d6d1c70b3394610",
    "random_image_04": "ee1ef4a7232a96e71cb1ff1b8c837af74977ba6e94d47f5016a355b08728946e",
    "beads_04": "463f06450082eab9ca1676b77a076886c499659a64a29f659b60ba3015f17f29",
    "texture_04": "c4d5bdf7ea757b9fe88021e96e4c5de65dbeec3da630a07384b779acc7f2f38f",
    "decon_psf_04": "9ab7fde416ee5b0e438bd53a2efa1605201d00893c44b70dca98384dbdb4cc92",
    "gaussian_filter_04": "f937c8f2b895aff397bb68b958ff8d1a468d2e8df6a4decf7a77d3a677750471",
    "uint32_05": "d43f6a0a5c7a47cdef5ab8c7af618c314072911c937ba7fca431e8fe09927d1a",
    "uniform64_05": "fa2a0ad417c17c19eacddfc8380970ec5ce22e7a9cdd4f08e0a6b748d14f29ae",
    "uniform32_05": "fc2faf52874745733cf2002f8679d732a3938899b295cff510f5f47798c1337e",
    "normal64_05": "8830ddd2211a5b37b755ab2ca72e753250c9a3dfe3cd73c8f3aa0ab2a1b8ab13",
    "normal32_05": "b84605138213d3f064cc5aee76b6aa5975ea4cb27be20cbc11b650c433ebecb3",
    "integers_05": "de95b168ff70a2afcc6f7328e9388ab70deb30f49d91ae7abe9dd80f5b989b66",
    "random_image_05": "732f83fd6143f9e99072b6a1aa254d24fe31682ea49c949f4ca755706a60b91a",
    "beads_05": "2c3fc97acb710e0567b2f70cf76d7fcd2f3a5d914c7def4d8706c96b330fea1e",
    "texture_05": "5a559e11e9b1c678f118ac3622e87a0de0e80616999d41250e1fb4fe4ad84d06",
    "decon_psf_05": "34215343e5aa90578c0a789ba07cdf1acfac73d958d23605b1b4c7956a58bcc3",
    "gaussian_filter_05": "1faec446a1ff98fdbff2c38c68ee5b49a72237ddb846befce877991b13fae11b",
    "uint32_06": "14c4a59e70ff8bcfc16c2461377442579654cefa2993479892c85d728070550a",
    "uniform64_06": "145837c3585adb0111eae365a31b86a3330a2ce9213e350eec7022ecdf739f8e",
    "uniform32_06": "5114578531fcabec375b933994512796266a1412abed1673fba78155612dc96c",
    "normal64_06": "b1cf866e57b1e2ba7ab871719684c0ddbd3433c3050209217a1db210c954ec40",
    "normal32_06": "e5dcab7e86257baf6cff371e4f990e2a7995fe0318568c8c2fc77be368a7a412",
    "integers_06": "424293186378941039cc8fc27a71d4b4e005116b5bce6e7536fc058c9bcc14bd",
    "random_image_06": "a64260332e6ef770ed9003287268cbbf112683e4792775ead13dc5b6f36b1dbc",
    "beads_06": "e29dc2be14d459a64e0bc1aae3a3b1a070230b09de9a6c2ad5f6f2cf6c39f15a",
    "texture_06": "fd14d3574265b7c99d4a5150c04fddf7226f9c4cfdd9bf33eee899ac164de471",
    "decon_psf_06": "587c9dd26217c4e008ff7f62ad3bdaed560ad9808a928755ca4cb0786a546891",
    "gaussian_filter_06": "c039e86b421e62f2d91b7f84d70ce8a3fff42649367797d4ec4f5fc181bbd5c5",
    "uint32_07": "acb6b1481edc5a199cbefac5a67392c5cf6cdd4c7a226bd49a82ab1bc01d2d4e",
    "uniform64_07": "c9fd8d3d6b4414c3c4b1b1e42555cc925348da898ea26927b6f56c2df5338d0f",
    "uniform32_07": "cf065c01b2f5fbab8c6dd85aac4ae91110d1bd7479c54a6ec56a2b02e8683d8c",
    "normal64_07": "ce48b2a4b1845f6a86480bacefc8718b04c95a7991cb1a2fb5c61f0ffef2236c",
    "normal32_07": "fec5121ea3a383772141a70429bdb76274edf8f65672983db2454f60538c3432",
    "integers_07": "33f108679229520e951569135d43388fd61514545bf1b27aed044242f9550fe1",
    "random_image_07": "20fb2719f65cd4b004fec20df9787dfd60ab259df5d59e99f7e2b47656d375f6",
    "beads_07": "29b286e825a07a5c68392b615a63a98a4fdd16356226875f3c2312883e4f82fc",
    "texture_07": "6445698872abe19324039aed09fcdee59d14e10bf6ad1931e803d76ef082ec13",
    "decon_psf_07": "2a75d50a56096cfedeb9596e321492efa34879f7cf50269a637b877470db761e",
    "gaussian_filter_07": "a68e352b9034665a73a39177a4c202159ab1948d935f70b9ba2e6bcb4e57f3de",
    "uint32_08": "d670c6bc55a61b1f3f646c4a39520722ee0dedbcf6c69f359a6fdfdbca14be17",
    "uniform64_08": "4fbe3675070166fc7cad4a9c983b7d7377bb0872a0ee98b7b55c728959747a77",
    "uniform32_08": "6168154379709981800111e00e5bb266a4fb3e4bed4f1c53ba9981139e9160d3",
    "normal64_08": "b4bf6bb5536b7395cdac29de2675f9b65913fb60cfa0e70315c1bd174d5aa62f",
    "normal32_08": "a0f537b309b9eb36acd2c138a8e3e39995c5f779ba0e881ec71cdca7c53c860e",
    "integers_08": "2951c4459414ca83fdafbefc56cdd1c7205be5532ebf3985790f1fd24a21abd5",
    "random_image_08": "6db47e1ac0b2a93993b3a2a1f41b28f64a3f2d548b32e5ee3ab9e48fbc3da9a4",
    "beads_08": "5166ed5e427a6cc68083e40aefd0dd43bc2882c40e5865ea46d194b1360b93d1",
    "texture_08": "f510aa08912783aa64c787bb8f1974afcaed7c00325ec8fd61a32421317f02c8",
    "decon_psf_08": "2fb83275ad925e8eda377a70acdb175785edb40e022d7f8468c7a7bd22d1a467",
    "gaussian_filter_08": "41b4ee6b507af72713c1ba080a5e2a6dd4aa2537ce9c0ad9b911fcc50864b857",
    "uint32_09": "4c545c4992c126bd39a6e2ddc3194acf0d4e823932d8eeaa6fd02ff91b0e86a0",
    "uniform64_09": "03456b05c2a498366b5c12ab42e9869538dba604a597d92b7d917371d9fb7b23",
    "uniform32_09": "283b7c30b25edb78c47a8fac7dadb0f68af1328001b32f2d8c26756cfa193100",
    "normal64_09": "8251e183536981548b088605601f7fbf58b00179d6799373a69cc091fa2db0eb",
    "normal32_09": "b30a36c141bbf7a4e8ef865b04cfa7a15970fc96378581e18f63ba43056252d3",
    "integers_09": "0dbbb8531f3b2f0aab7383e2d395ba42d0f3b9a507126f54c6541e05e62e763a",
    "random_image_09": "5003cd29098f01313f6884e4f3ffa8cbff079ecce048066bcdcbfda8a3286586",
    "beads_09": "75764edd9b7297b8ac09d72f871f815d4796769a3ca60095d9c75e4a8a7e870b",
    "texture_09": "ea5f7b0825680547cbddc1eb2aa0d32a582e32a6399f3e236c19bedb6e62a034",
    "decon_psf_09": "29dc3cda72332eeeab9fc595bed666acf981b77e3ead96059aa72a6445d07046",
    "gaussian_filter_09": "d622e1cd45f40a0e2b3d851aed8168fb1c4051da18876ef3d0ef7008db6bef2b",
    "uint32_10": "a2cc5f61541b7f3e0d3c20e9e763910a40b111fa4b23e4dce859ec4555ff3625",
    "uniform64_10": "0e844365fe1287478f4defd317a8836a95fee4a8f0ff5532d422b1d6c59d85db",
    "uniform32_10": "005ae6aa11e2a3c8c5ad007f683d7e911fbd9341c681f7d95c218a24051764ee",
    "normal64_10": "681443ec7f43fea572e08156375477679bf1b07566d6ba1761fab4635d2a4a6f",
    "normal32_10": "f2bc366c73da61d3cf020994aa77a6c7cc36879e04c0c650669c544f2d60587b",
    "integers_10": "8312d852cf26ed66e1500689ce479d984ca1b4433296f4991a20c4c18e608505",
    "random_image_10": "1e7e07c7e71cc14e83b0eddb1c8f48cdec9cffd8e9fe02beb904bf26b462fdc5",
    "beads_10": "b68a291dc106dd4135825c409e996fe1687f5736780a9321ae325ac2f17d43d9",
    "texture_10": "a362df68827f9c459859ce37cb3ef2c389335008c6a768e1c89e6391a6db679f",
    "decon_psf_10": "aa1cf927160313265e336b44e6e11f1cc25a7149802e69dcc0acaf156341e38d",
    "gaussian_filter_10": "372d929e33ce645b4c4c9601e2da01762f9226866e98618dee4fb17b9b3f9b02",
    "uint32_11": "e8f2c304cfa9aee9138d5e9ce9a18be64ad7bf9523e035632e7f8c0a4759ac55",
    "uniform64_11": "52db6c4a1edd9f0348a5c847032b769a698626acf102e7333f15824e9805801f",
    "uniform32_11": "f0112a6193a8ffe2d62b671cb7e0430268da1d2ae8cb9d843a7d44d0e0c8289d",
    "normal64_11": "4b56804595924d1551fa1b34fdcd66ec7b9d71fed8634c37d24dc1c7481e69de",
    "normal32_11": "3cdf3e6e66c927d558f2bb5f0a4b352eada8e55a26137f4e7f029661ccfad534",
    "integers_11": "e8f2c304cfa9aee9138d5e9ce9a18be64ad7bf9523e035632e7f8c0a4759ac55",
    "random_image_11": "b058b90876729a3c13d2be09f17423d97f601af524fa444a4ebc9812ec1edccf",
    "beads_11": "650992cc9a611748f228415a9177e76b7f52afb297bfb800355ba042e2a4cf5c",
    "texture_11": "8a918db9a7e0b966d9b56775d32f23e852014cb01dd9408e3b00efe4a49da1b9",
    "decon_psf_11": "017f74591a4ce434bb876721c7ae54af95217ec3274552e5fee05e9d1b3146b1",
    "gaussian_filter_11": "83cd324918bbbc36cd5e4e0520069ce31eed94e8e6fd42108922819c0e7153bf",
    "uint32_12": "88a2c4deb8ec7babf58fc46b3a8d95a08bfaaea4cc1d92a52f06256aee67620c",
    "uniform64_12": "a07664d02462bc0a10c0d97c2e3908acf11a2bf7001bf214bb49fd3da5ccfb4d",
    "uniform32_12": "7c5787281af49173b9428e4049f5613b2515ac0df0ee06152b058643c7e635ac",
    "normal64_12": "d367ff8b1e2faa25a901e116d003bf8af3af068d5f46e367c1aa4f05c90c3255",
    "normal32_12": "c1ee3fd85fee12e4af70c858538f814b613cbac2d07bd4ac45babada84c5afb6",
    "integers_12": "9adb41aa30f791b273efd6372aa094e3ec5719cbc6f0ef5e4410a1c829127b5e",
    "random_image_12": "5de8f0c9709148d73e9dc565dc28ad772e497bd6357ddb958c6648202ab8477a",
    "beads_12": "092cbcdb3a93e20fe893b17ef2e70d5b7cfb43b59eb4ac3ec0f8d43ac1364432",
    "texture_12": "749ba11f6a50dfc35680c5b0ae1f4e6c9efb5d87cf30d2793d9b0643edfcc433",
    "decon_psf_12": "6ffd76454b60977bf7fcc59e750b5c6200922c6bc73643a5104ba1210c57abbc",
    "gaussian_filter_12": "04fc7cff54f7eb93273b3d64cb08f5ede8b96ac006e466fe578ae939566939db",
    "uint32_13": "2108bcc3a992f58882f5943a541d2e98bf56b79b9d03b62449833fe2e875ddd9",
    "uniform64_13": "17a28918120f372e94a0d7871d934585f6c786fa6090f8cf5518ddff146b0d60",
    "uniform32_13": "6c6d844f5405df11b5ae563b2f0f28820f0b19c85e14a6ee363e94602d8c8ad5",
    "normal64_13": "322c32112a02b91469fe902cdaa72ce180d4e1b4b656de9dbccdbe3e2f62bc6d",
    "normal32_13": "add1fa0e6abdb9ff2eafa47d0c64804c6ca4226ed2078e1291000660869b035a",
    "integers_13": "ff066b194637009ab5900969543bffa7d61ab0a5012dcc0fb4c8adc2fe890c49",
    "random_image_13": "0994d35e104f27f6f11f136cc62c62ac0f4c07ced8250ea2963092a1275e7e39",
    "beads_13": "bb7d92e88ea51773729903d237bed6ddd1215c9538200d6e5a705745c4679321",
    "texture_13": "08f9de2f4185f3f80c547994d1b214bf2e5e4b6af5482f8c05c0c4c8c9d2161f",
    "decon_psf_13": "2c5f5f291ddb6c7b9c80dd64eb010fb67e03fc12def8e34ea4f2eb6951d8edfa",
    "gaussian_filter_13": "3558089783590867ff677df57555099b525260ed21265428553618cba606b8b1",
    "uint32_14": "8d25555fd33897bcd658c1a589487f4411249ba5a038615bd0b04894345e17c1",
    "uniform64_14": "4ff094687cd8ef4400c9494ef7e005192f94fc2ebac67b9df8fae50cd06391f7",
    "uniform32_14": "fcb5e3404175c4611cd06f778dc3fb95da4b18b65c930020f9cfcf56a52defc2",
    "normal64_14": "505c8e2f5fd32c8b3dcc068393d290a88105cc2119ffb23a1f4795a079a9dce3",
    "normal32_14": "309e13d0cd9a45b891b6787935290cf2b8eb417f16f1c1b4855da651ba0db9ad",
    "integers_14": "ae74b5e88397eb3cc1708a76096fb5430f6df052ac2c458b1d020b5121f881d0",
    "random_image_14": "5d02fc9024c782ab72249b258139b24e67b7be722de70c752d78a3b5249720ce",
    "beads_14": "a4f400d83d280ef42954dc91cac04faacf5e0f52f1612fa5d3035fb0bf80961a",
    "texture_14": "6a775f8d3b2126cf48ce7c7c4d2d424f07135b2e9b2fe985f7d8876e18c216d1",
    "decon_psf_14": "95fe82ce5ec77a166d291306915504b9f7368e483115dffb3d8fecbca3b97f00",
    "gaussian_filter_14": "6b1c01a5f8135b4dbb54525852f6cf924741e2c75e8d5b89d2c9dc16771698c4",
    "uint32_15": "2f73613602427de12d78c149d38eecb6759e3fda2b38d4efba5cc2046ebd20b4",
    "uniform64_15": "38ebc3a28a7c24ac72be298048d7db86a19bbabe67db5325c96af60738357015",
    "uniform32_15": "c5e0d395df05a9a303a208ebe90e4f6570007d2c1da6ed4479d84919c40d4ce0",
    "normal64_15": "24866ae8019dbcffff3b3e03166ccade081673058d355d3a472c4ea1e487c34a",
    "normal32_15": "0eb2be2b00f30e3f699b55204e4260b2f0a41789e568fea86cc8165b78d180fd",
    "integers_15": "89703fc4b21ea9a116f4e1ce3deed2a0c401b7d692a06ea660ee2d22a1addfa2",
    "random_image_15": "e5f58a4a7800b62d5e5417b977fce4dffde35dc49ca7c302f6ba775a4b7795d4",
    "beads_15": "59618dac658de31c819ef09865ba113170be1446cdab0281d45f427e7368dae6",
    "texture_15": "61e3de95853d530838d8fbf931d57f7a950ef5712bb10cb940ebd776b3fb026a",
    "decon_psf_15": "cb32bb1aaae721839f518941024b2d2e1f020df23dfbc62f1a4709753ffb95b1",
    "gaussian_filter_15": "1406f6398f0da5f743bf049318f45a3393550010796fd6ed472d96089956b8a3",
    "uint32_16": "d9516b5b214dc4a92514fcd517839447126b7e384c31d7310ba25d8e12b23eea",
    "uniform64_16": "23e6614178cc2af80fa41fe7cbe112530753f9c6309d4879e3fb7a0637e245a8",
    "uniform32_16": "f4868fedafa9eefe8a0425cc008b42ad55f5f11504126b062488365c5cdac018",
    "normal64_16": "3890056840a8dae3beb17c6c73e4b9667c4f10de211ff36c765f8a324331e4c6",
    "normal32_16": "31211d164d95d64055b18aaff95f85cc8a6fbc3e3e19a0f98c9cc8d16123a2ae",
    "integers_16": "cf4fa285b5f91e16c2281324b630b8a2382081a13f6cd5be2868bcef6b144686",
    "random_image_16": "44b09f5ba70e9a5f582864d263d92120ce5e5ff1042c405f0e3d451ffe80c389",
    "beads_16": "3d73638434d272cfc295adf2869391a99d66b6d7e9a3ff00899bb77595e3081e",
    "texture_16": "0ea6a865da1cf69f8ebd780cf91e1fae5b51ba41be63490432ece06693987aa8",
    "decon_psf_16": "6b9f8d97294814dfab7fb3df10956fcc2dcb30342c63827c03783f29d079cbbd",
    "gaussian_filter_16": "3c91c9cc4951c39aa5998e5fba45746acd9bd41af6df23891c8247cd74ee56e0",
    "uint32_17": "53b9b8c2c7d94ee4a766554662eaaddb73f3b14a0dd7e998adf12323746ab647",
    "uniform64_17": "9d6953b55516c80947eef2ee9b799c8e692c44dd84e0373359fd0cf0c12bde9e",
    "uniform32_17": "573c016e4345dbaedddbccf7d124edcc2d4d977e4cb5e959125747499b98e88e",
    "normal64_17": "8c5a383ddbe21990baa2c23de33bc9558b9d172f7929726206bad45c59ea08c5",
    "normal32_17": "66c2a5b86d0402663eef1e988279d5705af75874982b77ddf6ad32cc071f6c44",
    "integers_17": "185d88c6c4d370bf91f9f809c2e7bc5895c8197b29f387fe5cfe26e54d20d218",
    "random_image_17": "debd708eccc2d4bac9d11d86b83c56476189442011cf0a52d3a66361226ee6ff",
    "beads_17": "fb69b7cf643d2d633a3c4be7cad201b0e0074b161978e843198d2b0da4319760",
    "texture_17": "e6a945e573e39e8914a51e1fcf6404b1534b108c7ad9c3ecfe704f0f7893868e",
    "decon_psf_17": "c4825cd1a8bc06f30e49c8d2b06b20e50183bf13b82cb53d6f720786377e1fe6",
    "gaussian_filter_17": "3e2beba0aac52addafbf202e5d9f97639e08c17fa346a4f3d0905d11fde7d6a5",
    "uint32_18": "65f1a4b8ba6580ae6fdc6fb24512eb93fd1b7f4845bd09e494f9faab074a1ecf",
    "uniform64_18": "0cebb7228b1433e92d9ceeaff9d349d8f70218563313e78eef7cf446d2b717d1",
    "uniform32_18": "cd513cec8a97bea6719ecf8162dd7cc9a7dc2f47c9e345468ed3ed05f1082ad2",
    "normal64_18": "42d8cede2cd1b70ae3a2f8fae32e32b805be48d48b04e0070bba27d78d0b8e64",
    "normal32_18": "2ecf4b9ba2e22c4fa500306ac9dbd8ad0b447f394db1e987a72bb0ef2ba8d9ea",
    "integers_18": "65f1a4b8ba6580ae6fdc6fb24512eb93fd1b7f4845bd09e494f9faab074a1ecf",
    "random_image_18": "b3aa3fbe3ace1680930876d6877c3c5f51f7efcd3d9dcdbea3b65fa40f35f73d",
    "beads_18": "d3a0d1ef2a7e19c91cc1f96d9451811ee7917f33e0f333f0c281d35e9598ef10",
    "texture_18": "1d80b09ce23396e0f9205dc65c7e20935edad5cf7a169f626ed0da6c56775417",
    "decon_psf_18": "728be20764173bb8be72cbf7b1fe7fd8c687db86d5df63125f02d07954b1df03",
    "gaussian_filter_18": "44d5d68037bcd44f44d729ede06fd0a2818a9e508390006cfbfa60a681e040f6",
    "uint32_19": "4da9e76cbf863ba96c3be8fc7be2b4b02d850edd5fe2f43bf54d69656e9a880c",
    "uniform64_19": "4da67e990954e5b2bc0ffe621ed010599827c93f0a3b67e46c46c20c5ad9d561",
    "uniform32_19": "32c52de24e4d7547c2c774c9edcfc737ccc495515518139ab83a6084baef5448",
    "normal64_19": "f37376f8b097daa637ec499ab299d8f2007f674796135b8de8ed598ea51b3f80",
    "normal32_19": "a67c3f1ccc5f44699452cf2b45cc585bc683651cb2a6686e5e65a3a5d4a933be",
    "integers_19": "ac99ed4333247d52dd689d9f0c7cad4f09e39ab30dea9e0efb687e30808f4760",
    "random_image_19": "e738795d8e7e1e20471405471081a5e94081bad8c1a10891e4152605a1e30eee",
    "beads_19": "e9bbc9bddbbaf40c74658a7845f74aa5f3a893cfba93d20b7bb1f4d6fbe5201a",
    "texture_19": "f39cf85e3eca0072dca105509f304a779ae5b3deeb9c1b9e80f2282bfe8a93a8",
    "decon_psf_19": "d921789b57ac578876e6811ad9d2b45f7eac53b4281755e01ccd53f8c5785362",
    "gaussian_filter_19": "4b238037736a08088ff4248a6731a2e7dd9ee2fd1cc6ce2ee68b1f53b3878fcb",
}


@pytest.mark.parametrize("name, spec", PARAM_SETS, ids=[n for n, _ in PARAM_SETS])
def test_matlab_form_is_unchanged(name, spec):
    a = run_new(spec)
    assert a.flags.c_contiguous
    m = to_matlab(a)
    assert synth.array_hash(m) == EXPECTED_SETS[name], f"{name}: {a.dtype} {a.shape}"


def test_every_generator_has_20_parameter_sets():
    gens = {}
    for _, spec in PARAM_SETS:
        gens[spec["gen"]] = gens.get(spec["gen"], 0) + 1
    assert len(gens) == 11 and set(gens.values()) == {20}


# ------------------------------------------------------------------------------------------
# the generator and its distributions
# ------------------------------------------------------------------------------------------


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
    # MATLAB's linear order: element (z, y, x) is value y + ny (x + nx z), a C array
    v = CounterRNG(5, "s").uint32((20, 4, 5))
    assert v.flags.c_contiguous and v.shape == (20, 4, 5)
    assert np.array_equal(np.swapaxes(v, -1, -2).reshape(-1), whole)
    z, y, x = 7, 3, 2
    assert v[z, y, x] == whole[y + 4 * (x + 5 * z)]
    assert np.array_equal(to_matlab(v).reshape(-1, order="F"), whole)
    # 4D: the leading axis is MATLAB's 4th dimension
    w = CounterRNG(5, "s").uint32((2, 5, 4, 5))
    assert np.array_equal(to_matlab(w).reshape(-1, order="F"), whole[:200])


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


# ------------------------------------------------------------------------------------------
# the generators in the port's form
# ------------------------------------------------------------------------------------------


def test_beads_match_the_formula():
    pos = [[6.5, 9.3, 11.0], [2, 19, 4]]  # (z, y, x), 0-based
    vol = synth.beads((12, 24, 18), pos, [2.0, 1.5], [1000, 300], 100)
    z, y, x = np.meshgrid(np.arange(12), np.arange(24), np.arange(18), indexing="ij")
    ref = np.full(vol.shape, 100, np.float32)
    for (pz, py, px), a in zip(pos, [1000, 300], strict=True):
        ref = ref + (
            a * np.exp(-((y - py) ** 2 + (x - px) ** 2) / 4.5 - (z - pz) ** 2 / 8.0)
        ).astype(np.float32)
    assert vol.dtype == np.float32 and vol.flags.c_contiguous and vol.shape == (12, 24, 18)
    assert np.max(np.abs(vol - ref)) <= 4 * np.spacing(np.float32(1100))
    assert np.array_equal(
        synth.beads((7, 5, 6), np.zeros((0, 3)), 1.0, 1, 3), np.full((7, 5, 6), 3)
    )
    # an image: (ny, nx) and (y, x) positions
    im = synth.beads((20, 30), [[9.0, 14.0]], 1.5, 100)
    assert im.shape == (20, 30) and np.unravel_index(np.argmax(im), im.shape) == (9, 14)
    assert np.array_equal(im, synth.beads((1, 20, 30), [[0, 9.0, 14.0]], 1.5, 100)[0])
    with pytest.raises(ValueError):
        synth.beads((5, 5, 5), [[1, 1, 1]], [1.0, 2.0, 3.0])


def test_centroids_find_the_beads():
    pos = np.array([[8.2, 9.3, 11.6], [13, 21, 29]])  # (z, y, x), 0-based
    vol = synth.beads((24, 32, 40), pos, [1.6, 1.2], 1000, 50)
    c = synth.centroids(vol, pos + 0.4, 5)
    assert np.max(np.abs(c - pos)) < 0.02
    im = synth.beads((30, 40), [[12.3, 20.7]], 1.3, 500, 10)
    assert np.max(np.abs(synth.centroid(im, [12, 21], 4) - [0, 12.3, 20.7])) < 0.02


def test_texture_spans_its_range():
    t = synth.texture((12, 16, 20), 2, (100, 2000))
    assert t.dtype == np.float32 and t.flags.c_contiguous and t.shape == (12, 16, 20)
    assert t.min() == 100 and t.max() == 2000
    assert np.array_equal(t, synth.texture((12, 16, 20), 2, (100, 2000), CounterRNG(0)))
    assert not np.array_equal(t, synth.texture((12, 16, 20), 2, (100, 2000), CounterRNG(1)))
    # smooth along every axis: neighbors are closer than independent values
    for ax in range(3):
        assert np.std(np.diff(t, axis=ax)) < 0.5 * np.std(t)
    t2 = synth.texture((30, 25), 2)
    assert t2.shape == (30, 25) and t2.min() == 100 and t2.max() == 2000


def test_gaussian_filter_keeps_constants_and_sums():
    c = np.full((5, 9, 7), 3.0)
    assert np.allclose(synth.gaussian_filter(c, 1.7), 3.0, rtol=1e-15)
    for shape, at in (((41, 1, 1), (20, 0, 0)), ((1, 41, 1), (0, 20, 0)), ((1, 1, 41), (0, 0, 20))):
        d = np.zeros(shape)
        d[at] = 1
        f = synth.gaussian_filter(d, 2.0)
        assert f.flags.c_contiguous
        assert math.isclose(math.fsum(f.ravel().tolist()), 1.0, rel_tol=1e-14)
        assert np.count_nonzero(f) == 2 * math.ceil(2 * 2.0) + 1  # the kernel size of imgaussfilt


def test_decon_psf_peak_and_zero_slices():
    p = synth.decon_psf()
    assert p.shape == (31, 31, 31) and p.dtype == np.float32 and p.flags.c_contiguous
    assert np.unravel_index(np.argmax(p), p.shape) == (15, 15, 15)
    assert p.max() == 10000
    assert not p[:5].any() and not p[-5:].any() and p[5].any()
    e = synth.decon_psf((24, 32, 30))
    assert np.unravel_index(np.argmax(e), e.shape) == (12, 16, 15)
    # sigma (z, yx): wider along z
    assert np.count_nonzero(p[:, 15, 15] > 100) > np.count_nonzero(p[15, :, 15] > 100)


def test_random_image_classes():
    for name, dt in (
        ("logical", np.bool_),
        ("uint8", np.uint8),
        ("uint16", np.uint16),
        ("int16", np.int16),
        ("single", np.float32),
        ("double", np.float64),
    ):
        im = synth.random_image((4, 6, 5), name)
        assert im.dtype == dt and im.shape == (4, 6, 5) and im.flags.c_contiguous
        if dt not in (np.bool_, np.float32, np.float64):
            assert im.min() >= 1
    assert np.array_equal(
        synth.random_image((6, 5), np.uint16), synth.random_image((6, 5), "uint16")
    )
    with pytest.raises(TypeError):
        synth.random_image((2, 2), "int64")


def test_results_do_not_depend_on_threads(thread_count):  # noqa: F811
    # NumPy's element-wise loops are single threaded; this pins that nothing else is used
    assert synth.array_hash(to_matlab(GENERATORS["texture3d"]())) == EXPECTED["texture3d"]


# ------------------------------------------------------------------------------------------
# files
# ------------------------------------------------------------------------------------------


def test_temp_folder_style(pk_temp_folder):  # noqa: F811
    assert pk_temp_folder.endswith("/") and "\\" not in pk_temp_folder


def test_tiff_is_written_as_cpp_tiff_reads_it(pk_temp_folder):  # noqa: F811
    cpptiff = pytest.importorskip("cpptiff")
    vol = synth.random_image((3, 7, 5), "uint16")
    synth.write_tiff(pk_temp_folder + "a.tif", vol)
    assert np.array_equal(cpptiff.read_tiff(pk_temp_folder + "a.tif"), vol)
    # MATLAB's form gives the same file
    synth.write_tiff_matlab(pk_temp_folder + "b.tif", to_matlab(vol))
    assert np.array_equal(cpptiff.read_tiff(pk_temp_folder + "b.tif"), vol)
    im = synth.random_image((7, 5), "uint8")
    synth.write_tiff(pk_temp_folder + "c.tif", im)
    assert np.array_equal(np.squeeze(cpptiff.read_tiff(pk_temp_folder + "c.tif")), im)


def _zarray(path: str) -> dict:
    import json

    with open(f"{path}/.zarray") as f:
        return json.load(f)


def test_zarr_layouts(write_volume, zarr_layout):  # noqa: F811
    cppzarr = pytest.importorskip("cppzarr")
    vol = synth.random_image((9, 40, 33), "uint16")
    p = write_volume("v", vol, "zarr", (8, 16, 32), layout=zarr_layout)
    meta = _zarray(p)
    back = cppzarr.read_zarr(p)
    if zarr_layout == "native":
        assert meta["shape"] == [9, 40, 33] and meta["order"] == "C"
        assert meta["chunks"] == [8, 16, 32]
        assert np.array_equal(back, vol)
    else:
        # MATLAB's writezarr: (y, x, z), F order, chunks (y, x, z)
        assert meta["shape"] == [40, 33, 9] and meta["order"] == "F"
        assert meta["chunks"] == [16, 32, 8]
        assert np.array_equal(back, to_matlab(vol))


def test_matlab_layout_zarr_is_the_old_file(pk_temp_folder):  # noqa: F811
    cppzarr = pytest.importorskip("cppzarr")
    vol = synth.random_image((5, 12, 10), "single")
    synth.write_zarr(pk_temp_folder + "a.zarr", vol, (4, 8, 8), layout="matlab")
    synth.write_zarr_matlab(pk_temp_folder + "b.zarr", to_matlab(vol), (8, 8, 4))
    assert _zarray(pk_temp_folder + "a.zarr") == _zarray(pk_temp_folder + "b.zarr")
    a = cppzarr.read_zarr(pk_temp_folder + "a.zarr")
    assert np.array_equal(a, cppzarr.read_zarr(pk_temp_folder + "b.zarr"))
    with pytest.raises(ValueError):
        synth.write_zarr(pk_temp_folder + "c.zarr", vol, layout="yxz")


def test_paired_inputs_hold_the_same_data(paired_input):  # noqa: F811
    cppzarr = pytest.importorskip("cppzarr")
    cpptiff = pytest.importorskip("cpptiff")
    vol = synth.texture((7, 30, 21), 2, (0, 4000)).astype(np.uint16)
    z = paired_input("tile", vol, "zarr", (4, 16, 16))
    assert z.matlab != z.python and z.matlab.endswith("/matlab/tile.zarr")
    assert z.python.endswith("/python/tile.zarr")
    assert _zarray(z.matlab)["order"] == "F" and _zarray(z.python)["order"] == "C"
    m = cppzarr.read_zarr(z.matlab)
    p = cppzarr.read_zarr(z.python)
    # the same data under the mapping, each in its side's layout
    assert np.array_equal(p, vol) and np.array_equal(m, to_matlab(vol))
    assert m.shape == (30, 21, 7) and p.shape == (7, 30, 21)
    t = paired_input("tile", vol, "tif")
    assert t.matlab == t.python and t.python.endswith("/tile.tif")
    assert np.array_equal(cpptiff.read_tiff(t.python), vol)
    with pytest.raises(ValueError):
        paired_input("x", vol, "h5")


def test_stitch_tiles_names_and_image_list(pk_temp_folder):  # noqa: F811
    pytest.importorskip("cpptiff")
    gt = synth.texture((10, 24, 60), 2, (100, 3000)).astype(np.uint16)
    # offsets and voxel sizes (z, y, x); stage coordinates and tile indices are x, y, z by name
    tiles = synth.cut_tiles(gt, [[0, 0, 0], [2, 4, 30]], (8, 20, 30), (0.3, 0.1, 0.108))
    assert [t["idx"] for t in tiles] == [[0, 0, 0], [1, 1, 1]]
    assert tiles[1]["xyz"] == [30 * 0.108, 4 * 0.1, 2 * 0.3]
    assert np.array_equal(tiles[1]["im"], gt[2:10, 4:24, 30:60])
    assert tiles[1]["im"].flags.c_contiguous
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

    assert np.array_equal(cpptiff.read_tiff(f"{pk_temp_folder}{names[1]}.tif"), tiles[1]["im"])
    # 2D tiles: offsets and sizes (y, x)
    im = gt[0]
    t2 = synth.cut_tiles(im, [[0, 0], [4, 30]], (20, 30), (0.1, 0.108))
    assert np.array_equal(t2[1]["im"], im[4:24, 30:60]) and t2[1]["xyz"] == [30 * 0.108, 0.4, 0.0]


def test_zarr_tiles(pk_temp_folder, zarr_layout):  # noqa: F811
    cppzarr = pytest.importorskip("cppzarr")
    vol = synth.random_image((9, 40, 33), "uint16")
    names, _ = synth.write_stitch_tiles(
        pk_temp_folder,
        [{"im": vol, "idx": [0, 0, 0], "xyz": [0, 0, 0]}],
        save_format="zarr",
        zarr_layout=zarr_layout,
    )
    back = cppzarr.read_zarr(f"{pk_temp_folder}{names[0]}.zarr")
    assert np.array_equal(back, vol if zarr_layout == "native" else to_matlab(vol))


if __name__ == "__main__":
    # --old <synth.py before the change>: the hash table of the parameter sets (paste into
    # EXPECTED_SETS); without arguments: the hash table of GENERATORS
    if len(sys.argv) == 3 and sys.argv[1] == "--old":
        import importlib.util

        spec_ = importlib.util.spec_from_file_location("synth_old", sys.argv[2])
        old_ = importlib.util.module_from_spec(spec_)
        sys.modules["synth_old"] = old_
        spec_.loader.exec_module(old_)
        for n, s in PARAM_SETS:
            print(f'    "{n}": "{old_.array_hash(run_old(old_, s))}",')
    else:
        for k, f in GENERATORS.items():
            print(f'    "{k}": "{synth.array_hash(to_matlab(f()))}",')
    sys.exit(0)
