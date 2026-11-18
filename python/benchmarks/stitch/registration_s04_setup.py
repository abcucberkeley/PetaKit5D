"""Inputs of the registration_s04 benchmark, in the work folder.

- pair300: two uint16 tiles of 300 x 300 x 200 (MATLAB's y x z) of smooth texture, the second 20
  voxels on in x and listed 3, -2, 1 voxels off its position; zarr in each side's own layout
  (``matlab/`` and ``python/``, chunks 256), the cuboids and the overlap in ``pair300.json``;
- demo (with DEMO, a copy of the demo cell dataset stitched by s04_demo.m and the port): the
  stitching frame's 12 registration tasks of its 4 volumes of 4 tiles (``calls.json`` of the
  port's run), as MATLAB's calls (``demo_matlab.json``: 1-based pairs, data_order_mat, downSample
  [y x z], MaxOffset [x y z], MATLAB's zarr tiles) and the port's (``demo_python.json``: native
  copies of the tiles in ``native_zarr/``).

usage: registration_s04_setup.py WORK [DEMO]
"""

import json
import sys
from pathlib import Path

import numpy as np

from petakit5d.core.dtypes import cast
from petakit5d.stitch.geometry import check_cuboids_overlaps, data_order_axes_to_matlab
from petakit5d.testing.synth import CounterRNG, texture, write_paired_inputs


def pair300(work: Path) -> None:
    g = texture((200, 300, 320), 2.0, (100.0, 3000.0), rng=CounterRNG(12, "bench/pair300"))
    g = np.asarray(cast(g.astype(np.float64), "uint16"))
    write_paired_inputs(str(work), "t1", np.ascontiguousarray(g[:, :, :300]), "zarr", (256,) * 3)
    write_paired_inputs(str(work), "t2", np.ascontiguousarray(g[:, :, 20:]), "zarr", (256,) * 3)
    c1 = np.array([0.0, 0.0, 0.0, 299.0, 299.0, 199.0])
    p = np.array([23.0, -2.0, 1.0])
    c2 = np.concatenate([p, p + [299.0, 299.0, 199.0]])
    _, ov = check_cuboids_overlaps(c1, c2)
    spec = {"c1": c1.tolist(), "c2": c2.tolist(), "co": ov.ravel().tolist(), "vs": [1, 1, 1],
            "dom": [2, 1, 3], "ds": [1, 1, 1], "mo": [20, 20, 10]}  # fmt: skip
    (work / "pair300.json").write_text(json.dumps(spec), newline="\n")


def demo(work: Path, data: Path) -> None:
    calls = json.loads((data / "python_stitch_exact" / "calls.json").read_text())
    ml, py = [], []
    for k, c in enumerate(calls):
        a, kw = c["args"], c["kwargs"]
        out = f"xcorr_tile_{k + 1}.mat"
        ml.append({
            "img1": a[0], "img2": a[1], "out": out,
            "pairs": (np.asarray(a[3]) + 1).tolist(), "c1": a[4], "c2": a[5], "co": a[6],
            "vs": a[7], "dom": list(data_order_axes_to_matlab(a[8])),
            "ds": [kw["down_sample"][1], kw["down_sample"][2], kw["down_sample"][0]],
            "mo": kw["max_offset"][::-1],
        })  # fmt: skip
        b = list(a)
        b[0] = a[0].replace("/zarr/", "/native_zarr/")
        b[1] = [p.replace("/zarr/", "/native_zarr/") for p in a[1]]
        b[2] = out
        py.append({"args": b, "kwargs": {**kw, "parse_cluster": False}})
    (work / "demo_matlab.json").write_text(json.dumps(ml), newline="\n")
    (work / "demo_python.json").write_text(json.dumps(py), newline="\n")


if __name__ == "__main__":
    w = Path(sys.argv[1])
    pair300(w)
    if len(sys.argv) > 2 and sys.argv[2]:
        demo(w, Path(sys.argv[2]))
