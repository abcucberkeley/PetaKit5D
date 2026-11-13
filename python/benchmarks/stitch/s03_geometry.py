"""S03 benchmark: the stitching geometry for an n x n grid of DSR tiles (the layout of
ml/s03_geometry_ml.m): prepare_layout and place_tiles, with zero registration shifts and
feather blending (overlap type 'half')."""

import argparse

import numpy as np

from petakit5d.testing.bench import mark, record, timeit

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=40)
ap.add_argument("--repeats", type=int, default=5)
a = ap.parse_args()

from petakit5d.stitch import geometry as g  # noqa: E402

n = a.n
iy, ix = np.divmod(np.arange(n * n, dtype=np.float64), n)
coords = np.column_stack(
    (ix * 48.6 + np.mod(ix * 7 + iy * 3, 11) / 1000, iy * 48.6 + np.mod(ix * 5 + iy * 2, 13) / 1000,
     np.zeros(n * n))
)  # fmt: skip
sizes = np.tile([200.0, 500.0, 500.0], (n * n, 1))
tidx = np.column_stack((ix, iy, np.zeros(n * n), np.zeros(n * n)))
mark("ready")


def run():
    lay = g.prepare_layout(coords, sizes, tile_idx=tidx, dsr=True)
    return g.place_tiles(lay, lay.xyz + np.zeros((n * n, 3)) * lay.space.pixel_sizes,
                         blend_method="feather", xcorr_shift=True)  # fmt: skip


pl = run()
print("nv_size", pl.nv_size, "pairs", pl.pairs.shape[0], "batches", pl.num_batches)
record(f"geometry_{n * n}", timeit(run, repeats=a.repeats, warmup=1), threads=1)
