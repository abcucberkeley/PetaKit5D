"""The unmixing kernel benchmark (B02), Python side (MATLAB: bench_unmix_sum_b02.m): in a process
pinned to the CPUs of ``<work>/cpus_<t>.txt``, the sum of one demo tile pair in memory
(``unmix_channels_sum``: two uint16 channels, [1 -0.1], backgrounds [100 100], offset 1, uint16),
3 timed calls after one untimed: ``sum`` into a new output each call (its free not timed),
``sum_out`` into one output reused.

usage: python bench_unmix_sum_b02.py <work folder> <threads>
"""

import os
import sys

WORK, T = sys.argv[1], int(sys.argv[2])
with open(os.path.join(WORK, f"cpus_{T}.txt")) as f:
    os.sched_setaffinity(0, [int(c) for c in f.read().split()])

from petakit5d.io import read_tiff  # noqa: E402
from petakit5d.testing.bench import mark, record, timeit  # noqa: E402
from petakit5d.tools import unmix  # noqa: E402

mark("ready")
TILE = (
    "Scan_Iter_0000_0000_{}_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs"
    "_000x_003y_000z_0000t"
)
ims = [read_tiff(os.path.join(WORK, "demo", TILE.format(c) + ".tif")) for c in ("CamA", "CamB")]
mark("loaded")
args = (ims, (1, -0.1), (100, 100), 1, "uint16")
record("sum", timeit(lambda: unmix.unmix_channels_sum(*args), repeats=3), threads=T)
out = unmix.unmix_channels_sum(*args)
record("sum_out", timeit(lambda: unmix.unmix_channels_sum(*args, out=out), repeats=3), threads=T)
