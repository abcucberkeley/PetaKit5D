"""N06 benchmark setup: 900-tile assignment problems (the test module's generator) as MAT-files.

usage: n06_lsq_assign_setup.py <python/ tree> <work folder>
"""

import importlib.util
import sys

import numpy as np
import scipy.io as sio

src, work = sys.argv[1], sys.argv[2]
spec = importlib.util.spec_from_file_location("test_lsq", f"{src}/tests/matlabcompat/test_lsq.py")
tl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tl)
problems = {
    "grid": {"filt": 0.05},
    "free": {"filt": 0.1, "free_tiles": 12},
    "tight": {"filt": 0.05, "max_off": (6, 6, 3), "off": (12, 12, 4)},
}
for name, kw in problems.items():
    p = tl.grid_problem(1, 30, 30, **kw)
    ov = np.zeros((int(p["nF"]),) * 2, bool)
    pr = p["pairs"].astype(int) - 1
    ov[pr[:, 0], pr[:, 1]] = True
    sio.savemat(f"{work}/n06_p900_{name}.mat", {**p, "ov": ov})
print("problems written")
