"""W03 trial: every kernel case's digest, and the Anderson-Darling cases' outputs as hex, for the
comparison with Linux's (run from python/)

    python ../.github/w03_trial/dump_cases.py
"""

import sys

import numpy as np

sys.path.insert(0, "tests/platform")
import _kernel_cases as kc  # noqa: E402

for name in sorted(kc.CASES):
    print("DIGEST", name, kc.digest(kc.run(name)), flush=True)
for name in ("detection.ad_stat", "detection.fit_mex.ad"):
    for path, a in kc.leaves(kc.run(name)):
        a = np.ascontiguousarray(a)
        shape = "x".join(str(n) for n in a.shape) or "-"
        print("LEAF", name + path, a.dtype.str, shape, a.tobytes().hex(), flush=True)
