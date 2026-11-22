"""Track processing benchmark (C04b), the Python side: run_track_processing of FOLDER/Analysis
(the tracker's and the detector's files, the masks) for the movie of DATAMAT ('data', whose DS
frames are read in place), then rotate_tracks; the results removed before each run
(untimed).

usage: python processing_c04b.py FOLDER DATAMAT --repeats 1
"""

from __future__ import annotations

import argparse
import os
import shutil
import time

from petakit5d.testing.bench import mark, record

mark("ready")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("datamat")
    ap.add_argument("--repeats", type=int, default=1)
    a = ap.parse_args()
    from petakit5d.detection.processing import rotate_tracks, run_track_processing
    from petakit5d.io.matfile import read_mat

    data = read_mat(a.datamat, "data")["data"]
    data = data if isinstance(data, dict) else list(data)[0]
    rp = os.path.join(a.folder, "Analysis") + os.sep
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
    t_proc, t_rot = [], []
    for _ in range(a.repeats):
        for f in ("ProcessedTracks.mat", "RotatedTracks.mat"):
            if os.path.exists(rp + f):
                os.remove(rp + f)
        if os.path.isdir(rp + "Amira"):
            shutil.rmtree(rp + "Amira")
        t0 = time.perf_counter()
        run_track_processing(data, overwrite=True, results_path=[rp])
        t1 = time.perf_counter()
        rotate_tracks(data, overwrite=True, results_path=[rp], crop=False)
        t2 = time.perf_counter()
        t_proc.append(t1 - t0)
        t_rot.append(t2 - t1)
    record("process", t_proc, threads=cpus)
    record("rotate", t_rot, threads=cpus)


if __name__ == "__main__":
    main()
