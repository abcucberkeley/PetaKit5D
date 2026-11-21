"""Lifetime analysis benchmark (C05), the inputs: the demo condition as the lifetime analysis reads
it, once for each side (each writes lifetimeData.mat and Parameters.mat next to its tracks):

    WORK/<side>/<condition>/<experiment>/ch<k>/        links to the frames (ch1, ch2 and DS)
    WORK/<side>/<condition>/<experiment>/ch1/<analysis>/ProcessedTracks.mat   a copy

FRAMES: the condition folder of the frames (the demo dataset); TRACKS: the experiments of a
MATLAB run of the demo with their ch1/<analysis>/ProcessedTracks.mat (and DS frames, linked
where they exist), as in PetaKit5D_point_detection_demo_tracks.

usage: python lifetime_c05_data.py WORK FRAMES TRACKS ANALYSIS
"""

from __future__ import annotations

import os
import shutil
import sys


def main() -> None:
    work, frames, tracks, analysis = sys.argv[1:5]
    cond = os.path.basename(os.path.normpath(tracks))
    for side in ("matlab", "python"):
        for ex in sorted(os.listdir(tracks)):
            src = os.path.join(tracks, ex)
            if not os.path.isdir(src) or not ex.startswith("Ex"):
                continue
            for ch in sorted(c for c in os.listdir(os.path.join(frames, ex)) if c.startswith("ch")):
                dst = os.path.join(work, side, cond, ex, ch)
                os.makedirs(dst, exist_ok=True)
                fdir = os.path.join(frames, ex, ch)
                for f in sorted(os.listdir(fdir)):
                    if f.endswith(".tif"):
                        os.symlink(os.path.join(fdir, f), os.path.join(dst, f))
                ds = os.path.join(src, ch, "DS")
                if os.path.isdir(ds):
                    os.symlink(ds, os.path.join(dst, "DS"))
            a = os.path.join(work, side, cond, ex, "ch1", analysis)
            os.makedirs(a, exist_ok=True)
            shutil.copyfile(
                os.path.join(src, "ch1", analysis, "ProcessedTracks.mat"),
                os.path.join(a, "ProcessedTracks.mat"),
            )
    print(f"lifetime_c05 inputs: {cond} for matlab and python in {work}")


if __name__ == "__main__":
    main()
