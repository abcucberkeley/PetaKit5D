"""Inputs of the xcorr_k06 benchmark: template and image pairs as tiff files in the work folder.

- ovl300: an image of 300 x 300 x 200 (MATLAB's y x z) of smooth texture and a template of
  260 x 260 x 180 cut from it, shifts within +-20 (y, x) and +-10 (z) around the template's place;
- demo2p: the largest pair of the 2-photon stitching demo's registration (template 819 x 819 x
  11, image 1024 x 1024 x 32, a window of 401 x 401 x 21).
"""

import sys
from pathlib import Path

from petakit5d.testing.synth import CounterRNG, texture, write_tiff

PAIRS = {
    # name: image (z, y, x), template (z, y, x), the template's place in the image (z, y, x)
    "ovl300": ((200, 300, 300), (180, 260, 260), (10, 20, 20)),
    "demo2p": ((32, 1024, 1024), (11, 819, 819), (10, 102, 102)),
}


def main(work: str) -> None:
    for name, (sa, st, at) in PAIRS.items():
        im = texture(sa, 2.0, (100.0, 2000.0), rng=CounterRNG(11, f"bench/{name}"))
        t = im[tuple(slice(o, o + n) for o, n in zip(at, st, strict=True))].copy()
        write_tiff(Path(work) / f"{name}_A.tif", im)
        write_tiff(Path(work) / f"{name}_T.tif", t)


if __name__ == "__main__":
    main(sys.argv[1])
