"""The user function of the generic computing framework demo (MATLAB's
``demos/crop_deskew_rotate_demo_function.m``): crop a tile, deskew and rotate it, and write the
result under a temporary name that is renamed at the end, so that the framework sees either no
result or a complete one.

The framework runs it as a task (``Task('crop_deskew_rotate_demo_function:
crop_deskew_rotate_demo_function', args, path=[this folder])``), in this process or in a job.
"""

from __future__ import annotations

import os
from collections.abc import Sequence


def crop_deskew_rotate_demo_function(
    input_fullpath: str,
    output_fullpath: str,
    skew_angle: float,
    xy_pixel_size: float,
    dz: float,
    in_bbox: Sequence[int],
) -> None:
    """crops ``input_fullpath`` to ``in_bbox`` ``(z0, y0, x0, z1, y1, x1)`` (0-based,
    end-exclusive), deskews and rotates it and writes ``output_fullpath``"""
    from petakit5d.core.ids import get_uuid
    from petakit5d.deskew_rotate import deskew_rotate_volume
    from petakit5d.io import read_tiff, write_tiff
    from petakit5d.kernels.crop import crop_3d

    if not os.path.exists(input_fullpath):
        raise FileNotFoundError(f"The input file {input_fullpath} does not exist!")
    if os.path.exists(output_fullpath):
        # as in MATLAB's function, the message is printed and the file is made again
        print(f"The output file {output_fullpath} already exist, skip it!")

    im = read_tiff(input_fullpath)
    im = crop_3d(im, in_bbox)

    # the demo images' settings; the deskew/rotation runs in single, input and output are 16 bit
    dsr = deskew_rotate_volume(
        im,
        skew_angle,
        dz,
        xy_pixel_size,
        reverse=True,
        crop=True,
        objective_scan=False,
        resample_factor=None,
        interp_method="linear",
        save_16bit=True,
    )

    # write to a temporary file, then rename it to the final name
    tmp_fn = f"{output_fullpath[:-4]}_{get_uuid()}.tif"
    write_tiff(dsr, tmp_fn)
    os.replace(tmp_fn, output_fullpath)
