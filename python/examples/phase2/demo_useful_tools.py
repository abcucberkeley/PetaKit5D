"""The Phase 2 steps of the data tools demo (MATLAB: ``demos/demo_useful_tools.m``): the FFT
power spectra, the PSF analysis, the Fourier shell correlation (FSC) analysis and the chromatic
shift correction, in the demo's order; one function per step, with the demo's parameters.

The demo's other steps are ``examples/phase1/demo_useful_tools.py``: run it first, on the same
data (``PETAKIT5D_DEMO_DATA``). The FFT spectra take its deskewed decon results
(``matlab_decon_omw/DS_for_useful_tools/``) and the FSC analysis its decon with intermediate
results (``matlab_decon_omw_debug/``). These steps read nothing that the later steps of either
example write, so running them after the whole phase 1 example gives MATLAB's results.

Left out: the Imaris converter (``XR_imaris_conversion_data_wrapper`` is not ported yet: it waits
for the next release of the converter, Parallel_Imaris_Writer). The Bio-Formats conversion,
which MATLAB's demo shows after its ``return`` (no demo images for it), is
``bioformats_conversion`` below, not run, as in MATLAB. The PSF analysis needs matplotlib (the
``figures`` extra: its results are figures) and is left out without it; the FSC analysis then
writes no figures. Results: ``matlab_decon_omw/DS_for_useful_tools/FFT/``,
``PSF/DS/PSFAnalysis/``, ``matlab_decon_omw_debug/<tile>_debug/FSCs/`` (and ``figures/``) and
``Chromatic_Shift_Corrected/`` in the dataset folder.
"""

from __future__ import annotations

import os
import sys

import _phase1  # noqa: F401  (examples/phase1 on the import path)

# isort: split
import _common

DEMO = "demo_useful_tools"


def dz_deskewed() -> float:
    """the z step of the deskewed data: dz * sind(skew angle), with MATLAB's sind"""
    from petakit5d.matlabcompat.elementary import sind

    return 0.3 * float(sind(32.45))


def fft_analysis(data_path: str) -> None:
    """the power spectra of the deskewed decon results (zarr): their MIPs and central slices
    with gamma 0.5 (no 3D spectra), in ``DS_for_useful_tools/FFT/``"""
    from petakit5d.tools import fft_spectrum_dataset

    fft_spectrum_dataset(
        [data_path + "matlab_decon_omw/DS_for_useful_tools/"],
        result_dir_name="FFT",
        xy_pixel_size=0.108,
        dz=dz_deskewed(),
        zarr_file=True,
        out_size=(1001, 1001, 1001),
        # isotropic, a tenth of the wavelength in the medium
        out_pixel_size=0.1 * 0.488 / 1.33,
        save_3d_stack=False,
        parse_cluster=False,
        master_compute=True,
        cpus_per_task=24,
        mcc_mode=False,
        config_file="",
    )


def psf_analysis(data_path: str) -> None:
    """the PSF in skewed space deskewed and analyzed in real and Fourier space, against the
    theoretical widefield PSF of Richards & Wolf, in ``PSF/DS/PSFAnalysis/``"""
    from petakit5d.tools import psf_analysis_dataset

    psf_analysis_dataset(
        [data_path + "PSF/"],
        xy_pixel_size=0.108,
        dz=0.3,
        skew_angle=32.45,
        channel_patterns=["488"],
        channels=[488],
        deskew=True,
        objective_scan=False,
        source_str="Demo PSF",
        # the Richards & Wolf PSF of the channel
        rw_fn=[data_path + "PSF/RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif"],
        visible=True,
        parse_cluster=False,
        master_compute=True,
        cpus_per_task=24,
        mcc_mode=False,
        config_file="",
    )


def fsc_analysis(data_path: str) -> None:
    """the FSC of the decon's intermediate results (every 2 iterations), to choose the number
    of iterations, in ``<tile>_debug/FSCs/`` and ``figures/``"""
    import math

    from petakit5d.tools import fsc_analysis_dataset

    fsc_analysis_dataset(
        [data_path + f"matlab_decon_omw_debug/{_common.TILE_B}_debug/"],
        result_dir_name="FSCs",
        xy_pixel_size=0.108,
        dz=dz_deskewed(),
        dr=5,
        dtheta=math.pi / 6,
        res_thresh_method="one-bit",
        res_thresh=0.2,
        half_size=(201, 201, 201),
        input_bbox=None,
        channel_patterns=["Iter"],
        channels=[488],
        suffix="decon",
        iter_interval=2,
        parse_cluster=False,
        master_compute=True,
        cpus_per_task=24,
        mcc_mode=False,
        config_file="",
    )


def chromatic_shift_correction(data_path: str) -> None:
    """the CamB tiles shifted onto the CamA tiles by a given offset, both cropped to their
    overlap ('valid'), in ``Chromatic_Shift_Corrected/``"""
    from petakit5d.tools import chromatic_shift_correction_dataset

    chromatic_shift_correction_dataset(
        [data_path],
        # one row (z, y, x) per channel: none for CamA, CamB by +1 in y, -5 in x and -2 in z
        # (MATLAB's [0, 0, 0; 1, -5, -2] in [y x z])
        chromatic_offset=[(0, 0, 0), (-2, 1, -5)],
        result_dir_name="Chromatic_Shift_Corrected",
        mode="valid",
        pad_value=0,
        new_origin=None,
        channel_patterns=["CamA", "CamB"],
        psf_fullpaths=[],
        max_offset=(20, 20, 20),
        crop_length=(0, 0, 0),
        zarr_file=False,
        large_file=False,
        save_zarr=False,
        batch_size=(1024, 1024, 1024),
        block_size=(256, 256, 256),
        parse_cluster=True,
        master_compute=True,
        job_log_dir="../job_logs",
        cpus_per_task=3,
        config_file="",
        mcc_mode=False,
        uuid="",
        debug=False,
    )


def bioformats_conversion() -> None:
    """MATLAB's last section, after the demo's ``return`` (for illustration: the demo data has no
    such files): ND2 and CZI files to one tiff per series and channel. The port reads them
    without Bio-Formats, with the readers of the ``formats`` extra."""
    from petakit5d.io.bioformats import bioformats_to_tiff_or_zarr_dataset

    bioformats_to_tiff_or_zarr_dataset(
        ["/data/microscopy/sample_dataset_1/", "/data/microscopy/sample_dataset_2/"],
        result_dir_name="tiffs",
        channel_patterns=[".nd2", ".czi"],
        n_channels=2,
        data_formats=[".nd2", ".czi"],
        save_zarr=False,
        block_size=(256, 256, 256),
        over_write=False,
        uuid="",
    )


def main() -> None:
    steps = _common.Steps(DEMO)
    data_path = _common.data_path()
    needed = (
        data_path + "matlab_decon_omw/DS_for_useful_tools/",
        data_path + f"matlab_decon_omw_debug/{_common.TILE_B}_debug/",
    )
    missing = [p for p in needed if not os.path.isdir(p)]
    if missing:
        sys.exit(
            f"{missing[0]} is missing: run examples/phase1/demo_useful_tools.py first, on the "
            "same data (these steps use its decon and deskew results)."
        )

    with steps("FFT analysis"):
        fft_analysis(data_path)
    if _common.figures_available():
        with steps("PSF analysis"):
            psf_analysis(data_path)
    else:
        steps.skip("PSF analysis", "matplotlib is not installed (the figures extra)")
    with steps("FSC analysis"):
        fsc_analysis(data_path)
    steps.skip(
        "Imaris converter",
        "XR_imaris_conversion_data_wrapper is not ported yet: it waits for the next release of "
        "Parallel_Imaris_Writer",
    )
    with steps("chromatic shift correction"):
        chromatic_shift_correction(data_path)
    steps.skip("Bio-Formats conversion", "after the demo's return: bioformats_conversion()")


if __name__ == "__main__":
    _common.run_demo(main)
