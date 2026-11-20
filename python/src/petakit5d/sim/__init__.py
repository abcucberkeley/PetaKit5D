"""Structured illumination microscopy (SIM) reconstruction: MATLAB's ``third_parties/sim3D``.

========================================  ==========================================
``sim_recon_volume``                      simRecon
``sim_psf_to_otf``                        sim_PSFtoOTF_gen
``gaussian_fit_3d``                       Lsq_GaussianFit_3D
``sim_recon_file``                        simReconFrame
``deskew_phases_file``                    deskewPhasesFrame
``separate_phases_file``                  separatePhases
``combine_phases_file``                   combinePhases
``sim_recon_dataset``                     simRecon_data_wrapper
``deskew_phases_dataset``                 deskewPhases_data_wrapper
``rotate_3d_dataset``                     Rotate3D_data_wrapper
``sim_recon_automatic_processing``        simReconAutomaticProcessing
``perdecomp_3d``                          perdecomp_3D (from the published method)
``normxcorr2_masked``                     normxcorr2_masked (from the published method)
``masked_translation_registration``       MaskedTranslationRegistration2D_fit
``imtranslate``                           imtranslate_function (from its behaviour)
``fminsearch``                            MATLAB's fminsearch (Nelder-Mead)
========================================  ==========================================

The steps of the reconstruction (``separation_matrix``, ``resample_otf``, ``cyl_mask``,
``tukey_window``, ``fourier_shift``, ``apodize_ellipse``) are in ``sim.recon``. Nothing is
imported until a name is first used.
"""

_EXPORTS = {
    "SimError": "recon",
    "apodize_ellipse": "recon",
    "cyl_mask": "recon",
    "fourier_shift": "recon",
    "inverse_separation_matrix": "recon",
    "resample_otf": "recon",
    "separation_matrix": "recon",
    "sim_recon_volume": "recon",
    "tukey_window": "recon",
    "gaussian_fit_3d": "otf",
    "sim_psf_to_otf": "otf",
    "combine_phases_file": "file",
    "deskew_phases_file": "file",
    "separate_phases_file": "file",
    "sim_recon_file": "file",
    "deskew_phases_dataset": "dataset",
    "rotate_3d_dataset": "dataset",
    "sim_recon_dataset": "dataset",
    "sim_recon_automatic_processing": "automatic",
    "perdecomp_3d": "perdecomp",
    "masked_translation_registration": "masked_ncc",
    "normxcorr2_masked": "masked_ncc",
    "imtranslate": "translate",
    "FminsearchResult": "neldermead",
    "fminsearch": "neldermead",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    module = _EXPORTS.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    value = getattr(importlib.import_module(f"{__name__}.{module}"), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(_EXPORTS))
