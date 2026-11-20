"""3D point source detection, tracking and lifetime analysis (Phase 3).

New implementations from the published methods (Aguet et al. 2013 for the model-based detection,
Jaqaman et al. 2008 for the tracking) and from MATLAB's behaviour recorded as black-box
references; no code of the GPL sources is used. ``docs/detection_scope.md`` lists the stages,
their inputs and outputs, and the references.

Modules:

- ``stats``: MATLAB's statistics functions as the detection code calls them (norminv, tcdf, tinv,
  chi2inv, fcdf, ksdensity, ecdf, mad, kstest2).
"""
