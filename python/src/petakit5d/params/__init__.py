"""Parameter registry: one declaration per user-facing function (PLAN.md section 3.3).

    from petakit5d.params import Param, register

    @register("XR_crop_dataset", Param("dataPaths", "@(x) ischar(x) || iscell(x)", "..."), ...)
    def crop_dataset(data_paths, input_bbox, *, result_dir_name="Cropped", ...): ...

gives the snake_case API (the decorated function), ``petakit5d.matlab.XR_crop_dataset``
(MATLAB names in any case), the command-line conversion (``FunctionSpec.bind_strings``),
the saved parameters struct (``matlab_results``) and the mapping table. The spec format is
documented in docs/params.md. Importing this package loads neither NumPy nor any area.
"""

from petakit5d.params.literal import LiteralError, format_literal, matlab_eval, parse, str2num
from petakit5d.params.registry import (
    REGISTRY,
    NotPorted,
    Registry,
    find,
    matlab_results,
    not_ported,
    register,
    spec_of,
)
from petakit5d.params.spec import (
    FunctionSpec,
    MatlabArgumentError,
    Param,
    ParamSpec,
    RegistryError,
    snake_name,
)
from petakit5d.params.values import Cell, MArray, MString, from_python, to_python

__all__ = [
    "REGISTRY",
    "Cell",
    "FunctionSpec",
    "LiteralError",
    "MArray",
    "MString",
    "MatlabArgumentError",
    "NotPorted",
    "Param",
    "ParamSpec",
    "Registry",
    "RegistryError",
    "find",
    "format_literal",
    "from_python",
    "matlab_eval",
    "matlab_results",
    "not_ported",
    "parse",
    "register",
    "snake_name",
    "spec_of",
    "str2num",
    "to_python",
]
