"""Specs: one declaration, three ways in (snake_case, MATLAB names, command-line strings)."""

from __future__ import annotations

import pytest

from petakit5d.params import (
    Cell,
    MArray,
    MatlabArgumentError,
    MString,
    Param,
    Registry,
    RegistryError,
    snake_name,
)
from petakit5d.params.spec import validator_alternatives
from petakit5d.params.values import describe

ROOT = "/mnt/projects/embryo_2024"


@pytest.mark.parametrize(
    ("matlab", "py"),
    [
        ("xyPixelSize", "xy_pixel_size"),
        ("DSRDirName", "dsr_dir_name"),
        ("save16bit", "save_16bit"),
        ("save3DStack", "save_3d_stack"),
        ("parseCluster", "parse_cluster"),
        ("FFCorrection", "ff_correction"),
        ("FFImagePaths", "ff_image_paths"),
        ("GPUJob", "gpu_job"),
        ("useGPU", "use_gpu"),
        ("OTFCumThresh", "otf_cum_thresh"),
        ("dz", "dz"),
        ("dzPSF", "dz_psf"),
        ("RLMethod", "rl_method"),
        ("Overwrite", "overwrite"),
        ("zarrSubSize", "zarr_sub_size"),
        ("BKRemoval", "bk_removal"),
        ("IOScan", "io_scan"),
        ("stitchMIP", "stitch_mip"),
        ("tileOutBbox", "tile_out_bbox"),
        ("xcorrDownsample", "xcorr_downsample"),
        ("Decon", "decon"),
        ("XR_crop_dataset", "xr_crop_dataset"),
    ],
)
def test_snake_name(matlab, py):
    assert snake_name(matlab) == py


# -- the worked example: XR_crop_dataset three ways ------------------------------------------


def test_three_ways_give_the_same_call(sample_crop):
    spec = sample_crop.REG.find("XR_crop_dataset")
    box = [1, 1, 1, 100, 200, 50]
    a = sample_crop.crop_dataset(
        ROOT, box, batch_size=[512, 512, 512], pad=True, channel_patterns=["CamA_ch0"]
    )
    b = spec.call_matlab(
        ROOT, box, "BatchSize", [512, 512, 512], "pad", True, ChannelPatterns=["CamA_ch0"]
    )
    c = spec.call_strings(
        [ROOT, "[1,1,1,100,200,50]", "batchsize", "[512,512,512]", "Pad", "true",
         "channelPatterns", "{'CamA_ch0'}"]
    )  # fmt: skip
    import petakit5d.matlab as pm

    assert a == b == c
    assert a["batch_size"] == (512, 512, 512) and a["input_bbox"] == (1, 1, 1, 100, 200, 50)
    assert a["channel_patterns"] == ("CamA_ch0",)
    assert a["result_dir_name"] == "Cropped"
    # values the function gets are canonical: tuples, ints, str
    assert type(c["cpus_per_task"]) is int and type(c["data_paths"]) is str
    with pytest.raises(AttributeError):
        pm.XR_crop_dataset_not_there  # noqa: B018


def test_unknown_and_misspelled_names_are_errors_with_a_suggestion(sample_crop):
    spec = sample_crop.REG.find("XR_crop_dataset")
    with pytest.raises(TypeError, match=r"batch_sise.*Did you mean 'batch_size'\?"):
        sample_crop.crop_dataset(ROOT, [1] * 6, batch_sise=1)
    with pytest.raises(TypeError, match=r"Did you mean 'result_dir_name'\?"):
        sample_crop.crop_dataset(ROOT, [1] * 6, resultDirName="x")
    with pytest.raises(MatlabArgumentError, match=r"'bachSize' is not a recognized.*'batchSize'"):
        spec.call_matlab(ROOT, [1] * 6, "bachSize", [1, 1, 1])
    with pytest.raises(MatlabArgumentError, match="Did you mean 'parseCluster'") as e:
        spec.call_strings([ROOT, "[1]", "parseClustr", "true"])
    assert e.value.identifier == "MATLAB:InputParser:UnmatchedParameter"
    with pytest.raises(MatlabArgumentError) as e:
        spec.call_matlab(ROOT, [1] * 6, "save", True)
    assert e.value.identifier == "MATLAB:InputParser:AmbiguousParameter"
    assert "'save16bit', 'saveZarr'" in str(e.value)


def test_values_are_checked_as_matlab_checks_them(sample_crop):
    spec = sample_crop.REG.find("XR_crop_dataset")
    for bad in (
        lambda: sample_crop.crop_dataset(ROOT, [1] * 6, pad=1),
        lambda: spec.call_matlab(ROOT, [1] * 6, "pad", 1),
        lambda: spec.call_strings([ROOT, "[1]", "pad", "yes"]),
        lambda: spec.call_strings([ROOT, "[1]", "pad", "1"]),
        lambda: spec.call_matlab(ROOT, [1] * 6, "channelPatterns", "CamA_ch0"),
        lambda: spec.call_strings([ROOT, "[1]", "channelPatterns", "CamA_ch0"]),
        lambda: spec.call_matlab(ROOT, "x"),
        lambda: spec.call_matlab(ROOT, [1] * 6, "cpusPerTask", True),
    ):
        with pytest.raises(MatlabArgumentError) as e:
            bad()
        assert e.value.identifier == "MATLAB:InputParser:ArgumentFailedValidation"
    with pytest.raises(MatlabArgumentError, match="must satisfy the function: islogical"):
        spec.call_matlab(ROOT, [1] * 6, "pad", "true")


def test_string_conversions_follow_the_mcc_parser(sample_crop):
    spec = sample_crop.REG.find("XR_crop_dataset")
    got = spec.bind_strings(
        [ROOT, "1 1 1 10 10 10", "lastStartCoords", "[]", "uuid", "abc", "cropType", "{x}"]
    )
    assert got["input_bbox"] == (1, 1, 1, 10, 10, 10)
    assert got["last_start_coords"] is None
    assert got["uuid"] == "abc" and got["crop_type"] == "{x}"  # char: never converted
    got = spec.bind_strings(["{'/a','/b'}", "[1,2,3,4,5,6]"])
    assert got["data_paths"] == ("/a", "/b")
    got = spec.bind_strings(["[/a]", "[1,2,3,4,5,6]"])
    assert got["data_paths"] == "[/a]"  # ischar || iscell: only a leading { is evaluated
    with pytest.raises(MatlabArgumentError) as e:
        spec.bind_strings([ROOT])
    assert e.value.identifier == "MATLAB:badsubscript"
    with pytest.raises(MatlabArgumentError) as e:
        spec.bind_strings([ROOT, "[1]", "pad"])
    assert e.value.identifier == "MATLAB:InputParser:ParamMissingValue"


def test_matlab_layer_takes_structs_and_keywords(sample_crop):
    spec = sample_crop.REG.find("XR_crop_dataset")
    got = spec.bind_matlab([ROOT, [1] * 6, {"PAD": True, "uuid": "u"}, "uuid", "v"])
    assert got["pad"] is True and got["uuid"] == "v"
    got = spec.bind_matlab([ROOT, [1] * 6], {"cropTYPE": "center"})
    assert got["crop_type"] == "center"
    got = spec.bind_matlab([ROOT, [1] * 6, MString("uuid"), "w"])
    assert got["uuid"] == "w"
    with pytest.raises(MatlabArgumentError):  # a string scalar is not char (ischar fails)
        spec.bind_matlab([ROOT, [1] * 6, "uuid", MString("w")])


def test_matlab_results_is_ip_results(sample_crop):
    spec = sample_crop.REG.find("XR_crop_dataset")
    out = sample_crop.crop_dataset(ROOT, [1, 2, 3, 4, 5, 6], channel_patterns=["CamB_ch0"])
    res = spec.matlab_results(out)
    names = list(res)
    assert names == sorted(names, key=str.lower)
    assert names[0] == "batchSize" and "dataPaths" in names
    assert describe(res["batchSize"]) == describe(MArray("double", (1, 3), [1024.0] * 3))
    assert isinstance(res["channelPatterns"], Cell) and res["channelPatterns"].shape == (1, 1)
    assert res["lastStartCoords"].shape == (0, 0)
    assert res["uuid"].mclass == "char" and res["uuid"].shape == (0, 0)
    assert res["pad"].mclass == "logical"
    from petakit5d.params import matlab_results

    with pytest.raises(LookupError):
        matlab_results("XR_not_a_function_anywhere", {})


def test_matlab_name_function_object(sample_crop):
    from petakit5d.matlab import MatlabFunction

    spec = sample_crop.REG.find("XR_crop_dataset")
    f = MatlabFunction(spec)
    assert f.__name__ == "XR_crop_dataset"
    assert "batchSize (default: [1024,1024,1024])" in f.__doc__
    assert f(ROOT, [1] * 6, "uuid", "x")["uuid"] == "x"


# -- declarations that do not fit their function ---------------------------------------------


def test_spec_errors():
    reg = Registry()

    @reg.register("f_wrong_name", Param("fooBar", "@isnumeric"))
    def f1(foo=1):
        return foo

    with pytest.raises(RegistryError, match="foo_bar"):
        _ = reg.find("f_wrong_name").params

    @reg.register("f_count", Param("a", "@isnumeric"))
    def f2(a, b):
        return a

    with pytest.raises(RegistryError, match="2 parameters"):
        _ = reg.find("f_count").params

    @reg.register("f_override", Param("lambda", "@isnumeric", py="lambda_"))
    def f3(*, lambda_=1):
        return lambda_

    assert reg.find("f_override").params[0].py == "lambda_"
    with pytest.raises(RegistryError, match="registered twice"):

        @reg.register("f_override", Param("lambda", "@isnumeric", py="lambda_"))
        def f4(*, lambda_=1):
            return lambda_

    with pytest.raises(RegistryError):
        reg.register("bad name")
    with pytest.raises(RegistryError):
        reg.register("f5", "notAParam")  # type: ignore[arg-type]


def test_optional_arguments_take_matlabs_rule():
    reg = Registry()

    @reg.register(
        "f_opt",
        Param("frame", "@(x) ischar(x) || iscell(x)"),
        Param("deconPath", "@(x) ischar(x) || isempty(x)"),
        Param("Overwrite", "@islogical"),
    )
    def f(frame, decon_path="", *, overwrite=False):
        return locals()

    spec = reg.find("f_opt")
    assert spec.call_matlab("a.tif", "out")["decon_path"] == "out"
    # a char optional takes a parameter name too (inputParser does the same)
    with pytest.raises(MatlabArgumentError):
        spec.call_matlab("a.tif", "Overwrite", True)
    assert spec.call_matlab("a.tif", "", "Overwrite", True)["overwrite"] is True
    assert spec.call_strings(["a.tif", "out", "overwrite", "true"])["overwrite"] is True


@pytest.mark.parametrize(
    ("validator", "alts"),
    [
        ("@ischar", [{"ischar"}]),
        ("@(x) ischar(x) || iscell(x)", [{"ischar"}, {"iscell"}]),
        ("@(x) isnumeric(x) && numel(x) == 3", [{"isnumeric"}]),
        ("@(x) isempty(x) || isa(x, 'function_handle')", [{"isempty"}, {"function_handle"}]),
        ("@(x) any(strcmpi(x, {'a', 'b'}))", None),
        ("@myCheck", None),
        ("(1, :) {mustBeNumeric}", [{"isnumeric"}]),
        ("(1, 1) logical", [{"islogical"}, {"isnumeric"}]),
        ("char", [{"ischar"}, {"isstring"}]),
        (None, None),
    ],
)
def test_validator_alternatives(validator, alts):
    got = validator_alternatives(validator)
    assert (None if got is None else [set(a) for a in got]) == alts


def test_python_values_against_matlab_predicates():
    import numpy as np

    from petakit5d.params.spec import _pred

    assert _pred("isnumeric", 1) and not _pred("isnumeric", True)
    assert _pred("isnumeric", None) and _pred("isempty", None)
    assert _pred("isnumeric", (1, 2.5)) and not _pred("isnumeric", ("a",))
    assert _pred("islogical", (True, False)) and not _pred("islogical", (1, 0))
    assert _pred("iscell", ("a", 1)) and _pred("iscell", [])
    assert _pred("ischar", "") and not _pred("ischar", MString("a"))
    assert _pred("isnumeric", np.zeros((2, 2, 2), dtype=np.uint16))
    assert _pred("islogical", np.array([True]))
    assert _pred("isscalar", 5) and not _pred("isscalar", (1, 2))
    assert _pred("isvector", (1, 2)) and not _pred("isvector", ((1, 2), (3, 4)))
    assert _pred("isstruct", {"a": 1}) and _pred("function_handle", len)
    assert _pred("isnumeric", MArray("uint8", (1, 1), [1]))
