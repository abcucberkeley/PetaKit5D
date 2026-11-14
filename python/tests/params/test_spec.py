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
    # MATLAB's box [ymin xmin zmin ymax xmax zmax] and batch [y x z]; Python's (z, y, x) forms
    mbox = [2, 3, 4, 100, 200, 50]
    a = sample_crop.crop_dataset(
        ROOT,
        (3, 1, 2, 50, 100, 200),
        batch_size=[64, 512, 256],
        pad=True,
        channel_patterns=["CamA_ch0"],
    )
    b = spec.call_matlab(
        ROOT, mbox, "BatchSize", [512, 256, 64], "pad", True, ChannelPatterns=["CamA_ch0"]
    )
    c = spec.call_strings(
        [ROOT, "[2,3,4,100,200,50]", "batchsize", "[512,256,64]", "Pad", "true",
         "channelPatterns", "{'CamA_ch0'}"]
    )  # fmt: skip
    import petakit5d.matlab as pm

    assert a == b == c
    assert a["batch_size"] == (64, 512, 256) and a["input_bbox"] == (3, 1, 2, 50, 100, 200)
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


def test_snake_case_values_are_checked_in_matlabs_form(sample_crop, monkeypatch):
    """The validator sees the value as MATLAB has it: an open end (None) of a box or a range is
    MATLAB's Inf, so @isnumeric passes; the function gets the value as given."""
    from petakit5d.core import arrays

    transposed = []
    monkeypatch.setattr(arrays, "to_matlab", lambda v, **kw: transposed.append(v) or v)
    got = sample_crop.crop_dataset(ROOT, (0, 1, 2, None, 10, None))
    assert got["input_bbox"] == (0, 1, 2, None, 10, None)
    reg = Registry()
    seen = []

    def vol_in(v):
        return v

    def vol_out(v):
        seen.append(v)
        return v

    @reg.register(
        "f_open_ends",
        Param("zRange", "@isnumeric", convert="range1"),
        Param("boxes", "@isnumeric", convert="bboxes"),
        Param("vol", "@isnumeric", convert="volume"),
        Param("other", "@isnumeric", convert=(vol_in, vol_out)),
        Param("pad", "@islogical"),
    )
    def f(z_range=None, *, boxes=None, vol=None, other=None, pad=False):
        return dict(locals())

    got = f((4, None), boxes=((0, 0, 0, 5, None, 7), (1, 2, 3, 4, 5, 6)))
    assert got["z_range"] == (4, None) and got["boxes"][0] == (0, 0, 0, 5, None, 7)
    # image arrays and a function's own conversions are checked as given (not converted)
    import numpy as np

    v = np.zeros((2, 3, 4), np.uint16)
    assert f(vol=v, other=v)["vol"] is v and seen == [] and transposed == []
    # values the kind can't convert are checked as given; failing values still fail
    assert f((1, 2, 3))["z_range"] == (1, 2, 3)
    for bad in (lambda: f(("a", None)), lambda: f(pad=None), lambda: f(boxes=[[None, "x"]])):
        with pytest.raises(MatlabArgumentError) as e:
            bad()
        assert e.value.identifier == "MATLAB:InputParser:ArgumentFailedValidation"


def test_string_conversions_follow_the_mcc_parser(sample_crop):
    spec = sample_crop.REG.find("XR_crop_dataset")
    got = spec.bind_strings(
        [ROOT, "1 2 3 10 20 30", "lastStartCoords", "[]", "uuid", "abc", "cropType", "{x}"]
    )
    assert got["input_bbox"] == (2, 0, 1, 30, 10, 20)
    assert got["last_start_coords"] is None
    assert got["uuid"] == "abc" and got["crop_type"] == "{x}"  # char: never converted
    got = spec.bind_strings([ROOT, "[1 1 1 10 20 Inf]", "lastStartCoords", "[5 6 7]"])
    assert got["input_bbox"] == (0, 0, 0, None, 10, 20)
    assert got["last_start_coords"] == (6, 4, 5)
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


def test_integer_classes_keep_their_class_through_the_binding(sample_crop):
    """MATLAB passes an integer class on (str2num('[512, uint16(256)]') is uint16): the function
    gets a NumPy value of that class, and the saved parameters have it too."""
    import numpy as np

    spec = sample_crop.REG.find("XR_crop_dataset")
    box = "[1,1,1,10,20,30]"
    got = spec.bind_strings(
        [ROOT, box, "batchSize", "[512, uint16(256), 64]", "cpusPerTask", "uint8(4)"]
    )
    bs, cpus = got["batch_size"], got["cpus_per_task"]
    assert bs.dtype == np.uint16 and bs.tolist() == [64, 512, 256]  # (z, y, x)
    assert type(cpus) is np.uint8 and cpus == 4
    res = spec.matlab_results(got)
    assert describe(res["batchSize"]) == describe(MArray("uint16", (1, 3), [512, 256, 64]))
    assert describe(res["cpusPerTask"]) == describe(MArray("uint8", (1, 1), [4]))
    got = spec.bind_matlab([ROOT, [1] * 6, "blockSize", MArray("int16", (3, 1), [5, 6, 7])])
    assert got["block_size"].dtype == np.int16 and got["block_size"].tolist() == [7, 5, 6]
    m = spec.bind_matlab([ROOT, [1] * 6, "cpusPerTask", MArray("uint32", (2, 2), [1, 2, 3, 4])])
    assert m["cpus_per_task"].dtype == np.uint32 and m["cpus_per_task"].tolist() == [[1, 3], [2, 4]]
    v = spec.bind_matlab([ROOT, [1] * 6, "cpusPerTask", MArray("int8", (1, 2), [4, -5])])
    assert v["cpus_per_task"].shape == (2,)  # vectors are 1-D, as the tuples of doubles
    # doubles (and empty integer arrays) stay Python values
    got = spec.bind_strings([ROOT, box, "batchSize", "[512, 256, 64]", "blockSize", "uint8([])"])
    assert got["batch_size"] == (64, 512, 256) and got["block_size"] is None


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
    out = sample_crop.crop_dataset(
        ROOT, (2, 0, 1, 30, 10, 20), channel_patterns=["CamB_ch0"], batch_size=(8, 16, 32)
    )
    res = spec.matlab_results(out)
    names = list(res)
    assert names == sorted(names, key=str.lower)
    assert names[0] == "batchSize" and "dataPaths" in names
    # MATLAB's values: [y x z] sizes, [ymin xmin zmin ymax xmax zmax] boxes
    assert describe(res["batchSize"]) == describe(MArray("double", (1, 3), [16.0, 32.0, 8.0]))
    assert describe(res["inputBbox"]) == describe(
        MArray("double", (1, 6), [1.0, 2.0, 3.0, 10.0, 20.0, 30.0])
    )
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


def test_optional_without_validator_takes_no_text():
    """inputParser reads text in an optional's place as a parameter name when the optional
    has no validator (XR_estimateComputingMemory's steps), and fails for any other text."""
    reg = Registry()

    @reg.register(
        "f_steps",
        Param("filePath", None),
        Param("steps", None),
        Param("imSize", "@isnumeric", convert="size_yxz"),
    )
    def f(file_path, steps=("deskew", "rotate"), *, im_size=None):
        return locals()

    spec = reg.find("f_steps")
    got = spec.call_matlab("", "imSize", [100, 200, 600])
    assert got["steps"] == ("deskew", "rotate") and got["im_size"] == (600, 100, 200)
    assert spec.call_matlab("", "steps", ["deskew"])["steps"] == ("deskew",)
    assert spec.call_matlab("", ["rotate"])["steps"] == ("rotate",)
    with pytest.raises(MatlabArgumentError) as e:
        spec.call_matlab("", "deskew")
    assert e.value.identifier == "MATLAB:InputParser:ArgumentFailedDefaultValidation"
    assert str(e.value) == (
        "The argument 'deskew' is a string scalar or character vector and does not match any "
        "parameter names. It failed validation for the argument 'steps'."
    )


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
