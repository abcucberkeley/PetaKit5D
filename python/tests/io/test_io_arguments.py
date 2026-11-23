"""The io wrappers' arguments blocks through the MATLAB-name layer, against MATLAB's wrappers.

MATLAB converts readtiff's, writetiff's, readzarr's, writezarr's and createzarr's inputs to their
declared classes and sizes before the wrapper runs (``params.declarations``; each declaration's
rules on a battery of values: ``tests/params/test_declarations.py``). These cases call the lab's
wrappers themselves with values that meet the rules:

- conversions copied: a column or a 1x1xN vector as a range, a box or a size, numbers as flags,
  an array of flags as writezarr's ``create``; the files written and read back are MATLAB's;
- MATLAB's errors: an upper-case ``Compression`` or ``Mode``, a flag of two elements or NaN, a
  matrix as a range or a box, ``Inf`` in a range or a box (which the layer used to read as an open
  end), text as a range, a level of two elements: each MATLAB's identifier and message;
- conversions the port refuses (decided by the user on 2026-11-22): a number or a cell as a
  file name, a handle as a file name, text as a flag, numbers or a cell as the compression, a
  number as the chunk separator, a cell as the class. Each case records what MATLAB does (the
  files it writes: ``writetiff(img, 65)`` writes ``A``), and the test checks that the port
  refuses the call with its own error and writes nothing.
"""

from __future__ import annotations

import math
import os

import numpy as np
import pytest

pytest.importorskip("cpptiff")
pytest.importorskip("cppzarr")

from petakit5d.params import ConversionRefusedError, MatlabValidationError  # noqa: E402
from petakit5d.testing import oracle  # noqa: E402
from petakit5d.testing.cases import Case, In, Value, Work, mcell, mcol  # noqa: E402
from petakit5d.testing.compare import MFunc  # noqa: E402

SIZE = (12, 9, 6)  # MATLAB's (y, x, z), never cubic


def _volume() -> np.ndarray:
    a = np.random.default_rng(33).integers(0, 65535, size=SIZE, dtype=np.uint16)
    return np.asfortranarray(a)


def _region() -> np.ndarray:
    a = np.random.default_rng(34).integers(0, 65535, size=(10, 7, 5), dtype=np.uint16)
    return np.asfortranarray(a)


INPUTS = {"m": Value(_volume)}


def _matlab():
    from petakit5d import matlab

    return matlab


# --- the Python sides: MATLAB's calls through the MATLAB-name layer ------------------------------


def _read_tiff(m, opts, work):
    fn = work + "a.tif"
    _matlab().writetiff(m, fn)
    return _matlab().readtiff(fn, *opts)


def _write_tiff(m, opts, work):
    fn = work + "a.tif"
    _matlab().writetiff(m, fn, *opts)
    return _matlab().readtiff(fn)


def _read_zarr(m, opts, work):
    fn = work + "m.zarr"
    _matlab().writezarr(m, fn, "blockSize", [8, 8, 4])
    return _matlab().readzarr(fn, *opts)


def _write_zarr(m, opts, work):
    fn = work + "m.zarr"
    _matlab().writezarr(m, fn, *opts)
    return _matlab().readzarr(fn)


def _create_zarr(m, opts, work):
    fn = work + "m.zarr"
    _matlab().createzarr(fn, *opts)
    _matlab().writezarr(m, fn, "create", False)
    return _matlab().readzarr(fn)


def _region_zarr(r, opts, work):
    fn = work + "m.zarr"
    _matlab().createzarr(fn, "dataSize", [12, 9, 6], "blockSize", [8, 8, 4])
    _matlab().writezarr(r, fn, *opts)
    return _matlab().readzarr(fn)


_READ_TIFF = "fn = fullfile(work, 'a.tif'); writetiff(in{1}, fn); out{1} = readtiff(fn, in{2}{:});"
_WRITE_TIFF = "fn = fullfile(work, 'a.tif'); writetiff(in{1}, fn, in{2}{:}); out{1} = readtiff(fn);"
_READ_ZARR = (
    "fn = fullfile(work, 'm.zarr'); writezarr(in{1}, fn, 'blockSize', [8 8 4]); "
    "out{1} = readzarr(fn, in{2}{:});"
)
_WRITE_ZARR = (
    "fn = fullfile(work, 'm.zarr'); writezarr(in{1}, fn, in{2}{:}); out{1} = readzarr(fn);"
)
_CREATE_ZARR = (
    "fn = fullfile(work, 'm.zarr'); createzarr(fn, in{2}{:}); writezarr(in{1}, fn, 'create', "
    "false); out{1} = readzarr(fn);"
)
_REGION_ZARR = (
    "fn = fullfile(work, 'm.zarr'); createzarr(fn, 'dataSize', [12 9 6], 'blockSize', [8 8 4]); "
    "writezarr(in{1}, fn, in{2}{:}); out{1} = readzarr(fn);"
)
_SIDES = {
    "readtiff": (_READ_TIFF, _read_tiff, True),
    "writetiff": (_WRITE_TIFF, _write_tiff, True),
    "readzarr": (_READ_ZARR, _read_zarr, False),
    "writezarr": (_WRITE_ZARR, _write_zarr, False),
    "createzarr": (_CREATE_ZARR, _create_zarr, False),
    "region": (_REGION_ZARR, _region_zarr, False),
}


def _case(name: str, side: str, matlab_opts: list, py_opts: list, error: str | None = None) -> Case:
    """a call of the lab's wrapper with ``matlab_opts`` (MATLAB's values), and the MATLAB-name
    layer's with ``py_opts`` (the same values in the registry's forms)"""
    code, python, files = _SIDES[side]
    data = "r" if side == "region" else "m"
    return Case(
        f"args.{name}",
        code=code,
        args=[In(data), mcell(*matlab_opts), Work()],
        python=python,
        py_args=[In(data), py_opts, Work()],
        inputs={data: Value(_region if side == "region" else _volume)},
        files=files,
        **(
            {"expect_error": error, "py_error": (MatlabValidationError,), "error_message": True}
            if error
            else {}
        ),
    )


COL = np.array([[2.0], [5.0]])
CASES = [
    # conversions copied: vectors made rows, numbers made logical values
    _case("readtiff.range_column", "readtiff", ["range", mcol([2.0, 5.0])], ["range", COL]),
    _case(
        "readtiff.range_1x1x2",
        "readtiff",
        ["range", np.array([2.0, 5.0]).reshape(1, 1, 2)],
        ["range", np.array([2.0, 5.0]).reshape(1, 1, 2)],
    ),
    _case(
        "readzarr.box_column",
        "readzarr",
        ["inputBbox", mcol([2.0, 3.0, 1.0, 10.0, 8.0, 5.0])],
        ["inputBbox", np.array([[2.0], [3.0], [1.0], [10.0], [8.0], [5.0]])],
    ),
    _case(
        "writezarr.box_column",
        "region",
        ["bbox", mcol([2.0, 2.0, 2.0, 11.0, 8.0, 6.0]), "create", False],
        ["bbox", ((2,), (2,), (2,), (11,), (8,), (6,)), "create", False],
    ),
    _case(
        "createzarr.sizes_column",
        "createzarr",
        [
            "dataSize",
            mcol([12.0, 9.0, 6.0]),
            "blockSize",
            np.array([8.0, 8.0, 4.0]).reshape(1, 1, 3),
            "shardSize",
            None,
            "zarrSubSize",
            np.zeros((0, 1)),
        ],
        [
            "dataSize",
            np.array([[12.0], [9.0], [6.0]]),
            "blockSize",
            np.array([8.0, 8.0, 4.0]).reshape(1, 1, 3),
            "shardSize",
            None,
            "zarrSubSize",
            np.zeros((0, 1)),
        ],
    ),
    _case("writetiff.group_number", "writetiff", ["groupWrite", 0.0], ["groupWrite", 0.0]),
    _case("writezarr.uuid_number", "writezarr", ["uuid", 2.0], ["uuid", 2.0]),
    _case(
        "writezarr.create_flags",
        "writezarr",
        ["create", np.array([True, True])],
        ["create", np.array([True, True])],
    ),
    # MATLAB's errors
    _case(
        "writetiff.compression_upper",
        "writetiff",
        ["Compression", "LZW"],
        ["Compression", "LZW"],
        "MATLAB:validators:mustBeMember",
    ),
    _case(
        "writetiff.mode_upper",
        "writetiff",
        ["Mode", "Parallel"],
        ["Mode", "Parallel"],
        "MATLAB:validators:mustBeMember",
    ),
    _case(
        "writetiff.group_two",
        "writetiff",
        ["groupWrite", np.array([True, False])],
        ["groupWrite", (True, False)],
        "MATLAB:validation:IncompatibleSize",
    ),
    _case(
        "writetiff.group_nan",
        "writetiff",
        ["groupWrite", math.nan],
        ["groupWrite", math.nan],
        "MATLAB:validation:UnableToConvert",
    ),
    _case(
        "writezarr.overwrite_two",
        "writezarr",
        ["overwrite", [1.0, 1.0]],
        ["overwrite", [1.0, 1.0]],
        "MATLAB:validation:IncompatibleSize",
    ),
    _case(
        "readtiff.range_inf",
        "readtiff",
        ["range", [2.0, math.inf]],
        ["range", [2.0, math.inf]],
        "MATLAB:validators:mustBeFinite",
    ),
    _case(
        "readtiff.range_matrix",
        "readtiff",
        ["range", np.array([[1.0, 2.0], [3.0, 4.0]])],
        ["range", ((1.0, 2.0), (3.0, 4.0))],
        "MATLAB:validation:IncompatibleSize",
    ),
    _case(
        "readtiff.range_text",
        "readtiff",
        ["range", "ab"],
        ["range", "ab"],
        "MATLAB:validators:mustBeNumeric",
    ),
    _case(
        "readzarr.box_inf",
        "readzarr",
        ["inputBbox", [1.0, 1.0, 1.0, math.inf, 5.0, 5.0]],
        ["inputBbox", (1, 1, 1, math.inf, 5, 5)],
        "MATLAB:validators:mustBeFinite",
    ),
    _case(
        "createzarr.clevel_row",
        "createzarr",
        ["clevel", [1.0, 2.0]],
        ["clevel", (1, 2)],
        "MATLAB:validation:IncompatibleSize",
    ),
    _case(
        "writezarr.block_matrix",
        "writezarr",
        ["blockSize", np.ones((2, 3))],
        ["blockSize", np.ones((2, 3))],
        "MATLAB:validation:IncompatibleSize",
    ),
]


def _missing_create(m, flags, work):
    _matlab().writezarr(m, work + "m.zarr", "create", flags)


# writezarr's create of an array of flags, on a missing array: MATLAB's if (all true) is false
CASES.append(
    Case(
        "args.writezarr.create_flags_missing",
        code="writezarr(in{1}, fullfile(work, 'm.zarr'), 'create', in{2});",
        args=[In("m"), np.array([True, False]), Work()],
        python=_missing_create,
        py_args=[In("m"), np.array([True, False]), Work()],
        inputs=INPUTS,
        files=True,
        expect_error=True,
        py_error=(FileNotFoundError,),
    )
)


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# --- conversions the port refuses: MATLAB's result recorded, the port's refusal checked ----------

_FILES = "d = dir(work); out{1} = sort({d(~[d.isdir]).name});"
# case: (MATLAB's code, its arguments after the volume, the port's call of the same); the tiff
# files MATLAB writes are recorded with their tags (the compression in the codes case: zstd)
REFUSED = {
    "file_number": (
        "writetiff(in{1}, in{2}); " + _FILES,
        [65.0],
        lambda m, w: _matlab().writetiff(m, 65.0),
    ),
    "file_cell": (
        "writetiff(in{1}, in{2}); " + _FILES,
        [mcell("c.tif")],
        lambda m, w: _matlab().writetiff(m, ("c.tif",)),
    ),
    "file_handle": (
        "writetiff(in{1}, in{2}); " + _FILES,
        [MFunc("@sin")],
        lambda m, w: _matlab().writetiff(m, math.sin),
    ),
    "group_text": (
        "writetiff(in{1}, fullfile(work, 'b.tif'), 'groupWrite', in{2}); " + _FILES,
        ["a"],
        lambda m, w: _matlab().writetiff(m, w + "b.tif", "groupWrite", "a"),
    ),
    "compression_codes": (
        "writetiff(in{1}, fullfile(work, 'd.tif'), 'Compression', in{2}); " + _FILES,
        [[122.0, 115.0, 116.0, 100.0]],
        lambda m, w: _matlab().writetiff(m, w + "d.tif", "Compression", [122, 115, 116, 100]),
    ),
    "compression_cell": (
        "writetiff(in{1}, fullfile(work, 'e.tif'), 'Compression', in{2}); " + _FILES,
        [mcell("none")],
        lambda m, w: _matlab().writetiff(m, w + "e.tif", "Compression", ("none",)),
    ),
    "separator_number": (
        "fn = fullfile(work, 'm.zarr'); writezarr(in{1}, fn, 'dimSeparator', in{2}); "
        "z = jsondecode(fileread(fullfile(fn, '.zarray'))); out{1} = z.dimension_separator;",
        [47.0],
        lambda m, w: _matlab().writezarr(m, w + "m.zarr", "dimSeparator", 47),
    ),
    "class_cell": (
        "fn = fullfile(work, 'c.zarr'); createzarr(fn, 'dataSize', [4 5 3], 'blockSize', "
        "[4 5 3], 'dtype', in{2}); z = jsondecode(fileread(fullfile(fn, '.zarray'))); "
        "out{1} = z.dtype;",
        [mcell("uint8")],
        lambda m, w: _matlab().createzarr(
            w + "c.zarr", "dataSize", [4, 5, 3], "blockSize", [4, 5, 3], "dtype", ("uint8",)
        ),
    ),
}


def _refused_case(name: str) -> Case:
    code, args, _ = REFUSED[name]
    return Case(
        f"args.refused.{name}",
        code=code,
        args=[In("m"), *args, Work()],
        python=None,
        inputs=INPUTS,
        files=True,
    )


REFUSED_CASES = [_refused_case(n) for n in REFUSED]


@oracle.parametrize(REFUSED_CASES)
def test_matlab_refused(case, tmp_path):
    """MATLAB carries out these calls (its outputs and files recorded); the port refuses each with
    its own error, before it writes anything"""
    oracle.check(case, tmp_path)
    call = REFUSED[case.id.rsplit(".", 1)[1]][2]
    work = tmp_path / "port"
    work.mkdir()
    cwd = os.getcwd()
    os.chdir(work)
    try:
        with pytest.raises(ConversionRefusedError) as err:
            call(_volume(), work.as_posix() + "/")
    finally:
        os.chdir(cwd)
    assert err.value.identifier == "PETAKIT5D:params:conversionRefused"
    assert not os.listdir(work)


def test_the_refused_cases_record_matlabs_conversions():
    """what MATLAB did, as the references hold it (the files it wrote are compared with the
    store; the manifest keeps their names)"""
    want = {
        "file_number": ["A"],
        "file_cell": ["c.tif"],
        "file_handle": ["sin"],
        "group_text": ["b.tif"],
        "compression_codes": ["d.tif"],
        "compression_cell": ["e.tif"],
        "separator_number": "/",
        "class_cell": "|u1",
    }
    for case in REFUSED_CASES:
        ref = oracle.reference(case)
        out = ref.outputs[0]
        name = case.id.rsplit(".", 1)[1]
        if hasattr(out, "items") and not hasattr(out, "mclass"):  # a cell of file names
            assert out.items is not None, f"{name}: the manifest keeps the names in full"
            assert [x.text() for x in out.items] == want[name], name
        else:
            assert out.data is not None, f"{name}: the manifest keeps the text in full"
            assert out.text() == want[name], name


def test_write_zarr_takes_an_array_of_flags_for_create(tmp_path):
    """MATLAB's if: all true and not empty (the MATLAB-name layer gives writezarr's create, of
    any size, as MATLAB's logical array)"""
    from petakit5d.io import read_zarr, write_zarr

    vol = np.arange(24, dtype=np.uint16).reshape(2, 3, 4)
    for flags in ((True, False), np.array([], dtype=bool), np.array([[True], [False]])):
        with pytest.raises(FileNotFoundError):
            write_zarr(vol, tmp_path / "a.zarr", create=flags)
    write_zarr(vol, tmp_path / "a.zarr", create=np.array([True, True]))
    assert np.array_equal(read_zarr(tmp_path / "a.zarr"), vol)
    write_zarr(vol, tmp_path / "b.zarr", create=((True, True),))
    with pytest.raises(FileNotFoundError):
        _matlab().writezarr(vol.T, str(tmp_path / "c.zarr"), "create", [1, 0])
    _matlab().writezarr(vol.T, str(tmp_path / "c.zarr"), "create", np.ones((2, 2)))
    assert os.path.isdir(tmp_path / "c.zarr")
