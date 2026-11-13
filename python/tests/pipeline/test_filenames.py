"""File names of datasets: MATLAB's regular expressions, channel patterns, the tile name
grammar, parseImageFilenames, groupImageFilenamesByChannels and XR_parseImageFilenames.

Every name here is made up. The oracle cases compare with MATLAB; the unit tests pin the
translation rules and the streaming rule (file times can't go into references).
"""

from __future__ import annotations

import os
import re
import time
from pathlib import Path

import numpy as np
import pytest

from petakit5d.pipeline import filenames as fnm
from petakit5d.pipeline.filenames import (
    MatlabError,
    group_image_filenames_by_channels,
    matches_channel_pattern,
    matlab_regex,
    parse_dataset_filenames,
    parse_filename_fields,
    parse_image_filenames,
)
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, InPath, Tree, mcell

# ---------------------------------------------------------------------------------------------
# a corpus of made-up names
# ---------------------------------------------------------------------------------------------


def tile(
    cam: str | None = "A",
    ch: int = 0,
    x: int = 0,
    y: int = 0,
    z: int = 0,
    t: int = 0,
    it: str = "0000",
    sub: str = "_0000",
    laser: int = 488,
    prefix: str = "",
    msec: int = 0,
    abs_: int = 1234567,
    stack: int = 0,
    suffix: str = "",
) -> str:
    """an LLSM style tile name (without extension), as the scope software writes them"""
    cam_s = f"_Cam{cam}" if cam else ""
    return (
        f"{prefix}Scan_Iter_{it}{sub}{cam_s}_ch{ch}_CAM1_stack{stack:04d}_{laser}nm_"
        f"{msec:07d}msec_{abs_:010d}msecAbs_{x:03d}x_{y:03d}y_{z:03d}z_{t:04d}t{suffix}"
    )


NAMES = [
    tile(),
    tile("B", 1, laser=560),
    tile("A", 0, x=-1, y=2, z=-3),
    tile("A", 2, it="0003", sub="_0001_0002", laser=642, t=17),
    tile("B", 0, prefix="sample01_", msec=1500, abs_=987654321),
    tile(None, 0),
    tile(None, 1, x=1, prefix="emb_"),
    tile("A", 0, sub=""),
    tile("A", 0, suffix="_part0001"),
    tile("AB", 0),
    "RAW_exp01_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000558259msecAbs_000x_000y_000z_0000t",
    "ff_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0096911533msecAbs_000x_000y_000z_0017t",
    "cell1_488_ch0_stack0000",
    "Tile_Z001_Y002_X003_ch0",
    "scan_iter_0000_0000_cama_ch0_cam1_stack0000_488nm_0000000msec_0000000001msecabs_"
    "000x_000y_000z_0000t",
    "",
    "plain",
    "SCAN_ITER_0000_0000_CamA_ch0_cam1_STACK0000_488NM_0000000MSEC_0000000001MSECABS_"
    "001X_002Y_003Z_0004T",
    tile("A", 0, prefix="a") + "-" + tile("B", 1, x=5, it="0002", prefix="b"),
]

# MATLAB regular expressions with texts: each one pins a translation rule (module docstring)
REGEXP_CORPUS = [
    (r"a\bb", "a\bb"),
    (r"a\bb", "a b"),
    (r"\<b", "a b"),
    (r"\<b", "ab"),
    (r"b\>", "ab c"),
    (r"b\>", "abc"),
    (r"ab$", "ab\n"),
    (r"ab$", "ab"),
    (r"^ab", "\nab"),
    (r"a.b", "a\nb"),
    (r"a\wb", "aéb"),
    (r"a\db", "a٣b"),
    (r"a\sb", "a\x1cb"),
    (r"a\sb", "a\vb"),
    (r"a\Sb", "a\x1cb"),
    (r"a[\s]b", "a\tb"),
    (r"a\qb", "aqb"),
    (r"a\_b", "a_b"),
    (r"a\x5fb", "a_b"),
    (r"a\x5fb", "a׻"),
    (r"a\x{5f}b", "a_b"),
    (r"a\o137b", "a_b"),
    (r"a\o{137}b", "a_b"),
    (r"(?<n>a)\k<n>", "aa"),
    (r"(a)\1", "aa"),
    (r"(?i)a_b", "A_B"),
    (r"a(?i)_B", "a_b"),
    (r"a(?i)b|C", "xc"),
    (r"(a(?i)b)|C", "xc"),
    (r"(a(?i)b|C)", "xc"),
    (r"a)", "a)"),
    (r"*a", "*a"),
    (r"a**", "aa*"),
    (r"a*+b", "aab"),
    (r"a*?b", "aab"),
    (r"a+?", "aa"),
    (r"a??b", "ab"),
    (r"a{2}", "aaa"),
    (r"a{2,}", "aaaa"),
    (r"a{1,2}", "aaa"),
    (r"a{,1}b", "a{,1}b"),
    (r"a{,1}b", "ab"),
    (r"a{x}", "a{x}"),
    (r"{2}a", "{2}a"),
    (r"(?>a+)b", "aab"),
    (r"(?#note)ab", "ab"),
    (r"(?:ab)+", "abab"),
    (r"x(?=y)", "xy"),
    (r"x(?!y)", "xyxz"),
    (r"(?<=a_)b", "a_b"),
    (r"(?<!a)b", "abcb"),
    (r"[\w.]+", "a.b c"),
    (r"[^\w]", "ab-c"),
    (r"[]a]", "]"),
    (r"[^]a]", "a]b"),
    (r"[a-c]+", "abcd"),
    (r"[a&&b]", "&"),
    (r"[a|b]", "|"),
    (r"[a~b]", "~"),
    (r"[[a]", "["),
    (r"[[a]", "[a"),
    (r"[[a]]", "a"),
    (r"[x[a]]", "a"),
    (r"[x[a]]", "["),
    (r"[x[a]]", "]"),
    (r"[x[a]]", "x]"),
    (r"[[:alpha:]]", "b"),
    (r"[[:alpha:]]", ":"),
    (r"[[:digit:]]", "5"),
    (r"a[[b]c", "a[c"),
    (r"[x[a-c]y]", "b"),
    (r"[]", "a"),
    (r"[^]", "x"),
    (r"[^]a]", "xa]"),
    (r"[]a]", "a]"),
    (r"[]a]", "]a]"),
    (r"a[b]c[", "abc["),
    (r"a(b", "xab"),
    (r"(a|b", "(b"),
    (r"(?<n>a", "(?<n>a"),
    (r"a(?i)b(c)D", "aBcd"),
    (r"a|", "b"),
    (r"|a", "b"),
    (r"()", "b"),
    (r"a#b", "a#b"),
    (r"a b", "a b"),
    (r"a\.tif", "a.tif"),
    (r"a.tif", "aXtif"),
    (r"\\", "a\\b"),
    (r"a\\", "a\\"),
    (r"\t\d", "\t5"),
    ("a\\", "a\\"),
    ("a\\", "a"),
    (r"a(", "a("),
    (r"a[b", "a[b"),
    (r"a{2", "a{2"),
    (r"a{2,3", "a{2,3"),
    (r"((a)", "a"),
    (r"(?<na-me>a)", "a"),
    (r"_part\d+\.tif$", "x_part0012.tif"),
    (r"^(?!.*_part[0-9]*.tif).*$", "x_part0012.tif"),
    (r"^(?!.*_part[0-9]*.tif).*$", "x.tif"),
    (r"^(?!.*_part[0-9]*.tif).*$", ""),
    (r"Cam[AB]_ch[01]", "Scan_CamB_ch1_x"),
    (r"ch0_.*Iter_000[0-4]", "ch0_x_Iter_0003"),
]
# the forms MATLAB has and the translation refuses (variable-length lookbehind, dynamic
# expressions): the port fails with 'pk5d:regexp' instead of guessing
UNSUPPORTED = [r"(?<=a.?)b", r"(??a)", r"a(?@disp(1))"]

# channel patterns: as the tests, demos and wrappers use them, and patterns that are invalid
# expressions (MATLAB fails even where the text contains them)
CHANNEL_PATTERNS = [
    "CamA_ch0",
    "CamA",
    "CamB_ch0",
    "cama_ch0",
    "ch0",
    "Iter",
    "Tile",
    "Scan",
    "488",
    "000x_003y_000z",
    "0000t",
    "Scan_Iter_0000_0000_CamB_ch0",
    ".zarr",
    "Cam[AB]_ch[01]",
    "ch0_.*Iter_000[0-4]",
    "_ch[0-9]_",
    "",
    "Cam(A",
    "Cam)A",
    "Cam[A",
    "a{2",
    "*Cam",
]


def regexp_once(pats, texts, ignore_case=False):
    out = []
    for p, s in zip(pats, texts, strict=True):
        try:
            m = matlab_regex(p, ignore_case).search(s)
            out.append(m.group() if m else "")
        except MatlabError as err:
            out.append("ERR:" + err.identifier)
    return out


MATLAB_REGEXP_ONCE = (
    "p = in{1}; s = in{2}; r = cell(size(p)); "
    "for k = 1:numel(p); try; r{k} = %s(s{k}, p{k}, 'match', 'once'); "
    "catch ME; r{k} = ['ERR:' ME.identifier]; end; end; out = {r};"
)


def contains_matrix(names, patterns, ignore_case):
    out = []
    for p in patterns:
        try:
            out.append(np.array(matches_channel_pattern(names, p, ignore_case), bool))
        except MatlabError as err:
            out.append("ERR:" + err.identifier)
    return out


MATLAB_CONTAINS = (
    "n = in{1}; p = in{2}; ic = in{3}; r = cell(size(p)); "
    "for k = 1:numel(p); try; "
    "r{k} = contains(n, p{k}, 'IgnoreCase', ic) | contains(n, channel_regexp_pattern(p{k}), "
    "'IgnoreCase', ic); "
    "catch ME; r{k} = ['ERR:' ME.identifier]; end; end; out = {r};"
)

# the tile name grammar, as XR_parseSettingFiles_wrapper's lines
MATLAB_GRAMMAR = (
    "fsnames = in{1}; "
    "camExpression = '(?<prefix>\\w*)Scan_Iter_(?<Iter>(\\d+_)*\\d+?)_Cam(?<Cam>\\w+)_ch(?<ch>\\d+)"
    "_CAM1_stack(?<stack>\\d+)_(?<laser>\\d+)nm_(?<abstime>\\d+)msec_(?<fpgatime>\\d+)msecAbs_"
    "(?<x>-?\\d+)x_(?<y>-?\\d+)y_(?<z>-?\\d+)z_(?<t>\\d+)t'; "
    "expression = '(?<prefix>\\w*)Scan_Iter_(?<Iter>(\\d+_)*\\d+?)_ch(?<ch>\\d+)_CAM1_stack"
    "(?<stack>\\d+)_(?<laser>\\d+)nm_(?<abstime>\\d+)msec_(?<fpgatime>\\d+)msecAbs_(?<x>-?\\d+)x_"
    "(?<y>-?\\d+)y_(?<z>-?\\d+)z_(?<t>\\d+)t'; "
    "hasCam = ~cellfun(@isempty, regexp(fsnames, '_Cam\\w_ch', 'match')); "
    "tmp = cell(numel(fsnames), 1); "
    "tmp(hasCam) = regexpi(fsnames(hasCam), camExpression, 'names', 'once'); "
    "tmp(~hasCam) = regexpi(fsnames(~hasCam), expression, 'names', 'once'); out = {tmp};"
)


def grammar(fsnames):
    return parse_filename_fields(fsnames)


# ---------------------------------------------------------------------------------------------
# dataset folders
# ---------------------------------------------------------------------------------------------


def _touch(folder: Path, names) -> None:
    for n in names:
        (folder / n).write_bytes(b"")


def make_listing(root: Path) -> None:
    # parseImageFilenames lists names only: empty files are enough
    a = root / "d1"
    a.mkdir()
    _touch(a, [tile("A", 0, y=i) + ".tif" for i in range(3)])
    _touch(a, [tile("B", 0, y=i) + ".tif" for i in range(2)])
    _touch(a, [tile("A", 1) + ".tif", tile("A", 0, suffix="_part0001") + ".tif"])
    _touch(a, ["notes.txt", "Upper_CAMA_CH0.tif", "x_CamA_ch0.tiff", ".hidden_CamA_ch0.tif"])
    _touch(a, ["img_CamA_ch0.nd2", "img_CamB_ch0.CZI", "img_CamA_ch0.czi.bak"])
    (a / (tile("A", 0, t=1) + ".zarr")).mkdir()
    (a / (tile("B", 0, t=1) + ".zarr")).mkdir()
    b = root / "d2"
    b.mkdir()
    _touch(b, [tile(None, 0, x=i) + ".tif" for i in range(2)])
    (root / "empty").mkdir()


def listing_out(r: fnm.ImageFilenames):
    return [r.fnames, r.fsns, r.fd_inds + 1.0, r.filepaths, r.ch_inds + 1.0]


def py_listing(d, sub, zarr, pats, fmts=None):
    paths = [d + s for s in sub] if isinstance(sub, list) else d + sub
    if fmts is None:
        return parse_image_filenames(paths, zarr, pats)
    return parse_image_filenames(paths, zarr, pats, data_formats=fmts)


ML_LISTING = (
    "d = in{1}; sub = in{2}; if iscell(sub); p = cellfun(@(x) [d x], sub, 'unif', 0); "
    "else; p = [d sub]; end; "
    "if numel(in) > 4; [a, b, c, e, f] = parseImageFilenames(p, in{3}, in{4}, "
    "'dataFormats', in{5}); "
    "else; [a, b, c, e, f] = parseImageFilenames(p, in{3}, in{4}); end; out = {a, b, c, e, f};"
)

LISTING_OUTPUTS = ["fnames", "fsns", "fd_inds", "filepaths", "ch_inds"]


def _listing_case(cid, sub, zarr, pats, fmts=None, error=False):
    args = [InPath("root"), sub, zarr, pats] + ([fmts] if fmts is not None else [])
    return Case(
        cid,
        code=ML_LISTING,
        args=args,
        nargout=5,
        outputs=LISTING_OUTPUTS,
        python=py_listing,
        py_out=listing_out,
        inputs={"root": Tree(make_listing, "root")},
        warnings="ids",
        expect_error=error,
        py_error=MatlabError,
    )


EMPTY = mcell()


def make_groups_corpus():
    a = [tile("A", 0, y=i) for i in range(3)]
    b = [tile("B", 0, y=i) for i in range(3)]
    return a, b


def group_out(g):
    return [np.asarray(g, np.float64) + 1.0]


def py_group(names, pats, ch=None):
    if ch is None:
        return group_image_filenames_by_channels(names, pats)
    return group_image_filenames_by_channels(names, pats, [c - 1 for c in ch])


def _group_case(cid, names, pats, ch=None, error=False):
    args = [names, pats] + ([ch] if ch is not None else [])
    py_args = [[] if n is EMPTY else n for n in args]
    return Case(
        cid,
        matlab="groupImageFilenamesByChannels",
        args=args,
        py_args=py_args,
        python=py_group,
        py_out=group_out,
        expect_error=error,
        py_error=MatlabError,
    )


# XR_parseImageFilenames on small tiff and zarr datasets


def _tiffs(folder: Path, names, shape=(4, 6, 3), dtype=np.uint16) -> None:
    from petakit5d.testing import synth

    for i, n in enumerate(names):
        z = shape[2] + (i % 2)
        synth.write_tiff_matlab(folder / n, np.zeros((shape[0], shape[1], z), dtype))


def make_dataset(root: Path) -> None:
    a = root / "ds1"
    a.mkdir()
    _tiffs(a, [tile("A", 0, y=i) + ".tif" for i in range(3)])
    _tiffs(a, [tile("B", 0, y=i) + ".tif" for i in range(3)], dtype=np.uint8)
    _tiffs(a, [tile("A", 1, y=0) + ".tif"])
    b = root / "ds2"
    b.mkdir()
    _tiffs(b, [tile("A", 0, x=i, t=2) + ".tif" for i in range(2)], dtype=np.float32)
    (root / "ds3").mkdir()
    p = root / "part"
    p.mkdir()
    base = tile("A", 0)
    _tiffs(p, [base + ".tif", base + "_part0002.tif", base + "_part0001.tif"])
    _tiffs(p, [tile("A", 0, t=1) + ".tif", tile("B", 0, t=1) + "_part0001.tif"])
    zr = root / "zarr"
    zr.mkdir()
    from petakit5d.testing import synth

    for i in range(2):
        synth.write_zarr_matlab(
            zr / (tile("A", 0, y=i) + ".zarr"), np.zeros((5, 4, 3 + i), np.uint16)
        )
    # setting files next to the images (the stage line only: see test_settings.py)
    s = root / "set"
    s.mkdir()
    stage = "X Stage Offset, Interval (um), # of Pixels for Excitation (0) :\t0\t{}\t101"
    for i, interval in enumerate(["0.5", "-0.5", "-1"]):
        n = tile("A", 0, y=i)
        _tiffs(s, [n + ".tif"])
        text = "Settings\n" + stage.format(interval) + "\nEnd\n\n"
        (s / (n + "_Settings.txt")).write_bytes(text.encode())
    _tiffs(s, [tile("A", 0, y=1) + "_part0001.tif"])


def dataset_out(r: fnm.DatasetFilenames):
    return [
        r.fnames,
        r.fd_inds + 1.0,
        r.gfnames,
        r.partialvols,
        r.data_sizes,
        r.flip_zstack,
        r.latest_modify_times,
        r.ftp_inds + 1.0,
        r.mask_fullpaths,
    ]


DATASET_OUTPUTS = [
    "fnames",
    "fdinds",
    "gfnames",
    "partialvols",
    "dataSizes",
    "flipZstack_mat",
    "latest_modify_times",
    "FTP_inds",
    "maskFullpaths",
]

ML_DATASET = (
    "d = in{1}; p = cellfun(@(x) [d x '/'], in{2}, 'unif', 0); "
    "dp = cellfun(@(x) [d x], in{6}, 'unif', 0); "
    "[o1, o2, o3, o4, o5, o6, o7, o8, o9] = XR_parseImageFilenames(p, in{3}, in{4}, in{5}, "
    "in{7}, dp, false, 1, in{8}); out = {o1, o2, o3, o4, o5, o6, o7, o8, o9};"
)


def py_dataset(d, sub, pats, parse, flip, decon_sub, decon, zarr):
    return parse_dataset_filenames(
        [d + s + "/" for s in sub],
        pats,
        parse,
        flip,
        decon,
        [d + s for s in decon_sub],
        False,
        1,
        zarr,
    )


def _dataset_case(
    cid,
    sub,
    pats,
    flip=False,
    decon_sub=None,
    decon=False,
    zarr=False,
    parse=False,
    error=False,
    make=None,
    marks=(),
):
    decon_sub = decon_sub if decon_sub is not None else [s + "_decon" for s in sub]
    return Case(
        cid,
        code=ML_DATASET,
        args=[InPath("root"), sub, pats, parse, flip, decon_sub, decon, zarr],
        expect_error=error,
        py_error=MatlabError,
        nargout=9,
        outputs=DATASET_OUTPUTS,
        python=py_dataset,
        py_out=dataset_out,
        inputs={"root": Tree(make or make_dataset, "root")},
        warnings="ids",
        marks=marks,
    )


def make_special(root: Path) -> None:
    # main names are expressions in groupPartialVolumeFiles ([name, '_part\d+$']): a name
    # with '+' or '(1)' doesn't match its own part files, which are then left out
    d = root / "sp"
    d.mkdir()
    for base in ("x+y_CamA_ch0", "z(1)_CamA_ch0", "w.v_CamA_ch0", "W_V_CAMA_ch0"):
        _tiffs(d, [base + ".tif", base + "_part0001.tif"])
    _tiffs(d, ["wXv_CamA_ch0_part0001.tif", "W_V_CAMA_ch0_PART0002.tif"])


a3, b3 = make_groups_corpus()

CASES = [
    Case(
        "regexp_once_translation",
        code=MATLAB_REGEXP_ONCE % "regexp",
        args=[[p for p, _ in REGEXP_CORPUS], [s for _, s in REGEXP_CORPUS]],
        python=regexp_once,
    ),
    Case(
        "regexpi_once_translation",
        code=MATLAB_REGEXP_ONCE % "regexpi",
        args=[[p for p, _ in REGEXP_CORPUS], [s.upper() for _, s in REGEXP_CORPUS]],
        python=lambda p, s: regexp_once(p, s, True),
    ),
    Case(
        "channel_patterns_case_sensitive",
        code=MATLAB_CONTAINS,
        args=[NAMES, CHANNEL_PATTERNS, False],
        python=contains_matrix,
    ),
    Case(
        "channel_patterns_ignore_case",
        code=MATLAB_CONTAINS,
        args=[NAMES, CHANNEL_PATTERNS, True],
        python=contains_matrix,
    ),
    Case("grammar_cam", code=MATLAB_GRAMMAR, args=[NAMES[:5] + NAMES[17:19]], python=grammar),
    Case(
        "grammar_cam_odd_names",
        code=MATLAB_GRAMMAR,
        args=[[NAMES[i] for i in (7, 8, 9, 10, 11, 14, 17, 18)]],
        python=grammar,
    ),
    Case("grammar_no_cam", code=MATLAB_GRAMMAR, args=[[NAMES[5], NAMES[6]]], python=grammar),
    Case(
        "grammar_mixed", code=MATLAB_GRAMMAR, args=[[NAMES[0], NAMES[5], "plain"]], python=grammar
    ),
    Case("grammar_empty_list", code=MATLAB_GRAMMAR, args=[mcell()], py_args=[[]], python=grammar),
    # parseImageFilenames
    _listing_case("listing_tif_two_patterns", "d1", False, ["CamA_ch0", "CamB_ch0"]),
    _listing_case("listing_tif_first_pattern_wins", "d1", False, ["CamA", "ch0", "CamB"]),
    _listing_case("listing_tif_regexp_pattern", "d1", False, ["Cam[AB]_ch0", "ch1"]),
    _listing_case("listing_two_folders", ["d1", "d2", "empty"], False, ["ch0"]),
    _listing_case("listing_zarr", "d1", True, ["CamA", "CamB"]),
    _listing_case("listing_data_formats", "d1", False, ["CamA", "CamB"], [".nd2", ".czi"]),
    _listing_case("listing_empty_folder", "empty", False, ["CamA"]),
    _listing_case("listing_no_match", ["d1", "d2"], False, ["CamC"]),
    _listing_case("listing_invalid_pattern", "d1", False, ["Cam(A"]),
    # groupImageFilenamesByChannels
    _group_case("group_two_channels", a3 + b3, ["CamA", "CamB"]),
    _group_case("group_shuffled", [b3[2], a3[0], b3[0], a3[2], a3[1], b3[1]], ["CamA", "CamB"]),
    _group_case(
        "group_three_channels",
        [tile("A", c, y=i) for i in range(2) for c in (2, 0, 1)],
        ["ch0", "ch1", "ch2"],
    ),
    _group_case("group_given_channels", a3 + b3, ["x", "y"], [1, 1, 1, 2, 2, 2]),
    _group_case("group_given_channels_column", a3 + b3, ["x", "y"], [2, 2, 2, 1, 1, 1]),
    _group_case("group_tie_takes_first", ["a1_c0", "b1_c0", "x_c1", "y_c1"], ["c0", "c1"]),
    _group_case("group_ignore_case", ["s_CAMA_1", "s_camb_1"], ["CamA", "CamB"]),
    _group_case("group_unequal_counts", [a3[0], b3[0], a3[1]], ["CamA", "CamB"]),
    _group_case("group_unequal_counts_2", [b3[0], a3[0], b3[1]], ["CamA", "CamB"]),
    _group_case("group_unmatched_file", [a3[0], "other"], ["CamA", "CamB"], error=True),
    _group_case("group_channel_count", a3, ["CamA", "CamB"], error=True),
    _group_case("group_empty", EMPTY, EMPTY),
    # XR_parseImageFilenames
    _dataset_case("dataset_two_channels", ["ds1"], ["CamA_ch0", "CamB_ch0"]),
    _dataset_case("dataset_pattern_order", ["ds1"], ["CamB_ch0", "CamA_ch0"], flip=True),
    _dataset_case("dataset_setting_files", ["set"], ["CamA_ch0"], parse=True, flip=True),
    _dataset_case("dataset_setting_files_missing", ["ds1"], ["CamA_ch0"], parse=True),
    # images without a setting file keep flipZstack (1249929)
    _dataset_case(
        "dataset_setting_files_missing_flip", ["ds1"], ["CamA_ch0"], parse=True, flip=True
    ),
    _dataset_case("dataset_setting_files_no_file", ["ds1"], ["CamC_ch0"], parse=True),
    _dataset_case("dataset_three_folders_decon", ["ds1", "ds2", "ds3"], ["CamA_ch0"], decon=True),
    _dataset_case("dataset_partial_volumes", ["part"], ["CamA_ch0", "CamB_ch0"], decon=True),
    _dataset_case("dataset_zarr", ["zarr"], ["CamA_ch0"], zarr=True, decon=True),
    _dataset_case("dataset_none_found", ["ds3"], ["CamA_ch0"]),
    _dataset_case("dataset_no_pattern_match", ["ds1", "ds2"], ["CamC_ch0"]),
    _dataset_case("dataset_names_as_expressions", ["sp"], ["CamA_ch0"], make=make_special),
    _dataset_case("dataset_names_as_expressions_ic", ["sp"], ["_ch0"], make=make_special),
    _dataset_case("dataset_regexp_pattern", ["ds1", "ds2"], ["Cam[AB]_ch0_.*003x"]),
    _dataset_case("dataset_folder_in_pattern", ["ds1", "ds2"], ["ds2/+Scan.*001x", "ds1/Scan"]),
]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    if case.id == "dataset_folder_in_pattern" and os.name == "nt":
        # groupPartialVolumeFiles matches the patterns against [dataPath, filesep, name]: on
        # Windows the '/' of 'ds2/+Scan' meets '\', so MATLAB there, as the port, keeps nothing
        pytest.skip("MATLAB's filesep join makes a folder in a channel pattern miss on Windows")
    oracle.check(case, tmp_path)


# ---------------------------------------------------------------------------------------------
# unit tests
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("pattern", "text", "match"),
    [
        (r"a\bb", "a\bb", "a\bb"),
        (r"\<b", "a b", "b"),
        (r"ab$", "ab\n", None),
        (r"a.b", "a\nb", "a\nb"),
        (r"a\sb", "a\x1cb", None),
        (r"a\qb", "aqb", "aqb"),
        (r"a\x5fb", "a׻", "a׻"),
        (r"a(?i)_B", "a_b", "a_b"),
        (r"a)", "a)", "a)"),
        (r"*a", "*a", "*a"),
        (r"a**", "aa*", "aa*"),
        (r"a{,1}b", "a{,1}b", "a{,1}b"),
        (r"(?<n>a)\k<n>", "aa", "aa"),
        (r"[a&&b]", "&", "&"),
        (r"[]a]", "a]", None),
        (r"[^]", "x", "x"),
        (r"[x[a]]", "[", "["),
        (r"(?<n>a", "(?<n>a", "<n>a"),
        (r"a(?i)b|C", "xc", None),
    ],
)
def test_matlab_regex_rules(pattern, text, match):
    m = matlab_regex(pattern).search(text)
    assert (m.group() if m else None) == match


@pytest.mark.parametrize(
    ("pattern", "identifier"),
    [
        ("a(", "MATLAB:REGEXP:ErrorUnmatchedParen"),
        ("a[b", "MATLAB:REGEXP:ErrorUnmatchedBracket"),
        ("a{2", "MATLAB:REGEXP:ErrorUnmatchedCurly"),
    ]
    + [(p, "pk5d:regexp") for p in UNSUPPORTED],
)
def test_matlab_regex_errors(pattern, identifier):
    with pytest.raises(MatlabError) as e:
        matlab_regex(pattern, strict=True)
    assert e.value.identifier == identifier


@pytest.mark.parametrize(
    ("pattern", "text", "match"),
    [("a(", "a(", "a("), ("a[b", "a[b", "a[b"), ("a{2", "a{2", "a{2"), ("((a)", "(a", "(a")],
)
def test_regexp_reads_unmatched_brackets_as_text(pattern, text, match):
    assert matlab_regex(pattern).search(text).group() == match


def test_named_tokens_are_python_groups():
    rx = matlab_regex(fnm.FILENAME_EXPRESSION_CAM, True)
    m = rx.search(tile("B", 1, x=-2, it="0007", sub="_0001", laser=560, t=12, prefix="p_"))
    assert m.group("prefix") == "p_"
    assert m.group("Iter") == "0007_0001"
    assert (m.group("Cam"), m.group("ch"), m.group("laser")) == ("B", "1", "560")
    assert (m.group("x"), m.group("y"), m.group("t")) == ("-02", "000", "0012")


def test_channel_pattern_literal_or_regexp():
    names = ["a.zarr_x", "aXzarr", "a(b"]
    assert matches_channel_pattern(names, ".zarr") == [True, True, False]
    assert matches_channel_pattern(["CAMA"], "cama", True) == [True]
    assert matches_channel_pattern(["CAMA"], "cama") == [False]
    # not a valid expression: matched as text (channel_regexp_pattern)
    assert matches_channel_pattern(names, "a(b") == [False, False, True]


def test_parse_image_filenames_lists_in_dir_order(tmp_path):
    make_listing(tmp_path)
    r = parse_image_filenames(str(tmp_path / "d1"), False, ["CamB_ch0", "CamA"])
    assert r.fnames == sorted(r.fnames, key=fnm.matlab_sort_key)
    assert ".hidden_CamA_ch0.tif" in r.fnames
    assert "x_CamA_ch0.tiff" not in r.fnames
    assert r.ch_inds[r.fnames.index(tile("B", 0) + ".tif")] == 0
    assert r.ch_inds[r.fnames.index(tile("A", 1) + ".tif")] == 1
    assert all(f.startswith(str(tmp_path / "d1") + "/") for f in r.filepaths)
    fn, fs, fd, fp, ch = r
    assert fn is r.fnames and fd is r.fd_inds


def test_parse_image_filenames_warns_when_empty(tmp_path):
    with pytest.warns(UserWarning) as rec:
        r = parse_image_filenames([str(tmp_path)], True, ["CamA"])
    assert r.fnames == [] and r.fd_inds.size == 0
    assert len(rec) == 2 and "Zarr" in str(rec[0].message)


def test_group_rows_follow_names():
    a, b = make_groups_corpus()
    g = group_image_filenames_by_channels([b[1], a[0], b[0], a[1]], ["CamA", "CamB"])
    assert g.tolist() == [[3, 0], [1, 2]]
    g = group_image_filenames_by_channels(a + b[:1], ["CamA", "CamB"])
    assert g.tolist() == [[0, 3], [1, -1], [2, -1]]  # -1: no file of the channel left


def test_grammar_selects_expression():
    f = parse_filename_fields([NAMES[0], NAMES[5], "plain_CamA_ch0_x", NAMES[18]])
    assert "Cam" in f[0] and "Cam" not in f[1]
    assert f[2] is None
    assert f[3]["prefix"] == "a"  # the first match only


def _age(paths, minutes):
    t = time.time() - minutes * 60
    for p in paths:
        os.utime(p, (t, t))


def test_dataset_streaming_leaves_recent_files(tmp_path):
    # (MATLAB: AutomaticProcessingTest/streamingLeavesOutRecentFilesWithTheirSizes)
    names = [tile("A", 0, y=i) + ".tif" for i in range(3)]
    _tiffs(tmp_path, names[:1], shape=(16, 16, 16))
    _tiffs(tmp_path, names[1:2], shape=(16, 16, 20))
    _tiffs(tmp_path, names[2:], shape=(16, 16, 24))
    args = ([str(tmp_path) + "/"], ["CamA_ch0"], False, False, False, [])
    full = parse_dataset_filenames(*args, False, 1, False)
    assert full.latest_modify_times.tolist() == [0.0]

    r = parse_dataset_filenames(*args, True, 1, False)
    assert r.fnames == [] and r.gfnames == [] and r.data_sizes.size == 0
    assert r.partialvols.size == 0

    _age([tmp_path / names[2]], 5)
    r = parse_dataset_filenames(*args, True, 1, False)
    assert r.fnames == [names[2]] and r.gfnames == [[names[2]]]
    assert r.data_sizes.tolist() == [full.data_sizes[2]]
    assert r.latest_modify_times[0] < 0.5  # the newest file's age

    # all old: nothing is left out
    _age([tmp_path / n for n in names], 10)
    r = parse_dataset_filenames(*args, True, 1, False)
    assert r.fnames == names


def test_dataset_sizes_use_largest_of_first_and_last(tmp_path):
    names = [tile("A", 0, y=i) + ".tif" for i in range(3)]
    _tiffs(tmp_path, names, shape=(4, 6, 3))  # depths 3, 4, 3
    r = parse_dataset_filenames([str(tmp_path) + "/"], ["CamA"], False, False, False, [], False, 1)
    # MATLAB reads the image size of the first and the last file only, and gives every file
    # the largest of those and of the other files' sizes on disk
    expect = max(4 * 6 * 3 * 2.0, float(os.path.getsize(tmp_path / names[1])))
    assert r.data_sizes.tolist() == [expect] * 3
    assert r.ftp_inds.tolist() == [0] and r.mask_fullpaths == [None]


def test_dataset_partial_volume_groups(tmp_path):
    base = tile("A", 0)
    _tiffs(tmp_path, [base + ".tif", base + "_part0002.tif", base + "_part0001.tif"])
    r = parse_dataset_filenames(
        [str(tmp_path) + "/"], ["CamA"], False, True, True, [str(tmp_path)], False, 1
    )
    assert r.gfnames == [[base + ".tif", base + "_part0001.tif", base + "_part0002.tif"]]
    assert r.partialvols.tolist() == [True]
    assert r.flip_zstack.tolist() == [True]
    assert r.mask_fullpaths == [f"{tmp_path}/Masks/{base}_eroded.zarr"]


def test_regex_cache_is_per_flag():
    assert matlab_regex("a") is matlab_regex("a")
    assert matlab_regex("a", True).flags & re.IGNORECASE
    assert not matlab_regex("a").flags & re.IGNORECASE
