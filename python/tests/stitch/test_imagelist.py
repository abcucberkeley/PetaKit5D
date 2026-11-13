"""Image lists of stitching: readtable, the tile name grammar, row filtering, sortrows order,
the primary channel, several folders (task S01).

Oracle cases compare with MATLAB's readtable, stitch_parse_image_list_information,
stitch_parse_multi_loc_image_list_information, parseTileInfo and stitch_process_filenames on
made-up lists; the unit tests pin the rules found there.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

from petakit5d.stitch import imagelist as il
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, File, InPath, Tree, mcell

# ------------------------------------------------------------------------------------------
# helpers: MATLAB tables as values (names, descriptions, columns)
# ------------------------------------------------------------------------------------------

# T -> {names, descriptions, columns}
TCONV = (
    "tconv = @(T) {T.Properties.VariableNames, T.Properties.VariableDescriptions, "
    "cellfun(@(v) T.(v), T.Properties.VariableNames, 'UniformOutput', false)};\n"
)


def m_col(c: np.ndarray):
    kind = il.column_kind(c)
    if kind == "double":
        return np.asarray(c, dtype=np.float64).reshape(-1, 1)
    if kind == "char":
        return np.asarray(c, dtype="<U1").reshape(-1, 1)
    # a cell column (its size counts where a manifest keeps only the hash)
    out = np.empty((len(c), 1), dtype=object)
    out[:, 0] = list(c)
    return out


def m_table(t: il.Table | None):
    if t is None:
        return None
    return [t.names, list(t.descriptions), [m_col(t[n]) for n in t.names]]


def m_chars(s: str) -> np.ndarray:
    return np.array(list(s), dtype="<U1").reshape(-1, 1)


def write_text(text: str, name: str = "t.csv"):
    def make(path: Path) -> None:
        path.write_bytes(text.encode("utf-8") if isinstance(text, str) else text)

    return File(make, name)


# ------------------------------------------------------------------------------------------
# readtable
# ------------------------------------------------------------------------------------------

READTABLE_FILES = {
    "std": "Filepath,Filename,StageX (um),StageY(um),StageZ_um_,Is_CamB?,Data Image Index\n"
    "/a/b.tif,b.tif,1.5,-2,3e2,0,1\n/a/c.tif,c.tif,NaN,Inf,-Inf,1,2\n",
    "dupnames": "a,a,b,,1x,x y,for,a_1\n1,2,3,4,5,6,7,8\n",
    "dupnames_only": "a,a,b\n1,2,3\n",
    "emptyheader_only": "a,,c\n1,,3\n4,,6\n",
    "emptycells": "n,s,m\n1,,\n,x,\n3,y,\n",
    "quoted": 'n,s\n1,"a,b"\n2,"he said ""hi"""\n3,  spaced  \n',
    "quoted_spaces": 'n,s\n1,"  a b  "\n2," x"\n',
    "mixed": "n,s\n1,2\n2,x\n3,4\n",
    "m12": "n,s\n1,x\n2,y\nx,3\n",
    "m11": "n,s\n1,x\nx,2\n",
    "m23": "n\n1\n2\nx\ny\nz\n",
    "textfirst": "n\n" + "x\n" * 30 + "".join(f"{i}\n" for i in range(1, 41)),
    "textlate": "n\n" + "".join(f"{i}\n" for i in range(1, 41)) + "x\n" * 30,
    "textkeep": "a\nx\n01\ny\n",
    "allempty": "a,b,c\n1,,x\n2,,y\n",
    "trailcomma": "a,b,\n1,2,\n3,4,\n",
    "trail_hdr_empty_data": "a,b,\n1,2,3\n3,4,5\n",
    "trail_hdr_data_empty": "a,b,c\n1,2,\n3,4,\n",
    "alltext": "Filepath,Filename\n/a/b.tif,b.tif\n/a/c.tif,c.tif\n",
    "crlf": "a,b\r\n1,x\r\n2,y\r\n",
    "numfmt": "a,b,c,d,e,f,g\n.5,+3,1E+8,0x1F,-0,1.,1d3\n1,2,3,4,5,6,7\n",
    "upper_nan": "a\nnan\nNAN\n-inf\n+Inf\ninf\n1\n",
    "blanklines": "a,b\n\n1,x\n\n2,y\n",
    "onlyheader": "a,b\n",
    "spaces": "a , b\n 1 , 2 \n",
    "tabs": "a,b\n\t1\t,x\n",
    "unicode": "a,b\n1,été\n",
    "longname": "abcdefghij" * 7 + ",b\n1,2\n",
    "numlike": "a,b\n1,0001\n2,0002\n",
    "numheader": "1,2\n3,4\n",
    "numhdrcol": "Filename,1\nx,2\ny,3\n",
    "textnums": "a,b\nx,1\ny,2\n",
    "qnewline": 'a,b\n1,"x\ny"\n2,z\n',
    "qheader": '"StageX (um)","b,c"\n1,2\n',
    "na": "a,b\n1,N/A\n2,3\n3,4\n",
    "thousands": 'a\n"1,000"\n2\n',
    "quotednum": 'a,b\n"1",2\n"3",4\n',
    "onlycomma": "a,b\n,\n1,2\n",
    "semicolon": "a;b,c\n1;2,3\n",
    "truefalse": "a,b\n1,true\n2,false\n",
    "slicelist": "Filename,Data Image Index,Channel,Is_CamB?,X MCS2 Position (um),"
    "Encoder Timestamp (ms)\nScan_Iter_0000_CamA_ch0_x.tif,0,0,0,490.045755,1.060603E+8\n"
    "Scan_Iter_0000_CamA_ch0_x.tif,1,0,0,490.197084,1.060603E+8\n",
    "writetable_style": 'Filepath,Filename,StageX_um_\n"/a/x y.tif","x y.tif",0.333333333333333\n',
    "bom": b"\xef\xbb\xbfa,b\n1,2\n",
}

READTABLE_CODE = TCONV + "out = tconv(readtable(in{1}, 'Delimiter', 'comma'));"

READTABLE_CASES = [
    Case(
        f"readtable_{k}",
        code=READTABLE_CODE,
        nargout=3,
        args=[InPath("t.csv")],
        python=lambda p: m_table(il.read_table(p)),
        inputs={"t.csv": write_text(v)},
        files=False,
    )
    for k, v in READTABLE_FILES.items()
]

# ------------------------------------------------------------------------------------------
# image lists: tile folders and their lists
# ------------------------------------------------------------------------------------------

HEADER = (
    "Filepath,Filename,StageX (um),StageY (um),StageZ (um),ObjectiveX (um),ObjectiveY (um),"
    "ObjectiveZ (um)"
)


def tile_name(
    cam="A",
    ch=0,
    x=0,
    y=0,
    z=0,
    t=0,
    it="0000_0000",
    stack=0,
    laser=488,
    abst=1000,
    msec=0,
    prefix="",
    suffix="",
    ext=".tif",
):
    camstr = f"_Cam{cam}" if cam else ""
    return (
        f"{prefix}Scan_Iter_{it}{camstr}_ch{ch}_CAM1_stack{stack:04d}_{laser}nm_{msec:07d}msec_"
        f"{abst:010d}msecAbs_{x:03d}x_{y:03d}y_{z:03d}z_{t:04d}t{suffix}{ext}"
    )


def list_text(names, mapped=None, header=HEADER) -> str:
    rows = [header + (",mappedFilename" if mapped is not None else "")]
    for i, n in enumerate(names):
        v = [0.5 * i, -1.25 * i, 0.1 * (i % 3), 0.5 * i, -1.25 * i, 0.1 * (i % 3)]
        r = f"/mnt/projects/embryo_2024/{n},{n}," + ",".join(f"{x:.6f}" for x in v)
        if mapped is not None:
            r += f",{mapped[i]}"
        rows.append(r)
    return "\n".join(rows) + "\n"


def folder(listed, present=None, mapped=None, extra=(), zarr=False, sub="", header=HEADER):
    """A tile folder: list.csv with the listed names, (empty) tile files for present ones."""
    present = listed if present is None else present

    def make(root: Path) -> None:
        root.mkdir(parents=True, exist_ok=True)
        (root / "list.csv").write_text(list_text(listed, mapped, header), newline="\n")
        d = root / sub if sub else root
        d.mkdir(parents=True, exist_ok=True)
        for n in list(present) + list(extra):
            if zarr:
                (d / n.replace(".tif", ".zarr")).mkdir()
                (d / n.replace(".tif", ".zarr") / ".zarray").write_text("{}", newline="\n")
            else:
                (d / n).write_bytes(b"")
        (root / "stitchInfo.mat").write_bytes(b"x")

    return make


def grid(cams=("A",), chs=(0,), xs=(0, 1), ys=(0,), zs=(0,), its=("0000_0000",), **kw):
    out = []
    for it in its:
        for cam in cams:
            for ch in chs:
                for z in zs:
                    for y in ys:
                        for x in xs:
                            abst = 1000 + 37 * x + 101 * y + 7 * z + 3 * ch
                            out.append(tile_name(cam, ch, x, y, z, 0, it, abst=abst, **kw))
    return out


OPT_NAMES = {
    "Streaming": "streaming",
    "stitchInfoFullpath": "stitch_info_fullpath",
    "stitchInfoPath": "stitch_info_path",
    "zarrFile": "zarr_file",
    "channelPatterns": "channel_patterns",
    "useProcessedData": "use_processed_data",
    "ProcessedDirStr": "processed_dir_str",
    "timepoints": "timepoints",
    "subtimepoints": "subtimepoints",
    "xcorrMode": "xcorr_mode",
    "primaryCh": "primary_ch",
    "onlineStitch": "online_stitch",
}

PARSE_OUT = (
    "out = {tconv(tab), pt, fullIter, Ch, Cam, stackn, nz, specifyCam, prefix, zl, xm, sip};"
)
PARSE_CODE = (
    TCONV
    + "opts = in{2};\n"
    + "for k = 2 : 2 : numel(opts), if ischar(opts{k}) && startsWith(opts{k}, '<D>'), "
    + "opts{k} = [in{1}, opts{k}(4 : end)]; end, end\n"
    + "[tab, ptab, fullIter, Ch, Cam, stackn, nz, specifyCam, prefix, zl, xm, sip] = "
    + "stitch_parse_image_list_information(in{1}, [in{1}, 'list.csv'], opts{:});\n"
    + "if istable(ptab), pt = tconv(ptab); else, pt = []; end\n"
    + PARSE_OUT
)


def py_opts(d: str, opts) -> dict:
    kw = {}
    for k, v in zip(opts[0::2], opts[1::2], strict=True):
        if isinstance(v, str) and v.startswith("<D>"):
            v = d + v[3:]
        kw[OPT_NAMES[k]] = tuple(v) if isinstance(v, list) else v
    return kw


def info_out(info: il.ImageListInfo) -> list:
    return [
        m_table(info.tab),
        m_table(info.primary_tab),
        info.full_iter,
        info.ch,
        m_chars(info.cam),
        info.stackn,
        info.nz,
        info.specify_cam,
        info.prefix,
        info.zlayer_stitch,
        info.xcorr_mode,
        info.stitch_info_fullpath,
    ]


def py_parse(d: str, opts) -> list:
    return info_out(il.parse_image_list(d, d + "list.csv", **py_opts(d, opts)))


def parse_case(cid, make, opts=(), error=False, **kw):
    return Case(
        f"parse_{cid}",
        code=PARSE_CODE,
        nargout=12,
        args=[InPath("d"), mcell(*opts)],
        py_args=[InPath("d"), list(opts)],
        python=py_parse,
        inputs={"d": Tree(make, "d")},
        files=False,
        expect_error=error,
        py_error=(il.ImageListError, NotImplementedError),
        **kw,
    )


TWO_CAM = grid(cams=("A", "B"), chs=(0, 1), xs=(0, 1), ys=(0, 1))
NO_SUBITER = [tile_name("A", 0, x, 0, 0, 0, "0012", abst=1000 + x) for x in range(3)]
NOCAM = grid(cams=("",), chs=(0, 1), xs=(0, 1))
MULTI_T = grid(its=("0000_0000", "0001_0000", "0000_0001"), xs=(0, 1))
MIXED_FORMS = [
    tile_name("A", 0, 0, 0, 0, prefix="expA_"),
    tile_name("A", 0, 1, 0, 0, prefix="expA_", suffix="_decon"),
    tile_name("A", 0, 2, 0, 0, it="0000_0000_0003", prefix="expA_"),
    tile_name("A", 0, 3, 0, 0, prefix="expA_").replace("Scan_Iter", "scan_iter"),
    tile_name("A", 0, 0, 1, 0, prefix="expA_").replace("000x", "-01x"),
]
UNMATCHED = TWO_CAM[:4] + ["Scan_Iter_0000_CamA_ch0_bad_name.tif"] + TWO_CAM[4:8]
PARTIAL = grid(xs=(0, 1, 2))
PARTIAL_LIST = PARTIAL + [n.replace(".tif", "_part0001.tif") for n in PARTIAL[:2]]
SPLIT_ABST = [tile_name("A", 0, x, 0, 0, msec=10 * x, abst=2000 + x) for x in range(3)]

PARSE_CASES = [
    parse_case("two_cams_auto_primary", folder(TWO_CAM)),
    parse_case("two_cams_primary_camb_ch1", folder(TWO_CAM), ["primaryCh", "CamB_ch1"]),
    parse_case("primary_lowercase", folder(TWO_CAM), ["primaryCh", "camb_ch1"]),
    parse_case("channel_filter_camb", folder(TWO_CAM), ["channelPatterns", ["CamB"]]),
    parse_case("channel_filter_regexp", folder(TWO_CAM), ["channelPatterns", ["Cam[B]_ch1"]]),
    parse_case("iter_without_subiter", folder(NO_SUBITER)),
    parse_case("no_camera_in_names", folder(NOCAM), ["channelPatterns", ["ch"]]),
    parse_case(
        "no_camera_primary_ch1", folder(NOCAM), ["primaryCh", "ch1", "channelPatterns", ["ch"]]
    ),
    parse_case("time_points", folder(MULTI_T)),
    parse_case("time_points_filter", folder(MULTI_T), ["timepoints", [0]]),
    parse_case(
        "time_points_sub_mismatch",
        folder(MULTI_T),
        ["timepoints", [0], "subtimepoints", [0]],
        error=True,
    ),
    parse_case(
        "time_points_sub_all", folder(MULTI_T), ["timepoints", [0, 1], "subtimepoints", [0, 0]]
    ),
    parse_case("name_forms", folder(MIXED_FORMS)),
    parse_case("unmatched_rows", folder(UNMATCHED)),
    parse_case("partial_files", folder(PARTIAL_LIST, PARTIAL)),
    parse_case("missing_tiles", folder(TWO_CAM, TWO_CAM[2:])),
    parse_case("missing_tiles_streaming", folder(TWO_CAM, TWO_CAM[2:]), ["Streaming", True]),
    parse_case("substring_exists", folder(PARTIAL, [n.replace(".tif", "x.tif") for n in PARTIAL])),
    parse_case("xcorr_all", folder(TWO_CAM), ["xcorrMode", "all"]),
    parse_case("xcorr_primary", folder(TWO_CAM), ["xcorrMode", "primary"]),
    parse_case("stitch_info_given", folder(TWO_CAM), ["stitchInfoFullpath", "<D>stitchInfo.mat"]),
    parse_case(
        "stitch_info_missing", folder(TWO_CAM), ["stitchInfoFullpath", "<D>none.mat"], error=True
    ),
    parse_case("stitch_info_path", folder(TWO_CAM), ["stitchInfoPath", "/mnt/projects/stitchInfo"]),
    parse_case("zarr_tiles", folder(TWO_CAM, zarr=True), ["zarrFile", True]),
    parse_case(
        "online_stitch_zlayers",
        folder(grid(xs=(0, 1), zs=(0, 1, 2))),
        ["Streaming", True, "onlineStitch", True],
    ),
    parse_case(
        "processed_data",
        folder(TWO_CAM, sub="DSR"),
        ["useProcessedData", True, "ProcessedDirStr", "DSR"],
    ),
    parse_case("no_tiles", folder(TWO_CAM, []), error=True),
    parse_case("no_listed_tile_exists", folder(TWO_CAM, [], extra=["other_CamA.tif"]), error=True),
    parse_case("abstime_split", folder(SPLIT_ABST)),
    parse_case("laser_split", folder([tile_name("A", 0, x, laser=488 + 72 * x) for x in range(2)])),
    parse_case("single_tile", folder(TWO_CAM[:1])),
    parse_case("empty_list", folder([], extra=TWO_CAM[:1]), error=True),
    parse_case("unknown_grammar", folder(["tile_a.tif", "tile_b.tif"]), error=True),
    parse_case("primary_bad_format", folder(TWO_CAM), ["primaryCh", "A_0"], error=True),
    parse_case("primary_not_found", folder(TWO_CAM), ["primaryCh", "CamC_ch5"], error=True),
    parse_case(
        "mapped_names",
        folder(["tile_a.tif", "tile_b.tif", "tile_c.tif"], mapped=grid(xs=(0, 1, 2))),
    ),
    parse_case(
        "mapped_names_filter",
        folder(
            ["tile_a.tif", "tile_b.tif", "tile_c.tif"],
            mapped=["x_CamB_ch0_y.tif"] + grid(xs=(1, 0)),
        ),
        ["channelPatterns", ["CamA"]],
    ),
    parse_case("writetable_header", folder(TWO_CAM, header=HEADER.replace(" (um)", "_um_"))),
    parse_case("bad_xcorr_mode", folder(TWO_CAM), ["xcorrMode", "primaryfirst"], error=True),
    # 1,600 tiles (the benchmark's size), listed in reverse order
    parse_case(
        "tiles_1600",
        folder(grid(cams=("B", "A"), chs=(1, 0), xs=range(19, -1, -1), ys=range(19, -1, -1))),
    ),
]

# several folders

MULTI_CODE = (
    TCONV
    + "dps = {[in{1}, 'r1'], [in{1}, 'r2']};\n"
    + "lists = {[in{1}, 'r1/list.csv'], [in{1}, 'r2/list.csv']};\n"
    + "[tab, ptab, fullIter, Ch, Cam, stackn, nz, specifyCam, prefix, zl, xm, sip] = "
    + "stitch_parse_multi_loc_image_list_information(dps, lists, in{2}{:});\n"
    + "if istable(ptab), pt = tconv(ptab); else, pt = []; end\n"
    + PARSE_OUT
)


def two_folders(a_names, b_names, a_present=None, b_present=None):
    def make(root: Path) -> None:
        folder(a_names, a_present)(root / "r1")
        folder(b_names, b_present)(root / "r2")

    return make


def py_multi(d: str, opts) -> list:
    dps = [d + "r1", d + "r2"]
    lists = [d + "r1/list.csv", d + "r2/list.csv"]
    out = info_out(il.parse_multi_loc_image_list(dps, lists, **py_opts(d, opts)))
    # unique of the prefixes' char matrix is a column of characters
    out[8] = m_chars(out[8]) if out[8] else ""
    return out


def multi_case(cid, make, opts=(), error=False):
    return Case(
        f"multi_{cid}",
        code=MULTI_CODE,
        nargout=12,
        args=[InPath("d"), mcell(*opts)],
        py_args=[InPath("d"), list(opts)],
        python=py_multi,
        inputs={"d": Tree(make, "d")},
        files=False,
        expect_error=error,
        py_error=(il.ImageListError, NotImplementedError),
    )


R1 = grid(cams=("A", "B"), chs=(1,), xs=(0, 1))
R2 = grid(cams=("B",), chs=(0, 1), xs=(0, 1))
MULTI_CASES = [
    multi_case("two_folders", two_folders(R1, R2)),
    multi_case("two_folders_missing", two_folders(R1, R2, R1[1:], R2[:-1])),
    multi_case("two_folders_xcorr_all", two_folders(R1, R2), ["xcorrMode", "all"]),
    multi_case("prefix_chars", two_folders(grid(prefix="ab_"), grid(prefix="ba_"))),
    multi_case("mixed_cam_forms", two_folders(grid(), NOCAM), ["channelPatterns", ["ch"]]),
]

# parseTileInfo and stitch_process_filenames

TILE_NAMES = [
    "Tile_ch0_000x_001y_002z_0003t",
    "Tile_ch1_010x_000y_000z_0000t",
    "Tile_CamB_ch1_002x_003y_004z_0005t",
    "tile_camA_c2_001X_001Y_001Z_0001T",
    "Tile_ch0_noindex",
]


def py_tile_info(names, patterns, cam):
    info, cams = il.parse_tile_info(names, patterns)
    return [info, cams if cam else None]


CAM_NAMES = [
    "Tile_CamB_ch1_002x_003y_004z_0005t",
    "Tile_cama_ch0_000x_001y_002z_0003t",
    "Tile_CamB_ch0_001x_001y_001z_0001t",
]


def tile_info_case(cid, names, pats, cam=False, error=False):
    return Case(
        f"parse_tile_info_{cid}",
        matlab="parseTileInfo",
        args=[names, list(pats)],
        nargout=2,
        python=lambda n, p: py_tile_info(n, p, cam),
        expect_error=error,
        py_error=il.ImageListError,
    )


DEFAULT_PATTERNS = ("0000t", "ch0", "000x", "000y", "000z")
TILE_INFO_CASES = [
    tile_info_case("default", TILE_NAMES, DEFAULT_PATTERNS),
    tile_info_case("cam", CAM_NAMES, ("0000t", "CamA_ch0", "000x", "000y", "000z"), cam=True),
    tile_info_case(
        "cam_missing", TILE_NAMES, ("0000t", "CamA_ch0", "000x", "000y", "000z"), error=True
    ),
    tile_info_case("prefix_patterns", TILE_NAMES, ("t0", "c0", "x0", "y0", "z0")),
    tile_info_case("skip_time", TILE_NAMES, ("", "ch0", "000x", "000y", "000z")),
    tile_info_case("no_letters", TILE_NAMES, ("0000", "ch0", "000x", "000y", "000z"), error=True),
]

TILES = ["/mnt/projects/run1/a.tif", "/mnt/projects/run1//b.tif"]
PROCESS_CASES = [
    Case(
        f"process_filenames_{cid}",
        matlab="stitch_process_filenames",
        args=[TILES, *a],
        nargout=4,
        python=lambda *x: list(il.process_filenames(*x)),
    )
    for cid, a in [
        ("plain", ["", False, [], False, True]),
        ("dsr", ["DSR", False, [], False, True]),
        ("zarr_input", ["", False, [], True, True]),
        ("resample_1", ["", False, [2], False, True]),
        ("resample_2", ["", False, [2, 3], False, True]),
        ("resample_3", ["DSR", False, [1, 2, 0.5], False, True]),
        ("resample_ones", ["", False, [1, 1, 1], False, True]),
        ("mip", ["", [False, False, True], [], False, True]),
        ("no_processing", ["", False, [], False, False]),
    ]
]

ORACLE_CASES = READTABLE_CASES + PARSE_CASES + MULTI_CASES + TILE_INFO_CASES + PROCESS_CASES


@oracle.parametrize(ORACLE_CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ------------------------------------------------------------------------------------------
# unit tests
# ------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "header,name",
    [
        ("StageX (um)", "StageX_um_"),
        ("StageX(um)", "StageX_um_"),
        ("Is_CamB?", "Is_CamB_"),
        ("Data Image Index", "DataImageIndex"),
        ("x y", "xY"),
        ("1x", "x1x"),
        ("for", "xFor"),
        ("_a", "x_a"),
        ("a" * 70, "a" * 63),
    ],
)
def test_make_valid_name(header, name):
    assert il.make_valid_name(header) == name


def test_read_table_renames_and_warns(tmp_path):
    p = tmp_path / "l.csv"
    p.write_text("Filename,StageX (um)\na.tif,1.5\n", newline="\n")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        t = il.read_table(p)
    assert t.names == ["Filename", "StageX_um_"]
    assert t.descriptions == ["Filename", "StageX (um)"]
    assert [getattr(x.message, "identifier", "") for x in w] == [
        "MATLAB:table:ModifiedAndSavedVarnames"
    ]


def test_read_table_majority_and_leading_text(tmp_path):
    p = tmp_path / "l.csv"
    p.write_text("n,s\n1,x\nx,2\n", newline="\n")
    t = il.read_table(p)
    assert t.kind("n") == "double" and t.kind("s") == "double"  # ties are numbers
    p.write_text("a\nx\n01\ny\n", newline="\n")
    t = il.read_table(p)
    assert t.height == 2 and t["a"][0] == 1 and np.isnan(t["a"][1])


def test_read_table_nan_is_matlab_nan(tmp_path):
    p = tmp_path / "l.csv"
    p.write_text("a,b\n1,\nNaN,2\n", newline="\n")
    t = il.read_table(p)
    bits = t["a"].view(np.uint64)[1], t["b"].view(np.uint64)[0]
    assert all(int(b) == 0xFFF8000000000000 for b in bits)


def test_read_table_unsupported(tmp_path):
    p = tmp_path / "l.csv"
    p.write_text("a,b\n2024-01-02,x\n2024-01-03,y\n", newline="\n")
    with pytest.raises(NotImplementedError):
        il.read_table(p)
    p.write_text("a,b,c\n1,2\n3,4,5\n", newline="\n")
    with pytest.raises(NotImplementedError):
        il.read_table(p)


def test_write_table_format(tmp_path):
    t = il.Table(
        {
            "s": np.array(["plain", "with,comma", 'q"x', "", None], dtype=object),
            "x": np.array([0.1 + 0.2, 1 / 3, -0.0, float("nan"), 2.0**53 + 2]),
        }
    )
    p = tmp_path / "w.csv"
    il.write_table(t, p)
    assert p.read_text() == (
        's,x\n"plain",0.3\n"with,comma",0.333333333333333\n"q""x",-0\n"",NaN\n'
        ",9.00719925474099e+15\n"
    )


def test_sortrows_stable_nan_last():
    t = il.Table(
        {
            "fi": ["0001", "0000", "0001", "0000"],
            "camera": np.array(["B", "A", "A", "B"], dtype="<U1"),
            "ch": np.array([2.0, 1.0, 1.0, np.nan]),
        }
    )
    _, idx = il.sortrows(t, ["camera", "ch", "fi"])
    assert idx.tolist() == [1, 2, 0, 3]


def test_iter_without_subiter_quirk(tmp_path):
    folder(NO_SUBITER)(tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        info = il.parse_image_list(str(tmp_path), str(tmp_path / "list.csv"))
    # MATLAB's grammar puts the last digit into subIter when there is no sub-iteration
    assert info.tab["Iter"].tolist() == [1.0] * 3
    assert list(info.tab["subIter"]) == ["2"] * 3
    assert info.full_iter == ["0012"]


def test_primary_channel_first(tmp_path):
    folder(TWO_CAM)(tmp_path)
    info = il.parse_image_list(str(tmp_path), str(tmp_path / "list.csv"), primary_ch="CamB_ch1")
    assert info.ch.tolist() == [1.0, 0.0]
    assert info.cam == "BA"
    assert info.primary_tab.height == 4
    assert set(info.primary_tab["camera"].tolist()) == {"B"}


def test_missing_tiles_dropped_and_warned(tmp_path):
    folder(TWO_CAM, TWO_CAM[2:])(tmp_path)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        info = il.parse_image_list(str(tmp_path), str(tmp_path / "list.csv"))
    assert info.tab.height == len(TWO_CAM) - 2
    assert any("do not exist" in str(x.message) for x in w)


def test_multi_folder_did(tmp_path):
    two_folders(R1, R2)(tmp_path)
    d = str(tmp_path) + "/"
    info = il.parse_multi_loc_image_list(
        [d + "r1", d + "r2"], [d + "r1/list.csv", d + "r2/list.csv"]
    )
    assert sorted(set(info.tab["did"].tolist())) == [1.0, 2.0]
    assert info.ch.tolist() == [1.0, 0.0]


def test_table_selection():
    t = il.Table({"a": np.arange(4.0), "s": ["w", "x", "y", "z"]})
    assert t[t["a"] > 1]["s"].tolist() == ["y", "z"]
    assert t[[3, 0]]["a"].tolist() == [3.0, 0.0]
    assert t[1:3].height == 2
    with pytest.raises(il.ImageListError):
        t["missing"]
