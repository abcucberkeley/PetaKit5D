"""Image list generation (XR_generate_image_list_wrapper and its generators, task S01).

Ports ImageListGenerationTest in intent (the image list parts; the stitching that three of
its methods run on the generated lists belongs to the stitching wrapper's tests), and compares
the generated csv files, encoder_info files and step sizes with MATLAB's for made-up tile
folders, slice lists and acquisition databases.
"""

from __future__ import annotations

import sqlite3
import warnings
from pathlib import Path

import numpy as np
import pytest

from petakit5d.stitch import imagelist as il
from petakit5d.stitch import imagelist_generate as ig
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, InPath, Tree

XY_PIXEL_SIZE = 0.108


def matlab_wrapper(*args):
    from petakit5d import matlab

    matlab.XR_generate_image_list_wrapper(*args)
    return []


WRAPPER_CODE = "XR_generate_image_list_wrapper(in{:}); out = {};"


def wrapper_case(cid, make, *args, error=False, marks=()):
    return Case(
        cid,
        code=WRAPPER_CODE,
        nargout=0,
        args=[InPath("d"), *args],
        python=matlab_wrapper,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=(il.ImageListError, NotImplementedError, ValueError, TypeError),
        marks=marks,
    )


# ------------------------------------------------------------------------------------------
# inputs: tile folders, slice lists, databases
# ------------------------------------------------------------------------------------------


def write_tile(path: Path, shape, value=1) -> None:
    from petakit5d.io import write_tiff

    write_tiff(np.full(shape, value, dtype=np.uint16, order="F"), str(path), group_write=False)


def tiles(names, shape=(48, 64)):
    def make(root: Path) -> None:
        for n in names:
            sh = shape(n) if callable(shape) else shape
            write_tile(root / n, sh)

    return make


def llsm_name(cam, ch, x, y, z, t=0):
    camstr = f"_Cam{cam}" if cam else ""
    return (
        f"Scan_Iter_0000{camstr}_ch{ch}_CAM1_stack0000_642nm_0000000msec_0000000000msecAbs_"
        f"{x:03d}x_{y:03d}y_{z:03d}z_{t:04d}t.tif"
    )


GRID_2X2 = [llsm_name("B", 0, x, y, 0) for x in range(2) for y in range(2)]
GRID_3D = [llsm_name("A", c, x, 0, z) for c in range(2) for x in range(2) for z in range(2)]
TWO_PHOTON = [f"Tile_ch0_{x:03d}x_000y_{z:03d}z_0000t.tif" for x in range(2) for z in range(2)]
NO_CHANNEL = [f"Tile_c0_{x:03d}x_000y_000z_0000t.tif" for x in range(3)]
CAM_TILES = [f"Tile_Cam{c}_ch0_{x:03d}x_000y_000z_0000t.tif" for c in "AB" for x in range(2)]

SLICE_HEADER = (
    "Filename,Channel,Is_CamB_,XFPGAPosition_um_,YFPGAPosition_um_,ZFPGAPosition_um_,"
    "XMCS2Position_um_,YMCS2Position_um_,ZMCS2Position_um_"
)
SLICE_NAME = (
    "Scan_Iter_0000_0000_Cam{cam}_ch0_CAM1_stack0000_488nm_0000000msec_{abst:010d}msecAbs_"
    "000x_{y:03d}y_000z_0000t.tif"
)


def slice_name(i, cam="A"):
    return SLICE_NAME.format(cam=cam, abst=1000 * (i + 1), y=i)


def write_slice_lists(root: Path, mcs2, cams=None, nslice=3, step=0.2, dup=(0, 0), skip=()):
    """One SliceList csv per tile (the stage moves in x); dup: repeated first/last slices."""
    mcs2 = np.asarray(mcs2, dtype=np.float64)
    for i in range(mcs2.shape[0]):
        cam = cams[i] if cams else "A"
        fn = slice_name(i, cam)
        n = nslice - (1 if i in skip else 0)
        lines = [SLICE_HEADER]
        for s in range(1, n + 1):
            k = min(max(s - 1 - dup[0], 0), n - 1 - dup[0] - dup[1])
            x = mcs2[i, 0] + step * k
            lines.append(
                f"{fn},0,{1 if cam == 'B' else 0},{x + 0.2:f},{mcs2[i, 1]:f},{mcs2[i, 2]:f},"
                f"{x:f},{mcs2[i, 1]:f},{mcs2[i, 2]:f}"
            )
        (root / f"SliceList_tile{i + 1}_part{i:04d}.csv").write_text("\n".join(lines) + "\n")


def slice_tree(mcs2, **kw):
    def make(root: Path) -> None:
        write_slice_lists(root, mcs2, **kw)

    return make


MCS2 = [[490.05, -8052.12, -5473.28], [490.04, -7862.72, -5473.27]]

DB_COLUMNS = (
    ("ImageList_id", "INTEGER"),
    ("Filename", "TEXT"),
    ("Camera", "INTEGER"),
    ("Channel", "INTEGER"),
    ("Laser_nm", "INTEGER"),
    ("File_Stack_Index", "INTEGER"),
    ("Elapsed_ms", "REAL"),
    ("Absolute_ms", "INTEGER"),
    ("Tile_index_X", "INTEGER"),
    ("Tile_index_Y", "INTEGER"),
    ("Tile_index_Z", "INTEGER"),
    ("Tile_index_T", "INTEGER"),
    ("Frame_Type", "TEXT"),
) + tuple((f"Mean_{a}_{k}_Position_um", "REAL") for k in ("MCS2", "Target", "FPGA") for a in "XYZ")


def write_db(path: Path, rows) -> None:
    """A made-up acquisition database: the view's columns as a table of the same name."""
    con = sqlite3.connect(path)
    cols = ", ".join(f"{n} {t}" for n, t in DB_COLUMNS)
    con.execute(f"CREATE TABLE ImageListView_NoNull ({cols})")
    con.executemany(
        f"INSERT INTO ImageListView_NoNull VALUES ({', '.join('?' * len(DB_COLUMNS))})", rows
    )
    con.commit()
    con.close()
    # the library's version numbers in the header (offsets 92 to 99) made fixed
    b = bytearray(path.read_bytes())
    b[92:100] = bytes(8)
    path.write_bytes(bytes(b))


def db_rows(cams=(0, 1), chans=((0, 488), (1, 560)), xs=(0, 1), start=1):
    rows = []
    k = start
    for cam in cams:
        for ch, laser in chans:
            for x in xs:
                name = f"acq_Cam{'AB'[cam]}_ch{ch}_{laser}nm_{x:03d}x.tif"
                pos = [100.0 * x + 0.125, -20.5, 3.0 + 0.1 * k, 101.0 * x, -21.0, 3.5, 0.0, 0, 0]
                rows.append((k, name, cam, ch, laser, 0, 0.0, 5000 + 7 * k, x, 0, 0, 0, "", *pos))
                k += 1
    return rows


def db_tree(*dbs):
    def make(root: Path) -> None:
        for i, rows in enumerate(dbs):
            write_db(root / f"Scan_Iter_{i:04d}_slicelist.sqlite3", rows)

    return make


def ism_tree(root: Path) -> None:
    header = (
        "Filename,Data Image Index,X MCS2 Position (um),Y MCS2 Position (um),"
        "Z MCS2 Position (um),X Target Position (um),Y Target Position (um),"
        "Z Target Position (um)"
    )
    for i in range(2):
        lines = [header]
        for s in range(3):
            fn = f"Scan_Iter_0000_0000ISMslice_CamA_ch0_{i:03d}x.tif"
            v = [10.0 * i + 0.1 * s, -3.25, 7.0, 10.0 * i, -3.0, 7.5]
            lines.append(f"{fn},{s}," + ",".join(f"{x:.6f}" for x in v))
        (root / f"SliceList_{i}.csv").write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------------------------------------
# oracle cases
# ------------------------------------------------------------------------------------------

TILE_LIST_FNS = ["tile_a.tif", "tile_b.tif", "tile_c.tif", "tile_d.tif"]
TILE_LIST_IND = np.array(
    [[0, 0, 1, 0, 2], [0, 0, 2, 0, 2], [0, 0, 1, 1, 2], [0, 1, 2, 1, 3]], dtype=np.float64
)


def tile_list_args(ind=TILE_LIST_IND, interval=(50, 40, 10), fns=TILE_LIST_FNS):
    return [
        "tile_list",
        "tileFilenames",
        list(fns),
        "tileIndices",
        np.asarray(ind, dtype=np.float64),
        "tileInterval",
        list(interval),
        "uuid",
        "fixed",
    ]


def empty(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)


TP = ["xyPixelSize", XY_PIXEL_SIZE, "uuid", "fixed"]

GENERATE_CASES = [
    # tile_list
    wrapper_case("generate_tile_list", empty, *tile_list_args()),
    wrapper_case(
        "generate_tile_list_times_negative",
        empty,
        *tile_list_args(
            [[1, 0, -1, 0, 0], [0, 0, 0, 0, 0], [1, 1, 3, 2, 1]], (0.3, 0.7, 1.1), ["a", "b", "c"]
        ),
    ),
    wrapper_case(
        "generate_tile_list_rows_mismatch", empty, *tile_list_args(TILE_LIST_IND[:3]), error=True
    ),
    # tile_position
    wrapper_case(
        "generate_tile_position_ioscan_pixel",
        tiles(GRID_2X2),
        "tile_position",
        "channelPatterns",
        ["Scan"],
        "overlapSize",
        [8, 8, 1],
        "overlapSizeType",
        "pixel",
        "IOScan",
        True,
        *TP,
    ),
    wrapper_case(
        "generate_tile_position_objective_scan",
        tiles(GRID_3D, (16, 24, 10)),
        "tile_position",
        "channelPatterns",
        ["CamA"],
        "overlapSize",
        [4, 4, 2],
        "objectiveScan",
        True,
        "dz",
        0.5,
        *TP,
    ),
    wrapper_case(
        "generate_tile_position_ds",
        tiles(GRID_3D, (16, 24, 10)),
        "tile_position",
        "channelPatterns",
        ["CamA_ch1"],
        "overlapSize",
        [4, 4, 2],
        "DS",
        True,
        *TP,
    ),
    wrapper_case(
        "generate_tile_position_dsr",
        tiles(GRID_3D, (16, 24, 10)),
        "tile_position",
        "channelPatterns",
        ["CamA"],
        "DSR",
        True,
        "skewAngle",
        31.5,
        *TP,
    ),
    wrapper_case(
        "generate_tile_position_skewed",
        tiles(GRID_3D, (16, 24, 10)),
        "tile_position",
        "channelPatterns",
        ["CamA"],
        "overlapSize",
        [4, 2, 3],
        *TP,
    ),
    wrapper_case(
        "generate_tile_position_um_two_photon",
        tiles(TWO_PHOTON, (20, 20, 18)),
        "tile_position",
        "channelPatterns",
        ["Tile"],
        "tilePatterns",
        ["0000t", "ch0", "000x", "000y", "000z"],
        "overlapSize",
        [8 * XY_PIXEL_SIZE, 8 * XY_PIXEL_SIZE, 4 * 0.5],
        "overlapSizeType",
        "um",
        "IOScan",
        True,
        "dz",
        0.5,
        *TP,
    ),
    wrapper_case(
        "generate_tile_position_no_standard_channel",
        tiles(NO_CHANNEL, (20, 40, 8)),
        "tile_position",
        "channelPatterns",
        ["Tile"],
        "tilePatterns",
        ["0000t", "c0", "000x", "000y", "000z"],
        "overlapSize",
        [8, 8, 0],
        "IOScan",
        True,
        "dz",
        0.5,
        *TP,
    ),
    wrapper_case(
        "generate_tile_position_cam_pattern",
        tiles(CAM_TILES[:2], (20, 40, 8)),
        "tile_position",
        "channelPatterns",
        ["Tile"],
        "tilePatterns",
        ["0000t", "CamA_ch0", "000x", "000y", "000z"],
        "IOScan",
        True,
        *TP,
    ),
    # the channel pattern CamA_ch0 finds no channel in CamB names: rows MATLAB leaves unset
    wrapper_case(
        "generate_tile_position_cam_pattern_other_camera",
        tiles(CAM_TILES, (20, 40, 8)),
        "tile_position",
        "channelPatterns",
        ["Tile"],
        "tilePatterns",
        ["0000t", "CamA_ch0", "000x", "000y", "000z"],
        "IOScan",
        True,
        *TP,
        error=True,
    ),
    wrapper_case(
        "generate_tile_position_cam_without_index",
        tiles(CAM_TILES, (20, 40, 8)),
        "tile_position",
        "channelPatterns",
        ["Tile"],
        "IOScan",
        True,
        *TP,
        error=True,
    ),
    wrapper_case(
        "generate_tile_position_missing_index",
        tiles([*NO_CHANNEL, "Tile_c0_003x_000z_0000t.tif"], (20, 40, 8)),
        "tile_position",
        "channelPatterns",
        ["Tile"],
        "tilePatterns",
        ["0000t", "c0", "000x", "000y", "000z"],
        "IOScan",
        True,
        *TP,
        error=True,
    ),
    wrapper_case(
        "generate_tile_position_mode_size",
        tiles(GRID_2X2, lambda n: (40, 64) if "001x_001y" in n else (48, 64)),
        "tile_position",
        "channelPatterns",
        ["CamB"],
        "IOScan",
        True,
        *TP,
    ),
    # encoder
    wrapper_case(
        "generate_encoder", slice_tree(MCS2), "encoder", "dz", 0.2, "channelPatterns", ["CamA_ch0"]
    ),
    wrapper_case(
        "generate_encoder_two_cameras",
        slice_tree(MCS2 + [[1.5, 2.5, -3.5]], cams=["A", "B", "B"]),
        "encoder",
        "dz",
        0.2,
    ),
    # sqlite
    wrapper_case("generate_sqlite", db_tree(db_rows()), "sqlite"),
    wrapper_case("generate_sqlite_no_mapping", db_tree(db_rows()), "sqlite", "mapTilename", False),
    wrapper_case(
        "generate_sqlite_channel_filter",
        db_tree(db_rows(), db_rows(cams=(0,), xs=(2,), start=20)),
        "sqlite",
        "channelPatterns",
        ["CamA_ch1"],
    ),
    wrapper_case("generate_sqlite_empty_db", db_tree(db_rows(), []), "sqlite"),
    wrapper_case("generate_unknown_method", empty, "other", error=True),
]


def py_estimate(d, *args):
    kw = {il_name(k): v for k, v in zip(args[0::2], args[1::2], strict=True)}
    return ig.estimate_actual_step_size_from_encoder(d, **kw)


_EST_NAMES = {
    "dz": "dz",
    "CoordType": "coord_type",
    "channelPatterns": "channel_patterns",
    "skipMissingFile": "skip_missing_file",
    "onlyFirstRow": "only_first_row",
    "xThresh": "x_thresh",
}


def il_name(k):
    return _EST_NAMES[k]


def estimate_case(cid, make, *args, error=False):
    return Case(
        f"estimate_{cid}",
        matlab="XR_estimate_actual_step_size_from_encoder",
        args=[InPath("d"), *args],
        python=py_estimate,
        inputs={"d": Tree(make, "d")},
        expect_error=error,
        py_error=il.ImageListError,
    )


THREE = [[10.0, 0.0, 1.0], [10.0, 150.0, 1.0], [-4.0, 300.0, 1.0]]
ESTIMATE_CASES = [
    estimate_case("full", slice_tree(THREE, nslice=12, dup=(2, 1)), "dz", 0.2),
    estimate_case("mcs2", slice_tree(THREE, nslice=9, dup=(1, 2)), "dz", 0.2, "CoordType", "MCS2"),
    estimate_case("negative_dz", slice_tree(THREE, nslice=10, step=-0.25, dup=(1, 1)), "dz", -0.25),
    estimate_case("skipped_tile", slice_tree(THREE, nslice=8, skip=(1,)), "dz", 0.2),
    estimate_case(
        "missing_files_checked",
        slice_tree(THREE, nslice=6),
        "dz",
        0.2,
        "skipMissingFile",
        False,
        error=True,
    ),
    estimate_case("first_row", slice_tree(THREE, nslice=5), "dz", 0.2, "onlyFirstRow", True),
]

OTHER_CASES = [
    Case(
        f"ism_{pt}",
        code="stitch_generate_imagelist_from_encoder_ISM(in{:}); out = {};",
        nargout=0,
        args=[InPath("d"), pt],
        python=lambda d, p: (ig.generate_image_list_from_encoder_ism(d, p), [])[1],
        inputs={"d": Tree(ism_tree, "d")},
    )
    for pt in ("MCS2", "target")
] + [
    Case(
        "sqlite_old_format",
        code="stitch_generate_imagelist_from_sqlite_old_format(in{:}); out = {};",
        nargout=0,
        args=[InPath("d")],
        python=lambda d: (ig.generate_image_list_from_sqlite_old_format(d), [])[1],
        inputs={"d": Tree(db_tree(db_rows(), db_rows(cams=(1,), start=9)), "d")},
    )
]

ORACLE_CASES = GENERATE_CASES + ESTIMATE_CASES + OTHER_CASES


@oracle.parametrize(ORACLE_CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ------------------------------------------------------------------------------------------
# ImageListGenerationTest, in intent
# ------------------------------------------------------------------------------------------


def read_list(path) -> il.Table:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return il.read_table(path)


def test_tileListCoordinates(tmp_path):
    d = str(tmp_path) + "/"
    ig.generate_image_list(
        d,
        "tile_list",
        tile_filenames=TILE_LIST_FNS,
        tile_indices=TILE_LIST_IND,
        tile_interval=(50, 40, 10),
    )
    tab = read_list(d + "ImageList_from_tile_list.csv")
    assert tab.names == [
        "Filepath",
        "Filename",
        "mappedFilename",
        "StageX_um_",
        "StageY_um_",
        "StageZ_um_",
        "ObjectiveX_um_",
        "ObjectiveY_um_",
        "ObjectiveZ_um_",
    ]
    assert list(tab["Filename"]) == TILE_LIST_FNS
    # positions relative to the minimum tile index times the interval
    idx = TILE_LIST_IND[:, 2:5]
    expected = (idx - idx.min(axis=0)) * np.array([50, 40, 10])
    for prefix in ("Stage", "Objective"):
        got = np.column_stack([tab[f"{prefix}{a}_um_"] for a in "XYZ"])
        np.testing.assert_allclose(got, expected, atol=1e-9)
    assert tab["mappedFilename"][3] == (
        "Scan_Iter_0000_CamA_ch1_CAM1_stack0000_488nm_000000msec_000000msecAbs_002x_001y_003z_0000t.tif"
    )


def _check_encoder_list(folder: Path, mcs2) -> None:
    tab = read_list(folder / "ImageList_from_encoder.csv")
    names = [slice_name(i) for i in range(len(mcs2))]
    order = np.argsort(np.array(list(tab["Filename"]), dtype=object))
    assert [tab["Filename"][i] for i in order] == sorted(names)
    got = np.column_stack([tab[f"Stage{a}_um_"] for a in "XYZ"])[order]
    np.testing.assert_allclose(got, np.asarray(mcs2)[np.argsort(names)], atol=1e-6)


def test_encoderSliceListCoordinates(tmp_path):
    # stage coordinates come from the first MCS2 position of each tile in the SliceList files
    write_slice_lists(tmp_path, MCS2)
    ig.generate_image_list(str(tmp_path), "encoder", dz=0.2, channel_patterns=("CamA_ch0",))
    _check_encoder_list(tmp_path, MCS2)
    assert (tmp_path / "encoder_info" / "coordinates_info.mat").is_file()


def test_encoderListsForMultipleDatasets(tmp_path):
    # one image list per dataset folder (a cell of data paths was passed on as is and errored)
    mcs2 = [[[10, 20, 30], [10, 220, 30]], [[-5, 0, 2], [-5, 150, 2]]]
    folders = [tmp_path / "a", tmp_path / "b"]
    for f, m in zip(folders, mcs2, strict=True):
        f.mkdir()
        write_slice_lists(f, m)
    ig.generate_image_list([str(f) for f in folders], "encoder", channel_patterns=("CamA_ch0",))
    for f, m in zip(folders, mcs2, strict=True):
        _check_encoder_list(f, m)


def test_tilePositionStepsForNonSquareTiles(tmp_path):
    # the x and y steps used to be swapped for non-square tiles (48 x 64, y x)
    tiles(GRID_2X2)(tmp_path)
    ig.generate_image_list(
        str(tmp_path),
        "tile_position",
        xy_pixel_size=XY_PIXEL_SIZE,
        channel_patterns=("Scan",),
        overlap_size=(8, 8, 1),
        overlap_size_type="pixel",
        io_scan=True,
    )
    tab = read_list(tmp_path / "ImageList_from_tile_positions.csv")
    assert tab["StageX_um_"].max() == pytest.approx((64 - 8) * XY_PIXEL_SIZE, abs=1e-9)
    assert tab["StageY_um_"].max() == pytest.approx((48 - 8) * XY_PIXEL_SIZE, abs=1e-9)


def test_tilePositions2DPhase(tmp_path):
    # 2D tiles with standard names, IO scan and overlap in pixels (the list part; the
    # stitching of the tiles is the stitching wrapper's test)
    step = 56
    names = [llsm_name("B", 0, x, y, 0) for x in range(2) for y in range(2)]
    tiles(names, (64, 64))(tmp_path)
    ig.generate_image_list(
        str(tmp_path),
        "tile_position",
        xy_pixel_size=XY_PIXEL_SIZE,
        channel_patterns=("Scan",),
        overlap_size=(8, 8, 1),
        overlap_size_type="pixel",
        objective_scan=False,
        io_scan=True,
    )
    tab = read_list(tmp_path / "ImageList_from_tile_positions.csv")
    idx = np.array(
        [
            [int(n[n.index("x_") - 3 : n.index("x_")]), int(n[n.index("y_") - 3 : n.index("y_")])]
            for n in tab["Filename"]
        ]
    )
    got = np.column_stack([tab["StageX_um_"], tab["StageY_um_"], tab["StageZ_um_"]])
    expected = np.column_stack([idx * step * XY_PIXEL_SIZE, np.zeros(4)])
    np.testing.assert_allclose(got, expected, atol=1e-9)


def test_tilePositions3DNonStandardNames(tmp_path):
    # 2-photon style tiles with the overlap in um get mapped standard names
    tiles(TWO_PHOTON, (40, 40, 36))(tmp_path)
    ig.generate_image_list(
        str(tmp_path),
        "tile_position",
        xy_pixel_size=XY_PIXEL_SIZE,
        dz=0.5,
        channel_patterns=("Tile",),
        tile_patterns=("0000t", "ch0", "000x", "000y", "000z"),
        overlap_size=(8 * XY_PIXEL_SIZE, 8 * XY_PIXEL_SIZE, 4 * 0.5),
        overlap_size_type="um",
        objective_scan=False,
        io_scan=True,
    )
    tab = read_list(tmp_path / "ImageList_from_tile_positions.csv")
    assert "mappedFilename" in tab
    assert all(
        s.startswith("Scan_Iter_0000_ch0_CAM1_stack0000_488nm_") for s in tab["mappedFilename"]
    )
    # steps of 32 voxels in x and 32 slices in z
    assert sorted(set(np.round(tab["StageX_um_"] / XY_PIXEL_SIZE, 6))) == [0, 32]  # noqa: TID251
    assert sorted(set(np.round(tab["StageZ_um_"] / 0.5, 6))) == [0, 32]  # noqa: TID251


def test_tileNamesWithoutStandardChannel(tmp_path):
    # names that match neither LLSM pattern (no '_ch[0-9]_') still give a list (MATLAB once
    # failed with 'expression' undefined); the channel is in mappedFilename
    names = [f"Tile_c0_{x:03d}x_000y_000z_0000t.tif" for x in range(2)]
    tiles(names, (40, 40, 32))(tmp_path)
    ig.generate_image_list(
        str(tmp_path),
        "tile_position",
        xy_pixel_size=XY_PIXEL_SIZE,
        dz=0.5,
        channel_patterns=("Tile",),
        tile_patterns=("0000t", "c0", "000x", "000y", "000z"),
        overlap_size=(8, 8, 0),
        overlap_size_type="pixel",
        io_scan=True,
    )
    tab = read_list(tmp_path / "ImageList_from_tile_positions.csv")
    np.testing.assert_allclose(np.sort(tab["StageX_um_"]), [0, 32 * XY_PIXEL_SIZE], atol=1e-9)
    assert all("_ch0_" in s for s in tab["mappedFilename"])


# ------------------------------------------------------------------------------------------
# more unit tests
# ------------------------------------------------------------------------------------------


def test_publish_leaves_no_temporary_file(tmp_path):
    ig.generate_image_list_from_tile_list(
        str(tmp_path) + "/", TILE_LIST_FNS, TILE_LIST_IND, (1, 1, 1), uuid="u1"
    )
    assert sorted(p.name for p in tmp_path.iterdir()) == ["ImageList_from_tile_list.csv"]


def test_unknown_method_raises(tmp_path):
    with pytest.raises(il.ImageListError):
        ig.generate_image_list(str(tmp_path), "other")


def test_sind_exact_and_reduced():
    assert ig._sind(30.0) == 0.5
    assert ig._sind(-210.0) == 0.5
    assert ig._sind(32.45) == pytest.approx(np.sin(np.deg2rad(32.45)), rel=1e-15)


def test_estimate_step_size(tmp_path):
    write_slice_lists(tmp_path, THREE, nslice=12)
    dz = ig.estimate_actual_step_size_from_encoder(str(tmp_path), dz=0.2)
    assert dz == pytest.approx(0.2, abs=1e-6)
    assert (tmp_path / "encoder_info" / "coordinates_info.mat").is_file()
    assert any(p.name.startswith("dz_actual_") for p in (tmp_path / "encoder_info").iterdir())
