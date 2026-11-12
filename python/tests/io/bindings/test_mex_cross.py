"""The cross matrix against the MATLAB mex (cross_driver.py + matlab/io1_cross_mex.m).

Runs where a finished matrix folder is given (PETAKIT5D_IO1_CROSS_DIR, holding
results.json); stand-in for the oracle harness until it lands. Every check must
pass or be one of the known gaps below, each listed in the IO1 report.
"""

import json
import os
import re

import pytest

pytestmark = pytest.mark.oracle

# (check, case pattern) -> why it does not pass
KNOWN = [
    (
        "bytes_writetiff_eq_mex",
        r"^t_.*_none$",
        "writetiff 'none' goes through saveastiff (classic TIFF), not the mex",
    ),
    (
        "py_reads_writetiff",
        r"^t_u?int64_.*_none$",
        "saveastiff writes 64-bit integers as doubles (MATLAB quirk)",
    ),
    (
        "chunk_set_py_eq_mex",
        r"^z_(sparse_false|zeros|regions_zero)$",
        "Python skips all-zero chunks; writezarr writes them",
    ),
    (
        "chunk_set_py_eq_mex",
        r"^z_shard$",
        "Python shards only as Zarr v3; createzarr writes a v2 sharded layout",
    ),
    (
        "zarray_py_eq_mex",
        r"^z_shard$",
        "Python shards only as Zarr v3; createzarr writes a v2 sharded layout",
    ),
    ("zarr_python_reads_mex", r"^z_shard$", "zarr-python cannot read cpp-zarr's v2 sharded layout"),
    ("mex_size_of_py", r"^z_shard$", "getImageSize reads only .zarray (the Python file is v3)"),
    (
        "zarray_py2d_eq_mex",
        r"^z_2d$",
        "the package keeps 2D arrays 2D; PetaKit5D stores [ny, nx, 1]",
    ),
    (
        "py_regions_eq_expected",
        r"^z_convert$",
        "region writes convert with astype, not MATLAB's rounding",
    ),
    (
        "py_regions_on_mex_created",
        r"^z_convert$",
        "region writes convert with astype, not MATLAB's rounding",
    ),
    (
        "mex_reads_py_regions",
        r"^z_convert$",
        "region writes convert with astype, not MATLAB's rounding",
    ),
    (
        "chunk_bytes_py1_eq_mex1",
        r"^z_convert$",
        "region writes convert with astype, not MATLAB's rounding",
    ),
]


def _records():
    folder = os.environ.get("PETAKIT5D_IO1_CROSS_DIR")
    if not folder:
        return []
    out = []
    for name in sorted(os.listdir(folder)):
        p = os.path.join(folder, name, "results.json")
        if os.path.isfile(p):
            out += [dict(r, run=name) for r in json.load(open(p))["records"]]
    return out


RECORDS = _records()


def _known(r):
    return next(
        (why for chk, pat, why in KNOWN if chk == r["check"] and re.search(pat, r["case"])), None
    )


@pytest.mark.skipif(not RECORDS, reason="PETAKIT5D_IO1_CROSS_DIR is not set or holds no results")
def test_every_check_passes_or_is_a_known_gap():
    bad = [r for r in RECORDS if r["status"] in ("fail", "error") and not _known(r)]
    assert not bad, "\n".join(
        f"{r['run']} {r['side']} {r['case']} {r['check']}: {r['status']} {r['detail']}"
        for r in bad[:40]
    )


@pytest.mark.skipif(not RECORDS, reason="PETAKIT5D_IO1_CROSS_DIR is not set or holds no results")
def test_tiff_files_and_zarr_chunks_are_byte_identical_to_the_mex():
    for chk in (
        "bytes_py_eq_mex",
        "bytes_py1_eq_mex1",
        "bytes_pydirect_eq_py",
        "chunk_bytes_py_eq_mex",
        "zarray_py_eq_mex",
    ):
        rs = [r for r in RECORDS if r["check"] == chk]
        assert rs, chk
        assert all(r["status"] in ("pass", "n/a") or _known(r) for r in rs), chk


@pytest.mark.skipif(not RECORDS, reason="PETAKIT5D_IO1_CROSS_DIR is not set or holds no results")
def test_known_gaps_still_show():
    # a gap that no longer shows was fixed upstream: update KNOWN and the report
    for chk, pat, why in KNOWN:
        hits = [r for r in RECORDS if r["check"] == chk and re.search(pat, r["case"])]
        if hits:
            assert any(r["status"] in ("fail", "error") for r in hits), (
                f"{chk} {pat} passes now: {why}"
            )
