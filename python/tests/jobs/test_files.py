"""batch_file_exist and wait_for_files. Ported from tests/framework/BatchFileExistTest.m (16)
and WaitForFilesTest.m (7); the test names keep the MATLAB names in snake_case.

batch_file_exist's result must equal MATLAB's exist(..., 'file') > 0. The expected flags of
the edge cases are MATLAB's (recorded with the oracle tree on Linux); the one difference by
design is a relative name that MATLAB finds on its path.
"""

from __future__ import annotations

import os
import random
import sys
from pathlib import Path

import pytest

from petakit5d.jobs import files
from petakit5d.jobs.files import batch_file_exist, exist, fileparts, wait_for_files

linux_only = pytest.mark.skipif(
    not sys.platform.startswith("linux"), reason="The folder checks are used on Linux only."
)


def make_files(folder: Path, inds, pattern: str = "f_{:03d}.tif") -> list[str]:
    fns = [str(folder / pattern.format(i)) for i in inds]
    for fn in fns:
        Path(fn).touch()
    return fns


def queries_in(folder: Path, inds, pattern: str = "f_{:03d}.tif") -> list[str]:
    return [str(folder / pattern.format(i)) for i in inds]


def exist_flags(queries) -> list[bool]:
    return [exist(q) for q in queries]


# the edge cases of BatchFileExistTest: names, then queries with MATLAB's exist() result
EDGE_NAMES = [
    "a.tif",
    "Case.tif",
    " a.tif",
    "#a.tif",
    "-a.tif",
    "a b.tif",
    "a%b.tif",
    "a'b.tif",
    "a[1].tif",
    "a(1).tif",
    "a,b.tif",
    "a&b.tif",
    ".hidden",
    "a~1.tif",
    "café.tif",
    "a*b.tif",
    "a?b.tif",
    "a\\b.tif",
    'a"b.tif',
    "n\nl.tif",
    "trail .tif ",
]
EDGE_QUERIES = [  # (path relative to the folder, MATLAB's exist(...,'file') > 0)
    ("A.TIF", False),
    ("CASE.TIF", False),
    ("case.tif", False),
    ("b.zarr", True),
    ("b.zarr/", True),
    ("b.zarr//", True),
    ("a.tif/", True),
    ("b.zarr/0.0.0", True),
    ("b.zarr/.", True),
    ("b.zarr/..", True),
    ("sub//x.tif", True),
    ("sub/../a.tif", True),
    ("sub/./x.tif", True),
    ("sub", True),
    ("missing.tif", False),
    ("missing/x.tif", False),
    ("missing/deeper/x.tif", False),
    ("a.tif/x", False),
    ("a.tif", True),
    ("missing.tif", False),
    ("st*r/s.tif", True),
    ("st*r/t.tif", False),
    ("lnk_a", True),
    ("lnk_b", True),
    ("lnk_dangle", False),
    ("lnkdir/x.tif", True),
    ("lnkdir", True),
]
RELATIVE_QUERIES = [
    ("a.tif", True),
    ("sub/x.tif", True),
    ("nothing_here.tif", False),
    ("a.tif/", True),
    ("./a.tif", True),
    ("sub/../a.tif", True),
]


def make_edge_cases(folder: Path):
    """The edge-case tree in folder; returns (queries, expected, folders). Relative queries
    are relative to the current folder (MATLAB's test changes to the temp folder)."""
    for name in EDGE_NAMES:
        (folder / name).touch()
    (folder / "b.zarr").mkdir()
    (folder / "b.zarr" / "0.0.0").touch()
    (folder / "sub").mkdir()
    (folder / "sub" / "x.tif").touch()
    (folder / "st*r").mkdir()
    (folder / "st*r" / "s.tif").touch()
    for link, target in (
        ("lnk_a", "a.tif"),
        ("lnk_b", "b.zarr"),
        ("lnk_dangle", "nothing_here"),
        ("lnkdir", "sub"),
    ):
        os.symlink(target, folder / link)
    base = str(folder) + "/"
    queries = [base + n for n in EDGE_NAMES] + [base + q for q, _ in EDGE_QUERIES]
    expected = [True] * len(EDGE_NAMES) + [e for _, e in EDGE_QUERIES]
    rel = os.path.relpath(folder)
    queries += [os.path.join(rel, q) for q, _ in RELATIVE_QUERIES] + ["/", "~", "~/"]
    expected += [e for _, e in RELATIVE_QUERIES] + [True, True, True]
    # a relative name MATLAB finds on its path (Python has no such path)
    queries.append("batch_file_exist.m")
    expected.append(False)
    return queries, expected, [folder, folder / "sub", folder / "st*r"]


def mismatches(queries, expected, actual):
    return [q for q, e, a in zip(queries, expected, actual, strict=True) if e != a]


# --- BatchFileExistTest ---------------------------------------------------------------------


def test_empty_and_single_inputs(tmp_path):
    assert batch_file_exist([]) == []
    assert batch_file_exist([str(tmp_path)]) == [True], "A folder counts as existing."
    assert batch_file_exist([str(tmp_path / "missing.tif")]) == [False]


def test_checks_files_and_folders_individually(tmp_path):
    fns = make_files(tmp_path, range(1, 10, 2))
    (tmp_path / "a.zarr").mkdir()
    queries = queries_in(tmp_path, range(1, 10)) + [str(tmp_path / "a.zarr")]
    assert batch_file_exist(queries) == [q in fns for q in queries[:9]] + [True]


def test_lists_the_folder_for_many_files(tmp_path):
    fns = set(make_files(tmp_path, range(1, 151, 2)))
    queries = queries_in(tmp_path, range(1, 151))
    assert batch_file_exist(queries) == [q in fns for q in queries]


def test_missing_folder_with_many_files(tmp_path):
    queries = queries_in(tmp_path / "missing", range(1, 102))
    assert batch_file_exist(queries) == [False] * 101


def test_background_pool_for_many_files_in_many_folders(tmp_path):
    queries, expected = [], []
    for d in range(1, 4):
        folder = tmp_path / f"d{d}"
        folder.mkdir()
        fns = set(make_files(folder, range(1, 41, 3)))
        q = queries_in(folder, range(1, 41))
        queries += q
        expected += [x in fns for x in q]
    assert batch_file_exist(queries, None, True) == expected


def test_saves_result_to_output_file(tmp_path):
    matfile = pytest.importorskip("petakit5d.io.matfile", reason="MAT-files come with F07")
    fns = make_files(tmp_path, range(1, 4))
    queries = fns + [str(tmp_path / "missing.tif")]
    out_fn = str(tmp_path / "exist.mat")
    flags = batch_file_exist(queries, out_fn)
    assert flags == [True, True, True, False]
    a = matfile.read_mat(out_fn)
    assert [bool(x) for x in list(a["file_exist_mat"])] == flags
    assert list(a["fileFullpaths"]) == queries


@linux_only
def test_matches_exist_on_edge_cases(tmp_path):
    # few paths per folder (each folder is opened, then the paths checked by name), all at
    # once and one by one
    queries, expected, _ = make_edge_cases(tmp_path)
    for use_java in (True, False):
        flags = batch_file_exist(queries, None, True, use_java)
        assert mismatches(queries, expected, flags) == []
        flags = [batch_file_exist([q], None, True, use_java)[0] for q in queries]
        assert mismatches(queries, expected, flags) == []
    assert mismatches(queries, expected, exist_flags(queries)) == []


@linux_only
def test_matches_exist_on_edge_cases_in_listed_folders(tmp_path):
    # more than 100 paths in each folder: the folders are listed
    queries, expected, folders = make_edge_cases(tmp_path)
    for folder in folders:
        made = set(make_files(folder, range(1, 151, 2)))
        q = queries_in(folder, range(1, 151))
        queries += q
        expected += [x in made for x in q]
    for use_java in (True, False):
        flags = batch_file_exist(queries, None, True, use_java)
        assert mismatches(queries, expected, flags) == []


@linux_only
def test_folders_and_trailing_separators(tmp_path):
    # zarr outputs are folders; a trailing separator doesn't matter
    (tmp_path / "a.zarr").mkdir()
    (tmp_path / "b.tif").touch()
    p = str(tmp_path) + "/"
    queries = [
        p + "a.zarr",
        p + "a.zarr/",
        p + "a.zarr//",
        p + "b.tif",
        p + "b.tif/",
        p + "c.zarr/",
        p,
        p[:-1],
    ]
    expected = [True, True, True, True, True, False, True, True]
    assert exist_flags(queries) == expected
    for use_java in (True, False):
        assert batch_file_exist(queries, None, True, use_java) == expected
        assert batch_file_exist(queries * 20, None, True, use_java) == expected * 20


def test_missing_parent_folders(tmp_path):
    fns = make_files(tmp_path, range(1, 10, 2))
    (tmp_path / "file.tif").touch()
    for n in (1, 3, 150):
        missing = (
            queries_in(tmp_path / "missing", range(1, n + 1))
            + queries_in(tmp_path / "missing" / "deeper", range(1, n + 1))
            + queries_in(tmp_path / "file.tif", range(1, n + 1))
        )
        queries = fns + missing + [str(tmp_path / "f_002.tif")]
        expected = [True] * len(fns) + [False] * (len(missing) + 1)
        for use_java in (True, False):
            assert batch_file_exist(queries, None, True, use_java) == expected, n


def test_duplicate_paths(tmp_path):
    fns = set(make_files(tmp_path, range(1, 301, 2)))
    for n in (5, 300):
        inds = list(range(1, n + 1)) + list(range(n, 0, -1)) + [1, 2, 3]
        queries = queries_in(tmp_path, inds)
        for use_java in (True, False):
            flags = batch_file_exist(queries, None, True, use_java)
            assert flags == [q in fns for q in queries]


def test_mixed_found_and_missing_in_many_folders(tmp_path):
    # folders with 1 to 150 paths (opened or listed), some missing
    rng = random.Random(7)
    queries = []
    for d in range(1, 13):
        folder = tmp_path / f"d{d:02d}"
        if d % 4 != 0:
            folder.mkdir()
            make_files(folder, rng.sample(range(1, 151), rng.randint(1, 150)))
        queries += queries_in(folder, rng.sample(range(1, 151), rng.randint(1, 150)))
    rng.shuffle(queries)
    expected = [os.path.exists(q) for q in queries]
    assert any(expected) and not all(expected)
    for use_java in (True, False):
        assert batch_file_exist(queries, None, True, use_java) == expected


def test_finds_files_made_between_calls(tmp_path):
    for n in (1, 150):
        folder = tmp_path / f"n{n}"
        folder.mkdir()
        queries = queries_in(folder, range(1, n + 1)) + [str(folder / "out.zarr")]
        assert batch_file_exist(queries, None, True) == [False] * (n + 1)
        make_files(folder, [n])
        (folder / "out.zarr").mkdir()
        assert batch_file_exist(queries, None, True) == [False] * (n - 1) + [True, True]


def test_large_folder(tmp_path):
    # 20,000 files, 40,000 paths
    fns = set(make_files(tmp_path, range(1, 40001, 2), "f_{:05d}.tif"))
    queries = queries_in(tmp_path, range(40000, 0, -1), "f_{:05d}.tif")
    assert batch_file_exist(queries, None, True) == [q in fns for q in queries]


@pytest.mark.skipif(os.name == "nt", reason="Folder read permissions are POSIX only.")
def test_unreadable_folder(tmp_path):
    # a folder that can be searched but not listed: exist() finds its files, so
    # batch_file_exist must too
    if hasattr(os, "geteuid") and os.geteuid() == 0:
        pytest.skip("root reads any folder")
    folder = tmp_path / "noread"
    folder.mkdir()
    fns = set(make_files(folder, range(1, 151, 2)))
    folder.chmod(0o311)
    try:
        for n in (3, 150):
            queries = queries_in(folder, range(1, n + 1))
            expected = [q in fns for q in queries]
            assert exist_flags(queries) == expected, "exist() itself changed."
            assert batch_file_exist(queries, None, True) == expected
    finally:
        folder.chmod(0o755)


def test_existing_output_file_returns_cached_result(tmp_path):
    # KNOWN ISSUE of MATLAB, kept (batch_file_exist.m:43-45): with an existing result file the
    # call returns before it assigns its output, so a second call fails
    fns = make_files(tmp_path, range(1, 4))
    out_fn = tmp_path / "exist.mat"
    out_fn.touch()
    with pytest.raises(RuntimeError, match="not assigned"):
        batch_file_exist(fns + [str(tmp_path / "missing.tif")], str(out_fn))


def test_cells_of_paths_take_their_first_path(tmp_path):
    fns = make_files(tmp_path, range(1, 3))
    assert batch_file_exist([[fns[0], "x"], [str(tmp_path / "no"), fns[1]]]) == [True, False]
    with pytest.raises(TypeError):
        batch_file_exist([[fns[0], fns[1]]])


def test_fileparts_as_matlab():
    assert fileparts("/a/b/c.d.tif") == ("/a/b", "c.d", ".tif")
    assert fileparts("/a.tif") == ("/", "a", ".tif")
    assert fileparts("/a/.hidden") == ("/a", "", ".hidden")
    assert fileparts("c.zarr") == ("", "c", ".zarr")
    assert fileparts("/a/b/") == ("/a/b", "", "")
    assert fileparts("/a//b") == ("/a/", "b", "")
    assert fileparts("noext") == ("", "noext", "")


# --- WaitForFilesTest -----------------------------------------------------------------------


class PauseSpy:
    """Records each pause (without waiting) and runs actions at given calls (1-based)."""

    def __init__(self, actions=None):
        self.waits = []
        self.actions = actions or {}

    def __call__(self, t):
        self.waits.append(t)
        action = self.actions.get(len(self.waits))
        if action:
            action()


@pytest.fixture
def pauses(monkeypatch):
    def install(actions=None):
        spy = PauseSpy(actions)
        monkeypatch.setattr(files, "_pause", spy)
        return spy

    return install


def looks_at_least_a_second_apart(waits, n=None):
    # the first pause is a second, none but the last is shorter
    if n is not None:
        assert len(waits) == n
    assert waits[0] == 1
    assert all(w >= 1 - 1e-12 for w in waits[:-1])


def test_without_files_it_is_one_pause(pauses):
    spy = pauses()
    assert wait_for_files([], 4.5) is False
    assert wait_for_files([[], []], 3) is False
    assert spy.waits == [4.5, 3]


def test_files_that_do_not_come_take_the_whole_wait(tmp_path, pauses):
    spy = pauses()
    assert wait_for_files([[str(tmp_path / "a.tif")]], 3.5) is False
    looks_at_least_a_second_apart(spy.waits)
    assert abs(sum(spy.waits) - 3.5) < 1e-9


def test_wait_ends_once_the_files_of_a_group_exist(tmp_path, pauses):
    # groups {a, b} and {c}: b comes during the second second, a during the third
    fns = [str(tmp_path / n) for n in ("a.tif", "b.tif", "c.tif")]
    spy = pauses({2: lambda: Path(fns[1]).touch(), 3: lambda: Path(fns[0]).touch()})
    assert wait_for_files([fns[:2], fns[2:]], 30) is True
    looks_at_least_a_second_apart(spy.waits, 3)


def test_wait_ends_once_all_files_of_one_group_exist(tmp_path, pauses):
    fns = [str(tmp_path / n) for n in ("a.tif", "b.tif")]
    spy = pauses({2: lambda: Path(fns[0]).touch(), 3: lambda: Path(fns[1]).touch()})
    assert wait_for_files([fns], 30) is True
    looks_at_least_a_second_apart(spy.waits, 3)


def test_results_can_be_folders(tmp_path, pauses):
    result = tmp_path / "result.zarr"
    spy = pauses({2: result.mkdir})
    assert wait_for_files([[str(result)]], 30) is True
    looks_at_least_a_second_apart(spy.waits, 2)


def test_patterns_count_only_finished_names(tmp_path, pauses):
    # a pattern, as the pipeline looks for a stitched result: the result's temporary name
    # (with a suffix) doesn't count, its final name does
    result = tmp_path / "Scan_Iter_0000_msec_0000000001msecAbs.zarr"
    spy = pauses({2: lambda: Path(f"{result}_1234").mkdir(), 3: result.mkdir})
    assert wait_for_files([[str(tmp_path / "Scan_Iter_0000_msec*Abs.zarr")]], 30) is True
    looks_at_least_a_second_apart(spy.waits, 3)


def test_slow_looks_are_fewer(tmp_path, pauses, monkeypatch):
    # a look that takes 0.3 s (slow storage, many files): the next ones are 3 s or more
    # apart, so that looking takes at most a tenth of the wait, which stays the same
    now = [0.0]
    monkeypatch.setattr(files, "_clock", lambda: now[0])

    def slow_batch_file_exist(paths, *args):
        now[0] += 0.3
        return [False] * len(paths)

    monkeypatch.setattr(files, "batch_file_exist", slow_batch_file_exist)
    spy = pauses()
    assert wait_for_files([[str(tmp_path / "a.tif")]], 10) is False
    assert spy.waits[0] == 1
    assert all(w >= 3 - 1e-9 for w in spy.waits[1:-1])
    assert abs(sum(spy.waits) - 10) < 1e-9


def test_no_look_after_the_last_pause(tmp_path, pauses):
    # a file that comes during the last pause is for the caller's next round
    fn = tmp_path / "a.tif"
    spy = pauses({3: fn.touch})
    assert wait_for_files([[str(fn)]], 2.5) is False
    assert spy.waits == [1, 1, 0.5]


def test_patterns_need_the_whole_name(tmp_path, pauses):
    # (as MATLAB: start and end must not overlap)
    spy = pauses({1: lambda: (tmp_path / "aba").mkdir()})
    assert wait_for_files([[str(tmp_path / "ab*ba")]], 2) is False
    assert spy.waits == [1, 1]


def test_patterns_with_more_than_one_star(tmp_path, pauses):
    spy = pauses({1: lambda: (tmp_path / "ab_x_cd.zarr").mkdir()})
    assert wait_for_files([[str(tmp_path / "a*_x_*.zarr")]], 5) is True
    assert spy.waits == [1]
