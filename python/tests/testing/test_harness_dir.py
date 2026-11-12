"""The oracle harness finds the repo's MATLAB side also when petakit5d is installed (CI)."""

from pathlib import Path

from petakit5d.testing import oracle


def test_harness_dir_is_next_to_the_tests():
    want = Path(__file__).resolve().parents[1] / "oracle_matlab"
    assert oracle.harness_dir().resolve() == want
    assert (oracle.harness_dir() / "ORACLE_COMMIT").is_file()
