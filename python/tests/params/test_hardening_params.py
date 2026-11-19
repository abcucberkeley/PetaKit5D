"""Tests the mutation campaign asked for (E06, docs/hardening.md): the parameter registry, the
MATLAB-name layer and parameter literals."""

import subprocess
import sys

import pytest

import petakit5d.matlab
from petakit5d.params import Param
from petakit5d.params.check import check_signature_defaults
from petakit5d.params.kinds import KINDS
from petakit5d.params.matlab_source import _split_top_eq
from petakit5d.params.registry import Registry
from petakit5d.params.spec import _strip_parens


@pytest.mark.timeout(60)
def test_a_declaration_is_split_at_its_assignment_only():
    assert _split_top_eq("tol (1,1) double = 3") == ("tol (1,1) double", "3")
    assert _split_top_eq("a == b") == ("a == b", None)
    assert _split_top_eq("x ~= 1") == ("x ~= 1", None)
    assert _split_top_eq("y >= 0 = 2") == ("y >= 0", "2")


def test_two_flags_keep_their_order():
    k = KINDS["flags_yxz"]
    assert k.from_matlab([1, 0]) == (1, 0) and k.to_matlab([0, 1]) == (0, 1)
    assert k.from_matlab([0, 0, 1]) == (1, 0, 0)


def test_only_a_matching_pair_of_outer_parentheses_is_stripped():
    assert _strip_parens("(a) + (b)") == "(a) + (b)"
    assert _strip_parens("((a + b))") == "a + b"
    assert _strip_parens("(f(x))") == "f(x)"


def test_a_list_default_is_reported_as_mutable():
    reg = Registry()

    @reg.register("XR_hardening_fake", Param("sizes", "@isnumeric", "Sizes."))
    def fake(*, sizes=[1, 2]):  # noqa: B006 (the mutable default under test)
        return sizes

    probs = check_signature_defaults(reg.find("XR_hardening_fake"))
    assert any("mutable default [1, 2]" in p for p in probs)


def test_the_matlab_layer_lists_the_registered_names():
    names = dir(petakit5d.matlab)
    assert "XR_crop_dataset" in names and "XR_deskew_rotate_data_wrapper" in names


def test_a_dunder_lookup_on_the_matlab_layer_does_not_load_the_registry():
    code = (
        "import sys, petakit5d.matlab as m; print(hasattr(m, '__wrapped__')); "
        "print('petakit5d.params.registry' in sys.modules)"
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr
    assert r.stdout.split() == ["False", "False"]
