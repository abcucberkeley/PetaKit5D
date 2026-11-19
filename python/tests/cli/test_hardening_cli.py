"""Tests the mutation campaign asked for (E06, docs/hardening.md): the command line."""

import io

from petakit5d.cli import main as cli


def test_list_with_unported_functions_has_their_section():
    text = cli.list_text(include_unported=True)
    head, _, tail = text.partition("\nNot ported:\n")
    assert tail, "the 'Not ported:' section is missing"
    assert "XR_crop_dataset" in head
    assert all(line.startswith("  ") and ": " in line for line in tail.strip("\n").split("\n"))
    assert "Not ported:" not in cli.list_text()


def test_an_unknown_option_is_a_usage_error():
    out, err = io.StringIO(), io.StringIO()
    assert cli.main(["--no-such-option"], stdout=out, stderr=err) == cli.EXIT_USAGE
    assert err.getvalue() == "pk5d: error: unknown option --no-such-option\n"
