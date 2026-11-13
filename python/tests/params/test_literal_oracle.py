"""The literal parser against MATLAB R2024b (stored references, no MATLAB needed).

Every text of the corpus (mcc parser defaults, demo literals, edge cases) is evaluated with
``eval`` and ``str2num`` in MATLAB (literal_corpus.py, oracle_standin/); here the Python
parser must give the same class, size and bits, or fail where MATLAB fails. The formatter's
text of each value must evaluate in MATLAB to the same value (the round trip).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from petakit5d.params.literal import LiteralError, format_literal, parse, str2num
from petakit5d.params.values import describe

pytestmark = pytest.mark.oracle

MANIFEST = Path(__file__).parent / "data" / "literal_oracle.json"

# Valid MATLAB that the parser deliberately does not take (it raises LiteralError instead)
UNSUPPORTED = {
    "1i": "complex numbers",
    "0x10 + 1": "arithmetic on integer classes",
    "0x1F:0x22": "colon ranges of integer classes",
}


def _items() -> list[dict]:
    return json.loads(MANIFEST.read_text())["items"]


def _ids(items: list[dict]) -> list[str]:
    return [f"{i}:{it['text'][:40]!r}" for i, it in enumerate(items)]


ITEMS = _items()
EVAL_ITEMS = [it for it in ITEMS if it["kind"] != "str2num"]


def _matlab_ok(d: dict) -> bool:
    return "error" not in d


@pytest.mark.parametrize("it", EVAL_ITEMS, ids=_ids(EVAL_ITEMS))
def test_eval_agrees_with_matlab(it):
    m = it["eval"]
    if it["text"] in UNSUPPORTED:
        with pytest.raises(LiteralError):
            parse(it["text"])
        return
    if not _matlab_ok(m):
        with pytest.raises(LiteralError):
            parse(it["text"])
        return
    assert describe(parse(it["text"])) == m


@pytest.mark.parametrize("it", ITEMS, ids=_ids(ITEMS))
def test_str2num_agrees_with_matlab(it):
    if it["text"] in UNSUPPORTED:
        pytest.skip(UNSUPPORTED[it["text"]])
    x, ok, _ = str2num(it["text"])
    assert {"value": describe(x), "ok": ok} == it["str2num"]


@pytest.mark.parametrize("it", EVAL_ITEMS, ids=_ids(EVAL_ITEMS))
def test_round_trip_through_matlab(it):
    if not _matlab_ok(it["eval"]) or it["text"] in UNSUPPORTED:
        pytest.skip("not a literal")
    # the stored formatted text is what MATLAB evaluated; the formatter must still give it
    assert format_literal(parse(it["text"])) == it["formatted"]
    assert it["formatted_eval"] == it["eval"]
    assert describe(parse(it["formatted"])) == it["eval"]


def test_corpus_covers_parsers_and_demos():
    """Every default of the mcc parsers is in the corpus (or is site-specific text, kept out
    of the repository); the demos' literals are in it."""
    data = json.loads(MANIFEST.read_text())
    c = data["counts"]
    assert c["parser_files"] >= 55
    assert c["parser_defaults"] > 1000
    assert c["parser_defaults_not_literal"] == []
    assert c["parser_defaults_kept"] + c["parser_defaults_scan_dropped"] == c["parser_defaults"]
    demos = {it["source"].split(":")[0] for it in data["items"] if it["source"].startswith("demo")}
    assert c["demo_files"] >= 20 and len(demos) >= 10, sorted(demos)
