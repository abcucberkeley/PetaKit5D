"""Tests the mutation campaign asked for (E06, docs/hardening.md): the pipeline."""

import numpy as np

from petakit5d.pipeline import automatic
from petakit5d.pipeline.filenames import matches_channel_pattern


def test_a_pattern_that_is_no_regular_expression_is_matched_as_text():
    # (regexpPattern refuses the unmatched '(': '.' is then no wildcard either)
    texts = ["x_Cam(A.b_ch0", "x_Cam(AYb_ch0"]
    assert matches_channel_pattern(texts, "Cam(A.b") == [True, False]
    assert matches_channel_pattern(texts, "Cam.A") == [True, True]


def test_string_of_a_logical_and_of_empty():
    assert automatic._string([]) == ""
    assert automatic._string(True) == "true"
    assert automatic._string(np.array([False])) == "false"
