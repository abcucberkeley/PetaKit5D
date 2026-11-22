"""The probe helper of the port's checks of third-party code, on a fake MATLAB function.

``matlab_parse/pkpy_parse_probe.m`` calls a function with a battery of arguments and records,
for each call, 'accepted' (no error while the function's inputParser parsed it: the call
returned, or failed before or after the parse), the identifier of an error of inputParser's
own, or 'raised(<identifier>)' for an error that a check raised itself, which inputParser
passes on with the check's identifier ('raised()' for one without an identifier). inputParser
raises all of them at the function's call of parse, which a call it refuses (the calibration)
shows the helper. ``matlab_fake/pkpy_fake_checks.m`` (its own folder, so that the probe cases'
keys do not depend on it) has checks of each kind and errors before and after its parse; the
cases record what MATLAB R2024b makes of them through the helper. The helper in use before
2026-11-21 recorded an error of a check's own as accepted, and an error of a called function's
inputParser after the parse as a refusal.
"""

from __future__ import annotations

from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, mcell

HELPERS = ["tests/params/matlab_parse", "tests/params/matlab_fake"]
REFUSED = "MATLAB:InputParser:ArgumentFailedValidation"

# probe -> (the call's arguments, the helper's record)
CALLS = {
    "accepted": ([1.0], "accepted"),
    "refused": (["text"], REFUSED),
    "failed_before_the_parse": ([None], "accepted"),  # x(1) of an empty x
    "check_false": ([1.0, "Own", True], REFUSED),
    "unmatched": ([1.0, "Foo", 1.0], "MATLAB:InputParser:UnmatchedParameter"),
    "missing_value": ([1.0, "Own"], "MATLAB:InputParser:ParamMissingValue"),
    "raised": ([1.0, "Own", "c"], "raised(pkpyfake:ownCheck)"),
    "raised_without_identifier": ([1.0, "Own", mcell(1.0)], "raised()"),
    "raised_further_down": ([1.0, "Own", {"a": 1.0}], "raised(pkpyfake:deeper)"),
    "raised_through_anonymous": ([1.0, "Anon", "c"], "raised(pkpyfake:ownCheck)"),
    "raised_by_mathworks": ([1.0, "Positive", -1.0], "raised(MATLAB:expectedPositive)"),
    "failed_after_the_parse": ([1.0, "After", True], "accepted"),
    "nested_parser_after_the_parse": ([1.0, "Nested", True], "accepted"),
    # its own parse refuses a call of itself: the same place, one frame further down
    "own_parser_one_call_down": ([1.0, "Recursive", True], "accepted"),
}
NAMES = [f"p_{k}" for k in CALLS]


def _probe_args(calib: list, calls: dict) -> list:
    return [
        "pkpy_fake_checks",
        mcell(),
        mcell(*calib),
        mcell(*[f"p_{k}" for k in calls]),
        mcell(*[mcell(*a) for a, _ in calls.values()]),
    ]


def expected(*_args) -> dict[str, str]:
    """the helper's records of the calls, as stated in CALLS"""
    return {f"p_{k}": want for k, (_, want) in CALLS.items()}


CASES = [
    Case(
        "fake.outcomes",
        matlab="pkpy_parse_probe",
        args=_probe_args([1.0, "Own"], CALLS),
        python=expected,
        helpers=HELPERS,
        files=False,
    ),
    # a calibration call that inputParser takes: no parse site, the helper stops
    Case(
        "fake.no_parse_site",
        matlab="pkpy_parse_probe",
        args=_probe_args([1.0], {"accepted": CALLS["accepted"]}),
        python=None,
        expect_error="pkpy:probe:noParseSite",
        helpers=HELPERS,
        files=False,
    ),
    # one that fails after the parse, with an error of the function's own: none either
    Case(
        "fake.calibration_after_the_parse",
        matlab="pkpy_parse_probe",
        args=_probe_args([1.0, "After", True], {"accepted": CALLS["accepted"]}),
        python=None,
        expect_error="pkpy:probe:noParseSite",
        helpers=HELPERS,
        files=False,
    ),
]


@oracle.parametrize(CASES)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


def test_the_cases_state_every_outcome():
    """accepted (returned, failed before and after the parse), inputParser's own refusals and
    errors of a check's own, with and without identifiers"""
    wants = {w for _, w in CALLS.values()}
    assert "accepted" in wants and "raised()" in wants
    assert any(w.startswith("MATLAB:InputParser:") for w in wants)
    assert any(w.startswith("raised(") and w != "raised()" for w in wants)
    assert sorted(expected()) == sorted(NAMES)
