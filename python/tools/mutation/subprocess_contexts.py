"""pytest plugin for the campaign's coverage run: the coverage data of the subprocesses a test
starts (coverage's ``[run] patch = subprocess``) is labelled ``<node id>|sub``, so that lines run
only in child processes (job entry points, process pools, the isolation harness's children) map
to the tests that start them.

    PYTHONPATH=src:tools/mutation python -m pytest -p subprocess_contexts --cov \\
        --cov-config=covrc --cov-context=test --cov-report= --junitxml=junit.xml

(docs/hardening.md has the rc file). This module imports only the standard library and pytest.
"""

import base64
import json
import os

import pytest

_VAR = "COVERAGE_PROCESS_CONFIG"


def _with_context(raw: str, context: str) -> str:
    data = json.loads(base64.b64decode(raw.encode()).decode())
    data["context"] = context
    return base64.b64encode(json.dumps(data).encode()).decode()


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_setup(item):
    raw = os.environ.get(_VAR)
    if raw:
        item._pk5d_cov_raw = raw
        os.environ[_VAR] = _with_context(raw, item.nodeid + "|sub")


@pytest.hookimpl(trylast=True)
def pytest_runtest_teardown(item, nextitem):
    raw = getattr(item, "_pk5d_cov_raw", None)
    if raw is not None:
        os.environ[_VAR] = raw
