import random
import re
import sys

import numpy as np

from petakit5d.core import ids


def test_uuid_form():
    u = ids.get_uuid()
    if sys.platform.startswith("win"):
        assert len(u) == 4  # MATLAB keeps 4 characters on Windows
    else:
        # XRDeskewRotateFrameTest and TiffZarrConversionTest rely on numel(get_uuid()) >= 32
        assert re.fullmatch(
            r"[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}", u
        )
    assert ids.get_uuid() != u


def test_uuid_leaves_global_random_state():
    # TaskStateLeakTest: get_uuid must not change the global random numbers
    np.random.seed(5)  # noqa: NPY002 (the global state must not change)
    random.seed(5)
    before = (np.random.get_state()[1].copy(), random.getstate())  # noqa: NPY002 (the global state must not change)
    ids.get_uuid()
    after = (np.random.get_state()[1], random.getstate())  # noqa: NPY002 (the global state must not change)
    np.testing.assert_array_equal(before[0], after[0])
    assert before[1] == after[1]


def test_hostname(capsys, monkeypatch):
    monkeypatch.setenv("HOSTNAME", "node12")
    if not sys.platform.startswith("win"):
        assert ids.get_hostname() == "node12"
        assert capsys.readouterr().out == "Hostname: node12 \n"
    assert ids.get_hostname(verbose=False)
