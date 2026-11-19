"""Puts ``examples/phase1`` first on the import path: the Phase 2 examples use its helpers
(``_common``: data paths, step times) and run its examples, whose results they build on.

Modules of this folder that have a namesake there (``demo_useful_tools``,
``run_all_demos_for_the_paper``) are loaded by path under their own names (``load``), so
``import demo_useful_tools`` stays phase 1's.
"""

from __future__ import annotations

import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE1 = os.path.normpath(os.path.join(HERE, os.pardir, "phase1"))
if not sys.path or sys.path[0] != PHASE1:
    if PHASE1 in sys.path:
        sys.path.remove(PHASE1)
    sys.path.insert(0, PHASE1)


def load(name: str):
    """this folder's example ``name`` as the module ``phase2_<name>``"""
    key = "phase2_" + name
    if key in sys.modules:
        return sys.modules[key]
    spec = importlib.util.spec_from_file_location(key, os.path.join(HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[key] = mod
    try:
        spec.loader.exec_module(mod)
    except BaseException:
        del sys.modules[key]
        raise
    return mod
