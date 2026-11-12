"""The job framework: tasks run locally or as Slurm jobs (MATLAB's generic computing framework).

Modules: config (cluster configuration files), slurm (the Slurm commands), shapes (CPUs, time
limits, array ids of jobs), framework (run_tasks). Names are imported on first use, so importing
the package costs nothing.
"""

from __future__ import annotations

import importlib
from typing import Any

_EXPORTS = {
    "ClusterConfig": "config",
    "ConfigCache": "config",
    "cluster_config_value": "config",
    "load_config": "config",
    "Slurm": "slurm",
    "SlurmCheck": "slurm",
    "check_slurm_cluster": "slurm",
    "job_cpus": "shapes",
    "single_job_cpus": "shapes",
    "retry_cpus": "shapes",
    "time_str": "shapes",
    "run_tasks": "framework",
    "Task": "tasks",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str) -> Any:
    module = _EXPORTS.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(importlib.import_module(f"{__name__}.{module}"), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(_EXPORTS))
