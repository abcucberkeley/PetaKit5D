# End-to-end benchmarks against MATLAB

How the port's workflows are timed against MATLAB PetaKit5D as a whole: the demos and a real
dataset, in local and cluster mode, on CPU and GPU nodes, with the time of every stage apart
(compute, process start-up and teardown, framework waiting, queue wait), peak memory and CPU
hours. The per-function harness is in `benchmarks.md`; this page adds the files in
`benchmarks/e2e/` and the method for whole workflows.

| File | What |
|---|---|
| `benchmarks/e2e/task_costs.json` | one framework task as a process, three implementations on one node (harness spec) |
| `benchmarks/e2e/task_probe.py` | the Python side of it: a tile's deskew/rotation or a tiny crop, marked; or a call folder for `python -m petakit5d.jobs.exec` |
| `benchmarks/e2e/e2e_cluster_workflow.m` | a workflow in MATLAB, local or cluster mode (MATLAB jobs, or compiled jobs with `mccMode`) |
| `benchmarks/e2e/e2e_cluster_workflow.py` | the same workflow in the port, with the same steps and parameters |

## One task as a process

`task_costs.json` runs each task the way a job runs it, interleaved on one node:

| Variant | Process |
|---|---|
| `ml_dsr`, `ml_tiny` | `matlab -batch` (the framework's flags), `setup([])`, the function string |
| `mcc_dsr`, `mcc_tiny` | the compiled master (`run_mccMaster.sh <Runtime> <function> <arguments>`), when `PK5D_MCC_MASTER` and `PK5D_MCR` name a build and a Runtime |
| `py_dsr`, `py_tiny` | `task_probe.py`: the task's imports, `mark("ready")`, the task |
| `pyexec_dsr`, `pyexec_tiny` | `python -u -m petakit5d.jobs.exec <call folder> 1`, the job entry point on a call folder of the same task |

`dsr` is the deskew/rotation of one demo tile (500 frames, tiff in and out), `tiny` a crop of a
32³ box (nothing but the process). The harness gives start-up (launch to the first line of the
work), the task and teardown (the last line to the end of the process) for the marked
variants, and the total for the others; a compiled or entry-point task's overhead is its total
minus the marked task's time on the same node. The gate is start-up plus teardown at most 1 s
per CPU task (`process_costs_gpu.json` has the GPU side: at most 3 s).

## Workflows in local and cluster mode

`e2e_cluster_workflow.{m,py} WORKFLOW DATA CONFIG ...`:

| Workflow | Steps |
|---|---|
| `generic` | demo_generic_computing_framework: crop and deskew/rotate each of the 16 tiles, one task each, `memAllocate` 20 |
| `dsr` | demo_geometric_transformation's dataset step: the combined deskew/rotation of the 16 tiles, zarr results |
| `dsr_large` | the large-file deskew/rotation of zarr volumes (two channels), with the parameters of a real pipeline |

An empty config runs the workflow in this process (local mode); a config file runs it as the
framework's jobs. The three implementations get the same config apart from the launcher:
MATLAB jobs use `MatlabLaunchStr`, compiled jobs `MCCMasterStr` and `MCRParam` (with
`GNUparallel` and `wholeNodeJob` false, so each task is one job as in the others), Python jobs
`PythonLaunchStr` (the interpreter with its environment). `taskMetrics: true` in the config
gives every task a record (`jobs.metrics`; MATLAB writes the same format), from which the stage
times come:

| Stage | From |
|---|---|
| compute | the tasks' `wallSec` |
| process start-up | `start - processStart` of each task (the interpreter or MATLAB started to the task's first line) |
| teardown | the job's end (the job log's `End time`, or Slurm's) minus the task's `end` |
| queue wait | Slurm's `Start - Submit` per job; for the workflow, the time with jobs pending and none running |
| framework waiting | the rest of the workflow's wall time: polling, the grace periods, the master's own steps |
| peak memory, CPU hours | the tasks' `rssPeakGB` and `cpuSec`, and Slurm's `MaxRSS`, `TotalCPU` and allocated CPUs x elapsed |

## Method

- **Interleave on one node.** Local-mode comparisons run all implementations in one job on one
  node, in the order A B B A (or A B C, C B A), each on a fresh copy of the data; a cluster-mode
  comparison runs its masters from one job and pins the worker jobs to one node type (named in
  the job name), with the implementations rotated between rounds. Queue wait depends on the
  cluster's load at that moment and is reported apart, never in a ratio.
- **Fresh copies.** Each run gets its own copy of the inputs (hard links are enough when the
  inputs are only read) and its own result folders, so no run skips work another did.
- **Same inputs, each side's layout.** Tiff inputs are shared; zarr inputs written by MATLAB are
  read by both, as existing data would be.
- **Steps the port gets from MATLAB.** Where a step is not ported (until 2026-11-18 the stitching
  registration, whose results the demo examples were given in the runs of that date), the
  Python run is given MATLAB's results and skips that computation: report those steps with
  MATLAB's time for the skipped part beside them.
- **Memory state.** Record `/proc/buddyinfo` (free 2 MB blocks per NUMA node) at the start of
  each run. A fresh process on a node whose free memory is fragmented can stall for seconds in
  the kernel's compaction when NumPy asks for huge pages (`madvise`, NumPy's default on Linux):
  run Python with NumPy's default and with `NUMPY_MADVISE_HUGEPAGE=0` as two arms.
- **Shared storage.** A file read right after it was written can come from the server rather
  than the client's cache, at a third of the speed; when a demo reads a file it has just
  written, compare that step only with a side that did the same.
