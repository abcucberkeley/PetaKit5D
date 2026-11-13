# Benchmarks against MATLAB: the benchmark harness

The port must be no slower than MATLAB PetaKit5D, and it should lose MATLAB's start-up and
teardown. The harness (`petakit5d.testing.bench`) measures both the same way for every task:
MATLAB and Python processes interleaved on one node, the time a process spends before and
after its work measured apart from the work, everything about the run recorded, and the task's
speed gate checked.

| Part | Where | What |
|---|---|---|
| harness | `petakit5d.testing.bench` | runs a spec, marks and timings, summaries, comparisons, gates, history |
| MATLAB side | `benchmarks/common/pk5d_bench_*.m` | marks, timings and version facts from a MATLAB process |
| process costs | `benchmarks/common/process_costs*.json`, `startup_probe.py` | what a task process costs outside its work, Python and MATLAB |
| self-check | `benchmarks/common/self_check.json` | the same build twice: the harness's own noise |
| cluster helper | `brun.sh` (outside the repo, with the cluster's own names) | one comparison as one batch job on a pinned node type |

## Method

- **One job per comparison**, all variants on the same node: nodes of one type differ by up
  to 30%, so only processes that shared a node are compared.
- **Interleaved repeats**: every variant runs once to warm up, then in at least 5 rounds,
  forward and backward in turn (A B C, C B A, ...), so drift and the node's state hit every
  variant alike. Medians are reported with their range; operations timed inside a process
  pool all calls of all rounds and report the range of the per-round medians.
- **Node type pinned for timings only**, with the type in the job name (`..._c24_...`), and
  the job asks for the cores it measures with (24 for a gate at 1, 4 and 24 threads).
- **Recorded**: node name, CPU model, CPUs of the process, memory, load before and after, the
  busiest processes on the node (own or other), free memory per NUMA node, free 2 MB blocks
  per NUMA node (`/proc/buddyinfo`, order 9 and up) and the transparent huge page settings
  (`enabled`, `defrag`) at the start and the end, the NUMA nodes the processes may use, the
  Slurm job, GPUs (name, memory in use, load, persistence mode) before and after,
  Python and package versions (`pip freeze`), MATLAB's version and thread count, the source
  commit and a hash of uncommitted changes, every process's marks, exit status, peak memory
  (max RSS of the process and the processes it waited for), CPU time and log.

## Inputs: each side reads its own layout

The port's arrays are `(z, y, x)` C-ordered and its zarr files native (shape `[nz, ny, nx]`,
`"order": "C"`); MATLAB's are `(y, x, z)` Fortran-ordered and its zarr files `[ny, nx, nz]`,
`"order": "F"`. A comparison gives each side the same data in its own layout:

- **Tiff** inputs are shared: the file is the same for both.
- **Zarr** inputs are written twice from one `(z, y, x)` volume by
  `petakit5d.testing.synth.write_paired_inputs(folder, name, vol, 'zarr', block_size)` (the
  `paired_input` fixture in tests): `<folder>/matlab/<name>.zarr` in MATLAB's layout for the
  MATLAB variant and `<folder>/python/<name>.zarr` native for the Python variant, same name,
  same data, chunks `block_size` in `(z, y, x)` (MATLAB's file has them in `(y, x, z)`).
- **Reading MATLAB-layout zarr** in the port (one transposing pass) is a separate case: it is
  benchmarked and reported apart (its own op names, e.g. `zarr_read_matlab_layout`), never as
  the Python side of a MATLAB comparison.
- In-memory inputs of a Python benchmark are `(z, y, x)` C arrays (`np.zeros(shape)`); a
  MATLAB variant makes or reads its own.

The self-check (`self_check.py`) times an FFT pair and a C-order copy of a `(z, y, x)` volume
of the same size (it timed a copy into `(y, x, z)` order before the convention change, so its
history starts again at `copy_c_order`).

## Process costs: start-up and teardown apart from the work

Each process gets `PK5D_BENCH_FILE`, a JSON-lines file for marks and timings. The harness
writes `launch` just before it starts the process and `exit` when the process has been reaped.

| Mark | Python | MATLAB |
|---|---|---|
| `first` | the harness's bootstrap, before any import | the first statement of `-batch` (or `-r`) |
| `ready` | the script calls `mark("ready")` after its imports | after `setup([])` and the `addpath` of the variant's folders |
| `last` | when the script, module or code returns (or raises) | after the variant's code |

| Phase | From → to | Meaning |
|---|---|---|
| `start` | launch → first | the interpreter or MATLAB itself |
| `setup` | first → ready | imports; `setup([])` |
| `startup` | launch → ready | process start to the first line of the work |
| `task` | ready → last | the work |
| `teardown` | last → exit | from the last line to the end of the process |
| `overhead` | startup + teardown | what the process costs outside its work |
| `total` | launch → exit | |

Marks of a process's own (between `ready` and `last`, or instead of `ready`) give steps
(`seg:a->b`), as the GPU probe does for import, device count, context, kernel and FFT plan.
The first run of each variant (the warm-up) and the first MATLAB process of the job are
reported apart: MATLAB's first exit in a fresh job has been much slower than later ones.

**File cache.** A variant with `"cache": "cold"` has the page-cache pages of its files dropped
before each run (`posix_fadvise(DONTNEED)` on every file of the interpreter, the environment,
the source tree and the extension folder, or of its `evict` list). The kernel's file metadata
stays cached, and on network storage so does the server's cache: cold means the client's page
cache. A warm variant that runs after an eviction gets an untimed run first, so cold variants
don't make the warm ones after them cold.

**Environment on local disk.** `"env_copy": "local"` runs the variant with a copy of the
interpreter, its virtual environment, the source tree and the extension folder on the job's
temporary folder (node-local disk), made once per job.

## A spec

JSON (or a Python file that defines `SPEC`), in `benchmarks/<area>/`:

```json
{
 "name": "tiff_read",
 "repeats": 5,
 "setup": ["cp {data}/some_tile.tif {work}/in.tif"],
 "variants": [
  {"id": "ml", "impl": "matlab", "paths": ["benchmarks/io"],
   "code": "bench_read('{work}/in.tif', [1, 4, 24], 3)"},
  {"id": "py", "impl": "python", "script": "benchmarks/io/bench_read.py",
   "args": ["{work}/in.tif", "--threads", "1", "4", "24"]}
 ],
 "compare": [{"a": "py", "b": "ml"}],
 "gate": {"ratio_max": 1.05, "threads": [1, 4, 24]}
}
```

| Field | Meaning |
|---|---|
| `name` | the benchmark's name (results and history are grouped by it) |
| `repeats`, `warmup` | rounds (at least 5) and warm-up runs per variant (default 1) |
| `order` | `alternate` (default: forward, backward) or `rotate` |
| `setup` | shell commands run once in the work folder before the runs |
| `params` | more placeholders |
| `variants` | the processes (below) |
| `compare` | pairs `{a, b, op_map}`: candidate a against reference b; default: every Python variant against the only MATLAB variant |
| `gate` | overrides the task's gate from `tasks.json` (below) |
| `output_bytes` | the largest output's size, for the warning about free 2 MB blocks (default 1 GB) |

| Variant field | Meaning |
|---|---|
| `id`, `impl` | unique id; `python`, `matlab` or `command` |
| `script` / `module` / `code`, `args` | Python: a file (relative to `python/`), a module or code, with arguments |
| `python`, `python_flags` | another interpreter; flags before `-c` |
| `code`, `paths`, `setup`, `flags`, `launch`, `matlab` | MATLAB: statements, folders to add (relative to `python/`), `setup([])` first (default true), launcher flags (default the framework's `-nodisplay -nosplash -nodesktop -nojvm`), `batch` (default) or `r` (`-r` with `exit`, as the framework launched before), the MATLAB command (default `PETAKIT5D_MATLAB` or `matlab`) |
| `argv` | a command (`command`): only launch and exit are known |
| `env`, `cwd` | more environment variables; the folder to start in (MATLAB: the oracle tree; others: the work folder) |
| `cache`, `evict`, `env_copy` | `cold`, the files to drop from the cache, `local` (above) |
| `timeout` | seconds; a process past it is stopped with its process group (`rc` != 0, `timed_out`) |
| `when` | run the variant only when this environment variable is set (e.g. a compiled build) |

Placeholders: `{src}` (the synced `python/`), `{work}` (the job's work folder on local disk),
`{python}`, `{oracle}` (`PETAKIT5D_ORACLE_DIR`), `{data}` (`PETAKIT5D_DEMO_DATA`), `{common}`,
`{env[NAME]}`, and the spec's `params`. Other braces are left alone.

In a Python benchmark script (the import is cheap; nothing happens without `PK5D_BENCH_FILE`):

```python
from petakit5d.testing.bench import mark, record, timeit

import numpy as np

mark("ready")
x = np.zeros((500, 512, 512), np.float32)  # (z, y, x)
for t in (1, 4, 24):
    record("sum", timeit(lambda t=t: x.sum(), repeats=3), threads=t)
```

In MATLAB: `pk5d_bench_record('sum', t, times)` and `pk5d_bench_mark('label')`. The two sides
must use the same op names and thread counts (or an `op_map`).

**The free of an output is not timed** unless it belongs to the measured pipeline: `timeit`
keeps the result until the timer has stopped (`include_free=True` times it), and MATLAB
scripts `clear` the output after `toc`. On a node whose free memory has run out of 2 MB blocks,
large arrays get 4 KB pages, and freeing a 0.92 GB one costs 37-58 ms instead of about 3 ms.
The report shows the free 2 MB blocks per NUMA node and warns when a NUMA node the processes
may use had fewer than an output needs (the spec's `output_bytes`, default 1 GB) at the start
or the end; compare such a run only with runs on the same node.

## Gates

A task's gate comes from its `benchmark_gate` text in `tasks.json`; the harness checks the parts
it can read and shows the text for the rest:

| In the text | Checked |
|---|---|
| `(1.05x)`, or "no slower" (1.05, as PLAN.md defines it) | ratio of the medians, candidate / reference, of every compared operation |
| "at 1, 4 and 24 threads" | each compared operation was measured at all of them |
| "start-up + teardown at most 1 s" | the candidate's `overhead` median |
| "memory ... 1.1 x" | ratio of peak memory |
| `none ...` | nothing |

A spec's `gate` can set `ratio_max`, `threads`, `overhead_max`, `mem_ratio_max`, `metrics`
(glob patterns of metrics to check: `op:*`, `startup`, `task`, ...) and `within` (|ratio - 1|
at most this: a build against itself). An operation's metric carries its thread count,
`op:<name>@<threads>`: name an operation at every thread count as `op:<name>@*`. In a gate on
ratios (`ratio_max`, `within`, `threads`, `mem_ratio_max`), each pattern in `metrics` must
match a row of every compared pair, and the gate needs a compared
pair (and, without `metrics`, an operation timed on both sides); otherwise the gate fails with
an `unmatched:` check that lists the rows there are. With `metrics`, `threads` applies to the
operations they select. The verdict is met, missed or nothing to check; a met
gate whose text asks for more (data sizes, modes) still needs a reader. A missed gate is
reported, never averaged away. The median of the per-round ratios (`paired`) and whether the
round-median ranges overlap are shown next to each ratio.

## Command line

```
python -m petakit5d.testing.bench run SPEC [--out DIR] [--repeats N] [--keep] [--tasks tasks.json]
python -m petakit5d.testing.bench summary RESULT.json
python -m petakit5d.testing.bench gate RESULT.json [--task ID] [--text "gate text"]
python -m petakit5d.testing.bench history BENCH_DIR [--write] [--threshold 0.05]
python -m petakit5d.testing.bench node
```

`history` reads every result under a folder and writes `HISTORY.md` and `history.json`: per
benchmark and metric, each result's date, task, commit (`+` when it had local changes), node
type and node, the candidate's median, the reference's, their ratio, and the change from the
previous result of the same metric on the same node type (the ratio's change when there is a
reference, which cancels the node's speed; else the median's). Changes above 5% are flagged
when the median moved by more than 5 ms (`--threshold`, `--min-seconds`); results of
unpinned jobs (node type `any`) form their own line, and so do runs with interleaved memory
(`ilv`).

**Memory placement.** On two-socket nodes, a memory-bound kernel at 24 threads is timed with
both processes under `numactl --interleave=all` (the cluster helper's `-I` runs the driver that
way; the processes it starts inherit the policy). Without it, single-threaded first touch puts
a process's arrays on one socket. The result records the policy the driver ran under
(`node.numa_policy`).

## Environment variables

| Variable | Meaning |
|---|---|
| `PK5D_BENCH_FILE` | set by the harness for each process: its marks and timings |
| `PK5D_BENCH_ROUND` | set by the harness: the round (negative for warm-up runs) |
| `PK5D_SRC` | the `python/` tree specs are relative to (default: the current folder) |
| `PK5D_TASKS_JSON` | `tasks.json`, for the task's gate |
| `PK5D_BENCH_COMMIT`, `PK5D_BENCH_DIRTY`, `PK5D_BENCH_TASK`, `PK5D_BENCH_NODETYPE`, `PK5D_BENCH_WORKTREE` | recorded with the result (the cluster helper sets them) |
| `PETAKIT5D_MATLAB`, `PETAKIT5D_ORACLE_DIR`, `PETAKIT5D_DEMO_DATA` | as for the oracle harness |
