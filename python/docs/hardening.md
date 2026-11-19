# Hardening: mutation campaign, coverage, isolation and soak

Four checks find what the test suite misses, beyond the suite itself:

| Check | Where | What it shows |
|---|---|---|
| mutation campaign | `tools/mutation/` | single-line bugs that no test catches, per area |
| line coverage | `coverage` with per-test contexts | code no test runs; also the campaign's test selection |
| isolation | `tests/jobs/test_hardening_isolation.py` | every task function gives in one process, after any other tasks, what a fresh process gives |
| soak | `tests/cluster/test_hardening_soak.py` | the job framework on a real Slurm under crashes, memory kills, time limits and a second caller |

The tests the campaign asked for are in `tests/<area>/test_hardening_<area>.py`.

## Mutation campaign

A mutant is one changed span of one line of `src/petakit5d`, in a copy of `src/` (hard links of
the other modules), run against the tests that execute that line. `tools/mutation/mutants.py`
has the operators; `tools/mutation/campaign.py` the steps:

| Operator | Mutations |
|---|---|
| `cmp` | `<` `<=`, `>` `>=`, `==` `!=`, `is` `is not`, `in` `not in` |
| `arith` | `+` `-`, `*` `/`, `/` `*`, `//` `/`, `%` `//`, `**` `*`, `&` `\|`, `<<` `>>`, and the augmented forms |
| `num` | an int `n` to `n + 1`, a float `x` to `2 x` |
| `bool` | `True` `False`, `and` `or`, `not x` to `x` |
| `stmt` | a call, an assignment to an attribute or item, `return`, `raise`, `continue` or `break` alone on its line, replaced by `pass` |
| `cond` | the test of an `if`, `elif` or `while` negated |
| `swap` | two neighbouring items of a tuple or list of 2 to 4 numbers or names (axis orders) |
| `call` | `min` `max`, `floor` `ceil`, `any` `all`; a one-argument rounding or cast dropped (`mround(x)` to `(x)`) |

Strings, f-strings, docstrings, annotations, `assert`, `TYPE_CHECKING` blocks, `__all__`,
`__repr__`, the types of `isinstance` and `except`, `print` and logger calls and `flush()` are
never mutated. A mutant is kept only if the module still compiles and its syntax tree changed.

### 1. Coverage with per-test contexts

The selection comes from one run of the full suite with line coverage per test (it gives the
coverage numbers below too). Install `coverage>=7.10` and `pytest-cov>=7` next to the test
environment, and write an rc file:

```ini
[run]
source_pkgs = petakit5d
parallel = true
patch = subprocess
data_file = /local/scratch/cov/.coverage
```

`patch = subprocess` measures the child processes too (job entry points, process pools, the
isolation harness); the plugin `tools/mutation/subprocess_contexts.py` labels their data with
the test that started them. In the golden-store mode (`PETAKIT5D_GOLDEN_DIR`,
`PETAKIT5D_ORACLE_DIR` set):

```bash
PYTHONPATH=src:tools/mutation python -m pytest -p no:cacheprovider -p subprocess_contexts \
    --cov --cov-config=covrc --cov-context=test --cov-report= --junitxml=junit.xml tests
coverage combine --append --rcfile=covrc /local/scratch/cov
```

The suite can run in shards of test folders (one JUnit file each); `coverage combine` merges
their data files.

### 2. Generate, plan, run, report

```bash
python tools/mutation/campaign.py generate --root . --out mutants.jsonl --seed S \
    [--exclude 'src/petakit5d/jobs/exec.py' ...]          # ceil(lines / 100) per area, 12 to 50
python tools/mutation/campaign.py plan --root . --mutants mutants.jsonl --coverage .coverage \
    --junit junit_*.xml --out plan.jsonl --area-budget 600 --suite-budget 600 \
    --header-tests tests/params/test_check.py
python tools/mutation/campaign.py run --root . --plan plan.jsonl --out results_3.jsonl \
    --work /local/scratch/mut --shard 3/8 --parallel 4 --threads 2
python tools/mutation/campaign.py report --plan plan.jsonl --results results_*.jsonl \
    --rechecks recheck.jsonl --dispositions dispositions.json --out report.md
```

- **generate**: lines drawn uniformly per area, seeded (dense code isn't favoured). A pass over
  some files only (files that were changing during the first pass, say) takes `--only` globs
  and `--prefix p2/`, which keeps its ids apart from the first pass's.
- **plan**: per mutant, the area phase gets the tests of `tests/<area>/` that execute its line
  (fastest first, up to `--area-budget` seconds of the coverage run's durations), the suite
  phase the other tests that do. A line of a function's signature or decorators runs at import
  time and takes effect in calls: it gets the tests that run the function's body. A line at
  module or class level (constants, tables, dataclass fields) gets the whole `tests/<area>/`
  and the tests that run any line of its module. Both also get the tests of `--header-tests`
  files, which read signatures without calling them (the registry check compares every
  registered function's defaults with MATLAB's). Tests that failed or were skipped in the
  coverage run are left out. A line that no test runs is "not covered" and is not run.
- **run**: one pytest process per phase with `-x`, from the tree given by `--root` with
  `PYTHONPATH` at the mutant's `src/`, in a process group the tool ends at the phase's time
  limit (4 x the tests' durations + 2 minutes). Run it in batch jobs, not on a workstation:
  each shard is one job (`--shard i/n`), its `--work` folder on the job's local disk. A test
  that kills a mutant is run once on the unchanged tree too: if it fails there the status is
  `base_fails`. Collection errors and crashes count as killed; time-outs too. A rerun of a
  shard skips the ids already in its results file.
- **recheck**: `run --recheck map.json` (`{mutant id: [node ids]}`) runs each survivor against
  the tests written for it, from a tree that has them (the sources must be the same: a mutant
  holds the SHA-256 of its module and refuses another version).
- **report**: per area the mutants, those killed by the area's tests and by other tests, the
  suite's score, survivors, not covered, killed by new tests, those with a written reason, the
  ones still open and the score with the new tests; then every mutant the suite did not kill
  with its disposition (`{id: {"kind": "equivalent" | "unreachable" | "diagnostic" | "gpu" |
  "optional" | "owner", "why": ...}}`).

`campaign.py show` prints mutants, `campaign.py tree` writes one mutant's `src/` for debugging.

### The C++ extension modules

`tools/mutation/cpp.py` does the same for `cpp/<group>/` (not the vendored libraries, the
GPL add-on or the shared headers of `cpp/common`): sites by tokens (comments, strings and
preprocessor lines left alone; binary operators only with blanks around them, so templates
are not touched), the tests of a group from the coverage of the Python modules that load it,
and per job one copy of `python/` and one CMake build on local disk, then per mutant an
incremental build and pytest. A mutant that doesn't compile is reported as `uncompiled` and
left out of the scores; the results go through `campaign.py report` like the Python ones.

```bash
python tools/mutation/cpp.py generate --root . --out cpp_mutants.jsonl --seed S
python tools/mutation/cpp.py plan --root . --mutants cpp_mutants.jsonl --coverage .coverage \
    --junit junit_*.xml --out cpp_plan.jsonl --budget 900
python tools/mutation/cpp.py run --root . --plan cpp_plan.jsonl --out cpp_results_1.jsonl \
    --work /local/scratch/cppmut --shard 1/4 --jobs 8 --cmake-arg=-DCMAKE_CXX_COMPILER=g++
```

### What survives, and why

- **Equivalent mutants** are the usual survivors: a default that every caller passes, a
  `frozen=True` nothing relies on, `reshape(-2)` (NumPy reads any negative size as `-1`), a
  division by a sign of +-1, a fast path whose fallback computes the same, a cast of what is
  already a Python int, a swap of two names that the expression treats alike (`min(s, b)`).
- **GPU-only code** runs only in the `gpu`-marked tests, which a CPU campaign does not run; a GPU
  round needs its own job with per-test contexts on a GPU node.
- **Optional extras** (`figures`, `n5`, `formats`) are not in the test environments.
- **Platform-dependent bits**: NaN payloads equal MATLAB's on x86 where NumPy's operation is
  used directly; the tests that kill such mutants do so on ARM (macOS CI). In C++, branches
  for another compiler (MSVC) or CPU (without AVX2).
- **Memory-only effects** (C++): a write one element past a buffer that changes no result; only
  a sanitizer build (AddressSanitizer) sees it.

## Coverage

Line coverage of the full suite in the golden-store mode, test `069f3cb` (2026-11-18): 89.1% of
49,226 statements (at `ec71880`, with the hardening tests: 89.4% of 49,417).

| Area | Statements | Covered | Area | Statements | Covered |
|---|---|---|---|---|---|
| cli | 330 | 92.4% | kernels | 2995 | 84.8% |
| core | 1932 | 93.4% | matlab | 33 | 87.9% |
| decon | 4254 | 84.4% | matlabcompat | 4641 | 94.5% |
| deskew_rotate | 1654 | 96.3% | package | 62 | 93.5% |
| fft | 154 | 78.6% | params | 2801 | 89.9% |
| flatfield | 542 | 95.2% | pipeline | 1250 | 93.4% |
| gpu | 826 | 43.9% | stitch | 6258 | 95.5% |
| io | 3151 | 88.9% | testing | 5323 | 81.9% |
| jobs | 6970 | 92.2% | tools | 6050 | 87.6% |

Most of what is not covered is code the CPU suite cannot run (GPU paths, the exact FFT
add-on's loader, BaSiCPy and the figures extra) or developer command lines (the oracle's `run`,
`status` and `clean`, the benchmark runner).

## Isolation of every task function

`tests/jobs/test_hardening_isolation.py` (marked `slow`; Linux, the runners' platform: macOS
can't load cpp-tiff and cpp-zarr together yet; about a minute) builds a small synthetic dataset
per scenario and runs the port's wrapper until the framework call whose tasks it captures, as
MATLAB's `tests/isolation/isolation_prepare.m` does: `run_tasks` is replaced by a recorder that
lets the earlier calls run and stops the wrapper at the captured call; the frames that the decon
and stitching wrappers run themselves in local mode are recorded at `run_record`. The scenario
folders are the template that `jobs.isolation.run_harness` copies into a new work folder
before each run: every task in a process of its own, then all tasks of all scenarios in one
process in 5 random orders with the runners' restore between tasks; each run's outputs must
equal the fresh processes', bit for bit (`testing.compare`'s trees).

`test_every_task_function_has_a_scenario` compares the task functions the scenarios run with
the package's `Task.of(...)` and `Task("module:function")` call sites. A new task function needs
a scenario in `SCENARIOS` (a builder and the call to capture), or an entry in `NOT_RUN` with its
reason (an optional dependency, figures, cluster-only).

## Soak under failures

`tests/cluster/test_hardening_soak.py` (markers `cluster` and `slow`) runs with a cluster config
(`PETAKIT5D_TEST_SLURM_CONFIG`, giving jobs `--mem-per-cpu=1G` and at least 2 CPUs),
`PETAKIT5D_TEST_SLURM_SCRATCH` on storage the nodes share and
`PETAKIT5D_TEST_SLURM_HARDENING=<tasks>x<jobs>[x<max jobs>]`, once with the task runners and once
without. It injects crashes (SIGSEGV, SIGABRT), a task that needs 1.3 GiB (Slurm's memory limit
kills its 1-CPU job; the second try has 2 CPUs), and tasks that sleep past their job's time
limit; a second process calls `run_tasks` with the same tasks 20 s after the first. Both calls
must finish every task, each output must be published once, and no job or call folder may be
left. J11's `test_soak.py` covers raised errors, exits, `scancel` and outside cancels.
