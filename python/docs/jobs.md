# The job framework: run_tasks, task records, the job entry point, local, pool and cluster modes

`petakit5d.jobs` runs a list of tasks until their outputs exist: in this process (local
mode), on a local process pool, or as Slurm jobs. It is the port of MATLAB's
`generic_computing_frameworks_wrapper` and keeps its behavior; what changes is that a task is
a record, not MATLAB code text. Area wrappers build lists of inputs, outputs and tasks, call
`run_tasks`, and never call Slurm themselves.

## run_tasks

```python
from petakit5d.jobs.framework import run_tasks
from petakit5d.jobs.tasks import Task

tasks = [Task("petakit5d.tools.mip:mip_file", [inp, out], {"axis": [0, 0, 1]}) for inp, out in ...]
done = run_tasks(input_paths, output_paths, tasks, mem_allocate=20, config_file=config_file)
```

`run_tasks(input_paths, output_paths, tasks, **params) -> list[bool]`, one done flag per task.

- `input_paths[i]`: a path, or a list of paths. A task whose input is missing is not run
  and counts as done (MATLAB's behavior).
- `output_paths[i]`: the path whose existence means task i is done (a file or a folder;
  trailing separators are removed). Tasks publish it last, by one rename.
- `tasks[i]`: a `Task`, its record (a dict), a callable without arguments or a
  `functools.partial`; with `language='bash'` a shell command string.

The parameters are those of MATLAB's `generic_computing_frameworks_wrapper` in snake_case
with MATLAB's defaults (`framework.RUN_TASKS_PARAMS` lists the MATLAB name, Python name,
default and meaning of each):

| MATLAB | Python | default |
|---|---|---|
| finalOutFullpath | final_out_fullpath | `""` |
| clusterType | cluster_type | `"slurm"` (or `"parfor"`) |
| parseCluster | parse_cluster | `True` |
| masterCompute | master_compute | `True` |
| jobLogDir | job_log_dir | `"../job_logs"` (empty: the default) |
| tmpDir | tmp_dir | `""` |
| maxCPUNum, minCPUNum, cpusPerTask | max_cpu_num, min_cpu_num, cpus_per_task | 24, 1, 1 |
| memPerCPU, memAllocate | mem_per_cpu, mem_allocate | 20.9, `None` |
| wholeNodeJob | whole_node_job | `False` |
| uuid | uuid | `""` (a new one) |
| GPUJob, GPUNum | gpu_job, gpu_num | `False`, 0 |
| maxTrialNum | max_trial_num | 3 |
| unitWaitTime | unit_wait_time | 3 (minutes) |
| maxJobNum | max_job_num | `inf` |
| minTaskJobNum, taskBatchNum, minBatchNum, paraBatchNum | min_task_job_num, task_batch_num, min_batch_num, para_batch_num | 1, 1, 1, 1 |
| runExtraTasks | run_extra_tasks | `False` |
| PythonLaunchStr (new; was MatlabLaunchStr) | python_launch_str | `""` (this interpreter) |
| BashLaunchStr | bash_launch_str | `""` |
| SlurmParam | slurm_param | `""` (MATLAB's default named its site: it comes from the config) |
| SlurmConstraint | slurm_constraint | `""` |
| jobTimeLimit | job_time_limit | 24 (hours; `None`: no limit) |
| queryInterval | query_interval | 3 (seconds) |
| language | language | `"python"` (`"matlab"` is taken for it) or `"bash"` |
| GNUparallel, paraJobNum, masterParaFactor | gnu_parallel, para_job_num, master_para_factor | `False`, 1, 1 |
| ConfigFile | config_file | `""` |
| localWorkerNum | local_worker_num | 1 |
| taskMetrics | task_metrics | `False` |
| taskRunner | task_runner | `True` |
| SlurmAlternatives, gpuHeadroom | slurm_alternatives, gpu_headroom | `None`, `None` |
| nice | nice | 0 |

`mcc_mode`, `matlab_launch_str`, `mcr_cache_root`, `mcr_param` and `mcc_master_str` (and the
same keys in config files) are accepted and ignored. Types are checked as MATLAB's
inputParser checks them (`parse_cluster=1` is an error).

A config file (JSON; `.mat` through `petakit5d.io.matfile`) overrides the parameters: its keys
match the MATLAB names without regard to case (`MemPerCPU` sets `mem_per_cpu`), and `[]`
(MATLAB's empty) is `None` where a parameter takes `None`.

Order of decisions, as in MATLAB: config file; `jobLogDir` empty is the default; `nice > 0`
adds `--nice=<round(nice)>` to `slurm_param`; in cluster mode a master in a Slurm job runs
tasks only if they fit its memory; one task with `master_compute` runs here; `cluster_type
'parfor'` goes to the process pool; with `task_metrics`, the run's metrics folder and its
run_start record; cluster mode goes to the cluster loop, which runs the tasks locally (with
MATLAB's warning) when sinfo fails; without cluster and with `local_worker_num > 1` the
process pool runs first and the local loop finishes what it left; with `task_metrics`, the
run_end record after the backend returns (not after an exception, as in MATLAB).

## Local mode, as MATLAB runs it

`framework.run_local(call)`, the `parseCluster = false` path of
`slurm_cluster_generic_computing_wrapper`:

1. If `final_out_fullpath` exists, all tasks are done. If all outputs exist, all are done.
2. Inputs are looked at once (for tasks without output); a list of inputs is looked at in
   every round, each path.
3. Rounds until every task is done or tried `max_trial_num` times. Each round looks at the
   outputs again (from the second round), then goes through the batches in order
   (`task_batch_num` tasks each):
   - tasks with a missing input count as done; a task whose output exists is done and its
     marker removed;
   - the batch's marker is `<tmp_dir or output folder>/<output name without extension>.tmp`
     of its last task. If it exists and is younger than `unit_wait_time` minutes (whole
     seconds of its modification time), the batch is left to the other caller; else the
     marker is renewed and the batch runs. Without the output folder there is no marker;
   - the batch runs in this process: its tasks in order; if one raises, the whole batch
     runs once more from its first task, and a second exception propagates (markers stay).
     A bash batch is one shell command (`cmd1;cmd2`), its exit status ignored;
   - each task of the batch counts a try; the marker is removed if the output exists.
   - After the round: progress (`Time ... s: k / n (...%) are finished!`); if the final
     output exists, all tasks are done.
4. A last look at the outputs after the loop (a batch run in its last try counts).

Two callers on the same data share the work through the markers: each leaves the batches
whose markers are fresh to the other and waits for their outputs. The markers are not
atomic (as in MATLAB): two callers that start the same task within milliseconds both run it.

The one change in timing: a round that ran nothing (it waits for another caller's marker)
sleeps `LOCAL_IDLE_WAIT` (0.2 s); MATLAB loops without a pause.

## Local process pool

Without cluster and with `local_worker_num > 1`, `jobs.local.run_pool` runs the tasks on a pool
of worker processes before the sequential loop, as MATLAB's matlab_pool_generic_computing_wrapper
does; the loop then retries what the pool left (a task that failed there runs once more here).
Workers: `min(local_worker_num, floor(memory / 2 / (mem_allocate + 2)))` (the memory cap only
when mem_allocate is given; one or fewer: no pool), at most one per task to run. Each worker is
a fresh interpreter with this process's environment, current folder and sys.path; it gets
`floor(CPUs / workers)` threads (thread variables of OpenMP and BLAS; on Linux its own block of
CPUs) and runs one task at a time, in task order. A task's output is printed when it ends.
`jobs.local.worker_info()` tells a task which worker it runs on. `cluster_type='parfor'` runs the
tasks on the pool only (`max_trial_num` rounds); a GPU job gets one worker per usable GPU, each
with its own CUDA_VISIBLE_DEVICES. The pool runs Python and bash tasks.

## Cluster mode

`jobs.cluster.run_cluster(call)` is the parseCluster branch of
slurm_cluster_generic_computing_wrapper, a state machine with MATLAB's decisions and timings
(named constants in `jobs.cluster`). It checks Slurm once (sinfo, J01's `check_slurm_cluster`;
without it the tasks run here through `run_local`), shapes the jobs (`job_shape`: job_cpus,
whole-node slots, minBatchNum, paraBatchNum, minTaskJobNum, runExtraTasks, GPUNum), writes the
call folder in `tmp_dir` or the job log folder, and then goes round by round, in MATLAB's order:

- a query when `query_interval` has passed (one batched `squeue -j`); a failed or timed-out
  query is "unknown", never "gone", so nothing is submitted again on it;
- the outputs looked at with `batch_file_exist`;
- the batches in order: missing inputs done; `max_job_num` (a batch counts while its job is
  listed and not all its tasks are done); with `master_compute` the master takes the last
  waiting batch (scancel of its pending job; that try is not counted); running jobs left
  alone; a grace of 45 s (30 s with master compute) after a job left the queue; then one
  `sbatch --array=<b mod 5000> -o <logs>/job_%A_%a.out -e <logs>/job_%A_%a.err
  --cpus-per-task=N <SlurmParam> <SlurmConstraint> -t H:MM:00 --comment=<call>:<b> --wrap=...`
  per batch. CPU jobs get CPUs x (try + 1) on a retry, GPU jobs keep theirs;
- the master runs its batch with `run_batch_in_process(..., mode="master")` (or in slots in
  slot mode), after `min(max(30, t1 / 2), 180)` s when no job pends;
- a pause of `query_interval` minus the round's time (doubled when all jobs pend and the
  master doesn't compute); progress; the final output.

After the loop: a last look at the outputs, `scancel` of all the call's jobs (also after an
exception and on SIGTERM, which MATLAB leaves to orphans), and the call folder removed when all
tasks are done.

The wrap is MATLAB's job log layout: `printf "Start time: "; date ...; echo Python command:
'<cmd>'; <cmd>; echo; printf "End time: "; date ...`, with `<cmd>` =
`<PythonLaunchStr> -u -m petakit5d.jobs.exec <call folder> <b>`. No task text is in the
command line.

Slot mode (MATLAB's compiled mode: `whole_node_job`, or `gnu_parallel` with `para_job_num > 1`)
runs `python -m petakit5d.jobs.slots <call folder> <b> --slots P --tries T [--gpu]` in place of
GNU parallel: one process per task, up to P at once, starts 0.2 s apart, a failing task run
twice, tasks whose output (or the final output) exists skipped, `PARALLEL_JOBSLOT` and
`PARALLEL_SEQ` set, and in GPU jobs a GPU per slot (`CUDA_VISIBLE_DEVICES`). Tries 1 and 2 of a
CPU job keep its CPUs with P and then round(P / 2) slots, try 3 runs one task at a time; GPU
jobs keep P slots and get CPUs x (try + 1). The master runs its batch with P - 1 slots
(`master_para_factor`).

## Jobs on more than one place (SlurmAlternatives)

A cluster config may list other places for its jobs after its main settings, in order of
preference (`jobs.alternatives`, the port of MATLAB's slurm_alternatives; the module docstring
has the rules):

```json
{
  "SlurmParam": "-p cpu --qos normal -n1 --mem-per-cpu=21418M",
  "memPerCPU": 20.9,
  "SlurmAlternatives": [
    {"SlurmParam": "-p gpu --qos normal -n1 --mem-per-cpu=7800M", "memPerCPU": 7.6,
     "maxCPUNum": 32, "cpuJobs": true}
  ]
}
```

- An entry has any of SlurmParam, SlurmConstraint, memPerCPU, maxCPUNum, minCPUNum,
  gpuHeadroom and cpuJobs (any case); keys it leaves out are the main settings', except
  memPerCPU (its `--mem-per-cpu`, else the smallest memory per CPU of its nodes) and cpuJobs
  (false). Without the key, or with `[]` or `{}`, nothing changes: no Slurm command is added.
  Jobs of whole nodes use the main place only.
- A job goes to the first place where Slurm can start it now: the nodes and the user's pending
  jobs (read at most every 20 s; jobs placed since are taken from that view), then `sbatch
  --test-only`. It is shaped for that place: at a place with as much memory per CPU as the job's
  shape, the CPUs a config of only that place gives it; at one with less, its CPUs and memory.
  With no room anywhere it waits at the first place that isn't of GPU nodes only.
- A pending job that waited 2 minutes in its place moves to a place that can start it now, by
  one `scontrol update` (partition, QOS, account, CPUs, memory, GRES, features, node lists, job
  name). Only the call's own jobs are looked at, each at most once a minute.
- `nice` is the call's: its jobs get it at every place (in place of a place's own `--nice`), it
  doesn't keep jobs from moving between places, and a moved job keeps it.
- GPU nodes are for jobs that need a GPU: a CPU-only job goes to an alternative of GPU nodes
  only if its entry says `"cpuJobs": true`, leaves each idle GPU its share (gpuHeadroom: [CPUs,
  GB] per idle GPU, by default the node's share; 0 for none), never goes to a partition where
  GPU jobs of the user are pending (held ones too), and doesn't wait in the queue of GPU nodes.
- If Slurm doesn't answer the view's queries, jobs go to the main place and none moves until
  the next view.

## Task runners

On by default in cluster mode on Linux (`task_runner=False` turns them off); bash tasks and jobs
with in-job slots have none (`jobs.runner`, `jobs.claims`, `jobs.isolation`). call.json names
`petakit5d.jobs.runner:run` as the runner of the job entry point. A job's process runs its own
batch's tasks in order, then unclaimed tasks of the call's other batches, from the last task
backwards, while it has twice its longest task and a minute left (SLURM_JOB_END_TIME): only of
batches whose jobs get the same CPUs and GPU (call.json's keys), never of the batch the master
runs itself, of a batch the master submitted again, or a task whose output existed or whose
input was missing when the call started. A claimed task whose output (or the final output)
exists is skipped.

- Claims are symbolic links in `<call folder>/claims/`: `t<task>` (1-based) points to its
  owner `<job>_<array task>:<node>:<pid>:<slot>`; `f<task>` marks a task that failed there,
  `r<batch>` a batch the master submitted again, `s<job>_<array task>` a job that took tasks
  of other batches. A claim is live while its owner's job is running or pending (one
  `squeue -j` of the owners' jobs per round of the master), or, outside a job, while the
  owner's process runs on this node; when squeue doesn't answer, nothing is released.
- Between tasks the process state is given back (`jobs.isolation.ProcessState`: folder,
  environment, sys.path, warning filters, CPU affinity, floating-point environment, umask,
  locale, signal handlers, stdio, ...; a task runs in a copy of the context, so NumPy's print
  options and error state stay with it). The job goes on in a new process (`os.execv`, same
  pid) after a failed task (marked failed; the job takes no more tasks of other batches),
  when the process holds more than 0.5 GB above its memory after its first task
  (PETAKIT5D_RUNNER_MEM_GB), has 64 more open files or 256 more threads, or has run tasks for
  10 minutes (PETAKIT5D_RUNNER_MINUTES).
- The master (`jobs.runner.cluster_runners`, the cluster loop's runner object) leaves a batch
  while live runners hold its tasks (a first submission: all of them), removes the claims of
  owners that ended before it submits a batch again (and marks it resubmitted: its job runs
  it as without runners), claims the tasks of its own batch before it runs them, keeps the
  call's jobs in the queue to maxJobNum while jobs run other batches' tasks, and cancels the
  pending jobs of batches whose tasks are done or taken (`scancel --state=PENDING`).
- The isolation harness (`python -m petakit5d.jobs.isolation run <spec> <work> <results>`)
  runs task records in fresh processes and in one process in random orders and compares the
  outputs (`jobs.isolation.run_harness`).

## Task records

```python
Task("package.module:function", args=[...], kwargs={...}, path=[...])
Task.of(function, *args, **kwargs)  # import path from the function; partials too
```

A task record is JSON:

```json
{"func": "package.module:qualified.name", "args": [...], "kwargs": {...}, "path": ["/import/root"]}
{"bash": "command"}
```

- `func` names a callable that the job can import: module, a colon, the qualified name.
  `Task.of` fills `path` with the import root of a module outside petakit5d, which a job
  adds to `sys.path` only if the module can't be imported without it. A function of the
  script run as `__main__` is imported by the script's name (guard its top-level code with
  `if __name__ == "__main__":`, as for multiprocessing). Lambdas and nested functions can't
  be tasks.
- Arguments are JSON data: `None`, bool, int, float, str, lists, dicts with str keys.
  Tuples become lists, `os.PathLike` strings, NumPy values Python numbers (`tolist()`).
  Every record goes through a JSON round trip, so a task gets the same values in local mode
  and in a job; floats keep every bit, and inf and NaN are kept.
- Values that MATLAB passed as text keep MATLAB's rounding: the area wrapper rounds them
  (`%.10f`, `mat2str`) before it makes the task, with the helpers of `petakit5d.core`.

## Call folders and the job entry point

The cluster backend writes the call into a folder (`tasks.write_call`;
`tasks.new_call_folder(parent, uuid)` makes `<parent>/pk5d_call_<yyyyMMdd_HHmmss_SSS>_<uuid6>`):

```
<folder>/call.json     {"format": "petakit5d.jobs.call", "version": 1, "n_tasks", "language",
                        "task_batch_num", "final_out_fullpath", "bash_launch_str", "uuid",
                        "created", "python", "settings", "offsets": [byte offset of each line]}
<folder>/tasks.ndjson  one line per task: {"index": i (1-based), "input", "output", "task": record}
```

`call.json` is written last, by rename. A job runs

```
<PythonLaunchStr> -m petakit5d.jobs.exec <call folder> <batch>
```

with `batch` 1-based (the real batch number; Slurm's array index is `batch mod 5000`). The
entry point runs the batch's tasks in order and stops at the first exception (traceback on
stderr, exit status 1), as a MATLAB job stops at an error; a bash batch is one shell command,
`<BashLaunchStr>; cmd1;cmd2`, with the shell's exit status. Its last line is
`Elapsed time is <s> seconds.` Exit status 2: wrong arguments or no call folder. If
`call.json` has `"runner": "package.module:function"`, the batch goes to
`runner(call_folder, batch)` instead (the task runners plug in there). The entry point
imports only the standard library and the task's own module before the task.

Before any task, the entry point sets `OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`,
`MKL_NUM_THREADS` and `NUMEXPR_NUM_THREADS` to the CPUs the process may use (affinity and
cgroup quota), each only if it isn't set already: importing NumPy with 24 CPUs otherwise
spends 1.4 to 2.0 CPU-s on OpenBLAS's thread pool, which task processes packed on a node take
from each other. A slot manager that gives a process fewer CPUs than it may use sets them
itself.

## For the backends

- `tasks.Call(inputs, outputs, records, settings)`: what every backend gets; `settings` are the
  effective parameters by Python name (config applied, derived values set, `uuid` set).
  `call.batch(b)` gives the 0-based task indices of batch `b`, `call.batch_num` the batches.
- `framework.run_local(call)`: the local loop (also the fallback when there is no cluster).
- `framework.run_batch_in_process(call, indices, names)`: runs tasks in this process as the
  master does (retry once, then raise); returns the seconds.
- `framework.input_exist(call, out_exist)`: MATLAB's input check for the tasks without output.
- `files.batch_file_exist(paths)`, `files.wait_for_files(groups, wait_time)`.
- Cluster mode calls `petakit5d.jobs.cluster.run_cluster(call) -> list[bool]`; the local pool
  `petakit5d.jobs.local.run_pool(call, worker_num=...)` and, for `cluster_type='parfor'`,
  `local.run_parfor(call) -> list[bool]`.

## Task metrics

With `task_metrics=True` every task run leaves a record of what it used, in MATLAB's folders
and format, which MATLAB's read_task_metrics reads (`jobs.metrics`; the full record format is
in its docstring). Off by default, and then the metrics module isn't imported. A metrics error
never fails a task: it is one warning per process, and the task has no record.

- The run folder is `<jobLogDir>/task_metrics/<yyyyMMdd_HHmmss>_<function>_<id8>/` when
  jobLogDir is an existing absolute folder, else `<first input folder>/job_logs/task_metrics/...`.
  run_tasks makes it before any backend (`metrics.run_start(call, "python",
  poolWorkerNum=...)`, which sets `call.settings["task_metrics"]` to the folder, so call.json
  carries it into the jobs) and writes the run_end record after the backend returns
  (`metrics.run_end(call, done)`).
- Records: `task_<node>_<pid>_<process start ticks>.ndjson`, one JSON line per task (per batch
  for batches run in this process), one file per process; `run.ndjson` (the master's run_start,
  run_end and submit lines).
- Modes: `local` and `master` (`run_batch_in_process`), `pool` (pool workers), `slurm_python`
  (the job entry point: one record per Python task, one per bash batch with the shell's
  status). The hooks are `metrics.task_start(call, indices, mode) -> record` and
  `metrics.task_stop(record, outcome)`, outcome None, the exception or an exit status.
- Reading: `metrics.read_task_metrics(path, slurm_accounting=False) -> (tasks, runs)`, lists of
  dicts with MATLAB's columns.

## Publishing results and duplicate writers

A task is done when its output exists, so a task makes its secondary outputs first (MIPs,
flags, info files) and publishes its main result last, by one rename from a temporary name
(`<final>_<uuid>`). `petakit5d.jobs.publish` does that as MATLAB's publish_result and
keep_other_result do: `publish_result` moves the temporary result into place without ever
replacing or deleting a published one (`rename_no_replace`: a hard link for files, a rename
for folders, so that of two writers the first wins and the second keeps the other result),
`keep_other_result` removes a writer's own temporary result when another was published.
`petakit5d.jobs.writers` keeps records of who writes a result
(`zarr_flag/<name>_<uuid>/writer.txt`): `recorded_writer` (a context manager: registered,
released in `finally`), `wait` (a second writer waits while the first is alive, from squeue
and the process), and `orphaned` (block jobs of a writer that is gone skip their work). The
job entry point turns SIGTERM into an exception, so that the record of a cancelled job is
released at once. The module docstrings hold the rules.

## File checks

`files.batch_file_exist(paths, out_fullpath=None, use_parpool=False, use_java=True)` gives
MATLAB's `exist(path, 'file') > 0` for each path (files and folders, symbolic links followed,
trailing separators ignored, `~` expanded). On Linux each parent folder is opened (more than
100 paths in it: listed) before its paths are checked by name, so that files other nodes
made are seen at once on NFS; elsewhere the paths are checked one by one. A relative path is
relative to the current folder (MATLAB also finds names on its path). With a list of lists
the first path of each is checked. `out_fullpath` saves the result as a .mat file; MATLAB's
known issue is kept: a call with an existing result file fails.

`files.wait_for_files(groups, wait_time) -> bool` waits up to `wait_time` seconds until all
files of one group exist, looking after each second (less often when a look takes more than
a tenth of a second); a name with `*` is a pattern, matched as MATLAB's dir matches it
(`*` only; `?` is a plain character there).
