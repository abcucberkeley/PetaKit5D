# The job framework: run_tasks, task records, the job entry point, local mode

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
'parfor'` goes to the process pool; cluster mode without Slurm (sinfo fails) prints a warning
and runs locally; without cluster and with `local_worker_num > 1` the process pool runs
first and the local loop finishes what it left.

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
