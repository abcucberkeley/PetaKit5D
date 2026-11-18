# Running on a Slurm cluster

Every dataset-level function (deskew/rotation, decon, stitching, the tools, the pipeline) runs
its per-file or per-block tasks through the job framework (`petakit5d.jobs.run_tasks`, the port
of MATLAB's `generic_computing_frameworks_wrapper`), and takes its parameters:

| Python | MATLAB | Default | Meaning |
|---|---|---|---|
| `parse_cluster` | `parseCluster` | `True` | submit Slurm jobs; without Slurm (no `sinfo`) the tasks run in this process |
| `config_file` | `configFile` | `""` | a cluster configuration file (below) |
| `master_compute` | `masterCompute` | `True` | the calling process runs tasks too |
| `cpus_per_task`, `mem_allocate` | `cpusPerTask`, `memAllocate` | per function | what a task needs |
| `job_log_dir` | `jobLogDir` | `"../job_logs"` | the job logs, and the call folders of the jobs |
| `gpu_job` | `GPUJob` | `False` | GPU jobs (decon with `gpu_job=True`) |

`parse_cluster=False` runs the tasks in this process (local mode), as in MATLAB;
`local_worker_num` > 1 runs them on a pool of local processes first. The full list of the
framework's parameters, with MATLAB's names and defaults, and its exact behaviour are in
[jobs.md](../jobs.md).

## The configuration file

The cluster configuration is the JSON file MATLAB PetaKit5D uses (a `.mat` file also works). Its
keys are the framework's parameters in MATLAB's names, matched without regard to case, and
override the caller's. Only one key is new: **`PythonLaunchStr`**, what a job runs to start
Python, in place of MATLAB's `MatlabLaunchStr` and the compiled master.

```json
{
  "PythonLaunchStr": "source /opt/envs/petakit5d/bin/activate; python",
  "BashLaunchStr": "module load parallel",
  "SlurmParam": "-p compute --qos normal -n1 --mem-per-cpu=21418M",
  "jobLogDir": "/mnt/projects/job_logs",
  "maxCPUNum": 24,
  "minCPUNum": 1,
  "wholeNodeJob": true,
  "memPerCPU": 20.9,
  "GNUparallel": true,
  "jobTimeLimit": 12,
  "maxJobNum": 5000
}
```

- **`PythonLaunchStr`**: shell text before the interpreter (module loads, activating the
  environment, exports) followed by the interpreter. Each job runs
  `<PythonLaunchStr> -u -m petakit5d.jobs.exec <call folder> <batch>`. Empty or left out: the
  interpreter of the process that submits, which works when the environment is on storage every
  node sees.
- **`SlurmParam`**: partition, QOS, account, memory per CPU, as passed to `sbatch`. The port has
  no site default: it comes from the file.
- **`jobLogDir`**: a folder every node sees. The framework writes each call's task list there
  (`pk5d_call_*`, or in `tmpDir` when given) and removes it when the call is done.
- **`memPerCPU`, `maxCPUNum`, `minCPUNum`**: the memory per core in GB and the CPUs per job,
  which shape the jobs as in MATLAB.
- **`wholeNodeJob`, `GNUparallel`, `paraJobNum`**: several tasks per job, one process per task
  (the port's slot manager replaces GNU parallel; `BashLaunchStr` still runs first).
- **`jobTimeLimit`** (hours), **`maxJobNum`**, **`SlurmConstraint`**, **`nice`**: as in MATLAB.
- **`SlurmAlternatives`**: other partitions to place jobs on, in order of preference
  ([jobs.md](../jobs.md), "Jobs on more than one place").
- **Ignored keys:** `MatlabLaunchStr`, `MCRParam`, `MCCMasterStr`, `MCRCacheRoot` and `mccMode`
  are accepted and have no effect, so one file can serve MATLAB and Python. MATLAB, in turn,
  keeps `PythonLaunchStr` as an unused field.

Pass the file to any dataset function:

```python
from petakit5d.decon import decon_dataset

decon_dataset(
    ["/mnt/projects/embryo_2024/DSR/"],
    psf_fullpaths=["/mnt/projects/psfs/488.tif"],
    xy_pixel_size=0.108,
    dz=0.3,
    decon_iter=2,
    parse_cluster=True,
    config_file="/mnt/projects/configs/cluster_cpu.json",
)
```

or from the command line, `pk5d decon_dataset ... --config-file /mnt/projects/configs/cluster_cpu.json`,
or as MATLAB's string arguments through `pk5d-master` (the GUI's form).

## How the jobs run

As in MATLAB: tasks are grouped into batches, each batch is one array job (`sbatch --array`),
the caller watches the outputs and the queue, resubmits failed batches with more CPUs (up to
`maxTrialNum` tries), and runs batches itself when `masterCompute` is on. A task is done when
its output exists; outputs are published by one rename at the end, so a half-written result
never counts. The job logs keep MATLAB's layout (`job_<jobid>_<taskid>.out`), with
`Python command:` where MATLAB's said `Matlab command:`.

The differences from MATLAB's framework are all in speed: a task starts in about 0.1 s (the
compiled master took 3 to 4 s), a job's process runs several tasks, and on Linux idle jobs take
waiting tasks of other batches with the same resources (task runners; `task_runner=False` turns
them off). Details: [jobs.md](../jobs.md).

## GPU jobs

GPU work (decon with `gpu_job=True`) runs as GPU jobs with a configuration file of its own, as
in MATLAB: its `SlurmParam` names the GPU partition and the GPUs (`--gres=gpu:1`), and `GPUNum`
the GPUs per job. A job without GPUs of its own uses none. In jobs with several slots each task
gets its own GPU (`CUDA_VISIBLE_DEVICES`). CPU jobs may also use GPU nodes through
`SlurmAlternatives` entries with `"cpuJobs": true`, which leave each idle GPU its share of the
node.

## Testing a cluster setup

The test suite has tests that submit real jobs (`cluster` marker). They run when
`PETAKIT5D_TEST_SLURM_CONFIG` names a configuration file and `PETAKIT5D_TEST_SLURM_SCRATCH` a
folder every node sees:

```
PETAKIT5D_TEST_SLURM_CONFIG=/mnt/projects/configs/cluster_cpu.json \
PETAKIT5D_TEST_SLURM_SCRATCH=/mnt/projects/scratch \
python -m pytest python/tests -m cluster
```
