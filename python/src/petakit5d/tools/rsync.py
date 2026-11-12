"""Parallel rsync: copy a folder with rsync calls run as bash tasks of the job framework.

The port of ``XR_parallel_rsync_wrapper`` (and ``batch_file_exist_cluster``). The files of the
source are shuffled (with the permutation of MATLAB's ``rng(1)``; the global random state is
left alone), split into batches, and each batch is one bash task that rsyncs its files one by
one; a batch counts as done once its last file exists in the destination. Up to five rounds
copy what is still missing; a last ``rsync -avP`` of the whole folder finishes the copy.

As MATLAB since 03510da: every path is one single-quoted word (``shell_quote``, a leading ``~``
expanded), the batch script takes pairs of a file and its folder, and the final rsync gets
its paths unchanged. The one difference: the batch script is part of the task's command
(``bash -c <script> rsync_batch_files <in> <out folder> ...``, the text of MATLAB's
rsync_batch_files.sh), so the task works from any current folder and install location and
needs no file next to the module (MATLAB takes it from the wrapper's folder, fix 9047731).
"""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
import time
from collections.abc import Sequence

from petakit5d.params import Param, register

__all__ = [
    "BATCH_SCRIPT",
    "batch_file_exist_cluster",
    "batch_size_for",
    "shell_quote",
    "parallel_rsync",
    "randperm_seed1",
    "rsync_batch_command",
]

# rsync_batch_files.sh: its arguments are pairs of a source file and its destination folder
BATCH_SCRIPT = """n=$(( $# / 2 ))
i=0
while [ $# -ge 2 ]; do
    i=$(( i + 1 ))
    echo "$1"
    echo "$2"
    if [ ! -d "$2" ]; then
        mkdir -p "$2"
    fi
    rsync -aP "$1" "$2/"
    echo $i / $n
    shift 2
done"""

# file lists of this many files or more are checked by jobs (batch_file_exist_cluster)
CLUSTER_CHECK_FILES = 100000
# files per job of batch_file_exist_cluster
CLUSTER_CHECK_BATCH = 50000
# rounds of copying before the final rsync
ROUNDS = 5


def randperm_seed1(n: int) -> list[int]:
    """``randperm(RandStream('twister', 'Seed', 1), n)`` as 0-based indices.

    MATLAB's randperm(n) is the stable sort order of ``rand(1, n)``; MT19937 seeded with 1
    gives the same doubles in MATLAB and NumPy's legacy ``RandomState(1)``.
    """
    import numpy as np

    u = np.random.RandomState(1).random_sample(n)
    return np.argsort(u, kind="stable").tolist()


def batch_size_for(n_files: int) -> int:
    """Files per batch: ``max(1, round(min(nF / 500, min(max(5, nF / 10000), 200))))``."""
    from petakit5d.core.dtypes import mround

    x = min(n_files / 500, min(max(5.0, n_files / 10000), 200.0))
    return max(1, int(mround(x)))


def _expand_home(p: str) -> str:
    # a leading ~ as bash expands it in an unquoted path
    if p == "~" or p.startswith("~/"):
        return os.environ.get("HOME", "") + p[1:]
    return p


def shell_quote(p: str) -> str:
    """A path as one bash word, as MATLAB's shell_quote: in single quotes (a ``'`` as
    ``'\\''``), with a leading ``~`` expanded."""
    return "'" + _expand_home(p).replace("'", "'\\''") + "'"


def rsync_batch_command(inputs: Sequence[str], out_folders: Sequence[str]) -> str:
    """The bash task of one batch: rsync each input into its destination folder, in order."""
    args = []
    for i, o in zip(inputs, out_folders, strict=True):
        args += [shell_quote(i), shell_quote(o)]
    return " ".join(["bash -c", shell_quote(BATCH_SCRIPT), "rsync_batch_files", *args])


def _strrep(s: str, old: str, new: str) -> str:
    # MATLAB's strrep (non-overlapping here: old ends with '/')
    return s.replace(old, new) if old else s


def _list_files(source: str, include_subdir: bool) -> list[str]:
    from petakit5d.core.paths import dir_recursive, matlab_dir

    # MATLAB's dir expands a leading ~; the paths keep the folder as given
    folder = source[:-1] if source.endswith("/") else source
    real = _expand_home(folder)
    if include_subdir:
        return [folder + p[len(real) :] for p in dir_recursive(real, "file", True)]
    return [f"{source}/{e.name}" for e in matlab_dir(f"{real}/*") if not e.isdir]


def batch_file_exist_cluster(
    file_fullpaths: Sequence[str],
    tmp_dir: str,
    *,
    parse_cluster: bool = True,
    mcc_mode: bool = False,
    config_file: str = "",
) -> list[bool]:
    """``batch_file_exist_cluster``: whether each path exists, checked by jobs of 50000 paths.

    Each job runs ``batch_file_exist`` and saves its result to
    ``<tmp_dir>/file_exist_flags/file_exist_<s>_<t>.mat``; ``parse_cluster`` and ``config_file`` go
    to the framework. MATLAB returns a cell of empties when a job fails (its caller then
    fails); this one raises.
    """
    from petakit5d.core.paths import mkdir_recursive
    from petakit5d.io import matfile
    from petakit5d.jobs.framework import run_tasks
    from petakit5d.jobs.tasks import Task

    del mcc_mode
    paths = [os.fspath(p) for p in file_fullpaths]
    n = len(paths)
    num_batch = math.ceil(n / CLUSTER_CHECK_BATCH)
    flag_dir = f"{_expand_home(tmp_dir)}/file_exist_flags/"
    mkdir_recursive(flag_dir)
    param_fullpath = f"{flag_dir}/parameters.mat"
    matfile.write_mat(
        param_fullpath,
        {
            "pr": {
                "fileFullpaths": paths,
                "tmpDir": tmp_dir,
                "parseCluster": bool(parse_cluster),
                "mccMode": False,
                "configFile": config_file,
            }
        },
        version="7.3",
    )
    flag_paths, tasks = [], []
    for i in range(num_batch):
        s = i * CLUSTER_CHECK_BATCH + 1
        t = min((i + 1) * CLUSTER_CHECK_BATCH, n)
        flag = f"{flag_dir}/file_exist_{s}_{t}.mat"
        flag_paths.append(flag)
        tasks.append(Task("petakit5d.jobs.files:batch_file_exist", [paths[s - 1 : t], flag, True]))
    done = run_tasks(
        [param_fullpath] * num_batch,
        flag_paths,
        tasks,
        tmp_dir=flag_dir,
        language="python",
        cpus_per_task=1,
        mem_allocate=20,
        max_job_num=math.inf,
        master_compute=True,
        parse_cluster=parse_cluster,
        config_file=config_file,
    )
    if not all(done):
        raise RuntimeError("batch_file_exist_cluster: not every check job finished")
    out: list[bool] = []
    for flag in flag_paths:
        v = matfile.read_mat(flag, "file_exist_mat")["file_exist_mat"]
        out.extend(bool(x) for x in _flat(v))
    return out


def _flat(v) -> list:
    import numpy as np

    if isinstance(v, list | tuple):
        return [x for e in v for x in _flat(e)]
    return np.asarray(v).ravel(order="F").tolist()


@register(
    "XR_parallel_rsync_wrapper",
    Param("source", "@(x) ischar(x)", "Source folder (with a trailing / its contents)."),
    Param("dest", "@(x) ischar(x)", "Destination folder."),
    Param("cpusPerTask", "@isnumeric", "CPUs per rsync job."),
    Param("numStream", "@isnumeric", "Number of jobs at once (streams of the transfer)."),
    Param(
        "SlurmParam",
        "@ischar",
        "Slurm options of the jobs.",
        differs="site-specific MATLAB default (from the cluster config)",
    ),
    Param("niceFactor", "@isnumeric", "Slurm nice value of the jobs (--nice)."),
    Param("includeSubdir", "@islogical", "Also copy the subfolders."),
    Param("parseCluster", "@islogical", "Run the batches as Slurm jobs."),
    Param("masterCompute", "@islogical", "This process also runs batches."),
    Param("mccMode", "@islogical", "MATLAB's compiled mode (accepted and ignored)."),
    Param("configFile", "@ischar", "Cluster config file."),
    task="T05",
)
def parallel_rsync(
    source: str,
    dest: str,
    *,
    cpus_per_task: float = 1,
    num_stream: float = math.inf,
    slurm_param: str = "",
    nice_factor: float = 0,
    include_subdir: bool = True,
    parse_cluster: bool = True,
    master_compute: bool = False,
    mcc_mode: bool = False,
    config_file: str = "",
) -> None:
    """Copy a folder with many rsync calls run as jobs, then one rsync of the whole folder.

    Without a trailing ``/`` the source folder itself is copied into ``dest``, as rsync does.
    """
    from petakit5d.core.ids import get_uuid
    from petakit5d.core.paths import fileparts
    from petakit5d.core.text import sprintf
    from petakit5d.jobs.files import batch_file_exist
    from petakit5d.jobs.framework import run_tasks

    slurm_param = sprintf("%s --nice=%d", slurm_param, nice_factor)
    t0 = time.perf_counter()
    uuid = get_uuid()

    source = source.replace("//", "/")
    inputs = _list_files(source, include_subdir)
    n_files = len(inputs)
    dest = dest.replace("//", "/")
    if not dest.endswith("/"):
        dest += "/"
    if not source.endswith("/"):
        _, d, ext = fileparts(source)
        source += "/"
        dest = f"{dest}{d}{ext}/"

    outputs = [_strrep(p, source, dest) for p in inputs]
    out_folders = [fileparts(p)[0] for p in outputs]

    tmp_dir = f"{_expand_home(dest)}/tmp/{uuid}/"
    os.makedirs(tmp_dir, exist_ok=True)

    # shuffled so that large files spread over the workers (the order of rng(1))
    perm = randperm_seed1(n_files)
    inputs = [inputs[i] for i in perm]
    outputs = [outputs[i] for i in perm]
    out_folders = [out_folders[i] for i in perm]

    for _ in range(ROUNDS):
        if len(out_folders) < CLUSTER_CHECK_FILES:
            exist = batch_file_exist(outputs)
        else:
            stamp = time.strftime("%Y%m%d%H%M%S")
            exist = batch_file_exist_cluster(
                outputs,
                f"{tmp_dir}/{stamp}/",
                parse_cluster=parse_cluster,
                mcc_mode=mcc_mode,
                config_file=config_file,
            )
        if all(exist):
            break
        keep = [i for i, e in enumerate(exist) if not e]
        inputs = [inputs[i] for i in keep]
        outputs = [outputs[i] for i in keep]
        out_folders = [out_folders[i] for i in keep]
        n_files = len(inputs)

        batch_size = batch_size_for(n_files)
        num_batch = math.ceil(n_files / batch_size)
        commands = []
        for i in range(num_batch):
            s, t = i * batch_size, min((i + 1) * batch_size, n_files)
            commands.append(rsync_batch_command(inputs[s:t], out_folders[s:t]))
        # each batch is checked by its last file, which the batch copies last: the framework
        # counts a batch as done once its checked file exists
        last = [min((i + 1) * batch_size, n_files) - 1 for i in range(num_batch)]
        run_tasks(
            [inputs[i] for i in last],
            [outputs[i] for i in last],
            commands,
            tmp_dir=tmp_dir,
            language="bash",
            slurm_param=slurm_param,
            cpus_per_task=cpus_per_task,
            max_job_num=num_stream,
            master_compute=master_compute,
            parse_cluster=parse_cluster,
            mcc_mode=mcc_mode,
            config_file=config_file,
        )

    # final check with one rsync of the whole folder (its output shown, its status ignored)
    sys.stdout.flush()
    sys.stderr.flush()
    try:
        subprocess.run(["rsync", "-avP", _expand_home(source), _expand_home(dest)], check=False)
    except OSError as err:
        print(f"rsync: {err}", file=sys.stderr)

    shutil.rmtree(tmp_dir)
    print(f"Elapsed time is {time.perf_counter() - t0:f} seconds.")
