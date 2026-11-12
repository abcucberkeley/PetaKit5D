"""Cluster configuration files: reading, MATLAB's merge rules and cluster_config_value.

Ported from the MATLAB suite: MasterComputeFitTest/configValuesIgnoreCase, and the config tests
of GenericComputingFrameworkTest (configFileOverridesParameters, configFileKeysIgnoreCase,
emptyJobLogDirInAConfigFileIsTheDefaultFolder, missingConfigFileErrors) and
SingleJobSubmitWrapperTest (configFileKeysIgnoreCase, emptyJobLogPathsInAConfigFileAreTheDefaults,
logPathsStayWhenSinfoFails) at the level of this layer (the framework calls are J02's and J07's).
"""

import json
import os
import shlex
import sys

import pytest

from petakit5d.jobs import config as cfg
from petakit5d.jobs.config import (
    ClusterConfig,
    ConfigCache,
    cluster_config_value,
    decode_json,
    is_empty,
    load_config,
    matlab_valid_name,
)
from petakit5d.jobs.shapes import job_cpus, single_job_cpus
from petakit5d.jobs.slurm import (
    DEFAULT_JOB_LOG_DIR,
    check_slurm_cluster,
    job_log_dir_or_default,
    job_log_names,
)
from petakit5d.testing.slurm_fake import fake_slurm  # noqa: F401  (fixture)

# the parameters of generic_computing_frameworks_wrapper this file uses: MATLAB name -> Python
NAMES = {
    "maxTrialNum": "max_trial_num",
    "memPerCPU": "mem_per_cpu",
    "masterCompute": "master_compute",
    "jobLogDir": "job_log_dir",
    "maxCPUNum": "max_cpu_num",
}


def write_json(path, value):
    path.write_text(value if isinstance(value, str) else json.dumps(value))
    return path


def test_config_values_ignore_case(tmp_path):
    # MasterComputeFitTest/configValuesIgnoreCase
    fn = write_json(tmp_path / "config.json", '{"MasterCompute": false, "memPerCPU": 7.6}')
    assert cluster_config_value(fn, "masterCompute", True) is False
    assert cluster_config_value(fn, "MemPerCPU", 20.9) == 7.6
    assert cluster_config_value(fn, "maxCPUNum", 24) == 24
    assert cluster_config_value("", "masterCompute", True) is True
    assert cluster_config_value(None, "masterCompute", True) is True
    assert cluster_config_value(tmp_path / "none.json", "masterCompute", True) is True


def test_config_values_of_a_mat_file(tmp_path):
    # MasterComputeFitTest/configValuesIgnoreCase, the .mat part (reading needs F07's module)
    pytest.importorskip(
        "petakit5d.io.matfile", reason="petakit5d.io.matfile (F07) is not there yet"
    )
    import scipy.io

    fn = tmp_path / "config.mat"
    scipy.io.savemat(fn, {"masterCompute": False})
    assert cluster_config_value(fn, "masterCompute", True) in (False, 0)
    assert bool(ClusterConfig.load(fn).get("MASTERCOMPUTE")) is False


def test_other_extensions_give_the_value(tmp_path):
    fn = write_json(tmp_path / "config.txt", '{"masterCompute": false}')
    assert cluster_config_value(fn, "masterCompute", True) is True
    # (MATLAB's ext switch is case sensitive)
    fn = write_json(tmp_path / "config.JSON", '{"masterCompute": false}')
    assert cluster_config_value(fn, "masterCompute", True) is True
    with pytest.raises(ValueError):
        ClusterConfig.load(fn)


def test_config_file_overrides_parameters(tmp_path):
    # GenericComputingFrameworkTest/configFileOverridesParameters: maxTrialNum 3 here, 1 in the file
    fn = write_json(tmp_path / "config.json", {"maxTrialNum": 1})
    merged = load_config(fn, cache=None).merge({"max_trial_num": 3, "mem_per_cpu": 20.9}, NAMES)
    assert merged == {"max_trial_num": 1, "mem_per_cpu": 20.9}


def test_config_file_keys_ignore_case_batch_jobs(tmp_path):
    # GenericComputingFrameworkTest/configFileKeysIgnoreCase: 50 GB at 5 GB per CPU is 10 CPUs,
    # rounded up to 12 to fit 24-CPU nodes (the default 20.9 GB per CPU would give 3)
    fn = write_json(tmp_path / "config.json", {"MemPerCPU": 5})
    p = load_config(fn, cache=None).merge({"mem_per_cpu": 20.9, "max_cpu_num": 24}, NAMES)
    assert p["mem_per_cpu"] == 5
    assert job_cpus(1, 50, p["mem_per_cpu"], 1, p["max_cpu_num"], True) == 12
    assert job_cpus(1, 50, 20.9, 1, 24, True) == 3


def test_config_file_keys_ignore_case_single_jobs(tmp_path):
    # SingleJobSubmitWrapperTest/configFileKeysIgnoreCase: single jobs don't round (10 CPUs)
    fn = write_json(tmp_path / "config.json", {"MemPerCPU": 5})
    p = load_config(fn, cache=None).merge({"mem_per_cpu": 20.9}, NAMES)
    assert single_job_cpus(1, 50, p["mem_per_cpu"], 1, 24) == 10


def test_empty_job_log_dir_in_a_config_file_is_the_default_folder(
    tmp_path,
    fake_slurm,  # noqa: F811
    monkeypatch,
):
    # GenericComputingFrameworkTest/emptyJobLogDirInAConfigFileIsTheDefaultFolder: "jobLogDir": ""
    # gave the jobs "-o  -e"; it is the default folder, ../job_logs (from the current folder) if
    # it is there, else <data folder>/job_logs
    data = tmp_path / "work" / "data"
    data.mkdir(parents=True)
    monkeypatch.chdir(data)
    fn = write_json(tmp_path / "config.json", {"jobLogDir": ""})
    p = load_config(fn, cache=None).merge({"job_log_dir": DEFAULT_JOB_LOG_DIR}, NAMES)
    assert p["job_log_dir"] == ""
    log_dir = job_log_dir_or_default(p["job_log_dir"])
    assert log_dir == "../job_logs"
    with pytest.warns(UserWarning, match="job log directory does not exist"):
        check = check_slurm_cluster(str(data), log_dir)
    assert check.parse_cluster
    assert check.job_log_dir == str(data / "job_logs").replace(os.sep, "/")
    assert check.job_log_fname == f"{check.job_log_dir}/job_%A_%a.out"
    assert check.job_log_error_fname == f"{check.job_log_dir}/job_%A_%a.err"
    assert (data / "job_logs").is_dir()
    # ../job_logs, where it is there
    (tmp_path / "work" / "job_logs").mkdir()
    check = check_slurm_cluster(str(data), log_dir)
    assert check.job_log_dir == str(tmp_path / "work" / "job_logs").replace(os.sep, "/")


def test_empty_job_log_paths_in_a_config_file_are_the_defaults(tmp_path):
    # SingleJobSubmitWrapperTest/emptyJobLogPathsInAConfigFileAreTheDefaults
    fn = write_json(
        tmp_path / "config.json", {"jobLogDir": "", "jobLogFname": "", "jobErrorFname": ""}
    )
    conf = load_config(fn, cache=None)
    assert job_log_dir_or_default(conf.get("jobLogDir")) == "../job_logs"
    assert job_log_names(job_log_dir_or_default(conf.get("jobLogFname"))) == (
        "../job_logs/job_%A_%a.out",
        "../job_logs/job_%A_%a.err",
    )


def test_missing_config_file_errors(tmp_path):
    # GenericComputingFrameworkTest/missingConfigFileErrors
    with pytest.raises(FileNotFoundError, match="does not exist"):
        load_config(tmp_path / "missing.json")
    assert not load_config("")
    assert not load_config(None)


def test_matlab_only_keys_are_accepted_and_ignored(tmp_path):
    keys = dict.fromkeys(cfg.MATLAB_ONLY_KEYS, "x")
    fn = write_json(tmp_path / "config.json", {**keys, "mccmode": True, "maxTrialNum": 2})
    conf = load_config(fn, cache=None)
    names = {**NAMES, "MatlabLaunchStr": "matlab_launch_str", "mccMode": "mcc_mode"}
    merged = conf.merge({"max_trial_num": 3, "matlab_launch_str": "m", "mcc_mode": False}, names)
    assert merged == {"max_trial_num": 2, "matlab_launch_str": "m", "mcc_mode": False}
    assert sorted(conf.ignored_keys()) == sorted([*cfg.MATLAB_ONLY_KEYS, "mccmode"])
    assert conf.unknown_keys(names) == []


def test_unknown_keys_have_no_effect(tmp_path):
    fn = write_json(tmp_path / "config.json", {"noSuchKey": 1, "maxTrialNum": 2})
    conf = load_config(fn, cache=None)
    assert conf.merge({"max_trial_num": 3}, NAMES) == {"max_trial_num": 2}
    assert conf.unknown_keys(NAMES) == ["noSuchKey"]
    # without a name map the keys of params are MATLAB's names
    assert conf.merge({"maxTrialNum": 3}) == {"maxTrialNum": 2}


def test_the_last_of_two_spellings_wins(tmp_path):
    fn = write_json(tmp_path / "config.json", '{"memPerCPU": 7.6, "MemPerCPU": 5}')
    conf = load_config(fn, cache=None)
    assert conf.get("memPerCPU") == 5
    assert conf.merge({"mem_per_cpu": 20.9}, NAMES) == {"mem_per_cpu": 5}
    fn = write_json(tmp_path / "config2.json", '{"MemPerCPU": 5, "memPerCPU": 7.6}')
    assert load_config(fn, cache=None).merge({"mem_per_cpu": 20.9}, NAMES) == {"mem_per_cpu": 7.6}


def test_the_first_of_two_equal_keys_wins(tmp_path):
    # jsondecode names a repeated key a_1, which matches no parameter
    fn = write_json(tmp_path / "config.json", '{"maxTrialNum": 1, "maxTrialNum": 2}')
    conf = load_config(fn, cache=None)
    assert conf.keys() == ["maxTrialNum", "maxTrialNum_1"]
    assert conf.merge({"max_trial_num": 3}, NAMES) == {"max_trial_num": 1}


def test_python_launch_str(tmp_path):
    assert ClusterConfig.empty().python_launch_str == shlex.quote(sys.executable)
    fn = write_json(
        tmp_path / "config.json", {"PythonLaunchStr": "module load x; /opt/env/bin/python"}
    )
    assert load_config(fn, cache=None).python_launch_str == "module load x; /opt/env/bin/python"
    fn = write_json(tmp_path / "config2.json", {"pythonlaunchstr": ""})
    assert load_config(fn, cache=None).python_launch_str == shlex.quote(sys.executable)


def test_values_and_emptiness():
    conf = ClusterConfig.from_mapping(
        decode_json('{"a": null, "b": [], "c": "", "d": [4, 32], "e": {"f g": 1}}')
    )
    assert [is_empty(conf.get(k)) for k in "abcd"] == [True, True, True, False]
    assert conf.get("d") == [4, 32]
    assert conf.get("e") == {"fG": 1}
    assert is_empty(0) is False
    assert not is_empty("x") and not is_empty([4]) and not is_empty({"a": 1})


def test_cache_reads_a_file_again_when_it_changes(tmp_path):
    fn = write_json(tmp_path / "config.json", {"maxTrialNum": 1})
    cache = ConfigCache()
    conf, fresh = cache.load(fn)
    assert fresh and conf.get("maxTrialNum") == 1
    conf, fresh = cache.load(fn)
    assert not fresh
    write_json(fn, {"maxTrialNum": 22})
    st = os.stat(fn)
    os.utime(fn, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))
    conf, fresh = cache.load(fn)
    assert fresh and conf.get("maxTrialNum") == 22
    cache.clear()
    assert cache.load(fn)[1]
    cfg.clear_cache()
    assert load_config(fn).get("maxTrialNum") == 22


# jsondecode's field names (MATLAB R2024b, oracle run r1: out/jsondecode_names.json)
JSONDECODE_KEYS = [
    "MemPerCPU", "mem per cpu", " lead", "trail ", "Two  spaces", "tab\tkey", "1abc", "_x",
    "a-b", "a.b", "end", "for", "a", "a", "a_1", "b", "B", "x" * 70, "", "caf\u00e9",
    "x y-z", "A b", "a 1",
]  # fmt: skip
JSONDECODE_FIELDS = [
    "MemPerCPU", "memPerCpu", "lead", "trail", "TwoSpaces", "tabKey", "x1abc", "x_x", "a_b",
    "a_b_1", "end", "for", "a", "a_2", "a_1", "b", "B", "x" * 63, "x", "caf_", "xY_z", "AB",
    "a1",
]  # fmt: skip


def test_keys_get_matlabs_field_names():
    text = "{" + ", ".join(f"{json.dumps(k)}: {i + 1}" for i, k in enumerate(JSONDECODE_KEYS))
    decoded = decode_json(text + ', "nested": {"b c": 1, "end": 2}}')
    assert list(decoded) == [*JSONDECODE_FIELDS, "nested"]
    assert list(decoded["nested"]) == ["bC", "end"]
    assert list(decoded.values())[:-1] == list(range(1, len(JSONDECODE_KEYS) + 1))


@pytest.mark.parametrize(
    "key, name",
    [
        ("SlurmParam", "SlurmParam"),
        ("mem per cpu", "memPerCpu"),
        ("1abc", "x1abc"),
        ("a-b", "a_b"),
        ("end", "end"),
        ("a z", "aZ"),
        (" lead", "lead"),
    ],
)
def test_matlab_valid_name(key, name):
    assert matlab_valid_name(key) == name
