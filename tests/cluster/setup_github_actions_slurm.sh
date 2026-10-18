#!/bin/bash
# Set up a single-node Slurm cluster (controller and compute node on this
# machine) for the Cluster tests in GitHub Actions, see
# github_actions_slurm.json. Run as root on Ubuntu 24.04:
#   sudo bash tests/cluster/setup_github_actions_slurm.sh
# There is no accounting database: PetaKit5D only uses sbatch, squeue,
# scancel, sinfo, and scontrol.

set -euo pipefail

export DEBIAN_FRONTEND=noninteractive
apt-get update -q -o Acquire::Retries=3
apt-get install -y -q -o Acquire::Retries=3 --no-install-recommends slurmctld slurmd slurm-client munge

mkdir -p /etc/slurm /var/log/slurm /var/spool/slurmctld /tmp/slurmd
# The node has 4 CPUs and 8 GB whatever the runner has (config_overrides),
# the test jobs are small.
cat > /etc/slurm/slurm.conf << 'EOF'
ClusterName=ci
SlurmctldHost=localhost
SlurmUser=root
AuthType=auth/munge
ProctrackType=proctrack/pgid
SelectType=select/cons_tres
SelectTypeParameters=CR_Core_Memory
DefMemPerCPU=256
# Without an accounting database, the main scheduler of slurmctld 23.11
# rejects jobs ("invalid account") and only backfill starts them
SchedulerParameters=bf_interval=1
SlurmdParameters=config_overrides
ReturnToService=2
StateSaveLocation=/var/spool/slurmctld
SlurmdSpoolDir=/tmp/slurmd
SlurmctldLogFile=/var/log/slurm/slurmctld.log
SlurmdLogFile=/var/log/slurm/slurmd.log
NodeName=localhost CPUs=4 RealMemory=8000 State=UNKNOWN
PartitionName=ci Nodes=localhost Default=YES MaxTime=INFINITE State=UP
EOF
# slurmd 23.11 needs a cgroup plugin, but nothing here uses one (pgid process
# tracking, no task/cgroup). Unlike cgroup/v2, the v1 plugin doesn't need
# cgroup access to load, which slurmd lacks when not started by systemd.
echo 'CgroupPlugin=cgroup/v1' > /etc/slurm/cgroup.conf

# the munge package starts munged where systemd runs
if ! pgrep -x munged > /dev/null; then
    install -d -o munge -g munge /run/munge
    runuser -u munge -- munged
fi
slurmctld
slurmd

# wait until the node runs jobs
for i in $(seq 1 30); do
    timeout 60 srun --time=1 true 2> /dev/null && break
    if [ "$i" -eq 30 ]; then
        tail -n 50 /var/log/slurm/slurmctld.log /var/log/slurm/slurmd.log
        exit 1
    fi
    sleep 2
done
sinfo
