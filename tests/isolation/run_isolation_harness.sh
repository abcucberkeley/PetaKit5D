#!/bin/bash
# Isolation harness: run the tasks of tests/isolation/isolation_prepare.m in
# fresh MATLAB processes (one per task) and in single MATLAB processes in
# several orders, with and without a reset between tasks, and compare all
# outputs bit for bit with those of the fresh processes. Run it on a compute
# node (inside srun/sbatch); everything runs on that node, because results can
# differ between CPU models.
#
# usage: run_isolation_harness.sh -c <PetaKit5D folder> -o <output folder> [options]
#   -c  PetaKit5D folder (code; tests/isolation is used from there)
#   -o  output folder (shared storage): template/, results/<run>/, records/, logs/, summary.txt
#   -w  work folder (default /tmp/pk5d_isolation_work); node-local, the same absolute path
#       for all runs, because the tasks' files hold absolute paths
#   -r  runs (default none:none:fwd,proto:prototype:rev,full1:full:1,full2:full:2,full3:full:3):
#       comma separated <name>:<reset>:<order>, reset none|prototype|full|full_gpu or a
#       MATLAB function name called as fcn(st0) (see isolation_reset), order fwd|rev|<seed>
#   -s  scenarios (default: all), comma separated names (see isolation_prepare)
#   -g  1|0 include GPU scenarios (default: 1 if nvidia-smi sees a GPU)
#   -H  1 to include the heavy scenarios (PSF analysis, ~22 GB)
#   -D  folder that holds PetaKit5D_demo_cell_image_dataset: also the scenarios on two demo
#       tiles (copied into the work folder)
#   -F  0 to skip the fresh runs (reuse results/fresh and records/fresh_* in -o)
#   -j  fresh processes at once (default 1; MATLAB's start and exit take ~30 s per task);
#       GPU tasks always run one at a time
#   -P  0 to skip prepare (reuse template/ in -o)
#   -m  MATLAB command (default: matlab -nojvm -batch)
#   -x  MATLAB statement run in each single-process run instead of isolation_run_tasks,
#       with $TASKS (tasks .mat), $ORDER (order vector), $RUN (run name) and $RECORD
#       (record .mat) replaced, e.g. to run the tasks through the job framework:
#       "isolation_run_framework('$TASKS', 'order', $ORDER, 'frameworkArgs', {'parseCluster', false})"
#
# example (abc cluster):
#   srun -c8 --mem=64G --gres=gpu:1 -p abc --qos abc_high -A co_abc --time=04:00:00 \
#     bash tests/isolation/run_isolation_harness.sh -c $PWD -o /clusterfs/nvme/$USER/isolation_run1

CODE=""; OUT=""; WORK=/tmp/pk5d_isolation_work
RUNS="none:none:fwd,proto:prototype:rev,full1:full:1,full2:full:2,full3:full:3"
SCEN=""; GPU=""; HEAVY=0; DEMO=""; DOFRESH=1; DOPREP=1; MATLAB="matlab -nojvm -batch"; CUSTOM=""; JOBS=1
while getopts "c:o:w:r:s:g:H:D:F:P:m:x:j:" opt; do
    case $opt in
        c) CODE=$OPTARG;; o) OUT=$OPTARG;; w) WORK=$OPTARG;; r) RUNS=$OPTARG;; s) SCEN=$OPTARG;;
        g) GPU=$OPTARG;; H) HEAVY=$OPTARG;; D) DEMO=$OPTARG;; F) DOFRESH=$OPTARG;; P) DOPREP=$OPTARG;;
        m) MATLAB=$OPTARG;; x) CUSTOM=$OPTARG;; j) JOBS=$OPTARG;;
        *) echo "unknown option"; exit 1;;
    esac
done
if [ -z "$CODE" ] || [ -z "$OUT" ]; then
    echo "usage: $0 -c <PetaKit5D folder> -o <output folder> [options]"; exit 1
fi
if [ -z "$GPU" ]; then
    GPU=0
    if nvidia-smi -L > /dev/null 2>&1; then GPU=1; fi
fi
GPUBOOL=false; [ "$GPU" = "1" ] && GPUBOOL=true
HEAVYBOOL=false; [ "$HEAVY" = "1" ] && HEAVYBOOL=true
SCENCELL="{}"
if [ -n "$SCEN" ]; then SCENCELL="{'${SCEN//,/\',\'}'}"; fi

mkdir -p "$OUT/results" "$OUT/records" "$OUT/logs" "$OUT/template"
WORK=${WORK%/}
TASKS="$WORK/tasks.mat"
# the tasks start in the work folder, not in the code folder the framework
# changes to, so that a task that changes the current folder is seen
INIT="cd('$CODE'); setup([]); addpath('$CODE/tests'); addpath('$CODE/tests/isolation'); cd('$WORK');"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$OUT/logs/harness.log"; }

log "host $(hostname), $(nproc) CPUs, GPU $GPU, code $CODE, work $WORK"
nvidia-smi --query-gpu=index,name,memory.used,memory.total,utilization.gpu --format=csv 2>/dev/null | tee -a "$OUT/logs/harness.log"

restore() {
    mkdir -p "$WORK"
    rsync -a --delete "$OUT/template/" "$WORK/"
}
collect() {
    rm -rf "$OUT/results/$1"
    mkdir -p "$OUT/results/$1"
    rsync -a --exclude=tasks.mat "$WORK/" "$OUT/results/$1/"
}

if [ "$DOPREP" = "1" ]; then
    log "prepare"
    rm -rf "$WORK"
    mkdir -p "$WORK"
    $MATLAB "$INIT isolation_prepare('$WORK', 'scenarios', $SCENCELL, 'gpu', $GPUBOOL, 'heavy', $HEAVYBOOL, 'demoData', '$DEMO');" \
        > "$OUT/logs/prepare.log" 2>&1 || { log "prepare failed, see $OUT/logs/prepare.log"; exit 1; }
    rsync -a --delete "$WORK/" "$OUT/template/"
fi
NT=$($MATLAB "a = load('$OUT/template/tasks.mat', 'tasks'); fprintf('NTASKS=%d\n', numel(a.tasks));" 2>/dev/null | sed -n 's/^NTASKS=//p')
log "$NT tasks"

run_fresh() {
    $MATLAB "$INIT isolation_run_tasks('$TASKS', 'order', $1, 'recordFn', '$OUT/records/fresh_$(printf %03d $1).mat');" \
        > "$OUT/logs/fresh_$(printf %03d $1).log" 2>&1
}

if [ "$DOFRESH" = "1" ]; then
    # CPU tasks up to $JOBS at once (each process still has the node's thread
    # count, which is what the results depend on), then the GPU tasks one by one
    GPUTASKS=" $($MATLAB "a = load('$OUT/template/tasks.mat', 'tasks'); fprintf('GPUTASKS=%s\n', num2str(find([a.tasks.gpu])));" 2>/dev/null | sed -n 's/^GPUTASKS=//p') "
    log "fresh: one MATLAB process per task ($JOBS at once; GPU tasks:$GPUTASKS)"
    restore
    rm -f "$OUT"/records/fresh_*.mat
    for i in $(seq 1 "$NT"); do
        if [[ "$GPUTASKS" == *" $i "* ]]; then continue; fi
        while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
        run_fresh "$i" &
    done
    wait
    for i in $GPUTASKS; do
        run_fresh "$i"
    done
    collect fresh
fi

IFS=',' read -ra RUNLIST <<< "$RUNS"
for spec in "${RUNLIST[@]}"; do
    IFS=':' read -r name reset order <<< "$spec"
    case $order in
        fwd) ORDER="1 : $NT";;
        rev) ORDER="$NT : -1 : 1";;
        *) ORDER="randperm(RandStream('twister', 'Seed', $order), $NT)";;
    esac
    case $reset in
        none|prototype|full|full_gpu) RESETARG="'$reset'";;
        *) RESETARG="@$reset";;
    esac
    log "run $name: reset $reset, order $order"
    restore
    if [ -n "$CUSTOM" ]; then
        STMT=${CUSTOM//\$TASKS/$TASKS}
        STMT=${STMT//\$ORDER/[$ORDER]}
        STMT=${STMT//\$RUN/$name}
        STMT=${STMT//\$RECORD/$OUT/records/$name.mat}
    else
        STMT="isolation_run_tasks('$TASKS', 'order', [$ORDER], 'reset', $RESETARG, 'recordFn', '$OUT/records/$name.mat');"
    fi
    /usr/bin/time -v $MATLAB "$INIT $STMT" > "$OUT/logs/run_$name.log" 2>&1
    collect "$name"
done

log "compare"
$MATLAB "$INIT isolation_report('$OUT');" > "$OUT/logs/report.log" 2>&1
cat "$OUT/summary.txt"
log "done"
