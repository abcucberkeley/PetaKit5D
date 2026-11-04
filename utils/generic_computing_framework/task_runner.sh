# Job script of the job framework's task runners (see task_runner.m). The
# framework writes it for each call as <call folder>/run_tasks.sh, with the
# call's settings above this text (D, NT, BS, MASTER, MLAUNCH, LOCALRUN, UWTSEC,
# and WAITALL for the single jobs of generic_single_job_submit_wrapper).
#
# A runner is a MATLAB Runtime process (compiled jobs) or a MATLAB process
# (MATLAB jobs) that runs one task after another: the job script claims a
# task and starts a process for it with its task line (compiled) or the
# MATLAB command, and the process, after the task, resets itself to the state
# of a new process and claims the next task, first of its own batch, then
# (stealing) unstarted tasks of the call's other batches. The single jobs a
# MATLAB session submits with the same function and resources share a call
# folder (a pool, see task_runner_run): each job is a batch of one task, and
# tasks are added while jobs run.
#
# Claims are in $D/claims (shared by all jobs of the call): t<task> is a
# symbolic link that the claimer creates (atomic, also on NFS; it fails if
# the task is claimed) with its owner as the target
# (<job>_<array task>:<node>:<pid>:<slot>), f<task> marks a task that failed
# in that job (after its retry), r<batch> a batch the framework resubmitted
# (other jobs don't take its tasks). A task whose output exists is claimed
# and skipped.
#
# Usage:
#   bash run_tasks.sh job <own> <slots> <gpu 0|1> <steal 0|1|2|3> <retries>
#       a compiled job (or the compiled master): runs the tasks of <own>
#       with up to <slots> runners at once (GNU parallel), then steals (1:
#       from the last batch, 2: from the first one; 3 (a single job, one
#       slot): as 2, but first the job's own task runs by its task line in
#       a process of its own, as without runners, and a job whose task
#       another process has ends at once; it takes another batch's task
#       only while Slurm keeps that task's job pending (see waiting_tasks)
#       and cancels the job; only batches with the key of <own>, see
#       keys.txt); a task a runner doesn't finish runs
#       once more in a process of its own if retries > 1 (as GNU parallel
#       --retries 2 did), else it is marked failed. A slot with a failed
#       task stops stealing.
#   bash run_tasks.sh matlab <own> <steal>
#       a MATLAB job: the same with MATLAB processes; failed tasks are left
#       to the framework, which resubmits them with more CPUs
#   PK_PID=<pid> bash run_tasks.sh claim <own> <steal> <slot> <seconds> <status> <local> <out>
#       the next task of runner process <pid>: writes its index to <out>
#       (empty if none) and its arguments (compiled) to <local>/dump_<slot>;
#       or 'recycle <reason>' if the process has more open files (64) or
#       threads (256) than after its first task
# <own> is a batch number or the absolute path of a file of task indices.
# Files of the call folder:
#   tasks.txt     the task lines, one per task (empty: done when the call started)
#   tasks.idx     the byte offset of each line in tasks.txt (12 digits per line)
#   runnable.txt  the indices of the tasks with a line
#   keys.txt      per batch, the CPUs (threads) and GPU of the processes that
#                 run its tasks without runners ('master' for the batch the
#                 master runs); a runner takes only tasks of batches with the
#                 key of its own (MATLAB fixes how its FFTs use threads at a
#                 process's first FFT, and results depend on the thread count)
#   waiting.txt   the batches whose jobs have waited in the queue for 30 s
#                 (WAITALL=1: all, as in pools, whose tasks are added as
#                 their jobs are submitted)
#   jobs.txt      (pools) the Slurm job of each task: '<task> <job>_<array task>'
#   live/         (pools) a mark of the master (m) and of each task (t<task>)
#                 whose job may still use the pool; the last to leave
#                 removes the pool (see leave_pool)
#   tmps.txt      (local runs) the .tmp file and the output of each task
#   argdump/      a copy of the compiled master's run_mccMaster.sh with a
#                 stand-in master that writes the arguments it gets

C="$D/claims"
SELF="$D/run_tasks.sh"

owner() {
    # the claims' owner: job, node, process (PK_PID: the runner's) and slot
    echo "${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-0}}_${SLURM_ARRAY_TASK_ID:-0}:${SLURMD_NODENAME:-$HOSTNAME}:${PK_PID:-$$}:$1"
}

candidates() {
    # the tasks without a claim in the order to take them, after a line
    # 'n <own> <waiting> <other>' with their numbers: 'o <task>' for the tasks
    # of <own> (a batch, or a file of task indices), then (steal) 'a <task>'
    # for the runnable tasks of the batches not resubmitted whose processes
    # get the same threads and GPU as those of <own> (the same key in
    # keys.txt); <waiting> counts those of batches whose jobs wait in the
    # queue (waiting.txt, from the master; all with WAITALL=1)
    [ -d "$C" ] || return
    ls -f "$C" 2> /dev/null | awk -v bs="$BS" -v nt="$NT" -v own="$1" -v steal="$2" -v rf="$D/runnable.txt" \
        -v wf="$D/waiting.txt" -v kf="$D/keys.txt" -v allwait="${WAITALL:-0}" '
        { c[$0] = 1 }
        END {
            if (own ~ /^\//) {
                while ((getline x < own) > 0) { k++; ol[k] = x; ow[x] = 1 }
            } else {
                e = own * bs
                if (e > nt) e = nt
                for (x = (own - 1) * bs + 1; x <= e; x++) { k++; ol[k] = x; ow[x] = 1 }
            }
            for (i = 1; i <= k; i++) if (!(("t" ol[i]) in c)) { m++; om[m] = ol[i] }
            if (steal != 0) {
                while ((getline x < kf) > 0) { j++; key[j] = x }
                mykey = (own ~ /^\//) ? "" : key[own]
                while ((getline x < wf) > 0) w[x] = 1
                while ((getline x < rf) > 0) {
                    b = int((x - 1) / bs) + 1
                    if ((x in ow) || (("t" x) in c) || (("r" b) in c) || key[b] != mykey) continue
                    n++; o[n] = x
                    if (allwait == 1 || (b in w)) nw++
                }
            }
            print "n", m + 0, nw + 0, n + 0
            for (i = 1; i <= m; i++) print "o", om[i]
            if (steal == 1) { for (i = n; i >= 1; i--) print "a", o[i] } else { for (i = 1; i <= n; i++) print "a", o[i] }
        }'
}

task_line() {
    # the line of task $1 (tasks.idx has the byte offset of each line)
    local off
    if [ -f "$D/tasks.idx" ]; then
        off=$(dd if="$D/tasks.idx" bs=13 skip=$(( $1 - 1 )) count=1 status=none 2> /dev/null)
        [ -n "$off" ] && tail -c +$(( 10#$off + 1 )) "$D/tasks.txt" | head -n 1
    else
        sed -n "${1}{p;q}" "$D/tasks.txt"
    fi
}

task_runs() {
    # 0 if the task line of task $1 runs its task, 1 if it says its output
    # exists; the arguments it gives the compiled master go to file $2: the
    # line runs with a copy of the master's run_mccMaster.sh that starts a
    # stand-in, which writes them (no file: the task runs in a process of
    # its own)
    local line msg
    line=$(task_line "$1")
    [ -z "$line" ] && return 1
    rm -f "$2"
    [ -n "$MASTER" ] && line=${line//"$MASTER"/"$D/argdump/run_mccMaster.sh"}
    msg=$( ( export PK_DUMP="$2"; eval "$line" ) 2>&1 )
    [ -f "$2" ] && return 0
    [ -n "$msg" ] && echo "$msg"
    case "$msg" in
        *"already exists"*) return 1 ;;
    esac
    return 0
}

time_left() {
    # 0 if the job has the time for a task of about $1 seconds (0: none ran
    # yet; then if less than half of the job's time is used)
    [ -z "$SLURM_JOB_END_TIME" ] && return 0
    local now need
    now=$(date +%s)
    if [ "$1" -gt 0 ]; then
        need=$(( 2 * $1 + 60 ))
    else
        need=$(( (SLURM_JOB_END_TIME - ${SLURM_JOB_START_TIME:-$now}) / 2 ))
    fi
    [ $(( SLURM_JOB_END_TIME - now )) -ge "$need" ]
}

local_tmp() {
    # local runs: a .tmp file tells other framework calls that a task is
    # being processed; 1 if another one made it in the last UWTSEC seconds
    local tmp
    tmp=$(sed -n "${1}{p;q}" "$D/tmps.txt" | cut -f 1)
    [ -z "$tmp" ] && return 0
    if [ -f "$tmp" ] && [ $(( $(date +%s) - $(stat -c %Y "$tmp" 2> /dev/null || echo 0) )) -lt "$UWTSEC" ]; then
        return 1
    fi
    [ -d "$(dirname "$tmp")" ] && touch "$tmp"
    return 0
}

local_tmp_done() {
    # local runs: removes the .tmp file of a task whose output exists
    local tmp out
    tmp=$(sed -n "${1}{p;q}" "$D/tmps.txt" | cut -f 1)
    out=$(sed -n "${1}{p;q}" "$D/tmps.txt" | cut -f 2)
    [ -n "$tmp" ] && [ -e "$out" ] && rm -f "$tmp"
}

local_tmp_release() {
    # local runs: removes the .tmp file of a task given back unfinished
    local tmp out
    tmp=$(sed -n "${1}{p;q}" "$D/tmps.txt" | cut -f 1)
    out=$(sed -n "${1}{p;q}" "$D/tmps.txt" | cut -f 2)
    [ -n "$tmp" ] && [ ! -e "$out" ] && rm -f "$tmp"
}

claim_next() {
    # claims the next task: <own> <steal> <slot> <seconds> <status> <dump>
    # <out> [<own only>]; writes its index to <out> (empty if none), and to
    # <out>.more 'o' or 'a' (own or another batch's) and the numbers of
    # unclaimed own tasks, waiting batches' tasks and other batches' tasks
    # before it. With a dump file, tasks whose outputs exist are claimed and
    # skipped, and the arguments of the task go to the dump file; with '-'
    # the task's line checks its output itself. <own only> 1: a task of
    # <own> or none (the first task of a MATLAB single job)
    local own=$1 steal=$2 s=$3 est=$4 st=$5 dump=$6 out=$7 ownonly=${8:-0} kind f counts="0 0 0" wt=""
    local me="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-0}}_${SLURM_ARRAY_TASK_ID:-0}:${SLURMD_NODENAME:-$HOSTNAME}:${PK_PID:-$$}:$s"
    : > "$out"
    echo "o $counts" > "$out.more"
    while read -r kind f; do
        if [ "$kind" = n ]; then
            counts=$f
            continue
        fi
        [ -z "$f" ] && continue
        # tasks of other batches only if the job has the time for them
        if [ "$kind" = a ] && { [ "$ownonly" = 1 ] || ! time_left "$est"; }; then
            return 0
        fi
        # (pools) and while Slurm keeps their jobs in the queue
        if [ "$kind" = a ] && [ "$steal" = 3 ]; then
            [ -z "$wt" ] && wt=" $(waiting_tasks | tr '\n' ' ') "
            [[ "$wt" == *" $f "* ]] || continue
        fi
        ln -s "$me" "$C/t$f" 2> /dev/null || continue
        [ "$kind" = a ] && [ "$steal" = 3 ] && cancel_job_of "$f"
        if [ "$LOCALRUN" = 1 ] && ! local_tmp "$f"; then
            rm -f "$C/t$f"
            continue
        fi
        if [ "$dump" = - ] || task_runs "$f" "$dump"; then
            echo "$f" > "$out"
            echo "$kind $counts" > "$out.more"
            return 0
        fi
        echo "s $f" >> "$st"
        [ "$LOCALRUN" = 1 ] && local_tmp_done "$f"
    done < <(candidates "$own" "$steal")
    return 0
}

waiting_tasks() {
    # (steal 3, pools) the tasks whose jobs Slurm keeps in the queue: its
    # scheduler looked at them and could not start them (pending with a
    # reason: resources, priority, a limit), not jobs it has not looked at
    # yet (reason None), which may start at once, nor tasks whose jobs are
    # not recorded yet (being submitted)
    [ -s "$D/jobs.txt" ] || return 0
    timeout 10 squeue -r -h -u "${USER:-$(id -un)}" -t PENDING -o '%i %r' 2> /dev/null \
        | awk 'NR == FNR {if ($2 != "None") p[$1] = 1; next} ($2 in p) {print $1}' - "$D/jobs.txt"
}

cancel_job_of() {
    # (single jobs) the job of task $1, which this job took: cancelled while
    # it is pending (it would only start to end at once, at the front of the
    # queue); if it just started, its own claim failed and it ends anyway
    local job
    [ -f "$D/jobs.txt" ] || return 0
    job=$(awk -v k="$1" '$1 == k {j = $2} END {print j}' "$D/jobs.txt")
    [ -n "$job" ] && timeout 10 scancel --state=PENDING "$job" > /dev/null 2>&1
    return 0
}

recycle_reason() {
    # why runner process <pid> should end before its next task (nothing if
    # it shouldn't): it has 64 open files or 256 threads more than after its
    # first task (thread pools, loaded libraries); <local folder> keeps that
    # baseline. (The runner checks its memory itself, after its reset.)
    local pid=$1 L=$2 k v rest thr=0 fds fds0 thr0
    [ -r "/proc/$pid/status" ] || return 0
    while read -r k v rest; do
        [ "$k" = Threads: ] && thr=$v
    done < "/proc/$pid/status"
    fds=(/proc/"$pid"/fd/*)
    fds=${#fds[@]}
    if [ ! -f "$L/base_$pid" ]; then
        echo "$fds $thr" > "$L/base_$pid"
        return 0
    fi
    read -r fds0 thr0 < "$L/base_$pid"
    if [ "$fds" -gt $(( fds0 + 64 )) ]; then
        echo "$fds open files, $fds0 after its first task"
    elif [ "$thr" -gt $(( thr0 + 256 )) ]; then
        echo "$thr threads, $thr0 after its first task"
    fi
}

has_work() {
    # 0 if <own> has unclaimed tasks or (steal) the call has unstarted tasks
    candidates "$1" "$2" | grep -q '^[oa] '
}

unfinished() {
    # the tasks a process claimed and didn't finish (status file $1)
    awk '$1 == "c" {c[$2] = 1} $1 == "d" || $1 == "s" {delete c[$2]} END {for (u in c) print u}' "$1"
}

slot() {
    # one slot: <own> <gpu> <steal> <retries> <local folder> <slot> <slots>
    local own=$1 gpu=$2 steal=$3 retries=$4 L=$5 s=$6 P=$7 k=0 est=0 f u line st kind nown nwait nsteal later more
    [ "$gpu" = 1 ] && export CUDA_VISIBLE_DEVICES=$(( s - 1 ))
    while true; do
        k=$((k + 1))
        st="$L/s${s}_$k"
        : > "$st"
        if [ "$steal" = 3 ] && [ "$k" = 1 ]; then
            # (steal 3, a single job) its own task first, as before: claimed
            # without listing the pool and run by its task line in a process
            # of its own; a job whose task another process has ends at once.
            # Then the pool's unstarted tasks
            if ! ln -s "$(owner "$s")" "$C/t$own" 2> /dev/null; then
                echo "Task runner: no task left to run for $own"
                break
            fi
            f=$own
            nown=0
            nwait=0
            nsteal=1
        else
            claim_next "$own" "$steal" "$s" "$est" "$st" - "$L/next_$s"
            f=$(cat "$L/next_$s")
            [ -z "$f" ] && break
            read -r kind nown nwait nsteal < "$L/next_$s.more"
            [ "$kind" = o ] && nown=$(( nown - 1 ))
            [ "$kind" = a ] && nsteal=$(( nsteal - 1 ))
        fi
        echo "c $f" >> "$st"
        line=$(task_line "$f")
        # the process runs more tasks (a runner) only if the slot has more
        # to do: own tasks beyond one per slot that hasn't started, tasks of
        # batches that wait in the queue, or (local runs) any; else the task
        # line runs as before
        later=0
        [ "$k" = 1 ] && later=$(( P - s ))
        more=0
        if [ "$LOCALRUN" = 1 ] || [ "$nown" -gt "$later" ] || { [ "$steal" != 0 ] && [ "$nwait" -gt 0 ]; }; then
            echo "Task runner $(owner "$s"): task $f in a new process, which runs more tasks"
            ( export PETAKIT5D_RUNNER="$D|$own|$steal|$s|$st|$L|$f"; eval "$line" )
            more=1
        else
            local t0=$SECONDS
            ( eval "$line" ) && echo "d $f" >> "$st"
            [ $(( SECONDS - t0 )) -gt "$est" ] && est=$(( SECONDS - t0 ))
        fi
        # tasks the process didn't finish (an error, a crash): once more in a
        # process of their own (as before), then they count as failed and
        # the slot stops taking other batches' tasks
        for u in $(unfinished "$st"); do
            if [ "$retries" -gt 1 ]; then
                echo "Task runner $(owner "$s"): task $u did not finish, run it in a process of its own"
                line=$(task_line "$u")
                if ( eval "$line" ); then
                    echo "d $u" >> "$st"
                    continue
                fi
            fi
            ln -s "$(owner "$s")" "$C/f$u" 2> /dev/null
            steal=0
        done
        if [ "$LOCALRUN" = 1 ]; then
            for u in $(awk '$1 == "d" || $1 == "c" {print $2}' "$st"); do
                local_tmp_done "$u"
            done
            # a task the runner gave back (its claim removed, for a new
            # process) loses the .tmp of its claim, so it runs now
            for u in $(awk '$1 == "s" {print $2}' "$st"); do
                [ -L "$C/t$u" ] || local_tmp_release "$u"
            done
        fi
        # the runner found no task left, or (a task line run as before) there
        # was none left to claim
        if [ "$more" = 1 ]; then
            est=$(awk -v e="$est" '$1 == "t" && $2 > e {e = $2} END {printf "%d", e + 0.999}' "$st")
            grep -q '^x none' "$st" && break
        elif [ "$nown" -le 0 ] && { [ "$steal" = 0 ] || [ "$nsteal" -le 0 ]; }; then
            break
        fi
    done
}

job() {
    # a compiled job: <own> <slots> <gpu> <steal> <retries>
    local own=$1 P=$2 gpu=$3 steal=$4 retries=$5 L n
    # no more slots than tasks to take (GNU parallel started one process
    # per task line): the batch's unclaimed tasks and (steal) the others'
    # (steal 3: one, the single job's own task, which the slot claims first)
    if [ "$steal" = 3 ]; then
        n=1
    else
        n=$(candidates "$own" "$steal" | awk '$1 == "n" {print $2 + $4; exit}')
    fi
    if [ -z "$n" ] || [ "$n" -eq 0 ]; then
        echo "Task runner: no task left to run for $own"
        return 0
    fi
    [ "$n" -lt "$P" ] && P=$n
    L=$(mktemp -d "${TMPDIR:-/tmp}/pk5d_runner_XXXXXX")
    if [ "$P" -gt 1 ]; then
        parallel --ungroup --jobs "$P" --delay 0.2 bash "$SELF" slot "$own" "$gpu" "$steal" "$retries" "$L" {%} "$P" ::: $(seq "$P")
    else
        slot "$own" "$gpu" "$steal" "$retries" "$L" 1 1
    fi
    rm -rf "$L"
}

matlab_job() {
    # a MATLAB job: <own> <steal>
    local own=$1 steal=$2 L k=0 est=0 f u st first
    first=0
    [ "$steal" = 3 ] && first=1
    if ! has_work "$own" "$(( first ? 0 : steal ))"; then
        echo "Task runner: no task left to run for $own"
        return 0
    fi
    L=$(mktemp -d "${TMPDIR:-/tmp}/pk5d_runner_XXXXXX")
    while true; do
        k=$((k + 1))
        st="$L/s1_$k"
        : > "$st"
        # (steal 3) the job's first task is one of its own
        [ "$k" -gt 1 ] && first=0
        claim_next "$own" "$steal" 1 "$est" "$st" "$L/dump_1" "$L/next_1" "$first"
        f=$(cat "$L/next_1")
        [ -z "$f" ] && break
        echo "c $f" >> "$st"
        eval "$MLAUNCH \"setup([]);t0_=tic;task_runner('matlab','$D','$own',$steal,'$st','$L',$f);toc(t0_)\""
        # failed tasks: the framework resubmits them with more CPUs; the job
        # stops taking other batches' tasks
        for u in $(unfinished "$st"); do
            ln -s "$(owner 1)" "$C/f$u" 2> /dev/null
            steal=0
        done
        est=$(awk -v e="$est" '$1 == "t" && $2 > e {e = $2} END {printf "%d", e + 0.999}' "$st")
        grep -q '^x none' "$st" && break
    done
    rm -rf "$L"
}

leave_pool() {
    # (pools) a single job that ends removes the mark of its task $1 in
    # live/; the last to leave of the pool's jobs and its master (which
    # leaves at the end of its call) removes the pool: rmdir of live/
    # succeeds only if it is empty, which the file server decides. The job
    # ends with status $2
    [ -d "$D/live" ] || return 0
    rm -f "$D/live/t$1"
    case "$D" in */task_runner_single_*) ;; *) return 0 ;; esac
    if rmdir "$D/live" 2> /dev/null; then
        # (exec: no shell has the script open then, so all of it goes)
        exec sh -c 'cd / && rm -rf "$1"; exit "$2"' sh "$D" "$2"
    fi
    return 0
}

case "$1" in
    job)
        job "$2" "$3" "$4" "$5" "$6"
        rc=$?
        [ "$5" = 3 ] && leave_pool "$2" "$rc"
        exit "$rc"
        ;;
    slot)
        slot "$2" "$3" "$4" "$5" "$6" "$7" "$8"
        ;;
    matlab)
        matlab_job "$2" "$3"
        rc=$?
        [ "$3" = 3 ] && leave_pool "$2" "$rc"
        exit "$rc"
        ;;
    claim)
        why=$(recycle_reason "$PK_PID" "$7")
        if [ -n "$why" ]; then
            echo "recycle $why" > "$8"
        else
            claim_next "$2" "$3" "$4" "$5" "$6" "$7/dump_$4" "$8"
        fi
        ;;
    has_work)
        has_work "$2" "$3"
        ;;
esac
