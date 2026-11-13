#!/bin/bash
# For the CI job exit-probe: after a probe step whose result file is
# probe-results/$1.txt (outcome $2: success, MATLAB exited by itself; failure,
# the step's timeout or an error), adds the time, the MATLAB processes still
# running, the processes that they started and that still run, and which of the
# processes listed before still run to it. $1 is the result basename
# (suite_<part> for suite_exit_probe, the mode for mex_exit_probe).
f="probe-results/$1.txt"
mkdir -p probe-results
left=$(powershell -NoProfile -NonInteractive -Command \
    "(Get-Process MATLAB -ErrorAction SilentlyContinue | ForEach-Object { '{0} (threads {1}, CPU {2:F1} s)' -f \$_.Id, \$_.Threads.Count, \$_.CPU }) -join ', '")
{
    echo "step ended $(date -u +%H:%M:%S.%2N), outcome $2"
    echo "MATLAB processes now: ${left:-none}"
} >> "$f"
powershell -NoProfile -NonInteractive -ExecutionPolicy Bypass -File .github/scripts/child_processes.ps1 0 "$f" "when the step ended" ||
    echo "listing the processes failed" >> "$f"
