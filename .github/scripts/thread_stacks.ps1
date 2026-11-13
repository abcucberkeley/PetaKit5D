# Diagnostic for the CI job exit-probe, called by watch_matlab_exit.ps1 when
# MATLAB still runs 30 s after its tests ended: appends to File the stack of
# every thread of the process Id, each stack once with the threads that share
# it (cdb's !uniqstack; cdb attaches non-invasively, suspending the threads
# while it reads them; module!function from Microsoft's public symbols, or
# else from the modules' exports), or without cdb the thread list of a
# minidump that procdump writes (dump_threads.py), at most MaxLines lines and
# 90 s. (PROBE_CDB or PROBE_PROCDUMP: the debugger that the job found.)
param([int]$Id, [string]$File, [int]$MaxLines = 1500)

$out = Join-Path $env:RUNNER_TEMP "stacks_$Id.txt"
$t = [Diagnostics.Stopwatch]::StartNew()
if ($env:PROBE_CDB) {
    $sym = 'srv*{0}*https://msdl.microsoft.com/download/symbols' -f (Join-Path $env:RUNNER_TEMP 'symbols')
    $cmd = '-pv -p {0} -y "{1}" -c ".echo STACKS_BEGIN; !uniqstack; .echo STACKS_END; qd"' -f $Id, $sym
    $p = Start-Process -FilePath $env:PROBE_CDB -ArgumentList $cmd -RedirectStandardOutput $out -NoNewWindow -PassThru
    $what = 'cdb !uniqstack'
} elseif ($env:PROBE_PROCDUMP) {
    $dump = Join-Path $env:RUNNER_TEMP "matlab_$Id.dmp"
    & $env:PROBE_PROCDUMP -accepteula $Id $dump | Out-Null
    $cmd = '-I "{0}" "{1}"' -f (Join-Path $PSScriptRoot 'dump_threads.py'), $dump
    $p = Start-Process -FilePath python -ArgumentList $cmd -RedirectStandardOutput $out -NoNewWindow -PassThru
    $what = 'procdump and dump_threads.py'
} else {
    Add-Content -Path $File -Value 'no stacks: the job found no debugger'
    return
}
$state = ''
if (-not $p.WaitForExit(90000)) {
    $p.Kill()
    $state = ', stopped after 90 s'
}
$lines = @()
if (Test-Path $out) {
    $lines = @(Get-Content $out)
}
# (of cdb's output, the commands' between the markers, without the lines of
# its symbol loading; all of it if a marker is missing)
$b = -1
$e = $lines.Count
for ($k = 0; $k -lt $lines.Count; $k++) {
    if ($b -lt 0 -and $lines[$k].Trim() -eq 'STACKS_BEGIN') {
        $b = $k
    } elseif ($b -ge 0 -and $lines[$k].Trim() -eq 'STACKS_END') {
        $e = $k
        break
    }
}
if ($b -ge 0) {
    $lines = @($lines | Select-Object -Skip ($b + 1) -First ($e - $b - 1) | Where-Object { $_ -notmatch '^\*\*\* ' })
}
$text = @('== stacks of MATLAB {0} 30 s after the tests ended ({1}, {2:F0} s{3}): {4} lines' -f $Id, $what, $t.Elapsed.TotalSeconds, $state, $lines.Count)
$text += @($lines | Select-Object -First $MaxLines)
if ($lines.Count -gt $MaxLines) {
    $text += "== cut after $MaxLines lines"
}
$text += '== end of the stacks'
Add-Content -Path $File -Value $text
