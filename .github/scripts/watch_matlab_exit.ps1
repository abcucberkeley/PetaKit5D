# Diagnostic for the CI job exit-probe (started by suite_exit_probe.m): follows
# the MATLAB process Id after its tests ended (TestsEnd: Unix time in seconds)
# and appends to File a line while it is alive (its threads and CPU seconds,
# every 10 s for a minute, then every minute) and one when it is gone, and if it
# is alive 30 s after its tests ended (MATLAB that exits by itself is gone in
# seconds), the stacks of its threads (thread_stacks.ps1, up to 90 s, which
# delays the lines after them) and the processes it started that still run
# (child_processes.ps1).
param([int]$Id, [double]$TestsEnd, [string]$File)

$next = 0
$listed = $false
while ($true) {
    $p = Get-Process -Id $Id -ErrorAction SilentlyContinue
    $s = [DateTimeOffset]::UtcNow.ToUnixTimeMilliseconds() / 1000.0 - $TestsEnd
    if (-not $p) {
        Add-Content -Path $File -Value ("MATLAB gone {0:F1} s after the tests ended ({1:HH:mm:ss})" -f $s, (Get-Date).ToUniversalTime())
        break
    }
    if ($s -ge $next) {
        Add-Content -Path $File -Value ("MATLAB alive {0:F1} s after the tests ended: threads {1}, CPU {2:F1} s" -f $s, $p.Threads.Count, $p.CPU)
        if ($next -lt 60) {
            $next = $next + 10
        } else {
            $next = $next + 60
        }
    }
    if (-not $listed -and $s -ge 30) {
        $listed = $true
        try {
            & (Join-Path $PSScriptRoot 'thread_stacks.ps1') $Id $File
        } catch {
            Add-Content -Path $File -Value ('the stacks failed: {0}' -f $_)
        }
        try {
            & (Join-Path $PSScriptRoot 'child_processes.ps1') $Id $File '30 s after the tests ended'
        } catch {
            Add-Content -Path $File -Value ('listing the processes failed: {0}' -f $_)
        }
    }
    Start-Sleep -Milliseconds 250
}
