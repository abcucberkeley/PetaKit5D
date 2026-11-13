# Diagnostic for the CI job exit-probe (started by suite_exit_probe.m): follows
# the MATLAB process Id after its tests ended (TestsEnd: Unix time in seconds)
# and appends to File a line while it is alive (its threads and CPU seconds,
# every 10 s for a minute, then every minute) and one when it is gone.
param([int]$Id, [double]$TestsEnd, [string]$File)

$next = 0
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
    Start-Sleep -Milliseconds 250
}
