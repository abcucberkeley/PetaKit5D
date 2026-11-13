# Diagnostic for the CI job exit-probe: appends to File the processes that the
# process Id started and that still run, and those they started, recursively
# (Id 0: of each MATLAB process running now), one line each with its id, name,
# start time (UTC) and command line, after a line with Label; and which of the
# processes that File listed before still run. (A child is a process whose
# ParentProcessId is its parent's and that started after its parent, as ids
# are reused.)
param([int]$Id, [string]$File, [string]$Label)

$all = @(Get-CimInstance Win32_Process)
$key = { param($p) '{0} {1} {2}' -f $p.ProcessId, $p.Name, $p.CreationDate.ToUniversalTime().ToString('HH:mm:ss.ff') }
$now = @{}
foreach ($p in $all) {
    if ($p.CreationDate) {
        $now[(& $key $p)] = $true
    }
}
# (not this listing's process or the one that started it)
$skip = @($PID) + @($all | Where-Object { $_.ProcessId -eq $PID } | ForEach-Object { $_.ParentProcessId })

function Get-Descendants($parent, [int]$depth) {
    $kids = @($all | Where-Object { $_.ParentProcessId -eq $parent.ProcessId -and $_.CreationDate -and $_.CreationDate -ge $parent.CreationDate })
    foreach ($c in @($kids | Where-Object { $skip -notcontains $_.ProcessId })) {
        '{0}child {1}: {2}' -f ('  ' * $depth), (& $key $c), $c.CommandLine
        Get-Descendants $c ($depth + 1)
    }
}

$before = @()
if (Test-Path $File) {
    $before = @(Select-String -Path $File -Pattern '^\s*child (\d+ \S+ \S+): ' | ForEach-Object { $_.Matches[0].Groups[1].Value } |
        Sort-Object -Unique)
}
if ($Id -gt 0) {
    $roots = @($all | Where-Object { $_.ProcessId -eq $Id })
} else {
    $roots = @($all | Where-Object { $_.Name -eq 'MATLAB.exe' })
}
$out = @()
if ($before.Count -gt 0) {
    $alive = @($before | Where-Object { $now.ContainsKey($_) })
    $out += 'of the {0} processes listed before, running {1}: {2}' -f $before.Count, $Label, $(if ($alive.Count) { $alive -join ', ' } else { 'none' })
}
foreach ($r in $roots) {
    $lines = @(Get-Descendants $r 1)
    $out += 'processes started by MATLAB {0} (or by those) running {1}: {2}' -f $r.ProcessId, $Label, $lines.Count
    $out += $lines
}
if ($roots.Count -eq 0) {
    $out += "no MATLAB process running $Label"
}
Add-Content -Path $File -Value $out
