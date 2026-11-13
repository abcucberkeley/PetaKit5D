function [] = start_exit_watch(root, resultFn, tEnd)
% for the CI job exit-probe (suite_exit_probe, mex_exit_probe): on Windows,
% snapshot the processes this MATLAB started (child_processes.ps1) and start
% watch_matlab_exit.ps1 as a process of its own, which records whether this
% MATLAB exits by itself after the tests ended (tEnd, Unix seconds) and, if it
% is still alive 30 s later, its threads' stacks. Call it last, just before
% returning, so matlab -batch exits. Does nothing off Windows.

if ~ispc
    return;
end
scripts = fullfile(root, '.github', 'scripts');
system(sprintf('powershell -NoProfile -NonInteractive -ExecutionPolicy Bypass -File "%s" %d "%s" "when the tests ended"', ...
    fullfile(scripts, 'child_processes.ps1'), feature('getpid'), resultFn));
% (Start-Process: neither in this MATLAB's process tree nor holding its output
% pipes)
watcher = fullfile(scripts, 'watch_matlab_exit.ps1');
system(sprintf(['powershell -NoProfile -NonInteractive -Command "Start-Process powershell -WindowStyle Hidden ', ...
    '-ArgumentList ''-NoProfile'',''-ExecutionPolicy'',''Bypass'',''-File'',''%s'',%d,%.3f,''%s''"'], ...
    watcher, feature('getpid'), tEnd, resultFn));

end
