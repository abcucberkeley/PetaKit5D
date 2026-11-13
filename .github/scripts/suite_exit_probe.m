function [] = suite_exit_probe(part)
% diagnostic of MATLAB's exit on Windows after the PetaKit5D tests (CI job
% exit-probe; not a test). It runs the tests of one folder of tests/ (part;
% 'all' for all of them) as the test job does, but without pktest.finishCI,
% writes probe-results/suite_<part>.txt (the counts, the mex files in memory
% and the process's threads when the tests ended) and returns, so that matlab
% -batch exits. On Windows it then starts watch_matlab_exit.ps1, a process of
% its own, which adds to the file the seconds from the end of the tests until
% this MATLAB is gone, and while it is alive its threads and CPU time (every 10
% s for a minute, then every minute). The CI step after it adds the step's
% outcome (success: MATLAB exited by itself; failure: the step's timeout or an
% error) and the MATLAB processes left.

root = pwd;
outDir = fullfile(root, 'probe-results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
resultFn = fullfile(outDir, ['suite_', part, '.txt']);
note = @(varargin) writelines(string(sprintf(varargin{:})), resultFn, 'WriteMode', 'append');
% (times in UTC, as the watcher's and the CI step's)
utc = @() char(datetime('now', 'TimeZone', 'UTC', 'Format', 'HH:mm:ss.SS'));
note('suite %s start %s, MATLAB %s, maxNumCompThreads %d, process %d', part, utc(), version, ...
    maxNumCompThreads, feature('getpid'));

addpath(fullfile(root, 'tests'));
if strcmp(part, 'all')
    r = run_tests();
else
    r = run_tests('folder', part);
end
tEnd = posixtime(datetime('now', 'TimeZone', 'UTC'));
[~, mexFiles] = inmem('-completenames');
[~, mexNames] = cellfun(@fileparts, mexFiles, 'unif', 0);
note('tests ended %s: %d tests, %d passed, %d failed, %d incomplete; %d mex files in memory (%s); process threads %d', ...
    utc(), numel(r), nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete] & ~[r.Failed]), ...
    numel(mexNames), strjoin(mexNames, ' '), process_threads());
if ispc
    % (Start-Process: neither in this MATLAB's process tree nor holding its
    % output pipes)
    watcher = fullfile(root, '.github', 'scripts', 'watch_matlab_exit.ps1');
    system(sprintf(['powershell -NoProfile -NonInteractive -Command "Start-Process powershell -WindowStyle Hidden ', ...
        '-ArgumentList ''-NoProfile'',''-ExecutionPolicy'',''Bypass'',''-File'',''%s'',%d,%.3f,''%s''"'], ...
        watcher, feature('getpid'), tEnd, resultFn));
end
note('returns %s', utc());

end
