function [] = suite_exit_probe(part, names, first)
% diagnostic of MATLAB's exit on Windows after the PetaKit5D tests (CI job
% exit-probe; not a test). It runs the tests of one folder of tests/ or one
% test class (part; 'all' for all of them), or the tests of one class that
% names gives (run_tests' name patterns; part then names the results), or the
% first of these (as many as first, in the order that run_tests runs them), as
% the test job does, but without pktest.finishCI, writes
% probe-results/suite_<part>.txt (the counts, the mex files in memory, the
% process's threads, what could hold MATLAB open (open_handles) and, on
% Windows, the processes it started that still run (child_processes.ps1), all
% when the tests ended) and returns, so that matlab -batch exits. On Windows it
% then starts watch_matlab_exit.ps1, a process of its own, which adds to the
% file the seconds from the end of the tests until this MATLAB is gone, while
% it is alive its threads and CPU time (every 10 s for a minute, then every
% minute), and 30 s after the tests ended the stacks of its threads
% (thread_stacks.ps1) and the processes it started that still run. The CI
% step after it adds the step's outcome (success: MATLAB exited by itself;
% failure: the step's timeout or an error), the MATLAB processes left and the
% processes that they started and that still run.

if nargin < 2
    names = {[part, '/*']};
end
names = cellstr(names);
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
elseif exist(fullfile(root, 'tests', part), 'dir')
    r = run_tests('folder', part);
else
    % tests of one class: its folder's suite, of these tests only
    d = dir(fullfile(root, 'tests', '*', [extractBefore([names{1}, '/'], '/'), '.m']));
    [~, folder] = fileparts(d(1).folder);
    if nargin > 2
        % (the suite as run_tests makes it: the folder's, without the known
        % issues, of these names, in its order)
        s = matlab.unittest.TestSuite.fromFolder(fullfile(root, 'tests', folder), 'IncludingSubfolders', true);
        s = s.selectIf(~matlab.unittest.selectors.HasTag('KnownIssue'));
        expr = strjoin(cellfun(@(x) ['^', regexptranslate('wildcard', x), '$'], names, 'unif', 0), '|');
        s = s(~cellfun(@isempty, regexp({s.Name}, expr, 'once')));
        names = regexprep({s(1 : min(first, end)).Name}, '[()]', '?');
    end
    if nargin > 1
        note('tests %s', strjoin(names, ' '));
    end
    r = run_tests('folder', folder, 'name', names);
end
tEnd = posixtime(datetime('now', 'TimeZone', 'UTC'));
[~, mexFiles] = inmem('-completenames');
[~, mexNames] = cellfun(@fileparts, mexFiles, 'unif', 0);
note('tests ended %s: %d tests, %d passed, %d failed, %d incomplete; %d mex files in memory (%s); process threads %d', ...
    utc(), numel(r), nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete] & ~[r.Failed]), ...
    numel(mexNames), strjoin(mexNames, ' '), process_threads());
note('open when the tests ended: %s', open_handles());
if ispc
    scripts = fullfile(root, '.github', 'scripts');
    system(sprintf('powershell -NoProfile -NonInteractive -ExecutionPolicy Bypass -File "%s" %d "%s" "when the tests ended"', ...
        fullfile(scripts, 'child_processes.ps1'), feature('getpid'), resultFn));
    % (Start-Process: neither in this MATLAB's process tree nor holding its
    % output pipes)
    watcher = fullfile(scripts, 'watch_matlab_exit.ps1');
    system(sprintf(['powershell -NoProfile -NonInteractive -Command "Start-Process powershell -WindowStyle Hidden ', ...
        '-ArgumentList ''-NoProfile'',''-ExecutionPolicy'',''Bypass'',''-File'',''%s'',%d,%.3f,''%s''"'], ...
        watcher, feature('getpid'), tEnd, resultFn));
end
note('returns %s', utc());

end
