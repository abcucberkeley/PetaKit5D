function [summary] = isolation_report(outDir, varargin)
% compare the outputs of each harness run in outDir/results/<run> with those
% of the fresh processes (outDir/results/fresh), and summarize what the tasks
% changed in their processes (outDir/records/<run>.mat). Writes
% outDir/summary.txt and outDir/summary.mat.
%
% Parameters (as 'specifier'-value pairs):
%   baseline : char (default: 'fresh'). Run the others are compared with.
%     ignore : cell array (default: {'(^|/)logs/', 'job_logs'}). Relative paths to skip (regular expressions).


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('outDir', @ischar);
ip.addParameter('baseline', 'fresh', @ischar);
ip.addParameter('ignore', {'(^|/)logs/', 'job_logs'}, @iscell);
ip.parse(outDir, varargin{:});

pr = ip.Results;
resDir = [outDir, '/results/'];
d = dir(resDir);
runs = setdiff({d([d.isdir]).name}, {'.', '..', pr.baseline});
lines = {};
summary = struct('run', {}, 'nEqual', {}, 'nDiffer', {}, 'nOnlyRun', {}, 'nOnlyBaseline', {}, 'diffs', {}, ...
    'failedTasks', {}, 'changes', {});
for i = 1 : numel(runs)
    res = isolation_compare([resDir, pr.baseline], [resDir, runs{i}], 'ignore', pr.ignore);
    s.run = runs{i};
    s.nEqual = nnz(strcmp({res.status}, 'equal'));
    s.nDiffer = nnz(strcmp({res.status}, 'differ'));
    s.nOnlyRun = nnz(strcmp({res.status}, 'onlyB'));
    s.nOnlyBaseline = nnz(strcmp({res.status}, 'onlyA'));
    s.diffs = res(~strcmp({res.status}, 'equal'));
    [s.failedTasks, s.changes] = read_records(outDir, runs{i});
    summary(end + 1) = s; %#ok<AGROW>

    lines{end + 1} = sprintf('== %s vs %s: %d equal, %d differ, %d only in %s, %d only in %s; %d failed tasks', ...
        runs{i}, pr.baseline, s.nEqual, s.nDiffer, s.nOnlyRun, runs{i}, s.nOnlyBaseline, pr.baseline, numel(s.failedTasks)); %#ok<AGROW>
    for j = 1 : numel(s.diffs)
        lines{end + 1} = sprintf('   %s %s (%s) %s', s.diffs(j).status, s.diffs(j).path, s.diffs(j).kind, s.diffs(j).detail); %#ok<AGROW>
    end
    for j = 1 : numel(s.failedTasks)
        lines{end + 1} = sprintf('   failed task: %s', s.failedTasks{j}); %#ok<AGROW>
    end
end
[failed, changes] = read_records(outDir, pr.baseline);
lines{end + 1} = sprintf('== %s: %d failed tasks', pr.baseline, numel(failed));
for j = 1 : numel(failed)
    lines{end + 1} = sprintf('   failed task: %s', failed{j}); %#ok<AGROW>
end
lines{end + 1} = '== what the tasks changed in their processes (fresh runs):';
lines = [lines, changes];

fid = fopen([outDir, '/summary.txt'], 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
fprintf('%s\n', lines{:});
save([outDir, '/summary.mat'], 'summary');

end


function [failed, changes] = read_records(outDir, run)
% failed tasks of a run and, per task, what it changed in its process

failed = {};
changes = {};
% one file per run, or per task for the fresh processes (<run>_<task>.mat)
fns = dir(sprintf('%s/records/%s.mat', outDir, run));
if isempty(fns)
    fns = dir(sprintf('%s/records/%s_*.mat', outDir, run));
end
for i = 1 : numel(fns)
    a = load([fns(i).folder, '/', fns(i).name], 'records');
    for j = 1 : numel(a.records)
        r = a.records(j);
        if ~r.ok
            failed{end + 1} = sprintf('%s: %s', r.id, first_line(r.message)); %#ok<AGROW>
        end
        c = r.changes(~startsWith(r.changes, '(info)'));
        if ~isempty(c)
            changes{end + 1} = sprintf('   %s: %s', r.id, strjoin(c, '; ')); %#ok<AGROW>
        end
    end
end

end


function [s] = first_line(s)

s = strtrim(s);
k = find(s == newline, 1);
if ~isempty(k)
    s = s(1 : k - 1);
end

end
