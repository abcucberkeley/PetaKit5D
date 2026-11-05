function [st] = isolation_state()
% snapshot of the process-wide state that a task can change and a new process
% would start with: computational threads, warning states, random numbers,
% current folder, path, environment, global variables, open files and file
% descriptors, figures, timers, parallel pools, loaded mex files, memory, and
% the GPU this process uses. It doesn't change any of it (no GPU or graphics
% initialization: the GPU is read with nvidia-smi).
%
% Output: struct; compare two snapshots with isolation_state_diff.


st = struct();
st.time = datetime('now');
% the figures first: the first use of the graphics root turns some
% deprecation warnings off, which would look like a change of the next task
st.figures = NaN;
try
    st.figures = numel(findall(groot, 'Type', 'figure'));
catch
end
st.threads = maxNumCompThreads;
st.warnings = warning;
[st.lastwarnMsg, st.lastwarnId] = lastwarn;
st.rng = rng;
st.pwd = pwd;
st.path = path;
st.env = env_list();
st.globals = sort(who('global'));
st.baseVars = sort(evalin('base', 'who'));
st.files = open_files();
st.fds = NaN;
st.rssGB = NaN;
st.peakGB = NaN;
st.osThreads = NaN;
if isunix && ~ismac
    st.fds = file_descriptors();
    [st.rssGB, st.peakGB, st.osThreads] = proc_status();
end
st.timers = NaN;
try
    st.timers = numel(timerfindall);
catch
end
st.pool = false;
try
    st.pool = ~isempty(gcp('nocreate'));
catch
end
[~, mexs] = inmem('-completenames');
st.mex = sort(mexs(:));
st.gpuMemMB = gpu_memory();

end


function [env] = env_list()
% sorted {name, value} of the environment

env = cell(0, 2);
try
    e = getenv();
    keys = sort(e.keys);
    env = [cellstr(keys), cellstr(e(keys))];
catch
    if isunix
        [~, out] = system('env -0');
        lines = strsplit(out, char(0));
        lines = lines(~cellfun(@isempty, lines));
        env = cell(numel(lines), 2);
        for i = 1 : numel(lines)
            k = find(lines{i} == '=', 1);
            env(i, :) = {lines{i}(1 : k - 1), lines{i}(k + 1 : end)};
        end
        [~, inds] = sort(env(:, 1));
        env = env(inds, :);
    end
end

end


function [files] = open_files()
% {fid, name} of the files open with fopen

fids = isolation_open_fids();
files = cell(numel(fids), 2);
for i = 1 : numel(fids)
    files(i, :) = {fids(i), fopen(fids(i))};
end

end


function [rssGB, peakGB, nThreads] = proc_status()

rssGB = NaN;
peakGB = NaN;
nThreads = NaN;
try
    txt = fileread('/proc/self/status');
    rssGB = str2double(regexp(txt, 'VmRSS:\s*(\d+)', 'tokens', 'once')) / 2^20;
    peakGB = str2double(regexp(txt, 'VmHWM:\s*(\d+)', 'tokens', 'once')) / 2^20;
    nThreads = str2double(regexp(txt, 'Threads:\s*(\d+)', 'tokens', 'once'));
catch
end

end


function [memMB] = gpu_memory()
% GPU memory used by this process (MB, nvidia-smi; NaN without a GPU, 0
% without a CUDA context)

persistent hasSmi
memMB = NaN;
if ~isunix || ismac
    return;
end
if isempty(hasSmi)
    [status, ~] = system('nvidia-smi -L > /dev/null 2>&1');
    hasSmi = status == 0;
end
if ~hasSmi
    return;
end
[status, out] = system('nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader,nounits');
if status ~= 0
    return;
end
memMB = 0;
pid = feature('getpid');
lines = splitlines(strtrim(out));
for i = 1 : numel(lines)
    v = str2double(strsplit(lines{i}, ','));
    if numel(v) == 2 && v(1) == pid
        memMB = memMB + v(2);
    end
end

end


function [n] = file_descriptors()
% the open file descriptors of this process that point to files outside
% MATLAB, as a task can leave them open (a deleted file too). MATLAB's own
% descriptors (sockets, pipes, event and device descriptors, its files in
% matlabroot and prefdir) come and go in the background: their count changed
% by one between two looks on R2026b (IsolationHarnessTest on GitHub)

n = NaN;
[status, out] = system(sprintf('ls -l /proc/%d/fd 2>/dev/null', feature('getpid')));
if status ~= 0
    return;
end
targets = regexp(out, '-> ([^\n]+)', 'tokens');
targets = cellfun(@(x) strtrim(x{1}), targets, 'unif', 0);
n = nnz(startsWith(targets, '/') & ~startsWith(targets, {'/dev/', '/proc/', '/sys/', matlabroot, prefdir}));

end
