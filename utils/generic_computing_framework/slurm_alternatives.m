function [varargout] = slurm_alternatives(action, varargin)
% places on a Slurm cluster for the jobs of a job framework call besides its
% main settings (the SlurmAlternatives of a cluster config: entries with
% their own SlurmParam, memPerCPU, maxCPUNum, minCPUNum, gpuHeadroom and
% cpuJobs; keys an entry leaves out are the main settings', except memPerCPU:
% without it, the entry's --mem-per-cpu, else the memory per CPU of its
% nodes, the smallest RealMemory / CPUs of them, and cpuJobs: false). Each job
% is submitted to the first place (the main one first) where Slurm can start
% it now, and a job that waits in its place moves to another place where it
% can start now. GPU nodes are for jobs that need a GPU: a CPU-only job (one
% whose request has no GPUs) never goes to an alternative of GPU nodes only,
% unless the entry says cpuJobs: true (the main place is the user's choice,
% and used as before).
%
% A job's shape (CPUs and memory) at a place: a place with as much memory
% per CPU (memPerCPU) as the job's shape has gets the job a config with only
% that place submits (its CPUs from memAllocate and that place's memPerCPU,
% minCPUNum and maxCPUNum, its memory from that place's SlurmParam), so a
% node fills in CPUs and memory together; at a place with less memory per CPU
% the job keeps its CPUs and memory (more CPUs can mean more workers, and so
% more memory than estimated). A job is first shaped for the main place.
%
% Whether a place can start a job now comes from Slurm: the nodes (scontrol
% show node: free CPUs, memory and GPUs, without nodes the backfill
% scheduler planned for a pending job), the user's pending jobs (squeue: they
% start first, a QOS or association limit they wait for closes the place),
% and sbatch --test-only for the job's options, read at most every 20 s (a
% view; the jobs placed since are taken from it, so a burst of jobs doesn't
% go to the same free nodes). A CPU-only job on a node with idle GPUs leaves
% each of them its share of the node (gpuHeadroom: [CPUs, GB] per idle GPU,
% by default the node's CPUs and memory divided by its GPUs; 0 for none),
% doesn't go to a partition where GPU jobs of the user are pending (also held
% ones, or ones that wait for a dependency or a start time), and doesn't wait
% in the queue of a place of GPU nodes only: if no place has room, it goes to
% (or moves to) the first place that isn't one, and one waiting at a place of
% GPU nodes where GPU jobs of the user are pending moves to such a place.
%
% Usage:
%   places = slurm_alternatives('places', main, alternatives)
%       main: a struct with SlurmParam, SlurmConstraint, memPerCPU, maxCPUNum,
%       minCPUNum and gpuHeadroom; alternatives: a struct array, a cell array
%       of structs or empty. A struct array with the main place first (one
%       place without alternatives)
%   timeStrs = slurm_alternatives('time_strs', places, jobTimeLimit, fmt, mainStr)
%       the time limit options of the places as the framework writes them
%       (fmt, e.g. ' -t %d:%02d:00 '; none for a place whose SlurmParam has a
%       time limit; the main place's is mainStr)
%   [k, cpus, SlurmParam] = slurm_alternatives('choose', places, ownCpus, timeStrs)
%       the place for a job whose own CPUs at place k (as a config with only
%       place k gives them) are ownCpus(k), its CPUs there and the SlurmParam
%       of its shape (place k's, with the memory option of the shape)
%   slurm_alternatives('submitted', jobId, arrayId, places, ownCpus, k, timeStrs)
%       records a job submitted to place k, so that it can move later
%   slurm_alternatives('move', jobIds, arrayIds)
%       recorded jobs that the caller found pending, and that have waited in
%       their place for 2 minutes (since their submission or last move) and
%       weren't looked at in the last minute, move to the first other place
%       that has room for them now (scontrol update of partition, QOS,
%       account, CPUs, memory, features, node lists and job name), as many as
%       there is room for; between looks no Slurm command runs
%   slurm_alternatives('forget', jobIds, arrayIds)
%       drops the records of jobs that ended
%   slurm_alternatives('set', name, value, ...)
%       timing settings (s): viewAge (20), minWait (120), lookInterval (60),
%       startTolerance (60)


persistent S view reg disabled
if isempty(S)
    S = struct('viewAge', 20, 'minWait', 120, 'lookInterval', 60, 'startTolerance', 60);
    reg = containers.Map('KeyType', 'char', 'ValueType', 'any');
    disabled = {};
end

switch action
    case 'places'
        [places, view] = make_places(view, S, varargin{:});
        varargout{1} = places;
    case 'time_strs'
        varargout{1} = time_strings(varargin{:});
    case 'choose'
        [places, ownCpus, timeStrs] = varargin{1 : 3};
        view = current_view(view, S);
        [k, shape, view, disabled] = choose_place(view, places, ownCpus, timeStrs, S, disabled);
        if k > 1
            fprintf('The job goes to %s (%d CPUs): the main place has no room for it now.\n', ...
                place_name(places(k)), shape.cpus);
        end
        varargout = {k, shape.cpus, shape.param};
    case 'submitted'
        [jobId, arrayId, places, ownCpus, k, timeStrs] = varargin{1 : 6};
        t = now_s();
        reg(job_key(jobId, arrayId)) = struct('places', places, 'cpus', ownCpus, 'timeStrs', {timeStrs}, ...
            'k', k, 'shape', place_shape(places, ownCpus, k, main_shape(places, ownCpus)), 'placed', t, 'looked', t);
    case 'move'
        if reg.Count == 0
            return;
        end
        [view, reg, disabled] = move_jobs(view, reg, disabled, S, varargin{1 : 2});
    case 'forget'
        if reg.Count == 0
            return;
        end
        [jobIds, arrayIds] = varargin{1 : 2};
        for i = 1 : numel(jobIds)
            key = job_key(jobIds(i), arrayIds(i));
            if isKey(reg, key)
                remove(reg, key);
            end
        end
    case 'set'
        for i = 1 : 2 : numel(varargin)
            S.(varargin{i}) = varargin{i + 1};
        end
        view = [];
    otherwise
        error('slurm_alternatives:action', 'Unknown action %s.', action);
end

end


function [places, view] = make_places(view, S, main, alts)
% the main place and the alternatives, with the main settings for the keys
% an alternative leaves out (memPerCPU: its --mem-per-cpu, else its nodes')

names = {'SlurmParam', 'SlurmConstraint', 'memPerCPU', 'maxCPUNum', 'minCPUNum', 'gpuHeadroom', 'cpuJobs'};
place = struct();
for i = 1 : numel(names)
    place.(names{i}) = [];
    if isfield(main, names{i})
        place.(names{i}) = main.(names{i});
    end
end
if isempty(place.SlurmConstraint)
    place.SlurmConstraint = '';
end
% (the main place takes any job, as without alternatives)
place.cpuJobs = true;
place.opts = parse_options([place.SlurmParam, ' ', place.SlurmConstraint]);
places = place;
if nargin < 4 || isempty(alts)
    return;
end
if isstruct(alts)
    alts = num2cell(alts);
end
for a = 1 : numel(alts)
    alt = alts{a};
    if ~isstruct(alt)
        error('slurm_alternatives:entry', 'Each entry of SlurmAlternatives must be a struct (a JSON object).');
    end
    p = place;
    p.memPerCPU = [];
    p.cpuJobs = false;
    keys = fieldnames(alt);
    for i = 1 : numel(keys)
        name = names(strcmpi(names, keys{i}));
        if isempty(name)
            warning('slurm_alternatives:key', 'SlurmAlternatives entry %d: the key %s is ignored (an entry has %s).', ...
                a, keys{i}, strjoin(names, ', '));
            continue;
        end
        p.(name{1}) = alt.(keys{i});
    end
    if ~ischar(p.SlurmParam) && ~isstring(p.SlurmParam)
        error('slurm_alternatives:entry', 'SlurmAlternatives entry %d: SlurmParam must be text.', a);
    end
    p.SlurmParam = char(p.SlurmParam);
    p.SlurmConstraint = char(p.SlurmConstraint);
    p.cpuJobs = isequal(p.cpuJobs, true) || isequal(p.cpuJobs, 1);
    p.opts = parse_options([p.SlurmParam, ' ', p.SlurmConstraint]);
    if isempty(p.memPerCPU) && ~isempty(p.opts.memPerCPU)
        p.memPerCPU = p.opts.memPerCPU / 1024;
    end
    if isempty(p.memPerCPU)
        % the memory per CPU of the place's nodes (the smallest of them)
        view = current_view(view, S);
        [mask, view] = place_nodes(view, p);
        ratio = min(view.memTot(mask) ./ view.cpuTot(mask)) / 1024;
        p.memPerCPU = place.memPerCPU;
        if ~isempty(ratio) && ratio > 0
            p.memPerCPU = ratio;
        end
    end
    places(end + 1) = p; %#ok<AGROW>
end

end


function [timeStrs] = time_strings(places, jobTimeLimit, fmt, mainStr)
% each place's time limit option, as the framework makes it from jobTimeLimit

timeStrs = repmat({''}, numel(places), 1);
timeStrs{1} = mainStr;
for k = 2 : numel(places)
    SlurmParam = places(k).SlurmParam;
    if ~isempty(jobTimeLimit) && ~(contains(SlurmParam, ' -t ') || contains(SlurmParam, '--time'))
        jobTimeLimit_k = max(jobTimeLimit, 1 / 60);
        h = floor(jobTimeLimit_k);
        m = round((jobTimeLimit_k - h) * 60);
        timeStrs{k} = sprintf(fmt, h, m);
    end
end

end


function [opts] = parse_options(str)
% the options of a SlurmParam (and SlurmConstraint) that tell where a job
% may run and what it takes: partition, QOS, account, memory (MB per CPU or
% per node, and the option as given), GPUs, constraint, node lists, job name,
% tasks and time limit; other tokens are kept in order (places that differ
% in them don't trade jobs)

opts = struct('partition', '', 'qos', '', 'account', '', 'memPerCPU', [], 'mem', [], 'memToken', '', ...
    'gres', '', 'gpus', 0, 'constraint', '', 'exclude', '', 'nodelist', '', 'jobName', '', 'ntasks', 1, ...
    'time', '', 'other', {{}});
tokens = strsplit(strtrim(char(str)));
tokens = tokens(~cellfun(@isempty, tokens));
% option names: short (with a separate or attached value), long (with = or
% a separate value)
shortNames = {'-p', '-q', '-A', '-C', '-x', '-w', '-J', '-n', '-t', '-G'};
longNames = {'--partition', '--qos', '--account', '--constraint', '--exclude', '--nodelist', '--job-name', ...
    '--ntasks', '--time', '--gpus', '--mem', '--mem-per-cpu', '--gres', '--gpus-per-node'};
fields = {'partition', 'qos', 'account', 'constraint', 'exclude', 'nodelist', 'jobName', 'ntasks', 'time', 'gpus'};
longFields = [fields(1 : 8), {'time', 'gpus', 'mem', 'memPerCPU', 'gres', 'gpus'}];
i = 1;
while i <= numel(tokens)
    t = tokens{i};
    name = '';
    value = '';
    if startsWith(t, '--')
        eq = find(t == '=', 1);
        if isempty(eq)
            n = t;
        else
            n = t(1 : eq - 1);
        end
        ind = find(strcmp(longNames, n), 1);
        if ~isempty(ind)
            name = longFields{ind};
            if ~isempty(eq)
                value = t(eq + 1 : end);
            elseif i < numel(tokens)
                i = i + 1;
                value = tokens{i};
            end
        end
    elseif numel(t) >= 2 && t(1) == '-'
        ind = find(strcmp(shortNames, t(1 : 2)), 1);
        if ~isempty(ind)
            name = fields{ind};
            if numel(t) > 2
                value = t(3 : end);
            elseif i < numel(tokens)
                i = i + 1;
                value = tokens{i};
            end
        end
    end
    switch name
        case ''
            opts.other{end + 1} = t;
        case 'memPerCPU'
            opts.memPerCPU = memory_mb(value);
            opts.memToken = ['--mem-per-cpu=', value];
        case 'mem'
            opts.mem = memory_mb(value);
            opts.memToken = ['--mem=', value];
        case 'ntasks'
            opts.ntasks = str2double(value);
        case 'gpus'
            opts.gpus = gpu_count(value);
        case 'gres'
            opts.gres = value;
            % gpu, gpu:2, gpu:a100:1
            for g = strsplit(value, ',')
                parts = strsplit(g{1}, ':');
                if strcmp(parts{1}, 'gpu')
                    opts.gpus = str2double(parts{end});
                    if isscalar(parts) || isnan(opts.gpus)
                        opts.gpus = 1;
                    end
                end
            end
        otherwise
            opts.(name) = value;
    end
    i = i + 1;
end

end


function [n] = gpu_count(value)
% the count of a --gpus value ([type:]count)

n = 1;
tok = regexp(value, '(\d+)$', 'tokens', 'once');
if ~isempty(tok)
    n = str2double(tok{1});
end

end


function [mb] = memory_mb(value)
% a Slurm memory size in MB (K, M, G, T suffixes; M by default)

mb = [];
tok = regexp(value, '^([\d.]+)([KMGT]?)B?$', 'tokens', 'once', 'ignorecase');
if isempty(tok)
    return;
end
units = struct('K', 1 / 1024, 'M', 1, 'G', 1024, 'T', 1024 ^ 2);
u = upper(tok{2});
if isempty(u)
    u = 'M';
end
mb = str2double(tok{1}) * units.(u);

end


function [key] = job_key(jobId, arrayId)

key = sprintf('%d_%d', jobId, arrayId);

end


function [t] = now_s()

t = posixtime(datetime('now'));

end


function [name] = place_name(place)

name = place.opts.partition;
if isempty(name)
    name = place.SlurmParam;
end

end


function [shape] = main_shape(places, ownCpus)
% the shape of a job at the main place: its own CPUs and memory there (the
% framework's memPerCPU per CPU if its SlurmParam has no memory option)

o = places(1).opts;
memToken = o.memToken;
if isempty(memToken)
    memToken = sprintf('--mem-per-cpu=%dM', floor(places(1).memPerCPU * 1024));
end
shape = struct('cpus', ownCpus(1), 'memToken', memToken, 'ratio', places(1).memPerCPU, ...
    'param', places(1).SlurmParam);

end


function [shape] = place_shape(places, ownCpus, k, current)
% a job's shape at place k, for a job of shape current: the place's own
% shape if it has as much memory per CPU as the current one (1% less
% counts, as memPerCPU values are rounded), else the current CPUs and memory
% (more CPUs can mean more workers and more memory)

if places(k).memPerCPU >= current.ratio * 0.99
    shape = struct('cpus', ownCpus(k), 'memToken', places(k).opts.memToken, 'ratio', places(k).memPerCPU, ...
        'param', places(k).SlurmParam);
else
    shape = current;
    shape.param = with_memory(places(k).SlurmParam, current.memToken);
end

end


function [param] = with_memory(param, memToken)
% a SlurmParam with its memory option (--mem-per-cpu or --mem) replaced by
% memToken (added if it has none)

tokens = strsplit(strtrim(param));
keep = true(size(tokens));
at = 0;
for i = 1 : numel(tokens)
    if ~keep(i)
        continue;
    end
    t = tokens{i};
    if any(strcmp(t, {'--mem', '--mem-per-cpu'}))
        keep(i : min(i + 1, end)) = false;
    elseif startsWith(t, {'--mem=', '--mem-per-cpu='})
        keep(i) = false;
    else
        continue;
    end
    if at == 0
        at = i;
    end
end
if at == 0
    tokens = [tokens, {memToken}];
else
    rest = tokens(at : end);
    tokens = [tokens(1 : at - 1), {memToken}, rest(keep(at : end))];
end
param = strjoin(tokens, ' ');

end


function [mem] = shape_memory(place, shape)
% the memory (MB) of a job of a shape at a place

o = parse_options(['x ', shape.memToken]);
if ~isempty(o.memPerCPU)
    mem = o.memPerCPU * shape.cpus * place.opts.ntasks;
elseif ~isempty(o.mem)
    mem = o.mem;
else
    mem = place.memPerCPU * 1024 * shape.cpus * place.opts.ntasks;
end

end


function [view] = current_view(view, S)
% the cluster's nodes and the user's pending jobs, read again if older than
% S.viewAge s: free CPUs, memory (MB) and GPUs of each node after the jobs
% that run there and the user's pending jobs that would start first (each
% on the first node of its partition where it fits)

t = now_s();
if ~isempty(view) && t - view.time < S.viewAge
    return;
end
view = struct('time', t, 'testOnly', containers.Map('KeyType', 'char', 'ValueType', 'logical'), ...
    'masks', containers.Map('KeyType', 'char', 'ValueType', 'any'));

% nodes (a failing scontrol gives none: no place has room, and jobs go to
% the main place as without alternatives)
[status, out] = system('scontrol show node -o');
lines = {};
if status == 0
    lines = splitlines(strip(out));
    lines = lines(contains(lines, 'NodeName='));
end
val = @(key) cellfun(@(x) first_token(x, key), lines, 'unif', 0);
num = @(key) reshape(cellfun(@(x) str2double(first_token(x, key)), lines), [], 1);
view.names = reshape(val('NodeName'), [], 1);
view.parts = reshape(cellfun(@(x) strsplit(x, ','), val('Partitions'), 'unif', 0), [], 1);
view.feats = reshape(cellfun(@(x) strsplit(x, ','), val('AvailableFeatures'), 'unif', 0), [], 1);
states = upper(val('State'));
cpuTot = num('CPUEfctv');
tot = num('CPUTot');
cpuTot(isnan(cpuTot)) = tot(isnan(cpuTot));
memTot = num('RealMemory');
memSpec = num('MemSpecLimit');
memSpec(isnan(memSpec)) = 0;
memTot = memTot - memSpec;
gpuTot = reshape(cellfun(@(x) tres_count(first_token(x, 'CfgTRES'), 'gres/gpu'), lines), [], 1);
view.cpuTot = cpuTot;
view.memTot = memTot;
view.gpuTot = gpuTot;
view.freeCPU = cpuTot - num('CPUAlloc');
view.freeMem = memTot - num('AllocMem');
view.freeGPU = gpuTot - reshape(cellfun(@(x) tres_count(first_token(x, 'AllocTRES'), 'gres/gpu'), lines), [], 1);
bad = {'DOWN', 'DRAIN', 'FAIL', 'MAINT', 'RESERVED', 'NOT_RESPONDING', 'POWER', 'REBOOT', 'FUTURE', ...
    'INVALID', 'PLANNED', 'UNKNOWN', 'NO_RESPOND', 'BLOCKED'};
view.usable = reshape(~cellfun(@(x) any(contains(x, bad)), states), [], 1) & ~isnan(view.freeCPU) ...
    & ~isnan(view.freeMem);

% the user's pending jobs (job, partition, QOS, demand, reason), one entry
% per task of a pending array
lines = {};
[status, out] = system('squeue --me -t PD -h -O ''JobArrayID:|,Partition:|,QOS:|,NumCPUs:|,tres-alloc:|,Reason:|''');
if status == 0
    lines = splitlines(strip(out));
end
rows = cell(numel(lines), 1);
for l = 1 : numel(lines)
    f = strsplit(regexprep(lines{l}, '\|\s*$', ''), '|');
    if numel(f) < 6
        continue;
    end
    [keys, nTasks] = pending_keys(strtrim(f{1}));
    tres = strtrim(f{5});
    cpus = tres_count(tres, 'cpu');
    if isnan(cpus)
        cpus = str2double(f{4});
    end
    mem = [];
    tok = regexp(tres, '(?:^|,)mem=([^,]+)', 'tokens', 'once');
    if ~isempty(tok)
        mem = memory_mb(tok{1});
    end
    if isempty(mem)
        mem = 0;
    end
    nodes = max(1, tres_count(tres, 'node'));
    keys(end + 1 : max(1, nTasks)) = {''};
    rows{l} = struct('key', keys(:), 'partition', strtrim(f{2}), 'qos', strtrim(f{3}), 'cpus', cpus / nodes, ...
        'mem', mem / nodes, 'gpus', tres_count(tres, 'gres/gpu') / nodes, ...
        'reason', strtrim(strjoin(f(6 : end), '|')), 'node', 0);
end
view.pend = cat(1, rows{:}, struct('key', {}, 'partition', {}, 'qos', {}, 'cpus', {}, 'mem', {}, 'gpus', {}, ...
    'reason', {}, 'node', {}));
view.pendIndex = containers.Map('KeyType', 'char', 'ValueType', 'any');
keyed = find(~cellfun(@isempty, {view.pend.key}));
for i = keyed
    view.pendIndex(view.pend(i).key) = i;
end

% they start first: each takes its demand from the first usable node of its
% partition where it fits (node; one that fits nowhere takes nothing)
partMasks = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 1 : numel(view.pend)
    p = view.pend(i);
    if ~any(strcmp(p.reason, {'Resources', 'Priority', 'None'}))
        continue;
    end
    if ~isKey(partMasks, p.partition)
        partMasks(p.partition) = view.usable & cellfun(@(x) any(strcmp(x, p.partition)), view.parts);
    end
    j = find(partMasks(p.partition) & view.freeCPU >= p.cpus & view.freeMem >= p.mem & view.freeGPU >= p.gpus, 1);
    if ~isempty(j)
        view = take(view, j, p.cpus, p.mem, p.gpus);
        view.pend(i).node = j;
    end
end
% the QOS and association limits that pending jobs wait for (partition, QOS)
limited = ~cellfun(@isempty, regexp({view.pend.reason}, '^(QOS|Assoc)', 'once'));
view.limitParts = {view.pend(limited).partition};
view.limitQos = {view.pend(limited).qos};
% the partitions where jobs of the user wait for resources (they fit no
% node): a new job there waits behind them, also where it would fit (Slurm
% starts jobs by priority, and its backfill tests only so many of them)
waiting = ismember({view.pend.reason}, {'Resources', 'Priority', 'None'}) & [view.pend.node] == 0;
view.waitParts = unique({view.pend(waiting).partition});
% the partitions where GPU jobs of the user are pending, for any reason (held
% ones too): CPU-only jobs keep out of them
view.gpuParts = unique({view.pend([view.pend.gpus] > 0).partition});

end


function [view] = take(view, j, cpus, mem, gpus)

view.freeCPU(j) = view.freeCPU(j) - cpus;
view.freeMem(j) = view.freeMem(j) - mem;
view.freeGPU(j) = view.freeGPU(j) - gpus;

end


function [v] = first_token(line, key)
% the value of key=value in a line of scontrol's one-line output

v = '';
tok = regexp(line, ['(?:^|\s)', key, '=(\S*)'], 'tokens', 'once');
if ~isempty(tok) && ~strcmp(tok{1}, '(null)')
    v = tok{1};
end

end


function [n] = tres_count(tres, name)
% the count of name in a TRES list (cpu=4,mem=10G,gres/gpu=1); 0 for a GPU
% or node count that isn't there, NaN for others

tok = regexp(tres, ['(?:^|,)', regexptranslate('escape', name), '=([^,]+)'], 'tokens', 'once');
if isempty(tok)
    n = NaN;
    if strcmp(name, 'gres/gpu') || strcmp(name, 'node')
        n = 0;
    end
    return;
end
n = str2double(tok{1});

end


function [keys, n] = pending_keys(id)
% the job keys (job_task) of squeue's JobArrayID of pending jobs: 12_[3],
% 12_3, 12 (no array: never one of ours), 12_[1-4,7%2]

keys = {};
tok = regexp(id, '^(\d+)_\[?([^\]]*)\]?$', 'tokens', 'once');
if isempty(tok)
    n = 1;
    return;
end
spec = regexprep(tok{2}, '%\d+$', '');
tasks = [];
for part = strsplit(spec, ',')
    r = str2double(strsplit(part{1}, '-'));
    if numel(r) == 2
        tasks = [tasks, r(1) : r(2)]; %#ok<AGROW>
    elseif ~isnan(r)
        tasks(end + 1) = r; %#ok<AGROW>
    end
end
% (at most 1000 tasks of an array: more never start at once)
tasks = tasks(1 : min(end, 1000));
keys = arrayfun(@(x) sprintf('%s_%d', tok{1}, x), tasks, 'unif', 0);
n = numel(tasks);

end


function [mask, view] = place_nodes(view, place)
% the nodes of a place: in its partition, with its features, within its
% node lists (a place without a partition: none, its room comes from
% Slurm's check alone)

o = place.opts;
key = [o.partition, '|', o.constraint, '|', o.exclude, '|', o.nodelist];
if isKey(view.masks, key)
    mask = view.masks(key);
    return;
end
mask = false(numel(view.names), 1);
if ~isempty(o.partition)
    mask = cellfun(@(x) any(strcmp(x, o.partition)), view.parts);
end
if ~isempty(o.constraint)
    % a&b, a|b and single features (others are left to Slurm's check)
    if contains(o.constraint, '&') && ~contains(o.constraint, '|')
        need = strsplit(o.constraint, '&');
        mask = mask & cellfun(@(x) all(ismember(need, x)), view.feats);
    elseif ~any(ismember(o.constraint, '[]()*&'))
        either = strsplit(o.constraint, '|');
        mask = mask & cellfun(@(x) any(ismember(either, x)), view.feats);
    end
end
% (node names with or without their domain)
short = regexprep(view.names, '\..*$', '');
if ~isempty(o.exclude)
    mask = mask & ~ismember(short, regexprep(expand_hostlist(o.exclude), '\..*$', ''));
end
if ~isempty(o.nodelist)
    mask = mask & ismember(short, regexprep(expand_hostlist(o.nodelist), '\..*$', ''));
end
view.masks(key) = mask;

end


function [hosts] = expand_hostlist(str)
% the names of a Slurm host list (a,b, n00[01-03,07].x)

hosts = {};
parts = regexp(str, '[^,\[]+(\[[^\]]*\])?[^,]*', 'match');
for i = 1 : numel(parts)
    tok = regexp(parts{i}, '^([^\[]*)\[([^\]]*)\](.*)$', 'tokens', 'once');
    if isempty(tok)
        hosts{end + 1} = parts{i}; %#ok<AGROW>
        continue;
    end
    for r = strsplit(tok{2}, ',')
        b = strsplit(r{1}, '-');
        w = numel(b{1});
        for x = str2double(b{1}) : str2double(b{end})
            hosts{end + 1} = sprintf('%s%0*d%s', tok{1}, w, x, tok{3}); %#ok<AGROW>
        end
    end
end

end


function [availCPU, availMem] = cpu_job_room(view, place, gpus)
% each node's CPUs and memory a job may take: CPU-only jobs leave each idle
% GPU its share (gpuHeadroom, or the node's CPUs and memory per GPU)

availCPU = view.freeCPU;
availMem = view.freeMem;
if gpus > 0
    return;
end
idle = max(0, view.freeGPU);
h = place.gpuHeadroom;
if isempty(h)
    hCPU = view.cpuTot ./ max(1, view.gpuTot);
    hMem = view.memTot ./ max(1, view.gpuTot);
else
    h = double(h(:)');
    if isscalar(h)
        h = [h, h];
    end
    hCPU = h(1);
    hMem = h(2) * 1024;
end
availCPU = availCPU - idle .* hCPU;
availMem = availMem - idle .* hMem;

end


function [slots, j, view] = place_room(view, place, shape)
% the jobs of a shape that a place can start now by the view, and the first
% node with room for one

o = place.opts;
mem = shape_memory(place, shape);
need = shape.cpus * o.ntasks;
[mask, view] = place_nodes(view, place);
[availCPU, availMem] = cpu_job_room(view, place, o.gpus);
per = min(floor(availCPU / need), floor(availMem / max(mem, 1)));
if o.gpus > 0
    per = min(per, floor(view.freeGPU / o.gpus));
end
per(~mask | ~view.usable) = 0;
per = max(0, per);
slots = sum(per);
j = find(per > 0, 1);
% a QOS or association limit that a pending job of the user waits for here
% holds new jobs too, and so do jobs of the user that wait here for resources
if slots > 0 && (any(strcmp(view.limitParts, o.partition) & (isempty(o.qos) | strcmp(view.limitQos, o.qos))) ...
        || any(strcmp(view.waitParts, o.partition)))
    slots = 0;
    j = [];
end
% CPU-only jobs keep out of partitions where GPU jobs of the user are pending
if slots > 0 && o.gpus == 0 && any(strcmp(view.gpuParts, o.partition))
    slots = 0;
    j = [];
end

end


function [now_ok, view, disabled] = slurm_says_now(view, place, shape, timeStr, S, disabled, isMain)
% Slurm's own check (sbatch --test-only, once per view for a place and
% shape): whether the job would start now; a place Slurm refuses (partition,
% account, QOS) is not used again in this session

key = sprintf('%s|%s|%d|%s', shape.param, place.SlurmConstraint, shape.cpus, timeStr);
if isKey(view.testOnly, key)
    now_ok = view.testOnly(key);
    return;
end
cmd = sprintf('sbatch --test-only --cpus-per-task=%d %s %s %s --wrap="true" 2>&1', shape.cpus, ...
    shape.param, place.SlurmConstraint, timeStr);
[status, out] = system(cmd);
tok = regexp(out, 'to start at (\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d)', 'tokens', 'once');
now_ok = false;
if status ~= 0 || isempty(tok)
    % (a place Slurm refuses: not used again; an error of the moment, e.g.
    % a busy controller, or the main place: no room now)
    transient = contains(lower(out), {'timed out', 'unable to contact', 'try again', 'temporarily', 'busy'});
    if ~transient && ~isMain
        disabled{end + 1} = place.SlurmParam;
        warning('slurm_alternatives:testOnly', 'Slurm refuses jobs with "%s" (%s); this place is not used.', ...
            place.SlurmParam, strtrim(out));
    end
else
    start = datetime(tok{1}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
    now_ok = seconds(start - datetime('now')) <= S.startTolerance;
end
view.testOnly(key) = now_ok;

end


function [ok, j, view, disabled] = has_room(view, place, shape, timeStr, S, disabled, isMain)
% whether a place can start a job of a shape now, and the node it takes
% (isMain: the main place, which is never left out)

ok = false;
j = [];
if any(strcmp(disabled, place.SlurmParam))
    return;
end
[slots, j, view] = place_room(view, place, shape);
if slots < 1 && ~isempty(place.opts.partition)
    return;
end
[ok, view, disabled] = slurm_says_now(view, place, shape, timeStr, S, disabled, isMain);
% (a place without a partition: Slurm's check alone, one job per view)
if ok && isempty(place.opts.partition)
    view.testOnly(sprintf('%s|%s|%d|%s', shape.param, place.SlurmConstraint, shape.cpus, timeStr)) = false;
end

end


function [view] = take_room(view, place, shape, j)
% a job placed at a place takes its resources from node j of the view

if isempty(j)
    return;
end
view = take(view, j, shape.cpus * place.opts.ntasks, shape_memory(place, shape), place.opts.gpus);

end


function [tf] = gpu_nodes_only(view, place)
% whether all nodes of a place have GPUs

[mask, ~] = place_nodes(view, place);
tf = any(mask) && all(view.gpuTot(mask) > 0);

end


function [k, shape, view, disabled] = choose_place(view, places, ownCpus, timeStrs, S, disabled)
% the first place with room for the job now (shaped for the main place, or
% for each place, see place_shape); else the first place a CPU-only job may
% wait in (not of GPU nodes only), else the main place

current = main_shape(places, ownCpus);
for k = 1 : numel(places)
    if ~takes_job(view, places(k))
        continue;
    end
    shape = place_shape(places, ownCpus, k, current);
    [ok, j, view, disabled] = has_room(view, places(k), shape, timeStrs{k}, S, disabled, k == 1);
    if ok
        view = take_room(view, places(k), shape, j);
        return;
    end
end
k = wait_place(view, places, disabled);
shape = place_shape(places, ownCpus, k, current);

end


function [tf] = takes_job(view, place)
% whether a place takes the job: alternatives of GPU nodes only take CPU-only
% jobs only if they say so (cpuJobs)

tf = place.cpuJobs || place.opts.gpus > 0 || ~gpu_nodes_only(view, place);

end


function [k] = wait_place(view, places, disabled)
% the place where a job waits when no place has room: for CPU-only jobs the
% first place that isn't of GPU nodes only, else the main place

k = 1;
for i = 1 : numel(places)
    if any(strcmp(disabled, places(i).SlurmParam))
        continue;
    end
    if places(i).opts.gpus > 0 || ~gpu_nodes_only(view, places(i))
        k = i;
        return;
    end
end

end


function [fields] = move_fields(from, to, current, shape)
% the scontrol update fields that move a job of shape current from place
% from to place to with a shape ('' if the places differ in options a
% pending job can't change)

fields = '';
a = from.opts;
b = to.opts;
if ~isequal(a.other, b.other) || a.ntasks ~= b.ntasks || ~strcmp(a.time, b.time) || isempty(b.partition)
    return;
end
f = {sprintf('Partition=%s', b.partition)};
% Slurm keeps the QOS when only the partition changes, also one the new
% partition's association doesn't allow, and can't clear it: the QOS is
% always set (a place without one can't take jobs of another partition,
% account or QOS)
if ~strcmp(a.partition, b.partition) || ~strcmp(a.account, b.account) || ~strcmp(a.qos, b.qos)
    if isempty(b.qos)
        return;
    end
    f{end + 1} = sprintf('QOS=%s', b.qos);
end
if ~strcmp(a.account, b.account)
    if isempty(b.account)
        return;
    end
    f{end + 1} = sprintf('Account=%s', b.account);
end
f{end + 1} = sprintf('NumCPUs=%d CPUsPerTask=%d', shape.cpus * b.ntasks, shape.cpus);
m = parse_options(['x ', shape.memToken]);
if ~isempty(m.memPerCPU)
    f{end + 1} = sprintf('MinMemoryCPU=%d', ceil(m.memPerCPU));
elseif ~isempty(m.mem)
    f{end + 1} = sprintf('MinMemoryNode=%d', ceil(m.mem));
elseif ~isempty(current.memToken) || ~strcmp(a.partition, b.partition)
    % (the partition's default memory can't be set back)
    return;
end
if ~strcmp(a.gres, b.gres)
    if isempty(b.gres) || isempty(a.gres)
        return;
    end
    f{end + 1} = sprintf('Gres=%s', b.gres);
elseif a.gpus ~= b.gpus
    return;
end
pairs = {'constraint', 'Features'; 'exclude', 'ExcNodeList'; 'nodelist', 'ReqNodeList'; 'jobName', 'JobName'};
for i = 1 : size(pairs, 1)
    if ~strcmp(a.(pairs{i, 1}), b.(pairs{i, 1}))
        v = b.(pairs{i, 1});
        if strcmp(pairs{i, 1}, 'jobName') && isempty(v)
            v = 'wrap';
        end
        f{end + 1} = sprintf('%s=%s', pairs{i, 2}, v); %#ok<AGROW>
    end
end
fields = strjoin(f, ' ');

end


function [view, reg, disabled] = move_jobs(view, reg, disabled, S, jobIds, arrayIds)
% moves recorded pending jobs to places with room (see the usage above)

t = now_s();
cands = {};
for i = 1 : numel(jobIds)
    key = job_key(jobIds(i), arrayIds(i));
    if isKey(reg, key)
        e = reg(key);
        if t - e.placed >= S.minWait && t - e.looked >= S.lookInterval
            cands{end + 1} = key; %#ok<AGROW>
        end
    end
end
if isempty(cands)
    return;
end
view = current_view(view, S);
% (counts for the line about this look: room where they wait, moved, no room
% elsewhere)
nFit = 0;
nMoved = 0;
nStay = 0;
for c = 1 : numel(cands)
    key = cands{c};
    e = reg(key);
    e.looked = t;
    reg(key) = e;
    p = [];
    if isKey(view.pendIndex, key)
        p = view.pendIndex(key);
    end
    if isempty(p)
        continue;
    end
    % only jobs that wait for resources or a QOS or association limit
    % (not held, waiting for a dependency or start time, or not yet seen by
    % the scheduler), and that don't fit where they are (on a node of their
    % place; a CPU-only job at GPU nodes where GPU jobs of the user are
    % pending doesn't stay there)
    here = e.places(e.k);
    [mask, view] = place_nodes(view, here);
    gpuWork = here.opts.gpus == 0 && any(strcmp(view.gpuParts, here.opts.partition)) && gpu_nodes_only(view, here);
    if ~gpuWork && view.pend(p).node > 0 && (mask(view.pend(p).node) || isempty(here.opts.partition))
        nFit = nFit + 1;
        continue;
    end
    reason = view.pend(p).reason;
    if ~any(strcmp(reason, {'Resources', 'Priority'})) && isempty(regexp(reason, '^(QOS|Assoc)', 'once'))
        continue;
    end
    places = e.places;
    j = e.k;
    % a job that something else moved (e.g. a user's own script) is left to it
    if ~isempty(places(j).opts.partition) && ~strcmp(view.pend(p).partition, places(j).opts.partition)
        remove(reg, key);
        continue;
    end
    target = 0;
    for k = 1 : numel(places)
        if k == j || ~takes_job(view, places(k))
            continue;
        end
        shape = place_shape(places, e.cpus, k, e.shape);
        if isempty(move_fields(places(j), places(k), e.shape, shape))
            continue;
        end
        [ok, node, view, disabled] = has_room(view, places(k), shape, e.timeStrs{k}, S, disabled, k == 1);
        if ok
            target = k;
            break;
        end
    end
    % a CPU-only job doesn't wait in the queue of GPU nodes
    if target == 0 && places(j).opts.gpus == 0 && gpu_nodes_only(view, places(j))
        k = wait_place(view, places, disabled);
        shape = place_shape(places, e.cpus, k, e.shape);
        if k ~= j && ~isempty(move_fields(places(j), places(k), e.shape, shape))
            target = k;
            node = [];
        end
    end
    if target == 0
        nStay = nStay + 1;
        continue;
    end
    cmd = sprintf('scontrol update JobId=%s %s', key, move_fields(places(j), places(target), e.shape, shape));
    [status, out] = system(cmd);
    if status ~= 0
        fprintf('Job %s was not moved to %s: %s\n', key, place_name(places(target)), strtrim(out));
        nStay = nStay + 1;
        continue;
    end
    fprintf('Job %s moved from %s (%d CPUs) to %s (%d CPUs).\n', key, place_name(places(j)), e.shape.cpus, ...
        place_name(places(target)), shape.cpus);
    view = take_room(view, places(target), shape, node);
    e.k = target;
    e.shape = shape;
    e.placed = t;
    reg(key) = e;
    nMoved = nMoved + 1;
end
if nFit + nMoved + nStay > 0
    fprintf('%s: of %d waiting jobs, %d have room where they wait, %d moved, %d have no room elsewhere.\n', ...
        char(datetime('now', 'Format', 'HH:mm:ss')), nFit + nMoved + nStay, nFit, nMoved, nStay);
end

end
