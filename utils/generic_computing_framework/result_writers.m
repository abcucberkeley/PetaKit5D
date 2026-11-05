function [varargout] = result_writers(action, varargin)
% records of the tasks that make a result in parts with jobs of their own
% (large-file decon, deskew/rotate and crop), so that another task of the
% same result, e.g. of a call on the same data started while the first
% call's jobs still run, waits for the one at work instead of making the
% result a second time (and the two of them publishing it twice). A task
% records its process (Slurm job, node, process id and start) in its flag
% folder (<result folder>/zarr_flag/<name>_<uuid>/writer.txt), which goes
% with the folder once the result is published; a task that ends otherwise
% (an error, Ctrl+C) marks its record ended, and a task whose job was
% killed is seen by squeue. Tasks of the same uuid (one call and its
% retries) share the record.
%
%   result_writers('register', resultFullpath, uuid, flagPath)
%       records this process as the writer of uuid, in flagPath
%   result_writers('release', resultFullpath, uuid, flagPath)
%       marks the record in flagPath ended if it is this process's and still
%       there (for an onCleanup of the task: a published result's flags are
%       gone)
%   writers = result_writers('live', resultFullpath, uuid)
%       the owners of the records of other uuids that may still work on the
%       result: not ended, and their Slurm job (not this process's own) is
%       in the queue, or (outside Slurm, or in this process's own job) their
%       process runs on this node (Linux). Others (a process that ended, or
%       runs on another machine outside Slurm) are not waited for, as before.
%   done = result_writers('wait', resultFullpath, uuid)
%       after 'register': waits while a writer recorded before this one may
%       still work; true if the result exists then (another writer made it),
%       false if none works on it any more
%   tf = result_writers('orphaned', flagFullpath, maxAge)
%       for a block task of such a task (flagFullpath: the block's flag):
%       true if the writer recorded in the flag's folder has ended (its
%       record says so, or its Slurm job surely ended: squeue answered
%       without it; a failed squeue counts as running), e.g. one that was
%       cancelled: then nothing would use the block, and its job skips it. A
%       writer started again with the same uuid records itself again, and
%       uses the blocks. The record is read each time; squeue's answer about
%       a job is kept in this process for maxAge seconds (optional, default
%       60), so that a job that runs many blocks of one writer (task
%       runners) asks once a minute, and notices a cancel up to that late.


if strcmp(action, 'orphaned')
    varargout{1} = is_orphaned(varargin{:});
    return;
end
[resultFullpath, uuid] = varargin{1 : 2};
% (the tasks spell paths with repeated slashes; a leading // stays)
resultFullpath = strip(regexprep(resultFullpath, '(?<!^)/+', '/'), 'right', '/');
switch action
    case 'register'
        register_writer(resultFullpath, uuid, varargin{3});
    case 'release'
        release_writer(uuid, varargin{3});
    case 'live'
        varargout{1} = live_writers(resultFullpath, uuid, false);
    case 'wait'
        varargout{1} = wait_for_writers(resultFullpath, uuid);
end

end


function [tf] = is_orphaned(flagFullpath, maxAge)
% whether the writer recorded next to a block's flag has ended: its record
% says so, or it had a Slurm job (not this process's) that squeue says has
% ended (squeue's answers kept for maxAge seconds per job)

% (job, time of the answer in seconds, orphaned)
persistent answers
if isempty(answers)
    answers = struct('job', {}, 'time', {}, 'orphaned', {});
end
if nargin < 2 || isempty(maxAge)
    maxAge = 60;
end
tf = false;
fn = [fileparts(flagFullpath), '/writer.txt'];
if ~isfile(fn)
    return;
end
try
    txt = fileread(fn);
catch
    return;
end
if record_ended(txt)
    tf = true;
    owner = regexp(txt, 'owner=([^\n]*)', 'tokens', 'once');
    fprintf('The task that runs this block (%s) has ended, skip the block.\n', strjoin(owner, ''));
    return;
end
rec = regexp(txt, 'owner=(\d+)_(\d*):', 'tokens', 'once');
if isempty(rec) || str2double(rec{1}) <= 0
    return;
end
[myJob, myArr] = own_job();
if str2double(rec{1}) == myJob && strcmp(rec{2}, myArr)
    return;
end
jobStr = rec{1};
if ~isempty(rec{2})
    jobStr = sprintf('%s_%s', rec{1}, rec{2});
end
t = now * 86400;
i = find(strcmp({answers.job}, jobStr), 1);
if ~isempty(i) && t - answers(i).time < maxAge
    tf = answers(i).orphaned;
else
    [status, out] = system(sprintf('squeue -h -j %s -o %%T', jobStr));
    tf = (status == 0 && isempty(strtrim(out))) || contains(out, 'Invalid job id');
    if isempty(i)
        i = numel(answers) + 1;
    end
    answers(i) = struct('job', jobStr, 'time', t, 'orphaned', tf);
end
if tf
    fprintf('The task that runs this block (job %s) has ended, skip the block.\n', jobStr);
end

end


function [] = register_writer(resultFullpath, uuid, flagPath)
% writer.txt in the flag folder: the result, the uuid and the owner

write_record(flagPath, sprintf('result=%s\nuuid=%s\nowner=%s\n', resultFullpath, uuid, my_owner()));

end


function [] = release_writer(uuid, flagPath)
% the record in flagPath marked ended, if it is this process's (of uuid)
% and not ended yet; nothing if it is gone (the result was published)

fn = [strip(flagPath, 'right', '/'), '/writer.txt'];
try
    txt = fileread(fn);
catch
    return;
end
rec = regexp(txt, 'uuid=(?<uuid>[^\n]*)\nowner=(?<owner>[^\n]*)', 'names', 'once');
if isempty(rec) || ~strcmp(rec.uuid, uuid) || ~strcmp(rec.owner, my_owner()) || record_ended(txt)
    return;
end
% (it runs when the task ends, also by an error: it must not error itself)
try
    write_record(flagPath, sprintf('%s\nended=1\n', strtrim(txt)));
catch
end

end


function [] = write_record(flagPath, txt)
% writer.txt in flagPath written at once (a temporary file moved into place)

fn = [strip(flagPath, 'right', '/'), '/writer.txt'];
tmpFn = sprintf('%s_%s.txt', fn(1 : end - 4), get_uuid());
fid = fopen(tmpFn, 'w');
if fid < 0
    return;
end
fprintf(fid, '%s', txt);
fclose(fid);
movefile(tmpFn, fn);

end


function [tf] = record_ended(txt)
% whether a record says its writer has ended

tf = ~isempty(regexp(txt, '^ended=', 'once', 'lineanchors'));

end


function [writers, records] = live_writers(resultFullpath, uuid, earlierOnly)
% owners of the other uuids' records of the result that may still work (with
% earlierOnly, of those recorded before this uuid's record)

writers = {};
records = struct('owner', {}, 'uuid', {}, 'datenum', {});
[resultPath, fsname] = fileparts(resultFullpath);
flagRoot = [resultPath, '/zarr_flag/'];
if ~isfolder(flagRoot)
    return;
end
d = dir([flagRoot, fsname, '_*/writer.txt']);
mine = [];
for i = 1 : numel(d)
    try
        txt = fileread([d(i).folder, '/', d(i).name]);
    catch
        continue;
    end
    rec = regexp(txt, 'result=(?<result>[^\n]*)\nuuid=(?<uuid>[^\n]*)\nowner=(?<owner>[^\n]*)', 'names', 'once');
    if isempty(rec) || ~strcmp(rec.result, resultFullpath)
        continue;
    end
    if strcmp(rec.uuid, uuid)
        mine = d(i).datenum;
        continue;
    end
    % (a writer that has ended works on it no more)
    if record_ended(txt)
        continue;
    end
    records(end + 1) = struct('owner', rec.owner, 'uuid', rec.uuid, 'datenum', d(i).datenum);
end
if earlierOnly && ~isempty(mine)
    % (recorded at the same time: the smaller uuid first)
    earlier = [records.datenum] < mine | ([records.datenum] == mine & cellfun(@(x) string(x) < string(uuid), {records.uuid}));
    records = records(earlier);
end
records = records(owners_live({records.owner}));
writers = {records.owner};

end


function [done] = wait_for_writers(resultFullpath, uuid)
% waits while a writer recorded before this one may still work on the result

done = false;
announced = false;
while true
    if batch_file_exist({resultFullpath})
        done = true;
        fprintf('The result %s was made by another job.\n', resultFullpath);
        return;
    end
    writers = live_writers(resultFullpath, uuid, true);
    if isempty(writers)
        return;
    end
    if ~announced
        fprintf('Another job (%s) makes the result %s, wait for it...\n', writers{1}, resultFullpath);
        announced = true;
    end
    % (back at once when the result comes)
    wait_for_files({{resultFullpath}}, 60);
end

end


function [tf] = owners_live(owners)
% whether the owners <job>_<array task>:<node>:<pid>[:<start>] may still
% work: a Slurm job other than this process's own that is pending or
% running (one squeue for all array jobs), or else (outside Slurm, or in
% this process's own job, e.g. the record of a call this process made
% before, or of a job started again) a process other than this one that
% runs on this node

tf = false(size(owners));
tok = regexp(owners, '^(\d+)_(\d*):([^:]*):(\d+):?(\d*)$', 'tokens', 'once');
[myJob, myArr] = own_job();
jobs = zeros(size(owners));
for i = 1 : numel(owners)
    if isempty(tok{i})
        continue;
    end
    jobs(i) = str2double(tok{i}{1});
    if jobs(i) <= 0 || (jobs(i) == myJob && strcmp(tok{i}{2}, myArr))
        jobs(i) = 0;
        tf(i) = process_live(tok{i}{3}, tok{i}{4}, tok{i}{5});
    end
end
inJob = find(jobs > 0);
if isempty(inJob)
    return;
end
% array jobs (the jobs of the frameworks) with one squeue; others (a master's
% own job) one by one, without an array task
isArray = cellfun(@(x) ~isempty(x{2}), tok(inJob));
inds = inJob(isArray);
if ~isempty(inds)
    status = check_batch_slurm_jobs_status(jobs(inds), cellfun(@(x) str2double(x{2}), tok(inds)));
    tf(inds) = status >= 0;
end
for i = reshape(inJob(~isArray), 1, [])
    tf(i) = check_slurm_job_status(jobs(i)) >= 0;
end

end


function [tf] = process_live(node, pid, start)
% whether the process pid (char) runs on this node and is not this one: on
% Linux (/proc), and the one that started at start if it is known (a
% process id is used again once its process ends)

tf = strcmp(node, node_name()) && str2double(pid) ~= feature('getpid') && isfolder(['/proc/', pid]) ...
    && (isempty(start) || strcmp(start, process_start(pid)));

end


function [owner] = my_owner()
% this process as the owner of a record: <job>_<array task>:<node>:<pid>:<start>

[job, arr] = own_job();
pid = feature('getpid');
owner = sprintf('%d_%s:%s:%d:%s', job, arr, node_name(), pid, process_start(pid));

end


function [job, arr] = own_job()
% this process's Slurm job and array task ('' if not an array job; 0
% without a job)

job = str2double(getenv('SLURM_ARRAY_JOB_ID'));
arr = getenv('SLURM_ARRAY_TASK_ID');
if isnan(job)
    job = str2double(getenv('SLURM_JOB_ID'));
    arr = '';
end
if isnan(job)
    job = 0;
    arr = '';
end

end


function [start] = process_start(pid)
% the start time of a process (clock ticks after boot, the 22nd field of
% /proc/<pid>/stat on Linux; '' elsewhere or if it can't be read)

start = '';
if ~isunix || ismac
    return;
end
if isnumeric(pid)
    pid = sprintf('%d', pid);
end
try
    txt = fileread(['/proc/', pid, '/stat']);
catch
    return;
end
% (the fields after the name in parentheses, which may have spaces)
fields = strsplit(strtrim(txt(find(txt == ')', 1, 'last') + 1 : end)));
if numel(fields) >= 20
    start = fields{20};
end

end


function [node] = node_name()
% the node name (as the task runners record it)

node = getenv('SLURMD_NODENAME');
if isempty(node) && ispc
    node = getenv('COMPUTERNAME');
elseif isempty(node)
    [~, node] = system('uname -n');
    node = strtrim(node);
end

end

