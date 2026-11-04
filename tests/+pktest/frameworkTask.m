function [] = frameworkTask(inputFn, outputFn, logFn, mode)
% task function for generic computing framework tests. Every call appends a
% line to logFn (so tests can count runs), then depending on mode:
%   'copy'       : copy the input file to the output file
%   'error'      : throw pktest:frameworkTask:failed
%   'noop'       : do nothing (task "finishes" without writing its output)
%   'secondRun'  : as 'noop' in its first run, as 'copy' from its second on
%                  (a task done in its last try)
%   'where'      : write 'worker' or 'client' to the output file, whether it
%                  runs on a parallel pool worker or not
%   'clientonly' : as 'where', but throw on a parallel pool worker

if nargin < 4
    mode = 'copy';
end

append_line(logFn, outputFn);

switch mode
    case 'copy'
        copyfile(inputFn, outputFn);
    case 'error'
        error('pktest:frameworkTask:failed', 'Task failed on purpose for %s', inputFn);
    case 'noop'
    case 'secondRun'
        lines = splitlines(fileread(logFn));
        if nnz(strcmp(lines, outputFn)) >= 2
            copyfile(inputFn, outputFn);
        end
    case {'where', 'clientonly'}
        onWorker = ~isempty(which('getCurrentTask')) && ~isempty(getCurrentTask());
        if onWorker && strcmp(mode, 'clientonly')
            error('pktest:frameworkTask:failed', 'Task failed on purpose on a worker for %s', inputFn);
        end
        where = 'client';
        if onWorker
            where = 'worker';
        end
        fid = fopen(outputFn, 'w');
        fprintf(fid, '%s', where);
        fclose(fid);
end

end


function [] = append_line(logFn, line)
% appends a line to the log while holding a lock folder: tasks that run at
% once in other processes (pool workers, jobs) append to the same log, and on
% Windows two appends at the same time can overwrite each other. A lock left
% by a killed task is ignored after 30 s. (No pause: tests replace it.)

lockDir = [logFn, '.lock'];
t0 = tic;
locked = false;
while toc(t0) < 30
    % made by this task: the lock; already there, or (Windows) not made while
    % another task removes it: try again
    [status, ~, msgID] = mkdir(lockDir);
    if status && isempty(msgID)
        locked = true;
        break;
    end
    wait_briefly();
end
% (Windows can refuse an open for a moment after another process's close)
fid = fopen(logFn, 'a');
while fid < 0 && toc(t0) < 60
    wait_briefly();
    fid = fopen(logFn, 'a');
end
fprintf(fid, '%s\n', line);
fclose(fid);
if locked
    t1 = tic;
    while ~rmdir(lockDir) && toc(t1) < 10
        wait_briefly();
    end
end

end


function [] = wait_briefly()
% about 2 ms

t = tic;
while toc(t) < 0.002
end

end
