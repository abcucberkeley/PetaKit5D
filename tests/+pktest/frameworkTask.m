function [] = frameworkTask(inputFn, outputFn, logFn, mode)
% task function for generic computing framework tests. Every call appends a
% line to logFn (so tests can count runs), then depending on mode:
%   'copy'       : copy the input file to the output file
%   'error'      : throw pktest:frameworkTask:failed
%   'noop'       : do nothing (task "finishes" without writing its output)
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
    [~, ~, msgID] = mkdir(lockDir);
    if ~strcmp(msgID, 'MATLAB:MKDIR:DirectoryExists')
        locked = true;
        break;
    end
    t1 = tic;
    while toc(t1) < 0.002
    end
end
fid = fopen(logFn, 'a');
fprintf(fid, '%s\n', line);
fclose(fid);
if locked
    rmdir(lockDir);
end

end
