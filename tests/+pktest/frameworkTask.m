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

fid = fopen(logFn, 'a');
fprintf(fid, '%s\n', outputFn);
fclose(fid);

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
