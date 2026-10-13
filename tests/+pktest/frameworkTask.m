function [] = frameworkTask(inputFn, outputFn, logFn, mode)
% task function for generic computing framework tests. Every call appends a
% line to logFn (so tests can count runs), then depending on mode:
%   'copy'  : copy the input file to the output file
%   'error' : throw pktest:frameworkTask:failed
%   'noop'  : do nothing (task "finishes" without writing its output)

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
end

end
