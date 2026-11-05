function [kept] = keep_other_result(resultFullpath, tmpFullpath, flagPath)
% before a task moves its temporary result tmpFullpath into place as
% resultFullpath: whether another writer of the same result (a job of
% another call on the same data, or a task started again) already put its
% result there. Then the task keeps that result, as both are the same, and
% removes its own temporary result and its flag folder flagPath (optional),
% instead of deleting a result that a reader (the next step of a pipeline)
% may already use and doing the work twice. A task that overwrites results
% removes them when it starts, so a result there now came from another
% writer. batch_file_exist sees a result another node just made, which
% exist() can miss for 30 s or more.

if nargin < 3
    flagPath = '';
end

kept = batch_file_exist({resultFullpath});
if ~kept
    return;
end
fprintf('The result %s was made by another job, keep it.\n', resultFullpath);
if ~isempty(tmpFullpath) && (exist(tmpFullpath, 'dir') || exist(tmpFullpath, 'file'))
    if exist(tmpFullpath, 'dir')
        rmdir(tmpFullpath, 's');
    else
        delete(tmpFullpath);
    end
end
if ~isempty(flagPath) && exist(flagPath, 'dir')
    rmdir(flagPath, 's');
end

end
