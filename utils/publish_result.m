function [published] = publish_result(tmpFullpath, resultFullpath, flagPath)
% move the temporary result tmpFullpath (a folder, e.g. a zarr) of a task
% into place as resultFullpath, unless another writer of the same result
% already put its result there (see keep_other_result: that result is kept
% and this one removed), and remove the task's flag folder flagPath
% (optional). If the other writer's result came between the check and the
% move, movefile put this folder into it: it is removed from there.

if nargin < 3
    flagPath = '';
end

published = false;
if keep_other_result(resultFullpath, tmpFullpath, flagPath)
    return;
end
if exist(tmpFullpath, 'dir')
    movefile(tmpFullpath, resultFullpath);
    [~, fsn, ext] = fileparts(strip(tmpFullpath, 'right', '/'));
    nestedPath = sprintf('%s/%s%s', resultFullpath, fsn, ext);
    published = ~exist(nestedPath, 'dir');
    if ~published
        fprintf('The result %s was made by another job, keep it.\n', resultFullpath);
        rmdir(nestedPath, 's');
    end
end
if ~isempty(flagPath) && exist(flagPath, 'dir')
    rmdir(flagPath, 's');
end

end
