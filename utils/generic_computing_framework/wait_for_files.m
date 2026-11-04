function [done] = wait_for_files(fileGroups, waitTime)
% waits waitTime seconds, or until all files of one of the groups exist
% (then done is true). The callers of single jobs (XR_microscopeAutomaticProcessing
% and the decon and stitching wrappers) wait so between their rounds for the
% results that their jobs still make and for the files that tasks wait for,
% so that the next round follows them within about a second.
%
%   fileGroups : a cell array of groups, each a cell array of paths of files
%                or folders; a path whose name has * is a pattern (as for
%                dir). Empty groups are left out; without groups it is
%                pause(waitTime).
%   waitTime   : the longest wait in seconds.
%
% It looks for the files it hasn't found yet after each second of the wait,
% and less often when a look takes longer than a tenth of a second (many
% files, slow storage), so that looking takes at most a tenth of the wait.
% Callers pass files whose existence their next round acts on (it marks a
% result done or submits a task that waited for a file), so that a wait that
% ends early is followed by progress.


done = false;
if isempty(fileGroups)
    fileGroups = {};
end
fileGroups = fileGroups(~cellfun(@isempty, fileGroups));
if isempty(fileGroups)
    pause(waitTime);
    return;
end

% (the group of each path, a column also for one group)
paths = cellfun(@(x) reshape(x, [], 1), fileGroups(:), 'unif', 0);
groupInds = repelem((1 : numel(paths))', cellfun(@numel, paths));
groupInds = groupInds(:);
[paths, ~, pathInds] = unique(cat(1, paths{:}));
found = false(numel(paths), 1);
% patterns: a listing of their folder per look (names only), and dir for a
% name that matches, as the callers look for them
[pstrs, fsns, exts] = fileparts(paths);
names = strcat(fsns, exts);
isPattern = contains(names, '*');

interval = 1;
waited = 0;
while waitTime - waited > 1e-6
    t = min(interval, waitTime - waited);
    pause(t);
    waited = waited + t;
    % (the caller's next round looks right after the wait)
    if waitTime - waited <= 1e-6
        break;
    end

    t0 = tic;
    % batch_file_exist sees files another node just made, which exist()
    % can miss for 30 s or more; so does a listing of their folder
    inds = find(~found & ~isPattern);
    if ~isempty(inds)
        found(inds) = batch_file_exist(paths(inds));
    end
    inds = find(~found & isPattern);
    [folders, ~, folderInds] = unique(pstrs(inds));
    for d = 1 : numel(folders)
        listed = list_names(folders{d});
        for i = reshape(inds(folderInds == d), 1, [])
            found(i) = any(name_matches(listed, names{i})) && ~isempty(dir(paths{i}));
        end
    end
    if any(accumarray(groupInds, double(~found(pathInds)), [numel(fileGroups), 1]) == 0)
        done = true;
        return;
    end
    interval = max(1, 10 * toc(t0));
end

end


function [names] = list_names(folderPath)
% the names in a folder (with Java, without the file attributes dir reads)

names = {};
if usejava('jvm')
    % null (not listed) comes back as []
    fn_arr = java.io.File(folderPath).list();
    if isjava(fn_arr)
        names = cell(fn_arr);
    end
else
    dir_info = dir(folderPath);
    names = {dir_info.name};
end

end


function [tf] = name_matches(names, pattern)
% names that match pattern (* for any characters)

parts = strsplit(pattern, '*');
if numel(parts) == 2
    tf = startsWith(names, parts{1}) & endsWith(names, parts{2}) ...
        & cellfun(@numel, names) >= numel(parts{1}) + numel(parts{2});
else
    tf = ~cellfun(@isempty, regexp(names, ['^', regexptranslate('wildcard', pattern), '$'], 'once'));
end

end
