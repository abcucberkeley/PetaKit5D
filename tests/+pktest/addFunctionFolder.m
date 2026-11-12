function [] = addFunctionFolder(folder, functionName)
% puts folder first on the path, for the stub or spy functionName.m that a
% fixture just wrote there, and makes sure that MATLAB calls it: until which
% names that file, the path is read again (rehash). On Windows, MATLAB learns
% of the files in its path folders from the file system's change
% notifications, so a caller could still get the function of an earlier
% fixture or the real one for its first calls. (If which names the folder
% another way, e.g. with a short name, it gives up after 10 rehashes, as
% before.)

fn = strrep(fullfile(folder, [functionName, '.m']), '\', '/');
addpath(folder, '-begin');
for i = 1 : 10
    if strcmpi(strrep(which(functionName), '\', '/'), fn)
        return;
    end
    rehash;
end

end
