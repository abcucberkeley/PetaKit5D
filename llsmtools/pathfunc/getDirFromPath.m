%[dirName] = getDirFromPath(dpath) returns the last directory contained in the input path

% Francois Aguet, November 2010

function [dirName, dirPath] = getDirFromPath(dpath)

idx = regexp(dpath, filesep);
if idx(end) == length(dpath)
    dirName = dpath(idx(end-1)+1:end-1);
    dirPath = dpath(1:idx(end-1));
else
    dirName = dpath(idx(end)+1:end);
    dirPath = dpath(1:idx(end));
end