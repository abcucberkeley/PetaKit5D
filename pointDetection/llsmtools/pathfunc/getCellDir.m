function [dname, dpath] = getCellDir(data)

chParents = cellfun(@(c) getParentDir(c), data.channels, 'unif', 0);
sParent = getParentDir(data.source);

v = strcmp(sParent, chParents);

if numel(data.channels)==1 || ~all(v)
    dpath = data.source;
else % all channels at same level, below 'source'
    dpath = sParent;
end
dname = getDirFromPath(dpath);
dpath = getParentDir(dpath);