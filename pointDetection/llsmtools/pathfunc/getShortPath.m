%spath = getShortPath(data) returns the truncated path of a cell directory (3 levels) 

% Francois Aguet, 05/13/2011

function spath = getShortPath(data, level)
if nargin<2
    level = 3;
end

if numel(data.channels)>1
    mCh = strcmp(data.channels, data.source);
    sCh = setdiff(1:length(data.channels),mCh);
    spath = data.channels{sCh(1)};
    idx = regexp(spath, filesep);
    spath = spath(idx(end-level-1)+1:idx(end-1));
else
    spath = data.channels{1};
    idx = regexp(spath, filesep);
    spath = spath(idx(max(1,end-level))+1:idx(end)-1);
end
