function m = bwLargestObj(m,conn)
%BWLARGESTOBJ removes all but the largest object (connected component) in the input binary mask
%
% bw = bwLargestObj(m)
% bw = bwLargestObj(m,conn)
%
%  m - Input binary matrix
%  conn - neighborhood connectivity (4,8 etc) to use. Optional. Default is
%  bwconncomp default.
%

% Hunter Elliott
% 3/2013


if nargin < 2 || isempty(conn)
    cc=bwconncomp(m);
else    
    cc=bwconncomp(m,conn);
end
[~,iMax] = max(cellfun(@numel,cc.PixelIdxList));
m(:) = false;
m(cc.PixelIdxList{iMax}) = true;

