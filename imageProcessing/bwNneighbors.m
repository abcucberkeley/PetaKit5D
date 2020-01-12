function nn = bwNneighbors(bw,nHood)
%BWNNEIGHBORS counts the number of neighbors of each element in the input 2 or 3 dimensional binary matrix
% 
% nn = bwNneighbors(bw)
%
% nn = bwNneighbors(bw,nHood)
%
% This function returns a matrix the same size as the input BW matrix, but
% with each value representing the number of non-zero neighbors that
% element has. This function supports binary arrays of any dimension.
%
% NOTE: This is primarily of use in analyzing binary skeletons,
% where end-points are pixels with exactly 1 neighbor, line points have 2
% nieghbors and junctions have 3 or more neighbors. E.G., the command:
% 
% ep = (bwNneighbors(skel) == 1) & bw
%
% Will return all the end-points in the skeleton 'skel'
%
% Input:
% 
%   bw - A 3D binary matrix.
% 
%   nHood - N-Dimensional binary matrix specifying the neighborhood to use
%   to define neighbors. Optional. Default is ALL adjacent voxels (8
%   connectivity for 2D, 26 connectivity for 3D etc. )
% 
% 
% Output:
% 
%   nn - a uint8 matrix of the same size as bw, but with each element
%   replaced with the number of non-zero neighbors that element has, given
%   input neighborhood nHood
% 
% 
% Hunter Elliott
% 3/2010
%

nD = ndims(bw);

if nD < 2 || nD > 3
    error('Input matrix must be 2 or 3 dimensional!')
end

%Create default neighborhood if not input
if nargin < 2 || isempty(nHood)
    nHood = true(repmat(3,1,nD));
    nHood(ceil(numel(nHood)/2)) = false;
end

%Find the neighbors
iN = find(nHood(:));

%divide it into seperate neighborhoods with 1 neighbor each
iN = ((0:length(iN)-1)' * numel(nHood)) + iN;
nHood = false([size(nHood) nnz(nHood)]);
nHood(iN) = true;

%Initalize the output
nn = zeros(size(bw),'uint8');

%This can probably be sped-up, maybe using accumarray...HLE
%Go through each neighbor and sum up dilations
if ndims(bw) == 3
    for j = 1:nnz(nHood)    
        nn = nn + uint8(imdilate(bw,nHood(:,:,:,j)));                
    end
elseif ndims(bw) == 2
    for j = 1:nnz(nHood)    
        nn = nn + uint8(imdilate(bw,nHood(:,:,j)));                
    end    
end


