function distMat = bwMaxDirectDist(mask)
%BWMAXDIRECTDIST calculates the distance in each direction to a pixel in the mask
% 
% distMat = bwMaxDirectDist(mask)
% 
% This function calculates, for each 0/false pixel in the input mask, the
% distance to a 1/true pixel in a particular direction. This is analagous
% to a distance transform in that the minimum value over all 4 directions is
% the same as a city-block distance transform, while the maximum is a sort
% of maximum-distance transform.
% 
% 
% Input:
% 
%   mask - MxN 2D logical matrix
%
% Output:
%
%  distMat - MxNx4 single-precision matrix, where the 3rd dimension
%  corresponds to different directions. That is:
%       distMat(:,:,1) gives the distance to the nearest true pixel in
%       direction of increasing index in the first dimension,
%       distMat(:,:,2) gives the distance to the nearest true pixel in
%       direction of decreasing index in the first dimension,
%       distMat(:,:,3) gives the distance to the nearest true pixel in
%       direction of increasing index in the second dimension,
%       distMat(:,:,4) gives the distance to the nearest true pixel in
%       direction of decreasing index in the second dimension.
%
% Examples: 
% 
%   These commands would get the directional distance of a mask of the
%   peaks function, and show that the maximum along the 3rd dimension is
%   the same as the distance transform using the cityblock metric:
%   
%   distMat = bwMaxDirectDist(peaks>.5);
%   imagesc(min(distMat,[],3)),figure,imagesc(bwdist(~mask,'cityblock'))
% 
% 
% Hunter Elliott
% 9/2011 
%


[M,N] = size(mask);

distMat = zeros([size(mask) 4],'single');

distMat(1,:,1) = M;
distMat(M,:,2) = M;
distMat(:,1,3) = N;
distMat(:,N,4) = N;

for m = 2:M       
    distMat(m,~mask(m,:),1) = distMat(m-1,~mask(m,:),1) + 1;       
end
for m = M-1:-1:1   
    distMat(m,~mask(m,:),2) = distMat(m+1,~mask(m,:),2) + 1;       
end
for n = 2:N   
    distMat(~mask(:,n),n,3) = distMat(~mask(:,n),n-1,3) + 1;       
end
for n = N-1:-1:1   
    distMat(~mask(:,n),n,4) = distMat(~mask(:,n),n+1,4) + 1;       
end

%YOU NEED TO SET PIXELS WITH NO MASK OBJECT IN THEIR DIRECTION TO SOME
%OTHER VALUE THAN WHAT THEY ARE NOW, OTHERWISE THIS OUTPUT IS A PAIN TO
%USE...
