%SKELETON3D skeletonizes a 3D binary image
%
% skel = skeleton3D(mask)
%
% This function produces a 26-connected skeleton from the input 3D binary
% mask. The skeleton is produced by the iterative thinning algortith found
% in [1]. The function was converted to a MEX function from the C++ code
% found in the supplement of [2].
%
% Input: 
%
%   mask - A 3D binary (logical) matrix.
%
% Output:
%
%   skel - A 3D binary matrix containing the skeleton.
%
%
%
% References:
%
% [1] K. Palagyi & A. Kuba, A Parallel 3D 12-Subiteration Thinning
% Algorithm, Graphical Models and Image Processing, 61, p. 199, (1999)
%
% [2] N. Cornea et al, Curve Skeleton Properties, Applications and Algorithms, 
% IEEE Trans. on Vis. and Comp. Graphics, Vol 13, No 3 (2007)
%
%
% Compile Commands:
% 
%   To compile this function, simply change to the directory containing
%   both skeleton3D.h and skeleton3D.cpp and type this command at the
%   matlab prompt:
%
%       mex skeleton3D.cpp
%
% Hunter Elliott
% 6/2010
%