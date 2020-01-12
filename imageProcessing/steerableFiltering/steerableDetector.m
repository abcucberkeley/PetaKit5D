%[res, theta, nms, rotations] = steerableDetector(img, M, sigma) performs edge/ridge detection through a generalization of Canny's alorithm based on steerable filters
%
% Inputs: 
%         img : input image
%           M : order of the filter, between 1 and 5
%             : Odd orders: edge detectors, M = 1 is equivalent to Canny's detector
%             : Even orders: ridge detectors
%               Higher orders provide better orientation selectivity and are less sensitive to noise,
%               at a small trade-off in computational cost.
%       sigma : standard deviation of the Gaussian kernel on which the filters are based
%   {nAngles} : optional input specifying the number of angles computed for the filterbank output. Default: 36
%{bordercond} : optional character input specifying the border condition
%             : 'mirror' (default), 'replicate', 'periodic', 'zeros'
%
% Outputs: 
%         res : response to the filter
%       theta : orientation map
%         nms : non-maximum-suppressed response
%   rotations : response of the input to rotated versions of the filter, at 'nAngles' different angles.
%
% For more information, see:
% M. Jacob et al., IEEE Trans. Image Proc., 26(8) 1007-1019, 2004.

% Copyright (C) 2011-2012 Francois Aguet
% Adapted from the SteerableJ package, Copyright (C) 2005-2008 Francois Aguet, Biomedical Imaging Group, EPFL.

function [res, theta, nms, filterBank] = steerableDetector(img, M, sigma) %#ok<STOUT,INUSD>