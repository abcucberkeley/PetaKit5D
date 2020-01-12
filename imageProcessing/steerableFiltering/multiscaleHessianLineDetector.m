%MULTISCALEHESSIANLINEDETECTOR line / curvilinear object detection by hessian filtering
%
%[response, theta, nms, scaleindex] = multiscaleHessianLineDetector(input, sigmaVect)
%[response, theta, nms, scaleindex,eigenValues,eigenVectors] = multiscaleHessianLineDetector(...)
%[response, theta, nms, scaleindex] = multiscaleHessianLineDetector(input, sigmaVect,'ScaleSelectionValue','Value')
%
%   Calculates local image hessian at several scales by filtering with
%   gaussian second partial derivatives. Returns optimal values across
%   scales, and performs non-maximum suppression on the selected multiscale
%   response.
%
%       Input:
%
%           input - MxN image matrix
%
%           sigmaVect - 1xP vector of sigmas specifying scales to filter at
%
%           'ScaleSelectionValue' - specifies which value to use to
%           determine which scale has 'optimal' response. Options are:
%
%               'Response' {default} - scale is selected based on the
%               response magnitude (largest scale-normalized eigenvalue of hessian)
%
%               'Anisotropy' - scale is selected based on the maximum
%               anisotropy (ratio of largest to smalles hessian
%               eigenvalues, with negative values clipped to zero)
%
%           'NonMaximumSuppression' - true/false. If true {default true},
%           non-maximum suppression will be run.
%
%           'NormalizeEigenvalues' - true/false. If true {default false},
%           the output eigenvalues will be normalized by the gaussian-filtered
%           image intensity.
%
%       Output:
%
%           MxN matrices with the response, orientation, non-maximum
%           suppressed response, the index of the selected scale, and the
%           eigenvalues and eigenvectors of hessian. The scale,
%           eigenvalues and eigenvectors are scale-normalized (eigenvectors
%           are not unit length)
%
% Francois Aguet, Oct. 13, 2011
% Revised Hunter Elliott, Nov 2013;

function [maxResponse, maxTheta, nms, scaleIndex,maxEigVal,maxEigVec] = multiscaleHessianLineDetector(input, sigmaVect,varargin)

ip = inputParser;
ip.addParamValue('ScaleSelectionValue','Response',@(x)(ismember(x,{'Response','Anisotropy'})));
ip.addParamValue('NonMaximumSuppression',true,@islogical);
ip.addParamValue('NormalizeEigenvalues',false,@islogical);
ip.parse(varargin{:});
p = ip.Results;

[ny,nx] = size(input);
ns = numel(sigmaVect);


for si = 1:ns
    s = sigmaVect(si);

    w = ceil(4*s);
    x = -w:w;
    
    % 1-D components required for filtering
    g = exp(-x.^2/(2*s^2)) / (sqrt(2*pi)*s);
    gx = -x/s^2 .* g;
    gxx = x.^2 .* g / s^4; % -1/s^2 term subtracted below
    
    % compute 3 basis templates
    inputXT = padarray(input, [w w], 'symmetric');
    f_blur = conv2(g, g, inputXT, 'valid') / s^2; % col, row kernel
    f_xx = conv2(g, gxx, inputXT, 'valid') - f_blur;
    f_xy = conv2(gx, gx, inputXT, 'valid');
    f_yy = conv2(gxx, g, inputXT, 'valid') - f_blur;
    
    % eigenvalues -> response
    
    %Quadratic solution to eigenvalue problem for 2x2 symmetric matrix of H:
    % lambda1/2 = (f_xx + f_yy +/- sqrt((f_xx - f_yy) .^2 + 4*f_xy .^2)) ./ 2;
    %             ^----------^     ^------------------------------------^
    %                Alpha                   beta
    alpha = (f_xx + f_yy)/2;    
    beta = sqrt((f_xx - f_yy) .^2 + 4*f_xy .^2)/2; 
    eigVal(:,:,1) = -alpha - beta;%Flip sign because we want eigenvalues of -H
    eigVal(:,:,2) = -alpha + beta;            
    
    %Get non-unit eigenvectors - we only use the direction    
    eigVec(:,:,1,1) = eigVal(:,:,1) + f_yy;
    eigVec(:,:,2,1) = -f_xy;
            
    eigVec(:,:,1,2) = eigVal(:,:,2) + f_yy;
    eigVec(:,:,2,2) = -f_xy;
    
    %Scale-normalize eigenvalues and vectors
    eigVal = eigVal .* s^2;    
    eigVec = eigVec .* s^2;        
            
    %Second eigenvalue/vector will always be largest and is the response
    response = eigVal(:,:,2);
    theta = atan(eigVec(:,:,2,2) ./ eigVec(:,:,1,2));            
    
    %Set value used for scale selection
    switch p.ScaleSelectionValue        
        case 'Response'            
            ssVal = response;            
        case 'Anisotropy'            
            ssVal = eigVal;
            ssVal(ssVal<0) = 0;%Suppress valleys and the negative curvature component of saddles             
            ssVal = ssVal(:,:,2) - ssVal(:,:,1);
    end
        
    if p.NormalizeEigenvalues
        eigVal = eigVal ./ repmat(f_blur .* s^2,[1 1 2]);
    end
    
    if si == 1
        maxSsVal = ssVal;%Store scale selection value
        maxResponse = response;
        maxTheta = theta;        
        scaleIndex = ones(ny,nx,minIntClass(max(ns,2)));        
        maxEigVal = eigVal;
        maxEigVec = eigVec;        
    else
        idx = ssVal > maxSsVal;%Select pixels with larger scale selection value, store outputs        
        maxSsVal(idx) = ssVal(idx);
        maxResponse(idx) = response(idx);
        maxTheta(idx) = theta(idx);    
        scaleIndex(idx) = si;
        idx2 = repmat(idx,[1 1 2]);        
        maxEigVal(idx2) = eigVal(idx2);
        idx2 = repmat(idx,[1 1 2 2]);
        maxEigVec(idx2) = eigVec(idx2);        
    end
end

if p.NonMaximumSuppression && nargout > 2
    nms = nonMaximumSuppression(maxResponse, maxTheta);
else
    nms = [];
end
    
