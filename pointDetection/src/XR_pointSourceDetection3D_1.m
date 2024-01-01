%[pstruct, mask, imgLM, imgLoG] = pointSourceDetection3D(vol, sigma, varargin)
%
% Inputs :   
%                 vol : input volume
%               sigma : standard deviation of the Gaussian PSF
%                       If the PSF is anisotropic, 'sigma' should be a two-element vector:
%                       [sigma_xy sigma_z]
%
% Options (as 'specifier'-value pairs): 
%
%              'mode' : parameters to estimate. Default: 'xyzAc'.
%             'alpha' : alpha value used in the statistical tests. Default: 0.05.
%              'mask' : mask of pixels (i.e., cell mask) to include in the detection. Default: all.
%       'FitMixtures' : true|{false}. Toggles mixture-model fitting.
%       'MaxMixtures' : maximum number of mixtures to fit. Default: 5.
%   'RemoveRedundant' : {true}|false. Discard localizations that coincide within 'RedundancyRadius'.
%  'RedundancyRadius' : Radius for filtering out redundant localizatios. Default: 0.25
%     'RefineMaskLoG' : true|{false}. Apply threshold to LoG-filtered img to refine mask of significant pixels.
%   'RefineMaskValid' : true|{false}. Return only mask regions where a significant signal was localized.
%        'ConfRadius' : Confidence radius for positions, beyond which the fit is rejected. Default: 2*sigma
%        'WindowSize' : Window size for the fit. Default: 2*sigma, i.e., [-2*sigma ... 2*sigma]^2
%
% Outputs:  
%             pstruct : output structure with Gaussian parameters, standard deviations, p-values
%                mask : mask of significant (in amplitude) pixels
%               imgLM : image of local maxima
%              imgLoG : Laplacian of Gaussian-filtered image

% Francois Aguet, August 2013 (last modified: 11/05/2013)
% Xiongtao Ruan, July 2019 use matlab functions for convolution, which
% turns out to be much faster. 
% Xiongtao Ruan 09/09/2019 loose condition for initial detection using
% smaller region for local maximum

function [pstruct, mask, imgLM, imgLoG] = pointSourceDetection3D(vol, sigma, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('vol', @isnumeric);
ip.addRequired('sigma', @isnumeric);
ip.addParameter('Mode', 'xyzAc', @ischar);
ip.addParameter('Alpha', 0.05, @isscalar);
ip.addParameter('Mask', [], @(x) isnumeric(x) || islogical(x));
ip.addParameter('FitMixtures', false, @islogical);
ip.addParameter('MaxMixtures', 5, @(x) numel(x)==1 && x>0 && round(x)==x);
ip.addParameter('RemoveRedundant', true, @islogical);
ip.addParameter('RedundancyRadius', 0.25, @isscalar);
ip.addParameter('RefineMaskLoG', false, @islogical);
ip.addParameter('RefineMaskValid', false, @islogical);
ip.addParameter('ConfRadius', []); % Default: 2*sigma, see fitGaussians3D.
ip.addParameter('WindowSize', []); % Default: 2*sigma, see fitGaussians3D.
ip.addParameter('FitGaussianMethod', 'CD_Newton', @ischar);
ip.parse(vol, sigma, varargin{:});

if ~isa(vol, 'double')
    vol = double(vol);
end

if numel(sigma)==1
    sigma = [sigma sigma];
end

ws = ip.Results.WindowSize;
if isempty(ws)
    ws = ceil(2*sigma);
elseif numel(ws)==1
    ws = [ws ws];
end

%-------------------------------------------------------------------------------------------
% Convolutions
%-------------------------------XR_pointSourceDetection3D_1------------------------------------------------------------
% % right-hand side of symmetric kernels
% gx = exp(-(0:ws(1)).^2/(2*sigma(1)^2));
% gz = exp(-(0:ws(2)).^2/(2*sigma(2)^2));
% tic
% fg = conv3fast(vol, gx, gx, gz);
% toc
% tic
% fg = imgaussfilt3(vol, [sigma(1), sigma(1), sigma(2)], 'padding', 'symmetric');
% toc
% tic
% fg = imgaussfilt3(vol, [sigma(1), sigma(1), sigma(2)], 'FilterDomain', 'spatial', 'padding', 'symmetric');
% toc
% tic
gauss_kernel = fspecial3('gauss', [ws(1), ws(1), ws(2)] * 2 + 1, [sigma(1), sigma(1), sigma(2)]);
fg = imfilter(vol, gauss_kernel, 'conv', 'same', 'symmetric');
% tocoriginal


% tic
% fu =  conv3fast(vol,    ones(1,ws(1)+1), ones(1,ws(1)+1), ones(1,ws(2)+1));
% fu2 = conv3fast(vol.^2, ones(1,ws(1)+1), ones(1,ws(1)+1), ones(1,ws(2)+1));
% toc
% tic
% mean_kernel = fspecial3('average', ws(1) * 2 + 1);
% fu_1 =  convn(vol, mean_kernel, 'same');
% fu2_1 = convn(vol.^2, mean_kernel, 'same');
% toc
% tic
sum_kernel = ones(ws(1) * 2 + 1, ws(1) * 2 + 1, ws(2) * 2 + 1);
fu =  imfilter(vol, sum_kernel, 'conv', 'same', 'symmetric');
fu2 = imfilter(vol.^2, sum_kernel, 'conv', 'same', 'symmetric');
% toc

% tic
% fu_3 =  imboxfilt3(vol, [ws(1) ws(1) ws(2)] * 2 + 1, 'padding', 'symmetric');
% fu2_3 = imboxfilt3(vol.^2, [ws(1) ws(1) ws(2)] * 2 + 1, 'padding', 'symmetric');
% toc

% Laplacian of Gaussian-filtered input
% tic
% gx2 = (0:ws(1)).^2 .*gx;
% gz2 = (0:ws(2)).^2 .*gz;
% fgx2 = conv3fast(vol, gx2, gx, gz);
% fgy2 = conv3fast(vol, gx, gx2, gz);
% fgz2 = conv3fast(vol, gx, gx, gz2);
% imgLoG = (2/sigma(1)^2+1/sigma(2)^2)*fg - ((fgx2+fgy2)/sigma(1)^4 + fgz2/sigma(2)^4);
% 
% clear fgx2 fgy2 fgz2;
% toc

% tic 
log_kernel = fspecial3('log', [ws(1) ws(1) ws(2)] * 2 + 1, [sigma(1), sigma(1), sigma(2)]);
imgLoG = imfilter(vol, -log_kernel, 'conv', 'same', 'symmetric');
% toc

% Gaussian kernel (spatial)
[x,y,z] = meshgrid(-ws(1):ws(1),-ws(1):ws(1),-ws(2):ws(2));
% g = exp(-(x.^2+y.^2)/(2*sigma(1)^2)) .* exp(-z.^2/(2*sigma(2)^2));
g = exp(-(x.^2+y.^2)/(2*sigma(1)^2) -z.^2/(2*sigma(2)^2));
n = numel(g);
gsum = sum(g(:));
g2sum = sum(g(:).^2);

% xruan 07/25/2019 adapt to the old method, because of no normalization of
% the convolution
fg = fg * gsum;

% solution to linear system
A_est = (fg - gsum*fu/n) / (g2sum - gsum^2/n);
c_est = (fu - A_est*gsum)/n;

J = [g(:) ones(n,1)]; % g_dA g_dc
C = inv(J'*J);

f_c = fu2 - 2*c_est.*fu + n*c_est.^2; % f-c
RSS = A_est.^2*g2sum - 2*A_est.*(fg - c_est*gsum) + f_c;
% clear fg fu2;
RSS(RSS<0) = 0; % negative numbers may result from machine epsilon/roundoff precision
sigma_e2 = RSS/(n-3);

sigma_A = sqrt(sigma_e2*C(1,1));

% standard deviation of residuals
sigma_res = sqrt(RSS/(n-1));
clear fu;

kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1);
kLevel = norminv(1-ip.Results.Alpha, 0, 1);

SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
% df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
T = (A_est - sigma_res*kLevel) ./ scomb;

% mask of admissible positions for local maxima
% tic
% mask = tcdf(-T, df2) < 0.05;
% toc

% xruan 07/26/2019 in theory, df2 is constant, and it is possible only compute
% once to obtain the mask without perform T-test over each pixel
% tic
% tstat_thresh = tinv(0.05, df2_nu);
% save tstat_thresh to memory
persistent alpha tstat_thresh
if isempty(tstat_thresh) || alpha ~= ip.Results.Alpha * 16
    df2_nu = (n - 1) * (2 * (n - 1) ^ 2 * C(1, 1) + kLevel ^ 2 * (n - 3)) ^ 2 / (4 * (n - 1) ^ 4 * C(1, 1) ^ 2 + kLevel ^ 4 * (n - 3) ^ 2);
    alpha = ip.Results.Alpha * 16;
    tstat_thresh = tinv(alpha, df2_nu);
end
mask = -T < tstat_thresh;
% toc

% clear mask borders (change border conditions for conv3fast to 'zero')
mask([1 2 end-1 end],:,:) = 0;
mask(:,[1 2 end-1 end],:) = 0;
mask(:,:,[1 2 end-1 end]) = 0;

% 09/06/2019 change to use smaller kernel size for local maximum, which is
% better for two points very close or points with low SNR.
% all local max
% allMax = locmax3d(imgLoG, 2*ceil(sigma([1 1 2]))+1, 'ClearBorder', false);
allMax = locmax3d(imgLoG, [3, 3, 3], 'ClearBorder', false);

% local maxima above threshold in image domain
imgLM = allMax .* mask;

% log_kernel_1 = fspecial3('log', [3, 3, 3], [sigma(1), sigma(1), sigma(2)]);
% imgLoG_1 = imfilter(vol, -log_kernel_1, 'conv', 'same', 'symmetric');
% allMax_1 = locmax3d(imgLoG_1, [3, 3, 3], 'ClearBorder', false);

pstruct = [];
if sum(imgLM(:))~=0 % no local maxima found, likely a background image
    
    if ip.Results.RefineMaskLoG
        % -> set threshold in LoG domain
        logThreshold = min(imgLoG(imgLM~=0));
        logMask = imgLoG >= logThreshold;
        
        % combine masks
        mask = mask | logMask;
        % re-select local maxima
        imgLM = allMax .* mask;
    end
    
    allMax_1 = allMax;
    imgLM = imdilate(allMax_1 ~= 0, strel('sphere', 1)) .* mask;
    clear allMax;
    
    % apply exclusion mask
    if ~isempty(ip.Results.Mask)
        imgLM(ip.Results.Mask==0) = 0;
    end
    
    lmIdx = find(imgLM~=0);
    [lmy,lmx,lmz] = ind2sub(size(vol), lmIdx);
    
    if ~isempty(lmIdx)
        % run localization on local maxima
        if ~ip.Results.FitMixtures
            pstruct = XR_fitGaussians3D(vol, [lmx lmy lmz], A_est(lmIdx), sigma,...
                c_est(lmIdx), ip.Results.Mode, 'mask', mask, 'alpha', ip.Results.Alpha,...
                'ConfRadius', ip.Results.ConfRadius, 'WindowSize', ws, ...
                'FitGaussianMethod', ip.Results.FitGaussianMethod);
        else
            pstruct = fitGaussianMixtures3D(vol, [lmx lmy lmz], A_est(lmIdx), sigma,...
               c_est(lmIdx), 'mask', mask, 'alpha', ip.Results.Alpha,...
               'ConfRadius', ip.Results.ConfRadius, 'WindowSize', ws,...
               'maxM', ip.Results.MaxMixtures);
        end
    
        % remove NaN values
        idx = ~isnan([pstruct.x]);
        if sum(idx)~=0
            fnames = fieldnames(pstruct);           
            for k = 1:length(fnames)
                pstruct.(fnames{k}) = pstruct.(fnames{k})(:,idx);
            end
            
            % significant amplitudes
            idx = [pstruct.hval_Ar] == 1;
            
            % eliminate duplicate positions (resulting from localization)
            if ip.Results.RemoveRedundant
                pM = [pstruct.x' pstruct.y' pstruct.z'];
                idxKD = KDTreeBallQuery(pM, pM, ip.Results.RedundancyRadius*ones(numel(pstruct.x),1));
                idxKD = idxKD(cellfun(@numel, idxKD)>1);
                for k = 1:length(idxKD);
                    RSS = pstruct.RSS(idxKD{k});
                    idx(idxKD{k}(RSS ~= min(RSS))) = 0;
                end
            end
            
            if sum(idx)>0
                fnames = fieldnames(pstruct);
                for k = 1:length(fnames)
                    pstruct.(fnames{k}) = pstruct.(fnames{k})(:,idx);
                end
                pstruct.hval_Ar = logical(pstruct.hval_Ar);
                pstruct.hval_AD = logical(pstruct.hval_AD);
                pstruct.isPSF = ~pstruct.hval_AD;
                
                % adjust mixture index if only a single component remains
                if ip.Results.FitMixtures
                    mv = 0:max(pstruct.mixtureIndex);
                    multiplicity = getMultiplicity(pstruct.mixtureIndex);
                    pstruct.mixtureIndex(ismember(pstruct.mixtureIndex, mv(multiplicity==1))) = 0;
                end
            else
                pstruct = [];
            end
            else
            pstruct = [];
        end
    end
end

if ~isempty(pstruct) && ip.Results.RefineMaskValid
    CC = bwconncomp(mask);
    labels = labelmatrix(CC);
    loclabels = labels(sub2ind(size(vol), pstruct.y_init, pstruct.x_init));
    idx = setdiff(1:CC.NumObjects, loclabels);
    CC.PixelIdxList(idx) = [];
    CC.NumObjects = length(CC.PixelIdxList);
    mask = labelmatrix(CC)~=0;
end
