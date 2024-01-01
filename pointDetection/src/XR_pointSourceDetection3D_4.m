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
% Xiongtao Ruan 09/09/2019 loose condition for initial detection 
% combine original method and energy detection method
% 10/07/2019 estimate snr to adatively select threshold for the mask. 
% Xiongtao Ruan 10/09/2019 also add local maximum points if it is excluded.


function [pstruct, mask, imgLM, imgLoG] = XR_pointSourceDetection3D_4(vol, sigma, varargin)

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
ip.addParameter('RedundancyRadius', 0.25, @(x) numel(x) <= 2 && all(x > 0));
ip.addParameter('RefineMaskLoG', false, @islogical);
ip.addParameter('RefineMaskValid', false, @islogical);
ip.addParameter('ConfRadius', []); % Default: 2*sigma, see fitGaussians3D.
ip.addParameter('WindowSize', []); % Default: 2*sigma, see fitGaussians3D.
ip.addParameter('FitGaussianMethod', 'CD_Newton', @ischar);
ip.addParameter('snrThresh', 2.7, @isnumeric);
ip.addParameter('FilterByResidualSigma', true, @islogical);
ip.addParameter('ResidualSigmaThresh', 2.9, @isnumeric); % Typically, the camera noise std is ~4,
                                                        % if the residual std is very small, it 
                                                        % suggests the point is most probabily from noise. 
ip.addParameter('InitialTestFactor', 16, @isnumeric); % multiplier for the p-value for initial mask
ip.addParameter('EnergyDetectionTestFactor', 2, @isnumeric); % multiplier for the p-value in the energy detection
                                                        
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
%-------------------------------------------------------------------------------------------
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

% Gaussian kernel (spatial)
[x,y,z] = meshgrid(-ws(1):ws(1),-ws(1):ws(1),-ws(2):ws(2));
% g = exp(-(x.^2+y.^2)/(2*sigma(1)^2)) .* exp(-z.^2/(2*sigma(2)^2));
g = exp(-(x.^2+y.^2)/(2*sigma(1)^2) -z.^2/(2*sigma(2)^2));
n = numel(g);
% g(g < 0.1) = 0;
% n = sum(g > 0, 'all');
gsum = sum(g(:));
g2sum = sum(g(:).^2);

% gauss_kernel = fspecial3('gauss', [ws(1), ws(1), ws(2)] * 2 + 1, [sigma(1), sigma(1), sigma(2)]);
gauss_kernel = g;
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


% xruan 07/25/2019 adapt to the old method, because of no normalization of
% the convolution
% fg = fg * gsum;

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
% clear fu;

% kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1);
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
if isempty(tstat_thresh) || alpha ~= ip.Results.Alpha * ip.Results.InitialTestFactor
    df2_nu = (n - 1) * (2 * (n - 1) ^ 2 * C(1, 1) + kLevel ^ 2 * (n - 3)) ^ 2 / (4 * (n - 1) ^ 4 * C(1, 1) ^ 2 + kLevel ^ 4 * (n - 3) ^ 2);
    alpha = ip.Results.Alpha * ip.Results.InitialTestFactor;
    tstat_thresh = tinv(alpha, df2_nu);
end
mask = -T < tstat_thresh;
% toc


%%%%%%%%%%%%%%%%
% energy based detection
% [pstruct, mask_1, NLL] = pointSourceDetection3DMask_likelihood_ratio_test(vol, sigma);
noise_mu = c_est;
noise_sigma = sigma_res;

NLL = (fu2 - 2 * fu .* noise_mu + n * noise_mu .^ 2) ./ (2 * noise_sigma .^ 2);
% mask = NLL > (mean(NLL(:)) + 5.5 * std(NLL(:)));
% in theory, 2 * NLL ~ chi2(n)
chi2_thrsh = chi2inv(1 - ip.Results.Alpha * ip.Results.EnergyDetectionTestFactor, n);
mask_1 = 2 * NLL > chi2_thrsh;
%%%%%%%%%%%%%%%


% sum_mask = sum(mask, 'all')

mask = mask | mask_1;
% mask = mask_1;

% mask = true(size(mask));

% sum_mask_1 = sum(mask_1, 'all')
% sum_comb_mask = sum(mask, 'all')

% clear mask borders (change border conditions for conv3fast to 'zero')
mask([1 2 end-1 end],:,:) = 0;
mask(:,[1 2 end-1 end],:) = 0;
mask(:,:,[1 2 end-1 end]) = 0;

% 09/06/2019 change to use smaller kernel size for local maximum, which is
% better for two points very close or points with low SNR.
% all local max
% allMax = locmax3d(imgLoG, 2*ceil(sigma([1 1 2]))+1, 'ClearBorder', false);
% allMax_orig = locmax3d(imgLoG, [3, 3, 3], 'ClearBorder', false);
allMax_sigma = XR_locmax3d(imgLoG, 2*ceil(sigma([1 1 2]))+1, false);

allMax = true(size(mask));
% local maxima above threshold in image domain
imgLM = allMax .* mask;

% log_kernel_1 = fspecial3('log', [3, 3, 3], [sigma(1), sigma(1), sigma(2)]);
% imgLoG_1 = imfilter(vol, -log_kernel_1, 'conv', 'same', 'symmetric');
% allMax_1 = locmax3d(imgLoG_1, [3, 3, 3], 'ClearBorder', false);

% estimate snr, if snr is high, use strict condition for initial point estimation.
snr_est_inds = allMax_sigma .* mask > 0;
snr_est = A_est(snr_est_inds) ./ sigma_res(snr_est_inds);

median_snr_est = median(snr_est);
is_high_snr_case = median_snr_est > ip.Results.snrThresh;

% for high SNR situation, use the original method
if is_high_snr_case
    kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1);    
    df2_nu = (n - 1) * (2 * (n - 1) ^ 2 * C(1, 1) + kLevel ^ 2 * (n - 3)) ^ 2 / (4 * (n - 1) ^ 4 * C(1, 1) ^ 2 + kLevel ^ 4 * (n - 3) ^ 2);    
    tstat_thresh_high = tinv(ip.Results.Alpha, df2_nu);    
    
    SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
    % df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
    scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
    T = (A_est - sigma_res*kLevel) ./ scomb;
    
    mask = -T < tstat_thresh_high;
    
    mask([1 2 end-1 end],:,:) = 0;
    mask(:,[1 2 end-1 end],:) = 0;
    mask(:,:,[1 2 end-1 end]) = 0;

    allMax = allMax_sigma;
    imgLM = allMax .* mask;
end


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
        
    % apply exclusion mask
    if ~isempty(ip.Results.Mask)
        imgLM(ip.Results.Mask==0) = 0;
    end
    
%     lmIdx = find(imgLM~=0);
    % [lmy,lmx,lmz] = ind2sub(size(vol), lmIdx);

    if ~false && ~is_high_snr_case
        % 09/08/2019 enlarge local maximum mask
        allMax_1 = allMax;
        imgLM = imdilate(allMax_1 ~= 0, strel('sphere', 1)) .* mask;
        if ~isempty(ip.Results.Mask)
            imgLM(ip.Results.Mask==0) = 0;
        end
        
        % 10/09/2019
        allMax_orig = XR_locmax3d(imgLoG, [3, 3, 3], false);
        localMax_inds = find(allMax_orig .* imgLM ~= 0);
        [ly, lx, lz] = ind2sub(size(allMax_orig), localMax_inds); 
        to_include_flag = false(numel(localMax_inds), 1);
        
        % 09/09/2019 xruan only keep one or limited number of points in a 
        % single connected component
        CC = bwconncomp(imgLM ~= 0);
        lmy = zeros(CC.NumObjects, 1);
        lmx = zeros(CC.NumObjects, 1);
        lmz = zeros(CC.NumObjects, 1);
        
        for ind = 1 : CC.NumObjects
            cur_PixelIdx = CC.PixelIdxList{ind};
            [iy, ix, iz] = ind2sub(CC.ImageSize, cur_PixelIdx);
            if numel(cur_PixelIdx) > 1
                iy = mean(iy);
                ix = mean(ix);
                iz = mean(iz);
            end
            
            vol_thrsh = 4;
            if numel(cur_PixelIdx) > vol_thrsh
                cur_lm_ind = find(ismember(localMax_inds, cur_PixelIdx));
                % if there is one local maximum, check whether it is close
                % to the mean point, if not, include them.
                if any(cur_lm_ind)
                    sq_dist_mat = (lx(cur_lm_ind) - ix) .^ 2 + (ly(cur_lm_ind) - iy) .^ 2 + (lz(cur_lm_ind) - iz) .^ 2;
                    max_ind = find(sq_dist_mat > 2 ^ 2);
                    if any(max_ind)
                        to_include_flag(cur_lm_ind(max_ind)) = true;
                    end
                end
            end
            lmy(ind) = iy;
            lmx(ind) = ix;
            lmz(ind) = iz;
        end
        lmIdx = sub2ind(size(vol), round(lmy), round(lmx), round(lmz));
        lmIdx = unique([lmIdx; localMax_inds(to_include_flag)]);
        
        if numel(lmIdx) > 0 && ip.Results.FilterByResidualSigma
            lmIdx = lmIdx(sigma_res(lmIdx) > ip.Results.ResidualSigmaThresh * 0.9);
        end
        
        [lmy,lmx,lmz] = ind2sub(size(vol), lmIdx);
    else
        lmIdx = find(imgLM~=0);
        [lmy,lmx,lmz] = ind2sub(size(vol), lmIdx);        
    end
    
    clear allMax;
    
    % if ~isempty(lmIdx)
    if ~isempty(lmy)
        % run localization on local maxima
        if ~ip.Results.FitMixtures
            pstruct = XR_fitGaussians3D(vol, [lmx lmy lmz], A_est(lmIdx), sigma,...
                c_est(lmIdx), ip.Results.Mode, 'mask', mask, 'alpha', ip.Results.Alpha,...
                'ConfRadius', ip.Results.ConfRadius, 'WindowSize', ws, ...
                'FitGaussianMethod', ip.Results.FitGaussianMethod, 'kLevel', kLevel);
            
        else
            pstruct = XR_fitGaussianMixtures3D(vol, [lmx lmy lmz], A_est(lmIdx), sigma,...
               c_est(lmIdx), 'mask', mask, 'alpha', ip.Results.Alpha,...
               'ConfRadius', ip.Results.ConfRadius, 'WindowSize', ws,...
               'maxM', ip.Results.MaxMixtures, 'kLevel', kLevel);
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
                % xruan: change to allow radius for xy and z
                RedundancyRadius = ip.Results.RedundancyRadius;
                if numel(RedundancyRadius) == 1
                    idxKD = KDTreeBallQuery(pM, pM, RedundancyRadius*ones(numel(pstruct.x),1));
                elseif numel(RedundancyRadius) == 2
                    idxKD_xy = KDTreeBallQuery(pM(:, 1 : 2), pM(:, 1 : 2), RedundancyRadius(1)*ones(numel(pstruct.x),1));
                    idxKD_z = KDTreeBallQuery(pM(:, 3), pM(:, 3), RedundancyRadius(2)*ones(numel(pstruct.x),1));
                    
                    idxKD = cellfun(@(x, y) intersect(x, y), idxKD_xy, idxKD_z, 'uniformoutput', false);
                end
                
                idxKD = idxKD(cellfun(@numel, idxKD)>1);
                for k = 1:length(idxKD)
                    RSS = pstruct.RSS(idxKD{k});
                    idx(idxKD{k}(RSS ~= min(RSS))) = 0;
                end
            end
            
            % xruan 10/10/2019 add filtering based on residual sigma
            if sum(idx) > 0 && ip.Results.FilterByResidualSigma
                idx = idx & (pstruct.sigma_r > ip.Results.ResidualSigmaThresh);
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
