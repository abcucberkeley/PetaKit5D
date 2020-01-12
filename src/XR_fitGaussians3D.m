% pStruct = fitGaussians3D(vol, X, A, sigma, c, mode, varargin) fits
% 3D Gaussians to the input volume at the specified locations.
%
% Inputs:   vol : input volume
%             X : Nx3 matrix of initial (or fixed) (x,y,z)-positions
%             A : initial (or fixed) amplitudes
%         sigma : initial (or fixed) Gaussian PSF standard deviation. If z-sampling is different
%                 from x,y sampling, this should be a 2-element vector
%             c : initial (or fixed) local background intensities
%
% Options:
%          mode : selector for optimization parameters, any of 'xyzAsrc'. Default: 'xyzAc'
%                 's' selects the (x,y) s.d., and 'r' the (z) s.d. 
%
% Options ('specifier', value):
%        'Mask' : mask of spot locations
%
% Output: pStruct: structure with fields:
%                  x : estimated x-positions
%                  y : estimated y-positions
%                  z : estimated z-positions
%                  A : estimated amplitudes
%                  s : estimated standard deviations of the PSF
%                      2xN matrix, with the 1st row containing (x,y) s.d. and 
%                      the 2nd row containing the (z) s.d.
%                  c : estimated background intensities
%
%             x_pstd : standard deviations, estimated by error propagation
%             y_pstd : "
%             z_pstd : "
%             A_pstd : "
%             s_pstd : "
%             c_pstd : "
%            sigma_r : standard deviation of the background (residual)
%         SE_sigma_r : standard error of sigma_r
%            pval_Ar : p-value of an amplitude vs. background noise test (p < 0.05 -> significant amplitude)
%
% Usage for a volume with know spot locations (mask) and fixed sigma:
% fitGaussians3D(vol, X, A, sigma, c, 'xyzAc', 'mask', mask);
%
% See also fitGaussian3D

% Francois Aguet, August 2013 (last updated: 10/27/2013)
% Xiongtao Ruan, July 209 improve efficiency
% xruan, Oct, 2019 make KLevel as parameter

function pStruct = XR_fitGaussians3D(vol, X, A, sigma, c, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol', @isnumeric);
ip.addRequired('X');
ip.addRequired('A');
ip.addRequired('sigma');
ip.addRequired('c');
ip.addOptional('mode', 'xyzAc', @ischar);
ip.addParameter('Alpha', 0.05, @isscalar);
ip.addParameter('AlphaT', 0.05, @isscalar);
ip.addParameter('Mask', [], @islogical);
ip.addParameter('ConfRadius', []);
ip.addParameter('WindowSize', []);
ip.addParameter('FitGaussianMethod', 'CD_Newton');
ip.addParameter('kLevel', norminv(1-0.05/2.0, 0, 1), @isnumeric);
ip.parse(vol, X, A, sigma, c, varargin{:});

np = size(X,1);
sigma = ip.Results.sigma;
if numel(sigma)==1
    sigma = sigma*ones(np,2);
elseif numel(sigma)==2
    sigma = repmat(sigma(:)', [np 1]);
end
mode = ip.Results.mode;
if ~isempty(ip.Results.Mask)
    % labels = double(labelmatrix(bwconncomp(ip.Results.Mask)));
    labels = bwlabeln(ip.Results.Mask);
else
    labels = zeros(size(vol));
end

FitGaussianMethod = ip.Results.FitGaussianMethod;

pStruct = struct('x', [], 'y', [], 'z', [], 'A', [], 's', [], 'c', [],...
    'x_pstd', [], 'y_pstd', [], 'z_pstd', [], 'A_pstd', [], 's_pstd', [], 'c_pstd', [],...
    'x_init', [], 'y_init', [], 'z_init', [],...
    'sigma_r', [], 'SE_sigma_r', [], 'RSS', [], 'pval_Ar', [], 'hval_Ar', [], 'hval_AD', []);

[ny,nx,nz] = size(vol);
roundConstr = @(x,N) max(min(round(x),N),1);
xi = roundConstr(X(:,1), nx);
yi = roundConstr(X(:,2), ny);
zi = roundConstr(X(:,3), nz);

% kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1); % ~2 std above background
kLevel = ip.Results.kLevel;

iRange = [min(vol(:)) max(vol(:))];

estIdx = regexpi('xyzAsrc', ['[' mode ']']);


% initialize pStruct arrays
pStruct.x = NaN(1,np);
pStruct.y = NaN(1,np);
pStruct.z = NaN(1,np);
pStruct.A = NaN(1,np);
pStruct.s = NaN(2,np);
pStruct.c = NaN(1,np);

pStruct.x_pstd = NaN(1,np);
pStruct.y_pstd = NaN(1,np);
pStruct.z_pstd = NaN(1,np);
pStruct.A_pstd = NaN(1,np);
pStruct.s_pstd = NaN(2,np);
pStruct.c_pstd = NaN(1,np);

pStruct.x_init = reshape(xi, [1 np]);
pStruct.y_init = reshape(yi, [1 np]);
pStruct.z_init = reshape(zi, [1 np]);

pStruct.sigma_r = NaN(1,np);
pStruct.SE_sigma_r = NaN(1,np);
pStruct.RSS = NaN(1,np);

pStruct.pval_Ar = NaN(1,np);

pStruct.hval_AD = false(1,np);
pStruct.hval_Ar = false(1,np);

% if different sigma values are passed, use largest for filters
sigma_max = max(sigma,[],1);
w2 = ip.Results.ConfRadius;
if isempty(w2)
    w2 = ceil(2*sigma_max);
elseif numel(w2)==1
    w2 = [w2 w2];
end
ws = ip.Results.WindowSize;
if isempty(ws)
    ws = ceil(2*sigma_max);
elseif numel(ws)==1
    ws = [ws ws];
end

% xruan 07/29/2019 precompute the coordinates, so that there is no need to
% compute in each loop
[x_window_mat, y_window_mat, z_window_mat] = meshgrid(0 : ws(1) * 2, 0 : ws(1) * 2, 0 : ws(2) * 2);

T = zeros(1,np);
df2 = zeros(1,np);
res_mat = zeros(np, 3);
for p = 1:np
    
    % window boundaries
    xa = max(1,xi(p)-ws(1)):min(nx,xi(p)+ws(1));
    ya = max(1,yi(p)-ws(1)):min(ny,yi(p)+ws(1));
    za = max(1,zi(p)-ws(2)):min(nz,zi(p)+ws(2));
    
    % relative coordinates of (xi, yi, zi) in window. Origin at (0,0,0)
    ox = xi(p)-xa(1);
    oy = yi(p)-ya(1);
    oz = zi(p)-za(1);
    
    % label mask
    maskWindow = labels(ya, xa, za);
    maskWindow(maskWindow==maskWindow(oy+1,ox+1,oz+1)) = 0;
    
    window = vol(ya, xa, za);
    % window_orig = window;
    % set any other components to NaN
    window(maskWindow~=0) = NaN;
    npx = sum(isfinite(window(:)));    
    
    if npx >= 20 % only perform fit if window contains sufficient data points
        window_size = size(window);
        if any(window_size < [ws(1), ws(1), ws(2)] * 2 + 1)
            x_mat = x_window_mat(1 : window_size(1), 1 : window_size(2), 1 : window_size(3));
            y_mat = y_window_mat(1 : window_size(1), 1 : window_size(2), 1 : window_size(3));
            z_mat = z_window_mat(1 : window_size(1), 1 : window_size(2), 1 : window_size(3));
        else
            x_mat = x_window_mat;
            y_mat = y_window_mat;
            z_mat = z_window_mat;
        end

        % fit
        % tic
        switch FitGaussianMethod
            case 'original'
                [prm, prmStd, ~, res] = fitGaussian3D(window, [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)], mode);
        % toc;
            case 'CD_Newton'
                % tic
                % xruan 07/29/2019 if the window size is smaller than the precomputed
                % coordinates, change the precomputed coordinates
                % tic
                % [prm_1, prmStd_1, ~, res_1] = fitGaussian3D_test_type_2_refine(window, x_mat, y_mat, z_mat, [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)], mode);
                if strcmp(mode, 'xyzAc')
                	[prm, prmStd, ~, res] = XR_fitGaussian3D(window, x_mat, y_mat, z_mat, [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)], mode, 10);
                    % [prm_1, prmStd_1, ~, res_1] = XR_fitGaussian3D(window_orig, x_mat, y_mat, z_mat, [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)], mode, 10);
                elseif strcmp(mode, 'Ac')
                	[prm, prmStd, ~, res] = XR_fitGaussian3D_Ac(window, x_mat, y_mat, z_mat, [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)]);                    
                end
                % toc;
            case 'CD_LM'
                % tic
                [prm, prmStd, ~, res] = fitGaussian3D_test_type_3(window, x_mat, y_mat, z_mat, [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)], mode);
                % toc
        end
        
        if false && any(maskWindow ~= 0, 'all')
            prm
            prm_1
            sum(maskWindow~=0, 'all')

            flag = 1;
        end
        
        if false
            tic
            [prm_0, prmStd_0, ~, res_0] = fitGaussian3D(window, [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)], mode);
            toc;
            % tic
            % xruan 07/29/2019 if the window size is smaller than the precomputed
            % coordinates, change the precomputed coordinates
            tic
            % [prm_1, prmStd_1, ~, res_1] = fitGaussian3D_test_type_2_refine(window, x_mat, y_mat, z_mat, [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)], mode);
            [prm_1, prmStd_1, ~, res_1] = XR_fitGaussian3D(window, x_mat, y_mat, z_mat, [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)], mode, 8);
            toc;
            tic
            [prm_2, prmStd_2, ~, res_2] = XR_fitGaussian3D(window, x_mat, y_mat, z_mat, [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)], mode, 10);
            toc
            
            res_mat(p, :) = [res_0.RSS - res_1.RSS, res_0.RSS - res_2.RSS, res_1.RSS - res_2.RSS];
            if ~false && (res_mat(p, 1) < -1e-3 || res_mat(p, 2) > 1e-3 || res_mat(p, 3) > 1e-3)
                % res_mat(p, :)
                flag = 1;
                if any(res_mat(p, :) .* [-1, 1, 1] > 1e4)
                    flag = 1;
                end
            end
        end
        
        dx = prm(1)-ox;
        dy = prm(2)-oy;
        dz = prm(3)-oz;
        
        % exclude points where localization failed
        if (dx > -w2(1) && dx < w2(1) && dy > -w2(1) && dy < w2(1) && dz > -w2(2) && dz < w2(2) && prm(4)<2*diff(iRange))
            
            pStruct.x(p) = xi(p) + dx;
            pStruct.y(p) = yi(p) + dy;
            pStruct.z(p) = zi(p) + dz;
            pStruct.s(:,p) = prm(5:6);
            pStruct.A(p) = prm(4);
            pStruct.c(p) = prm(7);
            
            stdVect = zeros(1,7);
            stdVect(estIdx) = prmStd;
            
            pStruct.x_pstd(p) = stdVect(1);
            pStruct.y_pstd(p) = stdVect(2);
            pStruct.z_pstd(p) = stdVect(3);
            pStruct.A_pstd(p) = stdVect(4);
            pStruct.s_pstd(:,p) = stdVect(5:6);
            pStruct.c_pstd(p) = stdVect(7);
            
            pStruct.sigma_r(p) = res.std;
            pStruct.RSS(p) = res.RSS;
            
            pStruct.SE_sigma_r(p) = res.std/sqrt(2*(npx-1));
            SE_sigma_r = pStruct.SE_sigma_r(p) * kLevel;
            
            pStruct.hval_AD(p) = res.hAD;
            
            % H0: A <= k*sigma_r
            % H1: A > k*sigma_r
            sigma_A = stdVect(4);
            A_est = prm(4);
            df2(p) = (npx-1) * (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
            scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)/npx);
            T(p) = (A_est - res.std*kLevel) ./ scomb;
        end
    end
end
% 1-sided t-test: A_est must be greater than k*sigma_r
pStruct.pval_Ar = tcdf(-T, df2);
pStruct.hval_Ar = pStruct.pval_Ar < ip.Results.AlphaT;
