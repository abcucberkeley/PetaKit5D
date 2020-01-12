% pStruct = fitGaussianMixtures2D(vol, X, A, sigma, c, mode, varargin) fits 
% 3D Gaussian mixtures to the input volume at the specified locations.
%
% Inputs:   
%           vol : input volume
%             X : Nx3 matrix of initial (or fixed) (x,y,z)-positions
%             A : initial (or fixed) amplitudes
%         sigma : initial (or fixed) Gaussian PSF standard deviations
%             c : initial (or fixed) background intensities
%
% Optional inputs : ('Mask', mask) pair with a mask of spot locations
%
% Output: 
%         pStruct: structure with fields:
%                  x : estimated x-positions
%                  y : estimated y-positions
%                  z : estimated z-positions
%                  A : estimated amplitudes
%                  s : estimated standard deviations of the PSF
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
%            pval_Ar : p-value of an amplitude vs. background noise test (p > 0.05 -> significant amplitude)
%
%
% Usage for a volume with know spot locations (mask) and fixed sigma:
% fitGaussianMixtures3D(vol, X, A, sigma, 'mask', mask);
%
% See also fitGaussian3D, fitGaussians3D

% Francois Aguet, March 28 2011 (last modified: Feb 5 2013)
% updated -Gokul Upadhyayula, 2017 bug fix, line 219

function pStruct = XR_fitGaussianMixtures3D(vol, X, A, sigma, c, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('vol', @isnumeric);
ip.addRequired('X');
ip.addRequired('A');
ip.addRequired('sigma');
ip.addRequired('c');
ip.addParameter('Alpha', 0.05, @isscalar);
ip.addParameter('AlphaT', 0.05, @isscalar);
ip.addParameter('Mask', [], @islogical);
ip.addParameter('ConfRadius', []);
ip.addParameter('WindowSize', []);
ip.addParameter('maxM', 5, @isscalar);
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
mode = 'xyzAc'; % sigma values are fixed to improve stability of the fit
if ~isempty(ip.Results.Mask)
    labels = double(labelmatrix(bwconncomp(ip.Results.Mask)));
else
    labels = zeros(size(vol));
end

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

% initialize pStruct arrays
pStruct.x = cell(1,np);
pStruct.y = cell(1,np);
pStruct.z = cell(1,np);
pStruct.A = cell(1,np);
pStruct.s = cell(2,np);
pStruct.c = cell(1,np);

pStruct.x_pstd = cell(1,np);
pStruct.y_pstd = cell(1,np);
pStruct.z_pstd = cell(1,np);
pStruct.A_pstd = cell(1,np);
pStruct.s_pstd = cell(2,np);
pStruct.c_pstd = cell(1,np);

pStruct.x_init = cell(1,np);
pStruct.y_init = cell(1,np);
pStruct.z_init = cell(1,np);

pStruct.sigma_r = cell(1,np);
pStruct.SE_sigma_r = cell(1,np);
pStruct.RSS = cell(1,np);

pStruct.pval_Ar = cell(1,np);

pStruct.hval_AD = cell(1,np);
pStruct.hval_Ar = cell(1,np);

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

% xruan
[x_window_mat, y_window_mat, z_window_mat] = meshgrid(0 : ws(1) * 2, 0 : ws(1) * 2, 0 : ws(2) * 2);

n_count = zeros(1, 2);
ng_mat = zeros(np, 1);

mixtureIndex = 1;
% loop through initial points
for p = 1:np
    if p == 2074 
        flag = 1;
    end
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
    
    window = vol(ya,xa,za);
    % set any other mask components to NaN
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
        
        % initial fit with a single Gaussian 
        % Notation: reduced model: '_r', full model: '_f'                
        
        maxM = ip.Results.maxM;
        prm_f_cell = cell(maxM, 1);
        prmStd_f_cell = cell(maxM, 1);
        res_f_cell = cell(maxM, 1);
        
        for i = 1 : ip.Results.maxM
            if i == 1
                initV = [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)];
                method_type = 10; % cd2_newton_3               
                [prm_f, prmStd_f, ~, res_f] = XR_fitGaussian3D(window, x_mat, y_mat, z_mat, initV, mode, method_type);
                prm_f_cell{i} = prm_f;
                prmStd_f_cell{i} = [prmStd_f(1 : end - 1), 0, 0, prmStd_f(end)];
                res_f_cell{i} = res_f;                
            else
                [maxRes, idx] = max(res_f_cell{i - 1}.data(:));
                [y0, x0, z0] = ind2sub(size(window), idx);

                initV = [x0-1 y0-1 z0-1 maxRes prm_f_cell{i - 1}];
                method_type = 11; % block update of mus cd2_newton_4
                [prm_f, prmStd_f, ~, res_f] = XR_fitGaussianMixture3D(window, x_mat, y_mat, z_mat, initV, mode, method_type);
                % [prm_f_0, prmStd_f_0, ~, res_f_0] = fitGaussianMixture3D(window, initV, mode);
                % [prm_f_0, prmStd_f_0, ~, res_f_0] = XR_fitGaussianMixture3D(window, x_mat, y_mat, z_mat, initV, mode, 10);
                % [prm_f_2, prmStd_f_2, ~, res_f_2] = XR_fitGaussianMixture3D(window, x_mat, y_mat, z_mat, initV, mode, 10, 1e-5);
                prm_f_cell{i} = prm_f;
                prmStd_f_cell{i} = prmStd_f;
                res_f_cell{i} = res_f;    
                
                if false && (any(abs(res_f.RSS - res_f_0.RSS) > 0.01) || any(abs(res_f.RSS - res_f_0.RSS) > 0.01))
                    res_f.RSS - res_f_0.RSS
                    % res_f.RSS - res_f_1.RSS
                    n_count = n_count + ([res_f.RSS - res_f_0.RSS > 1, res_f.RSS - res_f_0.RSS < -1]); 
                    flag = 1;
                end
            end            
        end
        % use BIC to pick the optimal number of Gaussian
        BIC_mat = cellfun(@(x) x.BIC, res_f_cell);
        [~, optm_n_mix] = min(BIC_mat);
        prm_r = prm_f_cell{optm_n_mix};
        prmStd_r = prmStd_f_cell{optm_n_mix};
        res_r = res_f_cell{optm_n_mix};
        ng = optm_n_mix;
        ng_mat(p) = ng;
            
        % sigma, c are the same for each mixture
        x_est = prm_r(1:4:end-3)-ox;
        y_est = prm_r(2:4:end-3)-oy;
        z_est = prm_r(3:4:end-3)-oz;
        A_est = prm_r(4:4:end-3);
        
        % exclude points where localization failed
        if ng > 1 || (x_est > -w2(1) && x_est < w2(1) && y_est > -w2(1) && y_est < w2(1) &&...
                z_est > -w2(2) && z_est < w2(2) && A_est<2*diff(iRange))
            pStruct.x{p} = xi(p) + x_est;
            pStruct.y{p} = yi(p) + y_est;
            pStruct.z{p} = zi(p) + z_est;
            pStruct.A{p} = A_est;
            % sigma and background offset are identical for all mixture components
            pStruct.s{p} = repmat(prm_r([end-2 end-1])', [1 ng]); %%% edited by GU
            pStruct.c{p} = repmat(prm_r(end), [1 ng]);
            
            pStruct.x_pstd{p} = prmStd_r(1:4:end-3);
            pStruct.y_pstd{p} = prmStd_r(2:4:end-3);
            pStruct.z_pstd{p} = prmStd_r(3:4:end-3);
            pStruct.A_pstd{p} = prmStd_r(4:4:end-3);
            pStruct.s_pstd{p} = repmat(prmStd_r([end-2 end-1])', [1 ng]);
            pStruct.c_pstd{p} = repmat(prmStd_r(end), [1 ng]);
            
            pStruct.x_init{p} = repmat(xi(p), [1 ng]);
            pStruct.y_init{p} = repmat(yi(p), [1 ng]);
            pStruct.z_init{p} = repmat(zi(p), [1 ng]);
            
            pStruct.sigma_r{p} = repmat(res_r.std, [1 ng]);
            pStruct.RSS{p} = repmat(res_r.RSS, [1 ng]);
            
            SE_sigma_r = res_r.std/sqrt(2*(npx-1));
            pStruct.SE_sigma_r{p} = repmat(SE_sigma_r, [1 ng]);
            SE_sigma_r = SE_sigma_r * kLevel;
            
            pStruct.hval_AD{p} = repmat(res_r.hAD, [1 ng]);
            
            % H0: A <= k*sigma_r
            % H1: A > k*sigma_r
            for i = 1:ng
                sigma_A = pStruct.A_pstd{p}(i);
                A_est = pStruct.A{p}(i);
                df2 = (npx-1) * (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
                scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)/npx);
                T = (A_est - res_r.std*kLevel) ./ scomb;
                % 1-sided t-test: A_est must be greater than k*sigma_r
                pStruct.pval_Ar{p}(i) = tcdf(-T, df2);
                pStruct.hval_Ar{p}(i) = pStruct.pval_Ar{p}(i) < ip.Results.AlphaT;
                %pStruct.mask_Ar{p}(i) = sum(A_est*g>res_r.std*kLevel); % # significant pixels
            end
            if ng>1
                pStruct.mixtureIndex{p} = mixtureIndex*ones(1,ng);
                mixtureIndex = mixtureIndex+1;
            else
                pStruct.mixtureIndex{p} = 0;
            end
        end
    end
end

% concatenate cell arrays
fnames = fieldnames(pStruct);
for f = 1:numel(fnames)
    pStruct.(fnames{f}) = [pStruct.(fnames{f}){:}];
end
