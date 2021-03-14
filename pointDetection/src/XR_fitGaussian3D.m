function [prm, prmStd, C, res, J_mat] = XR_fitGaussian3D(window, x_mat, y_mat, z_mat, init_variable, mode, method_type)
% Fitting a single 3D Gaussian model for given window at given location use
% block coordinate descent method as the back optimization framework. 
%   Input:      window : 3D matrix of voxel values for a window
%                x_mat : x coordinates of the voxels
%                y_mat : y coordinates of the voxels
%                z_mat : z coordinates of the voxels
%        init_variable : initial values for the fitting: [x, y, z, A, sigma_xy, sigma_z, c]
%                 mode : parameters to estimate. Default: 'xyzAc'.
%          method_type : Gaussian fitting method. Default: 10 (method 10, fixed)
% 
% 
%  output:         prm : fitted parameters
%               prmStd : standarnd deviations of fitted parameters.
%                    C : Residual error covariance matrix of the parameters
%                  res : residual information
%                J_mat : Jacobian matrix of the fitted parameters. 
% 
% 
% xruan 11/06/2019, first check if the vol contains all the same values. 
% xruan 11/10/2019 change to use Hessian matrix to calculate parameter
% standard deviations. 
% xruan 11/13/2019 add early stop criterion in case it cannot improve the
% best solution for a long time. Also, check if the step is small enough
% with the combination of |F(i) - F(i + 1)|/F(1)
% xruan 11/15/2019 if there are too many nan values ( > 10%), skip it.  

if nargin < 6
    mode = 'xyzAc';
end
if nargin < 7 
    method_type = 10;
end

if nargin < 1
%     load('test_fitGaussian3D_function.mat');
%     init_variable = [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)];    
    load('test_fitGaussian3D_function_30.mat');
    method_type = 10;
    mode = 'xyzAc';
end

% default parameters
early_stop_iter_num = 30; % if it cannot be improved in given iterations, stop the update
t = 0.0001;
nIter = 2000;
F_epsilon = 1e-6;
step_epsilon = 1e-4; % threshold of max of absolute value of x, y, z step
F_ratio_epsilon = 0.1 * F_epsilon; % set |F(i) - F(i + 1)|/F(1) the same as F_epsilton
nan_num_thresh = 0.30; % if the ratio of nan value in the matrix is large, then skip the fitting (usually they are in the boundary). 

% tic
f_orig_mat = window(:);
% [x_mat, y_mat, z_mat] = meshgrid(1 : size(window, 2), 1 : size(window, 1), 1 : size(window, 3));
x_orig_mat = x_mat(:);
y_orig_mat = y_mat(:);
z_orig_mat = z_mat(:);
nonan_inds = ~isnan(f_orig_mat);
f_mat = f_orig_mat(nonan_inds);
x_mat = x_orig_mat(nonan_inds);
y_mat = y_orig_mat(nonan_inds);
z_mat = z_orig_mat(nonan_inds);
num_data = numel(f_mat);
max_f = max(f_mat);
max_x = max(x_mat);
max_y = max(y_mat);
max_z = max(z_mat);
% half_max_coords_vec = [max_x; max_y; max_z] / 2;
max_coords_vec = [max_x; max_y; max_z];
ones_mat = ones(num_data, 1);
persistent eye_2 eye_3 se_18
if isempty(eye_2)
    eye_2 = eye(2);
    eye_3 = eye(3);
    se_18 = true(3, 3, 3); % 18-conn matrix
    se_18([1, 3, 7, 9, 19, 21, 25, 27]) = false;
end
% toc

% xruan 11/06/2019 check if all values in the window are the same. 
if all(f_mat == f_mat(1)) % || sum(~nonan_inds) / num_data > nan_num_thresh
    prm = nan(1, 6);
    prmStd = nan(1, 5);
    C = nan;
    res = nan;
    J_mat = nan;
    return;
end

nan_param_inds = isnan(init_variable);
if any(nan_param_inds)
    if nan_param_inds(7)
        init_variable(7) = sum(f_mat) / num_data;
    end    
    if any(nan_param_inds(1 : 3))
        [~, max_ind] = max(f_mat);
        init_variable(1 : 3) = [x_mat(max_ind), y_mat(max_ind), z_mat(max_ind)];
        init_variable(4) = max_f - init_variable(7);
    elseif nan_param_inds(4)
        init_variable(4) = max(f_mat) - init_variable(7);
    end        
end

mu_x_0 = init_variable(1);
mu_y_0 = init_variable(2);
mu_z_0 = init_variable(3);
A_0 = init_variable(4);
sigma_xy = init_variable(5);
sigma_z = init_variable(6);
c_0 = init_variable(7);

beta_im1 = [mu_x_0; mu_y_0; mu_z_0; A_0; c_0];
h_mat = exp(-((x_mat - mu_x_0) .^ 2 + (y_mat - mu_y_0) .^ 2) / (2 * sigma_xy ^ 2) - (z_mat - mu_z_0) .^ 2 / (2 * sigma_z ^ 2));

% grad_F_epsilon = 1e-7;
F_mat = zeros(nIter + 1, 1);
g_mat = f_mat - beta_im1(4) * h_mat - beta_im1(5);
F_mat(1) = g_mat' * g_mat;

% methods. 
methods_set = {'cd', 'cd_1', 'cd_2', 'cd_2_accelerated', 'cd_2_backtracking', 'cd2_newton_raw', 'cd2_newton', 'cd2_newton_1', 'cd2_newton_2', 'cd2_newton_3', 'cd_3'};
method = methods_set{method_type};

if strcmp(mode, 'xyzAc')
    switch method
        case 'cd2_newton_1'
            alpha = 0.25;
            B = 0.5;
            t_min = 1e-8;
            counter = 0;
            max_count = 10;
            beta_i = beta_im1;  
            beta_xyz = beta_i(1 : 3);
            beta_A = beta_i(4);
            beta_c = beta_i(5);
            F_best = Inf;
            
            for i = 1 : nIter
                [grad_vec, hessian_mat] = compute_gradient_hessian(x_mat, y_mat, z_mat, beta_xyz, beta_A, g_mat, sigma_xy, sigma_z, h_mat);

                % repair hessian matrix
                if hessian_mat(1, 1) < 0 || hessian_mat(1, 1) * hessian_mat(2, 2) - hessian_mat(2, 1) * hessian_mat(1, 2) < 0 || det(hessian_mat) < 0
                    eig_vals = eig(hessian_mat);                
                    if sum(eig_vals < 0) == 3
                        hessian_mat = -hessian_mat;
                        hessian_mat = hessian_mat + 10 * eye_3;
                    else
                        min_eigval = min(eig_vals);
                        hessian_mat = hessian_mat - (min_eigval - max(0.5 * min(abs(eig_vals)), 10)) * eye_3;
                    end
                end
                v = -hessian_mat \ grad_vec;

                % beta_t = beta_im1;
                % beta_t(1 : 3) = beta_im1(1 : 3) + t * v;
                % t = 1;
                beta_t = beta_xyz + v;
                if any(beta_t < -1.5 | beta_t > max_coords_vec + 1.5) || beta_A < 0 
                    v_pos_ind = v > 0;
                    t = min([((max_coords_vec + 1.5) .* v_pos_ind - 1.5 * (1 - v_pos_ind) - beta_xyz) ./ v; 1]);
                    if t <= 0 || beta_A < 0
                        if grad_vec' * grad_vec > 1e-4 && counter <= max_count
                            if counter == 0
                                % record the best performance one
                                % F_best = F_mat(i);
                                % beta_i_best = [beta_xyz; beta_A; beta_c];

                                window_seg = window > prctile(f_mat, 90);
                                CC = bwconncomp(window_seg, 18);
                                volume_mat = cellfun(@numel, CC.PixelIdxList);
                                % only keep large region
                                volume_thrsh = 2;
                                obj_num_chosen = sum(volume_mat >= volume_thrsh);
                                if obj_num_chosen == 0
                                    obj_num_chosen = 3;
                                end
                                [~, inds] = maxk(volume_mat, obj_num_chosen);
                                centroids_mat = zeros(numel(inds), 3);
                                % min_f = min(f_mat);
                                for j = 1 : numel(inds)
                                    PixelIdxList_j = CC.PixelIdxList{inds(j)};
                                    f_mat_j = f_orig_mat(PixelIdxList_j);
                                    centroids_j = (f_mat_j .^ 2)' * [x_orig_mat(PixelIdxList_j), y_orig_mat(PixelIdxList_j), z_orig_mat(PixelIdxList_j)] / sum(f_mat_j .^ 2);
                                    centroids_mat(j, :) = centroids_j;
                                end
                            end
                            if counter + 1 <= size(centroids_mat, 1)
                                % if F_mat(i) < F_best
                                %     F_best = F_mat(i);
                                %    beta_i_best = [beta_xyz; beta_A; beta_c];
                                % end
                                beta_xyz = centroids_mat(counter + 1, :)';
                                counter = counter + 1;
                            end
                        end
                    end
                    beta_t = beta_xyz + t * v;
                end

                h_mat = exp(-((x_mat - beta_t(1)) .^ 2 + (y_mat - beta_t(2)) .^ 2) / (2 * sigma_xy ^ 2) - (z_mat - beta_t(3)) .^ 2 / (2 * sigma_z ^ 2));
                g_mat = f_mat - beta_A * h_mat - beta_c;                    

                alpha_grad_v = alpha * grad_vec' * v;            
                while g_mat' * g_mat > F_mat(i) + t * alpha_grad_v 
                    t = 0.4 * B * t;
                    if t < t_min
                        break;
                    end                    
                    % beta_t = beta_i(1 : 3) + t * v;
                    beta_t = beta_xyz + t * v;
                    h_mat = exp(-((x_mat - beta_t(1)) .^ 2 + (y_mat - beta_t(2)) .^ 2) / (2 * sigma_xy ^ 2) - (z_mat - beta_t(3)) .^ 2 / (2 * sigma_z ^ 2));  
                    g_mat = f_mat - beta_A * h_mat - beta_c;                                        
                end

                beta_xyz = beta_t;

                % update A and c
                % beta_i(4) = sum((f_mat - beta_i(5)) .* h_mat) / sum(h_mat .^ 2);
                beta_A = (f_mat - beta_c)' * h_mat / (h_mat' * h_mat);
                beta_c = ones_mat' * (f_mat - beta_A * h_mat) / num_data;

                % beta_im1 = beta_i;
                g_mat = f_mat - beta_A * h_mat - beta_c;
                F_mat(i + 1) = g_mat' * g_mat;

                if F_mat(i + 1) < F_best
                    F_best = F_mat(i + 1);
                    beta_xyz_best = beta_xyz; 
                    beta_A_best = beta_A; 
                    beta_c_best = beta_c;
                end

                if F_mat(i) - F_mat(i + 1) < F_epsilon && F_mat(i) - F_mat(i + 1) > -F_epsilon % || F_mat(i + 1) - F_mat(i) > 1e-2
                    break;
                end            
            end
            beta_i = [beta_xyz; beta_A; beta_c];
            if F_mat(i + 1) > F_best
                beta_i = [beta_xyz_best; beta_A_best; beta_c_best];           
            end
        case 'cd2_newton_3'
            % based on cd2_newton_1, use local maximum of log to pick candidate
            % points
            alpha = 0.25;
            B = 0.5;
            t_min = 1e-8;
            xyz_bound = 0;  % distance that xyz can be outside
            A_bound = max_f - min(f_mat);  % distance that xyz can be outside

            counter = 0;
            max_count = 15;
            beta_i = beta_im1;  
            beta_xyz = beta_i(1 : 3);
            beta_A = beta_i(4);
            beta_c = beta_i(5);
            i_best = 1;
            F_best = F_mat(1);
            persistent log_kernel sigma_mat

            h_A_mat = beta_A * h_mat;

            for i = 1 : nIter                
                x_c_mat = (x_mat - beta_xyz(1)) / sigma_xy ^ 2;
                y_c_mat = (y_mat - beta_xyz(2)) / sigma_xy ^ 2;
                z_c_mat = (z_mat - beta_xyz(3)) / sigma_z ^ 2;
                [grad_vec, hessian_mat] = compute_gradient_hessian_mex(x_c_mat, y_c_mat, z_c_mat, h_A_mat, g_mat, sigma_xy, sigma_z);
                % [grad_vec_0, hessian_mat_0] = compute_gradient_hessian(x_c_mat, y_c_mat, z_c_mat, h_A_mat, g_mat, sigma_xy, sigma_z);
                % repair hessian matrix
                if hessian_mat(1, 1) < 0 || ...
                        hessian_mat(1, 1) * hessian_mat(2, 2) - hessian_mat(2, 1) ^ 2 < 0 ...
                        || hessian_mat(1, 1) * hessian_mat(2, 2) * hessian_mat(3, 3) + 2 * hessian_mat(1, 2) * hessian_mat(1, 3) * hessian_mat(2, 3) - ...
                           hessian_mat(1, 1) * hessian_mat(2, 3) ^ 2 - hessian_mat(2, 2) * hessian_mat(1, 3) ^ 2 - hessian_mat(3, 3) * hessian_mat(1, 2) ^ 2 < 0
                        % || det(hessian_mat) < 0
                    eig_vals = eig(hessian_mat);                
                    if sum(eig_vals < 0) == 3
                        hessian_mat = -hessian_mat;
                        hessian_mat = hessian_mat + 10 * eye_3;
                    else
                        min_eigval = min(eig_vals);
                        hessian_mat = hessian_mat - (min_eigval - max(0.5 * min(abs(eig_vals)), 10)) * eye_3;
                    end
                end
                v = -hessian_mat \ grad_vec;
                
                t = 1;
                beta_t = beta_xyz + t * v;
                if any(beta_t < -xyz_bound | beta_t > max_coords_vec + xyz_bound) || beta_A < 0 || beta_A > A_bound
                    v_pos_ind = v > 0;
                    t = max(min([((max_coords_vec + xyz_bound) .* v_pos_ind - xyz_bound * (1 - v_pos_ind) - beta_xyz) ./ v; 1]), 0);
                    if t <= 0 || beta_A < 0 || beta_A > A_bound
                        if grad_vec' * grad_vec > 1e-4 && counter <= max_count
                            if counter == 0
                                % faster way to obtain the threshold
                                f_mat_sort = sort(f_mat);
                                f_thrsh = f_mat_sort(floor(numel(f_mat) * 0.8));
                                window_seg = window > f_thrsh;

                                % calculate log
                                if isempty(sigma_mat) || any(sigma_mat ~= [sigma_xy, sigma_xy, sigma_z])
                                    sigma_mat = [sigma_xy, sigma_xy, sigma_z];                                
                                    log_kernel = -fspecial3('log', [1, 1, 1] * 2 + 1, sigma_mat);
                                end
                                % imgLoG = imfilter(window, log_kernel, 'conv', 'same', 'symmetric');
                                imgLoG = imagesbuiltinImfilter(window, log_kernel, 'symmetric', 'same', 'conv');
                                
                                % allMax = XR_locmax3d(imgLoG, 2*ceil([sigma_xy, sigma_xy, sigma_z])+1, 'ClearBorder', false);
                                allMax = XR_locmax3d(imgLoG, 2*ceil([sigma_xy, sigma_xy, sigma_z])+1);

                                allMax_1 = allMax .* window_seg;

                                allMax_indices = find(allMax_1 ~= 0);

                                [~, sort_inds] = sort(allMax_1(allMax_indices), 'descend');
                                allMax_mat = [x_orig_mat(allMax_indices(sort_inds)), y_orig_mat(allMax_indices(sort_inds)), z_orig_mat(allMax_indices(sort_inds))];

                                % make sure there are at least a given number
                                % of points, if not, pick other local maximum
                                % points from high to low LoG values. 
                                min_num_chosen = 3;
                                if numel(allMax_indices) < min_num_chosen
                                    % if there is no local maximum in the
                                    % largest region, use the centroid.                                     
                                    allMax_indices_0 = find(allMax .* (~window_seg) > 0);
                                    [~, sort_inds] = sort(allMax(allMax_indices_0), 'descend');
                                    add_inds = allMax_indices_0(sort_inds(1 : min(numel(allMax_indices_0), min_num_chosen - numel(allMax_indices))));
                                    allMax_mat = [allMax_mat; x_orig_mat(add_inds), y_orig_mat(add_inds), z_orig_mat(add_inds)];
                                end

                                % CC = bwconncomp(window_seg, 18);
                                % use builtin function, which is much faster.
                                % CC.pixelIdxList = [];
                                if verLessThan('matlab','9.9')
                                    PixelIdxList = builtin('_pixelIdxListsn', window_seg, se_18); 
                                else
                                    PixelIdxList = images.internal.builtins.pixelIdxListsn(window_seg, se_18); 
                                end
                                volume_mat = cellfun(@numel, PixelIdxList);
                                % only keep large region
                                volume_thrsh = 2;
                                obj_num_chosen = sum(volume_mat >= volume_thrsh);
                                if obj_num_chosen == 0
                                    obj_num_chosen = 3;
                                end
                                [~, inds] = maxk(volume_mat, obj_num_chosen);

                                centroids_mat = zeros(numel(inds), 3);
                                for j = 1 : numel(inds)
                                    PixelIdxList_j = PixelIdxList{inds(j)};
                                    f_mat_j = f_orig_mat(PixelIdxList_j);
                                    centroids_j = f_mat_j' * [x_orig_mat(PixelIdxList_j), y_orig_mat(PixelIdxList_j), z_orig_mat(PixelIdxList_j)] / sum(f_mat_j);
                                    centroids_mat(j, :) = centroids_j;
                                end

                                % centroids_round_mat = sub2ind(size(window_seg), round(centroids_mat(:, 2)), round(centroids_mat(:, 1)), round(centroids_mat(:, 3)));
                                allMax_mat = [allMax_mat; centroids_mat];
                                % allMax_mat = [centroids_mat; allMax_mat];

                                allMax_mat(all(allMax_mat == init_variable(1 : 3), 2), :) = [];
                            end
                            if counter + 1 <= size(allMax_mat, 1)
                                beta_xyz = allMax_mat(counter + 1, :)';
                                t = 0;
                                counter = counter + 1;
                            end
                        end
                    end
                    beta_t = beta_xyz + t * v;
                end
                
                % h_mat_0 = h_mat_mex(x_mat, y_mat, z_mat, beta_t, sigma_xy, sigma_z);
                h_mat = exp(-((x_mat - beta_t(1)) .^ 2 + (y_mat - beta_t(2)) .^ 2) / (2 * sigma_xy ^ 2) - (z_mat - beta_t(3)) .^ 2 / (2 * sigma_z ^ 2));
                h_A_mat = beta_A * h_mat; 
                g_mat = f_mat - h_A_mat - beta_c;                    

                alpha_grad_v = alpha * (grad_vec' * v);
                while t > t_min && g_mat' * g_mat > F_mat(i) + t * alpha_grad_v 
                    % t = 0.4 * B * t;
                    t = B * t;
                    if t < t_min
                        break;
                    end                    
                    % beta_t = beta_i(1 : 3) + t * v;
                    beta_t = beta_xyz + t * v;
                    h_mat = exp(-((x_mat - beta_t(1)) .^ 2 + (y_mat - beta_t(2)) .^ 2) / (2 * sigma_xy ^ 2) - (z_mat - beta_t(3)) .^ 2 / (2 * sigma_z ^ 2));  
                    % h_mat = h_mat_fast_exp_mex(x_mat, y_mat, z_mat, beta_t, sigma_xy, sigma_z);
                    h_A_mat = beta_A * h_mat;                 
                    g_mat = f_mat - h_A_mat - beta_c;                                        
                end

                beta_xyz = beta_t;

                % update A and c
                % beta_i(4) = sum((f_mat - beta_i(5)) .* h_mat) / sum(h_mat .^ 2);
                beta_A = (f_mat - beta_c)' * h_mat / (h_mat' * h_mat);
                h_A_mat = beta_A * h_mat;                             
                beta_c = ones_mat' * (f_mat - h_A_mat) / num_data;
                % not fast by solving linear system
                % beta_A = (f_mat - beta_c - mean(f_mat - beta_c))' * (h_mat - mean(h_mat)) / ((h_mat - mean(h_mat))' * (h_mat - mean(h_mat)));
                % beta_c = ones_mat' * (f_mat - beta_A * h_mat) / num_data;

                % beta_im1 = beta_i;
                g_mat = f_mat - h_A_mat - beta_c;
                F_mat(i + 1) = g_mat' * g_mat;

                if F_mat(i + 1) < F_best && beta_A > 0
                    F_best = F_mat(i + 1);
                    beta_xyz_best = beta_xyz; 
                    beta_A_best = beta_A; 
                    beta_c_best = beta_c;
                    i_best = i;                
                end

                % xruan 11/13/2019
                if i - i_best > early_stop_iter_num
                    break;
                end

                % xruan 11/13/2019 add another criterion for stop
                delta_F = abs(F_mat(i) - F_mat(i + 1));
                if delta_F < F_epsilon || (delta_F / F_mat(1) < F_ratio_epsilon && max(abs(v)) < step_epsilon)
                    break;
                end            
            end
            beta_i = [beta_xyz; beta_A; beta_c];
            if (F_mat(i + 1) > F_best || beta_A < 0) && (exist('beta_A_best', 'var') && beta_A_best > 0)
                beta_i = [beta_xyz_best; beta_A_best; beta_c_best];
                h_mat = exp(-((x_mat - beta_i(1)) .^ 2 + (y_mat - beta_i(2)) .^ 2) / (2 * sigma_xy ^ 2) - (z_mat - beta_i(3)) .^ 2 / (2 * sigma_z ^ 2));            
                % h_mat = h_mat_fast_exp_mex(x_mat, y_mat, z_mat, beta_i, sigma_xy, sigma_z);
                h_A_mat = h_mat * beta_i(4);
                g_mat = f_mat - h_A_mat - beta_i(5);
            end
    end
elseif strcmp(mode, 'Ac')
    beta_i = beta_im1;
    beta_xyz = beta_i(1 : 3);
    % h_mat = exp(-((x_mat - beta_xyz(1)) .^ 2 + (y_mat - beta_xyz(2)) .^ 2) / (2 * sigma_xy ^ 2) - (z_mat - beta_xyz(3)) .^ 2 / (2 * sigma_z ^ 2));
    
    % directly obtain solution
    h_mat_sq = h_mat' * h_mat;
    h_sum = h_mat' * ones_mat;
    ss_xx = num_data * h_mat_sq - h_sum ^ 2;
    beta_A = (num_data * sum(f_mat .* h_mat) - h_sum * sum(f_mat)) / ss_xx;
    beta_c = ones_mat' * (f_mat - beta_A * h_mat) / num_data;
    beta_i = [beta_xyz; beta_A; beta_c];
    
    g_mat = f_mat - beta_A * h_mat - beta_c;
    F_mat(2) = g_mat' * g_mat;
    i = 1;
end

% F_mat(F_mat == 0) = [];
F_mat = F_mat(1 : i + 1);
prm = [beta_i(1 : 4); sigma_xy; sigma_z; beta_i(5)]';
res.RSS = F_mat(end);    

if i > 100
    flag = 1;
end

if sum(~nonan_inds) / num_data > nan_num_thresh 
    flag = 1;
end

% calculate parameter standard devation and covariance
if strcmp(mode, 'xyzAc')
    num_param = 5;    
    if (F_best < F_mat(end) || beta_A < 0)
        res.RSS = F_best;
    end    
    % tic, 
    % J_mat = [-beta_A * h_mat .* [(x_mat - beta_i(1)) / sigma_xy ^ 2, (y_mat - beta_i(2)) / sigma_xy ^ 2, (z_mat - beta_i(3)) / sigma_z ^ 2],  -h_mat, -ones(size(z_mat))];
    % C = eye(numel(beta_i)) / (J_mat' * J_mat);
    % JTJ = J_mat' * J_mat;
    % C = inv(JTJ); toc
 
    % xruan 11/10/2019 use hessian matrix to calculate parameter stds,
    % which is faster and more accuarate. 
    % [hessian_mat] = compute_full_hessian_given_partial_matrix(orig_hessian_mat, x_mat, y_mat, z_mat, beta_i, g_mat, sigma_xy, sigma_z, h_mat);
    x_c_mat = (x_mat - beta_i(1)) / sigma_xy ^ 2;
    y_c_mat = (y_mat - beta_i(2)) / sigma_xy ^ 2;
    z_c_mat = (z_mat - beta_i(3)) / sigma_z ^ 2;    
    [hessian_mat] = compute_full_hessian(x_c_mat, y_c_mat, z_c_mat, h_A_mat, g_mat, sigma_xy, sigma_z, h_mat);
    C = inv(hessian_mat); % toc
    prmStd = sqrt((res.RSS / (num_data - num_param - 1)) * abs(diag(C)))';
elseif strcmp(mode, 'Ac')
    num_param = 2;        
    % directly get the inverse of Hessian 2 * [h_mat_sq, h_sum; h_sum, num_data]
    C = [num_data, -h_sum; -h_sum, h_mat_sq] ./ ((h_mat_sq * num_data - h_sum ^ 2) * 2);
    prmStd = sqrt((res.RSS / (num_data - num_param - 1)) * diag(C))';
end    

% residual mean and standard deviation
res_mean = sum(g_mat) ./ num_data;
res_std = sqrt((sum(g_mat .^ 2) / num_data - res_mean .^ 2) * (num_data / (num_data - 1)));

BIC = num_data * log(res.RSS/num_data) + num_param * log(num_data);

% dist = makedist('normal', 'mu', res_mean, 'sigma', res_std);
% [hAD, p] = adtest(g_mat, 'Distribution', dist);

% [hAD_1, p_1] = adtest((g_mat - res_mean) / res_std);
% [hAD_2, p_2] = adtest(g_mat);
mu = 0;
alpha = 0.05;
mode = 2;
% [hAD, p] = adtest_XR(g_mat, mu, res_std, alpha, mode);
[hAD] = adtest_mex(sort(g_mat), mu, res_std, alpha, mode);
p = nan;

res.hAD = hAD;
res.pval = p;
res.mean = res_mean;
res.std = res_std;
res.data = zeros(size(window));
res.data(nonan_inds) = g_mat;
res.C = C;
res.BIC = BIC;

end



function [grad_vec, hessian_mat] = compute_gradient_hessian(x_c_mat, y_c_mat, z_c_mat, h_A_mat, g_mat, sigma_xy, sigma_z)
% calculate gradient and hessian matrix 

xyz_c = [x_c_mat, y_c_mat, z_c_mat];

w_1 = 2 * h_A_mat;
w_3 = (w_1 .* g_mat)';

grad_vec = -(w_3 * xyz_c)';

if nargout > 1    
    w_2_t = (w_1 .* (h_A_mat - g_mat));
    w_3_sum = sum(w_3);
    A_h_c = w_2_t .* xyz_c;
    hessian_mat = A_h_c' * xyz_c + diag([w_3_sum / sigma_xy ^ 2, w_3_sum / sigma_xy ^ 2, w_3_sum / sigma_z ^ 2]);
end


end


function [hessian_mat] = compute_full_hessian(x_c, y_c, z_c, h_A_mat, g_mat, sigma_xy, sigma_z, h_mat)
% calculate full hessian matrix 

w_1_t = (2 * h_A_mat)';

w_2_t = w_1_t .* (w_1_t / 2 - g_mat');
w_3_sum = w_1_t * g_mat;    
w_4_t = (2 * h_mat .* (h_mat - g_mat))';

H_11 = w_2_t * (x_c .^ 2) + w_3_sum / sigma_xy ^ 2;
H_21 = w_2_t * (x_c .* y_c); 
H_31 = w_2_t * (x_c .* z_c);
H_41 = w_4_t * x_c;
H_51 = w_1_t * x_c;
H_22 = w_2_t * (y_c .^ 2) + w_3_sum / sigma_xy ^ 2;
H_32 = w_2_t * (y_c .* z_c);
H_42 = w_4_t * y_c;
H_52 = w_1_t * y_c;
H_33 = w_2_t * (z_c .^ 2) + w_3_sum / sigma_z ^ 2;
H_43 = w_4_t * z_c;
H_53 = w_1_t * z_c;
H_44 = 2 * sum(h_mat .^ 2);
H_54 = 2 * sum(h_mat);
H_55 = 2 * numel(h_mat);

hessian_mat = [H_11, H_21, H_31, H_41, H_51; 
               H_21, H_22, H_32, H_42, H_52;
               H_31, H_32, H_33, H_43, H_53;
               H_41, H_42, H_43, H_44, H_54;
               H_51, H_52, H_53, H_54, H_55];

end



