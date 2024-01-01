function [prm, prmStd, C, res, J_mat] = XR_fitGaussian3D_Ac(window, x_mat, y_mat, z_mat, init_variable)
% copied from XR_fitGaussian3D.m and remove lots of unnesseary parts to
% accelarte the computing

if nargin < 1
%     load('test_fitGaussian3D_function.mat');
%     init_variable = [X(p,1)-xi(p)+ox X(p,2)-yi(p)+oy X(p,3)-zi(p)+oz A(p) sigma(p,:) c(p)];    
    load('test_fitGaussian3D_function_32.mat');
    method_type = 10;
    mode = '';
end

% tic
f_orig_mat = window(:);
nonan_inds = ~isnan(f_orig_mat);
f_mat = f_orig_mat(nonan_inds);
x_mat = x_mat(nonan_inds);
y_mat = y_mat(nonan_inds);
z_mat = z_mat(nonan_inds);
num_data = numel(f_mat);
ones_mat = ones(num_data, 1);

% xruan 11/06/2019 check if all values in the window are the same. 
if all(f_mat == f_mat(1))
    prm = nan(1, 6);
    prmStd = nan(1, 5);
    C = nan;
    res = nan;
    J_mat = nan;
    return;
end

% nan_param_inds = isnan(init_variable);
% if any(nan_param_inds)
%     if nan_param_inds(7)
%         init_variable(7) = sum(f_mat) / num_data;
%     end    
%     if any(nan_param_inds(1 : 3))
%         [~, max_ind] = max(f_mat);
%         init_variable(1 : 3) = [x_mat(max_ind), y_mat(max_ind), z_mat(max_ind)];
%         init_variable(4) = max_f - init_variable(7);
%     elseif nan_param_inds(4)
%         init_variable(4) = max(f_mat) - init_variable(7);
%     end        
% end

mu_x_0 = init_variable(1);
mu_y_0 = init_variable(2);
mu_z_0 = init_variable(3);
sigma_xy = init_variable(5);
sigma_z = init_variable(6);

h_mat = exp(-((x_mat - mu_x_0) .^ 2 + (y_mat - mu_y_0) .^ 2) / (2 * sigma_xy ^ 2) - (z_mat - mu_z_0) .^ 2 / (2 * sigma_z ^ 2));

beta_i = init_variable;

% directly obtain solution
h_mat_sq = h_mat' * h_mat;
h_sum = h_mat' * ones_mat;
ss_xx = num_data * h_mat_sq - h_sum ^ 2;
beta_A = (num_data * sum(f_mat .* h_mat) - h_sum * sum(f_mat)) / ss_xx;
A_h_mat = beta_A * h_mat;
beta_c = ones_mat' * (f_mat - A_h_mat) / num_data;
beta_i(4) = beta_A;
beta_i(7) = beta_c;

g_mat = f_mat - A_h_mat - beta_c;
F_mat = g_mat' * g_mat;

% F_mat(F_mat == 0) = [];
prm = beta_i';
res.RSS = F_mat;    

num_param = 2;        
C = [num_data, -h_sum; -h_sum, h_mat_sq] ./ ((h_mat_sq * num_data - h_sum ^ 2) * 2);
prmStd = sqrt((res.RSS / (num_data - num_param - 1)) * diag(C))';

% residual mean and standard deviation
res_mean = sum(g_mat) ./ num_data;
res_std = sqrt((sum(g_mat .^ 2) / num_data - res_mean .^ 2) * (num_data / (num_data - 1)));

% BIC = num_data * log(res.RSS/num_data) + num_param * log(num_data);

mu = 0;
alpha = 0.05;
mode = 2;
% [hAD, p] = adtest_XR(g_mat, mu, res_std, alpha, mode);
[hAD] = adtest_mex(g_mat, mu, res_std, alpha, mode);
p = nan;

res.hAD = hAD;
res.pval = p;
res.mean = res_mean;
res.std = res_std;
res.data = zeros(size(window));
res.data(nonan_inds) = g_mat;
res.C = C;
% res.BIC = BIC;

end


