function [h, p] = adtest_XR(X_mat, mu, sigma, alpha, mode)
% implementation of ADTEST Anderson-Darling goodness-of-fit hypothesis test
% for higher computing efficiency. 
%
% Based on the matlab function adtest 
% only return h, set p as NA so that it is more efficient
% 
% mode is the case number as defined in the findCriticalValues


num_X = numel(X_mat);

if nargin < 2
    mu = mean(X_mat);
    sigma = std(X_mat);
elseif nargin < 3
    sigma = sqrt(1 / num_X * (X_mat - mu) .^ 2);
end

if nargin < 4
    alpha = 0.05;
end

X_mat = sort(X_mat(:));
z_mat = (X_mat - mu) / sigma;

% calculate statstics
phi_mat = normcdf(z_mat);
phi_mat_r = flip(phi_mat);

A_sq = -num_X - 1 / num_X * (2 * (1 : num_X) - 1) * (log(phi_mat) + log(1 - phi_mat_r));

[alphas, CVs] = findCriticalValues(mode);

alpha_ind = find(alphas == alpha, 1, 'first');
h = A_sq > CVs(alpha_ind);
p = NaN;


end


function [alphas, CVs] = findCriticalValues(mode)
% From Francis's stats.h function
% Find rows of critical values at relevant significance levels.

alphas = [0.25,  0.15,  0.10,  0.05,  0.025, 0.01, 0.005];
CVs =  [1.248, 1.610, 1.933, 2.492, 3.070, 3.857, 4.500;  % /* case 0: normal, mu/sigma known */
       0.644, 0.782, 0.894, 1.087, 1.285, 1.551, 1.756;  % /* case 1: normal, mu unknown, sigma known */
       1.072, 1.430, 1.743, 2.308, 2.898, 3.702, 4.324;  % /* case 2: normal, mu known, sigma unknown */
       0.470, 0.561, 0.631, 0.752, 0.873, 1.035, 1.159;  % /* case 3: normal, mu/sigma unknown */
       0.736, 0.916, 1.062, 1.321, 1.591, 1.959, 2.244]; % /* case 4: exponential, µ unknown */

CVs = CVs(mode + 1, :);

end