function [varargout] = decon_lucy_accel(method, varargin)
% the acceleration step of decon_lucy_function's RL iterations on the CPU,
% with the sums in an order that doesn't depend on the number of threads
% (decon_sum_double: MATLAB's order at one thread):
%
% [lambda, J_4] = decon_lucy_accel('lambda', J_2, Y, J_4)
%   lambda = sum((J_2(:) - Y(:)) .* J_4, 'double') / (sum(J_4 .* J_4, 'double') + eps),
%   clamped to [0, 1] as max(min(lambda, 1), 0), and the new J_4 = J_2(:) - Y(:)
% Y = decon_lucy_accel('predict', J_2, J_3, lambda)
%   max(J_2 + lambda * (J_2 - J_3), 0)
%
% Arrays of 89000 or more elements (from there on MATLAB's sums depend on the
% number of threads) run in decon_lucy_accel_mex, each step in one pass over
% the data, with the results of the MATLAB code below; smaller ones, other
% classes, or without the mex, the MATLAB code.


switch method
    case 'lambda'
        [J_2, Y, J_4] = varargin{:};
        if use_mex(J_2, Y, J_4)
            try
                [lambda, J_4] = decon_lucy_accel_mex('lambda', J_2, Y, J_4, maxNumCompThreads);
                varargout = {lambda, J_4};
                return;
            catch ME
                mex_fallback('decon_lucy_accel_mex', ME);
            end
        end
        lambda = decon_sum_double((J_2(:) - Y(:)) .* J_4) / (decon_sum_double(J_4 .* J_4) + eps);
        lambda = max(min(lambda, 1), 0);
        J_4 = J_2(:) - Y(:);
        varargout = {lambda, J_4};
    case 'predict'
        [J_2, J_3, lambda] = varargin{:};
        if isa(lambda, 'double') && isscalar(lambda) && isreal(lambda) && use_mex(J_2, J_3)
            try
                varargout = {decon_lucy_accel_mex('predict', J_2, J_3, lambda, maxNumCompThreads)};
                return;
            catch ME
                mex_fallback('decon_lucy_accel_mex', ME);
            end
        end
        varargout = {max(J_2 + lambda * (J_2 - J_3), 0)};
    otherwise
        error('Unknown method %s, use ''lambda'' or ''predict''.', method);
end

end


function [tf] = use_mex(varargin)
% whether arrays of one class (real single or double, at least 89000 elements
% each) run in the mex; where the mex isn't available, the MATLAB code runs
% instead (recorded with mex_fallback)

a = varargin{1};
tf = (isa(a, 'single') || isa(a, 'double')) && isreal(a) && ~issparse(a) && numel(a) >= 89000;
for i = 2 : numel(varargin)
    b = varargin{i};
    tf = tf && strcmp(class(b), class(a)) && isreal(b) && ~issparse(b) && numel(b) == numel(a);
end
if tf && exist('decon_lucy_accel_mex', 'file') ~= 3
    mex_fallback('decon_lucy_accel_mex');
    tf = false;
end

end
