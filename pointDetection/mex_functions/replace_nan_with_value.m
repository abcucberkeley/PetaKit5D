function [x] = replace_nan_with_value(x, val)
% replace the NaNs of a single or double array x with val, in x's class: the
% MATLAB code of replace_nan_with_value.c, with its checks, which MATLAB runs
% where the mex isn't built (macOS); elsewhere the mex in this folder runs


if nargin ~= 2
    error('data:inputError', 'Two input required.');
end
if nargout ~= 1
    error('data:outputError', 'One output required.');
end
if ~isa(x, 'single') && ~isa(x, 'double')
    error('data:TypeError', 'Only single and double types are supported!');
end
if issparse(x)
    error('data:inputError', 'The input must be a full array.');
end
if ~(isnumeric(val) || islogical(val)) || issparse(val) || ~isreal(val) || isempty(val)
    error('data:inputError', 'The value must be a full real numeric scalar.');
end

% the first element as a double, as mxGetScalar gives it
x(isnan(x)) = cast(double(val(1)), class(x));

end
