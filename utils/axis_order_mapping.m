function [order_mat] = axis_order_mapping(inputAxisOrder, outputAxisOrder)
% map axis order based on the axis order strings for input and output
%
% Author: Xiongtao Ruan (10/08/2024)


if nargin < 2
    outputAxisOrder = 'yxz';
end

inputAxisOrder = lower(inputAxisOrder);
if numel(inputAxisOrder) ~= 3 || ~any(contains({'xyz', 'yxz', 'zyx', 'zxy', 'yzx', 'xzy'}, inputAxisOrder))
    error('Input axis order must contains x, y, and z.')
end

if numel(outputAxisOrder) ~= 3 || ~any(contains({'xyz', 'yxz', 'zyx', 'zxy', 'yzx', 'xzy'}, outputAxisOrder))
    error('Output axis order must contains x, y, and z.')
end

order_mat = zeros(1, 3); 

for i = 1 : 3
    order_mat(i) = strfind(outputAxisOrder, inputAxisOrder(i));
end

end