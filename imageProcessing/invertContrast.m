%[out] = invertContrast(in, rangeIn) inverts the input signal, preserving the dynamic range.
% If a dynamic range is provided, the input first gets truncated to that range.
% The built-in imcomplement.m function calculates a simpler inverse, i.e., [a b] gets mapped to [-b -a].

% Francois Aguet, 08/09/2012

function out = invertContrast(in, rangeIn)

if nargin<2
    rangeIn = [min(in(:)) max(in(:))];
end

if nargin>1
    in(in<rangeIn(1)) = rangeIn(1);
    in(in>rangeIn(2)) = rangeIn(2);
end

out = -in + sum(rangeIn);
