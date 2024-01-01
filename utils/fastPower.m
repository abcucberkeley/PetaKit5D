function result = fastPower(base, exponent)
% Only for integer exponent

if exponent < 0
    base = 1 ./ base;
    exponent = -exponent;
end

result = 1;
while exponent > 0
    % Check the least significant bit of the exponent
    if bitget(exponent, 1)
        result = result .* base;
    end
    % Square the base and right-shift the exponent
    base = base .* base;
    exponent = bitshift(exponent, -1);
end

end
    