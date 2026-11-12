function y = pkpy_f02_arith(a, b, op)
% element-wise operators in a function file
switch op
    case 'plus'
        y = a + b;
    case 'minus'
        y = a - b;
    case 'times'
        y = a .* b;
    case 'rdivide'
        y = a ./ b;
    case 'ldivide'
        y = a .\ b;
end
end
