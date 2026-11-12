function y = pkpy_f02_scalar_table(a, d, op, side)
% y(:, k) = a op d(k) (side 'r') or d(k) op a (side 'l'), with d(k) a double scalar
a = a(:);
y = zeros(numel(a), numel(d), class(a));
for k = 1:numel(d)
    if side == 'r'
        y(:, k) = pkpy_f02_arith(a, d(k), op);
    else
        y(:, k) = pkpy_f02_arith(d(k), a, op);
    end
end
end
