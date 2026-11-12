function y = pkpy_f02_table(a, op)
% every pair of a's values: y(i, j) = a(i) op a(j)
[A, B] = ndgrid(a(:), a(:));
y = pkpy_f02_arith(A, B, op);
end
