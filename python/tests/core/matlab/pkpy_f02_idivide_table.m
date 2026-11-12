function y = pkpy_f02_idivide_table(a, opt)
% every pair of a's values through idivide(a(i), a(j), opt)
[A, B] = ndgrid(a(:), a(:));
y = idivide(A, B, opt);
end
