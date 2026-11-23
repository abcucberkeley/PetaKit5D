function s = pkpy_decl_describe(o)
% The outcome of a call of a pkpy_decl_* fake function, as tests/params/test_declarations.py
% records it: 'accepted <class> <size>', then mat2str of a 2-D real numeric or logical value of
% at most 12 elements, or the codes of a 2-D char value of at most 24 characters.
v = o{3};
s = sprintf('accepted %s %s', o{1}, mat2str(o{2}));
if ismatrix(v) && numel(v) <= 12 && ((isnumeric(v) && isreal(v)) || islogical(v))
    s = [s ' ' mat2str(v)];
elseif ischar(v) && ismatrix(v) && numel(v) <= 24
    s = [s ' chars ' mat2str(double(v))];
end
end
