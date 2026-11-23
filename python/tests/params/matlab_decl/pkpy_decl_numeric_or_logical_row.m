function out = pkpy_decl_numeric_or_logical_row(options)
% A fake function with one of the io wrappers' arguments-block declarations (the lab's text
% and default), for the oracle cases of tests/params/test_declarations.py: {class, size,
% value} of the argument as MATLAB converts it.
arguments
    options.sparseData (1, :) {mustBeNumericOrLogical} = true
end
x = options.sparseData;
out = {class(x), size(x), x};
end
