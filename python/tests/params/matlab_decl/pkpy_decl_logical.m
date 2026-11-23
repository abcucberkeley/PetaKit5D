function out = pkpy_decl_logical(options)
% A fake function with one of the io wrappers' arguments-block declarations (the lab's text
% and default), for the oracle cases of tests/params/test_declarations.py: {class, size,
% value} of the argument as MATLAB converts it.
arguments
    options.create logical = true
end
x = options.create;
out = {class(x), size(x), x};
end
