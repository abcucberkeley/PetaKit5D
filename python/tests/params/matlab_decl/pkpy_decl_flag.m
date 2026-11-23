function out = pkpy_decl_flag(options)
% A fake function with one of the io wrappers' arguments-block declarations (the lab's text
% and default), for the oracle cases of tests/params/test_declarations.py: {class, size,
% value} of the argument as MATLAB converts it.
arguments
    options.groupWrite (1, 1) logical = true
end
x = options.groupWrite;
out = {class(x), size(x), x};
end
