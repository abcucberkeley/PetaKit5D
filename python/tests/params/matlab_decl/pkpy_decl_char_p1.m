function out = pkpy_decl_char_p1(filepath)
% A fake function with one of the io wrappers' arguments-block declarations (the lab's text
% and default), for the oracle cases of tests/params/test_declarations.py: {class, size,
% value} of the argument as MATLAB converts it.
arguments
    filepath char
end
out = {class(filepath), size(filepath), filepath};
end
