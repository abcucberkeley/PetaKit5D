function out = pkpy_decl_numeric(img)
% A fake function with one of the io wrappers' arguments-block declarations (the lab's text
% and default), for the oracle cases of tests/params/test_declarations.py: {class, size,
% value} of the argument as MATLAB converts it.
arguments
    img {mustBeNumeric}
end
out = {class(img), size(img), img};
end
