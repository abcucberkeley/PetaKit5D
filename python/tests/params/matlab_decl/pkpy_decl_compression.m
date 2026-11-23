function out = pkpy_decl_compression(options)
% A fake function with one of the io wrappers' arguments-block declarations (the lab's text
% and default), for the oracle cases of tests/params/test_declarations.py: {class, size,
% value} of the argument as MATLAB converts it.
arguments
    options.Compression (1, :) char {mustBeMember(options.Compression, {'none', 'lzw', 'zstd'})} = 'lzw'
end
x = options.Compression;
out = {class(x), size(x), x};
end
