function [x] = load_ref(fn)
% read an input value file, cached for the session (inputs are immutable:
% their names are their content hashes); the cache holds up to 2 GB. Each call
% returns arrays of their own (unshare): a case whose MATLAB code writes into
% an input in place through a mex must not change the input the pack's later
% cases read. PKPY_ORACLE_SHARE_INPUTS=1 returns the cached value itself, as
% the harness did before (for the harness's own test only).


persistent cache bytes
if isempty(cache)
    cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    bytes = 0;
end
if isKey(cache, fn)
    x = cache(fn);
else
    x = pkpy_oracle.read_value(fn);
    info = whos('x');
    if bytes + info.bytes > 2 ^ 31
        cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
        bytes = 0;
    end
    cache(fn) = x;
    bytes = bytes + info.bytes;
end
if ~strcmp(getenv('PKPY_ORACLE_SHARE_INPUTS'), '1')
    x = pkpy_oracle.unshare(x);
end

end
