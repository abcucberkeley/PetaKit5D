function [x] = load_ref(fn)
% read an input value file, cached for the session (inputs are immutable:
% their names are their content hashes); the cache holds up to 2 GB


persistent cache bytes
if isempty(cache)
    cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    bytes = 0;
end
if isKey(cache, fn)
    x = cache(fn);
    return;
end
x = pkpy_oracle.read_value(fn);
info = whos('x');
if bytes + info.bytes > 2 ^ 31
    cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    bytes = 0;
end
cache(fn) = x;
bytes = bytes + info.bytes;

end
