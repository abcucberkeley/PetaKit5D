function [] = write_value(fn, x)
% write a MATLAB value to a value file (format: python/docs/testing.md)
%
% Arrays keep their class, size and every bit (column-major, little-endian;
% complex: the real block, then the imaginary block). Cells, structs (field
% order kept), strings, sparse matrices, function handles (as text),
% containers.Map and inputParser objects (their Results) are kept; other
% objects are written as their class name and display text.


fid = fopen(fn, 'w', 'l');
if fid < 0
    error('pkpy:oracle:write', 'cannot write %s', fn);
end
cleaner = onCleanup(@() close_quietly(fid));
fwrite(fid, uint8(['PKV1', 0, 0, 0, 0]), 'uint8');
root = enc(x, fid);
js = unicode2native(jsonencode(struct('pkv', 1, 'root', root)), 'UTF-8');
fwrite(fid, js, 'uint8');
fwrite(fid, typecast(uint64(numel(js)), 'uint8'), 'uint8');
fwrite(fid, uint8(['PKV1END', 0]), 'uint8');
clear cleaner;

end


function [node] = enc(x, fid)

if isa(x, 'gpuArray')
    x = gather(x);
end
if issparse(x)
    [i, j, v] = find(x);
    node = struct('k', 'sparse', 'c', class(x), 's', size(x), 'i', enc(double(i), fid), ...
        'j', enc(double(j), fid), 'v', enc(v, fid));
elseif (isnumeric(x) || islogical(x) || ischar(x)) && ~isobject(x)
    cx = isnumeric(x) && ~isreal(x);
    if cx && isinteger(x)
        node = struct('k', 'object', 'c', ['complex ', class(x)], ...
            'v', enc(struct('re', real(x), 'im', imag(x)), fid), 't', '');
        return;
    end
    o = align(fid);
    n = write_raw(fid, x, cx);
    node = struct('k', 'a', 'c', class(x), 's', size(x), 'o', o, 'n', n, 'x', cx);
elseif iscell(x)
    v = cell(1, numel(x));
    for i = 1 : numel(x)
        v{i} = enc(x{i}, fid);
    end
    node = struct('k', 'cell', 's', size(x), 'v', {v});
elseif isstruct(x)
    f = reshape(fieldnames(x), 1, []);
    nf = numel(f);
    v = cell(1, numel(x) * nf);
    for e = 1 : numel(x)
        for j = 1 : nf
            v{(e - 1) * nf + j} = enc(x(e).(f{j}), fid);
        end
    end
    node = struct('k', 'struct', 's', size(x), 'f', {f}, 'v', {v});
elseif isstring(x)
    v = cell(1, numel(x));
    m = ismissing(x);
    for i = 1 : numel(x)
        if m(i)
            v{i} = [];
        else
            v{i} = char(x(i));
        end
    end
    node = struct('k', 'string', 's', size(x), 'v', {v});
elseif isa(x, 'function_handle')
    node = struct('k', 'func', 't', func2str(x));
elseif isa(x, 'containers.Map')
    node = struct('k', 'map', 'kt', x.KeyType, 'vt', x.ValueType, 'keys', enc(keys(x), fid), ...
        'values', enc(values(x), fid));
elseif isa(x, 'inputParser')
    node = struct('k', 'object', 'c', 'inputParser', 'v', enc(x.Results, fid), 't', '');
else
    try
        t = evalc('disp(x)');
    catch
        t = '';
    end
    node = struct('k', 'object', 'c', class(x), 'v', [], 't', t(1 : min(end, 10000)));
end

end


function [o] = align(fid)

o = ftell(fid);
pad = mod(-o, 8);
if pad
    fwrite(fid, zeros(1, pad, 'uint8'), 'uint8');
end
o = o + pad;

end


function [n] = write_raw(fid, x, cx)
% raw bits in column-major order, in blocks (typecast copies)

if ischar(x)
    x = uint16(x);
elseif islogical(x)
    x = uint8(x);
end
if cx
    parts = {real(x), imag(x)};
else
    parts = {x};
end
n = 0;
step = 2 ^ 24;
for p = 1 : numel(parts)
    v = parts{p};
    for s = 1 : step : numel(v)
        e = min(numel(v), s + step - 1);
        b = typecast(reshape(v(s : e), [], 1), 'uint8');
        fwrite(fid, b, 'uint8');
        n = n + numel(b);
    end
end

end


function [] = close_quietly(fid)

try
    fclose(fid);
catch
end

end
