function [x] = unshare(x)
% x with arrays of its own: MATLAB copies an array's data only when a variable
% that shares it is changed, but a mex that writes into an input in place
% (indexing4d_mex, as feather_blending_3d's fallback calls it) changes every
% variable that shares the data. Assigning an element of each array gives it
% its own copy (the same bits); cells and structs are unshared element by
% element; a containers.Map (a handle: a case that adds a key changes it for
% every case) is copied; other values (empty and sparse arrays, objects,
% function handles) are left as they are.


if isa(x, 'containers.Map')
    y = containers.Map('KeyType', x.KeyType, 'ValueType', x.ValueType);
    k = keys(x);
    for i = 1 : numel(k)
        y(k{i}) = pkpy_oracle.unshare(x(k{i}));
    end
    x = y;
elseif iscell(x)
    for i = 1 : numel(x)
        x{i} = pkpy_oracle.unshare(x{i});
    end
elseif isstruct(x)
    f = fieldnames(x);
    for e = 1 : numel(x)
        for j = 1 : numel(f)
            x(e).(f{j}) = pkpy_oracle.unshare(x(e).(f{j}));
        end
    end
elseif (isnumeric(x) || islogical(x) || ischar(x)) && ~isempty(x) && ~issparse(x) && ~isobject(x)
    x(1) = x(1);
end

end
