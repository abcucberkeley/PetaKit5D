function [x] = read_value(fn)
% read a value file written by petakit5d.testing.compare.write_pkv or write_value
%
% Placeholders come back as structs with field pkpy__ph (see
% fill_placeholders); references to other value files are read (and cached
% for the session, see load_ref).


fid = fopen(fn, 'r', 'l');
if fid < 0
    error('pkpy:oracle:read', 'cannot read %s', fn);
end
cleaner = onCleanup(@() fclose(fid));
m = fread(fid, 8, 'uint8=>uint8')';
if ~isequal(m, uint8(['PKV1', 0, 0, 0, 0]))
    error('pkpy:oracle:read', '%s is not a value file', fn);
end
fseek(fid, -16, 'eof');
n = double(typecast(fread(fid, 8, 'uint8=>uint8'), 'uint64'));
em = fread(fid, 8, 'uint8=>uint8')';
if ~isequal(em, uint8(['PKV1END', 0]))
    error('pkpy:oracle:read', '%s is incomplete', fn);
end
fseek(fid, -16 - n, 'eof');
hdr = jsondecode(native2unicode(fread(fid, n, 'uint8=>uint8')', 'UTF-8'));
x = dec(hdr.root, fid);

end


function [x] = dec(node, fid)

switch node.k
    case 'a'
        s = reshape(node.s, 1, []);
        cls = node.c;
        if node.n == 0
            x = empty_of(cls, s, node.x);
        else
            fseek(fid, node.o, 'bof');
            raw = fread(fid, node.n, 'uint8=>uint8');
            if node.x
                h = node.n / 2;
                x = complex(cast_raw(raw(1 : h), cls), cast_raw(raw(h + 1 : end), cls));
            else
                x = cast_raw(raw, cls);
            end
            x = reshape(x, s);
        end
    case 'cell'
        s = reshape(node.s, 1, []);
        x = cell(s);
        for i = 1 : prod(s)
            x{i} = dec(item(node.v, i), fid);
        end
    case 'struct'
        s = reshape(node.s, 1, []);
        f = node.f;
        if ~iscell(f)
            if ischar(f) && ~isempty(f)
                f = {f};
            else
                f = {};
            end
        end
        nf = numel(f);
        ne = prod(s);
        if nf == 0
            x = repmat(struct(), s);
        else
            C = cell(nf, ne);
            for e = 1 : ne
                for j = 1 : nf
                    C{j, e} = dec(item(node.v, (e - 1) * nf + j), fid);
                end
            end
            x = reshape(cell2struct(C, reshape(f, [], 1), 1), s);
        end
    case 'string'
        s = reshape(node.s, 1, []);
        x = strings(s);
        v = node.v;
        for i = 1 : prod(s)
            if iscell(v)
                e = v{i};
            elseif ischar(v)
                e = v;
            else
                e = [];
            end
            if ischar(e)
                x(i) = string(e);
            else
                x(i) = missing;
            end
        end
    case 'sparse'
        s = reshape(node.s, 1, []);
        x = sparse(dec(node.i, fid), dec(node.j, fid), dec(node.v, fid), s(1), s(2));
    case 'func'
        x = str2func(node.t);
    case 'map'
        k = dec(node.keys, fid);
        v = dec(node.values, fid);
        x = containers.Map('KeyType', node.kt, 'ValueType', node.vt);
        for i = 1 : numel(k)
            x(k{i}) = v{i};
        end
    case 'ph'
        x = struct('pkpy__ph', node.t, 'n', node.n, 'sub', node.sub);
    case 'ref'
        x = pkpy_oracle.load_ref(node.p);
    otherwise
        error('pkpy:oracle:read', 'value file: node kind %s cannot be read in MATLAB', node.k);
end

end


function [n] = item(v, i)
% jsondecode gives a struct array for nodes of one kind, else a cell array

if iscell(v)
    n = v{i};
else
    n = v(i);
end

end


function [x] = cast_raw(raw, cls)

switch cls
    case 'char'
        x = char(typecast(raw, 'uint16'));
    case 'logical'
        x = logical(raw);
    otherwise
        x = typecast(raw, cls);
end

end


function [x] = empty_of(cls, s, cx)

switch cls
    case 'char'
        x = char(zeros(s, 'uint16'));
    case 'logical'
        x = false(s);
    otherwise
        x = zeros(s, cls);
        if cx
            x = complex(x);
        end
end

end
