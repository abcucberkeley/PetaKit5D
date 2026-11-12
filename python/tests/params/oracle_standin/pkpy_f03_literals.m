function pkpy_f03_literals(corpusFile, outFile, probeFile)
% MATLAB side of the F03 literal corpus (a stand-in until the oracle harness of F04).
% For every corpus text: eval(text), str2num(text) and eval of the Python formatter's text,
% each described as {c: class, s: size, d: column-major data (hex for floats)}. With a
% third argument, also the inputParser probe (pkpy_f03_inputparser_probe).

c = jsondecode(fileread(corpusFile));
items = c.items;
if isstruct(items)
    items = num2cell(items);
end
out = cell(numel(items), 1);
for k = 1 : numel(items)
    it = items{k};
    r = struct();
    r.text = it.text;
    r.eval = pk_try_eval(it.text);
    [pk_x, pk_ok] = str2num(it.text); %#ok<ST2NM>
    r.str2num = struct('value', pk_describe(pk_x), 'ok', pk_ok);
    if isfield(it, 'formatted') && ~isempty(it.formatted)
        r.formatted_eval = pk_try_eval(it.formatted);
    end
    out{k} = r;
end
pk_write(outFile, out);
if nargin > 2
    pk_write(probeFile, pkpy_f03_inputparser_probe());
end
end


function pk_write(fn, val)
fid = fopen(fn, 'w');
fwrite(fid, jsonencode(val), 'char');
fclose(fid);
end


function d = pk_try_eval(pk_text__)
try
    pk_value__ = eval(pk_text__);
    d = pk_describe(pk_value__);
catch pk_err__
    d = struct('error', pk_err__.identifier, 'message', pk_err__.message);
end
end


function d = pk_describe(v)
d = struct();
if isstring(v)
    d.c = 'string';
    d.s = num2cell(size(v));
    d.d = cellstr(v(:)');
elseif iscell(v)
    d.c = 'cell';
    d.s = num2cell(size(v));
    d.d = cellfun(@pk_describe, v(:)', 'UniformOutput', false);
elseif isstruct(v)
    d.c = 'struct';
    d.s = num2cell(size(v));
    d.f = structfun(@pk_describe, v, 'UniformOutput', false);
elseif isnumeric(v) && ~isreal(v)
    d.c = [class(v), ' complex'];
    d.s = num2cell(size(v));
    d.d = {};
elseif isfloat(v)
    d.c = class(v);
    d.s = num2cell(size(v));
    d.d = arrayfun(@(t) num2hex(t), v(:)', 'UniformOutput', false);
elseif islogical(v) || isinteger(v)
    d.c = class(v);
    d.s = num2cell(size(v));
    d.d = num2cell(double(v(:)'));
elseif ischar(v)
    d.c = 'char';
    d.s = num2cell(size(v));
    d.d = {v(:)'};
else
    d.c = class(v);
    d.s = num2cell(size(v));
    d.d = {};
end
end
