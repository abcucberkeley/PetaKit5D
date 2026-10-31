function [res] = isolation_compare(dirA, dirB, varargin)
% compare the files under two folders bit for bit: tiff files byte for byte
% (and their pixels when the bytes differ), zarr arrays by voxel data and
% metadata (chunk files differ between identical runs), mat files by their
% variables (HDF5 files differ in their headers), json by content, png by
% bytes and pixels, pdf without its dates, other files byte for byte. Names
% are matched with uuids (8-4-4-4-12 hex digits) replaced by <uuid>.
%
% Parameters (as 'specifier'-value pairs):
%       ignore : cell array (default: {}). Regular expressions of relative paths to skip.
%   ignoreVars : cell array (default: {'t', 't1'}). Variables of mat files that are not compared:
%                the task times the block tasks save in their flag files. inputParser objects
%                (saved in stitchInfo files) are compared by their Results.
%
% Output: struct array with fields path (relative), kind, status ('equal',
%   'differ', 'onlyA' or 'onlyB') and detail.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dirA', @ischar);
ip.addRequired('dirB', @ischar);
ip.addParameter('ignore', {}, @iscell);
ip.addParameter('ignoreVars', {'t', 't1'}, @iscell);
ip.parse(dirA, dirB, varargin{:});

ignore = ip.Results.ignore;
ignoreVars = ip.Results.ignoreVars;
dirA = strip_slash(dirA);
dirB = strip_slash(dirB);
[filesA, zarrsA] = list_tree(dirA, ignore);
[filesB, zarrsB] = list_tree(dirB, ignore);

res = struct('path', {}, 'kind', {}, 'status', {}, 'detail', {});
res = compare_lists(res, dirA, dirB, filesA, filesB, 'file', ignoreVars);
res = compare_lists(res, dirA, dirB, zarrsA, zarrsB, 'zarr', ignoreVars);

end


function [res] = compare_lists(res, dirA, dirB, relA, relB, type, ignoreVars)

keyA = normalize(relA);
keyB = normalize(relB);
keys = unique([keyA; keyB]);
for i = 1 : numel(keys)
    ia = find(strcmp(keyA, keys{i}));
    ib = find(strcmp(keyB, keys{i}));
    n = max(numel(ia), numel(ib));
    for j = 1 : n
        r.path = keys{i};
        r.kind = type;
        r.detail = '';
        if j > numel(ia)
            r.status = 'onlyB';
            r.detail = relB{ib(j)};
        elseif j > numel(ib)
            r.status = 'onlyA';
            r.detail = relA{ia(j)};
        else
            fa = [dirA, '/', relA{ia(j)}];
            fb = [dirB, '/', relB{ib(j)}];
            if strcmp(type, 'zarr')
                [same, r.detail] = compare_zarr(fa, fb);
                r.kind = 'zarr';
            else
                [same, r.kind, r.detail] = compare_file(fa, fb, ignoreVars);
            end
            r.status = 'differ';
            if same
                r.status = 'equal';
            end
        end
        res(end + 1) = r; %#ok<AGROW>
    end
end

end


function [files, zarrs] = list_tree(root, ignore)
% relative paths of the files (not inside zarr arrays) and of the zarr arrays

files = {};
zarrs = {};
if ~exist(root, 'dir')
    return;
end
d = dir([root, '/**/*']);
d = d(~ismember({d.name}, {'.', '..'}));
rel = strrep(cellfun(@(f, n) [f, '/', n], {d.folder}, {d.name}, 'unif', 0), '\', '/');
rel = cellfun(@(x) x(numel(root) + 2 : end), rel, 'unif', 0);
isdir = [d.isdir];
inZarr = ~cellfun(@isempty, regexp(rel, '\.zarr/', 'once'));
zarrs = rel(isdir & endsWith(rel, '.zarr') & ~inZarr);
files = rel(~isdir & ~inZarr);
for i = 1 : numel(ignore)
    files = files(cellfun(@isempty, regexp(files, ignore{i}, 'once')));
    zarrs = zarrs(cellfun(@isempty, regexp(zarrs, ignore{i}, 'once')));
end
files = sort(files(:));
zarrs = sort(zarrs(:));

end


function [keys] = normalize(rel)

keys = regexprep(rel, '[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}', '<uuid>');

end


function [same, kind, detail] = compare_file(fa, fb, ignoreVars)

[~, ~, ext] = fileparts(fa);
kind = lower(ext);
detail = '';
ba = read_bytes(fa);
bb = read_bytes(fb);
switch kind
    case {'.tif', '.tiff'}
        same = isequal(ba, bb);
        if ~same
            try
                a = readtiff(fa);
                b = readtiff(fb);
                detail = pixel_detail(a, b, 'bytes differ');
            catch ME
                detail = ['bytes differ; ', ME.message];
            end
        end
    case '.mat'
        try
            a = load(fa);
            b = load(fb);
            [same, detail] = compare_structs(a, b, ignoreVars);
        catch ME
            same = isequal(ba, bb);
            detail = ME.message;
        end
    case '.json'
        try
            same = isequaln(jsondecode(char(ba)), jsondecode(char(bb)));
        catch
            same = isequal(ba, bb);
        end
    case '.png'
        same = isequal(ba, bb);
        if ~same
            try
                same = isequal(imread(fa), imread(fb));
                detail = 'bytes differ, pixels equal';
                if ~same
                    detail = 'pixels differ';
                end
            catch ME
                detail = ME.message;
            end
        end
    case '.pdf'
        pat = '/(CreationDate|ModDate) \([^)]*\)';
        same = strcmp(regexprep(char(ba), pat, ''), regexprep(char(bb), pat, ''));
    otherwise
        same = isequal(ba, bb);
end

end


function [same, detail] = compare_zarr(fa, fb)

detail = '';
try
    a = readzarr(fa);
    b = readzarr(fb);
    same = isequaln(a, b) && strcmp(class(a), class(b));
    if ~same
        detail = pixel_detail(a, b, '');
    end
    ma = zarr_meta(fa);
    mb = zarr_meta(fb);
    if ~isequaln(ma, mb)
        same = false;
        detail = strtrim([detail, ' metadata differs']);
    end
catch ME
    same = false;
    detail = ME.message;
end

end


function [meta] = zarr_meta(fn)

meta = [];
for name = {'.zarray', 'zarr.json'}
    if exist([fn, '/', name{1}], 'file')
        meta = jsondecode(fileread([fn, '/', name{1}]));
    end
end

end


function [same, detail] = compare_structs(a, b, ignoreVars)

detail = '';
fa = fieldnames(a);
fb = fieldnames(b);
same = isequal(sort(fa), sort(fb));
if ~same
    detail = sprintf('variables differ: %s', strjoin(setxor(fa, fb), ', '));
    return;
end
bad = {};
for i = 1 : numel(fa)
    va = a.(fa{i});
    vb = b.(fa{i});
    if any(strcmp(fa{i}, ignoreVars))
        continue;
    end
    % a handle object, equal only to itself: its parsed arguments
    if isa(va, 'inputParser') && isa(vb, 'inputParser')
        va = va.Results;
        vb = vb.Results;
    end
    if ~isequaln(va, vb)
        bad{end + 1} = fa{i}; %#ok<AGROW>
    end
end
same = isempty(bad);
if ~same
    detail = ['differ: ', strjoin(bad, ', ')];
end

end


function [detail] = pixel_detail(a, b, prefix)

if ~isequal(size(a), size(b)) || ~strcmp(class(a), class(b))
    detail = sprintf('%s; size or class differs (%s %s vs %s %s)', prefix, class(a), mat2str(size(a)), class(b), mat2str(size(b)));
else
    d = double(a(:)) - double(b(:));
    detail = sprintf('%s; %d of %d voxels differ, max abs %g', prefix, nnz(d), numel(d), max(abs(d)));
    if isequaln(a, b)
        detail = sprintf('%s; pixels equal', prefix);
    end
end
detail = regexprep(detail, '^; ', '');

end


function [b] = read_bytes(fn)

fid = fopen(fn, 'r');
b = fread(fid, inf, 'uint8=>uint8');
fclose(fid);

end


function [d] = strip_slash(d)

d = strrep(d, '\', '/');
while numel(d) > 1 && endsWith(d, '/')
    d = d(1 : end - 1);
end

end
