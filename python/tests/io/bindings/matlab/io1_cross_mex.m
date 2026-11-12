function io1_cross_mex(folder)
% MATLAB side of the IO1 binding matrix: writes every case of folder/manifest.json
% with the mex wrappers PetaKit5D uses (writetiff, parallelWriteTiff, writezarr,
% createzarr, parallelWriteZarr), reads the files the Python packages wrote, and
% saves the checks to folder/results_mex.json. Run from an oracle tree after setup.

man = jsondecode(fileread(fullfile(folder, 'manifest.json')));
cases = man.cases;
if ~iscell(cases)
    cases = num2cell(cases);
end
nThreads = maxNumCompThreads;
recs = {};
t0 = tic;
for i = 1 : numel(cases)
    c = cases{i};
    img = readRaw(fullfile(folder, 'inputs', [c.name, '.bin']), c.mclass, c.shape(:)');
    if strcmp(c.kind, 'tiff')
        recs = [recs, tiffCase(folder, c, img, nThreads)]; %#ok<AGROW>
    else
        recs = [recs, zarrCase(folder, c, img, nThreads)]; %#ok<AGROW>
    end
end
recs = [recs{:}];
meta = struct('case', '_session', 'kind', 'meta', 'check', 'session', 'status', 'pass', ...
    'detail', sprintf('MATLAB %s; maxNumCompThreads %d; %s; %s; %.1f s', version, nThreads, ...
    which('parallelWriteTiff'), which('parallelWriteZarr'), toc(t0)));
recs = [meta, recs];
fid = fopen(fullfile(folder, 'results_mex.json'), 'w');
fwrite(fid, jsonencode(recs));
fclose(fid);
fprintf('%d checks, %d not passing\n', numel(recs), sum(~strcmp({recs.status}, 'pass')));
end


function recs = tiffCase(folder, c, img, nThreads)
recs = {};
name = c.name;
comp = c.compression;
fnMex = fullfile(folder, 'mex', [name, '.tif']);
fnMex1 = fullfile(folder, 'mex1', [name, '.tif']);
fnWt = fullfile(folder, 'wt', [name, '.tif']);
mkdirs(folder, {'mex', 'mex1', 'wt'});
recs{end + 1} = run(c, 'mex_write', @() parallelWriteTiff(fnMex, img, 'w', comp));
if ~isBig(c)
    maxNumCompThreads(1);
    recs{end + 1} = run(c, 'mex_write_1thread', @() parallelWriteTiff(fnMex1, img, 'w', comp));
    maxNumCompThreads(nThreads);
end
recs{end + 1} = run(c, 'writetiff_write', @() writetiff(img, fnWt, 'Compression', comp));
recs{end + 1} = check(c, 'mex_reads_mex', @() bitsEqual(readtiff(fnMex), img));
fnPy = fullfile(folder, 'py', [name, '.tif']);
recs{end + 1} = check(c, 'mex_reads_py', @() bitsEqual(readtiff(fnPy), img));
recs{end + 1} = check(c, 'mex_reads_pydirect', @() bitsEqual(readtiff(fullfile(folder, 'pydirect', [name, '.tif'])), img));
if ~isempty(c.range)
    r = c.range(:)';
    recs{end + 1} = check(c, 'mex_reads_py_range', @() bitsEqual(readtiff(fnPy, 'range', r), img(:, :, r(1) : r(2))));
end
sz = [size(img), ones(1, 3 - ndims(img))];
recs{end + 1} = check(c, 'mex_size_of_py', @() isequal(getImageSize(fnPy), sz));
end


function recs = zarrCase(folder, c, img, nThreads)
recs = {};
name = c.name;
mkdirs(folder, {'mex', 'mex1', 'mexonpy'});
fn = fullfile(folder, 'mex', [name, '.zarr']);
fn1 = fullfile(folder, 'mex1', [name, '.zarr']);
chunks = c.chunks(:)';
sz = c.shape(:)';
switch c.writer
    case {'regions', 'convert'}
        regs = cell(1, size(c.regions, 1));
        for k = 1 : numel(regs)
            b = c.regions(k, :);
            cls = c.mclass;
            if strcmp(c.writer, 'convert')
                cls = 'double';
            end
            regs{k} = readRaw(fullfile(folder, 'inputs', sprintf('%s_r%d.bin', name, k - 1)), cls, b(4 : 6) - b(1 : 3) + 1);
        end
        expected = readRaw(fullfile(folder, 'inputs', [name, '_expected.bin']), c.mclass, sz);
        create = @(f) createzarr(f, 'dataSize', sz, 'blockSize', chunks, 'dtype', c.mclass);
        recs{end + 1} = run(c, 'mex_create_base', @() create(fullfile(folder, 'mex', [name, '_base.zarr'])));
        recs{end + 1} = run(c, 'mex_regions_write', @() writeRegions(fn, create, regs, c.regions));
        maxNumCompThreads(1);
        recs{end + 1} = run(c, 'mex_regions_write_1thread', @() writeRegions(fn1, create, regs, c.regions));
        maxNumCompThreads(nThreads);
        recs{end + 1} = check(c, 'mex_regions_eq_expected', @() bitsEqual(readzarr(fn), expected));
        onPy = fullfile(folder, 'mexonpy', [name, '.zarr']);
        if exist(onPy, 'dir')
            rmdir(onPy, 's');
        end
        copyfile(fullfile(folder, 'py', [name, '_base.zarr']), onPy);
        recs{end + 1} = run(c, 'mex_regions_on_py_created_write', @() writeRegions(onPy, [], regs, c.regions));
        recs{end + 1} = check(c, 'mex_regions_on_py_created', @() bitsEqual(readzarr(onPy), expected));
        recs{end + 1} = check(c, 'mex_reads_py_regions', @() bitsEqual(readzarr(fullfile(folder, 'py', [name, '.zarr'])), expected));
        return;
    case 'createonly'
        w = @(f) createzarr(f, 'dataSize', sz, 'blockSize', chunks, 'dtype', c.mclass);
    case 'writezarr'
        w = @(f) writezarr(img, f, 'blockSize', chunks, 'dimSeparator', c.sep, 'sparseData', logical(c.sparse));
    case 'createzarr'
        w = @(f) createAndWrite(f, img, chunks, c);
    case 'direct'
        w = @(f) parallelWriteZarr(f, img, 'chunks', min(chunks, size(img, 1 : numel(chunks))), ...
            'dimension_separator', c.sep, 'order', c.order, 'cname', c.cname, 'clevel', 1, 'sparse', false);
    case 'shard'
    otherwise
        error('unknown writer %s', c.writer);
end
if isfield(c, 'shard') && ~isempty(c.shard)
    w = @(f) writezarr(img, f, 'blockSize', chunks, 'shardSize', c.shard(:)');
end
recs{end + 1} = run(c, 'mex_write', @() w(fn));
if ~isBig(c)
    maxNumCompThreads(1);
    recs{end + 1} = run(c, 'mex_write_1thread', @() w(fn1));
    maxNumCompThreads(nThreads);
end
if strcmp(c.writer, 'createonly')
    return;
end
fnPy = fullfile(folder, 'py', [name, '.zarr']);
recs{end + 1} = check(c, 'mex_reads_mex', @() bitsEqual(readzarr(fn), img));
recs{end + 1} = check(c, 'mex_reads_py', @() bitsEqual(readzarr(fnPy), img));
if exist(fullfile(folder, 'py', [name, '_as2d.zarr']), 'dir')
    recs{end + 1} = check(c, 'mex_reads_py_as2d', @() bitsEqual(readzarr(fullfile(folder, 'py', [name, '_as2d.zarr'])), img));
end
for k = 1 : size(c.bboxes, 1)
    b = c.bboxes(k, :);
    recs{end + 1} = check(c, sprintf('mex_reads_py_bbox%d', k), ...
        @() bitsEqual(readzarr(fnPy, 'inputBbox', b), img(b(1) : b(4), b(2) : b(5), b(3) : b(6))));
end
recs{end + 1} = check(c, 'mex_size_of_py', @() isequal(getImageSize(fnPy), [size(img), ones(1, 3 - ndims(img))]));
end


function createAndWrite(f, img, chunks, c)
createzarr(f, 'dataSize', size(img), 'blockSize', chunks, 'dtype', class(img), 'compressor', c.cname, ...
    'order', c.order, 'dimSeparator', c.sep);
writezarr(img, f);
end


function writeRegions(fn, create, regs, boxes)
if ~isempty(create)
    create(fn);
end
for k = 1 : numel(regs)
    writezarr(regs{k}, fn, 'bbox', boxes(k, :), 'create', false);
end
end


function r = run(c, name, f)
try
    f();
    r = rec(c, name, 'pass', '');
catch ME
    r = rec(c, name, 'error', [ME.identifier, ': ', ME.message]);
end
end


function r = check(c, name, f)
try
    if f()
        r = rec(c, name, 'pass', '');
    else
        r = rec(c, name, 'fail', '');
    end
catch ME
    r = rec(c, name, 'error', [ME.identifier, ': ', ME.message]);
end
end


function r = rec(c, name, status, detail)
r = struct('case', c.name, 'kind', c.kind, 'check', name, 'status', status, 'detail', detail);
end


function tf = isBig(c)
tf = isfield(c, 'big') && ~isempty(c.big) && c.big;
end


function tf = bitsEqual(a, b)
% same class, size and bits (NaN payloads and -0 included)
tf = strcmp(class(a), class(b)) && isequal(size(a), size(b));
if ~tf
    return;
end
if isfloat(a)
    u = 'uint64';
    if isa(a, 'single')
        u = 'uint32';
    end
    tf = isequal(typecast(a(:), u), typecast(b(:), u));
else
    tf = isequal(a, b);
end
end


function img = readRaw(fn, cls, sz)
fid = fopen(fn, 'r', 'ieee-le');
if strcmp(cls, 'logical')
    img = fread(fid, Inf, 'uint8=>logical');
else
    img = fread(fid, Inf, ['*', cls]);
end
fclose(fid);
img = reshape(img, [sz, ones(1, 2 - numel(sz))]);
end


function mkdirs(folder, subs)
for k = 1 : numel(subs)
    d = fullfile(folder, subs{k});
    if ~exist(d, 'dir')
        mkdir(d);
    end
end
end
