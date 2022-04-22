function [dz_actual] = XR_estimate_actual_step_size_from_encoder(dataPath, varargin)
% estimate actual step size from encoder positions
% 
% xruan (01/27/2021): add support for bidirectional scan


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPath', @ischar);
ip.addParameter('CoordType', 'FPGA', @ischar); % FPGA or MCS2
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamB_ch0'}, @iscell);
ip.addParameter('pixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.2, @isnumeric);
ip.addParameter('angle', 32.45 , @isnumeric);
ip.addParameter('OnlyCoordinates', false, @islogical); % if true, only save encoder positions
ip.addParameter('xThresh', 0.3, @isnumeric); % percentage threshold to remove duplicate slices
ip.addParameter('skipMissingFile', true, @islogical); % skip the slice list of missing files
ip.addParameter('onlyFirstRow', false, @islogical); % only read the first row (coordinates)
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('cpusPerTask', 2, @isnumeric);
ip.addParameter('cpuOnlyNodes', true, @islogical);
ip.addParameter('uuid', '', @ischar);

ip.parse(dataPath, varargin{:});

tic
pr = ip.Results;
CoordType = pr.CoordType;
pixelSize = pr.pixelSize;
dz = pr.dz;
angle = pr.angle;
ChannelPatterns = pr.ChannelPatterns;
OnlyCoordinates = pr.OnlyCoordinates;
xThresh = pr.xThresh;
skipMissingFile = pr.skipMissingFile;
onlyFirstRow = pr.onlyFirstRow;

dataPath = [simplifyPath(dataPath), '/'];
resultPath = [dataPath, 'encoder_info/'];
if ~exist(resultPath, 'dir')
    mkdir(resultPath);
    fileattrib(resultPath, '+w', 'g');
end
save('-v7.3', [resultPath, '/parameters.mat'], 'pr');
resultPath = [simplifyPath(resultPath), '/'];


% load slice coordinates
fprintf('Process slice coordinates...\n');
dir_info = dir([dataPath, '/SliceList*csv']);

locFsns = {dir_info.name}';
% tokens = regexp(locFsns{1}, '(SliceList_.+part)\d+.csv', 'tokens');
tokens = regexp(locFsns, '(SliceList_.+part\d+).csv', 'tokens');
prefix = cellfun(@(x) x{1}{1}, tokens, 'unif', 0);
if onlyFirstRow
    locFn = sprintf('%s/%s.csv', dataPath, prefix{1});
    opts = detectImportOptions(locFn);
end

t_cell = cell(numel(locFsns), 1); 
for f = 1 : numel(locFsns)
    % locFn = sprintf('%s/%s%d.csv', dataPath, prefix, f - 1);
    locFn = sprintf('%s/%s.csv', dataPath, prefix{f});
    
    if onlyFirstRow
        fid = fopen(locFn, 'r');
        for i = 1 : 3
            fgetl(fid);
            tline = fgetl(fid);
            if ~isempty(tline) && (ischar(tline) || isstring(tline))          
                break;
            end
            pause(5);
            fid = fopen(locFn, 'r');
        end

        fclose(fid);

        tline = strsplit(tline, ',');
        tline(2 : end) = cellfun(@str2double, tline(2 : end), 'unif', 0);
        t = cell2table(tline, 'VariableName', opts.VariableNames);
    else
        t = readtable(locFn, 'Delimiter','comma');    
    end
    t_cell{f} = t;
end

t = cat(1, t_cell{:});

% save coordinate info to the disk 
t_fns = unique(t.Filename);
% channel patterns
include_flag = false(numel(t_fns), 1);
for c = 1 : numel(ChannelPatterns)
    include_flag = include_flag | contains(t_fns, ChannelPatterns{c}) | contains(t_fns, regexpPattern(ChannelPatterns{c}));
end

% exclude files that are not in the folder
if ~skipMissingFile
    for f = 1 : numel(t_fns)
        if include_flag(f) && ~exist([dataPath, t_fns{f}], 'file')
            include_flag(f) = false;
        end       
    end
end

t_fns = t_fns(include_flag);

% get time point/iteration, tile inds
tile_xyz = true;
if all(~cellfun(@isempty, regexp(t_fns{1}, '\d+x_\d+y_\d+z_\d+t', 'match')))
    expression = '_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t.tif';
elseif all(~cellfun(@isempty, regexp(t_fns{1}, 'Iter_\d+', 'match')))
    expression = '_Iter_(?<t>\d+)_';
    tile_xyz = false;
end
tmp = regexpi(t_fns, expression, 'names');

% skip those with incomplete z stacks
nz = sum(strcmp(t.Filename, t_fns{1}));
skipped_fns = {};

coords_struct = struct();
for f = 1 : numel(t_fns)
    fn = t_fns{f};
    inds_f = strcmp(t.Filename, fn);
    if sum(inds_f) < nz
        skipped_fns{end+1, 1} = fn;
        continue;
    end
    coords_struct(f).Filename = fn;
    coords_struct(f).ch = t.Channel(find(inds_f, 1));
    coords_struct(f).cam = t.Is_CamB_(find(inds_f, 1));
    
    coords_struct(f).xfpga = t.XFPGAPosition_um_(inds_f);
    coords_struct(f).yfpga = t.YFPGAPosition_um_(inds_f);
    coords_struct(f).zfpga = t.ZFPGAPosition_um_(inds_f);

    coords_struct(f).xmcs2 = t.XMCS2Position_um_(inds_f);
    coords_struct(f).ymcs2 = t.YMCS2Position_um_(inds_f);
    coords_struct(f).zmcs2 = t.ZMCS2Position_um_(inds_f);
    
    if strcmp(CoordType, 'FPGA')
        coords_struct(f).x = coords_struct(f).xfpga;
        coords_struct(f).y = coords_struct(f).yfpga;
        coords_struct(f).z = coords_struct(f).zfpga;
    else
        coords_struct(f).x = coords_struct(f).xmcs2;
        coords_struct(f).y = coords_struct(f).ymcs2;
        coords_struct(f).z = coords_struct(f).zmcs2;        
    end
    
    if tile_xyz
        coords_struct(f).tind = [str2double(tmp{f}.x), str2double(tmp{f}.y), str2double(tmp{f}.z)];
    else
        coords_struct(f).tind = [0, 0, 0];
    end
    coords_struct(f).t = str2double(tmp{f}.t);
end

coords_struct(arrayfun(@(x) isempty(coords_struct(x).x), 1 : numel(coords_struct))) = [];


% remove duplicate slices in the beginning and end
nF = numel(t_fns) - numel(skipped_fns);
duplicate_mat = true(nF, nz);
x_st = zeros(nF, 2);
for f = 1 : nF
    t_x = coords_struct(f).x;
    d_x = diff(t_x);
    d_x_inds = (dz > 0 & d_x < dz * xThresh) | (dz < 0 & d_x > dz * xThresh);
    s = find(d_x_inds(1 : round(nz / 4)), 1, 'last') + 1;
    if isempty(s)
        s = 1;
    end
    s = isempty(s) * 1 + (~isempty(s)) * s;
    t = find(d_x_inds(nz - round(nz / 4) + 1 : end), 1, 'first') + nz - round(nz / 4);
    if isempty(t)
        t = nz;
    end
    
    duplicate_mat(f, s : t) = false;
    x_st(f, :) = coords_struct(f).x([s, t]);
    coords_struct(f).duplicate = duplicate_mat(f, :);
end

% decide the common region to keep across time points
if dz > 0
    global_xrange = [max(x_st(:, 1)), min(x_st(:, 2))];
else
    global_xrange = [min(x_st(:, 1)), max(x_st(:, 2))];
end

% find the global minimum for each tile (computed with same tile indices)
% after removing the duplicates
encoded_tind = cat(1, coords_struct.tind) * [1e6, 1e3, 1]';
uniq_encoded_tind = unique(encoded_tind);

for i = 1 : numel(uniq_encoded_tind)
    tinds_i = find(encoded_tind == uniq_encoded_tind(i));
    
    x_i = cat(2, coords_struct(tinds_i).x);
    y_i = cat(2, coords_struct(tinds_i).y);
    z_i = cat(2, coords_struct(tinds_i).z);
    glbmin = [min(x_i(~duplicate_mat(tinds_i, :)')), min(y_i(~duplicate_mat(tinds_i, :)')), min(z_i(~duplicate_mat(tinds_i, :)'))];
    glbmax = [max(x_i(~duplicate_mat(tinds_i, :)')), max(y_i(~duplicate_mat(tinds_i, :)')), max(z_i(~duplicate_mat(tinds_i, :)'))];
    
    for j = 1 : numel(tinds_i)
        coords_struct(tinds_i(j)).glbmin = glbmin;
        coords_struct(tinds_i(j)).glbmax = glbmax;
    end
end

% estimate actual step size (dz) by excluding the top and bottome 10%
% coordinates
d_x_cell = cell(nF, 1);
for f = 1 : nF
    t_x = coords_struct(f).x;
    % duplicate = coords_struct(f).duplicate;
    % d_x = diff(t_x(~duplicate));
    d_x = diff(t_x);
    s = round(numel(d_x) * 0.10) + 1;
    t = round(numel(d_x) * 0.90);
    
    d_x_cell{f} = d_x(s : t);
end
d_x_all = cat(1, d_x_cell{:});
dz_actual = mean(abs(d_x_all));

% duplicate for all time points
global_duplicate = any(duplicate_mat, 1);
coords_info_fn = sprintf('%s%scoordinates_info.mat', resultPath, filesep);
coords_info_tmp_fn = sprintf('%s%scoordinates_info_%s.mat', resultPath, filesep, get_uuid());
save('-v7.3', coords_info_tmp_fn, 'coords_struct', 'skipped_fns', 'CoordType', ...
    'global_xrange', 'x_st', 'dz_actual', 'duplicate_mat', 'global_duplicate');
save(sprintf('%sdz_actual_%d.txt', resultPath, dz_actual), 'dz_actual', '-ASCII');
movefile(coords_info_tmp_fn, coords_info_fn);

end


