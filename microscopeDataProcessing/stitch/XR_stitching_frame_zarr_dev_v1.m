function XR_stitching_frame_zarr_dev_v1(tileFullpaths, coordinates, varargin)
% zarr-based stitching pipeline. 
% 
%
% Author: Xiongtao Ruan (10/03/2020)
% xruan (10/04/2020): add distributed computing for block processing
% xruan (10/05/2020): use tile num to decide pad size 
% xruan (10/17/2020): add blurred blend option and add for border size for the blocks in stiched image.
% xruan (10/24/2020): add support for user-defined processing on tiles
% xruan (10/24/2020): add bbox crop option
% xruan (11/13/2020): use absolute shift for primary channel instead of relative shift 
%                     in xcorr (sometimes the relative positions changes for some time points).
% xruan (11/18/2020): use primary channel distance map in feather blend for other time points/channels.
% xruan (11/19/2020): add support for saving tiles separately in the
%                     corresponding locations in the stitched image.
% xruan (12/06/2020): add support for flipped tiles
% xruan (12/09/2020): add support for using primary coordinates for secondary channels/tps
% xruan (02/24/2021): add support for user defined xy, z max offsets for xcorr registration
% xruan (07/05/2021): add support for user defined resample (arbitary factor)
% xruan (07/14/2021): add support for the runs using stitchInfo, but with
% different images sizes (may happen after decon).
% xruan (07/21/2021): add support for skewed space stitching
% xruan (09/23/2021): add support for including partial files in skewed space stitch
% xruan (10/13/2021): add support for cropping data; set option to add
% offset; add support for skewed space stitching with reference (for decon data)
% xruan (10/25/2021): add support for a single distance map for all tiles in 
% feather blending (save time for computing).
% xruan (10/28/2021): add support for IO scan
% xruan (11/08/2021): add support for saving as hierachical multi-resolution dataset
% xruan (12/17/2021): add support for 2d stitching (e.g., MIPs), only with dsr for now
% xruan (01/25/2022): add support for loading tileFullpaths and coordinates from file.
% xruan (01/27/2022): change block size to be equal to the median of tile sizes if it
% is larger, to reduce the workload for each stitching block. 
% xruan (03/03/2022): fix bug for skewed space stitch coordinate conversion for z coordinate 
% xruan (06/20/2022): add input variable axisWeight for user defined weights for optimization
% xruan (08/25/2022): change CropToSize to tileOutBbox (more generic)
% xruan (11/19/2022): compute feather power within the distance map to save
% time for stitching processing
% xruan (11/25/2022): use 100 GB as threshold for big/small data. Use zstd
% and node factor 1 for big data, and use lz4 and node factor 2 for small
% data. Also use stitch block info from primary channel for secondary channels.
% xruan (12/13/2022): change xcorr thresh as user defined parameter with default 0.25


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('tileFullpaths', @iscell);
ip.addRequired('coordinates', @isnumeric);
ip.addParameter('ResultPath', '', @ischar);
ip.addParameter('tileInfoFullpath', '', @ischar); % matfile contains tileFullpaths and coordinates for too many tiles
ip.addParameter('stitchInfoDir', 'stitchInfo', @ischar);
ip.addParameter('stitchInfoFullpath', '', @ischar); % filename that contain stitch information for secondrary channels
ip.addParameter('ProcessedDirstr', '', @ischar); % path str for processed data 
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('SkewAngle', 32.45, @isscalar);
ip.addParameter('axisOrder', 'x,y,z', @ischar);
ip.addParameter('dataOrder', 'y,x,z', @ischar);
ip.addParameter('flippedTile', [], @(x) isempty(x) || all(islogical(x) | isnumeric(x)));
ip.addParameter('dz', 0.5, @isscalar);
ip.addParameter('xyPixelSize', 0.108, @isscalar);
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('IOScan', false, @islogical);
ip.addParameter('sCMOSCameraFlip', false, @islogical);
ip.addParameter('Reverse', false, @islogical);
ip.addParameter('Crop', false, @islogical);
ip.addParameter('InputBbox', [], @isnumeric); % crop input tile before processing
ip.addParameter('tileOutBbox', [], @isnumeric); % crop tile after processing
ip.addParameter('TileOffset', 0, @isnumeric); % offset added to the tile
ip.addParameter('df', [], @isnumeric);
ip.addParameter('save16bit', false , @islogical); % saves deskewed data as 16 bit -- not for quantification
ip.addParameter('EdgeArtifacts', 0, @isnumeric);
ip.addParameter('Decon', false, @islogical);
ip.addParameter('DS', false, @islogical);
ip.addParameter('DSR', false, @islogical);
ip.addParameter('resampleType', 'xy_isotropic', @ischar);
ip.addParameter('resample', [], @isnumeric);
ip.addParameter('deconRotate', false, @islogical);
ip.addParameter('BlendMethod', 'none', @ischar);
ip.addParameter('blendWeightDegree', 10, @isnumeric);
ip.addParameter('halfOrder', [3, 2, 1], @isnumeric);
ip.addParameter('overlapType', '', @ischar); % '', 'none', 'half', or 'full'
ip.addParameter('xcorrShift', true, @islogical);
ip.addParameter('isPrimaryCh', true, @islogical);
ip.addParameter('usePrimaryCoords', false, @islogical); % use primary coordinates for secondary channels/tps
ip.addParameter('stitchPadSize', [2, 2, 1], @(x) isnumeric(x) && numel(x) == 3);
ip.addParameter('padSize', [], @(x) isnumeric(x) && (isempty(x) || numel(x) == 3));
ip.addParameter('outBbox', [], @(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6));
ip.addParameter('zNormalize', false, @islogical);
ip.addParameter('xcorrDownsample', [2, 2, 1], @isnumeric); % y,x,z
ip.addParameter('xcorrThresh', 0.25, @isnumeric); % threshold of of xcorr, ignore shift if xcorr below this threshold.
ip.addParameter('xyMaxOffset', 300, @isnumeric); % max offsets in xy axes
ip.addParameter('zMaxOffset', 50, @isnumeric); % max offsets in z axis
ip.addParameter('shiftMethod', 'grid', @ischar); % {'local', 'global', 'grid', 'test', 'group'}
ip.addParameter('axisWeight', [1, 0.1, 10], @isnumeric); % axis weight for optimization, y, x, z
ip.addParameter('groupFile', '', @ischar); % file to define tile groups
ip.addParameter('singleDistMap', ~false, @islogical); % compute distance map for the first tile and apply to all other tiles
ip.addParameter('distBboxes', [], @isnumeric); % bounding boxes for distance transform
ip.addParameter('zarrFile', false, @islogical); 
ip.addParameter('largeZarr', false, @islogical); 
ip.addParameter('poolSize', [], @isnumeric); % max pooling size for large zarr MIPs
ip.addParameter('blockSize', [500, 500, 500], @isnumeric); 
ip.addParameter('batchSize', [500, 500, 500], @isnumeric); 
ip.addParameter('shardSize', [], @isnumeric); 
ip.addParameter('saveMultires', false, @islogical); % save as multi resolution dataset
ip.addParameter('resLevel', 4, @isnumeric); % downsample to 2^1-2^resLevel
ip.addParameter('BorderSize', [0, 0, 0], @isnumeric);
ip.addParameter('saveMIP', true , @islogical); % save MIP-z for stitch. 
ip.addParameter('tileIdx', [] , @isnumeric); % tile indices 
ip.addParameter('processFunPath', '', @(x) isempty(x) || iscell(x) || ischar(x)); % path of user-defined process function handle
ip.addParameter('stitchMIP', [], @(x) isempty(x)  || (islogical(x) && (numel(x) == 1 || numel(x) == 3))); % 1x3 vector or vector, by default, stitch MIP-z
ip.addParameter('stitch2D', false, @(x)islogical(x));  
ip.addParameter('bigStitchData', false, @(x)islogical(x));  
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 30, @isnumeric);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(tileFullpaths, coordinates, varargin{:});

pr = ip.Results;
ResultPath = pr.ResultPath;
tileInfoFullpath = pr.tileInfoFullpath;
stitchInfoDir = pr.stitchInfoDir;
stitchInfoFullpath = pr.stitchInfoFullpath;
SkewAngle = pr.SkewAngle;
flippedTile = pr.flippedTile;
axisOrder = pr.axisOrder;
dataOrder = pr.dataOrder;
dz = pr.dz;
px = pr.xyPixelSize;
objectiveScan = pr.objectiveScan;
IOScan = pr.IOScan;
InputBbox = pr.InputBbox;
tileOutBbox = pr.tileOutBbox;
TileOffset = pr.TileOffset;
ProcessedDirstr = pr.ProcessedDirstr;
DS = pr.DS;
DSR = pr.DSR;
BlendMethod = pr.BlendMethod;
blendWeightDegree = pr.blendWeightDegree;
halfOrder = pr.halfOrder;
overlapType = pr.overlapType;
xcorrShift = pr.xcorrShift;
shiftMethod = pr.shiftMethod;
axisWeight = pr.axisWeight;
groupFile = pr.groupFile;
isPrimaryCh = pr.isPrimaryCh;
usePrimaryCoords = pr.usePrimaryCoords;
stitchPadSize = pr.stitchPadSize;
outBbox = pr.outBbox;
xcorrDownsample = pr.xcorrDownsample;
xcorrThresh = pr.xcorrThresh;
xyMaxOffset = pr.xyMaxOffset;
zMaxOffset = pr.zMaxOffset;
singleDistMap = pr.singleDistMap;
distBboxes = pr.distBboxes;
saveMultires = pr.saveMultires;
resLevel = pr.resLevel;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
save16bit = pr.save16bit;
EdgeArtifacts = pr.EdgeArtifacts;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
poolSize = pr.poolSize;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
shardSize = pr.shardSize;
BorderSize = pr.BorderSize;
saveMIP = pr.saveMIP;
tileIdx = pr.tileIdx;
processFunPath = pr.processFunPath;
stitchMIP = pr.stitchMIP;
stitch2D = pr.stitch2D;
bigStitchData = pr.bigStitchData;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;
debug = pr.debug;

%  load tile paths and coordinates from file.
if ~isempty(tileInfoFullpath)
    a = load(tileInfoFullpath);
    tileFullpaths = a.tile_fullpaths;
    coordinates = a.xyz;
    tileIdx = a.tileIdx;
    flippedTile = a.flippedTile;
    clear a;
end

[dataPath, fsname_first] = fileparts(tileFullpaths{1});
if isempty(ResultPath)
    ResultPath = sprintf('%s/%s/%s.zarr', dataPath, 'matlab_stitch', fsname_first(1:end-21));
end

nF = numel(tileFullpaths);

% for 2d stitch, use the same worker to do all computing, to reduce the
% overhead of lauching other workers. 
if any(stitchMIP)
    if nF < 50
        parseCluster = false;
    end
end        

% check if a slurm-based computing cluster exist
if parseCluster 
    [parseCluster] = checkSlurmCluster(dataPath);
end

if isempty(uuid)
    uuid = get_uuid();
end

nv_fullname = ResultPath;
[pstr, nv_fsname] = fileparts(ResultPath);
[~, ResultDir] = fileparts(pstr);

% first check if the result already exist
% nv_fullname = sprintf('%s/%s/%s.zarr', dataPath, ResultDir, fsname_first(1:end-21));
if exist(nv_fullname, 'dir')
    fprintf('Stitched result %s is already exist, skip it!\n', nv_fullname);
    return;
end

fprintf('Stitch for tiles:\n%s\n', strjoin(tileFullpaths, '\n'));

if IOScan
    DSR = false;
end

% expand stitchMIP to 3d array, if more than 1 axis true, use the first one.
if any(stitchMIP)
    if numel(stitchMIP) == 1
        stitchMIP = [false, false, stitchMIP];
    end
    if sum(stitchMIP) > 1
        aind = find(stitchMIP, 1, 'first');
        stitchMIP = false(1, 3);
        stitchMIP(aind) = true;        
    end
end

% check if axis order is valid
axisOrder = strrep(axisOrder, ' ', '');
pattern = '^(-?x,?-?y,?-?z|-?y,?-?x,?-?z|-?z,?-?y,?-?x|-?x,?-?z,?-?y|-?x,?-?z,?-?y|-?y,?-?z,?-?x)$';
if ~regexpi(axisOrder, pattern)
    error("The axisOrder is not right, it must has the form like 'y,x,z' or '-x,y,z' (flipped in x-axis)!");
end

order_sign_mat = zeros(2, 3); % first row: order for permute; second row: sign
axisOrder_pure = strrep(strrep(axisOrder, ',', ''), '-', '');
xyz_str = 'xyz';
for i = 1 : 3
    [~, order_sign_mat(1, :)] = sort(axisOrder_pure);
    order_sign_mat(2, i) = 1 - 2 * contains(axisOrder, ['-', xyz_str(i)]);
end
xyz = coordinates(:, order_sign_mat(1, :)) .* order_sign_mat(2, :);

% normalize coordinates to zero
xyz = xyz - min(xyz, [], 1);
xyz_orig = xyz;

% for secondary channel/time point, use flippedTile from the primary channel. 
if ~isPrimaryCh
    if ~exist(stitchInfoFullpath, 'file')
        error('The stitch information filename %s does not exist!', stitchInfoFullpaths);
    end
    
    if isempty(flippedTile)
        a = load(stitchInfoFullpath, 'flippedTile');
        flippedTile = a.flippedTile;
    end

    if usePrimaryCoords
        a = load(stitchInfoFullpath, 'xyz_orig');
        xyz = a.xyz_orig;
    end
end

if ~isempty(tileIdx)
    % sort tiles based on tileIdx (handle zigzag orders)
    [tileIdx, sinds] = sortrows(tileIdx, [4, 3, 2, 1]);
    xyz = xyz(sinds, :);
    tileFullpaths = tileFullpaths(sinds);
elseif isempty(tileIdx)
    tileIdx = zeros(nF, 4);
end
tileNum = [numel(unique(tileIdx(:,1))), numel(unique(tileIdx(:,2))), numel(unique(tileIdx(:,3)))];
if isempty(tileNum)
    tileNum = zeros(1, nF);
end
tileNum = tileNum(order_sign_mat(1, :));

if objectiveScan || IOScan
    zAniso = dz/px;    
else
    zAniso = sind(SkewAngle)*dz/px;
end
theta = SkewAngle * pi/180;
% dx = cos(theta)*dz/xyPixelSize;
resample_type = pr.resampleType;
resample = pr.resample;
if ~isempty(resample)
    resample_type = 'given';
end

switch resample_type
    case 'given'
        xf = resample(1);
        yf = resample(2);
        zf = resample(3);
    case 'method_1'  %% old one
        zf = cot(abs(theta));
        yf = 1;
        % xf = cos(abs(theta)) + tan(abs(theta))*sin(abs(theta));
        xf = 1 / cos(abs(theta));
    case 'isotropic'  %% isotropic
        zf = 1;
        yf = 1;
        xf = 1;
    case 'xy_isotropic' %% x, y isotropic and z: r_z / r_x
        zf = sqrt((sin(theta) ^ 2 + zAniso ^ 2 * cos(theta) ^ 2) / (cos(theta) ^ 2 + zAniso ^ 2 * sin(theta) ^ 2));
        yf = 1;
        xf = 1;
end

% for skewed space stitch, match the resample factor with actual aspect ratio [1, 1, zAniso]. 
if ~DSR
    zf = zAniso * zf;
end

% map xyz ratio with right data order: convert to the default yxz data order
switch strrep(dataOrder, ',', '')
    case 'yzx'
    case 'zyx'
        xyz_f = [xf, yf, zf];
        xf = xyz_f(2);
        yf = xyz_f(3);
        zf = xyz_f(1);
    case 'xyz'
end

% create a text file to indicate pixel size
pixelInfoFname = sprintf('px%0.5g_py%0.5g_pz%0.5g', px*xf, px*yf, px*zf);
pixelInfoFullpath = sprintf('%s/%s/%s', dataPath, ResultDir, pixelInfoFname);
% check if there is some other pixel size file
dir_info = dir(sprintf('%s/%s/px*_py*_pz*', dataPath, ResultDir));
pixelFnames = {dir_info.name}';
for i = 1 : numel(pixelFnames)
    if ~strcmp(pixelFnames{i}, pixelInfoFname)
        delete([dir_info(i).folder, '/', pixelFnames{i}])
    end
end

if ~exist(pixelInfoFullpath, 'file')
    fclose(fopen(pixelInfoFullpath, 'w'));
end

% use single distance map for 
if ~DS && ~DSR && ~any(stitchMIP)
    singleDistMap = true;
end

if isempty(processFunPath) || all(cellfun(@isempty, processFunPath))
    processFunPath = {''};
    if EdgeArtifacts > 0
        usrFun = '';
        if DSR 
            if ~any(stitchMIP)
                usrFun = sprintf('@(x)erodeVolumeBy2DProjection(x,%d)', EdgeArtifacts);
            end
        else
            % for skewed space stitch add 1 count to avoid 0 in the image
            if TileOffset ~= 0
                usrFun = sprintf('@(x)erodeVolumeBy2DProjection(x+%d,%d)', TileOffset, EdgeArtifacts);                
            else
                usrFun = sprintf('@(x)erodeVolumeBy2DProjection(x,%d)', EdgeArtifacts);
            end
        end
        dt = char(datetime('now', 'Format', 'yyyyMMddHHmmSS'));
        fn = sprintf('%s/tmp/processFunction_%s.mat', pstr, dt);
        save('-v7.3', fn, 'usrFun');  
        if ~isempty(usrFun)
            processFunPath = {fn};        
        end
    end
end

% change stitch resample to [1, 1, 1] for DSR (because we resample DSR for the future).
if isempty(resample)
    resample = [1, 1, 1];
end
zarr_flippedTile = false(size(flippedTile)); 
stitchResample = [1, 1, 1];
if ~DSR
    stitchResample = resample;
    if isempty(ProcessedDirstr)
        zarr_flippedTile = flippedTile > 0;
    end
end
% check if there are partial files when converting tiff to zarr
partialFile = ~DS && ~DSR;

% process tile filenames based on different processing for tiles
processTiles = ~zarrFile || (zarrFile && (~isempty(InputBbox) || ~isempty(tileOutBbox) ...
    || any(~cellfun(@isempty, processFunPath)) || any(zarr_flippedTile) || any(stitchResample ~= 1)));

[inputFullpaths, zarrFullpaths, fsnames, zarrPathstr] = stitch_process_filenames( ...
    tileFullpaths, ProcessedDirstr, stitchMIP, stitchResample, zarrFile, processTiles);

% first check if it is 2d stitch
for f = 1 : nF
    imSize = getImageSize(inputFullpaths{f});
    if imSize(3) > 1
        stitch2D = false;
        break;
    end
    stitch2D = true;
end

if stitch2D && nF < 50
    parseCluster = false;
end

% check if total input size is greater than 100 GB if bigStitchData is false
nodeFactor = 2;
compressor = 'lz4';
compressor = 'zstd';
if (~stitch2D && ~bigStitchData && nF > 4 && prod(imSize) * nF * 4 > (100 * 2^30)) || largeZarr
    bigStitchData = true;
end
if bigStitchData
    nodeFactor = 1;
    compressor = 'zstd';
end

% convert tiff to zarr (if inputs are tiff tiles), and process tiles
locIds = tileIdx(:, 4);
stitch_process_tiles(inputFullpaths, 'resultDirName', zarrPathstr, 'zarrFile', zarrFile, ...
    'locIds', locIds, 'blockSize', round(blockSize / 2), 'shardSize', round(shardSize / 2), ...
    'flippedTile', zarr_flippedTile, 'resample', stitchResample, 'partialFile', partialFile, ...
    'InputBbox', InputBbox, 'tileOutBbox', tileOutBbox, 'processFunPath', processFunPath, ...
    'parseCluster', parseCluster, 'masterCompute', masterCompute, 'bigData', bigStitchData, ...
    'mccMode', mccMode, 'configFile', configFile);

% load all zarr headers as a cell array and get image size for all tiles
imSizes = zeros(nF, 3);
for i = 1 : nF 
    zarrFullpath = zarrFullpaths{i};    
    sz_i = getImageSize(zarrFullpath);
    
    if any(stitchMIP) && numel(sz_i) == 2
        imSizes(i, :) = [sz_i, 1];        
    else
        imSizes(i, :) = sz_i;
    end
end

if all(imSizes(:, 3) == 1)
    stitch2D = true;
end

dtype = getImageDataType(zarrFullpaths{1});
if save16bit
    dtype = 'uint16';
end 

% adjust the x coordinates for flipped tiles.
if ~isempty(flippedTile)
    flippedTile = flippedTile > 0;
    if DSR    
        offsize = (imSizes(:, 2) - 1) * xf - imSizes(:, 3) * cot(theta) * zf;
        % offsize = - imSizes(:, 3) * cot(theta) * zf;
        % offsize = 0;
        xyz(:, 1) = xyz(:, 1) - offsize .* flippedTile(:) * px * order_sign_mat(2, 1);
    else
        offsize = (imSizes(:, 3) - 1) * zf / sin(theta);
        xyz(:, 1) = xyz(:, 1) - offsize .* flippedTile(:) * px * order_sign_mat(2, 1);        
    end
end

if ~IOScan && ~objectiveScan && ~DS && ~DSR
    % convert coordinates in DSR space to skewned space
   xyz = [xyz(:, 3) / sin(theta), xyz(:, 2), -xyz(:, 1) * sin(theta) + xyz(:, 3) * cos(theta)];
end

if objectiveScan || DS
    % convert coordinates in Objective Scan / DS space
    xyz = [xyz(:, 1) * cos(theta) + xyz(:, 3) * sin(theta), xyz(:, 2), -xyz(:, 1) * sin(theta) + xyz(:, 3) * cos(theta)];
end    

% identify pairs between pairs
% first slight shift xyz such that the distance between pairs are in
% integer folds of resolution. 
xyz = round((xyz - min(xyz, [], 1)) ./ ([xf, yf, zf] * px)) .* ([xf, yf, zf] * px);

cuboids = [xyz, xyz + (imSizes(:, [2, 1, 3]) - 1) .* [xf, yf, zf] * px];
overlap_matrix = false(nF);
overlap_regions = zeros(nF * (nF - 1) / 2, 6);
for i = 1 : nF - 1
    cuboid_i = cuboids(i, :);
    cuboid_js = cuboids(i + 1 : end, :);
    
    [is_overlap_mat, cuboid_overlap_mat] = check_cuboids_overlaps(cuboid_i, cuboid_js, stitch2D);
    for j = i + 1 : nF
        if is_overlap_mat(j - i)
            overlap_matrix(i, j) = true;
            % ind = 0.5 * i * (2 * j - i - 1);
            ind = (i - 1) * nF - i * (i + 1) / 2 + j;
            overlap_regions(ind, :) = cuboid_overlap_mat(j - i, :);
        end
    end
end

dxyz_shift = zeros(nF, 3);
xyz_shift = xyz;
% calculate relative/absolute shifts between tiles
if xcorrShift && isPrimaryCh
    xcorrDir = sprintf('%s/%s/xcorr/%s/', dataPath, ResultDir, nv_fsname);
    xcorrFinalFn = sprintf('%s/%s/xcorr/%s.mat', dataPath, ResultDir, nv_fsname);
    if ~exist(xcorrFinalFn, 'file')
        xcorTmpFn = sprintf('%s/%s/xcorr/%s_%s.mat', dataPath, ResultDir, nv_fsname, uuid);
        assign_method = shiftMethod;
        
        MaxOffset = [xyMaxOffset, xyMaxOffset, zMaxOffset];
        [xyz_shift, dxyz_shift] = stitch_shift_assignment(zarrFullpaths, xcorrDir, imSizes, xyz, ...
            px, [xf, yf, zf], overlap_matrix, overlap_regions, MaxOffset, xcorrDownsample, ...
            xcorrThresh, tileIdx, assign_method, stitch2D, axisWeight, groupFile, largeZarr, ...
            poolSize, parseCluster, nodeFactor, mccMode, configFile);
        save('-v7.3', xcorTmpFn, 'xyz_shift', 'dxyz_shift');
        movefile(xcorTmpFn, xcorrFinalFn)
    else
        a = load(xcorrFinalFn, 'xyz_shift', 'dxyz_shift');
        xyz_shift = a.xyz_shift;
        dxyz_shift = a.dxyz_shift;
    end
elseif ~isPrimaryCh
    if ~exist(stitchInfoFullpath, 'file')
        error('The stitch information filename %s does not exist!', stitchInfoFullpaths);
    end
    
    a = load(stitchInfoFullpath, 'overlap_matrix', 'dxyz_shift', 'pImSz', 'pTileSizes', 'imdistFullpaths', 'imdistFileIdx');
    if ~all(a.overlap_matrix == overlap_matrix, 'all')
        warning('The overlap matrix of current overlap matrix is different from the one in the stitch info, use the common ones!')
        overlap_matrix = overlap_matrix & a.overlap_matrix;
    end
    dxyz_shift = a.dxyz_shift;
    xyz_shift = xyz + dxyz_shift;

    pImSz = a.pImSz;
    pTileSizes = a.pTileSizes;
    imdistFullpaths = a.imdistFullpaths;
    imdistFileIdx = a.imdistFileIdx;
end

% use half of overlap region for pairs of tiles with overlap
if isempty(overlapType)
    if strcmp(BlendMethod, 'none')
        overlapType = 'none';
    else
        overlapType = 'half';
    end
end

% to do: change to dictionary in the future
cuboids = [xyz_shift, xyz_shift + (imSizes(:, [2, 1, 3]) - 1) .* [xf, yf, zf] * px];

half_ol_region_cell = cell(nF * 8, 1);
ol_region_cell = cell(nF * 8, 1);
overlap_map_mat = zeros(nF * 8, 2);
counter = 1;
for i = 1 : nF - 1
    cuboid_i = cuboids(i, :);
    js = find(overlap_matrix(i, i + 1 : end)) + i;
    if isempty(js)
        continue;
    end
    cuboid_js = cuboids(js, :);
    
    % recheck if the overlapped tiles become not overlapped after shift
    % (or with large resample factors)    
    [is_overlap_mat] = check_cuboids_overlaps(cuboid_i, cuboid_js, stitch2D);
    for j1 = 1 : numel(js)
        j = js(j1);
        if ~is_overlap_mat(j1) 
            overlap_matrix(i, j) = false;
            continue;
        end
        cuboid_j = cuboids(j, :);
        ind = (j - 1) * nF + i;
        
        [mregion_1, mregion_2] = compute_half_of_overlap_region(cuboid_i, cuboid_j, ...
            px, [xf, yf, zf]', 'overlapType', overlapType, 'halfOrder', halfOrder, 'stitch2D', stitch2D);
        % half_ol_region_cell{i, j} = {mregion_1, mregion_2};
        half_ol_region_cell{counter} = {mregion_1, mregion_2};

        [mregion_11, mregion_21] = compute_half_of_overlap_region(cuboid_i, cuboid_j, ...
            px, [xf, yf, zf]', 'overlapType', 'zero', 'halfOrder', halfOrder, 'stitch2D', stitch2D);
        % ol_region_cell{i, j} = {mregion_11, mregion_21};
        ol_region_cell{counter} = {mregion_11, mregion_21};
        overlap_map_mat(counter, :) = [counter, ind];
        counter = counter + 1;
    end
end
half_ol_region_cell(counter : end) = [];
ol_region_cell(counter : end) = [];
overlap_map_mat(counter : end, :) = [];

% put tiles into the stitched image
% also pad the tiles such that the images of all time points and different
% channels are the same
if nF > 1 && xcorrShift
    pad_size = [xf, yf, zf] .* stitchPadSize .* tileNum([2, 1, 3]) * px;
else
    pad_size = [0, 0, 0];
end
xyz_shift_orig = xyz_shift;
origin = min(xyz_shift_orig, [], 1) - pad_size;
xyz_shift = xyz_shift_orig - origin;

if isPrimaryCh
    dxyz = max(xyz_shift_orig + imSizes(:, [2, 1, 3]) .* [xf, yf, zf] .* px) + pad_size - origin;
    nxs = round(dxyz(1)/(px*xf) + 2);
    nys = round(dxyz(2)/(px*yf) + 2);
    nzs = round(dxyz(3)/(px*zf) + 2);    
    if any(stitchMIP) || max(imSizes(:, 3) == 1)
        nzs = 1;       
    end
else    
    nys = pImSz(1);
    nxs = pImSz(2);
    nzs = pImSz(3);
end

% bouding box crop by redefining coordinate system and image size
if ~isempty(outBbox)
    bbox = outBbox;
    if any(isinf(bbox(4 : 6)))
        stchSz = [nys, nxs, nzs];
        bbox_end = bbox(4 : 6);
        bbox_end(isinf(bbox_end)) = stchSz(isinf(bbox_end));
        bbox(4 : 6) = bbox_end;
    end
    
    pad_size_l = pad_size - (bbox([2, 1, 3]) - 1) .* [xf, yf, zf] * px; 
    % pad_size_r = pad_size + (bbox([5, 4, 6]) - [nys, nxs, nzs]) .* [xf, yf, zf] * px; 
    origin = min(xyz_shift_orig, [], 1) - pad_size_l;
    % dxyz = max(xyz_shift, [], 1) + pad_size_r - origin;
    xyz_shift = xyz_shift_orig - origin;
    
    nys = bbox(4) - bbox(1) + 1;
    nxs = bbox(5) - bbox(2) + 1;
    nzs = bbox(6) - bbox(3) + 1;
end

int_xyz_shift = round(xyz_shift ./ ([xf, yf, zf] * px));

% to increase the scalability of the code, we use block-based processing
% first obtain coordinate data structure for blocks
batchSize = min([nys, nxs, nzs], batchSize);
batchSize = min(median(imSizes), batchSize);
blockSize = min(batchSize, blockSize);
batchSize = ceil(batchSize ./ blockSize) .* blockSize;
bSubSz = ceil([nys, nxs, nzs] ./ batchSize);
numBatches = prod(bSubSz);

% save block info (for record and distributed computing in the future)
stichInfoPath = [dataPath, '/', ResultDir, '/', stitchInfoDir];
if ~exist(stichInfoPath, 'dir')
    mkdir(stichInfoPath);
    fileattrib(stichInfoPath, '+w', 'g');
end

if numBatches < 30
    parseCluster = false;
end

taskRatio = max(0.1, prod([512, 512, 512]) / prod(batchSize));
if parseCluster
    taskSize = 20; % the number of blocks a job should process for [512, 512, 512]
    % keep task size inversely propotional to block size
    taskSize = max(1, round(taskSize * taskRatio));
    taskSize = max(taskSize, min(round(50 * taskRatio), ceil(numBatches / 5000)));
else
    taskSize = min(numBatches, round(500 * taskRatio));
end

nvSize = [nys, nxs, nzs];

[block_info_fullname, PerBlockInfoFullpaths, block_info_bytes] = stitch_process_block_info(int_xyz_shift, ...
    imSizes, nvSize, batchSize, overlap_matrix, ol_region_cell, half_ol_region_cell, ...
    overlap_map_mat, BorderSize, zarrFullpaths, stichInfoPath, nv_fsname, isPrimaryCh, ...
    stitchInfoFullpath=stitchInfoFullpath, stitch2D=stitch2D, uuid=uuid, taskSize=taskSize, ...
    parseCluster=parseCluster, mccMode=mccMode, configFile=configFile);

% initial stitched block image and save header in the disk
if ispc && numel(uuid) > 4
    % for PC, the path length limit is 260, so make it shorter in case of beyond the limit
    nv_tmp_fullname = sprintf('%s/%s/%s_nv_%s', dataPath, ResultDir, nv_fsname, uuid(1:4));        
else
    nv_tmp_fullname = sprintf('%s/%s/%s_nv_%s', dataPath, ResultDir, nv_fsname, uuid);    
end
if saveMultires
    mkdir(nv_tmp_fullname);
    py.zarr.open_group(nv_tmp_fullname, pyargs('mode', 'w'));
    
    nv_tmp_raw_fullname = sprintf('%s/L_1_1_1', nv_tmp_fullname);
else
    nv_tmp_raw_fullname = nv_tmp_fullname;
end

fresh_stitch = true;
if exist(nv_tmp_raw_fullname, 'dir')
    bim = blockedImage(nv_tmp_raw_fullname, 'Adapter', CZarrAdapter);
    nv_tmp_size = bim.Size;
    nv_tmp_block_size = bim.BlockSize;
    if numel(nvSize) == numel(nv_tmp_size) && all(nvSize == nv_tmp_size) && ...
            all(blockSize == nv_tmp_block_size)
        fresh_stitch = false;
    else
        rmdir(nv_tmp_raw_fullname, 's');
    end
end

% do not overwrite the intermediate data if it exists, this is good for
% large data when the master job runs out of time and new master job can
% continue the stitching
if fresh_stitch
    dimSeparator = '.';
    if prod(ceil([nys, nxs, nzs] ./ blockSize)) > 10000
        dimSeparator = '/';
    end
    createzarr(nv_tmp_raw_fullname, dataSize=[nys, nxs, nzs], blockSize=blockSize, ...
        shardSize=shardSize, dtype=dtype, compressor=compressor, dimSeparator=dimSeparator);
end

% add support for feather blending
stitchPath = [dataPath, '/', ResultDir, '/'];
if strcmpi(BlendMethod, 'feather') 
    % xruan disable singleDistMap if some tiles have different image sizes
    [uniq_locIds, uniq_inds] = unique(locIds);
    for i = 1 : numel(uniq_locIds)
        singleDistMap = singleDistMap && all(imSizes(locIds == uniq_locIds(i), :) == imSizes(uniq_inds(i), :), 'all');
        if ~singleDistMap
            break;
        end
    end
    if isPrimaryCh 
        imdistPath = [dataPath, '/', ResultDir, '/imdist/'];
        mkdir(imdistPath);
        [imdistFullpaths, imdistFileIdx] = compute_tile_distance_transform(block_info_fullname, ...
            stitchPath, zarrFullpaths, 'blendWeightDegree', blendWeightDegree, ...
            'singleDistMap', singleDistMap, 'locIds', locIds, 'distBboxes', distBboxes, ...
            'blockSize', round(blockSize/2), 'shardSize', round(shardSize/2), ...
            'compressor', compressor, 'largeZarr', largeZarr, 'poolSize', poolSize, ...
            'parseCluster', parseCluster, 'mccMode', mccMode, 'configFile', configFile);
    else
        usePrimaryDist = all(imSizes == pTileSizes, 'all');
        if ~usePrimaryDist
            imdistPath = [dataPath, '/', ResultDir, '/imdist/'];
            mkdir(imdistPath);
            [imdistFullpaths, imdistFileIdx] = compute_tile_distance_transform(block_info_fullname, stitchPath, ...
                zarrFullpaths, 'blendWeightDegree', blendWeightDegree, 'singleDistMap', singleDistMap, ...
                'locIds', locIds, 'distBboxes', distBboxes, 'blockSize', round(blockSize/2), ...
                'shardSize', round(shardSize/2), 'compressor', compressor, 'largeZarr', largeZarr, ...
                'poolSize', poolSize, 'parseCluster', parseCluster, 'mccMode', mccMode, ...
                'configFile', configFile);
        end
    end
else
    imdistFullpaths = {};
    imdistFileIdx = [];
end

% save stitch information for primary channels if there is xcorr shift
if isPrimaryCh 
    stitch_info_tmp_fullname = sprintf('%s/%s/%s/%s_%s.mat', dataPath, ResultDir, stitchInfoDir, nv_fsname, uuid);
    stitchInfoFullpath =sprintf('%s/%s/%s/%s.mat', dataPath, ResultDir, stitchInfoDir, nv_fsname);
    pImSz = [nys, nxs, nzs];
    pTileSizes = imSizes;
    save('-v7.3', stitch_info_tmp_fullname, 'ip', 'flippedTile', 'overlap_regions', ...
        'overlap_matrix', 'dxyz_shift', 'shiftMethod', 'pImSz', 'pTileSizes', 'xf', 'yf', 'zf', ...
        'px', 'PerBlockInfoFullpaths', 'imdistFullpaths', 'imdistFileIdx', 'xyz_orig');
    movefile(stitch_info_tmp_fullname, stitchInfoFullpath);
end

% for separate blend, save each tile separately in the proposed location in
% stitched image
if strcmpi(BlendMethod, 'separate')
    st_indices = round(xyz_shift ./ ([xf, yf, zf] * px));    
    processStitchSeparteTiles(zarrFullpaths, stitchPath, st_indices, imSizes, pImSz);
    nv_fullname = sprintf('%s/%s/%s.zarr', dataPath, ResultDir, nv_fsname);
    if ~exist(nv_fullname, 'dir')
        movefile(nv_tmp_fullname, nv_fullname);
    else
        rmdir(nv_tmp_fullname, 's');
    end
    return;
end

% process for each block based on all BlockInfo use distributed computing
% framework (change to use the slurm generic computing framework)
fprintf('Process blocks for stitched result...\n')
numTasks = ceil(numBatches / taskSize);

% flag dir for files, make sure there is no flags from old runs.
% wait more time zarr file computing
zarrFlagPath = sprintf('%s/%s/zarr_flags/%s_%s/', dataPath, ResultDir, nv_fsname, uuid);
if exist(zarrFlagPath, 'dir') && fresh_stitch
    rmdir(zarrFlagPath, 's')
end
mkdir(zarrFlagPath);

zarrFlagFullpaths = cell(numTasks, 1);
funcStrs = cell(numTasks, 1);

imdistFullpaths_str = sprintf("{'%s'}", strjoin(imdistFullpaths, ''','''));    
for t = 1 : numTasks
    batchInds = (t - 1) * taskSize + 1 : min(t * taskSize, numBatches);

    % save block info searately for each task for faster access
    % PerBlockInfoFullpath = sprintf('%s/stitch_block_info_blocks_%d_%d.mat', PerBlockInfoPath, blockInds(1), blockInds(end));
    PerBlockInfoFullpath = PerBlockInfoFullpaths{t};
    
    zarrFlagFullpaths{t} = sprintf('%s/blocks_%d_%d.mat', zarrFlagPath, batchInds(1), batchInds(end));
    funcStrs{t} = sprintf(['processStitchBlock([%s],''%s'',''%s'',''%s'',''%s'',[],[],', ...
        '''imSize'',[%s],''batchSize'',[%s],''dtype'',''%s'',''BlendMethod'',''%s'',', ...
        '''BorderSize'',[%s],''imdistFullpaths'',%s,''imdistFileIdx'',%s,', ...
        '''poolSize'',%s,''weightDegree'',%d)'], strrep(mat2str(batchInds), ' ', ','), ...
        block_info_fullname, PerBlockInfoFullpath, zarrFlagFullpaths{t}, ...
        nv_tmp_fullname, strrep(mat2str(nvSize), ' ', ','), strrep(mat2str(batchSize), ' ', ','), ...
        dtype, BlendMethod, strrep(mat2str(BorderSize), ' ', ','), imdistFullpaths_str, ...
        strrep(mat2str(imdistFileIdx), ' ', ','), strrep(mat2str(poolSize), ' ', ','), blendWeightDegree);
end

inputFullpaths = PerBlockInfoFullpaths;
outputFullpaths = zarrFlagFullpaths;

if parseCluster
    % reorder tasks to make time consuming task first
    [~, sinds] = sort(block_info_bytes, 'descend');
    inputFullpaths = inputFullpaths(sinds);
    outputFullpaths = outputFullpaths(sinds);
    funcStrs = funcStrs(sinds);
end

byte_num = dataTypeToByteNumber(dtype);

% cluster setting
cpusPerTask = 1 * nodeFactor;
memAllocate = prod(batchSize) * byte_num / 1024^3 * 20;
maxTrialNum = 2;
jobTimeLimit = taskSize * (3 / 60);

if ~exist(nv_fullname, 'dir') 
    is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'finalOutFullpath', nv_fullname, 'cpusPerTask', cpusPerTask, ...
        'memAllocate', memAllocate, 'jobTimeLimit', jobTimeLimit, 'maxTrialNum', maxTrialNum, ...
        'masterCompute', masterCompute, 'parseCluster', parseCluster, 'mccMode', mccMode, ...
        'configFile', configFile);
end

% retry with more resources and longer time
for i = 1 : 3
    if ~exist(nv_fullname, 'dir') && ~all(is_done_flag)
        is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'finalOutFullpath', nv_fullname, 'cpusPerTask', cpusPerTask * 2^i, ...
            'jobTimeLimit', jobTimeLimit * 2^i, 'maxTrialNum', maxTrialNum, ...
            'masterCompute', masterCompute, 'parseCluster', parseCluster, 'mccMode', mccMode, ...
            'configFile', configFile);
    end
end

if ~exist(nv_fullname, 'dir') && ~all(is_done_flag)
    if ~debug
        rmdir(nv_tmp_fullname, 's');
        rmdir(zarrFlagPath, 's');
    end
    if exist(nv_fullname, 'dir')
        fprintf('Stitched result %s is finished by another job, skip it!\n', nv_fullname);
        return;
    else
        error('The block processing for blocks [%s] cannot be finished!', num2str(find(~is_done_flag)'));
    end
end

if ~debug && exist(zarrFlagPath, 'dir')
    rmdir(zarrFlagPath, 's');
end

% create multiresolution zarr file
if saveMultires
    zarrFullpath = nv_tmp_raw_fullname;
    outputFullpath = nv_tmp_raw_fullname;
    XR_multiresZarrGeneration(zarrFullpath, outputFullpath, varargin)    
end

% nv_fullname = sprintf('%s/%s/%s.zarr', dataPath, ResultDir, fsname_first(1:end-21));
if exist(nv_tmp_fullname, 'dir')
    if ~exist(nv_fullname, 'dir')
        movefile(nv_tmp_fullname, nv_fullname);
    else
        if any(stitchMIP)
            rmdir(nv_fullname, 's');
            movefile(nv_tmp_fullname, nv_fullname);        
        else
            rmdir(nv_tmp_fullname, 's');
        end
    end
end

% save MIP
if saveMIP
    t0 = tic;    
    stcMIPPath = sprintf('%s/%s/MIPs/', dataPath, ResultDir);
    if ~exist(stcMIPPath, 'dir')
        mkdir(stcMIPPath);
        fileattrib(stcMIPPath, '+w', 'g');
    end
    stcMIPname = sprintf('%s%s_MIP_z.tif', stcMIPPath, nv_fsname);
    % for data greater than 250 GB, use the cluster based MIP.
    if prod([nys, nxs, nzs]) * byte_num / 2^30 < 250
        saveMIP_zarr(nv_fullname, stcMIPname, dtype, [1, 1, 1]);
    else
        XR_MIP_zarr(nv_fullname, axis=[1, 1, 1], parseCluster=parseCluster, ...
            mccMode=mccMode, configFile=configFile);
    end
    toc(t0);
end

if exist(nv_tmp_fullname, 'dir')  
    if ~debug
        rmdir(nv_tmp_fullname, 's');
    end
end

% convert to tiff
if ~false
%     write()
end


end
