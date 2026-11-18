function [] = XR_chromatic_shift_correction_data_wrapper(dataPaths, varargin)
% Chromatic shift correction with integer offset voxels, with either user provided offsets or estimated from PSFs.
%
% Required inputs:
%           dataPaths : char or cell array. Path(s) to the dataset(s) to be processed. Can be a single path as a char string or multiple paths in a cell array.
%
% Parameters (as 'specifier'-value pairs):
%     chromaticOffset : numeric array, [] or #channelx3 (default: []). The chromatic shift to apply, specified as [y, x, z] offsets in voxels, each row corresponds to a channel defined in channelPatterns.
%       resultDirName : char (default: 'Chromatic_Shift_Corrected'). Name of the subdirectory within each dataPath where results will be saved.
%                mode : char (default: 'valid'). Determines the output size. Options are 'valid' (only the overlapping region), 'same' (same size as input), or 'full' (large enough to contain the entire shifted image), in any case.
%            padValue : numeric scalar (default: 0). Value used for padding areas with no image data after shifting.
%           newOrigin : numeric vector (default: []). Resets the origin coordinates, applicable only when 'mode' is 'same'.
%     channelPatterns : cell array (default: {'CamA', 'CamB'}). File extensions or patterns used to identify image files for each channel.
%        psfFullpaths : cell array (default: {}). Full paths to PSF files, if needed for the process, same order as channelPatterns.
%           maxOffset : 1x3 numeric vector (default: [20, 20, 20]). The maximum expected offset [y, x, z] for chromatic measurement.
%          cropLength : 1x3 numeric vector (default: [0, 0, 0]). The amount [y, x, z] to crop from the edges for chromatic measurement.
%            zarrFile : true|false (default: false). Specifies if the input data is in Zarr format.
%           largeFile : true|false (default: false). Flag to indicate if the input files are large and require special handling: Zarr files are cropped in batches (XR_crop_zarr, which pads with 0) and saved as Zarr; tiff files are cropped in memory.
%            saveZar  : true|false (default: false). If true, the output is saved in Zarr format; otherwise, saved in Tiff format.
%           batchSize : 1x3 numeric vector (default: [1024, 1024, 1024]). The size of data batches [y, x, z] for processing.
%           blockSize : 1x3 numeric vector (default: [256, 256, 256]). The chunk/block size [y, x, z] for writing data, especially for Zarr.
%        parseCluster : true|false (default: true). If true, sets up the parallel computing cluster environment.
%       masterCompute : true|false (default: true). If true, the master node will also participate in computation tasks.
%           jobLogDir : char (default: '../job_logs'). Directory for saving logs from cluster jobs.
%         cpusPerTask : numeric scalar (default: 3). The number of CPU cores to allocate for each parallel task.
%          configFile : char (default: ''). Path to an optional configuration file for slurm job submission settings.
%             mccMode : true|false (default: false). Set to true if running the code as a compiled application (using MATLAB Compiler).
%                uuid : char (default: ''). A unique identifier string for this processing run, often used for temporary files.
%               debug : true|false (default: false). If true, enables debug mode which may provide more verbose output.
% 
% Author: Xiongtao Ruan (09/19/2025)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('chromaticOffset', [], @(x) isnumeric(x)); % y, x, z in voxels
ip.addParameter('resultDirName', 'Chromatic_Shift_Corrected', @ischar); 
ip.addParameter('mode', 'valid', @(x) ischar(x) && any(strcmpi(x, {'valid', 'same', 'full'}))); % same, valid or full
ip.addParameter('padValue', 0, @isnumeric);
ip.addParameter('newOrigin', [], @isnumeric); % reset origin (only for same mode).
ip.addParameter('channelPatterns', {'CamA', 'CamB'}, @iscell);
ip.addParameter('psfFullpaths', {}, @iscell);
ip.addParameter('maxOffset', [20, 20, 20], @isnumeric); % max offset across channels in voxel in y, x, z
ip.addParameter('cropLength', [0, 0, 0], @isnumeric); % max offset across channels in voxel in y, x, z
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('saveZarr', false, @islogical);
ip.addParameter('batchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256] , @isvector); % in y, x, z
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 3, @isnumeric);
ip.addParameter('configFile', '', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
chromaticOffset = pr.chromaticOffset;
resultDirName = pr.resultDirName;
% (another mode cropped each channel's box shifted against its offset)
mode = lower(pr.mode);
padValue = pr.padValue;
newOrigin = pr.newOrigin;
channelPatterns = pr.channelPatterns;
psfFullpaths = pr.psfFullpaths;
maxOffset = pr.maxOffset;
cropLength = pr.cropLength;
zarrFile = pr.zarrFile;
largeFile = pr.largeFile;
saveZarr = pr.saveZarr;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
configFile = pr.configFile;
mccMode = pr.mccMode;
uuid = pr.uuid;

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

% (large Zarr files are cropped in batches, into a Zarr result)
largeZarr = zarrFile && largeFile;
if largeZarr
    saveZarr = true;
end

nd = numel(dataPaths);
resultPaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    if ~strcmp(dataPath(end), filesep)
        dataPaths{d} = [dataPath, filesep];
    end

    dataPaths{d} = [simplifyPath(dataPaths{d}), '/'];
    resultPath = [dataPaths{d}, '/' resultDirName, '/'];
    resultPaths{d} = resultPath;
    mkdir(resultPath);
    fileattrib(resultPath, '+w', 'g');
end

if isempty(uuid)
    uuid = get_uuid();
end

% parse image filenames
[fnames, fsns, fd_inds, filepaths, ch_inds] = parseImageFilenames(dataPaths, zarrFile, channelPatterns);
nF = numel(fnames);
nc = numel(channelPatterns);

uniq_ch = unique(ch_inds);
if numel(uniq_ch) ~= nc
    error('Number of channels of the images does not match the number of channel patterns!');
end

% measure chromatic shift using PSF files
if isempty(chromaticOffset) && ~isempty(psfFullpaths)
    if numel(psfFullpaths) ~= nc
        error('The number of PSF files does not match the number of the channels!')
    end
    fprintf('Measure chromatic offsets across channels with PSF files...\n');
    disp(psfFullpaths(:));

    measureResultPath = [resultPaths{1}, 'measurement/'];
    mkdir(measureResultPath);
    measurementFullpath = sprintf('%s/chromatic_shift_measurement_result.mat', measureResultPath);
    % (the PSFs are tiff files: zarrFile was set false for them and stayed
    % false for the data, so zarr data failed in the crops)
    [relative_shift_yxz_mat] = XR_chromatic_shift_estimation(psfFullpaths, ...
        resultFullpath=measurementFullpath, maxOffset=maxOffset, cropLength=cropLength, ...
        zarrFile=false, uuid=uuid);
    
    % unify shifts across channels and generate bboxes
    chromaticOffset = relative_shift_yxz_mat - min(relative_shift_yxz_mat, [], 1);
end

% (without offsets or PSFs the boxes failed on a size mismatch, with too few
% rows on an index; more rows changed the overlap)
if isempty(chromaticOffset)
    error('XR_chromatic_shift_correction_data_wrapper:noOffsets', ...
        'Give chromaticOffset (a [y, x, z] row per channel pattern) or psfFullpaths to measure it.');
end
if ~isequal(size(chromaticOffset), [nc, 3])
    error('XR_chromatic_shift_correction_data_wrapper:offsetSize', ...
        'chromaticOffset must have a [y, x, z] row per channel pattern (%d x 3), not %s.', nc, ...
        mat2str(size(chromaticOffset)));
end

disp('Chromatic offset across channels:');
disp(chromaticOffset);

% assume the data sizes across channels within each data folder are the same.
imSizes = zeros(nd, 3);
ch_bboxes_cell = cell(nd, 1);
pad = true;
for d = 1 : nd
    fd_ind_d = fd_inds(find(fd_inds == d, 1, 'first'));

    fn_d = [dataPaths{d}, fnames{fd_ind_d}];
    imSize = getImageSize(fn_d);
    imSizes(d, :) = imSize;

    ch_bboxes_d = [1, 1, 1, imSize] + [chromaticOffset, chromaticOffset];
    
    switch mode
        case 'valid'
            overlap_bbox = [max(ch_bboxes_d(:, 1 : 3), [], 1),  min(ch_bboxes_d(:, 4 : 6), [], 1)];
            ch_bboxes_d = overlap_bbox - [ch_bboxes_d(:, 1 : 3), ch_bboxes_d(:, 1 : 3)] + 1;
            pad = false;
        case 'same'
            origin_d = newOrigin;
            if isempty(origin_d)
                % center the output on the overlap of the channels
                overlap_bbox = [max(ch_bboxes_d(:, 1 : 3), [], 1),  min(ch_bboxes_d(:, 4 : 6), [], 1)];
                overlapSize = overlap_bbox(4 : 6) - overlap_bbox(1 : 3) + 1;
                origin_d = overlap_bbox(1 : 3) + floor((overlapSize - imSize) / 2);
            end
            % the output region in the coordinates of each channel
            ch_bboxes_d = [origin_d - ch_bboxes_d(:, 1 : 3) + 1, origin_d - ch_bboxes_d(:, 1 : 3) + imSize];
        case 'full'
            origin_d = min(ch_bboxes_d(:, 1 : 3), [], 1);
            fullSize = max(ch_bboxes_d(:, 4 : 6), [], 1) - origin_d + 1;
            ch_bboxes_d = [origin_d - ch_bboxes_d(:, 1 : 3) + 1, origin_d - ch_bboxes_d(:, 1 : 3) + fullSize];
    end
    ch_bboxes_cell{d} = ch_bboxes_d;
end

if largeZarr && pad && padValue ~= 0
    error('XR_chromatic_shift_correction_data_wrapper:padValue', ...
        'Large Zarr files are padded with 0 (XR_crop_zarr): padValue %g is not supported.', padValue);
end

ext = '.tif';
if saveZarr
    ext = '.zarr';
end

frameFullpaths = cell(nF, 1);
resultFullpaths = cell(nF, 1);
func_strs = cell(nF, 1);
for f = 1 : nF
    frameFullpath_f = sprintf('%s/%s', dataPaths{fd_inds(f)}, fnames{f});
    frameFullpaths{f} = frameFullpath_f;
    resultFullpaths{f} = sprintf('%s/%s%s', resultPaths{fd_inds(f)}, fsns{f}, ext);

    bbox_f = ch_bboxes_cell{fd_inds(f)}(ch_inds(f), :);

    func_strs{f} = sprintf(['XR_crop_frame(''%s'',''%s'',%s,''pad'',%s,''padValue'',%d,', ...
        '''zarrFile'',%s,''saveZarr'',%s,''batchSize'',%s,''blockSize'',%s,', ...
        '''uuid'',''%s'',''parseCluster'',%s,''masterCompute'',%s,''cpusPerTask'',%d,', ...
        '''mccMode'',%s,''configFile'',''%s'')'], frameFullpath_f, resultFullpaths{f}, ...
        mat2str_comma(bbox_f), string(pad), padValue, string(zarrFile), string(saveZarr), ...
        mat2str_comma(batchSize), mat2str_comma(blockSize), uuid, string(parseCluster), ...
        string(masterCompute), cpusPerTask, string(mccMode), configFile);
    % (largeFile only sized the jobs: the crops were made in memory with the
    % memory of a batch; passed only then, so that compiled masters built
    % before still run the other tasks)
    if largeZarr
        func_strs{f} = sprintf('%s,''largeFile'',true)', func_strs{f}(1 : end - 1));
    end
end

jobCPUs = cpusPerTask;
coordNice = 0;
if largeZarr
    % each file's job (XR_crop_zarr) runs its batches as jobs: it holds MATLAB
    % and the merge of the MIP slabs, and if it crops batches itself
    % (masterCompute) cpusPerTask and a batch's memory, as in XR_crop_dataset
    boxes = cat(1, ch_bboxes_cell{:});
    outSize = max(boxes(:, 4 : 6) - boxes(:, 1 : 3) + 1, [], 1);
    batchSize_b = min(outSize, batchSize);
    dtype = getImageDataType(frameFullpaths{1});
    mergeMem = MIP_zarr_merge_memory(outSize, batchSize_b, min(outSize, blockSize), dtype);
    [jobCPUs, memAllocate] = coordinator_resources(mergeMem, ...
        cluster_config_value(configFile, 'masterCompute', masterCompute), cpusPerTask, prod(batchSize_b) * 4 * 2 / 1024^3);
    % (and they wait for their batch jobs: these go first)
    if parseCluster
        coordNice = coordinator_nice();
    end
else
    % (a large tiff file is cropped in memory)
    memAllocate = prod(max(imSizes, [], 1)) / 1024^3 * 4 * 2.5;
end

is_done_flag = generic_computing_frameworks_wrapper(frameFullpaths, resultFullpaths, func_strs, ...
    parseCluster=parseCluster, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
    cpusPerTask=jobCPUs, memAllocate=memAllocate, mccMode=mccMode, configFile=configFile, nice=coordNice);

if ~all(is_done_flag)
    is_done_flag = generic_computing_frameworks_wrapper(frameFullpaths, resultFullpaths, func_strs, ...
        parseCluster=parseCluster, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
        cpusPerTask=jobCPUs, memAllocate=memAllocate * 2, mccMode=mccMode, configFile=configFile, nice=coordNice);
end

if ~all(is_done_flag)
    warning('Some files are not processed!');
end

end