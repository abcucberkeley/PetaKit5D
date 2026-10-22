function [] = XR_decon_data_wrapper(dataPaths, varargin)
% Dataset level deconvolution wrapper.
% 
%
% Required inputs :   
%          dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings and PSFs. 
%
% Parameters (as 'specifier'-value pairs): 
%      resultDirName : char (default: 'matlab_decon'). Result directory under data paths.
%          overwrite : true|false or 1x3 bool vector (default: false). Overwrite existing results
%    channelPatterns : a cell array (default: {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}).  Channel identifiers for included channels. 
%          skewAngle : a number (default: 32.45). Skew angle (in degree) of the stage.
%                 dz : a number (default: 0.5). Scan interval in um.
%        xyPixelSize : a number (default: 0.108). Pixel size in um.
%          save16bit : 1x1 or 1x2 bool vector (default: [true, true]). Save 16bit result for deskew/rotate and stitch. 
%   parseSettingFile : true|false (default: false). Use the setting file to decide whether filp z stacks or not.
%         flipZstack : true|false (default: false). Flip z stacks.
%         background : empty or a number or 1x#channel (default: []). Background subtraction before deconvolution. If empty, set as 100.
%              dzPSF : a number (default: 0.5). Scan interval in um for PSF.
%        edgeErosion : a number (default: 0). Number of pixels to erode after deconvolution to remove edge artifacts.
%         erodeByFTP : true|false (default: true). Apply the erosion mask from the first time point.
%       psfFullpaths : char or a cell array (default: {'','',''}). Paths of PSF with the same order as channelPatterns.
%          deconIter : a number or 1x#dataset (default: 15). Deconvolution iteration. If it is a vector with 1x#dateset, apply specific iterations for datasets.
%           RLMethod : 'original'|'simplified'|'omw' (default: 'simplified'). 'original': matlab deconvlucy.m (no GPU support); 'simplified': our implemention based on deconvlucy.m with GPU support; 'omw': OMW deconvolution.
%        wienerAlpha : a number or 1x#channel (default: 0.005). Wiener parameter. If it is a 1x#channel vector, apply a specific parameter per channel mapping to channelPatterns.
%       OTFCumThresh : a number or 1x#channel (default: 0.9). Accumulated energy quantile to segment OTF mask. Smaller number means tighter threshold with smaller masks. 
%      hannWinBounds : 1x2 vector (default: [0.8, 1]). Apodization range for distance matrix with Hann window function. 0 means the center of support, and 1 means the edge of support. 
%             skewed : empty or true|false (default: []). Parameter for OTF segmentation based on PSF space. If true, PSF is in the skewed space. If empty, automatically deteriming if PSF is in skewed space.
%              debug : true|false (default: false). Debug mode for simplified and OMW method. If true, save the intermediate steps every n iterations (n is defined by saveStep).
%           saveStep : a number (default: 5). For debug mode, determine the interval to save intermediate results.
%             psfGen : true|false (default: true). Generate a clean version of PSF by removing background or other isolated bright regions.
%             GPUJob : true|false (default: false). Use GPU for deconvolution
%        deconRotate : true|false (default: false). Run rotation after deconvolution.
%          batchSize : 1x3 vector (default: [1024, 1024, 1024]). Batch size per deconvolution task.
%          blockSize : 1x3 vector (default: [256, 256, 256]). Block/chunk size for zarr output.
%          largeFile : true|false (default: false). Use large scale deconvolution strategy.
%        largeMethod : 'inmemory'|'inplace' (default: 'inmemory'). Method for large scale deconvolution. 'inmemory': read data in memory and deconvolve region by region. 'inplace': only read the region for deconvolution in a batch, and write the result to disk, only support Zarr.
%           zarrFile : true|false (default: false). Use Zarr file as input.
%           saveZarr : true|false (default: false). Save results as Zarr files.
%         dampFactor : a number (default: 1). Damp factor for decon result by capping the intensity by the product of the raw data and the damp factor.
%        scaleFactor : empty or 1x#channel (default: []). Scale factor apply to deconvolution result. If it is a 1x#channel vector, apply a specific parameter per channel mapping to channelPatterns.
%        deconOffset : a number (default: 0). Offset added to the decon result.
%      maskFullpaths : empty or a cell array (default: {}). 2D masks to filter regions for decon in 'inplace' large decon method, in xy, xz, and yz order.
%       parseCluster : true|false (default: true). Use slurm cluster for the processing.
%        parseParfor : true|false (default: false). Use matlab parfor for paralle processing.
%      masterCompute : true|false (default: true). Master job node is involved in the processing.
%          jobLogDir : char (default: '../job_logs'). Path for the slurm job logs.
%        cpusPerTask : a number (default: 1). The number of cpu cores per task for slurm job submission. GPU jobs (with GPUConfigFile) use at most 4, so that every GPU of a node can work; the GPU config's minCPUNum can raise it.
%       taskBatchNum : empty or a number (default: []). The maximum number of files per slurm job. Files with the same settings and size share a job, which saves the start of MATLAB or the MATLAB Runtime and the GPU per file. If empty, GPU jobs take up to 8 files as long as there are at least 16 jobs, and other jobs one file.
%               uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%       unitWaitTime : a number (default: 1). The wait time per file in minutes to check whether the computing is done.
%        maxTrialNum : a number (default: 3). The max number of retries for a task.
%            mccMode : true|false (default: false). Use mcc mode.
%         configFile : empty or char (default: ''). Path for the config file for job submission in CPU.
%      GPUConfigFile : empty or char (default: ''). Path for the config file for job submission for GPU deconvolution (The conductor job can be in a CPU node).
%
%
% Author: Xiongtao Ruan
% Date: (08/26/2020)


ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x)); % data structure from loadConditionData
ip.addParameter('resultDirName', 'matlab_decon',  @(x) ischar(x));
ip.addParameter('overwrite', false,  @(x) islogical(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('skewAngle', 32.45, @isscalar);
ip.addParameter('dz', 0.5, @isscalar);
ip.addParameter('xyPixelSize', 0.108, @(x) isvector(x) && numel(x) <= 2);
ip.addParameter('save16bit', true, @(x) numel(x) == 1 && islogical(x));
ip.addParameter('parseSettingFile', false, @islogical); % use setting file to decide whether filp Z stack or not, it is  poirier over flipZstack
ip.addParameter('flipZstack', false, @islogical);
% decon parameters
ip.addParameter('background', [], @isnumeric);
ip.addParameter('dzPSF', 0.1, @isnumeric);
ip.addParameter('edgeErosion', 0, @isnumeric);
ip.addParameter('erodeByFTP', true, @islogical); % Edge erosion by the first time point (ranked the first in the inital file list for each dataset).
ip.addParameter('psfFullpaths', {'','',''}, @iscell);
ip.addParameter('deconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @isnumeric); 
ip.addParameter('OTFCumThresh', 0.9, @isnumeric); % OTF cumutative sum threshold
ip.addParameter('hannWinBounds', [0.8, 1.0], @isnumeric); % apodization range for distance matrix
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x)); % decon in skewed space
ip.addParameter('debug', false, @islogical); % debug mode for simplified code
ip.addParameter('saveStep', 5, @isnumeric); % save intermediate results every given iterations
ip.addParameter('psfGen', true, @islogical); % psf generation
ip.addParameter('GPUJob', false, @islogical); % use gpu for chuck deconvolution. 
% post decon
ip.addParameter('deconRotate', false, @islogical);
% job related parameters
ip.addParameter('batchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256], @isnumeric); % block size 
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('largeMethod', 'inmemory', @ischar); % inmemory, inplace. 
ip.addParameter('zarrFile', false, @islogical); % use zarr file as input
ip.addParameter('saveZarr', false, @islogical); % save as zarr
ip.addParameter('dampFactor', 1, @isnumeric); % damp factor for decon result
ip.addParameter('scaleFactor', [], @isnumeric); % scale factor for decon result
ip.addParameter('deconOffset', 0, @isnumeric); % offset for decon result
ip.addParameter('maskFullpaths', {}, @iscell); % 2d masks to filter regions to decon, in xy, xz, yz order
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @isnumeric);
ip.addParameter('taskBatchNum', [], @isnumeric); % files per job, [] for automatic
ip.addParameter('uuid', '', @ischar);
ip.addParameter('unitWaitTime', 1, @isnumeric);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);
ip.addParameter('GPUConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

% make sure the function is in the root PetaKit5D. 
mpath = fileparts(which(mfilename));
repo_rt = [mpath, '/../../'];
cd(repo_rt);
if ~exist([repo_rt, 'setup.m'], 'file')
    repo_rt = [mpath, '/../../'];
    cd(repo_rt);
end

pr = ip.Results;
overwrite = pr.overwrite;
resultDirName = pr.resultDirName;
% Resolution = pr.Resolution;
skewAngle = pr.skewAngle;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
channelPatterns = pr.channelPatterns;
save16bit = pr.save16bit;
parseSettingFile = pr.parseSettingFile;
flipZstack = pr.flipZstack;
% decon parameters
edgeErosion = pr.edgeErosion;
erodeByFTP = pr.erodeByFTP;
background = pr.background;
dzPSF = pr.dzPSF;
psfFullpaths = pr.psfFullpaths;
deconIter = pr.deconIter;
RLMethod = pr.RLMethod;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
hannWinBounds = pr.hannWinBounds;
skewed = pr.skewed;
GPUJob = pr.GPUJob;
% simplified version related options
debug = pr.debug;
saveStep = pr.saveStep;
psfGen = pr.psfGen;
% post decon
deconRotate = pr.deconRotate;
% job related
batchSize = pr.batchSize;
blockSize = pr.blockSize;
largeFile = pr.largeFile;
largeMethod = pr.largeMethod;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
dampFactor = pr.dampFactor;
scaleFactor = pr.scaleFactor;
deconOffset = pr.deconOffset;
maskFullpaths = pr.maskFullpaths;

jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
taskBatchNum = pr.taskBatchNum;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
masterCompute = pr.masterCompute;
uuid = pr.uuid;
unitWaitTime = pr.unitWaitTime;
maxTrialNum = pr.maxTrialNum;
mccMode = pr.mccMode;
configFile = pr.configFile;
GPUConfigFile = pr.GPUConfigFile;

if isempty(uuid)
    uuid = get_uuid();
end

% suppress directory exists warning
warning('off', 'MATLAB:MKDIR:DirectoryExists');

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
for d = 1 : nd
    dataPath = dataPaths{d};
    if ~strcmp(dataPath(end), filesep)
        dataPaths{d} = [dataPath, filesep];
    end
end

% check if decon iter is dataset specific
deconIter_mat = zeros(nd, 1);
if numel(deconIter) == nd
    deconIter_mat = deconIter;
else
    deconIter_mat = deconIter * ones(size(deconIter_mat));
end

if isscalar(overwrite)
    overwrite = repmat(overwrite, 1, 2);
end

if numel(channelPatterns) > 1 && isscalar(wienerAlpha)
    wienerAlpha = wienerAlpha * ones(1, numel(channelPatterns));
end

if numel(channelPatterns) > 1 && isscalar(OTFCumThresh)
    OTFCumThresh = OTFCumThresh * ones(1, numel(channelPatterns));
end

if isempty(background)
    background = 100;
end
if numel(channelPatterns) > 1 && isscalar(background)
    background = background * ones(1, numel(channelPatterns));
end

if isempty(scaleFactor)
    scaleFactor = 1;
end
if numel(channelPatterns) > 1 && isscalar(scaleFactor)
    scaleFactor = scaleFactor * ones(1, numel(channelPatterns));
end

if isempty(uuid)
    uuid = get_uuid();
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPath, jobLogDir);
end
    
% for deconvolution, check whether there is a gpu in the node. 
deconName = resultDirName;
if isempty(deconName)
    deconName = 'matlab_decon';
end

deconPaths = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};

    deconPath = [dataPath, deconName, filesep];
    if overwrite(1) && exist(deconPath, 'dir')
        rmdir(deconPath, 's');
    end
    if ~exist(deconPath, 'dir')
        mkdir(deconPath);
        fileattrib(deconPath, '+w', 'g');
    end
    deconPaths{d} = deconPath;

    % save decon parameters
    save('-v7.3', [deconPath, '/parameters.mat'], 'pr');
    writeJsonFile(pr, [deconPath, '/parameters.json']);
end

for f = 1 : numel(psfFullpaths)
    psfFn = psfFullpaths{f};
    if ~exist(psfFullpaths{f}, 'file')
        error('PSF file %s does not exist!', psfFn);
    end
    if psfGen 
        fprintf('PSF generation for %s ...\n', psfFn);
        [~, psfFsn] = fileparts(psfFn);

        medFactor = 1.5;
        PSFGenMethod = 'masked';
        psf = single(readtiff(psfFn));
        psf = psf_gen_new(psf, dzPSF, dz, medFactor, PSFGenMethod);
        
        if ~strcmp(RLMethod, 'omw')
            % crop psf to the bounding box (-/+ 1 pixel) and make sure the
            % center doesn't shift
            py = find(squeeze(sum(psf, [2, 3])));
            px = find(squeeze(sum(psf, [1, 3])));
            pz = find(squeeze(sum(psf, [1, 2])));
            cropSz = [min(py(1) - 1, size(psf, 1) - py(end)), min(px(1) - 1, size(psf, 2) - px(end)), min(pz(1) - 1, size(psf, 3) - pz(end))] - 1;
            cropSz = max(0, cropSz);
            bbox = [cropSz + 1, size(psf, 1 : 3) - cropSz];
            psf = psf(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
        end

        for d = 1 : nd 
            dataPath = dataPaths{d};
            psfgen_filename = sprintf('%s/%s/psfgen/%s_%s.tif', dataPath, deconName, psfFsn, RLMethod);
            tmp_filename = sprintf('%s/%s/psfgen/%s_%s.tif', dataPath, deconName, psfFsn, uuid);
            writetiff(psf, tmp_filename);
            movefile(tmp_filename, psfgen_filename)
        end
        fprintf('Done!\n\n');
    end
end

if largeFile && strcmp(largeMethod, 'inplace')
    saveZarr = true;
end

%% check existing files and parse channels
Streaming = false;
Decon = true;
minModifyTime = 1;
[fnames, fdinds, gfnames, partialvols, dataSizes, flipZstack_mat, latest_modify_times, FTP_inds, erodeMaskFullpaths] = ...
    XR_parseImageFilenames(dataPaths, channelPatterns, parseSettingFile, flipZstack, Decon, deconPaths, Streaming, minModifyTime, zarrFile);

nF = numel(fnames);

% for jobs on the cluster, make the OMW back projector of each PSF here once,
% instead of in every first job of a PSF; it also changes in the last bits
% with the number of threads, so this way it doesn't depend on the job sizes
if parseCluster && strcmp(RLMethod, 'omw') && ~largeFile
    % the parameters as the jobs get them from the function string
    hannWinBounds_j = str2num(mat2str_comma(hannWinBounds));
    for p = 1 : min(numel(psfFullpaths), numel(channelPatterns))
        psfFn = psfFullpaths{p};
        [~, psfFsn] = fileparts(psfFn);
        wienerAlpha_j = str2double(sprintf('%.20f', wienerAlpha(p)));
        OTFCumThresh_j = str2double(sprintf('%.20f', OTFCumThresh(p)));
        bpFsns = {sprintf('%s_back_projector_alpha_%0.6f_otf_cum_thresh_%0.6f.tif', psfFsn, wienerAlpha_j, OTFCumThresh_j), ...
            sprintf('%s_OTF.tif', psfFsn), sprintf('%s_OTF_mask_otf_cum_thresh_%0.6f.tif', psfFsn, OTFCumThresh_j)};
        srcPath = '';
        for d = 1 : nd
            % the first file of the channel in the dataset
            fullnames_d = strcat(dataPaths{d}, fnames(fdinds == d));
            ind = find(~cellfun(@isempty, regexpi(fullnames_d, channelPatterns{p})), 1);
            if isempty(ind)
                continue;
            end
            psfgenPath = sprintf('%s/%s/psfgen/', dataPaths{d}, deconName);
            if psfGen
                psf = single(readtiff(sprintf('%s/%s_%s.tif', psfgenPath, psfFsn, RLMethod)));
            else
                psf = single(readtiff(psfFn));
            end
            % RLdecon crops a PSF larger than the data first
            if any(size(psf, 1 : 3) > getImageSize(fullnames_d{ind}))
                continue;
            end
            if ~exist(psfgenPath, 'dir')
                mkdir(psfgenPath);
            end
            if isempty(srcPath)
                % without the figure, as in the cluster jobs (no Java)
                omw_backprojector_files(psf, psfgenPath, psfFsn, wienerAlpha_j, OTFCumThresh_j, skewed, hannWinBounds_j, psfFn, false);
                srcPath = psfgenPath;
            elseif ~exist([psfgenPath, bpFsns{1}], 'file')
                % same PSF and parameters as in the first dataset
                for i = 1 : numel(bpFsns)
                    if exist([srcPath, bpFsns{i}], 'file')
                        copyfile([srcPath, bpFsns{i}], [psfgenPath, bpFsns{i}]);
                    end
                end
            end
        end
    end
end

% flags: for thee: deskew w/o rotate, decon w/o rotate, rotate
is_done_flag = false(nF, 1);
trial_counter = zeros(nF, 1);

% For no-streaming computing, first skip the disabled options. 
if ~Decon
    is_done_flag(:, 1) = true;
end
if parseCluster
    job_ids = -ones(nF, 1);
    % array task id of each file's job (a job with several files has the
    % task id of its first file)
    job_task_ids = zeros(nF, 1);
end
% files whose job is pending or running: after the last try, the wrapper
% waits for its job before it returns
job_active = false(nF, 1);

% files per cluster job: a GPU deconvolves a tile in about as long as a job
% takes to start MATLAB or the MATLAB Runtime and the GPU, so by default GPU
% jobs take several files, while keeping enough jobs for the GPUs of several
% nodes (decided in each round from the files to submit)
autoBatch = isempty(taskBatchNum);
if ~autoBatch
    taskBatchNum = max(1, round(taskBatchNum));
end
batchFiles = parseCluster && ~largeFile && ((autoBatch && GPUJob) || (~autoBatch && taskBatchNum > 1));

nF_done = 0;
ts = tic;
% use while loop to perform computing for all images
while ~all(is_done_flag | (trial_counter >= maxTrialNum & ~job_active), 'all')
    % files to submit in jobs of several files in this round
    sub_flag = false(nF, 1);
    sub_func_strs = cell(nF, 1);
    sub_frames = cell(nF, 1);
    sub_sizes = cell(nF, 1);
    sub_outNames = cell(nF, 1);
    sub_mems = zeros(nF, 1);
    sub_cpus = zeros(nF, 1);
    sub_configFiles = cell(nF, 1);
    % results finished by cluster jobs: on NFS, exist sees a new file up to
    % ~30 s late, a listing of the folder at once (a finished job whose
    % result is not seen yet would be submitted again)
    if parseCluster
        deconListings = list_folders(deconPaths);
    end
    for f = 1 : nF
        if all(is_done_flag(f, :))
            continue;
        end

        % first deskew and rotation
        fname = fnames{f};
        [~, fsname] = fileparts(fname);
        fdind = fdinds(f);
        dataPath = dataPaths{fdind};
                
        frameFullpath = [dataPath, fname];
        % check wheter the file is deleted during the computing.
        if ~exist(frameFullpath, 'file') || (zarrFile && ~exist(frameFullpath, 'dir'))
            is_done_flag(f, :) = true;
            continue
        end
        
        task_id = rem(f, 5000);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % deconvolution
        % input chosen order is dsr, ds, raw (in decreased priority order)
        dcframeFullpath = frameFullpath;

        dc_dz = dz;
        dc_dzPSF = dzPSF;
        dc_psfFullpaths = psfFullpaths;
        
        if ~exist(dcframeFullpath, 'file')
            continue;
        end
                        
        deconPath = deconPaths{fdind};
        if saveZarr
            deconFullpath = sprintf('%s/%s.zarr', deconPath, fsname);
        else
            deconFullpath = sprintf('%s/%s.tif', deconPath, fsname);
        end
        dctmpFullpath = sprintf('%s.tmp', deconFullpath(1 : end - 4));

        [~, deconFsn, deconExt] = fileparts(deconFullpath);
        if exist(deconFullpath, 'file') || (saveZarr && exist(deconFullpath, 'dir')) ...
                || (parseCluster && any(strcmp(deconListings{fdind}, [deconFsn, deconExt])))
            is_done_flag(f, 1) = true;
            if exist(dctmpFullpath, 'file')
                delete(dctmpFullpath);
            end
        end
        
        % for erodeByFTP, check if the mask file exist
        erodeMaskFullpath = '';
        SaveMaskfile = false;
        % do not apply erode by first time point for cuda decon for now
        % (04/19/2020)
        if edgeErosion > 0 && erodeByFTP && ~largeFile
            FTP_ind = FTP_inds(fdind);
            if f == FTP_ind
                SaveMaskfile = true;
                % if decon result exist, but mask file not exist, rerun
                % it to save the mask.
                if is_done_flag(f, 1) && ~exist(erodeMaskFullpaths{fdind}, 'dir') ...
                        && ~(parseCluster && in_listing(erodeMaskFullpaths{fdind}))
                    is_done_flag(f, 1) = false;
                    delete(deconFullpath);
                end
            else
                % only check for the ones not finished. 
                if ~is_done_flag(f, 1)
                    erodeMaskFullpath = erodeMaskFullpaths{fdind};
                    if ~exist(erodeMaskFullpath, 'dir')
                        continue;
                    end

                    mask_sz = getImageSize(erodeMaskFullpath);
                    img_sz = getImageSize(dcframeFullpath);
                    if any(mask_sz ~= img_sz)
                        warning(['The image size [%s] does not match the defined ', ... 
                            'mask size [%s], use its own mask for edge erosion...'], ...
                            num2str(img_sz, '%d '), num2str(mask_sz, '%d '));
                        erodeMaskFullpath = '';
                    end
                end
            end
        end

        if ~is_done_flag(f, 1)     
            psfMapping =  ~cellfun(@isempty, regexpi(frameFullpath, channelPatterns));
            psfFullpath = dc_psfFullpaths{psfMapping};
            background_f = background(psfMapping);
            deconIter_f = deconIter_mat(fdind);
            if strcmp(RLMethod, 'omw')
                wienerAlpha_f = wienerAlpha(psfMapping);
                OTFCumThresh_f = OTFCumThresh(psfMapping);
            else
                wienerAlpha_f = wienerAlpha(1);
                OTFCumThresh_f = OTFCumThresh(1);
            end
            scaleFactor_f = scaleFactor(psfMapping);
            flipZstack = flipZstack_mat(f);
            maskFullpaths_str = sprintf('{''%s''}', strjoin(maskFullpaths, ''','''));

            func_str = sprintf(['XR_RLdeconFrame3D(''%s'',%.10f,%.10f,''%s'',''PSFfile'',''%s'',', ...
                '''dzPSF'',%.10f,''background'',[%d],''skewAngle'',%d,''flipZstack'',%s,', ...
                '''edgeErosion'',%d,''ErodeMaskfile'',''%s'',''SaveMaskfile'',%s,''Rotate'',%s,', ...
                '''deconIter'',%d,''RLMethod'',''%s'',''wienerAlpha'',%.20f,''OTFCumThresh'',%.20f,', ...
                '''hannWinBounds'',%s,''skewed'',[%s],''debug'',%s,''saveStep'',%d,''psfGen'',%s,''saveZarr'',%s,', ...
                '''parseCluster'',%s,''parseParfor'',%s,''GPUJob'',%s,''save16bit'',%s,', ...
                '''largeFile'',%s,''largeMethod'',''%s'',''batchSize'',%s,''blockSize'',%s,', ...
                '''dampFactor'',%d,''scaleFactor'',[%d],''deconOffset'',%d,''maskFullpaths'',%s,', ...
                '''uuid'',''%s'',''cpusPerTask'',%d,''mccMode'',%s,''configFile'',''%s'',''GPUConfigFile'',''%s'')'], ...
                dcframeFullpath, xyPixelSize, dc_dz, deconPath, psfFullpath, ...
                dc_dzPSF, background_f, skewAngle, string(flipZstack), edgeErosion, ...
                erodeMaskFullpath, string(SaveMaskfile), string(deconRotate), deconIter_f, ...
                RLMethod, wienerAlpha_f, OTFCumThresh_f, mat2str_comma(hannWinBounds), ...
                string(skewed), string(debug), saveStep, string(psfGen), string(saveZarr), ...
                string(parseCluster), string(parseParfor), string(GPUJob), string(save16bit), ...
                string(largeFile), largeMethod, strrep(mat2str(batchSize), ' ', ','), ...
                strrep(mat2str(blockSize), ' ', ','), dampFactor, scaleFactor_f, ...
                deconOffset, maskFullpaths_str, uuid, cpusPerTask, string(mccMode), ...
                configFile, GPUConfigFile);

            if exist(dctmpFullpath, 'file') || parseCluster
                if parseCluster
                    dsz = getImageSize(dcframeFullpath);
                    estMem = prod(dsz) * 4 / 2^30;
                    cur_configFile = configFile;
                    cpusPerTask_f = cpusPerTask;
                    if ~GPUJob
                        memAllocate = estMem * 10;
                        if largeFile && strcmp(largeMethod, 'inplace')
                            memAllocate = prod(batchSize * 2) * 4 / 2^30 * 10;
                        end
                    else
                        memAllocate = estMem * 10;
                        if largeFile && strcmp(largeMethod, 'inplace')
                            memAllocate = prod(batchSize) * 4 / 2^30 * 20;
                        else
                            cur_configFile = GPUConfigFile;
                            % the GPU deconvolves and the CPUs only read and write
                            % the file, so a GPU job takes at most 4 CPUs, which
                            % lets every GPU of a 16-core 4-GPU node work (the
                            % CPU-oriented cpusPerTask took a whole node for one GPU)
                            if ~largeFile && ~isempty(GPUConfigFile)
                                cpusPerTask_f = min(cpusPerTask, 4);
                            end
                        end
                    end

                    if batchFiles
                        % submitted after this round with other files of the
                        % same settings, unless its job is still pending or
                        % running or it was submitted maxTrialNum times
                        sub_flag(f) = true;
                        sub_func_strs{f} = func_str;
                        sub_frames{f} = dcframeFullpath;
                        sub_sizes{f} = mat2str(dsz);
                        sub_outNames{f} = [deconFsn, deconExt];
                        sub_mems(f) = memAllocate;
                        sub_cpus(f) = cpusPerTask_f;
                        sub_configFiles{f} = cur_configFile;
                    else
                        job_id = job_ids(f, 1);
                        if trial_counter(f, 1) < maxTrialNum
                            [job_id, job_status, submit_status] = generic_single_job_submit_wrapper(func_str, ...
                                job_id, task_id, jobLogFname=job_log_fname, jobErrorFname=job_log_error_fname, ...
                                lastFile=false, cpusPerTask=cpusPerTask_f, memAllocate=memAllocate, ...
                                mccMode=mccMode, configFile=cur_configFile);
                        else
                            % the last try: only check whether its job is still there
                            job_status = check_slurm_job_status(job_id, rem(task_id, 5000));
                            submit_status = 0;
                        end

                        job_ids(f, 1) = job_id;
                        job_task_ids(f) = task_id;
                        trial_counter(f, 1) = trial_counter(f, 1) + submit_status;
                        job_active(f) = submit_status > 0 || job_status >= 0;
                    end
                else
                    temp_file_info = dir(dctmpFullpath);
                    if minutes(datetime('now') - [temp_file_info.date]) < unitWaitTime
                        continue; 
                    else
                        fclose(fopen(dctmpFullpath, 'w'));
                    end
                end
            else
                fclose(fopen(dctmpFullpath, 'w'));
            end

            if ~parseCluster
                % fileInfo = imfinfo(dsFullpath);
                tic; feval(str2func(['@()', func_str])); toc;
                fprintf('\n');
                trial_counter(f, 1) = trial_counter(f, 1) + 1;
            end

            % check if computing is done
            if exist(deconFullpath, 'file') || (saveZarr && exist(deconFullpath, 'dir'))
                is_done_flag(f, 1) = true;
                if exist(dctmpFullpath, 'file')
                    delete(dctmpFullpath);
                end
            end
        end
    end

    % submit the collected files in their order, up to batchNum files with
    % the same settings per job; files whose job is still pending or running
    % wait for it
    if parseCluster && any(sub_flag)
        sub_inds = find(sub_flag);
        [jobs, ~, job_inds] = unique([job_ids(sub_inds), job_task_ids(sub_inds)], 'rows');
        job_live = false(size(jobs, 1), 1);
        for j = 1 : size(jobs, 1)
            if jobs(j, 1) > 0
                job_live(j) = check_slurm_job_status(jobs(j, 1), jobs(j, 2)) >= 0;
            end
        end
        % a file waits for its job while it is pending or running
        job_active(sub_inds) = job_live(job_inds);
        sub_inds = sub_inds(~job_live(job_inds));
        % a job may have finished since the listing at the start of the round
        rechecked = sub_inds(trial_counter(sub_inds) > 0);
        if ~isempty(rechecked)
            deconListings = list_folders(deconPaths);
            out_done = arrayfun(@(x) any(strcmp(deconListings{fdinds(x)}, sub_outNames{x})), rechecked);
            is_done_flag(rechecked(out_done), 1) = true;
            sub_inds = setdiff(sub_inds, rechecked(out_done), 'stable');
        end
        % files whose jobs ended maxTrialNum times are not submitted again
        sub_inds = sub_inds(trial_counter(sub_inds) < maxTrialNum);
        % files submitted before (their job failed) get a job each, so that
        % a file that fails can't use up the trials of others
        retry = trial_counter(sub_inds) > 0;
        % by default up to 8 files per job, as long as there are 16 jobs
        batchNum = taskBatchNum;
        if autoBatch
            batchNum = max(1, min(8, floor(sum(~retry) / 16)));
        end
        % a job takes files of the same settings (the function string after
        % the file path) and size, so the GPU keeps the OTFs of the first one
        param_strs = cellfun(@(x, y) x(numel(sprintf('XR_RLdeconFrame3D(''%s''', y)) + 1 : end), ...
            sub_func_strs(sub_inds), sub_frames(sub_inds), 'unif', 0);
        [~, ~, param_inds] = unique(strcat(param_strs, '|', sub_sizes(sub_inds)), 'stable');
        param_inds(retry) = max([param_inds; 0]) + (1 : nnz(retry))';
        for p = 1 : max([param_inds; 0])
            p_pos = find(param_inds == p);
            for b = 1 : ceil(numel(p_pos) / batchNum)
                b_pos = p_pos((b - 1) * batchNum + 1 : min(b * batchNum, numel(p_pos)));
                bfs = sub_inds(b_pos);
                func_str = sub_func_strs{bfs(1)};
                if numel(bfs) > 1
                    frame_strs = cellfun(@(x) sprintf('''%s''', x), sub_frames(bfs), 'unif', 0);
                    func_str = sprintf('XR_RLdeconFrame3D({%s}%s', strjoin(frame_strs(:)', ','), param_strs{b_pos(1)});
                end
                task_id = rem(bfs(1), 5000);
                [job_id, ~, submit_status] = generic_single_job_submit_wrapper(func_str, ...
                    -1, task_id, jobLogFname=job_log_fname, jobErrorFname=job_log_error_fname, ...
                    lastFile=false, cpusPerTask=sub_cpus(bfs(1)), memAllocate=max(sub_mems(bfs)), ...
                    mccMode=mccMode, configFile=sub_configFiles{bfs(1)});

                job_ids(bfs) = job_id;
                job_task_ids(bfs) = task_id;
                trial_counter(bfs) = trial_counter(bfs) + submit_status;
                job_active(bfs) = submit_status > 0;
            end
        end
    end

    % wait for running jobs finishing and checking for new coming images
    if ~all(is_done_flag | (trial_counter >= maxTrialNum & ~job_active), 'all')
        % with jobs on the cluster, check again soon when few files are left;
        % each check runs squeue per job, so wait longer for many files
        waitTime = 30;
        if parseCluster
            waitTime = min(30, max(3, 0.3 * sum(~all(is_done_flag, 2))));
        end
        pause(waitTime);
    end
    if nF_done < sum(is_done_flag)
        nF_done = sum(is_done_flag);
        fprintf('Time %0.2f s: %d / %d (%0.2f%%) are finished!\n', toc(ts), nF_done, nF, nF_done / nF * 100);
    end
end

if all(is_done_flag)
    fprintf('All output files (%d / %d) are finished!\n\n', nF, nF);
end

end


function [listings] = list_folders(folders)
% the names in each folder (results are moved into place under their final
% names, so a listed name is a finished result)

listings = cell(numel(folders), 1);
for i = 1 : numel(folders)
    dir_info = dir(folders{i});
    listings{i} = {dir_info.name};
end

end


function [tf] = in_listing(fullpath)
% true if the listing of the folder of fullpath has its name

[pstr, fsn, ext] = fileparts(fullpath);
listing = list_folders({pstr});
tf = any(strcmp(listing{1}, [fsn, ext]));

end

