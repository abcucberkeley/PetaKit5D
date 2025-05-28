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
%         background : empty or a number (default: []). Background subtraction before deconvolution.
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
%        cpusPerTask : a number (default: 1). The number of cpu cores per task for slurm job submission.
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

if numel(overwrite) == 1
    overwrite = repmat(overwrite, 1, 2);
end

if numel(channelPatterns) > 1 && numel(wienerAlpha) == 1
    wienerAlpha = wienerAlpha * ones(1, numel(channelPatterns));
end

if numel(channelPatterns) > 1 && numel(OTFCumThresh) == 1
    OTFCumThresh = OTFCumThresh * ones(1, numel(channelPatterns));
end

if isempty(scaleFactor)
    scaleFactor = 1;
end
if numel(channelPatterns) > 1 && numel(scaleFactor) == 1
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

% flags: for thee: deskew w/o rotate, decon w/o rotate, rotate
is_done_flag = false(nF, 1);
trial_counter = zeros(nF, 1);

% For no-streaming computing, first skip the disabled options. 
if ~Decon
    is_done_flag(:, 1) = true;
end
if parseCluster
    job_ids = -ones(nF, 1);
end

nF_done = 0;
ts = tic;
% use while loop to perform computing for all images
while ~all(is_done_flag | trial_counter >= maxTrialNum, 'all')
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

        if exist(deconFullpath, 'file') || (saveZarr && exist(deconFullpath, 'dir'))
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
        if edgeErosion > 0 && erodeByFTP && ~(largeFile && strcmp(largeMethod, 'inplace'))
            FTP_ind = FTP_inds(fdind);
            if f == FTP_ind
                SaveMaskfile = true;
                % if decon result exist, but mask file not exist, rerun
                % it to save the mask. 
                if is_done_flag(f, 1) && ~exist(erodeMaskFullpaths{fdind}, 'dir')
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
            deconIter_f = deconIter_mat(fdind);
            wienerAlpha_f = wienerAlpha(psfMapping);
            OTFCumThresh_f = OTFCumThresh(psfMapping);
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
                dc_dzPSF, background, skewAngle, string(flipZstack), edgeErosion, ...
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
                        end
                    end

                    job_id = job_ids(f, 1);
                    [job_id, ~, submit_status] = generic_single_job_submit_wrapper(func_str, ...
                        job_id, task_id, jobLogFname=job_log_fname, jobErrorFname=job_log_error_fname, ...
                        lastFile=false, cpusPerTask=cpusPerTask, memAllocate=memAllocate, ...
                        mccMode=mccMode, configFile=cur_configFile);

                    job_ids(f, 1) = job_id;
                    trial_counter(f, 1) = trial_counter(f, 1) + submit_status;
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

    % wait for running jobs finishing and checking for new coming images
    if ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') 
        pause(30);
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

