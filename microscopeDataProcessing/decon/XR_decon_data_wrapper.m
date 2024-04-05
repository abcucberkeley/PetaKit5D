function [] = XR_decon_data_wrapper(dataPaths, varargin)
% data-level deconvolution wrapper, support cuda decon, cpp decon and
% matlab decon. Adapted from the microscope pipeline. 
% Support of option to rotate PSF, also user-defined resolution.
% 
% Author: Xiongtao Ruan
% Date: (08/26/2020)
% 
% xruan (10/26/2020): add support for zarr
% xruan (03/24/2021): add support for normalization of intensity by the
% first time point when saving to 16bit. 
% xruan (03/25/2021): add options for different versions of rl method
% xruan (06/10/2021): add support for threshold and debug mode in simplified version. 
% xruan (06/11/2021): add support for gpu computing for chuck decon in matlab decon wrapper
% xruan (07/13/2021): add support for the processing of flipped files
% (currently only add support for matlab decon)
% xruan (10/19/2021): add support for dataset specific iteration


ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x)); % data structure from loadConditionData
ip.addParameter('deconPathstr', '',  @(x) ischar(x));
ip.addParameter('Overwrite', false,  @(x) (numel(x) == 1 || numel(x) == 2) && islogical(x));
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('Channels', [488, 560, 642], @isnumeric);
ip.addParameter('SkewAngle', 32.45, @isscalar);
ip.addParameter('dz', 0.5, @isscalar);
ip.addParameter('xyPixelSize', 0.108, @isscalar);
ip.addParameter('Reverse', true, @islogical);
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('sCMOSCameraFlip', false, @islogical);
ip.addParameter('Save16bit', false, @(x) numel(x) == 1 && islogical(x));
ip.addParameter('onlyFirstTP', false, @islogical);
ip.addParameter('parseSettingFile', false, @islogical); % use setting file to decide whether filp Z stack or not, it is  poirier over flipZstack
ip.addParameter('flipZstack', false, @islogical);
% pipeline steps
ip.addParameter('Decon', true, @islogical);
% decon parameters
ip.addParameter('cudaDecon', false, @islogical);
ip.addParameter('cppDecon', false, @islogical);
ip.addParameter('cppDeconPath', '/global/home/groups/software/sl-7.x86_64/modules/RLDecon_CPU/20200718/build-cluster/cpuDeconv', @ischar);
ip.addParameter('loadModules', 'module load gcc/4.8.5; module load fftw/3.3.6-gcc; module load boost/1.65.1-gcc; module load libtiff/4.1.0; ', @ischar);
ip.addParameter('cudaDeconPath', '/global/home/groups/software/sl-7.x86_64/modules/cudaDecon/bin/cudaDeconv' , @ischar);
ip.addParameter('OTFGENPath', '/global/home/groups/software/sl-7.x86_64/modules/cudaDecon/bin/radialft' , @ischar); % point to radialft file
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('dzPSF', 0.1, @isnumeric);
ip.addParameter('EdgeErosion', 8, @isnumeric);
ip.addParameter('ErodeByFTP', true, @islogical); % Edge erosion by the first time point (ranked the first in the inital file list for each dataset).
ip.addParameter('deconRotate', false, @islogical);
ip.addParameter('psfFullpaths', {'','',''}, @iscell);
ip.addParameter('rotatePSF', false, @islogical);
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @isnumeric); 
ip.addParameter('OTFCumThresh', 0.9, @isnumeric); % OTF cumutative sum threshold
ip.addParameter('hanWinBounds', [0.8, 1.0], @isnumeric); % apodization range for distance matrix
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x)); % decon in skewed space
ip.addParameter('fixIter', false, @islogical); 
ip.addParameter('errThresh', [], @isnumeric); % error threshold for simplified code
ip.addParameter('debug', false, @islogical); % debug mode for simplified code
ip.addParameter('saveStep', 5, @isnumeric); % save intermediate results every given iterations
ip.addParameter('psfGen', true, @islogical); % psf generation
ip.addParameter('GPUJob', false, @islogical); % use gpu for chuck deconvolution. 
% job related parameters
ip.addParameter('BatchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256], @isnumeric); % block size 
ip.addParameter('zarrSubSize', [], @isnumeric);
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('largeMethod', 'inmemory', @ischar); % inmemory, inplace. 
ip.addParameter('zarrFile', false, @islogical); % use zarr file as input
ip.addParameter('saveZarr', false, @islogical); % save as zarr
ip.addParameter('damper', 1, @isnumeric); % damp factor for decon result
ip.addParameter('scaleFactor', [], @isnumeric); % scale factor for decon result
ip.addParameter('deconOffset', 0, @isnumeric); % offset for decon result
ip.addParameter('deconMaskFns', {}, @iscell); % 2d masks to filter regions to decon, in xy, xz, yz order
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 10, @isnumeric);
ip.addParameter('maxWaitLoopNum', 10, @isnumeric); % the max number of loops the loop waits with all existing files processed. 
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('ConfigFile', '', @ischar);
ip.addParameter('GPUConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

% make sure the function is in the root of XR_Repository or LLSM5DTools. 
mpath = fileparts(which(mfilename));
repo_rt = [mpath, '/../../../'];
cd(repo_rt);
if ~exist([repo_rt, 'setup.m'], 'file')
    repo_rt = [mpath, '/../../'];
    cd(repo_rt);
end

pr = ip.Results;
Overwrite = pr.Overwrite;
deconPathstr = pr.deconPathstr;
% Resolution = pr.Resolution;
SkewAngle = pr.SkewAngle;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
ObjectiveScan = pr.ObjectiveScan;
Reverse = pr.Reverse;
ChannelPatterns = pr.ChannelPatterns;
Save16bit = pr.Save16bit;
onlyFirstTP = pr.onlyFirstTP;
parseSettingFile = pr.parseSettingFile;
flipZstack = pr.flipZstack;
% decon parameters
Decon = pr.Decon;
cppDecon = pr.cppDecon;
cudaDecon = pr.cudaDecon;
cppDeconPath = pr.cppDeconPath;
loadModules = pr.loadModules;
cudaDeconPath = pr.cudaDeconPath;
OTFGENPath = pr.OTFGENPath;
EdgeErosion = pr.EdgeErosion;
ErodeByFTP = pr.ErodeByFTP;
Background = pr.Background;
dzPSF = pr.dzPSF;
psfFullpaths = pr.psfFullpaths;
rotatePSF = pr.rotatePSF;
DeconIter = pr.DeconIter;
deconRotate = pr.deconRotate;
RLMethod = pr.RLMethod;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
skewed = pr.skewed;
GPUJob = pr.GPUJob;
% simplified version related options
fixIter = pr.fixIter;
errThresh = pr.errThresh;
debug = pr.debug;
saveStep = pr.saveStep;
psfGen = pr.psfGen;

% job related
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
zarrSubSize = pr.zarrSubSize;
largeFile = pr.largeFile;
largeMethod = pr.largeMethod;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
damper = pr.damper;
scaleFactor = pr.scaleFactor;
deconOffset = pr.deconOffset;
deconMaskFns = pr.deconMaskFns;

jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;
GPUConfigFile = pr.GPUConfigFile;

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
DeconIter_mat = zeros(nd, 1);
if numel(DeconIter) == nd
    DeconIter_mat = DeconIter;
else
    DeconIter_mat = DeconIter * ones(size(DeconIter_mat));
end

if numel(Overwrite) == 1
    Overwrite = repmat(Overwrite, 1, 2);
end

if numel(ChannelPatterns) > 1 && numel(wienerAlpha) == 1
    wienerAlpha = wienerAlpha * ones(1, numel(ChannelPatterns));
end

if numel(ChannelPatterns) > 1 && numel(OTFCumThresh) == 1
    OTFCumThresh = OTFCumThresh * ones(1, numel(ChannelPatterns));
end

if isempty(scaleFactor)
    scaleFactor = 1;
end
if numel(ChannelPatterns) > 1 && numel(scaleFactor) == 1
    scaleFactor = scaleFactor * ones(1, numel(ChannelPatterns));
end

if isempty(uuid)
    uuid = get_uuid();
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname] = checkSlurmCluster(dataPath, jobLogDir);
end
    
% for deconvolution, check whether there is a gpu in the node. if not, for
% cudaDecon, set parseCluster as true. 
if Decon
    % if both cudaDecon and cppDecon are true, use cppDecon
    if cudaDecon && cppDecon
        cudaDecon = false;
    end

    if cudaDecon && gpuDeviceCount() < 1 && ~parseCluster
        warning('There is no GPU in the node, and ther cluster is also not available. Set cudaDecon as false!');
        cudaDecon = false;
    end
        
    if cudaDecon
        deconName = 'GPUdecon';
    elseif cppDecon
        deconName = 'CPPdecon';
    else
        deconName = 'matlab_decon';
    end

    if ~isempty(deconPathstr)
        deconName = deconPathstr;
    end
        
    deconPaths = cell(nd, 1);
    for d = 1 : nd
        dataPath = dataPaths{d};

        deconPath = [dataPath, deconName, filesep];
        if Overwrite(1) && exist(deconPath, 'dir')
            rmdir(deconPath, 's');
        end
        if ~exist(deconPath, 'dir')
            mkdir(deconPath);
            fileattrib(deconPath, '+w', 'g');
        end
        deconPaths{d} = deconPath;
        
        % save decon parameters
        save('-v7.3', [deconPath, '/parameters.mat'], 'pr');
        writetable(struct2table(pr, 'AsArray', true), [deconPath, '/parameters.txt'])
    end
    if rotatePSF
        rotPSFFullpaths = cell(numel(psfFullpaths), 1);
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
        end
        
        %{
        if rotatePSF
            XR_rotate_PSF(psfFullpaths{f}, 'xyPixelSize', xyPixelSize, 'dz', dzPSF);
            [psfPath, fsname] = fileparts(psfFullpaths{f});
            rotPSFFullpaths{f} = [psfPath, '/Rotated/', fsname, '.tif'];
        end
        %}
    end
end

if largeFile && strcmp(largeMethod, 'inplace')
    saveZarr = true;
end

%% check existing files and parse channels
Streaming = false;
minModifyTime = 1;
[fnames, fdinds, gfnames, partialvols, dataSizes, flipZstack_mat, latest_modify_times, FTP_inds, maskFullpaths] = ...
    XR_parseImageFilenames(dataPaths, ChannelPatterns, parseSettingFile, flipZstack, Decon, deconPaths, Streaming, minModifyTime, zarrFile);

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
        if Decon
            % input chosen order is dsr, ds, raw (in decreased priority order)
            dcframeFullpath = frameFullpath;

            dc_dz = dz;
            dc_dzPSF = dzPSF;
            dc_psfFullpaths = psfFullpaths;
            if rotatePSF
                dc_dzPSF = xyPixelSize;
                dc_psfFullpaths = rotPSFFullpaths;
            end
            
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
            
            % for ErodeByFTP, check if the mask file exist
            maskFullpath = '';
            SaveMaskfile = false;
            % do not apply erode by first time point for cuda decon for now
            % (04/19/2020)
            if EdgeErosion > 0 && ErodeByFTP && ~cudaDecon && ~(largeFile && strcmp(largeMethod, 'inplace'))
                FTP_ind = FTP_inds(fdind);
                if f == FTP_ind
                    SaveMaskfile = true;
                    % if decon result exist, but mask file not exist, rerun
                    % it to save the mask. 
                    if is_done_flag(f, 1) && ~exist(maskFullpaths{fdind}, 'file')
                        is_done_flag(f, 1) = false;
                        delete(deconFullpath);
                    end
                else
                    % only check for the ones not finished. 
                    if ~is_done_flag(f, 1)
                        maskFullpath = maskFullpaths{fdind};
                        if ~exist(maskFullpath, 'file')
                            continue;
                        end

                        mask_sz = getImageSize(maskFullpath);
                        img_sz = getImageSize(dcframeFullpath);
                        if any(mask_sz ~= img_sz)
                            warning(['The image size [%s] does not match the defined ', ... 
                                'mask size [%s], use its own mask for edge erosion...'], ...
                                num2str(img_sz, '%d '), num2str(mask_sz, '%d '));
                            maskFullpath = '';
                        end
                    end
                end
            end            
        else
            is_done_flag(f, 1) = true;    
        end

        if ~is_done_flag(f, 1)     
            psfMapping =  ~cellfun(@isempty, regexpi(frameFullpath, ChannelPatterns));
            psfFullpath = dc_psfFullpaths{psfMapping};
            DeconIter_f = DeconIter_mat(fdind);
            wienerAlpha_f = wienerAlpha(psfMapping);
            OTFCumThresh_f = OTFCumThresh(psfMapping);
            scaleFactor_f = scaleFactor(psfMapping);
            flipZstack = flipZstack_mat(f);
            deconMaskFns_str = sprintf('{''%s''}', strjoin(deconMaskFns, ''','''));
            
            % do not use rotation in decon functions
            if cudaDecon
                func_str = sprintf(['XR_cudaDeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                    '''cudaDeconPath'',''%s'',''OTFGENPath'',''%s'',''dzPSF'',%.10f,''Background'',[%d],', ...
                    '''SkewAngle'',%d,''Rotate'',%s,''DeconIter'',%d,''Save16bit'',%s,''largeFile'',%s)'], ...
                    dcframeFullpath, xyPixelSize, dc_dz, psfFullpath, cudaDeconPath, OTFGENPath, dc_dzPSF, ...
                    Background, SkewAngle, string(deconRotate), DeconIter_f, string(Save16bit), string(largeFile));
            elseif cppDecon
                func_str = sprintf(['XR_cppDeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                    '''cppDeconPath'',''%s'',''loadModules'',''%s'',''dzPSF'',%.10f,''Background'',[%d],', ...
                    '''SkewAngle'',%d,''EdgeErosion'',%d,''ErodeMaskfile'',''%s'',''SaveMaskfile'',%s,', ...
                    '''Rotate'',%s,''DeconIter'',%d,''Save16bit'',%s,''largeFile'',%s)'], dcframeFullpath, ...
                    xyPixelSize, dc_dz, psfFullpath, cppDeconPath, loadModules, dc_dzPSF, Background, ...
                    SkewAngle, EdgeErosion, maskFullpath, string(SaveMaskfile), string(deconRotate), ...
                    DeconIter_f, string(Save16bit), string(largeFile));
            else
                func_str = sprintf(['XR_RLdeconFrame3D(''%s'',%.10f,%.10f,''%s'',''PSFfile'',''%s'',', ...
                    '''dzPSF'',%.10f,''Background'',[%d],''SkewAngle'',%d,''flipZstack'',%s,', ...
                    '''EdgeErosion'',%d,''ErodeMaskfile'',''%s'',''SaveMaskfile'',%s,''Rotate'',%s,', ...
                    '''DeconIter'',%d,''RLMethod'',''%s'',''wienerAlpha'',%.20f,''OTFCumThresh'',%.20f,', ...
                    '''skewed'',[%s],''fixIter'',%s,''errThresh'',[%0.20f],''debug'',%s,''saveStep'',%d,', ...
                    '''psfGen'',%s,''saveZarr'',%s,''parseCluster'',%s,''parseParfor'',%s,''GPUJob'',%s,', ...
                    '''Save16bit'',%s,''largeFile'',%s,''largeMethod'',''%s'',''BatchSize'',%s,''BlockSize'',%s,', ...
                    '''zarrSubSize'',%s,''damper'',%d,''scaleFactor'',[%d],''deconOffset'',%d,''deconMaskFns'',%s,', ...
                    '''uuid'',''%s'',''cpusPerTask'',%d,''mccMode'',%s,''ConfigFile'',''%s'',''GPUConfigFile'',''%s'')'], ...
                    dcframeFullpath, xyPixelSize, dc_dz, deconPath, psfFullpath, dc_dzPSF, Background, SkewAngle, ...
                    string(flipZstack), EdgeErosion, maskFullpath, string(SaveMaskfile), string(deconRotate), ...
                    DeconIter_f, RLMethod, wienerAlpha_f, OTFCumThresh_f, string(skewed), string(fixIter), errThresh, ...
                    string(debug), saveStep, string(psfGen), string(saveZarr), string(parseCluster),string(parseParfor), ... 
                    string(GPUJob), string(Save16bit), string(largeFile), largeMethod, strrep(mat2str(BatchSize), ' ', ','), ...
                    strrep(mat2str(BlockSize), ' ', ','), strrep(mat2str(zarrSubSize), ' ', ','), damper, ...
                    scaleFactor_f, deconOffset, deconMaskFns_str, uuid, cpusPerTask, string(mccMode), ...
                    ConfigFile, GPUConfigFile);
            end
            
            if exist(dctmpFullpath, 'file') || parseCluster
                if parseCluster
                    dsz = getImageSize(dcframeFullpath);
                    estMem = prod(dsz) * 4 / 2^30;
                    cur_ConfigFile = ConfigFile;
                    if ~GPUJob
                        memAllocate = estMem * 10;
                        if largeFile && strcmp(largeMethod, 'inplace')
                            memAllocate = prod(BatchSize * 2) * 4 / 2^30 * 10;
                        end
                    else
                        memAllocate = estMem * 10;
                        if largeFile && strcmp(largeMethod, 'inplace')
                            memAllocate = prod(BatchSize) * 4 / 2^30 * 20;
                        else
                            cur_ConfigFile = GPUConfigFile;
                        end                        
                    end

                    job_id = job_ids(f, 1);
                    [job_id, ~, submit_status] = generic_single_job_submit_wrapper(func_str, ...
                        job_id, task_id, jobLogFname=job_log_fname, jobErrorFname=job_log_error_fname, ...
                        lastFile=false, cpusPerTask=cpusPerTask, memAllocate=memAllocate, ...
                        mccMode=mccMode, ConfigFile=cur_ConfigFile);

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
        continue;
    end
end

end

