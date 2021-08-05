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


ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('dataPaths'); % data structure from loadConditionData
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
ip.addParameter('RotateAfterDecon', false, @islogical);
% decon parameters
ip.addParameter('cudaDecon', false, @islogical);
ip.addParameter('cppDecon', ~false, @islogical);
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
ip.addParameter('fixIter', false, @islogical); % CPU Memory in Gb
ip.addParameter('errThresh', [], @isnumeric); % error threshold for simplified code
ip.addParameter('debug', false, @islogical); % debug mode for simplified code
ip.addParameter('GPUJob', false, @islogical); % use gpu for chuck deconvolution. 
% job related parameters
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('zarrFile', false, @islogical); % use zarr file as input
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @isstr);
ip.addParameter('cpusPerTask', 2, @isnumeric);
ip.addParameter('cpuOnlyNodes', true, @islogical);
ip.addParameter('uuid', '', @isstr);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 1, @isnumeric);
ip.addParameter('maxWaitLoopNum', 10, @isnumeric); % the max number of loops the loop waits with all existing files processed. 
ip.addParameter('MatlabLaunchStr', 'module load matlab/r2020b; matlab -nodisplay -nosplash -nodesktop -nojvm -r', @ischar);
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar);

ip.parse(dataPaths, varargin{:});

% make sure the function is in the root of XR_Repository. 
mpath = fileparts(which(mfilename));
repo_rt = [mpath, '/../../'];
cd(repo_rt);

pr = ip.Results;
Overwrite = pr.Overwrite;
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
RotateAfterDecon = pr.RotateAfterDecon;
RLMethod = pr.RLMethod;
GPUJob = pr.GPUJob;
% simplified version related options
fixIter = pr.fixIter;
errThresh = pr.errThresh;
debug = pr.debug;
% job related
largeFile = pr.largeFile;
zarrFile = pr.zarrFile;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
cpusPerTask = pr.cpusPerTask;
cpuOnlyNodes = pr.cpuOnlyNodes;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
MatlabLaunchStr = pr.MatlabLaunchStr;
SlurmParam = pr.SlurmParam;

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

if numel(Overwrite) == 1
    Overwrite = repmat(Overwrite, 1, 2);
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str, jobLogDir] = checkSlurmCluster(dataPath, jobLogDir, cpuOnlyNodes);
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
        
    deconPaths = cell(nd, 1);
    if RotateAfterDecon
        rdcPaths = cell(nd, 1);
    end

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
        
        if RotateAfterDecon
            rdcPath = [deconPath filesep 'Rotated' filesep];
            if Overwrite(2) && exist(rdcPath, 'dir')
                rmdir(rdcPath, 's');
            end
            if ~exist(rdcPath, 'dir')
                mkdir(rdcPath);
                fileattrib(rdcPath, '+w', 'g');
            end
            rdcPaths{d} = rdcPath;

            if rotatePSF 
                warning('The rotation is already performed before deconvolution! Please check the setting to make sure there is no duplicate rotation.');
            end
        end
    end
    if rotatePSF
        rotPSFFullpaths = cell(numel(psfFullpaths), 1);
    end

    for f = 1 : numel(psfFullpaths)
        if ~exist(psfFullpaths{f}, 'file')
            error('PSF file %s does not exist!', psfFullpaths{f});
        end
        if rotatePSF
            XR_rotate_PSF(psfFullpaths{f}, 'xyPixelSize', xyPixelSize, 'dz', dzPSF);
            [psfPath, fsname] = fileparts(psfFullpaths{f});
            rotPSFFullpaths{f} = [psfPath, '/Rotated/', fsname, '.tif'];
        end
    end
end

%% check existing files and parse channels
Streaming = false;
[fnames, fdinds, gfnames, partialvols, dataSizes, flipZstack_mat, FTP_inds, maskFullpaths] = ...
    XR_parseImageFilenames(dataPaths, ChannelPatterns, parseSettingFile, flipZstack, Decon, deconPaths, Streaming);


nF = numel(fnames);

% flags: for thee: deskew w/o rotate, decon w/o rotate, rotate
is_done_flag = false(nF, 2);
trial_counter = zeros(nF, 2);

% For no-streaming computing, first skip the disabled options. 
if ~Decon
    is_done_flag(:, 1) = true;
end
if ~RotateAfterDecon
    is_done_flag(:, 2) = true;
end

if parseCluster
    job_ids = -ones(nF, 2);
end

matlab_setup_str = 'setup([],true)';

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
        if ~exist(frameFullpath, 'file')
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
            deconFullpath = sprintf('%s/%s_decon.tif', deconPath, fsname);
            dctmpFullpath = sprintf('%s.tmp', deconFullpath(1 : end - 4));

            if exist(deconFullpath, 'file')
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
            if ErodeByFTP && ~cudaDecon
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
            
            flipZstack = flipZstack_mat(f);
            
            % do not use rotation in decon functions
            if cudaDecon
                func_str = sprintf(['XR_cudaDeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                    '''cudaDeconPath'',''%s'',''OTFGENPath'',''%s'',''dzPSF'',%.10f,''Background'',[%d],', ...
                    '''SkewAngle'',%d,''Rotate'',%s,''DeconIter'',%d,''Save16bit'',%s,''largeFile'',%s)'], ...
                    dcframeFullpath, xyPixelSize, dc_dz, psfFullpath, cudaDeconPath, OTFGENPath, dc_dzPSF, ...
                    Background, SkewAngle, string(deconRotate), DeconIter, string(Save16bit), string(largeFile));
            elseif cppDecon
                func_str = sprintf(['XR_cppDeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                    '''cppDeconPath'',''%s'',''loadModules'',''%s'',''dzPSF'',%.10f,''Background'',[%d],', ...
                    '''SkewAngle'',%d,''EdgeErosion'',%d,''ErodeMaskfile'',''%s'',''SaveMaskfile'',%s,', ...
                    '''Rotate'',%s,''DeconIter'',%d,''Save16bit'',%s,''largeFile'',%s)'], dcframeFullpath, ...
                    xyPixelSize, dc_dz, psfFullpath, cppDeconPath, loadModules, dc_dzPSF, Background, ...
                    SkewAngle, EdgeErosion, maskFullpath, string(SaveMaskfile), string(deconRotate), ...
                    DeconIter, string(Save16bit), string(largeFile));
            else
                func_str = sprintf(['XR_RLdeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                    '''dzPSF'',%.10f,''Background'',[%d],''SkewAngle'',%d,''flipZstack'',%s,''EdgeErosion'',%d,', ...
                    '''ErodeMaskfile'',''%s'',''SaveMaskfile'',%s,''Rotate'',%s,''DeconIter'',%d,', ...
                    '''RLMethod'',''%s'',''fixIter'',%s,''errThresh'',[%0.20f],''debug'',%s,''GPUJob'',%s,', ...
                    '''Save16bit'',%s,''largeFile'',%s)'], dcframeFullpath, xyPixelSize, dc_dz, psfFullpath, ...
                    dc_dzPSF, Background, SkewAngle, string(flipZstack), EdgeErosion, maskFullpath, string(SaveMaskfile), ...
                    string(deconRotate), DeconIter, RLMethod, string(fixIter), errThresh, string(debug), ...
                    string(GPUJob), string(Save16bit), string(largeFile));
            end
            
            if exist(dctmpFullpath, 'file') || parseCluster
                if parseCluster
                    job_status = check_slurm_job_status(job_ids(f, 1), task_id);

                     % if the job is still running, skip it. 
                    if job_status == 1 
                        continue;
                    end

                    if job_status == -1
                        % for matlab decon,  decide how many cores. 
                        if ~cudaDecon
                            [estMem, estGPUMem] = XR_estimateComputingMemory(dcframeFullpath, {'deconvolution'}, ...
                                'cudaDecon', false);
                            
                            cpusPerTask_dc = cpusPerTask;
                            if cpusPerTask_dc * 20 < estMem
                                cpusPerTask_dc = min(24, ceil(estMem / 20));
                            end
                        end

                        % do not use rotation in decon functions
                        matlab_cmd = sprintf('%s;tic;%s;toc', matlab_setup_str, func_str);
                        process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);

                        if cudaDecon
                            cmd = sprintf(['sbatch --array=%d -o %s -e %s -p abc --gres=gpu:1 --qos ', ...
                                'abc_normal -n1 --mem-per-cpu=33G --cpus-per-task=%d --wrap="%s"'], ...
                                task_id, job_log_fname, job_log_error_fname, 5, process_cmd);
                        else
                            cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                                '--wrap="echo Matlab command:  \\\"%s\\\"; %s"'], task_id, job_log_fname, ...
                                job_log_error_fname, cpusPerTask_dc, SlurmParam, slurm_constraint_str, ...
                                matlab_cmd, process_cmd);
                        end
                        [status, cmdout] = system(cmd, '-echo');

                        job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                        job_id = str2double(job_id{1}{1});
                        job_ids(f, 1) = job_id;
                        trial_counter(f, 1) = trial_counter(f, 1) + 1;
                    end
                else
                    temp_file_info = dir(dctmpFullpath);
                    if (datenum(clock) - [temp_file_info.datenum]) * 24 * 60 < unitWaitTime
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
            if exist(deconFullpath, 'file')
                is_done_flag(f, 1) = true;
                if exist(dctmpFullpath, 'file')
                    delete(dctmpFullpath);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rotation after decon
        if RotateAfterDecon
            % input chosen order is dsr, ds, raw (in decreased priority order)
            if ~exist(deconFullpath, 'file')
                continue;
            end

            rdcPath = rdcPaths{fdind};
            rdcFullpath = sprintf('%s/%s_decon.tif', rdcPath, fsname);
            rdctmpFullpath = sprintf('%s.tmp', rdcFullpath(1 : end - 4));
            if exist(rdcFullpath, 'file')
                is_done_flag(f, 2) = true;
                if exist(rdctmpFullpath, 'file')
                    delete(rdctmpFullpath);
                end
            end
        else
            is_done_flag(f, 2) = true;    
        end

        if ~is_done_flag(f, 2)
            func_str = sprintf(['XR_RotateFrame3D(''%s'',%.20d,%.20d,''SkewAngle'',%.20d,', ...
                '''ObjectiveScan'',%s,''Reverse'',%s,''Save16bit'',%s)'], deconFullpath, ...
                xyPixelSize, dz, SkewAngle, string(ObjectiveScan), string(Reverse), string(Save16bit));
            
            if exist(rdctmpFullpath, 'file') || parseCluster
                if parseCluster
                    job_status = check_slurm_job_status(job_ids(f, 2), task_id);
                    
                    % if the job is still running, skip it. 
                    if job_status == 1 
                        continue;
                    end
                    
                    if job_status == -1
                        [estMem, estGPUMem] = XR_estimateComputingMemory(deconFullpath, {'deconvolution'}, 'cudaDecon', false);
                        cpusPerTask_dcr = cpusPerTask;

                        if cpusPerTask * 20 < estMem
                            cpusPerTask = min(24, ceil(estMem / 20));
                        end

                        matlab_cmd = sprintf('%s;tic;%s;toc', matlab_setup_str, func_str);
                        process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
                        cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                            '--wrap="echo Matlab command:  \\\"%s\\\"; %s"'], ...
                            task_id, job_log_fname, job_log_error_fname, cpusPerTask_dcr, SlurmParam, ...
                            slurm_constraint_str, matlab_cmd, process_cmd);
                        [status, cmdout] = system(cmd, '-echo');

                        job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                        job_id = str2double(job_id{1}{1});
                        job_ids(f, 2) = job_id;
                        trial_counter(f, 2) = trial_counter(f, 2) + 1;    
                    end
                else
                    temp_file_info = dir(rdctmpFullpath);
                    if (datenum(clock) - [temp_file_info.datenum]) * 24 * 60 < unitWaitTime
                    else
                        fclose(fopen(rdctmpFullpath, 'w'));
                    end
                end
            else
                fclose(fopen(rdctmpFullpath, 'w'));
            end
            
            if ~parseCluster
                tic; feval(str2func(['@()', func_str])); toc;
                trial_counter(f, 2) = trial_counter(f, 2) + 1;    
            end
            
            % check if computing is done
            if exist(rdcFullpath, 'file')
                is_done_flag(f, 2) = true;
                if exist(rdctmpFullpath, 'file')
                    delete(rdctmpFullpath);
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

