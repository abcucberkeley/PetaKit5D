function [] = XR_microscopeAutomaticProcessing(dataPaths, varargin)
% automatic image processing pipeline for microscopy data. It can perform
% deskew w/o rotation, deconvolution, rotation afte deconvolution. It supports 
% both existing data in the directory or real-time transferring of data from
% other sources, i.e., microscope. Once all computing done, it will wait another 
% 5 minutes for potential new coming files. 
%
% Inputs :   
%            dataPath : directory path for the data
%
%
% Options (as 'specifier'-value pairs): 
%
%         'Overwrite' : true|{false}. Overwrite existing results.
%         'SkewAngle' : skew angle of the stage. Default: 31.5.
%                'dz' : stage scan interval. Default: 0.5.
%       'xyPixelSize' : pixel size. Default: 0.108.
%     'ObjectiveScan' : true|{false}. Objective scan. Default: 5.
%   'sCMOSCameraFlip' : true|{false}. sCMOS camera flip. 
%           'Reverse' : true|{false}. Inverse direction of z axis. 
%            'Deskew' : {true}|false. Deskew the data.
%            'Rotate' : {true}|false. Rotate deskewed data.
%             'Decon' : {true}|false. Deconvolution on data.
%         'cudaDecon' : {true}|false. Use cudaDecon for deconvolution. If there is no GPU, false by default. 
%                'DS' : true|{false}. Use deskewed data for deconvolution.
%               'DSR' : {true}|false. Use deskewed rotated data for deconvolution. 
%        'Background' : Background intensity for deconvolution. Default: 99 if not provided.
%       'deconRotate' : true|{false}. Rotate deconvolution result within deconvolution steps.
%  'RotateAfterDecon' : true|{false}. Rotate deconvolution results using matlab rotation function.
%      'psfFullpaths' : Full paths of psf files. A file for a channel. The order is the same as ChannelPatterns.
%          'Channels' : Channel wavelengths. 
%   'ChannelPatterns' : String patterns to recognize each channel. Currently, use Cam[A/B]_ch[0-9] as patterns. 
%         'Save16bit' : true|{false}. Correcting backgrounds across z-slices. 
%      'parseCluster' : Use slurm-based cluster computing.
%         'jobLogDir' : Log directory for the slurm jobs.
%       'cpusPerTask' : Number of cpus for a job. Default: 12
%              'uuid' : unique string for a job for saving files. 
%       'maxTrialNum' : Max number of times to rerun failure cases. 
%      'unitWaitTime' : The wait time per file in minutes to check whether the computing is done.
%     'minModifyTime' : The minimum time in minutes for the latest modified file to decide whether it is fully transferred.
%     'maxModifyTime' :  The maximum time in minutes to check whether there are coming new files.
%    'maxWaitLoopNum' : Number of maximum loops without any computing. Default: 10
%
%
% Author: Xiongtao Ruan (03/2020)
% 
% xruan (03/18/2020): add option of rotation after deconvolution; add
% option for skipping deskew.
% xruan (04/16/2020): add option for CPP decon
% xruan (07/13/2020): add option for matlab stitching, and update framework
%                     structure for each component
% xruan (07/17/2020): add Overwrite option
% xruan (07/18/2020): add option for support of multiple datasets (with same settings)
% xruan (07/19/2020): add streaming option, set Channel pattern to filter channels to process
% xruan (07/23/2020): add CPU node only option, add mask file option for
% the the first time point for edge erosion in Decon
% xruan (07/28/2020): fix issue in decon, when ErodeByFTP true and the deconm
%                     result for the first time point exist but the mask file 
%                     doesn't exist, in which it will stuck. 
%                     Also check if the image size match the mask size, if
%                     not use its own mask for edge erosion.
% xruan (08/03/3030): add axis-order for the stitching
% xruan (08/07/3030): change file attributes for folders so that group users have write access
% xruan (08/08/3030): add more memory for per-cpu, fix issue for support of 
%                     multiple datasets in stitching; set xcorrShift as
%                     false by default



ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths'); % data structure from loadConditionData
ip.addParameter('Overwrite', false,  @(x) (numel(x) == 1 || numel(x) == 5) && islogical(x));
ip.addParameter('Streaming', true,  @islogical); % if true, check for new files. If false, assume all files transferred completely.
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('SkewAngle', 32.45, @isscalar);
ip.addParameter('dz', 0.5, @isscalar);
ip.addParameter('xyPixelSize', 0.108, @isscalar);
ip.addParameter('Reverse', true, @islogical);
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('sCMOSCameraFlip', false, @islogical);
% pipeline steps
ip.addParameter('Deskew', true, @islogical);
ip.addParameter('Rotate', true, @islogical);
ip.addParameter('Stitch', false, @islogical);
ip.addParameter('Decon', ~false, @islogical);
ip.addParameter('RotateAfterDecon', false, @islogical);
% deskew and rotation options
ip.addParameter('LLFFCorrection', false, @islogical);
ip.addParameter('LowerLimit', 0.4, @isnumeric); % this value is the lowest
ip.addParameter('LSImagePaths', {'','',''}, @iscell);
ip.addParameter('BackgroundPaths', {'','',''}, @iscell);
% stitch parameters
ip.addParameter('stitchResultDir', 'matlab_stitch', @islogical);
ip.addParameter('imageListFullpaths', '', @(x) ischar(x) || iscell(x));
ip.addParameter('axisOrder', 'xyz', @(x) ischar(x));
ip.addParameter('BlendMethod', 'none', @isstr);
ip.addParameter('xcorrShift', false, @islogical);
ip.addParameter('xcorrMode', 'primaryFirst', @(x) strcmpi(x, 'primary') || strcmpi(x, 'primaryFirst') || strcmpi(x, 'all')); % 'primary': choose one channel as primary channel, 
ip.addParameter('boundboxCrop', [], @(x) isnumeric(x) && (isempty(x) || all(size(x) == [3, 2]) || numel(x) == 6));
% decon parameters
ip.addParameter('cudaDecon', false, @islogical);
ip.addParameter('cppDecon', ~false, @islogical);
ip.addParameter('DS', true, @islogical);
ip.addParameter('DSR', false, @islogical);
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('dzPSF', 0.1, @isnumeric);
ip.addParameter('EdgeErosion', 8, @isnumeric);
ip.addParameter('ErodeByFTP', true, @islogical); % Edge erosion by the first time point (ranked the first in the inital file list for each dataset).
ip.addParameter('deconRotate', false, @islogical);
% ip.addParameter('PSFPath', '', @isstr);
ip.addParameter('psfFullpaths', {'','',''}, @iscell);
ip.addParameter('Channels', [488, 560, 642], @isnumeric);
ip.addParameter('Save16bit', [false, false, false, false], @(x) (numel(x) == 1 || numel(x) == 4) && islogical(x));
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @isstr);
ip.addParameter('cpusPerTask', 2, @isnumeric);
ip.addParameter('cpuOnlyNodes', true, @islogical);
ip.addParameter('uuid', '', @isstr);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 1, @isnumeric);
ip.addParameter('minModifyTime', 1, @isnumeric); % the minimum during of last modify time of a file, in minute.
ip.addParameter('maxModifyTime', 10, @isnumeric); % the maximum during of last modify time of a file, in minute.
ip.addParameter('maxWaitLoopNum', 10, @isnumeric); % the max number of loops the loop waits with all existing files processed. 

ip.parse(dataPaths, varargin{:});

% make sure the function is in the root of XR_Repository. 
repo_rt = fileparts(which(mfilename));
cd([repo_rt, '/..']);

pr = ip.Results;
Overwrite = pr.Overwrite;
Streaming = pr.Streaming;
% Resolution = pr.Resolution;
SkewAngle = pr.SkewAngle;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
ObjectiveScan = pr.ObjectiveScan;
LLFFCorrection = pr.LLFFCorrection;
LowerLimit = pr.LowerLimit;
LSImagePaths = pr.LSImagePaths;
BackgroundPaths = pr.BackgroundPaths;
Reverse = pr.Reverse;
Deskew = pr.Deskew;
Rotate = pr.Rotate;
% stitch parameters
Stitch = pr.Stitch;
stitchResultDir = pr.stitchResultDir;
imageListFullpaths = pr.imageListFullpaths;
axisOrder = pr.axisOrder;
BlendMethod = pr.BlendMethod;
xcorrShift = pr.xcorrShift;
xcorrMode = pr.xcorrMode;
boundboxCrop = pr.boundboxCrop;
% decon parameters
Decon = pr.Decon;
cppDecon = pr.cppDecon;
cudaDecon = pr.cudaDecon;
EdgeErosion = pr.EdgeErosion;
ErodeByFTP = pr.ErodeByFTP;
DS = pr.DS;
DSR = pr.DSR;
Background = pr.Background;
dzPSF = pr.dzPSF;
psfFullpaths = pr.psfFullpaths;
deconRotate = pr.deconRotate;
RotateAfterDecon = pr.RotateAfterDecon;
ChannelPatterns = pr.ChannelPatterns;
largeFile = pr.largeFile;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
cpusPerTask = pr.cpusPerTask;
cpuOnlyNodes = pr.cpuOnlyNodes;
Save16bit = pr.Save16bit;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
minModifyTime = pr.minModifyTime;
maxModifyTime = pr.maxModifyTime;
maxWaitLoopNum = pr.maxWaitLoopNum;

% suppress directory exists warning
warning('off', 'MATLAB:MKDIR:DirectoryExists');

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

if ischar(imageListFullpaths)
    imageListFullpaths = {imageListFullpaths};
end

nd = numel(dataPaths);
for d = 1 : nd
    dataPath = dataPaths{d};
    if ~strcmp(dataPath(end), filesep)
        dataPaths{d} = [dataPath, filesep];
    end
end

if numel(Overwrite) == 1
    Overwrite = repmat(Overwrite, 1, 5);
end

% check if a slurm-based computing cluster exists
if parseCluster 
    [status, ~] = system('sinfo');
    clusterAvailable = true;
    if status ~= 0
        warning('A slurm-based computing cluster is not exist. Set parseCluster as false.')
        parseCluster = false;
        clusterAvailable = false;
    end
    if parseCluster && ~exist(jobLogDir, 'dir')
        warning('The job log directory does not exist, use %s/job_logs as job log directory', dataPath)
        jobLogDir = [dataPath, '/tmp'];
        mkdir(jobLogDir);
        fileattrib(jobLogDir, '+w', 'g');        
    end
    job_log_fname = [jobLogDir, '/job_%A_%a.out'];
    job_log_error_fname = [jobLogDir, '/job_%A_%a.err'];
    
    % For cudaDecon, not applicable.
    if cpuOnlyNodes
        slurm_constraint_str = ' --constraint=c24 ';
    else
        slurm_constraint_str = '';
    end
end

% for stitching, enable DS and DSR
% For multiple datasets, if the image list is not provide for any dataset,
% set Stitch as false. 
if Stitch
    if numel(imageListFullpaths) ~= nd
        warning("The number of image list files does not match that of the data paths, please make sure image list files are provided for each dataset!")
        Stitch = false;
    else
        for d = 1 : nd
            if ~exist(imageListFullpaths{d}, 'file')
                warning('Image list filename %s does not exist, set Stitch option as False', imageListFullpaths{d});
                Stitch = false;
                break;
            end
        end
    end
    
    Deskew = true;
    Rotate = true;
    stchPaths = cell(nd, 1);
    for d = 1 : nd
        dataPath = dataPaths{d};
        stchPath = [dataPath, '/matlab_stitch/'];
        if Overwrite(3) && exist(stchPath, 'dir')
            rmdir(stchPath, 's');
        end
        mkdir(stchPath);
        fileattrib(stchPath, '+w', 'g');
        stchPaths{d} = stchPath;
    end
    % DSRDirstr = '/DSR/';
    if Decon && ErodeByFTP
        ErodeByFTP = false;
    end
    
    % check if axis order is valid
    axisOrder = strrep(axisOrder, ' ', '');
    pattern = '^(-?x,?-?y,?-?z|-?y,?-?x,?-?z|-?z,?-?y,?-?x|-?x,?-?z,?-?y|-?x,?-?z,?-?y|-?y,?-?z,?-?x)$';
    if ~regexpi(axisOrder, pattern)
        error("The axisOrder is not right, it must has the form like 'y,x,z' or '-x,y,z' (flipped in x-axis)!");
    end
end

% first check if DS and Deconvolution directories exist
if Deskew
    dsPaths = cell(nd, 1);
    for d = 1 : nd
        dataPath = dataPaths{d};
        dsPath = [dataPath, '/DS/'];
        if Overwrite(1) && exist(dsPath, 'dir')
            rmdir(dsPath, 's');
        end    
        mkdir(dsPath);
        fileattrib(dsPath, '+w', 'g');
        dsPaths{d} = dsPath;
    end
    
    % check LS and Background images for flat field correction
    if LLFFCorrection
        if numel(ChannelPatterns) ~= numel(LSImagePaths) || numel(ChannelPatterns) ~= numel(BackgroundPaths) 
            error('The number of channels in ChannelPatterns does not match the number of LS files or Background Files!')
        end
        for c = 1 : numel(LSImagePaths)
            if ~exist(LSImagePaths{c}, 'file')
                error('LS Image file %s does not exist!', LSImagePaths{c});
            end
            if ~exist(BackgroundPaths{c}, 'file')
                error('LS Image file %s does not exist!', BackgroundPaths{c});
            end
        end
    end
else
    Rotate = false;
    DS = false;
    DSR = false;
end
    
if Rotate
    dsrPaths = cell(nd, 1);
    for d = 1 : nd
        dataPath = dataPaths{d};
        dsrPath = [dataPath, '/DSR/'];
        if Overwrite(2) && exist(dsrPath, 'dir')
            rmdir(dsrPath, 's');
        end    
        mkdir(dsrPath);
        fileattrib(dsrPath, '+w', 'g');        
        dsrPaths{d} = dsrPath;
    end
end

% for deconvolution, check whether there is a gpu in the node. if not, for
% cudaDecon, set parseCluster as true. 
if Decon
    % if both cudaDecon and cppDecon are true, use cppDecon
    if cudaDecon && cppDecon
        cudaDecon = false;
    end

    if cudaDecon && gpuDeviceCount() < 1 && ~clusterAvailable
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
        if DS
            dsPath = dsPaths{d};
            deconPath = [dsPath, deconName, filesep];
        end
        if DSR
            dsrPath = dsrPaths{d};
            deconPath = [dsrPath, deconName, filesep];
        end
        if Stitch
            stchPath = stchPaths{d};            
            deconPath = [stchPath, deconName, filesep];
            RotateAfterDecon = false;
        end

        if Overwrite(4) && exist(deconPath, 'dir')
            rmdir(deconPath, 's');
        end    
        mkdir(deconPath);
        fileattrib(deconPath, '+w', 'g');                
        deconPaths{d} = deconPath;
        
        if RotateAfterDecon
            rdcPath = [deconPath filesep 'Rotated' filesep];
            if Overwrite(5) && exist(rdcPath, 'dir')
                rmdir(rdcPath, 's');
            end        
            mkdir(rdcPath);
            fileattrib(rdcPath, '+w', 'g');            
            rdcPaths{d} = rdcPath;

            if DSR 
                warning('The rotation is already performed before deconvolution! Please check the setting to make sure there is no duplicate rotation.');
            end
        end
    end
    
    % check whether a psf file exist, if not, show a warning
    if DSR || Stitch
        rotPSFFullpaths = cell(numel(psfFullpaths), 1);
    end
        
    for f = 1 : numel(psfFullpaths)
        if ~exist(psfFullpaths{f}, 'file')
            error('PSF file %s does not exist!', psfFullpaths{f});
        end
        if DSR || Stitch
            XR_rotate_PSF(psfFullpaths{f});
        end
        [psfPath, fsname] = fileparts(psfFullpaths{f});
        rotPSFFullpaths{f} = [psfPath, '/Rotated/', fsname, '.tif'];
    end
end

%% check existing files and parse channels
fnames_cell = cell(nd, 1);
for d = 1 : nd
    dataPath = dataPaths{d};
    dir_info = dir([dataPath, '*.tif']);
    fnames_d = {dir_info.name}';
    if Streaming
        last_modify_time = (datenum(clock) - [dir_info.datenum]) * 24 * 60;
        latest_modify_time = min(last_modify_time);

        % not include the lastest file if it is very recent
        if latest_modify_time < minModifyTime
            fnames_d(last_modify_time == latest_modify_time) = [];
        end
    end
    fnames_cell{d} = fnames_d;
end

fdinds = arrayfun(@(x) ones(numel(fnames_cell{x}), 1) * x, 1 : nd, 'unif', 0);
fnames = cat(1, fnames_cell{:});
fdinds = cat(1, fdinds{:});

% filter filenames by channel patterns
include_flag = false(numel(fnames), 1);
for c = 1 : numel(ChannelPatterns)
    include_flag = include_flag | contains(fnames, ChannelPatterns{c});
end
fnames = fnames(include_flag);
fdinds = fdinds(include_flag);

nF = numel(fnames);

% for ErodeByFTP, set up the frame numbers as first time point for each
% data
if Decon && ErodeByFTP
    FTP_inds = zeros(nd, 1);
    maskFullpaths = cell(nd, 1);
    for d = 1 : nd
        c = 1;
        FTPfname = '';
        while isempty(FTPfname)
            all_inds = contains(fnames_cell{d}, ChannelPatterns{c});
            FTPfname = fnames_cell{d}{find(all_inds, 1, 'first')};
            c = c + 1;
        end
        if ~isempty(FTPfname)
            ind_d = find(strcmp(fnames, FTPfname));
            FTP_inds(d) = ind_d;
            maskFullpaths{d} = sprintf('%s/Masks/%s_eroded.tif', deconPaths{d}, FTPfname(1 : end - 4));
        end
    end    
end

% flags: for thee: deskew w/o rotate, decon w/o rotate, rotate
is_done_flag = false(nF, 4);
trial_counter = zeros(nF, 4);
waitLoopCounter = 0;

% For no-streaming computing, first skip the disabled options. 
if ~Streaming
    if ~Stitch
        is_done_flag(:, 2) = true;
    end
    if ~Decon
        is_done_flag(:, 3) = true;
    end
    if ~RotateAfterDecon
        is_done_flag(:, 4) = true;
    end
end

if parseCluster
    job_ids = -ones(nF, 4);
end

% use while loop to perform computing for all images
while ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') || ...
        (Streaming && (nF == 0 || latest_modify_time < maxModifyTime || waitLoopCounter < maxWaitLoopNum))
    for f = 1 : nF
        if all(is_done_flag(f, :))
            continue;
        end
        
        % first deskew and rotation
        fname = fnames{f};
        fdind = fdinds(f);
        dataPath = dataPaths{fdind};
        
        frameFullpath = [dataPath, fname];
        % check wheter the file is deleted during the computing.
        if ~exist(frameFullpath, 'file')
            is_done_flag(f, :) = true;
            continue
        end
        
        task_id = f;
        
        %% deskew w/o rotate
        if Deskew
            dsPath = dsPaths{fdind};
            dsFullpath = [dsPath, fname];
            if Rotate
                dsrPath = dsrPaths{fdind};
                dsrFullpath = [dsrPath, fname];
            end
            tmpFullpath = sprintf('%s.tmp', dsFullpath(1 : end - 4));

            if exist(dsFullpath, 'file') && (~Rotate || exist(dsrFullpath, 'file'))
                is_done_flag(f, 1) = true;
                if exist(tmpFullpath, 'file')
                    delete(tmpFullpath);
                end
            end
        else
            is_done_flag(f, 1) = true;
        end
    
        if ~is_done_flag(f, 1)
            if LLFFCorrection
                LLFFMapping =  ~cellfun(@isempty, regexpi(fname, ChannelPatterns));
                LSImage = LSImagePaths{LLFFMapping};
                BackgroundImage = BackgroundPaths{LLFFMapping};
            else
                LSImage = '';
                BackgroundImage = '';
            end
                
            if exist(tmpFullpath, 'file') || parseCluster
                if parseCluster
                    job_status = check_slurm_job_status(job_ids(f, 1), task_id);
                    
                    % if the job is still running, skip it. 
                    if job_status == 1 
                        continue;
                    end
                    
                    if job_status == -1
                        % first estimate file size and decide whether cpusPerTask
                        % is enough
                        estRequiredMemory = XR_estimateComputingMemory(frameFullpath, {'deskew'});
                        if cpusPerTask * 20 < estRequiredMemory
                            cpusPerTask = min(24, ceil(estRequiredMemory / 20));
                        end
                        
                        matlab_cmd = sprintf(['addpath(genpath(''../XR_GU_Repository'')),addpath(genpath(pwd));', ...
                            'tic;XR_deskewRotateFrame(''%s'',%.20d,%.20d,''SkewAngle'',%.20d,''ObjectiveScan'',%s,', ...
                            '''Reverse'',%s,''LLFFCorrection'',%s,''LowerLimit'',%.20d,''LSImage'',''%s'',', ...
                            '''BackgroundImage'',''%s'',''Rotate'',%s,''Save16bit'',%s);toc;'], ...
                            frameFullpath, xyPixelSize, dz, SkewAngle, string(ObjectiveScan), string(Reverse), ...
                            string(LLFFCorrection), LowerLimit, LSImage, BackgroundImage, string(Rotate), ...
                            string(Save16bit(1)));
                        deskew_cmd = sprintf('module load matlab/r2020a; matlab -nodisplay -nosplash -nodesktop -nojvm -r \\"%s\\"', matlab_cmd);
                        cmd = sprintf('sbatch --array=%d %s -o %s -e %s -p abc --qos abc_normal -n1 --mem-per-cpu=20G --cpus-per-task=%d --wrap="%s"', ...
                            task_id, slurm_constraint_str, job_log_fname, job_log_error_fname, cpusPerTask, deskew_cmd);
                        [status, cmdout] = system(cmd, '-echo');

                        job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                        job_id = str2double(job_id{1}{1});
                        job_ids(f, 1) = job_id;
                        trial_counter(f, 1) = trial_counter(f, 1) + 1;
                    end
                else
                    temp_file_info = dir(tmpFullpath);
                    if (datenum(clock) - [temp_file_info.datenum]) * 24 * 60 < unitWaitTime
                        continue; 
                    else
                        fclose(fopen(cur_tmp_fname, 'w'));
                    end
                end
            else
                fclose(fopen(tmpFullpath, 'w'));
            end
            if ~parseCluster
                XR_deskewRotateFrame(frameFullpath, xyPixelSize, dz, 'SkewAngle', SkewAngle, ...
                    'ObjectiveScan', ObjectiveScan, 'Reverse', Reverse, 'LLFFCorrection', LLFFCorrection, ...
                    'LowerLimit', LowerLimit, 'LSImage', LSImage, 'BackgroundImage', BackgroundImage, ...
                    'Rotate', Rotate, 'Save16bit', Save16bit(1), 'uuid', uuid);
                trial_counter(f, 1) = trial_counter(f, 1) + 1;
            end

            % check if computing is done
            if exist(dsFullpath, 'file') && (~Rotate || exist(dsrFullpath, 'file'))
                is_done_flag(f, 1) = true;
                if exist(tmpFullpath, 'file')
                    delete(tmpFullpath);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% then stitching
        if Stitch
            stchPath = stchPaths{fdind};
            dir_info = dir(sprintf('%s/%s*Abs.tif', stchPath, fname(regexp(fname, 'Scan.*') : end - 43)));
            if ~isempty(dir_info) 
                stch_fname = dir_info(1).name;
                stchFullpath = [stchPath, stch_fname];
                stchTmpFullpath = sprintf('%s.tmp', stchFullpath(1 : end - 4));
                if exist(stchFullpath, 'file')
                    is_done_flag(f, 2) = true;
                    if exist(stchTmpFullpath, 'file')
                        delete(stchTmpFullpath);
                    end
                end
            end
        else
            is_done_flag(f, 2) = true;
        end
        
        if ~is_done_flag(f, 2)
            % if the deskew result is not available, wait for the deskew data
            % to be finished        
            if ~parseCluster && (~exist(frameFullpath, 'file') || ~exist(dsrFullpath, 'file'))
                % if DS or DSR, it means the deskew and rotation not done.
                continue;
            end
            
            useExistDSR = true;
            resampleType = 'isotropic';
            zNormalize = false;
            onlyFirstTP = false;
            bbox = [];
            imageListFullpath = imageListFullpaths{fdind};

            if parseCluster
                dfirst_ind = find(fdinds == fdind, 1, 'first');
                job_status = check_slurm_job_status(job_ids(dfirst_ind, 2), 1);

                if job_status == -1
                    % first estimate file size and decide whether cpusPerTask
                    % is enough
                    estRequiredMemory = XR_estimateComputingMemory(frameFullpath, {'deskew'});
                    if cpusPerTask * 20 < estRequiredMemory
                        cpusPerTask = min(24, ceil(estRequiredMemory / 20));
                    end

                    matlab_cmd = sprintf(['addpath(genpath(''../XR_GU_Repository''));', ...
                        'addpath(genpath(pwd));tic;XR_matlab_stitching_wrapper(''%s'',''%s'',' ...
                        '''ResultDir'',''%s'',''Streaming'',%s,''useExistDSR'',%s,''axisOrder'',''%s'',', ...
                        '''resampleType'',''%s'',''Reverse'',%s,''xcorrShift'',%s,''xcorrMode'',''%s'',', ...
                        '''BlendMethod'',''%s'',''zNormalize'',%s,''onlyFirstTP'',%s,''boundboxCrop'',[%s]);toc;'], ...
                        dataPath, imageListFullpath, stitchResultDir, string(Streaming), string(useExistDSR), ...
                        axisOrder, resampleType, string(Reverse), string(xcorrShift), xcorrMode, BlendMethod, ...
                        string(zNormalize), string(onlyFirstTP), strrep(num2str(bbox, '%d,'), ' ', ''));
                    stitch_cmd = sprintf('module load matlab/r2020a; matlab -nodisplay -nosplash -nodesktop -r \\"%s\\"', matlab_cmd);
                    cmd = sprintf(['sbatch --array=1 %s -o %s -e %s -p abc', ...
                        ' --qos abc_normal -n1 --mem-per-cpu=21418M --cpus-per-task=%d --wrap="%s"'], ...
                        slurm_constraint_str, job_log_fname, job_log_error_fname, cpusPerTask, stitch_cmd);            
                    [status, cmdout] = system(cmd, '-echo');

                    job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                    job_id = str2double(job_id{1}{1});
                    job_ids(dfirst_ind, 2) = job_id;
                end
            else
                XR_matlab_stitching_wrapper(dataPath, imageListFullpath, 'resultDir', stitchResultDir, ...
                    'Streaming', Streaming, 'useExistDSR', useExistDSR, 'axisOrder', axisOrder, ...
                    'resampleType', 'isotropic', 'Reverse', Reverse, 'BlendMethod', BlendMethod, ...
                    'xcorrShift', xcorrShift, 'xcorrMode', xcorrMode, 'zNormalize', zNormalize, ...
                    'onlyFirstTP', onlyFirstTP, 'boundboxCrop', bbox, 'parseCluster', false);
            end
            
            if exist('stchFullpath', 'var') && exist(stchFullpath, 'file')
                is_done_flag(f, 2) = true;
                if exist(stchTmpFullpath, 'file')
                    delete(stchTmpFullpath);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% then deconvolution
        if Decon
            % if the deskew result is not available, wait for the deskew data
            % to be finished

            % input chosen order is dsr, ds, raw (in decreased priority order)
            dcframeFullpath = frameFullpath;
            if DS
                dcframeFullpath = dsFullpath;
            end
            if DSR
                dcframeFullpath = dsrFullpath;
            end 

            dc_dz = dz;
            dc_dzPSF = dzPSF;
            dc_psfFullpaths = psfFullpaths;
            if Stitch
                % in case fname not start with Scan_
                dir_info = dir(sprintf('%s/%s*Abs.tif', stchPath, fname(regexp(fname, 'Scan.*') : end - 43)));
                if isempty(dir_info) 
                    continue;
                end
                stch_fname = dir_info(1).name;
                dcframeFullpath = [stchPath, stch_fname];
                dc_dz = xyPixelSize;
                dc_dzPSF = xyPixelSize;
                dc_psfFullpaths = rotPSFFullpaths;
                
                % if stitch, only do computing when it is for the first
                % tile, to avoid replicate computing
                if ~contains(fname, stch_fname(1 : end - 4))
                    is_done_flag(f, 3) = true;
                    continue;
                end
            end
            
            if ~exist(dcframeFullpath, 'file')
                % if DS or DSR, it means the deskew and rotation not done.
                if DS || DSR || Stitch
                    continue;
                end
            end
                            
            deconPath = deconPaths{fdind};
            deconFullpath = sprintf('%s/%s_decon.tif', deconPath, fname(1 : end - 4));
            if Stitch
                deconFullpath = sprintf('%s/%s_decon.tif', deconPath, stch_fname(1 : end - 4));
            end
            dctmpFullpath = sprintf('%s.tmp', deconFullpath(1 : end - 4));

            if exist(deconFullpath, 'file')
                is_done_flag(f, 3) = true;
                if exist(dctmpFullpath, 'file')
                    delete(dctmpFullpath);
                end
            end
            
            % for ErodeByFTP, check if the mask file exist
            maskFullpath = '';
            SaveMaskfile = false;
            if ErodeByFTP
                FTP_ind = FTP_inds(fdind);
                if f == FTP_ind
                    SaveMaskfile = true;
                    % if decon result exist, but mask file not exist, rerun
                    % it to save the mask. 
                    if is_done_flag(f, 3) && ~exist(maskFullpaths{fdind}, 'file')
                        is_done_flag(f, 3) = false;
                        delete(deconFullpath);
                    end
                else
                    % only check for the ones not finished. 
                    if ~is_done_flag(f, 3)
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
            is_done_flag(f, 3) = true;    
        end

        if ~is_done_flag(f, 3)            
            if exist(dctmpFullpath, 'file') || parseCluster
                psfMapping =  ~cellfun(@isempty, regexpi(fname, ChannelPatterns));
                psfFullpath = dc_psfFullpaths{psfMapping};

                if parseCluster
                    job_status = check_slurm_job_status(job_ids(f, 3), task_id);

                     % if the job is still running, skip it. 
                    if job_status == 1 
                        continue;
                    end

                    if job_status == -1
                        % for matlab decon,  decide how many cores. 
                        if ~cudaDecon
                            [estMem, estGPUMem] = XR_estimateComputingMemory(dsFullpath, {'deconvolution'}, 'cudaDecon', false);

                            if cpusPerTask * 20 < estMem
                                cpusPerTask = min(24, ceil(estMem / 20));
                            end
                        end

                        % do not use rotation in decon functions
                        if cudaDecon
                            matlab_cmd = sprintf(['addpath(genpath(''../XR_GU_Repository'')),addpath(genpath(pwd));tic;', ...
                                'XR_cudaDeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                                '''dzPSF'',%.10f,''Background'',[%d],''SkewAngle'',%d,', ...
                                '''Rotate'',%s,''Save16bit'',%s,''largeFile'',%s);toc;'], ...
                                dcframeFullpath, xyPixelSize, dc_dz, psfFullpath, ...
                                dc_dzPSF, Background, SkewAngle, string(deconRotate), string(Save16bit(4)), string(largeFile));
                            decon_cmd = sprintf('module load matlab/r2020a; matlab -nodisplay -nosplash -nodesktop -nojvm -r \\"%s\\"', matlab_cmd);
                            cmd = sprintf('sbatch --array=%d -o %s -e %s -p abc --gres=gpu:1 --qos abc_normal -n1 --mem-per-cpu=33G --cpus-per-task=%d --wrap="%s"', ...
                                task_id, job_log_fname, job_log_error_fname, 5, decon_cmd);
                        elseif cppDecon
                            matlab_cmd = sprintf(['addpath(genpath(''../XR_GU_Repository'')),addpath(genpath(pwd));tic;', ...
                                'XR_cppDeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                                '''dzPSF'',%.10f,''Background'',[%d],''SkewAngle'',%d,', ...
                                '''EdgeErosion'',%d,''ErodeMaskfile'',''%s'',''SaveMaskfile'',%s,', ...
                                '''Rotate'',%s,''Save16bit'',%s,''largeFile'',%s);toc;'], ...
                                dcframeFullpath, xyPixelSize, dc_dz, psfFullpath, ...
                                dc_dzPSF, Background, SkewAngle, EdgeErosion, maskFullpath, string(SaveMaskfile), ...
                                string(deconRotate), string(Save16bit(4)), string(largeFile));
                            decon_cmd = sprintf('module load matlab/r2020a; matlab -nodisplay -nosplash -nodesktop -nojvm -r \\"%s\\"', matlab_cmd);
                            cmd = sprintf('sbatch --array=%d %s -o %s -e %s -p abc --qos abc_normal -n1 --mem-per-cpu=21418M --cpus-per-task=%d --wrap="%s"', ...
                                task_id, slurm_constraint_str, job_log_fname, job_log_error_fname, cpusPerTask, decon_cmd);                        
                        else
                            matlab_cmd = sprintf(['addpath(genpath(''../GU_PrivateRepository'')),addpath(genpath(''../XR_GU_Repository'')),addpath(genpath(pwd));tic;', ...
                                'XR_RLdeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                                '''dzPSF'',%.10f,''Background'',[%d],''SkewAngle'',%d,', ...
                                '''Rotate'',%s,''Save16bit'',%s);toc;'], ...
                                dcframeFullpath, xyPixelSize, dc_dz, psfFullpath, ...
                                dc_dzPSF, Background, SkewAngle, string(deconRotate), string(Save16bit(4)));
                            decon_cmd = sprintf('module load matlab/r2020a; matlab -nodisplay -nosplash -nodesktop -r \\"%s\\"', matlab_cmd);
                            cmd = sprintf('sbatch --array=%d %s -o %s -e %s -p abc --qos abc_normal -n1 --mem-per-cpu=21418M --cpus-per-task=%d --wrap="%s"', ...
                                task_id, slurm_constraint_str, job_log_fname, job_log_error_fname, cpusPerTask, decon_cmd);
                        end
                        [status, cmdout] = system(cmd, '-echo');

                        job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                        job_id = str2double(job_id{1}{1});
                        job_ids(f, 3) = job_id;
                        trial_counter(f, 3) = trial_counter(f, 3) + 1;
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
                if cudaDecon
                    XR_cudaDeconFrame3D(dcframeFullpath, xyPixelSize, dc_dz, '', ...
                        'PSFfile', psfFullpath, 'dzPSF', dc_dzPSF, 'Background', Background, ...
                        'SkewAngle', SkewAngle, 'Rotate', deconRotate, 'Save16bit', Save16bit(4));
                elseif cppDecon
                    XR_cppDeconFrame3D(dcframeFullpath, xyPixelSize, dc_dz, '', ...
                        'PSFfile', psfFullpath, 'dzPSF', dc_dzPSF, 'Background', Background, ...
                        'SkewAngle', SkewAngle, 'EdgeErosion', EdgeErosion, 'ErodeMaskfile', maskFullpath, ...
                        'SaveMaskfile', SaveMaskfile, 'Rotate', deconRotate, 'Save16bit', Save16bit(4));
                else
                    XR_RLdeconFrame3D(dcframeFullpath, xyPixelSize, dc_dz, '', ...
                        'PSFfile', psfFullpath, 'dzPSF', dc_dzPSF, 'Background', Background, ...
                        'SkewAngle', SkewAngle, 'Rotate', deconRotate, 'Save16bit', Save16bit(4));
                end
                trial_counter(f, 3) = trial_counter(f, 3) + 1;
            end
            
            % check if computing is done
            if exist(deconFullpath, 'file')
                is_done_flag(f, 3) = true;
                if exist(dctmpFullpath, 'file')
                    delete(dctmpFullpath);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% rotation after decon
        if RotateAfterDecon
            % input chosen order is dsr, ds, raw (in decreased priority order)
            if ~exist(deconFullpath, 'file')
                if Decon
                    continue;
                end
            end

            rdcPath = rdcPaths{fdind};
            rdcFullpath = sprintf('%s/%s_decon.tif', rdcPath, fname(1 : end - 4));
            rdctmpFullpath = sprintf('%s.tmp', rdcFullpath(1 : end - 4));
            if exist(rdcFullpath, 'file')
                is_done_flag(f, 4) = true;
                if exist(rdctmpFullpath, 'file')
                    delete(rdctmpFullpath);
                end
            end
        else
            is_done_flag(f, 4) = true;    
        end

        if ~is_done_flag(f, 4)
            if exist(rdctmpFullpath, 'file') || parseCluster
                if parseCluster
                    job_status = check_slurm_job_status(job_ids(f, 4), task_id);
                    
                    % if the job is still running, skip it. 
                    if job_status == 1 
                        continue;
                    end
                    
                    if job_status == -1
                        [estMem, estGPUMem] = XR_estimateComputingMemory(deconFullpath, {'deconvolution'}, 'cudaDecon', false);

                        if cpusPerTask * 20 < estMem
                            cpusPerTask = min(24, ceil(estMem / 20));
                        end

                        matlab_cmd = sprintf('addpath(genpath(''../LLSM3DTools''));addpath(genpath(pwd));tic;XR_RotateFrame3D(''%s'',%.20d,%.20d,''SkewAngle'',%.20d,''ObjectiveScan'',%s,''Reverse'',%s,''Save16bit'',%s);toc;', ...
                            deconFullpath, xyPixelSize, dz, SkewAngle, string(ObjectiveScan), string(Reverse), string(Save16bit(4)));
                        decon_cmd = sprintf('module load matlab/r2020a; matlab -nodisplay -nosplash -nodesktop -nojvm -r \\"%s\\"', matlab_cmd);
                        cmd = sprintf('sbatch --array=%d %s -o %s -e %s -p abc --qos abc_normal -n1 --mem-per-cpu=21418M --cpus-per-task=%d --wrap="%s"', ...
                            task_id, slurm_constraint_str, job_log_fname, job_log_error_fname, cpusPerTask, decon_cmd);
                        [status, cmdout] = system(cmd, '-echo');

                        job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                        job_id = str2double(job_id{1}{1});
                        job_ids(f, 4) = job_id;
                        trial_counter(f, 4) = trial_counter(f, 4) + 1;    
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
                XR_RotateFrame3D(deconFullpath, xyPixelSize, dz, 'SkewAngle', SkewAngle, ...
                    'ObjectiveScan', ObjectiveScan, 'Reverse', Reverse, 'Save16bit', Save16bit(4), 'uuid', uuid);
                trial_counter(f, 4) = trial_counter(f, 4) + 1;    
            end
            
            % check if computing is done
            if exist(rdcFullpath, 'file')
                is_done_flag(f, 4) = true;
                if exist(rdctmpFullpath, 'file')
                    delete(rdctmpFullpath);
                end
            end
        end
    end
    
    %% wait for running jobs finishing and checking for new coming images
    if ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') 
        waitLoopCounter = 0;
        pause(30);
        % if not streaming, do not check for new files. 
        if ~Streaming
            continue;
        end
    else
         % if not streaming, exit when all existing files are processed. 
        if Streaming 
            if waitLoopCounter < maxWaitLoopNum
                waitLoopCounter = waitLoopCounter + 1;
                pause(30);
            end
        else
            break;
        end
    end
        
    % check whether there are new coming images (only for streaming option)
    cur_fnames_cell = cell(nd, 1);
    for d = 1 : nd
        dataPath = dataPaths{d};
        dir_info = dir([dataPath, '*.tif']);
        fnames_d = {dir_info.name}';
        
        last_modify_time = (datenum(clock) - [dir_info.datenum]) * 24 * 60;
        latest_modify_time = min(last_modify_time);

        % only check the modify time of the latest file
        inds = cellfun(@(x) ~any(contains(fnames_d, x)), fnames(fdinds == d));
        if latest_modify_time < minModifyTime
            inds = inds & (last_modify_time' ~= latest_modify_time);
        end
        cur_fnames_cell{d} = fnames_d(inds);
    end

    cur_fdinds = arrayfun(@(x) ones(numel(cur_fnames_cell{x}), 1) * x, 1 : nd, 'unif', 0);    
    cur_fnames = cat(1, cur_fnames_cell{:});
    cur_fdinds = cat(1, cur_fdinds{:});
    if isempty(cur_fnames)
        continue;
    end

    % filter filenames by channel patterns
    include_flag = false(numel(cur_fnames), 1);
    for c = 1 : numel(ChannelPatterns)
        include_flag = include_flag | contains(cur_fnames, ChannelPatterns{c});
    end
    cur_fnames = cur_fnames(include_flag);
    cur_fdinds = cur_fdinds(include_flag);
       
    % add new files and their computing flags
    nFnew = numel(cur_fnames);
    newfnames = cur_fnames;
    fnames(end + 1 : end + nFnew) = newfnames;
    fdinds = [fdinds; cur_fdinds];
    
    % For ErodeByFTP, if the folder is empty before, reassign first time
    % point if there are new coming files.
    if Decon && ErodeByFTP
        for d = 1 : nd
            if ~isempty(maskFullpaths{d})
                continue;
            end
            c = 1;
            FTPfname = '';
            while isempty(FTPfname) && numel(cur_fnames_cell{d}) > 0
                all_inds = contains(cur_fnames_cell{d}, ChannelPatterns{c});
                FTPfname = cur_fnames_cell{d}{find(all_inds, 1, 'first')};

                c = c + 1;
            end
            if ~isempty(FTPfname)
                ind_d = find(strcmp(fnames, FTPfname));
                FTP_inds(d) = ind_d;
                maskFullpaths{d} = sprintf('%s/Masks/%s_eroded.tif', deconPaths{d}, FTPfname(1 : end - 4));
            end
        end    
    end

    nF = numel(fnames);
    is_done_flag = [is_done_flag; false(nFnew, 4)];
    trial_counter = [trial_counter; zeros(nFnew, 4)];
    
    if parseCluster
        job_ids = [job_ids;  -ones(nFnew, 4)];
    end
end

end

