function [] = simRecon_data_wrapperBAK(dataPaths, varargin)
% Adapted from XR_microscopeAutomaticProcessing


ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('dataPaths'); % data structure from loadConditionData
ip.addParameter('psfFullpaths', {'','',''}, @iscell);
ip.addParameter('Overwrite', false,  @(x) (numel(x) == 1 || numel(x) == 5) && islogical(x));
ip.addParameter('Streaming', true,  @islogical); % if true, check for new files. If false, assume all files transferred completely.
ip.addParameter('ChannelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('Channels', [488, 560, 642], @isnumeric);
ip.addParameter('SkewAngle', 32.45, @isscalar);
ip.addParameter('dz', 0.5, @isscalar);
ip.addParameter('xyPixelSize', 0.108, @isscalar);
ip.addParameter('Reverse', true, @islogical);
ip.addParameter('ObjectiveScan', false, @islogical);
ip.addParameter('ZstageScan', false, @islogical);
ip.addParameter('sCMOSCameraFlip', false, @islogical);
ip.addParameter('Save16bit', [false, false, false, false], @(x) (numel(x) == 1 || numel(x) == 4) && islogical(x));
ip.addParameter('onlyFirstTP', false, @islogical);
ip.addParameter('dzFromEncoder', false, @islogical);

% pipeline parameters
ip.addParameter('Deskew', true, @islogical);
ip.addParameter('DeskewPSFs', false, @islogical);
ip.addParameter('Recon', true, @islogical);

% recon parameters
ip.addParameter('islattice', true, @islogical); %Flag to indicate if this is light sheet data
ip.addParameter('NA_det', 1.0, @isnumeric);
ip.addParameter('NA_ext', 0.55, @isnumeric);
ip.addParameter('nimm', 1.33, @isnumeric);
ip.addParameter('wvl_em', .605, @isnumeric);
ip.addParameter('wvl_ext', .560, @isnumeric);
ip.addParameter('w', 5e-3, @isnumeric); %Wiener coefficient for regularization
ip.addParameter('apodize', true, @islogical); %Flag to indicate whether or not to apodize the final data

ip.addParameter('DS', false, @islogical);
ip.addParameter('nphases', 5, @isnumeric);
ip.addParameter('norders', 5, @isnumeric);
ip.addParameter('norientations', 1, @isnumeric);
ip.addParameter('lattice_period', 1.2021, @isnumeric); %Lattice period in microns - this is the coarsest period
ip.addParameter('lattice_angle', [pi/2], @isnumeric); %Angle parellel to pattern modulation (assuming horizontal is zero)
ip.addParameter('phase_step', .232, @isnumeric); %Phase step in microns
ip.addParameter('pxl_dim_data', [0.11,0.11,0.3*sind(32.4)], @isnumeric); %Voxel dimensions of the image in microns - note, stored as [dy,dx,dz]
ip.addParameter('pxl_dim_PSF', [0.11,0.11,0.2*sind(32.4)], @isnumeric); %Voxel dimensions of the PSF in microns - note, stored as [dy,dx,dz]
ip.addParameter('Background', 105, @isnumeric);

ip.addParameter('normalize_orientations', false, @islogical); %Flag to indicate whether or not to normalize total intensity for each orientation
ip.addParameter('perdecomp', false, @islogical); %Flag to indicate whether or not to use periodic/smooth decomposition to reduce edge effects
ip.addParameter('edgeTaper', false, @islogical); %Flag to indicate whether or not to window the data to reduce edge effects
ip.addParameter('edgeTaperVal', 0.1, @isnumeric); %Roll-off parameter for Tukey windowing

ip.addParameter('useGPU', true, @islogical);

ip.addParameter('flipZstack', false, @islogical);
ip.addParameter('EdgeSoften', 5, @isnumeric); % # ofxy px to soften
ip.addParameter('zEdgeSoften', 2, @isnumeric); % # ofxy px to soften
ip.addParameter('Crop', [], @isnumeric); % requires lower and higher values for cropping
ip.addParameter('zFlip', false, @islogical);
ip.addParameter('GenMaxZproj', [0,0,1] , @isnumeric);
ip.addParameter('ResizeImages', [] , @isnumeric);
ip.addParameter('EdgeErosion', 0 , @isnumeric); % erode edges for certain size.
ip.addParameter('ErodeMaskfile', '', @ischar); % erode edges file
ip.addParameter('SaveMaskfile', false, @islogical); % save mask file for common eroded mask
ip.addParameter('ChunkSize', [250, 250, 250] , @isvector); % in y, x, z
ip.addParameter('Overlap', 10, @isnumeric); % block overlap
ip.addParameter('maxSubVolume', 3e8, @isnumeric);
ip.addParameter('CPUMaxMem', 500, @isnumeric); % CPU Memory in Gb

ip.addParameter('GPUJob', false, @islogical); % use gpu for chunk reconstruction. 
% job related parameters
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
ip.addParameter('MatlabLaunchStr', 'module load matlab/r2020b; matlab -nodisplay -nosplash -nodesktop -nojvm -r', @ischar);
ip.addParameter('SlurmParam', '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M', @ischar);

ip.parse(dataPaths, varargin{:});


pr = ip.Results;
psfFullpaths = pr.psfFullpaths;
Overwrite = pr.Overwrite;
Streaming = pr.Streaming;
% Resolution = pr.Resolution;
SkewAngle = pr.SkewAngle;
dz = pr.dz;
xyPixelSize = pr.xyPixelSize;
ObjectiveScan = pr.ObjectiveScan;
ZstageScan = pr.ZstageScan;
Reverse = pr.Reverse;
ChannelPatterns = pr.ChannelPatterns;
Save16bit = pr.Save16bit;
dzFromEncoder = pr.dzFromEncoder;

% pipeline parameters
Deskew = pr.Deskew;
DeskewPSFs = pr.DeskewPSFs;
Recon = pr.Recon;

% recon parameters
islattice = pr.islattice;
NA_det = pr.NA_det;
NA_ext = pr.NA_ext;
nimm = pr.nimm;
wvl_em = pr.wvl_em;
wvl_ext = pr.wvl_ext;
w = pr.w;
apodize = pr.apodize;

DS = pr.DS;
nphases = pr.nphases;
norders = pr.norders;
norientations = pr.norientations;
lattice_period = pr.lattice_period;
lattice_angle = pr.lattice_angle;
phase_step = pr.phase_step;
pxl_dim_data = pr.pxl_dim_data;
pxl_dim_PSF = pr.pxl_dim_PSF;
Background = pr.Background;

normalize_orientations = pr.normalize_orientations;
perdecomp = pr.perdecomp;
edgeTaper = pr.edgeTaper;
edgeTaperVal = pr.edgeTaperVal;

GPUJob = pr.GPUJob;

%frame job related
flipZstack = pr.flipZstack;
EdgeSoften = pr.EdgeSoften;
zEdgeSoften = pr.zEdgeSoften;
Crop = pr.Crop;
zFlip = pr.zFlip;
GenMaxZproj = pr.GenMaxZproj;
ResizeImages = pr.ResizeImages;
EdgeErosion = pr.EdgeErosion;
ErodeMaskfile = pr.ErodeMaskfile;
SaveMaskfile = pr.SaveMaskfile;
ChunkSize = pr.ChunkSize;
Overlap = pr.Overlap;
maxSubVolume = pr.maxSubVolume;
CPUMaxMem = pr.CPUMaxMem;

% job related
largeFile = pr.largeFile;
%zarrFile = pr.zarrFile;
jobLogDir = pr.jobLogDir;
parseCluster = pr.parseCluster;
cpusPerTask = pr.cpusPerTask;
cpuOnlyNodes = pr.cpuOnlyNodes;
uuid = pr.uuid;
maxTrialNum = pr.maxTrialNum;
unitWaitTime = pr.unitWaitTime;
MatlabLaunchStr = pr.MatlabLaunchStr;
SlurmParam = pr.SlurmParam;

% make sure the function is in the root of XR_Repository. 
mpath = fileparts(which(mfilename));
repo_rt = [mpath, '/../../../'];
cd(repo_rt);
%cd(/clusterfs/fiona/ABCcode/XR_Repository);

% suppress directory exists warning
warning('off', 'MATLAB:MKDIR:DirectoryExists');

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
for d = 1 : nd
    dataPath = dataPaths{d};
    if ~strcmp(dataPath(end), '/')
        dataPaths{d} = [dataPath, '/'];
    end
end

if ispc
    dataPaths = cellfun(@(x) strrep(x, '\', '/'), dataPaths, 'unif', 0);
    psfFullpaths = cellfun(@(x) strrep(x, '\', '/'), psfFullpaths, 'unif', 0);    
end

% check if a slurm-based computing cluster exists
if parseCluster
    [parseCluster, job_log_fname, job_log_error_fname, slurm_constraint_str, jobLogDir] = checkSlurmCluster(dataPath, jobLogDir, cpuOnlyNodes);
end

if Deskew
        dsPaths = cell(nd, 1);
        for d = 1 : nd
            dataPath = dataPaths{d};
            dsPath = [dataPath, '/DS/'];
            if Overwrite(1) && exist(dsPath, 'dir')
                rmdir(dsPath, 's');
            end
            if ~exist(dsPath, 'dir')
                mkdir(dsPath);
                fileattrib(dsPath, '+w', 'g');
            end
            dsPaths{d} = dsPath;
        end
end

% get actual dz from encoder positions
if dzFromEncoder && ~Streaming
    dz_all = zeros(nd, 1);
    for d = 1 : nd
        dataPath = dataPaths{d};        
        dz_actual = XR_estimate_actual_step_size_from_encoder(dataPath, 'dz', dz);
        dz_all(d) = dz_actual;
    end   
else
    dz_all = repmat(dz, nd, 1);
end

[fnames, fdinds, gfnames, partialvols, dataSizes, flipZstack_mat, FTP_inds, maskFullpaths] = ...
    XR_parseImageFilenames(dataPaths, ChannelPatterns, false, flipZstack, false, {}, false);

nF = numel(fnames);

% flags: for thee: deskew, recon
is_done_flag = false(nF, 3);
% ensure only one stitch wrapper is running of a dataset
stitch_running = false(nd, 1); 
trial_counter = zeros(nF, 3);
waitLoopCounter = 0;


if parseCluster
    job_ids = -ones(nF, 4);
    imSize_mat = zeros(nF, 3, 2);
    dataSize_mat = zeros(nF, 2);
    dataSize_mat(:, 1) = dataSizes;
end

matlab_setup_str = 'setup([],true)';

% gen otfs
if DeskewPSFs
    
else
    for i = 1:numel(psfFullpaths)
        cPSF = readtiff(psfFullpaths{i});
        cOTF = sim_PSFtoOTF_gen(cPSF,'useGPU',false);
    end
end

% use while loop to perform computing for all images
while ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') || ...
        (Streaming && (nF == 0 || latest_modify_time < maxModifyTime || waitLoopCounter < maxWaitLoopNum))
    for f = 1 : nF
        tic
        if all(is_done_flag(f, :))
            continue;
        end
        
        % first deskew
        fname = fnames{f};
        [~, fsname] = fileparts(fname);
        fdind = fdinds(f);
        partialvol = partialvols(f);
        gfname = gfnames{f};
        dataPath = dataPaths{fdind};
        dz_f = dz_all(fdind);
        
        frameFullpath = [dataPath, fname];
        % check wheter the file is deleted during the computing.
        if ~exist(frameFullpath, 'file')
            is_done_flag(f, :) = true;
            continue
        end
        
        task_id = rem(f, 5000);
        FTP_ind = FTP_inds(fdind);
        
        %% deskew
        if Deskew
            if true
                dsPath = dsPaths{fdind};
                dsFullpath = [dsPath, fname];
                tmpFullpath = sprintf('%s.tmp', dsFullpath(1 : end - 4));
            end

            if (exist(dsFullpath, 'file'))
                is_done_flag(f, 1) = true;
                if exist(tmpFullpath, 'file')
                    delete(tmpFullpath);
                end
            end
        else
            is_done_flag(f, 1) = true;
        end
        
        if ~is_done_flag(f, 1) 
            
            flipZstack = flipZstack_mat(f);
            
            % set up input file for either single volume file or a
            % group of files
            ds_input_path = {frameFullpath};
            if partialvol
                ds_input_path = cellfun(@(x) [dataPath, x], gfname, 'unif', 0);
            end
            ds_input_str = sprintf('{''%s''}', strjoin(ds_input_path, ''','''));
            
                    func_str = sprintf(['deskewPhasesFrame(%s,%.20d,%.20d,''SkewAngle'',%.20d,', ...
            '''Reverse'',%s,''nphases'',%.10d)'], ...
            ds_input_str, xyPixelSize, dz_f, SkewAngle, string(Reverse), nphases);
            
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
                        
                        %if ~true
                         %   estRequiredMemory = XR_estimateComputingMemory(frameFullpath, 'steps', {'deskew'});
                        %else
                            estRequiredMemory = dataSize_mat(f, 1) / 2^30 * 2 * (6 + 4 / prod([]) + trial_counter(f, 1) * 2);
                        %end
                        cpusPerTask_ds = cpusPerTask;
                        if cpusPerTask_ds * 20 < estRequiredMemory
                            cpusPerTask_ds = min(24, ceil(estRequiredMemory / 20));
                        else
                            cpusPerTask_ds = min(cpusPerTask_ds, ceil(estRequiredMemory / 20));
                        end
                            
                        matlab_cmd = sprintf('%s;tic;%s;toc', matlab_setup_str, func_str);
                        process_cmd = sprintf('%s \\"%s\\"', MatlabLaunchStr, matlab_cmd);
                        cmd = sprintf(['sbatch --array=%d -o %s -e %s --cpus-per-task=%d %s %s ', ...
                            '--wrap="echo Matlab command:  \\\"%s\\\"; %s"'], task_id, job_log_fname, ...
                            job_log_error_fname, cpusPerTask_ds, SlurmParam, slurm_constraint_str, ...
                            matlab_cmd, process_cmd);
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
                        fclose(fopen(tmpFullpath, 'w'));
                    end
                end
            else
                fclose(fopen(tmpFullpath, 'w'));
            end
            if ~parseCluster
                tic; feval(str2func(['@()', func_str])); toc;
                trial_counter(f, 1) = trial_counter(f, 1) + 1;
            end

            % check if computing is done
            if (exist(dsFullpath, 'file'))
                is_done_flag(f, 1) = true;
                if exist(tmpFullpath, 'file')
                    delete(tmpFullpath);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% then deconvolution
        if Recon
            % if the deskew result is not available, wait for the deskew data
            % to be finished

            % input chosen order is dsr, ds, raw (in decreased priority order)
            dcframeFullpath = frameFullpath;
            
            if Deskew
                dcframeFullpath = dsFullpath;
                dc_dz = dz_f;
                dc_dzPSF = dzPSF;
                dc_psfFullpaths = psfFullpaths;
            end
           
            
            if ~exist(dcframeFullpath, 'file')
                % if DS, it means the deskew and rotation not done.
                if Deskew
                    continue;
                end
            end
                            
            deconPath = deconPaths{fdind};
            deconFullpath = sprintf('%s/%s_decon.tif', deconPath, fsname);
            dctmpFullpath = sprintf('%s.tmp', deconFullpath(1 : end - 4));

            if exist(deconFullpath, 'file')
                is_done_flag(f, 2) = true;
                if exist(dctmpFullpath, 'file')
                    delete(dctmpFullpath);
                end
            end
            
            if parseCluster
                if dataSize_mat(f, 2) == 0 && (~all(is_done_flag(f, 2)) || f == FTP_ind)
                    dir_info = dir(dcframeFullpath);
                    dataSize_mat(f, 1) = dir_info.bytes;
                end
            end
            
            % for ErodeByFTP, check if the mask file exist
            maskFullpath = '';
            SaveMaskfile = false;
            %{
            if ErodeByFTP
                if f == FTP_ind
                    SaveMaskfile = true;
                    % if decon result exist, but mask file not exist, rerun
                    % it to save the mask. 
                    if Stitch
                        maskFullpaths{fdind} = sprintf('%s/Masks/%s_eroded.tif', deconPath, stch_fsname);
                    end
                    if is_done_flag(f, 3) && ~exist(maskFullpaths{fdind}, 'file')
                        is_done_flag(f, 3) = false;
                        fprintf('Mask file %s does not exist, delete the deconvolved result\n', maskFullpaths{fdind});
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
            %}
        else
            is_done_flag(f, 2) = true;    
        end

        if ~is_done_flag(f, 2) 
            % psfMapping =  ~cellfun(@isempty, regexpi(fname, ChannelPatterns));
            % change to contains.m to unify the matching
            psfMapping =  ~cellfun(@isempty, regexpi(frameFullpath, ChannelPatterns));
            
            psfFullpath = dc_psfFullpaths{psfMapping};
            
            % do not use rotation in decon functions
            if cudaDecon
                func_str = sprintf(['XR_cudaDeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                    '''cudaDeconPath'',''%s'',''OTFGENPath'',''%s'',''dzPSF'',%.10f,''Background'',[%d],', ...
                    '''SkewAngle'',%d,''Rotate'',%s,''DeconIter'',%d,''Save16bit'',%s,''largeFile'',%s)'], ...
                    dcframeFullpath, xyPixelSize, dc_dz, psfFullpath, cudaDeconPath, OTFGENPath, dc_dzPSF, ...
                    Background, SkewAngle, string(deconRotate), DeconIter, string(Save16bit(3)), string(largeFile));
            elseif cppDecon
                func_str = sprintf(['XR_cppDeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                    '''cppDeconPath'',''%s'',''loadModules'',''%s'',''dzPSF'',%.10f,''Background'',[%d],', ...
                    '''SkewAngle'',%d,''EdgeErosion'',%d,''ErodeMaskfile'',''%s'',''SaveMaskfile'',%s,', ...
                    '''Rotate'',%s,''DeconIter'',%d,''Save16bit'',%s,''largeFile'',%s)'], dcframeFullpath, ...
                    xyPixelSize, dc_dz, psfFullpath, cppDeconPath, loadModules, dc_dzPSF, Background, ...
                    SkewAngle, EdgeErosion, maskFullpath, string(SaveMaskfile), string(deconRotate), ...
                    DeconIter, string(Save16bit(3)), string(largeFile));
            else
                func_str = sprintf(['XR_RLdeconFrame3D(''%s'',%.10f,%.10f,'''',''PSFfile'',''%s'',', ...
                    '''dzPSF'',%.10f,''Background'',[%d],''SkewAngle'',%d,''EdgeErosion'',%d,''ErodeMaskfile'',''%s'',', ...
                    '''SaveMaskfile'',%s,''Rotate'',%s,''DeconIter'',%d,''RLMethod'',''%s'',''fixIter'',%s,', ...
                    '''errThresh'',[%0.20f],''debug'',%s,''GPUJob'',%s,''Save16bit'',%s,''largeFile'',%s)'], ...
                    dcframeFullpath, xyPixelSize, dc_dz, psfFullpath,  dc_dzPSF, Background, SkewAngle, ...
                    EdgeErosion, maskFullpath, string(SaveMaskfile), string(deconRotate), DeconIter, RLMethod, ...
                    string(fixIter), errThresh, string(debug), string(GPUJob), string(Save16bit(3)), string(largeFile));
            end
           
            if exist(dctmpFullpath, 'file') || parseCluster
                if parseCluster
                    job_status = check_slurm_job_status(job_ids(f, 3), task_id);

                     % if the job is still running, skip it. 
                    if job_status == 1 
                        continue;
                    end

                    if job_status == -1
                        % for matlab decon,  decide how many cores. 
                        if ~cudaDecon
                            [estMem, estGPUMem] = XR_estimateComputingMemory('', {'deconvolution'}, ...
                                'dataSize', dataSize_mat(f, 2), 'cudaDecon', false);
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
                tic; feval(str2func(['@()', func_str])); toc;
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
            func_str = sprintf(['XR_RotateFrame3D(''%s'',%.20d,%.20d,''SkewAngle'',%.20d,', ...
                '''ObjectiveScan'',%s,''Reverse'',%s,''Save16bit'',%s)'], deconFullpath, ...
                xyPixelSize, dz_f, SkewAngle, string(ObjectiveScan), string(Reverse), string(Save16bit(4)));

            if exist(rdctmpFullpath, 'file') || parseCluster
                if parseCluster
                    job_status = check_slurm_job_status(job_ids(f, 4), task_id);
                    
                    % if the job is still running, skip it. 
                    if job_status == 1 
                        continue;
                    end
                    
                    if job_status == -1
                        [estMem, estGPUMem] = XR_estimateComputingMemory('', {'deconvolution'}, ...
                            'dataSize', dataSize_mat(f, 2), 'cudaDecon', false);
                        cpusPerTask_dcr = cpusPerTask;
                        if cpusPerTask_dcr * 20 < estMem
                            cpusPerTask_dcr = min(24, ceil(estMem / 20));
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
                tic; feval(str2func(['@()', func_str])); toc;
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
        toc
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

    % [fnames, fdinds, gfnames, partialvols, dataSizes, flipZstack_mat, FTP_inds, maskFullpaths] = ...
    %     XR_parseImageFilenames(dataPaths, ChannelPatterns, flipZstack, Decon, Streaming);

    nF = numel(fnames);
    is_done_flag = [is_done_flag; false(nFnew, 4)];
    trial_counter = [trial_counter; zeros(nFnew, 4)];
    

    if parseCluster
        job_ids = [job_ids;  -ones(nFnew, 4)];
        imSize_mat = cat(1, imSize_mat, zeros(nFnew, 3, 2));
    end
end


end