function [] = XR_RLdeconFrame3D(frameFullpaths, pixelSize, dz, varargin)
% cuda Deconvolution for a single frame. It support both small file that
% can be fitted to a single GPU and also large files. For large files, it
% split files into chunks for deconvolution and the combine them together. It
% supports cluster computing for file chunks.  
%
% based on XR_matlabDeconFrame3D.m and GU_Decon.m 
% 
% Author: Xiongtao Ruan (03/15/2020)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('pixelSize', @isnumeric);
ip.addRequired('dz', @isnumeric);
ip.addOptional('deconPath', '', @(x) ischar(x) || isempty(x));
ip.addParameter('PSFfile', '', @ischar);
ip.addParameter('Overwrite', false , @islogical);
ip.addParameter('Save16bit', false , @islogical);
ip.addParameter('Rotate', false , @islogical);
ip.addParameter('Deskew', false , @islogical);
ip.addParameter('SkewAngle', -32.45 , @isnumeric);
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('EdgeSoften', 5, @isnumeric); % # ofxy px to soften
ip.addParameter('zEdgeSoften', 2, @isnumeric); % # ofxy px to soften
ip.addParameter('Crop', [], @isnumeric); % requires lower and higher values for cropping
ip.addParameter('zFlip', false, @islogical);
ip.addParameter('GenMaxZproj', [0,0,1] , @isnumeric);
ip.addParameter('ResizeImages', [] , @isnumeric);
ip.addParameter('dzPSF', 0.1 , @isnumeric); %in um
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
% ip.addParameter('DoNotAdjustResForFFT', true , @islogical); % not crop chunks for deconvolution
ip.addParameter('BlockSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('Overlap', 50, @isnumeric); % block overlap
ip.addParameter('CPUMaxMem', 500, @isnumeric); % GPU Memory in Gb
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('masterCPU', false, @islogical); % master node is a cpu node, which is just for large file deconvolution. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 5, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 2, @isnumeric);
ip.addParameter('debug', false, @islogical);

ip.parse(frameFullpaths, pixelSize, dz, varargin{:});

if ischar(frameFullpaths)
    frameFullpaths = {frameFullpaths};
end

pr = ip.Results;

PSF = pr.PSFfile;
if isempty(PSF)
    error('You should provide a PSF file for the frame...\n');
end
    
% parameters

dzPSF = pr.dzPSF;
DeconIter = pr.DeconIter;
Deskew = pr.Deskew;
SkewAngle = pr.SkewAngle;
Rotate = pr.Rotate;
Save16bit = pr.Save16bit;
Crop = pr.Crop;
zFlip = pr.zFlip;
GenMaxZproj = pr.GenMaxZproj;
ResizeImages = pr.ResizeImages;

% check if background information available, if not, estimate background
% info. Currently use 99. 
Background = pr.Background;
if isempty(Background)
    Background = 99;
end

tic
OL = pr.Overlap;
BlockSize = pr.BlockSize;

CPUMaxMem = pr.CPUMaxMem;
largeFile = pr.largeFile;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
maxTrialNum = pr.maxTrialNum;
uuid = pr.uuid;
if isempty(uuid)
    uuid = get_uuid();
end

% if master node is cpu node, msterCompute must be false, it is only for
% large file cluster computing
% if masterCPU
%     masterCompute = false;
%     largeFile = true;
%     parseCluster = true;
% end

if parseCluster 
    [status, ~] = system('sinfo');
    if status ~= 0
        warning('A slurm-based computing cluster is not exist. Set parseCluster as false.')
        parseCluster = false;
    end
end

nF = numel(frameFullpaths);
for f = 1 : nF
    frameFullpath = frameFullpaths{f};
    if ~exist(frameFullpath, 'file')
        warning('%s does not exist! Skip it!', frameFullpath);
        continue;
    end
    [pathstr, fname, ext] = fileparts(frameFullpath);
    deconPath = pr.deconPath;
    if isempty(deconPath)
        deconPath = [pathstr, filesep, 'matlab_decon'];
        mkdir(deconPath);
    end 

    % check file size
    [estMem] = XR_estimateComputingMemory(frameFullpath, {'deconvolution'}, 'cudaDecon', false);
    
    % The memory size is 500G
    if estMem > CPUMaxMem
        largeFile = true;
    end

    % first try single GPU deconvolution, if it fails split into multiple chunks
    deconFullPath = [deconPath filesep fname '_decon.tif'];
    % deconTempPath = [deconPath filesep fname '_decon.tif'];
    if exist(deconFullPath, 'file') && ~pr.Overwrite
        disp('Deconvolution results already exist, skip it!');
        continue;
    end

    if ~largeFile
        tic
        RLdecon(frameFullpath, PSF, Background, DeconIter, dzPSF, dz, Deskew, [], SkewAngle, ...
            pixelSize, Rotate, Save16bit, Crop, zFlip, GenMaxZproj, ResizeImages, [])
        toc
        if exist(deconFullPath, 'file')
            fprintf('Completed RL deconvolution of %s.\n', frameFullpath);
            continue;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Large-file deconvolution if the file is too large or single file
    % deconvolution fails.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('Start Large-file RL Decon for %s...\n', fname);

    tic
    fprintf(['reading ' fname '...\n'])
    im = readtiff(frameFullpath);
    toc
    % calculate number of chunks to break the image file
    % [my, mx, mz] = size(im);
    % [xmin,xmax,ymin,ymax,zmin,zmax,nn] = GU_extract_subVolCoordinates(mx,my,mz,csx,csy,csz,OL);
    imSize = size(im);
    dtype = class(im);
    if pr.debug
        [xmin,xmax,ymin,ymax,zmin,zmax,nn] = XR_subVolumeCoordinatesExtraction_test(imSize, 'BlockSize', BlockSize, 'overlapSize', OL);
    else
        [xmin,xmax,ymin,ymax,zmin,zmax,nn] = XR_subVolumeCoordinatesExtraction(imSize, 'BlockSize', BlockSize, 'overlapSize', OL, 'maxSubVolume', 1e11);
    end
    
    % create a folder for the file and write out the chunks
    chunkPath = [deconPath filesep fname];
    chunkDeconPath = [chunkPath filesep 'matlab_decon'];
    chunkDeconMIPPath = [chunkPath filesep 'matlab_decon' filesep 'MIPs'];
    mkdir(chunkPath);
    mkdir(chunkDeconPath);
    mkdir(chunkDeconMIPPath);
    
    if parseCluster
        if  ~exist(jobLogDir, 'dir')
            warning('The job log directory does not exist, use %s/job_logs as job log directory.', chunkPath)
            jobLogDir = sprintf('%s/job_logs', chunkPath);
            mkdir(jobLogDir);
        end
        job_log_fname = [jobLogDir, '/job_%A_%a.out'];
        job_log_error_fname = [jobLogDir, '/job_%A_%a.err'];
    end
    
    tic
    fprintf('Writing image chunks...\n')
    chunkFnames = cell(nn, 1);
    for ck = 1:nn
        chunkFnames{ck} = sprintf('chunk_%d_%d_%d_%d_%d_%d.tif', xmin(ck), ymin(ck), zmin(ck), xmax(ck), ymax(ck), zmax(ck));
        im_chunk = im(ymin(ck):ymax(ck), xmin(ck):xmax(ck), zmin(ck):zmax(ck));
        % for blank region, just skip it
        if all(im_chunk == 0, 'all')
            chunkFnames{ck} = '';
            continue;
        end
        if ~exist([chunkPath, filesep, chunkFnames{ck}], 'file') || pr.Overwrite
            writetiff(im_chunk, [chunkPath, filesep, chunkFnames{ck}(1 : end - 4), '_', uuid, '.tif']);
            movefile([chunkPath, filesep, chunkFnames{ck}(1 : end - 4), '_', uuid, '.tif'], [chunkPath, filesep, chunkFnames{ck}]);
        end        
    end
    
    clear im;
    toc

    tic
    fprintf('Deconvolving chunks...\n')

    is_done_flag = false(nn, 1);
    trial_counter = zeros(nn, 1);
    tmpFullnames = cell(nn, 1);
    
    if parseCluster
        job_ids = -ones(nn, 1);
        % job_status_flag = false(nn, 1);
    end

    while ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') % ...
            % || (parseCluster && any(job_status_flag & ~is_done_flag, 'all'))
        lastck = find(~is_done_flag & trial_counter < maxTrialNum, 1, 'last');
        for ck = 1 : nn
            if is_done_flag(ck) || trial_counter(ck) >= maxTrialNum 
                continue;
            end
            
            task_id = ck;
            
            chunkDeconFullpath = [chunkDeconPath filesep chunkFnames{ck}(1:end-4) '_decon.tif'];
            if exist(chunkDeconFullpath, 'file') || isempty(chunkFnames{ck}) % || pr.OverwriteDecon
                is_done_flag(ck) = true;
                continue;
            end
            
            % Create a softlink with uuid of the chunk and rename the file
            % after write to disk is done
            uuid = get_uuid();
            tmpChunkFullname = sprintf('%s_%s.tif', [chunkPath, filesep, chunkFnames{ck}(1:end-4)], uuid);
            softlink_cmd = sprintf('ln -s %s %s', [chunkPath, filesep, chunkFnames{ck}], tmpChunkFullname);
            matlab_cmd = sprintf(['addpath(genpath(pwd));tic;RLdecon(''%s'',''%s'',%.10f,%.10f,%.10f,%.10f,''%s'',[],', ...
                '%.10f,%.10f,''%s'',''%s'',''%s'',''%s'',[%s],[%s],[]);toc;'], ...
                tmpChunkFullname, PSF, Background, DeconIter, dzPSF, dz, string(Deskew), SkewAngle, pixelSize, ...
                string(Rotate), string(Save16bit), string(Crop), string(zFlip), num2str(GenMaxZproj, '%.10f,'), ...
                num2str(ResizeImages, '%.10f,'));
            DeconCommand = sprintf('module load matlab/r2020a; matlab -nodisplay -nosplash -nodesktop -r \\"%s\\"', matlab_cmd);
            rename_cmd = sprintf('mv %s_%s_decon.tif %s_decon.tif', [chunkDeconPath, filesep, chunkFnames{ck}(1:end-4)], uuid, [chunkDeconPath, filesep, chunkFnames{ck}(1:end-4)]);
            chunk_decon_cmd = sprintf('%s; %s; %s', softlink_cmd, DeconCommand, rename_cmd);
            if parseCluster
                job_status = check_slurm_job_status(job_ids(ck), task_id);

                % kill the first pending job and use master node do the computing.
                if job_status == 0.5 && (masterCompute && ck == lastck)
                    system(sprintf('scancel %d_%d', job_ids(ck), task_id), '-echo');
                    trial_counter(ck) = trial_counter(ck) - 1;
                end

                % if the job is still running, skip it. 
                if job_status == 1 
                    continue;
                end
                
                % for just finished job, pause a short time to make sure
                % the file is shown in the file system.
                if job_status == -1 && job_ids(ck) > 0
                    if exist(tmpFullnames{ck}, 'file')
                        pause(5);
                    end
                    pause(3);    
                    if exist(chunkDeconFullpath, 'file')
                        is_done_flag(ck) = true;
                        continue;
                    end
                end

                % If there is no job, submit a job
                if job_status == -1 && ~(masterCompute && ck == lastck)
                    cmd = sprintf('sbatch --array=%d -o %s -e %s -p abc --qos abc_normal -n1 --mem=500G --cpus-per-task=24 --wrap="%s"', ...
                        task_id, job_log_fname, job_log_error_fname, chunk_decon_cmd);

                    [status, cmdout] = system(cmd, '-echo');
                    tmpFullnames{ck} =  sprintf('%s/%s_%s_decon.tif', chunkDeconPath, chunkFnames{ck}(1:end-4), uuid);
                    
                    trial_counter(ck) = trial_counter(ck) + 1;    
    
                    job_id = regexp(cmdout, 'Submitted batch job (\d+)\n', 'tokens');
                    job_id = str2double(job_id{1}{1});
                    job_ids(ck) = job_id;
                end
            end

            % for non cluster, cluster but pending job or the last chunk, use
            % the master job for the computing. 
            if ~parseCluster || (parseCluster && masterCompute && ck == lastck)
                % [status, cmdout] = system(chunk_decon_cmd, '-echo');
                tic
                if ~ispc
                    system(softlink_cmd);
                    RLdecon(tmpChunkFullname, PSF, Background, DeconIter, dzPSF, dz, Deskew, [], SkewAngle, ...
                        pixelSize, Rotate, Save16bit, Crop, zFlip, GenMaxZproj, ResizeImages, []);
                    system(rename_cmd);
                else
                    RLdecon([chunkPath, filesep, chunkFnames{ck}], PSF, Background, DeconIter, dzPSF, dz, Deskew, [], SkewAngle, ...
                        pixelSize, Rotate, Save16bit, Crop, zFlip, GenMaxZproj, ResizeImages, []);
                end
                toc
                
                tmpFullnames{ck} =  sprintf('%s/%s_%s_decon.tif', chunkDeconPath, chunkFnames{ck}(1:end-4), uuid);

                trial_counter(ck) = trial_counter(ck) + 1;    
            end

            if exist(chunkDeconFullpath, 'file') % || pr.OverwriteDecon
                is_done_flag(ck) = true;
                continue;
            end
        end

        % wait 5 seconds if some tasks are still computing
        if ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') 
            pause(30);
        end    
    end

    % combine the deconvovled segments and write file
    tic
    fprintf('Combining deconvovled chunks...\n')
    im = zeros(imSize, dtype);
    lol = floor(OL / 2);
    rol = ceil(OL / 2); 
    for ck = 1:nn
        chunkDeconFullpath = [chunkDeconPath filesep chunkFnames{ck}(1:end-4) '_decon.tif'];
        if exist(chunkDeconFullpath, 'file')
            tim = readtiff(chunkDeconFullpath);
            [tsy, tsx, tsz] = size(tim);

            yrange = ymin(ck) + (ymin(ck) ~= 1) * lol : ymax(ck) - (ymax(ck) ~= imSize(1)) * rol;
            xrange = xmin(ck) + (xmin(ck) ~= 1) * lol : xmax(ck) - (xmax(ck) ~= imSize(2)) * rol;
            zrange = zmin(ck) + (zmin(ck) ~= 1) * lol : zmax(ck) - (zmax(ck) ~= imSize(3)) * rol;
            yrange = yrange(yrange - ymin(ck) + 1 <= tsy);
            xrange = xrange(xrange - xmin(ck) + 1 <= tsx);
            zrange = zrange(zrange - zmin(ck) + 1 <= tsz);
            tyrange = yrange - ymin(ck) + 1;
            txrange = xrange - xmin(ck) + 1;
            tzrange = zrange - zmin(ck) + 1;
            
            im(yrange, xrange, zrange) = tim(tyrange, txrange, tzrange);
            % im((ymin(ck)+OL/2):(ymin(ck)+tsy-1-OL/2), (xmin(ck)+OL/2):(xmin(ck)+tsx-1-OL/2), (zmin(ck)+OL/2):(zmin(ck)+tsz-1-OL/2)) = tim((OL/2+1):(end-OL/2),(OL/2+1):(end-OL/2),(OL/2+1):(end-OL/2));
        elseif ~isempty(chunkFnames{ck})
            error('Deconvolution of chunk %s is missing!', chunkFnames{ck});
        end
    end
    toc

    tic
    fprintf('Saving combined deconvolved file...\n')
    tmp_xy = max(im,[],3);

    if pr.Save16bit
        im = uint16(im);
    else
        im = single(im);
    end

    deconTmpPath = sprintf('%s_%s.tif', deconFullPath(1 : end - 4), uuid);
    writetiff(im, deconTmpPath);
    movefile(deconTmpPath, deconFullPath);
    
    deconMIPPath = sprintf('%s/MIPs/', deconPath);
    mkdir(deconMIPPath);
    deconMIPFullPath = sprintf('%s%s.tif', deconMIPPath, fname);
    writetiff(uint16(tmp_xy), deconMIPFullPath);
    toc

    % delete temporary files generated during the deconvolution
    rmdir(chunkPath, 's');

    fprintf('Completed RL deconvolution of %s.\n', frameFullpath);
end

toc

end
