function [] = XR_cudaDeconFrame3D(frameFullpaths, pixelSize, dz, varargin)
% cuda Deconvolution for a single frame. It support both small file that
% can be fitted to a single GPU and also large files. For large files, it
% split files into chunks for deconvolution and the combine them together. It
% supports cluster computing for file chunks.  
%
% It uses Lin Shao's GPU deconvolution function, and is based on GU_cudaDecon.m and
% GU_cudaDecon_multiGPU_LargeFiles.m
%
% 
% Author: Xiongtao Ruan (02/27/2020)
% xruan (03/06/2020): add support for a batch of files with same parameters
% and settings (reduce overhead of deconvolution for small files).
% xruan (03/11/2020): skip deconvolution of empty chunks.
% xruan (03/12/2020): add support for allow a cpu node as master node for
% large file deconvolution, which will save time for write image chunks and
% combine deconvolved chunks. 
% xruan (03/19/2020): use softlink for small files. 
% xruan (03/20/2020): throw error for small file if it fails, rather than go to split-based deconvolution.  
% xruan (08/21/2020): add psf2otf.m as backup function for otfgen


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('pixelSize', @isnumeric);
ip.addRequired('dz', @isnumeric);
ip.addOptional('deconPath', '', @(x) ischar(x) || isempty(x));
ip.addParameter('PSFfile', '' , @ischar);
ip.addParameter('cudaDeconPath', '/global/home/groups/software/sl-7.x86_64/modules/cudaDecon/bin/cudaDeconv' , @ischar);
ip.addParameter('OTFGENPath', '/global/home/groups/software/sl-7.x86_64/modules/cudaDecon/bin/radialft' , @ischar); % point to radialft file
ip.addParameter('Overwrite', false , @islogical);
ip.addParameter('OverwriteOTFs', ~true , @islogical);
ip.addParameter('OTFFileName', '', @isstr);
% ip.addParameter('Ch_search', '');
ip.addParameter('Verbose', true , @islogical);
ip.addParameter('SaveDeskew', false , @islogical);
ip.addParameter('Save16bit', false , @islogical);
ip.addParameter('lzwcompress', true, @islogical);
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
ip.addParameter('BlockSize', [1024,1024,1024] , @isvector); % in y, x, z
ip.addParameter('Overlap', 50, @isnumeric); % block overlap
ip.addParameter('GPUMaxMem', 12, @isnumeric); % GPU Memory in Gb
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

p = ip.Results;
masterCPU = p.masterCPU;
GPUMaxMem = p.GPUMaxMem;
gdc = gpuDeviceCount;
if gdc<1 && ~masterCPU
    error('To use cudaDecon, at least one GPU Device is needed!');
end

PSF = p.PSFfile;
if isempty(PSF)
    error('You should provide a PSF file for the frame...\n');
end

% write OTF file    
[OTFpathstr, psfFsname, ~] = fileparts(PSF);
OTFFileName = p.OTFFileName;
if isempty(OTFFileName)
    OTFFileName = sprintf('OTF_%s.tif', psfFsname);    
end
if p.OverwriteOTFs || ~exist([OTFpathstr '/' OTFFileName], 'file')
    otfCommand = [p.OTFGENPath ' ' PSF ' ' OTFpathstr '/' OTFFileName '  --nocleanup --fixorigin 10'];
    [stat, res] = system(otfCommand,'-echo');
    
    if ~exist([OTFpathstr '/' OTFFileName], 'file')
        psf = readtiff(PSF);
        otf = psf2otf(psf);
        writetiff(otf, [OTFpathstr '/' OTFFileName]);
    end
end

% parameters
Deskew_str = '';
if p.Deskew
    Deskew_str = sprintf('-D %.f', p.SkewAngle);
end

Rotate_str = '';
if p.Rotate
    Rotate_str = sprintf('-R %.f', p.SkewAngle);
end

Crop_str = '';
if ~isempty(p.Crop)
    Crop_str = sprintf('-C %s', num2str(p.Crop));
end

zFlip_str = '';
if p.zFlip
    zFlip_str = '-r';
end

SaveDeskew_str = '';
if p.SaveDeskew
    SaveDeskew_str = '-S';
end

EdgeSoften_str = '';
if ~isempty(p.EdgeSoften)
    EdgeSoften_str = sprintf('-e %d', p.EdgeSoften);
end

zEdgeSoften_str = '';
if ~isempty(p.zEdgeSoften)
    zEdgeSoften_str = sprintf('-E %d', p.zEdgeSoften);
end

GenMaxZproj_str = '';
if ~isempty(p.GenMaxZproj)
    GenMaxZproj_str = sprintf('-M %s', num2str(p.GenMaxZproj));
end

Save16bit_str = '';
if p.Save16bit
    Save16bit_str = '-u';
end

lzwcompress_str = '';
if p.lzwcompress
    lzwcompress_str = '--lzw';
end

Overwrite_str = '';
if p.Overwrite
    Overwrite_str = '--no_overwrite';
end

dzPSF = p.dzPSF;
DeconIter = p.DeconIter;

% check if background information available, if not, estimate background
% info. Currently use 99. 
Background = p.Background;
if isempty(Background)
%     frame = readtiff(frameFullpath);
%     frame(frame == 0) = nan;
%     [backgroundInfo] = XR_estimate_overall_background_information(frame); 
%     clear frame;
%     Background = backgroundInfo(1);
    Background = 99;
end

tic  
params = sprintf(' --drdata %.10f -Z %.10f -z %.10f -b %.10f -i %d %s %s %s %s %s %s %s %s %s %s %s', ...
    pixelSize, dzPSF, dz, Background, DeconIter, Deskew_str, Rotate_str, Crop_str, zFlip_str, SaveDeskew_str, ...
    EdgeSoften_str, GenMaxZproj_str, zEdgeSoften_str, Save16bit_str, lzwcompress_str, Overwrite_str);

OL = p.Overlap;
BlockSize = p.BlockSize;

largeFile = p.largeFile;
parseCluster = p.parseCluster;
jobLogDir = p.jobLogDir;
masterCompute = p.masterCompute;
maxTrialNum = p.maxTrialNum;
uuid = p.uuid;
if isempty(uuid)
    uuid = get_uuid();
end

% if master node is cpu node, msterCompute must be false, it is only for
% large file cluster computing
if masterCPU
    masterCompute = false;
    largeFile = true;
    parseCluster = true;
end

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
    [pathstr, fsname, ext] = fileparts(frameFullpath);
    fname = [fsname ext];
    deconPath = p.deconPath;
    if isempty(deconPath)
        deconPath = [pathstr, '/', 'GPUdecon'];
    end 

    % check file size
    [estMem, estGPUMem] = XR_estimateComputingMemory(frameFullpath, {'deconvolution'}, 'cudaDecon', true);
    
    % The gpu size is 12Gb
    if estGPUMem > GPUMaxMem
        largeFile = true;
    end

    % first try single GPU deconvolution, if it fails split into multiple chunks
    deconFullPath = [deconPath '/' fsname '_decon.tif'];
    % deconTempPath = [deconPath '/' fname '_decon.tif'];
    if exist(deconFullPath, 'file') && ~p.Overwrite
        disp('Deconvolution results already exist, skip it!');
        continue;
    end

    if ~largeFile
        if ispc
            DeconCommand = [p.cudaDeconPath params '--input-dir ' pathstr '/'  ' --filename-pattern ' fsname ' --otf-file ' OTFpathstr '/' OTFFileName ];
        else
            softlink_cmd = sprintf('ln -s %s %s_%s.tif', frameFullpath, frameFullpath(1:end-4), uuid);
            decon_cmd = [p.cudaDeconPath params pathstr '/' ' ' fsname '_' uuid ' ' OTFpathstr '/' OTFFileName ];
            rename_cmd = sprintf('mv %s/%s_%s_decon.tif %s/%s_decon.tif', deconPath, fsname, uuid, deconPath, fsname);
            rename_mip_cmd = sprintf('mv %s/MIPs/%s_%s_MIP_z.tif %s/MIPs/%s_MIP_z.tif', deconPath, fsname, uuid, deconPath, fsname);
            unlink_cmd = sprintf('rm %s_%s.tif', frameFullpath(1:end-4), uuid);
            DeconCommand = sprintf('%s; %s; %s; %s; %s', softlink_cmd, decon_cmd, rename_cmd, rename_mip_cmd, unlink_cmd);
        end

        if p.Verbose
            [~, ~] = system(DeconCommand, '-echo');
        else
            [~, ~] = system(DeconCommand);
        end

        if exist(deconFullPath, 'file')
            fprintf('Completed GPU deconvolution of %s.\n', frameFullpath);
            continue;
        else
            % For small file size, if it failed, just throw the error. 
            if estGPUMem < 0.1 * GPUMaxMem
                error('Deconvolution is failed with file size %d', estMem);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Large-file deconvolution if the file is too large or single file
    % deconvolution fails.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('Start Large-file Cuda Decon for %s...\n', fname);

    tic
    fprintf(['reading ' fname '...\n'])
    im = readtiff(frameFullpath);
    toc
    % calculate number of chunks to break the image file
    % [my, mx, mz] = size(im);
    % [xmin,xmax,ymin,ymax,zmin,zmax,nn] = GU_extract_subVolCoordinates(mx,my,mz,csx,csy,csz,OL);
    imSize = size(im);
    if p.Save16bit
        dtype = 'uint16';
    else
        dtype = 'single';
    end    
    if p.debug
        [xmin,xmax,ymin,ymax,zmin,zmax,nn] = XR_subVolumeCoordinatesExtraction_test(imSize, 'BlockSize', BlockSize, 'overlapSize', OL);
    else
        [xmin,xmax,ymin,ymax,zmin,zmax,nn] = XR_subVolumeCoordinatesExtraction(imSize, 'BlockSize', BlockSize, 'overlapSize', OL);
    end
    
    % create a folder for the file and write out the chunks
    chunkPath = [deconPath '/' fsname];
    chunkDeconPath = [chunkPath '/' 'GPUdecon'];
    chunkDeconMIPPath = [chunkPath '/' 'GPUdecon' '/' 'MIPs'];
    mkdir(deconPath);
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
        if ~exist(chunkFnames{ck}, 'file') || p.Overwrite
            writetiff(im_chunk, [chunkPath, '/', chunkFnames{ck}]);
        end
    end
    
    % clear im;
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
            
            chunkDeconFullpath = [chunkDeconPath '/' chunkFnames{ck}(1:end-4) '_decon.tif'];
            if exist(chunkDeconFullpath, 'file') || isempty(chunkFnames{ck}) % || p.OverwriteDecon
                is_done_flag(ck) = true;
                continue;
            end
            
            % Create a softlink with uuid of the chunk and rename the file
            % after write to disk is done
            uuid = get_uuid();
            softlink_cmd = sprintf('ln -s %s %s_%s.tif', [chunkPath, '/', chunkFnames{ck}], [chunkPath, '/', chunkFnames{ck}(1:end-4)], uuid);
            DeconCommand = [p.cudaDeconPath params '--input-dir ' chunkPath ' --filename-pattern ' chunkFnames{ck}(1:end-4) '_' uuid  ' --otf-file ' OTFpathstr '/' OTFFileName];
            rename_cmd = sprintf('mv %s_%s_decon.tif %s_decon.tif', [chunkDeconPath, '/', chunkFnames{ck}(1:end-4)], uuid, [chunkDeconPath, '/', chunkFnames{ck}(1:end-4)]);
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
                    cmd = sprintf('sbatch --array=%d -o %s -e %s -p abc --gres=gpu:1 --qos abc_normal -n1 --mem=167G --cpus-per-task=5 --wrap="%s"', ...
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
                [status, cmdout] = system(chunk_decon_cmd, '-echo');
                tmpFullnames{ck} =  sprintf('%s/%s_%s_decon.tif', chunkDeconPath, chunkFnames{ck}(1:end-4), uuid);

                trial_counter(ck) = trial_counter(ck) + 1;    
            end

            if exist(chunkDeconFullpath, 'file') % || p.OverwriteDecon
                is_done_flag(ck) = true;
                continue;
            end
        end

        % wait 5 seconds if some tasks are still computing
        if ~all(is_done_flag | trial_counter >= maxTrialNum, 'all') 
            pause(10);
        end    
    end

    % combine the deconvovled segments and write file
    tic
    fprintf('Combining deconvovled chunks...\n')
    im = zeros(imSize, dtype);
    lol = floor(OL / 2);
    rol = ceil(OL / 2); 
    for ck = 1:nn
        chunkDeconFullpath = [chunkDeconPath '/' chunkFnames{ck}(1:end-4) '_decon.tif'];
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

    if p.Save16bit
        im = uint16(im);
    else
        im = single(im);
    end

    deconTmpPath = sprintf('%s_%s.tif', deconFullPath(1 : end - 4), uuid);
    try
        writetiff(im, deconTmpPath);
    catch
        options = struct();
        if p.lzwcompress
            options.compress = 'lzw';
        end
        options.message = false;
        options.overwrite = true;
        saveastiff(im, deconTmpPath, options);
    end
    movefile(deconTmpPath, deconFullPath);
    
    deconMIPPath = sprintf('%s/MIPs/', deconPath);
    mkdir(deconMIPPath);
    deconMIPFullPath = sprintf('%s%s_MIP_z.tif', deconMIPPath, fsname);
    writetiff(uint16(tmp_xy), deconMIPFullPath);
    toc

    % delete temporary files generated during the deconvolution
    rmdir(chunkPath, 's');

    fprintf('Completed GPU deconvolution of %s.\n', frameFullpath);
end

toc

end
