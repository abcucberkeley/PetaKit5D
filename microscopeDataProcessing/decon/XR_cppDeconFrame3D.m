function [] = XR_cppDeconFrame3D(frameFullpaths, pixelSize, dz, varargin)
% CPP Deconvolution for a single frame. It support both small file that
% can be fitted to a single node and also large files. For large files, it
% split files into chunks for deconvolution and the combine them together. It
% supports cluster computing for file chunks.  
%
% It uses Lin Shao's CPP deconvolution function, and is based on
% XR_cudaDeconFrame3D.m
%
% 
% Author: Xiongtao Ruan (04/15/2020)
% xruan (06/25/2020): add pad option for a single file to avoid cropping.
%   To avoid artifects in the first and last slice, pad one slice before
%   first slice. 
% xruan (07/07/2020): add option for dxPSF
% xruan (07/08/2020): add option to remove edge region
% xruan (07/15/2020): change maxSubvolume to 2.0e9 (uint16) and 1.0e9
%   (single), also decide to use largeFile if the result file size is over 4Gb.
% xruan (07/18/2020): change Pad and Edge Erosion options to cppDecon itself. 
% xruan (07/20/2020): readd Edge Erosion options in Matlab as cpp erosion is very slow. 
% xruan (07/21/2020): add support for using existing eroded mask for edge erosion.  
% xruan (07/29/2020): Increase the size limit (from 4GB to 100GB) to use largeFile option, as cppDecon support bigTiff now.  
% xruan (10/15/2020): Add support for zarr format


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpaths', @(x) ischar(x) || iscell(x));
ip.addRequired('pixelSize', @isnumeric);
ip.addRequired('dz', @isnumeric);
ip.addOptional('deconPath', '', @(x) ischar(x) || isempty(x));
ip.addParameter('PSFfile', '' , @ischar);
ip.addParameter('cppDeconPath', '/global/home/groups/software/sl-7.x86_64/modules/RLDecon_CPU/20200718/build-cluster/cpuDeconv', @ischar);
ip.addParameter('loadModules', 'module load gcc/4.8.5; module load fftw/3.3.6-gcc; module load boost/1.65.1-gcc; module load libtiff/4.1.0; ', @ischar);
ip.addParameter('Overwrite', false , @islogical);
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
ip.addParameter('dxPSF', 0.108, @isnumeric); %in um
ip.addParameter('dzPSF', 0.1 , @isnumeric); %in um
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('padForVolume', true , @islogical); % pad for single volume to "good number"
ip.addParameter('EdgeErosion', 8 , @isnumeric); % erode edges for certain size.
ip.addParameter('ErodeMaskfile', '', @ischar); % erode edges file
ip.addParameter('SaveMaskfile', false, @islogical); % save mask file for common eroded mask
% ip.addParameter('DoNotAdjustResForFFT', true , @islogical); % not crop chunks for deconvolution
ip.addParameter('ChunkSize', [1000,1000,1000] , @isvector); % in y, x, z
ip.addParameter('Overlap', 200, @isnumeric); % block overlap
ip.addParameter('MaxMem', 500, @isnumeric); % GPU Memory in Gb
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('masterCPU', false, @islogical); % master node is a cpu node, which is just for large file deconvolution. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 2, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 2, @isnumeric);
ip.addParameter('debug', false, @islogical);

ip.parse(frameFullpaths, pixelSize, dz, varargin{:});

if ischar(frameFullpaths)
    frameFullpaths = {frameFullpaths};
end

p = ip.Results;

PSFfile = p.PSFfile;
if isempty(PSFfile)
    error('You should provide a PSF file for the frame...\n');
end

% parameters
Deskew_str = '';
if p.Deskew
    Deskew_str = sprintf('-D');
end

SkewAngle_str = '';
if p.SkewAngle
    SkewAngle_str = sprintf('-a %.f', p.SkewAngle);
end

Rotate_str = '';
if p.Rotate
    Rotate_str = sprintf('-R');
end

Crop_str = '';
if ~isempty(p.Crop)
    Crop_str = sprintf('-C %s', num2str(p.Crop));
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

Save16bit = p.Save16bit;
Save16bit_str = '';
if Save16bit
    Save16bit_str = '-u';
end

padForVolume = p.padForVolume;
% pad 2 by default
Pad_str = '';
if padForVolume
    Pad_str = '--Pad 2';
end

EdgeErosion = p.EdgeErosion;
EdgeErosion_str = '';
% if EdgeErosion > 0
%     EdgeErosion_str = sprintf('--EdgeErosion %d', EdgeErosion);
% end

ErodeMaskfile = p.ErodeMaskfile;
maskfile_str = '';
if ~isempty(ErodeMaskfile) && exist(ErodeMaskfile, 'file')
    EdgeErosion = 0; % set EdgeErosion length as 0.
    maskfile_str = sprintf('--mask-file %s ', ErodeMaskfile);
end

SaveMaskfile = p.SaveMaskfile;

system_lib_str = p.loadModules;
cpuDeconPath = p.cppDeconPath;
dxPSF = p.dxPSF;
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
params = sprintf(' --drdata %.10f --drpsf %.10f -Z %.10f -z %.10f -b %.10f -i %d %s  %s %s %s %s %s %s %s %s %s %s %s %s', ...
    pixelSize, dxPSF, dzPSF, dz, Background, DeconIter, Deskew_str, SkewAngle_str, Rotate_str, Crop_str, ...
    SaveDeskew_str, EdgeSoften_str, GenMaxZproj_str, zEdgeSoften_str, Save16bit_str, Pad_str, EdgeErosion_str);


OL = p.Overlap;
ChunkSize = p.ChunkSize;

largeFile = p.largeFile;
parseCluster = p.parseCluster;
cpusPerTask = p.cpusPerTask;
jobLogDir = p.jobLogDir;
masterCompute = p.masterCompute;
maxTrialNum = p.maxTrialNum;
uuid = p.uuid;
if isempty(uuid)
    uuid = get_uuid();
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
    if strcmp(ext, '.zarr')
        zarrFile = true;
    else
        zarrFile = false;
    end
    
    deconPath = p.deconPath;
    if isempty(deconPath)
        deconPath = [pathstr, '/', 'CPPdecon'];
    end 
    
    % first try single GPU deconvolution, if it fails split into multiple chunks
    deconFullPath = [deconPath '/' fsname '_decon.tif'];
    % if zarrFile
    %    deconFullPath_zarr = [deconPath '/' fsname '_decon.zarr'];
    % end
    % deconTempPath = [deconPath '/' fname '_decon.tif'];
    if exist(deconFullPath, 'file') && ~p.Overwrite
        disp('Deconvolution results already exist, skip it!');
        continue;
    end

    % check file size
    if ~zarrFile
        [estMem, estGPUMem, rawImageSize, sz] = XR_estimateComputingMemory(frameFullpath, {'deconvolution'}, 'cudaDecon', false);
        % due to 4Gb size limit in cppDecon program, we need to directly use
        % large file computing if the final size is greater than 2GB
        if rawImageSize / 2 > 2
            largeFile = true;
        end
    else
        largeFile = true;
    end
    if ~largeFile
        if EdgeErosion > 0
            fprintf('Create eroded masks using raw data...\n');
            im_raw = readtiff(frameFullpath);
            im_bw = im_raw > 0;
            % pad to avoid not erosion if a pixel touching the boundary
            im_bw_pad = false(size(im_raw) + 2);
            im_bw_pad(2 : end - 1, 2 : end - 1, 2 : end - 1) = im_bw;
            im_bw_erode = imerode(im_bw_pad, strel('sphere', EdgeErosion));
            im_bw_erode = im_bw_erode(2 : end - 1, 2 : end - 1, 2 : end - 1);
            clear im_raw im_bw im_bw_pad
            
            % save mask file as common one for other time points/channels
            if SaveMaskfile
                maskPath = [deconPath, '/', 'Masks'];
                if ~exist(maskPath, 'dir')
                    mkdir(maskPath);
                end
                maskFullPath = sprintf('%s/%s_eroded.tif', maskPath, fsname);
                maskTmpPath = sprintf('%s/%s_eroded_%s.tif', maskPath, fsname, uuid);
                writetiff(uint8(im_bw_erode), maskTmpPath);
                movefile(maskTmpPath, maskFullPath);
            end        
        end
        
        % perform mask-based edge erosion in cpp
        params_use = sprintf('%s %s', params, maskfile_str);
        
        if ispc
            DeconCommand = [system_lib_str cpuDeconPath params_use '--input-file ' frameFullpath ' --psf-file ' PSFfile ];
        else
            softlink_cmd = sprintf('ln -s %s %s_%s.tif', frameFullpath, frameFullpath(1:end-4), uuid);
            decon_cmd = [system_lib_str cpuDeconPath params_use '--input-file ' frameFullpath(1:end-4) '_' uuid '.tif' ' --psf-file ' PSFfile ];
            
            if EdgeErosion > 0
                rename_cmd = sprintf('mv %s/%s_%s_decon.tif %s/%s_decon_raw.tif', deconPath, fsname, uuid, deconPath, fsname);                
                rename_mip_cmd = sprintf('rm %s/MIPs/%s_%s_MIP_z.tif', deconPath, fsname, uuid);
            else
                rename_cmd = sprintf('mv %s/%s_%s_decon.tif %s/%s_decon.tif', deconPath, fsname, uuid, deconPath, fsname);                
                rename_mip_cmd = sprintf('mv %s/MIPs/%s_%s_MIP_z.tif %s/MIPs/%s_decon_MIP_z.tif', deconPath, fsname, uuid, deconPath, fsname);
            end
            unlink_cmd = sprintf('rm %s_%s.tif', frameFullpath(1:end-4), uuid);
            DeconCommand = sprintf('%s; %s; %s; %s; %s', softlink_cmd, decon_cmd, rename_cmd, rename_mip_cmd, unlink_cmd);
        end

        if p.Verbose
            [~, ~] = system(DeconCommand, '-echo');
        else
            [~, ~] = system(DeconCommand);
        end
        
        % check if the (intermediate) result exists.
        incompleteResult = false;
        if EdgeErosion > 0
            dcmFullpath = sprintf('%s/%s_decon_raw.tif', deconPath, fsname);
        else
            dcmFullpath = sprintf('%s/%s_decon.tif', deconPath, fsname);
        end
        if ~exist(dcmFullpath, 'file')
            incompleteResult = true;
        end
        
        if ~incompleteResult
            dcSize = getImageSize(dcmFullpath);
            if any(dcSize ~= sz)
                fprintf('The deconvolved file %s is not correctly saved! Will use large-file deconvolution method!\n', deconFrameFullpath_use)
                incompleteResult = true;
            end
        end

        if ~incompleteResult && EdgeErosion > 0    
            fprintf('Erode edges of deconvolved data w.r.t. raw data...\n');            
            im = readtiff(dcmFullpath);
            im = im .* cast(im_bw_erode, class(im));

            % save results
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

            % save MIP image
            deconMIPPath = sprintf('%s/MIPs/', deconPath);
            mkdir(deconMIPPath);
            deconMIPFullPath = sprintf('%s%s_decon_MIP_z.tif', deconMIPPath, fsname);
            writetiff(uint16(tmp_xy), deconMIPFullPath);
            delete(dcmFullpath);
        end
        if incompleteResult
            delete(dcmFullpath);
        end
        
        if exist(deconFullPath, 'file')
            fprintf('Completed CPP deconvolution of %s.\n', frameFullpath);
            continue;
        else
            % For small file size, if it failed, just throw the error. 
            if estMem < 0.1 * p.MaxMem
                error('Deconvolution is failed with file size %d', estMem);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Large-file deconvolution if the file is too large or single file
    % deconvolution fails.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('Start Large-file CPP Decon for %s...\n', fname);

    tic
    fprintf(['reading ' fname '...\n'])
    if zarrFile
        bim = blockedImage(frameFullpath, 'Adapter', ZarrAdapter);
        im = gather(bim);
    else
        im = readtiff(frameFullpath);
    end
    toc
    
    if EdgeErosion > 0
        if SaveMaskfile
            maskPath = [deconPath, '/', 'Masks'];
            if ~exist(maskPath, 'dir')
                mkdir(maskPath);
            end
            maskFullPath = sprintf('%s/%s_eroded.tif', maskPath, fsname);
            maskTmpPath = sprintf('%s/%s_eroded_%s.tif', maskPath, fsname, uuid);
        end
        if ~(SaveMaskfile && exist(maskFullPath, 'file'))
            im_bw = im > 0;
            % pad to avoid not erosion if a pixel touching the boundary
            im_bw_pad = false(size(im) + 2);
            im_bw_pad(2 : end - 1, 2 : end - 1, 2 : end - 1) = im_bw;
            im_bw_erode = imerode(im_bw_pad, strel('sphere', EdgeErosion));
            im_bw_erode = im_bw_erode(2 : end - 1, 2 : end - 1, 2 : end - 1);
            clear im_bw im_bw_pad
        else
            im_bw_erode = readtiff(maskFullPath) > 0;
        end
        
        % save mask file as common one for other time points/channels
        if SaveMaskfile
            writetiff(uint8(im_bw_erode), maskTmpPath);
            movefile(maskTmpPath, maskFullPath);
        end
    end

    % calculate number of chunks to break the image file    
    imSize = size(im);
    % dtype = class(im);
    if Save16bit
        dtype = 'uint16';
    else
        dtype = 'single';
    end
    
    if p.debug
        [xmin,xmax,ymin,ymax,zmin,zmax,nn] = XR_subVolumeCoordinatesExtraction_test(imSize, 'ChunkSize', ChunkSize, 'overlapSize', OL);
    else
        % maxSubVolume ~20G
        if p.Save16bit
            maxSubVolume = 2e10;
        else
            maxSubVolume = 1e10;
        end
            
        [xmin,xmax,ymin,ymax,zmin,zmax,nn] = XR_subVolumeCoordinatesExtraction(imSize, 'ChunkSize', ChunkSize, 'overlapSize', OL, 'maxSubVolume', maxSubVolume);        
    end
    
    % create a folder for the file and write out the chunks
    chunkPath = [deconPath '/' fsname];
    chunkDeconPath = [chunkPath '/' 'CPPdecon'];
    chunkDeconMIPPath = [chunkPath '/' 'CPPdecon' '/' 'MIPs'];
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
        if ~exist([chunkPath, '/', chunkFnames{ck}], 'file') || p.Overwrite
            writetiff(im_chunk, [chunkPath, '/', chunkFnames{ck}(1 : end - 4), '_', uuid, '.tif']);
            movefile([chunkPath, '/', chunkFnames{ck}(1 : end - 4), '_', uuid, '.tif'], [chunkPath, '/', chunkFnames{ck}]);
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
            
            chunkDeconFullpath = [chunkDeconPath '/' chunkFnames{ck}(1:end-4) '_decon.tif'];
            if exist(chunkDeconFullpath, 'file') || isempty(chunkFnames{ck}) % || p.OverwriteDecon
                is_done_flag(ck) = true;
                continue;
            end
            
            % Create a softlink with uuid of the chunk and rename the file
            % after write to disk is done
            uuid = get_uuid();
            
            softlink_cmd = sprintf('ln -s %s %s_%s.tif', [chunkPath, '/', chunkFnames{ck}], [chunkPath, '/', chunkFnames{ck}(1:end-4)], uuid);
            decon_cmd = [system_lib_str cpuDeconPath params '--input-file ' chunkPath, '/', chunkFnames{ck}(1:end-4), '_', uuid, '.tif' ' --psf-file ' PSFfile];
            rename_cmd = sprintf('mv %s_%s_decon.tif %s_decon.tif', [chunkDeconPath, '/', chunkFnames{ck}(1:end-4)], uuid, [chunkDeconPath, '/', chunkFnames{ck}(1:end-4)]);
            chunk_decon_cmd = sprintf('%s; %s; %s', softlink_cmd, decon_cmd, rename_cmd);
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
                    [estMem, estGPUMem] = XR_estimateComputingMemory([chunkPath, '/', chunkFnames{ck}], {'deconvolution'}, 'cudaDecon', false);

                    if cpusPerTask * 20 < estMem
                        cpusPerTask = min(24, ceil(estMem / 20));
                    end

                    cmd = sprintf('sbatch --array=%d --constraint=c24 -o %s -e %s -p abc --qos abc_normal -n1 --mem-per-cpu=21418M --cpus-per-task=%d --wrap="%s"', ...
                        task_id, job_log_fname, job_log_error_fname, cpusPerTask, chunk_decon_cmd);

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

        % wait 10 seconds if some tasks are still computing
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
    if EdgeErosion > 0
        fprintf('Erode edges of deconvolved data w.r.t. raw data...\n');        
        im = im .* cast(im_bw_erode, class(im));
    end
    
    if ~isempty(ErodeMaskfile) && exist(ErodeMaskfile, 'file')
        fprintf('Erode edges of deconvolved data using a predefined mask...\n');                
        im_bw_erode = readtiff(ErodeMaskfile);
        im = im .* cast(im_bw_erode, class(im));
    end
    
    tmp_xy = max(im,[],3);
    if p.Save16bit
        im = uint16(im);
    else
        im = single(im);
    end

    deconTmpPath = sprintf('%s_%s.tif', deconFullPath(1 : end - 4), uuid);
    writetiff(im, deconTmpPath);
    movefile(deconTmpPath, deconFullPath);
    
    deconMIPPath = sprintf('%s/MIPs/', deconPath);
    mkdir(deconMIPPath);
    deconMIPFullPath = sprintf('%s%s_decon_MIP_z.tif', deconMIPPath, fsname);
    writetiff(uint16(tmp_xy), deconMIPFullPath);
    toc

    % delete temporary files generated during the deconvolution
    if exist(chunkPath, 'dir')
        rmdir(chunkPath, 's');
    end

    fprintf('Completed CPP deconvolution of %s.\n', frameFullpath);
end

toc

end
