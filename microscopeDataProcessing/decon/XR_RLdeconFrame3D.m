function [] = XR_RLdeconFrame3D(frameFullpaths, pixelSize, dz, varargin)
% cuda Deconvolution for a single frame. It support both small file that
% can be fitted to a single GPU and also large files. For large files, it
% split files into chunks for deconvolution and the combine them together. It
% supports cluster computing for file chunks.  
%
% based on XR_matlabDeconFrame3D.m and GU_Decon.m 
% 
% Author: Xiongtao Ruan (03/15/2020)
% xruan (01/12/2021): add support for edge erosion and using existing eroded mask for edge erosion.  
% xruan (03/25/2021): add options for different versions of rl method
% xruan (06/10/2021): add support for threshold and debug mode in simplified version. 
% xruan (06/11/2021): add support for gpu computing for chuck decon in matlab decon wrapper
% xruan (07/13/2021): add support for the processing of flipped files
% (currently only add support for matlab decon)
% xruan (07/15/2021): add support for zarr input
% xruan (10/12/2021): setting overlap region based on the size of psf (and crop psf after psfgen)
% xruan (11/10/2021): add support for big zarr (may be tiff later) file chunk based decon


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
ip.addParameter('flipZstack', false, @islogical); 
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('EdgeSoften', 5, @isnumeric); % # ofxy px to soften
ip.addParameter('zEdgeSoften', 2, @isnumeric); % # ofxy px to soften
ip.addParameter('Crop', [], @isnumeric); % requires lower and higher values for cropping
ip.addParameter('zFlip', false, @islogical);
ip.addParameter('GenMaxZproj', [0,0,1] , @isnumeric);
ip.addParameter('ResizeImages', [] , @isnumeric);
ip.addParameter('dzPSF', 0.1 , @isnumeric); %in um
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('EdgeErosion', 8 , @isnumeric); % erode edges for certain size.
ip.addParameter('ErodeMaskfile', '', @ischar); % erode edges file
ip.addParameter('SaveMaskfile', false, @islogical); % save mask file for common eroded mask
% ip.addParameter('DoNotAdjustResForFFT', true , @islogical); % not crop chunks for deconvolution
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('fixIter', false, @islogical); % CPU Memory in Gb
ip.addParameter('errThresh', [], @isnumeric); % error threshold for simplified code
ip.addParameter('ChunkSize', [1600, 1600, 1600] , @isvector); % in y, x, z
ip.addParameter('Overlap', 200, @isnumeric); % block overlap
ip.addParameter('CPUMaxMem', 500, @isnumeric); % CPU Memory in Gb
ip.addParameter('BatchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256], @isnumeric); % block overlap
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('largeMethod', 'MemoryJobs', @ischar); % memory jobs, memory single, inplace. 
ip.addParameter('zarrFile', false, @islogical); % use zarr file as input
ip.addParameter('saveZarr', false, @islogical); % save as zarr
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('masterCPU', false, @islogical); % master node is a cpu node, which is just for large file deconvolution. 
ip.addParameter('GPUJob', false, @islogical); % use gpu for chuck deconvolution. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 5, @isnumeric);
ip.addParameter('cpuOnlyNodes', false, @islogical); % use gpu for chuck deconvolution. 
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 2, @isnumeric);
ip.addParameter('debug', false, @islogical);
ip.addParameter('psfGen', true, @islogical); % psf generation

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
flipZstack = pr.flipZstack;
Rotate = pr.Rotate;
Save16bit = pr.Save16bit;
Crop = pr.Crop;
zFlip = pr.zFlip;
GenMaxZproj = pr.GenMaxZproj;
ResizeImages = pr.ResizeImages;
RLMethod = pr.RLMethod;
GPUJob = pr.GPUJob;

EdgeErosion = pr.EdgeErosion;

ErodeMaskfile = pr.ErodeMaskfile;
if ~isempty(ErodeMaskfile) && exist(ErodeMaskfile, 'file')
    EdgeErosion = 0; % set EdgeErosion length as 0.
end
SaveMaskfile = pr.SaveMaskfile;

% check if background information available, if not, estimate background
% info. Currently use 99. 
Background = pr.Background;
if isempty(Background)
    Background = 99;
end

% simplified version related options
fixIter = pr.fixIter;
errThresh = pr.errThresh;
debug = pr.debug;
psfGen = pr.psfGen;

tic
OL = pr.Overlap;
ChunkSize = pr.ChunkSize;
CPUMaxMem = pr.CPUMaxMem;
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
largeFile = pr.largeFile;
largeMethod = pr.largeMethod;
zarrFile = pr.zarrFile;
saveZarr = pr.saveZarr;
parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
cpuOnlyNodes = pr.cpuOnlyNodes;
masterCompute = pr.masterCompute;
maxTrialNum = pr.maxTrialNum;
uuid = pr.uuid;

if isempty(uuid)
    uuid = get_uuid();
end

maxSubVolume = 1e11;
if GPUJob
    ChunkSize = round(pr.ChunkSize ./ [8, 8, 8]);
    OL = 200;
    maxSubVolume = 5e8;
    % masterCompute = gpuDeviceCount() > 0 & masterCompute;
end

OL_orig = OL;

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
    [pathstr, fsname, ext] = fileparts(frameFullpath);
    deconPath = pr.deconPath;
    if isempty(deconPath)
        deconPath = [pathstr, '/', 'matlab_decon'];
        mkdir(deconPath);
    end 

    % check file size
    [estMem] = XR_estimateComputingMemory(frameFullpath, {'deconvolution'}, 'cudaDecon', false);
    
    % The memory size is 500G
    if estMem > CPUMaxMem
        largeFile = true;
    end

    % first try single GPU deconvolution, if it fails split into multiple chunks
    if saveZarr
        deconFullPath = [deconPath '/' fsname '_decon.zarr'];
    else
        deconFullPath = [deconPath '/' fsname '_decon.tif'];
    end
    % deconTempPath = [deconPath '/' fname '_decon.tif'];
    if (exist(deconFullPath, 'file') || (saveZarr && exist(deconFullPath, 'dir'))) && ~pr.Overwrite
        disp('Deconvolution results already exist, skip it!');
        continue;
    end
    if ~largeFile
        if EdgeErosion > 0
            fprintf('Create eroded masks using raw data...\n');
            switch ext
                case {'.tif', '.tiff'}
                    try
                        im_raw = parallelReadTiff(frameFullpath);
                    catch
                        im_raw = readtiff(frameFullpath);
                    end
                case '.zarr'
                    im_raw = readzarr(frameFullpath);
            end
                    
            if flipZstack
                im_raw = flip(im_raw, 3);
            end
                
            im_bw = im_raw > 0;
            clear im_raw
            % pad to avoid not erosion if a pixel touching the boundary
            im_bw_pad = false(size(im_bw) + 2);
            im_bw_pad(2 : end - 1, 2 : end - 1, 2 : end - 1) = im_bw;
            im_bw_erode = imerode(im_bw_pad, strel('sphere', EdgeErosion));
            im_bw_erode = im_bw_erode(2 : end - 1, 2 : end - 1, 2 : end - 1);
            clear im_bw im_bw_pad
            
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
    end

    if ~largeFile
        tic
%         if ispc
%             softlink_cmd = sprintf('mklink %s_%s.tif %s', frameFullpath(1:end-4), uuid, frameFullpath); 
%             softlink_cmd = strrep(softlink_cmd, '/', '\');
%             unlink_cmd = sprintf('del %s_%s.tif', frameFullpath(1:end-4), uuid); 
%             unlink_cmd = strrep(unlink_cmd, '/', '\');            
%         else
%             softlink_cmd = sprintf('ln -s %s %s_%s.tif', frameFullpath, frameFullpath(1:end-4), uuid); 
%             unlink_cmd = sprintf('rm %s_%s.tif', frameFullpath(1:end-4), uuid);             
%         end
%         system(softlink_cmd);
        % frameTmpPath = sprintf('%s_%s.tif', frameFullpath(1:end-4), uuid); 
        deconTmpPath = sprintf('%s_%s_decon.tif', deconFullPath(1:end-10), uuid); 
        save3Dstack = false;
        im = RLdecon(frameFullpath, deconTmpPath, PSF, Background, DeconIter, dzPSF, dz, Deskew, [], SkewAngle, ...
            pixelSize, Rotate, Save16bit, Crop, zFlip, GenMaxZproj, ResizeImages, [], RLMethod, ...
            fixIter, errThresh, flipZstack, debug, 'save3Dstack', save3Dstack, 'psfGen', psfGen);
        toc
        % system(unlink_cmd);
            % im = readtiff(deconTmpPath);
        if EdgeErosion > 0
            fprintf('Erode edges of deconvolved data w.r.t. raw data...\n');        
            im = im .* cast(im_bw_erode, class(im));
        end

        if ~isempty(ErodeMaskfile) && exist(ErodeMaskfile, 'file')
            fprintf('Erode edges of deconvolved data using a predefined mask...\n');    
            try 
                im_bw_erode = parallelReadTiff(ErodeMaskfile);                    
            catch
                im_bw_erode = readtiff(ErodeMaskfile);
            end
            im = im .* cast(im_bw_erode, class(im));
        end
        deconTmpPath_eroded = sprintf('%s_%s_eroded.tif', deconFullPath(1:end-4), uuid);
        if Save16bit
            im = uint16(im);
        else
            im = single(im);
        end
        
        if saveZarr
            writezarr(im, deconTmpPath_eroded, 'blockSize', BlockSize);
        else
            writetiff(im, deconTmpPath_eroded);
        end
        movefile(deconTmpPath_eroded, deconFullPath);
        delete(deconTmpPath);
        deconTmpMIPPath = sprintf('%s/%s_%s_MIP_z.tif', deconPath, fsname, uuid); 
        delete(deconTmpMIPPath);

        deconMIPPath = sprintf('%s/MIPs/', deconPath);
        deconMIPFullPath = sprintf('%s%s_MIP_z.tif', deconMIPPath, fsname);
        writetiff(max(im,[],3), deconMIPFullPath);
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
    if strcmpi(largeMethod, 'inplace')
        RLdecon_large_in_place(frameFullpath, pixelSize, dz, deconPath, PSF, ...
            'Save16bit', Save16bit, 'Deskew', Deskew, 'SkewAngle', SkewAngle, ...
            'flipZstack', flipZstack, 'Background', Background, 'dzPSF', dzPSF, ...
            'DeconIter', DeconIter, 'fixIter', fixIter, 'BatchSize', BatchSize, ...
            'BlockSize', BlockSize, 'parseCluster', parseCluster, 'parseParfor', ...
            parseParfor, 'masterCompute', masterCompute, 'jobLogDir', jobLogDir, ...
            'cpuOnlyNodes', cpuOnlyNodes, 'GPUJob', GPUJob, 'uuid', uuid, 'debug', debug, ...
            'psfGen', psfGen);
        return;
    end
    
    % to do: put the code below to a function as in memory computing
    
    
    fprintf('Start Large-file RL Decon for %s...\n', fsname);
    
    % generate psf and use the cropped psf to decide the overlap region
    try
        pp = parallelReadTiff(PSF);
    catch
        pp = readtiff(PSF);
    end
    medFactor = 1.5;
    PSFGenMethod = 'masked';
    if Deskew
        dz_data_ratio = sind(SkewAngle);
    else
        dz_data_ratio = 1;
    end
    psf = psf_gen_new(pp, dzPSF, dz*dz_data_ratio, medFactor, PSFGenMethod);
    py = find(squeeze(sum(psf, [2, 3])));
    px = find(squeeze(sum(psf, [1, 3])));
    pz = find(squeeze(sum(psf, [1, 2])));
    cropSz = [min(py(1) - 1, size(psf, 1) - py(end)), min(px(1) - 1, size(psf, 2) - px(end)), min(pz(1) - 1, size(psf, 3) - pz(end))] - 1;
    cropSz = max(0, cropSz);
    bbox = [cropSz + 1, size(psf) - cropSz];
    
    % write generated psf to disk in case it will be deleted in the c
    psfgen_folder = sprintf('%s/psfgen/', deconPath);
    mkdir(psfgen_folder);
    [~, psf_fsn] = fileparts(PSF);
    psfgen_filename = sprintf('%s/%s.tif', psfgen_folder, psf_fsn);
    if ~exist(psfgen_filename, 'file')
        writetiff(psf, psfgen_filename);
    end
    
    % OL = min(OL, ceil((bbox(4 : 6) - bbox(1 : 3) + 1) / 2) * 3);
    OL = min(OL_orig, ceil((bbox(4 : 6) - bbox(1 : 3) + 1) / 2) * 2 + 10);
    
    tic
    fprintf(['reading ' fsname '...\n'])
    switch ext
        case {'.tif', '.tiff'}
            try
                im = parallelReadTiff(frameFullpath);
                % im = readtiff(frameFullpath);                
            catch
                im = readtiff(frameFullpath);
            end
        case '.zarr'
            bim = blockedImage(frameFullpath, 'Adapter', ZarrAdapter);
            im = gather(bim);
    end
    if flipZstack
        im = flip(im, 3);
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
            try
                im_bw_erode = parallelReadTiff(maskFullPath) > 0;
            catch
                im_bw_erode = readtiff(maskFullPath) > 0;                
            end
        end
        
        % save mask file as common one for other time points/channels
        if SaveMaskfile
            writetiff(uint8(im_bw_erode), maskTmpPath);
            movefile(maskTmpPath, maskFullPath);
        end
    end
    
    % calculate number of chunks to break the image file
    % [my, mx, mz] = size(im);
    % [xmin,xmax,ymin,ymax,zmin,zmax,nn] = GU_extract_subVolCoordinates(mx,my,mz,csx,csy,csz,OL);
    imSize = size(im);
    % dtype = class(im);
    if Save16bit
        dtype = 'uint16';
    else
        dtype = 'single';
    end
        
    if pr.debug
        [xmin,xmax,ymin,ymax,zmin,zmax,nn] = XR_subVolumeCoordinatesExtraction_test(imSize, 'ChunkSize', ChunkSize, 'overlapSize', OL);
    else
        [xmin,xmax,ymin,ymax,zmin,zmax,nn] = XR_subVolumeCoordinatesExtraction(imSize, 'ChunkSize', ChunkSize, 'overlapSize', OL, 'maxSubVolume', maxSubVolume);
    end
    
    % create a folder for the file and write out the chunks
    chunkPath = [deconPath '/' fsname];
    chunkDeconPath = [chunkPath '/' 'matlab_decon'];
    chunkDeconMIPPath = [chunkPath '/' 'matlab_decon' '/' 'MIPs'];
    mkdir(chunkPath);
    mkdir(chunkDeconPath);
    mkdir(chunkDeconMIPPath);
    
    if parseCluster
        if  ~exist(jobLogDir, 'dir')
            warning('The job log directory does not exist, use %s/job_logs as job log directory.', chunkPath)
            jobLogDir = sprintf('%s/job_logs', chunkPath);
            mkdir(jobLogDir);
        end
    end
    
    tic
    fprintf('Writing image chunks...\n')
    chunkFnames = cell(nn, 1);
    skip_flags = false(nn, 1);
    for ck = 1:nn
        chunkFnames{ck} = sprintf('chunk_%d_%d_%d_%d_%d_%d.tif', xmin(ck), ymin(ck), zmin(ck), xmax(ck), ymax(ck), zmax(ck));
        im_chunk = im(ymin(ck):ymax(ck), xmin(ck):xmax(ck), zmin(ck):zmax(ck));
        % for blank region, just skip it
        if all(im_chunk == 0, 'all')
            chunkFnames{ck} = '';
            skip_flags(ck) = true;
            continue;
        end
        if ~exist([chunkPath, '/', chunkFnames{ck}], 'file') || pr.Overwrite
            writetiff(im_chunk, [chunkPath, '/', chunkFnames{ck}(1 : end - 4), '_', uuid, '.tif']);
            movefile([chunkPath, '/', chunkFnames{ck}(1 : end - 4), '_', uuid, '.tif'], [chunkPath, '/', chunkFnames{ck}]);
        end        
    end
    
    clear im;
    toc

    tic
    fprintf('Deconvolving chunks...\n')
    funcStrs = cell(nn, 1);
    inputFullpaths = cell(nn, 1);
    outputFullpaths = cell(nn, 1);
    for ck = 1 : nn
        if skip_flags(ck)
            continue;
        end
        chunkFullpath = [chunkPath '/' chunkFnames{ck}];
        chunkDeconFullpath = [chunkDeconPath '/' chunkFnames{ck}(1:end-4) '_decon.tif'];
        inputFullpaths{ck} = chunkFullpath;
        outputFullpaths{ck} = chunkDeconFullpath;
    
        % no need to flip the chunks because the image is already flipped
        funcStrs{ck} =  sprintf(['RLdecon(''%s'',''%s'',''%s'',%.10f,%.10f,%.10f,%.10f,%s,[],', ...
                '%.10f,%.10f,%s,%s,[%s],%s,[%s],[%s],[],''%s'',%s,[%.20f],%s,%s,''psfGen'',%s)'], chunkFullpath, chunkDeconFullpath, ...
                PSF, Background, DeconIter, dzPSF, dz, string(Deskew), SkewAngle, pixelSize, string(Rotate), ...
                string(Save16bit),  strrep(num2str(Crop,'%d,'), ' ', ''), string(zFlip), num2str(GenMaxZproj, '%.10f,'), ...
                num2str(ResizeImages, '%.10f,'), RLMethod, string(fixIter), errThresh, string(false), string(debug), ...
                string(psfGen));
    end
    inputFullpaths(skip_flags) = [];
    outputFullpaths(skip_flags) = [];
    funcStrs(skip_flags) = [];
    
    if GPUJob
        maxJobNum = inf;
        cpusPerTask = 5;
        cpuOnlyNodes = false;
        taskBatchNum = 5;
        SlurmParam = '-p abc --qos abc_normal -n1 --mem=167G --gres=gpu:1';
    else
        maxJobNum = inf;
        cpusPerTask = 24;
        cpuOnlyNodes = true;
        taskBatchNum = 1;
        SlurmParam = '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M';
    end
        
    % submit jobs
    is_done_flag= slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', cpuOnlyNodes, 'SlurmParam', SlurmParam, ...
        'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster);

    if ~all(is_done_flag)
        slurm_cluster_generic_computing_wrapper(inputFullpaths, outputFullpaths, ...
            funcStrs, 'cpusPerTask', cpusPerTask, 'cpuOnlyNodes', cpuOnlyNodes, 'SlurmParam', SlurmParam, ...
            'maxJobNum', maxJobNum, 'taskBatchNum', taskBatchNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster);
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
            try
                tim = parallelReadTiff(chunkDeconFullpath);
            catch
                tim = readtiff(chunkDeconFullpath);
            end
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
        try
            im_bw_erode = parallelReadTiff(ErodeMaskfile);
        catch
            im_bw_erode = readtiff(ErodeMaskfile);
        end
        im = im .* cast(im_bw_erode, class(im));
    end
    

    if pr.Save16bit
        im = uint16(im);
    else
        im = single(im);
    end

    if savezarr 
        deconTmpPath = sprintf('%s_%s.zarr', deconFullPath(1 : end - 4), uuid);                
        writezarr(im, deconTmpPath, 'blockSize', BlockSize);
    else
        deconTmpPath = sprintf('%s_%s.tif', deconFullPath(1 : end - 4), uuid);        
        writetiff(im, deconTmpPath);
    end
    movefile(deconTmpPath, deconFullPath);
    
    tmp_xy = max(im,[],3);
    deconMIPPath = sprintf('%s/MIPs/', deconPath);
    mkdir(deconMIPPath);
    deconMIPFullPath = sprintf('%s%s_MIP_z.tif', deconMIPPath, fsname);
    writetiff(tmp_xy, deconMIPFullPath);
    toc

    % delete temporary files generated during the deconvolution
    rmdir(chunkPath, 's');

    fprintf('Completed RL deconvolution of %s.\n', frameFullpath);
end

toc

end
