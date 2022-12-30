function [] = XR_RLdeconFrame3D(frameFullpaths, xyPixelSize, dz, varargin)
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
ip.addRequired('xyPixelSize', @isnumeric);
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
ip.addParameter('dzPSF', 0.1 , @isnumeric); %in um
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('EdgeErosion', 8 , @isnumeric); % erode edges for certain size.
ip.addParameter('ErodeMaskfile', '', @ischar); % erode edges file
ip.addParameter('SaveMaskfile', false, @islogical); % save mask file for common eroded mask
% ip.addParameter('DoNotAdjustResForFFT', true , @islogical); % not crop chunks for deconvolution
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @isnumeric); 
ip.addParameter('OTFCumThresh', 0.9, @isnumeric); % OTF cumutative sum threshold
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x)); % decon in skewed space
ip.addParameter('fixIter', false, @islogical); % CPU Memory in Gb
ip.addParameter('errThresh', [], @isnumeric); % error threshold for simplified code
ip.addParameter('CPUMaxMem', 500, @isnumeric); % CPU Memory in Gb
ip.addParameter('BatchSize', [1024, 1024, 1024] , @isvector); % in y, x, z
ip.addParameter('BlockSize', [256, 256, 256], @isnumeric); % block overlap
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('largeMethod', 'MemoryJobs', @ischar); % memory jobs, memory single, inplace. 
ip.addParameter('saveZarr', false, @islogical); % save as zarr
ip.addParameter('deconMaskFns', {}, @iscell); % 2d masks to filter regions to decon, in xy, xz, yz order
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('masterCPU', false, @islogical); % master node is a cpu node, which is just for large file deconvolution. 
ip.addParameter('GPUJob', false, @islogical); % use gpu for chuck deconvolution. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 4, @isnumeric);
ip.addParameter('cpuOnlyNodes', false, @islogical); % use gpu for chuck deconvolution. 
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 2, @isnumeric);
ip.addParameter('debug', false, @islogical);
ip.addParameter('saveStep', 5, @isnumeric); % save intermediate results every given iterations
ip.addParameter('psfGen', true, @islogical); % psf generation

ip.parse(frameFullpaths, xyPixelSize, dz, varargin{:});

if ischar(frameFullpaths)
    frameFullpaths = {frameFullpaths};
end

pr = ip.Results;

PSFfile = pr.PSFfile;
if isempty(PSFfile)
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
RLMethod = pr.RLMethod;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
skewed = pr.skewed;
GPUJob = pr.GPUJob;

EdgeErosion = pr.EdgeErosion;

ErodeMaskfile = pr.ErodeMaskfile;
if ~isempty(ErodeMaskfile) && exist(ErodeMaskfile, 'file')
    EdgeErosion = 0; % set EdgeErosion length as 0.
end
SaveMaskfile = pr.SaveMaskfile;

% check if background information available. Currently use 100. 
Background = pr.Background;
if isempty(Background)
    Background = 100;
end

% simplified version related options
fixIter = pr.fixIter;
errThresh = pr.errThresh;
debug = pr.debug;
saveStep = pr.saveStep;
psfGen = pr.psfGen;

tic
CPUMaxMem = pr.CPUMaxMem;
BatchSize = pr.BatchSize;
BlockSize = pr.BlockSize;
largeFile = pr.largeFile;
largeMethod = pr.largeMethod;
saveZarr = pr.saveZarr;
deconMaskFns = pr.deconMaskFns;

parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
cpuOnlyNodes = pr.cpuOnlyNodes;
masterCompute = pr.masterCompute;
maxTrialNum = pr.maxTrialNum;
uuid = pr.uuid;

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
        deconFullPath = [deconPath '/' fsname '.zarr'];
    else
        deconFullPath = [deconPath '/' fsname '.tif'];
    end
    if (exist(deconFullPath, 'file') || (saveZarr && exist(deconFullPath, 'dir'))) && ~pr.Overwrite
        disp('Deconvolution results already exist, skip it!');
        continue;
    end
    if ~largeFile
        if EdgeErosion > 0
            fprintf('Create eroded masks using raw data...\n');
            switch ext
                case {'.tif', '.tiff'}
                    im_raw = readtiff(frameFullpath);
                case '.zarr'
                    im_raw = readzarr(frameFullpath);
            end
                    
            if flipZstack
                im_raw = flip(im_raw, 3);
            end
                
            im_bw_erode = im_raw > 0;
            clear im_raw;
            im_bw_erode([1, end], :, :) = false;
            im_bw_erode(:, [1, end], :) = false;
            im_bw_erode(:, :, [1, end]) = false;
            im_bw_erode = bwdist(~im_bw_erode) > EdgeErosion - 1;
            
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
        deconTmpPath = sprintf('%s_%s_decon.tif', deconFullPath(1:end-10), uuid); 
        inputFn = frameFullpath;
        outputFn = deconTmpPath;
        DSRCombined = true;
        Reverse = true;
        scaleFactor = 1;
        save3Dstack = [false, false, false];
        mipAxis = [0, 0, 0];
    
        im = RLdecon(inputFn, outputFn, PSFfile, xyPixelSize, dz, dzPSF, ...
            'rawdata', [], 'Save16bit', Save16bit, 'SkewAngle', SkewAngle, ...
            'Deskew', Deskew, 'Rotate', Rotate, 'DSRCombined', DSRCombined, 'Reverse', Reverse, ...
            'Background', Background, 'DeconIter', DeconIter, 'RLMethod', RLMethod, ...
            'skewed', skewed, 'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, ...
            'fixIter', fixIter, 'scaleFactor', scaleFactor, 'errThresh', errThresh, ...
            'saveStep', saveStep, 'useGPU', GPUJob, 'psfGen', psfGen, 'debug', debug, ...
            'save3Dstack', save3Dstack, 'mipAxis', mipAxis);
        toc

        if EdgeErosion > 0
            fprintf('Erode edges of deconvolved data w.r.t. raw data...\n');        
            im = im .* cast(im_bw_erode, class(im));
        end

        if ~isempty(ErodeMaskfile) && exist(ErodeMaskfile, 'file')
            fprintf('Erode edges of deconvolved data using a predefined mask...\n');    
            im_bw_erode = readtiff(ErodeMaskfile);
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
        if ~exist(deconMIPPath, 'dir')
            mkdir(deconMIPPath);
        end
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
        RLdecon_large_in_place(frameFullpath, xyPixelSize, dz, deconPath, PSFfile, ...
            'Save16bit', Save16bit, 'Deskew', Deskew, 'SkewAngle', SkewAngle, ...
            'flipZstack', flipZstack, 'Background', Background, 'dzPSF', dzPSF, ...
            'DeconIter', DeconIter, 'RLMethod', RLMethod, 'skewed', skewed, ...
            'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, 'fixIter', fixIter, ...
            'BatchSize', BatchSize, 'BlockSize', BlockSize, 'deconMaskFns', deconMaskFns, ...
            'parseCluster', parseCluster, 'parseParfor', parseParfor, 'masterCompute', masterCompute, ...
            'jobLogDir', jobLogDir,  'cpuOnlyNodes', cpuOnlyNodes, 'cpusPerTask', cpusPerTask, ...
            'GPUJob', GPUJob, 'uuid', uuid, 'debug', debug, 'psfGen', psfGen);
        return;
    end
    
    % to do: put the code below to a function as in memory computing
    if strcmpi(largeMethod, 'inmemory')
        RLdecon_large_in_memory(frameFullpath, PSFfile, deconFullPath, xyPixelSize, dz, ...
            'Save16bit', Save16bit, 'Deskew', Deskew, 'SkewAngle', SkewAngle, ...
            'flipZstack', flipZstack, 'Background', Background, 'dzPSF', dzPSF, ...
            'DeconIter', DeconIter, 'RLMethod', RLMethod, 'skewed', skewed, ...
            'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, 'EdgeErosion', EdgeErosion, ...
            'fixIter', fixIter,'BatchSize', BatchSize, 'saveZarr', saveZarr, ...
            'deconMaskFns', deconMaskFns, 'cpuOnlyNodes', cpuOnlyNodes, 'cpusPerTask', cpusPerTask, ...
            'useGPU', GPUJob, 'uuid', uuid, 'debug', debug, 'psfGen', psfGen);
        return;
    end
end

toc

end
