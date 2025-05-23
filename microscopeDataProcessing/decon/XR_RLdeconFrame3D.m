function [] = XR_RLdeconFrame3D(frameFullpaths, xyPixelSize, dz, deconPath, varargin)
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
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('Rotate', false, @islogical);
ip.addParameter('Deskew', false, @islogical);
ip.addParameter('SkewAngle', -32.45, @isnumeric);
ip.addParameter('flipZstack', false, @islogical); 
ip.addParameter('Reverse', true, @islogical);
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('EdgeSoften', 5, @isnumeric); % # ofxy px to soften
ip.addParameter('zEdgeSoften', 2, @isnumeric); % # ofxy px to soften
% ip.addParameter('Crop', [], @isnumeric); % requires lower and higher values for cropping
ip.addParameter('dzPSF', 0.1 , @isnumeric); %in um
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('EdgeErosion', 8 , @isnumeric); % erode edges for certain size.
ip.addParameter('ErodeMaskfile', '', @ischar); % erode edges file
ip.addParameter('SaveMaskfile', false, @islogical); % save mask file for common eroded mask
% ip.addParameter('DoNotAdjustResForFFT', true , @islogical); % not crop chunks for deconvolution
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @isnumeric); 
ip.addParameter('OTFCumThresh', 0.9, @isnumeric); % OTF cumutative sum threshold
ip.addParameter('hannWinBounds', [0.8, 1.0], @isnumeric); % apodization range for distance matrix
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x)); % decon in skewed space
ip.addParameter('fixIter', false, @islogical); % CPU Memory in Gb
ip.addParameter('errThresh', [], @isnumeric); % error threshold for simplified code
ip.addParameter('CPUMaxMem', 500, @isnumeric); % CPU Memory in Gb
ip.addParameter('batchSize', [1024, 1024, 1024] , @isnumeric); % in y, x, z
ip.addParameter('blockSize', [256, 256, 256], @isnumeric); % block overlap
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('largeMethod', 'inmemory', @ischar); % memory jobs, memory single, inplace. 
ip.addParameter('saveZarr', false, @islogical); % save as zarr
ip.addParameter('save3Dstack', [true, true, true], @islogical); % save 3d stack
ip.addParameter('mipAxis', [0, 0, 1], @isnumeric); % save MIPs, in y, x, z
ip.addParameter('dampFactor', 1, @isnumeric); % damp factor for decon result
ip.addParameter('scaleFactor', [], @isnumeric); % scale factor for decon result
ip.addParameter('deconOffset', 0, @isnumeric); % offset for decon result
ip.addParameter('maskFullpaths', {}, @iscell); % 2d masks to filter regions to decon, in xy, xz, yz order
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('parseParfor', false, @islogical);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('masterCPU', false, @islogical); % master node is a cpu node, which is just for large file deconvolution. 
ip.addParameter('GPUJob', false, @islogical); % use gpu for chuck deconvolution. 
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('cpusPerTask', 4, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('maxTrialNum', 3, @isnumeric);
ip.addParameter('unitWaitTime', 2, @isnumeric);
ip.addParameter('debug', false, @islogical);
ip.addParameter('saveStep', 5, @isnumeric); % save intermediate results every given iterations
ip.addParameter('psfGen', true, @islogical); % psf generation
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);
ip.addParameter('GPUConfigFile', '', @ischar);

ip.parse(frameFullpaths, xyPixelSize, dz, deconPath, varargin{:});

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
Reverse = pr.Reverse;
save16bit = pr.save16bit;
% Crop = pr.Crop;
RLMethod = lower(pr.RLMethod);
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
hannWinBounds = pr.hannWinBounds;
skewed = pr.skewed;
GPUJob = pr.GPUJob;

EdgeErosion = pr.EdgeErosion;

ErodeMaskfile = pr.ErodeMaskfile;
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

CPUMaxMem = pr.CPUMaxMem;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
largeFile = pr.largeFile;
largeMethod = pr.largeMethod;
saveZarr = pr.saveZarr;
save3Dstack = pr.save3Dstack;
mipAxis = pr.mipAxis;

dampFactor = pr.dampFactor;
scaleFactor = pr.scaleFactor;
deconOffset = pr.deconOffset;
maskFullpaths = pr.maskFullpaths;

parseCluster = pr.parseCluster;
parseParfor = pr.parseParfor;
jobLogDir = pr.jobLogDir;
cpusPerTask = pr.cpusPerTask;
masterCompute = pr.masterCompute;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;
GPUConfigFile = pr.GPUConfigFile;

if isempty(uuid)
    uuid = get_uuid();
end

% tic
if parseCluster 
    [status, ~] = system('sinfo');
    if status ~= 0
        warning('A slurm-based computing cluster is not exist. Set parseCluster as false.')
        parseCluster = false;
    end
end

save3Dstack(2) = save3Dstack(2) & Deskew;
save3Dstack(3) = save3Dstack(3) & Rotate;

if ~largeFile && ~any(save3Dstack) && ~any(mipAxis > 0)
    error('Either the 3D stack or the MIP for at least one axis needs to be saved!')
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
    if save3Dstack(3)
        DSRPath = sprintf('%s/DSR/', deconPath);
        mkdir(DSRPath);
    end
    % check file size
    % [estMem] = XR_estimateComputingMemory(frameFullpath, {'deconvolution'}, 'cudaDecon', false);
    imSize = getImageSize(frameFullpath);
    memRequired = prod(imSize) * 4 / 2^30 * 10;
    
    % The memory size is 500G
    if memRequired > CPUMaxMem
        largeFile = true;
    end

    % first try single GPU deconvolution, if it fails split into multiple chunks
    if saveZarr
        deconFullPath = [deconPath '/' fsname '.zarr'];
    else
        deconFullPath = [deconPath '/' fsname '.tif'];
    end
    if save3Dstack(3)
        if saveZarr
            dsrFullPath = [DSRPath '/' fsname '.zarr'];
        else
            dsrFullPath = [DSRPath '/' fsname '.tif'];
        end
    end    
    if (exist(deconFullPath, 'file') || (saveZarr && exist(deconFullPath, 'dir'))) && ~pr.Overwrite
        disp('Deconvolution results already exist, skip it!');
        continue;
    end
    if ~largeFile
        if EdgeErosion > 0 && SaveMaskfile
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
            
            im_bw_erode = decon_mask_edge_erosion(im_raw > 0, EdgeErosion);
            clear im_raw;

            % save mask file as common one for other time points/channels
            if SaveMaskfile
                maskPath = [deconPath, '/', 'Masks'];
                if ~exist(maskPath, 'dir')
                    mkdir(maskPath);
                end
                maskFullPath = sprintf('%s/%s_eroded.zarr', maskPath, fsname);
                maskTmpPath = sprintf('%s/%s_eroded_%s.zarr', maskPath, fsname, uuid);
                writezarr(uint8(im_bw_erode), maskTmpPath, 'blockSize', blockSize);
                if exist(maskFullPath, 'dir')
                    rmdirs(maskFullPath, 's');
                end
                movefile(maskTmpPath, maskFullPath);
            end
            ErodeMaskfile = maskFullPath;
        end

        % deconTmpPath = sprintf('%s_%s_decon.tif', deconFullPath(1:end-10), uuid); 
        inputFn = frameFullpath;
        outputFn = deconFullPath;
        DSRCombined = true;
        % Reverse = true;
        % save3Dstack = [false, false, false];
        % mipAxis = [0, 0, 0];

        t0 = tic;
        RLdecon(inputFn, outputFn, PSFfile, xyPixelSize, dz, dzPSF, 'rawdata', [], ...
            'save16bit', save16bit, 'SkewAngle', SkewAngle, 'Deskew', Deskew, ...
            'Rotate', Rotate, 'DSRCombined', DSRCombined, 'Reverse', Reverse, ...
            'Background', Background, 'DeconIter', DeconIter, 'RLMethod', RLMethod, ...
            'skewed', skewed, 'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, ...
            'hannWinBounds', hannWinBounds, 'saveZarr', saveZarr, 'blockSize', blockSize, ...
            'fixIter', fixIter, 'dampFactor', dampFactor, 'scaleFactor', scaleFactor, ...
            'deconOffset', deconOffset, 'EdgeErosion', EdgeErosion, 'ErodeMaskfile', ErodeMaskfile, ...
            'errThresh', errThresh, 'saveStep', saveStep, 'useGPU', GPUJob, ...
            'psfGen', psfGen, 'debug', debug, 'save3Dstack', save3Dstack, 'mipAxis', mipAxis);

        if (save3Dstack(1) && (exist(deconFullPath, 'file') || exist(deconFullPath, 'dir'))) ...
                || (save3Dstack(3) && (exist(dsrFullPath, 'file') || exist(dsrFullPath, 'dir')))
            fprintf('Done! ');
            toc(t0);
            continue;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Large-file deconvolution if the file is too large or single file
    % deconvolution fails.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(largeMethod, 'inplace')
        RLdecon_large_in_place(frameFullpath, xyPixelSize, dz, deconPath, PSFfile, ...
            'save16bit', save16bit, 'Deskew', Deskew, 'SkewAngle', SkewAngle, ...
            'flipZstack', flipZstack, 'Background', Background, 'dzPSF', dzPSF, ...
            'DeconIter', DeconIter, 'RLMethod', RLMethod, 'skewed', skewed, ...
            'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, 'hannWinBounds', hannWinBounds, ...
            'fixIter', fixIter, 'batchSize', batchSize, 'blockSize', blockSize, ...
            'dampFactor', dampFactor, 'scaleFactor', scaleFactor, 'deconOffset', deconOffset, ...
            'EdgeErosion', EdgeErosion, 'maskFullpaths', maskFullpaths, 'parseCluster', parseCluster, ...
            'parseParfor', parseParfor, 'masterCompute', masterCompute, 'jobLogDir', jobLogDir, ...
            'cpusPerTask', cpusPerTask, 'GPUJob', GPUJob, 'uuid', uuid, 'debug', debug, ...
            'psfGen', psfGen, 'mccMode', mccMode', 'configFile', configFile, ...
            'GPUConfigFile', GPUConfigFile);
        return;
    end

    % to do: put the code below to a function as in memory computing
    if strcmpi(largeMethod, 'inmemory')
        RLdecon_large_in_memory(frameFullpath, PSFfile, deconFullPath, xyPixelSize, dz, ...
            'save16bit', save16bit, 'Deskew', Deskew, 'SkewAngle', SkewAngle, ...
            'flipZstack', flipZstack, 'Background', Background, 'dzPSF', dzPSF, ...
            'DeconIter', DeconIter, 'RLMethod', RLMethod, 'skewed', skewed, ...
            'wienerAlpha', wienerAlpha, 'OTFCumThresh', OTFCumThresh, 'hannWinBounds', hannWinBounds, ...
            'EdgeErosion', EdgeErosion, 'fixIter', fixIter,'batchSize', batchSize, ...
            'saveZarr', saveZarr, 'dampFactor', dampFactor, 'scaleFactor', scaleFactor, ...
            'deconOffset', deconOffset, 'maskFullpaths', maskFullpaths, 'useGPU', GPUJob, ...
            'uuid', uuid, 'debug', debug, 'psfGen', psfGen);
        return;
    end
end

% toc

end
