function [memAllocate, arrayMem, outSize] = XR_estimateDeskewRotateMemory(imSize, xyPixelSize, dz, varargin)
% estimate the peak memory (in GB) of a deskew/rotate task, i.e., one
% XR_deskewRotateFrame call or one batch of XR_deskewRotateBlock, as the
% largest set of arrays the code path holds at the same time plus the memory
% of MATLAB (or the MATLAB Runtime) itself, with a margin
%
% The arrays per step follow XR_deskewRotateFrame, deskewRotateFrame3D and
% its mex files:
%   combined, linear, sample scan, no resampling (mex): input + interpolated
%     skewed-space slices (when the x step is above 2 pixels) + output
%   imwarp (objective scan, resampling, cubic): input + interpolated volume
%     + output (+ about 80 bytes per output voxel for cubic)
%   resampling of over 1000 slices: the input is binned and resized first
%   separate deskew and rotation: single input, deskewed and rotated
%     volumes and their 16-bit copies; over 1000 slices, deskewed blocks and
%     the DS write from them
%   reading partial files (concatenation), cropping, flat field, permuting
%     and zarr writes (buffers per thread)
% The estimate was at least the measured peak for every path and size
% measured (1800 x 512 x 100 to 5000 demo-based tiles, MATLAB and compiled
% tasks, 1 to 48 threads).
%
% Required inputs:
%              imSize : 1x3 vector. Size (yxz) of the volume deskewed/rotated by the task: all partial files of a group together and after inputBbox; for large files, one input batch.
%         xyPixelSize : a number. Pixel size in um.
%                  dz : a number. Scan interval in um.
%
% Parameters (as 'specifier'-value pairs):
%           skewAngle : a number (default: 32.45). Skew angle (in degree) of the stage.
%       objectiveScan : true|false (default: false). Objective scan.
%          zStageScan : true|false (default: false). Z stage scan (orthogonal to objective scan).
%         DSRCombined : true|false (default: true). Use combined processing for deskew and rotation.
%              rotate : true|false (default: true). Rotate after deskew.
%           save16bit : true|false (default: true). Save 16bit result.
%        interpMethod : 'linear'|'cubic' (default: 'linear'). Interpolation method.
%      resampleFactor : empty or 1x1, 1x2 or 1x3 vector (default: []). Resampling factor after rotation.
%            readSize : empty or 1x3 vector (default: []). Size read from disk before the crop by inputBbox. Empty: imSize.
%         partialFile : true|false (default: false). The input is read from several partial files and concatenated.
%          preprocess : true|false (default: false). Flat-field correction, flip or input axis permutation copy the input.
%             outSize : empty or 1x3 vector (default: []). Output size (yxz), e.g., of one batch in large-file mode. Empty: computed as in deskewRotateFrame3D.
%     outputAxisOrder : char (default: 'yxz'). Output axis order (other orders permute a copy of the output).
%            saveZarr : true|false (default: false). Save the output as zarr.
%           blockSize : 1x3 vector (default: [256, 256, 256]). Chunk size of the zarr output.
%           largeFile : true|false (default: false). One batch of large-file deskew/rotation (XR_deskewRotateBlock).
%            dataType : char (default: 'uint16'). Data type of the input.
%
% Outputs:
%         memAllocate : estimated peak memory of the task in GB.
%            arrayMem : the arrays held at the peak in GB (without MATLAB itself).
%             outSize : output size (yxz) of the deskewed/rotated volume.

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imSize', @(x) isnumeric(x) && numel(x) == 3);
ip.addRequired('xyPixelSize', @isscalar);
ip.addRequired('dz', @isscalar);
ip.addParameter('skewAngle', 32.45, @isscalar);
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('zStageScan', false, @islogical);
ip.addParameter('DSRCombined', true, @islogical);
ip.addParameter('rotate', true, @islogical);
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('interpMethod', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})));
ip.addParameter('resampleFactor', [], @(x) isempty(x) || isnumeric(x));
ip.addParameter('readSize', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 3));
ip.addParameter('partialFile', false, @islogical);
ip.addParameter('preprocess', false, @islogical);
ip.addParameter('outSize', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 3));
ip.addParameter('outputAxisOrder', 'yxz', @ischar);
ip.addParameter('saveZarr', false, @islogical);
ip.addParameter('blockSize', [256, 256, 256], @isnumeric);
ip.addParameter('largeFile', false, @islogical);
ip.addParameter('dataType', 'uint16', @(x) ischar(x) || isstring(x));
ip.parse(imSize, xyPixelSize, dz, varargin{:});

pr = ip.Results;
skewAngle = pr.skewAngle;
objectiveScan = pr.objectiveScan;
zStageScan = pr.zStageScan;
DSRCombined = pr.DSRCombined;
rotate = pr.rotate;
save16bit = pr.save16bit;
interpMethod = pr.interpMethod;
resampleFactor = pr.resampleFactor;
readSize = pr.readSize;
partialFile = pr.partialFile;
preprocess = pr.preprocess;
outSize = pr.outSize;
outputAxisOrder = pr.outputAxisOrder;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
largeFile = pr.largeFile;
dataType = char(pr.dataType);

% memory of MATLAB itself (measured up to 1.4 GB in a -nojvm task, 0.7 GB
% in a compiled task) and the margin on the arrays
baseMem = 2;
memFactor = 1.1;
% imwarp's own arrays per output voxel (bytes): none measurable for linear,
% about 80 for cubic
imwarpBytes = 0;
if strcmpi(interpMethod, 'cubic')
    imwarpBytes = 85;
end
% zarr writes buffer up to two chunks per thread (measured: about one chunk
% per thread), for up to 48 threads
zarrThreadNum = 48;

imSize = imSize(:)';
ny = imSize(1);
nx = imSize(2);
nz = imSize(3);
nIn = prod(imSize);
if isempty(readSize)
    readSize = imSize;
end
nRead = prod(readSize(:)');
byteNum = dataTypeToByteNumber(dataType);

% the angle as in XR_deskewRotateFrame
if zStageScan
    skewAngle = 90 - skewAngle;
end
theta = skewAngle * pi / 180;

rs = [];
if ~isempty(resampleFactor)
    rs = resampleFactor(:)';
    rs = [ones(1, 4 - numel(rs)) * rs(1), rs(2 : end)];
    if all(rs == 1)
        rs = [];
    end
end

% read (concatenation of partial files, crop) and preprocessing
arrays = byteNum * nRead * (1 + partialFile);
if nRead > nIn
    arrays(end + 1) = byteNum * (nRead + nIn);
end
if preprocess
    arrays(end + 1) = 2 * byteNum * nIn;
end

% XR_deskewRotateFrame bins (z) and resizes (xy) inputs of over 1000 slices
% before a combined deskew/rotation with resampling: a binned copy, then
% imresize (one dimension at a time), and the rest runs on the smaller volume
if DSRCombined && rotate && ~largeFile && ~isempty(rs) && nz > 1000
    outPixelSize = rs * xyPixelSize;
    preRs = min(outPixelSize) ./ [xyPixelSize, xyPixelSize, dz * sind(skewAngle)];
    preRs(3) = max(1, round(preRs(3)));
    if any(preRs ~= 1)
        nz = numel(round(preRs(3) / 2) : preRs(3) : nz);
        nBin = ny * nx * nz;
        ny_rs = round(ny / preRs(1));
        nx_rs = round(nx / preRs(2));
        arrays(end + 1) = byteNum * (nIn + nBin);
        arrays(end + 1) = byteNum * (nBin + nBin / min(preRs(1 : 2)) + ny_rs * nx_rs * nz);
        ny = ny_rs;
        nx = nx_rs;
        nIn = ny * nx * nz;
        xyPixelSize = xyPixelSize * preRs(1);
        dz = dz * preRs(3);
        rs = outPixelSize ./ min(outPixelSize);
        if all(rs == 1)
            rs = [];
        end
    end
end

% the x step as in deskewRotateFrame3D
dx = cos(theta) * dz / xyPixelSize;

if isempty(outSize)
    if ~objectiveScan
        zAniso = sin(abs(theta)) * dz / xyPixelSize;
        outSize = round([ny, (nx - 1) * cos(theta) + (nz - 1) * zAniso / sin(abs(theta)), (nx - 1) * sin(abs(theta)) - 4]);
    else
        zAniso = dz / xyPixelSize;
        outSize = round([ny, nx * cos(theta) + nz * zAniso * sin(abs(theta)), nz * zAniso * cos(theta) + nx * sin(abs(theta))]);
    end
    if ~isempty(rs)
        outSize = round(outSize ./ rs([2, 1, 3]));
    end
end
nOut = prod(max(outSize, 0));

if largeFile && ~save16bit
    % XR_deskewRotateBlock converts the batch to single
    arrays(end + 1) = (byteNum + 4) * nIn;
    byteNum = 4;
end

if ~DSRCombined || ~rotate
    % separate deskew (imwarp in single) and rotation (imwarp in single)
    nxDS = nx + ceil((nz - 1) * abs(dx));
    nDS = ny * nxDS * nz;
    arrays(end + 1) = (byteNum + 4) * nIn;
    if nz <= 1000
        arrays(end + 1) = 4 * nIn + (4 + imwarpBytes) * nDS;
        arrays(end + 1) = (4 + 2 * save16bit) * nDS;
    else
        % over 1000 slices, blocks of up to 150 rows are deskewed to disk
        % (blockedImage), DS is written from those blocks, and the rotation
        % reads DS back in single
        blockRows = min(ceil(ny / 8), 150) + 10;
        dsByteNum = 4 - 2 * save16bit;
        arrays(end + 1) = 8 * nIn + 3 * (4 + imwarpBytes) * blockRows * nxDS * nz;
        arrays(end + 1) = dsByteNum * (nDS + 3 * ny * nxDS * min(nz, 100));
        if rotate
            arrays(end + 1) = (4 + 2 * save16bit) * nDS;
        end
    end
    if rotate
        arrays(end + 1) = 4 * nDS + (4 + imwarpBytes) * nOut;
        arrays(end + 1) = (4 + 2 * save16bit) * nOut;
    end
else
    % skewed-space interpolation when the x step is above 2 pixels: k - 1 new
    % slices between neighboring slices
    nInterp = 1;
    if ~objectiveScan && abs(dx) > 2
        nInterp = ceil(abs(dx) / 2);
    end
    useMex = ~objectiveScan && isempty(rs) && strcmpi(interpMethod, 'linear');
    if save16bit && byteNum == 2
        intByteNum = 2;
    else
        intByteNum = 4;
    end
    if useMex
        nInt = ny * nx * (nz - 1) * (nInterp - 1);
        arrays(end + 1) = byteNum * nIn + intByteNum * (nInt + nOut);
        outByteNum = intByteNum;
    else
        % the interpolated volume includes the original slices; imwarp keeps
        % the class of its input
        nInt = 0;
        outByteNum = byteNum;
        if nInterp > 1
            nInt = ny * nx * ((nz - 1) * nInterp + 1);
            outByteNum = intByteNum;
        end
        arrays(end + 1) = byteNum * nIn + intByteNum * nInt + (outByteNum + imwarpBytes) * nOut;
        if ~save16bit && outByteNum < 4
            % single(dsr) copy of a 16-bit imwarp result
            arrays(end + 1) = (outByteNum + 4) * nOut;
            outByteNum = 4;
        end
    end
    if save16bit && outByteNum == 4
        % uint16(dsr) copy of a single result (single input)
        arrays(end + 1) = (outByteNum + 2) * nOut;
    end
    if ~strcmp(outputAxisOrder, 'yxz')
        arrays(end + 1) = 2 * outByteNum * nOut;
    end
    if saveZarr || largeFile
        chunkSize = min(blockSize(:)', max(outSize, 1));
        chunkNum = prod(ceil(max(outSize, 1) ./ chunkSize));
        arrays(end + 1) = outByteNum * nOut + 2 * min(zarrThreadNum, chunkNum) * prod(chunkSize) * outByteNum;
    end
end

arrayMem = max(arrays) / 2^30;
memAllocate = baseMem + memFactor * arrayMem;

end
