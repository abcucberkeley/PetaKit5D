function [done_flag] = XR_deskewRotateBlock_parser(batchInds, zarrFullpath, dsrFullpath, flagFullname, ...
    inBatchBBoxes, outBatchBBoxes, outRegionBBoxes, outLocalBboxes, xyPixelSize, dz, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('dsrFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('inBatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('outBatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('outRegionBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('outLocalBboxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('xyPixelSize', @(x) isnumeric(x) || ischar(x)); %in um
ip.addRequired('dz', @(x) isnumeric(x) || ischar(x)); %in um
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('SkewAngle', 32.45 , @(x) isnumeric(x) || ischar(x));
ip.addParameter('Reverse', true, @(x) islogical(x) || ischar(x)); 
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('interpMethod', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})) || ischar(x));
ip.addParameter('resampleFactor', [], @(x) isempty(x) || isnumeric(x) || ischar(x)); % resampling after rotation 
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, dsrFullpath, flagFullname, inBatchBBoxes, ...
    outBatchBBoxes, outRegionBBoxes, outLocalBboxes, xyPixelSize, dz, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
flipZstack = pr.flipZstack;
interpMethod = pr.interpMethod;
resampleFactor = pr.resampleFactor;
uuid = pr.uuid;
debug = pr.debug;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(inBatchBBoxes)
    inBatchBBoxes = str2num(inBatchBBoxes);
end
if ischar(outBatchBBoxes)
    outBatchBBoxes = str2num(outBatchBBoxes);
end
if ischar(outRegionBBoxes)
    outRegionBBoxes = str2num(outRegionBBoxes);
end
if ischar(outLocalBboxes)
    outLocalBboxes = str2num(outLocalBboxes);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(Reverse)
    Reverse = str2num(Reverse);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(resampleFactor)
    resampleFactor = str2num(resampleFactor);
end
if ischar(debug)
    debug = str2num(debug);
end

XR_deskewRotateBlock(batchInds, zarrFullpath, dsrFullpath, flagFullname, inBatchBBoxes, ...
    outBatchBBoxes, outRegionBBoxes, outLocalBboxes, xyPixelSize, dz, Overwrite=Overwrite, ...
    SkewAngle=SkewAngle, Reverse=Reverse, flipZstack=flipZstack, interpMethod=interpMethod, ...
    resampleFactor=resampleFactor, uuid=uuid, debug=debug);

end

