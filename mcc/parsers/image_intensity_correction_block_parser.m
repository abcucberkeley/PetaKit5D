function [done_flag] = image_intensity_correction_block_parser(batchInds, zarrFullpath, outFullpath, IntCorrFn, flagFullname, inBatchBBoxes, outBatchBBoxes, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('outFullpath', @(x) ischar(x));
ip.addRequired('IntCorrFn', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('inBatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('outBatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addParameter('constOffset', 0, @(x) isnumeric(x) || ischar(x));
ip.addParameter('constFactor', 1, @(x) isnumeric(x) || ischar(x));
ip.addParameter('lowerLimit', 0.4, @(x) isnumeric(x) || ischar(x));
ip.addParameter('upperLimit', 3, @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, outFullpath, IntCorrFn, flagFullname, inBatchBBoxes, outBatchBBoxes, varargin{:});

pr = ip.Results;
constOffset = pr.constOffset;
constFactor = pr.constFactor;
lowerLimit = pr.lowerLimit;
upperLimit = pr.upperLimit;
Overwrite = pr.Overwrite;

if ischar(blockInds)
    blockInds = str2num(blockInds);
end
if ischar(inBatchBBoxes)
    inBatchBBoxes = str2num(inBatchBBoxes);
end
if ischar(outBatchBBoxes)
    outBatchBBoxes = str2num(outBatchBBoxes);
end
if ischar(constOffset)
    constOffset = str2num(constOffset);
end
if ischar(constFactor)
    constFactor = str2num(constFactor);
end
if ischar(lowerLimit)
    lowerLimit = str2num(lowerLimit);
end
if ischar(upperLimit)
    upperLimit = str2num(upperLimit);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end

image_intensity_correction_block(blockInds, zarrFullpath, outFullpath, IntCorrFn, ...
    flagFullname, inBatchBBoxes, outBatchBBoxes, constOffset=constOffset, constFactor=constFactor, ...
    lowerLimit=lowerLimit, upperLimit=upperLimit, Overwrite=Overwrite);

end

