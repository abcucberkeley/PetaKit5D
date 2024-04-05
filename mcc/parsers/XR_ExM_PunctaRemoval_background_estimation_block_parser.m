function [done_flag] = XR_ExM_PunctaRemoval_background_estimation_block_parser(batchInds, zarrFullpath, ...
    resultFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('resultFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('RegionBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('localBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('Sigma', 2.5, @(x) isnumeric(x) || ischar(x)); % for 3D Gauss filtering
ip.addParameter('OTSUMaxPer', 99.9, @(x) isnumeric(x) || ischar(x)); % Max percentile of data for OTSU calculation
ip.addParameter('MinThreshold', 0, @(x) isnumeric(x) || ischar(x)); % if left empty, will use OTSU to calculate intensity to threshold;
ip.addParameter('MaxThreshold', [], @(x) isnumeric(x) || ischar(x)); % if left empty, will use OTSU to calculate intensity to threshold;
ip.addParameter('BaseThreshold', [], @(x) isnumeric(x) || ischar(x)); % if left empty, will use OTSU to calculate intensity to threshold;
ip.addParameter('volThrsh', 1000, @(x) isvector(x) || ischar(x)); % set to zere if no "pre-cleaning" is necessary to remove high-freq noise; second value is for the DAN channel
ip.addParameter('offset', 0, @(x) isnumeric(x) || ischar(x)); % offset to add to the cleaned image to make the background non-zero
ip.addParameter('localWinSize', [15, 15, 23], @(x) isvector(x) || ischar(x)); % local window size for cropping local region
ip.addParameter('SigmaThrsh', 4, @(x) isnumeric(x) || ischar(x)); % sigma threshold for point detection
ip.addParameter('intThrsh', 10000, @(x) isnumeric(x) || ischar(x)); % intensity threshold for the peak to be removed
ip.addParameter('initDetect', '2d', @ischar); % initial detection method, 2d mip or 3d stack
ip.addParameter('detVolThrsh', 5000, @(x) isnumeric(x) || ischar(x)); % volume threshold for the peak to be removed
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, resultFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
Sigma = pr.Sigma;
OTSUMaxPer = pr.OTSUMaxPer;
MinThreshold = pr.MinThreshold;
MaxThreshold = pr.MaxThreshold;
BaseThreshold = pr.BaseThreshold;
volThrsh = pr.volThrsh;
offset = pr.offset;
localWinSize = pr.localWinSize;
SigmaThrsh = pr.SigmaThrsh;
intThrsh = pr.intThrsh;
initDetect = pr.initDetect;
detVolThrsh = pr.detVolThrsh;
uuid = pr.uuid;
debug = pr.debug;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(BatchBBoxes)
    BatchBBoxes = str2num(BatchBBoxes);
end
if ischar(RegionBBoxes)
    RegionBBoxes = str2num(RegionBBoxes);
end
if ischar(localBBoxes)
    localBBoxes = str2num(localBBoxes);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(Sigma)
    Sigma = str2num(Sigma);
end
if ischar(OTSUMaxPer)
    OTSUMaxPer = str2num(OTSUMaxPer);
end
if ischar(MinThreshold)
    MinThreshold = str2num(MinThreshold);
end
if ischar(MaxThreshold)
    MaxThreshold = str2num(MaxThreshold);
end
if ischar(BaseThreshold)
    BaseThreshold = str2num(BaseThreshold);
end
if ischar(volThrsh)
    volThrsh = str2num(volThrsh);
end
if ischar(offset)
    offset = str2num(offset);
end
if ischar(localWinSize)
    localWinSize = str2num(localWinSize);
end
if ischar(SigmaThrsh)
    SigmaThrsh = str2num(SigmaThrsh);
end
if ischar(intThrsh)
    intThrsh = str2num(intThrsh);
end
if ischar(detVolThrsh)
    detVolThrsh = str2num(detVolThrsh);
end
if ischar(debug)
    debug = str2num(debug);
end

XR_ExM_PunctaRemoval_background_estimation_block(batchInds, zarrFullpath, resultFullname, ...
    BatchBBoxes, RegionBBoxes, localBBoxes, Overwrite=Overwrite, Sigma=Sigma, ...
    OTSUMaxPer=OTSUMaxPer, MinThreshold=MinThreshold, MaxThreshold=MaxThreshold, ...
    BaseThreshold=BaseThreshold, volThrsh=volThrsh, offset=offset, localWinSize=localWinSize, ...
    SigmaThrsh=SigmaThrsh, intThrsh=intThrsh, initDetect=initDetect, detVolThrsh=detVolThrsh, ...
    uuid=uuid, debug=debug);

end

