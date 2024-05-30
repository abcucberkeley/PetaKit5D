function [done_flag] = XR_ExM_Segmentation_interpolated_background_block_parser(batchInds, zarrFullpath, outFullpath, ...
    bgFullpath, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('outFullpath', @(x) ischar(x));
ip.addRequired('bgFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('RegionBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('localBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('bgFactor', 1.0, @(x) isnumeric(x) || ischar(x));
ip.addParameter('Sigma', 2.5, @(x) isnumeric(x) || ischar(x)); % for 3D Gauss filtering
ip.addParameter('volThrsh', 1000, @(x) isvector(x) || ischar(x)); % set to zere if no "pre-cleaning" is necessary to remove high-freq noise; second value is for the DAN channel
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, outFullpath, bgFullpath, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
bgFactor = pr.bgFactor;
Sigma = pr.Sigma;
volThrsh = pr.volThrsh;
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
if ischar(bgFactor)
    bgFactor = str2num(bgFactor);
end
if ischar(Sigma)
    Sigma = str2num(Sigma);
end
if ischar(volThrsh)
    volThrsh = str2num(volThrsh);
end
if ischar(debug)
    debug = str2num(debug);
end

XR_ExM_Segmentation_interpolated_background_block(batchInds, zarrFullpath, ...
    outFullpath, bgFullpath, flagFullname, BatchBBoxes, RegionBBoxes, localBBoxes, ...
    Overwrite=Overwrite, bgFactor=bgFactor, Sigma=Sigma, volThrsh=volThrsh, ...
    uuid=uuid, debug=debug);

end

