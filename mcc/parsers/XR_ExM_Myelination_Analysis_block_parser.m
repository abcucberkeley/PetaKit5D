function [done_flag] = XR_ExM_Myelination_Analysis_block_parser(batchInds, axonFullpath, myelinFullpath, ...
    maskFullpath, outFullpath, BatchBBoxes, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('axonFullpath', @(x) ischar(x));
ip.addRequired('myelinFullpath', @(x) ischar(x));
ip.addRequired('maskFullpath', @(x) ischar(x));
ip.addRequired('outFullpath', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('regionInds', [1], @(x) isnumeric(x) || ischar(x)); % #region
ip.addParameter('xyPixelSize', 0.108, @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('dz', 0.2, @(x) isnumeric(x) || ischar(x)); 
ip.addParameter('minDistThresh', 1, @(x) isnumeric(x) || ischar(x)); % um
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, axonFullpath, myelinFullpath, maskFullpath, outFullpath, BatchBBoxes, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
regionInds = pr.regionInds;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
minDistThresh = pr.minDistThresh;
uuid = pr.uuid;
debug = pr.debug;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(BatchBBoxes)
    BatchBBoxes = str2num(BatchBBoxes);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(regionInds)
    regionInds = str2num(regionInds);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(minDistThresh)
    minDistThresh = str2num(minDistThresh);
end
if ischar(debug)
    debug = str2num(debug);
end

XR_ExM_Myelination_Analysis_block(batchInds, axonFullpath, myelinFullpath, ...
    maskFullpath, outFullpath, BatchBBoxes, Overwrite=Overwrite, regionInds=regionInds, ...
    xyPixelSize=xyPixelSize, dz=dz, minDistThresh=minDistThresh, uuid=uuid, ...
    debug=debug);

end

