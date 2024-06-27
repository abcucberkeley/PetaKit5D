function [] = RLdecon_for_zarr_block_parser(batchInds, zarrFullpath, psfFullpath, deconFullpath, ...
    flagFullname, BatchBBoxes, RegionBBoxes, xyPixelSize, dz, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('psfFullpath', @(x) ischar(x));
ip.addRequired('deconFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('RegionBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('xyPixelSize', @(x) isnumeric(x) || ischar(x)); %in um
ip.addRequired('dz', @(x) isnumeric(x) || ischar(x)); %in um
ip.addParameter('save16bit', true , @(x) islogical(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('SkewAngle', -32.45 , @(x) isnumeric(x) || ischar(x));
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('Background', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('dzPSF', 0.1 , @(x) isnumeric(x) || ischar(x)); %in um
ip.addParameter('DeconIter', 15 , @(x) isnumeric(x) || ischar(x)); % number of iterations
ip.addParameter('dampFactor', 1, @(x) isnumeric(x) || ischar(x)); % damp factor for decon result
ip.addParameter('scaleFactor', 1.0, @(x) isnumeric(x) || ischar(x)); % scale factor for decon result
ip.addParameter('deconOffset', 0, @(x) isnumeric(x) || ischar(x)); % offset for decon result
ip.addParameter('EdgeErosion', 0, @(x) isnumeric(x) || ischar(x)); % edge erosion for decon result
ip.addParameter('maskFullpaths', {} , @(x) iscell(x) || ischar(x)); % Full paths of 2D mask zarr files, in xy, xz, yz order
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @(x) isnumeric(x) || ischar(x));
ip.addParameter('OTFCumThresh', 0.9, @(x) isnumeric(x) || ischar(x)); % OTF cumutative sum threshold
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x) || ischar(x)); % decon in skewed space
ip.addParameter('fixIter', false, @(x) islogical(x) || ischar(x)); % CPU Memory in Gb
ip.addParameter('useGPU', false, @(x) islogical(x) || ischar(x)); % use gpu for chuck deconvolution. 
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));
ip.addParameter('psfGen', true, @(x) islogical(x) || ischar(x)); % psf generation

ip.parse(batchInds, zarrFullpath, psfFullpath, deconFullpath, flagFullname, ...
    BatchBBoxes, RegionBBoxes, xyPixelSize, dz, varargin{:});

pr = ip.Results;
save16bit = pr.save16bit;
Overwrite = pr.Overwrite;
SkewAngle = pr.SkewAngle;
flipZstack = pr.flipZstack;
Background = pr.Background;
dzPSF = pr.dzPSF;
DeconIter = pr.DeconIter;
dampFactor = pr.dampFactor;
scaleFactor = pr.scaleFactor;
deconOffset = pr.deconOffset;
EdgeErosion = pr.EdgeErosion;
maskFullpaths = pr.maskFullpaths;
RLMethod = pr.RLMethod;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
skewed = pr.skewed;
fixIter = pr.fixIter;
useGPU = pr.useGPU;
uuid = pr.uuid;
debug = pr.debug;
psfGen = pr.psfGen;

if ischar(batchInds)
    batchInds = str2num(batchInds);
end
if ischar(BatchBBoxes)
    BatchBBoxes = str2num(BatchBBoxes);
end
if ischar(RegionBBoxes)
    RegionBBoxes = str2num(RegionBBoxes);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(save16bit)
    save16bit = str2num(save16bit);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(flipZstack)
    flipZstack = str2num(flipZstack);
end
if ischar(Background)
    Background = str2num(Background);
end
if ischar(dzPSF)
    dzPSF = str2num(dzPSF);
end
if ischar(DeconIter)
    DeconIter = str2num(DeconIter);
end
if ischar(dampFactor)
    dampFactor = str2num(dampFactor);
end
if ischar(scaleFactor)
    scaleFactor = str2num(scaleFactor);
end
if ischar(deconOffset)
    deconOffset = str2num(deconOffset);
end
if ischar(EdgeErosion)
    EdgeErosion = str2num(EdgeErosion);
end
if ischar(maskFullpaths) && ~isempty(maskFullpaths) && strcmp(maskFullpaths(1), '{')
    maskFullpaths = eval(maskFullpaths);
end
if ischar(wienerAlpha)
    wienerAlpha = str2num(wienerAlpha);
end
if ischar(OTFCumThresh)
    OTFCumThresh = str2num(OTFCumThresh);
end
if ischar(skewed)
    skewed = str2num(skewed);
end
if ischar(fixIter)
    fixIter = str2num(fixIter);
end
if ischar(useGPU)
    useGPU = str2num(useGPU);
end
if ischar(debug)
    debug = str2num(debug);
end
if ischar(psfGen)
    psfGen = str2num(psfGen);
end

RLdecon_for_zarr_block(batchInds, zarrFullpath, psfFullpath, deconFullpath, ...
    flagFullname, BatchBBoxes, RegionBBoxes, xyPixelSize, dz, save16bit=save16bit, ...
    Overwrite=Overwrite, SkewAngle=SkewAngle, flipZstack=flipZstack, Background=Background, ...
    dzPSF=dzPSF, DeconIter=DeconIter, dampFactor=dampFactor, scaleFactor=scaleFactor, ...
    deconOffset=deconOffset, EdgeErosion=EdgeErosion, maskFullpaths=maskFullpaths, ...
    RLMethod=RLMethod, wienerAlpha=wienerAlpha, OTFCumThresh=OTFCumThresh, ...
    skewed=skewed, fixIter=fixIter, useGPU=useGPU, uuid=uuid, debug=debug, ...
    psfGen=psfGen);

end

