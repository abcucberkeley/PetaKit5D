function [] = RLdecon_for_zarr_block_parser(batchInds, zarrFullpath, psfFullpath, deconFullpath, ...
    flagFullname, BatchBBoxes, RegionBBoxes, xyPixelSize, dz, varargin)
% RL decon for given zarr block/batches and write to the certain output
% location. 
% 
% xruan (05/24/2022): add support for masked decon, directly skip empty
% regions if the region is empty in the masks. Masks are 2D projections.

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('batchInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('psfFullpath', @(x) ischar(x));
ip.addRequired('deconFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('RegionBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('pixelSize', @(x) isnumeric(x) || ischar(x)); %in um
ip.addRequired('dz', @(x) isnumeric(x) || ischar(x)); %in um
% ip.addParameter('BlockSize', [], @isnumeric);
ip.addParameter('Save16bit', false , @(x) islogical(x) || ischar(x));
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('SkewAngle', -32.45 , @(x) isnumeric(x) || ischar(x));
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('Background', [], @(x) isnumeric(x) || ischar(x));
ip.addParameter('dzPSF', 0.1 , @(x) isnumeric(x) || ischar(x)); %in um
ip.addParameter('DeconIter', 15 , @(x) isnumeric(x) || ischar(x)); % number of iterations
ip.addParameter('scaleFactor', 1e8 , @(x) isnumeric(x) || ischar(x)); % scale factor for data
ip.addParameter('deconMaskFns', {} , @(x) iscell(x) || ischar(x)); % Full paths of 2D mask zarr files, in xy, xz, yz order
ip.addParameter('RLMethod', 'simplified' , @ischar); % rl method {'original', 'simplified', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @(x) isnumeric(x) || ischar(x));
ip.addParameter('OTFCumThresh', 0.9, @(x) isnumeric(x) || ischar(x)); % OTF cumutative sum threshold
ip.addParameter('skewed', [], @(x) (isempty(x) || islogical(x)) || ischar(x)); % decon in skewed space
ip.addParameter('fixIter', false, @(x) islogical(x) || ischar(x)); % CPU Memory in Gb
ip.addParameter('useGPU', false, @(x) islogical(x) || ischar(x)); % use gpu for chuck deconvolution. 
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));
ip.addParameter('psfGen', true, @(x) islogical(x) || ischar(x)); % psf generation

ip.parse(batchInds, zarrFullpath, psfFullpath, deconFullpath, flagFullname, ...
    BatchBBoxes, RegionBBoxes, xyPixelSize, dz, varargin{:});

pr = ip.Results;
Save16bit = pr.Save16bit;
Overwrite = pr.Overwrite;
SkewAngle = pr.SkewAngle;
flipZstack = pr.flipZstack;
Background = pr.Background;
dzPSF = pr.dzPSF;
DeconIter = pr.DeconIter;
scaleFactor = pr.scaleFactor;
deconMaskFns = pr.deconMaskFns;
RLMethod = pr.RLMethod;
skewed = pr.skewed;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
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
if ischar(Save16bit)
    Save16bit = strcmp(Save16bit,'true');
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite,'true');
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(flipZstack)
    flipZstack = strcmp(flipZstack,'true');
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
if ischar(scaleFactor)
    scaleFactor = str2num(scaleFactor);
end
if ischar(deconMaskFns)
    deconMaskFns = eval(deconMaskFns);
end
if ischar(wienerAlpha)
    wienerAlpha = str2num(wienerAlpha);
end
if ischar(OTFCumThresh)
    OTFCumThresh = str2num(OTFCumThresh);
end
if ischar(skewed)
    skewed = strcmp(skewed,'true');
end
if ischar(fixIter)
    fixIter = strcmp(fixIter,'true');
end
if ischar(useGPU)
    useGPU = strcmp(useGPU,'true');
end
if ischar(debug)
    debug = strcmp(debug,'true');
end
if ischar(psfGen)
    psfGen = strcmp(psfGen,'true');
end

RLdecon_for_zarr_block(batchInds, zarrFullpath, psfFullpath, deconFullpath, ...
    flagFullname, BatchBBoxes, RegionBBoxes, xyPixelSize, dz, 'Save16bit',Save16bit,...
    'Overwrite',Overwrite,'SkewAngle',SkewAngle,'flipZstack',flipZstack,'Background',Background,...
    'dzPSF',dzPSF,'DeconIter',DeconIter,'scaleFactor',scaleFactor,'deconMaskFns',deconMaskFns,...
    'RLMethod',RLMethod,'wienerAlpha',wienerAlpha,'OTFCumThresh',OTFCumThresh,...
    'skewed',skewed,'fixIter',fixIter,'useGPU',useGPU,'uuid',uuid,...
    'debug',debug,'psfGen',psfGen)