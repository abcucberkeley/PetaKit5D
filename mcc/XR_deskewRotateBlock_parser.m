function XR_deskewRotateBlock_parser(batchInds, zarrFullpath, dsrFullpath, flagFullname, ...
    BatchBBoxes, RegionBBoxes, borderSize, xyPixelSize, dz, varargin)
% Deskew and/or rotate data for given blocks

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('blockInds', @(x) isnumeric(x) || ischar(x));
ip.addRequired('zarrFullpath', @(x) ischar(x));
ip.addRequired('dsrFullpath', @(x) ischar(x));
ip.addRequired('flagFullname', @(x) ischar(x));
% ip.addParameter('ResultDir', 'matlab_stitch', @ischar);
ip.addRequired('BatchBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('RegionBBoxes', @(x) isnumeric(x) || ischar(x));
ip.addRequired('borderSize', @(x) isnumeric(x) || ischar(x));
ip.addRequired('pixelSize', @(x) isnumeric(x) || ischar(x)); %in um
ip.addRequired('dz', @(x) isnumeric(x) || ischar(x)); %in um
% ip.addParameter('BlockSize', [], @isnumeric);
ip.addParameter('Overwrite', false, @(x) islogical(x) || ischar(x));
ip.addParameter('SkewAngle', 32.45 , @(x) isnumeric(x) || ischar(x));
ip.addParameter('Reverse', true, @(x) islogical(x) || ischar(x)); 
ip.addParameter('flipZstack', false, @(x) islogical(x) || ischar(x)); 
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear'})) && ischar(x));
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));

ip.parse(batchInds, zarrFullpath, dsrFullpath, flagFullname, ...
    BatchBBoxes, RegionBBoxes, borderSize, xyPixelSize, dz, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
flipZstack = pr.flipZstack;
Interp = pr.Interp;
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
if ischar(borderSize)
    borderSize = str2num(borderSize);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(Overwrite)
    Overwrite = strcmp(Overwrite,'true');
end
if ischar(SkewAngle)
    SkewAngle = str2num(SkewAngle);
end
if ischar(Reverse)
    Reverse = strcmp(Reverse,'true');
end
if ischar(flipZstack)
    flipZstack = strcmp(flipZstack,'true');
end
if ischar(debug)
    debug = strcmp(debug,'true');
end

XR_deskewRotateBlock(batchInds, zarrFullpath, dsrFullpath, flagFullname, ...
    BatchBBoxes, RegionBBoxes, borderSize, xyPixelSize, dz,'Overwrite',Overwrite,...
    'SkewAngle',SkewAngle,'Reverse',Reverse,'flipZstack',flipZstack,'Interp',Interp,...
    'uuid',uuid,'debug',debug)
