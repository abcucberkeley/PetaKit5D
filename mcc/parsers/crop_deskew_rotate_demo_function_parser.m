function [] = crop_deskew_rotate_demo_function_parser(inputFullpath, outputFullpath, skewAngle, xyPixelSize, dz, inBbox)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFullpath', @ischar);
ip.addRequired('outputFullpath', @ischar);
ip.addRequired('skewAngle', @(x) isscalar(x) || ischar(x));
ip.addRequired('xyPixelSize', @(x) isscalar(x) || ischar(x));
ip.addRequired('dz', @(x) isscalar(x) || ischar(x));
ip.addRequired('inBbox', @(x) isvector(x) || ischar(x));

ip.parse(inputFullpath, outputFullpath, skewAngle, xyPixelSize, dz, inBbox);

pr = ip.Results;

if ischar(skewAngle)
    skewAngle = str2num(skewAngle);
end
if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(inBbox)
    inBbox = str2num(inBbox);
end

crop_deskew_rotate_demo_function(inputFullpath, outputFullpath, skewAngle, ...
    xyPixelSize, dz, inBbox);

end

