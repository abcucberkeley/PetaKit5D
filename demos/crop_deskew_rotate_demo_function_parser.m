function [] = crop_deskew_rotate_demo_function_parser(inputFullpath, outputFullpath, skewAngle, xyPixelSize, dz, inBbox)
% parser for crop_deskew_rotate_demo_function.m in mcc mode


%#function crop_deskew_rotate_demo_function

% for mcc mode, all input variables are char/string, we need to convert
% them accordingly for the actual input types.

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

% call the actual function with the correct types of inputs
crop_deskew_rotate_demo_function(inputFullpath, outputFullpath, skewAngle, xyPixelSize, dz, inBbox);

end