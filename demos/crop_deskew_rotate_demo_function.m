function [] = crop_deskew_rotate_crop_demo_function(inputFullpath, outputFullpath, skewAngle, xyPixelSize, dz, inBbox)
% the function needs the path for the input and output, along with
% parameters. The function needs to first check if input and output exist,
% it should skip the output file if it exists.

if ~exist(inputFullpath, 'file')
    error('The input file %s does not exist!', inputFullpath);
end

if exist(outputFullpath, 'file')
    fprintf('The output file %s already exist, skip it!\n', outputFullpath);
end

% read input file
im = readtiff(inputFullpath);

% crop the image
im = crop3d_mex(im, inBbox);

% deskew/rotate for the cropped image
% these are default parameters for the demo images
Reverse = true;
ObjectiveScan = false;
resample = [];
Interp = 'linear';

dsr = deskewRotateFrame3D(single(im), skewAngle, dz, xyPixelSize, 'reverse', Reverse, ...
    'Crop', true, 'ObjectiveScan', ObjectiveScan, 'resample', resample, 'Interp', Interp);

% convert image to uint16
dsr = uint16(dsr);

% write the image to disk with the given file name
% also to ensure the result file is complete, we first write to an
% temporary file, and then rename the temporary file to the final output file

% get uuid
uuid = get_uuid();
tmpFn = sprintf('%s_%s.tif', outputFullpath(1 : end - 4));

% write to temporary file
writetiff(dsr, tmpFn);

% rename to the final output
movefile(tmpFn, outputFullpath);

end
