function [] = crop_deskew_rotate_demo_function(inputFullpath, outputFullpath, skewAngle, xyPixelSize, dz, inBbox)
% the function needs the path for the input and output, along with
% parameters. The function needs to first check if input and output exist,
% it should skip the output file if it exists.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFullpath', @ischar);
ip.addRequired('outputFullpath', @ischar);
ip.addRequired('skewAngle', @isscalar);
ip.addRequired('xyPixelSize', @isscalar);
ip.addRequired('dz', @isscalar);
ip.addRequired('inBbox', @isvector);

ip.parse(inputFullpath, outputFullpath, skewAngle, xyPixelSize, dz, inBbox);

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
% direct input and output as 16 bit while actual deskew/rotation in single
save16bit = true;

dsr = deskewRotateFrame3D(im, skewAngle, dz, xyPixelSize, 'reverse', Reverse, ...
    'Crop', true, 'ObjectiveScan', ObjectiveScan, 'resample', resample, 'Interp', Interp, ...
    'save16bit', save16bit);

% write the image to disk with the given file nameedit demos/demo_RL_deconvolution.m
% also to ensure the result file is complete, we first write to an
% temporary file, and then rename the temporary file to the final output file

% get uuid
uuid = get_uuid();
tmpFn = sprintf('%s_%s.tif', outputFullpath(1 : end - 4), uuid);

% write to temporary file
writetiff(dsr, tmpFn);

% rename to the final output
movefile(tmpFn, outputFullpath);

end
