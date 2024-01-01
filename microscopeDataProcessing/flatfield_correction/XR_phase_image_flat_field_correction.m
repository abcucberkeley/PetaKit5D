function [] = XR_phase_image_flat_field_correction(fns, outPath, const)
% preprocessing for phase images for multiple tiles for a time point
%
% 1. subtract gaussian blurred image (sigma=100) from the raw data. 
% 2. add 1100 counts to the image from step 1 to make all pixels positive (this number may change when running all time points). 
% 3. use the BaSiC software to estimate flatfield using all 25 tiles from step 2, and apply the flatfield correction to the tiles. 
%
% Author: Xiongtao Ruan (08/22/2022)
%
% rename from prcs_janelia_20220428_phase_image.m to XR_phase_image_flat_field_correction.m


if nargin < 3
    const = 1100;
end

[~, fsns] = fileparts(fns);
if ~iscell(fsns)
    fsns = {fsns};
end

fnouts = cellfun(@(x) [outPath, x, '.tif'], fsns, 'unif', 0);
uuid = get_uuid();
% tmpouts = cellfun(@(x) [outPath, x, '_', uuid, '.tif'], fsns, 'unif', 0);
tmpouts = cellfun(@(x) ['/dev/shm/', x, '_', uuid, '.tif'], fsns, 'unif', 0);

nF = numel(fns);
is_done_flag = false(nF, 1);
for f = 1 : nF
    if exist(fnouts{f}, 'file')
        is_done_flag(f) = true;
    end
end

if all(is_done_flag)
    fprintf('All output files exist, skip it!\n')
    return;
end

% read images and subtract gaussian blurred image
ims = cell(nF, 1);

for f = 1 : nF
    imf = single(parallelReadTiff(fns{f}));
    g = imgaussfilt(imf, 100);

    imf = imf - g + const;

    ims{f} = imf;
end

im = cat(3, ims{:});
minval = min(im(:));

% flat field estimation and correction
[flatfield, darkfield] = BaSiC(im,'darkfield','false');

for f = 1 : nF
    imf = ims{f} ./ flatfield;

    writetiff(uint16(imf), tmpouts{f});
    movefile(tmpouts{f}, fnouts{f});
end

% save the minimum vals and flat field image in a mat file
ffPath = [outPath, 'BaSiC_flatfield/'];
if ~exist(ffPath, 'dir')
    mkdir(ffPath);
end

ffFn = sprintf('%s/%s_flatfield.mat', ffPath, fsns{1});
save('-v7.3', ffFn, 'minval', 'flatfield');



end

