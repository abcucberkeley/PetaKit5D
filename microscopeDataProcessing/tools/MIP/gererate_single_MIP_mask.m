function [] = gererate_single_MIP_mask(fn, fnout, intThresh, volThresh, dilateSize, blockSize, save16bit, uuid)
% generate MIP masks based on intensity thresholding and object size
%
% Author: Xiongtao Ruan (08/02/2023)


if nargin < 3
    intThresh = 100;
end

if nargin < 4
    volThresh = 1e5;
end

if nargin < 5
    dilateSize = [0, 0];
end

if nargin < 6
    blockSize = [256, 256, 1];
end

if nargin < 7
    save16bit = true;
end

if nargin < 8
    uuid = get_uuid();
end

if exist(fnout, 'file') || exist(fnout, 'dir')
    fprintf('The output file %s already exists, skip it!\n', fnout)
    return;
end    

if save16bit
    dtype = 'uint16';
else
    dtype = 'single';
end

[~, fsn, ext] = fileparts(fn);

switch ext
    case {'.tif', '.tiff'}
        im = readtiff(fn);
    case '.zarr'
        im = readzarr(fn);
end

sz = size(im);

bw = im > intThresh;

n_se = 5;
bw = imopen(bw, strel('disk', n_se));

bw = bwareaopen(bw, volThresh);

if dilateSize(1) == dilateSize(2)
    bw = imdilate(bw, strel('disk', n_se + dilateSize(1)));
else
    bw = imdilate(bw, strel('disk', n_se));
    if dilateSize(1) > 0
        bw = imdilate(bw, strel('line', dilateSize(1) * 2 + 1, 90));
    end
    if dilateSize(2) > 0
        bw = imdilate(bw, strel('line', dilateSize(2) * 2 + 1, 0));
    end
end

bw = imfill(bw, 'hole');

[outPath, fsn, ext] = fileparts(fnout);

% save bounding box
bbox = getImageBoundingBox(bw);

bboxFnout = sprintf('%s/%s_bbox.mat', outPath, fsn);
bboxTmpout = sprintf('%s/%s_bbox_%s.mat', outPath, fsn, uuid);
save('-v7.3', bboxTmpout, 'fn', 'fnout', 'intThresh', 'volThresh', 'sz', 'bbox');
movefile(bboxTmpout, bboxFnout);


% visualize input and output
if usejava('jvm')
    figure, 
    set(gcf, 'color', 'w', 'position', [1, 1, 1200, 600]);

    subplot(1, 2, 1);
    imagesc(im);
    clim([0, prctile(im(:), 99)]);
    axis equal;
    xlim([0, sz(2)]);
    ylim([0, sz(1)]);

    subplot(1, 2, 2);
    imagesc(bw);
    axis equal;
    xlim([0, sz(2)]);
    ylim([0, sz(1)]);

    figFn = sprintf('%s/%s_plot.png', outPath, fsn);
    print(gcf, figFn, '-dpng', '-r0');
    close(gcf);
end

tmpFnout = sprintf('%s/%s_%s%s', outPath, fsn, uuid, ext);

if save16bit
    bw = uint16(bw);
else
    bw = single(bw);
end

switch ext
    case {'.tif'}
        writetiff(bw, tmpFnout);
    case '.zarr'
        createzarr(tmpFnout, dataSize=size(im, 1:3), blockSize=blockSize, dtype=dtype);
        writezarr(bw, tmpFnout);
end

movefile(tmpFnout, fnout);

end

