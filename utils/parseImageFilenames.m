function [fnames, fsns, fd_inds, filepaths] = parseImageFilenames(dataPaths, zarrFile, ChannelPatterns)
% parse image filenames for given data paths for tiff or zarr for given
% channel pattens. This is a simpler version of XR_parseImageFilenames.m
% for only existing Tiff or Zarr files with given channel patterns.
%
% Author: Xiongtao Ruan (04/17/2024)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('zarrFile', @islogical);
ip.addRequired('ChannelPatterns', @iscell);

ip.parse(dataPaths, zarrFile, ChannelPatterns);

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

fprintf('Parse image filenames for data path(s):\n    %s\n', strjoin(dataPaths, '\n    '));

nd = numel(dataPaths);
nc = numel(ChannelPatterns);
fnames = cell(nd, 1);

for d = 1 : nd
    dataPath = dataPaths{d};
    if zarrFile
        dir_info = dir([dataPath, '/', '*.zarr']);        
    else
        dir_info = dir([dataPath, '/', '*.tif']);
    end
    fnames_d = {dir_info.name}';
    nF = numel(fnames_d);
    ch_inds = false(nF, nc);
    for c = 1 : nc
        ch_inds(:, c) = contains(fnames_d, ChannelPatterns{c}, 'IgnoreCase', true) | ...
            contains(fnames_d, regexpPattern(ChannelPatterns{c}), 'IgnoreCase', true);
    end
    fnames_d = fnames_d(any(ch_inds, 2));
    fnames{d} = fnames_d;
end

fd_inds = arrayfun(@(x) ones(numel(fnames{x}), 1) * x, 1 : nd, 'unif', 0);
fnames = cat(1, fnames{:});
[~, fsns] = fileparts(fnames);
fd_inds = cat(1, fd_inds{:});
nF = numel(fnames);
if nF == 0
    if zarrFile
        warning('No Zarr format (*.zarr) image files found in the input data folder(s)! Please check if dataPaths are correct!');
    else
        warning('No Tiff format (*.tif) image files found in the input data folder(s)! Please check if dataPaths are correct!');
    end
end
if nF == 1
    fsns = {fsns};
end

if nargout == 4
    filepaths = arrayfun(@(x) sprintf('%s/%s', dataPaths{fd_inds(x)}, fnames{x}), ...
        1 : nF, 'UniformOutput', false);
end

end
