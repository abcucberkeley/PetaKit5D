function [fsns, imageListFullpath] = writeStitchTiles(dataPath, tiles, varargin)
% write tiles with LLSM style filenames and add them to the image list csv
% used for stitching (Filename, StageX_um_, StageY_um_, StageZ_um_).
%
%   dataPath : char. Dataset folder (with trailing /).
%      tiles : struct array with fields
%                 im : tile image
%                idx : 1x3 [x, y, z] tile indices for the filename
%                xyz : 1x3 [x, y, z] stage coordinates in um
%
% Parameters (as 'specifier'-value pairs):
%                cam : char (default: 'A'). Camera letter in the filenames, '' for no camera.
%                 ch : a number (default: 0). Channel number in the filenames.
%              laser : a number (default: 488). Laser wavelength in the filenames.
%          imageList : char (default: 'ImageList.csv'). Image list filename in dataPath. Rows are appended if it exists.
%            tileDir : char (default: ''). Subfolder of dataPath for the tile files, e.g., 'DSR'.
%         saveFormat : 'tif'|'zarr' (default: 'tif'). Tile file format. The image list always has .tif names.

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPath', @ischar);
ip.addRequired('tiles', @isstruct);
ip.addParameter('cam', 'A', @ischar);
ip.addParameter('ch', 0, @isnumeric);
ip.addParameter('laser', 488, @isnumeric);
ip.addParameter('imageList', 'ImageList.csv', @ischar);
ip.addParameter('tileDir', '', @ischar);
ip.addParameter('saveFormat', 'tif', @(x) any(strcmpi(x, {'tif', 'zarr'})));
ip.parse(dataPath, tiles, varargin{:});

pr = ip.Results;
cam = pr.cam;
ch = pr.ch;
laser = pr.laser;
tileDir = pr.tileDir;
saveFormat = lower(pr.saveFormat);

if isempty(cam)
    camStr = '';
else
    camStr = ['_Cam', cam];
end

outPath = dataPath;
if ~isempty(tileDir)
    outPath = [dataPath, tileDir, '/'];
    if ~exist(outPath, 'dir')
        mkdir(outPath);
    end
end

nT = numel(tiles);
fsns = cell(nT, 1);
xyz = zeros(nT, 3);
for i = 1 : nT
    fsns{i} = sprintf('Scan_Iter_0000_0000%s_ch%d_CAM1_stack0000_%dnm_0000000msec_0000000000msecAbs_%03dx_%03dy_%03dz_0000t', ...
        camStr, ch, laser, tiles(i).idx);
    switch saveFormat
        case 'tif'
            writetiff(tiles(i).im, [outPath, fsns{i}, '.tif']);
        case 'zarr'
            writezarr(tiles(i).im, [outPath, fsns{i}, '.zarr'], blockSize=[32, 32, 32]);
    end
    xyz(i, :) = tiles(i).xyz;
end

imageListFullpath = [dataPath, pr.imageList];
tab = table(cellfun(@(x) [x, '.tif'], fsns, 'unif', 0), xyz(:, 1), xyz(:, 2), xyz(:, 3), ...
    'VariableNames', {'Filename', 'StageX_um_', 'StageY_um_', 'StageZ_um_'});
if exist(imageListFullpath, 'file')
    tab = [readtable(imageListFullpath, 'Delimiter', 'comma'); tab];
end
writetable(tab, imageListFullpath);

end
