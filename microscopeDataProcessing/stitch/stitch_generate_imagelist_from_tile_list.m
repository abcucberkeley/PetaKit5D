function [] = stitch_generate_imagelist_from_tile_list(dataPaths, tileFilenames, tileIndices, tileInterval, varargin)
% generate image list file from a given list of tile files.
% 
% Required inputs:
%           dataPaths : char or cell array. Directory paths for the datasets. Either a string for a single dataset or a cell array of paths for several datasets with same experimental settings.
%       tileFilenames : a #file x 1 cell array (default: {}). List of tile filenames
%         tileIndices : a #file x 5 array (default: []). Tile indices for corresponding tiles in tileFilenames, order: tcxyz.
%        tileInterval : a 1x3 array (default: []). Interval between adjancy tiles in um, order: xyz.
%
% Parameters (as 'specifier'-value pairs):
%     channelPatterns : a cell array (default: {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}).  Channel identifiers for included channels. 
%        tilePatterns : a 1x5 cell array (default: {'0000t', 'ch0', '000x', '000y', '000z'}). Patterns for time, channel, x, y and z to localize tiles. It should be the combination of word and numbers in the form of [a-zA-Z]*[0-9]* or [0-9]*[a-zA-Z]*.
%       tileFilenames : a #file x 1 cell array (default: {}). List of tile filenames
%         tileIndices : a #file x 5 array (default: []). Tile indices for corresponding tiles in tileFilenames, order: tcxyz.
%        tileInterval : a 1x3 array (default: []). Interval between adjancy tiles in um, order: xyz.
%                  DS : true|false (default: false). Data is in deskewed space.
%                 DSR : true|false (default: false). Data is in deskew/rotated space (with stage coordinates).
%         xyPixelSize : a number (default: 0.108). Pixel size in um.
%                  dz : a number (default: 0.5). Scan interval in um.
%           skewAngle : a number (default: 32.45). Skew angle (in degree) of the stage.
%           axisOrder : char (default: 'xyz'). Axis order mapping for coordinates in image list. With combinations of -, x, y, and z. '-yxz' means negative -y map to x, x maps to y, and z maps to z.
%           dataOrder : char (default: 'y,x,z'). Axis order mapping for data. 'y,x,z' means the first, second and third axes are y, x, and z, respectively.
%       objectiveScan : true|false (default: false). Objective scan.
%              IOScan : true|false (default: false). Inverted objective scan. This is the scan with the stage coordinates (DSR space). 
%            zarrFile : true|false (default: false). Use Zarr file as input.
%         overlapSize : empty or 1x3 vector (default: []). Overlap size between tiles in pixel or um. If in pixels, axis order is yxz; if in um, axis order is xyz.
%     overlapSizeType : 'pixel'|'um' (default: 'pixel'). The unit for the overlap size.
%                uuid : empty or a uuid string (default: ''). uuid string as part of the temporate result paths.
%
% Author: Xiongtao Ruan (11/18/2024)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('tileFilenames', @(x) ischar(x) || iscell(x));
ip.addRequired('tileIndices', @(x) isnumeric(x) && size(x, 2) == 5);
ip.addRequired('tileInterval', @(x) isnumeric(x) && size(x, 2) == 3);
ip.addParameter('uuid', '', @ischar);

ip.parse(dataPaths, tileFilenames, tileIndices, tileInterval, varargin{:});

pr = ip.Results;
uuid = pr.uuid;

if isempty(uuid)
    uuid = get_uuid();
end

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

fnames = tileFilenames;
nF = numel(fnames);

if size(tileIndices, 1) ~= nF
    error('The number of rows of tileIndices (%d) does not match the number of tiles (%d)!', size(tileIndices, 1), nF);
end

filepaths = cellfun(@(x) [dataPaths{1}, x], fnames, 'UniformOutput', false);

tab = cell2table(cat(2, filepaths(:), fnames(:)), 'VariableNames', ["Filepath", "Filename"]);

for f = 1 : nF
    tab.prefix{f} = 'Scan';
    tab.Iter(f) = tileIndices(f, 1);
    tab.subIter{f} = '';
    tab.fullIter{f} = sprintf('%04d', tileIndices(f, 1));
    tab.camera(f) = 'A';

    tab.ch(f) = tileIndices(f, 2);
    tab.stack(f) = 0;
    tab.x(f) = tileIndices(f, 3);
    tab.y(f) = tileIndices(f, 4);
    tab.z(f) = tileIndices(f, 5);
    tab.t(f) = tileIndices(f, 1);
    tab.mappedFilename(f) = {sprintf(['Scan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_', ...
        '488nm_000000msec_000000msecAbs_%03dx_%03dy_%03dz_%04dt.tif'], ...
        tab.fullIter{f}, tab.camera(f), tab.ch(f), tab.stack(f), tab.x(f), ...
        tab.y(f), tab.z(f), tab.t(f))};
end


% get image sizes for all tiles
nF = size(tab, 1);

% generate coordiantes 
tileIndices_xyz = [tab.x, tab.y, tab.z];
ts = min(tileIndices_xyz);
tt = max(tileIndices_xyz);

% calculate physical size in xyz order
fullIter = unique(tab.fullIter);
Ch = unique(tab.ch);
Cam = unique(tab.camera);
stackn = unique(tab.stack);

xyz = zeros(nF, 3);

for n = 1:numel(fullIter)
    for ncam = 1:numel(Cam)
        for s = 1:numel(stackn)
            for c = 1:numel(Ch)
                inds = tab.ch == Ch(c) & tab.camera == Cam(ncam) & strcmp(tab.fullIter, fullIter{n}) & tab.stack == stackn(s);
                cur_tab = tab(inds, :);                            

                tileIndices_xyz = [cur_tab.x, cur_tab.y, cur_tab.z];

                xyz(inds, :) = (tileIndices_xyz - ts) .* tileInterval;
            end
        end
    end
end

% visualize tile positions
if false
    figure, scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), '*')
end

% save image list csv file
counter = 0;
sz = [nF, 9];
varTypes = ["string", "string", "string", "double", "double", "double", "double", "double", "double"];
varNames = ["Filepath","Filename","mappedFilename","StageX_um_","StageY_um_","StageZ_um_","ObjectiveX_um_","ObjectiveY_um_","ObjectiveZ_um_"];

t = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);

for f = 1 : nF
    t(counter + 1, :) = {tab.Filepath{f}, tab.Filename{f}, tab.mappedFilename{f}, xyz(f, 1), xyz(f, 2), xyz(f, 3), xyz(f, 1), xyz(f, 2), xyz(f, 3)};
    counter = counter + 1;
end
t(counter + 1 : end, :) = [];

tmpout = sprintf('%s/ImageList_from_tile_list_%s.csv', dataPaths{1}, uuid);
fnout = sprintf('%s/ImageList_from_tile_list.csv', dataPaths{1});
writetable(t, tmpout, 'Delimiter', ',');
fileattrib(tmpout, '+w', 'g');
movefile(tmpout, fnout);

fprintf('Image list file %s is generated!\n', fnout);

end

