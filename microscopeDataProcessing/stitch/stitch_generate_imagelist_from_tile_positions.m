function [] = stitch_generate_imagelist_from_tile_positions(dataPaths, varargin)
% generate image list file from tile positions files. 
% It obtain related information from tile file names in the given dataPath, and generate image
% list based on image tile indices, image sizes and tile overlap size. 
% It saves as a csv file with name ImageList_from_tile_positions.csv 
% in the dataPath. The format is consistent with old csv files. 
% 
%
% Author: Xiongtao Ruan (06/07/2024)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('tilePatterns', {'0000t', 'ch0', '000x', '000y', '000z'}, @iscell);
ip.addParameter('DS', false, @islogical);
ip.addParameter('DSR', false, @islogical);
ip.addParameter('xyPixelSize', 0.108, @isnumeric); % in um
ip.addParameter('dz', 0.2, @isnumeric); % in um
ip.addParameter('skewAngle', 32.45, @isnumeric);
ip.addParameter('axisOrder', 'x,y,z', @ischar); % tile indices axis order, -x means scan in a reversed direction in x axis
ip.addParameter('dataOrder', 'y,x,z', @ischar); % data axis order
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('IOScan', false, @islogical);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('overlapSize', [], @isnumeric); % in stage coordinates
ip.addParameter('overlapSizeType', 'pixel', @(x) ischar(x) && ismember(lower(x), {'pixel', 'um'})); % yxz in pixel, xyz in um
ip.addParameter('uuid', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
% Overwrite = pr.Overwrite;
channelPatterns = pr.channelPatterns;
tilePatterns = pr.tilePatterns;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
skewAngle = pr.skewAngle;
axisOrder = pr.axisOrder;
dataOrder = pr.dataOrder;
DS = pr.DS;
DSR = pr.DSR;
objectiveScan = pr.objectiveScan;
IOScan =  pr.IOScan;
zarrFile = pr.zarrFile;
overlapSize = pr.overlapSize;
overlapSizeType = pr.overlapSizeType;
uuid = pr.uuid;

if isempty(uuid)
    uuid = get_uuid();
end

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

pixelSizes = [xyPixelSize(1), xyPixelSize(end), dz];

% parse image filenames
[fnames, fsns, fd_inds, filepaths] = parseImageFilenames(dataPaths, zarrFile, channelPatterns);
nF = numel(fnames);

% process filenames to obtain information for time point, channel, and tile indices
[tileInfo, camPatterns] = parseTileInfo(fsns, tilePatterns);

specifyCam = true;
if all(~cellfun(@isempty, regexp(fsns, '_Cam\w_ch', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_Cam(?<Cam>\w+)_ch(?<ch>\d+)_CAM\d_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t(?<suffix>_?\w*)';
elseif all(~cellfun(@isempty, regexp(fsns, '_ch[0-9]_', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_ch(?<ch>\d+)_CAM\d_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t(?<suffix>_?\w*)';
    specifyCam = false;
else
    expression = '';
    specifyCam = false;
end

tmp = regexp(fsns, expression, 'tokens');
empty_inds = cellfun(@(x) isempty(x), tmp);
standard_naming = ~any(empty_inds);

if standard_naming
    tmp = regexpi(fsns, expression, 'names');
    matched_inds = true(numel(tmp), 1);
    
    tab = cell2table(cat(2, filepaths, fnames), 'VariableNames', ["Filepath", "Filename"]);
    
    for f = 1:numel(tmp)
        if isempty(tmp{f})
            matched_inds(f) = false;
            continue;
        end
        
        tab.prefix{f} = tmp{f}.prefix;
        tab.Iter(f) = str2double(tmp{f}.Iter);
        tab.subIter{f} = tmp{f}.subIter;
        tab.fullIter{f} = [tmp{f}.Iter, tmp{f}.subIter];
        if specifyCam
            tab.camera(f) = (tmp{f}.Cam);
        else
            % use N to represent cam if it is not contained in the filename
            tab.camera(f) = 'N';
        end
        tab.ch(f) = str2double(tmp{f}.ch);
        tab.stack(f) = str2double(tmp{f}.stack);
        tab.laser(f) = str2double(tmp{f}.laser);
        tab.abstime(f) = str2double(tmp{f}.abstime);
        tab.fpgatime(f) = str2double(tmp{f}.fpgatime);
        tab.x(f) = str2double(tmp{f}.x);
        tab.y(f) = str2double(tmp{f}.y);
        tab.z(f) = str2double(tmp{f}.z);
        tab.t(f) = str2double(tmp{f}.t);
        tab.did(f) = 1; % data path id
    end
    
    tab = tab(matched_inds, :);
else
    matched_inds = true(nF, 1);
    tab = cell2table(cat(2, filepaths, fnames), 'VariableNames', ["Filepath", "Filename"]);
    
    for f = 1 : nF
        tileInfo_f = tileInfo(f, :);
        if any(isnan(tileInfo_f))
            matched_inds(f) = false;
            continue;
        end

        tab.prefix{f} = 'Scan';
        tab.Iter(f) = tileInfo_f(1);
        tab.subIter{f} = '';
        tab.fullIter{f} = sprintf('%04d', tileInfo_f(1));
        if specifyCam
            tab.camera(f) = camPatterns{tileInfo_f(6)};
        else
            % use A to represent cam if it is not contained in the filename
            tab.camera(f) = 'A';
        end
        tab.ch(f) = tileInfo_f(2);
        tab.stack(f) = 0;
        tab.x(f) = tileInfo_f(3);
        tab.y(f) = tileInfo_f(4);
        tab.z(f) = tileInfo_f(5);
        tab.t(f) = tileInfo_f(1);
        if specifyCam        
            tab.mappedFilename(f) = {sprintf(['Scan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_', ...
                '488nm_000000msec_000000msecAbs_%03dx_%03dy_%03dz_%04dt.tif'], ...
                tab.fullIter{f}, tab.camera(f), tab.ch(f), tab.stack(f), tab.x(f), ...
                tab.y(f), tab.z(f), tab.t(f))};
        else
            tab.mappedFilename(f) = {sprintf(['Scan_Iter_%s_ch%d_CAM1_stack%04d_', ...
                '488nm_000000msec_000000msecAbs_%03dx_%03dy_%03dz_%04dt.tif'], ...
                tab.fullIter{f}, tab.ch(f), tab.stack(f), tab.x(f), ...
                tab.y(f), tab.z(f), tab.t(f))};
        end
    end
end

% get image sizes for all tiles
nF = size(tab, 1);
imSizes = getImageSizeBatch(tab.Filepath);

% generate coordiantes 
tileIndices = [tab.x, tab.y, tab.z];
ts = min(tileIndices);
tt = max(tileIndices);

sz = mode(imSizes);

% coordinates conversion for skewed space and DS space data
if objectiveScan || IOScan
    zAniso = pixelSizes(3) / pixelSizes(1);
    DS = false;
    DSR = false;
else
    zAniso = sind(skewAngle) * pixelSizes(3) / pixelSizes(1);
end
theta = skewAngle * pi / 180;

switch lower(overlapSizeType)
    case 'pixel'
        overlapSize = overlapSize([2, 1, 3]) .* pixelSizes(1) .* [1, 1, zAniso];
    case 'um'
    otherwise
        error('Unknown overlap size type, must be pixel or um!');
end

% calculate physical size in xyz order
sz_phy = sz([2, 1, 3]) .* pixelSizes(1) .* [1, 1, zAniso];
sz_dsr = sz_phy;
if objectiveScan || DS
    % convert size in DS space to DSR space (stage coordinates)
    sz_dsr = [sz_phy(1) * cos(theta) + sz_phy(3) * sin(theta), sz_phy(2), sz_phy(1) * sin(theta) - sz_phy(3) * cos(theta)];
end

if ~IOScan && ~objectiveScan && ~DS && ~DSR
    % convert size in skewed space to DSR space (stage coordinates)
    sz_dsr = [sz_phy(1) * cos(theta) + sz_phy(3) / sin(theta), sz_phy(2), sz_phy(1) * sin(theta)];    
end

sz_1 = sz_dsr - overlapSize;

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

                tileIndices = [cur_tab.x, cur_tab.y, cur_tab.z];

                xyz(inds, :) = (tileIndices - ts) .* sz_1([2, 1, 3]);
            end
        end
    end
end

% visualize tile positions
if false
    figure, scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), '*')
end

% save image list csv file

% preallocate for 1000 rows, which should be enough for most siutations.
counter = 0;
if standard_naming
    sz = [nF, 8];
    varTypes = ["string", "string", "double", "double", "double", "double", "double", "double"];
    varNames = ["Filepath","Filename","StageX_um_","StageY_um_","StageZ_um_","ObjectiveX_um_","ObjectiveY_um_","ObjectiveZ_um_"];
    
    t = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
    
    for f = 1 : nF
        t(counter + 1, :) = {tab.Filepath{f}, tab.Filename{f}, xyz(f, 1), xyz(f, 2), xyz(f, 3), xyz(f, 1), xyz(f, 2), xyz(f, 3)};
        counter = counter + 1;
    end
else
    sz = [nF, 9];
    varTypes = ["string", "string", "string", "double", "double", "double", "double", "double", "double"];
    varNames = ["Filepath","Filename","mappedFilename","StageX_um_","StageY_um_","StageZ_um_","ObjectiveX_um_","ObjectiveY_um_","ObjectiveZ_um_"];
    
    t = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
    
    for f = 1 : nF
        t(counter + 1, :) = {tab.Filepath{f}, tab.Filename{f}, tab.mappedFilename{f}, xyz(f, 1), xyz(f, 2), xyz(f, 3), xyz(f, 1), xyz(f, 2), xyz(f, 3)};
        counter = counter + 1;
    end
end
t(counter + 1 : end, :) = [];

tmpout = sprintf('%s/ImageList_from_tile_positions_%s.csv', dataPaths{1}, uuid);
fnout = sprintf('%s/ImageList_from_tile_positions.csv', dataPaths{1});
writetable(t, tmpout, 'Delimiter', ',');
fileattrib(tmpout, '+w', 'g');
movefile(tmpout, fnout);

fprintf('Image list file %s is generated!\n', fnout);

end

