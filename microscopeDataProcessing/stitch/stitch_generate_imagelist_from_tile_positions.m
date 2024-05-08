function [] = stitch_generate_imagelist_from_tile_positions(dataPath, varargin)
% generate image list file from tile positions files. 
% It obtain related information from tile file names in the given dataPath, and generate image
% list based on image tile indices, image sizes and tile overlap size. 
% It saves as a csv file with name ImageList_from_tile_positions.csv 
% in the dataPath. The format is consistent with old csv files. 
% 
%
% Author: Xiongtao Ruan (08/19/2023)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPath', @(x) ischar(x) || iscell(x));
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0'}, @iscell);
ip.addParameter('DS', false, @islogical);
ip.addParameter('DSR', false, @islogical);
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.2, @isnumeric);
ip.addParameter('SkewAngle', 32.45, @isnumeric);
ip.addParameter('Reverse', false, @islogical);
ip.addParameter('axisOrder', 'x,y,z', @ischar); % stage coordinates axis order
ip.addParameter('dataOrder', 'y,x,z', @ischar); % data axis order, in case data in zyx order
ip.addParameter('objectiveScan', false, @islogical);
ip.addParameter('IOScan', false, @islogical);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('overlapSize', [100, 100, 0], @isnumeric);
ip.addParameter('onlyFirstTP', false, @islogical); % only compute first time point (for deciding cropping bouding box)
ip.addParameter('uuid', '', @ischar);

ip.parse(dataPath, varargin{:});

pr = ip.Results;
% Overwrite = pr.Overwrite;
channelPatterns = pr.channelPatterns;
DS = pr.DS;
DSR = pr.DSR;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
SkewAngle = pr.SkewAngle;
Reverse = pr.Reverse;
axisOrder = pr.axisOrder;
dataOrder = pr.dataOrder;
objectiveScan = pr.objectiveScan;
IOScan =  pr.IOScan;
zarrFile = pr.zarrFile;
overlapSize = pr.overlapSize;
uuid = pr.uuid;

if isempty(uuid)
    uuid = get_uuid();
end

pixelSizes = [xyPixelSize(1), xyPixelSize(end), dz];

% process filenames to obtain information for time point, channel, and tile indices

dir_info = dir([dataPath, '*.tif']);
fsn = {dir_info.name}';
fn = cellfun(@(x) [dataPath, x], fsn, 'UniformOutput', false);

specifyCam = true;
if all(~cellfun(@isempty, regexp(fsn, '_Cam\w_ch', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_Cam(?<Cam>\w+)_ch(?<ch>\d+)_CAM\d_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t(?<suffix>_?\w*)';
elseif all(~cellfun(@isempty, regexp(fsn, '_ch[0-9]_', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_ch(?<ch>\d+)_CAM\d_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t(?<suffix>_?\w*)';
    specifyCam = false;
end

tmp = regexpi(fsn, expression, 'names');
matched_inds = true(numel(tmp), 1);

tab = cell2table(cat(2, fn, fsn), 'VariableNames', ["Filepath", "Filename"]);

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

% get image sizes for all tiles
nF = size(tab, 1);
imSizes = zeros(nF, 3);

for f = 1 : nF
    fn_f = tab.Filepath{f};
    imSizes(f, :) = getImageSize(fn_f);
    if f > 1 && all(imSizes(f, :) == imSizes(1, :))
        imSizes = repmat(imSizes(1, :), nF, 1);
        break;
    end
end

% generate coordiantes 
tileIndices = [tab.x, tab.y, tab.z];
ts = min(tileIndices);
tt = max(tileIndices);

sz = mode(imSizes);
sz_1 = sz - overlapSize;

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

                xyz(inds, :) = (tileIndices - ts) .* sz_1([2, 1, 3]) .* pixelSizes([2, 1, 3]);
            end
        end
    end
end

if false
    figure, scatter3(xyz(:, 1), xyz(:, 2), xyz(:, 3), '*')
end

% preallocate for 1000 rows, which should be enough for most siutations.
sz = [nF, 8];
varTypes = ["string", "string", "double", "double", "double", "double", "double", "double"];
varNames = ["Filepath","Filename","StageX_um_","StageY_um_","StageZ_um_","ObjectiveX_um_","ObjectiveY_um_","ObjectiveZ_um_"];

t = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);

counter = 0;
for f = 1 : nF
    t(counter + 1, :) = {tab.Filepath{f}, tab.Filename{f}, xyz(f, 1), xyz(f, 2), xyz(f, 3), xyz(f, 1), xyz(f, 2), xyz(f, 3)};
    counter = counter + 1;
end
t(counter + 1 : end, :) = [];

tmpout = sprintf('%s/ImageList_from_tile_positions_%s.csv', dataPath, uuid);
fnout = sprintf('%s/ImageList_from_tile_positions.csv', dataPath);
writetable(t, tmpout, 'Delimiter', ',');
fileattrib(tmpout, '+w', 'g');
movefile(tmpout, fnout);

fprintf('Image list file %s is generated!\n', fnout);

end

