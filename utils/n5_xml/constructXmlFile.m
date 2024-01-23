function[s, t, voxelSize_xml] = constructXmlFile(dataPath, fileName,csvFile,voxelSize, ChannelPatterns, zarrFile, axisOrder, dataOrder, xcorrInfoFn)
% voxelSize: xyz order


if nargin < 5
    ChannelPatterns = {'CamA'};
end
if nargin < 6
    zarrFile = false;
end
if nargin < 7
    axisOrder = 'xyz';
end
if nargin < 8
    dataOrder = 'yxz'; % data order for tiff/zarr, opposite order for n5
end
if nargin < 9
    xcorrInfoFn = '';
end

tab = readtable(csvFile, 'Delimiter','comma');
fn = tab.Filename;
if zarrFile
    ext = '.zarr';
else
    ext = '.tif';
end
% remove records for partial volume files
[~, fsn] = fileparts(fn);

useMapFname = false;
if  any(ismember(tab.Properties.VariableNames, 'mappedFilename'))
    mapFsn = tab.mappedFilename;
    useMapFname = true;
    [~, fsn] = fileparts(mapFsn);
end

specifyCam = true;
if all(~cellfun(@isempty, regexp(fsn, '_Cam\w_ch', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_Cam(?<Cam>\w+)_ch(?<ch>\d+)_CAM\d_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t(?<suffix>_?\w*)';
elseif all(~cellfun(@isempty, regexp(fsn, '_ch[0-9]_', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_ch(?<ch>\d+)_CAM\d_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t(?<suffix>_?\w*)';
    specifyCam = false;
end

tmp = regexpi(fsn, expression, 'names');
matched_inds = true(numel(tmp), 1);

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
if zarrFile
    tab.Filename = strrep(tab.Filename, '.tif', ext);
end
% sort t by camera, ch, tile indices, time, Iter
[tab, tinds] = sortrows(tab, {'camera', 'ch', 'z', 'y', 'x', 't', 'fullIter'});
fsn = fsn(tinds);

prefix = unique(tab.prefix);
if ~isempty(prefix)
    prefix = prefix{1};
else
    prefix = '';
end
fullIter = unique(tab.fullIter);
Ch = unique(tab.ch);
Cam = unique(tab.camera);
stackn = unique(tab.stack);

% check whether the image files in the image list file exist 
dir_info = dir([dataPath, '/', '*', ext]);
imageFnames = {dir_info.name}';
if isempty(dir_info)
    % if useProcessedData
    %     dir_info = dir([dataPath, ProcessedDirStr, '/', '*', ext]);
    %     imageFnames = {dir_info.name}';
    % else
        error('The tiles do not exist!');
    % end
end

% filter filenames by channel patterns
if useMapFname
    include_flag = false(numel(fsn), 1);    
    imageFullnames = cellfun(@(x) sprintf('%s/%s%s', dataPath, x, ext), fsn, 'unif', 0);    
    for c = 1 : numel(ChannelPatterns)
        include_flag = include_flag | contains(imageFullnames, ChannelPatterns{c}) ...
            | contains(imageFullnames, regexpPattern(ChannelPatterns{c}));
    end
    if ~all(include_flag)
        tab = tab(include_flag, :);
    end
else
    include_flag = false(numel(imageFnames), 1);
    imageFullnames = cellfun(@(x) [dataPath, '/', x], imageFnames, 'unif', 0);
    for c = 1 : numel(ChannelPatterns)
        include_flag = include_flag | contains(imageFullnames, ChannelPatterns{c}) ...
            | contains(imageFullnames, regexpPattern(ChannelPatterns{c}));
    end
    imageFnames = imageFnames(include_flag);    
end

[~, tFsnames] = fileparts(tab.Filename);
if numel(tab.Filename) == 1
    tFsnames = {tFsnames};
end
image_file_exist_flag = true(numel(tFsnames), 1);
for f = 1 : numel(tab.Filename)
    if ~any(contains(imageFnames, tFsnames{f})) % || (useExistDecon && ~any(contains(imageFnames, tFsnames{f})))
        image_file_exist_flag(f) = false;
    end
end
if ~all(image_file_exist_flag)
    warning('Some files in the image list file do not exist! Ignore them in the stitching')
    disp(tab.Filename(~image_file_exist_flag));
    tab(~image_file_exist_flag, :) = [];
    
    fullIter = unique(tab.fullIter);
    Ch = unique(tab.ch);
    Cam = unique(tab.camera);
    % ntiles = numel(unique(t.x)) * numel(unique(t.y)) * numel(unique(t.z));
    stackn = unique(tab.stack);
end


% process coordinates with axis orders
coordinates = [tab.StageX_um_, tab.StageY_um_, tab.StageZ_um_];
axisOrder = strrep(axisOrder, ' ', '');
pattern = '^(-?x,?-?y,?-?z|-?y,?-?x,?-?z|-?z,?-?y,?-?x|-?x,?-?z,?-?y|-?x,?-?z,?-?y|-?y,?-?z,?-?x)$';
if ~regexpi(axisOrder, pattern)
    error("The axisOrder is not right, it must has the form like 'y,x,z' or '-x,y,z' (flipped in x-axis)!");
end

order_sign_mat = zeros(2, 3); % first row: order for permute; second row: sign
axisOrder_pure = strrep(strrep(axisOrder, ',', ''), '-', '');
xyz_str = 'xyz';
for i = 1 : 3
    [~, order_sign_mat(1, :)] = sort(axisOrder_pure);
    order_sign_mat(2, i) = 1 - 2 * contains(axisOrder, ['-', xyz_str(i)]);
end
xyz = coordinates(:, order_sign_mat(1, :)) .* order_sign_mat(2, :);

% normalize coordinates to zero
xyz = xyz - min(xyz, [], 1);

% get data order 
dataOrder_pure = strrep(strrep(dataOrder, ',', ''), '-', '');
[~, data_order_mat] = sort(dataOrder_pure);
switch dataOrder
    case 'yxz'
        dxyz_order_mat = [1, 2, 3];        
    case 'zyx'
        dxyz_order_mat = [3, 1, 2];
end

% coordinate for xml
xyz_xml = xyz(:, flip(data_order_mat));
xyz_xml = xyz_xml - min(xyz_xml, [], 1);
voxelSize_xml = voxelSize(flip(data_order_mat));

% load xcorr File if it is provided
genStitchInfo = false;
if ~isempty(xcorrInfoFn)
    genStitchInfo = true;
    a = load(xcorrInfoFn);
    dxyz_shift = a.dxyz_shift(:, dxyz_order_mat);
    dxyz_shift_xml = dxyz_shift(:, flip(data_order_mat));
end


%% generate xml file from the information in the image list

t = tab;
nFiles = size(t,1);
s = struct;
s.versionAttribute = 0.2;
s.BasePath.typeAttribute = "relative";
% s.BasePath.typeAttribute = "absolute";
s.BasePath.Text = ".";
s.SequenceDescription.ImageLoader.formatAttribute = "bdv.n5";
s.SequenceDescription.ImageLoader.versionAttribute = 1.0;
s.SequenceDescription.ImageLoader.n5.typeAttribute = "relative";
% s.SequenceDescription.ImageLoader.n5.typeAttribute = "absolute";
[~, fsn, ext] = fileparts(fileName);
s.SequenceDescription.ImageLoader.n5.Text = [fsn, ext];
% ViewSetups
for i = 1:nFiles
    sInner = struct;
    sInner.id = i-1;
    sInner.name = sInner.id;
    
    currFile = [dataPath, t.Filename{i}];
    imSize = getImageSize(currFile);

    sInner.size = sprintf("%d %d %d",imSize(3),imSize(2),imSize(1));
    sInnerVs = struct;
    sInnerVs.unit = "um";
    sInnerVs.size = sprintf("%f %f %f", voxelSize_xml);
    sInner.voxelSize = sInnerVs;
    s.SequenceDescription.ViewSetups.ViewSetup(i) = sInner;
end
% Attributes
sInnerAttr = struct;

sInnerAttr(1).nameAttribute = "illumination";
sInnerAttr(1).Illumination.id = 0;
sInnerAttr(1).Illumination.name = 0;
sInnerAttr(1).Channel = missing;
sInnerAttr(1).Tile = missing;
sInnerAttr(1).Angle = missing;

sInnerAttr(2).nameAttribute = "channel";
sInnerAttr(2).Illumination = missing;
sInnerAttr(2).Channel.id = 0;
sInnerAttr(2).Channel.name = 0;
sInnerAttr(2).Tile = missing;
sInnerAttr(2).Angle = missing;

sInnerAttr(3).nameAttribute = "tile";
sInnerAttr(3).Illumination = missing;
sInnerAttr(3).Channel = missing;
for i = 1:nFiles
    sInnerAttr(3).Tile(i).id = i-1;
    sInnerAttr(3).Tile(i).name = i-1;
end
sInnerAttr(3).Angle = missing;

sInnerAttr(4).nameAttribute = "angle";
sInnerAttr(4).Illumination = missing;
sInnerAttr(4).Channel = missing;
sInnerAttr(4).Tile = missing;
sInnerAttr(4).Angle.id = 0;
sInnerAttr(4).Angle.name = 0;

s.SequenceDescription.ViewSetups.Attributes = sInnerAttr;
% Timepoints
s.SequenceDescription.Timepoints.typeAttribute = "pattern";
s.SequenceDescription.Timepoints.integerpattern = 0;
s.SequenceDescription.MissingViews = "";

% ViewRegistration
for i = 1:nFiles
    sInner = struct;
    sInner.timepointAttribute = 0;
    sInner.setupAttribute = i-1;
    sInnerVt = struct;
    %TESTING
    %sInnerVt(1).typeAttribute = "affine";
    %sInnerVt(1).Name = "Stitching Transform";
    %sInnerVt(1).affine = sprintf("1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0");
    sInnerVt(1).typeAttribute = "affine";
    sInnerVt(1).Name = "Translation to Regular Grid";
    sInnerVt(1).affine = sprintf("1.0 0.0 0.0 %f 0.0 1.0 0.0 %f 0.0 0.0 1.0 %f", xyz_xml(i, 1) ./ voxelSize_xml(1), xyz_xml(i, 2) ./ voxelSize_xml(2), xyz_xml(i, 3) ./ voxelSize_xml(3));

    sInnerVt(2).typeAttribute = "affine";
    sInnerVt(2).Name = "calibration";
    sInnerVt(2).affine = sprintf("1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0");
    %sInnerVt(2).affine = sprintf("%f 0.0 0.0 0.0 0.0 %f 0.0 0.0 0.0 0.0 %f 0.0",voxelSize);

    if genStitchInfo
        sInnerVt(3).typeAttribute = "affine";
        sInnerVt(3).Name = "Stitching Transform";
        sInnerVt(3).affine = sprintf("1.0 0.0 0.0 %f 0.0 1.0 0.0 %f 0.0 0.0 1.0 %f", dxyz_shift_xml(i, 1) ./ voxelSize_xml(1), dxyz_shift_xml(i, 2) ./ voxelSize_xml(2), dxyz_shift_xml(i, 3) ./ voxelSize_xml(3));
    end

    sInner.ViewTransform = sInnerVt;
    s.ViewRegistrations.ViewRegistration(i) = sInner;
end
s.ViewInterestPoints = "";
s.BoundingBoxes = "";
s.PointSpreadFunctions = "";
%s.StitchingResults = "";
s.IntensityAdjustments = "";

end
