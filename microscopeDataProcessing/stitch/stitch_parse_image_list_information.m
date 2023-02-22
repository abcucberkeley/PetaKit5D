function [tab, primary_tab, fullIter, Ch, Cam, stackn, nz, specifyCam, prefix, zlayerStitch, ...
    stitchInfoFullpath] = stitch_parse_image_list_information(dataPath, imageListFileName, options)
% move the parase image list code from the main wrapper to this function 
%
% Author: Xiongtao Ruan (02/15/2023)


arguments
    dataPath char    
    imageListFileName char 
    options.Streaming (1, 1) {islogical} = false
    options.stitchInfoFullpath char = ''
    options.stitchInfoPath char = ''
    options.onlyFirstTP (1, 1) {islogical} = false
    options.ChannelPatterns {iscell} = {'CamA', 'CamB'}
    options.useProcessedData {islogical} = false
    options.ProcessedDirStr char = ''
    options.timepoints (:, 1) {mustBeNumeric} = []
    options.subtimepoints (:, 1) {mustBeNumeric} = []
    options.xcorrMode char {mustBeMember(options.xcorrMode, {'primaryFirst', 'primary', 'all'})} = 'primaryFirst'
    options.primaryCh char = ''
    options.onlineStitch (1, 1) {islogical} = false
end

Streaming = options.Streaming;
stitchInfoFullpath = options.stitchInfoFullpath;
stitch_info_path = options.stitchInfoPath;
onlyFirstTP = options.onlyFirstTP;
ChannelPatterns = options.ChannelPatterns;
useProcessedData = options.useProcessedData;
ProcessedDirStr = options.ProcessedDirStr;
timepoints = options.timepoints;
subtimepoints = options.subtimepoints;
xcorrMode = options.xcorrMode;
primaryCh = options.primaryCh;
onlineStitch = options.onlineStitch;

% read image list csv file
tab = readtable(imageListFileName, 'Delimiter','comma');
% t_column_name = t.Properties.VariableNames;
fn = tab.Filename;
% remove records for partial volume files
partialvol = ~cellfun(@isempty, regexpi(fn, '_part\d+.tif'));
if any(partialvol)
    tab(partialvol, :) = [];
    fn = tab.Filename;
end   

specifyCam = true;
if all(~cellfun(@isempty, regexp(fn, '_Cam\w_ch', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_Cam(?<Cam>\w+)_ch(?<ch>\d+)_CAM\d_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t(?<suffix>_?\w*).tif';
elseif all(~cellfun(@isempty, regexp(fn, '_ch[0-9]_', 'match')))
    expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>\d+)(?<subIter>_?(\d+_)*\d+?)_ch(?<ch>\d+)_CAM\d_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t(?<suffix>_?\w*).tif';
    specifyCam = false;
end

tmp = regexpi(fn, expression, 'names');

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
% sort t by tile indices, time, Iter
tab = sortrows(tab, {'z', 'y', 'x', 't', 'fullIter'});

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
if Streaming && onlineStitch
    nz = unique(tab.z);
    zlayerStitch = true;
    xcorrMode = 'all';
else
    nz = 1;
    zlayerStitch = false;    
end

% check whether the image files in the image list file exist 
dir_info = dir([dataPath, '/', '*.tif']);
imageFnames = {dir_info.name}';
if isempty(dir_info)
    if useProcessedData
        dir_info = dir([dataPath, ProcessedDirStr, '/', '*.tif']);
        imageFnames = {dir_info.name}';
    elseif useExistDecon
        dir_info = dir([dataPath, '/', DeconDirstr, '/', '*.tif']);
        imageFnames = cellfun(@(x) [x(1 : end - 10), '.tif'], {dir_info.name}', 'unif', 0);
    else
        error('The tiles do not exist!');
    end
end

% filter filenames by channel patterns
include_flag = false(numel(imageFnames), 1);
imageFullnames = cellfun(@(x) [dataPath, '/', x], imageFnames, 'unif', 0);
for c = 1 : numel(ChannelPatterns)
    include_flag = include_flag | contains(imageFullnames, ChannelPatterns{c}) | contains(imageFullnames, regexpPattern(ChannelPatterns{c}));
end
imageFnames = imageFnames(include_flag);

[~, tFsnames] = fileparts(tab.Filename);
image_file_exist_flag = true(numel(tFsnames), 1);
for f = 1 : numel(tab.Filename)
    if ~any(contains(imageFnames, tFsnames{f})) % || (useExistDecon && ~any(contains(imageFnames, tFsnames{f})))
        image_file_exist_flag(f) = false;
    end
end
if ~all(image_file_exist_flag)
    warning('Some files in the image list file do not exist! Ignore them in the stitching')
    disp(tab.Filename(~image_file_exist_flag));
    if ~Streaming
        tab(~image_file_exist_flag, :) = [];
        
        fullIter = unique(tab.fullIter);
        Ch = unique(tab.ch);
        Cam = unique(tab.camera);
        % ntiles = numel(unique(t.x)) * numel(unique(t.y)) * numel(unique(t.z));
        stackn = unique(tab.stack);
    end
end

if ~isempty(stitchInfoFullpath) 
    if exist(stitchInfoFullpath, 'file')
        % xcorrShift = true;
        xcorrMode = 'stitchInfo';
    else
        error('The user defined stitch info file %s does not exist!', stitchInfoFullpath);
    end
end

if (strcmpi(xcorrMode, 'primary') || strcmpi(xcorrMode, 'primaryFirst'))
    % first check whether the format is right
    if ~isempty(primaryCh)
        fprintf('The primary channel is %s.\n', primaryCh);
        if specifyCam && isempty(regexpi(primaryCh, 'Cam[a-z]_ch[0-9]', 'once'))
            error("primaryCh must be empty or with format 'Cam[a-z]_ch[0-9]'");
        elseif ~specifyCam && isempty(regexpi(primaryCh, 'ch[0-9]', 'once'))
            error("primaryCh must be empty or with format 'ch[0-9]'");            
        end
        pCam = primaryCh(4);
        pCh = str2double(primaryCh(end));
        if ~any(contains(fn, ['_', primaryCh, '_'], 'IgnoreCase', true))
            error('The given primary channel %s does not exist!', primaryCh);
        end
    else
        warning('The primary channel is not set, use the first available channel as primary channel...');
        pCam = [];
        for ncam = 1 : numel(Cam)
            for c = 1 : numel(Ch)
                if specifyCam
                    primaryCh = sprintf('Cam%s_ch%d', Cam(ncam), Ch(c));
                else
                    primaryCh = sprintf('ch%d', Ch(c));
                end
                if any(contains(fn, ['_', primaryCh, '_'], 'IgnoreCase', true))
                    pCam = Cam(ncam);
                    pCh =  Ch(c);
                    break;
                end
            end
            if ~isempty(pCam)
                break;
            end
        end       
        fprintf('Set primary channel as %s.\n', primaryCh);
    end

    % reorder Ch and Cam to make primary channel the first
    if Ch(1) ~= pCh
        Ch = [pCh; Ch(Ch ~= pCh)];
    end
    if ~strcmpi(Cam(1), pCam)
        ind = regexpi(Cam', pCam);
        Cam = Cam([ind, 1 : ind - 1, ind + 1 : end]);
    end
end

if onlyFirstTP
    fullIter = fullIter(1);
    timepoints = [];
end

% xruan (01/21/2021): modify to separate timepoints, subtimepoints, subsubtimepoints
if ~isempty(timepoints)
    Iter = cellfun(@(x) str2double(x(1 : 4)), fullIter);

    fullIter = fullIter(ismember(Iter, timepoints(:)));
    if ~isempty(subtimepoints)
        subIter = cellfun(@(x) str2double(x(6 : 9)), fullIter);
        fullIter = fullIter(ismember([Iter, subIter], [timepoints(:), subtimepoints(:)], 'rows'));
    end

%     if ~isempty(subsubtimepoints)
%         subIter = cellfun(@(x) str2double(x(6 : 9)), fullIter);
%         subSubIter = cellfun(@(x) str2double(x(11 : 14)), fullIter);
%         fullIter = fullIter(ismember([Iter, subIter, subSubIter], [timepoints(:), subtimepoints(:), subsubtimepoints(:)], 'rows'));
%     end
end

% predefine stitchInfo when xcorrMode is 'primaryFirst'
if strcmp(xcorrMode, 'primaryFirst')
    if zlayerStitch
        primary_tab = tab(tab.ch == Ch(1) & tab.camera == Cam(1) & strcmp(tab.fullIter, fullIter{1}) & tab.stack == stackn(1) & tab.z == nz(1), :);        
    else
        primary_tab = tab(tab.ch == Ch(1) & tab.camera == Cam(1) & strcmp(tab.fullIter, fullIter{1}) & tab.stack == stackn(1), :);
    end
    if isempty(primary_tab) 
        error('The Image List Info for the primary channel for the first time point does not exist!');
    end
    
    p_laser = unique(primary_tab.laser);
    p_abstime = unique(primary_tab.abstime);
    p_fpgatime = primary_tab.fpgatime(1);
    p_z = unique(primary_tab.z);
    p_z = p_z(1);

    if numel(p_laser) > 1
        p_laser = p_laser(1);
    end
    if zlayerStitch
        pz_str = sprintf('_%03dz', p_z);
    else
        pz_str = '';
    end

    if specifyCam
        stitchInfoFullpath = sprintf('%s/%sScan_Iter_%s_Cam%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s.mat', ...
            stitch_info_path, prefix, fullIter{1}, Cam(1), Ch(1), stackn(1), p_laser, p_abstime, p_fpgatime, pz_str);
    else
        stitchInfoFullpath = sprintf('%s/%sScan_Iter_%s_ch%d_CAM1_stack%04d_%dnm_%07dmsec_%010dmsecAbs%s.mat', ...
            stitch_info_path, prefix, fullIter{1}, Ch(1), stackn(1), p_laser, p_abstime, p_fpgatime, pz_str);        
    end
end

end

