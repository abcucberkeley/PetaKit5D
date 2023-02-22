function [tab, primary_tab, fullIter, Ch, Cam, stackn, nz, specifyCam, prefix, zlayerStitch, ...
    stitchInfoFullpath] = stitch_parse_multi_loc_image_list_information(dataPath, imageListFileName, options)
% parase image list for multiple location dataset
%
% Author: Xiongtao Ruan (02/16/2023)


arguments
    dataPath cell    
    imageListFileName cell 
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


% parase imagel list for each subregion
nd = numel(dataPath);

tab_cell = cell(nd, 1);
primary_tab_cell = cell(nd, 1);
fullIter_cell = cell(nd, 1);
Ch_cell = cell(nd, 1);
Cam_cell = cell(nd, 1);
stackn_cell = cell(nd, 1);
nz_cell = cell(nd, 1);
specifyCam_cell = cell(nd, 1);
prefix_cell = cell(nd, 1);
zlayerStitch_cell = cell(nd, 1);
stitchInfoFullpath_cell = cell(nd, 1);


for d = 1 : nd
    [tab, primary_tab, fullIter, Ch, Cam, stackn, nz, specifyCam, prefix, zlayerStitch, stitchInfoFullpath_d] = ...
        stitch_parse_image_list_information(dataPath{d}, imageListFileName{d}, Streaming=Streaming, ...
        onlineStitch=onlineStitch, stitchInfoFullpath=stitchInfoFullpath, stitchInfoPath=stitch_info_path, ...
        onlyFirstTP=onlyFirstTP, ChannelPatterns=ChannelPatterns, useProcessedData=useProcessedData, ...
        ProcessedDirStr=ProcessedDirStr, timepoints=timepoints, ...
        subtimepoints=subtimepoints, xcorrMode=xcorrMode, primaryCh=primaryCh);
    
    tab.did(:) = d;

    tab_cell{d} = tab;
    primary_tab_cell{d} = primary_tab;
    fullIter_cell{d} = fullIter;
    Ch_cell{d} = Ch;
    Cam_cell{d} = Cam;
    stackn_cell{d} = stackn;
    nz_cell{d} = nz;
    specifyCam_cell{d} = specifyCam;
    prefix_cell{d} = prefix;
    zlayerStitch_cell{d} = zlayerStitch;
    stitchInfoFullpath_cell{d} = stitchInfoFullpath_d;
end

% integrate image list information for all subregions
tab = cat(1, tab_cell{:});
fullIter = unique(cat(1, fullIter_cell{:}));
Ch = unique(cat(1, Ch_cell{:}));
Cam = unique(cat(1, Cam_cell{:}));
stackn = unique(cat(1, stackn_cell{:}));
nz = unique(cat(1, nz_cell{:}));
specifyCam = unique(cat(1, specifyCam_cell{:}));
prefix = unique(cat(1, prefix_cell{:}));
zlayerStitch = unique(cat(1, zlayerStitch_cell{:}));
stitchInfoFullpath = stitchInfoFullpath_cell{1};

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
end


end

