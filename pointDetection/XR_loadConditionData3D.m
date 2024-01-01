
% loadConditionData3D loads the relevant information for all the movies
% available for a specific condition; this requires a specific directory
% structure and nomenclature (see below)
%
% SYNOPSIS [data] = loadConditionData3D()
%
% INPUTS
%                      {condDir} : root directory where movies are located
%                      {chNames} : cell array of channel names
%                      {markers} : cell array of fluorescent markers
%
% OPTIONS ('Specifier', value)
%                        'Angle' : Scanning angle, in degrees. Default: 32.8
%                    'PixelSize' : (x,y) pixel size, in microns. Default: 0.114
%                'ObjectiveScan' : true|{false}. Specifies objective scan mode.
%                                  In this mode, the raw data is not skewed.
%                'MovieSelector' : selector string for movie folders, i.e., 'cell'
%          'IgnoreEmptyFolders'  : true | {false}; ignores cell folders that do not contain TIFF frames
%
% OUTPUT   data: structure with the fields
%                   .source      : path of the data/movie, location of master channel frames
%                   .channels    : cell array of paths for all channels
%                   .date        : date of the acquisition
%                   .framerate   : frame rate of the movie, in seconds
%                   .imagesize   : dimensions of the movie
%                   .movieLength : length of the movie, in frames
%                   .markers     : cell array of fluorescent marker names
%                   .pixelSize   : pixel size of the CCD, in meters

% Francois Aguet, 02/20/2014
% Xiongtao Ruan, 09/16/2019 make it work noninteractively
% xruan, 07/23/2020 add option for MovieSelector
% xruan (11/09/2020) add option to include experiments without experiment number
% xruan (02/08/2021) add option to use channel patterns to select files for
% a channel (in case of some other images, i.e., multiple channels in the same folder)
% xruan (03/29/2022) add support for user defined dz


function [data] = XR_loadConditionData3D(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.FunctionName = 'loadConditionData3D';
ip.addOptional('condDir', [], @(x) ischar(x) && ~any(strcmpi(x,...
    {'Parameters', 'MovieSelector', 'IgnoreEmptyFolders', 'FrameRate'})));
ip.addOptional('chNames', {}); %{'ch1','ch2','ch3'}
ip.addOptional('markers', {}); %{'alexa488','alexa568','alexa647'},
ip.addParameter('MovieSelector', 'Ex', @ischar);
ip.addParameter('IgnoreExpNumber', false, @islogical); % xruan
ip.addParameter('ChannelPatterns', {}, @iscell); % xruan
ip.addParameter('Angle', 31.5, @isscalar);
ip.addParameter('PixelSize', 0.104, @isscalar);
ip.addParameter('dz', [], @isnumeric);
ip.addParameter('ObjectiveScan', false, @isscalar);
ip.addParameter('Reversed', false, @isscalar);
ip.addParameter('FrameRate', [], @isscalar);
ip.addOptional('maxlevel', 3, @(x) isempty(x) || (isnumeric(x) && abs(round(x))==x));
ip.parse(varargin{:});


condDir = ip.Results.condDir;
chNames = ip.Results.chNames;
markers = ip.Results.markers;
MovieSelector = ip.Results.MovieSelector;
IgnoreExpNumber = ip.Results.IgnoreExpNumber;
FrameRate = ip.Results.FrameRate;
maxlevel = ip.Results.maxlevel;
dz = ip.Results.dz;

if isempty(FrameRate)
    data = XR_loadConditionData(condDir, chNames, markers, ...
        'MovieSelector', MovieSelector, 'IgnoreExpNumber', IgnoreExpNumber, ...
        'maxlevel', maxlevel);
else
    data = XR_loadConditionData(condDir, chNames, markers, ...
        'MovieSelector', MovieSelector, 'IgnoreExpNumber', IgnoreExpNumber, ...
        'FrameRate', FrameRate, 'maxlevel', maxlevel);
end

nd = numel(data);

if isempty(ip.Results.FrameRate)
    % frame rate from file names
    for i = 1:nd
        tmp = regexpi(data(i).framePaths{1}, '\d+(?=msec)', 'match', 'once');
        % in case there's only one z-stack
        if ischar(tmp) == 1
            tmp = {tmp};
            nFrame2Check = 1;
            framePaths = data(i).framePaths;
%             framePathsDS = data(i).framePathsDS;
%             framePathsDSR = data(i).framePathsDSR;
%             maskPaths = data(i).maskPaths;
            data(i).framePaths = cell(1,numel(data(i).channels));
            for c = 1:numel(data(i).channels)
                data(i).framePaths{c}{nFrame2Check} = framePaths{c};
%                 data(i).framePathsDS{c}{nFrame2Check} = framePathsDS{c};
%                 data(i).framePathsDSR{c}{nFrame2Check} = framePathsDSR{c};
%                 data(i).maskPaths{c}{nFrame2Check} = maskPaths{c};
            end
            data(i).imagesize = [data(i).imagesize data(i).movieLength];
            data(i).movieLength = 1;
        end
        % check whether all frames are indeed frames identified by 'msec'
        idx = ~cellfun(@isempty, tmp);
        % check whether last frame is complete
        
        info = imfinfo(data(i).framePaths{end}{end}, 'tif');
        
        if numel(info)~=data(i).imagesize(3)
            idx(end) = false;
            fprintf('Last frame of %s is truncated. Discarded from frame list.\n', getShortPath(data(i)));
        end
        for c = 1:numel(data(i).channels)
            data(i).framePaths{c} = data(i).framePaths{c}(idx);
        end
        data(i).movieLength = sum(idx);
        ms = cellfun(@str2double, tmp(idx));
        data(i).framerate = median(diff(ms))/1000;
    end
    if numel(unique([data.framerate]))~=1
        fprintf('Unequal frame rates detected: ');
        fprintf('%.3f ', [data.framerate]);
        r = input('\nIs this correct? [y/n]: ', 's');
        if strcmpi(r, 'n')
            [data.framerate] = deal(input('Enter correct frame rate: '));
        end
    end
else
    [data.framerate] = deal(ip.Results.FrameRate);
end

% pixel size
[data.pixelSize] = deal(ip.Results.PixelSize);
data = rmfield(data, {'M', 'NA'}); % not relevant for light sheet microscope

%parse z-step
if isempty(dz)
    dz = NaN(1,nd);
    for i = 1:nd
        % z-step is identified by prefix '_zp' or '_z0p' or 'zint_0p'
        dzSelect = '(?<=_zp)\d+|(?<=_z0p)\d+|(?<=zint_0p)\d+';
        idx = regexpi([getDirFromPath(data(i).source) data(i).framePaths{1}{1}], dzSelect, 'match', 'once');
        if ~isempty(idx)
            dz(i) = str2double(['0.' idx]);
        elseif i>1 % if not found, use previous value
            dz(i) = dz(i-1);
        end
        data(i).dz = dz(i);
    end
else
    if numel(dz) == 1
        dz = ones(nd, 1) * dz;
    end
    for i = 1:nd
        data(i).dz = dz(i);
    end
end


if ip.Results.ObjectiveScan
    f = dz./ip.Results.PixelSize;
else
    theta = ip.Results.Angle*pi/180;
    dz0 = sin(theta)*dz; % ~0.25 for dz0 = 0.45
    f = dz0./ip.Results.PixelSize;
end
f = num2cell(f);
[data.zAniso] = f{:};
[data.angle] = deal(ip.Results.Angle);

% time stamps of the frames (at resolution of seconds)
for i = 1:nd
    tmp = dir([data(i).source '*.tif']);
    t = [tmp.datenum];
    t = t-t(1);
%    data(i).timestamps = hour(t)*3600 + minute(t)*60 + second(t);
end

% sample/objective scan flag
objectiveScan = ip.Results.ObjectiveScan;
if numel(objectiveScan)==1 && nd>1
    objectiveScan = repmat(objectiveScan, [1 nd]);
end
objectiveScan = num2cell(objectiveScan);
[data.objectiveScan] = objectiveScan{:};

% reverse flag
reversed = ip.Results.Reversed;
if numel(reversed)==1 && nd>1
    reversed = repmat(reversed, [1 nd]);
end
reversed = num2cell(reversed);
[data.reversed] = reversed{:};

% add paths to de-skewed and rotated data
ext = '.tif';
for i = 1:nd
    % mask paths (loadConditionData assigns masks to Detection/Masks
    maskPath = [data(i).source 'Analysis' filesep 'Masks' filesep];
    data(i).maskPaths = arrayfun(@(i) sprintf('%sdmask_%.3d.tif',maskPath,i), 1:data(i).movieLength, 'unif', 0)';
    
    nCh = numel(data(i).channels);
    data(i).framePathsDS = cell(1,nCh);
    data(i).framePathsDSR = cell(1,nCh);
    for c = 1:nCh
        [~,fname] = cellfun(@fileparts, data(i).framePaths{c}, 'unif', 0);
        data(i).framePathsDS{c} = cellfun(@(f) [data(i).channels{c} 'DS' filesep f ext], fname, 'unif', 0);
        data(i).framePathsDSR{c} = cellfun(@(f) [data(i).channels{c} 'DSR' filesep f ext], fname, 'unif', 0);
    end
end

