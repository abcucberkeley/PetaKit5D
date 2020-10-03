function [data] = loadConditionData(varargin)
% loadConditionData loads the relevant information for all the data
% available for a specific data condition; this requires a specific
% directory structure and nomenclature (see below)
%
% SYNOPSIS [data] = loadConditionData()
%
% INPUT                   {condDir} : root directory where movies are located
%                         {chNames} : cell array of channel names
%                         {markers} : cell array of fluorescent markers
%             {'Parameters', value} : vector of microscope parameters: [NA M PixelSize]
%                                     NA: numerical aperture; M: magnification;
%                                     PixelSize: camera pixel size. Default [1.49 100 6.7e-6]
%          {'MovieSelector', value} : selector string for movie folders, i.e., 'cell'
%     {'IgnoreEmptyFolders', value} : true | {false}; ignores cell folders that do not contain TIFF frames
%
% OUTPUT   data: structure with the fields
%                   .source      : path of the data/movie, location of master channel frames
%                   .channels    : cell array of paths for all channels
%                   .date        : date of the acquisition
%                   .framerate   : frame rate of the movie, in seconds
%                   .imagesize   : dimensions of the movie
%                   .movieLength : length of the movie, in frames
%                   .markers     : cell array of fluorescent marker names
%                   .NA          : numerical aperture of the objective
%                   .M           : magnification of the objective
%                   .pixelSize   : pixel size of the CCD, in meters
%
%
% Francois Aguet, October 2010 (last modified: 05/27/2011)

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.FunctionName = 'loadConditionData';
ip.addOptional('condDir', [], @(x) ischar(x) && ~any(strcmpi(x,...
    {'Parameters', 'MovieSelector', 'IgnoreEmptyFolders', 'FrameRate'})));
ip.addOptional('chNames', [], @iscell);
ip.addOptional('markers', [], @iscell);
ip.addParamValue('Parameters', [1.49 100 6.45e-6], @(x) numel(x)==3);
ip.addParamValue('MovieSelector', 'cell', @ischar);
ip.addParamValue('StrictSelector', false, @islogical);
ip.addParamValue('IgnoreEmptyFolders', false, @islogical);
ip.addParamValue('FrameRate', [], @isscalar);
ip.parse(varargin{:});

condDir = ip.Results.condDir;
chNames = ip.Results.chNames;
markers = ip.Results.markers;
parameters = ip.Results.Parameters;

if isempty(condDir)
    condDir = uigetdir(pwd, 'Select the ''condition'' folder:');
    if condDir==0
        data = [];
        return
    else
        condDir = [condDir filesep];
    end
end
if ~strcmp(condDir(end), filesep)
    condDir = [condDir filesep];
end

fprintf('Root directory: %s\n', condDir);

cellPath = recursiveDir(condDir, 2);

% cell directories:
[cellDirs, cellPar] = cellfun(@(i) getDirFromPath(i), cellPath, 'unif', 0);

idx = regexpi(cellDirs, ip.Results.MovieSelector, 'once');
[idx{cellfun(@isempty, idx)}] = deal(NaN);
idx = vertcat(idx{:});
if ip.Results.StrictSelector
    idx(idx~=1) = NaN;
end
idx = ~isnan(idx);
cellPath = cellPath(idx);

if isempty(cellPath)
    error(['No movies found in: ' condDir]);
end

% sort by cell number (for each experiment)
[ui,~,ai] = unique(cellPar(idx));
for k = 1:numel(ui)
    cellPath(ai==k) = sortStringsByToken(cellPath(ai==k), ip.Results.MovieSelector, 'post');
end

nCells = length(cellPath);

data(1:nCells) = struct('source', [], 'channels', [], 'date', [], 'framerate', [],...
    'imagesize', [], 'movieLength', [], 'markers', [],...
    'framePaths', [], 'maskPaths', []);


% Load/determine channel names
if isempty(chNames)
    nCh = input('Enter the number of channels: ');
    chNames = cell(1,nCh);
    chPath = uigetdir(cellPath{1}, 'Select first (master) channel:');
    if chPath==0
        fprintf(2, 'LoadConditionData: cancelled.\n');
        return
    else
        chPath = [chPath filesep];
    end
    chNames{1} = chPath(length(cellPath{1})+1:end-1);
    for c = 2:nCh
        chPath = uigetdir(cellPath{1}, ['Select channel #' num2str(c) ':']);
        if chPath==0
            fprintf(2, 'LoadConditionData: cancelled.\n');
            return
        else
            chPath = [chPath filesep]; %#ok<AGROW>
        end
        chNames{c} = chPath(length(cellPath{1})+1:end-1);
    end
else
    nCh = length(chNames);
end
if isempty(markers)
    markers = cell(1,nCh);
    for c = 1:nCh
        markers{c} = input(['Enter the fluorescent marker for channel ' num2str(c) ': '], 's');
    end
end
% validate markers
s = getFluorPropStruct();
s = {s.name};
for c = 1:nCh
    if isnumeric(markers{c}) && markers{c}>=350 && markers{c}<=750
        markers{c} = 1e-9*markers{c};
    elseif ~any(strcmpi(markers{c}, s));
        markers{c} = 1e-9*input(['Marker ''' markers{c} ''' not recognized, enter max. emission wavelength in [nm]: ']);
    end
end
for c = 1:nCh
    fprintf('Channel %d name: "%s"\n', c, chNames{c});
end

channels = cell(1,nCh);
for k = 1:nCells
    
    % detect date
    data(k).date = cell2mat(regexp(cellPath{k}, '\d{6}+', 'match'));
    if isempty(date)
        data(k).date = '000000';
    end
    
    % detect frame
    fr = regexp(cellPath{k}, '_(\d+)?(\.)?\d+s', 'match');
    if ~isempty(ip.Results.FrameRate)
        data(k).framerate = ip.Results.FrameRate;
    elseif ~isempty(fr)
        data(k).framerate = str2double(fr{1}(2:end-1));
    else
        fr = regexp(cellPath{k}, '_\d+ms', 'match');
        if ~isempty(fr)
            data(k).framerate = str2double(fr{1}(2:end-2))/1000;
        else
            data(k).framerate = 2; % default: 2s
        end
    end
    
    % assign full channel paths
    framePaths = cell(1,nCh);
    for c = 1:nCh
        if ~isempty(chNames{c})
            channels{c} = [cellPath{k} chNames{c} filesep];
        else
            channels{c} = cellPath{k};
        end
        if ~(exist(channels{c}, 'dir')==7)
            channels{c} = uigetdir(cellPath{k}, ['Select channel #' num2str(c) ':']);
            if channels{c}==0
                fprintf(2, 'LoadConditionData: cancelled.\n');
                return
            else
                channels{c} = [channels{c} filesep];
            end
        end
 
        % list directory contents and select TIFFs
        tmp = dir(channels{c});
        tmp = tmp(cellfun(@(i) ~strcmpi(i(1), '.'), {tmp.name}));
        tmp = {tmp(~[tmp.isdir]).name}';
        tmp = tmp(~cellfun(@isempty, regexpi(tmp, '\.tif|\.stk')));
        % sort files in case leading zeros are missing
        idx = regexp(tmp','\d+(?=\.)', 'match', 'once');
        if numel(unique(cellfun(@numel, idx)))~=1
            idx = str2double(idx);
            [~,idx] = sort(idx);
            tmp = tmp(idx);
        end
        framePaths{c} = strcat(channels{c}, tmp);
    end
    data(k).channels = channels;
    data(k).source = channels{1}; % master channel default
    
    % only store frame paths if frames for all channels are found
    if all(~cellfun(@isempty, framePaths))
        if numel(framePaths{1})==1
            data(k).framePaths = cellfun(@(i) i{1}, framePaths, 'unif', 0);
            info = imfinfo(data(k).framePaths{1});
            data(k).imagesize = [info(1).Height info(1).Width];
            data(k).movieLength = numel(info);
        else
            data(k).framePaths = framePaths;
            info = imfinfo(framePaths{1}{1});
            data(k).imagesize = [info(1).Height info(1).Width];
            if numel(info)>1
                data(k).imagesize = [data(k).imagesize numel(info)];
            end
            data(k).movieLength = length(framePaths{1});
        end
    elseif exist([data(k).source 'Detection' filesep 'detection_v2.mat'], 'file')==2
        d = load([data(k).source 'Detection' filesep 'detection_v2.mat']);
        if isfield(d, 'data')
            data(k).framePaths = d.data.framePaths;
            data(k).imagesize = d.data.imagesize;
            data(k).movieLength = d.data.movieLength;
        end
    end
    
    % generate mask paths
    if iscell(data(k).framePaths{1})
        maskPath = [data(k).source 'Detection' filesep 'Masks' filesep];
        nf = data(k).movieLength;
        tmp = num2str((1:nf)');
        tmp(tmp==' ') = '0';
        tmp = [repmat([maskPath 'dmask_'], [nf 1]) tmp repmat('.tif', [nf 1])]; %#ok<AGROW>
        data(k).maskPaths = mat2cell(tmp, ones(1,nf), size(tmp,2));
    else
        % save masks as stack
        data(k).maskPaths = [data(k).source 'Detection' filesep 'dmasks.tif'];
    end    
    
    data(k).markers = markers;
    data(k).NA = parameters(1);
    data(k).M = parameters(2);
    data(k).pixelSize = parameters(3);
    
    fprintf('Loaded: %s\n', cellPath{k});
end

% remove empty cell folders
if ip.Results.IgnoreEmptyFolders
    data(arrayfun(@(x) isempty(x.framePaths), data)) = [];
end

if numel(data)==0
    data = [];
end
