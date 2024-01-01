%runTracking(data, varargin) tracks CCPs in the movies passed with the 'data' structure.
% This function generates a list of tracks in 'Tracking/trackedFeatures.mat' for each
% data set.
%
% Inputs   
%                  data : list of movies, using the structure returned by loadConditionData.m
%            {settings} : tracker settings. Default: 
%                         loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2);    
%
% Options ('specifier', value)
%           'Overwrite' : true|{false}. Overwrite previous tracking result.
%              'Frames' : Index array of frames to track (i.e., for downsampling). Default: all frames. 
%  'DownsamplingFactor' : Integer downsampling factor. Default: none.
%
%
% Example: runTracking(data, loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2), 'Overwrite', true) ;

% Francois Aguet, May 2010 (last modified 05/28/2013)

function [] = runTracking3D(data, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('settings', [], @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('FileName', 'trackedFeatures', @ischar); % default of the tracker
ip.addParamValue('Frames', [], @isvector);
ip.addParamValue('DetectionFile', 'Detection3D.mat', @ischar);
ip.addParamValue('ResultsPath', arrayfun(@(i) [i.source 'Analysis' filesep], data, 'unif', 0), @iscell);
ip.parse(data, varargin{:});
overwrite = ip.Results.Overwrite;
fileName = ip.Results.FileName;
idx = regexpi(fileName, '.mat');
if ~isempty(idx)
    fileName = fileName(1:idx-1);
end
frames = ip.Results.Frames;
detectionFile = ip.Results.DetectionFile;

% Determine file name
if ~isempty(frames)
    fileName = [fileName '_customFrames(' num2str(frames(1)) '_' num2str(frames(end)) ')'];
end

fileName = [fileName '.mat'];

% Load tracker settings
settings = ip.Results.settings;
if isempty(settings)
   settings = loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2);
end

% Run tracker on each data set
parfor i = 1:length(data)
    if ~(exist([ip.Results.ResultsPath{i} fileName], 'file')==2) || overwrite
        fprintf('Running tracker on %s\n', getShortPath(data(i)));
        main(data(i), settings, fileName, detectionFile, ip.Results.ResultsPath{i});
    else
        fprintf('Tracking has already been run for %s\n', getShortPath(data(i)));
    end
end


function main(data, settings, fileName, detectionFile, resultsPath)

dfile = [resultsPath detectionFile];
if exist(dfile, 'file')==2
    dfile = load(dfile);
    movieInfo = dfile.frameInfo;
else
    fprintf(['runTracking: no detection data found for ' getShortPath(data) '\n']);
    return;
end

[~,~] = mkdir(resultsPath);
saveResults.dir = resultsPath;
if ~isempty(fileName)
    saveResults.filename = fileName;
end
trackCloseGapsKalmanSparse(movieInfo, settings.costMatrices, settings.gapCloseParam,...
    settings.kalmanFunctions, 3, saveResults, 1);
