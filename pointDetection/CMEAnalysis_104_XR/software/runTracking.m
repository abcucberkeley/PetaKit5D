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
%
% Copyright (C) 2017, Danuser Lab - UTSouthwestern 
%
% This file is part of CMEAnalysis_Package.
% 
% CMEAnalysis_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% CMEAnalysis_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with CMEAnalysis_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Francois Aguet, May 2010 (last modified 05/28/2013)

function [] = runTracking(data, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('settings', [], @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('FileName', 'trackedFeatures', @ischar); % default of the tracker
ip.addParamValue('Frames', [], @isvector);
ip.addParamValue('DownsamplingFactor', [], @isscalar);
ip.addParamValue('DetectionFile', 'detection_v2.mat', @ischar);
ip.parse(data, varargin{:});
overwrite = ip.Results.Overwrite;
fileName = ip.Results.FileName;
idx = regexpi(fileName, '.mat');
if ~isempty(idx)
    fileName = fileName(1:idx-1);
end
frames = ip.Results.Frames;
dsfactor = ip.Results.DownsamplingFactor;
detectionFile = ip.Results.DetectionFile;

% Determine file name
if ~isempty(frames)
    fileName = [fileName '_customFrames(' num2str(frames(1)) '_' num2str(frames(end)) ')'];
end
if isempty(frames) && ~isempty(dsfactor)
    frames = 1:dsfactor:data.movieLength;
    fileName = [fileName '_' num2str(data.framerate*dsfactor) 's'];
end
fileName = [fileName '.mat'];

% Load tracker settings
settings = ip.Results.settings;
if isempty(settings)
   settings = loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2);
end

% Run tracker on each data set
parfor i = 1:length(data)
    if ~(exist([data(i).source 'Tracking' filesep 'trackedFeatures.mat'], 'file')==2) || overwrite  % Tracking is performed on master (in data.source) only
        fprintf('Running tracker on %s\n', getShortPath(data(i)));
        main(data(i), settings, fileName, detectionFile, frames);
    else
        fprintf('Tracking has already been run for %s\n', getShortPath(data(i)));
    end
end


function main(data, settings, fileName, detectionFile, frames)

dfile = [data.source 'Detection' filesep detectionFile];
if exist(dfile, 'file')==2
    dfile = load(dfile);
    movieInfo = dfile.frameInfo;
    if ~isempty(frames)
        movieInfo = movieInfo(frames);
    end
else
    fprintf(['runTracking: no detection data found for ' getShortPath(data) '\n']);
    return;
end

[~,~] = mkdir([data.source 'Tracking']);
saveResults.dir = [data.source 'Tracking' filesep];
if ~isempty(fileName)
    saveResults.filename = fileName;
end
trackCloseGapsKalmanSparse(movieInfo, settings.costMatrices, settings.gapCloseParam, settings.kalmanFunctions, 2, saveResults, 1);
