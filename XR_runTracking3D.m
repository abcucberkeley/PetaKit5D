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
%            'FileName' : tracking result filename.      
%       'DetectionFile' : detection filename as input.
%         'ResultsPath' : result directory. 
%      'ParallelMethod' : parallel computin method method. 
%                           'parfor' or 'lock_distributed' (default, and only available for Mac and Linux). 
%
% Example: runTracking(data, loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2), 'Overwrite', true) ;

% Francois Aguet, May 2010 (last modified 05/28/2013)
% Xiongtao Ruan, Sept 2019 skip failed cases
% Xiongtao Ruan, Jan 2020 add lock based parallel computing method

function [] = XR_runTracking3D(data, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('settings', [], @isstruct);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('FileName', 'trackedFeatures', @ischar); % default of the tracker
ip.addParameter('Frames', [], @isvector);
ip.addParameter('DetectionFile', 'Detection3D.mat', @ischar);
ip.addParameter('ResultsPath', arrayfun(@(i) [i.source 'Analysis' filesep], data, 'unif', 0), @iscell);
ip.addParameter('ParallelMethod', 'lock_distributed', @ischar); % parallization method

ip.parse(data, varargin{:});
overwrite = ip.Results.Overwrite;
fileName = ip.Results.FileName;
idx = regexpi(fileName, '.mat');
if ~isempty(idx)
    fileName = fileName(1:idx-1);
end
frames = ip.Results.Frames;
detectionFile = ip.Results.DetectionFile;

ParallelMethod = ip.Results.ParallelMethod;
if ispc
    ParallelMethod = 'parfor';
end
    
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

ResultsPath = ip.Results.ResultsPath;
FileName = ip.Results.FileName;

% Run tracker on each data set
switch ParallelMethod
    case 'parfor'
        parfor i = 1:length(data)
            if ~(exist([ResultsPath{i} fileName], 'file')==2) || overwrite
                fprintf('Running tracker on %s\n', getShortPath(data(i)));
                main(data(i), settings, fileName, detectionFile, ResultsPath{i});
            else
                fprintf('Tracking has already been run for %s\n', getShortPath(data(i)));
            end
        end
    case 'lock_distributed'
        for i = 1 : numel(data)
            chunk_lock_clean(ResultsPath{i}, 0.1 * data(i).movieLength);
            
            result_fullname = [ResultsPath{i} FileName];
            temp_filename = sprintf('%s.tmp', result_fullname);
            if exist(result_fullname, 'file')
                fprintf('Tracking has already been run for %s\n', getShortPath(data(i)));
                continue;
            end
            if exist(temp_filename, 'file')
                fprintf('Running tracker on %s\n', getShortPath(data(i)));
                continue;                
            else
                fclose(fopen(temp_filename, 'w'));
            end
            
            main(data(i), settings, fileName, detectionFile, ResultsPath{i});
            
            if exist(result_fullname, 'file') && exist(temp_filename, 'file')
                delete(temp_filename);
            end
        end
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

% try 
% trackCloseGapsKalmanSparse(movieInfo, settings.costMatrices, settings.gapCloseParam,...
%     settings.kalmanFunctions, 3, saveResults, 1);
% catch ME
%     disp(ME)
%     return;
% end
end