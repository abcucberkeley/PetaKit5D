%cmeAnalysis3D performs the analysis of clathrin-coated pit dynamics on data
% generated with a light sheet microscope.
% The analysis comprises detection, tracking, and selection CCP structures.
% The graphical output includes lifetime distribution and intensity cohort plots.
%
% Inputs (optional):
%           data : data structure returned by loadConditionData()
%
% Options:
%  'TrackingGapLength' : value defines the maximum number of consecutive missed frames in a trajectory. Default: 2
%     'TrackingRadius' : [minRadius maxRadius] search radii for frame-to-frame linking and gap closing. Default: [3 6]
%          'Overwrite' : true|[false] overwrites results of previous analyses. 

%        'RunAnalysis' : {true}|false. Disable to run de-skewing only.


% Author: Francois Aguet (2014)

function cmeAnalysis3D(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('data', [], @isstruct); % data structure from loadConditionData
ip.addOptional('apath', []); % analysis path
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('TrackingRadius', [3 6], @(x) numel(x)==2);
ip.addParamValue('TrackingGapLength', 2, @(x) numel(x)==1);
ip.addParamValue('Parameters', [], @(x) numel(x)==3);
% ip.addParamValue('PlotAll', false, @islogical);
% ip.addParamValue('ChannelNames', []);
% ip.addParamValue('FirstNFrames', [], @isposint);
% ip.addParamValue('DisplayMode', 'screen', @(x) any(strcmpi(x, {'print', 'screen'})));
ip.addParamValue('SkewAngle', 31.5, @isscalar);
ip.addParamValue('Sigma', []);
ip.addParamValue('WindowSize', []);
ip.addParamValue('Cutoff_f', 3, @isposint);
ip.addParamValue('Crop', true, @islogical);
ip.addParamValue('Rotate', false, @islogical);
ip.addParamValue('Buffer', [3 3]);
ip.addParamValue('BufferAll', false, @islogical);
ip.addParamValue('CheckMismatch', false, @islogical);
ip.addParamValue('RunAnalysis', true, @islogical);
ip.addParamValue('LoadSettings', true, @islogical);
ip.parse(varargin{:});
data = ip.Results.data;
apath = ip.Results.apath;

overwrite = ip.Results.Overwrite;
if numel(overwrite)==1
    overwrite = repmat(overwrite, [1 4]);
end

sigma = ip.Results.Sigma;
if isempty(sigma)
    fprintf('Enter the Gaussian s.d. for detection. If the z-value is different, \n');
    sigma = input('enter the two values as, e.g., "[1.5 1.3]": ');
end

% settings -> input options instead (check all, vs single)
% improve help for this function

%-------------------------------------------------------------------------------
% 1) Parse data, determine output structure
%-------------------------------------------------------------------------------
if isempty(data)
    data = loadConditionData3D();
end

if isempty(apath) % store results locally in cell/Analysis directory
    apath = arrayfun(@(i) [i.source 'Analysis' filesep], data, 'unif', 0);
else
    % expand results path for all data sets
    apath = arrayfun(@(i) [apath getShortPath(i,3) filesep 'Analysis' filesep], data, 'unif', 0);
end
[~,~] = cellfun(@mkdir, apath, 'unif', 0);

data = deskewData(data, 'Overwrite', overwrite(1), 'SkewAngle', ip.Results.SkewAngle,...
    'Rotate', ip.Results.Rotate, 'Crop', ip.Results.Crop, 'LoadSettings', ip.Results.LoadSettings);


if ip.Results.RunAnalysis
    
    %-------------------------------------------------------------------------------
    % 3) Detection
    %-------------------------------------------------------------------------------
    runDetection3D(data, 'Sigma', sigma, 'Overwrite', overwrite(2),...
        'WindowSize', ip.Results.WindowSize, 'ResultsPath', apath);
    
    %-------------------------------------------------------------------------------
    % 4) Tracking
    %-------------------------------------------------------------------------------
    runTracking3D(data, loadTrackSettings('Radius', ip.Results.TrackingRadius,...
        'MaxGapLength', ip.Results.TrackingGapLength),...
        'FileName', 'trackedFeatures.mat', 'Overwrite', overwrite(3), 'ResultsPath', apath);
    
    runTrackProcessing3D(data, 'Overwrite', overwrite(4),...
        'TrackerOutput', 'trackedFeatures.mat', 'FileName', 'ProcessedTracks.mat',...
        'Buffer', ip.Results.Buffer, 'BufferAll', ip.Results.BufferAll,...
        'WindowSize', ip.Results.WindowSize, 'ResultsPath', apath);
    
    %-------------------------------------------------------------------------------
    % 5) Rotate tracks
    %-------------------------------------------------------------------------------
    rotateTracks3D(data, 'Overwrite', overwrite(4), 'Crop', ip.Results.Crop);
    
    %-------------------------------------------------------------------------------
    % 6) Analysis functions
    %-------------------------------------------------------------------------------
    getLifetimeData(data, 'Overwrite', any(overwrite),...
        'Scale', false, 'Cutoff_f', ip.Results.Cutoff_f, 'ReturnValidOnly', true,...
        'ExcludeVisitors', false, 'Mask', false, 'AnalysisPath', 'Analysis');
    
    %lopts = {'Display', display, 'RemoveOutliers', true, 'Colormap', cmap, 'DisplayMode', ip.Results.DisplayMode,...
    %    'SlaveNames', chNames(2:end), 'FirstNFrames', ip.Results.FirstNFrames};
    lopts = {'Display', 'on', 'ExcludeVisitors', false};
    res.lftRes = runLifetimeAnalysis(data, [1 20 40 60 80 100 120], lopts{:});
    
    res.cohorts = plotIntensityCohorts(data, 'MaxIntensityThreshold', res.lftRes.MaxIntensityThreshold,...
        'ShowBackground', false, 'DisplayMode', 'screen', 'ScaleSlaveChannel', false,...
        'ShowPct', false);%, 'SlaveName', chNames(2:end));
end
