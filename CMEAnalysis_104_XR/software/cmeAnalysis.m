%[res, data] = cmeAnalysis(varargin) performs the analysis of clathrin-coated pit dynamics on a set of movies.
% The analysis comprises detection, tracking, and selection of bona fide CCP structures, and generates lifetime
% distribution and intensity cohort plots.
%
% Inputs (optional):
%           data : data structure returned by loadConditionData()
%
% Options:
%            'PlotAll' : true|[false] displays intermediary processing steps
%        'GaussianPSF' : ['model']|'data' toggles between a model-based or data-based estimation of the Gaussian PSF s.d.
%  'TrackingGapLength' : value defines the maximum number of consecutive missed frames in a trajectory. Default: 2
%     'TrackingRadius' : [minRadius maxRadius] search radii for frame-to-frame linking and gap closing. Default: [3 6]
%       'ChannelNames' : cell array of channel names (e.g., {'EGFP-CLCa'}). Default: fluorophore names
%          'Overwrite' : true|[false] overwrites results of previous analyses. 
%
% Outputs:
%            res : analysis results. Lifetime analysis in 'lftRes' field; intensity cohorts in 'cohorts' field
%           data : structure returned by loadConditionData()
%
% The function will ask for acquisition parameters and subsequently for the data location. The following acquisition
% parameters are required for the Gaussian point spread function (PSF) model used for CCP detection:
% numerical aperture (NA) and magnification of the objective used, and the physical pixel size of the camera (in �m).
%
%
% Notes:
%
% Gaussian PSF model parameters:
% ------------------------------
% The algorithm used for CCP detection estimates the Gaussian PSF s.d. based on an accurate model of a TIRFM PSF
% using the objective NA and magnification, camera pixel size, and fluorophores emission wavelength. In the case
% of misestimation of any of these parameters or non-ideal TIR conditions, the resulting PSF model may be sub-optimal
% for detection, and a data-derived parameterization of this model may be preferable.
%
% Tracking parameters:
% --------------------
% The two most important and sensitive parameters are the maximum number of consecutive missed detection, or ?gaps?,
% in a trajectory, and the search radii for frame-to-frame linking of the detections and for gap closing. These
% parameters are sensitive to the imaging frame rate and should be adjusted accordingly.
% The default values are recommended for data acquired at 0.5 frames/sec
%
% Comparing conditions:
% ---------------------
% cmeAnalysis() must be run separately on groups of movies from different experimental conditions
% (i.e., control vs. perturbation), since the automatic thresholds for identifying bona fide CCPs must be determined
% on control data. For such comparisons, first run, i.e.,
%  >> [resCtrl, dataCtrl] = cmeAnalysis;
% followed by
%  >> [resPert, dataPert] = cmeAnalysis('ControlData', resCtrl);
% In the first run, select the parent directory of the control data. In the second run, select the parent directory
% of the perturbation condition.
%
% Parallelization:
% ----------------
% This function takes advantage of Matlab?s parallelization capabilities. To enable this, enter
% >> matlabpool
% in the command prompt.
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

% Francois Aguet (last mod. 05/29/2013)

function [res, data] = cmeAnalysis(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('data', [], @isstruct);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('GaussianPSF', 'data', @(x) any(strcmpi(x, {'data', 'model'})));
ip.addParameter('Sigma', []);
ip.addParameter('TrackingRadius', [3 6], @(x) numel(x)==2);
ip.addParameter('TrackingGapLength', 2, @(x) numel(x)==1);
ip.addParameter('Parameters', [], @(x) numel(x)==3);
ip.addParameter('ControlData', [], @isstruct);
ip.addParameter('PlotAll', false, @islogical);
ip.addParameter('SlaveAmplitudeRatio', 0, @isnumeric);
ip.addParameter('ChannelNames', []);
ip.addParameter('FirstNFrames', [], @isposint);
ip.addParameter('MaxIntensityThreshold', [], @isscalar);
ip.addParameter('CompareHighLowSNR', false, @islogical);
ip.addParameter('DisplayMode', 'screen', @(x) any(strcmpi(x, {'print', 'screen'})));
ip.parse(varargin{:});
data = ip.Results.data;

overwrite = ip.Results.Overwrite;
if numel(overwrite)==1
    overwrite = repmat(overwrite, [1 4]);
end

if isempty(data)
    parameters = ip.Results.Parameters;
    if isempty(parameters)
        parameters = zeros(1,3);
        parameters(1) = input('Enter the N.A. of the objective: ');
        parameters(2) = input('Enter the magnification of the objective: ');
        parameters(3) = 1e-6*input('Enter the camera pixel size, in [um]: ');
    end
    data = loadConditionData('Parameters', parameters);
end

%-------------------------------------------------------------------------------
% 1) Detection
%-------------------------------------------------------------------------------
% Calculate cell masks and ask for user validation
getCellMask(data, 'Overwrite', overwrite(1), 'Validate', true);

runDetection(data, 'SigmaSource', ip.Results.GaussianPSF,...
    'Sigma', ip.Results.Sigma, 'Overwrite', overwrite(2));
[psnr, cmap] = plotPSNRDistribution(data, 'Pool', false);

%-------------------------------------------------------------------------------
% 2) Tracking & track processing
%-------------------------------------------------------------------------------
settings = loadTrackSettings('Radius', ip.Results.TrackingRadius, 'MaxGapLength', ip.Results.TrackingGapLength);
runTracking(data, settings, 'Overwrite', overwrite(3));

runTrackProcessing(data, 'Overwrite', overwrite(4),'SlaveAmplitudeRatio',ip.Results.SlaveAmplitudeRatio);

%-------------------------------------------------------------------------------
% 3) Analysis
%-------------------------------------------------------------------------------
chNames = ip.Results.ChannelNames;
if isempty(chNames)
    chNames = data(1).markers;
end

if ip.Results.PlotAll
    display = 'on';
else
    display = 'off';
end
lopts = {'Display', display, 'RemoveOutliers', true, 'Colormap', cmap, 'DisplayMode', ip.Results.DisplayMode,...
    'SlaveNames', chNames(2:end), 'FirstNFrames', ip.Results.FirstNFrames, 'Overwrite', any(overwrite)};
if isempty(ip.Results.ControlData)
    res.lftRes = runLifetimeAnalysis(data, lopts{:},...
        'MaxIntensityThreshold', ip.Results.MaxIntensityThreshold);
else
    res.lftRes = runLifetimeAnalysis(data, lopts{:},...
        'MaxIntensityThreshold', ip.Results.ControlData.lftRes.MaxIntensityThreshold);
end

% Graphical output
if ~ip.Results.CompareHighLowSNR
    plotLifetimes(res.lftRes, 'DisplayMode', ip.Results.DisplayMode, 'PlotAll', false,...
        'SlaveNames', chNames(2:end), 'SingleChannel', numel(data(1).channels)==1);
else  
%     nd = numel(data);
%     [~,idx] = sort(cellfun(@median, psnr), 'ascend');
%     id = idx(floor(nd/2)+1:end);
%     lftRes = runLifetimeAnalysis(data(id), lopts{:},...
%         'MaxIntensityThreshold', res.lftRes.MaxIntensityThreshold, 'Colormap', cmap(id,:));
%     plotLifetimes(lftRes, 'DisplayMode', ip.Results.DisplayMode, 'PlotAll', false,...
%         'SlaveNames', chNames(2:end), 'SingleChannel', numel(data(1).channels)==1);
%     
%     id = idx(1:floor(nd/2));
%     lftRes = runLifetimeAnalysis(data(id), lopts{:},...
%         'MaxIntensityThreshold', res.lftRes.MaxIntensityThreshold, 'Colormap', cmap(id,:));
%     plotLifetimes(lftRes, 'DisplayMode', ip.Results.DisplayMode, 'PlotAll', false,...
%         'SlaveNames', chNames(2:end), 'SingleChannel', numel(data(1).channels)==1);
end
    

res.cohorts = plotIntensityCohorts(data, 'MaxIntensityThreshold', res.lftRes.MaxIntensityThreshold,...
    'ShowBackground', false, 'DisplayMode', 'screen', 'ScaleSlaveChannel', false,...
    'ShowPct', false, 'SlaveName', chNames(2:end));
