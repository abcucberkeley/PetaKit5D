% useful only to perform deskew, detection and tracking
% based on cmeAnalysis
% Gokul Upadhyayula, December 2014
% updated Aug 2017 - enabled options for mixture model fitting
% Xiongtao Ruan, January 2019
% 
% 
% Inputs (optional):
%       data : structure returned by loadConditionData()
% 
% Options ('specifier', value):
%               'apath' : analysis root path
%               'aname' : the name of directory that stores the results (under the primary channel directory). 
%           'Overwrite' : overwrite if the result exists, defalt: false
%      'ParallelMethod' : use parfor or lock-based distributed computing. 
%                           'parfor' or 'lock-distributed' (only available in Mac and Linux, and by default)
%           'SkewAngle' : skew angle in microscopy
%               'Sigma' : sigma for the psf, N X 2 for N channels with form [sigma_xy, sigma_z] in each row. 
%          'WindowSize' : window size for detection, default: [], and it is ceiling of 4 folds of sigma. 
%                'Mode' : estimate x, y, z, A (intensity) and c (background).  
%         'FitMixtures' : fit mixture of Gaussians, default: false
%         'MaxMixtures' : number of mixture if choosing to fit mixtures. 
%     'DetectionMethod' : detection method chosen from 'original', 'speedup', 'lowSNR', 'adaptive'. 
%                           'original': original method in cmeAnalysis3D package
%                           'speedup': the speedup of the original method
%                           'lowSNR' (default): the low-SNR detection methods, 
%                           'adaptive' is for the adative method between speedup and lowSNR 
%'BackgroundCorrection' : adjust illumination across z-stacks. Default: false
%               'Track' : track after detection, default true
%        'RotateTracks' : rotate tracks, default true  
%      'TrackingRadius' : radius to link points
%   'TrackingGapLength' : number of frames for gap filling. 
%              'Buffer' : Length of buffer readout before/after each track
%           'BufferAll' : Buffer all tracks, default false
%     'sCMOSCameraFlip' : whether to flip image for sCMOS camera, default false
%                'Crop' : crop an ROI for deskew, default false
%              'Rotate' : rotate image to correct deskew angle, default false
%      'LLFFCorrection' : flat-field correction, default false
%          'LowerLimit' : Lower bound for mask in flat-field correction, default 0.4
%          'LSImageCh1' : path of light-sheet image for channel 1
%          'LSImageCh2' : path of light-sheet image for channel 2, leave as '' if there is no channel 2
%          'LSImageCh3' : path of light-sheet image for channel 3, leave as '' if there is no channel 3
%          'LSImageCh4' : path of light-sheet image for channel 4, leave as '' if there is no channel 4
%       'BackgroundCh1' : path of background image for channel 1
%       'BackgroundCh2' : path of background image for channel 2, leave as '' if there is no channel 2
%       'BackgroundCh3' : path of background image for channel 3, leave as '' if there is no channel 3
%       'BackgroundCh4' : path of background image for channel 4, leave as '' if there is no channel 4
%  'XZoffsetCorrection' : shift in X and/or Z directions to adjust potential
%                           offsets across different channels. Default: false. 


function XR_runDetTrack3d(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addOptional('data', [], @isstruct); % data structure from loadConditionData
ip.addOptional('apath', []); % analysis root path
ip.addOptional('aname', []); % analysis directory name
ip.addParameter('Overwrite', false, @islogical); % 1. deskew; 2. detection; 3. tracking; 4. track processing
ip.addParameter('ParallelMethod', 'lock_distributed', @ischar); % use parfor or lock-based distributed computing (default). 

% detection options
ip.addParameter('SkewAngle', 31.5, @isscalar);
ip.addParameter('Sigma', [1.26, 1.32; 1.41, 1.38;  1.58, 1.608;]); % for 0.5sample scan [1.41, 1.38; 1.26, 1.308; 1.58, 1.608;] %for z 0.6 [1.58, 1.34;1.41, 1.15; 1.26, 1.09;] based on zAniso = 2.5 sampling; hex [1.23, 1.88; 1.42, 2.2] // sq:[1.23, 2.54; 1.42, 3.15;]
ip.addParameter('WindowSize', []);
ip.addParameter('Mode', 'xyzAc', @ischar); % do not change unless necessary
ip.addParameter('FitMixtures', false, @islogical); % fairly robust
ip.addParameter('MaxMixtures', 2, @(x) numel(x)==1 && x>0 && round(x)==x);
ip.addParameter('DetectionMethod', 'update_3', @ischar);
ip.addParameter('BackgroundCorrection', false, @islogical); % background correction

% tracking options
ip.addParameter('Track', true, @islogical);
ip.addParameter('RotateTracks', true, @islogical);
ip.addParameter('TrackingRadius', [3 6], @(x) numel(x)==2);
ip.addParameter('TrackingGapLength', 2, @(x) numel(x)==1);
ip.addParameter('Buffer', [3 3]);
ip.addParameter('BufferAll', false, @islogical);

% deskew options
ip.addParameter('sCMOSCameraFlip', false, @islogical); % necessary when frame rotation is off during SLICE acq.
ip.addParameter('Crop', false, @islogical);
ip.addParameter('Rotate', false, @islogical);

% light sheet flat field corrections options below
ip.addParameter('LLFFCorrection', false, @islogical);
ip.addParameter('LowerLimit', 0.4, @isnumeric); % this value is the lowest
ip.addParameter('LSImageCh1', '' , @isstr);
ip.addParameter('LSImageCh2', '' , @isstr);
ip.addParameter('LSImageCh3', '' , @isstr);
ip.addParameter('LSImageCh4', '' , @isstr);
ip.addParameter('BackgroundCh1', '' , @isstr);
ip.addParameter('BackgroundCh2', '' , @isstr);
ip.addParameter('BackgroundCh3', '' , @isstr);
ip.addParameter('BackgroundCh4', '' , @isstr);

% z-offset correction option
ip.addParameter('XZoffsetCorrection', false , @islogical);


ip.parse(varargin{:});
data = ip.Results.data;
apath = ip.Results.apath;
aname = ip.Results.aname;
% xruan 01/01/2020 if it is pc use parfor as default parallel computing method
if ispc
    ParallelMethod = 'parfor';
else
    ParallelMethod = ip.Results.ParallelMethod;    
end

pr = ip.Results;
overwrite = ip.Results.Overwrite;
if numel(overwrite)==1
    overwrite = repmat(overwrite, [1 4]);
end

sigma = ip.Results.Sigma;
if isempty(sigma)
    fprintf('Enter the Gaussian s.d. for detection. If the z-value is different, \n');
    sigma = input('enter the two values as, e.g., "[1.5 1.3]": ');
end

LLCopts = {'LLFFCorrection',pr.LLFFCorrection, 'LowerLimit', pr.LowerLimit, 'LSImageCh1', pr.LSImageCh1, 'LSImageCh2',pr.LSImageCh2, 'LSImageCh3',pr.LSImageCh3,'LSImageCh4',pr.LSImageCh4,...
    'BackgroundCh1',pr.BackgroundCh1,'BackgroundCh2',pr.BackgroundCh2,'BackgroundCh3',pr.BackgroundCh3,'BackgroundCh4',pr.BackgroundCh4};
% settings -> input options instead (check all, vs single)
% improve help for this function

%-------------------------------------------------------------------------------
% 1) Parse data, determine output structure
%-------------------------------------------------------------------------------
if isempty(data)
    data = XR_loadConditionData3D();
end

if isempty(apath) % store results locally in cell/Analysis directory
    apath = arrayfun(@(i) [i.source], data, 'unif', 0);
else
    % expand results path for all data sets
    apath = arrayfun(@(i) [apath getShortPath(i,3) filesep], data, 'unif', 0);
end

if isempty(aname)
    apath = cellfun(@(x) [x, 'Analysis', filesep], apath, 'unif', 0);
else
    apath = cellfun(@(x) [x, aname, filesep], apath, 'unif', 0);
end 

[~,~] = cellfun(@mkdir, apath, 'unif', 0);

%-------------------------------------------------------------------------------
% 2) Deskew
%-------------------------------------------------------------------------------
% xruan 11/21/2019
% put z-offset correction as an option for deskew, and correct on the raw
% data before deskew
data = deskewData(data, 'Overwrite', overwrite(1), 'SkewAngle', ip.Results.SkewAngle, 'aname', aname, ...
    'Rotate', ip.Results.Rotate,'sCMOSCameraFlip', ip.Results.sCMOSCameraFlip, 'Crop', ip.Results.Crop, LLCopts{:}, ...
    'zoffsetCorrection', false);

% % xruan 11/19/2019 add option for zoffsetCorrection
% % z-offset correction
if ip.Results.XZoffsetCorrection
    % XR_correctZoffsetData3D(data, 'Overwrite', overwrite(1));
    XR_correctXZoffsetData3D(data, 'Overwrite', overwrite(1));
end

% return;

%-------------------------------------------------------------------------------
% 3) Detection
%-------------------------------------------------------------------------------
switch ParallelMethod
    case 'parfor'
        tic
        XR_runDetection3D(data, 'Sigma', sigma, 'Overwrite', overwrite(2),...
            'WindowSize', pr.WindowSize, 'ResultsPath', apath,'Mode', pr.Mode,...
            'FitMixtures',pr.FitMixtures, 'MaxMixtures',pr.MaxMixtures, ...
            'DetectionMethod', pr.DetectionMethod, 'BackgroundCorrection',pr.BackgroundCorrection);
        toc        
    case 'lock_distributed'
        tic
        XR_runDetection3D_frame_parallel(data, 'Sigma', sigma, 'Overwrite', overwrite(2),...
            'WindowSize', pr.WindowSize, 'ResultsPath', apath,'Mode', pr.Mode,...
            'FitMixtures',pr.FitMixtures, 'MaxMixtures',pr.MaxMixtures, ...
            'DetectionMethod', pr.DetectionMethod, 'BackgroundCorrection',pr.BackgroundCorrection);
        toc
end

%-------------------------------------------------------------------------------
% 4) Tracking
%-------------------------------------------------------------------------------
if ip.Results.Track
    data = data(~([data.movieLength]==1));
    XR_runTracking3D(data, loadTrackSettings('Radius', ip.Results.TrackingRadius,...
        'MaxGapLength', ip.Results.TrackingGapLength),...
        'FileName', 'trackedFeatures.mat', 'Overwrite', overwrite(3), ...
        'ResultsPath', apath, 'ParallelMethod', ParallelMethod);
    
%     GU_runTrackProcessing3D(data, 'Overwrite', overwrite(4),...
%         'TrackerOutput', 'trackedFeatures.mat', 'FileName', 'ProcessedTracks.mat',...
%         'Buffer', ip.Results.Buffer, 'BufferAll', ip.Results.BufferAll,...
%         'FitMixtures',pr.FitMixtures, 'WindowSize', ip.Results.WindowSize, 'ResultsPath', apath);

    XR_runTrackProcessing3D(data, 'Overwrite', overwrite(4),...
        'TrackerOutput', 'trackedFeatures.mat', 'FileName', 'ProcessedTracks.mat',...
        'Buffer', ip.Results.Buffer, 'BufferAll', ip.Results.BufferAll,...
        'FitMixtures',pr.FitMixtures, 'WindowSize', ip.Results.WindowSize, ...
        'ResultsPath', apath, 'ParallelMethod', ParallelMethod);
    
    %-------------------------------------------------------------------------------
    % 5) Rotate tracks
    %-------------------------------------------------------------------------------
    if ip.Results.RotateTracks
        rotateTracks3D(data, 'Overwrite', overwrite(4), 'Crop', ip.Results.Crop, 'ResultsPath', apath);
    end
end