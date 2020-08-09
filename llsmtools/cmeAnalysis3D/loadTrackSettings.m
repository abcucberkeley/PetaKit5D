% Adapted from 'ScriptTrackGeneral' in 'trackWithGapClosing/Kalman'

% Francois Aguet, November 2010

function trackSettings = loadTrackSettings(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Radius', []);
ip.addParamValue('GapRadius', [5 10]);
ip.addParamValue('LinkRadius', [5 10]);
ip.addParamValue('MaxGapLength', 2);
ip.parse(varargin{:});

gapRadius = ip.Results.GapRadius;
linkRadius = ip.Results.LinkRadius;
if ~isempty(ip.Results.Radius)
    gapRadius = ip.Results.Radius;
    linkRadius = ip.Results.Radius;
end

gapCloseParam.timeWindow = ip.Results.MaxGapLength+1;  % maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
gapCloseParam.mergeSplit = 1;  % 1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = 1; % minimum length of track segments from linking to be used in gap closing.
gapCloseParam.diagnostics = 0; % 1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.


% cost matrix for frame-to-frame linking
costMatrices(1).funcName = 'costMatLinearMotionLink2';
costMatrices(1).parameters.linearMotion = 0;     % use linear motion Kalman filter.
costMatrices(1).parameters.minSearchRadius = linkRadius(1); % minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
costMatrices(1).parameters.maxSearchRadius = linkRadius(2); % maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
costMatrices(1).parameters.brownStdMult = 3;     % multiplication factor to calculate search radius from standard deviation.

costMatrices(1).parameters.useLocalDensity = 1;  % 1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
costMatrices(1).parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
costMatrices(1).parameters.kalmanInitParam = []; % Kalman filter initialization parameters.
costMatrices(1).parameters.diagnostics = [];     % plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.


% cost matrix for gap closing
costMatrices(2).funcName = 'costMatLinearMotionCloseGaps2';
costMatrices(2).parameters.linearMotion = 0; % 1: linear motion Kalman filter.
costMatrices(2).parameters.minSearchRadius = gapRadius(1); % minimum allowed search radius.
costMatrices(2).parameters.maxSearchRadius = gapRadius(2); % maximum allowed search radius.
costMatrices(2).parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
costMatrices(2).parameters.brownScaling = [0.5 0.01];   % power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
costMatrices(2).parameters.timeReachConfB = gapCloseParam.timeWindow; %in the code, the search radius expands with the time gap (since a particle is expected to move further away in a longer gap than in a shorter one). This parameter controls how fast the search radius grows with time. timeReachConfB stands for time to reach confinement for the Brownian part of the motion. So before timeReachConfB, the search radius grows with the square root of time, after that it grows very, very slowly (it's almost fixed).

costMatrices(2).parameters.ampRatioLimit = [0 Inf]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.
costMatrices(2).parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
costMatrices(2).parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
costMatrices(2).parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
costMatrices(2).parameters.linStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.

costMatrices(2).parameters.linScaling = [1 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).

costMatrices(2).parameters.timeReachConfL = gapCloseParam.timeWindow; %same as timeReachConfB, but for the linear part of the motion.
costMatrices(2).parameters.maxAngleVV = 45; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.

%optional; if not input, 1 will be used (i.e. no penalty)
costMatrices(2).parameters.gapPenalty = []; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).
%optional; to calculate MS search radius
%if not input, MS search radius will be the same as gap closing search radius
costMatrices(2).parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.


% Kalman filter function names
kalmanFunctions.reserveMem  = 'kalmanResMemLM';
kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';


trackSettings.costMatrices = costMatrices;
trackSettings.gapCloseParam = gapCloseParam;
trackSettings.kalmanFunctions = kalmanFunctions;