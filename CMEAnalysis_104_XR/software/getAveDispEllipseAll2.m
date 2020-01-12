function [longVecS,longVecE,shortVecS,shortVecE,shortVecS3D,shortVecE3D,...
    longVecSMS,longVecEMS,shortVecSMS,shortVecEMS,shortVecS3DMS,shortVecE3DMS] = ...
    getAveDispEllipseAll2(xyzVel,brownStd,trackType,undetBrownStd,timeWindow,...
    brownStdMult,linStdMult,timeReachConfB,timeReachConfL,minSearchRadius,...
    maxSearchRadius,useLocalDensity,closestDistScale,maxStdMult,...
    nnDistLinkedFeat,nnWindow,trackStartTime,trackEndTime,probDim,resLimit,...
    brownScaling,linScaling)
%GETAVEDISPELLIPSEALL2 determines the search ellipse and expected displacement along x and y of a particle undergoing 2D diffusion with drift
%
%SYNOPSIS [longVecS,longVecE,shortVecS,shortVecE,shortVecS3D,shortVecE3D,...
%    longVecSMS,longVecEMS,shortVecSMS,shortVecEMS,shortVecS3DMS,shortVecE3DMS] = ...
%    getAveDispEllipseAll2(xyzVel,brownStd,trackType,undetBrownStd,timeWindow,...
%    brownStdMult,linStdMult,timeReachConfB,timeReachConfL,minSearchRadius,...
%    maxSearchRadius,useLocalDensity,closestDistScale,maxStdMult,...
%    nnDistLinkedFeat,nnWindow,trackStartTime,trackEndTime,probDim,resLimit)
%    brownScaling,linScaling)
%
%INPUT  xyzVel         : Velocity in x, y and z (if 3D).
%       brownStd       : Standard deviation of Brownian motion steps.
%       trackType      : Type of track. 1 for directed, 0 for Brownian, NaN for undetermined.
%       undetBrownStd  : Standard deviation of Brownian motion steps to be used
%                        for undetermined tracks.
%       timeWindow     : Maximum gap size.
%       brownStdMult   : Multiplication factor to go from average Brownian
%                        displacement to search radius.
%       linStdMult     : Multiplication factor to go from average linear
%                        displacement to search radius.
%       timeReachConfB : Time gap for Brownian motion to reach confinement.
%       timeReachConfL : Time gap for linear motion to reach confinement.
%       minSearchRadius: Minimum allowed search radius.
%       maxSearchRadius: Maximum allowed search radius for linking between
%                        two consecutive frames. It will be expanded for
%                        different gap lengths based on the time scaling of
%                        Brownian motion.
%       useLocalDensity: 1 if local density of features is used to expand 
%                        their search radius if possible, 0 otherwise.
%       closestDistScale:Scaling factor of nearest neighbor distance.
%       maxStdMult     : Maximum value of factor multiplying std to get
%                        search radius.
%       nnDistLinkedFeat:Matrix indicating the nearest neighbor
%                        distances of features linked together within
%                        tracks.
%       nnWindow       : Time window to be used in estimating the
%                        nearest neighbor distance of a track at its start
%                        and end.
%       trackStartTime : Starting time of all tracks.
%       trackEndTime   : Ending time of all tracks.
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%       resLimit       : Resolution limit, in whatever space units are
%                        being used.
%                        Optional. Default: 0 (i.e. don't use).
%       brownScaling   : Power by which Brownian part of search radius
%                        inceases with time.
%                        Optional. Default: 0.5.
%       linScaling     : Power by which linear part of search radius
%                        increases with time.
%                        Optional. Default: 0.5.
%
%OUTPUT longVecS  : Vector defining long radius of search ellipse/ellipsoid at the
%                   starts of tracks.
%       longVecE  : Vector defining long radius of search ellipse/ellipsoid at the
%                   ends of tracks.
%       shortVecS : Vector defining short radius of search ellipse/ellipsoid at the
%                   starts of tracks.
%       shortVecE : Vector defining short radius of search ellipse/ellipsoid at the
%                   ends of tracks.
%       shortVecS3D:Vector defining 2nd short radius of search ellipse/ellipsoid at the
%                   starts of tracks in case of 3D.
%       shortVecE3D:Vector defining 2nd short radius of search ellipse/ellipsoid at the
%                   ends of tracks in case of 3D.
%       longVecSMS, longVecEMS, shortVecSMS, shortVecEMS, shortVecS3DMS,
%       shortVecE3DMS: Same as above, but for merging ang splitting.
%
%REMARKS Drift is assumed to look more like 1D diffusion, i.e. the particle
%goes back and forth along a line
%
%Khuloud Jaqaman, April 2007
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dispDrift = [];
% dispBrown = [];
longVecS  = [];
longVecE  = [];
shortVecS = [];
shortVecE = [];
shortVecS3D = [];
shortVecE3D = [];
longVecSMS  = [];
longVecEMS  = [];
shortVecSMS = [];
shortVecEMS = [];
shortVecS3DMS = [];
shortVecE3DMS = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 19
    disp('--getAveDispEllipseAll2: Incorrect number of input arguments!');
    return
end

if nargin < 20 || isempty(resLimit)
    resLimit = 0;
end

if nargin < 21 || isempty(brownScaling)
    brownScaling = [0.5 0.01];
end

if nargin < 22 || isempty(linScaling)
    linScaling = [0.5 0.01];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine expected displacement and search ellipse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine number of tracks
numTracks = size(xyzVel,1);

%reserve memory for output
% dispDrift = zeros(probDim,timeWindow,numTracks);
% dispBrown = zeros(timeWindow,numTracks);
longVecS  = zeros(probDim,timeWindow,numTracks);
longVecE  = zeros(probDim,timeWindow,numTracks);
shortVecS = zeros(probDim,timeWindow,numTracks);
shortVecE = zeros(probDim,timeWindow,numTracks);
longVecSMS  = zeros(probDim,timeWindow,numTracks);
longVecEMS  = zeros(probDim,timeWindow,numTracks);
shortVecSMS = zeros(probDim,timeWindow,numTracks);
shortVecEMS = zeros(probDim,timeWindow,numTracks);
if probDim == 3
    shortVecS3D = zeros(probDim,timeWindow,numTracks);
    shortVecE3D = zeros(probDim,timeWindow,numTracks);
    shortVecS3DMS = zeros(probDim,timeWindow,numTracks);
    shortVecE3DMS = zeros(probDim,timeWindow,numTracks);
end

%define square root of "problem dimension" to avoid calculating it many times
sqrtDim = sqrt(probDim);

%put time scaling of linear motion in a vector
timeScalingLin = [(1:timeReachConfL).^linScaling(1) ...
    (timeReachConfL)^linScaling(1) * (2:timeWindow-timeReachConfL+1).^linScaling(2)];

%put time scaling of Brownian motion in a vector
timeScalingBrown = [(1:timeReachConfB).^brownScaling(1) ...
    (timeReachConfB)^brownScaling(1) * (2:timeWindow-timeReachConfB+1).^brownScaling(2)];

%scale maxSearchRadius like Brownian motion (it's only imposed on the
%Brownian aspect of tracks)
maxSearchRadius = maxSearchRadius * timeScalingBrown;

%calculate minimum and maximum search radii for merging and splitting,
%taking into account the point spread function width
minSearchRadiusMS = max(minSearchRadius,resLimit);
maxSearchRadiusMS = max(maxSearchRadius,resLimit);

%determine the nearest neighbor distances of tracks at their starts and ends
windowLimS = min([trackStartTime+nnWindow trackEndTime],[],2);
windowLimE = max([trackEndTime-nnWindow trackStartTime],[],2);
nnDistTracksS = zeros(numTracks,1);
nnDistTracksE = zeros(numTracks,1);
for iTrack = 1 : numTracks
    nnDistTracksS(iTrack) = min(nnDistLinkedFeat(iTrack,...
        trackStartTime(iTrack):windowLimS(iTrack)));
    nnDistTracksE(iTrack) = min(nnDistLinkedFeat(iTrack,...
        windowLimE(iTrack):trackEndTime(iTrack)));
end

for iTrack = 1 : numTracks
    
    switch trackType(iTrack)

        case 1
            
            %get velocity, its magnitude and direction of motion
            velDrift = xyzVel(iTrack,:)';
            velMag = sqrt(velDrift' * velDrift);
            directionMotion = velDrift / velMag;
            
            %obtain vector(s) perpendicular to direction of motion
            if probDim == 2 %in 2D case, 1 vector needed
                perpendicular = [-directionMotion(2) directionMotion(1)]';
            else %in 3D case, 2 vectors needed
                perpendicular = [-directionMotion(2) directionMotion(1) 0]';
                perpendicular = perpendicular / (sqrt(perpendicular'*perpendicular));
                perpendicular3D = cross(directionMotion,perpendicular);
            end

            %calculate the expected displacement due to drift for all time
            %gaps
            dispDrift1 = repmat(velDrift,1,timeWindow) .* repmat(timeScalingLin,probDim,1);

            %calculate the expected displacement along x (= along y, [z]) due to
            %brownian motion for all time gaps
            dispBrown1 = brownStd(iTrack) * timeScalingBrown;
            
            %copy brownStdMult into vector that might be modified using
            %local density
            brownStdMultModS = brownStdMult'; %for track start
            brownStdMultModE = brownStdMult'; %for track end

            %if local density information is used to expand search radius ...
            if useLocalDensity

                %divide the track's nearest neighbor distance at its start
                %/closestDistScale by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksS(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand search radius multiplication factor at track start
                %if possible
                brownStdMultModS = max([brownStdMultModS; ratioDist2Std]);

                %divide the track's nearest neighbor distance at its end
                %/closestDistScale by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksE(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand search radius multiplication factor at track end
                %if possible
                brownStdMultModE = max([brownStdMultModE; ratioDist2Std]);

            end

            %determine the "long vectors" of the search rectangles for all time
            %gaps
            longVec1 = repmat(linStdMult',probDim,1) .* dispDrift1 + ...
                repmat((brownStdMult' .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);
            longVecMag = sqrt((diag(longVec1' * longVec1))');  %magnitude
            longVecDir = longVec1 ./ repmat(longVecMag,probDim,1); %direction
            
            %determine the "short vectors" at track starts
            shortVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);
            shortVecSMag = sqrt((diag(shortVecS1' * shortVecS1))');  %magnitude
            shortVecSDir = shortVecS1 ./ repmat(shortVecSMag,probDim,1); %direction

            %determine the "short vectors" at track ends
            shortVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);
            shortVecEMag = sqrt((diag(shortVecE1' * shortVecE1))');  %magnitude
            shortVecEDir = shortVecE1 ./ repmat(shortVecEMag,probDim,1); %direction
            
            %             %output the absolute value of dispDrift1
            %             dispDrift1 = abs(dispDrift1);

            %make sure that "long vectors" are longer than minimum allowed
            longVecMagTmp = max([longVecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            longVec1 = repmat(longVecMagTmp,probDim,1) .* longVecDir; %new long vector
            %do the same for merging and splitting
            longVecMagTmp = max([longVecMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            longVec1MS = repmat(longVecMagTmp,probDim,1) .* longVecDir;

            %make sure that "short vectors" at track starts are within
            %allowed range
            shortVecSMagTmp = max([shortVecSMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            shortVecSMagTmp = min([shortVecSMagTmp;maxSearchRadius]); %compare to maximum
            shortVecS1 = repmat(shortVecSMagTmp,probDim,1) .* shortVecSDir; %new short vector
            %do the same for merging and splitting
            shortVecSMagTmpMS = max([shortVecSMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            shortVecSMagTmpMS = min([shortVecSMagTmpMS;maxSearchRadiusMS]);
            shortVecS1MS = repmat(shortVecSMagTmpMS,probDim,1) .* shortVecSDir;

            %make sure that "short vectors" at track ends are within allowed
            %range
            shortVecEMagTmp = max([shortVecEMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            shortVecEMagTmp = min([shortVecEMagTmp;maxSearchRadius]); %compare to maximum
            shortVecE1 = repmat(shortVecEMagTmp,probDim,1) .* shortVecEDir; %new short vector
            %do the same for merging and splitting
            shortVecEMagTmpMS = max([shortVecEMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            shortVecEMagTmpMS = min([shortVecEMagTmpMS;maxSearchRadiusMS]);
            shortVecE1MS = repmat(shortVecEMagTmpMS,probDim,1) .* shortVecEDir;
            
            %save values for this track
            %             dispDrift(:,:,iTrack) = dispDrift1;
            %             dispBrown(:,iTrack) = dispBrown1;
            longVecS(:,:,iTrack) = longVec1;
            longVecE(:,:,iTrack) = longVec1;
            shortVecS(:,:,iTrack) = shortVecS1;
            shortVecE(:,:,iTrack) = shortVecE1;
            longVecSMS(:,:,iTrack) = longVec1MS;
            longVecEMS(:,:,iTrack) = longVec1MS;
            shortVecSMS(:,:,iTrack) = shortVecS1MS;
            shortVecEMS(:,:,iTrack) = shortVecE1MS;
            
            %construct additional short vectors for 3D problems
            if probDim == 3
                shortVecS13D = repmat(shortVecSMagTmp,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE13D = repmat(shortVecEMagTmp,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3D(:,:,iTrack) = shortVecS13D;
                shortVecE3D(:,:,iTrack) = shortVecE13D;
                %do the same for merging and splitting
                shortVecS13DMS = repmat(shortVecSMagTmpMS,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE13DMS = repmat(shortVecEMagTmpMS,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3DMS(:,:,iTrack) = shortVecS13DMS;
                shortVecE3DMS(:,:,iTrack) = shortVecE13DMS;
            end

        case 0
            
            %take direction of motion to be along x and construct
            %perpendicular(s)
            if probDim == 2
                directionMotion = [1 0]';
                perpendicular = [0 1]';
            else
                directionMotion = [1 0 0]';
                perpendicular = [0 1 0]';
                perpendicular3D = [0 0 1]';
            end

            %calculate the expected displacement along x (= along y, [z]) due to
            %brownian motion for all time gaps
            dispBrown1 = brownStd(iTrack) * timeScalingBrown;

            %copy brownStdMult into vector that might be modified using
            %local density
            brownStdMultModS = brownStdMult'; %for track start
            brownStdMultModE = brownStdMult'; %for track end
            
            %if local density information is used to expand search radius ...
            if useLocalDensity

                %divide the track's nearest neighbor distance/closestDistScale
                %by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksS(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand start's search radius multiplication factor if possible
                brownStdMultModS = max([brownStdMultModS; ratioDist2Std]);

                %divide the track's nearest neighbor distance at its end
                %/closestDistScale by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksE(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand end's search radius multiplication factor if possible
                brownStdMultModE = max([brownStdMultModE; ratioDist2Std]);

            end

            %determine the long vectors of the search ellipses at track
            %starts for all time gaps
            longVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the long vectors of the search ellipses at track
            %ends for all time gaps
            longVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the short vectors at track starts
            shortVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);

            %determine the short vectors at track ends
            shortVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);

            %get magnitude and direction of both vectors at track starts
            vecMag = sqrt((diag(longVecS1' * longVecS1))'); %magnitude of both vectors
            longVecDir = longVecS1 ./ repmat(vecMag,probDim,1);   %direction of long vector
            shortVecDir = shortVecS1 ./ repmat(vecMag,probDim,1); %direction of short vector
            
            %make sure that magnitude is within allowed range
            vecMagTmp = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            vecMagTmp = min([vecMagTmp;maxSearchRadius]); %compare to maximum
            %repeat for merging and splitting
            vecMagTmpMS = max([vecMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            vecMagTmpMS = min([vecMagTmpMS;maxSearchRadiusMS]);
            
            %re-calculate both vectors based on modified magnitudes            
            longVecS1 = repmat(vecMagTmp,probDim,1) .* longVecDir; %new long vector
            shortVecS1 = repmat(vecMagTmp,probDim,1) .* shortVecDir; %new short vector
            %repeat for merging and splitting
            longVecS1MS = repmat(vecMagTmpMS,probDim,1) .* longVecDir;
            shortVecS1MS = repmat(vecMagTmpMS,probDim,1) .* shortVecDir;

            %construct additional short vectors for 3D problems and save
            %them
            if probDim == 3
                shortVecS13D = repmat(vecMagTmp,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3D(:,:,iTrack) = shortVecS13D;
                %repeat for merging and splitting
                shortVecS13DMS = repmat(vecMagTmpMS,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3DMS(:,:,iTrack) = shortVecS13DMS;
            end

            %get magnitude and direction of both vectors at track ends
            vecMag = sqrt((diag(longVecE1' * longVecE1))');  %magnitude of both vectors
            longVecDir = longVecE1 ./ repmat(vecMag,probDim,1);   %direction of long vector
            shortVecDir = shortVecE1 ./ repmat(vecMag,probDim,1); %direction of short vector
            
            %make sure that magnitude is larger than minimum allowed and
            %smaller than maximum allowed
            vecMagTmp = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            vecMagTmp = min([vecMagTmp;maxSearchRadius]); %compare to maximum
            %repeat for merging and splitting
            vecMagTmpMS = max([vecMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            vecMagTmpMS = min([vecMagTmpMS;maxSearchRadiusMS]);
            
            %re-calculate both vectors based on modified magnitudes            
            longVecE1 = repmat(vecMagTmp,probDim,1) .* longVecDir; %new long vector
            shortVecE1 = repmat(vecMagTmp,probDim,1) .* shortVecDir; %new short vector
            %repeat for merging and splitting
            longVecE1MS = repmat(vecMagTmpMS,probDim,1) .* longVecDir;
            shortVecE1MS = repmat(vecMagTmpMS,probDim,1) .* shortVecDir;

            %construct additional short vectors for 3D problems and save
            %them
            if probDim == 3
                shortVecE13D = repmat(vecMagTmp,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE3D(:,:,iTrack) = shortVecE13D;
                %repear for merging and splitting
                shortVecE13DMS = repmat(vecMagTmpMS,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE3DMS(:,:,iTrack) = shortVecE13DMS;
            end

            %save values for this track
            %             dispBrown(:,iTrack) = dispBrown1;
            longVecS(:,:,iTrack) = longVecS1;
            longVecE(:,:,iTrack) = longVecE1;
            shortVecS(:,:,iTrack) = shortVecS1;
            shortVecE(:,:,iTrack) = shortVecE1;
            longVecSMS(:,:,iTrack) = longVecS1MS;
            longVecEMS(:,:,iTrack) = longVecE1MS;
            shortVecSMS(:,:,iTrack) = shortVecS1MS;
            shortVecEMS(:,:,iTrack) = shortVecE1MS;

        otherwise

            %take direction of motion to be along x and construct
            %perpendicular(s)
            if probDim == 2
                directionMotion = [1 0]';
                perpendicular = [0 1]';
            else
                directionMotion = [1 0 0]';
                perpendicular = [0 1 0]';
                perpendicular3D = [0 0 1]';
            end

            %calculate the expected displacement along x (= along y) due to
            %brownian motion for all time gaps
            %             dispBrown1 = undetBrownStd * timeScalingBrown;
            if brownStd(iTrack)==1
                dispBrown1 = undetBrownStd * timeScalingBrown;
            else
                dispBrown1 = brownStd(iTrack) * timeScalingBrown;
            end

            %copy brownStdMult into vector that might be modified using
            %local density
            brownStdMultModS = brownStdMult'; %for track start
            brownStdMultModE = brownStdMult'; %for track end

            %if local density information is used to expand search radius ...
            if useLocalDensity

                %divide the track's nearest neighbor distance/closestDistScale
                %by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksS(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand start's search radius multiplication factor if possible
                brownStdMultModS = max([brownStdMultModS; ratioDist2Std]);

                %divide the track's nearest neighbor distance/closestDistScale
                %by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksE(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand end's search radius multiplication factor if possible
                brownStdMultModE = max([brownStdMultModE; ratioDist2Std]);

            end
            
            %determine the long vector of the search ellipse at track
            %starts for all time gaps
            longVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the long vector of the search ellipse at track
            %ends for all time gaps
            longVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the short vector at track starts
            shortVecS1 = repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);

            %determine the short vector at track ends
            shortVecE1 = repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);

            %get magnitude and direction of both vectors at track starts
            vecMag = sqrt((diag(longVecS1' * longVecS1))'); %magnitude of both vectors
            longVecDir = longVecS1 ./ repmat(vecMag,probDim,1);   %direction of long vector
            shortVecDir = shortVecS1 ./ repmat(vecMag,probDim,1); %direction of short vector

            %make sure that magnitude is larger than minimum allowed and
            %smaller than maximum allowed
            vecMagTmp = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            vecMagTmp = min([vecMagTmp;maxSearchRadius]); %compare to maximum
            %repeat for merging and spltting
            vecMagTmpMS = max([vecMag;repmat(minSearchRadiusMS,1,timeWindow)]); %compare to minimum
            vecMagTmpMS = min([vecMagTmpMS;maxSearchRadiusMS]); %compare to maximum

            %re-calculate both vectors based on modified magnitudes
            longVecS1 = repmat(vecMagTmp,probDim,1) .* longVecDir; %new long vector
            shortVecS1 = repmat(vecMagTmp,probDim,1) .* shortVecDir; %new short vector
            %repeat for merging and splitting
            longVecS1MS = repmat(vecMagTmpMS,probDim,1) .* longVecDir;
            shortVecS1MS = repmat(vecMagTmpMS,probDim,1) .* shortVecDir;

            %construct additional short vectors for 3D problems
            if probDim == 3
                shortVecS13D = repmat(vecMagTmp,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3D(:,:,iTrack) = shortVecS13D;
                %repeat for merging and splitting
                shortVecS13DMS = repmat(vecMagTmpMS,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3DMS(:,:,iTrack) = shortVecS13DMS;
            end

            %get magnitude and direction of both vectors at track ends
            vecMag = sqrt((diag(longVecE1' * longVecE1))'); %magnitude of both vectors
            longVecDir = longVecE1 ./ repmat(vecMag,probDim,1);   %direction of long vector
            shortVecDir = shortVecE1 ./ repmat(vecMag,probDim,1); %direction of short vector
            
            %make sure that magnitude is larger than minimum allowed and
            %smaller than maximum allowed
            vecMagTmp = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            vecMagTmp = min([vecMagTmp;maxSearchRadius]); %compare to maximum
            %repeat for merging and spltting
            vecMagTmpMS = max([vecMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            vecMagTmpMS = min([vecMagTmpMS;maxSearchRadiusMS]);
            
            %re-calculate both vectors based on modified magnitudes
            longVecE1 = repmat(vecMagTmp,probDim,1) .* longVecDir; %new long vector
            shortVecE1 = repmat(vecMagTmp,probDim,1) .* shortVecDir; %new short vector
            %repeat for merging and splitting
            longVecE1MS = repmat(vecMagTmpMS,probDim,1) .* longVecDir; %new long vector
            shortVecE1MS = repmat(vecMagTmpMS,probDim,1) .* shortVecDir; %new short vector

            %construct additional short vectors for 3D problems
            if probDim == 3
                shortVecE13D = repmat(vecMagTmp,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE3D(:,:,iTrack) = shortVecE13D;
                %repeat for merging and splitting
                shortVecE13DMS = repmat(vecMagTmpMS,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE3DMS(:,:,iTrack) = shortVecE13DMS;
            end
            
            %save values for this track
            %             dispBrown(:,iTrack) = dispBrown1;
            longVecS(:,:,iTrack) = longVecS1;
            longVecE(:,:,iTrack) = longVecE1;
            shortVecS(:,:,iTrack) = shortVecS1;
            shortVecE(:,:,iTrack) = shortVecE1;
            longVecSMS(:,:,iTrack) = longVecS1MS;
            longVecEMS(:,:,iTrack) = longVecE1MS;
            shortVecSMS(:,:,iTrack) = shortVecS1MS;
            shortVecEMS(:,:,iTrack) = shortVecE1MS;

    end %(switch trackType)

end %(for iTrack = 1 : numTracks)

%%%%% ~~ the end ~~ %%%%%

