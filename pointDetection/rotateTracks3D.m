%rotateTracks3D(data, varargin) rotates processed tracks to the same frame
% of reference as rotated stacks, with the coverslip horizontal
%
% See also rotateFrame3D
% 
% Inputs
%              data : list of movies, using the structure returned by loadConditionData.m
%
% Options
%       'Overwrite' : true|{false}. Overwrite previous processing result.
% 'ProcessedTracks' : filename for processed tracks.
%   'RotatedTracks' : filename for rotated tracks. 
%     'ResultsPath' : Result directory to store the results. 
%            'Crop' : crop a ROI for rotation. Default: true.  

% Francois Aguet, 2014
% updated angle sign - Gokul U
% xruan, 2019, update result directory and check whether processed track
% file exists.

function rotateTracks3D(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data');
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('ProcessedTracks', 'ProcessedTracks.mat', @ischar);
ip.addParameter('RotatedTracks', 'RotatedTracks.mat', @ischar);
ip.addParameter('ResultsPath', arrayfun(@(i) [i.source 'Analysis' filesep], data, 'unif', 0), @iscell);
ip.addParameter('Crop', true, @islogical);
ip.parse(data, varargin{:});

nd = numel(data);
for i = 1:nd
    % xruan first check if the processed track result exists. 
    tPath = [ip.Results.ResultsPath{i} filesep ip.Results.ProcessedTracks];
    if ~exist(tPath, 'file')
        fprintf('rotateTracks3D: no processed tracking data found for %s\n', getShortPath(data));
        continue;
    end
    
    if ~(exist([ip.Results.ResultsPath{i} filesep ip.Results.RotatedTracks], 'file')==2) || ip.Results.Overwrite
        
        vol = double(readtiff(data(i).framePathsDS{1}{1}));
        [ny,nx,nz] = size(vol);
        
        theta = -data(i).angle*pi/180;
        
        doCrop = false;
        if ~data(i).objectiveScan
            % calculate height of rotated volume
            % project and find left/right edges of deskewed rhomboid, take median width
            proj = squeeze(max(vol,[],1))'~=0;
            proj = diff(proj,1,2);
            
            startIndex = ones(nz,1);
            endIndex = nx*ones(nz,1);
            [srow,scol] = find(proj==1);
            [erow,ecol] = find(proj==-1);
            startIndex(srow) = scol;
            endIndex(erow) = ecol;
            
            doCrop = any(startIndex>1 & endIndex<nx);
            
            % calculate height; first & last 2 frames have interpolation artifacts
            w = median(diff([startIndex endIndex],[],2));
            h = round(abs(w*tan(theta)*cos(theta)))-4;
        end
        
        center = ([ny nx nz]+1)/2;
        T1 = [1 0 0 0
            0 1 0 0
            0 0 1 0
            -center([2 1 3]) 1];
        
        S = [1 0 0 0
            0 1 0 0
            0 0 data(i).zAniso 0
            0 0 0 1];
        
        % Rotate x,z
        R = [cos(theta) 0 -sin(theta) 0; % order for imwarp is x,y,z
            0 1 0 0;
            sin(theta) 0 cos(theta) 0;
            0 0 0 1];
        
        if ip.Results.Crop && doCrop
            outSize = round([ny nx/cos(theta) h]);
        else
            % exact proportions of rotated box
            outSize = round([ny nx*cos(theta)+nz*data(i).zAniso*sin(abs(theta)) nz*data(i).zAniso*cos(theta)+nx*sin(abs(theta))]);
        end
        
        T2 = [1 0 0 0
            0 1 0 0
            0 0 1 0
            (outSize([2 1 3])+1)/2 1];
        
        tform = affine3d(T1*S*R*T2);
        
        % tracks = load([data(i).source 'Analysis', filesep ip.Results.ProcessedTracks]);
        tracks = load([ip.Results.ResultsPath{i}, filesep ip.Results.ProcessedTracks]);
        tracks = tracks.tracks;
        
        % xruan check if there is tracks
        if isempty(fieldnames(tracks))
            fprintf('rotateTracks3D: no track is found for %s\n', getShortPath(data));
            continue;
        end
            
        nc = numel(data(i).channels);
        for t = 1:numel(tracks)
            for c = 1:nc
                [tracks(t).x(c,:),tracks(t).y(c,:),tracks(t).z(c,:)] ...
                    = transformPointsForward(tform, tracks(t).x(c,:), tracks(t).y(c,:), tracks(t).z(c,:));
                
                if ~isempty(tracks(t).startBuffer)
                    [tracks(t).startBuffer.x(c,:), tracks(t).startBuffer.y(c,:), tracks(t).startBuffer.z(c,:)] ...
                        = transformPointsForward(tform, tracks(t).startBuffer.x(c,:),...
                        tracks(t).startBuffer.y(c,:), tracks(t).startBuffer.z(c,:));
                end
                
                if ~isempty(tracks(t).endBuffer)
                    [tracks(t).endBuffer.x(c,:), tracks(t).endBuffer.Y(c,:), tracks(t).endBuffer.z(c,:)] ...
                        = transformPointsForward(tform, tracks(t).endBuffer.x(c,:),...
                        tracks(t).endBuffer.y(c,:), tracks(t).endBuffer.z(c,:));
                end
            end
        end
        % save([data(i).source 'Analysis' filesep ip.Results.RotatedTracks], 'tracks');
        save('-v7.3', [ip.Results.ResultsPath{i} filesep ip.Results.RotatedTracks], 'tracks');
    end
end
