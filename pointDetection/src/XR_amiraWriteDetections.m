function XR_amiraWriteDetections(filename,frameInfo,data, varargin)
% Write an Amira Mesh file with name [<filename>_%04d.am] representing tracks. 
% Options
%    - <scales>: [x y z] defines relative pixel size (must be synced to
%    amira stack opening) 
%    - <vertexProp>: {{'name',{NVertex x 1, ...}},...} is the vertex-associated 
%    properties  each frame must be described
%    in a cell. 
%    - <edgeProp>: {{'name',NTrack x 1}, ...} is the
%    edge-associated properties
%    Based on Filipe's script, modified to use Francois Aguet's tracking data structure
% usage %
% ch = 1;
% load([data.channels{ch} 'Analysis' filesep 'ProcessedTracks.mat']);
% tracks = tracks(ismember([tracks.catIdx], [1:4]));
% filename = [data.channels{ch} 'AmiraTracks' filesep 'am'];
% %%%%remove nans from amplitudes%%%%
% for i = 1:numel(tracks)
%    tracks(i).A(isnan(tracks(i).A)) = 0;
% end
%  GU_amiraWriteTracks(filename,tracks,data)

% Gokul Upadhyayula, Dec 2015
% 
% copied from GU_amiraWriteTracks.m


ip=inputParser();
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParamValue('scales', [], @isnumeric);
ip.addParamValue('vertexProp',{}, @iscell);
ip.addParamValue('FlipZ',true, @islogical);
ip.addParamValue('edgeProp',{}, @iscell);
ip.addParamValue('AmpMultiplier',1, @isnumeric);
ip.addParamValue('Ch',1, @isnumeric);
ip.parse( varargin{:});
p=ip.Results;

[pathstr,name,ext] = fileparts(filename); 
basename=[pathstr filesep name];

% find the volume dimensions
if p.FlipZ
    im = readtiff(data.framePathsDS{1}{1});
    [~, ~, sz] = size(im);
    if ~isempty(ip.Results.scales);
        s=ip.Results.scales;
    else
        s = [1 1 -data.zAniso];
    end
    zOffset = sz*data.zAniso;
else
    if ~isempty(ip.Results.scales);
        s=ip.Results.scales;
    else
        s = [1 1 data.zAniso];
    end
    zOffset = 0;
end

if(~exist(fileparts(filename)))
    mkdir(fileparts(filename)); 
end;

for fIdx=1:data.movieLength
    %% Indx of tracks on the current frame
%     tracksOn=([tracks.end]>=fIdx)&(fIdx>[tracks.start]);
%     nbTracsOn=sum(tracksOn);
    frameInfoOn = frameInfo(fIdx);
    %% tracks extremity
%     tracksEnds=size(nbTracsOn*2,3);nbTracsOn
    % relative Idx of the end of each tracks (e.g. for use in tracks(i).x)
%     endRelIdx=fIdx-[tracks.start]+1; 
%     count=1;
%     for tIdx=find(tracksOn)
%         tr=tracks(tIdx);
%         % Handling gap
%         tracksEnds(count+1,1)=(tr.x(endRelIdx(tIdx))-1)*s(1);
%         tracksEnds(count+1,2)=(tr.y(endRelIdx(tIdx))-1)*s(2);
%         tracksEnds(count+1,3)=((tr.z(endRelIdx(tIdx))-1))*s(3) + zOffset ;
%         tracksEnds(count,1)=(tr.x(1)-1)*s(1);
%         tracksEnds(count,2)=(tr.y(1)-1)*s(2);
%         tracksEnds(count,3)=(tr.z(1)-1)*s(3) + zOffset;                        
%         count=count+2;
%     end
% tracksEnds(:,3) = tracksEnds(:,3)+120;
    
    %% tracks extremity property (start or end)
    % startEnd=ones(nbTracsOn*2,1);
    % startEnd(1:2:end)=0;

    %% tracks edgeConnectivity
    % tracksEdge=reshape(1:2*nbTracsOn,[2,nbTracsOn])'-1;    
    
    %% tracks point    
    % numPointsPerEdge=max(2,endRelIdx(tracksOn)'); %if one single time point Amira still needs two time points...
    numPoints=size(frameInfoOn.x, 2);
    numPointsPerEdge=ones(numPoints, 1);    
    tracksPoints=zeros(numPoints,3);

    tracksPoints(:, 1) = frameInfoOn.x(1, :)' * s(1);
    tracksPoints(:, 2) = frameInfoOn.y(1, :)' * s(2);
    tracksPoints(:, 3) = frameInfoOn.z(1, :)' * s(3) + zOffset;
    tracksA = frameInfoOn.A(1, :)' * p.AmpMultiplier;   % extract amplitude
    
% tracksPoints(:,3) = tracksPoints(:,3)+120;
    %% Track id (edge property)
%     tracksId=find(tracksOn)';
    
    %% Track lifetime (edge property)
%     tracksLft=[tracks(tracksOn).lifetime_s]';
%     tracksA=[tracks(tracksOn).A(1,:)]'; 
    % Write 
    paramCount=0;
    frameFilename=[basename '_t_' num2str(fIdx,'%04.0f'),'.am'];
    fid = fopen(frameFilename, 'w');
    fprintf(fid,['# AmiraMesh 3D ASCII 2.0\n\n']);
    fprintf(fid,['define VERTEX ' num2str(numPoints) '\n']);
    fprintf(fid,['define EDGE ' num2str(numPoints) ' \n']);    
    fprintf(fid,['define POINT ' num2str(numPoints) '\n']);
    fprintf(fid,'\n');
    fprintf(fid,'Parameters { ContentType "HxSpatialGraph"}\n\n');
    fprintf(fid,'VERTEX { float[3] VertexCoordinates } @1\n');
    fprintf(fid,'EDGE { int[2] EdgeConnectivity } @2\n');
    fprintf(fid,'EDGE { int NumEdgePoints } @3\n');
    fprintf(fid,'POINT { float[3] EdgePointCoordinates } @4\n');
    fprintf(fid,'POINT { int Amplitude } @5\n');
%     fprintf(fid,'VERTEX {int Amplitude } @9\n');
    fclose(fid);
                
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@1\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksPoints, '-append', 'delimiter',' ','precision', 16);
                
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@2\n');
        fclose(fid);
        
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@3\n');
        fclose(fid);
        dlmwrite(frameFilename, numPointsPerEdge, '-append', 'delimiter',' ','precision', 16);
        
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@4\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksPoints, '-append', 'delimiter',' ','precision', 16);
        
        paramCount=paramCount+1;
        fid = fopen(frameFilename, 'a');
        fprintf(fid,'\n@5\n');
        fclose(fid);
        dlmwrite(frameFilename, tracksA, '-append', 'delimiter',' ','precision', 16);
        
%         paramCount=paramCount+1;
%         fid = fopen(frameFilename, 'a');
%         fprintf(fid,'\n@9\n');
%         fclose(fid);
%         dlmwrite(frameFilename, tracksA, '-append', 'delimiter',' ','precision', 16);
        
%         for propIdx=1:length(p.vertexProp)
%             fid = fopen(frameFilename, 'a');
%             paramCount=paramCount+1;
%             fprintf(fid,['\nVERTEX { float ' p.vertexProp{propIdx}{1} ' } @' num2str(paramCount+propIdx-1) '\n']);
%             fprintf(fid,['@' num2str(paramCount+propIdx-1) '\n']);
%             fclose(fid);
%             dlmwrite(frameFilename, p.vertexProp{propIdx}{2}{fIdx}, '-append', 'delimiter',' ','precision', 16)           
%         end
%         
%         for propIdx=1:length(p.edgeProp)
%             fid = fopen(frameFilename, 'a');
%             paramCount=paramCount+1;
%             fprintf(fid,['\nEDGE { float ' p.edgeProp{propIdx}{1} ' } @' num2str(paramCount+propIdx-1) '\n']);
%             fprintf(fid,['@' num2str(paramCount+propIdx-1) '\n']);
%             fclose(fid);
%             dlmwrite(frameFilename, p.edgeProp{propIdx}{2}(tracksOn), '-append', 'delimiter',' ','precision', 16)
%         end
    end
end
    
    
