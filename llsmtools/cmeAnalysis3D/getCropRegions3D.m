%[reg, hasDS] = getCropRegions3D(data, varargin) enables cropping movies prior to de-skewing in order to reduce file sizes.
% The function loops through the movies in the input structure and displays projections
% of the first and last frame to facilitate cropping.
%
% Inputs:
%    data : structure returned by loadConditionData()
%
% Outputs:
%     roi : cell array of selected regions of interest for cropping
%   hasDS : flag indicating whether the inputs have already been de-skewed

% Author: Francois Aguet

function [roi, hasDS] = getCropRegions3D(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct); % 
ip.addParamValue('Overwrite', false, @islogical);
ip.parse(data, varargin{:});

nd = numel(data);
roi = cell(1,nd);

hasDS = false(1,nd);
for i = 1:nd
    % check whether de-skewing has been run
    tmp = dir([data(i).source 'DS']);
    tmp = {tmp(~[tmp.isdir]).name}';
    tmp = tmp(~cellfun(@isempty, regexpi(tmp, '\.tif|\.mat')));
    hasDS(i) = (numel(tmp)==data(i).movieLength);
    
    if ~hasDS(i) || ip.Results.Overwrite
        
        % load first frame
        frame1 = double(readtiff(data(i).framePaths{1}{1}));
        frame2 = double(readtiff(data(i).framePaths{1}{end}));
        proj = log(max(cat(3,frame1,frame2),[],3));
        
        nCh = numel(data(i).channels);
        if nCh>1
            [ny,nx] = size(proj);
            rframe = zeros(ny,nx,3,'uint8');
            idxRGB = getRGBindex(data(i).markers);
            rframe(:,:,idxRGB(1)) = uint8(scaleContrast(proj));
                for ci = 2:nCh
                    frame1 = double(readtiff(data(i).framePaths{ci}{1}));
                    frame2 = double(readtiff(data(i).framePaths{ci}{end}));
                    rframe(:,:,idxRGB(ci)) = uint8(scaleContrast(log(max(cat(3,frame1,frame2),[],3))));
                end
            proj = rframe;
        end         
        
        % crop x,y (opens UI)
        hf = figure('Name', getShortPath(data(i),2), 'Position', [100 100 500 500]);
        imshow(proj,[]); axis image;
        nx = data(i).imagesize(2);
        ny = data(i).imagesize(1);
        b = max(nx,ny)/10;
        hr = imrect(gca, [b b nx-2*b ny-2*b]);
        fcn = makeConstrainToRectFcn('imrect', get(gca,'XLim'), get(gca,'YLim'));
        setPositionConstraintFcn(hr, fcn);
        pos = round(wait(hr));
        close(hf);
        drawnow; % forces window to close
        roi{i} = pos;

    end
end
