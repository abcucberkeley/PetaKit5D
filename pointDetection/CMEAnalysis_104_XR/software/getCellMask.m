%[mask] = getCellMask(data, varargin)
%
% Notes: for the 'bgproj' method, the detection must be run before
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

% Francois Aguet, 02/17/2012

function [mask, proj] = getCellMask(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('Channel', 1, @isposint);
ip.addParameter('Connect', true, @islogical);
ip.addParameter('Display', false, @islogical);
ip.addParameter('ShowHistogram', false, @islogical);
ip.addParameter('ModeRatio', 0.8, @isscalar);
ip.addParameter('Mode', 'maxproj', @(x) any(strcmpi(x, {'maxproj', 'bgproj', 'density'})));
ip.addParameter('Validate', false, @islogical);
ip.parse(data, varargin{:});

nd = numel(data);

mask = cell(1,nd);
proj = cell(1,nd);

ch = ip.Results.Channel;

se = strel('disk', 5);

hasStoredMask = true(1,nd);
for i = 1:nd
    maskPath = [data(i).source 'Detection' filesep 'cellmask.tif'];
    projPath = [data(i).source 'Detection' filesep 'maxproj.tif'];
    
    if ~(exist(maskPath, 'file')==2) || ~(exist(projPath, 'file')==2) || ip.Results.Overwrite
        hasStoredMask(i) = false;
        
        switch ip.Results.Mode
            case 'maxproj'                
                if ~iscell(data(i).framePaths{ch})
                    stack = readtiff(data(i).framePaths{ch});
                else
                    stack = zeros([data(i).imagesize data(i).movieLength], 'uint16');
                    for f = 1:data(i).movieLength
                        stack(:,:,f) = imread(data(i).framePaths{ch}{f});
                    end
                end
                
                proj{i} = max(stack, [], 3);
                mask{i} = maskFromFirstMode(double(proj{i}), 'Connect', ip.Results.Connect,...
                    'Display', ip.Results.ShowHistogram, 'ModeRatio', ip.Results.ModeRatio);
                
            case 'bgproj'
                % load max. 100 frames
                frameRange = unique(round(linspace(1, data(i).movieLength, 100)));
                aip = zeros(data(i).imagesize);
                mproj = zeros(data(i).imagesize);
                
                if iscell(data(i).framePaths{ch})
                    parfor f = 1:numel(frameRange)
                        frame = double(imread(data(i).framePaths{ch}{frameRange(f)})); %#ok<PFBNS>
                        dmask = 0~=double(imread(data(i).maskPaths{frameRange(f)}));
                        
                        dmask = imdilate(dmask, se);
                        frame(dmask) = 0;
                        aip = aip + frame;
                        mproj = mproj + dmask;
                    end
                else
                    info = imfinfo(data(i).framePaths{ch});
                    minfo = imfinfo(data(i).maskPaths);
                    parfor f = 1:numel(frameRange)
                        frame = double(readtiff(data(i).framePaths{ch}, frameRange(f), info)); %#ok<PFBNS>
                        dmask = 0~=double(readtiff(data(i).maskPaths, frameRange(f), minfo));
                        dmask = imdilate(dmask, se);
                        frame(dmask) = 0;
                        aip = aip + frame;
                        mproj = mproj + dmask;
                    end
                end
                
                aip = aip./(numel(frameRange)-mproj);
                % fill to enable filtering
                aip(isnan(aip)) = prctile(aip(:), 95);                  
                
                mask{i} = maskFromFirstMode(aip, 'Connect', ip.Results.Connect,...
                    'Display', ip.Results.ShowHistogram, 'ModeRatio', ip.Results.ModeRatio);
                
                proj{i} = aip;
                
            case 'density'
                error('Density-based mode not implemented.');
    
        end
        mask{i} = imfill(mask{i}, 'holes');
    else
        mask{i} = double(imread(maskPath));
        proj{i} = double(imread(projPath));
    end
end

% Ask user to inspect and potentially manually adjust masks
if ip.Results.Validate && any(~hasStoredMask)
    vmask = manualSegmentationTweakGUI(cellfun(@double, proj(~hasStoredMask), 'unif', 0), mask(~hasStoredMask));
    mask(~hasStoredMask) = vmask;
end

% save masks
for i = 1:numel(data)
    if ~hasStoredMask(i)
        [~,~] = mkdir([data(i).source 'Detection']);
        maskPath = [data(i).source 'Detection' filesep 'cellmask.tif'];
        imwrite(uint8(mask{i}), maskPath, 'tif', 'compression' , 'lzw');
        
        projPath = [data(i).source 'Detection' filesep 'maxproj.tif'];
        imwrite(proj{i}, projPath, 'tif', 'compression' , 'lzw');
    end
end


if ip.Results.Display
    for i = 1:nd
        if ~isempty(mask{i})
            [ny,nx] = size(mask{i});
            B = bwboundaries(mask{i});
            B = cellfun(@(i) sub2ind([ny nx], i(:,1), i(:,2)), B, 'unif', 0);
            B = cell2mat(B);
            bmask = zeros([ny nx]);
            bmask(B) = 1;
            bmask = bwmorph(bmask, 'dilate');
            
            iproj = imread([data(i).source 'Detection' filesep 'maxproj.tif']);
            iproj = scaleContrast(log(double(iproj)));            
            iproj(bmask==1) = 0;
            overlay = iproj;
            overlay(bmask==1) = 255;
            overlay = uint8(cat(3, overlay, iproj, iproj));
            figure; imagesc(overlay); axis image;
        end
    end
end

if nd==1
    mask = mask{1};
    proj = proj{1};
end
