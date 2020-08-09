%getCellVolume(data, varargin) calculates cell volume & area based on thresholding
% of the smoothened input data
%
% Inputs: 
%    data : structure returned by loadConditionData3D
% 
% Parameters (specifier/value pairs):
%  'SmoothingSigma' : s.d. for Gaussian smoothing of the data
%       'MinVolume' : minimum accepted size of connected components, in voxels.
%                     This is used for eliminating spurious detections, for example other
%                     cells at the borders of the input volume
%            'Mode' : {'Average'}|'Timeseries' specifies whether to calculate volume based
%                      on an average (max. int. projection) of the data, or for all time points
%       'FillHoles' : {true}|false specifies whether to fill holes at the interior of the
%                      thresholded volume automatically
%        'Filetype' : 'framePathsDS' | {'framePathsDSR'} selects whether to use 
%                     deskewed (i.e. with objective-scan data) or rotated volumes
%
% Outputs:
%   In 'Average' mode, the files 'maxProjRotated.tif' and 'volmaskRotated.tif' are 
%   written to [data(i).source 'Analysis']
%
%   In 'TimeSeries' mode, volume masks are written to [data(i).source 'Analysis/VolMasks']
%   and volume/area statistics are saved in [data(i).source 'Analysis/voldata.mat']

% Author: Francois Aguet (2014)

function getCellVolume(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data');
ip.addParamValue('SmoothingSigma', 2, @isnumeric);
ip.addParamValue('Display', true, @islogical); % for time series only
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('FillHoles', true, @islogical);
ip.addParamValue('ExcludeIndex', [], @isposint);
ip.addParamValue('MinVolume', 200000, @isscalar);
ip.addParamValue('Mode', 'Average', @(x) any(strcmpi(x, {'Average', 'TimeSeries'})));
ip.addParamValue('Filepath', 'framePathsDSR', @(x) any(strcmpi(x, {'framePathsDS', 'framePathsDSR'})));
ip.addParamValue('ThresholdMode', 'Otsu', @(x) any(strcmpi(x, {'Background', 'Otsu'})));
ip.parse(data, varargin{:});

nd = numel(data);

if strcmpi(ip.Results.Mode, 'TimeSeries')
    for i = 1:nd
        nf = data(i).movieLength;
        fmt = ['%.' num2str(ceil(log10(nf))) 'd'];
        
        mpath = [data(i).source 'Analysis' filesep 'VolMasks' filesep];
        rfile = [data(i).source 'Analysis' filesep 'voldata.mat'];
        
        if ~(exist(rfile,'file')==2) || ip.Results.Overwrite
            [~,~] = mkdir(mpath);
            
            % define interpolation grids
            info = imfinfo(data(i).(ip.Results.Filepath){1}{1}, 'tif');
            nz = numel(info);
            nx = info(1).Width;
            ny = info(1).Height;
            [y,x,z] = ndgrid(1:ny,1:nx,1:nz);
            [yi,xi,zi] = ndgrid(1:ny,1:nx,1:1/data(i).zAniso:nz);
            
            px = data(i).pixelSize; % [um]
            
            vol = zeros(1,nf);
            area = zeros(1,nf);
            parfor f = 1:nf
                % determine background level from volume borders
                stack = readtiff(data(i).(ip.Results.Filepath){1}{f}); %#ok<PFBNS>
                gstack = filterGauss3D(stack, ip.Results.SmoothingSigma);
                
                if strcmpi(ip.Results.ThresholdMode, 'Otsu');
                    % threshold smoothened volume
                    T = thresholdOtsu(gstack);
                else
                    % quick hack, not an optimal method
                    b = 10; % distance from border
                    z0 = 3; % first and last 3 slices have interpolation noise
                    tmp = stack([1:1+b end-b:end], [1:1+b end-b:end], 1+z0:end-z0);
                    T = prctile(tmp(:), 99);
                end
                volMask = gstack>T;
                
                % fill potential holes inside the cell mask (slow)
                % this version is saved; not used for actual volume calculation
                if ip.Results.FillHoles
                    volMask = fillHoles2D(volMask);
                end
                CC = bwconncomp(volMask, 26);
                % discard small components (assumed to be noise or debris on glass slide)
                csize = cellfun(@numel, CC.PixelIdxList);
                idx = csize>=ip.Results.MinVolume;
                CC.NumObjects = sum(idx);
                CC.PixelIdxList = CC.PixelIdxList(idx);
                volMask = labelmatrix(CC)~=0;
                
                % store mask
                writetiff(uint8(volMask), [mpath 'volmask_' num2str(f,fmt) '.tif']);
                
                % interpolate for surface area (and volume) calculation
                volMask2 = interpn(y, x, z, gstack, yi, xi, zi, 'linear')>T;
                if ip.Results.FillHoles
                    volMask2 = fillHoles2D(volMask2);
                end
                CC = bwconncomp(volMask2, 26);
                % discard small components (assumed to be noise or debris on glass slide)
                csize = cellfun(@numel, CC.PixelIdxList);
                idx = csize>=2*ip.Results.MinVolume;
                CC.NumObjects = sum(idx);
                CC.PixelIdxList = CC.PixelIdxList(idx);
                volMask2 = labelmatrix(CC)~=0;
                                
                % calculate volume & surface area
                vol(f) = sum(cellfun(@numel, CC.PixelIdxList));
                perim = bwperim(volMask2);
                area(f) = sum(perim(:));                
            end
            vol = vol*px^3; %#ok<NASGU>
            area = area*px^2; %#ok<NASGU>
            save(rfile, 'vol', 'area');
        else
            load(rfile);
        end
        
        if ip.Results.Display
            plotMitosisStatistics(data(i), 'ExcludeIndex', ip.Results.ExcludeIndex);
        end
    end
else
    for i = 1:nd
        if ~(exist([data(i).source 'Analysis' filesep 'volmaskRotated.tif'], 'file')==2) || ip.Results.Overwrite
            
            % max int. proj. of all frames
            maxProj = readtiff(data(i).(ip.Results.Filepath){1}{1});
            for f = 2:data(i).movieLength
                maxProj = max(maxProj, readtiff(data(i).(ip.Results.Filepath){1}{f}));
            end
            gstack = filterGauss3D(double(maxProj), ip.Results.SmoothingSigma);
            
            T = thresholdOtsu(gstack);
            volMask = gstack>T;
            if ip.Results.FillHoles
                volMask = fillHoles2D(volMask);
            end
            CC = bwconncomp(volMask, 26);
            % discard small components (assumed to be noise or debris on glass slide)
            csize = cellfun(@numel, CC.PixelIdxList);
            idx = csize>=ip.Results.MinVolume;
            CC.NumObjects = sum(idx);
            CC.PixelIdxList = CC.PixelIdxList(idx);
            volMask = labelmatrix(CC)~=0;
            
            % store mask
            writetiff(maxProj, [data(i).source 'Analysis' filesep 'maxProjRotated.tif']);
            writetiff(uint8(volMask), [data(i).source 'Analysis' filesep 'volmaskRotated.tif']);
            
            % interpolate for surface area (and volume) calculation
            if ~data(i).objectiveScan
                vol = sum(cellfun(@numel, CC.PixelIdxList));
                perim = bwperim(volMask);
                area = sum(perim(:));
            else % assumes DS input!
                % define interpolation grids
                info = imfinfo(data(i).(ip.Results.Filepath){1}{1}, 'tif');
                nz = numel(info);
                nx = info(1).Width;
                ny = info(1).Height;
                [y,x,z] = ndgrid(1:ny,1:nx,1:nz);
                [yi,xi,zi] = ndgrid(1:ny,1:nx,1:1/data(i).zAniso:nz);
                
                volMask2 = interpn(y, x, z, gstack, yi, xi, zi, 'linear')>T;
                if ip.Results.FillHoles
                    volMask2 = fillHoles2D(volMask2);
                end
                CC = bwconncomp(volMask2, 26);
                % discard small components (assumed to be noise or debris on glass slide)
                csize = cellfun(@numel, CC.PixelIdxList);
                idx = csize>=2*ip.Results.MinVolume;
                CC.NumObjects = sum(idx);
                CC.PixelIdxList = CC.PixelIdxList(idx);
                volMask2 = labelmatrix(CC)~=0;
                
                % calculate volume & surface area
                vol = sum(cellfun(@numel, CC.PixelIdxList));
                perim = bwperim(volMask2);
                area = sum(perim(:));
            end
            px = data(i).pixelSize; % [um]
            vol = vol*px^3; %#ok<NASGU>
            area = area*px^2; %#ok<NASGU>
            rfile = [data(i).source 'Analysis' filesep 'voldata.mat'];
            save(rfile, 'vol', 'area');
        end
    end
end


function vol = fillHoles2D(vol)

for z = 1:size(vol,3)
    vol(:,:,z) = imfill(vol(:,:,z), 8, 'holes');
end
for x = 1:size(vol,2)
    vol(:,x,:) = imfill(squeeze(vol(:,x,:)), 8, 'holes');
end
for y = 1:size(vol,1)
    vol(y,:,:) = imfill(squeeze(vol(y,:,:)), 8, 'holes');
end
for z = 1:size(vol,3)
    vol(:,:,z) = imfill(vol(:,:,z), 8, 'holes');
end
for x = 1:size(vol,2)
    vol(:,x,:) = imfill(squeeze(vol(:,x,:)), 8, 'holes');
end
for y = 1:size(vol,1)
    vol(y,:,:) = imfill(squeeze(vol(y,:,:)), 8, 'holes');
end
