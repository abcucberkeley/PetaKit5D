function [tiffFullpaths, zarrFullpaths, fsnames, zarrPathstr, overall_z_median] = stitch_process_filenames(tileFullpaths, DS, DSR, Decon, DSRDirstr, DeconDirstr, stitchMIP, ObjectiveScan, zNormalize, px, xf, zf, resample)
% return tile paths and other info based on the processing required for the
% data (i.e., dsr, dsr/decon, decon, raw etc)
% 
% The code is taken from the main stitching function (for simplify), and
% with further improvement
%
% Author: Xiongtao Ruan (07/23/2021)
%
% xruan (07/27/2021): if resample not empty, indicate resample factor in
% the zarr folder
% xruan (10/13/2021): add support for stitching with reference for decon data
% xruan (12/20/2021): add support for mip stitching


[dataPath, fsname] = fileparts(tileFullpaths{1});

if DSR
    % create DSR or Rotated folder
    if ~isempty(DSRDirstr)
        dsrpath = [dataPath filesep DSRDirstr];
    elseif ~isempty(DeconDirstr)
        dsrpath = [dataPath filesep DeconDirstr];
        Decon = true;
    else
        if ObjectiveScan
            dsrpath = [dataPath filesep 'Rotated_dx' num2str(px*xf) '_dz' num2str(px*zf)];        
        else
            dsrpath = [dataPath filesep 'DSR_dx' num2str(px*xf) '_dz' num2str(px*zf)];
        end
    end
    if any(stitchMIP)
        dsrpath = [dsrpath, filesep, 'MIPs', filesep];
        suffix_strs = {'x', 'y', 'z'};
        if numel(stitchMIP) == 1
            suffix_str = 'z';
        else
            suffix_str = suffix_strs{stitchMIP};
        end
    end
    tilePath = dsrpath;
else
    if ~isempty(DeconDirstr)
        tilePath = [dataPath filesep DeconDirstr];
        Decon = true;  
    else
        tilePath = dataPath;
    end
end

nF = numel(tileFullpaths);

if zNormalize
    zmed_cell = cell(nF, 1);
end

if DSR
    for k = 1:nF
        tic
        tileFullpath = tileFullpaths{k};
        [dataPath, fsname] = fileparts(tileFullpath);

        dsrFullpath = [dsrpath filesep fsname '.tif'];
        if Decon && ~isempty(DeconDirstr)
            dsrFullpath = [dsrpath filesep fsname '_decon.tif'];
        end
        if any(stitchMIP)
            dsrFullpath = sprintf('%s/%s_MIP_%s.tif', dsrpath, fsname, suffix_str);         
        end
        
        if ~exist(dsrFullpath, 'file')
            im = readtiff(tileFullpath);

            if ~ObjectiveScan
                fprintf('deskewing Data...')
                tic
                im = deskewFrame3D(im, SkewAngle, dz, xyPixelSize, Reverse,...
                    'Crop', Crop);
                toc
            end

            % deconvovle
            if Decon && isempty(DeconDirstr)
                deconpath = [rt filesep 'matlab_decon'];
                if ~exist(deconpath, 'dir')
                    mkdir(deconpath);
                end

                if ~exist([deconpath filesep fsname '.tif'], 'file')
                    fprintf('Deconvolving Data...')
                    tic
                    volpath = [rt fsname];
                    mask = logical(im);
                    im = RLdecon_for_ExMPipeline(im, volpath, PSFpath, Background, nIter, dzPSF, zAniso, 0, [], 0, ...
                        px, 0, 0, 0, 0, [0,0,1],0, []) ;toc

                    % remove edge artifacts
                    fprintf('Removing Edges...'); tic
                    r = EdgeArtifacts;
                    if r > 0
                        se = strel('disk', r);
                        mask = imerode(mask,se);
                        mask(:,:,1:r) = 0;
                        mask(:,:,end-r:end) = 0;
                        im(~mask) = 0;
                        clear mask
                        toc
                    end

                    if Save16bit
                        writetiff(uint16(im), [deconpath filesep fsname '_' uuid '.tif']);
                    else
                        writetiff(im, [deconpath filesep fsname '_' uuid '.tif']);
                    end
                    movefile([deconpath filesep fsname '_' uuid '.tif'], [deconpath filesep fsname '.tif']);
                end
            end
            fprintf('Rotating Data...')
            tic
            dsr = rotateFrame3D(im, SkewAngle, zAniso, Reverse,...
                'Crop', true, 'ObjectiveScan', ObjectiveScan);
            toc

            fprintf('resampling Rotated Data...')
            if any([xf, yf, zf] ~= 1)
                tic
                dsr = GU_resampleStack3D(dsr, xf, yf, zf,'Interp', 'linear');toc
            end

            if ~exist(dsrpath, 'dir')
                mkdir(dsrpath);
            end

            fprintf('saving resampled + Rotated Data...')
            tic
            writetiff(uint16(dsr), dsrFullpath); toc
        else
            % only load data when zNormalize is true
            if zNormalize
                fprintf('loading resampling Rotated Data...')
                tic
                dsr = readtiff(dsrFullpath); toc
            end
        end

        if zNormalize
            zmed = zeros(size(dsr, 3), 1);
            for z = 1 : size(dsr, 3)
                dsr_z = dsr(:, :, z);
                if any(dsr_z ~= 0, 'all')
                    zmed(z) = median(dsr_z(dsr_z ~= 0));
                end
            end
            zmed_cell{k} = zmed;
        end    
    end
end

overall_z_median = 0;
if zNormalize
    zmed_mat = cat(2, zmed_cell{:});
    overall_z_median = median(zmed_mat(zmed_mat ~= 0));
end

% convert tif to zarr files (first assume the files exist)
tiffFullpaths = cell(nF, 1);
zarrFullpaths = cell(nF, 1);
fsnames = cell(nF, 1);
for i = 1 : nF
    tileFullpath = tileFullpaths{i};
    [dataPath, fsname] = fileparts(tileFullpath);
    fsnames{i} = fsname;
    
    if ~isempty(resample)
        % complete resample to 3d
        rs = [ones(1, 4 - numel(resample)) * resample(1), resample(2 : end)];    
        zarrPathstr = sprintf('zarr_%d_%d_%d', rs(1), rs(2), rs(3));
    else
        zarrPathstr = 'zarr';
    end
    if Decon && ~isempty(DeconDirstr)
        tiffFullpaths{i} = sprintf('%s/%s_decon.tif', tilePath, fsname);
        zarrFullpaths{i} = sprintf('%s/%s/%s_decon.zarr', tilePath, zarrPathstr, fsname);        
    elseif any(stitchMIP)
        tiffFullpaths{i} = sprintf('%s/%s_MIP_%s.tif', tilePath, fsname, suffix_str);
        zarrFullpaths{i} = sprintf('%s/%s/%s_MIP_%s.zarr', tilePath, zarrPathstr, fsname, suffix_str);                
    else
        tiffFullpaths{i} = sprintf('%s/%s.tif', tilePath, fsname);
        zarrFullpaths{i} = sprintf('%s/%s/%s.zarr', tilePath, zarrPathstr, fsname);
    end
end


end

