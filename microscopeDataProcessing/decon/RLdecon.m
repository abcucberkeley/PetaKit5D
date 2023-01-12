function [deconvolved, ds, dsr] = RLdecon(inputFn, outputFn, PSFfn, xyPixelSize, dz, dzPSF, varargin)
% new RL decon framework for multiple decon methods, and includes deskew/rotate after decon
% 
% Author: Xiongtao Ruan (11/12/2022)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFn');
ip.addRequired('outputFn');
ip.addRequired('PSFfn');
ip.addRequired('xyPixelSize', @isnumeric);
ip.addRequired('dz', @isnumeric);
ip.addRequired('dzPSF', @isnumeric);
ip.addParameter('rawdata', [], @isnumeric); 
ip.addParameter('nTapering', 0, @isnumeric); 
ip.addParameter('Save16bit', false , @islogical);
ip.addParameter('ObjectiveScan', false, @islogical);
% deskew and rotation options
ip.addParameter('Deskew', false , @islogical);
ip.addParameter('Rotate', false , @islogical);
ip.addParameter('DSRCombined', true , @islogical);
ip.addParameter('Reverse', false, @islogical);
ip.addParameter('SkewAngle', 32.45 , @isnumeric);
ip.addParameter('flipZstack', false, @islogical); 
ip.addParameter('Interp', 'linear', @isnumeric); 
ip.addParameter('Crop', false, @islogical); 
ip.addParameter('xStepThresh', 2, @isnumeric); 
% decon parameters
ip.addParameter('Background', [], @isnumeric);
ip.addParameter('Resample', [] , @isnumeric);
ip.addParameter('DeconIter', 15 , @isnumeric); % number of iterations
ip.addParameter('RLMethod', 'omw' , @ischar); % rl method {'original', 'simplified', 'omw', 'cudagen'}
ip.addParameter('wienerAlpha', 0.005, @isnumeric); % alpha for wiener in OMW method
ip.addParameter('OTFCumThresh', 0.9, @isnumeric); % OTF cumutative sum threshold
ip.addParameter('skewed', [], @(x) isempty(x) || islogical(x)); % decon in skewed space
ip.addParameter('fixIter', true, @islogical);
ip.addParameter('errThresh', [], @isnumeric); % error threshold for simplified code
ip.addParameter('saveZarr', false, @islogical); % save as zarr
ip.addParameter('BlockSize', [256, 256, 256], @isnumeric); % block overlap
ip.addParameter('scaleFactor', [], @isnumeric); % scale factor for result
ip.addParameter('deconBbox', [], @isnumeric); % bounding box to crop data after decon
ip.addParameter('debug', false, @islogical);
ip.addParameter('saveStep', 5, @isnumeric); % save intermediate results every given iterations
ip.addParameter('useGPU', true, @islogical); % use GPU processing
ip.addParameter('save3Dstack', [true, true, true], @(x) islogical(x) && (numel(x) == 1 || numel(x) == 3)); % save 3d stack to disk
ip.addParameter('mipAxis', [0,0,1] , @(x) isnumeric(x) && numel(x) == 3);
ip.addParameter('psfGen', true, @islogical); % psf generation
ip.addParameter('uuid', '', @ischar);

ip.parse(inputFn, outputFn, PSFfn, xyPixelSize, dz, dzPSF, varargin{:});

pr = ip.Results;
rawdata = pr.rawdata;
pr.rawdata = [];
nTapering = pr.nTapering;
Save16bit = pr.Save16bit;
ObjectiveScan = pr.ObjectiveScan;
% deskew and rotation
Deskew = pr.Deskew;
DSRCombined = pr.DSRCombined;
Rotate = pr.Rotate;
Reverse = pr.Reverse;
SkewAngle = pr.SkewAngle;
flipZstack = pr.flipZstack;
Interp = pr.Interp;
xStepThresh = pr.xStepThresh;
% decon
Background = pr.Background;
Crop = pr.Crop;
mipAxis = pr.mipAxis;
Resample = pr.Resample;
nIter = pr.DeconIter;
RLMethod = pr.RLMethod;
wienerAlpha = pr.wienerAlpha;
OTFCumThresh = pr.OTFCumThresh;
skewed = pr.skewed;
fixIter = pr.fixIter;
errThresh = pr.errThresh;
saveZarr = pr.saveZarr;
BlockSize = pr.BlockSize;
scaleFactor = pr.scaleFactor;
deconBbox = pr.deconBbox;
debug = pr.debug;
saveStep = pr.saveStep;
useGPU = pr.useGPU;
save3Dstack = pr.save3Dstack;
psfGen = pr.psfGen;
uuid = pr.uuid;

SkewAngle = abs(SkewAngle);
if ~ObjectiveScan
    dz_ratio = sind(SkewAngle);
else
    dz_ratio = 1;
end
zAniso = dz_ratio * dz / xyPixelSize;

if isempty(uuid)
    uuid = get_uuid();
end
loc_uuid = get_uuid();

if ~isempty(inputFn)
    [dataPath, fsname] = fileparts(inputFn);
end

if ~isempty(outputFn)
    [deconPath, outFsn] = fileparts(outputFn);    
else
    if ~isempty(inputFn)
        deconPath = [dataPath, '/matlab_decon/'];    
    end
end

if saveZarr
    ext = '.zarr';
else
    ext = '.tif';
end

% order: decon, ds, dsr/rotated
file_exist_mat = false(3, 1);
if numel(save3Dstack) == 1
    save3Dstack = repmat(save3Dstack, 3, 1);
end
if Deskew 
    if ~DSRCombined
        dsPath = sprintf('%s/DS/', deconPath);
        dsFn = sprintf('%s/%s%s', dsPath, outFsn, ext);
    end
    
    if Rotate
        dsrPath = sprintf('%s/DSR/', deconPath);
        dsrFn = sprintf('%s/%s%s', dsrPath, outFsn, ext);
    end
end

if ObjectiveScan
    Deskew = false;
    if Rotate
        dsrPath = sprintf('%s/Rotated/', deconPath);    
        dsrFn = sprintf('%s/%s%s', dsrPath, outFsn, ext);
    end
end

if ~Deskew && Rotate && ~ObjectiveScan
    Deskew = true;
    DSRCombined = true;
end

if ~Deskew || DSRCombined
    save3Dstack(2) = false;
    file_exist_mat(2) = true;
end

if ~Rotate
    save3Dstack(3) = false;
    file_exist_mat(3) = true;    
end

% check if the files need to be saved and exists
if save3Dstack(1)
    if exist(outputFn, 'file')
        file_exist_mat(1) = true;
    end
end
if save3Dstack(2)
    if exist(dsFn, 'file')
        file_exist_mat(2) = true;
    end
end
if save3Dstack(3)
    if exist(dsrFn, 'file')
        file_exist_mat(3) = true;
    end
end

deconvolved = [];
ds = [];
dsr = [];

% if all file exist, or not save decon, but save ds or dsr, and they
% already exist, then skip the computing.
if all(file_exist_mat) || (~save3Dstack(1) && (Deskew || Rotate) && all(file_exist_mat(2 : 3)))
    fprintf('The output file(s) %s and related files already exist, skip the processing!\n', outputFn);
    return;
end

if ~isempty(inputFn)
    imSize = getImageSize(inputFn);
else
    imSize = size(rawdata);
end

% prepare for PSF
if psfGen || ~isempty(inputFn) || strcmp(RLMethod, 'omw')
    psfgenPath = sprintf('%s/psfgen/', deconPath);
    if ~exist(psfgenPath, 'dir')
        mkdir(psfgenPath);
    end
end

if ischar(PSFfn)
    [~, psfFsn, suffix]=fileparts(PSFfn);
    if strcmp(suffix, '.mat')
        load(PSFfn, 'psf');
    elseif strcmp(suffix, '.tif')
        % xruan (01/12/2021)
        if psfGen || ~isempty(inputFn)
            psfgen_filename = sprintf('%s/%s.tif', psfgenPath, psfFsn);
        end

        if psfGen && exist(psfgen_filename, 'file')
            fprintf('Load existing generated PSF %s for %s ...\n', psfgen_filename, PSFfn);   
            try 
                psf = double(readtiff(psfgen_filename));
                psfGen = false;
            catch ME
                disp(ME);
                fprintf('Load PSF %s ...\n', PSFfn);                                            
                psf = double(readtiff(PSFfn));
            end
        else
            fprintf('Load PSF %s ...\n', PSFfn);                                            
            psf = double(readtiff(PSFfn));
        end

        try 
            % xruan (05/05/2021) change to psf_gen_new
            if psfGen
                fprintf('PSF generation for %s ...\n', PSFfn);

                medFactor = 1.5;
                PSFGenMethod = 'masked';
                psf = psf_gen_new(psf, dzPSF, dz * dz_ratio, medFactor, PSFGenMethod);
                
                % crop psf to the bounding box (-/+ 1 pixel) and make sure the
                % center doesn't shift
                py = find(squeeze(sum(psf, [2, 3])));
                px = find(squeeze(sum(psf, [1, 3])));
                pz = find(squeeze(sum(psf, [1, 2])));
                cropSz = [min(py(1) - 1, size(psf, 1) - py(end)), min(px(1) - 1, size(psf, 2) - px(end)), min(pz(1) - 1, size(psf, 3) - pz(end))] - 1;
                cropSz = max(0, cropSz);
                bbox = [cropSz + 1, size(psf) - cropSz];
                psf = psf(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
            end
            
            % crop psf if it is larger than data in any dimension
            if any(size(psf) > imSize)
                warning('The psf size %s is larger than the data size %s, crop it!', mat2str(size(psf)), mat2str(imSize));
                s = max(0, floor((size(psf) - imSize) / 2)) + 1;
                t = s + min(size(psf), imSize) - 1;
                psf = psf(s(1) : t(1), s(2) : t(2), s(3) : t(3));
            end
            
            if ~isempty(inputFn)
                if ~exist(psfgen_filename, 'file')
                    tmp_filename = sprintf('%s/%s_%s.tif', psfgenPath, psfFsn, loc_uuid);
                    writetiff(psf, tmp_filename);
                    if ~exist(psfgen_filename, 'file')
                        movefile(tmp_filename, psfgen_filename);
                    end
                end
            end
        catch ME
            disp(ME)
        end
    else
        psf = [];  % dummy PSF
    end
end

% generate or load back projectors
switch RLMethod
    case 'omw'
        bpFn = sprintf('%s/%s_back_projector_alpha_%0.6f_otf_cum_thresh_%0.6f.tif', psfgenPath, psfFsn, wienerAlpha, OTFCumThresh);
        if exist(bpFn, 'file')
            fprintf('Load existing OMW back projector %s for %s ...\n', bpFn, PSFfn);                    
            psf_b = readtiff(bpFn);
        else
            fprintf('OMW back projector generation for %s ...\n', PSFfn);        
            
            bpTmpFn = sprintf('%s/%s_back_projector_alpha_%0.6f_otf_cum_thresh_%0.6f_%s.tif', psfgenPath, psfFsn, wienerAlpha, OTFCumThresh, loc_uuid);
            [psf_b, OTF_bp_omw, abs_OTF_c, OTF_mask] = omw_backprojector_generation(psf, wienerAlpha, skewed, 'OTFCumThresh', OTFCumThresh);

            if usejava('jvm')
                fig = visualize_OTF_and_mask_outline(abs_OTF_c, OTF_mask);
                figFn = sprintf('%s/%s_back_projector_alpha_%0.6f_otf_cum_thresh_%0.6f_figure.png', psfgenPath, psfFsn, wienerAlpha, OTFCumThresh);
                print(fig, figFn, '-dpng', '-r0');
                close(fig);
            end

            % write back projector
            writetiff(psf_b, bpTmpFn);
            movefile(bpTmpFn, bpFn);
            
            % write OTF
            bpFn = sprintf('%s/%s_OTF.tif', psfgenPath, psfFsn);
            bpTmpFn = sprintf('%s/%s_OTF_%s.tif', psfgenPath, psfFsn, loc_uuid);
            writetiff(abs_OTF_c, bpTmpFn);
            movefile(bpTmpFn, bpFn);

            % write OTF mask 
            bpFn = sprintf('%s/%s_OTF_mask_otf_cum_thresh_%0.6f.tif', psfgenPath, psfFsn, OTFCumThresh);
            bpTmpFn = sprintf('%s/%s_OTF_mask_otf_cum_thresh_%0.6f_%s.tif', psfgenPath, psfFsn, OTFCumThresh, loc_uuid);
            writetiff(uint8(OTF_mask), bpTmpFn);
            movefile(bpTmpFn, bpFn);
        end
        if any(size(psf_b) > imSize)
            warning('The psf size %s is larger than the data size %s, crop it!', mat2str(size(psf_b)), mat2str(imSize));
            s = max(0, floor((size(psf_b) - imSize) / 2)) + 1;
            t = s + min(size(psf_b), imSize) - 1;
            psf_b = psf_b(s(1) : t(1), s(2) : t(2), s(3) : t(3));
        end
end

if isempty(rawdata)
    [~, ~, ext] = fileparts(inputFn);
    switch ext
        case {'.tif', '.tiff'}
            rawdata = readtiff(inputFn);
        case '.zarr'
            rawdata = readzarr(inputFn);
    end
end

% add support for flip z stack
if flipZstack
    rawdata = flip(rawdata, 3);
end

nz = size(rawdata, 3);
% soften the x-y borders if requested
if nTapering > 0
    taperKernel = fspecial('gaussian', nTapering+30, nTapering);
    for iz=1:nz
        rawdata(:,:,iz) = edgetaper(rawdata(:,:,iz), taperKernel);
    end
end
% rawdata(rawdata<0.0) = 0.0;
% for omw and simplified RL method, subtrack background within the decon function 
if ~strcmp(RLMethod, 'omw') && ~strcmp(RLMethod, 'simplified')
    rawdata = max(single(rawdata) - Background, 0);        
end

% call Richardson-Lucy
if ~file_exist_mat(1)
    if nIter>0 && any(rawdata, 'all')
        fprintf('Deconvolution for frame %s with %s Method ...\n', inputFn, RLMethod);            
        if isempty(scaleFactor)
            % scaleFactor = numel(rawdata);
            scaleFactor = 1;
        end
        psf = psf ./ sum(psf(:));
        if debug
            debug_folder = sprintf('%s/%s_debug/', deconPath, fsname);
            mkdir(debug_folder)
            saveStep = min(saveStep, nIter);
        else
            debug_folder = '/tmp/';
        end
        
        switch RLMethod 
            case 'original'
                deconvolved = deconvlucy(rawdata, psf, nIter);
            case 'simplified'
                [deconvolved, err_mat, iter_run] = decon_lucy_function(rawdata, psf, nIter, Background, useGPU, Save16bit, deconBbox, debug, debug_folder, saveStep);
            case 'omw'
                [deconvolved, err_mat] = decon_lucy_omw_function(rawdata, psf, psf_b, nIter, Background, useGPU, Save16bit, deconBbox, debug, debug_folder, saveStep);          
            case 'cudagen'
                deconvolved = decon_lucy_cuda_function(single(rawdata), single(psf), nIter);            
        end
        if scaleFactor ~= 1
            deconvolved = deconvolved * scaleFactor;
        end
    else
        deconvolved = rawdata;
    end
elseif Deskew || Rotate
    % load existing decon results
    [~, ~, ext] = fileparts(outputFn);
    switch ext
        case {'.tif', '.tiff'}
            deconvolved = single(readtiff(outputFn));
        case '.zarr'
            deconvolved = single(readzarr(outputFn));
    end
end

clear rawdata;

if ObjectiveScan
    Deskew = false;
end

% Deskew and Rotate the decon result if requested
if Deskew && (~Rotate || ~DSRCombined)
    if ~file_exist_mat(2)
        fprintf('Deskew deconvolved frame %s...\n', fsname);
        if ~exist(dsPath, 'dir')
            mkdir(dsPath);
            if ~ispc                
                fileattrib(dsPath, '+w', 'g');
            end        
        end
    
        ds = deskewFrame3D(deconvolved, SkewAngle, dz, xyPixelSize, Reverse, ...
            'Crop', Crop, 'Interp', Interp); 
            
        MIPFn = sprintf('%s/MIPs/%s_MIP_z.tif', dsPath, fsname);
        saveMIP_frame(ds, MIPFn, 'axis', mipAxis);
        
        if save3Dstack(2)
            if saveZarr 
                dsTmpPath = sprintf('%s%s_%s.zarr', dsPath, fsname, uuid);
                if Save16bit
                    writezarr(uint16(ds), dsTmpPath, 'blockSize', BlockSize);            
                else
                    writezarr(ds, dsTmpPath, 'blockSize', BlockSize);
                end
            else
                dsTmpPath = sprintf('%s%s_%s.tif', dsPath, fsname, uuid);                
                writetiff(ds, dsTmpPath);
                if Save16bit
                    writetiff(uint16(ds), dsTmpPath);
                else
                    writetiff(ds, dsTmpPath);
                end        
            end
            movefile(dsTmpPath, dsFn);
        end
    elseif Rotate && ~file_exist_mat(3)
        [~, ~, ext] = fileparts(outputFn);
        switch ext
            case {'.tif', '.tiff'}
                ds = single(readtiff(dsFn));
            case '.zarr'
                ds = single(readzarr(dsFn));
        end
    end

    if Rotate
        fprintf('Rotate deskewed deconvolved frame %s...\n', fsname);        
        dsr = rotateFrame3D(ds, SkewAngle, zAniso, Reverse,...
            'Crop', true, 'ObjectiveScan', ObjectiveScan, 'Interp', Interp);
        ds = [];
    end
elseif Deskew && Rotate && DSRCombined
    fprintf('Deskew, rotate and resample for deconvolved frame %s...\n', fsname);                
    dsr = deskewRotateFrame3D(deconvolved, SkewAngle, dz, xyPixelSize, ...
        'reverse', Reverse, 'Crop', true, 'ObjectiveScan', ObjectiveScan, ...
        'resample', Resample, 'Interp', Interp, 'xStepThresh', xStepThresh);
end

if ObjectiveScan && Rotate    
    fprintf('Rotate deconvolved frame %s...\n', fsname);                
    dsr = rotateFrame3D(deconvolved, SkewAngle, zAniso, Reverse,...
        'Crop', true, 'ObjectiveScan', ObjectiveScan, 'Interp', Interp);
end

if Rotate 
    if Save16bit
        dsr = uint16(dsr);
    else
        dsr = single(dsr);
    end

    MIPFn = sprintf('%s/MIPs/%s_MIP_z.tif', dsrPath, fsname);
    saveMIP_frame(dsr, MIPFn, 'axis', mipAxis);

    if save3Dstack(3)
        if saveZarr 
            dsrTmpPath = sprintf('%s%s_%s.zarr', dsrPath, fsname, uuid);
            writezarr(dsr, dsrTmpPath, 'blockSize', BlockSize);
        else
            dsrTmpPath = sprintf('%s%s_%s.tif', dsrPath, fsname, uuid);
            writetiff(dsr, dsrTmpPath);
        end
        movefile(dsrTmpPath, dsrFn);
    end
end

% if the output file is empty, directly return the deconvolved results. 
if isempty(outputFn)
    return;
end

% not save 3D stack
if ~any(save3Dstack)
    return;
end

% save err mat for simplified method
if nIter > 0 && strcmp(RLMethod, 'simplified')
    infoPath = sprintf('%s/%s/decon_info/', deconPath);
    if ~exist(infoPath, 'dir')
        mkdir(infoPath);
    end
    infoFn = sprintf('%s/%s_info.mat', infoPath, outputFsn);
    save('-v7.3', infoFn, 'err_mat', 'nIter', 'fixIter', 'errThresh', 'debug');
    save(sprintf('%sactual_iterations_%d.txt', infoFn(1 : end - 8), iter_run), 'iter_run', '-ASCII');
end

if Save16bit
    deconvolved = uint16(deconvolved);
end

if length(mipAxis) ~= 3
    warning('Exactly 3 parameters are needed for -M')
else
    if any(mipAxis)
        MIPFn = sprintf('%s/MIPs/%s_MIP_z.tif', deconPath, fsname);
        saveMIP_frame(deconvolved, MIPFn, 'axis', mipAxis);
    end
end

if save3Dstack(1)
    if saveZarr 
        deconTmpPath = sprintf('%s_%s.zarr', outputFn(1 : end - 4), uuid);                
        writezarr(deconvolved, deconTmpPath, 'blockSize', BlockSize);
    else
        deconTmpPath = sprintf('%s_%s.tif', outputFn(1 : end - 4), uuid);        
        writetiff(deconvolved, deconTmpPath);
    end
    movefile(deconTmpPath, outputFn);
end

end

