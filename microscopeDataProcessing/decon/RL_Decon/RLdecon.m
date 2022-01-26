function [deconvolved] = RLdecon(input_tiff, output_filename, psf, background, nIter, dz_psf, dz_data, ...
    zalign, zalignParams, rotateByAngle, xypixelsize, bRotFinal, ...
    bSaveUint16, cropFinal, bFlipZ, axesMaxIntProj, resizeFactor, scalingThresh, ...
    RLMethod, fixIter, errThresh, flipZstack, debug, varargin)
%RLdecon Summary of this function goes here
%   input_tiff: input TIFF file name
%   psf: psf array in 'double'
%   background: backgroud to subtract
%   nIter: number of iterations
%   zalign: is 1 when de-skewing is needed for sample-scan data
%   zalignParams: parameters related to de-skewing [extrashift reverse mean saveAligned]
%   rotateByAngle: angle used in deskewing and rotating final results
%   bRotFinal: 1 to rotate the final result to make horizontal plane
%   correspond to coverslip plane
%   CropFinal: [lower upper] to crop the final result in Z; if
%   upper==lower then no cropping will be done
% 
% xruan (01/10/2021): for conversion to uint16, first rescale intensity if the max intensity is over 2^16 -1. 
% xruan (01/12/2021): for conversion to uint16, add support for smaller psf
% xruan (03/25/2021): add options for different versions of rl method
% xruan (03/26/2021): change loadtiff to readtiff
% xruan (06/10/2021): add support for threshold and debug mode in simplified version. 
% xruan (06/16/2021): add support for saving generated psf, and change default psfgen method as masked
% xruan (06/22/2021): add support to save err mat for simplified version.
% xruan (07/13/2021): add support for the processing of flipped files;
% first save to intermediate file and then move the final file; crop psf if
% it is larger than the data in any dimension.
% xruan (07/15/2021): add support for zarr input
% xruan (11/11/2021): add support for using the data as input instead of
% filename and also user defined normalization factor for the result
% xruan (11/23/2021): normalize psf by the sum so the decon scale remain
% the normal range. 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('input_tiff');
ip.addRequired('output_filename');
ip.addRequired('psf');
ip.addRequired('background');
ip.addRequired('nIter');
ip.addRequired('dz_psf');
ip.addRequired('dz_data');
ip.addRequired('zalign');
ip.addRequired('zalignParams');
ip.addRequired('rotateByAngle');
ip.addRequired('xypixelsize');
ip.addRequired('bRotFinal');
ip.addRequired('bSaveUint16');
ip.addRequired('cropFinal');
ip.addRequired('bFlipZ');
ip.addRequired('axesMaxIntProj');
ip.addRequired('resizeFactor');
ip.addRequired('scalingThresh');
ip.addRequired('RLMethod');
ip.addRequired('fixIter');
ip.addRequired('errThresh');
ip.addRequired('flipZstack');
ip.addRequired('debug');
ip.addParameter('nTapering', [], @isnumeric); 
ip.addParameter('rawdata', [], @isnumeric); 
ip.addParameter('scaleFactor', [], @isnumeric); % scale factor for result
ip.addParameter('useGPU', true, @islogical); % use GPU processing
ip.addParameter('save3Dstack', true, @islogical); % use GPU processing


ip.parse(input_tiff, output_filename, psf, background, nIter, dz_psf, dz_data, ...
    zalign, zalignParams, rotateByAngle, xypixelsize, bRotFinal, ...
    bSaveUint16, cropFinal, bFlipZ, axesMaxIntProj, resizeFactor, scalingThresh, ...
    RLMethod, fixIter, errThresh, flipZstack, debug, varargin{:});

pr = ip.Results;

if ischar(dz_psf)
    dz_psf=str2double(dz_psf);
end
if ischar(dz_data)
    dz_data=str2double(dz_data);
end

if ischar(zalign)
    bZAlign = logical(str2double(zalign));
else
    bZAlign = logical(zalign);
end

if ischar(rotateByAngle)
    rotateByAngle = str2double(rotateByAngle);
end

if ischar(xypixelsize)
    xypixelsize = str2double(xypixelsize);
end

bReverse = false;
if bZAlign
    if ischar(zalignParams)
        zalignParams = str2num(zalignParams);
    end
    outputWidth = zalignParams(1);
    extraShift = zalignParams(2);
    bReverse = logical(zalignParams(3));
    fillval = zalignParams(4);
    nphases = zalignParams(5);
    bSaveRawT = logical(zalignParams(6));
    
    dz_data_ratio = sin(rotateByAngle*pi/180);
else
    dz_data_ratio = 1;
end

if ~isempty(input_tiff)
    [datafolder, inputfile, suffix] = fileparts(input_tiff);
end

if ~isempty(output_filename) && exist(output_filename, 'file')
    return;
end

if ischar(psf)
    [a,b,suffix]=fileparts(psf);
    if strcmp(suffix, '.mat')
        load(psf, 'psf');
    elseif strcmp(suffix, '.tif')
        % xruan (01/12/2021)
        try 
            % psf=psf_gen(psf, dz_psf, dz_data*dz_data_ratio, 48);
            % xruan (05/05/2021) change to psf_gen_new
            try
                pp = parallelReadTiff(psf);
                % pp = readtiff(psf);                
            catch
                pp = readtiff(psf);
            end
            medFactor = 1.5;
            PSFGenMethod = 'masked';
            psf = psf_gen_new(pp, dz_psf, dz_data*dz_data_ratio, medFactor, PSFGenMethod);
            
            % crop psf to the bounding box (-/+ 1 pixel) and make sure the
            % center doesn't shift
            py = find(squeeze(sum(psf, [2, 3])));
            px = find(squeeze(sum(psf, [1, 3])));
            pz = find(squeeze(sum(psf, [1, 2])));
            cropSz = [min(py(1) - 1, size(psf, 1) - py(end)), min(px(1) - 1, size(psf, 2) - px(end)), min(pz(1) - 1, size(psf, 3) - pz(end))] - 1;
            cropSz = max(0, cropSz);
            bbox = [cropSz + 1, size(psf) - cropSz];
            psf = psf(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
            
            % test decon without psf gen
            if false
                psf = double(pp);
            end
            
            % crop psf if it is larger than data in any dimension
            if ~isempty(input_tiff)
                imSize = getImageSize(input_tiff);
            else
                imSize = size(pr.rawdata);
            end
            if any(size(psf) > imSize)
                warning('The psf size %s is larger than the data size %s, crop it!', mat2str(size(psf)), mat2str(imSize));
                s = max(0, floor((size(psf) - imSize) / 2)) + 1;
                t = s + min(size(psf), imSize) - 1;
                psf = psf(s(1) : t(1), s(2) : t(2), s(3) : t(3));
            end
            
            if ~isempty(input_tiff)
                psfgen_folder = sprintf('%s/%s/psfgen/', datafolder, 'matlab_decon');
                mkdir(psfgen_folder);
                psfgen_filename = sprintf('%s/%s.tif', psfgen_folder, b);
                if ~exist(psfgen_filename, 'file')
                    writetiff(psf, psfgen_filename);
                end
            end
        catch ME
            disp(ME)
            % psf=psf_gen(psf, dz_psf, dz_data*dz_data_ratio, 24);
        end
    else
        psf = [];  % dummy PSF
    end
end

if ischar(background)
    background=str2double(background);
end
if ischar(nIter)
    nIter=str2double(nIter);
end

if isempty(RLMethod)
    RLMethod = 'simplied';
end

nTapering = pr.nTapering;
% nTapering = 0;
% for k = 1 : length(varargin)
%     switch k
%         case 1
%             % number of pixel for x-y tapering
%             if ischar(varargin{k})
%                 nTapering = str2num(varargin{k});
%             else
%                 nTapering = varargin{k};
%             end
%         otherwise
%             disp('Unknown varargin index')
%     end
% end

% rawdata = loadtiff(input_tiff);
rawdata = pr.rawdata;
pr.rawdata = [];
if isempty(rawdata)
    [~, ~, ext] = fileparts(input_tiff);
    switch ext
        case {'.tif', '.tiff'}
            try
                rawdata = parallelReadTiff(input_tiff);
                % rawdata = readtiff(input_tiff);                
            catch
                rawdata = readtiff(input_tiff);
            end
        case '.zarr'
            rawdata = readzarr(input_tiff);
    end
end

scaleFactor = pr.scaleFactor;
useGPU = pr.useGPU;
save3Dstack = pr.save3Dstack;

% add support for flip z stack
if flipZstack
    rawdata = flip(rawdata, 3);
end

% rawdata = hpbuster(rawdata, background, 2, 1);
% tic
nz = size(rawdata, 3);

rawdata = single(rawdata) - background;

% soften the x-y borders if requested
if nTapering > 0
    taperKernel = fspecial('gaussian', nTapering+30, nTapering);
    for iz=1:nz
        rawdata(:,:,iz) = edgetaper(rawdata(:,:,iz), taperKernel);
    end
end
rawdata(rawdata<0.0) = 0.0;

% de-skew data if requested
if bZAlign
    size(rawdata)
    rawdata = shear3DinDim2(rawdata, rotateByAngle, bReverse, ...
                dz_data, xypixelsize, 0, outputWidth, extraShift, nphases);
    size(rawdata)
    if bSaveRawT
        if isempty(datafolder)
            if ~exist('deskewed', 'dir')
                mkdir('deskewed');
            end
            write3Dtiff(rawdata, ['deskewed/',inputfile,'_deskewed.tif']);
        else
            if ~exist([datafolder, '/deskewed'], 'dir')
                mkdir([datafolder, '/deskewed']);
            end
            write3Dtiff(rawdata, [datafolder,'/deskewed/',inputfile,'_deskewed.tif']);
        end
%        rawdata = rawdata - background;
        return
    end
end

% call Richardson-Lucy
if nIter>0
    if isempty(scaleFactor)
        % scaleFactor = numel(rawdata);
        scaleFactor = 1;
    end
    psf = psf ./ sum(psf(:));
    switch RLMethod 
        case 'original'
            deconvolved = deconvlucy(rawdata, psf, nIter) * scaleFactor;
        case 'simplified'
            % psf = psf ./ sqrt(mean(psf .^ 2, 'all'));
            if debug
                decon_folder = [datafolder, '/matlab_decon' '/'];
                [~, fsn] = fileparts(input_tiff);
                debug_folder = sprintf('%s/%s_debug/', decon_folder, fsn);
                mkdir(debug_folder)
            else
                debug_folder = '/tmp/'; 
            end                
                
            [deconvolved, err_mat, iter_run] = decon_lucy_function(rawdata, psf, nIter, fixIter, errThresh, debug, debug_folder, useGPU);
            % [deconvolved, err_mat, iter_run] = decon_lucy_function_test(rawdata, psf, nIter, fixIter, errThresh, debug, debug_folder, useGPU);
            % [deconvolved, err_mat, iter_run] = decon_lucy_function_test_1(rawdata, psf, nIter, fixIter, errThresh, debug, debug_folder, useGPU);
            deconvolved = deconvolved * scaleFactor;
        case 'cudagen'
            deconvolved = decon_lucy_cuda_function(single(rawdata), single(psf), nIter) * scaleFactor;            
    end
else
    deconvolved = rawdata;
end

% Rotate the decon result if requested
if ischar(bRotFinal)
    bRotFinal = logical(str2double(bRotFinal));
end

if bRotFinal
    if bReverse
        sign = -1;
    else
        sign = 1;
    end
    deconvolved = rotate3D(deconvolved, sign * rotateByAngle, ...
        dz_data*dz_data_ratio/xypixelsize);
end

if ischar(bFlipZ)
    bFlipZ = logical(str2double(bFlipZ));
end

if bFlipZ
    deconvolved = flipdim(deconvolved, 3);
end

% Z-crop the decon result if requested
if ischar(cropFinal)
    cropFinal = str2num(cropFinal);
end

if ~isempty(cropFinal)
    if cropFinal(2) > cropFinal(1)
        deconvolved = deconvolved(cropFinal(1):cropFinal(2),:,:);
    end
    if cropFinal(4) > cropFinal(3)
        deconvolved = deconvolved(:,cropFinal(3):cropFinal(4),:);
    end
    if cropFinal(6) > cropFinal(5)
        deconvolved = deconvolved(:,:,cropFinal(5):cropFinal(6));
    end
end
%toc

% if the output file is empty, directly return the deconvolved results. 
if isempty(output_filename)
    return;
end

% construct output file name

% decon_filename_tail = '_decon.tif';
decon_folder = ['matlab_decon' '/'];
thumnail_folder = ['downsampled_data' '/'];
MIPs_folder = ['MIPs' '/'];
[~, output_tiff] = fileparts(output_filename);

% output_tiff = strrep(input_tiff, '.tif', decon_filename_tail);
output_tiff = [output_tiff, '.tif'];

if isempty(datafolder)
    output_tiff1 = strcat(decon_folder, output_tiff);
    output_tiff2 = strcat([decon_folder,thumnail_folder], output_tiff);
    MIPs_tiff = strcat([decon_folder,MIPs_folder], inputfile);
else
    output_tiff1 = strcat(datafolder, '/', decon_folder, output_tiff);
    output_tiff2 = strcat(datafolder, '/', [decon_folder, thumnail_folder], output_tiff);
    MIPs_tiff = strcat(datafolder, '/', [decon_folder, MIPs_folder], inputfile);
end

% save err mat for simplified method
if nIter > 0 && strcmp(RLMethod, 'simplified')
    info_folder = sprintf('%s/%s/decon_info/', datafolder, decon_folder);
    if ~exist(info_folder, 'dir')
        mkdir(info_folder);
    end
    info_fn = sprintf('%s/%s_info.mat', info_folder, output_tiff(1 : end - 4));
    save('-v7.3', info_fn, 'err_mat', 'nIter', 'fixIter', 'errThresh', 'debug');
    save(sprintf('%sactual_iterations_%d.txt', info_fn(1 : end - 8), iter_run), 'iter_run', '-ASCII');
end

% not save 3D stack
if ~save3Dstack
    return;
end

if ischar(bSaveUint16)
    bSaveUint16 = logical(str2double(bSaveUint16));
end

if bSaveUint16
    max_val = max(deconvolved(:));
    if max_val > 65535
        % deconvolved = deconvolved * (65535 / max_val);
        warning('%0.2d% voxels are saturated!', mean(deconvolved(:) > 65535));
    end
    deconvolved = uint16(deconvolved);
end

% write3Dtiff(deconvolved, output_tiff1);
uuid = get_uuid();
output_tmp_tiff = [output_tiff1(1 : end - 4), '_', uuid, '.tif'];
writetiff(deconvolved, output_tmp_tiff);
movefile(output_tmp_tiff, output_tiff1);

% generate max-intensity projections if requested:
if ischar(axesMaxIntProj)
    axesMaxIntProj = str2num(axesMaxIntProj);
end

if length(axesMaxIntProj) ~= 3
    warning('Exactly 3 parameters are needed for -M')
else
    if any(axesMaxIntProj)
        mkdir([datafolder,'/',decon_folder,MIPs_folder]);
        axislabel = 'yxz';
        for d=1:3
            if axesMaxIntProj(d)
                % write3Dtiff(max(deconvolved, [], d), strcat(MIPs_tiff, '_MIP_', axislabel(d), '.tif'))
                writetiff(max(deconvolved, [], d), strcat(MIPs_tiff, '_MIP_', axislabel(d), '.tif'))
            end
        end
    end
end

% generate thumbnails in "downsampled_data"
if ischar(resizeFactor)
    resizeFactor = str2double(resizeFactor);
end

if ischar(scalingThresh)
    scalingThresh = str2double(scalingThresh);
end

if resizeFactor ~= 1
    deconvolved = rescale_resample(deconvolved, resizeFactor, scalingThresh);
    mkdir([datafolder, '/' ,decon_folder,thumnail_folder]);
    % write3Dtiff(deconvolved, output_tiff2);
    writetiff(deconvolved, output_tiff2);
end

end
