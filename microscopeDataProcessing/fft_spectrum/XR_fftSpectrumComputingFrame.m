function XR_fftSpectrumComputingFrame(frameFullpath, spectrumFullname, varargin)
% compute fft spectrum in log scale for a given image
% 
% Author: Xiongtao Ruan (11/25/2020)
%
% xruan (07/09/2021): change fft image to original scale (instead of log scale), add writing of central slices. 
% xruan (07/20/2021): change to normalize spectrum to [0, 1] and apply
% gamma (0.5) to central slices and mips
% xruan (11/11/2021): add support for rescale to isotropic and crop/pad to
% given size (1001 in each dimension)
% xruan (11/12/2021): for big data, change to first predefine the region to load
% xruan (12/13/2021): add background subtraction
% xruan (01/03/2022): add support for output voxel size
% xruan (07/18/2022): add support for different interpolation methods for resampling


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpath', @ischar);
ip.addRequired('spectrumFullname', @ischar);
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.1, @isnumeric); % actual pixel size in z
ip.addParameter('zarrFile', false, @islogical); % input as zarr
ip.addParameter('outPixelSize', [], @(x) isnumeric(x) || isempty(x)); % output pixel size
ip.addParameter('outSize', [1001, 1001, 1001], @isnumeric);
ip.addParameter('save3DStack', false, @islogical);
ip.addParameter('background', 0, @isnumeric);
ip.addParameter('interpMethod', 'linear', @ischar);
ip.parse(frameFullpath, spectrumFullname, varargin{:});

pr = ip.Results;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
zarrFile = pr.zarrFile;
outPixelSize = pr.outPixelSize;
N = pr.outSize;
save3DStack = pr.save3DStack;
background = pr.background;
interpMethod = pr.interpMethod;

uuid = get_uuid();

if exist(spectrumFullname, 'file')
    return;
end

sz = getImageSize(frameFullpath);

% resize image to make it isotropic
if isempty(outPixelSize)
    px = min(xyPixelSize, dz);
else
    px = outPixelSize;
end
% first crop the image if it is too larger (add 10% buffer)
keep_size = min(sz, ceil(N ./ ([xyPixelSize, xyPixelSize, dz] ./ px)) .* 1.1);
s = max(1, ceil((sz - keep_size) / 2));
t = min(sz, s + keep_size - 1);

if zarrFile
    im = readzarr(frameFullpath, inputBbox=[s, t]);
else
    im = readtiff(frameFullpath, 'range', [s(3), t(3)]);
    im = im(s(1) : t(1), s(2) : t(2), :);
end

if background ~= 0
    im = single(im) - background;
    im = im .* (im > 0);
end

% resize to the output pixel size, and pad the dimensions smaller than N and
% crop the center N voxels: only that region is resized (from the input
% region it depends on)
rsz = round(size(im) .* [xyPixelSize, xyPixelSize, dz] ./ px);

hfN = (N - 1) / 2;

if any(rsz ~= N)
    if numel(N) ~= 3 || any(mod(N, 2) == 0)
        error('The output size %s must be three odd numbers to crop or pad the resized image (size %s) to it.', ...
            mat2str(N), mat2str(rsz));
    end
    % the center N voxels, or all voxels of a dimension smaller than N
    c = round((rsz + 1) / 2);
    s = max(1, c - hfN);
    t = min(rsz, c + hfN);
    im = imresize3_fast(single(im), rsz, interpMethod, [s, t]);
    if any(rsz < N)
        im = padarray(im, max(0, floor((N - rsz) / 2)), 0, 'pre');
        im = padarray(im, max(0, ceil((N - rsz) / 2)), 0, 'post');
    end
else
    im = imresize3_fast(single(im), rsz, interpMethod);
end

% |FFT| without fftshift: the MIPs and central slices (the DC planes) are
% taken from it and then shifted, and the normalization and gamma (monotonic)
% are applied to them, which gives the same values as on the whole spectrum
im = abs(fftn(im));
imMax = max(im(:));

[FFTPath, fsname] = fileparts(spectrumFullname);

spectrumTmpname = [FFTPath, filesep, fsname, '_', uuid, '.tif'];
if save3DStack
    writetiff(log10(single(fftshift(im) ./ imMax)), spectrumTmpname);
else
    fclose(fopen(spectrumTmpname, 'w'));
end

FFTMIPPath = [FFTPath, '/MIPs/'];
% mkdir(rotMIPPath);

% apply gamma
MIP = normalize_gamma(max(im, [], 1), imMax);
FFTMIPFullpath = [FFTMIPPath, fsname, '_MIP_y.tif'];
writetiff(squeeze(MIP), FFTMIPFullpath);    

MIP = normalize_gamma(max(im, [], 2), imMax);
FFTMIPFullpath = [FFTMIPPath, fsname, '_MIP_x.tif'];
writetiff(squeeze(MIP), FFTMIPFullpath);

MIP = normalize_gamma(max(im, [], 3), imMax);
FFTMIPFullpath = [FFTMIPPath, fsname, '_MIP_z.tif'];
writetiff(squeeze(MIP), FFTMIPFullpath);


% write central slices
FFTCSPath = [FFTPath, '/central_slices/'];
if ~exist(FFTCSPath, 'dir')
    mkdir(FFTCSPath);
end

MIP = normalize_gamma(im(:, :, 1), imMax);
FFTMIPFullpath = [FFTCSPath, fsname, '_center_xy.tif'];
writetiff(squeeze(MIP), FFTMIPFullpath);    

MIP = normalize_gamma(im(:, 1, :), imMax);
FFTMIPFullpath = [FFTCSPath, fsname, '_center_yz.tif'];
writetiff(squeeze(MIP), FFTMIPFullpath);

MIP = normalize_gamma(im(1, :, :), imMax);
FFTMIPFullpath = [FFTCSPath, fsname, '_center_xz.tif'];
writetiff(squeeze(MIP), FFTMIPFullpath);

% move to final filename
movefile(spectrumTmpname, spectrumFullname);


end


function [im] = normalize_gamma(im, imMax)
% normalize to the maximum of the spectrum, apply gamma 0.5, and fftshift

im = fftshift(single((im ./ imMax) .^ 0.5));

end
