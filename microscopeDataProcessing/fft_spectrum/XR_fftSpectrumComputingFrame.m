function XR_fftSpectrumComputingFrame(FrameFullpath, SpectrumFullname, varargin)
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
ip.addRequired('FrameFullpath');
ip.addRequired('SpectrumFullname');
ip.addParameter('xyPixelSize', 0.108, @isnumeric);
ip.addParameter('dz', 0.1, @isnumeric); % actual pixel size in z
ip.addParameter('outPixelSize', [], @(x) isnumeric(x) || isempty(x)); % output pixel size
ip.addParameter('N', [1001, 1001, 1001], @isnumeric);
ip.addParameter('save3DStack', false, @islogical);
ip.addParameter('background', 0, @isnumeric);
ip.addParameter('Interp', 'linear', @ischar);
ip.parse(FrameFullpath, SpectrumFullname, varargin{:});

pr = ip.Results;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
outPixelSize = pr.outPixelSize;
N = pr.N;
save3DStack = pr.save3DStack;
background = pr.background;
Interp = pr.Interp;

uuid = get_uuid();

if exist(SpectrumFullname, 'file')
    return;
end

sz = getImageSize(FrameFullpath);

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

try
    im = parallelReadTiff(FrameFullpath, [s(3), t(3)]);
catch
    im = readtiff(FrameFullpath, 'range', s(3) : t(3));
end
im = im(s(1) : t(1), s(2) : t(2), :);

if background ~= 0
    im = double(im) - background;
    im = im .* (im > 0);
end

im = imresize3(double(im), round(size(im) .* [xyPixelSize, xyPixelSize, dz] ./ px), 'method', Interp);
sz = size(im);

hfN = (N - 1) / 2;

if any(size(im) ~= N)
    if any(sz < N)
        im = padarray(im, max(0, floor((N - sz) / 2)), 0, 'pre');
        im = padarray(im, max(0, ceil((N - sz) / 2)), 0, 'post');
    end
    
    c = round((size(im) + 1) / 2);
    im = im(c(1) - hfN(1) : c(1) + hfN(1), c(2) - hfN(2) : c(2) + hfN(2), c(3) - hfN(3) : c(3) + hfN(3));
end

im = abs(fftshift(fftn(im)));
im = im ./ max(im(:));

[FFTPath, fsname] = fileparts(SpectrumFullname);

spectrumTmpname = [FFTPath, filesep, fsname, '_', uuid, '.tif'];
if save3DStack
    writetiff(log10(single(im)), spectrumTmpname);
else
    fclose(fopen(spectrumTmpname, 'w'));
end

FFTMIPPath = [FFTPath, '/MIPs/'];
% mkdir(rotMIPPath);

% apply gamma
im = single(im .^ 0.5);
MIP = max(im, [], 1);
FFTMIPFullpath = [FFTMIPPath, fsname, '_MIP_y.tif'];
writetiff(squeeze(MIP), FFTMIPFullpath);    

MIP = max(im, [], 2);
FFTMIPFullpath = [FFTMIPPath, fsname, '_MIP_x.tif'];
writetiff(squeeze(MIP), FFTMIPFullpath);

MIP = max(im, [], 3);
FFTMIPFullpath = [FFTMIPPath, fsname, '_MIP_z.tif'];
writetiff(squeeze(MIP), FFTMIPFullpath);


% write central slices
FFTCSPath = [FFTPath, '/central_slices/'];
if ~exist(FFTCSPath, 'dir')
    mkdir(FFTCSPath);
end

[ny, nx, nz] = size(im);
MIP = im(:, :, round((nz + 1) / 2));
FFTMIPFullpath = [FFTCSPath, fsname, '_center_xy.tif'];
writetiff(squeeze(MIP), FFTMIPFullpath);    

MIP = im(:, round((nx + 1) / 2), :);
FFTMIPFullpath = [FFTCSPath, fsname, '_center_yz.tif'];
writetiff(squeeze(MIP), FFTMIPFullpath);

MIP = im(round((ny + 1) / 2), :, :);
FFTMIPFullpath = [FFTCSPath, fsname, '_center_xz.tif'];
writetiff(squeeze(MIP), FFTMIPFullpath);

% move to final filename
movefile(spectrumTmpname, SpectrumFullname);


end


