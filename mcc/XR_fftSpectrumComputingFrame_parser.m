function XR_fftSpectrumComputingFrame_parser(FrameFullpath, SpectrumFullname, varargin)
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
ip.addParameter('xyPixelSize', 0.108, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dz', 0.1, @(x) isnumeric(x) || ischar(x)); % actual pixel size in z
ip.addParameter('outPixelSize', [], @(x) isnumeric(x) || isempty(x) || ischar(x)); % output pixel size
ip.addParameter('N', [1001, 1001, 1001], @(x) isnumeric(x) || ischar(x));
ip.addParameter('save3DStack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('background', 0, @(x) isnumeric(x) || ischar(x));
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

if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(outPixelSize)
    outPixelSize = str2num(outPixelSize);
end
if ischar(N)
    N = str2num(N);
end
if ischar(save3DStack)
    save3DStack = strcmp(save3DStack,'true');
end
if ischar(background)
    background = str2num(background);
end

XR_fftSpectrumComputingFrame(FrameFullpath,SpectrumFullname,'xyPixelSize',xyPixelSize,...
    'dz',dz,'outPixelSize',outPixelSize,'N',N,'save3DStack',save3DStack,...
    'background',background,'Interp',Interp);