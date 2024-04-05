function XR_fftSpectrumComputingFrame_parser(FrameFullpath, SpectrumFullname, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('FrameFullpath', @ischar);
ip.addRequired('SpectrumFullname', @ischar);
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
    save3DStack = str2num(save3DStack);
end
if ischar(background)
    background = str2num(background);
end

XR_fftSpectrumComputingFrame(FrameFullpath, SpectrumFullname, xyPixelSize=xyPixelSize, ...
    dz=dz, outPixelSize=outPixelSize, N=N, save3DStack=save3DStack, background=background, ...
    Interp=Interp);

end

