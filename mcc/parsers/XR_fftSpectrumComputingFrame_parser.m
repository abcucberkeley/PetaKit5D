function XR_fftSpectrumComputingFrame_parser(frameFullpath, spectrumFullname, varargin)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('frameFullpath', @ischar);
ip.addRequired('spectrumFullname', @ischar);
ip.addParameter('xyPixelSize', 0.108, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dz', 0.1, @(x) isnumeric(x) || ischar(x)); % actual pixel size in z
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x)); % input as zarr
ip.addParameter('outPixelSize', [], @(x) isnumeric(x) || isempty(x) || ischar(x)); % output pixel size
ip.addParameter('outSize', [1001, 1001, 1001], @(x) isnumeric(x) || ischar(x));
ip.addParameter('save3DStack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('background', 0, @(x) isnumeric(x) || ischar(x));
ip.addParameter('interpMethod', 'linear', @ischar);

ip.parse(frameFullpath, spectrumFullname, varargin{:});

pr = ip.Results;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
zarrFile = pr.zarrFile;
outPixelSize = pr.outPixelSize;
outSize = pr.outSize;
save3DStack = pr.save3DStack;
background = pr.background;
interpMethod = pr.interpMethod;

if ischar(xyPixelSize)
    xyPixelSize = str2num(xyPixelSize);
end
if ischar(dz)
    dz = str2num(dz);
end
if ischar(zarrFile)
    zarrFile = str2num(zarrFile);
end
if ischar(outPixelSize)
    outPixelSize = str2num(outPixelSize);
end
if ischar(outSize)
    outSize = str2num(outSize);
end
if ischar(save3DStack)
    save3DStack = str2num(save3DStack);
end
if ischar(background)
    background = str2num(background);
end

XR_fftSpectrumComputingFrame(frameFullpath, spectrumFullname, xyPixelSize=xyPixelSize, ...
    dz=dz, zarrFile=zarrFile, outPixelSize=outPixelSize, outSize=outSize, save3DStack=save3DStack, ...
    background=background, interpMethod=interpMethod);

end

