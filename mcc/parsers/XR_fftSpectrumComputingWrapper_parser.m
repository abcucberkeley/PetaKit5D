function XR_fftSpectrumComputingWrapper_parser(dataPaths, varargin)


%#function XR_fftSpectrumComputingFrame

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x)); % data structure from loadConditionData
ip.addParameter('Overwrite', false,  @(x) (numel(x) == 1 || numel(x) == 5) && islogical(x) || ischar(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dz', 0.1, @(x) isnumeric(x) || ischar(x)); % actual pixel size in z
ip.addParameter('outPixelSize', [], @(x) isnumeric(x) || isempty(x) || ischar(x)); % output pixel size
ip.addParameter('N', [1001, 1001, 1001], @(x) isnumeric(x) || ischar(x));
ip.addParameter('ChannelPatterns', {}, @(x) ischar(x) || iscell(x));
ip.addParameter('save3DStack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('background', 0, @(x) isnumeric(x) || ischar(x));
ip.addParameter('Interp', 'linear', @ischar);
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('ConfigFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
Overwrite = pr.Overwrite;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
outPixelSize = pr.outPixelSize;
N = pr.N;
ChannelPatterns = pr.ChannelPatterns;
save3DStack = pr.save3DStack;
background = pr.background;
Interp = pr.Interp;
mccMode = pr.mccMode;
ConfigFile = pr.ConfigFile;

if ischar(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(Overwrite)
    Overwrite = str2num(Overwrite);
end
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
if ischar(ChannelPatterns) && strcmp(ChannelPatterns(1), '{')
    ChannelPatterns = eval(ChannelPatterns);
end
if ischar(save3DStack)
    save3DStack = str2num(save3DStack);
end
if ischar(background)
    background = str2num(background);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_fftSpectrumComputingWrapper(dataPaths, Overwrite=Overwrite, xyPixelSize=xyPixelSize, ...
    dz=dz, outPixelSize=outPixelSize, N=N, ChannelPatterns=ChannelPatterns, ...
    save3DStack=save3DStack, background=background, Interp=Interp, mccMode=mccMode, ...
    ConfigFile=ConfigFile);

end

