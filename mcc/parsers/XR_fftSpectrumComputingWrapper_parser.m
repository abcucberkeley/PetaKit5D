function XR_fftSpectrumComputingWrapper_parser(dataPaths, varargin)


%#function XR_fftSpectrumComputingFrame

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x)); % data structure from loadConditionData
ip.addParameter('resultDirName', 'FFT',  @ischar);
ip.addParameter('overwrite', false,  @(x) (numel(x) == 1) && islogical(x) || ischar(x));
ip.addParameter('xyPixelSize', 0.108, @(x) isnumeric(x) || ischar(x));
ip.addParameter('dz', 0.1, @(x) isnumeric(x) || ischar(x)); % actual pixel size in z
ip.addParameter('zarrFile', false, @(x) islogical(x) || ischar(x)); % input as zarr
ip.addParameter('outPixelSize', [], @(x) isnumeric(x) || isempty(x) || ischar(x)); % output pixel size
ip.addParameter('outSize', [1001, 1001, 1001], @(x) isnumeric(x) || ischar(x));
ip.addParameter('channelPatterns', {}, @(x) ischar(x) || iscell(x));
ip.addParameter('save3DStack', false, @(x) islogical(x) || ischar(x));
ip.addParameter('background', 0, @(x) isnumeric(x) || ischar(x));
ip.addParameter('interpMethod', 'linear', @ischar);
ip.addParameter('parseCluster', true, @(x) islogical(x) || ischar(x));
ip.addParameter('masterCompute', true, @(x) islogical(x) || ischar(x)); % master node participate in the task computing. 
ip.addParameter('cpusPerTask', 3, @(x) isscalar(x) || ischar(x));
ip.addParameter('debug', false, @(x) islogical(x) || ischar(x));
ip.addParameter('mccMode', false, @(x) islogical(x) || ischar(x));
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, varargin{:});

pr = ip.Results;
resultDirName = pr.resultDirName;
overwrite = pr.overwrite;
xyPixelSize = pr.xyPixelSize;
dz = pr.dz;
zarrFile = pr.zarrFile;
outPixelSize = pr.outPixelSize;
outSize = pr.outSize;
channelPatterns = pr.channelPatterns;
save3DStack = pr.save3DStack;
background = pr.background;
interpMethod = pr.interpMethod;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
debug = pr.debug;
mccMode = pr.mccMode;
configFile = pr.configFile;

if ischar(dataPaths) && ~isempty(dataPaths) && strcmp(dataPaths(1), '{')
    dataPaths = eval(dataPaths);
end
if ischar(overwrite)
    overwrite = str2num(overwrite);
end
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
if ischar(channelPatterns) && ~isempty(channelPatterns) && strcmp(channelPatterns(1), '{')
    channelPatterns = eval(channelPatterns);
end
if ischar(save3DStack)
    save3DStack = str2num(save3DStack);
end
if ischar(background)
    background = str2num(background);
end
if ischar(parseCluster)
    parseCluster = str2num(parseCluster);
end
if ischar(masterCompute)
    masterCompute = str2num(masterCompute);
end
if ischar(cpusPerTask)
    cpusPerTask = str2num(cpusPerTask);
end
if ischar(debug)
    debug = str2num(debug);
end
if ischar(mccMode)
    mccMode = str2num(mccMode);
end

XR_fftSpectrumComputingWrapper(dataPaths, resultDirName=resultDirName, overwrite=overwrite, ...
    xyPixelSize=xyPixelSize, dz=dz, zarrFile=zarrFile, outPixelSize=outPixelSize, ...
    outSize=outSize, channelPatterns=channelPatterns, save3DStack=save3DStack, ...
    background=background, interpMethod=interpMethod, parseCluster=parseCluster, ...
    masterCompute=masterCompute, cpusPerTask=cpusPerTask, debug=debug, mccMode=mccMode, ...
    configFile=configFile);

end

