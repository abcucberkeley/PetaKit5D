function [] = XR_resample_dataset(dataPaths, resampleFactor, varargin)
% resample a dataset via cluster computing
% resampleFactor [5,5,5] (yxz) means downsample data by 5x5x5
% 
%
% Author: Xiongtao Ruan (01/15/2021)
% 
% xruan (03/09/2022): add support for zarr
% xruan (04/17/2024): add support for large zarr


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPaths', @(x) ischar(x) || iscell(x));
ip.addRequired('resampleFactor', @isnumeric);
ip.addParameter('resultDirName', 'resampled', @ischar);
ip.addParameter('channelPatterns', {'CamA_ch0', 'CamA_ch1', 'CamB_ch0', 'CamB_ch1'}, @iscell);
ip.addParameter('inputBbox', [], @isnumeric); % bbox for input
ip.addParameter('interpMethod', 'linear', @(x) ischar(x) && any(strcmpi(x, {'cubic', 'linear', 'nearest'})));
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('zarrFile', false, @islogical);
ip.addParameter('largeZarr', false, @islogical);
ip.addParameter('saveZarr', false, @islogical); % use zarr file as output
ip.addParameter('blockSize', [256, 256, 256], @isnumeric); % blcoksize
ip.addParameter('batchSize', [512, 512, 512], @isnumeric); % size to process in one batch
ip.addParameter('borderSize', [5, 5, 5], @isnumeric); % padded boarder for each batch
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('jobLogDir', '../job_logs', @ischar);
ip.addParameter('masterCompute', true, @islogical); % master node participate in the task computing. 
ip.addParameter('cpusPerTask', 2, @isscalar);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(dataPaths, resampleFactor, varargin{:});

warning('off', 'MATLAB:MKDIR:DirectoryExists');

pr = ip.Results;
resultDirName = pr.resultDirName;
channelPatterns = pr.channelPatterns;
inputBbox = pr.inputBbox;
interpMethod = pr.interpMethod;
save16bit = pr.save16bit;
zarrFile = pr.zarrFile;
largeZarr = pr.largeZarr;
saveZarr = pr.saveZarr;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
borderSize = pr.borderSize;
parseCluster = pr.parseCluster;
jobLogDir = pr.jobLogDir;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

if isempty(uuid)
    uuid = get_uuid();
end

if ischar(dataPaths)
    dataPaths = {dataPaths};
end

nd = numel(dataPaths);
resultPaths = cell(nd, 1);
for d = 1 : nd
    dataPaths{d} = simplifyPath(dataPaths{d});
    resultPath = [dataPaths{d}, '/' resultDirName, '/'];
    resultPaths{d} = resultPath;
    mkdir(resultPath);
    fileattrib(resultPath, '+w', 'g');

    save('-v7.3', [resultPath, '/parameters.mat'], 'pr');
    s = jsonencode(pr, PrettyPrint=true);
    fid = fopen([resultPath, '/parameters.json'], 'w');
    fprintf(fid, s);
    fclose(fid);
end

% parse image filenames
[fnames, fsns, fd_inds, filepaths] = parseImageFilenames(dataPaths, zarrFile, channelPatterns);
nF = numel(fnames);
sz = getImageSize(filepaths{1});
dtype = getImageDataType(filepaths{1});
byteNum = dataTypeToByteNumber(dtype);

%% use generic framework for the cropping do computing

frameFullpaths = filepaths;

ext = '.tif';
if saveZarr
    ext = '.zarr';
end
resultFullpaths = arrayfun(@(x) sprintf('%s/%s%s', resultPaths{fd_inds(x)}, fsns{x}, ext), ...
    1 : nF, 'UniformOutput', false);

if largeZarr
    func_strs = arrayfun(@(x) sprintf(['XR_resampleSingleZarr(''%s'',''%s'',%s,', ...
        '''inputBbox'',%s,''blockSize'',%s,''batchSize'',%s,''borderSize'',%s,', ...
        '''interpMethod'',''%s'',''parseCluster'',%s,''cpusPerTask'',%d,''uuid'',''%s'',', ...
        '''mccMode'',%s,''configFile'',''%s'')'], frameFullpaths{x}, ...
        resultFullpaths{x}, mat2str_comma(resampleFactor, 10), mat2str_comma(inputBbox), ...
        mat2str_comma(blockSize), mat2str_comma(batchSize), mat2str_comma(borderSize), ...
        interpMethod, string(parseCluster),cpusPerTask, uuid, string(mccMode), string(configFile)), ...
        1 : nF, 'unif', false);

    memAllocate = prod(batchSize) * byteNum / 2^30 * (prod(resampleFactor) + 1) * 3;
else
    func_strs = arrayfun(@(x) sprintf(['XR_resampleFrame(''%s'',''%s'',[%s],', ...
        '''inputBbox'',%s,''zarrFile'',%s,''saveZarr'',%s,''interpMethod'',''%s'',''save16bit'',%s,', ...
        '''uuid'',''%s'')'], frameFullpaths{x}, resultFullpaths{x}, mat2str_comma(resampleFactor, 10), ...
        mat2str_comma(inputBbox), string(zarrFile), string(saveZarr), interpMethod, string(save16bit), ...
        uuid), 1 : nF, 'unif', false);

    memAllocate = prod(sz) * byteNum / 1024^3 * (2 + 2 / prod(resampleFactor));
end

generic_computing_frameworks_wrapper(frameFullpaths, resultFullpaths, func_strs, ...
    parseCluster=parseCluster, jobLogDir=jobLogDir, masterCompute=masterCompute, ...
    cpusPerTask=cpusPerTask, memAllocate=memAllocate, mccMode=mccMode, configFile=configFile);

end
