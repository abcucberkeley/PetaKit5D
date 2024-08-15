function [] = processStitchSeparteTiles(zarrFullpaths, stitchPath, stIndices, pImSz, varargin)
% process the tiles in the corresponding location in the stitched image as
% separte images. 
% 
% 
% Author: Xiongtao Ruan (11/19/2020)
%
% (07/08/2024): use crop function to pad the tiles to the size of the stitched image.


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpaths', @iscell);
ip.addRequired('stitchPath', @(x) ischar(x));
ip.addRequired('stIndices', @isnumeric);
ip.addRequired('pImSz', @isnumeric);
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('saveMIP', true, @islogical);
ip.addParameter('batchSize', [1024, 1024, 1024] , @isvector);
ip.addParameter('blockSize', [256, 256, 256] , @isnumeric);
ip.addParameter('parseCluster', true, @islogical);
ip.addParameter('masterCompute', true, @islogical);
ip.addParameter('cpusPerTask', 1, @isnumeric);
ip.addParameter('uuid', '', @ischar);
ip.addParameter('debug', false, @islogical);
ip.addParameter('mccMode', false, @islogical);
ip.addParameter('configFile', '', @ischar);

ip.parse(zarrFullpaths, stitchPath, stIndices, pImSz, varargin{:});

pr = ip.Results;
saveMIP = pr.saveMIP;
batchSize = pr.batchSize;
blockSize = pr.blockSize;
parseCluster = pr.parseCluster;
masterCompute = pr.masterCompute;
cpusPerTask = pr.cpusPerTask;
uuid = pr.uuid;
mccMode = pr.mccMode;
configFile = pr.configFile;

if isempty(uuid)
    uuid = get_uuid();
end

nF = numel(zarrFullpaths);
outputFullpaths = cell(nF, 1);
funcStrs = cell(nF, 1);

for f = 1 : nF
    zarrFullpath = zarrFullpaths{f};
    [~, fsname, ext] = fileparts(zarrFullpath);
    outputFullpath = [stitchPath, '/', fsname, ext];
    outputFullpaths{f} = outputFullpath;
    
    st_idx = stIndices(f, :);
    bbox = [2 - st_idx([2, 1, 3]), 2 - st_idx([2, 1, 3]) + pImSz - 1];
    pad = true;
    parseParfor = false;

    funcStrs{f} = sprintf(['XR_crop_zarr(''%s'',''%s'',%s,''pad'',%s,''batchSize'',%s,', ...
        '''blockSize'',%s,''saveMIP'',%s,''parseCluster'',%s,''parseParfor'',%s,', ...
        '''masterCompute'',%s,''cpusPerTask'',%d,''uuid'',''%s'',''mccMode'',%s,''configFile'',''%s'')'], ...
        zarrFullpath, outputFullpath, mat2str_comma(bbox), string(pad), mat2str_comma(batchSize), ...
        mat2str_comma(blockSize), string(saveMIP), string(parseCluster), string(parseParfor), ...
        string(masterCompute), cpusPerTask, uuid, string(mccMode), configFile);
end

dtype = getImageDataType(zarrFullpaths{1});
byte_num = dataTypeToByteNumber(dtype);

% cluster setting
memAllocate = prod(batchSize) * byte_num / 1024^3 * 10;
maxTrialNum = 2;
is_done_flag = false;

% retry with more resources and longer time
for i = 1 : 3
    if ~all(is_done_flag)
        is_done_flag = generic_computing_frameworks_wrapper(zarrFullpaths, outputFullpaths, ...
            funcStrs, 'cpusPerTask', cpusPerTask * 2^(i-1), 'memAllocate', memAllocate * 2^(i-1), ...
            'maxTrialNum', maxTrialNum, 'masterCompute', masterCompute, 'parseCluster', parseCluster, ...
            'mccMode', mccMode, 'configFile', configFile);
    end
end

end

