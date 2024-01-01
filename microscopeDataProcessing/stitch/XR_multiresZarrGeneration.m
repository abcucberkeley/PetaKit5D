function [] = XR_multiresZarrGeneration(zarrFullpath, outputFullpath, varargin)
% create a multiresolution zarr file from a given zarr file
% 
% Author: Xiongtao Ruan (12/14/2020)
% xruan (11/08/2021): add support for multires within the same file and
% refactor the code


if nargin < 1
    zarrFullpath = '/clusterfs/fiona/Data/BetzigLab/Janelia_MOSAIC/20211022_ISM_AO/Fish3_XYTile_Processed/Linescan_Processed_AO/_0000t/matlab_stitch_xcorr_feather_zarr/Scan_Iter_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0113378311msecAbs_local_assignment_test.zarr';
    outputFullpath = [zarrFullpath, '_multires'];
    % outputFullpath = [zarrFullpath];
end

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('zarrFullpath'); 
ip.addRequired('outputFullpath'); 
ip.addParameter('dsFactors', [2, 2; 4, 4; 8, 8; 16, 16; 32, 32], @isnumeric); % downsample factors
ip.addParameter('blockSize', [256, 256, 256], @isnumeric); % blockSize
ip.addParameter('batchSize', [512, 512, 512], @isnumeric); % batchSize
ip.addParameter('Interp', 'linear', @(x) any(strcmpi(x, {'cubic', 'linear', 'nearest'})));
ip.addParameter('includeRaw', ~false, @islogical); % include raw data
% ip.addParameter('suffix', '', @ischar); % suffix for the folder
ip.addParameter('uuid', '', @ischar);

ip.parse(zarrFullpath, outputFullpath, varargin{:});


pr = ip.Results;
dsFactors = pr.dsFactors;
blockSize = pr.blockSize;
batchSize = pr.batchSize;
Interp = pr.Interp;
includeRaw = pr.includeRaw;
uuid = pr.uuid;

if isempty(uuid)
    uuid = get_uuid();
end

% check if the group folder exist

resLevel = size(dsFactors, 1); 
dataFullpaths = cell(resLevel + includeRaw, 1);
for i = 1 : resLevel + includeRaw
    if i == 1 && includeRaw
        dataFullpaths{i} = sprintf('%s/L_1_1_1', outputFullpath);
    else
        dataFullpaths{i} = sprintf('%s/L_%d_%d_%d', outputFullpath, dsFactors(i - includeRaw, 1), dsFactors(i - includeRaw, 1), dsFactors(i - includeRaw, 2));
    end
end

if exist(outputFullpath, 'dir')
    exist_flags = false(numel(dataFullpaths), 1);
    for i = 1 : numel(dataFullpaths)
        exist_flags(i) = exist(dataFullpaths{i}, 'dir');
    end
    if all(exist_flags)
        disp('The output result exists, skip it!');
        return;
    end
end

% create the group file
outputTmppath = sprintf('%s_%s', outputFullpath, uuid);
py.zarr.open_group(outputTmppath, pyargs('mode', 'w'));

dsTmpFullpath_im1 = zarrFullpath;
for i = 1 : size(dsFactors, 1)
    if i == 1
        dsFactor_i = dsFactors(i, :);
    else
        dsFactor_i = dsFactors(i, :) ./ dsFactors(i - 1, :);
    end
    blockSize_i = round(blockSize ./ [dsFactors(i, 1), dsFactors(i, 1), dsFactors(i, 2)]);
    batchSize_i = round(batchSize ./ [dsFactors(i, 1), dsFactors(i, 1), dsFactors(i, 2)]);
    
    dsFullpath = sprintf('%s/L_%d_%d_%d', outputFullpath, dsFactors(i, 1), dsFactors(i, 1), dsFactors(i, 2));
    dsTmpFullpath = sprintf('%s/L_%d_%d_%d', outputTmppath, dsFactors(i, 1), dsFactors(i, 1), dsFactors(i, 2));
    if exist(dsFullpath, 'dir')
        dsTmpFullpath_im1 = dsFullpath;
        continue;
    end        
        
    XR_resampleSingleZarr(dsTmpFullpath_im1, dsTmpFullpath, dsFactor_i, 'blockSize', blockSize_i, ...
        'batchSize', batchSize_i, 'Interp', Interp, 'uuid', uuid);
    
    if exist(dsTmpFullpath, 'dir')
        movefile(dsTmpFullpath, dsFullpath);
    end
    
    if exist(dsFullpath, 'dir')
        dsTmpFullpath_im1 = dsFullpath;
    end
end

% if the raw data is directly inside the group folder, move it first
if includeRaw
    dsTmpFullpath = sprintf('%s/L_1_1_1', outputFullpath);    
    movefile(zarrFullpath, dsTmpFullpath);

    system(sprintf('mv %s/{.,}* %s', zarrFullpath, dataFullpaths{1}));
end

movefile(outputTmppath, outputFullpath);

end


