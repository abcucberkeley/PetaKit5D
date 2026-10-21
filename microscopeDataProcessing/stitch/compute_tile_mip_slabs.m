function [mipDirStr] = compute_tile_mip_slabs(zarrFullpaths, poolSize, parseCluster, mccMode, configFile)
% compute the max pooled MIP slabs of large tiles, used by the registration
% (xcorrShift) and the feather blending distance maps of large file
% stitching; existing slabs are skipped
%
% moved from stitch_shift_assignment.m


nF = numel(zarrFullpaths);
axis = [1, 1, 1];
batchSize = [1024, 1024, 1024];
% poolSize = [5, 5, 5];
if numel(poolSize) == 3
    mipDirStr = sprintf('MIP_Slabs_pooling_%d_%d_%d', poolSize(1), poolSize(2), poolSize(3));
else
    mipDirStr = sprintf('MIP_Slabs_pooling_%d_%d_%d_%d_%d_%d', poolSize(1), ...
        poolSize(2), poolSize(3), poolSize(4), poolSize(5), poolSize(6));
end
mipSlab = true;
parseParfor = false;

[dataPaths, fsns] = fileparts(zarrFullpaths);
uniq_dataPaths = unique(dataPaths);
for d = 1 : numel(uniq_dataPaths)
    mipPath = [uniq_dataPaths{d}, '/', mipDirStr];
    if ~exist(mipPath, 'dir')
        mkdir(mipPath);
    end
end

funcStrs = cell(nF, 1);
outputFullpaths = cell(nF, 1);
for f = 1 : nF
    funcStrs{f} = sprintf(['XR_MIP_zarr(''%s'',''resultDirName'',''%s'',''axis'',%s,', ...
        '''batchSize'',%s,''poolSize'',%s,''mipSlab'',%s,''parseCluster'',%s,', ...
        '''parseParfor'',%s,''mccMode'',%s,''configFile'',''%s'')'], ...
        zarrFullpaths{f}, mipDirStr, strrep(mat2str(axis), ' ', ','), strrep(mat2str(batchSize), ' ', ','), ...
        strrep(mat2str(poolSize), ' ', ','), string(mipSlab), string(parseCluster), ...
        string(parseParfor), string(mccMode), configFile);
    outputFullpaths{f} = sprintf('%s/%s/%s_z.zarr', dataPaths{f}, mipDirStr, fsns{f});
end

inputFullpaths = zarrFullpaths;
memAllocate = prod(batchSize) * 4 / 2^30 * 3;
minTaskJobNum = round(nF / 2);
for i = 1 : 3
    is_done_flag = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, ...
        funcStrs, 'cpusPerTask', 1, 'maxTrialNum', 2, 'parseCluster', parseCluster, ...
        'memAllocate', memAllocate * 2^(i - 1), 'minTaskJobNum', minTaskJobNum, ...
        'mccMode', mccMode, 'configFile', configFile);
    if all(is_done_flag)
        break;
    end
end

end
