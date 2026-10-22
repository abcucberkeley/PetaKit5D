function [] = verifyLargeFileMIPs(testCase, resultFn)
% verify the MIPs next to a large-file result (a zarr) as XR_MIP_zarr makes
% them: for each axis, the slab zarr MIPs/<name>_{y,x,z}.zarr, whose
% element k along the axis is the max over MIP batch k (MIP_zarr_batch_size)
% of the result, with chunks of a batch and one element along the axis and
% the result's data type, and the MIP tiff MIPs/<name>_MIP_{y,x,z}.tif

[resultPath, fsn] = fileparts(resultFn);
zInfo = getZarrInfo(resultFn);
im = readzarr(resultFn);
batchSize = MIP_zarr_batch_size(zInfo.size, zInfo.blockSize, zInfo.dtype);
axis_strs = {'y', 'x', 'z'};
for j = 1 : 3
    starts = 1 : batchSize(j) : zInfo.size(j);
    ends = [starts(2 : end) - 1, zInfo.size(j)];
    idx = {':', ':', ':'};
    parts = cell(numel(starts), 1);
    for k = 1 : numel(starts)
        idx{j} = starts(k) : ends(k);
        parts{k} = max(im(idx{:}), [], j);
    end
    slab = cat(j, parts{:});

    slabFn = sprintf('%s/MIPs/%s_%s.zarr', resultPath, fsn, axis_strs{j});
    testCase.assertTrue(exist(slabFn, 'dir') > 0, sprintf('No %s slab zarr.', axis_strs{j}));
    sInfo = getZarrInfo(slabFn);
    slabSize = zInfo.size;
    slabSize(j) = numel(starts);
    chunks = batchSize;
    chunks(j) = 1;
    testCase.verifyEqual(sInfo.size, slabSize, sprintf('Size of the %s slab.', axis_strs{j}));
    testCase.verifyEqual(sInfo.blockSize, min(chunks, slabSize), sprintf('Chunks of the %s slab.', axis_strs{j}));
    testCase.verifyEqual(sInfo.dtype, zInfo.dtype, sprintf('Data type of the %s slab.', axis_strs{j}));
    testCase.verifyEqual(readzarr(slabFn), slab, sprintf('The %s slab.', axis_strs{j}));

    mip = readtiff(sprintf('%s/MIPs/%s_MIP_%s.tif', resultPath, fsn, axis_strs{j}));
    testCase.verifyEqual(mip, squeeze(max(im, [], j)), sprintf('The %s MIP.', axis_strs{j}));
end

end
