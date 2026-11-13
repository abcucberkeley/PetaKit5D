classdef MIPTest < pktest.TestCase
    % maximum intensity projections: XR_MIP_wrapper (tiff, zarr and large
    % zarr paths), saveMIP_tiff, and max_pooling_3d; and min_bbox_3d

    properties (Constant)
        ImSize = [40, 36, 12]
        Channels = {'CamA_ch0', 'CamB_ch0'}
    end

    properties (TestParameter)
        poolSize = struct('cube', [2, 2, 2], 'anisotropic', [3, 2, 1], 'nonDivisible', [4, 5, 5])
        poolType = {'uint16', 'single'}
        zarrType = {'uint16', 'single'}
        % min_bbox_3d_mex takes these but int16
        minType = {'logical', 'uint8', 'uint16', 'single', 'double', 'int16'}
    end

    methods (Test, TestTags = {'Integration'})
        function allAxesOfTiff(testCase)
            [dataPath, ims, names] = testCase.writeDataset('tif');
            XR_MIP_wrapper(dataPath, 'axis', [1, 1, 1], 'channelPatterns', testCase.Channels, ...
                'parseCluster', false, 'masterCompute', true);
            for i = 1 : 2
                testCase.verifyMIPs([dataPath, 'MIPs/', names{i}], ims{i}, 'uint16');
            end
            testCase.verifyFalse(exist([dataPath, 'MIPs/', names{3}, '_MIP_z.tif'], 'file') > 0, ...
                'A file that matches no channel pattern was projected.');
        end

        function defaultIsZOnly(testCase)
            [dataPath, ims, names] = testCase.writeDataset('tif');
            XR_MIP_wrapper(dataPath, 'channelPatterns', testCase.Channels(1), 'parseCluster', false, 'masterCompute', true);
            base = [dataPath, 'MIPs/', names{1}];
            testCase.verifyEqual(readtiff([base, '_MIP_z.tif']), max(ims{1}, [], 3));
            testCase.verifyFalse(exist([base, '_MIP_x.tif'], 'file') > 0);
            testCase.verifyFalse(exist([base, '_MIP_y.tif'], 'file') > 0);
        end

        function singleOutputWithoutSave16bit(testCase)
            [dataPath, ims, names] = testCase.writeDataset('tif', 'single');
            XR_MIP_wrapper(dataPath, 'axis', [1, 1, 1], 'channelPatterns', testCase.Channels(1), ...
                'save16bit', false, 'parseCluster', false, 'masterCompute', true);
            testCase.verifyMIPs([dataPath, 'MIPs/', names{1}], ims{1}, 'single');
        end

        function allAxesOfZarr(testCase)
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            XR_MIP_wrapper(dataPath, 'axis', [1, 1, 1], 'channelPatterns', testCase.Channels, 'zarrFile', true, ...
                'parseCluster', false, 'masterCompute', true);
            for i = 1 : 2
                testCase.verifyMIPs([dataPath, 'MIPs/', names{i}], ims{i}, 'uint16');
            end
        end

        function allAxesOfLargeZarr(testCase)
            % batch-wise MIPs of a zarr file, merged across batches
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            XR_MIP_wrapper(dataPath, 'axis', [1, 1, 1], 'channelPatterns', testCase.Channels(1), 'zarrFile', true, ...
                'largeFile', true, 'batchSize', [16, 16, 4], 'parseCluster', false, 'masterCompute', true);
            testCase.verifyMIPs([dataPath, 'MIPs/', names{1}], ims{1}, 'uint16');
        end

        function inputBboxProjectsTheCroppedRegion(testCase)
            [dataPath, ims, names] = testCase.writeDataset('tif');
            b = [5, 7, 3, 30, 33, 10];
            XR_MIP_wrapper(dataPath, 'axis', [1, 1, 1], 'inputBbox', b, 'channelPatterns', testCase.Channels(1), ...
                'parseCluster', false, 'masterCompute', true);
            testCase.verifyMIPs([dataPath, 'MIPs/', names{1}], ims{1}(b(1) : b(4), b(2) : b(5), b(3) : b(6)), 'uint16');
        end

        function largeZarrInputBboxThatEndsBeforeTheImageEnd(testCase)
            % (XR_MIP_zarr) the last batches of a box whose size isn't a
            % multiple of the batch size ran to the end of the image instead
            % of the box, so their MIPs went past the end of the MIP slabs
            % and writezarr raised "Upper bound is invalid"
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            b = [5, 7, 3, 30, 33, 10];
            XR_MIP_wrapper(dataPath, 'axis', [1, 1, 1], 'inputBbox', b, 'channelPatterns', testCase.Channels(1), ...
                'zarrFile', true, 'largeFile', true, 'batchSize', [16, 16, 4], 'parseCluster', false, 'masterCompute', true);
            testCase.verifyMIPs([dataPath, 'MIPs/', names{1}], ims{1}(b(1) : b(4), b(2) : b(5), b(3) : b(6)), 'uint16');
        end
    end

    methods (Test, TestTags = {'Unit'})
        function saveMIPTiffWritesRequestedAxes(testCase)
            dataPath = testCase.createTempFolder();
            im = uint16(randi(60000, testCase.ImSize));
            writetiff(im, [dataPath, 'im.tif']);
            saveMIP_tiff([dataPath, 'im.tif'], [dataPath, 'out_MIP_z.tif'], 'axis', [1, 0, 1]);
            testCase.verifyEqual(readtiff([dataPath, 'out_MIP_y.tif']), squeeze(max(im, [], 1)));
            testCase.verifyEqual(readtiff([dataPath, 'out_MIP_z.tif']), max(im, [], 3));
            testCase.verifyFalse(exist([dataPath, 'out_MIP_x.tif'], 'file') > 0);
        end

        function slabAnotherWriterMadeIsKept(testCase)
            % a MIP slab that another writer (of the same result) had put in
            % place was deleted and replaced; it is kept, and the MIP of its
            % axis is made from it
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            mkdir([dataPath, 'MIPs']);
            slabFn = [dataPath, 'MIPs/', names{1}, '_y.zarr'];
            other = uint16(7 * ones(1, testCase.ImSize(2), testCase.ImSize(3)));
            writezarr(other, slabFn);
            evalc('XR_MIP_zarr([dataPath, names{1}, ''.zarr''], ''axis'', [1, 1, 1], ''parseCluster'', false)');
            testCase.verifyEqual(readzarr(slabFn), other, 'The other writer''s slab was replaced.');
            testCase.verifyEqual(readtiff([dataPath, 'MIPs/', names{1}, '_MIP_y.tif']), squeeze(other));
            testCase.verifyEqual(readtiff([dataPath, 'MIPs/', names{1}, '_MIP_z.tif']), max(ims{1}, [], 3));
            testCase.verifyEmpty(dir([dataPath, 'MIPs/', names{1}, '_y_*']), 'The temporary slab was not removed.');
            testCase.verifyEmpty(dir([slabFn, '/*.zarr']), 'The temporary slab was moved into the other''s.');
        end

        function largeZarrJobsGoAfterTheirBatches(testCase)
            % the job of a large file waits for its MIP batch jobs: it asks
            % to go after them (nice, coordinator_nice); the jobs of other
            % files don't (the framework is a spy)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            createzarr([dataPath, 'scan_CamA_ch0_t0000.zarr'], dataSize=[512, 512, 64], blockSize=[64, 64, 64], ...
                dtype='uint16');
            testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.applyFixture(pktest.FunctionStubFixture('coordinator_nice', 1000));
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture());
            for largeFile = [true, false]
                XR_MIP_wrapper(dataPath, 'axis', [1, 1, 1], 'channelPatterns', {'CamA_ch0'}, 'zarrFile', true, ...
                    'largeFile', largeFile, 'parseCluster', true, 'resultDirName', sprintf('MIPs_%d', largeFile), ...
                    'jobLogDir', [dataPath, 'logs']);
            end
            testCase.verifyEqual(spy.getOption(1, 'nice'), 1000);
            testCase.verifyEqual(spy.getOption(2, 'nice'), 0);
        end

        function largeZarrJobsAskForWhatTheyHold(testCase)
            % the job of a large file (XR_MIP_zarr) runs its MIP batches as
            % jobs (1 CPU each). It asked for 16 B per voxel of the whole
            % input: whole nodes for large data. Now it asks for MATLAB with
            % one CPU, and if it computes batches itself (masterCompute, the
            % config's if it sets it), for cpusPerTask (default 3) and a
            % batch's estimate (2048^3 uint16, 64 GB in MATLAB jobs). The
            % input is a 34 G-voxel header (the framework is a spy).
            dataPath = testCase.createTempFolder();
            createzarr([dataPath, 'scan_CamA_ch0_t0000.zarr'], dataSize=[4096, 4096, 2048], blockSize=[256, 256, 256], ...
                dtype='uint16');
            offConfig = [dataPath, 'cpu_config.json'];
            fid = fopen(offConfig, 'w');
            fprintf(fid, '{"masterCompute": false}');
            fclose(fid);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture());
            cases = {true, ''; false, ''; true, offConfig};
            for i = 1 : size(cases, 1)
                XR_MIP_wrapper(dataPath, 'axis', [1, 1, 1], 'channelPatterns', {'CamA_ch0'}, 'zarrFile', true, ...
                    'largeFile', true, 'parseCluster', false, 'masterCompute', cases{i, 1}, 'configFile', cases{i, 2}, ...
                    'resultDirName', sprintf('MIPs_%d', i));
            end
            testCase.verifyEqual(spy.getOption(1, 'cpusPerTask'), 3);
            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), 64, 'AbsTol', 1e-9);
            for i = 2 : 3
                testCase.verifyEqual(spy.getOption(i, 'cpusPerTask'), 1, sprintf('Case %d.', i));
                testCase.verifyEqual(spy.getOption(i, 'memAllocate'), 2, sprintf('Case %d.', i));
            end
            calls = spy.getCalls();
            testCase.verifyFalse(any(contains(calls{1}.funcStrs, '''masterCompute'',false')));
            testCase.verifyTrue(all(contains(calls{2}.funcStrs, '''masterCompute'',false')));
        end

        function mipsAreMovedIntoPlace(testCase)
            % XR_MIP_wrapper's tasks count as done once their MIPs exist, so a
            % MIP written in place could be seen, and its job cancelled,
            % while half written
            dataPath = testCase.createTempFolder();
            im = uint16(randi(60000, testCase.ImSize));
            writetiff(im, [dataPath, 'im.tif']);
            writezarr(im, [dataPath, 'im.zarr']);
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            saveMIP_tiff([dataPath, 'im.tif'], [dataPath, 'tif_MIP_z.tif'], 'axis', [1, 1, 1]);
            saveMIP_zarr([dataPath, 'im.zarr'], [dataPath, 'zarr_MIP_z.tif'], 'axis', [1, 1, 1]);
            for f = {'tif', 'zarr'}
                testCase.verifyMIPs([dataPath, f{1}], im, 'uint16');
                for a = {'y', 'x', 'z'}
                    fn = sprintf('%s%s_MIP_%s.tif', dataPath, f{1}, a{1});
                    testCase.verifyTrue(moves.wasMoved(fn), sprintf('%s was written in place.', fn));
                end
            end
            testCase.verifyEqual(numel(dir([dataPath, '*.tif'])), 7, 'Temporary files were left behind.');
        end

        function zarrIsReadInSlabsOfWholeChunks(testCase, zarrType)
            % saveMIP_zarr reads zarr files in z slabs of whole chunk layers
            % with at least 64 chunks (3 slabs of one layer here, 2 slabs of
            % two layers in the box) and keeps a running max along z. For
            % single, lines of NaN and of -0 and 0, where max takes the first
            % (0 first in the second and third slab here).
            dataPath = testCase.createTempFolder();
            im = cast(randi(60000, [128, 128, 12]), zarrType);
            if strcmp(zarrType, 'single')
                im(rand(size(im)) < 0.05) = NaN;
                im(3, :, :) = NaN;
                im(:, 5, :) = -0;
                im(:, 5, [2, 9]) = 0;
                im(7, 9, :) = -0;
                im(7, 9, [5, end]) = 0;
            end
            writezarr(im, [dataPath, 'im.zarr'], 'blockSize', [16, 16, 4]);
            saveMIP_zarr([dataPath, 'im.zarr'], [dataPath, 'all_MIP_z.tif'], 'axis', [1, 1, 1]);
            saveMIP_zarr([dataPath, 'im.zarr'], [dataPath, 'z_MIP_z.tif']);
            b = [5, 7, 3, 120, 110, 10];
            saveMIP_zarr([dataPath, 'im.zarr'], [dataPath, 'bbox_MIP_z.tif'], 'axis', [1, 1, 1], 'inputBbox', b);
            names = {'y', 'x', 'z'};
            for i = 1 : 3
                testCase.verifySameBits(readtiff(sprintf('%sall_MIP_%s.tif', dataPath, names{i})), squeeze(max(im, [], i)), ...
                    sprintf('MIP along %s', names{i}));
                testCase.verifySameBits(readtiff(sprintf('%sbbox_MIP_%s.tif', dataPath, names{i})), ...
                    squeeze(max(im(b(1) : b(4), b(2) : b(5), b(3) : b(6)), [], i)), sprintf('MIP along %s in the box', names{i}));
            end
            testCase.verifySameBits(readtiff([dataPath, 'z_MIP_z.tif']), max(im, [], 3), 'MIP along z alone');
            testCase.verifyFalse(exist([dataPath, 'z_MIP_y.tif'], 'file') > 0);
        end

        function maxPoolingMatchesBlockwiseMax(testCase, poolSize, poolType)
            im = cast(rand([23, 19, 11]) * 1000, poolType);
            expected = testCase.blockMax(im, poolSize);
            testCase.verifyEqual(max_pooling_3d_mex(im, poolSize), expected, 'max_pooling_3d_mex');
            testCase.verifyEqual(max_pooling_3d(im, poolSize), expected, 'max_pooling_3d');
        end

        function minimumOfARegion(testCase, minType)
            % min_bbox_3d: the minimum of a region, by min_bbox_3d_mex (Linux
            % only) or by the MATLAB code (elsewhere, for int16, or where the
            % mex fails), which returned the minimum along y instead
            if strcmp(minType, 'logical')
                im = rand([20, 30, 10]) > 0.05;
            else
                im = cast(rand([20, 30, 10]) * 200 + 1, minType);
            end
            bbox = [3, 4, 2, 15, 25, 9];
            expected = min(im(3 : 15, 4 : 25, 2 : 9), [], 'all');
            % a region without voxels: the largest value of the class
            emptyBbox = [5, 1, 1, 4, 30, 10];
            if strcmp(minType, 'logical')
                largest = true;
            elseif any(strcmp(minType, {'single', 'double'}))
                largest = realmax(minType);
            else
                largest = intmax(minType);
            end
            testCase.verifyEqual(min_bbox_3d(im, bbox), expected);
            testCase.verifyEqual(min_bbox_3d(im, emptyBbox), largest);
            testCase.applyFixture(pktest.FunctionStubFixture('min_bbox_3d_mex'));
            testCase.expectMexFallback('min_bbox_3d_mex');
            testCase.verifyEqual(min_bbox_3d(im, bbox), expected, 'without the mex');
            testCase.verifyEqual(min_bbox_3d(im, emptyBbox), largest, 'without the mex');
        end

        function minimumOfARegionSkipsNaN(testCase)
            % as the mex: NaN is skipped, and a region of NaN only or without
            % voxels gives the largest value of the class
            im = single(rand([20, 30, 10]) + 1);
            im(2 : 9, :, :) = NaN;
            im(12, 7, 3) = 0.5;
            boxes = {[1, 1, 1, 20, 30, 10], [2, 1, 1, 9, 30, 10], [5, 1, 1, 4, 30, 10]};
            expected = {single(0.5), realmax('single'), realmax('single')};
            for useMex = [true, false]
                if ~useMex
                    testCase.applyFixture(pktest.FunctionStubFixture('min_bbox_3d_mex'));
                    testCase.expectMexFallback('min_bbox_3d_mex');
                end
                for i = 1 : numel(boxes)
                    testCase.verifyEqual(min_bbox_3d(im, boxes{i}), expected{i}, sprintf('bbox %s, mex %d', mat2str(boxes{i}), useMex));
                end
            end
        end
    end

    methods (Test, TestTags = {'Unit', 'KnownIssue'})
        function maxPoolingOfNegativeValues(testCase)
            % KNOWN ISSUE (may be intended): max_pooling_3d_mex initializes
            % each block maximum with 0, so blocks with only negative values
            % return 0 instead of their maximum. The MATLAB fallback
            % (max_pooling_3d.m) zero-pads edge blocks, so it does the same
            % there. Non-negative input may be an intended precondition (the
            % mex has no signed integer types, and stitching treats 0 as no
            % data); if not, initialize with the first voxel of the block.
            im = -single(1 : 8)';
            im = reshape(im, 2, 2, 2);
            testCase.verifyEqual(max_pooling_3d_mex(im, [2, 2, 2]), single(-1));
        end
    end

    methods (Access = private)
        function [dataPath, ims, names] = writeDataset(testCase, format, dtype)
            if nargin < 3
                dtype = 'uint16';
            end
            dataPath = testCase.createTempFolder();
            names = {'scan_CamA_ch0_t0000', 'scan_CamB_ch0_t0000', 'scan_CamA_ch1_t0000'};
            ims = cell(1, 3);
            for i = 1 : 3
                ims{i} = cast(randi(60000, testCase.ImSize), dtype);
                if strcmp(format, 'zarr')
                    writezarr(ims{i}, [dataPath, names{i}, '.zarr'], 'blockSize', [16, 16, 4]);
                else
                    writetiff(ims{i}, [dataPath, names{i}, '.tif']);
                end
            end
        end

        function [] = verifySameBits(testCase, actual, expected, name)
            % equal values, and for floats equal bits (NaN payloads, -0)
            testCase.verifyClass(actual, class(expected), name);
            testCase.verifySize(actual, size(expected), name);
            if isa(actual, class(expected)) && isequal(size(actual), size(expected))
                if isfloat(expected)
                    testCase.verifyEqual(typecast(actual(:), 'uint8'), typecast(expected(:), 'uint8'), name);
                else
                    testCase.verifyEqual(actual, expected, name);
                end
            end
        end

        function [] = verifyMIPs(testCase, base, im, dtype)
            mipY = readtiff([base, '_MIP_y.tif']);
            mipX = readtiff([base, '_MIP_x.tif']);
            mipZ = readtiff([base, '_MIP_z.tif']);
            testCase.verifyClass(mipZ, dtype);
            testCase.verifyEqual(mipY, cast(squeeze(max(im, [], 1)), dtype), 'MIP along y');
            testCase.verifyEqual(mipX, cast(squeeze(max(im, [], 2)), dtype), 'MIP along x');
            testCase.verifyEqual(mipZ, cast(max(im, [], 3), dtype), 'MIP along z');
        end
    end

    methods (Static, Access = private)
        function [out] = blockMax(im, pool)
            sz = size(im, 1 : 3);
            outSize = ceil(sz ./ pool);
            out = zeros(outSize, class(im));
            for k = 1 : outSize(3)
                for j = 1 : outSize(2)
                    for i = 1 : outSize(1)
                        ys = (i - 1) * pool(1) + 1 : min(i * pool(1), sz(1));
                        xs = (j - 1) * pool(2) + 1 : min(j * pool(2), sz(2));
                        zs = (k - 1) * pool(3) + 1 : min(k * pool(3), sz(3));
                        out(i, j, k) = max(im(ys, xs, zs), [], 'all');
                    end
                end
            end
        end
    end
end
