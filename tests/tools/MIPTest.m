classdef MIPTest < pktest.TestCase
    % maximum intensity projections: XR_MIP_wrapper (tiff, zarr and large
    % zarr paths), saveMIP_tiff, and max_pooling_3d

    properties (Constant)
        ImSize = [40, 36, 12]
        Channels = {'CamA_ch0', 'CamB_ch0'}
    end

    properties (TestParameter)
        poolSize = struct('cube', [2, 2, 2], 'anisotropic', [3, 2, 1], 'nonDivisible', [4, 5, 5])
        poolType = {'uint16', 'single'}
        zarrType = {'uint16', 'single'}
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
