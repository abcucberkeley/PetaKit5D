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

        function maxPoolingMatchesBlockwiseMax(testCase, poolSize, poolType)
            im = cast(rand([23, 19, 11]) * 1000, poolType);
            expected = testCase.blockMax(im, poolSize);
            testCase.verifyEqual(max_pooling_3d_mex(im, poolSize), expected, 'max_pooling_3d_mex');
            testCase.verifyEqual(max_pooling_3d(im, poolSize), expected, 'max_pooling_3d');
        end
    end

    methods (Test, TestTags = {'Unit', 'KnownIssue'})
        function maxPoolingOfNegativeValues(testCase)
            % KNOWN ISSUE: max_pooling_3d_mex initializes each block maximum
            % with 0, so blocks with only negative values return 0 instead of
            % their maximum (the MATLAB fallback in max_pooling_3d.m is
            % correct). Affects 'max' resampling of data with negative values,
            % e.g., after background subtraction. Fix: initialize with the
            % first voxel of the block (or -inf).
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
