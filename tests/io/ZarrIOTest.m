classdef ZarrIOTest < pktest.TestCase
    % readzarr/writezarr/createzarr (Cpp-Zarr mex) round trips, partial
    % reads and writes, compressors, chunk layouts, and .zarray metadata

    properties (TestParameter)
        dtype = {'uint8', 'uint16', 'single', 'double'}
        % non-divisible, a single chunk, and chunks larger than the image
        blockSize = struct('nonDivisible', [16, 16, 5], 'singleChunk', [37, 41, 13], 'larger', [64, 64, 64])
        compressor = {'zstd', 'lz4', 'lz4hc', 'blosclz', 'zlib', 'gzip'}
        order = {'F', 'C'}
        % the C order bug fixed in 53c4e20 only shows if the z chunk is larger
        % than the y chunk
        chunkShape = struct('zSmaller', [16, 8, 5], 'zLarger', [8, 16, 13])
    end

    properties (Constant)
        ImSize = [37, 41, 13]
    end

    methods (Test, TestTags = {'Unit'})
        function roundTripPreservesTypeAndValues(testCase, dtype, blockSize)
            im = testCase.randomImage(testCase.ImSize, dtype);
            fn = [testCase.createTempFolder(), 'im.zarr'];
            writezarr(im, fn, 'blockSize', blockSize);
            out = readzarr(fn);
            testCase.verifyClass(out, dtype);
            testCase.verifyEqual(out, im);
        end

        function compressorsRoundTrip(testCase, compressor)
            im = testCase.randomImage(testCase.ImSize, 'uint16');
            fn = [testCase.createTempFolder(), 'im.zarr'];
            createzarr(fn, 'dataSize', size(im), 'blockSize', [16, 16, 5], 'dtype', 'uint16', 'compressor', compressor);
            writezarr(im, fn);
            testCase.verifyEqual(readzarr(fn), im);
            meta = jsondecode(fileread([fn, '/.zarray']));
            if strcmp(compressor, 'gzip')
                testCase.verifyEqual(meta.compressor.id, 'gzip');
            else
                testCase.verifyEqual(meta.compressor.cname, compressor);
            end
        end

        function chunkOrderRoundTrip(testCase, order, chunkShape)
            % regression test for the C order chunk size issue (53c4e20)
            im = testCase.randomImage(testCase.ImSize, 'uint16');
            fn = [testCase.createTempFolder(), 'im.zarr'];
            createzarr(fn, 'dataSize', size(im), 'blockSize', chunkShape, 'dtype', 'uint16', 'order', order);
            writezarr(im, fn);
            testCase.verifyEqual(readzarr(fn), im);
            meta = jsondecode(fileread([fn, '/.zarray']));
            testCase.verifyEqual(meta.order, order);
        end

        function defaultMetadata(testCase)
            im = testCase.randomImage(testCase.ImSize, 'uint16');
            fn = [testCase.createTempFolder(), 'im.zarr'];
            writezarr(im, fn, 'blockSize', [16, 16, 5]);
            meta = jsondecode(fileread([fn, '/.zarray']));
            testCase.verifyEqual(meta.zarr_format, 2);
            testCase.verifyEqual(meta.shape(:)', testCase.ImSize);
            testCase.verifyEqual(meta.chunks(:)', [16, 16, 5]);
            testCase.verifyEqual(meta.dtype, '<u2');
            testCase.verifyEqual(meta.order, 'F');
            testCase.verifyEqual(meta.fill_value, 0);
            testCase.verifyEqual(meta.compressor.id, 'blosc');
            testCase.verifyEqual(meta.compressor.cname, 'zstd');
            testCase.verifyEqual(meta.compressor.clevel, 1);
            % chunk files are named by 0-based chunk index in y.x.z
            testCase.verifyTrue(exist([fn, '/2.2.2'], 'file') == 2);
            testCase.verifyFalse(exist([fn, '/3.0.0'], 'file') == 2);
        end

        function blockSizeIsCappedAtImageSize(testCase)
            fn = [testCase.createTempFolder(), 'im.zarr'];
            writezarr(testCase.randomImage(testCase.ImSize, 'uint16'), fn, 'blockSize', [64, 64, 64]);
            meta = jsondecode(fileread([fn, '/.zarray']));
            testCase.verifyEqual(meta.chunks(:)', testCase.ImSize);
        end

        function nestedDimensionSeparator(testCase)
            im = testCase.randomImage(testCase.ImSize, 'uint16');
            fn = [testCase.createTempFolder(), 'im.zarr'];
            writezarr(im, fn, 'blockSize', [16, 16, 5], 'dimSeparator', '/');
            testCase.verifyEqual(readzarr(fn), im);
            testCase.verifyTrue(exist([fn, '/2/2/2'], 'file') == 2, 'Chunks are not stored in nested folders.');
        end

        function shardedRoundTrip(testCase)
            im = testCase.randomImage(testCase.ImSize, 'uint16');
            fn = [testCase.createTempFolder(), 'im.zarr'];
            writezarr(im, fn, 'blockSize', [32, 32, 10], 'shardSize', [16, 16, 5]);
            testCase.verifyEqual(readzarr(fn), im);
        end

        function image2dIsStoredWithSingletonZ(testCase)
            im = testCase.randomImage([21, 17], 'uint16');
            fn = [testCase.createTempFolder(), 'im.zarr'];
            writezarr(im, fn);
            testCase.verifyEqual(readzarr(fn), im);
            testCase.verifyEqual(getImageSize(fn), [21, 17, 1]);
        end

        function readBoundingBox(testCase)
            im = testCase.randomImage(testCase.ImSize, 'uint16');
            fn = [testCase.createTempFolder(), 'im.zarr'];
            writezarr(im, fn, 'blockSize', [16, 16, 5]);
            % chunk-crossing, single voxel, single slice, and the whole image
            bboxes = [3, 5, 2, 30, 33, 11; 17, 17, 6, 17, 17, 6; 1, 1, 7, 37, 41, 7; 1, 1, 1, testCase.ImSize];
            for i = 1 : size(bboxes, 1)
                b = bboxes(i, :);
                testCase.verifyEqual(readzarr(fn, 'inputBbox', b), im(b(1) : b(4), b(2) : b(5), b(3) : b(6)), ...
                    sprintf('bbox %s', mat2str(b)));
            end
        end

        function readBoundingBoxPastTheEndErrors(testCase)
            fn = [testCase.createTempFolder(), 'im.zarr'];
            writezarr(testCase.randomImage(testCase.ImSize, 'uint16'), fn);
            testCase.verifyError(@() readzarr(fn, 'inputBbox', [1, 1, 1, testCase.ImSize + [1, 0, 0]]), ?MException);
        end

        function readMissingFileErrors(testCase)
            testCase.verifyError(@() readzarr([testCase.createTempFolder(), 'missing.zarr']), ?MException);
        end

        function partialWritesUpdateOnlyTheRegion(testCase)
            % also rewrites existing chunks, which takes a different code
            % path (rename over an existing file) than the first write. On
            % Windows, rewrites kept the old data before Cpp-Zarr 1.6.0.
            im = testCase.randomImage(testCase.ImSize, 'uint16');
            fn = [testCase.createTempFolder(), 'im.zarr'];
            createzarr(fn, 'dataSize', size(im), 'blockSize', [16, 16, 5], 'dtype', 'uint16');
            expected = zeros(size(im), 'uint16');
            testCase.verifyEqual(readzarr(fn), expected, 'A new zarr file must read as the fill value.');

            b = [5, 7, 2, 30, 20, 9];
            writezarr(im(b(1) : b(4), b(2) : b(5), b(3) : b(6)), fn, 'bbox', b, 'create', false);
            expected(b(1) : b(4), b(2) : b(5), b(3) : b(6)) = im(b(1) : b(4), b(2) : b(5), b(3) : b(6));
            testCase.verifyEqual(readzarr(fn), expected, 'First partial write.');

            b = [10, 1, 1, 37, 15, 13];
            region = testCase.randomImage(b(4 : 6) - b(1 : 3) + 1, 'uint16');
            writezarr(region, fn, 'bbox', b, 'create', false);
            expected(b(1) : b(4), b(2) : b(5), b(3) : b(6)) = region;
            testCase.verifyEqual(readzarr(fn), expected, 'Overlapping partial write did not replace existing chunks.');

            names = testCase.chunkFiles(fn);
            testCase.verifyTrue(all(~cellfun(@isempty, regexp(names, '^\d+\.\d+\.\d+$', 'once'))), ...
                sprintf('Unexpected (temporary?) files in the zarr folder: %s', strjoin(names, ', ')));
        end

        function overwriteReplacesExistingArray(testCase)
            fn = [testCase.createTempFolder(), 'im.zarr'];
            writezarr(testCase.randomImage([20, 20, 20], 'single'), fn);
            im = testCase.randomImage(testCase.ImSize, 'uint16');
            writezarr(im, fn, 'overwrite', true);
            testCase.verifyEqual(readzarr(fn), im);
        end

        function sparseWriteSkipsEmptyChunks(testCase)
            im = zeros(32, 32, 10, 'uint16');
            im(1 : 16, 1 : 16, 1 : 5) = 7;
            fn = [testCase.createTempFolder(), 'im.zarr'];
            writezarr(im, fn, 'blockSize', [16, 16, 5], 'sparseData', true);
            testCase.verifyEqual(readzarr(fn), im);
            testCase.verifyEqual(testCase.chunkFiles(fn), {'0.0.0'}, 'Only the nonzero chunk should be written.');
        end

        function createzarrRejectsUnsupportedType(testCase)
            fn = [testCase.createTempFolder(), 'im.zarr'];
            testCase.verifyError(@() createzarr(fn, 'dataSize', [4, 4, 4], 'dtype', 'int16'), ?MException);
        end
    end

    methods (Static, Access = private)
        function [names] = chunkFiles(fn)
            dir_info = dir(fn);
            names = {dir_info(~[dir_info.isdir]).name};
            names = names(~startsWith(names, '.'));
        end

        function [im] = randomImage(sz, dtype)
            if any(strcmp(dtype, {'single', 'double'}))
                im = cast(randn(sz) * 1000, dtype);
            else
                im = cast(randi(double(intmax(dtype)), sz), dtype);
            end
        end
    end
end
