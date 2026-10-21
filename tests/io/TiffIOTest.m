classdef TiffIOTest < pktest.TestCase
    % readtiff/writetiff (Cpp-Tiff mex with libtiff/saveastiff fallbacks),
    % readtiff_matlab, and interoperability with MATLAB's Tiff reader/writer

    properties (TestParameter)
        dtype = {'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32', 'single', 'double', 'uint64', 'int64'}
        imSize = struct('image2d', [37, 53], 'stack', [37, 53, 11], 'singleRow', [1, 64, 5])
        writeMode = struct('parallelLzw', {{'Compression', 'lzw', 'Mode', 'parallel'}}, ...
            'parallelNone', {{'Compression', 'none', 'Mode', 'parallel'}}, ...
            'libtiff', {{'Mode', 'libtiff'}})
    end

    methods (Test, TestTags = {'Unit'})
        function roundTripPreservesTypeAndValues(testCase, dtype)
            im = testCase.randomImage([37, 53, 11], dtype);
            fn = [testCase.createTempFolder(), 'im.tif'];
            writetiff(im, fn);
            out = readtiff(fn);
            testCase.verifyClass(out, dtype);
            testCase.verifyEqual(out, im);
        end

        function roundTripPreservesShape(testCase, imSize)
            im = testCase.randomImage(imSize, 'uint16');
            fn = [testCase.createTempFolder(), 'im.tif'];
            writetiff(im, fn);
            testCase.verifyEqual(readtiff(fn), im);
            testCase.verifyEqual(getImageSize(fn), [imSize, ones(1, 3 - numel(imSize))]);
        end

        function writeModesAreReadableByAllReaders(testCase, writeMode)
            im = testCase.randomImage([37, 53, 11], 'uint16');
            fn = [testCase.createTempFolder(), 'im.tif'];
            writetiff(im, fn, writeMode{:});
            testCase.verifyEqual(readtiff(fn), im, 'readtiff');
            testCase.verifyEqual(readtiff_matlab(fn), im, 'readtiff_matlab');
            testCase.verifyEqual(testCase.readWithMatlabTiff(fn), im, 'MATLAB Tiff');
        end

        function readsTiffWrittenByMatlab(testCase)
            % multi-page stacks written by imwrite (libtiff), with and
            % without compression
            im = testCase.randomImage([37, 53, 11], 'uint16');
            dataPath = testCase.createTempFolder();
            compressions = {'none', 'lzw', 'packbits', 'deflate'};
            for i = 1 : numel(compressions)
                fn = sprintf('%s%s.tif', dataPath, compressions{i});
                for z = 1 : size(im, 3)
                    imwrite(im(:, :, z), fn, 'Compression', compressions{i}, 'WriteMode', 'append');
                end
                testCase.verifyEqual(readtiff(fn), im, compressions{i});
                testCase.verifyEqual(readtiff(fn, 'range', [4, 9]), im(:, :, 4 : 9), compressions{i});
            end
        end

        function rangeReadsInclusiveSlices(testCase)
            im = testCase.randomImage([37, 53, 11], 'uint16');
            fn = [testCase.createTempFolder(), 'im.tif'];
            writetiff(im, fn);
            testCase.verifyEqual(readtiff(fn, 'range', [3, 7]), im(:, :, 3 : 7));
            testCase.verifyEqual(readtiff(fn, 'range', [5, 5]), im(:, :, 5));
            testCase.verifyEqual(readtiff(fn, 'range', [1, 11]), im);
            testCase.verifyEqual(readtiff_matlab(fn, 3 : 7), im(:, :, 3 : 7));
        end

        function moreThan65535Pages(testCase)
            % regression test for the page number limit (d2a8d7d)
            im = testCase.randomImage([4, 4, 70000], 'uint8');
            fn = [testCase.createTempFolder(), 'pages.tif'];
            writetiff(im, fn);
            testCase.verifyEqual(getImageSize(fn), [4, 4, 70000]);
            testCase.verifyEqual(readtiff(fn), im);
            testCase.verifyEqual(readtiff(fn, 'range', [69990, 70000]), im(:, :, 69990 : 70000));
        end

        function overwritesExistingFile(testCase)
            fn = [testCase.createTempFolder(), 'im.tif'];
            writetiff(testCase.randomImage([20, 30, 6], 'uint16'), fn);
            im = testCase.randomImage([10, 12, 3], 'single');
            writetiff(im, fn);
            testCase.verifyEqual(readtiff(fn), im);
        end

        function createsMissingOutputFolders(testCase)
            fn = [testCase.createTempFolder(), 'a/b/c/im.tif'];
            im = testCase.randomImage([20, 30, 6], 'uint16');
            writetiff(im, fn);
            testCase.verifyEqual(readtiff(fn), im);
        end

        function readsImageJStack(testCase)
            % ImageJ/Fiji stacks are big-endian with contiguous pixel data and
            % one IFD per slice; also a regression test for range reads of
            % ImageJ files (f369d36)
            im = testCase.randomImage([20, 30, 9], 'uint16');
            fn = [testCase.createTempFolder(), 'imagej.tif'];
            testCase.writeImageJStack(im, fn, false);
            testCase.verifyEqual(testCase.readWithMatlabTiff(fn), im, 'Test fixture is not a valid TIFF stack.');

            testCase.verifyEqual(readtiff(fn), im);
            testCase.verifyEqual(readtiff(fn, 'range', [3, 6]), im(:, :, 3 : 6));
            testCase.verifyEqual(getImageSize(fn), size(im));
        end

        function readsSingleIfdImageJStack(testCase)
            % ImageJ saves stacks over 4 GB with only the first IFD; the
            % number of slices is in the 'images=' description entry. The
            % ARM macOS binaries read only the first slice before Cpp-Tiff
            % 1.6.1.
            im = testCase.randomImage([20, 30, 9], 'uint16');
            fn = [testCase.createTempFolder(), 'imagej.tif'];
            testCase.writeImageJStack(im, fn, true);
            testCase.verifyEqual(testCase.readWithMatlabTiff(fn), im(:, :, 1), 'Test fixture is not a valid TIFF.');

            testCase.verifyEqual(readtiff(fn), im);
            testCase.verifyEqual(getImageSize(fn), size(im));
            % regression test for range reads of such stacks (f369d36), which
            % used to fail the bounds check; the content of range reads is
            % checked by rangeReadOfImageJStack
            testCase.verifyEqual(size(readtiff(fn, 'range', [3, 6])), [20, 30, 4]);
        end

        function rangeReadOfImageJStack(testCase)
            % before Cpp-Tiff 1.6.1, range reads of single-IFD ImageJ stacks
            % ignored the start of the range (slices 1-4 for [3, 6]), so
            % cropping, MIPs and resampling with inputBbox got the wrong
            % slices of large Fiji-saved stacks
            im = testCase.randomImage([20, 30, 9], 'uint16');
            fn = [testCase.createTempFolder(), 'imagej.tif'];
            testCase.writeImageJStack(im, fn, true);
            testCase.verifyEqual(readtiff(fn, 'range', [3, 6]), im(:, :, 3 : 6), ...
                'readtiff with a range returned the wrong slices of an ImageJ stack.');
        end

        function unsupportedInputIsNotSilentlyDropped(testCase)
            % writetiff used to ignore the error code of its fallback writer
            dataPath = testCase.createTempFolder();
            inputs = {zeros(4, 4, 2, 2, 'uint16'), zeros(0, 0, 'uint16')};
            names = {'4-D', 'empty'};
            for i = 1 : numel(inputs)
                fn = sprintf('%s%d.tif', dataPath, i);
                written = false;
                try
                    writetiff(inputs{i}, fn);
                    written = exist(fn, 'file') == 2;
                catch
                    % raising an error is acceptable
                    written = true;
                end
                testCase.verifyTrue(written, sprintf('writetiff of %s data neither wrote a file nor raised an error.', names{i}));
            end
        end
    end

    methods (Static, Access = private)
        function [im] = randomImage(sz, dtype)
            if any(strcmp(dtype, {'single', 'double'}))
                im = cast(randn(sz) * 1000, dtype);
            else
                lo = double(intmin(dtype));
                hi = min(double(intmax(dtype)), 2 ^ 52);
                im = cast(lo + floor(rand(sz) * (hi - lo)), dtype);
            end
        end

        function [im] = readWithMatlabTiff(fn)
            t = Tiff(fn, 'r');
            cleanupObj = onCleanup(@() close(t));
            slices = {};
            while true
                slices{end + 1} = t.read(); %#ok<AGROW>
                if t.lastDirectory()
                    break;
                end
                t.nextDirectory();
            end
            im = cat(3, slices{:});
        end

        function [] = writeImageJStack(im, fn, singleIFD)
            % uint16 stack in the layout ImageJ uses: big-endian, the
            % description 'ImageJ=...\nimages=N' in the first IFD, pixel data
            % of all slices contiguous, one IFD per slice after the data (or
            % only the first IFD, as ImageJ does for stacks over 4 GB)
            [ny, nx, nz] = size(im);
            nIFD = nz;
            if singleIFD
                nIFD = 1;
            end
            desc = [sprintf('ImageJ=1.54f\nimages=%d\nslices=%d\nloop=false\n', nz, nz), char(0)];
            nEntries = [11, 10 * ones(1, nIFD - 1)];
            ifdBytes = 2 + nEntries * 12 + 4;
            descOffset = 8 + ifdBytes(1);
            dataOffset = descOffset + numel(desc);
            sliceBytes = ny * nx * 2;
            ifdOffsets = [8, dataOffset + nz * sliceBytes + [0, cumsum(ifdBytes(2 : end - 1))]];

            fid = fopen(fn, 'w', 'ieee-be');
            cleanupObj = onCleanup(@() fclose(fid));
            fwrite(fid, 'MM', 'char');
            fwrite(fid, 42, 'uint16');
            fwrite(fid, 8, 'uint32');
            for z = 1 : nIFD
                if z == 2
                    fwrite(fid, desc, 'char');
                    fwrite(fid, permute(im, [2, 1, 3]), 'uint16');
                end
                % tag, type (3: short, 4: long, 2: ascii), count, value
                entries = {254, 4, 1, 0; 256, 4, 1, nx; 257, 4, 1, ny; 258, 3, 1, 16; 259, 3, 1, 1; ...
                    262, 3, 1, 1; 270, 2, numel(desc), descOffset; ...
                    273, 4, 1, dataOffset + (z - 1) * sliceBytes; 277, 3, 1, 1; 278, 4, 1, ny; 279, 4, 1, sliceBytes};
                if z > 1
                    entries(7, :) = [];
                end
                fwrite(fid, size(entries, 1), 'uint16');
                for i = 1 : size(entries, 1)
                    fwrite(fid, entries{i, 1}, 'uint16');
                    fwrite(fid, entries{i, 2}, 'uint16');
                    fwrite(fid, entries{i, 3}, 'uint32');
                    if entries{i, 2} == 3
                        fwrite(fid, [entries{i, 4}, 0], 'uint16');
                    else
                        fwrite(fid, entries{i, 4}, 'uint32');
                    end
                end
                if z < nIFD
                    fwrite(fid, ifdOffsets(z + 1), 'uint32');
                else
                    fwrite(fid, 0, 'uint32');
                end
            end
            if nIFD == 1
                fwrite(fid, desc, 'char');
                fwrite(fid, permute(im, [2, 1, 3]), 'uint16');
            end
        end
    end
end
