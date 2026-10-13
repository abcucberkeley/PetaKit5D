classdef ImageInfoTest < pktest.TestCase
    % image metadata without loading the data: getImageSize,
    % getImageSizeBatch, getImageDataType, getZarrInfo

    properties (TestParameter)
        dtype = {'uint8', 'uint16', 'single', 'double'}
    end

    methods (Test, TestTags = {'Unit'})
        function tiffSizeAndType(testCase, dtype)
            dataPath = testCase.createTempFolder();
            fn = [dataPath, 'im.tif'];
            writetiff(zeros(23, 31, 7, dtype), fn);
            testCase.verifyEqual(getImageSize(fn), [23, 31, 7]);
            testCase.verifyEqual(char(getImageDataType(fn)), dtype);

            fn2d = [dataPath, 'im2d.tif'];
            writetiff(zeros(23, 31, dtype), fn2d);
            testCase.verifyEqual(getImageSize(fn2d), [23, 31, 1]);
        end

        function zarrSizeAndType(testCase, dtype)
            dataPath = testCase.createTempFolder();
            fn = [dataPath, 'im.zarr'];
            writezarr(zeros(23, 31, 7, dtype), fn, 'blockSize', [16, 16, 4]);
            testCase.verifyEqual(getImageSize(fn), [23, 31, 7]);
            testCase.verifyEqual(char(getImageDataType(fn)), dtype);

            info = getZarrInfo(fn);
            testCase.verifyEqual(info.size, [23, 31, 7]);
            testCase.verifyEqual(info.blockSize, [16, 16, 4]);
            testCase.verifyEqual(info.dtype, dtype);
            testCase.verifyEqual(info.order, 'F');
            testCase.verifyEqual(info.zarr_format, 2);
        end

        function zarrFolderWithoutExtension(testCase)
            % zarr detection also works by the .zarray file
            fn = [testCase.createTempFolder(), 'im_zarr_folder'];
            createzarr(fn, 'dataSize', [23, 31, 7], 'blockSize', [16, 16, 4], 'dtype', 'single');
            testCase.verifyEqual(getImageSize(fn), [23, 31, 7]);
            testCase.verifyEqual(char(getImageDataType(fn, true)), 'single');
        end

        function batchSizesOfMixedFiles(testCase)
            dataPath = testCase.createTempFolder();
            fns = {[dataPath, 'a.tif'], [dataPath, 'b.tif'], [dataPath, 'c.zarr']};
            writetiff(zeros(10, 20, 3, 'uint16'), fns{1});
            writetiff(zeros(11, 21, 4, 'uint16'), fns{2});
            writezarr(zeros(12, 22, 5, 'uint16'), fns{3});
            testCase.verifyEqual(getImageSizeBatch(fns, false), [10, 20, 3; 11, 21, 4; 12, 22, 5]);
        end

        function getZarrInfoOfMissingFileErrors(testCase)
            testCase.verifyError(@() getZarrInfo([testCase.createTempFolder(), 'missing.zarr']), ?MException);
        end

        function unknownFormatErrors(testCase)
            fn = [testCase.createTempFolder(), 'im.png'];
            imwrite(zeros(4, 4, 'uint8'), fn);
            testCase.verifyError(@() getImageDataType(fn), ?MException);
        end
    end

    methods (Test, TestTags = {'Unit', 'KnownIssue'})
        function queryingMetadataKeepsWarningState(testCase)
            % KNOWN ISSUE: getImageSize.m (line 15) and getZarrInfo.m (line 8)
            % call warning('off', 'all') and never restore it. Nearly every
            % wrapper calls getImageSize, so after the first wrapper call all
            % warnings stay disabled for the rest of the MATLAB session.
            % Fix: ws = warning('off', 'all'); c = onCleanup(@() warning(ws));
            % or remove the call.
            dataPath = testCase.createTempFolder();
            writetiff(zeros(4, 4, 2, 'uint16'), [dataPath, 'im.tif']);
            writezarr(zeros(4, 4, 2, 'uint16'), [dataPath, 'im.zarr']);

            warning('on', 'all');
            getImageSize([dataPath, 'im.tif']);
            testCase.verifyEqual(warning('query', 'MATLAB:singularMatrix').state, 'on', ...
                'getImageSize turned all warnings off.');
            warning('on', 'all');
            getZarrInfo([dataPath, 'im.zarr']);
            testCase.verifyEqual(warning('query', 'MATLAB:singularMatrix').state, 'on', ...
                'getZarrInfo turned all warnings off.');
        end
    end
end
