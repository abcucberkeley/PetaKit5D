classdef CropMexTest < pktest.TestCase
    % crop/indexing mex functions and their .m wrappers against plain MATLAB
    % indexing. bbox convention: [ymin, xmin, zmin, ymax, xmax, zmax]
    % (4-D: [ymin, xmin, zmin, tmin, ymax, xmax, zmax, tmax]), 1-based and
    % inclusive. The indexing*_mex functions write into their first input in
    % place, so expected values are built from independent copies.

    properties (TestParameter)
        dtype = {'logical', 'uint8', 'uint16', 'single', 'double'}
        bbox3d = struct('interior', [3, 5, 2, 20, 17, 9], 'full', [1, 1, 1, 23, 19, 11], ...
            'singleVoxel', [7, 8, 9, 7, 8, 9], 'singleSlice', [1, 1, 11, 23, 19, 11], 'lastRow', [23, 1, 1, 23, 19, 11])
    end

    properties (Constant)
        ImSize = [23, 19, 11]
        ImSize4d = [13, 11, 7, 3]
    end

    methods (Test, TestTags = {'Unit'})
        function crop3dMatchesIndexing(testCase, dtype, bbox3d)
            im = testCase.randomImage(testCase.ImSize, dtype);
            b = bbox3d;
            expected = im(b(1) : b(4), b(2) : b(5), b(3) : b(6));
            testCase.verifyEqual(crop3d_mex(im, b), expected, 'crop3d_mex');
            testCase.verifyEqual(crop3d(im, b), expected, 'crop3d');
        end

        function crop3dRejectsInvalidBbox(testCase)
            im = testCase.randomImage(testCase.ImSize, 'uint16');
            testCase.verifyError(@() crop3d_mex(im, [1, 1, 1, 5, 5]), 'crop:inputError');
            testCase.verifyError(@() crop3d_mex(im, [0, 1, 1, 5, 5, 5]), 'crop:inputError');
            testCase.verifyError(@() crop3d_mex(im, [1, 1, 1, testCase.ImSize + [0, 0, 1]]), 'crop:inputError');
            testCase.verifyError(@() crop3d_mex(int16(im), [1, 1, 1, 5, 5, 5]), 'crop:dataTypeError');
        end

        function crop4dMatchesIndexing(testCase, dtype)
            im = testCase.randomImage(testCase.ImSize4d, dtype);
            bboxes = [2, 3, 2, 1, 12, 9, 6, 3; 1, 1, 1, 1, testCase.ImSize4d; 5, 5, 5, 2, 5, 5, 5, 2];
            for i = 1 : size(bboxes, 1)
                b = bboxes(i, :);
                expected = im(b(1) : b(5), b(2) : b(6), b(3) : b(7), b(4) : b(8));
                testCase.verifyEqual(crop4d_mex(im, b), expected, sprintf('crop4d_mex %s', mat2str(b)));
                testCase.verifyEqual(crop4d(im, b), expected, sprintf('crop4d %s', mat2str(b)));
            end
        end

        function crop4dRejectsInvalidBbox(testCase)
            im = testCase.randomImage(testCase.ImSize4d, 'uint16');
            testCase.verifyError(@() testCase.callWithOutput(@crop4d_mex, im, [1, 1, 1, 1, 2, 2, 2]), 'crop:inputError');
            testCase.verifyError(@() testCase.callWithOutput(@crop4d_mex, im, [1, 1, 1, 1, testCase.ImSize4d + [1, 0, 0, 0]]), 'crop:inputError');
            testCase.verifyError(@() testCase.callWithOutput(@crop4d_mex, im, [3, 1, 1, 1, 2, 2, 2, 1]), 'crop:inputError');
        end

        function indexing3dAssignsRegionInPlace(testCase, dtype)
            if strcmp(dtype, 'logical')
                A = false(testCase.ImSize);
                expected = false(testCase.ImSize);
            else
                A = zeros(testCase.ImSize, dtype);
                expected = zeros(testCase.ImSize, dtype);
            end
            b = [3, 5, 2, 20, 17, 9];
            R = testCase.randomImage(b(4 : 6) - b(1 : 3) + 1, dtype);
            indexing3d_mex(A, b, R);
            expected(b(1) : b(4), b(2) : b(5), b(3) : b(6)) = R;
            testCase.verifyEqual(A, expected);
        end

        function indexing3dRejectsMismatchedRegion(testCase)
            A = zeros(testCase.ImSize, 'uint16');
            b = [3, 5, 2, 20, 17, 9];
            testCase.verifyError(@() indexing3d_mex(A, b, zeros(3, 3, 3, 'uint16')), 'indexing:inputError');
            testCase.verifyError(@() indexing3d_mex(A, b, zeros(b(4 : 6) - b(1 : 3) + 1, 'single')), 'indexing:inputError');
        end

        function indexing4dWrapperMatchesIndexing(testCase, dtype)
            % the wrapper falls back to MATLAB indexing if the mex fails
            [A, expected, b, R] = testCase.indexing4dInputs(dtype);
            expected(b(1) : b(5), b(2) : b(6), b(3) : b(7), b(4) : b(8)) = R;
            testCase.verifyEqual(indexing4d(A, R, b), expected);
        end

        function indexing4dCropAssignsSubregion(testCase)
            [A, expected, b] = testCase.indexing4dInputs('uint16');
            R = testCase.randomImage([15, 15, 9, 3], 'uint16');
            rb = [2, 3, 1, 1, b(5 : 8) - b(1 : 4) + [2, 3, 1, 1]];
            indexing4d_crop_mex(A, b, R, rb);
            expected(b(1) : b(5), b(2) : b(6), b(3) : b(7), b(4) : b(8)) = R(rb(1) : rb(5), rb(2) : rb(6), rb(3) : rb(7), rb(4) : rb(8));
            testCase.verifyEqual(A, expected);
        end

        function trimBorderModes(testCase)
            im = testCase.randomImage(testCase.ImSize, 'uint16');
            bs = [2, 3, 1];
            testCase.verifyEqual(trimBorder(im, bs, 'both'), im(3 : end - 2, 4 : end - 3, 2 : end - 1));
            testCase.verifyEqual(trimBorder(im, bs, 'pre'), im(3 : end, 4 : end, 2 : end));
            testCase.verifyEqual(trimBorder(im, bs, 'post'), im(1 : end - 2, 1 : end - 3, 1 : end - 1));
        end

        function indexing4dMexSupportsLogical(testCase)
            % logical input used to fall through to the uint8 branch and be rejected
            [A, expected, b, R] = testCase.indexing4dInputs('logical');
            expected(b(1) : b(5), b(2) : b(6), b(3) : b(7), b(4) : b(8)) = R;
            try
                indexing4d_mex(A, b, R);
            catch ME
                testCase.verifyFail(sprintf('indexing4d_mex failed for logical input: %s', ME.message));
                return;
            end
            testCase.verifyEqual(A, expected);
        end
    end

    methods (Access = private)
        function [A, expected, b, R] = indexing4dInputs(testCase, dtype)
            sz = testCase.ImSize4d;
            if strcmp(dtype, 'logical')
                A = false(sz);
                expected = false(sz);
            else
                A = zeros(sz, dtype);
                expected = zeros(sz, dtype);
            end
            b = [2, 3, 2, 1, 12, 9, 6, 3];
            R = testCase.randomImage(b(5 : 8) - b(1 : 4) + 1, dtype);
        end
    end

    methods (Static, Access = private)
        function [im] = randomImage(sz, dtype)
            switch dtype
                case 'logical'
                    im = rand(sz) > 0.5;
                case {'single', 'double'}
                    im = cast(randn(sz) * 1000, dtype);
                otherwise
                    im = cast(randi(double(intmax(dtype)), sz), dtype);
            end
        end

        function [] = callWithOutput(fun, varargin)
            % call a function that requires an output argument
            out = fun(varargin{:}); %#ok<NASGU>
        end
    end
end
