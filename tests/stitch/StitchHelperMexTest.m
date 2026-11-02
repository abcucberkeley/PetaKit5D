classdef StitchHelperMexTest < pktest.TestCase
    % stitching helpers that run a mex with MATLAB code as the fallback and
    % that the stitching tests don't reach: the crop of a long registration
    % region to its brightest part (crop_subregion_by_intensity, crop3d_mex)
    % and the check that the major tiles of a batch have data at every voxel
    % (check_major_tile_valid, any_4th_dim_mex); the edge erosion of tiles
    % (erodeVolumeBy2DProjection, indexing3d_mex) and the integral images
    % of the registration (integral_image_3d_mex). pktest.TestCase fails a
    % test in which the mex call fails.

    properties (TestParameter)
        % uint16 and single tiles run in the mex, others in the MATLAB code
        tileClass = {'uint16', 'single', 'uint8', 'double'}
    end

    methods (Test, TestTags = {'Unit'})
        function longRegionIsCroppedAroundItsBrightestPart(testCase)
            region = 100 + uint16(randi(20, [30, 40, 12]));
            region(20 : 22, 31 : 33, 5 : 7) = 5000;
            % axes longer than 10 are cropped to 10 voxels around the first
            % brightest position (or from the start if it is within 5)
            [cropped, bbox] = crop_subregion_by_intensity(region, 10);
            testCase.verifyEqual(bbox, [16, 27, 1, 25, 36, 10]);
            testCase.verifyEqual(cropped, region(16 : 25, 27 : 36, 1 : 10));
            % regions within the limit are kept
            [kept, bbox] = crop_subregion_by_intensity(region, 40);
            testCase.verifyEqual(kept, region);
            testCase.verifyEqual(bbox, [1, 1, 1, 30, 40, 12]);
        end

        function majorTilesWithDataEverywhereAreValid(testCase, tileClass)
            % tile 1 and 3 are major; tile 2 has no data (0) and tile 4 only
            % some, which doesn't matter
            fmat = zeros(8, 9, 5, 4, tileClass);
            fmat(:, :, :, 1) = 7;
            fmat(:, :, 1 : 2, 4) = 3;
            major = logical([1, 0, 1, 0]);
            [valid, bbox] = check_major_tile_valid(fmat, major);
            testCase.verifyEqual(valid, true);
            testCase.verifyEqual(bbox, []);
        end

        function voxelWithoutDataInTheMajorTilesIsFound(testCase, tileClass)
            % one voxel (many threads: boundsOfVoxelsWithoutDataWithManyThreads)
            fmat = zeros(8, 9, 5, 3, tileClass);
            fmat(:, :, :, 1) = 7;
            fmat(:, :, :, 3) = 5;
            fmat(3, 6, 2, 1) = 0;
            fmat(3, 6, 2, 3) = 0;
            [valid, bbox] = check_major_tile_valid(fmat, logical([1, 0, 1]));
            testCase.verifyEqual(valid, false);
            testCase.verifyEqual(bbox, [3, 6, 2, 3, 6, 2]);
            % data in one major tile is enough
            fmat(3, 6, 2, 3) = 1;
            testCase.verifyEqual(check_major_tile_valid(fmat, logical([1, 0, 1])), true);
        end

        function boundsOfVoxelsWithoutDataWithManyThreads(testCase)
            % a block without data that all threads find and single voxels at
            % the bounds that one thread each finds: the mex updated the bounds
            % from all threads without synchronization and lost some of them
            % (wrong in most calls with 8 threads); the threads keep their own
            % now, combined at the end, as one thread would find them
            testCase.assumeMexAvailable('any_4th_dim_mex');
            fmat = 7 * ones(256, 256, 64, 2, 'uint16');
            fmat(100 : 150, 100 : 150, 10 : 50, :) = 0;
            fmat(1, 120, 30, :) = 0;
            fmat(256, 130, 33, :) = 0;
            fmat(120, 1, 40, :) = 0;
            fmat(125, 256, 20, :) = 0;
            fmat(110, 115, 1, :) = 0;
            fmat(111, 116, 64, :) = 0;
            nWrong = 0;
            for k = 1 : 20
                [valid, bbox] = check_major_tile_valid(fmat, true(1, 2));
                nWrong = nWrong + (valid || ~isequal(bbox, [1, 1, 1, 256, 256, 64]));
            end
            testCase.verifyEqual(nWrong, 0, sprintf('%d of 20 calls gave other bounds, the last %s.', nWrong, mat2str(bbox)));
        end

        function edgeErosionWithoutTheMexMatches(testCase)
            % erodeVolumeBy2DProjection zeroes esize slices at both ends of y
            % with indexing3d_mex; the MATLAB code that runs where the mex
            % fails zeroed one more at the end
            vol = uint16(randi([1, 1000], [30, 40, 20]));
            vol(:, 1 : 3, :) = 0;
            esizes = [1, 3, 15];
            expected = cell(size(esizes));
            for i = 1 : numel(esizes)
                e = esizes(i);
                expected{i} = erodeVolumeBy2DProjection(vol, e);
                testCase.verifyEqual(nnz(expected{i}([1 : e, end - e + 1 : end], :, :)), 0);
                if e < 15
                    testCase.verifyGreaterThan(nnz(expected{i}([e + 1, end - e], :, :)), 0);
                end
            end
            testCase.applyFixture(pktest.FunctionStubFixture('indexing3d_mex'));
            testCase.expectMexFallback('indexing3d_mex');
            for i = 1 : numel(esizes)
                testCase.verifyEqual(erodeVolumeBy2DProjection(vol, esizes(i)), expected{i}, sprintf('esize %d', esizes(i)));
            end
        end

        function integralImageMexErrorsHaveIdentifiers(testCase)
            % the mex's error identifiers had spaces, which MATLAB rejects
            % (MATLAB:mexapi:invalidErrorIdentifier, without the message)
            testCase.assumeMexAvailable('integral_image_3d_mex');
            A = single(rand(10, 12, 8));
            import matlab.unittest.constraints.Throws
            testCase.verifyThat(@() integral_image_3d_mex(double(A), [3, 3, 3]), ...
                Throws('integral_image_3d:dataTypeError', 'WhenNargoutIs', 1));
            testCase.verifyThat(@() integral_image_3d_mex(A, [3, 3]), Throws('integral_image_3d:inputError', 'WhenNargoutIs', 1));
            testCase.verifyThat(@() integral_image_3d_mex(A, [3, 0, 3]), Throws('integral_image_3d:inputError', 'WhenNargoutIs', 1));
            testCase.verifyError(@() integral_image_3d_mex(A, [3, 3, 3]), 'integral_image_3d:outputError');

            % other types are rejected by design: the MATLAB code runs without
            % a warning
            expected = integral_image_3d(A, [3, 4, 2]);
            mex_fallback('-reset');
            out = integral_image_3d(double(A), [3, 4, 2]);
            records = mex_fallback('-log');
            testCase.verifyEqual({records.mex, records.kind, records.identifier}, ...
                {'integral_image_3d_mex', 'designed', 'integral_image_3d:dataTypeError'});
            testCase.verifyClass(out, 'double');
            testCase.verifyEqual(out, double(expected), 'AbsTol', 1e-4);
        end

        function integralImageMexRejectsMoreThan3D(testCase)
            % the mex copied the dimensions into an array of 3 (a 4-D input
            % wrote past it)
            testCase.assumeMexAvailable('integral_image_3d_mex');
            A = single(rand(10, 12, 8, 2));
            testCase.verifyThat(@() integral_image_3d_mex(A, [3, 3, 3]), ...
                matlab.unittest.constraints.Throws('integral_image_3d:inputError', 'WhenNargoutIs', 1));
        end
    end
end
