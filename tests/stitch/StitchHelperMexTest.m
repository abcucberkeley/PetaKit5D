classdef StitchHelperMexTest < pktest.TestCase
    % stitching helpers that run a mex with MATLAB code as the fallback and
    % that the stitching tests don't reach: the crop of a long registration
    % region to its brightest part (crop_subregion_by_intensity, crop3d_mex)
    % and the check that the major tiles of a batch have data at every voxel
    % (check_major_tile_valid, any_4th_dim_mex). pktest.TestCase fails a
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
            % (one voxel: the mex updates the bounding box from several
            % threads without synchronization)
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
    end
end
