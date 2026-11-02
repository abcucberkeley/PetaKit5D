classdef DeskewRotateFrame3DTest < pktest.TestCase
    % deskewRotateFrame3D (combined deskew/rotation in memory) on a grid of
    % synthetic beads with known raw positions.
    %
    % The raw -> output map is fitted from the bead centroids and compared
    % with the acquisition geometry. For sample scan, one raw x pixel lies
    % along the light sheet and maps to (cos(theta), -/+sin(theta)) in
    % output (x, z); one raw z slice is a stage step of dz and maps to
    % +/-dz/xyPixelSize in output x; y is unchanged. Output voxels are
    % isotropic with size xyPixelSize. Fitting many beads averages out the
    % +/-0.3 voxel centroid jitter from interpolation.

    properties (Constant)
        SkewAngle = 32.45
        XYPixelSize = 0.108
        RawSize = [64, 96, 80]
    end

    properties (TestParameter)
        % 0.3 exceeds the x step threshold and uses skewed space interpolation
        dz = struct('dz02', 0.2, 'dz03', 0.3)
        reverse = struct('reverse', true, 'forward', false)
    end

    methods (Test, TestTags = {'Unit'})
        function outputSizeMatchesGeometry(testCase, dz, reverse)
            vol = zeros(testCase.RawSize, 'single');
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, 'reverse', reverse);
            testCase.verifyEqual(size(out), testCase.sampleScanOutSize(dz));
        end

        function geometryMatchesAcquisition(testCase, dz, reverse)
            [raw, vol] = testCase.beadGrid();
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, 'reverse', reverse);

            [expected, expectedPos] = testCase.expectedSampleScanMap(dz, reverse, raw);
            [actual, residual] = testCase.fitMap(out, raw, expectedPos);
            testCase.verifyEqual(actual(:, 1 : 3), expected(:, 1 : 3), 'AbsTol', 0.005 * dz / testCase.XYPixelSize, ...
                sprintf('Fitted raw -> output map does not match the acquisition geometry.\nFitted:\n%s\nExpected:\n%s', ...
                mat2str(actual(:, 1 : 3), 4), mat2str(expected(:, 1 : 3), 4)));
            testCase.verifyLessThan(residual, 0.6, 'Bead positions are not consistent with an affine map.');
        end

        function canvasOriginMatchesConvention(testCase, dz, reverse)
            % the output canvas origin matters downstream (e.g., stitching
            % in DSR space): the fitted map must put the bead grid where the
            % centered deskew/rotation convention puts it. The fast (mex)
            % path used to be one shear step off for reverse=false.
            [raw, vol] = testCase.beadGrid();
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, 'reverse', reverse);

            [~, expectedPos] = testCase.expectedSampleScanMap(dz, reverse, raw);
            actual = testCase.fitMap(out, raw, expectedPos);
            testCase.verifyEqual([raw, ones(size(raw, 1), 1)] * actual', expectedPos, 'AbsTol', 0.25);
        end

        function totalIntensityScalesWithVoxelVolume(testCase, dz, reverse)
            % raw voxels are dz * sin(theta) / xyPixelSize times larger than
            % output voxels
            [~, vol] = testCase.beadGrid();
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, 'reverse', reverse);
            zAniso = sind(testCase.SkewAngle) * dz / testCase.XYPixelSize;
            testCase.verifyEqual(double(sum(out, 'all') / sum(vol, 'all')), zAniso, 'RelTol', 0.02);
        end

        function uint16InputWithSave16bitGivesUint16(testCase)
            [~, vol] = testCase.beadGrid();
            vol = uint16(vol + 100);
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, 0.2, testCase.XYPixelSize, 'reverse', true, 'save16bit', true);
            testCase.verifyClass(out, 'uint16');
            outSingle = deskewRotateFrame3D(single(vol), testCase.SkewAngle, 0.2, testCase.XYPixelSize, 'reverse', true);
            testCase.verifyClass(outSingle, 'single');
            % +/-1 for rounding differences between the uint16 and single paths
            testCase.verifyLessThanOrEqual(max(abs(double(out(:)) - double(round(outSingle(:))))), 1);
        end

        function resampleFactorScalesOutputSize(testCase)
            vol = zeros(testCase.RawSize, 'single');
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, 0.2, testCase.XYPixelSize, 'reverse', true, ...
                'resampleFactor', [2, 2, 2]);
            testCase.verifyEqual(size(out), round(testCase.sampleScanOutSize(0.2) ./ [2, 2, 2]));
        end

        function bboxCropsTheFullOutput(testCase)
            [~, vol] = testCase.beadGrid();
            full = deskewRotateFrame3D(vol, testCase.SkewAngle, 0.2, testCase.XYPixelSize, 'reverse', true);
            % a y extent that isn't a multiple of 8 also runs the scalar tail
            % of the SIMD loops (see DeskewMexSourceTest)
            bbox = [9, 21, 5, 53, 150, 40];
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, 0.2, testCase.XYPixelSize, 'reverse', true, 'bbox', bbox);
            testCase.verifyImagesClose(out, full(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6)), 'relTol', 1e-6);
        end

        function separateMexStepsReplaceTheCombinedMex(testCase)
            % with skewed space interpolation (dz 0.3), where the combined mex
            % fails (here hidden by a stub that errors), its two steps run in
            % their own mex functions (the test fails if they fall back too)
            [~, vol] = testCase.beadGrid();
            bbox = [9, 21, 5, 53, 150, 40];
            combined = deskewRotateFrame3D(vol, testCase.SkewAngle, 0.3, testCase.XYPixelSize, 'reverse', true, 'bbox', bbox);
            testCase.applyFixture(pktest.FunctionStubFixture('skewed_space_interp_volume_deskew_rotate_warp_mex'));
            testCase.expectMexFallback('skewed_space_interp_volume_deskew_rotate_warp_mex');
            separate = deskewRotateFrame3D(vol, testCase.SkewAngle, 0.3, testCase.XYPixelSize, 'reverse', true, 'bbox', bbox);
            testCase.verifyEqual({mex_fallback('-log').mex}, {'skewed_space_interp_volume_deskew_rotate_warp_mex'});
            testCase.verifyImagesClose(separate, combined, 'relTol', 1e-6);
        end

        function bboxPastTheInputRowsMatchesImwarp(testCase, dz)
            % an output region whose y range goes past the input's rows is 0
            % there, as imwarp gives (the mex functions hidden by stubs that
            % error); the warp mex used to read past the rows
            [~, vol] = testCase.beadGrid();
            bbox = [9, 21, 5, testCase.RawSize(1) + 12, 150, 40];
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, 'reverse', true, 'bbox', bbox);
            testCase.applyFixture(pktest.FunctionStubFixture('volume_deskew_rotate_warp_mex'));
            testCase.applyFixture(pktest.FunctionStubFixture('skewed_space_interp_volume_deskew_rotate_warp_mex'));
            testCase.expectMexFallback('volume_deskew_rotate_warp_mex', 'skewed_space_interp_volume_deskew_rotate_warp_mex');
            ref = deskewRotateFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, 'reverse', true, 'bbox', bbox);

            nIn = testCase.RawSize(1) - bbox(1) + 1;
            testCase.verifyEqual(size(out), bbox(4 : 6) - bbox(1 : 3) + 1);
            testCase.verifyEqual(nnz(ref(nIn + 1 : end, :, :)), 0);
            testCase.verifyEqual(nnz(out(nIn + 1 : end, :, :)), 0, 'Values past the input''s rows.');
            testCase.verifyImagesClose(out, ref, 'relTol', 1e-6);
        end

        function forwardInterpolationEndsAtTheLastRawSlice(testCase)
            % the separate skewed space interpolation (with resampling, on the
            % GPU, or where the combined mex fails) at dz 0.6 (step 1/3,
            % 0.33333334 in single precision): the last interpolated slice
            % passes the last raw slice by a rounding error, and forward scans
            % interpolated it with a slice past the input instead of copying it
            testCase.assumeMexAvailable('skewed_space_interp_defined_stepsize_mex');
            xstep = cosd(testCase.SkewAngle) * 0.6 / testCase.XYPixelSize;
            v = single(rand(16, 40, 30) * 1000);
            for save16bit = [true, false]
                for input = {v, uint16(v)}
                    out = skewed_space_interp_defined_stepsize_mex(input{1}, xstep, 1 / 3, false, save16bit);
                    testCase.verifyEqual(size(out, 3), 88);
                    testCase.verifyEqual(out(:, :, end), cast(input{1}(:, :, end), class(out)), ...
                        sprintf('%s input, save16bit %d', class(input{1}), save16bit));
                end
            end
        end

        function inputsNarrowerThanTheShiftRunTheMatlabCode(testCase)
            % the rows of an interpolated slice come from rows up to the x step
            % apart: the mex functions that interpolate read past narrower
            % inputs; they reject them, so the MATLAB interpolation runs
            testCase.assumeMexAvailable('skewed_space_interp_volume_deskew_rotate_warp_mex', ...
                'skewed_space_interp_defined_stepsize_mex');
            rejected = matlab.unittest.constraints.Throws('volume_warp:inputError', 'WhenNargoutIs', 1);
            v = single(rand(8, 6, 12));
            testCase.verifyThat(@() skewed_space_interp_volume_deskew_rotate_warp_mex(v, 4.69, 1 / 3, true, eye(4), ...
                [1, 1, 1, 8, 5, 5], false), rejected);
            testCase.verifyThat(@() skewed_space_interp_defined_stepsize_mex(v, 4.69, 1 / 3, true, false), rejected);

            % dz 1.5 (x step 11.7), 13 columns
            v = single(rand(8, 13, 12) * 1000);
            testCase.expectMexFallback('skewed_space_interp_volume_deskew_rotate_warp_mex', 'skewed_space_interp_defined_stepsize_mex');
            out = deskewRotateFrame3D(v, testCase.SkewAngle, 1.5, testCase.XYPixelSize, 'reverse', true);
            testCase.applyFixture(pktest.FunctionStubFixture('skewed_space_interp_volume_deskew_rotate_warp_mex'));
            testCase.applyFixture(pktest.FunctionStubFixture('skewed_space_interp_defined_stepsize_mex'));
            ref = deskewRotateFrame3D(v, testCase.SkewAngle, 1.5, testCase.XYPixelSize, 'reverse', true);
            testCase.verifyGreaterThan(nnz(ref), 0);
            testCase.verifyEqual(out, ref);
        end

        function singleSliceOrColumnReadsNothingOutside(testCase)
            % output points interpolate from two slices (columns); for a point
            % exactly on the only slice (column) the warp mex stepped to the one
            % before the array. They are outside now, as all points are there.
            testCase.assumeMexAvailable('volume_deskew_rotate_warp_mex', 'skewed_space_interp_volume_deskew_rotate_warp_mex');
            v = single(rand(8, 40, 6) * 1000 + 1);
            % (the identity puts the first output slice on the first input slice)
            bbox = [1, 1, 1, 8, 5, 5];
            testCase.verifyEqual(volume_deskew_rotate_warp_mex(v(:, :, 1), eye(4), bbox, false), zeros(8, 5, 5, 'single'));
            testCase.verifyEqual(volume_deskew_rotate_warp_mex(v(:, 1, :), eye(4), bbox, false), zeros(8, 5, 5, 'single'));
            testCase.verifyEqual(volume_deskew_rotate_warp_mex(uint16(v(:, :, 1)), eye(4), bbox, true), zeros(8, 5, 5, 'uint16'));
            testCase.verifyEqual(skewed_space_interp_volume_deskew_rotate_warp_mex(v(:, :, 1), 4.69, 1 / 3, true, eye(4), bbox, false), ...
                zeros(8, 5, 5, 'single'));
            % with two slices the points on them have values
            testCase.verifyGreaterThan(nnz(volume_deskew_rotate_warp_mex(v(:, :, 1 : 2), eye(4), bbox, false)), 0);
        end

        function combinedMatchesSeparateDeskewThenRotate(testCase)
            % the combined transform must equal deskew followed by rotation
            % (up to the extra interpolation of the separate path)
            dz = 0.2;
            [raw, vol] = testCase.beadGrid();
            combined = deskewRotateFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, 'reverse', true);

            theta = testCase.SkewAngle * pi / 180;
            ds = deskewFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, true, 'Crop', false);
            separate = rotateFrame3D(ds, testCase.SkewAngle, sin(theta) * dz / testCase.XYPixelSize, true, ...
                'outSize', size(combined));

            [~, expectedPos] = testCase.expectedSampleScanMap(dz, true, raw);
            testCase.verifyEqual(testCase.fitMap(separate, raw, expectedPos), testCase.fitMap(combined, raw, expectedPos), ...
                'AbsTol', 0.1);
            % the separate path interpolates twice, which blurs the beads
            testCase.verifyImagesClose(separate, combined, 'relTol', 0.1);
        end

        function objectiveScanGeometryWithCubicInterpolation(testCase, reverse)
            % objective scan is a pure rotation (no deskew); the cubic
            % (imwarp) path doesn't have the mex issue below
            dz = 0.2;
            q = dz / testCase.XYPixelSize;
            [raw, vol] = testCase.beadGrid();
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, 'objectiveScan', true, ...
                'reverse', reverse, 'interpMethod', 'cubic');

            sz = testCase.RawSize;
            theta = testCase.SkewAngle * pi / 180;
            outSize = round([sz(1), sz(2) * cos(theta) + sz(3) * q * sin(theta), sz(3) * q * cos(theta) + sz(2) * sin(theta)]);
            testCase.verifyEqual(size(out), outSize);

            if reverse
                theta = -theta;
            end
            expected = [1, 0, 0, 0; 0, cos(theta), sin(theta) * q, 0; 0, -sin(theta), cos(theta) * q, 0];
            expected(:, 4) = [0; (outSize(2 : 3)' + 1) / 2] - expected(:, 1 : 3) * [0; (sz(2 : 3)' + 1) / 2];
            expectedPos = [raw, ones(size(raw, 1), 1)] * expected';
            actual = testCase.fitMap(out, raw, expectedPos);
            testCase.verifyEqual(actual(:, 1 : 3), expected(:, 1 : 3), 'AbsTol', 0.005 * q);
            testCase.verifyEqual([raw, ones(size(raw, 1), 1)] * actual', expectedPos, 'AbsTol', 0.25);
            testCase.verifyEqual(double(sum(out, 'all') / sum(vol, 'all')), q, 'RelTol', 0.02);
        end

        function objectiveScanPreservesContent(testCase)
            % objective scan used to go through volume_deskew_rotate_warp_mex,
            % which assumes a sample scan and kept ~6.5% of the content
            dz = 0.2;
            [~, vol] = testCase.beadGrid();
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, 'objectiveScan', true);
            testCase.verifyEqual(double(sum(out, 'all') / sum(vol, 'all')), dz / testCase.XYPixelSize, 'RelTol', 0.02, ...
                'Objective scan deskew/rotation lost the image content.');
        end

        function save16bitRoundsToNearest(testCase)
            % the AVX uint16 kernels used to add 0.5 and then round, which is ceil
            [~, vol] = testCase.beadGrid();
            vol = uint16(vol + 100);
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, 0.2, testCase.XYPixelSize, 'reverse', true, 'save16bit', true);
            ref = deskewRotateFrame3D(single(vol), testCase.SkewAngle, 0.2, testCase.XYPixelSize, 'reverse', true);
            % values within 1e-3 of a tie can round either way
            notTie = abs(ref - floor(ref) - 0.5) > 1e-3;
            mismatched = nnz(out(notTie) ~= uint16(round(ref(notTie)))) / nnz(notTie);
            testCase.verifyLessThan(mismatched, 1e-3, ...
                sprintf('%.1f%% of the uint16 voxels are not rounded to the nearest integer.', 100 * mismatched));
        end
    end

    methods (Access = private)
        function [raw, vol] = beadGrid(testCase)
            % 15 beads on a raw x-z grid, each at its own y so they stay
            % separated after the transform. y is scrambled so it is not a
            % linear function of x and z (that would make the fit singular).
            [xr, zr] = ndgrid([30, 48, 66], [12, 26, 40, 54, 68]);
            raw = [8 + 3.4 * mod((0 : 14)' * 7, 15), xr(:), zr(:)];
            vol = pktest.beads(testCase.RawSize, raw, 2, 1000);
        end

        function [outSize] = sampleScanOutSize(testCase, dz)
            sz = testCase.RawSize;
            theta = testCase.SkewAngle * pi / 180;
            outSize = round([sz(1), (sz(2) - 1) * cos(theta) + (sz(3) - 1) * dz / testCase.XYPixelSize, ...
                (sz(2) - 1) * sin(theta) - 4]);
        end

        function [map, pos] = expectedSampleScanMap(testCase, dz, reverse, raw)
            % affine map [y; x; z]_out = map * [y; x; z; 1]_raw of the
            % centered deskew/rotation convention used by deskewRotateFrame3D
            sz = testCase.RawSize;
            theta = testCase.SkewAngle * pi / 180;
            q = dz / testCase.XYPixelSize;
            dx = cos(theta) * q;
            f = ceil((sz(3) - 1) * dx) - (sz(3) - 1) * dx;
            c = (testCase.sampleScanOutSize(dz) + 1) / 2;
            s = 1 - 2 * reverse;

            map = [1, 0, 0, 0; 0, cos(theta), s * q, 0; 0, -s * sin(theta), 0, 0];
            map(:, 4) = [0; c(2) - s * cos(theta) * f / 2; c(3) + sin(theta) * f / 2] - map(:, 1 : 3) * [0; (sz(2 : 3)' + 1) / 2];
            pos = [raw, ones(size(raw, 1), 1)] * map';
        end

        function [map, residual] = fitMap(~, out, raw, approxPos)
            % least squares affine map from bead centroids measured near
            % their approximate output positions
            c = zeros(size(raw));
            for i = 1 : size(raw, 1)
                c(i, :) = pktest.centroid(out, approxPos(i, :), 6);
            end
            A = [raw, ones(size(raw, 1), 1)];
            map = (A \ c)';
            residual = max(abs(A * map' - c), [], 'all');
        end
    end
end
