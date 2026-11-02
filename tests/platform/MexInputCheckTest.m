classdef MexInputCheckTest < pktest.TestCase
    % the input checks of the PetaKit5D mex functions: inputs they used to
    % read or write out of bounds with (more dimensions than they copy into
    % fixed arrays, sparse arrays, whose memory has only the nonzero values,
    % bounding boxes that aren't double, start before 1 or end past the
    % output, NaN coordinates, zero or negative steps, empty inputs, too few
    % outputs) raise an error with the mex's identifier, and the valid edge
    % cases next to them still work as before. Inputs that went out of bounds
    % crash or corrupt MATLAB with the old binaries, so most of these tests
    % only pass or fail on binaries with the checks. An int64 array with the
    % bits of a valid double array (asInt64) is read as that array without
    % the class checks (also by the old binaries, safely).

    methods (Test, TestTags = {'Unit'})
        function crop3dChecksItsInputs(testCase)
            testCase.assumeMexAvailable('crop3d_mex');
            im = reshape(uint16(1 : 6 * 5 * 4), 6, 5, 4);
            id = 'crop:inputError';
            testCase.verifyRejects(@() crop3d_mex(sparse(6, 5), [1, 1, 1, 2, 2, 1]), 0, id, 'sparse input');
            testCase.verifyRejects(@() crop3d_mex(im, single([1, 1, 1, 2, 2, 2])), 0, id, 'single bbox');
            testCase.verifyRejects(@() crop3d_mex(im, int32([1, 1, 1, 2, 2, 2])), 0, id, 'int32 bbox');
            testCase.verifyRejects(@() crop3d_mex(im, testCase.asInt64([1, 1, 1, 2, 2, 2])), 0, id, 'int64 bbox');
            testCase.verifyRejects(@() crop3d_mex(im, zeros(0, 6)), 0, id, 'empty bbox');
            testCase.verifyRejects(@() crop3d_mex(im, [-1, 1, 1, 2, 2, 2]), 0, id, 'start -1');
            testCase.verifyRejects(@() crop3d_mex(im, [NaN, 1, 1, 2, 2, 2]), 0, id, 'NaN start');
            testCase.verifyRejects(@() crop3d_mex(im, [1, 1, 1, Inf, 2, 2]), 0, id, 'Inf end');
            testCase.verifyRejects(@() crop3d_mex(im, [4, 1, 1, 2, 2, 2]), 0, id, 'end < start - 1');
            % empty ranges and non-integer bounds (truncated) as before
            testCase.verifyEqual(crop3d_mex(im, [4, 1, 1, 3, 5, 4]), im(4 : 3, :, :));
            testCase.verifyEqual(crop3d_mex(im, [1.5, 1, 2, 6.9, 5, 4]), im(1 : 6, :, 2 : 4));
        end

        function crop4dChecksItsInputs(testCase)
            testCase.assumeMexAvailable('crop4d_mex');
            im = reshape(uint16(1 : 4 * 3 * 2 * 2), 4, 3, 2, 2);
            id = 'crop:inputError';
            testCase.verifyRejects(@() crop4d_mex(zeros(2, 2, 2, 2, 2, 'uint16'), [1, 1, 1, 1, 2, 2, 2, 2]), 1, id, '5-D input');
            testCase.verifyRejects(@() crop4d_mex(sparse(4, 3), [1, 1, 1, 1, 2, 2, 1, 1]), 1, id, 'sparse input');
            testCase.verifyRejects(@() crop4d_mex(im, single([1, 1, 1, 1, 2, 2, 2, 2])), 1, id, 'single bbox');
            testCase.verifyRejects(@() crop4d_mex(im, testCase.asInt64([1, 1, 1, 1, 2, 2, 2, 2])), 1, id, 'int64 bbox');
            testCase.verifyRejects(@() crop4d_mex(im, zeros(0, 8)), 1, id, 'empty bbox');
            testCase.verifyRejects(@() crop4d_mex(im, [NaN, 1, 1, 1, 2, 2, 2, 2]), 1, id, 'NaN start');
            testCase.verifyRejects(@() crop4d_mex(im, [1, 1, 1, 1, 2, NaN, 2, 2]), 1, id, 'NaN end');
            testCase.verifyEqual(crop4d_mex(im, [2, 1, 1, 2, 3.5, 3, 2, 2]), im(2 : 3, :, :, 2));
        end

        function indexing3dChecksItsInputs(testCase)
            testCase.assumeMexAvailable('indexing3d_mex');
            A = zeros(6, 5, 4, 'uint16');
            R = ones(2, 2, 2, 'uint16');
            id = 'indexing:inputError';
            testCase.verifyRejects(@() indexing3d_mex(zeros(6, 5, 4, 2, 'uint16'), [1, 1, 1, 2, 2, 2], R), 0, id, '4-D data');
            testCase.verifyRejects(@() indexing3d_mex(A, [1, 1, 1, 2, 2, 2], ones(2, 2, 2, 2, 'uint16')), 0, id, '4-D region');
            testCase.verifyRejects(@() indexing3d_mex(sparse(6, 5), [1, 1, 1, 2, 2, 1], sparse(2, 2)), 0, id, 'sparse data');
            testCase.verifyRejects(@() indexing3d_mex(zeros(6, 5), [1, 1, 1, 2, 2, 1], sparse(2, 2)), 0, id, 'sparse region');
            testCase.verifyRejects(@() indexing3d_mex(A, single([1, 1, 1, 2, 2, 2]), R), 0, id, 'single bbox');
            testCase.verifyRejects(@() indexing3d_mex(A, testCase.asInt64([1, 1, 1, 2, 2, 2]), R), 0, id, 'int64 bbox');
            testCase.verifyRejects(@() indexing3d_mex(A, zeros(0, 6), R), 0, id, 'empty bbox');
            % a start before 1 wrote before the array
            testCase.verifyRejects(@() indexing3d_mex(A, [-1, 1, 1, 0, 2, 2], R), 0, id, 'start -1');
            testCase.verifyRejects(@() indexing3d_mex(A, [NaN, 1, 1, 2, 2, 2], R), 0, id, 'NaN start');
            testCase.verifyEqual(A, zeros(6, 5, 4, 'uint16'));
            % non-integer bounds are truncated as before
            indexing3d_mex(A, [2.5, 1, 1, 3.5, 2, 2], R);
            expected = zeros(6, 5, 4, 'uint16');
            expected(2 : 3, 1 : 2, 1 : 2) = 1;
            testCase.verifyEqual(A, expected);
        end

        function indexing4dChecksItsInputs(testCase)
            testCase.assumeMexAvailable('indexing4d_mex', 'indexing4d_crop_mex');
            A = zeros(4, 3, 2, 2, 'uint16');
            R = ones(2, 2, 1, 1, 'uint16');
            b = [1, 1, 1, 1, 2, 2, 1, 1];
            id = 'indexing:inputError';
            for f = {@indexing4d_mex, @indexing4d_crop_mex}
                fn = func2str(f{1});
                testCase.verifyRejects(@() f{1}(zeros(4, 3, 2, 2, 2, 'uint16'), b, R), 0, id, [fn, ' 5-D data']);
                testCase.verifyRejects(@() f{1}(A, b, ones(2, 2, 1, 1, 2, 'uint16')), 0, id, [fn, ' 5-D region']);
                testCase.verifyRejects(@() f{1}(sparse(4, 3), b, sparse(2, 2)), 0, id, [fn, ' sparse data']);
                testCase.verifyRejects(@() f{1}(A, single(b), R), 0, id, [fn, ' single bbox']);
                testCase.verifyRejects(@() f{1}(A, testCase.asInt64(b), R), 0, id, [fn, ' int64 bbox']);
                testCase.verifyRejects(@() f{1}(A, zeros(0, 8), R), 0, id, [fn, ' empty bbox']);
                testCase.verifyRejects(@() f{1}(A, [NaN, 1, 1, 1, 2, 2, 1, 1], R), 0, id, [fn, ' NaN start']);
            end
            testCase.verifyRejects(@() indexing4d_mex(A, [-1, 1, 1, 1, 0, 2, 1, 1], R), 0, id, 'indexing4d_mex start -1');
            testCase.verifyRejects(@() indexing4d_crop_mex(A, b, R, single([1, 1, 1, 1, 2, 2, 1, 1])), 0, id, 'single region bbox');
            testCase.verifyRejects(@() indexing4d_crop_mex(A, b, R, testCase.asInt64([1, 1, 1, 1, 2, 2, 1, 1])), 0, id, 'int64 region bbox');
            testCase.verifyRejects(@() indexing4d_crop_mex(A, b, R, [NaN, 1, 1, 1, 2, 2, 1, 1]), 0, id, 'NaN region bbox');
            testCase.verifyEqual(A, zeros(4, 3, 2, 2, 'uint16'));
        end

        function deskewRotateWarpChecksItsInputs(testCase)
            testCase.assumeMexAvailable('volume_deskew_rotate_warp_mex');
            vol = ones(16, 12, 10, 'single');
            tmat = eye(4);
            bbox = [1, 1, 1, 16, 12, 10];
            id = 'volume_warp:inputError';
            f = @volume_deskew_rotate_warp_mex;
            testCase.verifyRejects(@() f(ones(16, 12, 10, 2, 'single'), tmat, bbox, false), 1, id, '4-D input');
            testCase.verifyRejects(@() f(vol, single(tmat), bbox, false), 1, id, 'single transform');
            testCase.verifyRejects(@() f(vol, testCase.asInt64(tmat), bbox, false), 1, id, 'int64 transform');
            testCase.verifyRejects(@() f(vol, tmat, testCase.asInt64(bbox), false), 1, id, 'int64 bbox');
            % an entry the kernel doesn't use; NaN coordinates read outside the input
            t = tmat;
            t(4, 1) = NaN;
            testCase.verifyRejects(@() f(vol, t, bbox, false), 1, id, 'NaN in the transform');
            t = tmat;
            t(3, 4) = Inf;
            testCase.verifyRejects(@() f(vol, t, bbox, false), 1, id, 'Inf in the transform');
            testCase.verifyRejects(@() f(vol, tmat, single(bbox), false), 1, id, 'single bbox');
            testCase.verifyRejects(@() f(vol, tmat, zeros(0, 6), false), 1, id, 'empty bbox');
            testCase.verifyRejects(@() f(vol, tmat, [NaN, 1, 1, 16, 12, 10], false), 1, id, 'NaN bbox');
            testCase.verifyEqual(size(f(vol, tmat, bbox, false)), [16, 12, 10]);
        end

        function skewedSpaceInterpWarpChecksItsInputs(testCase)
            testCase.assumeMexAvailable('skewed_space_interp_volume_deskew_rotate_warp_mex');
            vol = ones(16, 12, 10, 'single');
            tmat = eye(4);
            bbox = [1, 1, 1, 16, 12, 10];
            id = 'volume_warp:inputError';
            f = @skewed_space_interp_volume_deskew_rotate_warp_mex;
            testCase.verifyRejects(@() f(ones(16, 12, 10, 2, 'single'), 2.5, 0.5, true, tmat, bbox, false), 1, id, '4-D input');
            % xstep 0 read a row past the slices
            testCase.verifyRejects(@() f(vol, 0, 0.5, true, tmat, bbox, false), 1, id, 'xstep 0');
            testCase.verifyRejects(@() f(vol, [2.5, 1], 0.5, true, tmat, bbox, false), 1, id, 'xstep vector');
            testCase.verifyRejects(@() f(vol, single(2.5), 0.5, true, tmat, bbox, false), 1, id, 'single xstep');
            testCase.verifyRejects(@() f(vol, testCase.asInt64(2.5), 0.5, true, tmat, bbox, false), 1, id, 'int64 xstep');
            testCase.verifyRejects(@() f(vol, 2.5, testCase.asInt64(0.5), true, tmat, bbox, false), 1, id, 'int64 step');
            testCase.verifyRejects(@() f(vol, 2.5, 0.5, true, testCase.asInt64(tmat), bbox, false), 1, id, 'int64 transform');
            testCase.verifyRejects(@() f(vol, 2.5, 0.5, true, tmat, testCase.asInt64(bbox), false), 1, id, 'int64 bbox');
            testCase.verifyRejects(@() f(vol, 2.5, 0, true, tmat, bbox, false), 1, id, 'step 0');
            testCase.verifyRejects(@() f(vol, 2.5, -0.5, true, tmat, bbox, false), 1, id, 'negative step');
            testCase.verifyRejects(@() f(vol, 2.5, NaN, true, tmat, bbox, false), 1, id, 'NaN step');
            testCase.verifyRejects(@() f(vol, 2.5, 1e-20, true, tmat, bbox, false), 1, id, 'tiny step');
            t = tmat;
            t(4, 1) = NaN;
            testCase.verifyRejects(@() f(vol, 2.5, 0.5, true, t, bbox, false), 1, id, 'NaN in the transform');
            testCase.verifyRejects(@() f(vol, 2.5, 0.5, true, tmat, [1, NaN, 1, 16, 12, 10], false), 1, id, 'NaN bbox');
            testCase.verifyRejects(@() f(vol, 2.5, 0.5, true, tmat, int32(bbox), false), 1, id, 'int32 bbox');
            % an input without slices gives zeros, as before (on x86)
            testCase.verifyEqual(f(ones(16, 12, 0, 'single'), 2.5, 0.5, true, tmat, bbox, false), zeros(16, 12, 10, 'single'));
        end

        function skewedSpaceInterpChecksItsInputs(testCase)
            testCase.assumeMexAvailable('skewed_space_interp_defined_stepsize_mex');
            vol = ones(16, 12, 10, 'single');
            id = 'volume_warp:inputError';
            f = @skewed_space_interp_defined_stepsize_mex;
            testCase.verifyRejects(@() f(ones(16, 12, 10, 2, 'single'), 2.5, 0.5, true, false), 1, id, '4-D input');
            testCase.verifyRejects(@() f(vol, 0, 0.5, true, false), 1, id, 'xstep 0');
            testCase.verifyRejects(@() f(vol, -0.5, 0.5, true, false), 1, id, 'xstep -0.5');
            testCase.verifyRejects(@() f(vol, [2.5, 1], 0.5, true, false), 1, id, 'xstep vector');
            testCase.verifyRejects(@() f(vol, testCase.asInt64(2.5), 0.5, true, false), 1, id, 'int64 xstep');
            testCase.verifyRejects(@() f(vol, 2.5, testCase.asInt64(0.5), true, false), 1, id, 'int64 step');
            testCase.verifyRejects(@() f(vol, 2.5, 0, true, false), 1, id, 'step 0');
            testCase.verifyRejects(@() f(vol, 2.5, Inf, true, false), 1, id, 'Inf step');
            testCase.verifyRejects(@() f(vol, 2.5, 1e-20, true, false), 1, id, 'tiny step');
            testCase.verifyEqual(size(f(vol, 2.5, 0.5, true, false)), [16, 12, 19]);
            testCase.verifyEqual(size(f(ones(16, 12, 0, 'single'), 2.5, 0.5, true, false)), [16, 12, 0]);
        end

        function stitchMexCheckTheirInputs(testCase)
            testCase.assumeMexAvailable('any_4th_dim_mex', 'feather_blending_3d_mex', 'feather_blending_3d_with_indexing_mex', ...
                'feather_distance_map_resize_3d_mex', 'integral_image_3d_mex');
            % any_4th_dim_mex
            testCase.verifyRejects(@() any_4th_dim_mex(ones(4, 3, 2, 2, 2, 'uint16'), true(1, 2)), 2, 'any_4th_dim_mex:inputError', '5-D tiles');
            testCase.verifyRejects(@() any_4th_dim_mex(ones(4, 3, 2, 2, 'uint16'), sparse([true, false])), 2, 'any_4th_dim_mex:inputError', ...
                'sparse tile flags');
            % feather_blending_3d_mex
            id = 'feather_blending:inputError';
            f = ones(4, 3, 2, 2, 'uint16');
            d = ones(4, 3, 2, 2, 'single');
            testCase.verifyRejects(@() feather_blending_3d_mex(ones(4, 3, 2, 2, 2, 'uint16'), ones(4, 3, 2, 2, 2, 'single')), 1, id, '5-D tiles');
            testCase.verifyRejects(@() feather_blending_3d_mex(f(:, :, :, 1), d), 1, id, 'distance map with more dimensions');
            testCase.verifyRejects(@() feather_blending_3d_mex(f, d(:, :, :, 1)), 1, id, 'distance map with fewer dimensions');
            testCase.verifyRejects(@() feather_blending_3d_mex(f, double(d)), 1, id, 'double distance map');
            testCase.verifyRejects(@() feather_blending_3d_mex(f, uint16(d)), 1, id, 'uint16 distance map');
            testCase.verifyEqual(feather_blending_3d_mex(f, d), ones(4, 3, 2, 'uint16'));
            % feather_blending_3d_with_indexing_mex: the region is written at bbox of the output
            out = zeros(10, 10, 10, 'uint16');
            testCase.verifyRejects(@() feather_blending_3d_with_indexing_mex(ones(4, 3, 2, 2, 2, 'uint16'), ones(4, 3, 2, 2, 2, 'single'), out, ...
                [1, 1, 1, 4, 3, 2]), 0, id, '5-D tiles');
            testCase.verifyRejects(@() feather_blending_3d_with_indexing_mex(f, d(:, :, :, 1), out, [1, 1, 1, 4, 3, 2]), 0, id, ...
                'distance map with fewer dimensions');
            testCase.verifyRejects(@() feather_blending_3d_with_indexing_mex(f, d, zeros(10, 10, 10, 2, 'uint16'), [1, 1, 1, 4, 3, 2]), 0, id, ...
                '4-D output');
            testCase.verifyRejects(@() feather_blending_3d_with_indexing_mex(f, d, out, single([1, 1, 1, 4, 3, 2])), 0, id, 'single bbox');
            testCase.verifyRejects(@() feather_blending_3d_with_indexing_mex(f, d, out, testCase.asInt64([1, 1, 1, 4, 3, 2])), 0, id, 'int64 bbox');
            % past the end of the rows (it wrote into the next rows)
            testCase.verifyRejects(@() feather_blending_3d_with_indexing_mex(f, d, out, [8, 1, 1, 11, 3, 2]), 0, id, 'bbox past the rows');
            testCase.verifyRejects(@() feather_blending_3d_with_indexing_mex(f, d, out, [NaN, 1, 1, 4, 3, 2]), 0, id, 'NaN bbox');
            testCase.verifyEqual(out, zeros(10, 10, 10, 'uint16'));
            feather_blending_3d_with_indexing_mex(f, d, out, [7, 8, 9, 10, 10, 10]);
            expected = zeros(10, 10, 10, 'uint16');
            expected(7 : 10, 8 : 10, 9 : 10) = 1;
            testCase.verifyEqual(out, expected);
            % feather_distance_map_resize_3d_mex
            testCase.verifyRejects(@() feather_distance_map_resize_3d_mex(ones(4, 3, 2, 2, 'single'), [4, 3, 2], 2), 1, id, '4-D distance map');
            testCase.verifyRejects(@() feather_distance_map_resize_3d_mex(ones(4, 3, 2, 'single'), single([4, 3, 2]), 2), 1, ...
                'distmap_resize:inputError', 'single size');
            testCase.verifyRejects(@() feather_distance_map_resize_3d_mex(ones(4, 3, 2, 'single'), testCase.asInt64([4, 3, 2]), 2), 1, ...
                'distmap_resize:inputError', 'int64 size');
            testCase.verifyRejects(@() feather_distance_map_resize_3d_mex(ones(4, 3, 2, 'single'), [4, NaN, 2], 2), 1, ...
                'distmap_resize:inputError', 'NaN size');
            testCase.verifyRejects(@() feather_distance_map_resize_3d_mex(ones(0, 3, 2, 'single'), [4, 3, 2], 2), 1, ...
                'distmap_resize:inputError', 'empty distance map');
            testCase.verifyEqual(feather_distance_map_resize_3d_mex(ones(4, 3, 2, 'single'), [8, 3, 2.5], 2), ones(8, 3, 2, 'single'));
            testCase.verifyEqual(size(feather_distance_map_resize_3d_mex(ones(0, 3, 2, 'single'), [0, 3, 2], 2)), [0, 3, 2]);
            % integral_image_3d_mex: the sums read the input up to size T - 1 back
            id = 'integral_image_3d:inputError';
            testCase.verifyRejects(@() integral_image_3d_mex(ones(4, 4, 4, 'single'), [5, 1, 1]), 1, id, 'input smaller than T');
            testCase.verifyRejects(@() integral_image_3d_mex(ones(4, 0, 4, 'single'), [1, 1, 1]), 1, id, 'empty input');
            testCase.verifyRejects(@() integral_image_3d_mex(ones(4, 4, 4, 'single'), single([2, 2, 2])), 1, id, 'single size T');
            testCase.verifyRejects(@() integral_image_3d_mex(ones(4, 4, 4, 'single'), testCase.asInt64([2, 2, 2])), 1, id, 'int64 size T');
            testCase.verifyRejects(@() integral_image_3d_mex(ones(4, 4, 4, 'single'), [2, NaN, 2]), 1, id, 'NaN size T');
            out = integral_image_3d_mex(ones(4, 4, 4, 'single'), [4, 4, 4]);
            testCase.verifyEqual(size(out), [7, 7, 7]);
            testCase.verifyEqual(out(4, 4, 4), single(64));
        end

        function flatFieldCorrectionChecksItsInputs(testCase)
            testCase.assumeMexAvailable('flat_field_correction_3d_mex');
            id = 'flat_field_correction_3d_mex:inputError';
            im = ones(4, 3, 2, 'uint16');
            ff = ones(4, 3, 'single');
            bg = zeros(4, 3, 'single');
            testCase.verifyRejects(@() flat_field_correction_3d_mex(ones(4, 3, 2, 2, 'uint16'), ff, bg, true), 1, id, '4-D data');
            % the identifier had spaces (MATLAB:mexapi:invalidErrorIdentifier)
            testCase.verifyRejects(@() flat_field_correction_3d_mex(im, ff, zeros(4, 3, 2, 'single'), true), 1, id, '3-D background');
            testCase.verifyEqual(flat_field_correction_3d_mex(im, ff, bg, true), im);
        end

        function mipMexCheckTheirInputs(testCase)
            testCase.assumeMexAvailable('max_pooling_3d_mex', 'range_max_3d_mex');
            id = 'pooling:inputError';
            im = reshape(single(1 : 4 * 3 * 2), 4, 3, 2);
            testCase.verifyRejects(@() max_pooling_3d_mex(im), 1, id, 'one input');
            testCase.verifyRejects(@() max_pooling_3d_mex(sparse(4, 3), [2, 2, 1]), 1, id, 'sparse input');
            testCase.verifyRejects(@() max_pooling_3d_mex(im, single([2, 2, 1])), 1, id, 'single pooling size');
            testCase.verifyRejects(@() max_pooling_3d_mex(im, testCase.asInt64([2, 2, 1])), 1, id, 'int64 pooling size');
            testCase.verifyRejects(@() max_pooling_3d_mex(im, [2, NaN, 1]), 1, id, 'NaN pooling size');
            testCase.verifyRejects(@() max_pooling_3d_mex(im, [0.5, 1, 1]), 1, id, 'pooling size 0.5');
            testCase.verifyEqual(max_pooling_3d_mex(im, [2, 3, 2.5]), single([22; 24]));
            id = 'range_max_3d:inputError';
            testCase.verifyRejects(@() range_max_3d_mex(sparse(4, 3), 1, [1, 3]), 1, id, 'sparse input');
            testCase.verifyRejects(@() range_max_3d_mex(im, 1, sparse([1, 3])), 1, id, 'sparse starts');
        end

        function omwBackProjectorWindowNeedsBothOutputs(testCase)
            testCase.assumeMexAvailable('omw_backprojector_mex');
            id = 'omw_backprojector_mex:invalidInput';
            E = [1, 0, 0; 0, 1, 0; 0, 0, 1];
            R = [1; 2; 3];
            % it set both outputs
            testCase.verifyRejects(@() omw_backprojector_mex('window', [8, 8, 8], E, R, [0.8, 1], 1), 1, id, 'one output');
            testCase.verifyRejects(@() omw_backprojector_mex('window', [8, 8, 8], sparse(E), R, [0.8, 1], 1), 2, id, 'sparse E');
            testCase.verifyRejects(@() omw_backprojector_mex('convhull', true(4, 4, 4), sparse(1), [], [], 1), 1, id, 'sparse indices');
            [mask, inds] = omw_backprojector_mex('window', [8, 8, 8], E, R, [0.8, 1], 1);
            testCase.verifyEqual(size(mask), [8, 8, 8]);
            testCase.verifyEqual(size(inds, 2), 1);
        end

        function pointDetectionMexCheckTheirInputs(testCase)
            testCase.assumeMexAvailable('replace_nan_with_value', 'adtest_mex', 'compute_gradient_hessian_mex');
            % replace_nan_with_value: the value is read as an element of the input's class
            id = 'data:inputError';
            testCase.verifyRejects(@() replace_nan_with_value(single([1, 2])), 1, id, 'one input');
            testCase.verifyRejects(@() replace_nan_with_value(sparse([1, 0, NaN]), 0), 1, id, 'sparse input');
            testCase.verifyRejects(@() replace_nan_with_value([1, NaN], single(0)), 1, id, 'single value for double input');
            testCase.verifyRejects(@() replace_nan_with_value([1, 2], []), 1, id, 'empty value');
            testCase.verifyRejects(@() replace_nan_with_value(single([1, NaN]), {0}), 1, id, 'cell value');
            testCase.verifyRejects(@() replace_nan_with_value(uint16([1, 2]), []), 1, 'data:TypeError', 'uint16 input');
            testCase.verifyEqual(replace_nan_with_value(single([1, NaN]), 0), single([1, 0]));
            testCase.verifyEqual(replace_nan_with_value([1, NaN], 5), [1, 5]);
            testCase.verifyEqual(replace_nan_with_value(single([1, NaN]), single(5)), single([1, 5]));
            % adtest_mex
            id = 'adtest_mex:inputError';
            x = sort(randn(50, 1));
            testCase.verifyRejects(@() adtest_mex(x, 0, 1, 0.05), 1, id, 'four inputs');
            testCase.verifyRejects(@() adtest_mex(single(x), 0, 1, 0.05, 1), 1, id, 'single data');
            testCase.verifyRejects(@() adtest_mex(testCase.asInt64(x), 0, 1, 0.05, 1), 1, id, 'int64 data');
            testCase.verifyRejects(@() adtest_mex(x, 0, 1, testCase.asInt64(0.05), 1), 1, id, 'int64 alpha');
            testCase.verifyRejects(@() adtest_mex(x, [], 1, 0.05, 1), 1, id, 'empty mean');
            testCase.verifyRejects(@() adtest_mex(x, 0, 1, single(0.05), 1), 1, id, 'single alpha');
            testCase.verifyClass(adtest_mex(x, 0, 1, 0.05, 1), 'double');
            % compute_gradient_hessian_mex: five vectors of the first's rows and two scalars
            id = 'compute_gradient_hessian_mex:inputError';
            v = ones(10, 1);
            testCase.verifyRejects(@() compute_gradient_hessian_mex(v, v, v, v, v, 1), 2, id, 'six inputs');
            testCase.verifyRejects(@() compute_gradient_hessian_mex(v, v(1 : 5), v, v, v, 1, 1), 2, id, 'shorter vector');
            testCase.verifyRejects(@() compute_gradient_hessian_mex(v, v, v, v, single(v), 1, 1), 2, id, 'single vector');
            testCase.verifyRejects(@() compute_gradient_hessian_mex(v, v, testCase.asInt64(v), v, v, 1, 1), 2, id, 'int64 vector');
            testCase.verifyRejects(@() compute_gradient_hessian_mex(v, v, v, v, v, testCase.asInt64(1), 1), 2, id, 'int64 sigma');
            testCase.verifyRejects(@() compute_gradient_hessian_mex(v, v, v, v, v, [], 1), 2, id, 'empty sigma');
            testCase.verifyRejects(@() compute_gradient_hessian_mex(ones(10, 0), v, v, v, v, 1, 1), 2, id, 'empty first vector');
            [g, h] = compute_gradient_hessian_mex(v, v, v, v, v, 1, 1);
            testCase.verifyEqual(size(g), [3, 1]);
            testCase.verifyEqual(size(h), [3, 3]);
        end
    end

    methods (Access = private)
        function [] = verifyRejects(testCase, fcn, nOut, id, what)
            % fcn errors with identifier id when called with nOut outputs
            testCase.verifyThat(fcn, matlab.unittest.constraints.Throws(id, 'WhenNargoutIs', nOut), what);
        end
    end

    methods (Static, Access = private)
        function [y] = asInt64(x)
            % an int64 array of x's size with the bits of double(x)
            y = reshape(typecast(double(x(:)), 'int64'), size(x));
        end
    end
end
