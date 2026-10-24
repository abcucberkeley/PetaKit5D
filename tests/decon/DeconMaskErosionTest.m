classdef DeconMaskErosionTest < pktest.TestCase
    % the edge erosion of the decon masks (decon_mask_edge_erosion). Masks
    % that aren't full (e.g., skewed data with zero borders) are eroded on the
    % CPU by decon_mask_edge_erosion_mex, which applies the decomposition of
    % strel('disk', EdgeErosion - 1) and the line along z as imerode does,
    % with the same masks as the previous imerode code
    % (tests/decon/private/decon_mask_edge_erosion_before.m).

    properties (TestParameter)
        shape = {'random', 'skewed', 'thin', 'single', 'singleFalse', 'blobs'}
        maskSize = struct('odd3D', [37, 41, 23], 'small', [5, 5, 5], 'twoSlices', [17, 19, 2], ...
            'oneRow', [1, 30, 20], 'oneColumn', [30, 1, 20], 'image2D', [45, 47])
    end

    methods (Test, TestTags = {'Unit'})
        function masksEqualImerode(testCase, shape, maskSize)
            % disks of radius 0 to 40 (no decomposition up to 2, then periodic
            % lines with morphop's padding, with or without extra lines)
            % (full masks thinner than the erosion error in both)
            testCase.assumeMexAvailable('decon_mask_edge_erosion_mex');
            A = testCase.makeMask(shape, maskSize);
            for E = 1 : 41
                testCase.verifyEqual(testCase.outcome(@() decon_mask_edge_erosion(A, E)), ...
                    testCase.outcome(@() decon_mask_edge_erosion_before(A, E)), sprintf('EdgeErosion %d', E));
            end
        end

        function demoSettingOnASkewedVolume(testCase)
            % EdgeErosion 20 as in the large-scale demo, zero borders as in
            % skewed data, more slices than threads and a few isolated holes
            testCase.assumeMexAvailable('decon_mask_edge_erosion_mex');
            sz = [150, 200, 90];
            A = testCase.makeMask('skewed', sz);
            A(sub2ind(sz, [70, 20, 100], [70, 90, 110], [40, 60, 75])) = false;
            expected = decon_mask_edge_erosion_before(A, 20);
            testCase.assertTrue(any(expected(:)) && ~all(expected(:)));
            testCase.verifyEqual(decon_mask_edge_erosion(A, 20), expected);
            % numeric input as decon passes it (e.g., the raw data)
            raw = uint16(A) .* uint16(randi(1000, size(A)));
            testCase.verifyEqual(decon_mask_edge_erosion(raw, 20), expected);
        end

        function mexEqualsImerodeForOtherNeighborhoods(testCase)
            % sequences of random (asymmetric, even-sized) neighborhoods, with
            % the padding morphop uses for them, and single neighborhoods
            testCase.assumeMexAvailable('decon_mask_edge_erosion_mex');
            for trial = 1 : 40
                ns = randi(3);
                nh = cell(1, ns);
                for k = 1 : ns
                    e = rand(randi(6, 1, 2)) > 0.5;
                    e(randi(numel(e))) = true;
                    nh{k} = e;
                end
                A = rand([randi(30), randi(30), randi(4)]) > 0.2;
                se = strel(nh{1});
                for k = 2 : ns
                    se(k) = strel(nh{k});
                end
                padUL = [0, 0];
                padLR = [0, 0];
                if ns > 1 && ~isdecompositionorthogonal(se)
                    [padBefore, padAfter] = getpadsize(se);
                    padUL = padAfter(1 : 2);
                    padLR = padBefore(1 : 2);
                end
                for nt = [1, 3]
                    testCase.verifyEqual(decon_mask_edge_erosion_mex(A, nh, padUL, padLR, 0, false, nt), imerode(A, se), ...
                        sprintf('Trial %d, %d threads.', trial, nt));
                end
            end
            % a line along z as imerode with true(2 * r + 1, 1) along the first dimension
            A = rand(20, 30, 40) > 0.1;
            expected = permute(imerode(permute(A, [3, 1, 2]), true(9, 1)), [2, 3, 1]);
            testCase.verifyEqual(decon_mask_edge_erosion_mex(A, {true}, [0, 0], [0, 0], 4, false, 2), expected);
        end

        function resultDoesNotDependOnTheThreadCount(testCase)
            % whole slices per thread, or the columns of each slice split among
            % the threads when there are few slices
            testCase.assumeMexAvailable('decon_mask_edge_erosion_mex');
            nt = maxNumCompThreads;
            testCase.addTeardown(@() maxNumCompThreads(nt));
            for sz = {[60, 70, 40], [60, 70, 3], [80, 90]}
                A = testCase.makeMask('skewed', sz{1}) & rand(sz{1}) > 0.001;
                expected = decon_mask_edge_erosion_before(A, 9);
                for threads = [1, 2, 5]
                    maxNumCompThreads(threads);
                    testCase.verifyEqual(decon_mask_edge_erosion(A, 9), expected, sprintf('%d threads, size %s', threads, mat2str(sz{1})));
                end
            end
        end

        function irregularMasksDoNotUseImerode(testCase)
            % (imerode on large masks took a fifth of the large-file decon)
            testCase.assumeMexAvailable('decon_mask_edge_erosion_mex');
            A = testCase.makeMask('skewed', [40, 50, 30]);
            expected = decon_mask_edge_erosion_before(A, 7);
            testCase.applyFixture(pktest.FunctionStubFixture('imerode', false));
            testCase.verifyEqual(decon_mask_edge_erosion(A, 7), expected);
        end

        function withoutTheMexImerodeRuns(testCase)
            % a .m file of the same name first on the path hides the mex, as on
            % a platform without it
            A = testCase.makeMask('skewed', [40, 50, 30]);
            expected = decon_mask_edge_erosion_before(A, 7);
            testCase.applyFixture(pktest.FunctionStubFixture('decon_mask_edge_erosion_mex'));
            testCase.verifyEqual(exist('decon_mask_edge_erosion_mex', 'file'), 2);
            testCase.verifyEqual(decon_mask_edge_erosion(A, 7), expected);
        end

        function fullMasksAndOtherErosionValues(testCase)
            % full masks get their edges set directly; erosion values the mex
            % doesn't take run the previous code (with its errors)
            A = true(20, 22, 12);
            for E = [1, 3, 6]
                testCase.verifyEqual(decon_mask_edge_erosion(A, E), decon_mask_edge_erosion_before(A, E));
            end
            testCase.verifyEqual(decon_mask_edge_erosion(uint8(A) * 7, 2), decon_mask_edge_erosion_before(A, 2));
            S = sparse(testCase.makeMask('skewed', [20, 22]));
            testCase.verifyEqual(testCase.outcome(@() decon_mask_edge_erosion(S, 3)), ...
                testCase.outcome(@() decon_mask_edge_erosion_before(S, 3)), 'sparse mask');
            B = testCase.makeMask('skewed', [20, 22, 12]);
            for E = {0, single(3), 2.5}
                testCase.verifyEqual(testCase.outcome(@() decon_mask_edge_erosion(B, E{1})), ...
                    testCase.outcome(@() decon_mask_edge_erosion_before(B, E{1})), sprintf('EdgeErosion %g', E{1}));
            end
        end

        function mexRejectsInvalidInputs(testCase)
            testCase.assumeMexAvailable('decon_mask_edge_erosion_mex');
            A = rand(10, 12, 5) > 0.5;
            id = 'decon_mask_edge_erosion_mex:inputError';
            testCase.verifyError(@() decon_mask_edge_erosion_mex(A, {true}, [0, 0], [0, 0], 1, true), id);
            testCase.verifyError(@() decon_mask_edge_erosion_mex(double(A), {true}, [0, 0], [0, 0], 1, true, 1), id);
            testCase.verifyError(@() decon_mask_edge_erosion_mex(true(3, 3, 3, 2), {true}, [0, 0], [0, 0], 1, true, 1), id);
            testCase.verifyError(@() decon_mask_edge_erosion_mex(sparse(A(:, :, 1)), {true}, [0, 0], [0, 0], 1, true, 1), id);
            testCase.verifyError(@() decon_mask_edge_erosion_mex(A, true, [0, 0], [0, 0], 1, true, 1), id);
            testCase.verifyError(@() decon_mask_edge_erosion_mex(A, {false(3)}, [0, 0], [0, 0], 1, true, 1), id);
            testCase.verifyError(@() decon_mask_edge_erosion_mex(A, {true(2, 2, 2)}, [0, 0], [0, 0], 1, true, 1), id);
            testCase.verifyError(@() decon_mask_edge_erosion_mex(A, {true}, [0, -1], [0, 0], 1, true, 1), id);
            testCase.verifyError(@() decon_mask_edge_erosion_mex(A, {true}, [0, 0], 0, 1, true, 1), id);
            testCase.verifyError(@() decon_mask_edge_erosion_mex(A, {true}, [0, 0], [0, 0], -1, true, 1), id);
            testCase.verifyError(@() decon_mask_edge_erosion_mex(A, {true}, [0, 0], [0, 0], 1.5, true, 1), id);
            % numeric neighborhoods count their nonzero elements
            testCase.verifyEqual(decon_mask_edge_erosion_mex(A, {[0, 2, 1]}, [0, 0], [0, 0], 0, false, 1), ...
                imerode(A, strel([false, true, true])));
            testCase.verifyEqual(decon_mask_edge_erosion_mex(false(0, 3, 2), {true}, [0, 0], [0, 0], 1, true, 1), false(0, 3, 2));
        end
    end

    methods (Test, TestTags = {'Unit', 'GPU'})
        function gpuArrayMasksEqualTheCPUMasks(testCase)
            % the GPU path erodes gpuArray masks with imerode on the GPU
            testCase.assumeGPUAvailable();
            testCase.assumeMexAvailable('decon_mask_edge_erosion_mex');
            A = testCase.makeMask('skewed', [61, 57, 33]) & rand([61, 57, 33]) > 0.001;
            for E = [1, 2, 3, 4, 9, 20, 36]
                out = decon_mask_edge_erosion(gpuArray(A), E);
                testCase.verifyEqual(gather(out), decon_mask_edge_erosion(A, E), sprintf('EdgeErosion %d', E));
            end
            % masks of one row or column (the previous GPU code
            % errored: gpuArray assignments take no repeated subscripts)
            for sz = {[1, 30, 20], [30, 1, 20], [30, 20, 1]}
                A = testCase.makeMask('random', sz{1});
                testCase.verifyEqual(gather(decon_mask_edge_erosion(gpuArray(A), 3)), decon_mask_edge_erosion_before(A, 3), ...
                    sprintf('size %s', mat2str(sz{1})));
            end
        end
    end

    methods (Static, Access = private)
        function [out] = outcome(f)
            % the output of f or the identifier of its error
            try
                out = f();
            catch ME
                out = ME.identifier;
            end
        end

        function [A] = makeMask(shape, sz)
            % 2D or 3D masks: random, skewed data's (zero triangles along
            % x-z), a thin plane, one true or one false voxel, blobs
            sz3 = [sz, ones(1, 3 - numel(sz))];
            switch shape
                case 'random'
                    A = rand(sz) > 0.3;
                case 'skewed'
                    [~, x, z] = ndgrid(1 : sz3(1), 1 : sz3(2), 1 : sz3(3));
                    A = reshape((x > 0.6 * z + sz3(2) * 0.05) & (x < 0.6 * z + sz3(2) * 0.5), sz);
                case 'thin'
                    A = false(sz);
                    A(:, max(1, round(end / 2)) + (0 : min(1, sz(2) - 1)), :) = true;
                case 'single'
                    A = false(sz);
                    A(ceil(numel(A) / 2)) = true;
                case 'singleFalse'
                    A = true(sz);
                    A(ceil(numel(A) / 2)) = false;
                case 'blobs'
                    A = imdilate(rand(sz) > 0.99, ones(3, 3));
            end
        end
    end
end
