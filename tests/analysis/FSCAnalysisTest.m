classdef FSCAnalysisTest < pktest.TestCase
    % XR_FSC_analysis_wrapper ("FSC analysis" in demo_useful_tools.m), which
    % runs Fourier shell correlation on decon iterations (Iter_*.tif) to
    % pick the number of iterations. The single-image FSC correlates
    % interleaved sub-lattices of the image, so it measures the noise-limited
    % resolution: the same structure with more noise (as with too many decon
    % iterations) must get a worse (larger) resolution value.
    %
    % XR_FSC_resolution uses fit() from the Curve Fitting Toolbox (listed in
    % the README for the FSC analysis).

    properties (TestParameter)
        % {image size, class, xyPixelSize, dz, halfSize, bbox, resAxis,
        % clipPer}: resized up in z or in x and y, isotropic (no resize),
        % cropped in some dimensions and padded in others, padded in all,
        % an even half size, and the options applied before the resize
        fscCase = struct( ...
            'zUp', {{[60, 57, 23], 'uint16', 0.108, 0.3, [21, 21, 21], [], 'xz', []}}, ...
            'xyUp', {{[30, 28, 50], 'uint16', 0.108, 0.05, [21, 21, 21], [], 'xz', []}}, ...
            'padZ', {{[60, 57, 9], 'uint16', 0.108, 0.3, [21, 21, 21], [], 'xz', []}}, ...
            'padAll', {{[20, 22, 9], 'uint16', 0.108, 0.3, [21, 21, 21], [], 'xz', []}}, ...
            'isotropic', {{[50, 47, 45], 'single', 0.1, 0.1, [21, 21, 21], [], 'xz', []}}, ...
            'evenHalfSize', {{[60, 57, 23], 'uint16', 0.108, 0.3, [20, 20, 20], [], 'xz', []}}, ...
            'yzBbox', {{[60, 57, 23], 'uint16', 0.108, 0.2, [15, 15, 15], [5, 3, 2, 55, 50, 21], 'yz', []}}, ...
            'xyClip', {{[40, 37, 35], 'double', 0.108, 0.15, [17, 17, 17], [], 'xy', 99}})
    end

    methods (TestMethodSetup)
        function requireCurveFitting(testCase)
            testCase.assumeTrue(~isempty(ver('curvefit')) && license('test', 'Curve_Fitting_Toolbox'), ...
                'FSC resolution estimation needs the Curve Fitting Toolbox.');
            % XR_one_image_FSC_analysis caches the last image by file name
            clear XR_one_image_FSC_analysis
        end
    end

    methods (Test, TestTags = {'Integration'})
        function noisierIterationHasWorseResolution(testCase)
            dataPath = testCase.createTempFolder();
            b = imgaussfilt3(randn(64, 64, 64, 'single'), 2);
            signal = 1000 + 300 * b / std(b, 0, 'all');
            noise = randn(64, 64, 64, 'single');
            writetiff(signal + 20 * noise, [dataPath, 'Iter_0002.tif']);
            writetiff(signal + 600 * noise, [dataPath, 'Iter_0004.tif']);
            writetiff(1000 + 600 * noise, [dataPath, 'Iter_0006.tif']);

            XR_FSC_analysis_wrapper({dataPath}, 'resultDirName', 'FSCs', 'xyPixelSize', 0.1, 'dz', 0.1, ...
                'dr', 2, 'dtheta', pi / 6, 'resThreshMethod', 'fixed', 'resThresh', 0.2, ...
                'halfSize', [32, 32, 32], 'channelPatterns', {'Iter'}, 'channels', 488, ...
                'iterInterval', 2, 'parseCluster', false, 'masterCompute', true);

            for fsn = {'Iter_0002', 'Iter_0004', 'Iter_0006'}
                a = load([dataPath, 'FSCs/', fsn{1}, '.mat']);
                testCase.verifyTrue(all(isfield(a, {'fsc_mu', 'res_mu', 'fsc_cell', 'res_cell'})));
                testCase.verifyNumElements(a.fsc_cell, 4);
            end

            info = load([dataPath, 'FSCs/fsc_res_info_c1.mat']);
            ntheta = round(2 * pi / (pi / 6)) + 1;
            testCase.verifySize(info.res_mat, [3, ntheta]);
            res = info.res_mat;
            res(~isfinite(res)) = nan;
            meanRes = mean(res, 2, 'omitnan');
            testCase.verifyLessThan(meanRes(1), meanRes(2), 'More noise should give a worse resolution.');
            testCase.verifyLessThan(meanRes(2), meanRes(3), 'Pure noise should give the worst resolution.');
            % no better than the voxel size allows (Nyquist on the diagonal)
            testCase.verifyGreaterThanOrEqual(meanRes, 2 * 0.1 / sqrt(2) - 1e-6);

            figs = dir([dataPath, 'figures/*.png']);
            testCase.verifyGreaterThanOrEqual(numel(figs), 4);
        end

        function existingResultsAreNotRecomputed(testCase)
            dataPath = testCase.createTempFolder();
            writetiff(1000 + 300 * imgaussfilt3(randn(48, 48, 48, 'single'), 1), [dataPath, 'Iter_0002.tif']);
            args = {'xyPixelSize', 0.1, 'dz', 0.1, 'dr', 2, 'halfSize', [24, 24, 24], 'channelPatterns', {'Iter'}, ...
                'channels', 488, 'parseCluster', false, 'masterCompute', true};

            XR_FSC_analysis_wrapper({dataPath}, args{:});
            resultFn = [dataPath, 'FSCs/Iter_0002.mat'];
            first = dir(resultFn);
            XR_FSC_analysis_wrapper({dataPath}, args{:});
            second = dir(resultFn);
            testCase.verifyEqual(second.datenum, first.datenum);
        end

        function regionsAreCentredOnEachFile(testCase)
            % the region centres came from cs(1), cs(2) and cs(3) of the files'
            % centres (the first files' y centres when there are several
            % files), and with one region per image the wrapper failed when
            % it read the results (MATLAB:cellRefFromNonCell): each file's
            % regions are around its own centre, here one region each
            dataPath = testCase.createTempFolder();
            sizes = [64, 72, 40; 80, 56, 44];
            for i = 1 : 2
                writetiff(single(1000 + 300 * pktest.texture(sizes(i, :), 1.5, [0, 1]) + 30 * randn(sizes(i, :))), ...
                    sprintf('%sIter_%04d.tif', dataPath, 2 * i));
            end
            h = [16, 16, 16];
            XR_FSC_analysis_wrapper({dataPath}, 'xyPixelSize', 0.1, 'dz', 0.1, 'dr', 2, 'dtheta', pi / 6, ...
                'halfSize', h, 'channelPatterns', {'Iter'}, 'channels', 488, 'iterInterval', 2, ...
                'multiRegions', true, 'regionInterval', [-1, -1, -1], 'parseCluster', false, 'masterCompute', true);
            for i = 1 : 2
                a = load(sprintf('%sFSCs/Iter_%04d.mat', dataPath, 2 * i), 'bbox');
                c = round((sizes(i, :) + 1) / 2);
                testCase.verifyEqual(a.bbox, [c - h, c + h - 1], sprintf('file %d', i));
            end
            info = load([dataPath, 'FSCs/fsc_res_info_c1.mat']);
            testCase.verifySize(info.res_mat, [2, 13]);

            % a grid of regions: those of the second file are around its centre
            rmdir([dataPath, 'FSCs'], 's');
            XR_FSC_analysis_wrapper({dataPath}, 'xyPixelSize', 0.1, 'dz', 0.1, 'dr', 2, 'dtheta', pi / 6, ...
                'halfSize', h, 'channelPatterns', {'Iter'}, 'channels', 488, 'iterInterval', 2, ...
                'multiRegions', true, 'regionInterval', [10, 10, -1], 'parseCluster', false, 'masterCompute', true);
            a = load([dataPath, 'FSCs/Iter_0004.mat'], 'bbox');
            c = round((sizes(2, :) + 1) / 2);
            centres = round((a.bbox(:, 1 : 3) + a.bbox(:, 4 : 6)) / 2);
            % (steps of 10 from its centre, at least the half size from its edges)
            ys = [fliplr(c(1) - 10 : -10 : h(1) + 1), c(1) : 10 : sizes(2, 1) - h(1) + 1];
            xs = [fliplr(c(2) - 10 : -10 : h(2) + 1), c(2) : 10 : sizes(2, 2) - h(2) + 1];
            testCase.verifyEqual(unique(centres(:, 1))', ys);
            testCase.verifyEqual(unique(centres(:, 2))', xs);
            testCase.verifyEqual(unique(centres(:, 3))', c(3));
        end

        function resultIsMovedIntoPlace(testCase)
            % the job framework counts a task as done once its result
            % exists, and the wrapper loads each result right after it
            % returns; a result saved in place could be read half written
            dataPath = testCase.createTempFolder();
            writetiff(1000 + 300 * imgaussfilt3(randn(48, 48, 48, 'single'), 1), [dataPath, 'Iter_0002.tif']);
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            XR_FSC_analysis_wrapper({dataPath}, 'xyPixelSize', 0.1, 'dz', 0.1, 'dr', 2, 'halfSize', [24, 24, 24], ...
                'channelPatterns', {'Iter'}, 'channels', 488, 'parseCluster', false, 'masterCompute', true);
            testCase.verifyTrue(moves.wasMoved([dataPath, 'FSCs/Iter_0002.mat']), 'The result was saved in place.');
            testCase.verifyEmpty(dir([dataPath, 'FSCs/Iter_0002_*.mat']), 'A temporary file was left behind.');
        end
    end

    methods (Test, TestTags = {'Unit'})
        function oneImageMatchesResizingTheWholeImage(testCase, fscCase)
            % only the region kept after the isotropic resize is resized
            % (from the input region it depends on); the FSC must equal that
            % of resizing the whole image and cropping it
            [sz, cls, xy, dz, h, bbox, resAxis, clipPer] = fscCase{:};
            im = cast(1000 + 300 * pktest.texture(sz, 1.5, [0, 1]) + 30 * randn(sz), cls);
            [fsc_ref, res_ref] = FSCAnalysisTest.referenceFSC(im, xy, dz, h, bbox, resAxis, clipPer);

            % the resize must not fall back to imresize3 of the whole image
            testCase.applyFixture(pktest.FunctionSpyFixture('imresize3'));
            [fsc_mu, res_mu] = XR_one_image_FSC_analysis(im, 'xyPixelSize', xy, 'dz', dz, 'dr', 2, 'dtheta', pi / 6, ...
                'halfSize', h, 'bbox', bbox, 'resAxis', resAxis, 'clipPer', clipPer);
            testCase.verifyEqual(fsc_mu, fsc_ref);
            testCase.verifyEqual(res_mu, res_ref);
        end

        function twoImageAnalysisRuns(testCase)
            % XR_two_image_FSC_analysis passed 'N' to XR_FSC, which has no such
            % parameter (its halfSize), so it always failed; its crop kept the
            % second half of an image larger than 2 * N. Two copies of one
            % image correlate fully (the FSC never falls below a fixed 0.2: the
            % resolution is that of the largest frequency, sqrt(2) of the
            % Nyquist frequency), two images of independent noise much less
            im = 1000 + 300 * pktest.texture([60, 56, 52], 1.5, [0, 1]) + 30 * randn(60, 56, 52);
            args = {'xyPixelSize', 0.1, 'dz', 0.1, 'dr', 2, 'N', [24, 24, 24], 'resThreshMethod', 'fixed'};
            [~, res] = XR_two_image_FSC_analysis(im, im, args{:});
            testCase.verifySize(res.res, [25, 1]);
            testCase.verifyEqual(res.res, repmat(2 * 0.1 / sqrt(2), 25, 1), 'AbsTol', 1e-12);
            [~, resNoise] = XR_two_image_FSC_analysis(randn(60, 56, 52), randn(60, 56, 52), args{:});
            testCase.verifyGreaterThan(min(resNoise.res), 1.5 * res.res(1));

            % the crop of a larger image is its centre (as a centred crop of
            % 2 * N given to it)
            big = 1000 + 300 * pktest.texture([60, 120, 52], 1.5, [0, 1]);
            other = big + 100 * randn(size(big));
            fscC = XR_two_image_FSC_analysis(big, other, args{:});
            s = max(1, floor((size(big) - 48) / 2));
            fscR = XR_two_image_FSC_analysis(big(:, s(2) : s(2) + 47, :), other(:, s(2) : s(2) + 47, :), args{:});
            testCase.verifyEqual(fscC.fsc, fscR.fsc);
        end

        function sectorsInsideTheConesAreNaN(testCase)
            % with dtheta = pi/24 the sectors at 0 and pi/2 lie inside the cones
            % that are skipped: XR_FSC gives the 22 others, which took their
            % places in the full circle (45 values for 49 angles). Each angle
            % phi takes the sector at phi + pi/2 (0 is x, pi/2 z), and the
            % missing sectors are NaN
            a = randn(41, 41, 41);
            [fsc, res] = XR_FSC(a, a + 0.5 * randn(41, 41, 41), 'halfSize', [41, 41, 41], 'dr', 2, 'dtheta', pi / 24);
            sectors = round(fsc.thetas / (pi / 24)) + 1;
            testCase.assertEqual(sectors', [2 : 12, 14 : 24]);
            testCase.verifySize(res.res, [49, 1]);
            testCase.verifySize(res.thetas, [49, 1]);
            values = 2 * 0.108 ./ res.points;
            for j = 1 : 49
                k = find(sectors == mod(j - 1 + 12, 24) + 1);
                if isempty(k)
                    testCase.verifyTrue(isnan(res.res(j)), sprintf('angle %d', j));
                else
                    testCase.verifyEqual(res.res(j), values(k), sprintf('angle %d', j));
                end
            end
            testCase.verifyEqual(find(isnan(res.res))', [1, 13, 25, 37, 49]);
        end

        function oddNumberOfSectorsIsRefused(testCase)
            % with an odd number of sectors (pi / dtheta = 5) no sector is at
            % pi/2 (x): the full circle had 2n values for 2n + 1 angles, with
            % two MATLAB:colon:nonIntegerIndex warnings
            a = randn(41, 41, 41);
            testCase.verifyError(@() XR_FSC(a, a + randn(41, 41, 41), 'halfSize', [41, 41, 41], 'dr', 2, ...
                'dtheta', pi / 5), 'XR_FSC_resolution:oddSectorNumber');
        end

        function imagesWithoutSignalHaveNoFSC(testCase)
            % all-zero images failed with an index error (xi(end) = [] of an
            % empty xi)
            testCase.verifyError(@() XR_FSC(zeros(41, 41, 41), zeros(41, 41, 41), 'halfSize', [41, 41, 41], 'dr', 2), ...
                'XR_FSC:noSignal');
        end

        function shellCacheFollowsParameters(testCase)
            % XR_FSC keeps the window and the frequency shells between calls:
            % the same call twice and calls with another dr, dtheta,
            % skipConeRegion or halfSize must give the results of fresh calls
            a = randn(21, 21, 21);
            b = a + randn(21, 21, 21);
            c = randn(25, 25, 25);
            d = c + randn(25, 25, 25);
            h = {'halfSize', [21, 21, 21]};
            calls = {[{a, b, 'dr', 2}, h], [{a, b, 'dr', 2}, h], [{a, b, 'dr', 3}, h], ...
                [{a, b, 'dr', 2, 'dtheta', pi / 4}, h], [{a, b, 'dr', 2, 'skipConeRegion', false}, h], ...
                {c, d, 'dr', 2, 'halfSize', [25, 25, 25]}, [{a, b, 'dr', 2}, h]};
            fscs = cell(size(calls));
            ress = cell(size(calls));
            clear XR_FSC
            for i = 1 : numel(calls)
                [fscs{i}, ress{i}] = XR_FSC(calls{i}{:});
            end
            for i = 1 : numel(calls)
                clear XR_FSC
                [fsc, res] = XR_FSC(calls{i}{:});
                testCase.verifyEqual(fscs{i}, fsc, sprintf('Call %d', i));
                testCase.verifyEqual(ress{i}, res, sprintf('Call %d', i));
            end
        end
    end

    methods (Static, Access = private)
        function [fsc_mu, res_mu] = referenceFSC(im, xy, dz, h, bbox, resAxis, clipPer)
            % XR_one_image_FSC_analysis before only the kept region of the
            % resize was computed: resize the whole image, crop, pad
            im = double(im);
            if ~isempty(clipPer)
                clip_val = prctile(im(im > 0), clipPer);
                im = im .* (im <= clip_val) + clip_val .* (im > clip_val);
            end
            if ~isempty(bbox)
                im = im(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
            end
            % ('xy': the pixel sizes of z, x and y, and its plane, x and y, has
            % the lateral pixel size; they were [xy xy dz] and dz)
            pixelSizes = [xy, xy, dz];
            fscDz = dz;
            switch resAxis
                case 'yz'
                    im = permute(im, [2, 1, 3]);
                case 'xy'
                    im = permute(im, [3, 2, 1]);
                    pixelSizes = [dz, xy, xy];
                    fscDz = xy;
            end
            if dz ~= xy
                im = imresize3(im, round(size(im) .* pixelSizes / min(xy, dz)), 'linear');
            end
            sz = size(im);
            if any(sz > h * 2)
                s = max(1, floor((sz - h * 2) / 2));
                t = min(sz, s + h * 2 - 1);
                im = im(s(1) : t(1), s(2) : t(2), s(3) : t(3));
                sz = size(im);
            end
            if any(sz < h * 2)
                lsz = floor((h * 2 - sz) / 2);
                im = padarray(padarray(im, lsz, 0, 'pre'), h * 2 - sz - lsz, 0, 'post');
            end
            im = permute(reshape(im, 2, h(1), 2, h(2), 2, h(3)), [2, 4, 6, 1, 3, 5]);
            pairs = {[1, 1, 1, 2, 2, 2], [2, 1, 1, 1, 2, 2], [1, 2, 1, 2, 1, 2], [1, 1, 2, 2, 2, 1]};
            fscs = cell(1, 4);
            for i = 1 : 4
                p = pairs{i};
                fscs{i} = XR_FSC(im(:, :, :, p(1), p(2), p(3)), im(:, :, :, p(4), p(5), p(6)), 'xyPixelSize', xy, ...
                    'dz', fscDz, 'dr', 2, 'dtheta', pi / 6, 'halfSize', h);
            end
            fsc_mu = fscs{1};
            fsc_mu.fsc = (fscs{1}.fsc + fscs{2}.fsc + fscs{3}.fsc + fscs{4}.fsc) / 4;
            res_mu = XR_FSC_resolution(fsc_mu, 'xyPixelSize', xy, 'dz', fscDz, 'resThreshMethod', 'fixed', 'resThresh', 0.2);
        end
    end
end
