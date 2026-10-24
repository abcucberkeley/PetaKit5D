classdef FSCAnalysisTest < pktest.TestCase
    % XR_FSC_analysis_wrapper ("FSC analysis" in demo_useful_tools.m), which
    % runs Fourier shell correlation on decon iterations (Iter_*.tif) to
    % pick the number of iterations. The single-image FSC correlates
    % interleaved sub-lattices of the image, so it measures the noise-limited
    % resolution: the same structure with more noise (as with too many decon
    % iterations) must get a worse (larger) resolution value.
    %
    % XR_FSC_resolution uses fit() from the Curve Fitting Toolbox, which is
    % not in the README's list of required toolboxes.

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
            switch resAxis
                case 'yz'
                    im = permute(im, [2, 1, 3]);
                case 'xy'
                    im = permute(im, [3, 2, 1]);
            end
            if dz ~= xy
                im = imresize3(im, round(size(im) .* [xy, xy, dz] / min(xy, dz)), 'linear');
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
                    'dz', dz, 'dr', 2, 'dtheta', pi / 6, 'halfSize', h);
            end
            fsc_mu = fscs{1};
            fsc_mu.fsc = (fscs{1}.fsc + fscs{2}.fsc + fscs{3}.fsc + fscs{4}.fsc) / 4;
            res_mu = XR_FSC_resolution(fsc_mu, 'xyPixelSize', xy, 'dz', dz, 'resThreshMethod', 'fixed', 'resThresh', 0.2);
        end
    end
end
