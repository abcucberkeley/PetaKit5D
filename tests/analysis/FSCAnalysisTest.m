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
    end
end
