classdef FlatfieldBackgroundEstimationTest < pktest.TestCase
    % XR_flatfield_background_wrapper: estimates the camera background from
    % dark (sample free) stacks. Slices noisier than the median slice and
    % slices with bright structures are excluded, so the estimate is the
    % dark level even with some sample content in the stack.

    properties (Constant)
        Fsn = 'Scan_Iter_0000_0000_%s_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t.tif'
        DarkLevel = 100
    end

    methods (Test, TestTags = {'Integration'})
        function estimatesDarkLevelPerCamera(testCase)
            dataPath = testCase.createTempFolder();
            pattern = 100 + 10 * sin((1 : 48)' / 6) * cos((1 : 40) / 5);
            for cam = {'CamA', 'CamB'}
                writetiff(testCase.darkStack(pattern), [dataPath, sprintf(testCase.Fsn, cam{1})]);
            end

            XR_flatfield_background_wrapper({dataPath}, 'channelPatterns', {'CamA_ch0', 'CamB_ch0'}, ...
                'parseCluster', false, 'masterCompute', true);

            resultPath = [dataPath, 'matlab_flat_field_estimation/'];
            for cam = {'A', 'B'}
                bgFn = sprintf('%sScan_Cam%s_ch0_CAM1_background.tif', resultPath, cam{1});
                testCase.assertTrue(exist(bgFn, 'file') == 2, sprintf('No background image for Cam%s.', cam{1}));
                bg = readtiff(bgFn);
                testCase.verifySize(bg, [48, 40]);
                % the per-pixel dark pattern is recovered within the noise of
                % a mean over ~20 slices with sigma 2-4
                err = double(bg(:)) - pattern(:);
                testCase.verifyLessThan(abs(mean(err)), 0.3);
                testCase.verifyLessThan(max(abs(err)), 5);
            end
            a = load(sprintf('%s%s.mat', resultPath, sprintf(testCase.Fsn(1 : end - 4), 'CamA')));
            testCase.verifyGreaterThan(a.numSlice, 0);
            testCase.verifyLessThan(a.numSlice, 40, 'Noisy or bright slices should be excluded.');
        end

        function brightSlicesAreExcluded(testCase)
            % a bright object in some slices must not raise the estimate
            dataPath = testCase.createTempFolder();
            pattern = testCase.DarkLevel * ones(48, 40);
            stack = single(testCase.darkStack(pattern));
            stack(10 : 30, 10 : 30, 1 : 6) = stack(10 : 30, 10 : 30, 1 : 6) + 2000;
            writetiff(uint16(stack), [dataPath, sprintf(testCase.Fsn, 'CamA')]);

            XR_flatfield_background_wrapper({dataPath}, 'channelPatterns', {'CamA_ch0'}, ...
                'parseCluster', false, 'masterCompute', true);

            bg = readtiff([dataPath, 'matlab_flat_field_estimation/Scan_CamA_ch0_CAM1_background.tif']);
            err = double(bg(:)) - testCase.DarkLevel;
            testCase.verifyLessThan(abs(mean(err)), 0.3);
            testCase.verifyLessThan(max(abs(err)), 5);
        end

        function acceptsMccModeFalse(testCase)
            % mccMode used to be validated with @ischar
            dataPath = testCase.createTempFolder();
            writetiff(testCase.darkStack(testCase.DarkLevel * ones(48, 40)), [dataPath, sprintf(testCase.Fsn, 'CamA')]);

            XR_flatfield_background_wrapper({dataPath}, 'channelPatterns', {'CamA_ch0'}, ...
                'parseCluster', false, 'masterCompute', true, 'mccMode', false);

            testCase.verifyTrue(exist([dataPath, 'matlab_flat_field_estimation/Scan_CamA_ch0_CAM1_background.tif'], 'file') == 2);
        end
    end

    methods (Access = private)
        function [stack] = darkStack(~, pattern)
            % 40 dark slices whose noise level increases with z (the
            % estimator keeps slices quieter than the median slice)
            sigmas = permute(linspace(2, 6, 40), [1, 3, 2]);
            stack = uint16(pattern + sigmas .* randn([size(pattern), 40]));
        end
    end
end
