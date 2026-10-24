classdef MexInventoryTest < pktest.TestCase
    % every mex function built for Linux must also be built for Windows and
    % ARM macOS, except the known gaps below. This catches a mex function
    % that is added or rebuilt for one platform but not committed for the
    % others. Linux binaries are the reference set.
    %
    % When you build a missing binary, remove it from the gap list.

    properties (Constant)
        % not called by any current code path (only in comments, *_bak.m
        % files, or not at all), and the decon mex files, which are built on
        % Linux only (the CPU one equals MATLAB's FFTs there; mexcuda needs
        % MSVC on Windows, and macOS has no CUDA), so decon runs the MATLAB
        % code there (the CPU erosion mex is built on all platforms)
        WindowsGaps = {'min_bbox_3d_mex', 'h_mat_mex', 'h_mat_fast_exp_mex', ...
            'h_mat_mthreads_mex', 'replace_nan_inf_with_value', 'decon_lucy_rfft_mex', 'decon_lucy_rfft_gpu_mex', ...
            'decon_mask_edge_erosion_gpu_mex'}

        % unused ones as above, plus the point detection and imageProcessing
        % mex functions, which have never been built for ARM macOS. The .m
        % files with the same names are help-text stubs, so point detection
        % does not work on ARM macOS.
        MacArmGaps = [MexInventoryTest.WindowsGaps, {'feather_distance_map_resize_3d_w_space_mex', ...
            'conv3fast', 'binterp', 'skeleton3D', 'steerableDetector', ...
            'steerableDetector3D', 'createDistanceMatrix', 'createDiffMatrix', 'fitGaussian2D', ...
            'fitGaussianMixture2D', 'fitGaussian3D', 'fitGaussianMixture3D', 'KDTreeBallQuery', 'mexLap', ...
            'vectorialPSF', 'adtest_mex', 'compute_cross_hessian_component_mex', ...
            'compute_gradient_hessian_mex', 'replace_nan_with_value'}]
    end

    methods (Test, TestTags = {'Unit'})
        function everyLinuxMexIsBuiltForThisPlatform(testCase)
            ext = mexext;
            switch ext
                case 'mexw64'
                    knownGaps = testCase.WindowsGaps;
                case 'mexmaca64'
                    knownGaps = testCase.MacArmGaps;
                otherwise
                    testCase.assumeFail(sprintf('Only checked on Windows and ARM macOS (this is %s; Linux is the reference).', ext));
            end

            linuxNames = testCase.mexNames('mexa64');
            platformNames = testCase.mexNames(ext);

            missing = setdiff(setdiff(linuxNames, platformNames), knownGaps);
            testCase.verifyEmpty(missing, sprintf(['These mex functions have a Linux build but no .%s build: %s. ', ...
                'Build and commit them, or add them to the known gaps with a reason.'], ext, strjoin(missing, ', ')));

            nowBuilt = intersect(knownGaps, platformNames);
            if ~isempty(nowBuilt)
                testCase.log(sprintf('Listed as known gaps but now built for %s (remove from the gap list): %s', ...
                    ext, strjoin(nowBuilt, ', ')));
            end
        end
    end

    methods (Static, Access = private)
        function [names] = mexNames(ext)
            dir_info = dir(fullfile(pktest.repoRoot(), '**', ['*.', ext]));
            [~, names] = cellfun(@fileparts, {dir_info.name}, 'unif', 0);
            names = unique(names);
        end
    end
end
