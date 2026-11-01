function [gaps, tracked] = mexPlatformGaps(ext)
% the mex functions built for Linux but not for the platform with mex
% extension ext (default: this one), as known gaps: MexInventoryTest allows
% them, and a test that runs the MATLAB code instead of one of them doesn't
% fail (pktest.TestCase). tracked is false for platforms other than Linux,
% Windows and ARM macOS, whose builds aren't tracked.
%
% When you build a missing binary, remove it from the list.

if nargin < 1
    ext = mexext;
end

% not called by any current code path (only in comments, *_bak.m files, or
% not at all), and the decon mex files, which are built on Linux only (the
% CPU one equals MATLAB's FFTs there; mexcuda needs MSVC on Windows, and
% macOS has no CUDA), so decon runs the MATLAB code there (the CPU erosion
% mex is built on all platforms)
windowsGaps = {'min_bbox_3d_mex', 'h_mat_mex', 'h_mat_fast_exp_mex', ...
    'h_mat_mthreads_mex', 'replace_nan_inf_with_value', 'decon_lucy_rfft_mex', 'decon_lucy_rfft_gpu_mex', ...
    'decon_mask_edge_erosion_gpu_mex'};

% unused ones as above, plus the point detection and imageProcessing mex
% functions, which have never been built for ARM macOS. The .m files with
% the same names are help-text stubs, so point detection does not work on
% ARM macOS.
macArmGaps = [windowsGaps, {'feather_distance_map_resize_3d_w_space_mex', ...
    'conv3fast', 'binterp', 'skeleton3D', 'steerableDetector', ...
    'steerableDetector3D', 'createDistanceMatrix', 'createDiffMatrix', 'fitGaussian2D', ...
    'fitGaussianMixture2D', 'fitGaussian3D', 'fitGaussianMixture3D', 'KDTreeBallQuery', 'mexLap', ...
    'vectorialPSF', 'adtest_mex', 'compute_cross_hessian_component_mex', ...
    'compute_gradient_hessian_mex', 'replace_nan_with_value'}];

tracked = true;
switch ext
    case 'mexa64'
        gaps = {};
    case 'mexw64'
        gaps = windowsGaps;
    case 'mexmaca64'
        gaps = macArmGaps;
    otherwise
        gaps = {};
        tracked = false;
end

end
