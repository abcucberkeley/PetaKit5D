function [] = mex_exit_probe(mode, mexDir)
% diagnostic of MATLAB's exit on Windows after the PetaKit5D OpenMP mex files
% ran (CI job exit-probe; not a test). It calls each of them once with an
% input big enough to start its OpenMP threads, writes
% probe-results/<mode>.txt with a line per call (and whether the mex is
% locked), and returns, so that matlab -batch exits; the CI step's timeout
% shows a MATLAB that does not exit by itself.
%   mode 'committed': the repository's mex files
%   mode 'locked': the mex files in mexDir, built to lock themselves on
%     Windows too (shadowing the repository's)
%   mode 'unload': the repository's, then clear mex (which unloads them) with
%     a line before and after it
%   mode 'io': only the cpp-tiff and cpp-zarr mex files (writetiff, readtiff,
%     writezarr, readzarr)

if nargin < 2
    mexDir = '';
end
root = pwd;
outDir = fullfile(root, 'probe-results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
resultFn = fullfile(outDir, [mode, '.txt']);
note = @(varargin) writelines(string(sprintf(varargin{:})), resultFn, 'WriteMode', 'append');
note('%s start %s, MATLAB %s, %d threads', mode, datestr(now, 'HH:MM:SS'), version, maxNumCompThreads);

if strcmp(mode, 'io')
    addpath(root);
    evalc('setup(root, false, '''')');
    try
        im = uint16(randi(60000, 128, 128, 64));
        fn = fullfile(tempdir, 'mex_exit_probe');
        writetiff(im, [fn, '.tif']);
        a = readtiff([fn, '.tif']);
        writezarr(im, [fn, '.zarr']);
        b = readzarr([fn, '.zarr']);
        note('io: tiff %d (%s), zarr %d', isequal(a, im), which('parallelReadTiff'), isequal(b, im));
    catch ME
        note('io: error %s: %s', ME.identifier, ME.message);
    end
    note('done %s', datestr(now, 'HH:MM:SS'));
    return;
end

for d = {'crop', 'decon', 'deskew_rotate', 'flatfield_correction', 'stitch', 'tools/MIP', 'unmix'}
    addpath(fullfile(root, 'microscopeDataProcessing', d{1}, 'mex'));
end
addpath(fullfile(root, 'pointDetection', 'mex_functions'));
if ~isempty(mexDir)
    addpath(mexDir, '-begin');
end

rng(1);
nt = 4;
A16 = uint16(randi(60000, 128, 128, 64));
L = rand(96, 96, 48) > 0.3;
S = single(rand(128, 128, 64));
names = {'crop3d_mex', 'crop4d_mex', 'indexing3d_mex', 'indexing4d_mex', 'indexing4d_crop_mex', 'decon_lucy_accel_mex', ...
    'decon_mask_edge_erosion_mex', 'omw_backprojector_mex', 'skewed_space_interp_defined_stepsize_mex', ...
    'skewed_space_interp_volume_deskew_rotate_warp_mex', 'volume_deskew_rotate_warp_mex', 'flat_field_correction_3d_mex', ...
    'any_4th_dim_mex', 'feather_blending_3d_mex', 'feather_blending_3d_with_indexing_mex', 'feather_distance_map_resize_3d_mex', ...
    'integral_image_3d_mex', 'tile_blending_3d_mex', 'max_pooling_3d_mex', 'range_max_3d_mex', 'unmix_channels_sum_mex', ...
    'replace_nan_with_value'};
for i = 1 : numel(names)
    name = names{i};
    t0 = tic;
    try
        switch name
            case 'crop3d_mex'
                out = crop3d_mex(A16, [1, 1, 1, 100, 100, 50]);
            case 'crop4d_mex'
                out = crop4d_mex(repmat(A16, 1, 1, 1, 2), [1, 1, 1, 1, 100, 100, 50, 2]);
            case 'indexing3d_mex'
                indexing3d_mex(zeros(128, 128, 64, 'uint16'), [1, 1, 1, 100, 100, 50], A16(1 : 100, 1 : 100, 1 : 50));
            case 'indexing4d_mex'
                indexing4d_mex(zeros(128, 128, 64, 2, 'uint16'), [1, 1, 1, 1, 100, 100, 50, 1], A16(1 : 100, 1 : 100, 1 : 50));
            case 'indexing4d_crop_mex'
                indexing4d_crop_mex(zeros(128, 128, 64, 2, 'uint16'), [1, 1, 1, 1, 100, 100, 50, 1], A16, [1, 1, 1, 1, 100, 100, 50, 1]);
            case 'decon_lucy_accel_mex'
                out = decon_lucy_accel_mex('sum', S(:), nt);
            case 'decon_mask_edge_erosion_mex'
                out = decon_mask_edge_erosion_mex(L, {true}, [0, 0], [0, 0], 4, false, nt);
            case 'omw_backprojector_mex'
                out = omw_backprojector_mex('medfilt3', S(1 : 64, 1 : 64, 1 : 32), nt);
            case 'skewed_space_interp_defined_stepsize_mex'
                out = skewed_space_interp_defined_stepsize_mex(S, 2.5, 0.5, true, false);
            case 'skewed_space_interp_volume_deskew_rotate_warp_mex'
                out = skewed_space_interp_volume_deskew_rotate_warp_mex(S, 2.5, 0.5, true, eye(4), [1, 1, 1, 128, 128, 64], false);
            case 'volume_deskew_rotate_warp_mex'
                out = volume_deskew_rotate_warp_mex(S, eye(4), [1, 1, 1, 128, 128, 64], false);
            case 'flat_field_correction_3d_mex'
                out = flat_field_correction_3d_mex(A16, ones(128, 128, 'single'), zeros(128, 128, 'single'), true);
            case 'any_4th_dim_mex'
                [out, ~] = any_4th_dim_mex(repmat(A16, 1, 1, 1, 3), [true, false, true]);
            case 'feather_blending_3d_mex'
                out = feather_blending_3d_mex(repmat(A16, 1, 1, 1, 2), ones(128, 128, 64, 2, 'single'));
            case 'feather_blending_3d_with_indexing_mex'
                feather_blending_3d_with_indexing_mex(repmat(A16, 1, 1, 1, 2), ones(128, 128, 64, 2, 'single'), ...
                    zeros(128, 128, 64, 'uint16'), [1, 1, 1, 128, 128, 64]);
            case 'feather_distance_map_resize_3d_mex'
                out = feather_distance_map_resize_3d_mex(S(1 : 64, 1 : 64, 1 : 32), [128, 128, 64], 2);
            case 'integral_image_3d_mex'
                out = integral_image_3d_mex(S, [4, 4, 4]);
            case 'tile_blending_3d_mex'
                out = tile_blending_3d_mex(repmat(A16, 1, 1, 1, 2), [1, 1, 1, 1, 100, 100, 50; 2, 1, 1, 1, 100, 100, 50], 'mean');
            case 'max_pooling_3d_mex'
                out = max_pooling_3d_mex(A16, [2, 2, 2]);
            case 'range_max_3d_mex'
                out = range_max_3d_mex(A16, 3, 1);
            case 'unmix_channels_sum_mex'
                out = unmix_channels_sum_mex({A16, A16}, [1, -0.5], [100, 100], 0, 'uint16');
            case 'replace_nan_with_value'
                out = replace_nan_with_value(S, single(0));
        end
        status = 'ok';
    catch ME
        status = ['error ', ME.identifier, ': ', ME.message];
    end
    note('%-50s locked %d %5.2f s %s (%s)', name, mislocked(name), toc(t0), status, which(name));
end
clear out
if strcmp(mode, 'unload')
    note('before clear mex %s', datestr(now, 'HH:MM:SS'));
    clear mex
    note('after clear mex %s', datestr(now, 'HH:MM:SS'));
end
note('done %s', datestr(now, 'HH:MM:SS'));

end
