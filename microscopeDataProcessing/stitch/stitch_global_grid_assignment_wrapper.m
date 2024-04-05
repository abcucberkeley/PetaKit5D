function [] = stitch_global_grid_assignment_wrapper(inputFullpath, outputFullpath, xcorr_thresh, axisWeight, uuid)
% wrapper for grid assigment optimization used for parallel optimization
% for multiple groups


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inputFullpath', @(x) ischar(x));
ip.addRequired('outputFullpath', @(x) ischar(x));
ip.addRequired('xcorr_thresh', @(x) isscalar(x));
ip.addRequired('axisWeight', @(x) isvector(x));
ip.addRequired('uuid', @(x) ischar(x));

ip.parse(inputFullpath, outputFullpath, xcorr_thresh, axisWeight, uuid);

if exist(outputFullpath, 'file')
    fprintf('output %s already exists, skip it!\n', outputFullpath);
    return;
end

a = load(inputFullpath, 'nF_g', 'max_xcorr_mat_g', 'absolute_shift_mat_g', ...
    'overlap_matrix_g', 'max_allow_shift_g', 'tileIdx_g');
nF_g = a.nF_g;
max_xcorr_mat_g = a.max_xcorr_mat_g;
absolute_shift_mat_g = a.absolute_shift_mat_g;
overlap_matrix_g = a.overlap_matrix_g;
max_allow_shift_g = a.max_allow_shift_g;
tileIdx_g = a.tileIdx_g;

[d_shift_g] = stitch_global_grid_assignment(nF_g, max_xcorr_mat_g, absolute_shift_mat_g, ...
    overlap_matrix_g, max_allow_shift_g, xcorr_thresh, tileIdx_g, axisWeight);

outputTmppath = sprintf('%s_%s.mat', outputFullpath(1 : end - 4), uuid);
save('-v7.3', outputTmppath, 'd_shift_g');
movefile(outputTmppath, outputFullpath);

end
