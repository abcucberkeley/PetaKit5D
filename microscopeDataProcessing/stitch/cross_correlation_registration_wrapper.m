function [] = cross_correlation_registration_wrapper(imgFullpath_1, imgFullpath_2, xcorrFullpath, pair_indices, cuboid_1, cuboid_2, cuboid_overlap_12, px, xyz_factors, varargin)
% wrapper for 2d/3d cross correlation computing, for giving tile indices,
% boundled by leading tile indices. 

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imgFullpath_1', @ischar);
ip.addRequired('imgFullpath_2', @iscell);
ip.addRequired('xcorrFullpath', @ischar);
ip.addRequired('pair_indices', @isnumeric);
ip.addRequired('cuboid_1', @isnumeric);
ip.addRequired('cuboid_2', @isnumeric);
ip.addRequired('cuboid_overlap_12', @isnumeric);
ip.addRequired('px', @isnumeric);
ip.addRequired('xyz_factors', @isnumeric);
ip.addParameter('Stitch2D', false, @islogical);
ip.addParameter('downSample', [1, 1, 1], @isnumeric);
ip.addParameter('MaxOffset', [300, 300, 50], @isnumeric);
ip.addParameter('dimNumThrsh', 10000, @isnumeric);
% ip.addParameter('xyz_factor', 0.5, @isnumeric);

ip.parse(imgFullpath_1, imgFullpath_2, xcorrFullpath, pair_indices, cuboid_1, cuboid_2, cuboid_overlap_12, px, xyz_factors, varargin{:});

pr = ip.Results;
Stitch2D = pr.Stitch2D;
downSample = pr.downSample;
MaxOffset = pr.MaxOffset;
dimNumThrsh = pr.dimNumThrsh;

if exist(xcorrFullpath, 'file')
    fprintf('The xcorr result %s already exists!\n', xcorrFullpath);
    return;
end

imgFullpath_2_str = sprintf('{''%s''}', strjoin(imgFullpath_2, ''',\n'''));
fprintf('Compute xcorr shifts for tile:\n    %s\n  and tiles:\n   %s\n', imgFullpath_1, imgFullpath_2_str);

nF = size(pair_indices, 1);

relative_shift_mat = zeros(nF, 3);
max_xcorr_mat = zeros(nF, 1);
for i = 1 : nF
    tic
    imgFullpath_2i = imgFullpath_2{i};
    cuboid_1i = cuboid_1;
    cuboid_2i = cuboid_2(i, :);
    cuboid_overlap_12i = cuboid_overlap_12(i, :);

    if Stitch2D
        [relative_shift, max_xcorr] = cross_correlation_registration_2d(imgFullpath_1, ...
            imgFullpath_2i, '', cuboid_1i, cuboid_2i, cuboid_overlap_12i, px, xyz_factors, ...
            downSample=downSample, MaxOffset=MaxOffset, dimNumThrsh=dimNumThrsh);
    else
        [relative_shift, max_xcorr] = cross_correlation_registration_3d(imgFullpath_1, ...
            imgFullpath_2i, '', cuboid_1i, cuboid_2i, cuboid_overlap_12i, px, xyz_factors, ...
            downSample=downSample, MaxOffset=MaxOffset, dimNumThrsh=dimNumThrsh);
    end
    relative_shift_mat(i, :) = relative_shift;
    max_xcorr_mat(i) = max_xcorr;
    toc
end

fprintf('Save results ...\n');

uuid = get_uuid();
xcorrTmppath = sprintf('%s_%s.mat', xcorrFullpath(1 : end - 4), uuid);
save('-v7.3', xcorrTmppath, 'relative_shift_mat', 'max_xcorr_mat', 'pair_indices');
fileattrib(xcorrTmppath, '+w', 'g');
movefile(xcorrTmppath, xcorrFullpath);

fprintf('Done!\n');

end
