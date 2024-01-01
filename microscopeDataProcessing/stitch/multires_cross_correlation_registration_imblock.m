function [relative_shift, max_xcorr] = multires_cross_correlation_registration_imblock(imgFullpath_1, imgFullpath_2, xcorrFullpath, cuboid_1, cuboid_2, cuboid_overlap_12, px, xyz_factors, varargin)
% compute shift between a pair of tiles in stitching based on
% cross-correlation using image block framework in a multi-resolution way.
% 
% relative_shift: negative for move toward; postive for move away.
% 
% Author: Xiongtao Ruan (10/19/2020)
% 
% xruan (12/07/2020): update constraints for the maximum shift 
% xruan (02/24/2021): add support for user defined xy, z max offsets for xcorr registration
% xruan (07/13/2021): fix issue if none part of region 2 is selected. 
% xruan (08/20/2021): make overlap region smaller when calculating xcorr
% xruan (10/08/2021): fix shape issue for cuboid coordinates after change to the generic framework


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imgFullpath_1', @ischar);
ip.addRequired('imgFullpath_2', @ischar);
ip.addRequired('xcorrFullpath', @ischar);
ip.addRequired('cuboid_1', @isnumeric);
ip.addRequired('cuboid_2', @isnumeric);
ip.addRequired('cuboid_overlap_12', @isnumeric);
ip.addRequired('px', @isnumeric);
ip.addRequired('xyz_factors', @isnumeric);
ip.addParameter('downSample', [1, 1, 1], @isnumeric);
ip.addParameter('xyMaxOffset', 300, @isnumeric);
ip.addParameter('zMaxOffset', 50, @isnumeric);
% ip.addParameter('xyz_factor', 0.5, @isnumeric);

ip.parse(imgFullpath_1, imgFullpath_2, xcorrFullpath, cuboid_1, cuboid_2, cuboid_overlap_12, px, xyz_factors, varargin{:});

pr = ip.Results;
downSample = pr.downSample;
xy_offset_mat = pr.xyMaxOffset;
z_offset_mat = pr.zMaxOffset;

% xy_offset_mat = -10 : 2 : 10;
% z_offset_mat = -3 : 3;
% xy_offset_mat = -5 : 5; 
% z_offset_mat = -3 : 3;

if exist(xcorrFullpath, 'file')
    fprintf('The xcorr result %s already exists!\n', xcorrFullpath);
    return;
end

fprintf('Compute xcorr shifts for tiles:\n    %s\n    %s\n', imgFullpath_1, imgFullpath_2);
fprintf('Load overlap regions ...\n');

% img_1(img_1 == 0) = nan;
% img_2(img_2 == 0) = nan;
[~, ~, ext] = fileparts(imgFullpath_1);
switch ext
    case '.tif'
        bim_1 = blockedImage(imgFullpath_1, "Adapter", MPageTiffAdapter);
        bim_2 = blockedImage(imgFullpath_2, "Adapter", MPageTiffAdapter);
    case '.zarr'
        bim_1 = blockedImage(imgFullpath_1, "Adapter", ZarrAdapter);
        bim_2 = blockedImage(imgFullpath_2, "Adapter", ZarrAdapter);
end

xy_num_off = numel(xy_offset_mat);
z_num_off = numel(z_offset_mat);
% xcorr_mat = zeros(xy_num_off, xy_num_off, z_num_off);

% [ny_1, nx_1, nz_1] = bim_1.Size;
% [ny_2, nx_2, nz_2] = bim_2.Size;
sz_1 = bim_1.Size;
% sz_2 = bim_2.Size;

% ny_1 = sz_1(1); nx_1 = sz_1(2); nz_1 = sz_1(3);
% ny_2 = sz_2(1); nx_2 = sz_2(2); nz_2 = sz_2(3);

% note: cuboid_overlap_12, cuboid_1 and xyz_factors are in xyz order.
% s1 = round((cuboid_overlap_12(:, 1) - cuboid_1(:, 1)) ./ (px * xyz_factors)) + 1;
s1 = round((cuboid_overlap_12(1 : 3) - cuboid_1(1 : 3))' ./ (px * xyz_factors)) + 1;
s2 = round((cuboid_overlap_12(1 : 3) - cuboid_2(1 : 3))' ./ (px * xyz_factors)) + 1;

t1 = round((cuboid_overlap_12(4 : 6) - cuboid_1(1 : 3))' ./ (px * xyz_factors)) + 1;
t2 = round((cuboid_overlap_12(4 : 6) - cuboid_2(1 : 3))' ./ (px * xyz_factors)) + 1;

% calculate bbox for overlap regions
maxoff_xy = max(abs(xy_offset_mat));
maxoff_z = max(abs(z_offset_mat));
boarder = [maxoff_xy; maxoff_xy; maxoff_z];
sr1 = max(s1 - boarder, 1);
tr1 = min(t1 + boarder, sz_1([2, 1, 3])');

% region_1 = bim_1.getRegion(sr1([2, 1, 3]), tr1([2, 1, 3]));
region_1 = bim_1.Adapter.getIORegion(sr1([2, 1, 3]), tr1([2, 1, 3]));
region_1 = single(region_1);

% region_2 = bim_2.getRegion(s2([2, 1, 3]), t2([2, 1, 3]));
region_2 = bim_2.Adapter.getIORegion(s2([2, 1, 3]), t2([2, 1, 3]));
region_2 = single(region_2);

% first resize to downsampled ones
if any(downSample ~= 1)
%     region_1 = double(imresize3(region_1, 1 ./ downSample .* size(region_1), 'nearest'));
%     region_2 = double(imresize3(region_2, 1 ./ downSample .* size(region_2), 'nearest'));
    region_1 = double(imresize3(region_1, 1 ./ downSample .* size(region_1), 'nearest'));
    region_2 = double(imresize3(region_2, 1 ./ downSample .* size(region_2), 'nearest'));
end

% region_1 = (region_1 / 100) .^ 5;
% region_2 = (region_2 / 100) .^ 5;

% crop region 2
sz_2 = size(region_2, [1, 2, 3]);
% region_2 = region_2(s2(2) : t2(2), s2(1) : t2(1), s2(3) : t2(3));
region_2_nozeros = region_2 ~= 0;
% select the region with at least 1/4 areas in the selected dimensions.
z_inds = sum(region_2_nozeros, [1, 2]) > prod(sz_2([1, 2])) / 4;
x_inds = sum(region_2_nozeros(:, :, z_inds), [1, 3]) > sz_2(1) * sum(z_inds) / 4;
y_inds = sum(region_2_nozeros(:, x_inds, z_inds), [2, 3]) > sum(x_inds) * sum(z_inds) / 4;

% reduce area threshold if the inds are empty
if ~any(z_inds) || ~any(x_inds) || ~any(y_inds)  
    z_inds = sum(region_2_nozeros, [1, 2]) > prod(sz_2([1, 2])) / 8;
    x_inds = sum(region_2_nozeros(:, :, z_inds), [1, 3]) > sz_2(1) * sum(z_inds) / 8;
    y_inds = sum(region_2_nozeros(:, x_inds, z_inds), [2, 3]) > sum(x_inds) * sum(z_inds) / 8;
end

% (08/20/2021): check if the indices are full range, if so, exclude the
% first and last several slices/planes

    function [inds] = cap_region_2_inds(inds, sz, max_off, bounds)
        inds = find(inds);
        inds = inds(inds > max(sz * bounds(1), min(sz * bounds(2), max_off / 2)) & inds < min(sz * (1 - bounds(1)), max(sz * (1 - bounds(2)), sz - max_off / 2)));
    end

bounds = [0.1, 0.25];
maxoff = [maxoff_xy, maxoff_xy, maxoff_z] ./ downSample;
y_inds = cap_region_2_inds(y_inds, sz_2(1), maxoff(1), bounds);
x_inds = cap_region_2_inds(x_inds, sz_2(2), maxoff(2), bounds);
z_inds = cap_region_2_inds(z_inds, sz_2(3), maxoff(3), bounds);


% if empty again, just return none shift
if isempty(z_inds) || isempty(x_inds) || isempty(y_inds)  
    fprintf('The non-empty area of the region 2 is too small, skip the xcorr computing, set relative_shift as [0, 0, 0]\n. Save results ...\n');
    relative_shift = [0, 0, 0];
    max_xcorr = 1;
    
    uuid = get_uuid();
    xcorrTmppath = sprintf('%s_%s.mat', xcorrFullpath(1 : end - 4), uuid);
    save('-v7.3', xcorrTmppath, 'relative_shift', 'max_xcorr')
    fileattrib(xcorrTmppath, '+w', 'g');
    movefile(xcorrTmppath, xcorrFullpath);

    fprintf('Done!\n');
    return;
end

region_2 = region_2(y_inds, x_inds, z_inds);
% src2 = [find(x_inds, 1, 'first'); find(y_inds, 1, 'first'); find(z_inds, 1, 'first')];
src2 = [x_inds(1); y_inds(1); z_inds(1)];
    
% set lower and upper bound of the maxShifts
sp2_down = (sr1 - s1) ./ downSample(:) - (src2 - 1);
maxoff = [maxoff_xy, maxoff_xy, maxoff_z] ./ downSample;
maxShifts_lb = -sp2_down' - maxoff;
maxShifts_ub = -sp2_down' + maxoff;
maxShifts = [maxShifts_lb; maxShifts_ub];
maxShifts = maxShifts(:, [2, 1, 3]);
[offset_yxz, max_xcorr, C] = normxcorr3_max_shift(single(region_2), single(region_1), maxShifts);
% [offset_yxz, max_xcorr, C] = normxcorr3_max_shift_crop(single(region_2), single(region_1), maxShifts);
sp2 = (sr1 - s1) - (src2 - 1) .* downSample(:);
relative_shift = offset_yxz([2, 1, 3]) .* downSample([2, 1, 3]) + sp2'; % - [maxoff_xy, maxoff_xy, maxoff_z];
relative_shift = relative_shift .* ((s1 >= s2)' - 0.5) * 2;

fprintf('Save results ...\n');

uuid = get_uuid();
xcorrTmppath = sprintf('%s_%s.mat', xcorrFullpath(1 : end - 4), uuid);
save('-v7.3', xcorrTmppath, 'relative_shift', 'max_xcorr')
fileattrib(xcorrTmppath, '+w', 'g');
movefile(xcorrTmppath, xcorrFullpath);

fprintf('Done!\n');

end


