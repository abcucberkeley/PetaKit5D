function [relative_shift, max_xcorr] = cross_correlation_registration_2d(imgFullpath_1, imgFullpath_2, xcorrFullpath, cuboid_1, cuboid_2, cuboid_overlap_12, xyz_voxelsizes, data_order_mat, varargin)
% compute shift between a pair of tiles in stitching based on
% cross-correlation for 2d images.
% based on multires_cross_correlation_registration_imblock_test.m
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
% xruan (10/08/2021): accelearte maxoff search; add support to crop region
% with rich signal for large space sample. 


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imgFullpath_1', @ischar);
ip.addRequired('imgFullpath_2', @ischar);
ip.addRequired('xcorrFullpath', @ischar);
ip.addRequired('cuboid_1', @isnumeric);
ip.addRequired('cuboid_2', @isnumeric);
ip.addRequired('cuboid_overlap_12', @isnumeric);
ip.addRequired('xyz_voxelsizes', @isnumeric);
ip.addRequired('data_order_mat', @isnumeric);
ip.addParameter('downSample', [1, 1, 1], @isnumeric);
ip.addParameter('MaxOffset', [300, 300, 50], @isnumeric);
ip.addParameter('dimNumThrsh', 10000, @isnumeric);
% ip.addParameter('xyz_factor', 0.5, @isnumeric);

ip.parse(imgFullpath_1, imgFullpath_2, xcorrFullpath, cuboid_1, cuboid_2, cuboid_overlap_12, xyz_voxelsizes, data_order_mat, varargin{:});

pr = ip.Results;
downSample = pr.downSample;
MaxOffset = pr.MaxOffset;
% dimNumThrsh = pr.dimNumThrsh;

downSample = downSample(:)';
xdo = data_order_mat == 1;
ydo = data_order_mat == 2;
[~, data_order_reverse_mat] = sort(data_order_mat);
maxXOffset = MaxOffset(xdo);
maxYOffset = MaxOffset(ydo);

saveResult = ~isempty(xcorrFullpath);
if saveResult && exist(xcorrFullpath, 'file')
    fprintf('The xcorr result %s already exists!\n', xcorrFullpath);
    return;
end

fprintf('Compute xcorr shifts for tiles:\n    %s\n    %s\n', imgFullpath_1, imgFullpath_2);
% fprintf('Load overlap regions ...\n');

% img_1(img_1 == 0) = nan;
% img_2(img_2 == 0) = nan;
[~, ~, ext] = fileparts(imgFullpath_1);
switch ext
    case '.tif'
        bim_1 = blockedImage(imgFullpath_1, "Adapter", MPageTiffAdapter);
        bim_2 = blockedImage(imgFullpath_2, "Adapter", MPageTiffAdapter);
    case '.zarr'
        try
            bim_1 = blockedImage(imgFullpath_1, "Adapter", CZarrAdapter);
            bim_2 = blockedImage(imgFullpath_2, "Adapter", CZarrAdapter);
        catch ME
            disp(ME)
            bim_1 = blockedImage(imgFullpath_1, "Adapter", ZarrAdapter);
            bim_2 = blockedImage(imgFullpath_2, "Adapter", ZarrAdapter);
        end
end

sz_1 = bim_1.Size;

% note: cuboid_overlap_12, cuboid_1 and xyz_factors are in xyz order.
xyz_voxelsizes = xyz_voxelsizes(:);
s1 = round((cuboid_overlap_12(1 : 3) - cuboid_1(1 : 3))' ./ xyz_voxelsizes) + 1;
s2 = round((cuboid_overlap_12(1 : 3) - cuboid_2(1 : 3))' ./ xyz_voxelsizes) + 1;

t1 = round((cuboid_overlap_12(4 : 6) - cuboid_1(1 : 3))' ./ xyz_voxelsizes) + 1;
t2 = round((cuboid_overlap_12(4 : 6) - cuboid_2(1 : 3))' ./ xyz_voxelsizes) + 1;

region_2 = bim_2.Adapter.getIORegion(s2(data_order_reverse_mat), t2(data_order_reverse_mat));
region_2 = single(region_2);

% crop region 2 if one or more dimension is too large, if so, find a small
% region with rich signal. 
if false    
    [region_2, crop_bbox] = crop_subregion_by_intensity(region_2, dimNumThrsh);
    
    s1 = s1 + (crop_bbox([2, 1, 3]) - 1)';
    t1 = s1 + (crop_bbox([5, 4, 6]) - crop_bbox([2, 1, 3]))';
    
    s2 = s2 + (crop_bbox([2, 1, 3]) - 1)';
    % t2 = s2 + (crop_bbox([5, 4, 6]) - crop_bbox([2, 1, 3]))';    
end

% first resize to downsampled ones
if any(downSample ~= 1)
    % resize by max pooling
    region_2 = max_pooling_3d(region_2, downSample);
end

% calculate bbox for overlap regions
maxoff_x = ceil(maxXOffset ./ downSample(xdo)) * downSample(xdo);
maxoff_y = ceil(maxYOffset ./ downSample(ydo)) * downSample(ydo);
% maxoff_z = ceil(maxZOffset ./ downSample(1)) * downSample(3);

boarder = [maxoff_x; maxoff_y; 0];
sr1 = max(s1 - boarder, 1);
tr1 = min(t1 + boarder, sz_1(data_order_mat)');

% region_1 = bim_1.getRegion(sr1([2, 1, 3]), tr1([2, 1, 3]));
region_1 = bim_1.Adapter.getIORegion(sr1(data_order_reverse_mat), tr1(data_order_reverse_mat));
region_1 = single(region_1);

% first resize to downsampled ones
if any(downSample ~= 1)
    % resize by max pooling
    region_1 = max_pooling_3d(region_1, downSample);
end

% crop region 2
sz_2 = size(region_2, [1, 2, 3]);
% region_2 = region_2(s2(2) : t2(2), s2(1) : t2(1), s2(3) : t2(3));
region_2_nozeros = region_2 ~= 0;
% select the region with at least 1/4 areas in the selected dimensions.
z_inds = 1;
x_inds = sum(region_2_nozeros(:, :, z_inds), [1, 3]) > sz_2(1) * sum(z_inds) / 4;
y_inds = sum(region_2_nozeros(:, x_inds, z_inds), [2, 3]) > sum(x_inds) * sum(z_inds) / 4;

% reduce area threshold if the inds are empty
if ~any(z_inds) || ~any(x_inds) || ~any(y_inds)  
    z_inds = 1;
    x_inds = sum(region_2_nozeros(:, :, z_inds), [1, 3]) > sz_2(1) * sum(z_inds) / 8;
    y_inds = sum(region_2_nozeros(:, x_inds, z_inds), [2, 3]) > sum(x_inds) * sum(z_inds) / 8;
end

% if empty again, just return none shift
if ~any(z_inds) || ~any(x_inds) || ~any(y_inds)  
    fprintf('The non-empty area of the region 2 is too small, skip the xcorr computing, set relative_shift as [0, 0, 0]\n. Save results ...\n');
    relative_shift = [0, 0, 0];
    max_xcorr = 1;
    
    if saveResult
        uuid = get_uuid();
        xcorrTmppath = sprintf('%s_%s.mat', xcorrFullpath(1 : end - 4), uuid);
        save('-v7.3', xcorrTmppath, 'relative_shift', 'max_xcorr')
        fileattrib(xcorrTmppath, '+w', 'g');
        movefile(xcorrTmppath, xcorrFullpath);
    end

    fprintf('Done!\n');
    return;
end

if ~all(y_inds([1, end]) & x_inds([1, end]))
    region_2 = region_2(y_inds, x_inds, z_inds);
end
src2 = [find(y_inds, 1, 'first'); find(x_inds, 1, 'first'); 1];
src2 = src2(data_order_mat);
% src2 = [x_inds(1); y_inds(1); z_inds(1)];
    
% set lower and upper bound of the maxShifts
sp2_down = (sr1 - s1) ./ downSample(data_order_mat)' - (src2 - 1);
maxoff = [maxoff_x, maxoff_y, 1] ./ downSample(data_order_mat);
maxShifts_lb = -sp2_down' - maxoff;
maxShifts_ub = -sp2_down' + maxoff;
maxShifts = [maxShifts_lb; maxShifts_ub];
maxShifts = maxShifts(:, data_order_reverse_mat);
if data_order_mat(3) ~= 3
    region_1 = permute(region_1, data_order_mat);
    region_2 = permute(region_2, data_order_mat);
    maxShifts = maxShifts(:, data_order_mat);
end
[offset_yxz, max_xcorr, C] = normxcorr2_max_shift(single(region_2), single(region_1), maxShifts);
if data_order_mat(3) ~= 3
    offset_yxz = offset_yxz(data_order_reverse_mat);
end
% [offset_yxz, max_xcorr, C] = normxcorr3_max_shift_crop(single(region_2), single(region_1), maxShifts);
sp2 = (sr1 - s1) - (src2 - 1) .* downSample(data_order_mat)';
relative_shift = offset_yxz(data_order_mat) .* downSample(data_order_mat) + sp2'; % - [maxoff_xy, maxoff_xy, maxoff_z];
relative_shift = relative_shift .* ((s1 >= s2)' - 0.5) * 2;

if saveResult
    fprintf('Save results ...\n');
    
    uuid = get_uuid();
    xcorrTmppath = sprintf('%s_%s.mat', xcorrFullpath(1 : end - 4), uuid);
    save('-v7.3', xcorrTmppath, 'relative_shift', 'max_xcorr')
    fileattrib(xcorrTmppath, '+w', 'g');
    movefile(xcorrTmppath, xcorrFullpath);
end

fprintf('Done!\n');

end


