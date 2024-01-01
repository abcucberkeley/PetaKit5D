function [relative_shift, max_xcorr] = cross_correlation_stitching_matching(imgFullpath_1, imgFullpath_2, xcorrFullpath, cuboid_1, cuboid_2, cuboid_overlap_12, px, xyz_factors, varargin)
% decide shift between a pair of tiles in stitching based on
% cross-correlation
% 
% 
% Author: Xiongtao Ruan
% 
% xruan (04/08/2020): add option for downsampling of xcorr calculation.
% xruan (08/14/2020): change input as filenames and also save the result in file
% xruan (08/28/2020): fix issue for different image size for the tiles


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
ip.addParameter('downSample', [2, 2, 1], @isnumeric);
ip.addParameter('xy_offset_mat', -10 : 2 : 10, @isnumeric);
ip.addParameter('z_offset_mat', -5 : 5, @isnumeric);
% ip.addParameter('xyz_factor', 0.5, @isnumeric);

ip.parse(imgFullpath_1, imgFullpath_2, xcorrFullpath, cuboid_1, cuboid_2, cuboid_overlap_12, px, xyz_factors, varargin{:});

pr = ip.Results;
downSample = pr.downSample;
xy_offset_mat = pr.xy_offset_mat;
z_offset_mat = pr.z_offset_mat;

% xy_offset_mat = -10 : 2 : 10;
% z_offset_mat = -3 : 3;
% xy_offset_mat = -5 : 5;
% z_offset_mat = -3 : 3;

if exist(xcorrFullpath, 'file')
    fprintf('The xcorr result %s already exists!\n', xcorrFullpath);
    return;
end

fprintf('Compute xcorr shifts for tiles:\n    %s\n    %s\n', imgFullpath_1, imgFullpath_2);
fprintf('Load the tiles ...\n');
tic
img_1 = double(readtiff(imgFullpath_1));
toc
% img_1(img_1 < prctile(img_1(:), 80)) = nan;
tic
img_2 = double(readtiff(imgFullpath_2));
toc
% img_2(img_2 < prctile(img_2(:), 80)) = nan;

img_1(img_1 == 0) = nan;
img_2(img_2 == 0) = nan;

xy_num_off = numel(xy_offset_mat);
z_num_off = numel(z_offset_mat);
xcorr_mat = zeros(xy_num_off, xy_num_off, z_num_off);

[ny_1, nx_1, nz_1] = size(img_1);
[ny_2, nx_2, nz_2] = size(img_2);

s1 = round((cuboid_overlap_12(:, 1) - cuboid_1(:, 1)) ./ (px * xyz_factors)) + 1;
s2 = round((cuboid_overlap_12(:, 1) - cuboid_2(:, 1)) ./ (px * xyz_factors)) + 1;

t1 = round((cuboid_overlap_12(:, 2) - cuboid_1(:, 1)) ./ (px * xyz_factors)) + 1;
t2 = round((cuboid_overlap_12(:, 2) - cuboid_2(:, 1)) ./ (px * xyz_factors)) + 1;

for i = 1 : xy_num_off
    xoff = xy_offset_mat(i);
    [xa_1, xa_2] = calc_cross_correlation_indices([s1(1), t1(1)], [s2(1), t2(1)], nx_1, nx_2, xoff, downSample(1));
    
    for j = 1 : xy_num_off
        yoff = xy_offset_mat(j);
        [ya_1, ya_2] = calc_cross_correlation_indices([s1(2), t1(2)], [s2(2), t2(2)], ny_1, ny_2, yoff, downSample(2));

        for k = 1 : z_num_off
            zoff = z_offset_mat(k);
            [za_1, za_2] = calc_cross_correlation_indices([s1(3), t1(3)], [s2(3), t2(3)], nz_1, nz_2, zoff, downSample(3));
            
            frame_1_xyz = img_1(ya_1, xa_1, za_1);
            frame_2_xyz = img_2(ya_2, xa_2, za_2);

            nonnan_inds = ~isnan(frame_1_xyz + frame_2_xyz);

            inds = nonnan_inds;
            if ~any(inds, 'all')
                continue;
            end

            xcorr_xyz = corr(frame_1_xyz(inds) , frame_2_xyz(inds));
            xcorr_mat(j, i, k) = xcorr_xyz;
        end
    end
end

BW = imregionalmax(xcorr_mat);
max_xcorr = max(xcorr_mat, [], 'all');
% remove boundary values
% BW(:, :, [1, end]) = false;
% BW(:, [1, end], :) = false;
% BW([1, end], :, :) = false;

[lmy, lmx, lmz] = ind2sub(size(BW), find(BW));

% if there are multiple ones, use the one mostly close the center
% (conservative correction)
if numel(lmz) > 1
    [~, min_ind] = min((z_offset_mat(lmz) * xyz_factors(3)) .^ 2 + ...
        (xy_offset_mat(lmx) * xyz_factors(3)) .^ 2 + (xy_offset_mat(lmy) * xyz_factors(3)) .^ 2);
    lmy = lmy(min_ind);      
    lmx = lmx(min_ind);      
    lmz = lmz(min_ind);
end

if ~isempty(lmz)
    yoff = xy_offset_mat(lmy);    
    xoff = xy_offset_mat(lmx);    
    zoff = z_offset_mat(lmz);
else
    yoff = 0;
    xoff = 0;
    zoff = 0;
end

relative_shift = [xoff, yoff, zoff];

fprintf('Save results ...\n');

uuid = get_uuid();
xcorrTmppath = sprintf('%s_%s.mat', xcorrFullpath(1 : end - 4), uuid);
save('-v7.3', xcorrTmppath, 'relative_shift', 'max_xcorr')
fileattrib(xcorrTmppath, '+w', 'g');
movefile(xcorrTmppath, xcorrFullpath);

fprintf('Done!\n');

end


function [inds_1, inds_2] = calc_cross_correlation_indices(st_1, st_2, nsize_1, nsize_2, n_shift, dsf)
% determine the indices of the two arrays if there is a shift of n_shift.
% The idea is to find the longest mapping array for given mappings with a
% shift. 
%
% st_1 start and end indice of the mapping in array 1
% 

if nargin < 5
    dsf = 1;
end

inds_1 = min(1, st_1(1) + n_shift) : max(st_1(2) + n_shift, nsize_1);
inds_2 = min(1, st_2(1) - n_shift) : max(st_2(2) - n_shift, nsize_2);

% decide which way generates longer mapping, and change the other array.
if sum(inds_1 >= 1 & inds_1 <= nsize_1) >= sum(inds_2 >= 1 & inds_2 <= nsize_2)
    inds_2 = inds_1 - st_1(1) + st_2(1) - n_shift;  
else
    inds_1 = inds_2 - st_2(1) + st_1(1) + n_shift;  
end

kinds = inds_1 >=1 & inds_1 <= nsize_1 & inds_2 >=1 & inds_2 <= nsize_2;
inds_1 = inds_1(kinds);
inds_2 = inds_2(kinds);

if dsf > 1
    inds_1 = inds_1(1 : dsf : numel(inds_1));
    inds_2 = inds_2(1 : dsf : numel(inds_2));
end
    
end
