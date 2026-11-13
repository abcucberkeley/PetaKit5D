function out = s03_frame_geometry_excerpt(in)
% the geometry lines of XR_stitching_frame_zarr_dev_v1 (FRAME_CODE of tests/stitch/test_geometry.py)

P = in{1};
res = [];
for kk = 1 : numel(P)
p = P{kk};
coordinates = p.coordinates; imSizes0 = p.imSizes; tileIdx = p.tileIdx;
flippedTile = p.flippedTile; axisOrder = p.axisOrder; dataOrder = p.dataOrder;
xyPixelSize = p.xyPixelSize; dz = p.dz; SkewAngle = p.SkewAngle;
objectiveScan = p.objectiveScan; IOScan = p.IOScan; DS = p.DS; DSR = p.DSR;
xcorrShift = p.xcorrShift; isPrimaryCh = p.isPrimaryCh; BlendMethod = p.BlendMethod;
overlapType = p.overlapType; halfOrder = p.halfOrder; stitchPadSize = p.stitchPadSize;
stitchMIP = p.stitchMIP; outBbox = p.outBbox; batchSize = p.batchSize; blockSize = p.blockSize;
nF = size(coordinates, 1);
if IOScan
	DSR = false;
end
axisOrder = strrep(axisOrder, ' ', '');
order_sign_mat = zeros(2, 3); % first row: order for permute; second row: sign
axisOrder_pure = strrep(strrep(axisOrder, ',', ''), '-', '');
xyz_str = 'xyz';
for i = 1 : 3
	[~, order_sign_mat(1, :)] = sort(axisOrder_pure);
	order_sign_mat(2, i) = 1 - 2 * contains(axisOrder, ['-', xyz_str(i)]);
end
dataOrder_pure = strrep(strrep(dataOrder, ',', ''), '-', '');
[~, data_order_mat] = sort(dataOrder_pure);
xyz = coordinates(:, order_sign_mat(1, :)) .* order_sign_mat(2, :);
xyz = xyz - min(xyz, [], 1);
xyz_orig = xyz;
sinds = (1 : nF)';
if ~isempty(tileIdx)
	[tileIdx, sinds] = sortrows(tileIdx, [4, 3, 2, 1]);
	xyz = xyz(sinds, :);
elseif isempty(tileIdx)
	tileIdx = zeros(nF, 4);
end
tileNum = [numel(unique(tileIdx(:,1))), numel(unique(tileIdx(:,2))), numel(unique(tileIdx(:,3)))];
if isempty(tileNum)
	tileNum = zeros(1, nF);
end
tileNum = tileNum(order_sign_mat(1, :));
if objectiveScan || IOScan
	zAniso = dz/xyPixelSize;
else
	zAniso = sind(SkewAngle)*dz/xyPixelSize;
end
theta = SkewAngle * pi/180;
resample_type = p.resampleType;
resampleFactor = p.resampleFactor;
if ~isempty(resampleFactor)
	resample_type = 'given';
end
switch resample_type
	case 'given'
		xf = resampleFactor(1);
		yf = resampleFactor(2);
		zf = resampleFactor(3);
	case 'method_1'  %% old one
		zf = cot(abs(theta));
		yf = 1;
		xf = 1 / cos(abs(theta));
	case 'isotropic'  %% isotropic
		zf = 1;
		yf = 1;
		xf = 1;
	case 'xy_isotropic' %% x, y isotropic and z: r_z / r_x
		zf = sqrt((sin(theta) ^ 2 + zAniso ^ 2 * cos(theta) ^ 2) / (cos(theta) ^ 2 + zAniso ^ 2 * sin(theta) ^ 2));
		yf = 1;
		xf = 1;
end
if ~DSR
	zf = zAniso * zf;
end
xp = xyPixelSize * xf;
yp = xyPixelSize * yf;
zp = xyPixelSize * zf;
% the tiles as read (the input tiles' sizes are those of the processed ones here)
imSizes = imSizes0(sinds, :);
for f = 1 : nF
	imSize = imSizes(f, :);
	if imSize(data_order_mat == 3) > 1
		stitch2D = false;
		break;
	end
	stitch2D = true;
end
if all(imSizes(:, data_order_mat(3)) == 1)
	stitch2D = true;
end
xyz = stitch_process_coordinates(xyz, imSizes, flippedTile, order_sign_mat, data_order_mat, xp, yp, zp, theta, IOScan, objectiveScan, DS, DSR);
[~, data_order_reverse_mat] = sort(data_order_mat);
cuboids = [xyz, xyz + (imSizes(:, data_order_mat) - 1) .* [xp, yp, zp]];
cuboids0 = cuboids;
overlap_matrix = false(nF);
overlap_regions = zeros(nF * (nF - 1) / 2, 6);
for i = 1 : nF - 1
	cuboid_i = cuboids(i, :);
	cuboid_js = cuboids(i + 1 : end, :);
	[is_overlap_mat, cuboid_overlap_mat] = check_cuboids_overlaps(cuboid_i, cuboid_js, stitch2D);
	for j = i + 1 : nF
		if is_overlap_mat(j - i)
			overlap_matrix(i, j) = true;
			ind = (i - 1) * nF - i * (i + 1) / 2 + j;
			overlap_regions(ind, :) = cuboid_overlap_mat(j - i, :);
		end
	end
end
overlap_matrix0 = overlap_matrix;
xyz_shift = xyz;
if xcorrShift && isPrimaryCh
	xyz_shift = xyz + p.shiftVox .* [xp, yp, zp];
elseif ~isPrimaryCh
	if ~all(p.pOverlap == overlap_matrix, 'all')
		overlap_matrix = overlap_matrix & p.pOverlap;
	end
	dxyz_shift = p.shiftVox .* [xp, yp, zp];
	xyz_shift = xyz + dxyz_shift;
	pImSz = p.pImSz;
end
if isempty(overlapType)
	if strcmp(BlendMethod, 'none')
		overlapType = 'none';
	else
		overlapType = 'half';
	end
end
cuboids = [xyz_shift, xyz_shift + (imSizes(:, data_order_mat) - 1) .* [xp, yp, zp]];
half_ol_region_cell = cell(nF * 8, 1);
ol_region_cell = cell(nF * 8, 1);
overlap_map_mat = zeros(nF * 8, 2);
counter = 1;
for i = 1 : nF - 1
	cuboid_i = cuboids(i, :);
	js = find(overlap_matrix(i, i + 1 : end)) + i;
	if isempty(js)
		continue;
	end
	cuboid_js = cuboids(js, :);
	[is_overlap_mat] = check_cuboids_overlaps(cuboid_i, cuboid_js, stitch2D);
	for j1 = 1 : numel(js)
		j = js(j1);
		if ~is_overlap_mat(j1)
			overlap_matrix(i, j) = false;
			continue;
		end
		cuboid_j = cuboids(j, :);
		ind = (j - 1) * nF + i;
		[mregion_1, mregion_2] = compute_half_of_overlap_region(cuboid_i, cuboid_j, ...
			[xp, yp, zp]', 'overlapType', overlapType, 'halfOrder', halfOrder, 'stitch2D', stitch2D);
		half_ol_region_cell{counter} = {mregion_1, mregion_2};
		[mregion_11, mregion_21] = compute_half_of_overlap_region(cuboid_i, cuboid_j, ...
			[xp, yp, zp]', 'overlapType', 'zero', 'halfOrder', halfOrder, 'stitch2D', stitch2D);
		ol_region_cell{counter} = {mregion_11, mregion_21};
		overlap_map_mat(counter, :) = [counter, ind];
		counter = counter + 1;
	end
end
half_ol_region_cell(counter : end) = [];
ol_region_cell(counter : end) = [];
overlap_map_mat(counter : end, :) = [];
if nF > 1 && xcorrShift
	pad_size = [xp, yp, zp] .* stitchPadSize(data_order_mat) .* tileNum;
else
	pad_size = [0, 0, 0];
end
xyz_shift_orig = xyz_shift;
origin = min(xyz_shift_orig, [], 1) - pad_size;
xyz_shift = xyz_shift_orig - origin;
if isPrimaryCh
	dxyz = max(xyz_shift_orig + imSizes(:, data_order_mat) .* [xp, yp, zp]) + pad_size - origin;
	vxyz = round(dxyz ./ [xp, yp, zp] + 2);
	nvSize = vxyz(data_order_reverse_mat);
	if any(stitchMIP) || max(imSizes(:, data_order_mat(3)) == 1)
		nvSize(data_order_mat(3)) = 1;
	end
else
	nvSize = pImSz;
end
if ~isempty(outBbox)
	bbox = outBbox;
	if any(isinf(bbox(4 : 6)))
		stchSz = nvSize;
		bbox_end = bbox(4 : 6);
		bbox_end(isinf(bbox_end)) = stchSz(isinf(bbox_end));
		bbox(4 : 6) = bbox_end;
	end
	pad_size_l = pad_size - (bbox(data_order_mat) - 1) .* [xp, yp, zp];
	origin = min(xyz_shift_orig, [], 1) - pad_size_l;
	xyz_shift = xyz_shift_orig - origin;
	nvSize = bbox(4 : 6) - bbox(1 : 3) + 1;
end
int_xyz_shift = round(xyz_shift ./ [xp, yp, zp]);
batchSize = min(nvSize, batchSize);
batchSize = min(median(imSizes) * 2, batchSize);
blockSize = min(batchSize, blockSize);
batchSize = ceil(batchSize ./ blockSize) .* blockSize;
bSubSz = ceil(nvSize ./ batchSize);
numBatches = prod(bSubSz);
% stitch_process_block_info.m 49-84
bboxes = zeros(nF, 6);
bboxStart_mat = zeros(nF, 3);
bboxEnd_mat = zeros(nF, 3);
zind = data_order_mat(3);
int_xyz_shift_d = int_xyz_shift(:, data_order_reverse_mat);
for i = 1 : nF
    st_idx = int_xyz_shift_d(i, :);
    imSize = imSizes(i, :);
    for j = 1 : 3
        ridx = max(1, 1 - st_idx(j)) : min(imSize(j), nvSize(j) - st_idx(j));
        idx = st_idx(j) + ridx;
        if stitch2D && zind == j
            idx = 1;
            ridx = 1;
        end
        if isempty(idx)
            bboxes(i, :) = 0;
            bboxStart_mat(i, :) = 0;
            bboxEnd_mat(i, :) = 0;
            break;
        end
        bboxes(i, j) = idx(1);
        bboxes(i, j + 3) = idx(end);
        bboxStart_mat(i, j) = ridx(1);
        bboxEnd_mat(i, j) = ridx(end);
    end
end
K = numel(half_ol_region_cell);
half = zeros(K, 12);
ol = zeros(K, 12);
for c = 1 : K
	half(c, :) = [half_ol_region_cell{c}{1}(:)', half_ol_region_cell{c}{2}(:)'];
	ol(c, :) = [ol_region_cell{c}{1}(:)', ol_region_cell{c}{2}(:)'];
end
r = struct('order_sign_mat', order_sign_mat, 'data_order_mat', data_order_mat, ...
	'xyz_orig', xyz_orig, 'sinds', sinds, 'tileNum', tileNum, 'xp', xp, 'yp', yp, 'zp', zp, ...
	'stitch2D', stitch2D, 'xyz', xyz, 'cuboids0', cuboids0, 'overlap_matrix0', overlap_matrix0, ...
	'overlap_regions', overlap_regions, 'cuboids', cuboids, 'overlap_matrix', overlap_matrix, ...
	'half', half, 'ol', ol, 'overlap_map_mat', overlap_map_mat, 'pad_size', pad_size, ...
	'origin', origin, 'xyz_shift', xyz_shift, 'nvSize', nvSize, 'int_xyz_shift', int_xyz_shift, ...
	'batchSize', batchSize, 'blockSize', blockSize, 'bSubSz', bSubSz, 'numBatches', numBatches, ...
	'bboxes', bboxes, 'bboxStart_mat', bboxStart_mat, 'bboxEnd_mat', bboxEnd_mat);
if isempty(res)
	res = r;
else
	res(end + 1) = r;
end
end
out = {res};

end
