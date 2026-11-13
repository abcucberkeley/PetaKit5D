"""Stitching geometry (task S03): axis and data orders, voxel sizes, coordinates in stitching
space, cuboids and overlaps, half-overlap regions, canvas, tile boxes, batches, the size
estimate and the major-tile cover check.

Oracle cases:

- ``elem.*``: MATLAB's scalar ``sind``, ``sin``, ``cos`` and ``cot`` of skew angles (the
  geometry's values depend on their last bits);
- ``cuboids.*``, ``half.*``, ``coords.*``, ``estimate.*``, ``cover.*``: the MATLAB functions
  themselves (``check_cuboids_overlaps``, ``compute_half_of_overlap_region``,
  ``stitch_process_coordinates``, ``stitch_volume_size_estimate``, ``check_major_tile_cover``)
  on generated inputs;
- ``frame.*``: 2,000 generated tile layouts (2D and 3D, every scan type, axis and data order,
  flipped tiles, shifts, secondaries, output boxes) through an excerpt of the geometry of
  ``XR_stitching_frame_zarr_dev_v1`` (its lines, with the registration's shifts given) and of
  the tile boxes of ``stitch_process_block_info``: every derived integer and box;
- ``framerun.*``: the stitching function itself on small tiles (shifts from a prepared xcorr
  file), its stitch information and block information files read back: checks the excerpt.

Unit tests port the geometry assertions of StitchingTest (canvas convention, placements,
skewed space, the overlap split) and the size-estimate test in intent.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from petakit5d.core import axes as ax
from petakit5d.core import bbox
from petakit5d.stitch import geometry as g
from petakit5d.testing import oracle
from petakit5d.testing.cases import Case, In, Value

# ------------------------------------------------------------------------------------------
# generated layouts (MATLAB's form; integers and k / 1000 only, the same bytes everywhere)
# ------------------------------------------------------------------------------------------

SCANS = ("dsr", "dsr", "dsr", "skewed", "skewed", "ds", "objective", "io")
AXIS_ORDERS = (
    "x,y,z",
    "x,y,z",
    "-x,y,z",
    "-x,y,z",
    "y,x,z",
    "x,-y,z",
    "-x,-y,z",
    "y,-x,z",
    "x,y,-z",
    " -x, y, z",
    "-y,x,z",
    "z,y,x",
)
DATA_ORDERS = ("y,x,z",) * 6 + ("x,y,z", "z,y,x", "y,z,x", "z,x,y", "x,z,y")
SKEWS = (32.45, 32.45, 31.5, 30.0, 45.0, 32.8, -32.45)
RESAMPLES = ((1, 1, 1), (2, 2, 1), (1, 1, 2), (0.5, 0.5, 1), (1.5, 1.5, 1.5), (1, 2, 3))
HALF_ORDERS = ((3, 2, 1), (3, 2, 1), (2, 3, 1), (1, 2, 3), (2, 1, 3))
SIZES = ((500, 500, 500),) * 3 + ((64, 64, 32), (100, 80, 50), (32, 48, 16), (40, 40, 40))


def _dom(data_order: str) -> list[int]:
    """MATLAB's data_order_mat (1-based dimensions of x, y, z)"""
    pure = data_order.replace(",", "").replace("-", "")
    return [k + 1 for k in sorted(range(3), key=lambda k: pure[k])]


def make_layout(seed: int, *, uniform: bool = False) -> dict:
    """a random tile layout in MATLAB's form (sizes [d1 d2 d3], stage x, y, z in um)"""
    r = np.random.default_rng(seed)

    def pick(seq):
        return seq[int(r.integers(len(seq)))]

    def chance(k: int) -> bool:
        return int(r.integers(k)) == 0

    scan = pick(SCANS)
    two_d = chance(8)
    data_order = "y,x,z" if two_d and chance(2) else pick(DATA_ORDERS)
    dom = _dom(data_order)
    if chance(14):
        gx = gy = gz = 1
    else:
        gx, gy = int(r.integers(1, 6)), int(r.integers(1, 5))
        gz = 1 if two_d else int(r.integers(1, 3))
    n = gx * gy * gz
    sxyz = [
        int(r.integers(12, 90)),
        int(r.integers(12, 90)),
        1 if two_d else int(r.integers(5, 40)),
    ]
    msize = [0, 0, 0]
    for k in range(3):
        msize[dom[k] - 1] = sxyz[k]
    sizes = np.tile(np.array(msize, dtype=np.float64), (n, 1))
    if not uniform and n > 1 and chance(4):
        t, d = int(r.integers(n)), int(r.integers(3))
        if sizes[t, d] > 5:
            sizes[t, d] -= int(r.integers(1, 5))
    xy = pick((0.108, 0.104, 0.147, 0.1)) if chance(2) else int(r.integers(50, 300)) / 1000
    dz = pick((0.5, 0.3, 0.25, 0.2, 1.0)) if chance(2) else int(r.integers(100, 1000)) / 1000
    skew = pick(SKEWS) if chance(2) else int(r.integers(2000, 4500)) / 100
    # stage steps (um, 3 decimals): a tile's extent less an overlap, with jitter
    steps = []
    for k, px in enumerate((xy, xy, dz)):
        ov = int(r.integers(5, 45)) / 100
        steps.append(int(sxyz[k] * px * (1 - ov) * 1000))
    idx, coords = [], []
    for iz in range(gz):
        for iy in range(gy):
            for ix in range(gx):
                idx.append((ix, iy, iz, 0))
                c = [ix * steps[0], iy * steps[1], iz * steps[2]]
                c = [v + int(r.integers(-400, 401)) for v in c]
                coords.append([v / 1000 for v in c])
    tidx = np.array(idx, dtype=np.float64)
    coords = np.array(coords, dtype=np.float64)
    if n > 1 and chance(4):
        perm = r.permutation(n)
        tidx, coords, sizes = tidx[perm], coords[perm], sizes[perm]
    if n > 1 and chance(10):
        tidx[:, 3] = r.integers(0, 2, n)
    p = {
        "coordinates": coords,
        "imSizes": sizes,
        "tileIdx": None if chance(6) else tidx,
        "flippedTile": r.integers(0, 2, n).astype(np.float64).reshape(-1, 1) if chance(4) else None,
        "axisOrder": pick(AXIS_ORDERS),
        "dataOrder": data_order,
        "xyPixelSize": xy,
        "dz": dz,
        "SkewAngle": skew,
        "objectiveScan": scan == "objective",
        "IOScan": scan == "io",
        "DS": scan == "ds",
        "DSR": scan in ("dsr", "io") and not (scan == "io" and chance(2)),
        "resampleType": pick(("xy_isotropic", "xy_isotropic", "isotropic", "method_1")),
        "resampleFactor": list(pick(RESAMPLES)) if chance(5) else None,
        "xcorrShift": chance(2),
        "shiftVox": r.integers(-4, 5, (n, 3)).astype(np.float64),
        "isPrimaryCh": not chance(8),
        "BlendMethod": pick(("none", "feather", "max", "mean")),
        "overlapType": pick(("", "", "", "none", "half", "full", "zero")),
        "halfOrder": list(pick(HALF_ORDERS)),
        "stitchPadSize": [2, 2, 2] if chance(2) else [int(v) for v in r.integers(0, 5, 3)],
        "stitchMIP": [False, False, True] if chance(10) else None,
        "outBbox": None,
        "batchSize": list(pick(SIZES)),
        "blockSize": list(pick(SIZES)),
    }
    if not p["isPrimaryCh"]:
        p["pImSz"] = [float(v) for v in r.integers(20, 300, 3)]
        p["pOverlap"] = np.triu(r.integers(0, 2, (n, n)).astype(bool), 1)
    else:
        p["pImSz"] = None
        p["pOverlap"] = None
    if chance(6):
        st = [int(v) for v in r.integers(1, 6, 3)]
        en = [float(s + int(r.integers(5, 60))) if not chance(3) else math.inf for s in st]
        p["outBbox"] = [float(v) for v in st] + en
    return p


# ------------------------------------------------------------------------------------------
# the frame's geometry: MATLAB's lines and the port
# ------------------------------------------------------------------------------------------

# the lines of XR_stitching_frame_zarr_dev_v1.m (c49fc95) from 'if IOScan' (240-242), 256-351,
# 410-420 (the 2D check on the input tiles), 470-473, 529-545, 568-585 (a secondary: the
# primary's overlaps and shifts), 588-692, and of stitch_process_block_info.m 49-84; the
# registration's xyz_shift is xyz + shiftVox .* [xp, yp, zp]
FRAME_CODE = r"""
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
pattern = '^(-?x,?-?y,?-?z|-?y,?-?x,?-?z|-?z,?-?y,?-?x|-?x,?-?z,?-?y|-?z,?-?x,?-?y|-?y,?-?z,?-?x)$';
if isempty(regexpi(axisOrder, pattern, 'once'))
	error('XR_stitching_frame_zarr_dev_v1:axisOrder', 'The axisOrder is not right');
end
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
flippedTile_sorted = flippedTile;
if numel(flippedTile) == nF
	flippedTile_sorted = flippedTile(sinds);
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
		xyzf = resampleFactor(data_order_mat);
		xf = xyzf(1);
		yf = xyzf(2);
		zf = xyzf(3);
	case 'method_1'  %% old one
		zf = cot(abs(theta));
		yf = 1;
		xf = 1 / cos(abs(theta));
	case 'isotropic'  %% isotropic
		zf = 1;
		yf = 1;
		xf = 1;
	case 'xy_isotropic' %% x, y isotropic and z: r_z / r_x
		zf = sqrt((sin(theta) ^ 2 + zAniso ^ 2 * cos(theta) ^ 2) / ...
			(cos(theta) ^ 2 + zAniso ^ 2 * sin(theta) ^ 2));
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
	if imSize(data_order_mat(3)) > 1
		stitch2D = false;
		break;
	end
	stitch2D = true;
end
if all(imSizes(:, data_order_mat(3)) == 1)
	stitch2D = true;
end
xyz = stitch_process_coordinates(xyz, imSizes, flippedTile_sorted, order_sign_mat, ...
	data_order_mat, ...
	xp, yp, zp, theta, IOScan, objectiveScan, DS, DSR);
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
			[xp, yp, zp]', 'overlapType', overlapType, 'halfOrder', halfOrder, ...
			'stitch2D', stitch2D);
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
	dxyz = max(xyz_shift_orig + imSizes(:, data_order_mat) .* [xp, yp, zp], [], 1) ...
		+ pad_size - origin;
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
batchSize = min(median(imSizes, 1) * 2, batchSize);
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
"""


def _none(v):
    return None if v is None or np.size(v) == 0 else v


def py_size(m):
    """MATLAB's sizes [d1 d2 d3] (rows) as (z, y, x) rows"""
    return np.asarray(m, dtype=np.float64).reshape(-1, 3)[:, [2, 0, 1]]


def _regions(pairs) -> np.ndarray:
    rows = [
        np.concatenate(
            (g.region_to_matlab(b1).ravel(order="F"), g.region_to_matlab(b2).ravel(order="F"))
        )
        for b1, b2 in pairs
    ]
    return np.asarray(rows, dtype=np.float64).reshape(-1, 12)


def py_frame(p: dict) -> dict:
    """the port's geometry of one layout, in the form of FRAME_CODE's struct"""
    data_order = ax.axis_string_from_matlab(p["dataOrder"].replace(" ", ""))
    sizes = py_size(p["imSizes"])
    tidx = _none(p["tileIdx"])
    n = sizes.shape[0]
    sinds = g.sort_tiles(tidx)[1] if tidx is not None else np.arange(n)
    lay = g.prepare_layout(
        p["coordinates"],
        sizes,
        tile_idx=tidx,
        flipped_tile=_none(p["flippedTile"]),
        axis_order=p["axisOrder"],
        data_order=data_order,
        xy_pixel_size=p["xyPixelSize"],
        dz=p["dz"],
        skew_angle=p["SkewAngle"],
        objective_scan=p["objectiveScan"],
        io_scan=p["IOScan"],
        ds=p["DS"],
        dsr=p["DSR"],
        resample_type=p["resampleType"],
        resample_factor=None if p["resampleFactor"] is None else ax.from_yxz(p["resampleFactor"]),
        stitch_2d=g.stitch_2d_from_input_sizes(sizes[sinds], data_order),
    )
    om0 = lay.overlap_matrix.copy()
    cub0 = lay.cuboids.copy()
    psz = lay.space.pixel_sizes
    xyz_shift = lay.xyz
    p_im_sz = None
    if p["xcorrShift"] and p["isPrimaryCh"]:
        xyz_shift = lay.xyz + np.asarray(p["shiftVox"]) * psz
    elif not p["isPrimaryCh"]:
        lay.overlap_matrix = lay.overlap_matrix & np.asarray(p["pOverlap"], dtype=bool)
        xyz_shift = lay.xyz + np.asarray(p["shiftVox"]) * psz
        p_im_sz = ax.from_yxz(p["pImSz"])
    out_bbox = None if p["outBbox"] is None else bbox.box_from_matlab(p["outBbox"])
    pl = g.place_tiles(
        lay,
        xyz_shift,
        blend_method=p["BlendMethod"],
        overlap_type=p["overlapType"],
        half_order=g.half_order_from_matlab(p["halfOrder"]),
        xcorr_shift=p["xcorrShift"],
        stitch_pad_size=ax.from_yxz(p["stitchPadSize"]),
        stitch_mip=None if p["stitchMIP"] is None else ax.from_yxz(p["stitchMIP"]),
        out_bbox=out_bbox,
        p_im_sz=p_im_sz,
        batch_size=ax.from_yxz(p["batchSize"]),
        block_size=ax.from_yxz(p["blockSize"]),
    )
    perm, signs = g.parse_axis_order(p["axisOrder"])
    bxs, bst, ben = g.tile_canvas_boxes_to_matlab(pl.canvas_boxes, pl.tile_boxes, pl.tile_valid)
    cub = g.tile_cuboids(xyz_shift, lay.im_sizes, data_order, lay.space)
    return {
        "order_sign_mat": np.array([[k + 1 for k in perm], list(signs)], dtype=np.float64),
        "data_order_mat": np.array(
            g.data_order_axes_to_matlab(g.parse_data_order(data_order)), dtype=np.float64
        ),
        "xyz_orig": lay.xyz_orig,
        "sinds": lay.sinds.astype(np.float64) + 1,
        "tileNum": lay.tile_num,
        "xp": lay.space.xp,
        "yp": lay.space.yp,
        "zp": lay.space.zp,
        "stitch2D": lay.stitch_2d,
        "xyz": lay.xyz,
        "cuboids0": cub0,
        "overlap_matrix0": om0,
        "overlap_regions": lay.overlap_regions,
        "cuboids": cub,
        "overlap_matrix": pl.overlap_matrix,
        "half": _regions(pl.half_regions),
        "ol": _regions(pl.full_regions),
        "overlap_map_mat": g.overlap_map_to_matlab(pl.pairs, n),
        "pad_size": pl.pad_size,
        "origin": pl.origin,
        "xyz_shift": pl.xyz_shift,
        "nvSize": np.asarray(ax.to_yxz(pl.nv_size), dtype=np.float64),
        "int_xyz_shift": pl.int_xyz_shift,
        "batchSize": np.asarray(ax.to_yxz(pl.batch_size), dtype=np.float64),
        "blockSize": np.asarray(ax.to_yxz(pl.block_size), dtype=np.float64),
        "bSubSz": np.asarray(ax.to_yxz(pl.b_sub_sz), dtype=np.float64),
        "numBatches": float(pl.num_batches),
        "bboxes": bxs,
        "bboxStart_mat": bst,
        "bboxEnd_mat": ben,
    }


FRAME_GROUP = 50


def _frame_case(k: int) -> Case:
    seeds = range(1000 + k * FRAME_GROUP, 1000 + (k + 1) * FRAME_GROUP)

    def make():
        return [make_layout(s) for s in seeds]

    return Case(
        f"frame.{k:02d}",
        code=FRAME_CODE,
        args=[In("p")],
        python=lambda P: [py_frame(p) for p in P],
        inputs={"p": Value(make)},
    )


FRAME_CASES = [_frame_case(k) for k in range(2000 // FRAME_GROUP)]


# the stitching function itself on small tiles of ones: its stitch information and block
# information files. With registration, its xcorr file holds the shifts (the excerpt's
# coordinates plus shiftVox voxels), so that no registration runs.
_PREFIX_END = "[~, data_order_reverse_mat] = sort(data_order_mat);"
FRAMERUN_CODE = (
    r"""
p0 = in{1};
P = {p0};
tdir = [work 'tiles/'];
mkdir(tdir);
fns = cell(size(p0.coordinates, 1), 1);
for k = 1 : numel(fns)
	fns{k} = sprintf('%stile_%02d.tif', tdir, k);
	writetiff(ones(p0.imSizes(k, :), 'uint16'), fns{k});
end
rdir = [tdir 'matlab_stitch/'];
mkdir(rdir);
"""
    + FRAME_CODE[FRAME_CODE.index("res = [];") : FRAME_CODE.index(_PREFIX_END) + len(_PREFIX_END)]
    + r"""
end
if p.xcorrShift
	mkdir([rdir 'xcorr']);
	xyz_shift = xyz + p.shiftVox .* [xp, yp, zp];
	dxyz_shift = p.shiftVox .* [xp, yp, zp];
	save([rdir 'xcorr/vol.mat'], 'xyz_shift', 'dxyz_shift');
end
args = {'ResultPath', [rdir 'vol.zarr'], 'axisOrder', p.axisOrder, 'dataOrder', p.dataOrder, ...
	'flippedTile', p.flippedTile, 'dz', p.dz, 'xyPixelSize', p.xyPixelSize, ...
	'SkewAngle', p.SkewAngle, ...
	'objectiveScan', p.objectiveScan, 'IOScan', p.IOScan, 'DS', p.DS, 'DSR', p.DSR, ...
	'resampleType', p.resampleType, 'resampleFactor', p.resampleFactor, ...
	'BlendMethod', p.BlendMethod, 'halfOrder', p.halfOrder, ...
	'overlapType', p.overlapType, 'xcorrShift', p.xcorrShift, 'stitchPadSize', p.stitchPadSize, ...
	'outBbox', p.outBbox, 'batchSize', p.batchSize, 'blockSize', p.blockSize, ...
	'tileIdx', p.tileIdx, ...
	'parseCluster', false, 'jobLogDir', [work 'job_logs'], 'saveMIP', false};
evalc('XR_stitching_frame_zarr_dev_v1(fns, p.coordinates, args{:})');
si = load([rdir 'stitchInfo/vol.mat']);
d = dir([rdir 'stitchInfo/vol_block_info_task_size_*.mat']);
bi = load([d(1).folder '/' d(1).name]);
K = numel(bi.half_ol_region_cell);
half = zeros(K, 12);
ol = zeros(K, 12);
for c = 1 : K
	half(c, :) = [bi.half_ol_region_cell{c}{1}(:)', bi.half_ol_region_cell{c}{2}(:)'];
	ol(c, :) = [bi.ol_region_cell{c}{1}(:)', bi.ol_region_cell{c}{2}(:)'];
end
out = {struct('overlap_regions', si.overlap_regions, 'overlap_matrix', si.overlap_matrix, ...
	'pImSz', si.pImSz, 'pTileSizes', si.pTileSizes, 'xp', si.xp, 'yp', si.yp, 'zp', si.zp, ...
	'data_order_mat', si.data_order_mat, 'xyz_orig', si.xyz_orig, 'dxyz_shift', si.dxyz_shift, ...
	'int_xyz_shift', bi.int_xyz_shift, 'half', half, 'ol', ol, ...
	'overlap_map_mat', bi.overlap_map_mat, ...
	'bboxes', bi.bboxes, 'bboxStart_mat', bi.bboxStart_mat, 'bboxEnd_mat', bi.bboxEnd_mat, ...
	'nvSize', bi.nvSize, 'blockSize', bi.blockSize, 'bSubSz', bi.bSubSz)};
"""
)


def make_run_layout(seed: int) -> dict:
    """a layout the stitching function runs in seconds: 2 to 6 tiles of at most 28 voxels per
    axis, a primary, no resampling, no MIP and no output box"""
    for k in range(200):
        p = make_layout(20000 + 1000 * seed + k)
        n = p["coordinates"].shape[0]
        if 2 <= n <= 6 and p["isPrimaryCh"]:
            break
    sizes = np.minimum(p["imSizes"], 28)
    # the stage steps shrink with the tiles (3 decimals)
    scale = float((sizes / p["imSizes"]).min())
    p["imSizes"] = sizes
    p["coordinates"] = np.floor(p["coordinates"] * (scale * 1000)) / 1000
    p["resampleFactor"] = None
    p["stitchMIP"] = None
    p["outBbox"] = None
    return p


def py_framerun(p: dict) -> dict:
    r = py_frame(p)
    dom = [int(v) for v in r["data_order_mat"]]
    rev = [dom.index(d) for d in (1, 2, 3)]
    n = p["coordinates"].shape[0]
    sinds = (r["sinds"] - 1).astype(int)
    shift = np.asarray(p["shiftVox"]) * np.array([r["xp"], r["yp"], r["zp"]])
    return {
        "overlap_regions": r["overlap_regions"],
        "overlap_matrix": r["overlap_matrix"],
        "pImSz": r["nvSize"],
        "pTileSizes": np.asarray(p["imSizes"], dtype=np.float64)[sinds],
        "xp": r["xp"],
        "yp": r["yp"],
        "zp": r["zp"],
        "data_order_mat": r["data_order_mat"],
        "xyz_orig": r["xyz_orig"],
        "dxyz_shift": shift if p["xcorrShift"] else np.zeros((n, 3)),
        "int_xyz_shift": r["int_xyz_shift"][:, rev],
        "half": r["half"],
        "ol": r["ol"],
        "overlap_map_mat": r["overlap_map_mat"],
        "bboxes": r["bboxes"],
        "bboxStart_mat": r["bboxStart_mat"],
        "bboxEnd_mat": r["bboxEnd_mat"],
        "nvSize": r["nvSize"],
        "blockSize": r["batchSize"],
        "bSubSz": r["bSubSz"],
    }


FRAMERUN_CASES = [
    Case(
        f"framerun.{seed:02d}",
        code=FRAMERUN_CODE,
        args=[In("p")],
        python=py_framerun,
        inputs={"p": Value(lambda seed=seed: make_run_layout(seed))},
        files=False,
        marks=("slow",),
    )
    for seed in range(16)
]


def _run_base(coords, sizes, **kw) -> dict:
    """a layout for the stitching function: DSR tiles, no registration, defaults otherwise"""
    n = len(coords)
    p = {
        "coordinates": np.asarray(coords, dtype=np.float64),
        "imSizes": np.asarray(sizes, dtype=np.float64)
        .reshape(-1, 3)
        .repeat(n if np.size(sizes) == 3 else 1, axis=0),
        "tileIdx": None,
        "flippedTile": None,
        "axisOrder": "x,y,z",
        "dataOrder": "y,x,z",
        "xyPixelSize": 0.108,
        "dz": 0.3,
        "SkewAngle": 32.45,
        "objectiveScan": False,
        "IOScan": False,
        "DS": False,
        "DSR": True,
        "resampleType": "isotropic",
        "resampleFactor": None,
        "xcorrShift": False,
        "shiftVox": np.zeros((n, 3)),
        "isPrimaryCh": True,
        "BlendMethod": "none",
        "overlapType": "",
        "halfOrder": [3, 2, 1],
        "stitchPadSize": [2, 2, 2],
        "stitchMIP": None,
        "outBbox": None,
        "batchSize": [500, 500, 500],
        "blockSize": [500, 500, 500],
        "pImSz": None,
        "pOverlap": None,
    }
    p.update(kw)
    return p


# the fixes of c49fc95 that the random layouts rarely reach, on the stitching function itself
_REVERSED = np.array([[1, 0, 0, 0], [0, 0, 0, 0]], dtype=np.float64)
FRAMERUN_FIXES = {
    # flippedTile sorted with the tiles (a listing against the tileIdx order), DSR and skewed
    "flip_dsr": lambda: _run_base(
        [[1.5, 0, 0], [0, 0, 0]], [20, 24, 10], tileIdx=_REVERSED,
        flippedTile=np.array([[1.0], [0.0]]),
    ),
    "flip_skewed": lambda: _run_base(
        [[0.9, 0, 0], [0, 0, 0]], [20, 24, 16], tileIdx=_REVERSED,
        flippedTile=np.array([[1.0], [0.0]]), DSR=False,
    ),
    # 'z,x,y' tiles of one x plane are 3D (the first 2D check read the size at z's place)
    "zxy_one_x_plane": lambda: _run_base(
        [[0, 0, 0], [0, 1.5, 0]], [8, 1, 20], dataOrder="z,x,y",
        tileIdx=np.array([[0, 0, 0, 0], [0, 1, 0, 0]], dtype=np.float64),
    ),
    # one tile: the canvas is the tile plus two voxels
    "single_tile": lambda: _run_base([[0, 0, 0]], [24, 28, 10]),
    # a resample factor in the tiles' axis order ([y x z]), y and x differing
    "resample_dsr": lambda: _run_base(
        [[0, 0, 0], [1.5, 0, 0]], [20, 24, 10], resampleFactor=[2.0, 1.0, 1.0],
        tileIdx=np.array([[0, 0, 0, 0], [1, 0, 0, 0]], dtype=np.float64),
    ),
}  # fmt: skip

FRAMERUN_CASES += [
    Case(
        f"framerun.fix.{name}",
        code=FRAMERUN_CODE,
        args=[In("p")],
        python=py_framerun,
        inputs={"p": Value(make)},
        files=False,
        marks=("slow",),
    )
    for name, make in FRAMERUN_FIXES.items()
]

# the frame's axisOrder check (before the tiles are read: these tiles don't exist)
AXIS_ORDER_CODE = r"""
fns = {[work 't1.tif'], [work 't2.tif']};
evalc(['XR_stitching_frame_zarr_dev_v1(fns, [0, 0, 0; 1, 0, 0], ''axisOrder'', in{1}, ', ...
	'''parseCluster'', false, ''ResultPath'', [work ''res/out.zarr''])']);
out = {};
"""
AXIS_ORDER_CASES = [
    Case(
        f"axisorder.error.{k}",
        code=AXIS_ORDER_CODE,
        args=[ao],
        nargout=0,
        python=lambda ao: g.prepare_layout(
            np.array([[0.0, 0, 0], [1.0, 0, 0]]), np.ones((2, 3)), axis_order=ao
        ),
        expect_error="XR_stitching_frame_zarr_dev_v1:axisOrder",
        py_error=g.GeometryError,
        files=False,
    )
    for k, ao in enumerate(("x,x,z", "a,b,c", "x,y", "x,y,z,y", "x,y;z"))
]


# ------------------------------------------------------------------------------------------
# MATLAB's elementary functions of the geometry
# ------------------------------------------------------------------------------------------


def _angles() -> np.ndarray:
    a = [k / 100 for k in range(-36000, 36001, 137)] + [k / 4 for k in range(-720, 721)]
    a += [32.45, 31.5, 30.0, 45.0, 60.0, 90.0, -30.0, 32.8, 0.0, 359.95, 0.01, 1e-9]
    return np.array(a, dtype=np.float64)


def _elem_py(a: np.ndarray):
    th = a * math.pi / 180
    s = np.array([g._sind(float(v)) for v in a])
    sn = np.array([g._sin(float(v)) for v in th])
    cs = np.array([g._cos(float(v)) for v in th])
    ct = np.array([g._cot(float(v)) for v in np.abs(th)[np.abs(th) > 0]])
    return s, sn, cs, ct


ELEM_CASES = [
    Case(
        "elem.trig",
        code=(
            "a = in{1}; th = a * pi/180; ath = abs(th); ath = ath(ath > 0);\n"
            "out = {arrayfun(@sind, a), arrayfun(@sin, th), arrayfun(@cos, th), "
            "arrayfun(@cot, ath)};"
        ),
        args=[In("a")],
        nargout=4,
        python=_elem_py,
        inputs={"a": Value(_angles)},
    )
]

# ------------------------------------------------------------------------------------------
# check_cuboids_overlaps, compute_half_of_overlap_region
# ------------------------------------------------------------------------------------------


def _cuboid_rows(r, n: int, *, grid: bool) -> np.ndarray:
    """closed cuboids: starts and ends on a coarse grid (touching faces and single points)
    or k / 1000"""
    if grid:
        st = r.integers(0, 6, (n, 3)).astype(np.float64)
        ln = r.integers(0, 4, (n, 3)).astype(np.float64)
    else:
        st = r.integers(0, 6000, (n, 3)) / 1000
        ln = r.integers(0, 3000, (n, 3)) / 1000
    return np.hstack((st, st + ln))


def _cuboid_inputs(seed: int):
    r = np.random.default_rng(seed)
    grid = seed % 2 == 0
    n2 = int(r.integers(1, 12))
    n1 = 1 if seed % 3 else n2
    return {
        "c1": _cuboid_rows(r, n1, grid=grid),
        "c2": _cuboid_rows(r, n2, grid=grid),
        "two_d": bool(seed % 5 == 0),
    }


def _cuboids_case(seed: int, nargout: int) -> Case:
    def py(p):
        is_ov, ov = g.check_cuboids_overlaps(p["c1"], p["c2"], p["two_d"], regions=nargout > 1)
        return (is_ov, ov) if nargout > 1 else is_ov

    return Case(
        f"cuboids.n{nargout}.{seed:03d}",
        code=(
            "p = in{1}; o = cell(1, nargout_); [o{:}] = check_cuboids_overlaps(p.c1, p.c2, "
            "p.two_d); out = o;"
        ).replace("nargout_", str(nargout)),
        args=[In("p")],
        nargout=nargout,
        python=py,
        inputs={"p": Value(lambda: _cuboid_inputs(seed))},
    )


CUBOID_ERRORS = [
    Case(
        f"cuboids.error.{name}",
        matlab="check_cuboids_overlaps",
        args=[np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64), False],
        python=lambda a, b, c: g.check_cuboids_overlaps(a, b, c),
        expect_error=True,
        py_error=g.GeometryError,
        error_message=True,
    )
    for name, a, b in (
        ("cols", np.zeros((2, 5)), np.zeros((2, 6))),
        ("counts", np.zeros((2, 6)), np.zeros((3, 6))),
        ("counts_rev", np.zeros((3, 6)), np.zeros((2, 6))),
    )
]

OVERLAP_TYPE_LIST = ("half", "none", "full", "zero")


def _half_inputs(seed: int):
    r = np.random.default_rng(seed)
    vs = [int(v) / 1000 for v in r.integers(50, 600, 3)]
    sz = r.integers(2, 40, (2, 3))
    st1 = r.integers(0, 20, 3)
    st2 = st1 + r.integers(-10, 25, 3)
    if seed % 7 == 0:  # a single plane of overlap along one axis
        a = int(r.integers(3))
        st2[a] = st1[a] + sz[0, a] - 1
    c1 = np.concatenate((st1 * vs, (st1 + sz[0] - 1) * vs))
    c2 = np.concatenate((st2 * vs, (st2 + sz[1] - 1) * vs))
    return {
        "c1": c1,
        "c2": c2,
        "vs": np.array(vs).reshape(-1, 1),
        "type": OVERLAP_TYPE_LIST[seed % 4],
        "order": list(HALF_ORDERS[seed % len(HALF_ORDERS)]),
        "two_d": bool(seed % 9 == 0),
    }


def _half_py(p):
    b1, b2 = g.compute_half_of_overlap_region(
        p["c1"],
        p["c2"],
        np.ravel(p["vs"]),
        overlap_type=p["type"],
        half_order=g.half_order_from_matlab(p["order"]),
        stitch_2d=p["two_d"],
    )
    return g.region_to_matlab(b1), g.region_to_matlab(b2)


def _half_case(seed: int) -> Case:
    p = _half_inputs(seed)
    _, ov = g.check_cuboids_overlaps(p["c1"], p["c2"], p["two_d"])
    return Case(
        f"half.{seed:03d}",
        code=(
            "p = in{1}; [a, b] = compute_half_of_overlap_region(p.c1, p.c2, p.vs, "
            "'overlapType', p.type, 'halfOrder', p.order, 'stitch2D', p.two_d); out = {a, b};"
        ),
        args=[In("p")],
        nargout=2,
        python=_half_py,
        inputs={"p": Value(lambda: _half_inputs(seed))},
        expect_error=ov is None,
        py_error=g.GeometryError,
    )


HALF_DEFAULT_CASE = Case(
    "half.default_order",
    matlab="compute_half_of_overlap_region",
    args=[
        np.array([0, 0, 0, 6.3, 4.1, 2.0]),
        np.array([5.0, 3.0, 1.0, 9.0, 7.0, 4.0]),
        np.array([0.1, 0.1, 0.2]),
    ],
    nargout=2,
    python=lambda a, b, v: tuple(
        g.region_to_matlab(r) for r in g.compute_half_of_overlap_region(a, b, v)
    ),
)

HALF_ERROR_CASE = Case(
    "half.error.type",
    matlab="compute_half_of_overlap_region",
    args=[
        np.array([0, 0, 0, 6.3, 4.1, 2.0]),
        np.array([5.0, 3.0, 1.0, 9.0, 7.0, 4.0]),
        np.array([0.1, 0.1, 0.2]),
        "overlapType",
        "quarter",
    ],
    nargout=2,
    python=lambda a, b, v, _n, t: g.compute_half_of_overlap_region(a, b, v, overlap_type=t),
    expect_error="MATLAB:UndefinedFunction",
    py_error=g.GeometryError,
)

# ------------------------------------------------------------------------------------------
# stitch_process_coordinates
# ------------------------------------------------------------------------------------------


def _coords_inputs(seed: int):
    p = make_layout(5000 + seed)
    r = np.random.default_rng(seed)
    n = p["coordinates"].shape[0]
    p["xyz"] = r.integers(0, 50000, (n, 3)) / 1000
    p["xp"], p["yp"], p["zp"] = (int(v) / 1000 for v in r.integers(40, 900, 3))
    if p["flippedTile"] is None and seed % 3 == 0:
        p["flippedTile"] = r.integers(0, 2, n).astype(np.float64).reshape(-1, 1)
    return p


COORDS_CODE = r"""
p = in{1};
axisOrder = strrep(p.axisOrder, ' ', '');
order_sign_mat = zeros(2, 3);
axisOrder_pure = strrep(strrep(axisOrder, ',', ''), '-', '');
xyz_str = 'xyz';
for i = 1 : 3
	[~, order_sign_mat(1, :)] = sort(axisOrder_pure);
	order_sign_mat(2, i) = 1 - 2 * contains(axisOrder, ['-', xyz_str(i)]);
end
[~, data_order_mat] = sort(strrep(strrep(p.dataOrder, ',', ''), '-', ''));
DSR = p.DSR && ~p.IOScan;
theta = p.SkewAngle * pi/180;
out = {stitch_process_coordinates(p.xyz, p.imSizes, p.flippedTile, order_sign_mat, ...
	data_order_mat, p.xp, p.yp, p.zp, theta, p.IOScan, p.objectiveScan, p.DS, DSR)};
"""


def _coords_py(p):
    space = g.StitchSpace(
        p["xp"],
        p["yp"],
        p["zp"],
        1.0,
        1.0,
        1.0,
        p["SkewAngle"] * math.pi / 180,
        1.0,
        bool(p["DSR"]) and not bool(p["IOScan"]),
    )
    return g.process_coordinates(
        p["xyz"],
        py_size(p["imSizes"]),
        _none(p["flippedTile"]),
        p["axisOrder"],
        ax.axis_string_from_matlab(p["dataOrder"]),
        space,
        io_scan=p["IOScan"],
        objective_scan=p["objectiveScan"],
        ds=p["DS"],
    )


COORDS_CASES = [
    Case(
        f"coords.{seed:03d}",
        code=COORDS_CODE,
        args=[In("p")],
        python=_coords_py,
        inputs={"p": Value(lambda seed=seed: _coords_inputs(seed))},
    )
    for seed in range(120)
]

# ------------------------------------------------------------------------------------------
# stitch_volume_size_estimate
# ------------------------------------------------------------------------------------------


def _estimate_inputs(seed: int) -> dict:
    p = make_layout(9000 + seed, uniform=True)
    r = np.random.default_rng(seed)
    sz = [float(v) for v in p["imSizes"][0]]
    q = {
        "coordinates": p["coordinates"],
        "imSize": sz[:2] if sz[2] == 1 and p["dataOrder"] == "y,x,z" and seed % 2 else sz,
        "tileIdx": p["tileIdx"],
        "flippedTile": p["flippedTile"],
        "axisOrder": p["axisOrder"],
        "dataOrder": p["dataOrder"],
        "xyPixelSize": p["xyPixelSize"],
        "dz": p["dz"],
        "SkewAngle": p["SkewAngle"],
        "objectiveScan": p["objectiveScan"],
        "IOScan": p["IOScan"],
        "DS": p["DS"],
        "DSR": p["DSR"],
        "resampleType": p["resampleType"],
        "resampleFactor": p["resampleFactor"],
        "InputBbox": None,
        "tileOutBbox": None,
        "xcorrShift": p["xcorrShift"],
        "xyMaxOffset": float(r.integers(1, 400)) if seed % 3 == 0 else 300.0,
        "zMaxOffset": float(r.integers(1, 80)) if seed % 3 == 0 else 50.0,
        "stitchPadSize": p["stitchPadSize"],
        "outBbox": p["outBbox"],
        "stitchMIP": p["stitchMIP"],
        "batchSize": p["batchSize"],
        "blockSize": p["blockSize"],
    }
    for key, k in (("InputBbox", 4), ("tileOutBbox", 5)):
        if seed % k == 0 and len(q["imSize"]) == 3:
            st = [int(r.integers(1, max(1, int(v) // 3) + 1)) for v in sz]
            en = [int(r.integers(s, int(v) + 1)) for s, v in zip(st, sz, strict=True)]
            q[key] = [float(v) for v in st + en]
    return q


ESTIMATE_CODE = r"""
p = in{1};
nv = {'tileIdx', 'flippedTile', 'axisOrder', 'dataOrder', 'xyPixelSize', 'dz', 'SkewAngle', ...
	'objectiveScan', 'IOScan', 'DS', 'DSR', 'resampleType', 'resampleFactor', 'InputBbox', ...
	'tileOutBbox', 'xcorrShift', 'xyMaxOffset', 'zMaxOffset', 'stitchPadSize', 'outBbox', ...
	'stitchMIP', 'batchSize', 'blockSize'};
args = cell(1, 2 * numel(nv));
for i = 1 : numel(nv)
	args{2 * i - 1} = nv{i};
	args{2 * i} = p.(nv{i});
end
o = cell(1, 4);
[o{:}] = stitch_volume_size_estimate(p.coordinates, p.imSize, args{:});
out = o;
"""


def _box(v):
    return None if v is None else bbox.box_from_matlab(v)


def _estimate_py(p):
    nv, nb, rng, nbr = g.stitch_volume_size_estimate(
        p["coordinates"],
        ax.from_yxz(p["imSize"]),
        tile_idx=_none(p["tileIdx"]),
        flipped_tile=_none(p["flippedTile"]),
        axis_order=p["axisOrder"],
        data_order=ax.axis_string_from_matlab(p["dataOrder"]),
        xy_pixel_size=p["xyPixelSize"],
        dz=p["dz"],
        skew_angle=p["SkewAngle"],
        objective_scan=p["objectiveScan"],
        io_scan=p["IOScan"],
        ds=p["DS"],
        dsr=p["DSR"],
        resample_type=p["resampleType"],
        resample_factor=None if p["resampleFactor"] is None else ax.from_yxz(p["resampleFactor"]),
        input_bbox=_box(p["InputBbox"]),
        tile_out_bbox=_box(p["tileOutBbox"]),
        xcorr_shift=p["xcorrShift"],
        xy_max_offset=p["xyMaxOffset"],
        z_max_offset=p["zMaxOffset"],
        stitch_pad_size=ax.from_yxz(p["stitchPadSize"]),
        out_bbox=_box(p["outBbox"]),
        stitch_mip=None if p["stitchMIP"] is None else ax.from_yxz(p["stitchMIP"]),
        batch_size=ax.from_yxz(p["batchSize"]),
        block_size=ax.from_yxz(p["blockSize"]),
    )
    rows = np.array([ax.to_yxz(rng[0]), ax.to_yxz(rng[1])], dtype=np.float64)
    return np.asarray(ax.to_yxz(nv), dtype=np.float64), float(nb), rows, nbr


ESTIMATE_CASES = [
    Case(
        f"estimate.{seed:03d}",
        code=ESTIMATE_CODE,
        args=[In("p")],
        nargout=4,
        python=_estimate_py,
        inputs={"p": Value(lambda seed=seed: _estimate_inputs(seed))},
    )
    for seed in range(300)
]

# ------------------------------------------------------------------------------------------
# check_major_tile_cover
# ------------------------------------------------------------------------------------------


def _cover_inputs(seed: int):
    r = np.random.default_rng(seed)
    n = int(r.integers(2, 7))
    st = r.integers(1, 30, (n, 3))
    en = st + r.integers(0, 25, (n, 3))
    if seed % 3 == 0:  # minor tiles inside a major one
        st[1:] = np.maximum(st[1:], st[0])
        en[1:] = np.minimum(en[1:], en[0])
        en[1:] = np.maximum(en[1:], st[1:])
    major = r.integers(0, 2, n).astype(bool)
    major[int(r.integers(n))] = True
    if seed % 3 == 0:
        major[0] = True
    return {"coords": np.hstack((st, en)).astype(np.float64), "major": major}


def _cover_py(p):
    cover, unc = g.check_major_tile_cover(bbox.boxes_from_matlab(p["coords"]), p["major"])
    m = bbox.boxes_to_matlab(unc).astype(np.float64)
    m[cover] = 0
    return cover.reshape(-1, 1), m


COVER_CASES = [
    Case(
        f"cover.{seed:03d}",
        code="p = in{1}; [a, b] = check_major_tile_cover(p.coords, p.major); out = {a, b};",
        args=[In("p")],
        nargout=2,
        python=_cover_py,
        inputs={"p": Value(lambda seed=seed: _cover_inputs(seed))},
    )
    for seed in range(80)
]

CUBOID_CASES = [_cuboids_case(s, k) for s in range(60) for k in (1, 2)]
HALF_CASES = [_half_case(s) for s in range(150)]


@oracle.parametrize(
    ELEM_CASES
    + CUBOID_CASES
    + CUBOID_ERRORS
    + HALF_CASES
    + [HALF_DEFAULT_CASE, HALF_ERROR_CASE]
    + COORDS_CASES
    + ESTIMATE_CASES
    + COVER_CASES
    + FRAME_CASES
    + FRAMERUN_CASES
    + AXIS_ORDER_CASES
)
def test_matlab(case, tmp_path):
    oracle.check(case, tmp_path)


# ------------------------------------------------------------------------------------------
# StitchingTest's geometry in intent (canvas convention, placements, skewed space, splits)
# ------------------------------------------------------------------------------------------

PX = 0.108


def _grid(starts_yxz, tile_size_yxz, pixel_sizes=(PX, PX, PX), **kw):
    """StitchingTest's tiles: starts [y0 x0 z0] (voxels) at stage x, y, z = starts * pixels"""
    st = np.asarray(starts_yxz, dtype=np.float64)
    ps = np.asarray(pixel_sizes, dtype=np.float64)
    stage = st[:, [1, 0, 2]] * ps[[1, 0, 2]]
    n = st.shape[0]
    tidx = np.zeros((n, 4))
    for d in range(3):
        tidx[:, d] = np.unique(stage[:, d], return_inverse=True)[1]
    sizes = np.tile(np.asarray(ax.from_yxz(list(tile_size_yxz)), dtype=np.float64), (n, 1))
    opts = {"dsr": True, "resample_type": "isotropic", "xy_pixel_size": PX, "xcorr_shift": False}
    opts.update(kw)
    return stage, sizes, tidx, opts


def _stitch(stage, sizes, tidx, opts, **place):
    lay_kw = {k: v for k, v in opts.items() if k != "xcorr_shift"}
    lay = g.prepare_layout(stage, sizes, tile_idx=tidx, **lay_kw)
    return lay, g.place_tiles(lay, xcorr_shift=opts["xcorr_shift"], **place)


GRID = [[0, 0, 0], [36, 0, 0], [0, 48, 0], [36, 48, 0]]


def test_canvas_has_two_voxels_at_the_end_and_tiles_at_their_starts():
    # exactPositionsReproduceGroundTruth: 100 x 120 x 32 from 2 x 2 tiles of 64 x 72 x 32
    lay, pl = _stitch(*_grid(GRID, [64, 72, 32]))
    assert tuple(pl.nv_size) == (34, 102, 122)
    # the minimum corner at voxel 0; each tile at its start (z, y, x); tiles in tileIdx order
    zyx = [(s[2], s[0], s[1]) for s in GRID]
    assert list(lay.sinds) == [0, 2, 1, 3]
    for k, t in enumerate(lay.sinds):
        z, y, x = zyx[t]
        assert tuple(pl.canvas_boxes[k]) == (z, y, x, z + 32, y + 64, x + 72)
        assert tuple(pl.tile_boxes[k]) == (0, 0, 0, 32, 64, 72)
    assert pl.tile_valid.all()
    assert np.array_equal(
        pl.int_xyz_shift, np.array([[0, 0, 0], [48, 0, 0], [0, 36, 0], [48, 36, 0]])
    )


def test_canvas_with_registration_is_padded_by_two_voxels_per_tile_on_both_sides():
    # StitchingTest's canvas convention with xcorrShift: 2 * #tiles per axis on both sides
    stage, sizes, tidx, opts = _grid(GRID, [64, 72, 32], xcorr_shift=True)
    lay, pl = _stitch(stage, sizes, tidx, opts)
    pad = np.array([2 * 2, 2 * 2, 2 * 1])  # x, y, z: 2 tiles in x and y, 1 in z
    assert np.allclose(pl.pad_size, pad * PX)
    assert tuple(pl.nv_size) == (34 + 2 * 2, 102 + 2 * 4, 122 + 2 * 4)
    assert np.array_equal(pl.int_xyz_shift[0], pad)


def test_tiles_2d():
    lay, pl = _stitch(*_grid(GRID, [64, 72, 1]))
    assert lay.stitch_2d
    assert tuple(pl.nv_size) == (1, 102, 122)
    assert tuple(pl.canvas_boxes[3]) == (0, 36, 48, 1, 100, 120)


def test_io_scan_uses_stage_positions_as_is():
    starts = [[0, 0, 0], [0, 40, 0], [0, 0, 16], [0, 40, 16]]
    stage, sizes, tidx, opts = _grid(
        starts, [64, 64, 24], (PX, PX, 0.5), dsr=False, io_scan=True, dz=0.5
    )
    lay, pl = _stitch(stage, sizes, tidx, opts)
    assert tuple(pl.nv_size) == (42, 66, 106)
    assert tuple(pl.canvas_boxes[3][:3]) == (16, 0, 40)


@pytest.mark.parametrize("axis_order", ["x,y,z", "-x,y,z"])
def test_skewed_space_placement(axis_order):
    # skewedSpacePlacement: a stage x step of n * dz is n raw frames in z; stage y is y
    dz, n = 0.3, 10
    flipped = axis_order.startswith("-")
    z0 = (0, n) if flipped else (n, 0)
    starts = [[0, 0, z0[0]], [36, 0, z0[0]], [0, 0, z0[1]], [36, 0, z0[1]]]
    stage = np.array([[0, 0, 0], [0, 36 * PX, 0], [n * dz, 0, 0], [n * dz, 36 * PX, 0]])
    tidx = np.array([[0, 0, 0, 0], [0, 1, 0, 0], [1, 0, 0, 0], [1, 1, 0, 0]], dtype=np.float64)
    sizes = np.tile([40.0, 64.0, 64.0], (4, 1))
    lay = g.prepare_layout(
        stage,
        sizes,
        tile_idx=tidx,
        axis_order=axis_order,
        dz=dz,
        skew_angle=32.45,
        dsr=False,
        xy_pixel_size=PX,
        resample_type="isotropic",
    )
    pl = g.place_tiles(lay, xcorr_shift=False, blend_method="feather")
    assert tuple(pl.nv_size) == (52, 102, 66)
    for k, t in enumerate(lay.sinds):
        y, x, z = starts[t]
        assert tuple(pl.canvas_boxes[k][:3]) == (z, y, x)


def _two_tiles():
    # noneBlendTakesHalfOverlapFromEachTile: 64 x 64 x 32 tiles 40 voxels apart in x
    stage = np.array([[0, 0, 0], [40 * PX, 0, 0]])
    sizes = np.tile([32.0, 64.0, 64.0], (2, 1))
    tidx = np.array([[0, 0, 0, 0], [1, 0, 0, 0]], dtype=np.float64)
    return g.prepare_layout(stage, sizes, tile_idx=tidx, dsr=True, resample_type="isotropic")


def test_none_blend_splits_the_overlap_in_halves_along_z():
    pl = g.place_tiles(_two_tiles(), xcorr_shift=False, blend_method="none")
    assert pl.overlap_type == "none"
    (b1, b2), (f1, f2) = pl.half_regions[0], pl.full_regions[0]
    # the full overlap: x 40..63 of tile 1, 0..23 of tile 2, all of y and z
    assert f1 == (0, 0, 40, 32, 64, 64) and f2 == (0, 0, 0, 32, 64, 24)
    # z is split (halfOrder [3, 2, 1]): half of the 32 slices taken off each region
    assert b1 == (0, 0, 40, 16, 64, 64) and b2 == (16, 0, 0, 32, 64, 24)


def test_blending_takes_a_quarter_from_each_side():
    # blendingSplitsTheOverlapInQuarters: overlap type 'half' takes off round(0.75 * 32) = 24
    pl = g.place_tiles(_two_tiles(), xcorr_shift=False, blend_method="mean")
    assert pl.overlap_type == "half"
    b1, b2 = pl.half_regions[0]
    assert b1 == (0, 0, 40, 8, 64, 64) and b2 == (24, 0, 0, 32, 64, 24)


SIZE_SETUPS = {
    "grid": {"starts": GRID, "tile": [64, 72, 32], "kw": {}},
    "gridPadded": {"starts": GRID, "tile": [64, 72, 32], "kw": {"xcorr_shift": True}},
    "gridOutBbox": {
        "starts": GRID,
        "tile": [64, 72, 32],
        "kw": {},
        "out_bbox": bbox.box_from_matlab([3, 5, 2, 90, 100, 30]),
    },
    "tiles2D": {"starts": GRID, "tile": [64, 72, 1], "kw": {}},
    "ioScan": {
        "starts": [[0, 0, 0], [0, 40, 0], [0, 0, 16], [0, 40, 16]],
        "tile": [64, 64, 24],
        "ps": (PX, PX, 0.5),
        "kw": {"dsr": False, "io_scan": True, "dz": 0.5},
    },
}


@pytest.mark.parametrize("name", sorted(SIZE_SETUPS))
def test_volume_size_estimate_matches_the_stitched_volume(name):
    s = SIZE_SETUPS[name]
    stage, sizes, tidx, opts = _grid(s["starts"], s["tile"], s.get("ps", (PX, PX, PX)), **s["kw"])
    out_bbox = s.get("out_bbox")
    lay, pl = _stitch(stage, sizes, tidx, opts, out_bbox=out_bbox)
    est, nb, rng, nbr = g.stitch_volume_size_estimate(
        stage, sizes[0], tile_idx=tidx, out_bbox=out_bbox, **opts
    )
    assert np.array_equal(est, pl.nv_size)
    margin = np.array([50, 300, 300]) * 2 * (opts["xcorr_shift"] and out_bbox is None)
    margin[est == 1] = 0
    assert np.array_equal(rng, np.vstack((np.maximum(est - margin, 1), est + margin)))


# ------------------------------------------------------------------------------------------
# rules and quirks
# ------------------------------------------------------------------------------------------


def test_axis_and_data_orders():
    assert g.parse_axis_order("x,y,z") == ((0, 1, 2), (1.0, 1.0, 1.0))
    assert g.parse_axis_order(" -x, y,z") == ((0, 1, 2), (-1.0, 1.0, 1.0))
    assert g.parse_axis_order("y,-x,z") == ((1, 0, 2), (-1.0, 1.0, 1.0))
    assert g.parse_axis_order("z,x,y") == ((1, 2, 0), (1.0, 1.0, 1.0))
    with pytest.raises(g.GeometryError):
        g.parse_axis_order("x,y")
    assert g.parse_data_order("z,y,x") == (2, 1, 0)
    for m in ("y,x,z", "x,y,z", "z,y,x", "y,z,x", "z,x,y", "x,z,y"):
        py = ax.axis_string_from_matlab(m)
        assert g.data_order_axes_to_matlab(g.parse_data_order(py)) == tuple(_dom(m))
        assert g.data_order_axes_from_matlab(_dom(m)) == g.parse_data_order(py)


def test_axis_order_must_be_a_permutation_of_xyz():
    # StitchGeometryTest/axisOrderMustBeAPermutationOfXYZ: the check never raised, and 'z,x,y'
    # was missing from its pattern
    for ao in ("x,x,z", "a,b,c", "x,y", "x,y,z,y", "x,y;z"):
        with pytest.raises(g.GeometryError) as e:
            g.check_axis_order(ao)
        assert e.value.identifier == "XR_stitching_frame_zarr_dev_v1:axisOrder"
        with pytest.raises(g.GeometryError):
            g.prepare_layout(np.zeros((2, 3)), np.ones((2, 3)), axis_order=ao)
    for ao in ("x,y,z", "y,x,z", "z,y,x", "x,z,y", "z,x,y", "y,z,x", "-x,y,z", "x,-y,z",
               "-z,x,-y", "yxz", " y, x, z", "X,Y,Z"):  # fmt: skip
        g.check_axis_order(ao)


def test_half_order_conversion():
    assert g.half_order_from_matlab([3, 2, 1]) == (0, 1, 2)
    assert g.half_order_from_matlab([2, 3, 1]) == (1, 0, 2)
    assert g.half_order_to_matlab((1, 0, 2)) == (2, 3, 1)


def test_pair_index_is_matlabs_minus_one():
    n = 7
    for i in range(n - 1):
        for j in range(i + 1, n):
            i1, j1 = i + 1, j + 1
            assert g.pair_index(i, j, n) == (i1 - 1) * n - i1 * (i1 + 1) // 2 + j1 - 1


def test_region_round_trip():
    m = np.array([[3, 9], [1, 64], [5, 5]], dtype=np.float64)
    b = g.region_from_matlab(m)
    assert b == (4, 0, 2, 5, 64, 9)
    assert np.array_equal(g.region_to_matlab(b), m)


def test_first_2d_check_reads_zs_dimension():
    # StitchingTest/zxyTilesOfOneXPlaneAreNot2D: the check read the size at z's place in the
    # data order, another axis for 'y,z,x' and 'z,x,y'
    assert g.stitch_2d_from_input_sizes(np.array([[1.0, 30.0, 40.0]]), "z,y,x")
    # MATLAB 'y,z,x' is ours 'x,y,z': z is the last axis
    zyx = ax.axis_string_from_matlab("y,z,x")
    assert g.stitch_2d_from_input_sizes(np.array([[40.0, 30.0, 1.0]]), zyx)
    assert not g.stitch_2d_from_input_sizes(np.array([[40.0, 1.0, 9.0]]), zyx)
    # MATLAB 'z,x,y' (ours 'y,z,x'): one x plane is not 2D
    zxy = ax.axis_string_from_matlab("z,x,y")
    assert not g.stitch_2d_from_input_sizes(np.array([[30.0, 8.0, 1.0]]), zxy)
    assert g.stitch_2d_from_input_sizes(np.array([[30.0, 1.0, 9.0]]), zxy)


def test_single_tile_canvas_is_the_tile_plus_two_voxels():
    # StitchingTest/singleTileCanvasIsTheTilePlusTwoVoxels: max and median of one tile took the
    # row (its largest extent on every axis); now per axis, as the size estimate
    stage = np.zeros((1, 3))
    sizes = np.array([[10.0, 30.0, 40.0]])
    lay = g.prepare_layout(stage, sizes, dsr=True, resample_type="isotropic")
    pl = g.place_tiles(lay, xcorr_shift=False)
    assert tuple(pl.nv_size) == (12, 32, 42)
    est = g.stitch_volume_size_estimate(stage, sizes[0], dsr=True, resample_type="isotropic")[0]
    assert tuple(est) == (12, 32, 42)
    assert tuple(pl.batch_size) == (12, 32, 42)
    bs, bl, sub, nb = g.batch_grid([12, 32, 42], sizes, (20, 20, 20), (20, 20, 20))
    assert tuple(bs) == (12, 20, 20) and tuple(sub) == (1, 2, 3) and nb == 6


def test_batch_grid_rounds_batches_up_to_blocks():
    sizes = np.array([[30.0, 100.0, 90.0], [30.0, 100.0, 91.0]])
    bs, bl, sub, nb = g.batch_grid([70, 300, 250], sizes, (64, 128, 100), (25, 50, 40))
    # min(nv, batch) = (64, 128, 100), min(2 * median, .) = (60, 128, 100) -> blocks (25, 50, 40)
    assert tuple(bl) == (25, 50, 40)
    assert tuple(bs) == (75, 150, 120)
    assert tuple(sub) == (1, 2, 3) and nb == 6


def test_flipped_tiles_stay_with_their_tiles_when_sorted():
    # StitchingTest/flippedTilesStayWithTheirTilesWhenSorted: the flags stayed in the listed
    # order when the tiles were sorted by tileIdx
    stage = np.array([[3.0, 0, 0], [0.0, 0, 0]])
    sizes = np.tile([20.0, 30.0, 30.0], (2, 1))
    tidx = np.array([[1, 0, 0, 0], [0, 0, 0, 0]], dtype=np.float64)
    for dsr in (True, False):
        a = g.prepare_layout(stage, sizes, tile_idx=tidx, flipped_tile=[1, 0], dsr=dsr)
        b = g.prepare_layout(stage[::-1], sizes, tile_idx=tidx[::-1], flipped_tile=[0, 1], dsr=dsr)
        assert list(a.sinds) == [1, 0] and list(b.sinds) == [0, 1]
        assert np.array_equal(a.xyz, b.xyz)
        # the flag still matters
        c = g.prepare_layout(stage, sizes, tile_idx=tidx, flipped_tile=[0, 1], dsr=dsr)
        assert not np.array_equal(a.xyz, c.xyz)


def test_resample_factor_is_in_the_tiles_axis_order():
    # StitchingTest/resampleFactorIsInTheTilesAxisOrder: the voxel sizes took the factor as
    # x, y, z while the tiles are resampled with it in their axis order
    sp = g.stitch_pixel_sizes(0.1, 0.2, 30, resample_factor=(1, 2, 1), dsr=True)
    assert (sp.xf, sp.yf, sp.zf) == (1.0, 2.0, 1.0)
    # (tiles named 'y,z,x': y is axis 0, z axis 1, x axis 2)
    sp = g.stitch_pixel_sizes(0.1, 0.2, 30, resample_factor=(2, 3, 1), dsr=True, data_order="y,z,x")
    assert (sp.xf, sp.yf, sp.zf) == (1.0, 2.0, 3.0)


def test_cuboid_overlaps():
    c = np.array([[0, 0, 0, 2, 2, 2]], dtype=np.float64)
    others = np.array(
        [[1, 1, 1, 3, 3, 3], [2, 0, 0, 4, 2, 2], [3, 0, 0, 4, 1, 1]], dtype=np.float64
    )
    is_ov, ov = g.check_cuboids_overlaps(c, others)
    # a shared face is no overlap, with one output too
    assert list(is_ov) == [True, False, False]
    assert np.array_equal(ov[0], [1, 1, 1, 2, 2, 2]) and not ov[1].any()
    assert list(g.check_cuboids_overlaps(c, others, regions=False)[0]) == [True, False, False]
    # 2D: only the first two axes count for a single point
    flat = np.array([[1, 1, 2, 3, 3, 2]], dtype=np.float64)
    assert g.check_cuboids_overlaps(c, flat, True)[0][0]
    assert not g.check_cuboids_overlaps(c, flat, False)[0][0]
    assert g.check_cuboids_overlaps(c, np.array([[5.0, 5, 5, 6, 6, 6]]))[1] is None
    with pytest.raises(g.GeometryError, match="n X 6"):
        g.check_cuboids_overlaps(np.zeros((1, 5)), c)


def test_cuboids_that_only_touch_do_not_overlap():
    # StitchGeometryTest/cuboidsThatOnlyTouchDoNotOverlap
    c = np.array([0, 0, 0, 1, 1, 1], dtype=np.float64)
    js = np.array(
        [[1, 0, 0, 2, 1, 1], [0, 1, 0, 1, 2, 1], [0, 0, 1, 1, 1, 2], [1, 1, 0, 2, 2, 1],
         [1, 1, 1, 2, 2, 2], [0.5, 0.5, 0.5, 2, 2, 2], [1.5, 0, 0, 2, 1, 1]], dtype=np.float64,
    )  # fmt: skip
    expected = [False] * 5 + [True, False]
    assert list(g.check_cuboids_overlaps(c, js, False, regions=False)[0]) == expected
    is_ov, regions = g.check_cuboids_overlaps(c, js, False)
    assert list(is_ov) == expected
    assert regions[5].tolist() == [0.5, 0.5, 0.5, 1, 1, 1]
    c = np.array([0, 0, 0, 1, 1, 0], dtype=np.float64)
    js = np.array([[0.5, 0.5, 0, 2, 2, 0], [1, 0, 0, 2, 1, 0]], dtype=np.float64)
    assert list(g.check_cuboids_overlaps(c, js, True)[0]) == [True, False]
    assert list(g.check_cuboids_overlaps(c, js, False)[0]) == [False, False]


def test_major_tile_cover_box_holds_every_uncovered_piece():
    # StitchGeometryTest/majorTileCoverBoxHoldsEveryUncoveredPiece (MATLAB boxes there)
    def cover(rows):
        c, unc = g.check_major_tile_cover(bbox.boxes_from_matlab(np.array(rows)), [False, True])
        m = bbox.boxes_to_matlab(unc)
        m[c] = 0
        return c.tolist(), m[0].tolist()

    assert cover([[1, 1, 1, 10, 10, 10], [1, 4, 1, 10, 7, 10]]) == (
        [False, True],
        [1, 1, 1, 10, 10, 10],
    )
    assert cover([[1, 1, 1, 10, 10, 10], [1, 4, 1, 10, 10, 10]]) == (
        [False, True],
        [1, 1, 1, 10, 4, 10],
    )
    assert cover([[2, 2, 2, 9, 9, 9], [1, 1, 1, 10, 10, 10]]) == ([True, True], [0] * 6)


def test_cover_of_a_minor_tile():
    major = [True, False]
    inside = np.array([[0, 0, 0, 10, 10, 10], [2, 2, 2, 8, 8, 8]])
    cover, unc = g.check_major_tile_cover(inside, major)
    assert list(cover) == [True, True] and not unc.any()
    outside = np.array([[0, 0, 0, 10, 10, 10], [2, 2, 2, 8, 8, 14]])
    cover, unc = g.check_major_tile_cover(outside, major)
    assert list(cover) == [True, False]
    assert unc[0].tolist() == [0] * 6 and unc[1][2] >= 9


def test_matlab_trig():
    # fdlibm's values (MATLAB's bits are checked by elem.trig): within an ulp or two of the C
    # library's
    for x in (0.1, 0.5663, 0.785398, 1.2, 3.0, 10.0, -2.5):
        assert abs(g._tan(x) - math.tan(x)) <= 2 * math.ulp(math.tan(x))
        assert abs(g._cot(x) - 1 / math.tan(x)) <= 4 * math.ulp(1 / math.tan(x))
    assert g._sind(30.0) == 0.5 and g._sind(-150.0) == -0.5 and g._sind(90.0) == 1.0
    assert abs(g._sind(32.45) - math.sin(math.radians(32.45))) <= math.ulp(0.5)


def test_pixel_sizes():
    sp = g.stitch_pixel_sizes(0.108, 0.5, 32.45, dsr=True)
    assert sp.xp == sp.yp == 0.108 and sp.zf == pytest.approx(
        math.sqrt(
            (math.sin(sp.theta) ** 2 + sp.z_aniso**2 * math.cos(sp.theta) ** 2)
            / (math.cos(sp.theta) ** 2 + sp.z_aniso**2 * math.sin(sp.theta) ** 2)
        )
    )
    sk = g.stitch_pixel_sizes(0.108, 0.5, 32.45, resample_type="isotropic")
    assert sk.zp == pytest.approx(0.108 * math.sin(math.radians(32.45)) * 0.5 / 0.108)
    # IO scans are never DSR
    assert not g.stitch_pixel_sizes(io_scan=True, dsr=True).dsr
    # resample factors in the tiles' axis order: (z, y, x)
    rs = g.stitch_pixel_sizes(0.1, 0.2, 30, resample_factor=(3, 2, 1), dsr=True)
    assert (rs.xf, rs.yf, rs.zf) == (1.0, 2.0, 3.0)
    with pytest.raises(g.GeometryError):
        g.stitch_pixel_sizes(resample_factor=(2, 1))
    with pytest.raises(g.GeometryError):
        g.stitch_pixel_sizes(resample_type="other")


def test_frame_layouts_run():
    # the generated layouts run through the port (the oracle compares their values)
    for s in range(1000, 1040):
        r = py_frame(make_layout(s))
        assert r["numBatches"] >= 1
        assert r["half"].shape[0] == r["overlap_map_mat"].shape[0]
