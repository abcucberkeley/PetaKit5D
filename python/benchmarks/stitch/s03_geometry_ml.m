function s03_geometry_ml(n, repeats)
% S03 benchmark: the stitching geometry of the frame (s03_frame_geometry_excerpt: axis and
% data orders, voxel sizes, coordinates, cuboids, both overlap passes with the half-overlap
% regions, canvas, batch grid, tile boxes) for an n x n grid of DSR tiles with registration
% shifts (all zero) and feather blending, as s03_geometry.py builds it
p = s03_layout(n);
for t = 1
    s03_frame_geometry_excerpt({{p}});  % warm-up
    times = zeros(1, repeats);
    for r = 1 : repeats
        t0 = tic;
        s03_frame_geometry_excerpt({{p}});
        times(r) = toc(t0);
    end
    pk5d_bench_record(sprintf('geometry_%d', n * n), t, times);
end
end

function p = s03_layout(n)
[ix, iy] = ndgrid(0 : n - 1, 0 : n - 1);
ix = ix(:);
iy = iy(:);
nF = n * n;
coords = [ix * 48.6 + mod(ix * 7 + iy * 3, 11) / 1000, iy * 48.6 + mod(ix * 5 + iy * 2, 13) / 1000, zeros(nF, 1)];
p = struct('coordinates', coords, 'imSizes', repmat([500, 500, 200], nF, 1), ...
    'tileIdx', [ix, iy, zeros(nF, 2)], 'flippedTile', [], 'axisOrder', 'x,y,z', 'dataOrder', 'y,x,z', ...
    'xyPixelSize', 0.108, 'dz', 0.5, 'SkewAngle', 32.45, 'objectiveScan', false, 'IOScan', false, ...
    'DS', false, 'DSR', true, 'resampleType', 'xy_isotropic', 'resampleFactor', [], 'xcorrShift', true, ...
    'shiftVox', zeros(nF, 3), 'isPrimaryCh', true, 'BlendMethod', 'feather', 'overlapType', '', ...
    'halfOrder', [3, 2, 1], 'stitchPadSize', [2, 2, 2], 'stitchMIP', [], 'outBbox', [], ...
    'batchSize', [500, 500, 500], 'blockSize', [500, 500, 500], 'pImSz', [], 'pOverlap', []);
end
