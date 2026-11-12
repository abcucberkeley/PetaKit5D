function out = pkpy_f02_subvol(P)
% XR_subVolumeCoordinatesExtraction for each row [vol(3), chunk(3), overlap, maxSubVolume, useGoodNum];
% per row: [row; nn; xmin; xmax; ymin; ymax; zmin; zmax], or [row; -1] on an error
out = cell(size(P, 1), 1);
for i = 1:size(P, 1)
    p = P(i, :);
    try
        [xmin, xmax, ymin, ymax, zmin, zmax, nn] = XR_subVolumeCoordinatesExtraction(p(1:3), 'ChunkSize', p(4:6), ...
            'overlapSize', p(7), 'maxSubVolume', p(8), 'useGoodNum', p(9));
        out{i} = [i; nn; xmin; xmax; ymin; ymax; zmin; zmax];
    catch
        out{i} = [i; -1];
    end
end
out = cat(1, out{:});
end
