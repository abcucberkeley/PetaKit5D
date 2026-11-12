function out = pkpy_f02_zarr_chunks(P)
% XR_zarrChunkCoordinatesExtraction for each parameter row
% [vol(3), batch(3), block(3), border(3), sameBatchSize, useBbox, bbox(6)];
% per row: [row; nBatch; batch(:); region(:); local(:)], or [row; -1] on an error
out = cell(size(P, 1), 1);
for i = 1:size(P, 1)
    p = P(i, :);
    bbox = [];
    if p(14)
        bbox = p(15:20);
    end
    try
        [b, r, l] = XR_zarrChunkCoordinatesExtraction(p(1:3), 'batchSize', p(4:6), 'BlockSize', p(7:9), ...
            'BorderSize', p(10:12), 'sameBatchSize', logical(p(13)), 'bbox', bbox);
        out{i} = [i; size(b, 1); b(:); r(:); l(:)];
    catch
        out{i} = [i; -1];
    end
end
out = cat(1, out{:});
end
