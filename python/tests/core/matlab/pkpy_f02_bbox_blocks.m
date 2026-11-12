function out = pkpy_f02_bbox_blocks(P)
% bboxToBlocks for each row [bbox(6), bboxStart(3), blockSize(3), imageSize(3), borderSize(3)];
% per row: [row; n; blockSub(:); blockInd(:); bCoords(:); wCoords(:); bboxCoords(:); bSubSize(:)]
out = cell(size(P, 1), 1);
for i = 1:size(P, 1)
    p = P(i, :);
    try
        s = bboxToBlocks(p(1:6), p(7:9), p(10:12), p(13:15), p(16:18));
        out{i} = [i; size(s.blockSub, 1); s.blockSub(:); s.blockInd(:); s.bCoords(:); s.wCoords(:); ...
            s.bboxCoords(:); s.bSubSize(:)];
    catch
        out{i} = [i; -1];
    end
end
out = cat(1, out{:});
end
