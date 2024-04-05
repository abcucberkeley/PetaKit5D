function [] = writeTextFile(text_lines, filename, batchSize)
% wrapper to write a text file for cells of strings or a single string
%
% Author: Xiongtao Ruan (04/04/2024)

if nargin < 3
    batchSize = 10000;
end

if ischar(text_lines) || (iscell(text_lines) && numel(text_lines) <= batchSize)
    fid = fopen(filename, 'w');
    if iscell(text_lines)
        text_lines = strjoin(text_lines, '\n');
    end

    fwrite(fid, text_lines);
elseif iscell(text_lines) && numel(text_lines) > batchSize
    if exist(filename, 'file')
        delete(filename)
    end
    fid = fopen(filename, 'a');
    nL = numel(text_lines);
    nB = ceil(nL / batchSize);
    for b = 1 : nB
        s = (b - 1) * batchSize + 1;
        t = min(b * batchSize, nL);
        text_line_b = sprintf('%s\n', strjoin(text_lines(s : t), '\n'));
        fwrite(fid, text_line_b);
    end
end
fclose(fid);

end
