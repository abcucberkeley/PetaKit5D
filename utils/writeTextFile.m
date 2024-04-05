function [] = writeTextFile(text_lines, filename)
% wrapper to write a text file for cells of strings or a single string
%
% Author: Xiongtao Ruan (04/04/2024)


fid = fopen(filename, 'w');
if iscell(text_lines)
    text_lines = strjoin(text_lines, '\n');
end
fwrite(fid, text_lines);
fclose(fid);

end