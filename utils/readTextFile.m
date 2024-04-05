function [file_lines] = readTextFile(filename)
% wrapper to read a text file line by line to a cell array
%
% Author: Xiongtao Ruan (04/04/2024)

if ~exist(filename, 'file')
    error('%s does not exist, please check the path!', filename);
end

fid = fopen(filename, 'r');
file_lines = cell(10000, 1);
ind = 1;
while ~feof(fid)
    cur_line = fgetl(fid);
    file_lines{ind} = cur_line;
    ind = ind + 1;
end
fclose(fid);
file_lines(ind : end) = [];


end