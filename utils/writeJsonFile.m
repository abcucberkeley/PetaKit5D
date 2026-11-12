function [] = writeJsonFile(pr, filename)
% wrapper to write a Json file for struct variable
%
% Author: Xiongtao Ruan (05/01/2025)


s = jsonencode(pr, PrettyPrint=true);
fid = fopen(filename, 'w');
% the text as it is (as the format, a % in a value cut the file off there
% and the doubled backslashes of Windows paths became single ones)
fprintf(fid, '%s', s);
fclose(fid);

end
