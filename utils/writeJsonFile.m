function [] = writeJsonFile(pr, filename)
% wrapper to write a Json file for struct variable
%
% Author: Xiongtao Ruan (05/01/2025)


s = jsonencode(pr, PrettyPrint=true);
fid = fopen(filename, 'w');
fprintf(fid, s);
fclose(fid);

end
