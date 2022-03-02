function [simplePath] = simplifyPath(Path)
% simplify a path: make it absolute path, for folder, there is no filesep in the end
% 
% Author: Xiongtao Ruan (03/11/2020)

dir_info = dir(Path);
if numel(dir_info) >= 2
    simplePath = dir_info(1).folder;
else
    simplePath = [dir_info.folder, filesep, dir_info.name];
end

end