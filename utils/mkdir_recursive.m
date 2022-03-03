function [status] = mkdir_recursive(full_path, group_write)
% xruan 07/25/2015 
% mkdir recursively
%
% xruan (08/22/2020): add option for group write
% xruan (03/02/2022): add support for windows

if nargin < 2
    group_write = false;
end

path_list = strsplit(full_path, '/');

path_depth_i = '';
for i = 1 : length(path_list)
    if strcmp(path_list{i}, '~')
        path_depth_i = '~';
        continue;
    end
    if ispc && isempty(path_depth_i)
        path_depth_i = path_list{i};
    else
        path_depth_i = strjoin({path_depth_i, path_list{i}}, '/');
    end
    if ~exist(path_depth_i, 'dir')
        mkdir(path_depth_i);
        if group_write
            fileattrib(path_depth_i, '+w', 'g');
        end
    end

end
status = 1;

end