function [cameraClass] = get_camera_type_from_setting_file(data, varargin)
% get camera type
% 
% Author: Xiongtao Ruan 11/08/2019
% 

ip = inputParser;
ip.addParameter('CameraTypeNames', {'Orca', 'Andor'}, @(x) iscell(x) & all(cellfun(@isstr, x)));
ip.addParameter('CameraClasses', {'SCMOS', 'EMCCD'}, @(x) iscell(x) & all(cellfun(@isstr, x)));
ip.parse(varargin{:});

source = data.source;

pattern = 'settings.txt$';

% from the current folder to the parent and higher level and find the setting file
parts = strsplit(source, filesep);
% It may include an empty character in the end because of filesep in the end
if isempty(parts{end})
    parts = parts(1 : end - 1);
end
for i = 1 : numel(parts)
    cur_path = strjoin(parts(1 : numel(parts) - i + 1), filesep);
    
    dir_info = dir(cur_path);
    dir_info(1 : 2) = [];

    filenames = {dir_info.name};
    
    setting_file_ind = cellfun(@(x) ~isempty(x), regexpi(filenames, pattern));
    
    if any(setting_file_ind)
        exp_rt = cur_path;
        setting_filename = filenames{setting_file_ind};
        break;
    end
end
    
cameraClass = '';
if isempty(setting_filename)
    return;
end

fid = fopen([exp_rt, filesep, setting_filename], 'r');
setting_str = textscan(fid, '%s', 'delimiter', '\n');
setting_str = setting_str{1};

camera_type_pattern = '^Camera Type';
camera_type_str = setting_str{cellfun(@(x) ~isempty(x), regexpi(setting_str, camera_type_pattern))};

if isempty(camera_type_str)
    return;
end

% camera_type_names = {'Orca', 'Andor'};
camera_type_names = ip.Results.CameraTypeNames;
CameraClasses = ip.Results.CameraClasses;

for i = 1 : numel(camera_type_names)
    if contains(camera_type_str, camera_type_names{i}, 'IgnoreCase', true)
        cameraClass = CameraClasses{i};
    end
end


end

