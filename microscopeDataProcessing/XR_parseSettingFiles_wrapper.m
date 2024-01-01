function [settingInfo] = XR_parseSettingFiles_wrapper(imageFilenames, varargin)
% parse a given list of settting files to get useful setting info for the
% image processing
% 
% Author: Xiongtao Ruan (12/05/2020)
% 
% xruan (06/16/2021): exclude partial files
% xruan (07/13/2021): fix bug for image filenames order not match for input
% and output
% xruan (01/31/2021): add support for negative tile number


if nargin < 1
    imageFilenames = {'/Users/xruan/Images/20201201_p35_p4_LLS_Calibrations/488_totalPSF_Cropping_0p1_env_0_DOscan_z0p1.tif', ...
                      '/Users/xruan/Images/20201201_p35_p4_LLS_Calibrations/488_totalPSF_Cropping_0p1_env_0_DOscan_z0p1.tif', ...
                      '/Users/xruan/Images/20201201_p35_p4_LLS_Calibrations/488_totalPSF_Cropping_0p2_env_0_DOscan_z0p1.tif'};
end

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('imageFilenames'); %
ip.parse(imageFilenames, varargin{:});

% exclude partial files
imageFilenames = imageFilenames(~cellfun(@isempty, regexp(imageFilenames, '^(?!.*_part[0-9]*.tif).*$')));
% imageFilenames = unique(imageFilenames);

% first find the corresponding setting file for each image
[pathstrs, fsnames] = fileparts(imageFilenames);
sfns = arrayfun(@(x) [pathstrs{x}, filesep, fsnames{x}, '_Settings.txt'], 1 : numel(fsnames), 'unif', 0);
sfns_exist_mat = batch_file_exist(sfns);

if ~all(sfns_exist_mat)
    % mapping with the iter and xyz tile number    
    if all(~cellfun(@isempty, regexp(fsnames, '_Cam\w_ch', 'match')))
        expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>(\d+_)*\d+?)_Cam(?<Cam>\w+)_ch(?<ch>\d+)_CAM1_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t';
    elseif all(~cellfun(@isempty, regexp(fsnames, '_ch[0-9]_', 'match')))
        expression = '(?<prefix>\w*)Scan_Iter_(?<Iter>(\d+_)*\d+?)_ch(?<ch>\d+)_CAM1_stack(?<stack>\d+)_(?<laser>\d+)nm_(?<abstime>\d+)msec_(?<fpgatime>\d+)msecAbs_(?<x>-?\d+)x_(?<y>-?\d+)y_(?<z>-?\d+)z_(?<t>\d+)t';
        % specifyCam = false;
    end
    tmp = regexpi(fsnames, expression, 'names');
    
    sfns_1 = arrayfun(@(x) sprintf('%s/%sScan_Iter_%s_%sx_%sy_%sz_%st_Settings.txt', ...
        pathstrs{x}, tmp{x}.prefix, tmp{x}.Iter, tmp{x}.x, tmp{x}.y, tmp{x}.z, tmp{x}.t), ...
        1 : numel(fsnames), 'unif', 0);
    sfns_1_exist_mat = batch_file_exist(sfns_1);
end

sfn_cell = cell(numel(imageFilenames), 1);
for i = 1 : numel(imageFilenames)
    [pathstr, fsname] = fileparts(imageFilenames{i});
    
    % mapping with the filename
    if sfns_exist_mat(i)
        sfn_cell{i} = sfns{i};
        continue;
    end
        
    if sfns_1_exist_mat(i)
        sfn_cell{i} = sfns_1{i};
        continue;
    end
    
    % add another pattern
    sfn_2 = sprintf('%s/%sScan_Iter_%s_Settings.txt', pathstr, tmp{i}.prefix, tmp{i}.Iter);
    if exist(sfn_2, 'file')
        sfn_cell{i} = sfn_2;
        continue;
    end
end

% get the unique setting filenames
[uniq_sfn_cell, ia, ic] = unique(sfn_cell);

% then get setting info for each setting file
uniq_settingInfo = XR_parseSettingFiles(uniq_sfn_cell);
fieldnms = fieldnames(uniq_settingInfo);

settingInfo = struct();
for i = 1 : numel(imageFilenames)
    for j = 1 : numel(fieldnms)
        settingInfo(i).(fieldnms{j}) = uniq_settingInfo(ic(i)).(fieldnms{j});  
    end
end

end

