function [settingInfo] = XR_parseSettingFiles(settingFilenames, varargin)
% parse a given list of settting files to get useful setting info for the
% image processing
% 
% Author: Xiongtao Ruan (12/05/2020)


if nargin < 1
    settingFilenames = {'/Users/xruan/Images/20201201_p35_p4_LLS_Calibrations/488_totalPSF_Cropping_0p1_env_0_DOscan_z0p1_Settings.txt', ...
                        '/Users/xruan/Images/20201201_p35_p4_LLS_Calibrations/488_totalPSF_Cropping_0p1_env_0_DOscan_z0p1_Settings_test.txt'};
end

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('settingFilenames'); %
ip.parse(settingFilenames, varargin{:});


settingInfo = struct();
for i = 1 : numel(settingFilenames)
    
    fn = settingFilenames{i};
    
    if ~exist(fn, 'file')
        warning('Setting file %s does not exist, skip it.', fn)
        continue;
    end
    
    sInfo_i = parseSingleSettingFile(fn);
    fieldnms = fieldnames(sInfo_i);
    for j = 1 : numel(fieldnms)
        settingInfo(i).(fieldnms{j}) = sInfo_i.(fieldnms{j});  
    end
end

end


function [sInfo] = parseSingleSettingFile(fn)
% parse a single setting file


fid = fopen(fn, 'r');

line_cell = cell(10000, 1);
tline = fgetl(fid);
counter = 1;
line_cell{counter} = tline;
while ischar(tline)
    tline = fgetl(fid);
    counter = counter + 1;
    line_cell{counter} = tline;
end
fclose(fid);

% tline in the last loop is -1, not a string
line_cell(counter : end) = [];
allline = strjoin(line_cell, '\n');

% get the sign of the X Stage interval
pattern = '\n(XZ? +Stage Offset, Interval \(um\), # of Pixels for Excitation \(0\) :\t\d\t-?\d*?.?\d+\t\d+)\n';
str = regexpi(allline, pattern, 'tokens');
str = str{1}{1};

expression = '(?<stageOff>\d+)\t(?<interval>-?\d*?.?\d+)\t(?<numExiciation>\d+)';

tmp = regexpi(str, expression, 'names');

sInfo.stageOff = str2double(tmp.stageOff);
sInfo.StageInterval = str2double(tmp.interval);
sInfo.numExiciation = str2double(tmp.numExiciation);

end

