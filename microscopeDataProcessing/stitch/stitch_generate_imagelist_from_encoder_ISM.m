function stitch_generate_imagelist_from_encoder_ISM(dataPath, positionType)
% generate image list file from encoder file for the raw scan data. 
% It reads all encoder files in the given dataPath, and generate image
% list from all the files. It saves as a csv file with name ImageList_from_sqlite.csv 
% in the dataPath. The format is consistent with old csv files. 
% 
%
% Author: Xiongtao Ruan (11/19/2021)
% xruan (08/08/2022): add support for user defined position type

if nargin < 2
    positionType = 'MCS2'
end

% get the slice list files and load the entries of the first image indices
dir_info = dir([dataPath, 'SliceList_*.csv']);
fsns = {dir_info.name}';
t_cell = cell(numel(fsns));
for i = 1 : numel(fsns)
    t = readtable([dataPath, fsns{i}]);
    
    t_cell{i} = t(t.DataImageIndex == 0, :);
end

tin = cat(1, t_cell{:});
clear t t_cell;

% generate image list file from encoder. 
sz = [size(tin, 1), 8];
varTypes = ["string", "string", "double", "double", "double", "double", "double", "double"];
varNames = ["Filepath","Filename","StageX_um_","StageY_um_","StageZ_um_","ObjectiveX_um_","ObjectiveY_um_","ObjectiveZ_um_"];

t = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
t.Filename = strrep(tin.Filename, '_0000ISMslice', '');
t.Filepath = t.Filename;
switch lower(positionType)
    case 'mcs2'
        t.StageX_um_ = tin.XMCS2Position_um_;
        t.StageY_um_ = tin.YMCS2Position_um_;
        t.StageZ_um_ = tin.ZMCS2Position_um_;
        t.ObjectiveX_um_ = tin.XMCS2Position_um_;
        t.ObjectiveY_um_ = tin.YMCS2Position_um_;
        t.ObjectiveZ_um_ = tin.ZMCS2Position_um_;
    case 'target'
        t.StageX_um_ = tin.XTargetPosition_um_;
        t.StageY_um_ = tin.YTargetPosition_um_;
        t.StageZ_um_ = tin.ZTargetPosition_um_;
        t.ObjectiveX_um_ = tin.XTargetPosition_um_;
        t.ObjectiveY_um_ = tin.YTargetPosition_um_;
        t.ObjectiveZ_um_ = tin.ZTargetPosition_um_;
end
fnout = sprintf('%s/ImageList_ISM_from_encoder.csv', dataPath);
writetable(t, fnout, 'Delimiter', ',');


end