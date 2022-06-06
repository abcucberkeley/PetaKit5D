function [fnout] = stitch_generate_imagelist_from_encoder(dataPath, dz, ChannelPatterns)
% generate image list file from encoder files. 
% It reads all encoder files in the given dataPath, and generate image
% list from all the files. It saves as a csv file with name ImageList_from_sqlite.csv 
% in the dataPath. The format is consistent with old csv files. 
% 
%
% Author: Xiongtao Ruan (11/19/2021)
% 
% xruan (06/06/2022): add support for user define channel patterns

if nargin < 3
    ChannelPatterns = {'CamA_ch0', 'CamB_ch0'};
end

% generate coordinate info
dz_actual = XR_estimate_actual_step_size_from_encoder(dataPath, 'dz', dz, ...
    'ChannelPatterns', ChannelPatterns, 'onlyFirstRow', true, 'skipMissingFile', true);
encoder_fn = [dataPath, '/encoder_info/coordinates_info.mat'];

% generate image list file from encoder. 
a = load(encoder_fn);
coords_struct = a.coords_struct;

sz = [1000, 8];
varTypes = ["string", "string", "double", "double", "double", "double", "double", "double"];
varNames = ["Filepath","Filename","StageX_um_","StageY_um_","StageZ_um_","ObjectiveX_um_","ObjectiveY_um_","ObjectiveZ_um_"];

t = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);

counter = 0;
for i = 1 : numel(coords_struct)
    fpth = sprintf('%s/%s', dataPath, coords_struct(i).Filename);
    fpth_a = strrep(fpth, 'CamB', 'CamA');
    fsn_a = strrep(coords_struct(i).Filename, 'CamB', 'CamA');

    % t(counter + 1, :) = {fpth_a, fsn_a, coords_struct(i).xmcs2(1), coords_struct(i).ymcs2(1), coords_struct(i).zmcs2(1), coords_struct(i).xmcs2(1), coords_struct(i).ymcs2(1), coords_struct(i).zmcs2(1)};
    t(counter + 1, :) = {fpth, coords_struct(i).Filename, coords_struct(i).xmcs2(1), coords_struct(i).ymcs2(1), coords_struct(i).zmcs2(1), coords_struct(i).xmcs2(1), coords_struct(i).ymcs2(1), coords_struct(i).zmcs2(1)};
    counter = counter + 1;
end
t(counter + 1 : end, :) = [];
t = unique(t, 'rows');

uuid = get_uuid();
tmpout = sprintf('%s/ImageList_from_encoder_%s.csv', dataPath, uuid);
fnout = sprintf('%s/ImageList_from_encoder.csv', dataPath);
writetable(t, tmpout, 'Delimiter', ',');
movefile(tmpout, fnout);

end

