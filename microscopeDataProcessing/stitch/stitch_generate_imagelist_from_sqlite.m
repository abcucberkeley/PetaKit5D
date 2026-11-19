function [fnout] = stitch_generate_imagelist_from_sqlite(dataPath, varargin)
% generate image list file from sqlite file. 
% It reads all *.sqlite3 files in the given dataPath and generate image
% list from all the files. It saves as a csv file with name ImageList_from_sqlite.csv 
% in the dataPath. The format is consistent with old csv files. 
% 
%
% Author: Xiongtao Ruan (11/12/2021)


ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dataPath', @(x) ischar(x) || iscell(x));
ip.addParameter('channelPatterns', {'CamA', 'CamB'}, @iscell);
ip.addParameter('mapTilename', true, @islogical);
ip.parse(dataPath, varargin{:});

pr = ip.Results;
channelPatterns = pr.channelPatterns;
mapTilename = pr.mapTilename;

dir_info = dir([dataPath, '/*.sqlite3']);
fsns = {dir_info.name}';

% preallocate for 1000 rows, which should be enough for most siutations.
if mapTilename
    sz = [1000, 9];
    varTypes = ["string", "string", "string", "double", "double", "double", "double", "double", "double"];
    varNames = ["Filepath","Filename","mappedFilename","StageX_um_","StageY_um_","StageZ_um_","ObjectiveX_um_","ObjectiveY_um_","ObjectiveZ_um_"];

    tile_t_sz = [1000, 11];
    tile_t_varTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    tile_t_varNames = ["Filename","Camera","Channel","Laser_nm","File_Stack_Index","Elapsed_ms","Absolute_ms","Tile_index_X","Tile_index_Y","Tile_index_Z","Tile_index_T"];

    tile_info_tab = table('Size', tile_t_sz, 'VariableTypes', tile_t_varTypes, 'VariableNames', tile_t_varNames);
else
    sz = [1000, 8];
    varTypes = ["string", "string", "double", "double", "double", "double", "double", "double"];
    varNames = ["Filepath","Filename","StageX_um_","StageY_um_","StageZ_um_","ObjectiveX_um_","ObjectiveY_um_","ObjectiveZ_um_"];
end

t = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);

counter = 0;
for f = 1 : numel(fsns)
    dbfile = [dataPath, '/', fsns{f}];
    conn = sqlite(dbfile, 'readonly');
    connectionCleanup = onCleanup(@() close(conn));

    sqlquery = 'SELECT * FROM ImageListView_NoNull';

    try
        results = fetch(conn, sqlquery);
    catch ME
        if isempty(regexp(ME.message, 'no such column: JSONsettings_MD5\>', 'once'))
            rethrow(ME);
        end
        warning('stitch_generate_imagelist_from_sqlite:missingMetadataHashInView', ...
            'ImageListView_NoNull cannot resolve JSONsettings_MD5 in %s. Reading image and slice tables instead.', dbfile);
        sqlquery = [ ...
            'SELECT IFNULL(ImageList.Filename, ''N/A'') AS Filename, ' ...
            'IFNULL(ImageList.Camera, ''N/A'') AS Camera, ' ...
            'IFNULL(ImageList.Channel, ''N/A'') AS Channel, ' ...
            'IFNULL(ImageList.Laser_nm, ''N/A'') AS Laser_nm, ' ...
            'IFNULL(ImageList.File_Stack_Index, ''N/A'') AS File_Stack_Index, ' ...
            'IFNULL(ImageList.Elapsed_ms, ''N/A'') AS Elapsed_ms, ' ...
            'IFNULL(ImageList.Absolute_ms, ''N/A'') AS Absolute_ms, ' ...
            'IFNULL(ImageList.Tile_index_X, ''N/A'') AS Tile_index_X, ' ...
            'IFNULL(ImageList.Tile_index_Y, ''N/A'') AS Tile_index_Y, ' ...
            'IFNULL(ImageList.Tile_index_Z, ''N/A'') AS Tile_index_Z, ' ...
            'IFNULL(ImageList.Tile_index_T, ''N/A'') AS Tile_index_T, ' ...
            'AVG(IFNULL(SliceList.X_Target_Position_um, ''N/A'')) AS Mean_X_Target_Position_um, ' ...
            'AVG(IFNULL(SliceList.Y_Target_Position_um, ''N/A'')) AS Mean_Y_Target_Position_um, ' ...
            'AVG(IFNULL(SliceList.Z_Target_Position_um, ''N/A'')) AS Mean_Z_Target_Position_um ' ...
            'FROM SliceList ' ...
            'INNER JOIN ImageList ON ImageList.ImageList_id = SliceList.ImageList_id ' ...
            'INNER JOIN MetadataList ON MetadataList.JSONsettings_MD5 = ImageList.JSONsettings_MD5 ' ...
            'GROUP BY ImageList.ImageList_id ORDER BY ImageList.Filename'];
        results = fetch(conn, sqlquery);
    end
    clear connectionCleanup;

    if isempty(results)
        warning('%s is empty, skip it!', dbfile);
        continue;
    end

    ntile = size(results, 1);
    if mapTilename
        t(counter + 1 : counter + ntile, :) = results(:, {'Filename', 'Filename', 'Filename', ...
            'Mean_X_Target_Position_um', 'Mean_Y_Target_Position_um', 'Mean_Z_Target_Position_um', ...
            'Mean_X_Target_Position_um', 'Mean_Y_Target_Position_um', 'Mean_Z_Target_Position_um'});
        tile_info_tab(counter + 1 : counter + ntile, :) = results(:, {'Filename', ...
            'Camera', 'Channel', 'Laser_nm', 'File_Stack_Index', 'Elapsed_ms', ...
            'Absolute_ms', 'Tile_index_X', 'Tile_index_Y', 'Tile_index_Z', 'Tile_index_T'});
    else
        t(counter + 1 : counter + ntile, :) = results(:, {'Filename', 'Filename', ...
            'Mean_X_Target_Position_um', 'Mean_Y_Target_Position_um', 'Mean_Z_Target_Position_um', ...
            'Mean_X_Target_Position_um', 'Mean_Y_Target_Position_um', 'Mean_Z_Target_Position_um'});
    end
    counter = counter + ntile;
end
t(counter + 1 : end, :) = [];

fsn = t.Filename;
include_flag = false(numel(fsn), 1);
imageFullnames = cellfun(@(x) sprintf('%s/%s%s', dataPath, x), fsn, 'unif', 0);
for c = 1 : numel(channelPatterns)
    include_flag = include_flag | contains(imageFullnames, channelPatterns{c}) ...
        | contains(imageFullnames, channel_regexp_pattern(channelPatterns{c}));
end
if ~all(include_flag)
    t = t(include_flag, :);
end

if mapTilename
    tile_info_tab = tile_info_tab(include_flag, :);

    % update channel id for mapped channels
    cam_ch_laser_nm = [tile_info_tab.Camera, tile_info_tab.Channel, tile_info_tab.Laser_nm];
    uniq_cam_ch_laser_nm = unique(cam_ch_laser_nm, 'rows');

    uniq_cam = unique(uniq_cam_ch_laser_nm(:, 1));

    mapped_cam_ch_laser_nm_cell = cell(numel(uniq_cam), 1);
    for i = 1 : numel(uniq_cam)
        uniq_cam_ch_laser_nm_i = uniq_cam_ch_laser_nm(uniq_cam_ch_laser_nm(:, 1) == uniq_cam(i), :);
        mapped_ch_i = 0 : size(uniq_cam_ch_laser_nm_i, 1) - 1;
        mapped_cam_ch_laser_nm_cell{i} = [uniq_cam_ch_laser_nm_i, mapped_ch_i'];
    end
    mapped_uniq_cam_ch_laser_nm = cat(1, mapped_cam_ch_laser_nm_cell{:});

    mapped_cam_ch_laser_nm = cam_ch_laser_nm;
    for i = 1 : size(mapped_uniq_cam_ch_laser_nm, 1)
        mapped_ch_i = mapped_uniq_cam_ch_laser_nm(i, 4);
        mapped_cam_ch_laser_nm(all(cam_ch_laser_nm == mapped_uniq_cam_ch_laser_nm(i, 1 : 3), 2), 2) = mapped_ch_i;
    end
    tile_info_tab.Channel = mapped_cam_ch_laser_nm(:, 2);

    mappedFilenames = cell(size(tile_info_tab, 1), 1);
    cam_strs = {'A', 'B', 'C', 'D'};
    for i = 1 : numel(mappedFilenames)
        tinfo_i = table2array(tile_info_tab(i, 2 : end));
        mappedFilenames{i} = sprintf('Scan_Iter_%04d_Cam%s_ch%d_CAM1_stack%04d_%03dnm_%06dmsec_%08dmsecAbs_%03dx_%03dy_%03dz_%04dt.tif', ...
            tinfo_i(10), cam_strs{tinfo_i(1) + 1}, tinfo_i(2), tinfo_i(4), tinfo_i(3), ...
            tinfo_i(5), tinfo_i(6), tinfo_i(7), tinfo_i(8), tinfo_i(9), tinfo_i(10));
    end

    t(:, 3) = mappedFilenames(:);
end

% deduplicate tilenames if there are
[~, uniq_inds] = unique(t.Filename);
t = t(uniq_inds, :);

uuid = get_uuid();
tmpout = sprintf('%s/ImageList_from_sqlite_%s.csv', dataPath, uuid);
fnout = sprintf('%s/ImageList_from_sqlite.csv', dataPath);
writetable(t, tmpout, 'Delimiter', ',');
fileattrib(tmpout, '+w', 'g');
movefile(tmpout, fnout);

end

