function [fnout] = stitch_generate_imagelist_from_sqlite(dataPath)
% generate image list file from sqlite file. 
% It reads all *.sqlite3 files in the given dataPath, and generate image
% list from all the files. It saves as a csv file with name ImageList_from_sqlite.csv 
% in the dataPath. The format is consistent with old csv files. 
% 
%
% Author: Xiongtao Ruan (11/12/2021)


dir_info = dir([dataPath, '*.sqlite3']);
fsns = {dir_info.name}';

% preallocate for 1000 rows, which should be enough for most siutations.
sz = [1000, 8];
varTypes = ["string", "string", "double", "double", "double", "double", "double", "double"];
varNames = ["Filepath","Filename","StageX_um_","StageY_um_","StageZ_um_","ObjectiveX_um_","ObjectiveY_um_","ObjectiveZ_um_"];

t = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);

counter = 0;
for f = 1 : numel(fsns)
    dbfile = [dataPath, fsns{f}];
    conn = sqlite(dbfile, 'readonly');

    sqlquery = 'SELECT * FROM ImageListView_NoNull';

    results = fetch(conn, sqlquery);
    close(conn);
    
    if isempty(results)
        warning('%s is empty, skip it!', dbfile);
        continue;
    end
    
    ntile = size(results, 1);
    
    t(counter + 1 : counter + ntile, :) = results(:, [2, 2, 14 : 19]);
    counter = counter + ntile;
end
t(counter + 1 : end, :) = [];

uuid = get_uuid();
tmpout = sprintf('%s/ImageList_from_sqlite_%s.csv', dataPath, uuid);
fnout = sprintf('%s/ImageList_from_sqlite.csv', dataPath);
writetable(t, tmpout, 'Delimiter', ',');
movefile(tmpout, fnout);

end

