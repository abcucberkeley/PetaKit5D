function [file_exist_mat] = batch_file_exist(fileFullpaths)
% check a batch of files by dir() for the files in the same directory to
% reduce the io for the check per file. 
% 
% Author: Xiongtao Ruan (02/06/2022)


if numel(fileFullpaths) == 1
    file_exist_mat = exist(fileFullpaths{1}, 'file') || exist(fileFullpaths{1}, 'dir');
    return;
end

[pstrs, fsns, exts] = fileparts(fileFullpaths);

uniq_pstrs = unique(pstrs);
nd = numel(uniq_pstrs);

nF = numel(fsns);
fns = arrayfun(@(x) [fsns{x}, exts{x}], 1 : nF, 'unif', 0);

file_exist_mat = false(nF, 1);
for d = 1 : nd
    pstr_d = uniq_pstrs{d};
    
    dir_info = dir(pstr_d);
    dir_info(1 : 2) = [];
    
    fns_d = {dir_info.name};
    
    inds_d = find(strcmp(pstrs, pstr_d));
    
    for f = 1 : numel(inds_d)
        file_exist_mat(inds_d(f)) = any(strcmp(fns_d, fns{inds_d(f)}));
    end    
end

end
