function [file_exist_mat] = batch_file_exist(fileFullpaths, outFullpath, useParpool)
% check a batch of files by dir() for the files in the same directory to
% reduce the io for the check per file. 
% 
% Author: Xiongtao Ruan (02/06/2022)
% 
% xruan (03/15/2022): if nd is large, directly check the input files
% xruan (09/06/2022): add support for thread based parallel computing 

if nargin < 2
    outFullpath = [];
end
if nargin < 3
    useParpool = false;
end

if numel(fileFullpaths) == 0
    file_exist_mat = [];
    return;
end
    
if numel(fileFullpaths) == 1
    file_exist_mat = exist(fileFullpaths{1}, 'file') || exist(fileFullpaths{1}, 'dir');
    return;
end

if ~isempty(outFullpath) && exist(outFullpath, 'file')
    return;
end

if iscell(fileFullpaths) && iscell(fileFullpaths{1})
    fileFullpaths = cellfun(@(x) x{1}, fileFullpaths, 'unif', 0);
end

[pstrs, fsns, exts] = fileparts(fileFullpaths);

uniq_pstrs = unique(pstrs);
nd = numel(uniq_pstrs);

nF = numel(fsns);
fns = arrayfun(@(x) [fsns{x}, exts{x}], 1 : nF, 'unif', 0);

file_exist_mat = false(nF, 1);

% if there are less than 3 folders, and nF / nD > 100
check_folder = false;
if (nd < 3) && (nF / nd > 100)
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
    check_folder = true;
end

% go through each file
if useParpool && nF > 100 && ~check_folder
    % if ~isempty(gcp('nocreate'))
    %     delete(gcp('nocreate'));
    % end
    p = backgroundPool;
    nworker = p.NumWorkers;

    fs = parallel.FevalFuture;
    for f = 1 : nF
        fs(f) = parfeval(p, @exist, 1, fileFullpaths{f}, 'file');            
    end
    
    wait(fs, 'finished', nF / nworker * 0.1);
    file_exist_mat = fetchOutputs(fs) > 0;
else
    if usejava('jvm')
        for f = 1 : nF
            file_exist_mat(f) = java.io.File(fileFullpaths{f}).exists;
        end
    else
        for f = 1 : nF
            file_exist_mat(f) = exist(fileFullpaths{f}, 'file');
        end
    end
end

if ~isempty(outFullpath) && ~exist(outFullpath, 'file')
    uuid = get_uuid();
    tmpFullpath = sprintf('%s_uuid_%s.mat', outFullpath(1 : end - 4), uuid);
    save('-v7.3', tmpFullpath, 'file_exist_mat', 'fileFullpaths');
    movefile(tmpFullpath, outFullpath);
end

end
