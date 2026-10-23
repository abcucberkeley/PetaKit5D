function [file_exist_mat] = batch_file_exist(fileFullpaths, outFullpath, useParpool, useJava)
% check whether a batch of files or folders exist, with the same result as
% exist(..., 'file'). On Linux, each parent folder is opened (or listed for
% many files) before its files are checked: on NFS, exist() can miss a file
% another node made for up to the attribute cache time (30-60 s), while a
% file in a folder opened just before is seen at once. Elsewhere, the files
% are checked one by one (useParpool: with the background pool).
%
%   useJava : true|false (default: true). On Linux, open and list folders
%             with Java if the JVM runs; false uses dir(), as without it.
% 
% Author: Xiongtao Ruan (02/06/2022)
% 
% xruan (03/15/2022): if nd is large, directly check the input files
% xruan (09/06/2022): add support for thread based parallel computing 
% xruan (10/05/2023): java.io seems not working for paths with ~ 

if nargin < 2
    outFullpath = [];
end
if nargin < 3
    useParpool = false;
end
if nargin < 4
    useJava = true;
end
useJava = useJava && usejava('jvm');

if numel(fileFullpaths) == 0
    file_exist_mat = [];
    return;
end
    
if numel(fileFullpaths) == 1
    if isunix && ~ismac
        file_exist_mat = folder_file_exist(fileFullpaths, useJava);
    else
        file_exist_mat = exist(fileFullpaths{1}, 'file') || exist(fileFullpaths{1}, 'dir');
    end
    return;
end

if ~isempty(outFullpath) && exist(outFullpath, 'file')
    return;
end

if iscell(fileFullpaths) && iscell(fileFullpaths{1})
    fileFullpaths = cellfun(@(x) x{1}, fileFullpaths, 'unif', 0);
end

if isunix && ~ismac
    file_exist_mat = folder_file_exist(fileFullpaths, useJava);
else
    file_exist_mat = per_file_exist(fileFullpaths, useParpool);
end

if ~isempty(outFullpath) && ~exist(outFullpath, 'file')
    uuid = get_uuid();
    tmpFullpath = sprintf('%s_uuid_%s.mat', outFullpath(1 : end - 4), uuid);
    save('-v7.3', tmpFullpath, 'file_exist_mat', 'fileFullpaths');
    movefile(tmpFullpath, outFullpath);
end

end


function [file_exist_mat] = folder_file_exist(fileFullpaths, useJava)
% check each unique path once. Absolute paths are grouped by folder: a folder
% with few paths is opened and the paths checked by name; a larger one is
% listed, and the paths whose names the listing has (ignoring case, for
% case-insensitive file systems) are checked by name (a listing also shows
% broken symlinks), as are names a listing may spell differently (non-ASCII,
% short names with ~). Other paths (relative, ~, '/', '.', '..') go to
% exist(), which also finds relative names on the MATLAB path.

% more paths than this in a folder: list it instead of checking each path
maxNameCheckNum = 100;

[uniq_fullpaths, ~, ic] = unique(fileFullpaths(:));
nU = numel(uniq_fullpaths);
uniq_exist_mat = false(nU, 1);

fullpaths = regexprep(uniq_fullpaths, '/+$', '');
[pstrs, fsns, exts] = fileparts(fullpaths);
fns = strcat(fsns, exts);
by_folder = startsWith(fullpaths, '/') & ~ismember(fns, {'', '.', '..'});

for f = reshape(find(~by_folder), 1, [])
    uniq_exist_mat(f) = exist(uniq_fullpaths{f}, 'file') > 0;
end

odd_name = contains(fns, '~');
if any([fns{:}] > 127)
    odd_name = odd_name | cellfun(@(x) any(x > 127), fns);
end

folder_inds = reshape(find(by_folder), [], 1);
[uniq_pstrs, ~, pinds] = unique(pstrs(folder_inds));
nd = numel(uniq_pstrs);
if nd > 0
    inds_cell = accumarray(pinds(:), folder_inds, [nd, 1], @(x) {x});
end
for d = 1 : nd
    pstr_d = uniq_pstrs{d};
    inds_d = inds_cell{d};
    if ~isfolder(pstr_d)
        open_parent_folders(pstr_d, useJava);
        if ~isfolder(pstr_d)
            continue;
        end
    end

    if numel(inds_d) <= maxNameCheckNum
        open_folder(pstr_d, useJava);
        uniq_exist_mat(inds_d) = path_exist(fullpaths(inds_d));
        continue;
    end

    [fns_d, is_listed] = list_folder(pstr_d, useJava);
    if is_listed
        inds_d = inds_d(ismember(lower(fns(inds_d)), lower(fns_d)) | odd_name(inds_d));
    end
    uniq_exist_mat(inds_d) = path_exist(fullpaths(inds_d));
end

file_exist_mat = uniq_exist_mat(ic);

end


function [] = open_folder(folderPath, useJava)
% opening a folder revalidates its cached attributes on NFS, so a check by
% name right after it sees the files other nodes made in it

if useJava
    try
        ds = java.nio.file.Files.newDirectoryStream(java.io.File(folderPath).toPath());
        ds.close();
    catch
    end
else
    % without Java, read it with a pattern that matches nothing (a full
    % dir() costs ~20 us per file)
    dir_info = dir([folderPath, '/*.batch_file_exist_no_match']); %#ok<NASGU>
end

end


function [] = open_parent_folders(folderPath, useJava)
% a folder another node made is hidden to a check by name until its parent
% is revalidated: open the parents, from the first one that exists

pstr = fileparts(folderPath);
if isempty(pstr) || strcmp(pstr, folderPath)
    return;
end
if ~isfolder(pstr)
    open_parent_folders(pstr, useJava);
end
open_folder(pstr, useJava);

end


function [fns, is_listed] = list_folder(folderPath, useJava)
% file and folder names in a folder; is_listed is false if it can't be
% listed (e.g., no read permission)

fns = {};
if useJava
    % null (not listed) comes back as []
    fn_arr = java.io.File(folderPath).list();
    is_listed = isjava(fn_arr);
    if is_listed
        fns = cell(fn_arr);
    end
else
    dir_info = dir(folderPath);
    fns = {dir_info.name};
    is_listed = any(strcmp(fns, '.'));
    fns(strcmp(fns, '.') | strcmp(fns, '..')) = [];
end

end


function [is_exist] = path_exist(fullpaths)
% files and folders by name (follows symlinks, as exist())

fullpaths = fullpaths(:);
is_exist = isfile(fullpaths);
is_exist(~is_exist) = isfolder(fullpaths(~is_exist));

end


function [file_exist_mat] = per_file_exist(fileFullpaths, useParpool)
% check a batch of files by dir() for the files in the same directory to
% reduce the io for the check per file. 

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
        if ~exist(pstr_d, 'dir')
            continue;
        end
        
        dir_info = dir(pstr_d);
        fns_d = {dir_info.name};
        % without . and .. (a folder that can't be listed has neither)
        fns_d(strcmp(fns_d, '.') | strcmp(fns_d, '..')) = [];

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
    if usejava('jvm') && ~any(contains(fileFullpaths, '~'))
        for f = 1 : nF
            file_exist_mat(f) = java.io.File(fileFullpaths{f}).exists;
        end
    else
        for f = 1 : nF
            file_exist_mat(f) = exist(fileFullpaths{f}, 'file');
        end
    end
end

end
