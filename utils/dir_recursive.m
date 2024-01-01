function [filenames] = dir_recursive(given_dir, fmode, useParpool)

% xruan 08/13/2015
% get all filenames under the current dir and its subdiretories 
% use DFS to recursively get all filenames (no directory)
% 
% xruan (03/16/2022): replace dir with jave methods for acceleration,
% especially for direcotries with millions or more files. 
% the idea is from https://stackoverflow.com/questions/6385531/very-slow-dir
% 
% xruan (03/18/2022): add support to only get folders
% xruan (04/28/2022): fix bug for files/folder with leading non-letter/digit symbols
% xruan (09/06/2022): add support for thread based parallel computing 

if nargin < 2
    fmode = 'file';
end  

if nargin < 3
    useParpool = false;
end

check_dir = false;
switch fmode
    case 'file'
        check_dir = false;
    case 'folder'
        check_dir = true;
end

if strcmp(given_dir(end), filesep)
    given_dir(end) = [];
end

filenames = {};

stack = {};

if isfolder(given_dir)
    stack{end + 1} = given_dir;
end

if useParpool 
    if ~isempty(gcp('nocreate'))
        delete(gcp('nocreate'));
    end
    p = backgroundPool;
    nworker = p.NumWorkers;
end

loop_counter = 0;
while ~isempty(stack)
    if ~useParpool || numel(stack) == 1
        % pop one element
        current_dir = stack{end};
        stack(end) = [];
        
        % tic
        fileinfo = dir(current_dir);
        % fileinfo(1 : 2) = [];    
        if numel(fileinfo) == 2
            continue;
        end
        subdir_files = {fileinfo.name}';
        % remove . and ..
        idx = ~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), subdir_files);
        subdir_files = subdir_files(idx);
        sub_isdir = [fileinfo(idx).isdir];
    
        % toc
    %     tic
    %     jFile = java.io.File(current_dir); 
    %     jPaths = jFile.listFiles;
    %     jNames = jFile.list;
    %     isFolder = arrayfun(@isDirectory,jPaths);
    %     
    %     subdir_files = cellstr(char(jNames));
    %     sub_isdir = isFolder;
    %     toc
    
        subdirs = cellfun(@(x) [current_dir, filesep, x], subdir_files(sub_isdir), 'UniformOutput', false);
        if ~isempty(subdirs)
            stack(end + 1 : end + length(subdirs)) = subdirs;
        end
        
        if check_dir
            subfiles = cellfun(@(x) [current_dir, filesep, x], subdir_files(sub_isdir), 'UniformOutput', false);
        else
            subfiles = cellfun(@(x) [current_dir, filesep, x], subdir_files(~sub_isdir), 'UniformOutput', false);        
        end
        if ~isempty(subfiles)
            filenames(end + 1 : end + length(subfiles)) = subfiles;
        end
    else
        fs = parallel.FevalFuture;
        for s = 1 : numel(stack)
            current_dir = stack{end - s + 1};
            fs(s) = parfeval(p, @dir, 1, current_dir);            
        end
        
        wait(fs, 'finished', numel(stack) / nworker * 0.1);
        stack = {};

        % retrieve results
        fileinfo = fetchOutputs(fs);
        if numel(fileinfo) == 2
            continue;
        end
        dir_names = {fileinfo.folder}';
        subdir_files = {fileinfo.name}';
        % remove . and ..
        idx = ~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), subdir_files);
        dir_names = dir_names(idx);
        subdir_files = subdir_files(idx);
        sub_isdir = [fileinfo(idx).isdir];
        sub_fns = cellfun(@(x, y) sprintf('%s/%s', x, y), dir_names, subdir_files, 'UniformOutput', false);

        subdirs = sub_fns(sub_isdir);
        if ~isempty(subdirs)
            % stack(end + 1 : end + length(subdirs)) = subdirs;
            % remove any subdir that is already in the stack (some symbolic
            % directories may point to the upper level of directories).
            stack = subdirs;
            if loop_counter > 20
                base_depth = numel(regexp(given_dir, '/'));
                dir_depths = cellfun(@numel, regexp(subdirs, '/'));
                dir_included = dir_depths < base_depth + loop_counter;
                stack = subdirs(~dir_included);
            end
        end
        
        if check_dir
            subfiles = subdirs;
        else
            subfiles = sub_fns(~sub_isdir);
        end
        if ~isempty(subfiles)
            filenames(end + 1 : end + length(subfiles)) = subfiles;
        end
    end

    loop_counter = loop_counter + 1;
end

filenames = unique(filenames);

end

