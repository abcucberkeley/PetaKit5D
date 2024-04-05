function [func_names, dpnt_lists] = parse_mcc_function_list_file(filename)
% parse the mcc function list file to extract the function list to generate
% mcc parser and their dependent functions


if ~exist(filename, 'file')
    error('The function list file %s does not exist, please check the path!', filename);
end

func_lists = readTextFile(filename);

% remove the comment lines with leading % or #
tmp = regexp(func_lists, '^ *[%#]');
finds = cellfun(@(x) isempty(x), tmp);
func_lists = func_lists(finds);

% get actual function list
tmp = regexp(func_lists, ' *\w.*');
inds = ~cellfun(@(x) isempty(x), tmp);
func_lists = func_lists(inds);

% remove potential .m
func_lists = strrep(func_lists, '.m', '');

% strip spaces
func_lists = strip(func_lists);

% remove duplicate functions
func_lists = unique(func_lists, 'stable');

func_names = cell(numel(func_lists), 1);
dpnt_lists = cell(numel(func_lists), 1);

% separate function name and dependency functions
for i = 1 : numel(func_lists)
    func_list_i = func_lists{i};
    tmp = strsplit(strip(func_list_i), ' ');
    empty_inds = cellfun(@(x) contains(x, ' '), tmp);
    tmp(empty_inds) = [];

    func_names{i} = tmp{1};
    if numel(tmp) > 1
        dpnt_lists{i} = tmp(2 : end);
    end
end

end
