function [] = mcc_autogen(funcListFn, outPath, compileMCC, nojvm)
% automatically generate the mcc parsers and mccMaster for a given function
% list file. The funciton list file should be a text file that list the
% candidate functions to generate mcc parser. Each line contains a
% function. If a function has dependency functions, the dependency
% functions should be list after the function separated by space. 
% 
% Author: Xiongtao Ruan (04/04/2024)


if nargin < 1
    [fpath, fname] = fileparts(which('mcc_autogen.m'));    
    funcListFn = [fpath, '/mcc_function_list.txt'];
end
if nargin < 2 || isempty(outPath)
    pstr = fileparts(funcListFn);
    outPath = [pstr, '/', 'parsers/'];
    mkdir(outPath);
end
if nargin < 3
    compileMCC = false;
end
if nargin < 4
    nojvm = true;
end

if ~exist(funcListFn, 'file')
    error('The function list file %s does not exist, please check the path!', funcListFn);
end

% parse the function list file
[func_names, dpnt_lists] = parse_mcc_function_list_file(funcListFn);

% generate parser functions for the func names 
pos_arg_num_mat = zeros(numel(func_names), 1);
param_num_mat = zeros(numel(func_names), 1);
for i = 1 : numel(func_names)
    func_name_i = func_names{i};
    dpnt_list_i = dpnt_lists{i};
    fprintf('Generate parser for %s.m...\n', func_name_i);
    [~, ~, arg_types] = generate_function_mcc_parser_file(func_name_i, dpnt_list_i, outPath);
    param_inds = strcmp(arg_types, 'parameter');
    pos_arg_num_mat(i) = sum(~param_inds);
    param_num_mat(i) = sum(param_inds);
end

% generate mccMaster
fprintf('Generate mccMaster.m...\n')
generate_mcc_master_file(func_names, pos_arg_num_mat, param_num_mat, outPath)
fprintf('Done!\n')

if compileMCC
    compile_mccMaster(nojvm);
end

end
