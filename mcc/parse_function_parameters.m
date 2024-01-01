function [param_names, char_types] = parse_function_parameters(func_name)
% parse function names and whether it contains char type in the validators
% 
% Author: Xiongtao Ruan (03/16/2023)


if nargin < 1 
    func_name = 'XR_microscopeAutomaticProcessing.m';
end

func_fullpath = which(func_name);

fid = fopen(func_fullpath, 'r');
func_lines = cell(10000, 1);
ind = 1;
while ~feof(fid)
    cur_line = fgetl(fid);
    func_lines{ind} = cur_line;
    ind = ind + 1;
end
fclose(fid);
func_lines(ind : end) = [];

% get lines for parameters
tmp = regexp(func_lines, '^ip.addParameter');
pinds = ~cellfun(@(x) isempty(x), tmp);

param_lines = func_lines(pinds);

param_names = regexp(param_lines, 'ip.addParameter\(''(\w+)'',', 'tokens');
param_names = cellfun(@(x) x{1}{1}, param_names, 'UniformOutput', false);

char_types = regexp(param_lines, 'ischar');
char_types = ~cellfun(@(x) isempty(x), char_types);

end

