function [func_lines, arg_lines, arg_types] = generate_function_mcc_parser_file(funcName, dpntFuncNames, outPath)
% automatically parse the given function and generate function mcc parser
% The input function must use the input parser to define the input arguments
%
% Author: Xiongtao Ruan (04/04/2024)

if nargin < 1
    funcName = 'XR_psf_analysis_plot';
end
if nargin < 2
    dpntFuncNames = '';
end
if nargin < 3
    outPath = '~/Downloads/';
end

func_fullpath = which(funcName);

% read the function as a cell of strings
func_lines = readTextFile(func_fullpath);

% get function definition line and function name
tmp = regexp(func_lines, '^ *function');
finds = find(~cellfun(@(x) isempty(x), tmp), 1, 'first');

func_define_line = func_lines{finds};
if ~contains(func_define_line, ')')
    func_lines_1 = strjoin(func_lines, '\n');
    tmp = regexp(func_lines_1, '^ *(function [^)]*\))', 'tokens'); 
    func_define_line = tmp{1}{1};
end
if contains(func_define_line, '=')
    tmp = regexp(func_define_line, '^ *function *.* *= *([^ ]*)(', 'tokens');
else
    tmp = regexp(func_define_line, '^ *function *([^ ]*)(', 'tokens');
end
func_name = tmp{1}{1};

% get the required and optional arguments from the function line
func_define_line_1 = func_define_line;
if contains(func_define_line_1, '...')
    func_define_line_1 = strrep(func_define_line_1, '...', '');
    func_define_line_1 = strrep(func_define_line_1, newline, '');
end
tmp = regexp(func_define_line_1, '\((.*)\)', 'tokens');
func_arg_names = strip(strsplit(tmp{1}{1}, ','));
if strcmp(func_arg_names{end}, 'varargin')
    func_arg_names = func_arg_names(1 : end - 1);
end

parser_func_name = sprintf('%s_parser', func_name);

% remove comment lines 
tmp = regexp(func_lines, '^ *%');
ninds = cellfun(@(x) isempty(x), tmp);
func_lines = func_lines(ninds);

% only keep lines start with ip
s = find(contains(func_lines, 'inputParser'));
t = find(contains(func_lines, 'ip.parse'));
func_lines = func_lines(s : t + 5);

% get lines for required arguments
[required_lines, required_names, required_func_handles] = parse_parameter_structure(func_lines, 'required');

% get lines for optional arguments
[optional_lines, optional_names, optional_func_handles] = parse_parameter_structure(func_lines, 'optional');

% get lines for parameters
[param_lines, param_names, param_func_handles] = parse_parameter_structure(func_lines, 'parameter');

% aggragte together
arg_lines = cat(1, required_lines, optional_lines, param_lines);
arg_names = cat(1, required_names, optional_names, param_names);
arg_types = cat(1, repmat({'required'}, numel(required_names), 1), repmat({'optional'}, numel(optional_names), 1), ...
    repmat({'parameter'}, numel(param_names), 1));
arg_func_handles = cat(1, required_func_handles, optional_func_handles, param_func_handles);

% get the line for parameter parser
tmp = regexp(func_lines, 'ip.parse(');
ipinds = ~cellfun(@(x) isempty(x), tmp);
if isempty(ipinds)
    input_parser_func_line = '';
else
    input_parser_func_line = func_lines{ipinds};
    if ~contains(input_parser_func_line, ')')
        func_lines_1 = strjoin(func_lines, '\n');
        tmp = regexp(func_lines_1, ' *(ip.parse\([^)]*\);)', 'tokens');
        input_parser_func_line = tmp{1}{1};
    end
end

% check if required and optional argument missing and if they follow the same orders
ro_arg_names = arg_names(~contains(arg_types, 'parameter'));
if numel(func_arg_names) ~= numel(ro_arg_names)
    error('The number of required/optional arguments does not match that defined in the input parser!');
else
    for i = 1 : numel(func_arg_names) 
        if ~strcmpi(func_arg_names{i}, ro_arg_names{i})
            error('Argument %s does not match the one defined in the same order %s in the input parser!\n', func_arg_names{i}, ro_arg_names{i});
        end
    end

    % sort the order of the argument to see if they match
    func_arg_names_s = sort(func_arg_names);
    ro_arg_names_s = sort(ro_arg_names);
    for i = 1 : numel(func_arg_names_s) 
        if ~strcmpi(func_arg_names_s{i}, ro_arg_names_s{i})
            warning('Argument %s is not defined in the input parser!\n', func_arg_names_s{i}, ro_arg_names_s{i});
        end
    end
end

%% construct parser function

parser_cell = cell(numel(func_lines), 1);

% first line for the function definition
parser_func_define_line = strrep(func_define_line, func_name, parser_func_name);

parser_cell{1} = parser_func_define_line;

% add empty lines
parser_cell(2 : 3) = repmat({''}, 2, 1);
t = 3;

% add declaration of dependency
for i = 1 : numel(dpntFuncNames)
    parser_cell{t+i} = sprintf('%%#function %s', dpntFuncNames{i});
end
t = t + numel(dpntFuncNames);
if numel(dpntFuncNames) > 0
    parser_cell{t + 1} = '';
    t = t + 1;
end

parser_cell{t + 1} = 'ip = inputParser;';
parser_cell{t + 2} = 'ip.CaseSensitive = false;';
t = t + 2;

% construct parser function parser lines
parser_arg_lines = construct_parser_function_parser_lines(arg_lines, arg_func_handles);
parser_cell(t + 1 : t + numel(parser_arg_lines)) = parser_arg_lines;
t = t + numel(parser_arg_lines);

parser_cell{t + 1} = '';
parser_cell{t + 2} = input_parser_func_line;
parser_cell{t + 3} = '';
parser_cell{t + 4} = 'pr = ip.Results;';
t = t + 4;

% define the parameters
for i = 1 : numel(param_lines)
    param_name_i = param_names{i};
    parser_cell{t+i} = sprintf('%s = pr.%s;', param_name_i, param_name_i);
end
t = t + numel(param_lines);
parser_cell{t + 1} = '';
t = t + 1;

% add lines to check the actual data type and do the conversion for the parameters
% do that for required inputs
arg_converter_lines = check_and_update_actual_argment_type(arg_names, arg_func_handles);
parser_cell(t + 1 : t + numel(arg_converter_lines)) = arg_converter_lines;
t = t + numel(arg_converter_lines);

parser_cell{t+1} = '';

%  construct function all
max_char_num = 75;
function_call_line = construct_function_call_line(func_name, arg_names, arg_types, max_char_num);

parser_cell{t+2} = function_call_line;
parser_cell{t+3} = '';
parser_cell{t+4} = 'end';
parser_cell{t+5} = '';
parser_cell(t+6 : end) = [];

% write out to a file
% if the parser function exist, check if the content is the same
fnout = sprintf('%s/%s.m',  outPath, parser_func_name);
if exist(fnout, 'file')
    exist_parser_cell = readTextFile(fnout);
    tmpout = sprintf('%s/%s_tmp.m',  outPath, parser_func_name);
    writeTextFile(parser_cell, tmpout);
    tmp_parser_cell = readTextFile(tmpout);

    if strcmp(strjoin(tmp_parser_cell, '\n'), strjoin(exist_parser_cell, '\n'))
        fprintf('    The input parser remains the same, skip the update of the parser function.\n');
        delete(tmpout);
    else
        fprintf('    The function has changed, update the parser function.\n');
        movefile(tmpout, fnout);
    end
else
    writeTextFile(parser_cell, fnout);
end

end


function [arg_lines, arg_names, arg_func_handles] = parse_parameter_structure(func_lines, type)
% function to parse input parser, get argument name and arg function handle
% based on the argument type

switch type
    case 'required'
        pstr = 'ip.addRequired';
    case 'optional'
        pstr = 'ip.addOptional';
    case 'parameter'
        pstr = 'ip.addParameter';
end

tmp = regexp(func_lines, pstr);
rinds = ~cellfun(@(x) isempty(x), tmp);

arg_lines = func_lines(rinds);

tmp = regexp(arg_lines, [pstr, '\(''(\w+)'',?'], 'tokens');
arg_names = cellfun(@(x) x{1}{1}, tmp, 'UniformOutput', false);

if strcmp(type, 'parameter')
    tmp = regexp(arg_lines, [pstr, '\(''\w+'', *.+, *(@[^;]+)\);'], 'tokens');
    arg_func_handles = cellfun(@(x) x{1}{1}, tmp, 'UniformOutput', false);    
else
    tmp = regexp(arg_lines, [pstr, '\(''\w+'', *(@[^;]+)\);'], 'tokens');
    tmp_1 = regexp(arg_lines, [pstr, '\(''\w+'', *.+, *(@[^;]+)\);'], 'tokens');

    arg_func_handles = cell(numel(arg_names), 1);
    for i = 1 : numel(tmp)
        if ~isempty(tmp{i})
            arg_func_handles{i} = tmp{i}{1}{1};
        end
        if ~isempty(tmp_1{i})
            arg_func_handles{i} = tmp_1{i}{1}{1};
        end
    end
end

end


function [parser_lines] = construct_parser_function_parser_lines(arg_lines, arg_func_handles)
% construct the parser lines for parser function by adding ischar check if
% the original function handle does not contain char type check

n = numel(arg_lines);
parser_lines = cell(n, 1);

for i = 1 : n
    arg_line_i = arg_lines{i};
    func_handle = arg_func_handles{i};
    if ~isempty(func_handle) && ~contains(func_handle, 'ischar')
        if contains(func_handle, '@(x)')
            func_handle_updated = sprintf('%s || ischar(x)', func_handle);
        else
            func_handle_updated = sprintf('@(x) %s(x) || ischar(x)', strip(strrep(func_handle, '@', ''), ' '));
        end
        arg_line_i = strrep(arg_line_i, func_handle, func_handle_updated);
    end
    parser_lines{i} = arg_line_i;
end

end

function [converter_lines] = check_and_update_actual_argment_type(arg_names, arg_func_handles)
% check if the actual argument data type, if it is not char, convert back
% to the actual type. 


n = numel(arg_names);
converter_lines = cell(n, 1);
t = 0;

for i = 1 : numel(arg_names)
    arg_name_i = arg_names{i};
    func_handle = arg_func_handles{i};
    if isempty(func_handle) || (contains(func_handle, 'ischar') && ~contains(func_handle, 'iscell') ...
            && ~contains(func_handle, 'function_handle') && ~contains(func_handle, 'isstruct') && ~contains(func_handle, 'isempty'))
        continue;
    end

    if contains(func_handle, 'isnumeric') || contains(func_handle, 'islogical') || contains(func_handle, 'isvector') || contains(func_handle, 'isscalar')
        converter_lines{t+1} = sprintf('if ischar(%s)\n    %s = str2num(%s);\nend', arg_name_i, arg_name_i, arg_name_i); 
        t = t + 1;
        continue
    end
    
    if contains(func_handle, 'iscell')
        converter_lines{t+1} = sprintf('if ischar(%s) && ~isempty(%s) && strcmp(%s(1), ''{'')\n    %s = eval(%s);\nend', ...
            arg_name_i, arg_name_i, arg_name_i, arg_name_i, arg_name_i); 
        t = t + 1;
        continue
    end

    if contains(func_handle, 'isstruct')
        converter_lines{t+1} = sprintf('if ischar(%s) && ~isempty(%s) && strcmp(%s(1), ''['')\n    %s = eval(%s);\nend', ...
            arg_name_i, arg_name_i, arg_name_i, arg_name_i, arg_name_i); 
        t = t + 1;
        continue
    end
    
    if contains(func_handle, 'function_handle')
        converter_lines{t+1} = sprintf('if ischar(%s) && ~isempty(%s) && strcmp(%s(1), ''@'')\n    %s = eval(%s);\nend', ...
            arg_name_i, arg_name_i, arg_name_i, arg_name_i, arg_name_i);
        t = t + 1;
        continue
    end
    
    if contains(func_handle, 'isempty')
        converter_lines{t+1} = sprintf('if ischar(%s) && ~isempty(%s) && (strcmp(%s(1), ''['') || strcmp(%s(1), ''{''))\n    %s = eval(%s);\nend', ...
            arg_name_i, arg_name_i, arg_name_i, arg_name_i, arg_name_i, arg_name_i); 
        t = t + 1;
        continue
    end    
end

converter_lines(t + 1 : end) = [];

end


function [function_call_line] = construct_function_call_line(func_name, arg_names, arg_types, max_char_num)
% construct the function call line based on the argument type
% 
% arg_types: required, optional, parameter


if nargin < 3
    max_char_num = 75;
end

function_call_line = sprintf('%s(', func_name);

% required inputs
counter = numel(function_call_line);
for i = 1 : numel(arg_names)
    arg_name_i = arg_names{i};
    arg_type_i = arg_types{i};
    if counter >= max_char_num
        function_call_line = sprintf('%s...\n    ', function_call_line);
        counter = 4;
    end
    if strcmp(arg_type_i, 'parameter')
        cur_str = sprintf('%s=%s, ', arg_name_i, arg_name_i);
        if i == numel(arg_names)
            cur_str = sprintf('%s=%s);', arg_name_i, arg_name_i);        
        end
    else
        cur_str = sprintf('%s, ', arg_name_i);
        if i == numel(arg_names)
            cur_str = sprintf('%s);', arg_name_i);        
        end        
    end
    counter = counter + numel(cur_str);    
    function_call_line = sprintf('%s%s', function_call_line, cur_str);
end

end

