function [] = generate_mcc_master_file(funcNames, posArgNums, paramArgNums, outPath)
% automatically generate mcc master function
%
% Author: Xiongtao Ruan (04/04/2024)


parser_cell = cell(numel(funcNames) * 2 + 20, 1);
parser_cell{1} = sprintf('function mccMaster(functionName, varargin)');
parser_cell{2} = sprintf('');
parser_cell{3} = sprintf('%%#function setup');
parser_cell{4} = sprintf('');
parser_cell{5} = sprintf('get_hostname();');
parser_cell{6} = sprintf('fprintf(''Function name: %%s\\nVariables: %%s \\n'', functionName, strjoin(varargin, '' ''));');
parser_cell{7} = sprintf('t0 = tic;');
parser_cell{8} = sprintf('');
parser_cell{9} = sprintf('switch functionName');

t = 9;
for i = 1 : numel(funcNames)
    func_name_i = funcNames{i};
    pos_arg_num = posArgNums(i);
    param_arg_num = paramArgNums(i);
    case_line = sprintf('    case ''%s''', func_name_i);
    
    func_call_line = sprintf('        %s_parser(', func_name_i);
    for j = 1 : pos_arg_num
        if j ~= pos_arg_num
            func_call_line = sprintf('%svarargin{%d}, ', func_call_line, j);
        else
            func_call_line = sprintf('%svarargin{%d}', func_call_line, j);
        end
    end

    if pos_arg_num > 0
        if param_arg_num > 0
            func_call_line = sprintf('%s, varargin{%d:end});', func_call_line, pos_arg_num + 1);
        else
            func_call_line = sprintf('%s);', func_call_line);
        end
    else
        func_call_line = sprintf('%svarargin{:});', func_call_line);
    end

    parser_cell{t+i*2-1} = case_line;
    parser_cell{t+i*2} = func_call_line;   
end

t = t + 2 * numel(funcNames);
parser_cell{t + 1} = sprintf('    otherwise\n        error(''The parser for function %%s does not exist!'', functionName);');
parser_cell{t + 2} = sprintf('end');
parser_cell{t + 3} = sprintf('');
parser_cell{t + 4} = sprintf('toc(t0);');
parser_cell{t + 5} = sprintf('');
parser_cell{t + 6} = sprintf('exit;');
parser_cell{t + 7} = sprintf('');
parser_cell{t + 8} = sprintf('end');
parser_cell{t + 9} = sprintf('');
parser_cell(t + 10 : end) = [];

% write out to a file
% if the parser function exist, check if the content is the same
fnout = sprintf('%s/mccMaster.m',  outPath);
if exist(fnout, 'file')
    exist_parser_cell = readTextFile(fnout);
    tmpout = sprintf('%s/mccMaster_tmp.m',  outPath);
    writeTextFile(parser_cell, tmpout);
    tmp_parser_cell = readTextFile(tmpout);

    if strcmp(strjoin(tmp_parser_cell, '\n'), strjoin(exist_parser_cell, '\n'))
        fprintf('    The function mccMaster.m remains the same, skip the update of the function.\n');
        delete(tmpout);
    else
        fprintf('    The function mccMaster.m has changed, update the parser function.\n');
        movefile(tmpout, fnout);
    end
else
    writeTextFile(parser_cell, fnout);
end

end
