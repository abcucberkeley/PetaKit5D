function [] = generate_mcc_master_file(funcNames, posArgNums, paramArgNums, outPath, outFuncName)
% automatically generate mcc master function
%
% Author: Xiongtao Ruan (04/04/2024)


if nargin < 5
    outFuncName = 'mccMaster';
end

parser_cell = cell(numel(funcNames) * 2 + 30, 1);
parser_cell{1} = sprintf('function %s(functionName, varargin)', outFuncName);
parser_cell{2} = sprintf('');
parser_cell{3} = sprintf('%%#function setup');
parser_cell{4} = sprintf('');
parser_cell{5} = sprintf('get_hostname();');
parser_cell{6} = sprintf('fprintf(''Function name: %%s\\nVariables: %%s \\n'', functionName, strjoin(varargin, '' ''));');
parser_cell{7} = sprintf('t0 = tic;');
% the task's resources when the job framework asks for them (in the
% environment of the task line, so builds without it run the same lines)
parser_cell{8} = sprintf('%% the task''s resources, when the job framework asks for them (task_metrics)');
parser_cell{9} = sprintf('tm = [];');
parser_cell{10} = sprintf('if ~isempty(getenv(''PETAKIT5D_TASK_METRICS''))\n    tm = task_metrics(''start'', functionName);\nend');
parser_cell{11} = sprintf('');
parser_cell{12} = sprintf('try');
parser_cell{13} = sprintf('    switch functionName');

t = 13;
for i = 1 : numel(funcNames)
    func_name_i = funcNames{i};
    pos_arg_num = posArgNums(i);
    param_arg_num = paramArgNums(i);
    case_line = sprintf('        case ''%s''', func_name_i);
    
    func_call_line = sprintf('            %s_parser(', func_name_i);
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
parser_cell{t + 1} = sprintf('        otherwise\n            error(''The parser for function %%s does not exist!'', functionName);');
parser_cell{t + 2} = sprintf('    end');
parser_cell{t + 3} = sprintf('catch ME');
parser_cell{t + 4} = sprintf('    if ~isempty(tm)\n        task_metrics(''stop'', tm, ME);\n    end');
parser_cell{t + 5} = sprintf('    rethrow(ME);');
parser_cell{t + 6} = sprintf('end');
parser_cell{t + 7} = sprintf('if ~isempty(tm)\n    task_metrics(''stop'', tm);\nend');
parser_cell{t + 8} = sprintf('');
parser_cell{t + 9} = sprintf('toc(t0);');
parser_cell{t + 10} = sprintf('fprintf(''\\n'');');
parser_cell{t + 11} = sprintf('');
parser_cell{t + 12} = sprintf('if ismcc\n    exit;\nend');
parser_cell{t + 13} = sprintf('');
parser_cell{t + 14} = sprintf('end');
parser_cell{t + 15} = sprintf('');
parser_cell(t + 16 : end) = [];

% write out to a file
% if the parser function exist, check if the content is the same
fnout = sprintf('%s/%s.m',  outPath, outFuncName);
if exist(fnout, 'file')
    exist_parser_cell = readTextFile(fnout);
    tmpout = sprintf('%s/%s_tmp.m',  outPath, outFuncName);
    writeTextFile(parser_cell, tmpout);
    tmp_parser_cell = readTextFile(tmpout);

    if strcmp(strjoin(tmp_parser_cell, '\n'), strjoin(exist_parser_cell, '\n'))
        fprintf('    The function %s.m remains the same, skip the update of the function.\n', outFuncName);
        delete(tmpout);
    else
        fprintf('    The function %s.m has changed, update the parser function.\n', outFuncName);
        movefile(tmpout, fnout);
    end
else
    writeTextFile(parser_cell, fnout);
end

end
