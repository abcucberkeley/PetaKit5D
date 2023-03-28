%% check whether mcc functions in this folder have the exactly same parameters as in their functions

cd(fileparts(which(mfilename)));

dir_info = dir('*_parser.m');

parser_func_names = {dir_info.name}';

%% check if all the variables in the parser functions have char validators

nF = numel(parser_func_names);

for i = 1 : nF
    cur_func_name = parser_func_names{i};
    fprintf('Function name: %s\n', cur_func_name)
    [parser_param_names, parser_char_types] = parse_function_parameters(cur_func_name);

    % check parser's main function
    main_func_name = cur_func_name(1 : end - 9);
    [main_param_names, main_char_types] = parse_function_parameters(main_func_name);

    % check if all variables contains ischar validator
    if ~all(parser_char_types)
        miss_char_param_names = parser_param_names(~parser_char_types);
        if ~isempty(miss_char_param_names)
            disp('  Parameters that are missing ischar validator:')
            disp(miss_char_param_names);
        end
    end

    % check if the parameters between parser and main functions match
    common_param_names = intersect(parser_param_names, main_param_names);
    extra_param_names = setdiff(parser_param_names, common_param_names);
    missing_param_names = setdiff(main_param_names, common_param_names);
    if ~isempty(extra_param_names)
        disp('  Extra parameters:')        
        disp(extra_param_names);
    end
    if ~isempty(missing_param_names)
        disp('  Missing parameters:')        
        disp(missing_param_names);
    end
    disp(' ');
end
