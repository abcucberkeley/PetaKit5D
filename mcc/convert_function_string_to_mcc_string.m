function [func_name, var_str] = convert_function_string_to_mcc_string(func_str)

[func_name, var_cell] = separate_function_string_variables(func_str);

for j = 1 : numel(var_cell)
    var_j = var_cell{j};
    if strcmp(var_j(1), '{')
        var_j = '"' + string(var_j) + '"';
    end
    if strcmp(var_j(1), '[') && contains(var_j, ';')
        var_j = '"' + string(var_j) + '"';
    end
    if j == 1
        var_str = string(var_j);
    else
        var_str = var_str + ' ' + string(var_j);
    end
end

end