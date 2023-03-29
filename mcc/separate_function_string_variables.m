function [func_name, var_cell] = separate_function_string_variables(funcStr)
% separate function name, variables from a single function string. The
% function string must be the standard function string as input for slurm
% wrapper
%
% TO DO: add support for "keyword=variable" type separation
%
% Author: Xiongtao Ruan (03/04/2023)
% 
% xruan (03/11/2022): fix issue for function handles with nested strings


% find out function name 
func_name = regexp(funcStr, '^(\w+)(', 'tokens');
func_name = func_name{1}{1};

% find out variable string
var_str = regexp(funcStr, sprintf('^%s\\((.+)\\)$', func_name), 'tokens');
var_str = var_str{1}{1};

% group strings in {}, [] and '' as same variables (for now, we don't
% consider escape cases)
p_11 = regexp(var_str, '{');
p_12 = regexp(var_str, '}');

p_21 = regexp(var_str, '[');
p_22 = regexp(var_str, ']');

p_31 = regexp(var_str, '''');

p_41 = regexp(var_str, '''''');
p_51 = regexp(var_str, ','''',');
p_52 = regexp(var_str, '{''''}');

p_61 = regexp(var_str, '"');

% remove ,'', from p_41
p_41(ismember(p_41, p_51 + 1)) = [];
p_41(ismember(p_41, p_52 + 1)) = [];

% remove '''' from p_31
p_31(ismember(p_31, [p_41, p_41 + 1])) = [];

if numel(p_11) ~= numel(p_12)
    error('{ and } are not fully paired.');
end

if numel(p_21) ~= numel(p_21)
    error('[ and ] are not fully paired.');
end

if rem(numel(p_31), 2) == 1
    error(' '' are not fully paired.');
end

if rem(numel(p_41), 2) == 1
    error(' '''' are not fully paired.');
end

if rem(numel(p_61), 2) == 1
    error(' " are not fully paired.');
end

p_pairs = cell(5, 1);
p_pairs{1} = pair_parenthesis(p_11', p_12');
p_pairs{2} = pair_parenthesis(p_21', p_22');
p_pairs{3} = reshape(p_31, 2, [])';
p_pairs{4} = reshape(p_41, 2, [])';
p_pairs{5} = reshape(p_61, 2, [])';

p_pairs = sortrows(cat(1, p_pairs{:}), [1, 2]);

% remove nested {}, [], '' '''' pairs
keep_flag = true(size(p_pairs, 1), 1);
for i = 1 : size(p_pairs, 1)
    if ~keep_flag(i)
        continue;
    end
    pr_i = p_pairs(i, :);
    ind_i = p_pairs(:, 1) > pr_i(1) & p_pairs(:, 2) < pr_i(2) & keep_flag;
    if any(ind_i)
        keep_flag(ind_i) = false;
    end
end

p_pairs = p_pairs(keep_flag, :);

% split strings based on the pairs
var_cell = cell(size(p_pairs, 1) * 2, 1);
ind = 1;
kpts = unique([1; p_pairs(:, 1); p_pairs(:, 2) + 1]);
for i = 1 : numel(kpts)
    s = kpts(i);
    % before first pair
    if all(s < p_pairs(:, 1))
        var_cell(ind) = var_str(s : min(s, p_pairs(:)) - 1);
        ind = ind + 1;
    end
    if any(p_pairs(:, 1) <= s & p_pairs(:, 2) > s)
        % within a pair
        ind_i = p_pairs(:, 1) <= s & p_pairs(:, 2) > s;
        var_cell{ind} = var_str(p_pairs(ind_i, 1) : p_pairs(ind_i, 2));  
        ind = ind + 1;
    elseif any(s < p_pairs(:, 1) & p_pairs(:, 1) - s ~= 1)
        % between two pairs
        ind_i = find(s < p_pairs(:, 1), 1, 'first');
        if p_pairs(ind_i, 1) - s == 1
            continue;
        end
        var_i = var_str(s : p_pairs(ind_i, 1) - 1);
        if strcmp(var_i(1) , ',')
            var_i = var_i(2 : end);
        end
        if strcmp(var_i(end) , ',')
            var_i = var_i(1 : end - 1);
        end        
        var_cell{ind} = strsplit(var_i, ',');
        ind = ind + 1;
    end
    % last pair to the end
    if all(s > p_pairs(:, 2)) && s < numel(var_str)
        var_i = var_str(s : end);
        if strcmp(var_i(1) , ',')
            var_i = var_i(2 : end);
        end
        if strcmp(var_i(end) , ',')
            var_i = var_i(1 : end - 1);
        end        
        var_cell(ind) = strsplit(var_i, ',');
        ind = ind + 1;
    end
end

var_cell(ind : end) = [];
var_cell = var_cell';
if any(cellfun(@(x) iscell(x), var_cell))
    var_cell = cat(2, var_cell{:});
end

end



function [pr] = pair_parenthesis(p_1, p_2)
% p_1: forward indices, p_2: backward indices

if isempty(p_1) && isempty(p_2)
    pr = [];
    return;
end
if numel(p_1) ~= numel(p_2) || p_1(1) >= p_2(1)
    error('The forward and backward parenthesis cannot be fully paired!')
end
pr = zeros(numel(p_1), 2);

p = sort([p_1; p_2]);
f = ismember(p, p_1);

s = zeros(numel(p_1), 1);
ind = 1;
pind = 1;

for i = 1 : numel(p)
    if f(i)
        s(ind) = p(i);
        ind = ind + 1;
    else
        ind = ind - 1;        
        pr(pind, :) = [s(ind), p(i)];
        pind = pind + 1;
    end
end

end

