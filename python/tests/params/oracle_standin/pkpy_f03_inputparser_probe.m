function out = pkpy_f03_inputparser_probe()
% How inputParser binds arguments (CaseSensitive false, as every PetaKit5D function sets):
% name matching, prefixes, structs, duplicates, optional positionals, the order of Results,
% and its error identifiers and messages. test_inputparser_probe.py replays the same calls
% through the registry. Parser 'num': the optional 'opt' is numeric; parser 'char': it is
% char; parser 'none': no optional. (F03 stand-in until the oracle harness of F04.)

calls = {
    'num', 'base', {1}
    'num', 'opt_given', {1, 7}
    'num', 'opt_then_pair', {1, 7, 'zeta', 5}
    'num', 'pair', {1, 'zeta', 5}
    'num', 'case', {1, 'ZETA', 5}
    'num', 'prefix', {1, 'zet', 5}
    'num', 'prefix_mixed_case', {1, 'al', 7}
    'num', 'ambiguous', {1, 'save', true}
    'num', 'exact_over_prefix', {1, 'deskew', true}
    'num', 'prefix_of_two_with_exact', {1, 'deskewR', true}
    'num', 'unknown', {1, 'foo', 1}
    'num', 'unknown_close', {1, 'zetta', 1}
    'num', 'odd', {1, 'zeta'}
    'num', 'dup', {1, 'zeta', 1, 'zeta', 2}
    'num', 'dup_case', {1, 'zeta', 1, 'ZETA', 2}
    'num', 'struct', {1, struct('zeta', 3, 'Alpha', 4)}
    'num', 'struct_then_pair', {1, struct('zeta', 3), 'zeta', 9}
    'num', 'pair_then_struct', {1, 'zeta', 9, struct('zeta', 3)}
    'num', 'struct_unknown', {1, struct('foo', 3)}
    'num', 'struct_case', {1, struct('ZETA', 3)}
    'num', 'struct_prefix', {1, struct('zet', 3)}
    'num', 'invalid', {1, 'zeta', 'abc'}
    'num', 'invalid_anon', {1, 'beta', 1}
    'num', 'valid_anon_char', {1, 'beta', 'x'}
    'num', 'name_is_required', {1, 'req', 5}
    'num', 'name_is_optional', {1, 'opt', 5}
    'num', 'opt_invalid_char', {1, 'hello', 5}
    'num', 'opt_invalid_struct', {1, struct('zeta', 3)}
    'num', 'nonchar_name', {1, 7, 5, 5}
    'num', 'string_name', {1, "zeta", 5}
    'num', 'empty_name', {1, 7, '', 5}
    'num', 'required_missing', {}
    'num', 'required_invalid', {'x'}
    'num', 'cell_param', {1, 'Gamma', {'a', 'b'}}
    'num', 'cell_param_invalid', {1, 'Gamma', 'a'}
    'num', 'logical_param_numeric', {1, 'save16bit', 1}
    'num', 'numeric_param_logical', {1, 'zeta', true}
    'num', 'numeric_param_empty', {1, 'zeta', []}
    'num', 'numeric_param_empty_char', {1, 'zeta', ''}
    'char', 'base', {1}
    'char', 'opt_given', {1, 'hello'}
    'char', 'opt_is_param_name', {1, 'zeta', 5}
    'char', 'opt_then_pair', {1, 'hello', 'zeta', 5}
    'char', 'opt_string', {1, "zeta", 5}
    'char', 'opt_nonchar', {1, 5}
    'none', 'pair', {1, 'zeta', 5}
    'none', 'unknown', {1, 'foo', 1}
    'none', 'unknown_close', {1, 'zetta', 1}
    'none', 'odd', {1, 'zeta'}
    'none', 'ambiguous', {1, 'save', true}
    'none', 'nonchar_name', {1, 5, 5}
    'none', 'struct_unknown', {1, struct('foo', 3)}
    'none', 'struct_mixed', {1, struct('zeta', 3, 'foo', 4)}
    'none', 'struct_prefix', {1, struct('zet', 3)}
    'none', 'pair_unknown_after', {1, 'zeta', 1, 'foo', 2}
    'none', 'empty_name', {1, '', 5}
    'none', 'name_is_required', {1, 'req', 5}
    'none', 'prefix_of_required', {1, 're', 5}
    'none', 'string_value', {1, 'Gamma', {"a"}}
    'none', 'too_many_positional', {1, 2}
};
out = cell(size(calls, 1), 1);
for k = 1 : size(calls, 1)
    r = struct('parser', calls{k, 1}, 'name', calls{k, 2});
    ip = make_parser(calls{k, 1});
    try
        ip.parse(calls{k, 3}{:});
        res = ip.Results;
        r.fields = fieldnames(res)';
        vals = struct();
        for f = fieldnames(res)'
            vals.(f{1}) = describe(res.(f{1}));
        end
        r.values = vals;
        r.using_defaults = ip.UsingDefaults;
    catch err
        r.error = err.identifier;
        r.message = err.message;
    end
    out{k} = r;
end
end


function ip = make_parser(kind)
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('req', @isnumeric);
if strcmp(kind, 'num')
    ip.addOptional('opt', 0, @isnumeric);
elseif strcmp(kind, 'char')
    ip.addOptional('opt', 'dflt', @ischar);
end
ip.addParameter('zeta', 0, @isnumeric);
ip.addParameter('Alpha', 0, @isnumeric);
ip.addParameter('beta', false, @(x) islogical(x) || ischar(x));
ip.addParameter('save16bit', true, @islogical);
ip.addParameter('saveZarr', false, @islogical);
ip.addParameter('deskew', false, @islogical);
ip.addParameter('deskewRotate', false, @islogical);
ip.addParameter('Gamma', {}, @iscell);
end


function d = describe(v)
d = struct('c', class(v));
d.s = num2cell(size(v));
if isfloat(v)
    d.d = arrayfun(@(t) num2hex(t), v(:)', 'UniformOutput', false);
elseif islogical(v)
    d.d = num2cell(double(v(:)'));
elseif ischar(v)
    d.d = {v(:)'};
elseif iscell(v)
    d.d = cellfun(@describe, v(:)', 'UniformOutput', false);
elseif isstring(v)
    d.d = cellstr(v(:)');
end
end
