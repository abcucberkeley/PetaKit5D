function [values, header] = pkpy_j12_master(ipcode, drop, functionName, varargin)
% mccMaster with command-line strings, the task function replaced by a recorder: the values
% its mcc parser hands the function, and the "Function name / Variables" lines it printed.
% The recorded arguments are then checked by the function's own inputParser (ipcode: MATLAB
% code that makes it as ip, from the function's inputParser block), as the function's first
% step would check them.
%
% values: a struct with pos1..posN (the arguments the parser passes by position) and one
% field per name-value pair the parser passes, named as the parser names them; names in
% drop (site-specific defaults) are left out. Vectors and cells are rows (the port's Python
% values are tuples, which have no orientation).

d = tempname;
mkdir(d);
fid = fopen(fullfile(d, [functionName, '.m']), 'w');
fprintf(fid, 'function varargout = %s(varargin)\npkpy_j12_record(''set'', varargin);\nvarargout = cell(1, nargout);\nend\n', functionName);
fclose(fid);
addpath(d, '-begin');
clear(functionName);
c = onCleanup(@() cleanup(d, functionName));
pkpy_j12_record('set', []);

txt = evalc('mccMaster(functionName, varargin{:})');
got = pkpy_j12_record('get');
eval(ipcode);
ip.parse(got{:});

n = nargin([functionName, '_parser']);
if n < 0
    n = -n - 1;
end
values = struct();
for i = 1 : n
    values.(sprintf('pos%d', i)) = rows(got{i});
end
for i = n + 1 : 2 : numel(got)
    name = char(got{i});
    if ~any(strcmpi(name, drop))
        values.(name) = rows(got{i + 1});
    end
end
header = regexp(txt, 'Function name:[^\n]*\nVariables:[^\n]*\n', 'match', 'once');

end


function v = rows(v)
if iscell(v)
    for i = 1 : numel(v)
        v{i} = rows(v{i});
    end
end
if (isnumeric(v) || islogical(v) || iscell(v)) && iscolumn(v) && ~isscalar(v)
    v = v.';
end
end


function cleanup(d, functionName)
rmpath(d);
rmdir(d, 's');
clear(functionName);
pkpy_j12_record('set', []);
end
