function [x] = fill_placeholders(x, work)
% replace the placeholders in a case's arguments: Work(sub) by the work
% folder with a trailing slash and sub, InPath(name) by the path of the input
% copied there (with a trailing slash for folders)


if iscell(x)
    for i = 1 : numel(x)
        x{i} = pkpy_oracle.fill_placeholders(x{i}, work);
    end
elseif isstruct(x) && isscalar(x) && isfield(x, 'pkpy__ph')
    root = strrep(work, '\', '/');
    switch x.pkpy__ph
        case 'work'
            x = [root, '/', x.sub];
        case 'path'
            p = [root, '/', x.n];
            if isfolder(p)
                p = [p, '/'];
            end
            x = p;
        otherwise
            error('pkpy:oracle:placeholder', 'unknown placeholder %s', x.pkpy__ph);
    end
elseif isstruct(x)
    f = fieldnames(x);
    for e = 1 : numel(x)
        for j = 1 : numel(f)
            x(e).(f{j}) = pkpy_oracle.fill_placeholders(x(e).(f{j}), work);
        end
    end
end

end
