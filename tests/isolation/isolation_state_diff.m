function [changes] = isolation_state_diff(a, b)
% what changed between two isolation_state snapshots: a cell array of short
% descriptions (empty if nothing a new process would differ in changed).
% Memory and newly loaded mex files are reported as well, as information.


changes = {};
if a.threads ~= b.threads
    changes{end + 1} = sprintf('threads %d -> %d', a.threads, b.threads);
end
if ~isequal(a.warnings, b.warnings)
    changes{end + 1} = ['warnings: ', warning_diff(a.warnings, b.warnings)];
end
if ~strcmp(a.lastwarnId, b.lastwarnId) || ~strcmp(a.lastwarnMsg, b.lastwarnMsg)
    changes{end + 1} = sprintf('lastwarn: %s', b.lastwarnId);
end
if ~isequal(a.rng, b.rng)
    changes{end + 1} = sprintf('rng: %s seed %d -> %s seed %d (state %s)', a.rng.Type, a.rng.Seed, b.rng.Type, ...
        b.rng.Seed, string(isequal(a.rng.State, b.rng.State)));
end
if ~strcmp(a.pwd, b.pwd)
    changes{end + 1} = sprintf('pwd %s -> %s', a.pwd, b.pwd);
end
if ~strcmp(a.path, b.path)
    pa = strsplit(a.path, pathsep);
    pb = strsplit(b.path, pathsep);
    changes{end + 1} = sprintf('path: +%d -%d entries (%s)', numel(setdiff(pb, pa)), numel(setdiff(pa, pb)), ...
        strjoin([setdiff(pb, pa), strcat('-', setdiff(pa, pb))], ', '));
end
if ~isequal(a.env, b.env)
    ka = strcat(a.env(:, 1), '=', a.env(:, 2));
    kb = strcat(b.env(:, 1), '=', b.env(:, 2));
    names = unique([a.env(ismember(ka, setdiff(ka, kb)), 1); b.env(ismember(kb, setdiff(kb, ka)), 1)]);
    changes{end + 1} = ['env: ', strjoin(names, ', ')];
end
if ~isequal(a.globals, b.globals)
    changes{end + 1} = ['globals: ', strjoin(setxor(a.globals, b.globals), ', ')];
end
if ~isequal(a.baseVars, b.baseVars)
    changes{end + 1} = ['base workspace: ', strjoin(setxor(a.baseVars, b.baseVars), ', ')];
end
if ~isequal(a.files, b.files)
    changes{end + 1} = sprintf('open files %d -> %d', size(a.files, 1), size(b.files, 1));
end
if ~isequaln(a.fds, b.fds)
    changes{end + 1} = sprintf('file descriptors %d -> %d', a.fds, b.fds);
end
if ~isequaln(a.figures, b.figures)
    changes{end + 1} = sprintf('figures %d -> %d', a.figures, b.figures);
end
if ~isequaln(a.timers, b.timers)
    changes{end + 1} = sprintf('timers %d -> %d', a.timers, b.timers);
end
if a.pool ~= b.pool
    changes{end + 1} = sprintf('parallel pool %d -> %d', a.pool, b.pool);
end
newMex = setdiff(b.mex, a.mex);
if ~isempty(newMex)
    [~, names] = cellfun(@fileparts, newMex, 'unif', 0);
    changes{end + 1} = ['(info) mex loaded: ', strjoin(names, ', ')];
end
if ~isequaln(a.gpuMemMB, b.gpuMemMB)
    changes{end + 1} = sprintf('GPU memory %g -> %g MB', a.gpuMemMB, b.gpuMemMB);
end
if isfinite(a.rssGB) && isfinite(b.rssGB) && abs(b.rssGB - a.rssGB) > 0.05
    changes{end + 1} = sprintf('(info) RSS %.2f -> %.2f GB', a.rssGB, b.rssGB);
end

end


function [str] = warning_diff(a, b)

ids = unique([{a.identifier}, {b.identifier}]);
parts = {};
for i = 1 : numel(ids)
    sa = a(strcmp({a.identifier}, ids{i}));
    sb = b(strcmp({b.identifier}, ids{i}));
    va = 'default';
    vb = 'default';
    if ~isempty(sa)
        va = sa(1).state;
    end
    if ~isempty(sb)
        vb = sb(1).state;
    end
    if ~strcmp(va, vb)
        parts{end + 1} = sprintf('%s %s -> %s', ids{i}, va, vb); %#ok<AGROW>
    end
end
str = strjoin(parts, ', ');

end
