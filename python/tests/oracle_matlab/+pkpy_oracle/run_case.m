function [r] = run_case(c, pack, info, ws0, home)
% run one oracle case and publish its result to the golden store (see run_pack);
% home: the folder the session returns to after the call (the local work root)
%
% Status: 'ok' or 'error' (published; 'error' when the call raised),
% 'cached' (already in the store), 'lost_race' (another session published
% it first), 'harness_error' (nothing published).


tCase = tic;
r = struct('key', c.key, 'id', c.id, 'status', '', 't_call', NaN, 't_case', NaN, 'message', '');
final = fullfile(pack.store, c.key(1 : 2), c.key);
if isfolder(final)
    r.status = 'cached';
    r.t_case = toc(tCase);
    return;
end
stage = fullfile(pack.staging, c.key);
work = fullfile(pack.work_root, c.key(1 : 16));

try
    if isfolder(work)
        remove_dir(work);
    end
    mkdir_p(stage);
    mkdir_p(work);
    copies = as_cell(c.copy);
    for j = 1 : numel(copies)
        cp = copies{j};
        dst = fullfile(work, cp.dst);
        if strcmp(cp.kind, 'tree')
            mkdir_p(dst);
            copyfile(cp.src, dst);
        else
            mkdir_p(fileparts(dst));
            copyfile(cp.src, dst);
        end
    end
    if startsWith(c.device, 'gpu:')
        g = gpuDevice();
        if ~contains(lower(g.Name), lower(c.device(5 : end)))
            error('pkpy:oracle:device', 'the case needs a %s GPU; this session has %s', c.device(5 : end), g.Name);
        end
    end
    args = pkpy_oracle.read_value(c.args);
    args = pkpy_oracle.fill_placeholders(args, work);
    before = pkpy_oracle.snapshot(work);
catch ME
    r = harness_error(r, ME, stage, work, tCase);
    return;
end

% the state a test starts in (as pktest.TestCase sets it)
reset_rng();
warning(ws0);
warning('off', 'backtrace');
warning('on', 'verbose');
lastwarn('');
clear('decon_lucy_function', 'decon_lucy_omw_function');
figs = get(groot, 'Children');
cd(work);
try
    [outs, err, txt, tCall] = pkpy_oracle.call(c, args, work);
catch ME
    cd(home);
    warning(ws0);
    r = harness_error(r, ME, stage, work, tCase);
    return;
end
[lastMsg, lastId] = lastwarn;
cd(home);
warning(ws0);
fclose('all');
figsNow = get(groot, 'Children');
if ~isempty(figsNow)
    newFigs = setdiff(figsNow, figs);
    if ~isempty(newFigs)
        close(newFigs, 'force');
    end
end
r.t_call = tCall;

try
    files = {};
    deleted = {};
    if c.capture
        [files, deleted] = pkpy_oracle.changed_files(before, pkpy_oracle.snapshot(work));
    end
    funcFile = '';
    if ~isempty(c.func)
        funcFile = func_file(c.func, info.oracle_dir);
    end
    % warnings as displayed (MATLAB wraps long ones at the window width); the last one exactly
    warns = pkpy_oracle.parse_warnings(txt, lastMsg, lastId);
    meta = struct('format', pack.format, 'key', c.key, 'id', c.id, 'func', c.func, 'func_file', funcFile, ...
        'nargout', c.nargout, 'status', 'ok', 'error', [], 'warnings', {warns}, ...
        'capture', logical(c.capture), 'files', {files}, 'deleted', {deleted}, 'decoded', {{}}, ...
        'decode_errors', {{}}, 'work', strrep(work, '\', '/'), 'host', info.host, ...
        'matlab_version', info.matlab_version, 'threads', info.threads, 'job', info.job, 'pack', info.pack, ...
        'oracle_commit', info.oracle_commit, 'device', c.device, 'gpu', '', ...
        'date', char(java.time.LocalDateTime.now().withNano(0).toString()), 't_call', tCall, ...
        'last_warning', struct('identifier', lastId, 'message', lastMsg));
    if startsWith(c.device, 'gpu')
        try
            g = gpuDevice();
            meta.gpu = g.Name;
        catch
        end
    end
    if isempty(err)
        pkpy_oracle.write_value(fullfile(stage, 'outputs.pkv'), outs);
    else
        meta.status = 'error';
        meta.error = struct('identifier', err.identifier, 'message', err.message, 'stack', {stack_text(err, info)});
    end
    for j = 1 : numel(files)
        src = fullfile(work, files{j});
        dst = fullfile(stage, 'files', files{j});
        if isfolder(src)
            mkdir_p(dst);
            copyfile(src, dst);
        else
            mkdir_p(fileparts(dst));
            copyfile(src, dst);
            if endsWith(lower(files{j}), '.mat')
                try
                    S = load(src);
                    dfn = fullfile(stage, 'decoded', [files{j}, '.pkv']);
                    mkdir_p(fileparts(dfn));
                    pkpy_oracle.write_value(dfn, S);
                    meta.decoded{end + 1} = files{j};
                catch MEd
                    meta.decode_errors{end + 1} = sprintf('%s: %s', files{j}, MEd.message);
                end
            end
        end
    end
    if ~isempty(txt)
        write_text(fullfile(stage, 'stdout.txt'), txt);
    end
    write_text(fullfile(stage, 'meta.json'), jsonencode(meta));

    % publish: one rename; a second writer finds the result there and drops its own
    mkdir_p(fileparts(final));
    if java.io.File(stage).renameTo(java.io.File(final))
        r.status = meta.status;
    elseif exist(fullfile(final, 'meta.json'), 'file')
        r.status = 'lost_race';
        remove_dir(stage);
    else
        error('pkpy:oracle:publish', 'cannot rename %s to %s', stage, final);
    end
catch ME
    r = harness_error(r, ME, stage, work, tCase);
    return;
end
remove_dir(work);
r.t_case = toc(tCase);

end


function [f] = func_file(name, oracleDir)
% where a function comes from (relative to the oracle tree), cached for the session

persistent cache
if isempty(cache)
    cache = containers.Map('KeyType', 'char', 'ValueType', 'char');
end
if ~isKey(cache, name)
    f = strrep(which(name), '\', '/');
    cache(name) = strrep(f, [strrep(oracleDir, '\', '/'), '/'], '<ORACLE>/');
end
f = cache(name);

end


function [] = reset_rng()
% rng(0, 'twister'), without building a new stream when the global one is a
% default twister already

s = RandStream.getGlobalStream();
if strcmp(s.Type, 'mt19937ar') && strcmp(s.NormalTransform, 'Ziggurat') && ~s.Antithetic && s.FullPrecision
    reset(s, 0);
else
    rng(0, 'twister');
end

end


function [r] = harness_error(r, ME, stage, work, tCase)

r.status = 'harness_error';
r.message = ME.message;
fprintf('PKPY harness error in %s:\n%s\n', r.id, getReport(ME, 'extended', 'hyperlinks', 'off'));
remove_dir(stage);
remove_dir(work);
r.t_case = toc(tCase);

end


function [s] = stack_text(err, info)
% the error's stack without the harness frames, paths relative to the oracle tree

s = {};
root = [strrep(info.oracle_dir, '\', '/'), '/'];
for i = 1 : numel(err.stack)
    f = strrep(err.stack(i).file, '\', '/');
    if contains(f, '+pkpy_oracle/')
        break;
    end
    s{end + 1} = sprintf('%s:%d %s', strrep(f, root, '<ORACLE>/'), err.stack(i).line, err.stack(i).name); %#ok<AGROW>
end

end


function [] = write_text(fn, txt)

fid = fopen(fn, 'w');
if fid < 0
    error('pkpy:oracle:write', 'cannot write %s', fn);
end
fwrite(fid, unicode2native(txt, 'UTF-8'), 'uint8');
fclose(fid);

end


function [] = mkdir_p(d)
% (MATLAB's mkdir costs about 2.5 ms on network storage, Java's 0.8 ms)

if ~isempty(d) && ~isfolder(d)
    java.io.File(d).mkdirs();
    if ~isfolder(d)
        error('pkpy:oracle:mkdir', 'cannot create %s', d);
    end
end

end


function [] = remove_dir(d)
% (an empty folder is removed through Java: rmdir costs about 2 ms)

if isfolder(d)
    if numel(dir(d)) <= 2 && java.io.File(d).delete()
        return;
    end
    rmdir(d, 's');
end

end


function [c] = as_cell(x)

if isempty(x)
    c = {};
elseif iscell(x)
    c = x;
else
    c = num2cell(x);
end

end
