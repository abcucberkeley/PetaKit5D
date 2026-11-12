function [summary] = run_pack(packFn)
% run a pack of oracle cases in this MATLAB session (python/docs/testing.md)
%
% packFn: pack.json written by petakit5d.testing.oracle. For each case the
% outputs, the files it wrote, its warnings and its error go to a staging
% folder, which is published to the golden store with one rename: a stopped
% session leaves no partial case, and of two sessions that run the same case
% the first rename wins.
%
% Call it after setup([]) of the oracle tree, with the folder holding
% +pkpy_oracle on the path (the package shadows nothing).


tStart = tic;
pack = jsondecode(fileread(packFn));
cases = pack.cases;
if isstruct(cases)
    cases = num2cell(cases);
end
nc = numel(cases);

try
    host = char(java.net.InetAddress.getLocalHost().getHostName());
catch
    host = getenv('HOSTNAME');
end
info = struct('host', host, 'matlab_version', version, 'threads', maxNumCompThreads, ...
    'job', getenv('SLURM_JOB_ID'), 'pack', pack.pack, 'oracle_commit', pack.oracle_commit, ...
    'oracle_dir', pack.oracle_dir);
try
    feature('hotlinks', false);
catch
end
ws0 = warning;
pwd0 = pwd;
% work from the local work root: returning to a large folder on network storage after
% each case costs MATLAB a rescan of it
if ~isfolder(pack.work_root)
    mkdir(pack.work_root);
end
cd(pack.work_root);

fprintf('PKPY pack %s: %d cases, MATLAB %s, %d threads\n', pack.pack, nc, version, maxNumCompThreads);
rows = cell(1, nc);
for i = 1 : nc
    fprintf('PKPY start %s\n', cases{i}.id);
    rows{i} = pkpy_oracle.run_case(cases{i}, pack, info, ws0, pack.work_root);
    r = rows{i};
    fprintf('PKPY %-13s %s (%.1f ms, call %.1f ms)%s\n', r.status, r.id, 1000 * r.t_case, 1000 * r.t_call, ...
        pkpy_oracle.ternary(isempty(r.message), '', [': ', r.message]));
end
cd(pwd0);
if isfolder(pack.staging) && numel(dir(pack.staging)) <= 2
    rmdir(pack.staging);
end

summary = struct('pack', pack.pack, 'n', nc, 't_total', toc(tStart), 'matlab_version', version, ...
    'threads', maxNumCompThreads, 'host', host);
summary.cases = rows;
fid = fopen(pack.summary, 'w');
if fid >= 0
    fwrite(fid, unicode2native(jsonencode(summary), 'UTF-8'), 'uint8');
    fclose(fid);
end
st = cellfun(@(r) r.status, rows, 'unif', 0);
fprintf('PKPY done: %d published (%d MATLAB errors), %d cached, %d lost races, %d harness errors in %.1f s\n', ...
    sum(ismember(st, {'ok', 'error'})), sum(strcmp(st, 'error')), sum(strcmp(st, 'cached')), ...
    sum(strcmp(st, 'lost_race')), sum(strcmp(st, 'harness_error')), toc(tStart));

end
