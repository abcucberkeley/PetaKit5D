function [s] = open_handles()
% what could hold this MATLAB open, as one line (exit probe of the CI job
% exit-probe; not a test): the parallel pool (gcp('nocreate')) and the
% background pool with their futures not done, the timers, and Java's live
% threads that are not daemons. Nothing here starts a pool: the background
% pool is looked up in MATLAB's pool manager, as backgroundPool would start
% one.

s = {};
try
    p = gcp('nocreate');
    if isempty(p)
        s{end + 1} = 'no parallel pool';
    else
        s{end + 1} = sprintf('parallel pool (%s) of %d workers, connected %d, futures queued %d, running %d', ...
            class(p), p.NumWorkers, p.Connected, numel(p.FevalQueue.QueuedFutures), numel(p.FevalQueue.RunningFutures));
    end
catch ME
    s{end + 1} = ['parallel pool: ', ME.identifier];
end
try
    pm = parallel.internal.pool.PoolManager.getInstance();
    b = pm.getAllPools(parallel.internal.pool.PoolApiTag.Background);
    if isempty(b)
        s{end + 1} = 'no background pool';
    end
    for k = 1 : numel(b)
        s{end + 1} = sprintf('background pool of %d workers, busy %d, futures queued %d, running %d', b(k).NumWorkers, ...
            b(k).Busy, numel(b(k).FevalQueue.QueuedFutures), numel(b(k).FevalQueue.RunningFutures)); %#ok<AGROW>
    end
catch ME
    s{end + 1} = ['background pool: ', ME.identifier];
end
t = timerfindall;
s{end + 1} = sprintf('%d timers', numel(t));
for k = 1 : numel(t)
    s{end} = sprintf('%s, %s (running %s)', s{end}, t(k).Name, t(k).Running);
end
try
    th = java.lang.Thread.getAllStackTraces().keySet().toArray();
    names = {};
    for k = 1 : numel(th)
        if th(k).isAlive() && ~th(k).isDaemon()
            names{end + 1} = char(th(k).getName()); %#ok<AGROW>
        end
    end
    s{end + 1} = sprintf('%d live non-daemon Java threads (%s)', numel(names), strjoin(names, ', '));
catch ME
    s{end + 1} = ['Java threads: ', ME.identifier];
end
s = strjoin(s, '; ');

end
