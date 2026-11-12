function pk5d_bench_mark(label)
% appends a time mark (seconds since 1970, UTC) to this process's benchmark
% file (PK5D_BENCH_FILE, set by petakit5d.testing.bench); does nothing
% without it. The harness puts 'first', 'ready' and 'last' around the code of
% a MATLAB variant; code may add marks of its own in between.
%
% Usage: pk5d_bench_mark('label')

f = getenv('PK5D_BENCH_FILE');
if isempty(f)
    return;
end
% MATLAB's tic counts microseconds since 1970 on Linux (as task_metrics uses it)
t = double(tic) * 1e-6;
if t < 1e9 || t > 5e9
    t = posixtime(datetime('now', 'TimeZone', 'UTC'));
end
fid = fopen(f, 'a');
fprintf(fid, '{"mark": "%s", "t": %.6f}\n', label, t);
fclose(fid);

end
