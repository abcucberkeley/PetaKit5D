function pk5d_bench_record(op, threads, times, varargin)
% appends the timings (seconds per call) of an operation measured inside this
% MATLAB process to its benchmark file (PK5D_BENCH_FILE), as
% petakit5d.testing.bench.record does in Python; does nothing without it.
% Extra name-value pairs (numbers or text) are kept with the record.
%
% Usage: pk5d_bench_record('fftn', 4, [0.51, 0.50, 0.52], 'size', '512x512x500')

f = getenv('PK5D_BENCH_FILE');
if isempty(f)
    return;
end
if isempty(threads)
    th = 'null';
else
    th = sprintf('%d', threads);
end
ts = sprintf('%.9g, ', times);
s = sprintf('{"op": "%s", "threads": %s, "times": [%s]', op, th, ts(1 : end - 2));
for i = 1 : 2 : numel(varargin)
    v = varargin{i + 1};
    if ischar(v) || isstring(v)
        s = sprintf('%s, "%s": "%s"', s, varargin{i}, char(v));
    else
        s = sprintf('%s, "%s": %.17g', s, varargin{i}, v);
    end
end
fid = fopen(f, 'a');
fprintf(fid, '%s}\n', s);
fclose(fid);

end
