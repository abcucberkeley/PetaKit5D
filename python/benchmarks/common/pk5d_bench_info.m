function pk5d_bench_info()
% appends the MATLAB version and thread settings of this process to its
% benchmark file (PK5D_BENCH_FILE); does nothing without it.

f = getenv('PK5D_BENCH_FILE');
if isempty(f)
    return;
end
fid = fopen(f, 'a');
fprintf(fid, '{"info": {"matlab": "%s", "comp_threads": %d, "cores": %d}}\n', ...
    version, maxNumCompThreads, feature('numcores'));
fclose(fid);

end
