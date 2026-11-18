function n5_b09_ml(folder, threads, repeats)
% N5 benchmark (B09), the MATLAB side, with maxNumCompThreads set to each thread count:
% zarrToN5 of folder/t.zarr and N5ToZarr of folder/s.n5 (defaults, parseCluster false), each
% result folder removed before each run (untimed). N5Adapter needs zarr-python 2 and numcodecs:
% the interpreter of the active virtual environment.

pe = pyenv;
if pe.Status ~= "Loaded"
    pyenv('Version', fullfile(getenv('VIRTUAL_ENV'), 'bin', 'python'));
end
insert(py.sys.path, int64(0), fullfile(fileparts(which('setup')), 'microscopeDataProcessing', 'python'));
py.importlib.import_module('zarrAPI');

nt0 = maxNumCompThreads;
for t = threads
    maxNumCompThreads(t);
    tz = zeros(1, repeats);
    tn = zeros(1, repeats);
    for r = 1 : repeats
        if exist([folder, '/n5'], 'dir')
            rmdir([folder, '/n5'], 's');
        end
        t0 = tic;
        zarrToN5([folder, '/t.zarr'], 'parseCluster', false);
        tz(r) = toc(t0);
        if exist([folder, '/zarr'], 'dir')
            rmdir([folder, '/zarr'], 's');
        end
        t0 = tic;
        N5ToZarr([folder, '/s.n5'], 'parseCluster', false);
        tn(r) = toc(t0);
    end
    pk5d_bench_record('zarr_to_n5', t, tz);
    pk5d_bench_record('n5_to_zarr', t, tn);
end
maxNumCompThreads(nt0);
end
