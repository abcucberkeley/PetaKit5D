function puncta_b08_ml(folder, threads, repeats)
% Puncta removal benchmark (B08), the MATLAB side: XR_ExM_PunctaRemoval_zarr of
% folder/expanded.zarr (defaults, parseCluster false) with maxNumCompThreads set to each
% thread count; the result folder removed before each run (untimed).

nt0 = maxNumCompThreads;
fn = [folder, '/expanded.zarr'];
out = [folder, '/Puncta_Removed'];
for t = threads
    maxNumCompThreads(t);
    times = zeros(1, repeats);
    for r = 1 : repeats
        if exist(out, 'dir')
            rmdir(out, 's');
        end
        t0 = tic;
        XR_ExM_PunctaRemoval_zarr(fn, 'parseCluster', false);
        times(r) = toc(t0);
    end
    pk5d_bench_record('puncta_removal_zarr', t, times);
end
maxNumCompThreads(nt0);
end
