function bench_resample_t04(work, t)
% The resampling benchmark (T04), MATLAB side (Python: bench_resample_t04.py): this process
% pins itself (all its threads) to the CPUs of <work>/cpus_<t>.txt before any mex runs and
% sets maxNumCompThreads; then each operation once untimed and then timed:
%   frame   XR_resampleFrame of the tile by [2 2 2] (linear, uint16 tiff to uint16 tiff)
%   zarr    XR_resampleSingleZarr of the tile as zarr in MATLAB's layout by [2 2 2] (batches
%           512^3, chunks 256^3, local)

cpus = strtrim(fileread(fullfile(work, sprintf('cpus_%d.txt', t))));
[s, msg] = system(sprintf('taskset -a -p -c %s %d', strrep(cpus, ' ', ','), feature('getpid')));
if s ~= 0
    error('taskset failed: %s', msg);
end
maxNumCompThreads(t);
pk5d_bench_mark('pinned');

tile = [work, '/tile.tif'];
zarrIn = [work, '/zarr/matlab/tile.zarr'];
timed('frame', @(o) XR_resampleFrame(tile, o, [2, 2, 2]), work, '.tif', t);
timed('zarr', @(o) XR_resampleSingleZarr(zarrIn, o, [2, 2, 2], 'parseCluster', false), work, '.zarr', t);
end


function timed(op, fn, work, ext, t)
warm = [work, '/w_ml_', op, ext];
fn(warm);
remove(warm);
out = [work, '/t_ml_', op, ext];
t0 = tic;
fn(out);
dt = toc(t0);
remove(out);
pk5d_bench_record(op, t, dt);
end


function remove(p)
if isfolder(p)
    rmdir(p, 's');
elseif isfile(p)
    delete(p);
end
end
