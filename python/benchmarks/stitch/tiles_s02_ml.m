function tiles_s02_ml(folder, threads, repeats)
% Tile preparation benchmark (S02), the MATLAB side: stitch_process_tiles of the 16 tiles with
% the flat-field and erosion function file, tiff tiles to zarr (tiff_tiles) and zarr tiles
% (zarr_tiles), parseCluster false, with maxNumCompThreads set to each thread count; the
% result folders removed before each run (untimed). Times are per tile.

nt0 = maxNumCompThreads;
d = dir([folder, '/tif/*.tif']);
tifs = cellfun(@(x) [folder, '/tif/', x], sort({d.name}), 'unif', 0);
d = dir([folder, '/zarr/matlab/*.zarr']);
zarrs = cellfun(@(x) [folder, '/zarr/matlab/', x], sort({d.name}), 'unif', 0);
common = {'blockSize', [250, 250, 125], 'processFunPath', {[folder, '/fun.mat']}, 'parseCluster', false};
runs = {'tiff_tiles', tifs, [folder, '/tif/zarr'], {'zarrFile', false, 'resultDirName', 'zarr'}; ...
    'zarr_tiles', zarrs, [folder, '/zarr/matlab/proc'], {'zarrFile', true, 'resultDirName', 'proc'}};
for t = threads
    maxNumCompThreads(t);
    for k = 1 : size(runs, 1)
        times = zeros(1, repeats);
        for r = 1 : repeats
            if exist(runs{k, 3}, 'dir')
                rmdir(runs{k, 3}, 's');
            end
            t0 = tic;
            stitch_process_tiles(runs{k, 2}, common{:}, runs{k, 4}{:});
            times(r) = toc(t0) / numel(runs{k, 2});
        end
        pk5d_bench_record(runs{k, 1}, t, times);
    end
end
maxNumCompThreads(nt0);
end
