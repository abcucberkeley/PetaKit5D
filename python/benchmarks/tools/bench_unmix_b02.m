function bench_unmix_b02(work, t, doDataset)
% The unmixing benchmark (B02), MATLAB side (Python: bench_unmix_b02.py): this process pins
% itself (all its threads) to the CPUs of <work>/cpus_<t>.txt before any mex runs and sets
% maxNumCompThreads; then each operation once untimed and then timed:
%   frame         XR_unmix_channels_frame of one demo tile pair (tiff, [1 -0.1], backgrounds
%                 [100 100], offset 1, uint16)
%   zarr          XR_unmix_channels_zarr of the same pair as zarr in MATLAB's layout (batches
%                 1024^3, chunks 256^3, local)
%   dataset_demo  (doDataset, timed once) XR_unmix_channels_data_wrapper of the demo's 8 tile
%                 pairs with the demo's settings (local)

cpus = strtrim(fileread(fullfile(work, sprintf('cpus_%d.txt', t))));
[s, msg] = system(sprintf('taskset -a -p -c %s %d', strrep(cpus, ' ', ','), feature('getpid')));
if s ~= 0
    error('taskset failed: %s', msg);
end
maxNumCompThreads(t);
pk5d_bench_mark('pinned');

demo = [work, '/demo/'];
tile = 'Scan_Iter_0000_0000_%s_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_003y_000z_0000t';
frames = cellfun(@(c) [demo, sprintf(tile, c), '.tif'], {'CamA', 'CamB'}, 'unif', 0);
zarrDir = [work, '/zarr/matlab/'];
zarrs = cellfun(@(c) [zarrDir, sprintf(tile, c), '.zarr'], {'CamA', 'CamB'}, 'unif', 0);

frame = @(name) XR_unmix_channels_frame(frames, [1, -0.1], 'resultDirName', name, ...
    'constBackground', [100, 100], 'constOffset', 1, 'save16bit', true);
zarrf = @(name) XR_unmix_channels_zarr(zarrs, [1, -0.1], 'resultDirName', name, 'parseCluster', false);
dataset = @(name) XR_unmix_channels_data_wrapper(demo, 'unmixFactors', [1, -0.1], 'mode', 'linear', ...
    'resultDirName', name, 'channelPatterns', {'CamA', 'CamB'}, 'channelInd', 1, ...
    'constBackground', [100, 100], 'constOffset', 1, 'save16bit', true, 'parseCluster', false, ...
    'masterCompute', true);

timed('frame', frame, demo, t, true);
timed('zarr', zarrf, zarrDir, t, true);
if doDataset
    timed('dataset_demo', dataset, demo, t, false);
end
end


function timed(op, fn, folder, t, warm)
if warm
    fn('Uw_ml');
    rmdir([folder, 'Uw_ml'], 's');
end
t0 = tic;
fn('Ut_ml');
dt = toc(t0);
rmdir([folder, 'Ut_ml'], 's');
pk5d_bench_record(op, t, dt);
end
