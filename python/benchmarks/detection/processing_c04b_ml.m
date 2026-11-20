function processing_c04b_ml(folder, datamat, repeats)
% Track processing benchmark (C04b), the MATLAB side: XR_runTrackProcessing3D of folder/Analysis
% for the movie of datamat ('data', whose DS frames are read in place) with the demos' options,
% then rotateTracks3D; the results removed before each run (untimed).

try
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false;  % the Amira writer's parfor runs in this process, as Python's
catch
end
rp = [folder, '/Analysis/'];
S = load(datamat, 'data');
data = S.data;
clear S
t_proc = zeros(1, repeats);
t_rot = zeros(1, repeats);
for r = 1 : repeats
    for f = {'ProcessedTracks.mat', 'RotatedTracks.mat'}
        if exist([rp, f{1}], 'file')
            delete([rp, f{1}]);
        end
    end
    if exist([rp, 'Amira'], 'dir')
        rmdir([rp, 'Amira'], 's');
    end
    t0 = tic;
    XR_runTrackProcessing3D(data, 'Overwrite', true, 'TrackerOutput', 'trackedFeatures.mat', ...
        'FileName', 'ProcessedTracks.mat', 'Buffer', [3, 3], 'BufferAll', false, ...
        'FitMixtures', false, 'WindowSize', [], 'ResultsPath', {rp});
    t_proc(r) = toc(t0);
    t1 = tic;
    rotateTracks3D(data, 'Overwrite', true, 'Crop', false, 'ResultsPath', {rp});
    t_rot(r) = toc(t1);
end
pk5d_bench_record('process', maxNumCompThreads, t_proc);
pk5d_bench_record('rotate', maxNumCompThreads, t_rot);
end
