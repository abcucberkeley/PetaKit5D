function bench_c04a(work)
% The C04a benchmark, MATLAB side (Python: bench_c04a.py): on a fresh copy of the master
% condition, timed one after the other: XR_loadConditionData3D, deskewData with
% XR_runDetTrack3d's options, XR_correctXZoffsetData3D.

run = sprintf('%s/run_ml_%d', work, feature('getpid'));
copyfile([work '/master'], run);
movie = [run '/Ex1_demo_z0p5/'];
LLC = {'LLFFCorrection', false, 'LowerLimit', 0.4, 'LSImageCh1', '', 'LSImageCh2', '', ...
    'LSImageCh3', '', 'LSImageCh4', '', 'BackgroundCh1', '', 'BackgroundCh2', '', ...
    'BackgroundCh3', '', 'BackgroundCh4', ''};
t0 = tic;
data = XR_loadConditionData3D(movie, {'ch1', 'ch2'}, {'gfp', 'rfp'});
t1 = toc(t0);
t0 = tic;
data = deskewData(data, 'Overwrite', false, 'SkewAngle', 31.5, 'aname', 'Analysis', ...
    'Rotate', false, 'sCMOSCameraFlip', false, 'Crop', false, LLC{:}, 'zoffsetCorrection', false);
t2 = toc(t0);
t0 = tic;
XR_correctXZoffsetData3D(data, 'Overwrite', false);
t3 = toc(t0);
pk5d_bench_record('load', [], t1);
pk5d_bench_record('deskew', [], t2);
pk5d_bench_record('offsets', [], t3);
rmdir(run, 's');
end
