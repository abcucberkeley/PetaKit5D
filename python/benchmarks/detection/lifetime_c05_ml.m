function lifetime_c05_ml(folder, aname, repeats)
% Lifetime analysis benchmark (C05), the MATLAB side: runLifetimeAnalysis3D with 'overwrite'
% true (the tracks read again each run) and plotIntensityCohorts on the demo condition prepared
% by lifetime_c05_data.py, each run timed whole, as the multiple-experiments demo calls them
% (the condition loaded first, untimed). Figures are closed after each run (untimed).

d = dir(folder);
d = d([d.isdir] & ~startsWith({d.name}, '.'));
cond = [fullfile(folder, d(1).name), filesep];
data = XR_loadConditionData3D(cond, {'ch1', 'ch2'}, {'gfp', 'rfp'});
nt = 1;
t1 = zeros(1, repeats);
t2 = zeros(1, repeats);
for r = 1 : repeats
    t0 = tic;
    lftRes = runLifetimeAnalysis3D(data, 'AnalysisPath', aname, 'overwrite', true); %#ok<NASGU>
    t1(r) = toc(t0);
    close all force
    t0 = tic;
    cohorts = plotIntensityCohorts(data, 'AnalysisPath', aname); %#ok<NASGU>
    t2(r) = toc(t0);
    close all force
end
pk5d_bench_record('lifetime_analysis', nt, t1);
pk5d_bench_record('intensity_cohorts', nt, t2);
end
