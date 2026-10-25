function [] = nestedFrameworkTask(inputFn, outputFn, logFn)
% task function for job framework tests that runs two framework calls of
% its own (locally, with task metrics; each copies inputFn to a file, the
% second to outputFn), as large-file tasks run their steps

[outPath, name, ext] = fileparts(outputFn);
outFns = {sprintf('%s/%s_step1%s', outPath, name, ext), outputFn};
for i = 1 : 2
    funcStr = sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''copy'')', inputFn, outFns{i}, logFn);
    generic_computing_frameworks_wrapper({inputFn}, outFns(i), {funcStr}, 'parseCluster', false, ...
        'masterCompute', true, 'taskMetrics', true);
end

end
