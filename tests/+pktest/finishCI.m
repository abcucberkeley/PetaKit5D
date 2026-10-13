function [] = finishCI(results, outputDir)
% record the test outcome for CI and end the MATLAB session. On Windows and
% macOS, MATLAB hangs (Windows) or crashes (macOS) during its shutdown after
% the PetaKit5D mex files have been used (see tests/README.md), so there the
% process is stopped right away, and the CI job reads the outcome from
% <outputDir>/status.txt instead of MATLAB's exit code.
%
%    results : matlab.unittest.TestResult array from run_tests.
%  outputDir : char (default: 'test-results'). Folder for status.txt.

if nargin < 2
    outputDir = 'test-results';
end
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

nFailed = nnz([results.Failed]);
fid = fopen(fullfile(outputDir, 'status.txt'), 'w');
fprintf(fid, 'tests %d\npassed %d\nfailed %d\nincomplete %d\n', numel(results), nnz([results.Passed]), ...
    nFailed, nnz([results.Incomplete] & ~[results.Failed]));
fclose(fid);

if ispc
    system(sprintf('taskkill /F /PID %d', feature('getpid')));
elseif ismac
    system(sprintf('kill -9 %d', feature('getpid')));
end

% Linux shuts down normally
assertSuccess(results);

end
