function [] = update_golden_data(scenarios)
% regenerate the reference outputs compared by GoldenOutputTest. Only do
% this when an output change is intended, preferably on Linux, and commit
% the new tests/golden/data/*.mat files together with that change.
%
%   update_golden_data                         % all scenarios
%   update_golden_data('decon_omw')            % one scenario
%   update_golden_data({'decon_omw', 'psf_gen_masked'})

goldenDir = fileparts(mfilename('fullpath'));
testDir = fileparts(goldenDir);
addpath(testDir, goldenDir);
setup(pktest.repoRoot());

mc = ?GoldenOutputTest;
allScenarios = mc.PropertyList(strcmp({mc.PropertyList.Name}, 'scenario')).DefaultValue;
if nargin < 1 || isempty(scenarios)
    scenarios = allScenarios;
end
scenarios = cellstr(scenarios);
unknown = setdiff(scenarios, allScenarios);
if ~isempty(unknown)
    error('Unknown scenario(s): %s', strjoin(unknown, ', '));
end

dataDir = fullfile(goldenDir, 'data');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

info = struct();
info.platform = computer;
info.version = version;
info.date = char(datetime('now', 'Format', 'yyyy-MM-dd'));
[status, commit] = system(sprintf('git -C "%s" rev-parse --short HEAD', pktest.repoRoot()));
info.commit = '';
if status == 0
    info.commit = strtrim(commit);
end

for i = 1 : numel(scenarios)
    tmpDir = [tempname, '/'];
    mkdir(tmpDir);
    cleanupObj = onCleanup(@() rmdir(tmpDir, 's'));
    fprintf('Computing %s...\n', scenarios{i});
    expected = GoldenOutputTest.computeScenario(scenarios{i}, strrep(tmpDir, '\', '/'));
    save(fullfile(dataDir, [scenarios{i}, '.mat']), 'expected', 'info', '-v7');
    clear cleanupObj;
end

end
