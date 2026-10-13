function [results] = run_tests(varargin)
% run the PetaKit5D test suite. Tests that need something this machine
% lacks (GPU, demo dataset, or a mex not built for this platform) are
% reported as incomplete instead of failed.
%
%
% Parameters (as 'specifier'-value pairs):
%                 tag : char or cell array (default: {}). Only run tests with any of these tags: 'Unit', 'Integration', 'Golden', 'Slow', 'GPU', 'DemoData', 'KnownIssue'.
%          excludeTag : char or cell array (default: {}). Skip tests with any of these tags.
%  includeKnownIssues : true|false (default: false). Also run 'KnownIssue' tests, which check the correct behavior for bugs that are not fixed yet (see tests/README.md). They always run if 'KnownIssue' is in tag.
%                name : char or cell array (default: {}). Only run tests whose names match any of these patterns (* is a wildcard), e.g., '*Decon*'.
%              folder : char (default: ''). Only run tests under this subfolder of tests/, e.g., 'io'.
%           outputDir : char (default: ''). Write JUnit results (test-results.xml) and coverage (coverage.xml) to this folder.
%            coverage : true|false (default: false). Collect code coverage of the PetaKit5D sources (slow).
%           verbosity : 'Terse'|'Concise'|'Detailed'|'Verbose' (default: 'Concise'). Amount of output.
%        summaryTitle : char (default: ''). Heading for the GitHub Actions job summary (default: MATLAB release and platform).
%
% Examples:
%   run_tests                                     % all tests
%   run_tests('tag', 'Unit')                      % fast unit tests only
%   run_tests('name', '*Decon*')                  % deconvolution tests only
%   assertSuccess(run_tests('outputDir', 'test-results'))   % CI


ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('tag', {}, @(x) ischar(x) || iscell(x) || isstring(x));
ip.addParameter('excludeTag', {}, @(x) ischar(x) || iscell(x) || isstring(x));
ip.addParameter('includeKnownIssues', false, @islogical);
ip.addParameter('name', {}, @(x) ischar(x) || iscell(x) || isstring(x));
ip.addParameter('folder', '', @(x) ischar(x) || isstring(x));
ip.addParameter('outputDir', '', @(x) ischar(x) || isstring(x));
ip.addParameter('coverage', false, @islogical);
ip.addParameter('verbosity', 'Concise', @(x) any(strcmpi(x, {'Terse', 'Concise', 'Detailed', 'Verbose'})));
ip.addParameter('summaryTitle', '', @(x) ischar(x) || isstring(x));

ip.parse(varargin{:});

pr = ip.Results;
tags = cellstr(pr.tag);
excludeTags = cellstr(pr.excludeTag);
if ~pr.includeKnownIssues && ~any(strcmp(tags, 'KnownIssue'))
    excludeTags{end + 1} = 'KnownIssue';
end
namePatterns = cellstr(pr.name);
folder = char(pr.folder);
outputDir = char(pr.outputDir);
coverage = pr.coverage;
verbosity = matlab.unittest.Verbosity.(char(pr.verbosity));

testDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(testDir);
% for the pktest helper package
addpath(testDir);

searchDir = fullfile(testDir, folder);
suite = matlab.unittest.TestSuite.fromFolder(searchDir, 'IncludingSubfolders', true);

if ~isempty(tags)
    selector = matlab.unittest.selectors.HasTag(tags{1});
    for i = 2 : numel(tags)
        selector = selector | matlab.unittest.selectors.HasTag(tags{i});
    end
    suite = suite.selectIf(selector);
end
for i = 1 : numel(excludeTags)
    suite = suite.selectIf(~matlab.unittest.selectors.HasTag(excludeTags{i}));
end
if ~isempty(namePatterns)
    expr = strjoin(cellfun(@(x) ['^', regexptranslate('wildcard', x), '$'], namePatterns, 'unif', 0), '|');
    suite = suite(~cellfun(@isempty, regexp({suite.Name}, expr, 'once')));
end

fprintf('Running %d tests from %s on %s (MATLAB %s)\n', numel(suite), searchDir, computer, version);

runner = matlab.unittest.TestRunner.withTextOutput('OutputDetail', verbosity, 'LoggingLevel', verbosity);

if ~isempty(outputDir)
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    runner.addPlugin(matlab.unittest.plugins.XMLPlugin.producingJUnitFormat(fullfile(outputDir, 'test-results.xml')));
end
if coverage
    % only the PetaKit5D sources, not tests, demos, mcc parsers or third party code
    dir_info = dir(repoDir);
    srcDirs = setdiff({dir_info([dir_info.isdir]).name}, {'.', '..', '.git', '.github', 'tests', 'demos', 'mcc', 'third_parties'});
    srcDirs = cellfun(@(x) fullfile(repoDir, x), srcDirs, 'unif', 0);
    if isempty(outputDir)
        outputDir = pwd;
    end
    runner.addPlugin(matlab.unittest.plugins.CodeCoveragePlugin.forFolder(srcDirs, 'IncludingSubfolders', true, ...
        'Producing', matlab.unittest.plugins.codecoverage.CoberturaFormat(fullfile(outputDir, 'coverage.xml'))));
end

results = runner.run(suite);

summary = sprintf('%d passed, %d failed, %d incomplete (%.1f s)', nnz([results.Passed]), ...
    nnz([results.Failed]), nnz([results.Incomplete] & ~[results.Failed]), sum([results.Duration]));
fprintf('\n%s\n', summary);
failed = results([results.Failed]);
if ~isempty(failed)
    fprintf('Failed:\n');
    fprintf('  %s\n', failed.Name);
end

write_github_step_summary(results, summary, char(pr.summaryTitle));

end


function [] = write_github_step_summary(results, summary, title)
% on GitHub Actions, list the results on the job summary page

summaryFn = getenv('GITHUB_STEP_SUMMARY');
if isempty(summaryFn)
    return;
end
fid = fopen(summaryFn, 'a');
if fid < 0
    return;
end
cleanupObj = onCleanup(@() fclose(fid));

if isempty(title)
    title = sprintf('MATLAB %s on %s', version('-release'), computer);
end
fprintf(fid, '### %s\n\n%s\n\n', title, summary);
failed = results([results.Failed]);
if ~isempty(failed)
    fprintf(fid, '**Failed**\n\n');
    fprintf(fid, '- `%s`\n', failed.Name);
    fprintf(fid, '\n');
end
skipped = results([results.Incomplete] & ~[results.Failed]);
if ~isempty(skipped)
    fprintf(fid, '<details><summary>Skipped (%d)</summary>\n\n', numel(skipped));
    fprintf(fid, '- `%s`\n', skipped.Name);
    fprintf(fid, '\n</details>\n\n');
end

end
