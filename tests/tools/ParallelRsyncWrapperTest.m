classdef ParallelRsyncWrapperTest < pktest.TestCase
    % XR_parallel_rsync_wrapper, which copies a folder in batches of rsync
    % calls run by the job framework (a spy here), then with one rsync of
    % the whole folder

    properties (TestParameter)
        % the folder of another setup.m first on the path
        decoyFolder = struct('otherName', 'apps', 'nameWithPetaKit5D', 'PetaKit5D_app')
    end

    methods (Test, TestTags = {'Unit'})
        function batchScriptIsFoundWhicheverSetupIsFirst(testCase, decoyFolder)
            % the wrapper found the code by the first setup.m on the path and
            % a folder named PetaKit5D: with another setup.m first (e.g., of
            % an app that bundles PetaKit5D), or with the code in a folder of
            % another name, it changed to a folder that doesn't exist, or
            % gave the jobs a batch script that doesn't exist
            testCase.assumeFalse(ispc, 'rsync is not available on Windows.');
            src = [testCase.createTempFolder(), 'src'];
            mkdir(src);
            for i = 1 : 4
                fclose(fopen(sprintf('%s/f%d.txt', src, i), 'w'));
            end
            dest = testCase.createTempFolder();
            decoy = [testCase.createTempFolder(), decoyFolder];
            mkdir(decoy);
            fid = fopen([decoy, '/setup.m'], 'w');
            fprintf(fid, '%s\n', 'function [] = setup(varargin)', '% another setup.m (test decoy)', 'end');
            fclose(fid);
            % the wrapper changes to the repository folder
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(pwd));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(decoy));
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            evalc('XR_parallel_rsync_wrapper(src, dest, ''parseCluster'', false)');

            calls = spy.getCalls();
            testCase.assertNotEmpty(calls, 'The job framework was not called.');
            script = regexp(calls{1}.funcStrs{1}, '^bash (\S+/rsync_batch_files\.sh) ', 'tokens', 'once');
            testCase.assertNotEmpty(script, 'No batch script in the task.');
            testCase.verifyTrue(isfile(script{1}), sprintf('The batch script %s does not exist.', script{1}));
            testCase.verifyEqual(numel(dir([dest, 'src/*.txt'])), 4, 'The final rsync did not copy every file.');
        end

        function eachBatchIsCheckedByItsLastFile(testCase)
            % the job framework counts a batch as done once its checked file
            % exists, and cancels the remaining jobs when all are; the second
            % to last file of each batch was checked, so the last one could
            % be cut short, and a last batch of fewer files had none
            testCase.assumeFalse(ispc, 'rsync is not available on Windows.');
            src = [testCase.createTempFolder(), 'src'];
            mkdir(src);
            % 1252 files: batches of 3, the last with 1 file
            for i = 1 : 1252
                fclose(fopen(sprintf('%s/f%04d.txt', src, i), 'w'));
            end
            dest = testCase.createTempFolder();
            % the wrapper changes to the repository folder
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(pwd));
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            evalc('XR_parallel_rsync_wrapper(src, dest, ''parseCluster'', false)');

            calls = spy.getCalls();
            testCase.assertNotEmpty(calls, 'The job framework was not called.');
            c = calls{1};
            testCase.assertEqual(numel(c.inputFullpaths), numel(c.funcStrs), 'Each batch needs a checked file.');
            lastIns = cell(numel(c.funcStrs), 1);
            for b = 1 : numel(c.funcStrs)
                ins = regexp(c.funcStrs{b}, '''''([^'']*)''''', 'tokens', 'once');
                ins = strsplit(ins{1}, ',');
                lastIns{b} = ins{end};
            end
            testCase.verifyEqual(c.inputFullpaths(:), lastIns, 'Batches are not checked by their last file.');
            testCase.verifyEqual(numel(dir([dest, 'src/*.txt'])), 1252, 'The final rsync did not copy every file.');
        end
    end
end
