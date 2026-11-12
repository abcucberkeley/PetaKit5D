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
            words = ParallelRsyncWrapperTest.shellWords(calls{1}.funcStrs{1});
            testCase.assertTrue(numel(words) > 1 && strcmp(words{1}, 'bash') && endsWith(words{2}, '/rsync_batch_files.sh'), ...
                'No batch script in the task.');
            testCase.verifyTrue(isfile(words{2}), sprintf('The batch script %s does not exist.', words{2}));
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
                [ins, outs] = ParallelRsyncWrapperTest.batchPairs(c.funcStrs{b});
                lastIns{b} = ins{end};
                % each file goes to its folder in the copy
                testCase.verifyEqual(outs, fileparts(strrep(ins, [src, '/'], [dest, 'src/'])), ...
                    'A file goes to another folder.');
            end
            testCase.verifyEqual(c.inputFullpaths(:), lastIns, 'Batches are not checked by their last file.');
            testCase.verifyEqual(numel(dir([dest, 'src/*.txt'])), 1252, 'The final rsync did not copy every file.');
        end

        function namesWithSpacesCommasAndQuotesAreCopied(testCase)
            % paths had their spaces escaped as '\ ' (also in the checked
            % files, so files with spaces were never found, and a batch
            % checked by one was skipped), the batch script split its lists
            % at commas (a name with a comma), and the final rsync had its
            % paths unquoted (a folder with a space)
            testCase.assumeFalse(ispc, 'rsync is not available on Windows.');
            src = [testCase.createTempFolder(), 'my src'];
            mkdir([src, '/sub dir']);
            special = {'a b.txt', 'c,d.txt', 'it''s.txt', 'sub dir/e,f g.txt'};
            % 750 files: batches of 2
            names = [special, arrayfun(@(x) sprintf('f%04d.txt', x), 1 : 746, 'unif', 0)];
            for i = 1 : numel(names)
                fid = fopen([src, '/', names{i}], 'w');
                fprintf(fid, '%s', names{i});
                fclose(fid);
            end
            dest = [testCase.createTempFolder(), 'out dir'];
            copy = [dest, '/my src/'];
            % the wrapper changes to the repository folder
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(pwd));
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            evalc('XR_parallel_rsync_wrapper(src, dest, ''parseCluster'', false)');

            calls = spy.getCalls();
            testCase.assertNotEmpty(calls, 'The job framework was not called.');
            c = calls{1};
            testCase.verifyTrue(all(isfile(c.inputFullpaths)), 'A batch is checked by a file that does not exist.');
            testCase.verifyEqual(sort(strrep(c.outputFullpaths(:), copy, '')), ...
                sort(strrep(c.inputFullpaths(:), [src, '/'], '')), 'A checked file is not the copy of its file.');
            testCase.verifyCopies(copy, names, 'The final rsync');

            % the batches with the special names, run as the framework runs
            % them: in this MATLAB, and in the script of a Slurm job
            % (sbatch --wrap="...")
            bs = find(contains(c.funcStrs, {' b.txt', ',d.txt', 's.txt', ' g.txt'}));
            testCase.assertNotEmpty(bs);
            files = {};
            for b = bs(:)'
                ins = ParallelRsyncWrapperTest.batchPairs(c.funcStrs{b});
                testCase.verifyNumElements(ins, 2);
                files = [files, strrep(ins(:)', [src, '/'], '')];
            end
            testCase.verifyEmpty(setdiff(special, files), 'The batches leave out a file.');
            % (in an empty folder, which must stay empty)
            cwd = testCase.createTempFolder();
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(cwd));
            for mode = {'local', 'job'}
                rmdir(copy, 's');
                for b = bs(:)'
                    cmd = c.funcStrs{b};
                    if strcmp(mode{1}, 'job')
                        cmd = sprintf('bash -c "%s"', cmd);
                    end
                    [status, out] = system(cmd);
                    testCase.verifyEqual(status, 0, out);
                    % (no failed command in the batch, which ends with an echo)
                    testCase.verifyEmpty(regexpi(out, 'error|cannot', 'once'), out);
                end
                testCase.verifyCopies(copy, files, sprintf('The batches (%s)', mode{1}));
                left = dir(cwd);
                testCase.verifyEmpty(setdiff({left.name}, {'.', '..'}), 'The batches made files elsewhere.');
            end
        end

        function clusterFileCheckFindsTheFiles(testCase)
            % batch_file_exist_cluster (the check of 100000 files or more)
            % used undefined MatlabLaunchStr, mccMode and configFile, so it
            % always failed
            root = testCase.createTempFolder();
            names = {'a b.txt', 'c,d.txt', 'it''s.txt', 'plain.txt', 'missing.txt'};
            for i = 1 : 4
                fclose(fopen([root, names{i}], 'w'));
            end
            mkdir([root, 'cwd']);
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture([root, 'cwd']));
            paths = strcat(root, names(:));
            % (no Slurm here: the framework runs the checks in this MATLAB)
            tf = [];
            evalc('tf = batch_file_exist_cluster(paths, [root, ''check1'']);');
            testCase.verifyEqual(tf, [true; true; true; true; false]);
            tf = [];
            evalc('tf = batch_file_exist_cluster(paths, [root, ''check2''], ''parseCluster'', false);');
            testCase.verifyEqual(tf, [true; true; true; true; false]);
        end

        function clusterFileCheckOfABatchOfOnePath(testCase)
            % the task of a batch of one path (e.g., the last batch of
            % 50001 files) wrote no result, so its check never finished
            root = testCase.createTempFolder();
            fclose(fopen([root, 'a.txt'], 'w'));
            mkdir([root, 'cwd']);
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture([root, 'cwd']));
            % (no Slurm here: the framework runs the check in this MATLAB)
            tf = [];
            evalc('tf = batch_file_exist_cluster({[root, ''a.txt'']}, [root, ''check'']);');
            testCase.verifyEqual(tf, true);
            testCase.verifyNumElements(dir([root, 'check/file_exist_flags/file_exist_1_1.mat']), 1, ...
                'The task wrote no result.');
        end

        function failedClusterFileChecksAreCheckedHere(testCase)
            % a batch whose check failed in its jobs (here every batch: a
            % framework spy that runs nothing) made the result a cell of
            % empties, on which the rsync wrapper failed; such batches are
            % checked in this MATLAB, with a warning
            root = testCase.createTempFolder();
            names = {'a b.txt', 'it''s.txt', 'missing.txt'};
            for i = 1 : 2
                fclose(fopen([root, names{i}], 'w'));
            end
            paths = strcat(root, names(:));
            testCase.applyFixture(pktest.FrameworkSpyFixture([], false));
            tf = testCase.verifyWarning(@() batch_file_exist_cluster(paths, [root, 'check']), ...
                'batch_file_exist_cluster:checkFailed');
            testCase.verifyEqual(tf, [true; true; false]);
        end
    end

    methods (Test, TestTags = {'Unit'})
        function largeCopiesPassTheirSettingsToTheClusterCheck(testCase)
            % copies of 100000 files or more are checked by jobs: the
            % wrapper's settings go to the check (a stub here, given 100000
            % file names that report all files as copied)
            testCase.assumeFalse(ispc, 'rsync is not available on Windows.');
            src = [testCase.createTempFolder(), 'src'];
            mkdir(src);
            dest = testCase.createTempFolder();
            paths = arrayfun(@(x) sprintf('%s/f%06d.txt', src, x), (1 : 100000)', 'unif', 0);
            testCase.applyFixture(pktest.FunctionStubFixture('dir_recursive', paths));
            stub = testCase.createTempFolder();
            fid = fopen([stub, 'batch_file_exist_cluster.m'], 'w');
            fprintf(fid, '%s\n', 'function [tf] = batch_file_exist_cluster(fileFullpaths, tmpDir, varargin)', ...
                'save([fileparts(mfilename(''fullpath'')), ''/call.mat''], ''fileFullpaths'', ''tmpDir'', ''varargin'');', ...
                'tf = true(numel(fileFullpaths), 1);', 'end');
            fclose(fid);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(stub));
            testCase.addTeardown(@() clear('batch_file_exist_cluster'));
            % the wrapper changes to the repository folder
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(pwd));
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            evalc('XR_parallel_rsync_wrapper(src, dest, ''parseCluster'', false, ''configFile'', ''/no/config.json'')');

            testCase.assertTrue(isfile([stub, 'call.mat']), 'The cluster check was not called.');
            a = load([stub, 'call.mat']);
            testCase.verifyNumElements(a.fileFullpaths, 100000);
            opts = cell2struct(a.varargin(2 : 2 : end), cellfun(@char, a.varargin(1 : 2 : end), 'unif', 0), 2);
            testCase.verifyEqual(opts, struct('parseCluster', false, 'mccMode', false, 'configFile', '/no/config.json'));
            testCase.verifyEmpty(spy.getCalls(), 'The files were copied again.');
        end
    end

    methods (Access = private)
        function [] = verifyCopies(testCase, copy, names, what)
            % each file in copy has its name as its content
            for i = 1 : numel(names)
                fn = [copy, names{i}];
                testCase.verifyTrue(isfile(fn), sprintf('%s did not copy %s.', what, names{i}));
                if isfile(fn)
                    testCase.verifyEqual(fileread(fn), names{i});
                end
            end
        end
    end

    methods (Static, Access = private)
        function [words] = shellWords(str)
            % the words of a bash command of plain and single-quoted words
            words = regexp(str, '(''[^'']*''(\\''''[^'']*'')*|[^ '']+)', 'match');
            for i = 1 : numel(words)
                if startsWith(words{i}, '''')
                    words{i} = strrep(words{i}(2 : end - 1), '''\''''', '''');
                end
            end
        end

        function [ins, outs] = batchPairs(funcStr)
            % the input files and output folders of a batch task
            words = ParallelRsyncWrapperTest.shellWords(funcStr);
            ins = words(3 : 2 : end)';
            outs = words(4 : 2 : end)';
        end
    end
end
