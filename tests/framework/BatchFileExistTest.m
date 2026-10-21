classdef BatchFileExistTest < pktest.TestCase
    % batch_file_exist, which the job framework uses to check task outputs.
    % It has three strategies: per-file checks, one folder listing for many
    % files in few folders, and a background pool for many files.

    methods (Test, TestTags = {'Unit'})
        function emptyAndSingleInputs(testCase)
            dataPath = testCase.createTempFolder();
            testCase.verifyEmpty(batch_file_exist({}));
            testCase.verifyTrue(logical(batch_file_exist({dataPath(1 : end - 1)})), 'A folder counts as existing.');
            testCase.verifyFalse(logical(batch_file_exist({[dataPath, 'missing.tif']})));
        end

        function checksFilesAndFoldersIndividually(testCase)
            dataPath = testCase.createTempFolder();
            fns = testCase.makeFiles(dataPath, 1 : 2 : 9);
            mkdir([dataPath, 'a.zarr']);
            queries = [arrayfun(@(x) sprintf('%sf_%03d.tif', dataPath, x), 1 : 9, 'unif', 0), {[dataPath, 'a.zarr']}];
            expected = [ismember(queries(1 : 9), fns), true]';
            testCase.verifyEqual(logical(batch_file_exist(queries)), expected);
        end

        function listsTheFolderForManyFiles(testCase)
            % more than 100 files in fewer than 3 folders: one dir() call
            dataPath = testCase.createTempFolder();
            fns = testCase.makeFiles(dataPath, 1 : 2 : 150);
            queries = arrayfun(@(x) sprintf('%sf_%03d.tif', dataPath, x), 1 : 150, 'unif', 0);
            testCase.verifyEqual(logical(batch_file_exist(queries)), ismember(queries, fns)');
        end

        function missingFolderWithManyFiles(testCase)
            % eb1e3b8 regression: listing a folder that doesn't exist errored
            dataPath = testCase.createTempFolder();
            queries = arrayfun(@(x) sprintf('%smissing/f_%03d.tif', dataPath, x), 1 : 101, 'unif', 0);
            testCase.verifyEqual(logical(batch_file_exist(queries)), false(101, 1));
        end

        function backgroundPoolForManyFilesInManyFolders(testCase)
            % more than 100 files in 3 or more folders with useParpool: parfeval
            dataPath = testCase.createTempFolder();
            queries = {};
            expected = false(0, 1);
            for d = 1 : 3
                folder = sprintf('%sd%d/', dataPath, d);
                mkdir(folder);
                fns = testCase.makeFiles(folder, 1 : 3 : 40);
                q = arrayfun(@(x) sprintf('%sf_%03d.tif', folder, x), 1 : 40, 'unif', 0);
                queries = [queries, q]; %#ok<AGROW>
                expected = [expected; ismember(q, fns)']; %#ok<AGROW>
            end
            testCase.verifyEqual(logical(batch_file_exist(queries, [], true)), expected);
        end

        function savesResultToOutputFile(testCase)
            dataPath = testCase.createTempFolder();
            fns = testCase.makeFiles(dataPath, 1 : 3);
            queries = [fns, {[dataPath, 'missing.tif']}];
            outFn = [dataPath, 'exist.mat'];
            flags = batch_file_exist(queries, outFn);
            testCase.verifyEqual(logical(flags), [true; true; true; false]);
            a = load(outFn);
            testCase.verifyEqual(logical(a.file_exist_mat), logical(flags));
            testCase.verifyEqual(a.fileFullpaths, queries);
        end
    end

    methods (Test, TestTags = {'Unit', 'KnownIssue'})
        function existingOutputFileReturnsCachedResult(testCase)
            % KNOWN ISSUE (may be intended): batch_file_exist.m:28-30 returns early when the
            % result file outFullpath already exists, before file_exist_mat
            % is assigned, so a second call errors with "Output argument
            % file_exist_mat not assigned". Fix: load file_exist_mat from
            % outFullpath before returning.
            dataPath = testCase.createTempFolder();
            fns = testCase.makeFiles(dataPath, 1 : 3);
            queries = [fns, {[dataPath, 'missing.tif']}];
            outFn = [dataPath, 'exist.mat'];
            first = batch_file_exist(queries, outFn);
            second = batch_file_exist(queries, outFn);
            testCase.verifyEqual(logical(second), logical(first));
        end
    end

    methods (Access = private)
        function [fns] = makeFiles(~, folder, inds)
            fns = arrayfun(@(x) sprintf('%sf_%03d.tif', folder, x), inds, 'unif', 0);
            for i = 1 : numel(fns)
                fclose(fopen(fns{i}, 'w'));
            end
        end
    end
end
