classdef BatchFileExistTest < pktest.TestCase
    % batch_file_exist, which the job framework uses to check task inputs
    % and outputs. Its result must equal exist(..., 'file'). On Linux it
    % opens each folder before checking its files by name (up to 100 paths
    % in a folder) or looks the names up in a listing of the folder (more),
    % with Java or, without the JVM, dir(); on NFS this shows files other
    % nodes made at once. Elsewhere it checks the files one by one, with a
    % background pool for many files.

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
            % more than 100 files in a folder: one listing
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
            % on macOS and Windows
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

        function matchesExistOnEdgeCases(testCase)
            % few paths per folder (each folder is opened, then the paths
            % checked by name), all at once and one by one
            testCase.assumeLinux();
            queries = testCase.makeEdgeCases();
            expected = testCase.existFlags(queries);
            for useJava = [true, false]
                flags = logical(batch_file_exist(queries, [], true, useJava));
                testCase.verifyEqual(flags, expected, testCase.mismatchDiag(queries, expected, flags, useJava));
                flags = cellfun(@(x) logical(batch_file_exist({x}, [], true, useJava)), queries(:));
                testCase.verifyEqual(flags, expected, testCase.mismatchDiag(queries, expected, flags, useJava));
            end
        end

        function matchesExistOnEdgeCasesInListedFolders(testCase)
            % more than 100 paths in each folder: the folders are listed
            testCase.assumeLinux();
            [queries, folders] = testCase.makeEdgeCases();
            for d = 1 : numel(folders)
                testCase.makeFiles(folders{d}, 1 : 2 : 150);
                queries = [queries, arrayfun(@(x) sprintf('%sf_%03d.tif', folders{d}, x), 1 : 150, 'unif', 0)]; %#ok<AGROW>
            end
            expected = testCase.existFlags(queries);
            for useJava = [true, false]
                flags = logical(batch_file_exist(queries, [], true, useJava));
                testCase.verifyEqual(flags, expected, testCase.mismatchDiag(queries, expected, flags, useJava));
            end
        end

        function foldersAndTrailingSeparators(testCase)
            % zarr outputs are folders; a trailing separator doesn't matter
            testCase.assumeLinux();
            dataPath = testCase.createTempFolder();
            mkdir([dataPath, 'a.zarr']);
            fclose(fopen([dataPath, 'b.tif'], 'w'));
            queries = {[dataPath, 'a.zarr'], [dataPath, 'a.zarr/'], [dataPath, 'a.zarr//'], ...
                [dataPath, 'b.tif'], [dataPath, 'b.tif/'], [dataPath, 'c.zarr/'], dataPath, dataPath(1 : end - 1)};
            expected = [true; true; true; true; true; false; true; true];
            testCase.verifyEqual(testCase.existFlags(queries), expected, 'exist() itself changed.');
            for useJava = [true, false]
                testCase.verifyEqual(logical(batch_file_exist(queries, [], true, useJava)), expected);
                testCase.verifyEqual(logical(batch_file_exist(repmat(queries, 1, 20), [], true, useJava)), ...
                    repmat(expected, 20, 1));
            end
        end

        function missingParentFolders(testCase)
            dataPath = testCase.createTempFolder();
            fns = testCase.makeFiles(dataPath, 1 : 2 : 9);
            fclose(fopen([dataPath, 'file.tif'], 'w'));
            for n = [1, 3, 150]
                missingFns = [arrayfun(@(x) sprintf('%smissing/f_%03d.tif', dataPath, x), 1 : n, 'unif', 0), ...
                    arrayfun(@(x) sprintf('%smissing/deeper/f_%03d.tif', dataPath, x), 1 : n, 'unif', 0), ...
                    arrayfun(@(x) sprintf('%sfile.tif/f_%03d.tif', dataPath, x), 1 : n, 'unif', 0)];
                queries = [fns, missingFns, {[dataPath, 'f_002.tif']}];
                expected = [true(numel(fns), 1); false(numel(missingFns) + 1, 1)];
                for useJava = [true, false]
                    testCase.verifyEqual(logical(batch_file_exist(queries, [], true, useJava)), expected, ...
                        sprintf('%d paths per missing folder, useJava %d', n, useJava));
                end
            end
        end

        function duplicatePaths(testCase)
            dataPath = testCase.createTempFolder();
            fns = testCase.makeFiles(dataPath, 1 : 2 : 300);
            for n = [5, 300]
                queries = arrayfun(@(x) sprintf('%sf_%03d.tif', dataPath, x), [1 : n, n : -1 : 1, 1 : 3], 'unif', 0);
                for useJava = [true, false]
                    flags = batch_file_exist(queries, [], true, useJava);
                    testCase.verifyEqual(size(flags), [numel(queries), 1]);
                    testCase.verifyEqual(logical(flags), ismember(queries, fns)');
                end
            end
        end

        function mixedFoundAndMissingInManyFolders(testCase)
            % folders with 1 to 150 paths (opened or listed), some missing
            dataPath = testCase.createTempFolder();
            queries = {};
            for d = 1 : 12
                folder = sprintf('%sd%02d/', dataPath, d);
                if rem(d, 4) ~= 0
                    mkdir(folder);
                    testCase.makeFiles(folder, randperm(150, randi(150)));
                end
                n = randi(150);
                queries = [queries, arrayfun(@(x) sprintf('%sf_%03d.tif', folder, x), randperm(150, n), 'unif', 0)]; %#ok<AGROW>
            end
            queries = queries(randperm(numel(queries)));
            expected = testCase.existFlags(queries);
            testCase.assertTrue(any(expected) && ~all(expected));
            for useJava = [true, false]
                testCase.verifyEqual(logical(batch_file_exist(queries, [], true, useJava)), expected);
            end
        end

        function findsFilesMadeBetweenCalls(testCase)
            dataPath = testCase.createTempFolder();
            for n = [1, 150]
                folder = sprintf('%sn%d/', dataPath, n);
                mkdir(folder);
                queries = [arrayfun(@(x) sprintf('%sf_%03d.tif', folder, x), 1 : n, 'unif', 0), {[folder, 'out.zarr']}];
                for useJava = [true, false]
                    testCase.verifyEqual(logical(batch_file_exist(queries, [], true, useJava)), false(n + 1, 1));
                end
                testCase.makeFiles(folder, n);
                mkdir([folder, 'out.zarr']);
                expected = [false(n - 1, 1); true; true];
                for useJava = [true, false]
                    testCase.verifyEqual(logical(batch_file_exist(queries, [], true, useJava)), expected);
                end
            end
        end

        function largeFolder(testCase)
            % 20,000 files, 40,000 paths
            dataPath = testCase.createTempFolder();
            fns = testCase.makeFiles(dataPath, 1 : 2 : 40000, 'f_%05d.tif');
            queries = arrayfun(@(x) sprintf('%sf_%05d.tif', dataPath, x), 40000 : -1 : 1, 'unif', 0);
            expected = ismember(queries, fns)';
            for useJava = [true, false]
                testCase.verifyEqual(logical(batch_file_exist(queries, [], true, useJava)), expected);
            end
        end

        function unreadableFolder(testCase)
            % a folder that can be searched but not listed: exist() finds
            % its files, so batch_file_exist must too
            testCase.assumeTrue(isunix, 'Folder read permissions are POSIX only.');
            testCase.assumeTrue(usejava('jvm'), 'Removing read permission needs Java.');
            dataPath = testCase.createTempFolder();
            folder = [dataPath, 'noread/'];
            mkdir(folder);
            fns = testCase.makeFiles(folder, 1 : 2 : 150);
            % fileattrib can't change read permission
            java.io.File(folder).setReadable(false, false);
            testCase.addTeardown(@() java.io.File(folder).setReadable(true, false));
            for n = [3, 150]
                queries = arrayfun(@(x) sprintf('%sf_%03d.tif', folder, x), 1 : n, 'unif', 0);
                expected = ismember(queries, fns)';
                testCase.verifyEqual(testCase.existFlags(queries), expected, 'exist() itself changed.');
                for useJava = [true, false]
                    testCase.verifyEqual(logical(batch_file_exist(queries, [], true, useJava)), expected);
                end
            end
        end
    end

    methods (Test, TestTags = {'Unit', 'KnownIssue'})
        function existingOutputFileReturnsCachedResult(testCase)
            % KNOWN ISSUE (may be intended): batch_file_exist.m:43-45 returns early when the
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
        function [fns] = makeFiles(~, folder, inds, pattern)
            if nargin < 4
                pattern = 'f_%03d.tif';
            end
            fns = arrayfun(@(x) sprintf(['%s', pattern], folder, x), inds, 'unif', 0);
            for i = 1 : numel(fns)
                fclose(fopen(fns{i}, 'w'));
            end
        end

        function [queries, folders] = makeEdgeCases(testCase)
            % paths that exist() answers in unusual ways (Linux), and the
            % existing folders they are in. The current folder is the temp
            % folder.
            dataPath = testCase.createTempFolder();
            cd(dataPath);
            names = {'a.tif', 'Case.tif', ' a.tif', '#a.tif', '-a.tif', 'a b.tif', 'a%b.tif', 'a''b.tif', ...
                'a[1].tif', 'a(1).tif', 'a,b.tif', 'a&b.tif', '.hidden', 'a~1.tif', ['caf', char(233), '.tif'], ...
                'a*b.tif', 'a?b.tif', 'a\b.tif', 'a"b.tif', ['n', newline, 'l.tif'], 'trail .tif '};
            for i = 1 : numel(names)
                fclose(fopen([dataPath, names{i}], 'w'));
            end
            mkdir([dataPath, 'b.zarr']);
            fclose(fopen([dataPath, 'b.zarr/0.0.0'], 'w'));
            mkdir([dataPath, 'sub']);
            fclose(fopen([dataPath, 'sub/x.tif'], 'w'));
            folders = {dataPath, [dataPath, 'sub/']};

            queries = [strcat(dataPath, names), strcat(dataPath, {'A.TIF', 'CASE.TIF', 'case.tif', ...
                'b.zarr', 'b.zarr/', 'b.zarr//', 'a.tif/', 'b.zarr/0.0.0', 'b.zarr/.', 'b.zarr/..', ...
                'sub//x.tif', 'sub/../a.tif', 'sub/./x.tif', 'sub', 'missing.tif', 'missing/x.tif', ...
                'missing/deeper/x.tif', 'a.tif/x', 'a.tif', 'missing.tif'}), ...
                {'a.tif', 'sub/x.tif', 'nothing_here.tif', '/'}];
            mkdir([dataPath, 'st*r']);
            fclose(fopen([dataPath, 'st*r/s.tif'], 'w'));
            folders{end + 1} = [dataPath, 'st*r/'];
            % '~' is the home folder; exist() finds relative names on the MATLAB path
            queries = [queries, strcat(dataPath, {'st*r/s.tif', 'st*r/t.tif'}), {'~', '~/', 'batch_file_exist.m'}];
            if usejava('jvm')
                % symlinks to a file, a folder and nothing
                links = {'lnk_a', 'a.tif'; 'lnk_b', 'b.zarr'; 'lnk_dangle', 'nothing_here'; 'lnkdir', 'sub'};
                for i = 1 : size(links, 1)
                    java.nio.file.Files.createSymbolicLink(java.io.File([dataPath, links{i, 1}]).toPath(), ...
                        java.io.File(links{i, 2}).toPath(), javaArray('java.nio.file.attribute.FileAttribute', 0));
                end
                queries = [queries, strcat(dataPath, {'lnk_a', 'lnk_b', 'lnk_dangle', 'lnkdir/x.tif', 'lnkdir'})];
            end
        end

        function [] = assumeLinux(testCase)
            % macOS and Windows keep the per-file checks
            testCase.assumeTrue(isunix && ~ismac, 'The folder checks are used on Linux only.');
        end

        function [flags] = existFlags(~, queries)
            flags = cellfun(@(x) exist(x, 'file') > 0, queries(:));
        end

        function [diagStr] = mismatchDiag(~, queries, expected, actual, useJava)
            bad = find(logical(actual(:)) ~= expected);
            diagStr = sprintf('useJava %d, %d mismatches: %s', useJava, numel(bad), ...
                strjoin(cellfun(@(x) sprintf('[%s]', x), queries(bad), 'unif', 0), ' '));
        end
    end
end
