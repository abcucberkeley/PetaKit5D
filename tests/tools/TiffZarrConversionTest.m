classdef TiffZarrConversionTest < pktest.TestCase
    % tiff <-> zarr conversion (XR_tiffToZarr_wrapper, XR_zarrToTiff_wrapper,
    % tiffToZarr, zarrToTiff) and zarr resaving (XR_resave_zarr_wrapper).
    % XR_zarrToTiff_wrapper doesn't forward parseCluster to the job
    % framework, so on a machine with Slurm it submits jobs; the shared
    % pktest.NoSlurmFixture keeps it local.

    properties (Constant)
        ImSize = [40, 36, 12]
        Channels = {'CamA_ch0', 'CamB_ch0'}
    end

    methods (Test, TestTags = {'Integration'})
        function tiffToZarrWrapper(testCase)
            [dataPath, ims, names] = testCase.writeDataset('tif');
            XR_tiffToZarr_wrapper(dataPath, 'channelPatterns', testCase.Channels, 'blockSize', [16, 16, 4], ...
                'parseCluster', false, 'masterCompute', true);
            resultPath = [dataPath, 'zarr/'];
            for i = 1 : 2
                fn = [resultPath, names{i}, '.zarr'];
                testCase.verifyEqual(readzarr(fn), ims{i}, names{i});
                info = getZarrInfo(fn);
                testCase.verifyEqual(info.blockSize, [16, 16, 4]);
                testCase.verifyEqual(info.compressor.cname, 'zstd');
            end
            testCase.verifyFalse(exist([resultPath, names{3}, '.zarr'], 'dir') > 0, ...
                'A file that matches no channel pattern was converted.');
        end

        function tiffToZarrWrapperWithInputBbox(testCase)
            [dataPath, ims, names] = testCase.writeDataset('tif');
            b = [5, 7, 3, 30, 33, 10];
            XR_tiffToZarr_wrapper(dataPath, 'channelPatterns', testCase.Channels(1), 'inputBbox', b, ...
                'resultDirName', 'tiff_to_zarr', 'parseCluster', false, 'masterCompute', true);
            out = readzarr([dataPath, 'tiff_to_zarr/', names{1}, '.zarr']);
            testCase.verifyEqual(out, ims{1}(b(1) : b(4), b(2) : b(5), b(3) : b(6)));
        end

        function tiffToZarrWrapperWithoutReadingWholeTiff(testCase)
            % reads the tiff in blocks with blockedImage instead of at once
            [dataPath, ims, names] = testCase.writeDataset('tif');
            b = [1, 1, 1, testCase.ImSize];
            XR_tiffToZarr_wrapper(dataPath, 'channelPatterns', testCase.Channels(1), 'readWholeTiff', false, ...
                'inputBbox', b, 'blockSize', [16, 16, 4], 'parseCluster', false, 'masterCompute', true);
            testCase.verifyEqual(readzarr([dataPath, 'zarr/', names{1}, '.zarr']), ims{1});
        end

        function zarrToTiffWrapper(testCase)
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            XR_zarrToTiff_wrapper(dataPath, 'channelPatterns', testCase.Channels, ...
                'parseCluster', false, 'masterCompute', true);
            resultPath = [dataPath, 'tiffs/'];
            for i = 1 : 2
                testCase.verifyEqual(readtiff([resultPath, names{i}, '.tif']), ims{i}, names{i});
            end
            testCase.verifyFalse(exist([resultPath, names{3}, '.tif'], 'file') > 0, ...
                'A file that matches no channel pattern was converted.');
        end

        function zarrToTiffWrapperPassesJobSettings(testCase)
            % each task converts one file, so the memory per task is for one
            % file (it was the sum over all files), and parseCluster,
            % masterCompute and cpusPerTask were not passed on
            dataPath = testCase.writeDataset('zarr');
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            XR_zarrToTiff_wrapper(dataPath, 'channelPatterns', testCase.Channels, ...
                'parseCluster', false, 'masterCompute', false, 'cpusPerTask', 3);
            testCase.verifyEqual(spy.getOption(1, 'parseCluster'), false);
            testCase.verifyEqual(spy.getOption(1, 'masterCompute'), false);
            testCase.verifyEqual(spy.getOption(1, 'cpusPerTask'), 3);
            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), prod(testCase.ImSize) * 4 / 1024 ^ 3 * 2.5, ...
                'RelTol', 1e-12);
        end

        function resaveZarrChangesChunks(testCase)
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            XR_resave_zarr_wrapper(dataPath, 'blockSize', [8, 8, 4], 'channelPatterns', {'.zarr'}, ...
                'largeFile', false, 'parseCluster', false, 'masterCompute', true);
            for i = 1 : 3
                fn = [dataPath, 'resaved_zarr/', names{i}, '.zarr'];
                testCase.verifyEqual(readzarr(fn), ims{i}, names{i});
                testCase.verifyEqual(getZarrInfo(fn).blockSize, [8, 8, 4]);
            end
        end

        function resaveKeepsResultOfAnotherCopy(testCase)
            % two copies of a task (e.g. a retry next to a slow first try)
            % share the temporary result; the copy that finished second
            % deleted the result the first one had moved into place
            [dataPath, ~, names] = testCase.writeDataset('zarr');
            outFn = [dataPath, 'resaved.zarr'];
            % the other copy moves the result into place while this one runs
            testCase.applyFixture(pktest.FrameworkSpyFixture(@() movefile([outFn, '_test'], outFn)));
            XR_resaveSingleZarr([dataPath, names{1}, '.zarr'], outFn, 'blockSize', [8, 8, 4], ...
                'parseCluster', false, 'uuid', 'test');
            testCase.verifyTrue(isfolder(outFn), 'The result of the other copy was deleted.');
        end

        function resaveKeepsResultOfAnotherRun(testCase)
            % a copy from another run (its own uuid) that finds the result
            % published replaced it with its own, incomplete result
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            outFn = [dataPath, 'resaved.zarr'];
            % the other run publishes a complete result; this copy's blocks
            % never run (the framework reports them done once the result exists)
            testCase.applyFixture(pktest.FrameworkSpyFixture(@() writezarr(ims{1}, outFn, 'blockSize', [8, 8, 4])));
            XR_resaveSingleZarr([dataPath, names{1}, '.zarr'], outFn, 'blockSize', [8, 8, 4], ...
                'parseCluster', false, 'uuid', 'test');
            testCase.verifyEqual(readzarr(outFn), ims{1}, 'The published result was replaced.');
            testCase.verifyFalse(isfolder([outFn, '_test']), 'The temporary result was left behind.');
        end

        function resaveLargeZarrInBatches(testCase)
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            b = [5, 7, 3, 30, 33, 10];
            XR_resave_zarr_wrapper(dataPath, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 4], 'inputBbox', b, ...
                'channelPatterns', testCase.Channels(1), 'largeFile', true, 'parseCluster', false, 'masterCompute', true);
            out = readzarr([dataPath, 'resaved_zarr/', names{1}, '.zarr']);
            testCase.verifyEqual(out, ims{1}(b(1) : b(4), b(2) : b(5), b(3) : b(6)));
        end

        function resaveGroupsBatchesWithoutCluster(testCase)
            % without a cluster, neighboring batches are read and written
            % together (along y, then x and z, up to one chunk per thread),
            % so that cpp-zarr compresses their chunks in parallel (writes of
            % one chunk at a time kept a few cores busy); cluster tasks, and
            % too little or unknown memory, keep one batch per read and write
            [dataPath, ~, names] = testCase.writeDataset('zarr');
            src = [dataPath, names{1}, '.zarr'];

            % 12 threads, 3 x 5 x 3 batches of two chunks ([16, 8, 4]): up
            % to 6 batches, 3 along y and 2 along x
            testCase.applyFixture(pktest.FunctionStubFixture('maxNumCompThreads', 12));
            regions = testCase.resaveRegions(src, [dataPath, 'local2.zarr'], [16, 8, 4], false);
            testCase.verifySize(regions, [9, 6]);
            testCase.verifyEqual(unique(regions(:, 4 : 6) - regions(:, 1 : 3) + 1, 'rows'), [40, 4, 4; 40, 16, 4]);
            testCase.verifyRegionsTile(regions, [16, 8, 4]);

            % 6 threads, 5 x 5 x 3 batches of [8, 8, 4]: 15 groups of the 5
            % batches along y
            testCase.applyFixture(pktest.FunctionStubFixture('maxNumCompThreads', 6));
            regions = testCase.resaveRegions(src, [dataPath, 'local.zarr'], [8, 8, 4], false);
            testCase.verifySize(regions, [15, 6]);
            testCase.verifyEqual(unique(regions(:, 4 : 6) - regions(:, 1 : 3) + 1, 'rows'), [40, 4, 4; 40, 8, 4]);
            testCase.verifyRegionsTile(regions, [8, 8, 4]);

            regions = testCase.resaveRegions(src, [dataPath, 'cluster.zarr'], [8, 8, 4], true);
            testCase.verifySize(regions, [75, 6]);
            testCase.verifyRegionsTile(regions, [8, 8, 4]);

            for mem = [NaN, 1e-6]
                testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', mem));
                regions = testCase.resaveRegions(src, sprintf('%slow_memory_%g.zarr', dataPath, mem), [8, 8, 4], false);
                testCase.verifySize(regions, [75, 6], sprintf('memory %g GB', mem));
            end
        end

        function resaveGroupsWithinTheShareOfTheMemory(testCase)
            % the batches read and written together take at most 1/16 of this
            % process's share of the memory: the memory times its share of
            % the CPUs (a local process pool divides the threads among its
            % workers), so that processes that resave at once stay within it
            [dataPath, ~, names] = testCase.writeDataset('zarr');
            src = [dataPath, names{1}, '.zarr'];
            cores = getenv('SLURM_CPUS_ON_NODE');
            if isempty(cores)
                testCase.addTeardown(@unsetenv, 'SLURM_CPUS_ON_NODE');
            else
                testCase.addTeardown(@setenv, 'SLURM_CPUS_ON_NODE', cores);
            end
            setenv('SLURM_CPUS_ON_NODE', '24');
            % 64 kB of memory: 1/16 is 8 batches of [8, 8, 4] uint16
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 2^16 / 2^30));

            % all 24 CPUs: 8 batches, i.e. the 5 along y
            testCase.applyFixture(pktest.FunctionStubFixture('maxNumCompThreads', 24));
            regions = testCase.resaveRegions(src, [dataPath, 'all.zarr'], [8, 8, 4], false);
            testCase.verifySize(regions, [15, 6]);

            % 6 of 24 CPUs: a quarter of the memory, 2 batches
            testCase.applyFixture(pktest.FunctionStubFixture('maxNumCompThreads', 6));
            regions = testCase.resaveRegions(src, [dataPath, 'quarter.zarr'], [8, 8, 4], false);
            testCase.verifySize(regions, [45, 6]);
            testCase.verifyEqual(unique(regions(:, 4 : 6) - regions(:, 1 : 3) + 1, 'rows'), ...
                [8, 4, 4; 8, 8, 4; 16, 4, 4; 16, 8, 4]);
            testCase.verifyRegionsTile(regions, [8, 8, 4]);
        end

        function resaveGroupedBatchesDoNotTakeOtherFlags(testCase)
            % a run that groups batches, with the uuid of an interrupted run
            % that wrote one batch at a time (and finished its first task,
            % batches 1-5), must not take that task's flag for its own first
            % task (groups 1-5)
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            out = [dataPath, 'resumed.zarr'];
            createzarr([out, '_u'], dataSize=testCase.ImSize, blockSize=[8, 8, 4], dtype='uint16');
            writezarr(ims{1}(:, 1 : 8, 1 : 4), [out, '_u'], bbox=[1, 1, 1, 40, 8, 4], create=false);
            flagPath = [dataPath, 'zarr_flag/resumed_u/'];
            mkdir(flagPath);
            fclose(fopen([flagPath, 'batches_1_5.mat'], 'w'));

            testCase.applyFixture(pktest.FunctionStubFixture('maxNumCompThreads', 6));
            XR_resaveSingleZarr([dataPath, names{1}, '.zarr'], out, 'blockSize', [8, 8, 4], 'batchSize', [8, 8, 4], ...
                'parseCluster', false, 'uuid', 'u');
            testCase.verifyEqual(readzarr(out), ims{1});
        end

        function resaveGroupedBatchesWriteTheSameChunks(testCase)
            % batches read and written together give the same chunk files
            % as one batch at a time (cpp-zarr writes deterministic chunks),
            % also with a crop and batches of two chunks that don't divide
            % the image
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            b = [5, 7, 3, 30, 33, 10];
            threads = [1, 6];
            for k = 1 : 2
                testCase.applyFixture(pktest.FunctionStubFixture('maxNumCompThreads', threads(k)));
                XR_resave_zarr_wrapper(dataPath, 'resultDirName', sprintf('resaved_%d', k), 'batchSize', [16, 8, 4], ...
                    'blockSize', [8, 8, 4], 'inputBbox', b, 'channelPatterns', testCase.Channels(1), ...
                    'largeFile', true, 'parseCluster', false, 'masterCompute', true);
            end
            fns = {[dataPath, 'resaved_1/', names{1}, '.zarr'], [dataPath, 'resaved_2/', names{1}, '.zarr']};
            testCase.verifyEqual(readzarr(fns{2}), ims{1}(b(1) : b(4), b(2) : b(5), b(3) : b(6)));
            files = {dir(fns{1}), dir(fns{2})};
            files = cellfun(@(x) sort({x(~[x.isdir]).name}), files, 'unif', 0);
            testCase.verifyEqual(files{2}, files{1});
            % 4 x 4 x 2 chunks and the metadata
            testCase.verifyNumElements(files{1}, 33);
            for i = 1 : numel(files{1})
                bytes = cell(1, 2);
                for k = 1 : 2
                    fid = fopen([fns{k}, '/', files{1}{i}]);
                    bytes{k} = fread(fid, inf, '*uint8');
                    fclose(fid);
                end
                testCase.verifyEqual(bytes{2}, bytes{1}, files{1}{i});
            end
        end
    end

    methods (Test, TestTags = {'Integration', 'Slow'})
        function resaveOnLocalPoolGroupsWithinTheWorkersThreads(testCase)
            % files resaved on a local process pool (localWorkerNum in the
            % configuration file): each worker has its share of the threads,
            % and groups no more batches than that
            testCase.assumeTrue(license('test', 'Distrib_Computing_Toolbox') && ~isempty(which('gcp')), ...
                'Needs the Parallel Computing Toolbox.');
            % the framework starts a pool only if half of the memory holds the
            % workers (2 GB each besides the task's; GitHub's macOS runners
            % have 7 GB)
            testCase.assumeEqual(min(2, floor(getSystemMemory() / 2 / 2.1)), 2, ...
                'Too little memory for the framework to start a local pool of 2 workers.');
            if isempty(gcp('nocreate'))
                testCase.addTeardown(@() delete(gcp('nocreate')));
            end
            [dataPath, ims, names] = testCase.writeDataset('zarr');
            cfg = [dataPath, 'config.json'];
            fid = fopen(cfg, 'w');
            fprintf(fid, '{"localWorkerNum": 2}');
            fclose(fid);
            threadNum = max(1, floor(maxNumCompThreads / 2));

            out = evalc(['XR_resave_zarr_wrapper(dataPath, ''batchSize'', [8, 8, 4], ''blockSize'', [8, 8, 4], ', ...
                '''largeFile'', true, ''parseCluster'', false, ''masterCompute'', true, ''configFile'', cfg)']);

            testCase.verifySubstring(out, 'Run 3 tasks on the local pool');
            for i = 1 : 3
                testCase.verifyEqual(readzarr([dataPath, 'resaved_zarr/', names{i}, '.zarr']), ims{i}, names{i});
            end
            % the flags of grouped batches name the group size
            groups = regexp(out, 'batches_(\d+)_(\d+)_(\d+)_\d+_\d+', 'tokens');
            if threadNum > 1
                testCase.verifyNotEmpty(groups, 'The workers did not group batches.');
            end
            for i = 1 : numel(groups)
                batchNum = prod(str2double(groups{i}) ./ [8, 8, 4]);
                testCase.verifyLessThanOrEqual(batchNum, threadNum, sprintf('A worker grouped %d batches.', batchNum));
            end
        end
    end

    methods (Test, TestTags = {'Unit'})
        function singleFileRoundTrip(testCase)
            dataPath = testCase.createTempFolder();
            im = uint16(randi(60000, testCase.ImSize));
            writetiff(im, [dataPath, 'im.tif']);
            % the third input (in-memory frame) must be given, see the known
            % issue optionalFrameCanBeOmitted
            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'im.zarr'], [], 'blockSize', [16, 16, 4]);
            testCase.verifyEqual(readzarr([dataPath, 'im.zarr']), im);
            testCase.verifyEqual(getZarrInfo([dataPath, 'im.zarr']).blockSize, [16, 16, 4]);
            zarrToTiff([dataPath, 'im.zarr'], [dataPath, 'back.tif']);
            testCase.verifyEqual(readtiff([dataPath, 'back.tif']), im);
        end

        function tiffToZarrFromMemory(testCase)
            % the optional third input converts data that is already loaded
            dataPath = testCase.createTempFolder();
            im = single(randn(testCase.ImSize));
            tiffToZarr('', [dataPath, 'im.zarr'], im, 'blockSize', [16, 16, 4]);
            testCase.verifyEqual(readzarr([dataPath, 'im.zarr']), im);
        end

        function tiffToZarrProcessingOptions(testCase)
            % flip in z, then crop (tileOutBbox), resample, and the user
            % function, in that order; the result has the input's type
            dataPath = testCase.createTempFolder();
            im = uint16(randi(60000, testCase.ImSize));
            writetiff(im, [dataPath, 'im.tif']);
            b = [3, 4, 2, 38, 30, 11];
            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'flip.zarr'], [], 'blockSize', [16, 16, 4], ...
                'flipZstack', true, 'tileOutBbox', b, 'usrFcn', '@(x) single(x) / 3');
            expected = flip(im, 3);
            expected = uint16(single(expected(b(1) : b(4), b(2) : b(5), b(3) : b(6))) / 3);
            testCase.verifyEqual(readzarr([dataPath, 'flip.zarr']), expected);

            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'rs.zarr'], [], 'blockSize', [16, 16, 4], ...
                'resampleFactor', [1, 1, 2]);
            testCase.verifyEqual(readzarr([dataPath, 'rs.zarr']), imresize3(im, round(testCase.ImSize ./ [1, 1, 2])));
        end

        function tiffToZarrReadsRegions(testCase)
            % readWholeTiff false reads the tiff by regions; without
            % inputBbox it failed on an undefined variable (outBbox)
            dataPath = testCase.createTempFolder();
            im = uint16(randi(60000, testCase.ImSize));
            writetiff(im, [dataPath, 'im.tif']);
            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'whole.zarr'], [], 'blockSize', [16, 16, 4], ...
                'readWholeTiff', false);
            testCase.verifyEqual(readzarr([dataPath, 'whole.zarr']), im);
            b = [3, 4, 2, 38, 30, 11];
            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'crop.zarr'], [], 'blockSize', [16, 16, 4], ...
                'readWholeTiff', false, 'tileOutBbox', b);
            testCase.verifyEqual(readzarr([dataPath, 'crop.zarr']), im(b(1) : b(4), b(2) : b(5), b(3) : b(6)));
        end

        function tiffToZarrPartialFilesAndFrames(testCase)
            % partial files are joined along z; an in-memory frame gets the
            % same processing as a file; a 2D image becomes one z slice
            dataPath = testCase.createTempFolder();
            im = uint16(randi(60000, testCase.ImSize));
            writetiff(im(:, :, 1 : 5), [dataPath, 'part1.tif']);
            writetiff(im(:, :, 6 : end), [dataPath, 'part2.tif']);
            tiffToZarr({[dataPath, 'part1.tif'], [dataPath, 'part2.tif']}, [dataPath, 'parts.zarr'], [], ...
                'blockSize', [16, 16, 4]);
            testCase.verifyEqual(readzarr([dataPath, 'parts.zarr']), im);

            tiffToZarr('', [dataPath, 'frame.zarr'], im, 'blockSize', [16, 16, 4], 'flipZstack', true);
            testCase.verifyEqual(readzarr([dataPath, 'frame.zarr']), flip(im, 3));

            writetiff(im(:, :, 1), [dataPath, 'im2d.tif']);
            tiffToZarr([dataPath, 'im2d.tif'], [dataPath, 'im2d.zarr'], [], 'blockSize', [16, 16, 4]);
            testCase.verifyEqual(readzarr([dataPath, 'im2d.zarr']), im(:, :, 1));
            testCase.verifyEqual(getZarrInfo([dataPath, 'im2d.zarr']).size, [testCase.ImSize(1 : 2), 1]);
        end

        function tiffToZarrWritesPrivateResultInPlace(testCase)
            % with a full uuid made by the call (not on Windows, where
            % get_uuid shortens it) no other process writes the temporary
            % result, so its chunks are written without cpp-zarr's temporary
            % files ('uuid' false); a given uuid can be shared by copies of a
            % task, which keep them
            dataPath = testCase.createTempFolder();
            writetiff(uint16(randi(60000, testCase.ImSize)), [dataPath, 'im.tif']);
            givenUuid = '01234567-89ab-cdef-0123-456789abcdef';
            spy = testCase.applyFixture(pktest.FunctionSpyFixture('parallelWriteZarr'));
            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'whole.zarr'], [], 'blockSize', [16, 16, 4]);
            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'regions.zarr'], [], 'blockSize', [16, 16, 4], ...
                'readWholeTiff', false);
            tiffToZarr([dataPath, 'im.tif'], [dataPath, 'given.zarr'], [], 'blockSize', [16, 16, 4], 'uuid', givenUuid);
            calls = spy.getCalls();
            % one write of the whole image, one per z block of 4 slices, one with the given uuid
            testCase.assertNumElements(calls, 5);
            for i = 1 : 4
                if numel(get_uuid()) >= 32
                    testCase.verifyEqual(spy.getOption(i, 'uuid'), false, sprintf('call %d (%s)', i, calls{i}{1}));
                else
                    testCase.verifyEmpty(spy.getOption(i, 'uuid'), sprintf('call %d (%s)', i, calls{i}{1}));
                end
            end
            testCase.verifyEqual(calls{5}{1}, [dataPath, 'given.zarr_', givenUuid]);
            testCase.verifyEmpty(spy.getOption(5, 'uuid'), 'A temporary result named with a given uuid was written in place.');
        end
    end

    methods (Access = private)
        function [regions] = resaveRegions(testCase, src, out, batchSize, parseCluster)
            % the output regions of the tasks XR_resaveSingleZarr hands to the
            % framework (with chunks of [8, 8, 4]), one row per read and write
            frameworkSpy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            XR_resaveSingleZarr(src, out, 'blockSize', [8, 8, 4], 'batchSize', batchSize, ...
                'parseCluster', parseCluster);
            calls = frameworkSpy.getCalls();
            blockSpy = testCase.applyFixture(pktest.FunctionSpyFixture('resaveZarrBlock'));
            for i = 1 : numel(calls{1}.funcStrs)
                eval(calls{1}.funcStrs{i});
            end
            tasks = blockSpy.getCalls();
            regions = cell2mat(cellfun(@(x) x{6}, tasks(:), 'unif', 0));
        end

        function [] = verifyRegionsTile(testCase, regions, batchSize)
            % the regions are whole batches that cover the image once
            counts = zeros(testCase.ImSize);
            for i = 1 : size(regions, 1)
                r = regions(i, :);
                counts(r(1) : r(4), r(2) : r(5), r(3) : r(6)) = counts(r(1) : r(4), r(2) : r(5), r(3) : r(6)) + 1;
            end
            testCase.verifyTrue(all(counts(:) == 1), 'The regions do not cover the image once.');
            testCase.verifyTrue(all(mod(regions(:, 1 : 3) - 1, batchSize) == 0, 'all'), 'A region starts inside a batch.');
            testCase.verifyTrue(all(mod(regions(:, 4 : 6), batchSize) == 0 | regions(:, 4 : 6) == testCase.ImSize, 'all'), ...
                'A region ends inside a batch.');
        end

        function [dataPath, ims, names] = writeDataset(testCase, format)
            dataPath = testCase.createTempFolder();
            names = {'scan_CamA_ch0_t0000', 'scan_CamB_ch0_t0000', 'scan_CamA_ch1_t0000'};
            ims = cell(1, 3);
            for i = 1 : 3
                ims{i} = uint16(randi(60000, testCase.ImSize));
                if strcmp(format, 'zarr')
                    writezarr(ims{i}, [dataPath, names{i}, '.zarr'], 'blockSize', [16, 16, 4]);
                else
                    writetiff(ims{i}, [dataPath, names{i}, '.tif']);
                end
            end
        end
    end
end
