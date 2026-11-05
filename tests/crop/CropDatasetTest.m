classdef CropDatasetTest < pktest.TestCase
    % XR_crop_dataset / XR_crop_frame / XR_crop_zarr on small synthetic
    % datasets (two channels plus a channel that must be ignored)

    properties (Constant)
        ImSize = [40, 36, 12]
        Bbox = [5, 7, 3, 30, 33, 10]
        Channels = {'CamA_ch0', 'CamB_ch0'}
    end

    properties
        DataPath
        Images
        Names
    end

    methods (Test, TestTags = {'Integration'})
        function tiffToTiff(testCase)
            testCase.writeDataset('tif');
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels, ...
                'parseCluster', false, 'masterCompute', true);

            resultPath = [testCase.DataPath, 'Cropped/'];
            for i = 1 : 2
                out = readtiff([resultPath, testCase.Names{i}, '.tif']);
                testCase.verifyEqual(out, testCase.crop(testCase.Images{i}, testCase.Bbox), testCase.Names{i});
            end
            testCase.verifyFalse(exist([resultPath, testCase.Names{3}, '.tif'], 'file') > 0, ...
                'A file that matches no channel pattern was cropped.');
            testCase.verifyTrue(exist([resultPath, 'parameters.mat'], 'file') > 0);
        end

        function zarrToZarr(testCase)
            testCase.writeDataset('zarr');
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels, ...
                'zarrFile', true, 'saveZarr', true, 'parseCluster', false, 'masterCompute', true);
            for i = 1 : 2
                out = readzarr([testCase.DataPath, 'Cropped/', testCase.Names{i}, '.zarr']);
                testCase.verifyEqual(out, testCase.crop(testCase.Images{i}, testCase.Bbox), testCase.Names{i});
            end
        end

        function zarrToTiff(testCase)
            testCase.writeDataset('zarr');
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels, ...
                'zarrFile', true, 'saveZarr', false, 'parseCluster', false, 'masterCompute', true);
            for i = 1 : 2
                out = readtiff([testCase.DataPath, 'Cropped/', testCase.Names{i}, '.tif']);
                testCase.verifyEqual(out, testCase.crop(testCase.Images{i}, testCase.Bbox), testCase.Names{i});
            end
        end

        function largeZarrMatchesInMemoryCrop(testCase)
            % small batches so the crop is assembled from several blocks
            testCase.writeDataset('zarr');
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels, ...
                'zarrFile', true, 'largeFile', true, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 4], ...
                'parseCluster', false, 'masterCompute', true);
            for i = 1 : 2
                out = readzarr([testCase.DataPath, 'Cropped/', testCase.Names{i}, '.zarr']);
                testCase.verifyEqual(out, testCase.crop(testCase.Images{i}, testCase.Bbox), testCase.Names{i});
            end
        end

        function largeZarrMakesMIPsBeforePublishing(testCase)
            % the job framework counts a crop as done, and may cancel its
            % job, once the result exists; the large file path moved the
            % result into place and then made its MIPs
            testCase.writeDataset('zarr');
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels, ...
                'zarrFile', true, 'largeFile', true, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 4], ...
                'parseCluster', false, 'masterCompute', true);
            for i = 1 : 2
                cropFn = [testCase.DataPath, 'Cropped/', testCase.Names{i}, '.zarr'];
                testCase.assertTrue(moves.wasMoved(cropFn), sprintf('%s was not moved into place.', testCase.Names{i}));
                mips = strcat('MIPs/', testCase.Names{i}, {'_MIP_y.tif'; '_MIP_x.tif'; '_MIP_z.tif'});
                testCase.verifyEqual(ismember(mips, moves.filesWhenMoved(cropFn)), true(3, 1), ...
                    sprintf('MIPs of %s were made after it was moved into place.', testCase.Names{i}));
                mip = readtiff([testCase.DataPath, 'Cropped/MIPs/', testCase.Names{i}, '_MIP_z.tif']);
                testCase.verifyEqual(mip, max(testCase.crop(testCase.Images{i}, testCase.Bbox), [], 3), testCase.Names{i});
            end
        end

        function largeZarrMIPsAreMergedFromTheBlocks(testCase)
            % the crop tasks save the MIPs of their blocks in their flags, and
            % XR_MIP_zarr merges them into the MIP slabs instead of reading
            % the result again in MIP_block tasks
            testCase.writeDataset('zarr');
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            % MIP_block tasks would read the result; the spy errors at once, as the
            % framework would retry a task without output for minutes
            mipTasks = testCase.applyFixture(pktest.FunctionSpyFixture('MIP_block', [], 0));
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels(1), ...
                'zarrFile', true, 'largeFile', true, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 4], ...
                'parseCluster', false, 'masterCompute', true);
            testCase.verifyEmpty(mipTasks.getCalls(), 'The MIPs were computed from the result.');
            calls = moves.getCalls();
            testCase.verifyTrue(any(contains({calls.dst}, '/zarr_flag/') & endsWith({calls.dst}, '.mat')), ...
                'The block flags were not moved into place.');
            pktest.verifyLargeFileMIPs(testCase, [testCase.DataPath, 'Cropped/', testCase.Names{1}, '.zarr']);
        end

        function largeZarrBlocksGetTheMIPOptionOnlyForMIPs(testCase)
            % compiled masters built before the block tasks' mipBatchSize
            % option still run the tasks of crops without MIPs
            testCase.writeDataset('zarr');
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture());
            zarrFn = [testCase.DataPath, testCase.Names{1}, '.zarr'];
            for saveMIP = [false, true]
                evalc(['XR_crop_zarr(zarrFn, sprintf(''%scrop_%d.zarr'', testCase.DataPath, saveMIP), testCase.Bbox, ', ...
                    '''batchSize'', [16, 16, 4], ''blockSize'', [8, 8, 4], ''saveMIP'', saveMIP, ''parseCluster'', false)']);
            end
            calls = spy.getCalls();
            testCase.verifyFalse(any(contains(calls{1}.funcStrs, 'mipBatchSize')), 'Tasks without MIPs got mipBatchSize.');
            testCase.verifyTrue(all(contains(calls{2}.funcStrs, '''mipBatchSize'',[')), 'Tasks with MIPs got no mipBatchSize.');
        end

        function largeZarrKeepsAResultAnotherWriterPublished(testCase)
            % two job sets on the same data: the one that finished last
            % deleted the crop the other had published and moved its own
            % into place; now it keeps the other's and removes its own
            % (incomplete here: the framework is a spy that reports the
            % blocks done while the other crop comes)
            testCase.writeDataset('zarr');
            cropFn = [testCase.DataPath, 'crop.zarr'];
            other = uint16(7 * ones(5, 5, 5));
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture(@() pktest.atCall([testCase.DataPath, 'calls.txt'], 1, ...
                @() writezarr(other, cropFn))));
            evalc(['XR_crop_zarr([testCase.DataPath, testCase.Names{1}, ''.zarr''], cropFn, testCase.Bbox, ', ...
                '''batchSize'', [16, 16, 4], ''blockSize'', [8, 8, 4], ''parseCluster'', false)']);
            testCase.verifyEqual(readzarr(cropFn), other, 'The other writer''s crop was replaced.');
            testCase.verifyEmpty(dir([cropFn, '_*']), 'The temporary crop was not removed.');
            testCase.verifyEmpty(dir([testCase.DataPath, 'zarr_flag/crop_*']), 'The flags were not removed.');
            % (no MIPs of this one's incomplete crop, named after the crop)
            testCase.verifyEmpty(dir([testCase.DataPath, 'MIPs/*.tif']), 'MIPs were made from the incomplete crop.');
            % (and the blocks stop once the crop is there)
            testCase.verifyEqual(spy.getOption(1, 'finalOutFullpath'), cropFn);
        end

        function largeZarrKeepsAResultPublishedDuringItsMIPs(testCase)
            % as above, with the other crop coming while the MIPs are made (a
            % spy of XR_MIP_zarr makes it)
            testCase.writeDataset('zarr');
            cropFn = [testCase.DataPath, 'crop.zarr'];
            other = uint16(7 * ones(5, 5, 5));
            testCase.applyFixture(pktest.FunctionSpyFixture('XR_MIP_zarr', @() writezarr(other, cropFn)));
            evalc(['XR_crop_zarr([testCase.DataPath, testCase.Names{1}, ''.zarr''], cropFn, testCase.Bbox, ', ...
                '''batchSize'', [16, 16, 4], ''blockSize'', [8, 8, 4], ''parseCluster'', false)']);
            testCase.verifyEqual(readzarr(cropFn), other, 'The other writer''s crop was replaced.');
            testCase.verifyEmpty(dir([cropFn, '_*']), 'The temporary crop was not removed.');
            testCase.verifyEmpty(dir([cropFn, '/*.zarr_*']), 'The crop was moved into the other''s.');
        end

        function largeZarrWaitsForAWriterOfAnotherCall(testCase)
            % a crop of another call (another uuid) that may still make the
            % result: the second one made it as well; now it waits for it.
            % The other is a process on this node, whose crop comes during
            % the wait (a spy of wait_for_files); no block may run.
            testCase.assumeTrue(isunix && ~ismac, 'Writers outside Slurm are checked on Linux only.');
            testCase.writeDataset('zarr');
            cropFn = [testCase.DataPath, 'crop.zarr'];
            pktest.writeWriterRecord(testCase, cropFn, 'otheruuid');
            other = uint16(7 * ones(5, 5, 5));
            testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', @() pktest.atCall( ...
                [testCase.DataPath, 'calls.txt'], 1, @() writezarr(other, cropFn)), 10));
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture());
            evalc(['XR_crop_zarr([testCase.DataPath, testCase.Names{1}, ''.zarr''], cropFn, testCase.Bbox, ', ...
                '''batchSize'', [16, 16, 4], ''blockSize'', [8, 8, 4], ''parseCluster'', false)']);
            testCase.verifyEmpty(spy.getCalls(), 'The crop was made a second time.');
            testCase.verifyEqual(readzarr(cropFn), other);
        end

        function largeZarrBlocksOfAnEndedTaskAreSkipped(testCase)
            % the block jobs of a large crop whose job was cancelled ran on
            % (the queued ones when they started), for a result nothing
            % would publish; a block now skips its work when its task's job
            % has ended (its writer record; a fake Slurm that no longer
            % lists the job)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.writeDataset('zarr');
            cropFn = [testCase.DataPath, 'crop.zarr'];
            tmpFn = [cropFn, '_u1'];
            outSize = testCase.Bbox(4 : 6) - testCase.Bbox(1 : 3) + 1;
            createzarr(tmpFn, dataSize=outSize, blockSize=[8, 8, 4], dtype='uint16');
            pktest.writeWriterRecord(testCase, cropFn, 'u1', '777_3:othernode:5');
            flagFn = [testCase.DataPath, 'zarr_flag/crop_u1/blocks_1_1.mat'];
            testCase.applyFixture(pktest.FakeSlurmFixture);
            % (squeue's answers are kept per job for a minute; fake job ids repeat)
            clear('result_writers');
            box = [testCase.Bbox(1 : 3), testCase.Bbox(1 : 3) + [7, 7, 3]];
            evalc(['XR_crop_block(1, [testCase.DataPath, testCase.Names{1}, ''.zarr''], tmpFn, flagFn, box, ', ...
                '[1, 1, 1, 8, 8, 4], ''uuid'', ''u1'')']);
            testCase.verifyFalse(isfile(flagFn), 'The block of an ended task was cropped.');
            testCase.verifyEqual(nnz(readzarr(tmpFn)), 0);
        end

        function largeZarrTaskThatFailsMarksItsRecordEnded(testCase)
            % a crop that ends by an error (its blocks failed: a framework
            % spy) leaves its writer record, and its process may stay (an
            % interactive session): the record says it has ended, so that
            % no crop of another call waits for it and its blocks are skipped
            testCase.writeDataset('zarr');
            cropFn = [testCase.DataPath, 'crop.zarr'];
            testCase.applyFixture(pktest.FrameworkSpyFixture([], false));
            err = [];
            try
                evalc(['XR_crop_zarr([testCase.DataPath, testCase.Names{1}, ''.zarr''], cropFn, testCase.Bbox, ', ...
                    '''batchSize'', [16, 16, 4], ''blockSize'', [8, 8, 4], ''parseCluster'', false)']);
            catch err
            end
            testCase.assertNotEmpty(err, 'The crop did not fail.');
            testCase.assertSubstring(err.message, 'not finished');
            records = dir([testCase.DataPath, 'zarr_flag/crop_*/writer.txt']);
            testCase.assertNumElements(records, 1, 'The crop did not record itself.');
            testCase.verifySubstring(fileread([records.folder, '/', records.name]), sprintf('\nended=1'), ...
                'The crop that failed did not mark its record ended.');
        end

        function largeZarrJobsAskForWhatTheyHold(testCase)
            % the job of a large crop (XR_crop_zarr) runs its blocks as jobs
            % (1 CPU each). It asked for 32 B per voxel of the whole input:
            % whole nodes for large data. Now it asks for MATLAB and the
            % merge of the MIP slabs with one CPU, and if it crops blocks
            % itself (masterCompute, the config's if it sets it), for
            % cpusPerTask (default 2) and a block's estimate (8 B per voxel
            % of a batch, 1024 x 1024 x 512 here: 4 GB). The input is an 8
            % G-voxel header (the framework is a spy).
            dataPath = testCase.createTempFolder();
            createzarr([dataPath, 'scan_CamA_ch0_t0000.zarr'], dataSize=[4096, 4096, 512], blockSize=[256, 256, 256], ...
                dtype='uint16');
            offConfig = [dataPath, 'cpu_config.json'];
            fid = fopen(offConfig, 'w');
            fprintf(fid, '{"masterCompute": false}');
            fclose(fid);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture());
            cases = {true, ''; false, ''; true, offConfig};
            for i = 1 : size(cases, 1)
                XR_crop_dataset(dataPath, [1, 1, 1, 4096, 4096, 512], 'channelPatterns', {'CamA_ch0'}, 'zarrFile', true, ...
                    'largeFile', true, 'parseCluster', false, 'masterCompute', cases{i, 1}, 'configFile', cases{i, 2}, ...
                    'resultDirName', sprintf('Cropped_%d', i));
            end
            testCase.verifyEqual(spy.getOption(1, 'cpusPerTask'), 2);
            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), 4, 'AbsTol', 1e-9);
            for i = 2 : 3
                testCase.verifyEqual(spy.getOption(i, 'cpusPerTask'), 1, sprintf('Case %d.', i));
                testCase.verifyLessThan(spy.getOption(i, 'memAllocate'), 4, sprintf('Case %d.', i));
            end
            % (the wrapper's masterCompute reaches the crops)
            calls = spy.getCalls();
            testCase.verifyFalse(any(contains(calls{1}.funcStrs, '''masterCompute'',false')));
            testCase.verifyTrue(all(contains(calls{2}.funcStrs, '''masterCompute'',false')));
        end

        function largeZarrPassesJobSettingsToMIPs(testCase)
            % XR_crop_zarr passed parseCluster to XR_MIP_zarr, but not
            % masterCompute or jobLogDir
            testCase.writeDataset('zarr');
            spy = testCase.applyFixture(pktest.FunctionSpyFixture('XR_MIP_zarr'));
            XR_crop_zarr([testCase.DataPath, testCase.Names{1}, '.zarr'], [testCase.DataPath, 'crop.zarr'], ...
                testCase.Bbox, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 4], 'parseCluster', false, ...
                'masterCompute', false, 'jobLogDir', [testCase.DataPath, 'logs']);
            testCase.assertNotEmpty(spy.getCalls(), 'XR_MIP_zarr was not called.');
            testCase.verifyEqual(spy.getOption(1, 'parseCluster'), false);
            testCase.verifyEqual(spy.getOption(1, 'masterCompute'), false);
            testCase.verifyEqual(spy.getOption(1, 'jobLogDir'), [testCase.DataPath, 'logs']);
        end

        function padPastImageEdge(testCase)
            % bbox starts before and ends after the image in every axis
            testCase.writeDataset('tif');
            bbox = [-2, -1, 0, 45, 40, 15];
            XR_crop_dataset(testCase.DataPath, bbox, 'channelPatterns', testCase.Channels, 'pad', true, ...
                'parseCluster', false, 'masterCompute', true);
            out = readtiff([testCase.DataPath, 'Cropped/', testCase.Names{1}, '.tif']);
            testCase.verifyEqual(out, testCase.padCrop(testCase.Images{1}, bbox));
        end

        function padPastImageEdgeLargeZarr(testCase)
            % regression test for padding in large zarr cropping (81f3341)
            testCase.writeDataset('zarr');
            bbox = [-2, -1, 0, 45, 40, 15];
            XR_crop_dataset(testCase.DataPath, bbox, 'channelPatterns', testCase.Channels, 'pad', true, ...
                'zarrFile', true, 'largeFile', true, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 4], ...
                'parseCluster', false, 'masterCompute', true);
            out = readzarr([testCase.DataPath, 'Cropped/', testCase.Names{1}, '.zarr']);
            testCase.verifyEqual(out, testCase.padCrop(testCase.Images{1}, bbox));
        end

        function clipPastImageEdgeWithoutPad(testCase)
            testCase.writeDataset('tif');
            bbox = [30, 20, 8, 60, 50, 20];
            XR_crop_dataset(testCase.DataPath, bbox, 'channelPatterns', testCase.Channels, 'pad', false, ...
                'parseCluster', false, 'masterCompute', true);
            out = readtiff([testCase.DataPath, 'Cropped/', testCase.Names{1}, '.tif']);
            testCase.verifyEqual(out, testCase.Images{1}(30 : end, 20 : end, 8 : end));
        end

        function singleSliceZarrCrop(testCase)
            % regression test for single slice crops of zarr files (d513b73)
            testCase.writeDataset('zarr');
            bbox = [5, 7, 6, 30, 33, 6];
            XR_crop_dataset(testCase.DataPath, bbox, 'channelPatterns', testCase.Channels(1), ...
                'zarrFile', true, 'saveZarr', true, 'parseCluster', false, 'masterCompute', true, ...
                'resultDirName', 'Cropped_small');
            XR_crop_dataset(testCase.DataPath, bbox, 'channelPatterns', testCase.Channels(1), ...
                'zarrFile', true, 'largeFile', true, 'batchSize', [16, 16, 4], 'blockSize', [8, 8, 4], ...
                'parseCluster', false, 'masterCompute', true, 'resultDirName', 'Cropped_large');
            expected = testCase.crop(testCase.Images{1}, bbox);
            testCase.verifyEqual(readzarr([testCase.DataPath, 'Cropped_small/', testCase.Names{1}, '.zarr']), expected);
            testCase.verifyEqual(readzarr([testCase.DataPath, 'Cropped_large/', testCase.Names{1}, '.zarr']), expected);
        end

        function centerCrop(testCase)
            testCase.writeDataset('tif');
            XR_crop_dataset(testCase.DataPath, [1, 1, 1, 10, 9, 4], 'cropType', 'center', ...
                'channelPatterns', testCase.Channels(1), 'parseCluster', false, 'masterCompute', true);
            % start = round((imSize + 1) / 2) - floor(bboxSize / 2)
            start = round((testCase.ImSize + 1) / 2) - floor([10, 9, 4] / 2);
            out = readtiff([testCase.DataPath, 'Cropped/', testCase.Names{1}, '.tif']);
            testCase.verifyEqual(out, testCase.crop(testCase.Images{1}, [start, start + [10, 9, 4] - 1]));
        end

        function movingCrop(testCase)
            % the bbox moves linearly from inputBbox to lastStartCoords over
            % the time points of a channel
            dataPath = testCase.createTempFolder();
            nt = 3;
            ims = cell(nt, 1);
            for t = 1 : nt
                ims{t} = uint16(randi(60000, testCase.ImSize));
                writetiff(ims{t}, sprintf('%sscan_CamA_ch0_t%04d.tif', dataPath, t - 1));
            end
            bbox = [3, 3, 2, 12, 12, 6];
            lastStart = [21, 13, 6];
            XR_crop_dataset(dataPath, bbox, 'cropType', 'moving', 'lastStartCoords', lastStart, ...
                'channelPatterns', testCase.Channels(1), 'parseCluster', false, 'masterCompute', true);
            bboxSize = bbox(4 : 6) - bbox(1 : 3) + 1;
            for t = 1 : nt
                s = bbox(1 : 3) + round((lastStart - bbox(1 : 3)) / (nt - 1) * (t - 1));
                out = readtiff(sprintf('%sCropped/scan_CamA_ch0_t%04d.tif', dataPath, t - 1));
                testCase.verifyEqual(out, testCase.crop(ims{t}, [s, s + bboxSize - 1]), sprintf('time point %d', t));
            end
        end

        function existingOutputIsNotOverwritten(testCase)
            testCase.writeDataset('tif');
            resultPath = [testCase.DataPath, 'Cropped/'];
            mkdir(resultPath);
            sentinel = uint16(magic(4));
            writetiff(sentinel, [resultPath, testCase.Names{1}, '.tif']);
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels, ...
                'parseCluster', false, 'masterCompute', true);
            testCase.verifyEqual(readtiff([resultPath, testCase.Names{1}, '.tif']), sentinel);
            testCase.verifyEqual(readtiff([resultPath, testCase.Names{2}, '.tif']), testCase.crop(testCase.Images{2}, testCase.Bbox));
        end

        function tiffInputSavedAsZarr(testCase)
            % tif input saved as zarr used to get a zarr folder named *.tif
            testCase.writeDataset('tif');
            XR_crop_dataset(testCase.DataPath, testCase.Bbox, 'channelPatterns', testCase.Channels(1), ...
                'saveZarr', true, 'parseCluster', false, 'masterCompute', true);
            fn = [testCase.DataPath, 'Cropped/', testCase.Names{1}, '.zarr'];
            testCase.verifyTrue(exist(fn, 'dir') > 0, 'Cropped zarr output was not saved as <name>.zarr.');
            if exist(fn, 'dir')
                testCase.verifyEqual(readzarr(fn), testCase.crop(testCase.Images{1}, testCase.Bbox));
            end
        end
    end

    methods (Access = private)
        function [] = writeDataset(testCase, format)
            testCase.DataPath = testCase.createTempFolder();
            testCase.Names = {'scan_CamA_ch0_t0000', 'scan_CamB_ch0_t0000', 'scan_CamA_ch1_t0000'};
            testCase.Images = cell(1, 3);
            for i = 1 : 3
                testCase.Images{i} = uint16(randi(60000, testCase.ImSize));
                fn = [testCase.DataPath, testCase.Names{i}];
                if strcmp(format, 'zarr')
                    writezarr(testCase.Images{i}, [fn, '.zarr'], 'blockSize', [16, 16, 4]);
                else
                    writetiff(testCase.Images{i}, [fn, '.tif']);
                end
            end
        end
    end

    methods (Static, Access = private)
        function [out] = crop(im, b)
            out = im(b(1) : b(4), b(2) : b(5), b(3) : b(6));
        end

        function [out] = padCrop(im, b)
            % crop with zero padding outside the image
            out = zeros(b(4 : 6) - b(1 : 3) + 1, class(im));
            s = max(b(1 : 3), 1);
            t = min(b(4 : 6), size(im));
            o = s - b(1 : 3) + 1;
            out(o(1) : o(1) + t(1) - s(1), o(2) : o(2) + t(2) - s(2), o(3) : o(3) + t(3) - s(3)) = ...
                im(s(1) : t(1), s(2) : t(2), s(3) : t(3));
        end
    end
end
