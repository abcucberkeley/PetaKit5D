classdef DeskewRotateDataWrapperTest < pktest.TestCase
    % XR_deskew_rotate_data_wrapper on small synthetic datasets, run locally
    % (parseCluster=false). Outputs must equal the in-memory functions.

    properties (Constant)
        SkewAngle = 32.45
        XYPixelSize = 0.108
        Dz = 0.2
        RawSize = [64, 96, 80]
        FsnA0 = 'Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t'
        FsnB0 = 'Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t'
        FsnA1 = 'Scan_Iter_0001_0000_CamA_ch0_CAM1_stack0000_488nm_0011000msec_0106071251msecAbs_000x_000y_000z_0001t'
    end

    methods (Test, TestTags = {'Integration'})
        function processesEveryMatchingFile(testCase)
            [dataPath, frames] = testCase.writeDataset();
            testCase.runWrapper(dataPath);

            fsns = {testCase.FsnA0, testCase.FsnB0, testCase.FsnA1};
            for i = 1 : numel(fsns)
                dsr = readtiff([dataPath, 'DSR/', fsns{i}, '.tif']);
                testCase.verifyEqual(dsr, testCase.expectedDSR(frames{i}), fsns{i});
                testCase.verifyEqual(readtiff([dataPath, 'DSR/MIPs/', fsns{i}, '_MIP_z.tif']), max(dsr, [], 3));
            end
            testCase.verifyTrue(isfile([dataPath, 'DSR/parameters.json']));
            testCase.verifyTrue(isfile([dataPath, 'DSR/parameters.mat']));
            testCase.verifyEmpty(dir([dataPath, 'DSR/*.tmp']), 'Task flag files were left behind.');
            pr = jsondecode(fileread([dataPath, 'DSR/parameters.json']));
            testCase.verifyEqual(pr.dz, testCase.Dz);
            testCase.verifyEqual(pr.skewAngle, testCase.SkewAngle);
        end

        function channelPatternsSelectFiles(testCase)
            dataPath = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'channelPatterns', {'CamB'});
            dir_info = dir([dataPath, 'DSR/*.tif']);
            testCase.verifyEqual({dir_info.name}, {[testCase.FsnB0, '.tif']});
        end

        function multipleDataPathsAreProcessed(testCase)
            rootPath = testCase.createTempFolder();
            dataPaths = {[rootPath, 'tile1/'], [rootPath, 'tile2/']};
            frame = testCase.makeFrame(1);
            for i = 1 : 2
                writetiff(frame, [dataPaths{i}, testCase.FsnA0, '.tif']);
            end
            testCase.runWrapper(dataPaths);
            for i = 1 : 2
                testCase.verifyEqual(readtiff([dataPaths{i}, 'DSR/', testCase.FsnA0, '.tif']), testCase.expectedDSR(frame));
            end
        end

        function deskewOnlyWritesDS(testCase)
            [dataPath, frames] = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'rotate', false, 'DSRCombined', false);
            ds = readtiff([dataPath, 'DS/', testCase.FsnA0, '.tif']);
            testCase.verifyEqual(ds, uint16(deskewFrame3D(single(frames{1}), testCase.SkewAngle, testCase.Dz, ...
                testCase.XYPixelSize, true, 'Crop', false)));
            testCase.verifyFalse(isfolder([dataPath, 'DSR']));
        end

        function deskewOnlyWithCombinedDefaultWritesDS(testCase)
            % DSRCombined=true (the default) used to send deskew only to the
            % combined path, which writes DSR, and then wait minutes for DS
            % before failing
            [dataPath, frames] = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'rotate', false, 'DSRCombined', true);
            ds = readtiff([dataPath, 'DS/', testCase.FsnA0, '.tif']);
            testCase.verifyEqual(ds, uint16(deskewFrame3D(single(frames{1}), testCase.SkewAngle, testCase.Dz, ...
                testCase.XYPixelSize, true, 'Crop', false)));
            testCase.verifyFalse(isfolder([dataPath, 'DSR']));
        end

        function deskewOnlyIsCroppedWithCrop(testCase)
            % the wrapper didn't pass crop to XR_deskewRotateFrame: the
            % deskewed data was never cropped
            [dataPath, frames] = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'rotate', false, 'DSRCombined', false, 'crop', true);
            ds = readtiff([dataPath, 'DS/', testCase.FsnA0, '.tif']);
            testCase.verifyEqual(ds, uint16(deskewFrame3D(single(frames{1}), testCase.SkewAngle, testCase.Dz, ...
                testCase.XYPixelSize, true, 'Crop', true)));
        end

        function separateDeskewRotateWritesBoth(testCase)
            dataPath = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'DSRCombined', false);
            testCase.verifyTrue(isfile([dataPath, 'DS/', testCase.FsnA0, '.tif']));
            dsrSize = getImageSize([dataPath, 'DSR/', testCase.FsnA0, '.tif']);
            testCase.verifyEqual(dsrSize, size(testCase.expectedDSR(testCase.makeFrame(1))));
        end

        function existingOutputsAreSkippedUnlessOverwrite(testCase)
            [dataPath, frames] = testCase.writeDataset();
            testCase.runWrapper(dataPath);
            dsrFn = [dataPath, 'DSR/', testCase.FsnA0, '.tif'];

            writetiff(frames{1} * 2, [dataPath, testCase.FsnA0, '.tif']);
            testCase.runWrapper(dataPath);
            testCase.verifyEqual(readtiff(dsrFn), testCase.expectedDSR(frames{1}));

            testCase.runWrapper(dataPath, 'overwrite', true);
            testCase.verifyEqual(readtiff(dsrFn), testCase.expectedDSR(frames{1} * 2));
        end

        function saveZarrMatchesTiff(testCase)
            [dataPath, frames] = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'saveZarr', true);
            testCase.verifyEqual(readzarr([dataPath, 'DSR/', testCase.FsnA0, '.zarr']), testCase.expectedDSR(frames{1}));
        end

        function zarrInputMatchesTiffInput(testCase)
            dataPath = testCase.createTempFolder();
            frame = testCase.makeFrame(1);
            writezarr(frame, [dataPath, testCase.FsnA0, '.zarr']);
            testCase.runWrapper(dataPath, 'zarrFile', true);
            testCase.verifyEqual(readtiff([dataPath, 'DSR/', testCase.FsnA0, '.tif']), testCase.expectedDSR(frame));
        end

        function largeFileMatchesRegular(testCase)
            % the large file path splits the output along y into batches
            % (4 here) and must reproduce the in-memory result exactly
            dataPath = testCase.createTempFolder();
            frame = testCase.makeFrame(1);
            writezarr(frame, [dataPath, testCase.FsnA0, '.zarr']);
            testCase.runWrapper(dataPath, 'zarrFile', true, 'largeFile', true, 'batchSize', [16, 1024, 1024], ...
                'blockSize', [16, 64, 64]);
            dsr = readzarr([dataPath, 'DSR/', testCase.FsnA0, '.zarr']);
            testCase.verifyEqual(dsr, testCase.expectedDSR(frame));
            testCase.verifyEqual(readtiff([dataPath, 'DSR/MIPs/', testCase.FsnA0, '_MIP_z.tif']), max(dsr, [], 3));
        end

        function largeFileReadsTheResampleFactorAsRegularFiles(testCase)
            % XR_deskewRotateZarr read the factor as [y x z] for the result's
            % size and its batches, while its blocks (deskewRotateFrame3D) and
            % the regular path read [x y z], and 1 or 2 values failed with an
            % index error. The result is the regular path's, completed to 3
            % values as XR_deskewRotateFrame does; a fractional factor along y
            % (the batches' axis) only gives its size here
            dataPath = testCase.createTempFolder();
            frame = testCase.makeFrame(1);
            writezarr(frame, [dataPath, testCase.FsnA0, '.zarr']);
            % factors, the [x y z] factor of the regular path, and whether the
            % values are compared (integer factors and 1 along y give the
            % regular path's values)
            factors = {[1.5, 1, 1.25], [1.5, 1, 1.25], true; [1, 2, 1], [1, 2, 1], true; 2, [2, 2, 2], true; ...
                [2, 1], [2, 2, 1], true; [1, 1.5, 1.25], [1, 1.5, 1.25], false};
            for i = 1 : size(factors, 1)
                name = sprintf('DSR_%d', i);
                testCase.runWrapper(dataPath, 'zarrFile', true, 'largeFile', true, 'batchSize', [16, 1024, 1024], ...
                    'blockSize', [16, 64, 64], 'saveMIP', false, 'resampleFactor', factors{i, 1}, 'DSRDirName', name);
                dsr = readzarr([dataPath, name, '/', testCase.FsnA0, '.zarr']);
                expected = deskewRotateFrame3D(frame, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                    'reverse', true, 'save16bit', true, 'resampleFactor', factors{i, 2});
                testCase.verifySize(dsr, size(expected), mat2str(factors{i, 1}));
                if factors{i, 3}
                    testCase.verifyEqual(dsr, expected, mat2str(factors{i, 1}));
                end
            end
        end

        function largeFileIsFlippedAsRegularFiles(testCase)
            % the large-file path took flipZstack (the wrapper's, not the one
            % of each file's setting file) and never flipped the data
            dataPath = testCase.createTempFolder();
            frame = testCase.makeFrame(1);
            writezarr(frame, [dataPath, testCase.FsnA0, '.zarr']);
            args = {'zarrFile', true, 'largeFile', true, 'batchSize', [16, 1024, 1024], 'blockSize', [16, 64, 64], ...
                'saveMIP', false};
            testCase.runWrapper(dataPath, args{:}, 'flipZstack', true);
            testCase.verifyEqual(readzarr([dataPath, 'DSR/', testCase.FsnA0, '.zarr']), testCase.expectedDSR(flip(frame, 3)));

            % a setting file with a negative stage interval flips its file
            testCase.writeStageInterval([dataPath, testCase.FsnA0, '_Settings.txt'], '-0.2');
            testCase.runWrapper(dataPath, args{:}, 'parseSettingFile', true, 'DSRDirName', 'DSR_settings');
            testCase.verifyEqual(readzarr([dataPath, 'DSR_settings/', testCase.FsnA0, '.zarr']), ...
                testCase.expectedDSR(flip(frame, 3)));
        end

        function largeFileMIPMaskBoxesFollowTheFactorAndTheFlip(testCase)
            % XR_getDeskeRotateBoxesFromMasks took x from the factor's second
            % value and the output rows from its first ([y x z]; one value or
            % none failed), and had no flip. The output box is the box of the
            % regular path's result of the masked region, and with the flip
            % that of the flipped masks (the input box stays in the
            % coordinates of the raw data: it is flipped after the crop)
            root = testCase.createTempFolder();
            m = false(testCase.RawSize);
            m(10 : 50, 20 : 40, 15 : 30) = true;
            m(10 : 50, 50 : 70, 40 : 60) = true;
            maskFns = testCase.writeMIPMasks([root, 'raw_'], m);
            flippedFns = testCase.writeMIPMasks([root, 'flipped_'], flip(m, 3));
            args = {[0, 0, 0], testCase.XYPixelSize, testCase.Dz, testCase.SkewAngle, true};
            factors = {[1, 1, 1], [1.5, 1, 1.25], [1, 1.5, 1.25], 1.5, [1.5, 1]};
            outBboxes = cell(size(factors));
            for i = 1 : numel(factors)
                rs = factors{i};
                msg = mat2str(rs);
                r = [ones(1, 4 - numel(rs)) * rs(1), rs(2 : end)];
                [inBbox, outBboxes{i}] = XR_getDeskeRotateBoxesFromMasks(maskFns, args{:}, rs);
                testCase.verifyEqual(inBbox, [10, 20, 15, 50, 70, 60], msg);
                % (as the masks' MIP: 16 bits without resampling along x and z)
                dsr = deskewRotateFrame3D(uint16(m(10 : 50, 20 : 70, 15 : 60)), testCase.SkewAngle, testCase.Dz, ...
                    testCase.XYPixelSize, 'reverse', true, 'resampleFactor', r, 'save16bit', all(r([1, 3]) == 1));
                xz = squeeze(max(dsr, [], 1)) > 0;
                xs = find(any(xz, 2));
                zs = find(any(xz, 1));
                testCase.verifyEqual(outBboxes{i}, [1, xs(1), zs(1), size(dsr, 1), xs(end), zs(end)], msg);
            end
            for i = 1 : numel(factors)
                rs = factors{i};
                msg = mat2str(rs);
                [inFlip, outFlip, mipFlip] = XR_getDeskeRotateBoxesFromMasks(maskFns, args{:}, rs, [], true);
                [inFlipped, outFlipped, mipFlipped] = XR_getDeskeRotateBoxesFromMasks(flippedFns, args{:}, rs);
                testCase.verifyEqual(inFlip, inBbox, msg);
                testCase.verifyEqual(inFlipped([3, 6]), testCase.RawSize(3) + 1 - inBbox([6, 3]), msg);
                testCase.verifyEqual(outFlip, outFlipped, msg);
                % (the deskewed and rotated MIP the box is taken from)
                testCase.verifyEqual(mipFlip, mipFlipped, msg);
                testCase.verifyNotEqual(outFlip, outBboxes{i}, msg);
            end
        end

        function largeFileParforRunsTheBlocks(testCase)
            % XR_deskewRotateZarr with parseParfor and no cluster failed on an
            % undefined GPUJob. The parfor wrapper is a stand-in that runs the
            % tasks in this process (GenericComputingFrameworkTest runs the
            % real one on a pool)
            dataPath = testCase.createTempFolder();
            frame = testCase.makeFrame(1);
            zarrFn = [dataPath, testCase.FsnA0, '.zarr'];
            writezarr(frame, zarrFn);
            optionsFn = testCase.writeParforStandIn();
            evalc(['XR_deskewRotateZarr(zarrFn, testCase.XYPixelSize, testCase.Dz, ''skewAngle'', testCase.SkewAngle, ', ...
                '''reverse'', true, ''batchSize'', [16, 1024, 1024], ''blockSize'', [16, 64, 64], ''saveMIP'', false, ', ...
                '''parseCluster'', false, ''parseParfor'', true)']);
            testCase.verifyEqual(readzarr([dataPath, 'DSR/', testCase.FsnA0, '.zarr']), testCase.expectedDSR(frame));
            a = load(optionsFn);
            testCase.verifyEqual(a.GPUJob, false);
        end

        function largeFilePassesJobSettingsToMIPs(testCase)
            % XR_deskewRotateZarr didn't pass parseCluster, masterCompute or
            % jobLogDir to XR_MIP_zarr, so on a machine with Slurm the MIPs
            % of a local run were made in Slurm jobs
            dataPath = testCase.createTempFolder();
            frame = testCase.makeFrame(1);
            writezarr(frame, [dataPath, testCase.FsnA0, '.zarr']);
            spy = testCase.applyFixture(pktest.FunctionSpyFixture('XR_MIP_zarr'));
            testCase.runWrapper(dataPath, 'zarrFile', true, 'largeFile', true, 'batchSize', [16, 1024, 1024], ...
                'blockSize', [16, 64, 64], 'jobLogDir', [dataPath, 'logs']);
            testCase.assertNotEmpty(spy.getCalls(), 'XR_MIP_zarr was not called.');
            testCase.verifyEqual(spy.getOption(1, 'parseCluster'), false);
            testCase.verifyEqual(spy.getOption(1, 'masterCompute'), true);
            testCase.verifyEqual(spy.getOption(1, 'jobLogDir'), [dataPath, 'logs']);
        end

        function largeFileMakesMIPsBeforePublishing(testCase)
            % the job framework counts a task as done, and may cancel its
            % job, once the result exists; the large file path moved the
            % result into place and then made its MIPs
            dataPath = testCase.createTempFolder();
            writezarr(testCase.makeFrame(1), [dataPath, testCase.FsnA0, '.zarr']);
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            testCase.runWrapper(dataPath, 'zarrFile', true, 'largeFile', true, 'batchSize', [16, 1024, 1024], ...
                'blockSize', [16, 64, 64]);
            dsrFn = [dataPath, 'DSR/', testCase.FsnA0, '.zarr'];
            testCase.assertTrue(moves.wasMoved(dsrFn), 'The result was not moved into place.');
            published = moves.filesWhenMoved(dsrFn);
            mips = strcat('MIPs/', testCase.FsnA0, {'_MIP_y.tif'; '_MIP_x.tif'; '_MIP_z.tif'});
            testCase.verifyEqual(ismember(mips, published), true(3, 1), 'MIPs were made after the result was moved into place.');
            added = setdiff(pktest.MovefileSpyFixture.listFiles([dataPath, 'DSR']), [published; {[testCase.FsnA0, '.zarr']}]);
            testCase.verifyEmpty(added, 'Files were written after the result was moved into place.');
        end

        function largeFileMIPsAreMergedFromTheBlocks(testCase)
            % the block tasks save the MIPs of their blocks in their flags,
            % and XR_MIP_zarr merges them into the MIP slabs instead of
            % reading the result again in MIP_block tasks
            dataPath = testCase.createTempFolder();
            writezarr(testCase.makeFrame(1), [dataPath, testCase.FsnA0, '.zarr']);
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            % MIP_block tasks would read the result; the spy errors at once, as the
            % framework would retry a task without output for minutes
            mipTasks = testCase.applyFixture(pktest.FunctionSpyFixture('MIP_block', [], 0));
            testCase.runWrapper(dataPath, 'zarrFile', true, 'largeFile', true, 'batchSize', [16, 1024, 1024], ...
                'blockSize', [16, 64, 64]);
            testCase.verifyEmpty(mipTasks.getCalls(), 'The MIPs were computed from the result.');
            % a task's flag holds its MIPs and marks it done: it is written
            % under a temporary name first
            calls = moves.getCalls();
            testCase.verifyTrue(any(contains({calls.dst}, '/zarr_flag/') & endsWith({calls.dst}, '.mat')), ...
                'The block flags were not moved into place.');
            pktest.verifyLargeFileMIPs(testCase, [dataPath, 'DSR/', testCase.FsnA0, '.zarr']);
        end

        function largeFileBlocksGetTheMIPOptionOnlyForMIPs(testCase)
            % compiled masters built before the block tasks' mipBatchSize
            % option still run the tasks of results without MIPs
            dataPath = testCase.createTempFolder();
            frameFn = [dataPath, testCase.FsnA0, '.zarr'];
            writezarr(testCase.makeFrame(1), frameFn);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture());
            for saveMIP = [false, true]
                evalc(['XR_deskewRotateZarr(frameFn, testCase.XYPixelSize, testCase.Dz, ''skewAngle'', testCase.SkewAngle, ', ...
                    '''reverse'', true, ''batchSize'', [16, 1024, 1024], ''blockSize'', [16, 64, 64], ''saveMIP'', saveMIP, ', ...
                    '''resultDirStr'', sprintf(''DSR_%d/'', saveMIP), ''parseCluster'', false)']);
            end
            calls = spy.getCalls();
            testCase.verifyFalse(any(contains(calls{1}.funcStrs, 'mipBatchSize')), 'Tasks without MIPs got mipBatchSize.');
            testCase.verifyTrue(all(contains(calls{2}.funcStrs, '''mipBatchSize'',[')), 'Tasks with MIPs got no mipBatchSize.');
        end

        function largeFileWithFailedBlocksIsNotPublished(testCase)
            % when blocks still failed after the retries, the large file path
            % moved the incomplete result into place, so it counted as done
            dataPath = testCase.createTempFolder();
            frameFn = [dataPath, testCase.FsnA0, '.zarr'];
            writezarr(testCase.makeFrame(1), frameFn);
            % a block framework that reports every block as failed; without
            % MIPs, which would fail on the missing blocks too
            testCase.applyFixture(pktest.FrameworkSpyFixture([], false));
            testCase.verifyError(@() XR_deskewRotateZarr(frameFn, testCase.XYPixelSize, testCase.Dz, ...
                'skewAngle', testCase.SkewAngle, 'reverse', true, 'batchSize', [16, 1024, 1024], ...
                'blockSize', [16, 64, 64], 'saveMIP', false, 'parseCluster', false), ?MException);
            testCase.verifyFalse(isfolder([dataPath, 'DSR/', testCase.FsnA0, '.zarr']), ...
                'An incomplete result was moved into place.');
        end

        function largeFileKeepsAResultAnotherWriterPublished(testCase)
            % two job sets on the same data: the one that finished last
            % deleted the result the other had published (which a reader
            % may use) and moved its own into place. Now it keeps the
            % other's and removes its own (incomplete here: the framework
            % is a spy that reports the blocks done while the other result
            % comes)
            dataPath = testCase.createTempFolder();
            frameFn = [dataPath, testCase.FsnA0, '.zarr'];
            writezarr(testCase.makeFrame(1), frameFn);
            dsrFn = [dataPath, 'DSR/', testCase.FsnA0, '.zarr'];
            other = uint16(7 * ones(10, 10, 10));
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture(@() pktest.atCall([dataPath, 'calls.txt'], 1, ...
                @() writezarr(other, dsrFn))));
            evalc(['XR_deskewRotateZarr(frameFn, testCase.XYPixelSize, testCase.Dz, ''skewAngle'', testCase.SkewAngle, ', ...
                '''reverse'', true, ''batchSize'', [16, 1024, 1024], ''blockSize'', [16, 64, 64], ''parseCluster'', false)']);
            testCase.verifyEqual(readzarr(dsrFn), other, 'The other writer''s result was replaced.');
            testCase.verifyEmpty(dir([dsrFn, '_*']), 'The temporary result was not removed.');
            testCase.verifyEmpty(dir([dataPath, 'DSR/zarr_flag/', testCase.FsnA0, '_*']), 'The flags were not removed.');
            % (no MIPs of this one's incomplete result, named after the result)
            testCase.verifyEmpty(dir([dataPath, 'DSR/MIPs/*.tif']), 'MIPs were made from the incomplete result.');
            % (and the blocks stop once the result is there)
            testCase.verifyEqual(regexprep(spy.getOption(1, 'finalOutFullpath'), '/+', '/'), dsrFn);
        end

        function largeFileKeepsAResultPublishedDuringItsMIPs(testCase)
            % as above, with the other result coming while the MIPs are made
            % (a spy of XR_MIP_zarr makes it)
            dataPath = testCase.createTempFolder();
            frameFn = [dataPath, testCase.FsnA0, '.zarr'];
            writezarr(testCase.makeFrame(1), frameFn);
            dsrFn = [dataPath, 'DSR/', testCase.FsnA0, '.zarr'];
            other = uint16(7 * ones(10, 10, 10));
            testCase.applyFixture(pktest.FunctionSpyFixture('XR_MIP_zarr', @() writezarr(other, dsrFn)));
            evalc(['XR_deskewRotateZarr(frameFn, testCase.XYPixelSize, testCase.Dz, ''skewAngle'', testCase.SkewAngle, ', ...
                '''reverse'', true, ''batchSize'', [16, 1024, 1024], ''blockSize'', [16, 64, 64], ''parseCluster'', false)']);
            testCase.verifyEqual(readzarr(dsrFn), other, 'The other writer''s result was replaced.');
            testCase.verifyEmpty(dir([dsrFn, '_*']), 'The temporary result was not removed.');
            testCase.verifyEmpty(dir([dsrFn, '/*.zarr_*']), 'The result was moved into the other''s.');
        end

        function largeFileBlocksOfAnEndedTaskAreSkipped(testCase)
            % the block jobs of a large file whose job was cancelled ran on,
            % for a result nothing would publish; a block now skips its work
            % when its task's job has ended (its writer record; a fake Slurm
            % that no longer lists the job). The block task is the one the
            % task made (a framework spy that reports it failed).
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            frameFn = [dataPath, testCase.FsnA0, '.zarr'];
            writezarr(testCase.makeFrame(1), frameFn);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture([], false));
            testCase.verifyError(@() XR_deskewRotateZarr(frameFn, testCase.XYPixelSize, testCase.Dz, ...
                'skewAngle', testCase.SkewAngle, 'reverse', true, 'batchSize', [16, 1024, 1024], ...
                'blockSize', [16, 64, 64], 'saveMIP', false, 'parseCluster', false), ?MException);
            calls = spy.getCalls();
            flagFn = calls{1}.outputFullpaths{1};
            recordFn = [fileparts(flagFn), '/writer.txt'];
            testCase.assertTrue(isfile(recordFn), 'The task did not record itself.');
            % (it ended by an error, and says so: its process may stay, e.g.
            % an interactive session, and no one may wait for it)
            testCase.verifySubstring(fileread(recordFn), sprintf('\nended=1'), ...
                'The task that failed did not mark its record ended.');
            % a task whose job was cancelled can't say so: its job is gone
            txt = regexprep(fileread(recordFn), {'owner=[^\n]*', 'ended=[^\n]*\n'}, {'owner=777_3:othernode:5', ''});
            fid = fopen(recordFn, 'w');
            fprintf(fid, '%s', txt);
            fclose(fid);
            testCase.applyFixture(pktest.FakeSlurmFixture);
            % (squeue's answers are kept per job for a minute; fake job ids repeat)
            clear('result_writers');
            funcStr = calls{1}.funcStrs{1};
            evalc('feval(str2func([''@()'', funcStr]))');
            testCase.verifyFalse(isfile(flagFn), 'The block of an ended task was deskewed.');
        end

        function largeFileWaitsForAWriterOfAnotherCall(testCase)
            % a task of another call (another uuid) that may still make the
            % result: the second one made it as well; now it waits for it.
            % The other is a process on this node, whose result comes during
            % the wait (a spy of wait_for_files); no block may run.
            testCase.assumeTrue(isunix && ~ismac, 'Writers outside Slurm are checked on Linux only.');
            dataPath = testCase.createTempFolder();
            frameFn = [dataPath, testCase.FsnA0, '.zarr'];
            writezarr(testCase.makeFrame(1), frameFn);
            dsrFn = [dataPath, 'DSR/', testCase.FsnA0, '.zarr'];
            pktest.writeWriterRecord(testCase, dsrFn, 'otheruuid');
            other = uint16(7 * ones(10, 10, 10));
            testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', @() pktest.atCall([dataPath, 'calls.txt'], 1, ...
                @() writezarr(other, dsrFn)), 10));
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture());
            evalc(['XR_deskewRotateZarr(frameFn, testCase.XYPixelSize, testCase.Dz, ''skewAngle'', testCase.SkewAngle, ', ...
                '''reverse'', true, ''batchSize'', [16, 1024, 1024], ''blockSize'', [16, 64, 64], ''parseCluster'', false)']);
            testCase.verifyEmpty(spy.getCalls(), 'The result was made a second time.');
            testCase.verifyEqual(readzarr(dsrFn), other);
        end

        function largeFileJobsGoAfterTheirBlocks(testCase)
            % the job of a large file waits for its block jobs: it asks to
            % go after them (nice, coordinator_nice); the jobs of other
            % files don't (the framework is a spy)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            createzarr([dataPath, testCase.FsnA0, '.zarr'], dataSize=[512, 256, 64], blockSize=[64, 64, 64], ...
                dtype='uint16');
            testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.applyFixture(pktest.FunctionStubFixture('coordinator_nice', 1000));
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture());
            for largeFile = [true, false]
                testCase.runWrapper(dataPath, 'zarrFile', true, 'largeFile', largeFile, 'parseCluster', true, ...
                    'DSRDirName', sprintf('DSR_%d', largeFile), 'jobLogDir', [dataPath, 'logs']);
            end
            testCase.verifyEqual(spy.getOption(1, 'nice'), 1000);
            testCase.verifyEqual(spy.getOption(2, 'nice'), 0);
        end

        function largeFileJobsAskForWhatTheyHold(testCase)
            % the job of a large file (XR_deskewRotateZarr) runs its y
            % batches as jobs. It asked for 40 B per voxel of the whole
            % input (masterCompute) or of a y slab (without, on one CPU):
            % whole nodes for large data. Now it asks for MATLAB and the
            % merge of the MIP slabs on one CPU, and if it computes batches
            % itself (masterCompute, the config's if it sets it), for a
            % batch's estimate on the CPUs of a batch's job. The input is an
            % 8 G-voxel header (the framework is a spy).
            dataPath = testCase.createTempFolder();
            createzarr([dataPath, testCase.FsnA0, '.zarr'], dataSize=[8192, 2048, 512], blockSize=[256, 256, 256], ...
                dtype='uint16');
            offConfig = [dataPath, 'cpu_config.json'];
            fid = fopen(offConfig, 'w');
            fprintf(fid, '{"masterCompute": false}');
            fclose(fid);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture());
            cases = {true, ''; false, ''; true, offConfig};
            for i = 1 : size(cases, 1)
                testCase.runWrapper(dataPath, 'zarrFile', true, 'largeFile', true, 'masterCompute', cases{i, 1}, ...
                    'configFile', cases{i, 2}, 'cpusPerTask', 8, 'DSRDirName', sprintf('DSR_%d', i));
            end
            % a batch of 1024 rows across x and z, less than half of the
            % input in single (32 GB)
            testCase.verifyEqual(spy.getOption(1, 'cpusPerTask'), 8);
            testCase.verifyGreaterThan(spy.getOption(1, 'memAllocate'), 2);
            testCase.verifyLessThan(spy.getOption(1, 'memAllocate'), 16);
            for i = 2 : 3
                testCase.verifyEqual(spy.getOption(i, 'cpusPerTask'), 1, sprintf('Case %d.', i));
                testCase.verifyLessThan(spy.getOption(i, 'memAllocate'), 4, sprintf('Case %d.', i));
            end
        end

        function flatFieldUsesEachChannelsImages(testCase)
            [dataPath, frames] = testCase.writeDataset();
            bg = zeros(testCase.RawSize(1 : 2), 'single') + 100;
            ffA = bg + single(linspace(500, 1000, testCase.RawSize(2)));
            ffB = bg + single(linspace(1000, 500, testCase.RawSize(2)));
            writetiff(ffA, [dataPath, 'ff/ffA.tif']);
            writetiff(ffB, [dataPath, 'ff/ffB.tif']);
            writetiff(bg, [dataPath, 'ff/bg.tif']);
            testCase.runWrapper(dataPath, 'FFCorrection', true, 'lowerLimit', 0.4, ...
                'FFImagePaths', {[dataPath, 'ff/ffA.tif'], [dataPath, 'ff/ffB.tif']}, ...
                'backgroundPaths', {[dataPath, 'ff/bg.tif'], [dataPath, 'ff/bg.tif']});

            expectedA = testCase.expectedDSR(XR_LSFlatFieldCorrection(frames{1}, ffA, bg, 'LowerLimit', 0.4));
            expectedB = testCase.expectedDSR(XR_LSFlatFieldCorrection(frames{2}, ffB, bg, 'LowerLimit', 0.4));
            testCase.verifyEqual(readtiff([dataPath, 'DSR/', testCase.FsnA0, '.tif']), expectedA);
            testCase.verifyEqual(readtiff([dataPath, 'DSR/', testCase.FsnB0, '.tif']), expectedB);
        end

        function partialFilesAreJoinedAlongZ(testCase)
            % acquisitions split into <name>.tif + <name>_part0001.tif are
            % processed as one volume
            dataPath = testCase.createTempFolder();
            frame = testCase.makeFrame(1);
            writetiff(frame(:, :, 1 : 40), [dataPath, testCase.FsnA0, '.tif']);
            writetiff(frame(:, :, 41 : end), [dataPath, testCase.FsnA0, '_part0001.tif']);
            testCase.runWrapper(dataPath);
            testCase.verifyEqual(readtiff([dataPath, 'DSR/', testCase.FsnA0, '.tif']), testCase.expectedDSR(frame));
            testCase.verifyFalse(isfile([dataPath, 'DSR/', testCase.FsnA0, '_part0001.tif']));
        end

        function saveMIPFalseWritesNoMIP(testCase)
            % saveMIP and blockSize used to be dropped from the task string
            % of regular (not large) files
            dataPath = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'saveMIP', false);
            testCase.verifyFalse(isfolder([dataPath, 'DSR/MIPs']), 'MIPs were written although saveMIP is false.');
        end

        function blockSizeIsUsedForZarrOutput(testCase)
            dataPath = testCase.writeDataset();
            testCase.runWrapper(dataPath, 'saveZarr', true, 'blockSize', [16, 32, 16]);
            info = getZarrInfo([dataPath, 'DSR/', testCase.FsnA0, '.zarr']);
            testCase.verifyEqual(info.blockSize(:)', [16, 32, 16]);
        end

        function objectiveScanLargeFileIsRejected(testCase)
            % the large-file block processing only supports sample scan; it
            % used to fail inside the tasks with an undefined variable
            dataPath = testCase.createTempFolder();
            writezarr(testCase.makeFrame(1), [dataPath, testCase.FsnA0, '.zarr']);
            try
                testCase.runWrapper(dataPath, 'zarrFile', true, 'largeFile', true, 'objectiveScan', true, ...
                    'batchSize', [16, 1024, 1024]);
                testCase.verifyFail('Large-file objective scan did not raise an error.');
            catch ME
                testCase.verifySubstring(ME.message, 'objective scan is not supported');
            end
        end
    end

    methods (Access = private)
        function [] = runWrapper(testCase, dataPaths, varargin)
            if ischar(dataPaths)
                dataPaths = {dataPaths};
            end
            XR_deskew_rotate_data_wrapper(dataPaths, 'deskew', true, 'rotate', true, 'DSRCombined', true, ...
                'xyPixelSize', testCase.XYPixelSize, 'dz', testCase.Dz, 'skewAngle', testCase.SkewAngle, ...
                'reverse', true, 'channelPatterns', {'CamA', 'CamB'}, 'save16bit', true, ...
                'parseCluster', false, 'masterCompute', true, varargin{:});
        end

        function [frame] = makeFrame(testCase, seed)
            % uint16 beads on a background of 100, different per seed
            rs = RandStream('twister', 'Seed', seed);
            raw = [8 + 48 * rand(rs, 15, 1), 26 + 44 * rand(rs, 15, 1), 10 + 60 * rand(rs, 15, 1)];
            frame = uint16(pktest.beads(testCase.RawSize, raw, 2, 1000, 100));
        end

        function [dataPath, frames] = writeDataset(testCase)
            dataPath = testCase.createTempFolder();
            fsns = {testCase.FsnA0, testCase.FsnB0, testCase.FsnA1};
            frames = cell(1, numel(fsns));
            for i = 1 : numel(fsns)
                frames{i} = testCase.makeFrame(i);
                writetiff(frames{i}, [dataPath, fsns{i}, '.tif']);
            end
        end

        function [dsr] = expectedDSR(testCase, frame)
            dsr = deskewRotateFrame3D(frame, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                'reverse', true, 'save16bit', true);
        end

        function [] = writeStageInterval(~, fn, interval)
            % a setting file with the stage interval (negative: the stack is
            % flipped)
            lines = {'***** ***** *****   General   ***** ***** *****', 'Date :  1/1/2024 10:00:00 AM', '', ...
                sprintf('X Stage Offset, Interval (um), # of Pixels for Excitation (0) :\t0\t%s\t101', interval), ''};
            fid = fopen(fn, 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
        end

        function [fns] = writeMIPMasks(~, prefix, m)
            % the MIP masks of a mask volume, in xy, xz and yz order
            fns = {[prefix, 'xy.tif'], [prefix, 'xz.tif'], [prefix, 'yz.tif']};
            writetiff(uint8(max(m, [], 3)), fns{1});
            writetiff(uint8(squeeze(max(m, [], 1))), fns{2});
            writetiff(uint8(squeeze(max(m, [], 2))), fns{3});
        end

        function [optionsFn] = writeParforStandIn(testCase)
            % a matlab_parfor_generic_computing_wrapper first on the path that
            % saves its GPUJob and nworker and runs the tasks in this process
            folder = [testCase.createTempFolder(), 'parfor/'];
            mkdir(folder);
            optionsFn = [folder, 'options.mat'];
            lines = {
                'function [is_done_flag] = matlab_parfor_generic_computing_wrapper(inputFullpaths, outputFullpaths, funcStrs, varargin)'
                '% test stand-in written by DeskewRotateDataWrapperTest'
                'ip = inputParser;'
                'ip.KeepUnmatched = true;'
                'ip.addParameter(''GPUJob'', true);'
                'ip.addParameter(''nworker'', 24);'
                'ip.parse(varargin{:});'
                'GPUJob = ip.Results.GPUJob;'
                'nworker = ip.Results.nworker;'
                sprintf('save(''%s'', ''GPUJob'', ''nworker'');', optionsFn)
                'for i = 1 : numel(funcStrs)'
                '    eval(funcStrs{i});'
                'end'
                'is_done_flag = true(numel(funcStrs), 1);'
                'end'};
            fid = fopen([folder, 'matlab_parfor_generic_computing_wrapper.m'], 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
            testCase.addTeardown(@() clear('matlab_parfor_generic_computing_wrapper'));
            testCase.addTeardown(@() rmpath(folder));
            pktest.addFunctionFolder(folder, 'matlab_parfor_generic_computing_wrapper');
        end
    end
end
