classdef AutomaticProcessingTest < pktest.TestCase
    % XR_microscopeAutomaticProcessing on existing data (streaming false),
    % following demos/demo_microscope_pipeline_single_dataset*.m and
    % demos/demo_microscope_pipeline_multiple_dataset.m: flat field
    % correction, combined deskew/rotation, and stitching in DSR space.
    %
    % Expected results come from the underlying functions run directly.
    % Deskew/rotation works row by row in y, so tiles cut from a raw volume
    % along y deskew/rotate and stitch back to the DSR of the whole volume.

    properties (Constant)
        XYPixelSize = 0.108
        Dz = 0.3
        SkewAngle = 32.45
    end

    methods (Test, TestTags = {'Integration'})
        function flatFieldAndDeskewRotate(testCase)
            % demo_microscope_pipeline_single_dataset_nostitch.m: two
            % channels with per channel flat field and background images
            dataPath = testCase.createTempFolder();
            chs = {'CamB_ch0', 'CamA_ch0'};
            raw = cell(1, 2);
            ffFns = cell(1, 2);
            bgFns = cell(1, 2);
            for c = 1 : 2
                raw{c} = uint16(pktest.beads([64, 64, 40], [32, 30, 20; 20, 40, 12], 2, 2000, 300));
                writetiff(raw{c}, testCase.tileFilename(dataPath, chs{c}, 0));
                % flat field falling off to 50% at the right edge
                ffFns{c} = sprintf('%sFF/ff_%s.tif', dataPath, chs{c});
                bgFns{c} = sprintf('%sFF/bg_%s.tif', dataPath, chs{c});
                writetiff(uint16(100 + repmat(linspace(1000, 500, 64) * c, 64, 1)), ffFns{c});
                writetiff(100 * ones(64, 64, 'uint16'), bgFns{c});
            end

            XR_microscopeAutomaticProcessing(dataPath, 'xyPixelSize', testCase.XYPixelSize, 'dz', testCase.Dz, ...
                'reverse', true, 'channelPatterns', chs, 'save16bit', [true, true], 'overwrite', false, ...
                'streaming', false, 'deskew', true, 'rotate', true, 'stitch', false, 'DSRCombined', true, ...
                'FFCorrection', true, 'lowerLimit', 0.4, 'constOffset', 1, 'FFImagePaths', ffFns, ...
                'backgroundPaths', bgFns, 'minModifyTime', 1, 'maxModifyTime', 10, 'maxWaitLoopNum', 10, ...
                'parseCluster', false, 'cpusPerTask', 4);

            for c = 1 : 2
                corrected = XR_LSFlatFieldCorrection(raw{c}, readtiff(ffFns{c}), readtiff(bgFns{c}), ...
                    'LowerLimit', 0.4, 'constOffset', 1);
                expected = deskewRotateFrame3D(corrected, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                    'reverse', true, 'save16bit', true);
                [~, fsn] = fileparts(testCase.tileFilename(dataPath, chs{c}, 0));
                dsr = readtiff(sprintf('%sDSR/%s.tif', dataPath, fsn));
                testCase.verifyClass(dsr, 'uint16');
                testCase.verifyImagesClose(dsr, expected, 'relTol', 1e-6, 'absTol', 1, 'name', [chs{c}, ' DSR']);
                testCase.verifyTrue(exist(sprintf('%sDSR/MIPs/%s_MIP_z.tif', dataPath, fsn), 'file') == 2);
            end
        end

        function demoParameterSetIsAccepted(testCase)
            % the demos used to pass a 4-element Save16bit (for steps that were
            % removed from the pipeline in abd53f9), which the pipeline rejects
            demos = {'demo_microscope_pipeline_single_dataset', 'demo_microscope_pipeline_single_dataset_nostitch', ...
                'demo_microscope_pipeline_multiple_dataset'};
            save16bits = cell(1, numel(demos));
            for i = 1 : numel(demos)
                tokens = regexp(fileread(fullfile(pktest.repoRoot(), 'demos', [demos{i}, '.m'])), ...
                    '(?m)^\s*Save16bit\s*=\s*([^;\r\n]+);', 'tokens', 'once');
                testCase.assertNotEmpty(tokens, sprintf('%s no longer sets Save16bit.', demos{i}));
                save16bits{i} = eval(tokens{1});
            end

            dataPath = testCase.createTempFolder();
            writetiff(uint16(pktest.beads([64, 64, 40], [32, 30, 20], 2, 2000, 300)), testCase.tileFilename(dataPath, 'CamA_ch0', 0));
            for i = 1 : numel(demos)
                try
                    XR_microscopeAutomaticProcessing(dataPath, 'xyPixelSize', testCase.XYPixelSize, 'dz', testCase.Dz, ...
                        'Reverse', true, 'ChannelPatterns', {'CamA_ch0'}, 'Save16bit', save16bits{i}, ...
                        'Overwrite', false, 'Streaming', false, 'Deskew', true, 'Rotate', true, 'Stitch', false, ...
                        'DSRCombined', true, 'minModifyTime', 1, 'maxModifyTime', 10, 'maxWaitLoopNum', 10, ...
                        'parseCluster', false, 'cpusPerTask', 4, 'configFile', '', 'mccMode', false);
                catch ME
                    testCase.verifyFail(sprintf('The Save16bit of %s (%s) is rejected: %s', demos{i}, ...
                        mat2str(save16bits{i}), ME.message));
                end
            end
            testCase.verifyTrue(exist(strrep(testCase.tileFilename(dataPath, 'CamA_ch0', 0), dataPath, [dataPath, 'DSR/']), 'file') == 2);
        end

        function clusterModeWaitGrowsWithUnfinishedFiles(testCase)
            % with jobs on the cluster the pipeline always waited 30 s between
            % rounds, so the last results were seen up to 30 s late; it now
            % waits 0.3 s per unfinished file, at least 3 s and at most 30 s
            % (or until the results come, see clusterModeEndsSoonAfterTheLastResults)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            fns = arrayfun(@(y) testCase.tileFilename(dataPath, 'CamA_ch0', y), 0 : 19, 'unif', 0);
            for i = 1 : 20
                writetiff(zeros(16, 16, 8, 'uint16'), fns{i});
            end
            logPath = [dataPath, 'logs/'];
            mkdir(logPath);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            waitSpy = testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', [], 20));

            % the fake jobs never run: one round submits, the next resubmits
            args = {'xyPixelSize', testCase.XYPixelSize, 'dz', testCase.Dz, 'reverse', true, ...
                'channelPatterns', {'CamA_ch0'}, 'streaming', false, 'deskew', true, 'rotate', true, ...
                'stitch', false, 'DSRCombined', true, 'parseCluster', true, 'maxTrialNum', 2, 'jobLogDir', logPath};
            evalc('XR_microscopeAutomaticProcessing(dataPath, args{:})');

            testCase.verifyEqual(slurm.callCount('sbatch'), 40);
            calls = waitSpy.getCalls();
            waits = cellfun(@(x) x{2}, calls);
            % the second round also waits for the last tries' jobs to end
            testCase.verifyEqual(waits, [6, 6], 'AbsTol', 1e-9, 'The wait for 20 unfinished files should be 6 s.');
            % for their results
            awaited = regexprep(calls{1}{1}{1}(:), '/+', '/');
            testCase.verifyEqual(sort(awaited), sort(strrep(fns(:), dataPath, [dataPath, 'DSR/'])));
        end

        function clusterModeEndsSoonAfterTheLastResults(testCase)
            % the pipeline waited the whole wait between its rounds (0.3 s
            % per unfinished file, at least 3 s), so it saw the last
            % results up to that late; it now looks for the results still
            % to come each second and goes on once they are all there. The
            % fake jobs don't run; their results come during the second
            % second.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            fns = arrayfun(@(y) testCase.tileFilename(dataPath, 'CamA_ch0', y), 0 : 1, 'unif', 0);
            for i = 1 : 2
                writetiff(zeros(16, 16, 8, 'uint16'), fns{i});
            end
            logPath = [dataPath, 'logs/'];
            mkdir(logPath);
            dsrFns = strrep(fns, dataPath, [dataPath, 'DSR/']);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [dataPath, 'calls.txt'], 2, @() cellfun(@(x) fclose(fopen(x, 'w')), dsrFns)), 20));
            args = {'xyPixelSize', testCase.XYPixelSize, 'dz', testCase.Dz, 'reverse', true, ...
                'channelPatterns', {'CamA_ch0'}, 'streaming', false, 'deskew', true, 'rotate', true, ...
                'stitch', false, 'DSRCombined', true, 'parseCluster', true, 'maxTrialNum', 2, 'jobLogDir', logPath};
            evalc('XR_microscopeAutomaticProcessing(dataPath, args{:})');
            testCase.verifyEqual(slurm.callCount('sbatch'), 2, 'A job was submitted again.');
            testCase.verifyEqual(cellfun(@(x) x{1}, pauses.getCalls()), [1, 1], ...
                'The pipeline did not go on within a second of the results.');
        end

        function clusterModeDoesNotWaitForStepsGivenUp(testCase)
            % a step given up (the job of its last try ended without its
            % result) doesn't hold the wait for the results still to come:
            % tile 1's job is gone at the second round, tile 2's runs, and
            % tile 2's result comes during the first second of the second
            % wait, after which the pipeline goes on (and ends)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            fns = arrayfun(@(y) testCase.tileFilename(dataPath, 'CamA_ch0', y), 0 : 1, 'unif', 0);
            for i = 1 : 2
                writetiff(zeros(16, 16, 8, 'uint16'), fns{i});
            end
            logPath = [dataPath, 'logs/'];
            mkdir(logPath);
            dsrFns = strrep(fns, dataPath, [dataPath, 'DSR/']);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('1001_2', {'R'});
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [dataPath, 'calls.txt'], 4, @() fclose(fopen(dsrFns{2}, 'w'))), 20));
            args = {'xyPixelSize', testCase.XYPixelSize, 'dz', testCase.Dz, 'reverse', true, ...
                'channelPatterns', {'CamA_ch0'}, 'streaming', false, 'deskew', true, 'rotate', true, ...
                'stitch', false, 'DSRCombined', true, 'parseCluster', true, 'maxTrialNum', 1, 'jobLogDir', logPath};
            evalc('XR_microscopeAutomaticProcessing(dataPath, args{:})');
            testCase.verifyEqual(slurm.callCount('sbatch'), 2);
            testCase.verifyEqual(cellfun(@(x) x{1}, pauses.getCalls()), [1, 1, 1, 1], ...
                'The pipeline waited for the result of a step it gave up.');
        end

        function clusterModeSeesResultsOtherNodesJustMade(testCase)
            % on NFS, exist() can miss a file another node just made for 30
            % s or more, so the pipeline saw finished deskew results a round
            % late; it now checks them with batch_file_exist, which sees them
            % at once. The stub reports them made, so no job is submitted.
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            writetiff(zeros(16, 16, 8, 'uint16'), testCase.tileFilename(dataPath, 'CamA_ch0', 0));
            logPath = [dataPath, 'logs/'];
            mkdir(logPath);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.applyFixture(pktest.FunctionStubFixture('batch_file_exist', true));
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            args = {'xyPixelSize', testCase.XYPixelSize, 'dz', testCase.Dz, 'reverse', true, ...
                'channelPatterns', {'CamA_ch0'}, 'streaming', false, 'deskew', true, 'rotate', true, ...
                'stitch', false, 'DSRCombined', true, 'parseCluster', true, 'maxTrialNum', 2, 'jobLogDir', logPath};
            evalc('XR_microscopeAutomaticProcessing(dataPath, args{:})');
            testCase.verifyEqual(slurm.callCount('sbatch'), 0, 'A deskew job was submitted for a finished result.');
        end
    end

    methods (Test, TestTags = {'Integration', 'Slow'})
        function deskewRotateAndStitchMultipleDatasets(testCase)
            % demo_microscope_pipeline_multiple_dataset.m: two datasets, each
            % with its own image list; tiles overlap in y by 16 rows
            px = testCase.XYPixelSize;
            dataPaths = {testCase.createTempFolder(), testCase.createTempFolder()};
            csvFns = cell(1, 2);
            expected = cell(1, 2);
            for d = 1 : 2
                raw = uint16(pktest.beads([80, 64, 40], [20, 30, 15; 40, 34, 25; 64, 28, 20] + [d, 0, 0], 2, 3000, 200));
                expected{d} = deskewRotateFrame3D(raw, testCase.SkewAngle, testCase.Dz, px, 'reverse', true, 'save16bit', true);
                tiles = struct('im', {raw(1 : 48, :, :), raw(33 : 80, :, :)}, 'idx', {[0, 0, 0], [0, 1, 0]}, ...
                    'xyz', {[0, 0, 0], [0, 32 * px, 0]});
                [~, csvFns{d}] = pktest.writeStitchTiles(dataPaths{d}, tiles);
            end

            % edgeArtifacts 0: the default (2) erodes tile borders, which
            % would change the outer voxels of the stitched volume
            XR_microscopeAutomaticProcessing(dataPaths, 'xyPixelSize', px, 'dz', testCase.Dz, 'reverse', true, ...
                'channelPatterns', {'CamA_ch0'}, 'save16bit', [true, true], 'streaming', false, 'deskew', true, ...
                'rotate', true, 'stitch', true, 'DSRCombined', true, 'imageListFullpaths', csvFns, ...
                'axisOrder', '-x,y,z', 'edgeArtifacts', 0, 'parseCluster', false, 'cpusPerTask', 4);

            for d = 1 : 2
                stitchFn = sprintf('%smatlab_stitch/Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs.zarr', dataPaths{d});
                testCase.assertTrue(exist(stitchFn, 'dir') == 7, sprintf('Stitched result %s is missing.', stitchFn));
                out = readzarr(stitchFn);
                sz = size(expected{d});
                testCase.verifyEqual(size(out), sz + 2);
                testCase.verifyEqual(out(1 : sz(1), 1 : sz(2), 1 : sz(3)), expected{d}, ...
                    sprintf('Dataset %d: stitched DSR tiles differ from the DSR of the whole volume.', d));
            end
        end

        function localStitchingFinishesWithoutWaiting(testCase)
            % the stitching step looked for its result before computing it,
            % so a local run always idled 30 s before seeing that it was done
            px = testCase.XYPixelSize;
            dataPath = testCase.createTempFolder();
            raw = uint16(pktest.beads([80, 64, 40], [20, 30, 15; 64, 28, 20], 2, 3000, 200));
            tiles = struct('im', {raw(1 : 48, :, :), raw(33 : 80, :, :)}, 'idx', {[0, 0, 0], [0, 1, 0]}, ...
                'xyz', {[0, 0, 0], [0, 32 * px, 0]});
            [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles);
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause'));

            XR_microscopeAutomaticProcessing(dataPath, 'xyPixelSize', px, 'dz', testCase.Dz, 'reverse', true, ...
                'channelPatterns', {'CamA_ch0'}, 'save16bit', [true, true], 'streaming', false, 'deskew', true, ...
                'rotate', true, 'stitch', true, 'DSRCombined', true, 'imageListFullpaths', {csvFn}, ...
                'axisOrder', '-x,y,z', 'parseCluster', false, 'cpusPerTask', 4);

            stitchFn = [dataPath, 'matlab_stitch/Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs.zarr'];
            testCase.verifyTrue(exist(stitchFn, 'dir') == 7, 'The stitched result is missing.');
            waits = cellfun(@(x) x{1}, pauses.getCalls());
            testCase.verifyFalse(any(waits >= 30), sprintf('The pipeline waited %s s.', mat2str(waits)));
        end

        function clusterModeWaitsForTheLastStitchingTry(testCase)
            % after submitting the last try of a dataset's stitching job, the
            % pipeline returned without waiting for it; it now checks the job
            % until it ends, without submitting again
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            px = testCase.XYPixelSize;
            dataPath = testCase.createTempFolder();
            raw = uint16(pktest.beads([80, 64, 40], [20, 30, 15; 64, 28, 20], 2, 3000, 200));
            tiles = struct('im', {raw(1 : 48, :, :), raw(33 : 80, :, :)}, 'idx', {[0, 0, 0], [0, 1, 0]}, ...
                'xyz', {[0, 0, 0], [0, 32 * px, 0]});
            [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            % jobs in submission order: deskew of tile 1 (1000) and tile 2
            % (1001), stitching of the dataset (1002); the deskew jobs end at
            % once without results, the stitching job runs for longer, so
            % that only it keeps the run going. Each file of the dataset
            % checks it once per round.
            slurm.setJobStates('1002_1', {'PD', 'R', 'R', 'R', 'done'});
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            evalc(['XR_microscopeAutomaticProcessing(dataPath, ''xyPixelSize'', px, ''dz'', testCase.Dz, ', ...
                '''reverse'', true, ''channelPatterns'', {''CamA_ch0''}, ''streaming'', false, ', ...
                '''stitch'', true, ''imageListFullpaths'', {csvFn}, ''parseCluster'', true, ', ...
                '''maxTrialNum'', 1, ''jobLogDir'', [dataPath, ''logs/''])']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 3, 'A last try was submitted again.');
            checks = slurm.history();
            checks = checks(startsWith(checks, 'squeue') & contains(checks, '1002_1'));
            testCase.verifyGreaterThanOrEqual(numel(checks), 5, 'The stitching job was not checked until it ended.');
        end

        function clusterModeSubmitsStitchingAfterTheDeskewJobs(testCase)
            % a dataset's stitching job was submitted right after the deskew
            % job of its first file, so on a busy cluster it took a slot
            % ahead of the dataset's other deskew jobs, only to wait for
            % their results; it is now submitted after them
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            [dataPath, csvFn] = testCase.writeTwoTiles();
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            args = testCase.clusterStitchArgs(dataPath, csvFn, 1);
            evalc('XR_microscopeAutomaticProcessing(dataPath, args{:})');
            jobs = slurm.calls('sbatch');
            testCase.assertNumElements(jobs, 3);
            testCase.verifyEqual(find(contains(jobs, 'XR_matlab_stitching_wrapper')), 3, ...
                'The stitching job was not submitted after the deskew jobs.');
        end

        function clusterModeSubmitsStitchingWhenTilesWereDeskewedBefore(testCase)
            % (with the stitching job after the deskew jobs) the deskewed
            % tile of a file done before needs no job: the stitching job
            % follows the deskew job of the other file
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            [dataPath, csvFn] = testCase.writeTwoTiles();
            [~, fsn] = fileparts(testCase.tileFilename(dataPath, 'CamA_ch0', 0));
            mkdir([dataPath, 'DSR']);
            fclose(fopen(sprintf('%sDSR/%s.tif', dataPath, fsn), 'w'));
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            args = testCase.clusterStitchArgs(dataPath, csvFn, 1);
            evalc('XR_microscopeAutomaticProcessing(dataPath, args{:})');
            jobs = slurm.calls('sbatch');
            testCase.assertNumElements(jobs, 2);
            testCase.verifyTrue(contains(jobs{2}, 'XR_matlab_stitching_wrapper'), ...
                'The stitching job was not submitted after the deskew job.');
        end

        function clusterModeLooksForTheStitchedResultEachSecond(testCase)
            % (as clusterModeEndsSoonAfterTheLastResults, for the stitched
            % result, which the pipeline finds by a pattern) the tiles were
            % deskewed before; the stitching job doesn't run, its result
            % comes during the third second, after its temporary folder
            % (which doesn't count) during the second
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            [dataPath, csvFn] = testCase.writeTwoTiles();
            mkdir([dataPath, 'DSR']);
            for y = 0 : 1
                [~, fsn] = fileparts(testCase.tileFilename(dataPath, 'CamA_ch0', y));
                fclose(fopen(sprintf('%sDSR/%s.tif', dataPath, fsn), 'w'));
            end
            stitchFn = [dataPath, 'matlab_stitch/Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs.zarr'];
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() pktest.atCall( ...
                [dataPath, 'calls.txt'], [2, 3], {@() mkdir([stitchFn, '_1234']), @() mkdir(stitchFn)}), 20));
            args = testCase.clusterStitchArgs(dataPath, csvFn, 1);
            evalc('XR_microscopeAutomaticProcessing(dataPath, args{:})');
            testCase.verifyEqual(slurm.callCount('sbatch'), 1);
            testCase.verifyEqual(cellfun(@(x) x{1}, pauses.getCalls()), [1, 1, 1], ...
                'The pipeline did not go on within a second of the stitched result.');
        end

        function clusterModeTellsStitchingThatTilesAreComing(testCase)
            % the stitching job gave up a time point whose deskewed tiles
            % were not all there after 3 checks 5 s apart, and the pipeline
            % then needed another stitching job for it (or lost it, after
            % the tries of its stitching, while deskewing took long). While
            % deskew jobs of a dataset are pending or running, a marker in
            % its DSR folder now tells the stitching job that the tiles are
            % coming; it is written again once it is a minute old, so that
            % it stays recent, and removed when they are done (at the start
            % of a round). The fake deskew jobs never run: their two tries
            % end in the third round, while the stitching job (1002) runs
            % for a round more (each file checks it once per round).
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.assumeTrue(usejava('jvm'), 'Setting file times needs Java.');
            [dataPath, csvFn] = testCase.writeTwoTiles();
            markerFn = [dataPath, 'DSR/deskew_in_progress.tmp'];
            logFn = [dataPath, 'marker_ages.txt'];
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('1002_1', {'R', 'R', 'R', 'R', 'R', 'R', 'done'});
            % at each wait: the marker's age, then 5 minutes added to it
            testCase.applyFixture(pktest.FunctionSpyFixture('wait_for_files', ...
                @() testCase.recordAndAgeFile(markerFn, logFn, false), 20));
            args = testCase.clusterStitchArgs(dataPath, csvFn, 2);
            evalc('XR_microscopeAutomaticProcessing(dataPath, args{:})');
            ages = str2double(splitlines(strip(fileread(logFn))));
            testCase.assertGreaterThanOrEqual(numel(ages), 4);
            testCase.verifyTrue(all(ages(1 : 3) < 60), sprintf(['The marker was missing (NaN) or ', ...
                'not written again (age in s): %s.'], mat2str(ages')));
            testCase.verifyTrue(all(isnan(ages(4 : end))), sprintf(['The marker was kept after the deskew ', ...
                'tries (age in s): %s.'], mat2str(ages')));
            testCase.verifyFalse(isfile(markerFn), 'The marker was left after the run.');
        end

        function deskewMarkerIsRemovedWhenThePipelineStops(testCase)
            % an error (or Ctrl-C) in the pipeline must not leave its marker
            % for a stitching job to wait for
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            testCase.assumeTrue(usejava('jvm'), 'Setting file times needs Java.');
            [dataPath, csvFn] = testCase.writeTwoTiles();
            markerFn = [dataPath, 'DSR/deskew_in_progress.tmp'];
            logFn = [dataPath, 'marker_ages.txt'];
            testCase.applyFixture(pktest.FakeSlurmFixture);
            % the first wait stops the run
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', ...
                @() testCase.recordAndAgeFile(markerFn, logFn, true), 20));
            args = testCase.clusterStitchArgs(dataPath, csvFn, 2);
            stopped = false;
            try
                evalc('XR_microscopeAutomaticProcessing(dataPath, args{:})');
            catch ME
                stopped = strcmp(ME.identifier, 'pktest:stop');
            end
            testCase.assertTrue(stopped, 'The pipeline did not stop at its first wait.');
            testCase.verifyLessThan(str2double(strip(fileread(logFn))), 60, ...
                'There was no marker while the deskew jobs were pending.');
            testCase.verifyFalse(isfile(markerFn), 'The marker was left after the pipeline stopped.');
        end

        function streamingLeavesOutRecentFilesWithTheirSizes(testCase)
            % in streaming mode, files changed in the last minModifyTime
            % minutes (still being written) are left for later, but their
            % sizes stayed in the list: a cluster run that started with
            % only such files failed (3 sizes for no file), and the files
            % after such a file got the sizes of the files before them
            testCase.assumeTrue(usejava('jvm'), 'Setting file times needs Java.');
            dataPath = testCase.createTempFolder();
            fns = arrayfun(@(y) testCase.tileFilename(dataPath, 'CamA_ch0', y), 0 : 2, 'unif', 0);
            depths = [16, 20, 24];
            for i = 1 : 3
                writetiff(zeros(16, 16, depths(i), 'uint16'), fns{i});
            end
            args = {{dataPath}, {'CamA_ch0'}, false, false, false, {}};
            [~, ~, ~, ~, allSizes] = XR_parseImageFilenames(args{:}, false, 1, false);

            % all files new: none of them
            [fnames, ~, gfnames, partialvols, dataSizes] = XR_parseImageFilenames(args{:}, true, 1, false);
            testCase.verifyEmpty(fnames);
            testCase.verifyEmpty(gfnames);
            testCase.verifyEmpty(partialvols);
            testCase.verifyEmpty(dataSizes, 'The sizes of the files left out stayed in the list.');

            % the last file 5 minutes old, the others new: only it, with its size
            java.io.File(fns{3}).setLastModified(java.lang.System.currentTimeMillis() - 5 * 60 * 1000);
            [fnames, ~, gfnames, partialvols, dataSizes] = XR_parseImageFilenames(args{:}, true, 1, false);
            [~, fsn, ext] = fileparts(fns{3});
            testCase.verifyEqual(fnames, {[fsn, ext]});
            testCase.verifyEqual(gfnames, {{[fsn, ext]}}, 'The file got the group of another file.');
            testCase.verifyNumElements(partialvols, 1);
            testCase.verifyEqual(dataSizes, allSizes(3), 'The file got the size of another file.');
        end

        function streamingProcessesNewFiles(testCase)
            % with streaming on, files that arrive during the run used to
            % make the pipeline fail (it added 4 flag columns to 2)
            testCase.assumeTrue(usejava('jvm'), 'Setting file times needs Java.');
            px = testCase.XYPixelSize;
            dataPath = testCase.createTempFolder();
            stagePath = testCase.createTempFolder();
            raw = uint16(pktest.beads([48, 64, 40], [20, 30, 15], 2, 3000, 200));
            writetiff(raw, testCase.tileFilename(dataPath, 'CamA_ch0', 0));
            newFn = testCase.tileFilename(dataPath, 'CamA_ch0', 1);
            stagedFn = testCase.tileFilename(stagePath, 'CamA_ch0', 1);
            writetiff(raw, stagedFn);
            % files changed in the last minute count as still being written
            oldTime = java.lang.System.currentTimeMillis() - 5 * 60 * 1000;
            java.io.File(testCase.tileFilename(dataPath, 'CamA_ch0', 0)).setLastModified(oldTime);
            java.io.File(stagedFn).setLastModified(oldTime);
            % the second file arrives (with its old time) at the first wait
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', @() testCase.moveIfPresent(stagedFn, newFn)));
            evalc(['XR_microscopeAutomaticProcessing(dataPath, ''xyPixelSize'', px, ''dz'', testCase.Dz, ', ...
                '''reverse'', true, ''channelPatterns'', {''CamA_ch0''}, ''streaming'', true, ''deskew'', true, ', ...
                '''rotate'', true, ''DSRCombined'', true, ''stitch'', false, ''parseCluster'', false, ', ...
                '''maxModifyTime'', 1, ''maxWaitLoopNum'', 2, ''cpusPerTask'', 4)']);
            testCase.verifyTrue(isfile(strrep(newFn, dataPath, [dataPath, 'DSR/'])), 'The new file was not processed.');
            testCase.verifyTrue(isfile(strrep(testCase.tileFilename(dataPath, 'CamA_ch0', 0), dataPath, [dataPath, 'DSR/'])));
        end

        function localRunSubmitsNoStitchingJobs(testCase)
            % the stitching task didn't get parseCluster, so a local run on a
            % machine with Slurm stitched in Slurm jobs
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            px = testCase.XYPixelSize;
            dataPath = testCase.createTempFolder();
            raw = uint16(pktest.beads([80, 64, 40], [20, 30, 15; 64, 28, 20], 2, 3000, 200));
            tiles = struct('im', {raw(1 : 48, :, :), raw(33 : 80, :, :)}, 'idx', {[0, 0, 0], [0, 1, 0]}, ...
                'xyz', {[0, 0, 0], [0, 32 * px, 0]});
            [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);

            XR_microscopeAutomaticProcessing(dataPath, 'xyPixelSize', px, 'dz', testCase.Dz, 'reverse', true, ...
                'channelPatterns', {'CamA_ch0'}, 'save16bit', [true, true], 'streaming', false, 'deskew', true, ...
                'rotate', true, 'stitch', true, 'DSRCombined', true, 'imageListFullpaths', {csvFn}, ...
                'axisOrder', '-x,y,z', 'parseCluster', false, 'cpusPerTask', 4);

            testCase.verifyEqual(slurm.callCount('sbatch'), 0, 'The local run submitted Slurm jobs.');
            stitchFn = [dataPath, 'matlab_stitch/Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs.zarr'];
            testCase.verifyTrue(exist(stitchFn, 'dir') == 7, 'The stitched result is missing.');
        end
    end

    methods (Access = private)
        function [fn] = tileFilename(~, dataPath, ch, yi)
            fn = sprintf(['%sScan_Iter_0000_0000_%s_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs_', ...
                '000x_%03dy_000z_0000t.tif'], dataPath, ch, yi);
        end

        function [] = moveIfPresent(~, src, dst)
            % move a file (keeping its time) unless it was moved already
            if isfile(src)
                movefile(src, dst);
            end
        end

        function [dataPath, csvFn] = writeTwoTiles(testCase)
            % two tiles of one volume, overlapping in y, and their image list
            px = testCase.XYPixelSize;
            dataPath = testCase.createTempFolder();
            raw = uint16(pktest.beads([80, 64, 40], [20, 30, 15; 64, 28, 20], 2, 3000, 200));
            tiles = struct('im', {raw(1 : 48, :, :), raw(33 : 80, :, :)}, 'idx', {[0, 0, 0], [0, 1, 0]}, ...
                'xyz', {[0, 0, 0], [0, 32 * px, 0]});
            [~, csvFn] = pktest.writeStitchTiles(dataPath, tiles);
        end

        function [args] = clusterStitchArgs(testCase, dataPath, csvFn, maxTrialNum)
            % deskew/rotation and stitching with (fake) Slurm jobs
            args = {'xyPixelSize', testCase.XYPixelSize, 'dz', testCase.Dz, 'reverse', true, ...
                'channelPatterns', {'CamA_ch0'}, 'streaming', false, 'stitch', true, ...
                'imageListFullpaths', {csvFn}, 'parseCluster', true, 'maxTrialNum', maxTrialNum, ...
                'jobLogDir', [dataPath, 'logs/']};
        end

        function [] = recordAndAgeFile(~, fn, logFn, stop)
            % append the file's age in s (NaN: no file) to logFn, then set
            % its time 5 minutes back; stop: then error
            age = NaN;
            dir_info = dir(fn);
            if ~isempty(dir_info)
                age = (now - dir_info(1).datenum) * 24 * 3600;
                java.io.File(fn).setLastModified(java.lang.System.currentTimeMillis() - 5 * 60 * 1000);
            end
            fid = fopen(logFn, 'a');
            fprintf(fid, '%.3f\n', age);
            fclose(fid);
            if stop
                error('pktest:stop', 'Stopped by the test.');
            end
        end
    end
end
