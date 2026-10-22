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
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            dataPath = testCase.createTempFolder();
            for y = 0 : 19
                writetiff(zeros(16, 16, 8, 'uint16'), testCase.tileFilename(dataPath, 'CamA_ch0', y));
            end
            logPath = [dataPath, 'logs/'];
            mkdir(logPath);
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            pauses = testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));

            % the fake jobs never run: one round submits, the next resubmits
            args = {'xyPixelSize', testCase.XYPixelSize, 'dz', testCase.Dz, 'reverse', true, ...
                'channelPatterns', {'CamA_ch0'}, 'streaming', false, 'deskew', true, 'rotate', true, ...
                'stitch', false, 'DSRCombined', true, 'parseCluster', true, 'maxTrialNum', 2, 'jobLogDir', logPath};
            evalc('XR_microscopeAutomaticProcessing(dataPath, args{:})');

            testCase.verifyEqual(slurm.callCount('sbatch'), 40);
            waits = cellfun(@(x) x{1}, pauses.getCalls());
            % the second round also waits for the last tries' jobs to end
            testCase.verifyEqual(waits, [6, 6], 'AbsTol', 1e-9, 'The wait for 20 unfinished files should be 6 s.');
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
            % jobs in submission order: deskew of tile 1 (1000), stitching of
            % the dataset (1001), deskew of tile 2 (1002); the deskew jobs end
            % at once without results, the stitching job runs for longer, so
            % that only it keeps the run going. Each file of the dataset
            % checks it once per round.
            slurm.setJobStates('1001_1', {'PD', 'R', 'R', 'R', 'done'});
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', [], 20));
            evalc(['XR_microscopeAutomaticProcessing(dataPath, ''xyPixelSize'', px, ''dz'', testCase.Dz, ', ...
                '''reverse'', true, ''channelPatterns'', {''CamA_ch0''}, ''streaming'', false, ', ...
                '''stitch'', true, ''imageListFullpaths'', {csvFn}, ''parseCluster'', true, ', ...
                '''maxTrialNum'', 1, ''jobLogDir'', [dataPath, ''logs/''])']);
            testCase.verifyEqual(slurm.callCount('sbatch'), 3, 'A last try was submitted again.');
            checks = slurm.history();
            checks = checks(startsWith(checks, 'squeue') & contains(checks, '1001_1'));
            testCase.verifyGreaterThanOrEqual(numel(checks), 5, 'The stitching job was not checked until it ended.');
        end

        function streamingProcessesNewFiles(testCase)
            % with streaming on, files that arrive during the run used to
            % make the pipeline fail (it added 4 flag columns to 2)
            px = testCase.XYPixelSize;
            dataPath = testCase.createTempFolder();
            stagePath = testCase.createTempFolder();
            raw = uint16(pktest.beads([48, 64, 40], [20, 30, 15], 2, 3000, 200));
            writetiff(raw, testCase.tileFilename(dataPath, 'CamA_ch0', 0));
            newFn = testCase.tileFilename(dataPath, 'CamA_ch0', 1);
            stagedFn = testCase.tileFilename(stagePath, 'CamA_ch0', 1);
            writetiff(raw, stagedFn);
            % files changed in the last minute count as still being written
            system(sprintf('touch -d "-5 min" %s %s', testCase.tileFilename(dataPath, 'CamA_ch0', 0), stagedFn));
            % the second file arrives (with its old time) at the first wait
            testCase.applyFixture(pktest.FunctionSpyFixture('pause', ...
                @() system(sprintf('mv -n %s %s 2>/dev/null; true', stagedFn, newFn))));
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
    end
end
