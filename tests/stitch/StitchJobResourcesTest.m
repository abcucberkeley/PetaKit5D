classdef StitchJobResourcesTest < pktest.TestCase
    % a stitching job of the cluster is sized for what it holds (see
    % XR_matlab_stitching_wrapper), so the stitching function decides with
    % its job's memory: it stitches a small volume in its own process only if
    % its blocks and distance maps fit the job (stitch_in_process), and it
    % measures the volume for its MIPs against at least the largest job of
    % the config (stitch_mip_memory), as its jobs used to have memory for all
    % their tiles.

    methods (Test, TestTags = {'Unit'})
        function smallVolumesAreStitchedInTheirJobIfTheyFitIt(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_JOB_ID', '4242'));
            % 512^3 batches: 5 GB (20 B per voxel of uint16); tiles of
            % 1024^3: 40 GB distance maps
            batch = [512, 512, 512];
            tiles = repmat([1024, 1024, 1024], 4, 1);
            testCase.verifyFalse(stitch_in_process(30, batch, 'uint16', 'none', tiles, false, ''), ...
                'A volume of 30 batches was stitched in the job.');
            testCase.verifyTrue(stitch_in_process(29, batch, 'uint16', 'feather', tiles, false, ''), ...
                'Without the cluster, a small volume was not stitched here.');
            stub = testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 8)); %#ok<NASGU>
            testCase.verifyTrue(stitch_in_process(29, batch, 'uint16', 'none', tiles, true, ''));
            testCase.verifyFalse(stitch_in_process(29, batch, 'uint16', 'feather', tiles, true, ''), ...
                'Distance maps that do not fit the job were made in it.');
            testCase.verifyFalse(stitch_in_process(29, batch * 2, 'uint16', 'none', tiles, true, ''), ...
                'Blocks that do not fit the job were stitched in it.');
            % (outside Slurm: as before)
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_JOB_ID', ''));
            testCase.verifyTrue(stitch_in_process(29, batch * 2, 'uint16', 'feather', tiles, true, ''));
        end

        function mipMemoryOfAClusterJobIsAtLeastTheLargestJob(testCase)
            folder = testCase.createTempFolder();
            configFn = [folder, 'config.json'];
            fid = fopen(configFn, 'w');
            fprintf(fid, '{"maxCPUNum": 32, "memPerCPU": 7.6}');
            fclose(fid);
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 8));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_JOB_ID', '4242'));
            testCase.verifyEqual(stitch_mip_memory(true, configFn), 32 * 7.6, 'AbsTol', 1e-9);
            % (the job framework's defaults without a config)
            testCase.verifyEqual(stitch_mip_memory(true, ''), 24 * 20.9, 'AbsTol', 1e-9);
            % not in cluster mode, or outside Slurm: this process's memory
            testCase.verifyEqual(stitch_mip_memory(false, configFn), 8);
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_JOB_ID', ''));
            testCase.verifyEqual(stitch_mip_memory(true, configFn), 8);
        end
    end
end
