classdef GenericComputingFrameworkDemoTest < pktest.TestCase
    % the workflow of demos/demo_generic_computing_framework.m on small
    % synthetic data: a user-defined function (crop + deskew/rotation) run
    % through the generic computing framework, and the slurm configuration
    % file written by the demo

    properties (Constant)
        SkewAngle = 32.45
        XYPixelSize = 0.108
        Dz = 0.3
        InBbox = [9, 1, 1, 56, 48, 30]
    end

    methods (Test, TestTags = {'Integration'})
        function userFunctionRunsThroughFramework(testCase)
            dataPath = testCase.createTempFolder();
            outPath = [dataPath, 'Cropped_DSR/'];
            mkdir(outPath);
            fsns = {'Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t.tif', ...
                'Scan_Iter_0000_0000_CamB_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t.tif'};
            inputFullnames = cellfun(@(x) [dataPath, x], fsns, 'unif', 0);
            outputFullnames = cellfun(@(x) [outPath, x], fsns, 'unif', 0);
            for f = 1 : 2
                writetiff(uint16(pktest.texture([64, 48, 32], 2)), inputFullnames{f});
            end

            % task strings built exactly as in the demo (%d prints 32.45 as 3.245000e+01)
            funcStrs = cell(2, 1);
            for f = 1 : 2
                funcStrs{f} = sprintf(['crop_deskew_rotate_demo_function(''%s'',''%s'',%d,', ...
                    '%d,%d,%s)'], inputFullnames{f}, outputFullnames{f}, testCase.SkewAngle, ...
                    testCase.XYPixelSize, testCase.Dz, strrep(mat2str(testCase.InBbox), ' ', ','));
            end
            generic_computing_frameworks_wrapper(inputFullnames, outputFullnames, funcStrs, ...
                parseCluster=false, masterCompute=true, memAllocate=20, configFile='', mccMode=false);

            for f = 1 : 2
                im = crop3d_mex(readtiff(inputFullnames{f}), testCase.InBbox);
                expected = deskewRotateFrame3D(im, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                    'reverse', true, 'objectiveScan', false, 'interpMethod', 'linear', 'save16bit', true);
                testCase.verifyEqual(readtiff(outputFullnames{f}), expected);
            end
            testCase.verifyEmpty(dir([outPath, '*.tmp']));
        end

        function demoConfigFileRoundTripsAndRunsLocally(testCase)
            % the slurm configuration from the demo, written the same way
            dataPath = testCase.createTempFolder();
            a = struct();
            a.MCRParam = '/global/home/users/user/bin/MATLAB_Runtime/R2024b';
            a.MCCMasterStr = '/clusterfs/fiona/user/PetaKit5D/mcc/linux/run_mccMaster.sh';
            a.BashLaunchStr = 'module load parallel/20190822';
            a.SlurmParam = '-p abc --qos abc_normal -n1 --mem-per-cpu=21418M';
            a.jobLogDir = '/clusterfs/fiona/user/Projects/job_logs';
            a.maxCPUNum = 24;
            a.minCPUNum = 1;
            a.wholeNodeJob = true;
            a.memPerCPU = 20.9;
            a.GNUparallel = true;
            a.jobTimeLimit = 12;
            a.maxJobNum = 5000;
            s = jsonencode(a, PrettyPrint=true);
            configFile = [dataPath, 'test_slurm_cpu_config_file.json'];
            fid = fopen(configFile, 'w');
            fprintf(fid, s);
            fclose(fid);

            testCase.verifyEqual(jsondecode(fileread(configFile)), a);

            % with the config file and no cluster, tasks still run locally
            clear generic_computing_frameworks_wrapper
            inFn = [dataPath, 'in.txt'];
            outFns = {[dataPath, 'out_1.txt'], [dataPath, 'out_2.txt']};
            fclose(fopen(inFn, 'w'));
            logFn = [dataPath, 'calls.log'];
            funcStrs = cellfun(@(x) sprintf('pktest.frameworkTask(''%s'',''%s'',''%s'',''copy'')', inFn, x, logFn), ...
                outFns, 'unif', 0);
            flags = generic_computing_frameworks_wrapper({inFn, inFn}, outFns, funcStrs, ...
                parseCluster=false, masterCompute=true, configFile=configFile, mccMode=false);
            testCase.verifyEqual(flags, true(2, 1));
        end
    end
end
