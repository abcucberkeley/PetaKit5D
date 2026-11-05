classdef NoSlurmFixture < matlab.unittest.fixtures.Fixture
    % hide Slurm from the code under test by putting failing sinfo, sbatch,
    % squeue and scancel stubs first on the PATH. A wrapper that doesn't
    % forward parseCluster to the job framework (XR_zarrToTiff_wrapper
    % didn't) would submit real jobs on a machine with Slurm. With the
    % stubs, checkSlurmCluster sees no cluster and everything runs locally.
    % Windows has no Slurm, so nothing is changed there. The stub folder is
    % in PETAKIT5D_TEST_SLURM_STUBS so pktest.RealSlurmFixture can lift it
    % for the opt-in cluster tests.
    %
    % The Slurm job of the tests (when they run in one, e.g. with srun or
    % sbatch) is hidden too: its environment variables are unset while the
    % fixture is set up, so the code under test acts as outside a job
    % wherever the tests run (a master in a Slurm job computes only the tasks
    % that fit the job's memory, so in a job of 16 GB the tests of tasks the
    % master computes failed). Tests that need a job set these themselves.

    properties (Constant, Access = private)
        % the variables of a job that PetaKit5D and its scripts read
        JobVariables = {'SLURM_JOB_ID', 'SLURM_JOBID', 'SLURM_ARRAY_JOB_ID', 'SLURM_ARRAY_TASK_ID', ...
            'SLURMD_NODENAME', 'SLURM_MEM_PER_NODE', 'SLURM_MEM_PER_CPU', 'SLURM_CPUS_ON_NODE', ...
            'SLURM_JOB_START_TIME', 'SLURM_JOB_END_TIME'}
    end
    properties (Access = private)
        OriginalPath
        OriginalStubs
        OriginalJobValues
        StubFolder
    end

    methods
        function setup(fixture)
            if ispc
                fixture.SetupDescription = 'No Slurm stubs needed on Windows';
                return;
            end
            fixture.StubFolder = tempname;
            mkdir(fixture.StubFolder);
            cmds = {'sinfo', 'sbatch', 'squeue', 'scancel'};
            for i = 1 : numel(cmds)
                fn = [fixture.StubFolder, '/', cmds{i}];
                fid = fopen(fn, 'w');
                fprintf(fid, '#!/bin/sh\nexit 1\n');
                fclose(fid);
                fileattrib(fn, '+x');
            end
            fixture.OriginalPath = getenv('PATH');
            % (a test that runs tests of its own sets the fixture up again)
            fixture.OriginalStubs = getenv('PETAKIT5D_TEST_SLURM_STUBS');
            setenv('PATH', [fixture.StubFolder, pathsep, fixture.OriginalPath]);
            setenv('PETAKIT5D_TEST_SLURM_STUBS', fixture.StubFolder);
            fixture.OriginalJobValues = cellfun(@getenv, fixture.JobVariables, 'unif', 0);
            cellfun(@unsetenv, fixture.JobVariables);
            fixture.SetupDescription = 'Hid Slurm commands behind failing stubs, and the Slurm job of the tests';
        end

        function teardown(fixture)
            if ispc
                return;
            end
            setenv('PATH', fixture.OriginalPath);
            setenv('PETAKIT5D_TEST_SLURM_STUBS', fixture.OriginalStubs);
            for i = 1 : numel(fixture.JobVariables)
                if isempty(fixture.OriginalJobValues{i})
                    unsetenv(fixture.JobVariables{i});
                else
                    setenv(fixture.JobVariables{i}, fixture.OriginalJobValues{i});
                end
            end
            if exist(fixture.StubFolder, 'dir')
                rmdir(fixture.StubFolder, 's');
            end
        end
    end
end
