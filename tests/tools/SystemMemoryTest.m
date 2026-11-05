classdef SystemMemoryTest < pktest.TestCase
    % getSystemMemory, which stitching uses to choose how to make the MIP of
    % the result

    methods (Test, TestTags = {'Unit'})
        function slurmJobMemoryIsInGB(testCase)
            % in a job started with --mem, the memory was returned as the
            % text of SLURM_MEM_PER_NODE (MB)
            testCase.assumeTrue(isunix && ~ismac, 'Slurm memory is only read on Linux.');
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_JOB_ID', '123'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_MEM_PER_NODE', '32768'));
            % (a job started with --mem has no SLURM_MEM_PER_CPU of its own)
            perCPU = getenv('SLURM_MEM_PER_CPU');
            if ~isempty(perCPU)
                testCase.addTeardown(@setenv, 'SLURM_MEM_PER_CPU', perCPU);
                unsetenv('SLURM_MEM_PER_CPU');
            end
            testCase.verifyEqual(getSystemMemory(), 32);
        end

        function memoryIsNotTheValueOfAParentJob(testCase)
            % a job submitted from another job (as the jobs of a job that
            % coordinates jobs are) inherits the other's SLURM_MEM_PER_NODE
            % (a 32-CPU controller got 16 GB, its master's): the memory is
            % the limit of the job's own cgroup. In a Slurm job only: the job
            % of the tests, from its cgroup (pktest.NoSlurmFixture hides its
            % variables)
            testCase.assumeTrue(isunix && ~ismac, 'Slurm memory is only read on Linux.');
            job = regexp(fileread('/proc/self/cgroup'), '/job_(\d+)(/|\s)', 'tokens', 'once');
            testCase.assumeNotEmpty(job, 'Not in a Slurm job with a cgroup of its own.');
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_JOB_ID', job{1}));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_MEM_PER_NODE', '1'));
            testCase.verifyGreaterThan(getSystemMemory(), 0.5);
        end

        function memoryPerNodeNextToMemoryPerCPUIsNotTaken(testCase)
            % without the job's cgroup (here a job id that isn't this
            % process's), SLURM_MEM_PER_NODE next to SLURM_MEM_PER_CPU is
            % one job's value in another's environment, so scontrol tells
            % (it knows no job 123: unknown)
            testCase.assumeTrue(isunix && ~ismac, 'Slurm memory is only read on Linux.');
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_JOB_ID', '123'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_MEM_PER_NODE', '32768'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_MEM_PER_CPU', '1024'));
            mem = [];
            evalc('mem = getSystemMemory();');
            testCase.verifyNotEqual(mem, 32);
        end

        function memoryIsANumber(testCase)
            totalMem = getSystemMemory();
            testCase.verifyClass(totalMem, 'double');
            testCase.verifySize(totalMem, [1, 1]);
            testCase.verifyGreaterThan(totalMem, 0);
        end
    end
end
