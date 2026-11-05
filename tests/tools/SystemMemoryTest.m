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
            % the limit of the job's own cgroup. In a Slurm job only.
            testCase.assumeTrue(isunix && ~ismac && ~isempty(getenv('SLURM_JOB_ID')), 'Not in a Slurm job.');
            txt = fileread('/proc/self/cgroup');
            testCase.assumeNotEmpty(regexp(txt, ['/job_', getenv('SLURM_JOB_ID'), '(/|\s)'], 'once'), ...
                'The job has no cgroup of its own.');
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
