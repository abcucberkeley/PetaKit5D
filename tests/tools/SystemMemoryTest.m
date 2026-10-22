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
            testCase.verifyEqual(getSystemMemory(), 32);
        end

        function memoryIsANumber(testCase)
            totalMem = getSystemMemory();
            testCase.verifyClass(totalMem, 'double');
            testCase.verifySize(totalMem, [1, 1]);
            testCase.verifyGreaterThan(totalMem, 0);
        end
    end
end
