classdef SingleJobSubmitWrapperTest < pktest.TestCase
    % generic_single_job_submit_wrapper builds the sbatch command for one
    % task. pktest.NoSlurmFixture (shared by all tests) makes squeue and
    % sbatch fail, so nothing is submitted; the command is read from the
    % console output.

    properties (TestParameter)
        % name-value arguments and the expected --cpus-per-task
        cpuCase = struct( ...
            'cpusAboveMax', {{{'cpusPerTask', 40}, 24}}, ...
            'minCPUNum', {{{'minCPUNum', 4}, 4}}, ...
            'memAllocate', {{{'memAllocate', 100}, 5}}, ...
            'memAllocateAboveMax', {{{'memAllocate', 1000}, 24}})
    end

    methods (Test, TestTags = {'Unit'})
        function cpusPerTaskIsClampedToLimits(testCase, cpuCase)
            % fb6a98b regression: user-defined cpusPerTask/minCPUNum were
            % only clamped to [minCPUNum, maxCPUNum] when memAllocate was set
            cmd = testCase.submitCommand(cpuCase{1}{:});
            cpus = regexp(cmd, '--cpus-per-task=(\d+)', 'tokens', 'once');
            testCase.assertNotEmpty(cpus, 'No sbatch command was printed.');
            testCase.verifyEqual(str2double(cpus{1}), cpuCase{2});
        end

        function commandWrapsTheTaskString(testCase)
            cmd = testCase.submitCommand('SlurmParam', '-p test_partition -n1');
            testCase.verifySubstring(cmd, 'sbatch --array=7 ');
            testCase.verifySubstring(cmd, '-p test_partition -n1');
            testCase.verifySubstring(cmd, 'disp(12345)');
            testCase.verifySubstring(cmd, 'matlab/r2024b');
        end
    end

    methods (Access = private)
        function [cmd] = submitCommand(~, varargin)
            out = evalc('generic_single_job_submit_wrapper(''disp(12345)'', 1, 7, ''maxTrialNum'', 1, varargin{:})');
            cmd = regexp(out, 'sbatch [^\n]*', 'match', 'once');
        end
    end
end
