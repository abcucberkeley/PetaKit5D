classdef NoSlurmFixtureTest < pktest.TestCase
    % pktest.NoSlurmFixture, which all tests share: besides hiding the Slurm
    % commands, it hides the Slurm job the tests run in, so the code under
    % test acts as outside a job wherever the tests run. A master in a Slurm
    % job computes only the tasks that fit the job's memory, so in a job of
    % 16 GB TaskMetricsTest/masterComputedTasksAreRecordedAsMaster failed,
    % while it passed in larger jobs and outside Slurm.

    methods (Test, TestTags = {'Unit'})
        function slurmJobOfTheTestsIsHiddenAndRestored(testCase)
            testCase.assumeFalse(ispc, 'There is no Slurm on Windows.');
            % the job the tests run in
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_JOB_ID', '123'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURM_MEM_PER_NODE', '16384'));
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture('SLURMD_NODENAME', 'node12'));
            fixture = pktest.NoSlurmFixture;
            fixture.setup();
            hidden = {getenv('SLURM_JOB_ID'), getenv('SLURM_MEM_PER_NODE'), getenv('SLURMD_NODENAME')};
            % (the master of a job framework call acts as outside Slurm)
            [~, ~, masterMem] = master_task_fits(100, 1, 20.9);
            [~, child] = system('printf "%s" "${SLURM_JOB_ID-unset}"');
            fixture.teardown();

            testCase.verifyEqual(hidden, {'', '', ''}, 'The Slurm job is visible to the code under test.');
            testCase.verifyEqual(masterMem, inf, 'The master acts as in the Slurm job.');
            testCase.verifyEqual(child, 'unset', 'A child process sees the Slurm job.');
            testCase.verifyEqual({getenv('SLURM_JOB_ID'), getenv('SLURM_MEM_PER_NODE'), getenv('SLURMD_NODENAME')}, ...
                {'123', '16384', 'node12'}, 'The Slurm job was not restored.');
        end
    end
end
