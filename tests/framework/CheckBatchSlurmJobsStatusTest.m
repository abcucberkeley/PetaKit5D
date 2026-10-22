classdef CheckBatchSlurmJobsStatusTest < pktest.TestCase
    % check_batch_slurm_jobs_status, which the job framework uses to poll
    % its jobs, with pktest.FakeSlurmFixture: 1 running, 0 pending, -1 any
    % other state or not in the queue

    methods (Test, TestTags = {'Unit'})
        function mapsQueueStatesToJobs(testCase)
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('100_1', {'R'});
            slurm.setJobStates('100_2', {'PD'});
            slurm.setJobStates('100_3', {'CG'});
            slurm.setJobStates('101_1', {'R'});
            % 100_4 isn't in the queue, and job id -1 was never submitted
            status = check_batch_slurm_jobs_status([100; 100; 100; 100; 101; -1], [1; 2; 3; 4; 1; 1]);
            testCase.verifyEqual(status, [1; 0; -1; -1; 1; -1]);
        end

        function tasksSharingAJobGetItsStatus(testCase)
            % tasks batched into one job have the same job and array id
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('200_1', {'R'});
            slurm.setJobStates('200_2', {'PD'});
            status = check_batch_slurm_jobs_status([200; 200; 200; 200], [1; 1; 2; 2]);
            testCase.verifyEqual(status, [1; 1; 0; 0]);
            testCase.verifyEqual(slurm.callCount('squeue'), 1, 'All jobs should be checked with one squeue call.');
        end
    end
end
