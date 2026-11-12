classdef CheckBatchSlurmJobsStatusTest < pktest.TestCase
    % check_batch_slurm_jobs_status, which the job framework uses to poll
    % its jobs, with pktest.FakeSlurmFixture: 1 running, 0 pending, -1 any
    % other state or not in the queue (or unknown, the second output, if
    % squeue failed)

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

        function jobsInAnyOrderGetTheirOwnStatus(testCase)
            % distinct jobs whose ids aren't in order (e.g., the job of a
            % retry has a newer id than the jobs after it) got the states in
            % the order of their sorted ids, so a pending job could get the
            % state of a running one and a running job that of a gone one
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('101_1', {'PD'});
            slurm.setJobStates('100_1', {'R'});
            % 102_1 isn't in the queue
            status = check_batch_slurm_jobs_status([101; 100; 102], [1; 1; 1]);
            testCase.verifyEqual(status, [0; 1; -1]);
            % with tasks sharing a job, as before
            status = check_batch_slurm_jobs_status([101; 100; 101], [1; 1; 1]);
            testCase.verifyEqual(status, [0; 1; 0]);
        end

        function failedQueryLeavesStatesUnknown(testCase)
            % a failed squeue (a timeout of a busy controller) read as all
            % jobs gone, which the framework then submitted again; their
            % states are unknown now. "Invalid job id" means gone, as before
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('100_1', {'R'});
            slurm.failQueries(1);
            [status, unknown] = check_batch_slurm_jobs_status([100; -1], [1; 1]);
            testCase.verifyEqual([status, unknown], [-1, 1; -1, 0]);
            [status, unknown] = check_batch_slurm_jobs_status([100; -1], [1; 1]);
            testCase.verifyEqual([status, unknown], [1, 0; -1, 0]);
            slurm.failQueries(1, 'slurm_load_jobs error: Invalid job id specified');
            [status, unknown] = check_batch_slurm_jobs_status(100, 1);
            testCase.verifyEqual([status, unknown], [-1, 0]);
        end

        function singleJobOfAFailedQueryMayStillRun(testCase)
            % check_slurm_job_status (single jobs) read a failed squeue as
            % the job gone, so its callers submitted or ran the task again;
            % the job may still run now (1, and unknown). "Invalid job id"
            % means gone, as before
            testCase.assumeFalse(ispc, 'The fake Slurm needs a POSIX shell.');
            slurm = testCase.applyFixture(pktest.FakeSlurmFixture);
            slurm.setJobStates('100_1', {'PD'});
            slurm.failQueries(1);
            testCase.verifyEqual(check_slurm_job_status(101, 1), 1);
            slurm.failQueries(1);
            [status, unknown] = check_slurm_job_status(100, 1);
            testCase.verifyEqual([status, unknown], [1, 1]);
            [status, unknown] = check_slurm_job_status(100, 1);
            testCase.verifyEqual([status, unknown], [0, 0]);
            testCase.verifyEqual(check_slurm_job_status(101, 1), -1);
            slurm.failQueries(1, 'slurm_load_jobs error: Invalid job id specified');
            [status, unknown] = check_slurm_job_status(100, 1);
            testCase.verifyEqual([status, unknown], [-1, 0]);
        end
    end
end
