classdef FakeSlurmFixture < matlab.unittest.fixtures.Fixture
    % a stand-in Slurm for testing how the job framework talks to Slurm
    % without a cluster. Stubs first on the PATH: sinfo succeeds, sbatch
    % prints a new job id (1000, 1001, ...), squeue lists the jobs whose
    % states a test set with setJobStates, and scancel does nothing. Every
    % call is logged, so tests can count them, and history lists all calls
    % in order. Not available on Windows.
    %
    % A job's states are used one per squeue call, and the last one stays:
    % setJobStates('123_7', {'PD', 'PD', 'R'}) reports the job pending
    % twice and then running. The state 'done' means the job is gone. The
    % jobs run nothing, but a job can write (empty) output files when it
    % is first reported running, and rejectSubmissions makes sbatch fail.

    properties (SetAccess = private)
        Folder
    end

    properties (Access = private)
        OriginalPath
    end

    methods
        function setup(fixture)
            if ispc
                error('pktest:FakeSlurmFixture', 'The fake Slurm needs a POSIX shell.');
            end
            fixture.Folder = tempname;
            mkdir([fixture.Folder, '/state']);
            fixture.writeStub('sinfo', {'echo "PARTITION AVAIL TIMELIMIT NODES STATE NODELIST"'});
            fixture.writeStub('sbatch', { ...
                'r=$(cat "$D/reject" 2>/dev/null || echo 0)', ...
                'if [ "$r" -gt 0 ]; then', ...
                '  echo $((r - 1)) > "$D/reject"', ...
                '  echo "sbatch: error: Batch job submission failed: Socket timed out on send/recv operation" >&2', ...
                '  exit 1', ...
                'fi', ...
                'n=$(cat "$D/next_id" 2>/dev/null || echo 1000)', ...
                'echo $((n + 1)) > "$D/next_id"', ...
                'echo "Submitted batch job $n"'});
            fixture.writeStub('scancel', {});
            fixture.writeStub('squeue', { ...
                'ids=""', ...
                'while [ $# -gt 0 ]; do', ...
                '  if [ "$1" = "-j" ]; then ids="$2"; shift; fi', ...
                '  shift', ...
                'done', ...
                'echo "             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)"', ...
                'for id in $(echo "$ids" | tr "," " "); do', ...
                '  f="$D/state/$id"', ...
                '  [ -f "$f" ] || continue', ...
                '  st=$(head -n 1 "$f")', ...
                '  if [ $(wc -l < "$f") -gt 1 ]; then tail -n +2 "$f" > "$f.next" && mv "$f.next" "$f"; fi', ...
                '  if [ "$st" = "R" ] && [ -f "$f.outputs" ]; then', ...
                '    while IFS= read -r o; do : > "$o"; done < "$f.outputs"', ...
                '    rm "$f.outputs"', ...
                '  fi', ...
                '  [ "$st" = "done" ] && continue', ...
                '  echo "             $id       abc     wrap     user  $st       0:01      1 n0001"', ...
                'done'});
            fixture.OriginalPath = getenv('PATH');
            setenv('PATH', [fixture.Folder, pathsep, fixture.OriginalPath]);
            fixture.SetupDescription = 'Put fake Slurm commands first on the PATH';
        end

        function teardown(fixture)
            setenv('PATH', fixture.OriginalPath);
            if exist(fixture.Folder, 'dir')
                rmdir(fixture.Folder, 's');
            end
        end

        function setJobStates(fixture, jobId, states, outputs)
            % jobId as squeue shows it, e.g., '123_7' for task 7 of job 123.
            % outputs (optional): files the job writes (empty) when squeue
            % first reports it running, e.g., the result of its first task,
            % or all its results while it is still exiting
            fid = fopen(sprintf('%s/state/%s', fixture.Folder, jobId), 'w');
            fprintf(fid, '%s\n', states{:});
            fclose(fid);
            if nargin > 3
                fid = fopen(sprintf('%s/state/%s.outputs', fixture.Folder, jobId), 'w');
                fprintf(fid, '%s\n', outputs{:});
                fclose(fid);
            end
        end

        function rejectSubmissions(fixture, n)
            % the next n sbatch calls fail without submitting a job, as
            % when the Slurm controller is busy
            fid = fopen([fixture.Folder, '/reject'], 'w');
            fprintf(fid, '%d\n', n);
            fclose(fid);
        end

        function [n] = callCount(fixture, cmd)
            fn = sprintf('%s/%s.log', fixture.Folder, cmd);
            n = 0;
            if exist(fn, 'file')
                n = numel(splitlines(strip(fileread(fn))));
            end
        end

        function [args] = calls(fixture, cmd)
            % the arguments of every call, one cell per call
            fn = sprintf('%s/%s.log', fixture.Folder, cmd);
            args = {};
            if exist(fn, 'file')
                args = splitlines(strip(fileread(fn)));
            end
        end

        function [cmds] = history(fixture)
            % every call of all commands in call order, as 'command arguments'
            fn = [fixture.Folder, '/history.log'];
            cmds = {};
            if exist(fn, 'file')
                cmds = splitlines(strip(fileread(fn)));
            end
        end
    end

    methods (Access = protected)
        function [tf] = isCompatible(~, ~)
            % each test gets its own fake Slurm
            tf = false;
        end
    end

    methods (Access = private)
        function writeStub(fixture, cmd, body)
            fn = sprintf('%s/%s', fixture.Folder, cmd);
            fid = fopen(fn, 'w');
            fprintf(fid, '#!/bin/sh\nD="$(dirname "$0")"\necho "$*" >> "$D/%s.log"\necho "%s $*" >> "$D/history.log"\n', cmd, cmd);
            fprintf(fid, '%s\n', body{:});
            fclose(fid);
            fileattrib(fn, '+x');
        end
    end
end
