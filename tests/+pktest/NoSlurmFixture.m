classdef NoSlurmFixture < matlab.unittest.fixtures.Fixture
    % hide Slurm from the code under test by putting failing sinfo, sbatch,
    % squeue and scancel stubs first on the PATH. Some wrappers don't
    % forward parseCluster to the job framework (e.g., XR_zarrToTiff_wrapper),
    % so on a machine with Slurm they would submit real jobs. With the
    % stubs, checkSlurmCluster sees no cluster and everything runs locally.
    % Windows has no Slurm, so nothing is changed there.

    properties (Access = private)
        OriginalPath
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
            setenv('PATH', [fixture.StubFolder, pathsep, fixture.OriginalPath]);
            fixture.SetupDescription = 'Hid Slurm commands behind failing stubs';
        end

        function teardown(fixture)
            if ispc
                return;
            end
            setenv('PATH', fixture.OriginalPath);
            if exist(fixture.StubFolder, 'dir')
                rmdir(fixture.StubFolder, 's');
            end
        end
    end
end
