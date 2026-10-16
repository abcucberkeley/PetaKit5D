classdef RealSlurmFixture < matlab.unittest.fixtures.Fixture
    % make the real Slurm commands visible again by taking the stub folder
    % of pktest.NoSlurmFixture off the PATH, for the opt-in cluster tests.
    % The original PATH is restored on teardown.

    properties (Access = private)
        OriginalPath
    end

    methods
        function setup(fixture)
            fixture.OriginalPath = getenv('PATH');
            stubFolder = getenv('PETAKIT5D_TEST_SLURM_STUBS');
            if ~isempty(stubFolder)
                entries = strsplit(fixture.OriginalPath, pathsep);
                setenv('PATH', strjoin(entries(~strcmp(entries, stubFolder)), pathsep));
            end
            fixture.SetupDescription = 'Made the real Slurm commands visible';
        end

        function teardown(fixture)
            setenv('PATH', fixture.OriginalPath);
        end
    end
end
