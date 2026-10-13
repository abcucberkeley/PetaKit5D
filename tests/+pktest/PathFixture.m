classdef PathFixture < matlab.unittest.fixtures.Fixture
    % add PetaKit5D to the path with setup.m, as users do, and restore the
    % original path afterwards. Shared by all tests through pktest.TestCase.

    properties (Access = private)
        OriginalPath
    end

    methods
        function setup(fixture)
            fixture.OriginalPath = path;
            root = pktest.repoRoot();
            addpath(root);
            % setup.m prints the hostname and path messages, keep logs clean
            evalc('feval(''setup'', root)');
            fixture.SetupDescription = sprintf('Added PetaKit5D (%s) to the path with setup.m', root);
        end

        function teardown(fixture)
            path(fixture.OriginalPath);
            fixture.TeardownDescription = 'Restored the original MATLAB path';
        end
    end
end
