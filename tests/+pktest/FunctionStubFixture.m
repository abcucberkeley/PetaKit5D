classdef FunctionStubFixture < matlab.unittest.fixtures.Fixture
    % replace a function with a stub that returns fixed values, e.g. to
    % make getSystemMemory report little memory so the large-data branch
    % of the code under test runs on small test data. The stub is first on
    % the path while the fixture is set up.

    properties (SetAccess = private)
        FunctionName
        Outputs
        Folder
    end

    properties (Access = private)
        WarningState
    end

    methods
        function fixture = FunctionStubFixture(functionName, varargin)
            % varargin: the stub's outputs, in order
            fixture.FunctionName = functionName;
            fixture.Outputs = varargin;
        end

        function setup(fixture)
            fixture.Folder = strrep(tempname, '\', '/');
            mkdir(fixture.Folder);
            outputs = fixture.Outputs;
            save([fixture.Folder, '/outputs.mat'], 'outputs');
            lines = {
                sprintf('function varargout = %s(varargin)', fixture.FunctionName)
                '% test stub written by pktest.FunctionStubFixture'
                'a = load([fileparts(mfilename(''fullpath'')), ''/outputs.mat'']);'
                'varargout = a.outputs(1 : max(nargout, 1));'
                'end'};
            fid = fopen(sprintf('%s/%s.m', fixture.Folder, fixture.FunctionName), 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
            % a stub of a built-in (e.g., pause) would warn at each use
            fixture.WarningState = warning('off', 'MATLAB:dispatcher:nameConflict');
            addpath(fixture.Folder, '-begin');
            fixture.SetupDescription = sprintf('Replaced %s with a stub', fixture.FunctionName);
        end

        function teardown(fixture)
            rmpath(fixture.Folder);
            clear(fixture.FunctionName);
            warning(fixture.WarningState);
            if exist(fixture.Folder, 'dir')
                rmdir(fixture.Folder, 's');
            end
        end
    end

    methods (Access = protected)
        function [tf] = isCompatible(fixture, other)
            tf = strcmp(fixture.FunctionName, other.FunctionName) && isequal(fixture.Outputs, other.Outputs);
        end
    end
end
