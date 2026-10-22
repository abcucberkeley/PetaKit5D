classdef FunctionSpyFixture < matlab.unittest.fixtures.Fixture
    % replace a function that returns nothing with a spy that records its
    % arguments and does nothing else, to test what callers pass to it
    % (e.g., the job settings a wrapper forwards). The spy is first on the
    % path while the fixture is set up.

    properties (SetAccess = private)
        FunctionName
        Folder
    end

    methods
        function fixture = FunctionSpyFixture(functionName)
            fixture.FunctionName = functionName;
        end

        function setup(fixture)
            fixture.Folder = strrep(tempname, '\', '/');
            mkdir(fixture.Folder);
            lines = {
                sprintf('function %s(varargin)', fixture.FunctionName)
                '% test spy written by pktest.FunctionSpyFixture'
                'spyFn = [fileparts(mfilename(''fullpath'')), ''/calls.mat''];'
                'calls = {};'
                'if exist(spyFn, ''file'')'
                '    a = load(spyFn);'
                '    calls = a.calls;'
                'end'
                'calls{end + 1} = varargin;'
                'save(spyFn, ''calls'');'
                'end'};
            fid = fopen(sprintf('%s/%s.m', fixture.Folder, fixture.FunctionName), 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
            addpath(fixture.Folder, '-begin');
            fixture.SetupDescription = sprintf('Replaced %s with a spy', fixture.FunctionName);
        end

        function teardown(fixture)
            rmpath(fixture.Folder);
            clear(fixture.FunctionName);
            if exist(fixture.Folder, 'dir')
                rmdir(fixture.Folder, 's');
            end
        end

        function [calls] = getCalls(fixture)
            % the arguments of each recorded call
            spyFn = [fixture.Folder, '/calls.mat'];
            calls = {};
            if exist(spyFn, 'file')
                a = load(spyFn);
                calls = a.calls;
            end
        end

        function [value] = getOption(fixture, callInd, name)
            % value of a name-value argument of a recorded call ([] if absent);
            % the arguments before the first name are skipped
            calls = fixture.getCalls();
            args = calls{callInd};
            value = [];
            for i = 1 : numel(args) - 1
                if (ischar(args{i}) || isstring(args{i})) && strcmpi(char(args{i}), name)
                    value = args{i + 1};
                end
            end
        end
    end

    methods (Access = protected)
        function [tf] = isCompatible(fixture, other)
            tf = strcmp(fixture.FunctionName, other.FunctionName);
        end
    end
end
