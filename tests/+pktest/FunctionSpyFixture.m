classdef FunctionSpyFixture < matlab.unittest.fixtures.Fixture
    % replace a function that returns nothing with a spy that records its
    % arguments and does nothing else, to test what callers pass to it
    % (e.g., the job settings a wrapper forwards). The spy is first on the
    % path while the fixture is set up. An optional function given to the
    % constructor runs at each call (e.g., to make files another job would).
    % An optional maximum number of calls makes the spy error once it is
    % exceeded, so that a caller stuck in a loop (e.g., of pause) fails the
    % test instead of hanging it.

    properties (SetAccess = private)
        FunctionName
        Folder
        OnCall
        MaxCalls = Inf
    end

    methods
        function fixture = FunctionSpyFixture(functionName, onCall, maxCalls)
            fixture.FunctionName = functionName;
            if nargin > 1
                fixture.OnCall = onCall;
            end
            if nargin > 2
                fixture.MaxCalls = maxCalls;
            end
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
                'save(spyFn, ''calls'');'};
            if isfinite(fixture.MaxCalls)
                lines = [lines; {
                    sprintf('if numel(calls) > %d', fixture.MaxCalls)
                    sprintf('    error(''pktest:FunctionSpyFixture:TooManyCalls'', ''%s was called more than %d times.'');', ...
                        fixture.FunctionName, fixture.MaxCalls)
                    'end'}];
            end
            lines = [lines; {
                'onCallFn = [fileparts(mfilename(''fullpath'')), ''/onCall.mat''];'
                'if exist(onCallFn, ''file'')'
                '    a = load(onCallFn);'
                '    a.onCall();'
                'end'
                'end'}];
            fid = fopen(sprintf('%s/%s.m', fixture.Folder, fixture.FunctionName), 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
            if ~isempty(fixture.OnCall)
                onCall = fixture.OnCall;
                save([fixture.Folder, '/onCall.mat'], 'onCall');
            end
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
            tf = strcmp(fixture.FunctionName, other.FunctionName) && isempty(fixture.OnCall) && isempty(other.OnCall) ...
                && isequal(fixture.MaxCalls, other.MaxCalls);
        end
    end
end
