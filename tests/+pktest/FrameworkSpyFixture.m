classdef FrameworkSpyFixture < matlab.unittest.fixtures.Fixture
    % replace generic_computing_frameworks_wrapper with a spy that records
    % its arguments, runs nothing, and reports every task as done. Used to
    % test what the data wrappers pass to the job framework (e.g.,
    % memAllocate). The spy is first on the path while the fixture is set up.
    % Optional constructor arguments: a function that runs inside each call,
    % e.g. to change files as other jobs would while the framework runs, and
    % whether the calls report their tasks as done (default true).

    properties (SetAccess = private)
        Folder
        OnCall
        Done = true
    end

    methods
        function fixture = FrameworkSpyFixture(onCall, done)
            if nargin > 0
                fixture.OnCall = onCall;
            end
            if nargin > 1
                fixture.Done = done;
            end
        end

        function setup(fixture)
            fixture.Folder = strrep(tempname, '\', '/');
            mkdir(fixture.Folder);
            lines = {
                'function [is_done_flag] = generic_computing_frameworks_wrapper(inputFullpaths, outputFullpaths, funcStrs, varargin)'
                '% test spy written by pktest.FrameworkSpyFixture'
                'spyFn = [fileparts(mfilename(''fullpath'')), ''/calls.mat''];'
                'calls = {};'
                'if exist(spyFn, ''file'')'
                '    a = load(spyFn);'
                '    calls = a.calls;'
                'end'
                'calls{end + 1} = struct(''inputFullpaths'', {inputFullpaths}, ''outputFullpaths'', {outputFullpaths}, ...'
                '    ''funcStrs'', {funcStrs}, ''options'', {varargin});'
                'save(spyFn, ''calls'');'
                'onCallFn = [fileparts(mfilename(''fullpath'')), ''/onCall.mat''];'
                'if exist(onCallFn, ''file'')'
                '    a = load(onCallFn);'
                '    a.onCall();'
                'end'
                'a = load([fileparts(mfilename(''fullpath'')), ''/done.mat'']);'
                'is_done_flag = repmat(a.done, numel(inputFullpaths), 1);'
                'end'};
            fid = fopen([fixture.Folder, '/generic_computing_frameworks_wrapper.m'], 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
            if ~isempty(fixture.OnCall)
                onCall = fixture.OnCall;
                save([fixture.Folder, '/onCall.mat'], 'onCall');
            end
            done = fixture.Done;
            save([fixture.Folder, '/done.mat'], 'done');
            addpath(fixture.Folder, '-begin');
            fixture.SetupDescription = 'Replaced generic_computing_frameworks_wrapper with a spy';
        end

        function teardown(fixture)
            rmpath(fixture.Folder);
            clear generic_computing_frameworks_wrapper
            if exist(fixture.Folder, 'dir')
                rmdir(fixture.Folder, 's');
            end
        end

        function [calls] = getCalls(fixture)
            % recorded calls, each with inputFullpaths, outputFullpaths,
            % funcStrs and options (the name-value arguments)
            spyFn = [fixture.Folder, '/calls.mat'];
            calls = {};
            if exist(spyFn, 'file')
                a = load(spyFn);
                calls = a.calls;
            end
        end

        function [value] = getOption(fixture, callInd, name)
            % value of a name-value argument of a recorded call ([] if absent)
            calls = fixture.getCalls();
            options = calls{callInd}.options;
            value = [];
            for i = 1 : 2 : numel(options) - 1
                if strcmpi(char(options{i}), name)
                    value = options{i + 1};
                end
            end
        end
    end

    methods (Access = protected)
        function [tf] = isCompatible(fixture, other)
            tf = isempty(fixture.OnCall) && isempty(other.OnCall) && isequal(fixture.Done, other.Done);
        end
    end
end
