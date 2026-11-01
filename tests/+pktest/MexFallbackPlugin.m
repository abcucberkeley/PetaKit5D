classdef MexFallbackPlugin < matlab.unittest.plugins.TestRunnerPlugin
    % collects the mex fallbacks each test recorded (mex_fallback: the MATLAB
    % code ran instead of a mex) for the summary of run_tests

    properties (SetAccess = private)
        % struct array: test, mex, kind, identifier, message, count
        Records = struct('test', {}, 'mex', {}, 'kind', {}, 'identifier', {}, 'message', {}, 'count', {})
    end

    methods (Access = protected)
        function runTest(plugin, pluginData)
            % fallbacks since the last test (in a class setup or teardown),
            % which the setup of the test (pktest.TestCase) checks and clears
            plugin.collect(sprintf('(class setup or teardown before %s)', pluginData.Name));
            runTest@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
        end

        function teardownTestMethod(plugin, pluginData)
            % the fallbacks of the test, before its teardown checks and clears them
            plugin.collect(pluginData.Name);
            teardownTestMethod@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
        end
    end

    methods
        function [lines] = summary(plugin)
            % one line per mex and kind: the number of fallbacks and the tests
            lines = {};
            r = plugin.Records;
            if isempty(r)
                return;
            end
            keys = strcat({r.mex}, {' '}, {r.kind});
            [ukeys, inds] = unique(keys, 'stable');
            for i = 1 : numel(ukeys)
                sel = strcmp(keys, ukeys{i});
                tests = unique({r(sel).test}, 'stable');
                ids = unique({r(sel).identifier}, 'stable');
                ids = ids(~cellfun(@isempty, ids));
                idStr = '';
                if ~isempty(ids)
                    idStr = sprintf(' (%s)', strjoin(ids, ', '));
                end
                testStr = strjoin(tests(1 : min(5, end)), ', ');
                if numel(tests) > 5
                    testStr = sprintf('%s and %d more', testStr, numel(tests) - 5);
                end
                lines{end + 1} = sprintf('%s %s%s: %d times in %d tests: %s', r(inds(i)).mex, r(inds(i)).kind, ...
                    idStr, sum([r(sel).count]), numel(tests), testStr); %#ok<AGROW>
            end
        end
    end

    methods (Access = private)
        function [] = collect(plugin, testName)
            % (PetaKit5D is on the path only while the tests' shared fixture is set up)
            try
                records = mex_fallback('-log');
            catch
                return;
            end
            for i = 1 : numel(records)
                rec = records(i);
                rec.test = testName;
                plugin.Records(end + 1) = orderfields(rec, plugin.Records);
            end
        end
    end
end
