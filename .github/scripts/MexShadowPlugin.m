classdef MexShadowPlugin < matlab.unittest.plugins.TestRunnerPlugin
    % for the CI job exit-probe (suite_exit_probe with a mex variant): put a mex
    % folder at the front of the path after the shared fixtures run (the path
    % fixture's setup.m does addpath(genpath(repo)), which would otherwise be
    % ahead), and again before each test class, so the tests load the variant's
    % mex files in place of the committed ones. Not a test helper.

    properties
        MexDir char
    end

    methods
        function plugin = MexShadowPlugin(mexDir)
            plugin.MexDir = mexDir;
        end
    end

    methods (Access = protected)
        function setupSharedTestFixture(plugin, pluginData)
            setupSharedTestFixture@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
            addpath(plugin.MexDir, '-begin');
        end

        function setupTestClass(plugin, pluginData)
            setupTestClass@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
            addpath(plugin.MexDir, '-begin');
        end
    end
end
