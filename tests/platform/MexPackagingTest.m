classdef MexPackagingTest < pktest.TestCase
    % the ARM macOS binaries must only load libraries that every Mac has or
    % that ship next to them. This reads the Mach-O load commands, so it
    % runs on any OS and doesn't depend on what is installed on the runner.

    properties (Constant)
        AllowedPrefixes = {'@rpath/', '@loader_path/', '@executable_path/', '/usr/lib/', '/System/Library/'}
    end

    methods (Test, TestTags = {'Unit', 'KnownIssue'})
        function macArmBinariesOnlyLoadBundledOrSystemLibraries(testCase)
            % KNOWN ISSUE: the 15 crop, deskew/rotate, flat-field, stitch
            % and MIP .mexmaca64 files and the bundled cpp-zarr
            % libcppZarr.dylib load /opt/homebrew/opt/gcc@13/lib/gcc/13/libgomp.1.dylib
            % by absolute path. On a Mac without Homebrew gcc@13 they only
            % load if a cpp-tiff/cpp-zarr call happened to load the bundled
            % libgomp (which has that install name) first; otherwise the
            % code silently falls back to slower MATLAB code, or errors
            % where there is no fallback. Fix: point the dependency at the
            % bundled copy, e.g. install_name_tool -change <path>
            % @loader_path/libgomp.1.dylib, and ship libgomp next to them.
            files = [testCase.findFiles('*.mexmaca64'), testCase.findFiles('*.dylib')];
            files = files(contains(files, {'.mexmaca64', [filesep, 'macArm', filesep]}));
            testCase.assertNotEmpty(files);

            bad = {};
            for i = 1 : numel(files)
                deps = pktest.machoDylibs(files{i});
                for j = 1 : numel(deps)
                    if ~any(startsWith(deps{j}, testCase.AllowedPrefixes))
                        bad{end + 1} = sprintf('%s -> %s', strrep(files{i}, [pktest.repoRoot(), filesep], ''), deps{j});
                    end
                end
            end
            testCase.verifyEmpty(bad, sprintf('ARM macOS binaries with absolute library paths:\n  %s', strjoin(bad, '\n  ')));
        end
    end

    methods (Test, TestTags = {'Unit'})
        function cppTiffAndZarrLoadBundledLibraries(testCase)
            % the cpp-tiff/cpp-zarr mex files themselves are packaged
            % correctly; keep them that way
            files = [testCase.findFiles('getImageSizeMex.mexmaca64'), testCase.findFiles('parallel*Tiff.mexmaca64'), ...
                testCase.findFiles('*Zarr*.mexmaca64')];
            testCase.assertNotEmpty(files);
            for i = 1 : numel(files)
                deps = pktest.machoDylibs(files{i});
                testCase.verifyTrue(all(startsWith(deps, testCase.AllowedPrefixes)), ...
                    sprintf('%s loads %s', files{i}, strjoin(deps, ', ')));
            end
        end
    end

    methods (Static, Access = private)
        function [files] = findFiles(pattern)
            dir_info = dir(fullfile(pktest.repoRoot(), '**', pattern));
            files = fullfile({dir_info.folder}, {dir_info.name});
        end
    end
end
