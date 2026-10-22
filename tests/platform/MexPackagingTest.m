classdef MexPackagingTest < pktest.TestCase
    % the ARM macOS binaries must only load libraries that every Mac has or
    % that ship next to them. This reads the Mach-O load commands, so it
    % runs on any OS and doesn't depend on what is installed on the runner.

    properties (Constant)
        AllowedPrefixes = {'@rpath/', '@loader_path/', '@executable_path/', '/usr/lib/', '/System/Library/'}
    end

    methods (Test, TestTags = {'Unit'})
        function macArmBinariesOnlyLoadBundledOrSystemLibraries(testCase)
            % the PetaKit5D .mexmaca64 files used to load Homebrew's libgomp
            % by absolute path, so on a Mac without Homebrew gcc@13 they only
            % loaded if cpp-tiff/cpp-zarr had loaded the bundled libgomp
            % first. compileAllMexFiles now links the GCC runtime statically;
            % cpp-tiff/cpp-zarr load theirs from next to them (@loader_path).
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

    methods (Static, Access = private)
        function [files] = findFiles(pattern)
            dir_info = dir(fullfile(pktest.repoRoot(), '**', pattern));
            files = fullfile({dir_info.folder}, {dir_info.name});
        end
    end
end
