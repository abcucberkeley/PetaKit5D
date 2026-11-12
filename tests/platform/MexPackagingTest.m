classdef MexPackagingTest < pktest.TestCase
    % the ARM macOS binaries must only load libraries that every Mac has or
    % that ship next to them, the Linux cpp-tiff/cpp-zarr binaries must use
    % their own codecs, and the other Linux mex files must not name folders
    % of the machine they were built on. This reads the Mach-O load
    % commands and ELF dynamic sections, so it runs on any OS and doesn't
    % depend on what is installed on the runner.

    properties (Constant)
        AllowedPrefixes = {'@rpath/', '@loader_path/', '@executable_path/', '/usr/lib/', '/System/Library/'}
        % functions of libtiff, zstd, zlib and libdeflate, which MATLAB has too
        CodecPrefixes = {'TIFF', '_TIFF', 'ZSTD_', 'ZDICT_', 'deflate', 'inflate', 'compress', 'uncompress', 'crc32', 'adler32', 'libdeflate_'}
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

        function linuxTiffZarrBinariesUseTheirOwnCodecs(testCase)
            % on Linux, once Tiff or imread had loaded MATLAB's libtiff, the
            % cpp-tiff mex files called it instead of their own (before
            % cpp-tiff/cpp-zarr v1.7.3), so readtiff couldn't read zstd tiffs
            % in that session. Their libraries must bind their own symbols
            % first (-Bsymbolic), and the mex files must not call a codec
            % themselves (that call would go to MATLAB's).
            ioDir = fullfile(pktest.repoRoot(), 'microscopeDataProcessing', 'io');
            libs = {fullfile(ioDir, 'cpp-tiff', 'linux', 'libcppTiff.so'), fullfile(ioDir, 'cpp-zarr', 'linux', 'libcppZarr.so')};
            for i = 1 : numel(libs)
                testCase.verifyTrue(pktest.elfDynamic(libs{i}), sprintf('%s is not linked with -Bsymbolic.', libs{i}));
            end

            dir_info = [dir(fullfile(ioDir, 'cpp-tiff', 'linux', '*.mexa64')); dir(fullfile(ioDir, 'cpp-zarr', 'linux', '*.mexa64'))];
            testCase.assertNotEmpty(dir_info);
            bad = {};
            for i = 1 : numel(dir_info)
                [~, imports] = pktest.elfDynamic(fullfile(dir_info(i).folder, dir_info(i).name));
                codecImports = imports(startsWith(imports, testCase.CodecPrefixes));
                if ~isempty(codecImports)
                    bad{end + 1} = sprintf('%s -> %s', dir_info(i).name, strjoin(codecImports, ', '));
                end
            end
            testCase.verifyEmpty(bad, sprintf('Linux mex files calling codecs directly:\n  %s', strjoin(bad, '\n  ')));
        end

        function linuxMexFilesNameNoFoldersOfTheBuildMachine(testCase)
            % Linux mex files built with a compiler outside the system
            % folders kept its library folder, an absolute path on the
            % build machine, as their RPATH, and patchelf --remove-rpath
            % left that path in their dynamic string table. They need
            % neither: they use the runtime libraries that MATLAB or the
            % system loads. Their RPATH and RUNPATH may only name folders
            % relative to the file ($ORIGIN), and no string of their
            % dynamic string table may be an absolute path. The vendored
            % cpp-tiff and cpp-zarr binaries come from those projects'
            % releases and are left out. (The messages don't repeat the
            % paths.)
            files = testCase.findFiles('*.mexa64');
            files = files(~contains(files, fullfile('io', {'cpp-tiff', 'cpp-zarr'})));
            testCase.assertNotEmpty(files);

            bad = {};
            for i = 1 : numel(files)
                [~, ~, searchPaths, dynStrings] = pktest.elfDynamic(files{i});
                name = strrep(files{i}, [pktest.repoRoot(), filesep], '');
                nOutside = nnz(~startsWith(searchPaths, {'$ORIGIN', '${ORIGIN}'}));
                if nOutside > 0
                    bad{end + 1} = sprintf('%s: %d RPATH/RUNPATH folders not relative to $ORIGIN', name, nOutside);
                end
                nAbsolute = nnz(startsWith(dynStrings, '/'));
                if nAbsolute > 0
                    bad{end + 1} = sprintf('%s: %d absolute paths in its dynamic string table', name, nAbsolute);
                end
            end
            testCase.verifyEmpty(bad, sprintf('Linux mex files naming folders of the build machine:\n  %s', strjoin(bad, '\n  ')));
        end
    end

    methods (Static, Access = private)
        function [files] = findFiles(pattern)
            % (not in folders whose names start with a dot, such as other
            % checkouts of the repository inside it)
            dir_info = dir(fullfile(pktest.repoRoot(), '**', pattern));
            files = fullfile({dir_info.folder}, {dir_info.name});
            rel = strrep(files, [pktest.repoRoot(), filesep], '');
            files = files(~contains(strcat(filesep, rel), [filesep, '.']));
        end
    end
end
