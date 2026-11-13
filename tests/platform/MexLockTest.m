classdef MexLockTest < pktest.TestCase
    % on macOS, MATLAB crashed on exit (segmentation fault, exit status
    % 137) after a session had called a mex file built with GCC's OpenMP:
    % its runtime leaves thread exit handlers in the mex file, which MATLAB
    % unloads at exit before the threads end. A mex file that uses OpenMP
    % keeps itself loaded (mexLock, on macOS) from its first call;
    % omw_backprojector_mex, min_bbox_3d_mex, replace_nan_with_value and
    % replace_nan_inf_with_value did not. The vendored cpp-tiff and
    % cpp-zarr mex files lock themselves in their releases and are left out.

    methods (Test, TestTags = {'Unit'})
        function openMPMexSourcesLockThemselves(testCase)
            % (reads the sources, so it runs on any OS)
            srcs = testCase.openMPMexSources();
            testCase.assertNotEmpty(srcs);
            bad = {};
            for i = 1 : numel(srcs)
                if ~contains(fileread(srcs{i}), 'mexLock()')
                    bad{end + 1} = strrep(srcs{i}, [pktest.repoRoot(), filesep], '');
                end
            end
            testCase.verifyEmpty(bad, sprintf('OpenMP mex sources that don''t call mexLock:\n  %s', strjoin(bad, '\n  ')));
        end

        function macOSMexFilesStayLoadedAfterTheirFirstCall(testCase)
            % the binaries of those sources are locked after a call. A call
            % without inputs is enough: the lock comes before the input
            % checks, whose error the test ignores.
            testCase.assumeTrue(ismac, 'The mex files lock themselves on macOS only.');
            srcs = testCase.openMPMexSources();
            names = {};
            notLocked = {};
            for i = 1 : numel(srcs)
                [folder, name] = fileparts(srcs{i});
                mexFn = fullfile(folder, [name, '.', mexext]);
                if ~exist(mexFn, 'file') || ~strcmp(which(name), mexFn)
                    continue;
                end
                names{end + 1} = name;
                try
                    feval(name);
                catch
                end
                if ~mislocked(name)
                    notLocked{end + 1} = name;
                end
            end
            testCase.assertNotEmpty(names);
            testCase.verifyEmpty(notLocked, sprintf('%s mex files not locked after a call: %s', mexext, ...
                strjoin(notLocked, ', ')));
        end
    end

    methods (Static, Access = private)
        function [srcs] = openMPMexSources()
            % the C and C++ mex sources that include omp.h (not in folders
            % whose names start with a dot, such as other checkouts of the
            % repository inside it, the Python port or the vendored projects)
            root = pktest.repoRoot();
            dir_info = [dir(fullfile(root, '**', '*.c')); dir(fullfile(root, '**', '*.cpp'))];
            srcs = fullfile({dir_info.folder}, {dir_info.name});
            rel = strcat(filesep, strrep(srcs, [root, filesep], ''));
            skip = contains(rel, [filesep, '.']) | startsWith(rel, [filesep, 'python', filesep]) | ...
                contains(rel, fullfile('io', {'cpp-tiff', 'cpp-zarr'})) | contains(rel, 'Parallel_Imaris_Writer');
            srcs = srcs(~skip);
            keep = false(size(srcs));
            for i = 1 : numel(srcs)
                txt = fileread(srcs{i});
                keep(i) = contains(txt, 'mexFunction') && ~isempty(regexp(txt, '#\s*include\s*[<"]omp\.h[>"]', 'once'));
            end
            srcs = srcs(keep);
        end
    end
end
