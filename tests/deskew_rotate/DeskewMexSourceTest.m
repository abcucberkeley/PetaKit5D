classdef DeskewMexSourceTest < pktest.TestCase
    % static checks of the deskew/rotation mex sources for bugs that can't
    % be tested by running the code without crashing MATLAB

    properties (TestParameter)
        sourceFile = {'volume_deskew_rotate_warp_mex.cpp', ...
            'skewed_space_interp_volume_deskew_rotate_warp_mex.cpp', ...
            'skewed_space_interp_defined_stepsize_mex.c'}
    end

    methods (Test, TestTags = {'Unit', 'KnownIssue'})
        function simdLoopsStayInsideTheRow(testCase, sourceFile)
            % KNOWN ISSUE: the AVX float loops are written as
            %   for (size_t i = 0; i < len; i += 8) { ... _mm256_storeu_ps(&result[i], ...) }
            % so the last iteration loads and stores up to 7 floats past the
            % end of the row when len is not a multiple of 8. For
            % single-precision input whose first dimension isn't a multiple
            % of 8, this overwrites the next row and finally runs past the
            % output buffer (heap corruption, MATLAB crashes with
            % "malloc(): corrupted top size"). Linux and Windows only (x86
            % SIMD path). Fix: use i + 8 <= len, as the 16-bit loops do.
            src = fileread(fullfile(pktest.repoRoot(), 'microscopeDataProcessing', 'deskew_rotate', 'mex', sourceFile));
            lines = splitlines(src);
            unsafe = find(~cellfun(@isempty, regexp(lines, 'for\s*\(\s*size_t\s+i\s*=\s*0\s*;\s*i\s*<\s*len\s*;\s*i\s*\+=\s*(8|16)\s*\)', 'once')));
            testCase.verifyEmpty(unsafe, sprintf('%s: SIMD loop can run past the end of the row at line(s) %s.', ...
                sourceFile, mat2str(unsafe)));
        end
    end
end
