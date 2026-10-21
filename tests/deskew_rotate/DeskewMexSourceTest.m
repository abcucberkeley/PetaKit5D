classdef DeskewMexSourceTest < pktest.TestCase
    % static checks of the deskew/rotation mex sources for bugs that can't
    % be tested by running the code without crashing MATLAB

    properties (TestParameter)
        sourceFile = {'volume_deskew_rotate_warp_mex.cpp', ...
            'skewed_space_interp_volume_deskew_rotate_warp_mex.cpp', ...
            'skewed_space_interp_defined_stepsize_mex.c'}
    end

    methods (Test, TestTags = {'Unit'})
        function simdLoopsStayInsideTheRow(testCase, sourceFile)
            % the AVX float loops used to run up to 7 floats past each row
            % (heap corruption for single-precision input whose first
            % dimension isn't a multiple of 8)
            src = fileread(fullfile(pktest.repoRoot(), 'microscopeDataProcessing', 'deskew_rotate', 'mex', sourceFile));
            lines = splitlines(src);
            unsafe = find(~cellfun(@isempty, regexp(lines, 'for\s*\(\s*size_t\s+i\s*=\s*0\s*;\s*i\s*<\s*len\s*;\s*i\s*\+=\s*(8|16)\s*\)', 'once')));
            testCase.verifyEmpty(unsafe, sprintf('%s: SIMD loop can run past the end of the row at line(s) %s.', ...
                sourceFile, mat2str(unsafe)));
        end
    end
end
