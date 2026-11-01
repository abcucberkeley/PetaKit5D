classdef MexInventoryTest < pktest.TestCase
    % every mex function built for Linux must also be built for Windows and
    % ARM macOS, except the known gaps (pktest.mexPlatformGaps, which also
    % tells pktest.TestCase which missing mex functions are expected). This
    % catches a mex function that is added or rebuilt for one platform but
    % not committed for the others. Linux binaries are the reference set.
    %
    % When you build a missing binary, remove it from the gap list.

    methods (Test, TestTags = {'Unit'})
        function everyLinuxMexIsBuiltForThisPlatform(testCase)
            ext = mexext;
            if ~any(strcmp(ext, {'mexw64', 'mexmaca64'}))
                testCase.assumeFail(sprintf('Only checked on Windows and ARM macOS (this is %s; Linux is the reference).', ext));
            end
            knownGaps = pktest.mexPlatformGaps(ext);

            linuxNames = testCase.mexNames('mexa64');
            platformNames = testCase.mexNames(ext);

            missing = setdiff(setdiff(linuxNames, platformNames), knownGaps);
            testCase.verifyEmpty(missing, sprintf(['These mex functions have a Linux build but no .%s build: %s. ', ...
                'Build and commit them, or add them to the known gaps with a reason.'], ext, strjoin(missing, ', ')));

            nowBuilt = intersect(knownGaps, platformNames);
            if ~isempty(nowBuilt)
                testCase.log(sprintf('Listed as known gaps but now built for %s (remove from the gap list): %s', ...
                    ext, strjoin(nowBuilt, ', ')));
            end
        end
    end

    methods (Static, Access = private)
        function [names] = mexNames(ext)
            dir_info = dir(fullfile(pktest.repoRoot(), '**', ['*.', ext]));
            [~, names] = cellfun(@fileparts, {dir_info.name}, 'unif', 0);
            names = unique(names);
        end
    end
end
