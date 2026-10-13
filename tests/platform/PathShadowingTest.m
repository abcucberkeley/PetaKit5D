classdef PathShadowingTest < pktest.TestCase
    % setup.m puts every folder of the repo on the path, so a function
    % defined in two places silently resolves to whichever folder comes
    % first, and a repo function with the name of a MATLAB function hides
    % it. Both make fixes "not take effect" and break with path changes or
    % new MATLAB releases.

    properties (Constant)
        % names already defined more than once (mostly vendored point
        % detection code). Don't add to this list; rename the new function.
        KnownDuplicates = {'binterp', 'cmeAnalysis3D', 'conv3fast', 'createDistanceMatrix', ...
            'crop_deskew_rotate_demo_function_parser', 'cropMovie3D', 'deskewData', 'deskewFrame3D', ...
            'estGaussianAmplitude3D', 'filterGauss1D', 'filterGauss2D', 'filterGauss3D', ...
            'fitGaussianModeToPDF', 'getCellDir', 'getCellVolume', 'getCropRegions3D', 'getDirFromPath', ...
            'getFluorPropStruct', 'getLifetimeData', 'getMovieName', 'getRGBindex', 'getShortPath', ...
            'getVisitorIndex', 'interpln', 'loadConditionData', 'loadConditionData3D', 'loadFigureSettings', ...
            'loadTracks', 'loadTrackSettings', 'locmax3d', 'plotIntensityCohorts', 'plotLifetimes', ...
            'pointSourceDetection3D', 'recursiveDir', 'rgbOverlay', 'robustMean', 'rotateFrame3D', ...
            'runDetection3D', 'runLifetimeAnalysis', 'runTracking3D', 'runTrackProcessing3D', 'scaleContrast', ...
            'scaleEDFs', 'setupFigure', 'sortStringsByToken', 'stackviewer', 'tests', 'writetiff_obsolete'}

        % duplicated functions that the main processing code calls, and the
        % copy it currently uses
        PipelineCopies = {'deskewFrame3D', 'pointDetection/llsmtools/transforms/deskewFrame3D.m'; ...
            'rotateFrame3D', 'pointDetection/llsmtools/transforms/rotateFrame3D.m'; ...
            'scaleContrast', 'imageProcessing/scaleContrast.m'}
    end

    methods (Test, TestTags = {'Unit'})
        function noNewDuplicateFunctionNames(testCase)
            names = testCase.repoFunctionNames();
            [uniqNames, ~, ic] = unique(names);
            duplicates = uniqNames(accumarray(ic, 1) > 1);

            newDuplicates = setdiff(duplicates, testCase.KnownDuplicates);
            testCase.verifyEmpty(newDuplicates, sprintf(['These function names are now defined in more ', ...
                'than one place: %s. Only one of them runs, depending on the path order.'], strjoin(newDuplicates, ', ')));

            resolved = setdiff(testCase.KnownDuplicates, duplicates);
            if ~isempty(resolved)
                testCase.log(sprintf('No longer duplicated (remove from KnownDuplicates): %s', strjoin(resolved, ', ')));
            end
        end

        function pipelineUsesExpectedCopies(testCase)
            root = [strrep(pktest.repoRoot(), '\', '/'), '/'];
            for i = 1 : size(testCase.PipelineCopies, 1)
                fn = strrep(which(testCase.PipelineCopies{i, 1}), '\', '/');
                testCase.verifyEqual(strrep(fn, root, ''), testCase.PipelineCopies{i, 2}, ...
                    sprintf('%s now resolves to a different copy, which changes the processing results.', ...
                    testCase.PipelineCopies{i, 1}));
            end
        end

        function noFunctionShadowsMatlab(testCase)
            % e.g., a new MATLAB release adding a function with the same
            % name as a repo function
            names = unique(testCase.repoFunctionNames());
            mr = [matlabroot, filesep];
            shadowed = {};
            for i = 1 : numel(names)
                w = which('-all', names{i});
                % only functions: class methods (in @class folders or
                % classdef files) and private functions can't be shadowed
                [~, fsns] = fileparts(w);
                inMatlab = w((startsWith(w, mr) & strcmp(fsns, names{i}) | contains(w, 'built-in')) & ...
                    ~contains(w, [filesep, '@']) & ~contains(w, [filesep, 'private', filesep]));
                if ~isempty(inMatlab)
                    shadowed{end + 1} = sprintf('%s (MATLAB: %s)', names{i}, inMatlab{1});
                end
            end
            testCase.verifyEmpty(shadowed, sprintf('Repo functions with the same name as MATLAB functions:\n  %s', ...
                strjoin(shadowed, '\n  ')));
        end
    end

    methods (Static, Access = private)
        function [names] = repoFunctionNames()
            root = pktest.repoRoot();
            dir_info = dir(fullfile(root, '**', '*.m'));
            dir_info = dir_info(~startsWith({dir_info.folder}, fullfile(root, 'tests')));
            [~, names] = cellfun(@fileparts, {dir_info.name}, 'unif', 0);
        end
    end
end
