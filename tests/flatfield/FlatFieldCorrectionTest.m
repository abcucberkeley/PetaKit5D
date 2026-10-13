classdef FlatFieldCorrectionTest < pktest.TestCase
    % XR_LSFlatFieldCorrection and processFFCorrectionFrame with values
    % chosen so the expected result is exact:
    %   ff = (FF - BG) / max(FF - BG), capped below at LowerLimit
    %   out = max(raw - BG, 0) ./ ff + BG   (or + constOffset if given)
    % With BG = 100, FF - BG = 1000/500/200 in three column blocks gives
    % ff = 1/0.5/0.4 (0.2 capped at 0.4), so raw = 300 gives 300/500/600.

    properties (Constant)
        Background = 100
    end

    methods (Test, TestTags = {'Unit'})
        function correctsWithExpectedFormula(testCase)
            [raw, ff, bg] = testCase.exactInputs();
            out = XR_LSFlatFieldCorrection(raw, ff, bg, 'LowerLimit', 0.4);
            testCase.verifyClass(out, 'uint16');
            testCase.verifyEqual(out, testCase.expectedOutput(raw, 100));
        end

        function constOffsetReplacesBackground(testCase)
            [raw, ff, bg] = testCase.exactInputs();
            out = XR_LSFlatFieldCorrection(raw, ff, bg, 'LowerLimit', 0.4, 'constOffset', 5);
            testCase.verifyEqual(out, testCase.expectedOutput(raw, 5));
        end

        function singleInputGivesSingle(testCase)
            [raw, ff, bg] = testCase.exactInputs();
            out = XR_LSFlatFieldCorrection(single(raw), ff, bg, 'LowerLimit', 0.4);
            testCase.verifyClass(out, 'single');
            testCase.verifyEqual(out, single(testCase.expectedOutput(raw, 100)), 'AbsTol', single(1e-3));

            out = XR_LSFlatFieldCorrection(raw, ff, bg, 'LowerLimit', 0.4, 'castDataType', false);
            testCase.verifyClass(out, 'single');
            testCase.verifyEqual(out, single(testCase.expectedOutput(raw, 100)), 'AbsTol', single(1e-3));
        end

        function matlabFallbackMatchesMex(testCase)
            % double input isn't supported by the mex, so this runs the
            % MATLAB implementation
            [raw, ff, bg] = testCase.exactInputs();
            evalc('out = XR_LSFlatFieldCorrection(double(raw), ff, bg, ''LowerLimit'', 0.4);');
            testCase.verifyClass(out, 'double');
            testCase.verifyEqual(out, double(testCase.expectedOutput(raw, 100)), 'AbsTol', 1e-3);
        end

        function flatFieldStackIsAveragedOverZ(testCase)
            [raw, ff, bg] = testCase.exactInputs();
            % mean over z is ff; the deviations cancel
            ff3 = cat(3, ff - 50, ff + 50);
            out = XR_LSFlatFieldCorrection(raw, ff3, bg, 'LowerLimit', 0.4);
            testCase.verifyEqual(out, testCase.expectedOutput(raw, 100));
        end

        function largerImagesAreCenterCropped(testCase)
            [raw, ff, bg] = testCase.exactInputs();
            % 2 extra pixels on each side with values that would change
            % the max normalization if they were used
            ffLarge = padarray(ff, [2, 2], 5000);
            bgLarge = padarray(bg, [2, 2], 0);
            out = XR_LSFlatFieldCorrection(raw, ffLarge, bgLarge, 'LowerLimit', 0.4);
            testCase.verifyEqual(out, testCase.expectedOutput(raw, 100));
        end

        function nonSingleFlatFieldUsesFastPath(testCase)
            % regression for 824f476: uint16 flat field/background images
            % made the mex error and silently fall back to MATLAB
            testCase.assumeTrue(testCase.flatFieldMexLoads(), ...
                'flat_field_correction_3d_mex does not load on this platform (see MexPackagingTest).');
            [raw, ~, bg] = testCase.exactInputs();
            % without rescaling the flat field is used as is: 1 or 2
            ff = ones(size(bg), 'uint16');
            ff(:, 13 : end) = 2;
            log = evalc(['out = XR_LSFlatFieldCorrection(raw, ff, uint16(bg), ''LowerLimit'', 0.4, ', ...
                '''removeFFImBackground'', false, ''LSRescale'', false);']);
            testCase.verifyFalse(contains(log, 'MException'), ...
                sprintf('Fell back to the MATLAB implementation:\n%s', log));
            expected = uint16(max(double(raw) - 100, 0) ./ double(ff) + 100);
            testCase.verifyEqual(out, expected);
        end
    end

    methods (Test, TestTags = {'Integration'})
        function processFrameReadsImagesFromDisk(testCase)
            [raw, ff, bg] = testCase.exactInputs();
            dataPath = testCase.createTempFolder();
            writetiff(single(ff), [dataPath, 'ff.tif']);
            writetiff(single(bg), [dataPath, 'bg.tif']);
            out = processFFCorrectionFrame(raw, [dataPath, 'ff.tif'], [dataPath, 'bg.tif'], [], 0.4);
            testCase.verifyEqual(out, testCase.expectedOutput(raw, 100));
        end
    end

    methods (Access = private)
        function [raw, ff, bg] = exactInputs(testCase)
            bg = zeros(16, 24, 'single') + testCase.Background;
            ff = bg + repelem(single([1000, 500, 200]), 16, 8);
            raw = zeros(16, 24, 4, 'uint16') + 300;
            % below background, clamped to 0 before adding the offset back
            raw(3, :, 2) = 50;
        end

        function [out] = expectedOutput(~, raw, offset)
            ffNorm = repelem([1, 0.5, 0.4], 16, 8);
            out = uint16(max(double(raw) - 100, 0) ./ ffNorm + offset);
        end

        function [tf] = flatFieldMexLoads(~)
            try
                out = flat_field_correction_3d_mex(zeros(2, 2, 2, 'uint16'), ones(2, 2, 'single'), zeros(2, 2, 'single'), true); %#ok<NASGU>
                tf = true;
            catch
                tf = false;
            end
        end
    end
end
