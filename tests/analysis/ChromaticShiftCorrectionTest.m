classdef ChromaticShiftCorrectionTest < pktest.TestCase
    % chromatic shift correction (demo_useful_tools.m) with integer voxel
    % offsets. Both channels image the same synthetic scene; channel c voxel
    % i shows scene voxel i + offset(c, :), so after correction the channels
    % must be identical where they overlap.

    properties (Constant)
        ImSize = [48, 40, 24]
        % offsets from the demo: no shift for CamA; [1, -5, -2] for CamB
        Offset = [0, 0, 0; 1, -5, -2]
        Fsn = 'Scan_Iter_0000_0000_%s_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t'
    end

    methods (Test, TestTags = {'Integration'})
        function validModeAlignsChannels(testCase)
            [dataPath, scene, base] = testCase.shiftedChannels(testCase.Offset);

            testCase.runCorrection(dataPath, 'valid', 'chromaticOffset', testCase.Offset);

            [outA, outB] = testCase.readOutputs(dataPath);
            % overlap of the channels in scene coordinates
            s = max(1 + testCase.Offset, [], 1);
            t = min(testCase.ImSize + testCase.Offset, [], 1);
            testCase.verifySize(outA, t - s + 1);
            testCase.verifyEqual(outB, outA, 'Corrected channels are not aligned.');
            expected = scene(s(1) + base(1) : t(1) + base(1), s(2) + base(2) : t(2) + base(2), s(3) + base(3) : t(3) + base(3));
            testCase.verifyEqual(outA, expected);
        end

        function offsetsEstimatedFromPSFs(testCase)
            % without chromaticOffset, offsets are measured from bead images
            % (one per channel) of the same beads
            [dataPath, ~, ~] = testCase.shiftedChannels(testCase.Offset);
            psfPath = [dataPath, 'PSF/'];
            mkdir(psfPath);
            psfFns = {[psfPath, 'psf_CamA.tif'], [psfPath, 'psf_CamB.tif']};
            bead = [20, 20, 12];
            for c = 1 : 2
                writetiff(uint16(pktest.beads([40, 40, 24], bead - testCase.Offset(c, :), [1.5, 2], 2000, 100)), psfFns{c});
            end

            testCase.runCorrection(dataPath, 'valid', 'chromaticOffset', [], 'psfFullpaths', psfFns, ...
                'maxOffset', [8, 8, 6]);

            [outA, outB] = testCase.readOutputs(dataPath);
            testCase.verifyEqual(outB, outA, 'Channels corrected with the measured offsets are not aligned.');
            testCase.verifyTrue(exist([dataPath, 'Chromatic_Shift_Corrected/measurement/chromatic_shift_measurement_result.mat'], 'file') == 2);
        end

        function memAllocateIsInGigabytes(testCase)
            % 979fa55 regression: memAllocate (GB) was computed in bytes, so
            % every Slurm job requested the maximum number of CPUs
            dataPath = testCase.shiftedChannels(testCase.Offset);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);

            testCase.runCorrection(dataPath, 'valid', 'chromaticOffset', testCase.Offset);

            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), prod(testCase.ImSize) * 4 * 2.5 / 1024 ^ 3, 'RelTol', 1e-12);
        end

        function zarrOutput(testCase)
            [dataPath, ~, ~] = testCase.shiftedChannels(testCase.Offset);
            testCase.runCorrection(dataPath, 'valid', 'chromaticOffset', testCase.Offset, 'saveZarr', true);
            outA = readzarr([dataPath, 'Chromatic_Shift_Corrected/', sprintf(testCase.Fsn, 'CamA'), '.zarr']);
            outB = readzarr([dataPath, 'Chromatic_Shift_Corrected/', sprintf(testCase.Fsn, 'CamB'), '.zarr']);
            testCase.verifyEqual(outB, outA);
        end
    end

    methods (Test, TestTags = {'Integration', 'KnownIssue'})
        function sameModeKeepsImageSize(testCase)
            % KNOWN ISSUE: XR_chromatic_shift_correction_data_wrapper.m:160
            % tests isempty(neworigin) (lowercase), an undefined variable,
            % so 'same' mode always errors. Fix: use newOrigin.
            [dataPath, ~, ~] = testCase.shiftedChannels(testCase.Offset);

            testCase.runCorrection(dataPath, 'same', 'chromaticOffset', testCase.Offset);

            [outA, outB] = testCase.readOutputs(dataPath);
            testCase.verifySize(outA, testCase.ImSize);
            testCase.verifySize(outB, testCase.ImSize);
            both = outA > 0 & outB > 0;
            testCase.verifyEqual(outB(both), outA(both), 'Corrected channels are not aligned.');
        end

        function fullModeContainsBothChannels(testCase)
            % KNOWN ISSUE: in 'full' mode (wrapper lines 166-168), fullSize
            % is computed from the bbox starts only (columns 1:3 instead of
            % 4:6), so the output is offset-sized instead of the union of
            % the channels, and the bbox is relative to the canvas instead
            % of each channel (sign flipped).
            [dataPath, ~, ~] = testCase.shiftedChannels(testCase.Offset);

            testCase.runCorrection(dataPath, 'full', 'chromaticOffset', testCase.Offset);

            [outA, outB] = testCase.readOutputs(dataPath);
            unionSize = max(testCase.ImSize + testCase.Offset, [], 1) - min(1 + testCase.Offset, [], 1) + 1;
            testCase.verifySize(outA, unionSize);
            testCase.verifySize(outB, unionSize);
            both = outA > 0 & outB > 0;
            testCase.verifyEqual(outB(both), outA(both), 'Corrected channels are not aligned.');
            testCase.verifyEqual(nnz(outA), prod(testCase.ImSize), 'Channel A content is missing from the output.');
        end
    end

    methods (Access = private)
        function [dataPath, scene, base] = shiftedChannels(testCase, offsets)
            % channel c voxel i = scene voxel i + offsets(c, :) + base
            dataPath = testCase.createTempFolder();
            base = max(0, -min(offsets, [], 1));
            sceneSize = testCase.ImSize + max(offsets, [], 1) + base;
            scene = uint16(pktest.texture(sceneSize, 1.5, [500, 5000]));
            cams = {'CamA', 'CamB'};
            for c = 1 : size(offsets, 1)
                s = 1 + offsets(c, :) + base;
                t = s + testCase.ImSize - 1;
                writetiff(scene(s(1) : t(1), s(2) : t(2), s(3) : t(3)), [dataPath, sprintf(testCase.Fsn, cams{c}), '.tif']);
            end
        end

        function [] = runCorrection(~, dataPath, mode, varargin)
            XR_chromatic_shift_correction_data_wrapper({dataPath}, 'mode', mode, ...
                'channelPatterns', {'CamA', 'CamB'}, 'parseCluster', false, 'masterCompute', true, varargin{:});
        end

        function [outA, outB] = readOutputs(testCase, dataPath)
            outA = readtiff([dataPath, 'Chromatic_Shift_Corrected/', sprintf(testCase.Fsn, 'CamA'), '.tif']);
            outB = readtiff([dataPath, 'Chromatic_Shift_Corrected/', sprintf(testCase.Fsn, 'CamB'), '.tif']);
        end
    end
end
