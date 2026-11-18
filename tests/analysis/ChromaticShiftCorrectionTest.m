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

        function psfBoundsAreInYXZ(testCase)
            % the bounds (cropLength +- maxOffset, in y, x, z) were turned into
            % x, y, z for normxcorr3_max_shift, which takes y, x, z: y was
            % bounded by the x values and x by the y values. A bead 5 voxels
            % off in y and 1 in x is found within bounds of 6 in y and 2 in x
            % as within 6 in all axes
            psfPath = testCase.createTempFolder();
            psfFns = {[psfPath, 'psf_CamA.tif'], [psfPath, 'psf_CamB.tif']};
            beads = [20, 20, 12; 25, 21, 12];
            for c = 1 : 2
                writetiff(uint16(pktest.beads([40, 40, 24], beads(c, :), [1.5, 2], 2000, 100)), psfFns{c});
            end
            wide = XR_chromatic_shift_estimation(psfFns, 'maxOffset', [6, 6, 6]);
            testCase.assertEqual(abs(wide(2, :)), [5, 1, 0]);
            testCase.verifyEqual(XR_chromatic_shift_estimation(psfFns, 'maxOffset', [6, 2, 2]), wide);
            % (a char is one path: its characters were taken as the files)
            testCase.verifyError(@() XR_chromatic_shift_estimation(psfFns{1}), ?MException);
            try
                XR_chromatic_shift_estimation(psfFns{1});
            catch ME
                testCase.verifyEqual(ME.message, 'At least two PSF files are needed!');
            end
        end

        function zarrDataWithOffsetsFromPSFs(testCase)
            % the offsets measured from (tiff) PSFs set zarrFile false for
            % them, and it stayed false for the crops of the zarr data, which
            % failed
            [tifPath, ~, ~] = testCase.shiftedChannels(testCase.Offset);
            dataPath = testCase.zarrCopy(tifPath);
            psfFns = testCase.writePSFs(dataPath);

            testCase.runCorrection(dataPath, 'valid', 'chromaticOffset', [], 'psfFullpaths', psfFns, ...
                'maxOffset', [8, 8, 6], 'zarrFile', true);

            [outA, outB] = testCase.readOutputs(dataPath);
            testCase.verifyEqual(outB, outA, 'Channels corrected with the measured offsets are not aligned.');
        end

        function largeZarrFilesAreCroppedInBatches(testCase)
            % largeFile only sized the jobs (with a batch's memory), and the
            % files were cropped in memory; large zarr files are cropped in
            % batches (XR_crop_zarr) into zarr results, equal to the crops of
            % the regular path
            [tifPath, ~, ~] = testCase.shiftedChannels(testCase.Offset);
            dataPath = testCase.zarrCopy(tifPath);
            args = {'chromaticOffset', testCase.Offset, 'zarrFile', true, 'batchSize', [16, 16, 8], 'blockSize', [8, 8, 8]};
            for mode = {'valid', 'full'}
                testCase.runCorrection(dataPath, mode{1}, args{:}, 'saveZarr', true, 'resultDirName', ['regular_', mode{1}]);
                testCase.runCorrection(dataPath, mode{1}, args{:}, 'largeFile', true, 'resultDirName', ['large_', mode{1}]);
                for cam = {'CamA', 'CamB'}
                    fn = [sprintf(testCase.Fsn, cam{1}), '.zarr'];
                    large = [dataPath, 'large_', mode{1}, '/', fn];
                    testCase.verifyEqual(readzarr(large), readzarr([dataPath, 'regular_', mode{1}, '/', fn]), ...
                        sprintf('%s, %s', mode{1}, cam{1}));
                    % (XR_crop_zarr publishes its result with its MIPs)
                    testCase.verifyTrue(isfile([dataPath, 'large_', mode{1}, '/MIPs/', sprintf(testCase.Fsn, cam{1}), ...
                        '_MIP_z.tif']), sprintf('%s, %s: not cropped in batches.', mode{1}, cam{1}));
                end
            end
            % (XR_crop_zarr pads with 0)
            testCase.verifyError(@() testCase.runCorrection(dataPath, 'full', args{:}, 'largeFile', true, ...
                'padValue', 7, 'resultDirName', 'large_pad'), 'XR_chromatic_shift_correction_data_wrapper:padValue');
        end

        function offsetsHaveARowPerChannel(testCase)
            % without offsets (or PSFs) the boxes failed on a size mismatch,
            % with a row too few on an index, and a row too many changed the
            % overlap
            dataPath = testCase.shiftedChannels(testCase.Offset);
            id = 'XR_chromatic_shift_correction_data_wrapper:';
            testCase.verifyError(@() testCase.runCorrection(dataPath, 'valid', 'chromaticOffset', []), [id, 'noOffsets']);
            testCase.verifyError(@() testCase.runCorrection(dataPath, 'valid', 'chromaticOffset', testCase.Offset(2, :)), ...
                [id, 'offsetSize']);
            testCase.verifyError(@() testCase.runCorrection(dataPath, 'valid', 'chromaticOffset', [testCase.Offset; 0, 9, 0]), ...
                [id, 'offsetSize']);
        end

        function modeIsOneOfThreeInAnyCase(testCase)
            % another mode cropped each channel's box shifted against its
            % offset (without error), and the switch took the case into
            % account: 'Valid' is 'valid'
            dataPath = testCase.shiftedChannels(testCase.Offset);
            testCase.verifyError(@() testCase.runCorrection(dataPath, 'none', 'chromaticOffset', testCase.Offset), ...
                'MATLAB:InputParser:ArgumentFailedValidation');
            testCase.runCorrection(dataPath, 'Valid', 'chromaticOffset', testCase.Offset, 'resultDirName', 'upper');
            testCase.runCorrection(dataPath, 'valid', 'chromaticOffset', testCase.Offset);
            fn = [sprintf(testCase.Fsn, 'CamB'), '.tif'];
            testCase.verifyEqual(readtiff([dataPath, 'upper/', fn]), readtiff([dataPath, 'Chromatic_Shift_Corrected/', fn]));
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

        function sameModeKeepsImageSize(testCase)
            % 'same' mode used to error (neworigin undefined), flip the offsets and
            % shift the output by half the image by default
            [dataPath, ~, ~] = testCase.shiftedChannels(testCase.Offset);

            testCase.runCorrection(dataPath, 'same', 'chromaticOffset', testCase.Offset);

            [outA, outB] = testCase.readOutputs(dataPath);
            testCase.verifySize(outA, testCase.ImSize);
            testCase.verifySize(outB, testCase.ImSize);
            both = outA > 0 & outB > 0;
            testCase.verifyEqual(outB(both), outA(both), 'Corrected channels are not aligned.');
        end

        function fullModeContainsBothChannels(testCase)
            % 'full' mode used to size the output from the bounding box starts and
            % flip the offsets
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

        function [zarrPath] = zarrCopy(testCase, tifPath)
            % the channels' tiff files of tifPath as zarr files in a folder of
            % their own
            zarrPath = testCase.createTempFolder();
            for cam = {'CamA', 'CamB'}
                fsn = sprintf(testCase.Fsn, cam{1});
                writezarr(readtiff([tifPath, fsn, '.tif']), [zarrPath, fsn, '.zarr'], 'blockSize', [16, 16, 8]);
            end
        end

        function [psfFns] = writePSFs(testCase, dataPath)
            % bead images (tiff) of the same bead for the two channels
            psfPath = [dataPath, 'PSF/'];
            mkdir(psfPath);
            psfFns = {[psfPath, 'psf_CamA.tif'], [psfPath, 'psf_CamB.tif']};
            bead = [20, 20, 12];
            for c = 1 : 2
                writetiff(uint16(pktest.beads([40, 40, 24], bead - testCase.Offset(c, :), [1.5, 2], 2000, 100)), psfFns{c});
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
