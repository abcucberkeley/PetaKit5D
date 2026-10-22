classdef XRDeskewRotateFrameTest < pktest.TestCase
    % XR_deskewRotateFrame: single file deskew/rotation with outputs on disk.
    % The in-memory functions are tested in DeskewRotateFrame3DTest; here
    % the files must contain exactly what those functions return.

    properties (Constant)
        SkewAngle = 32.45
        XYPixelSize = 0.108
        Dz = 0.2
        RawSize = [64, 96, 80]
        Fsn = 'Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t'
    end

    methods (Test, TestTags = {'Integration'})
        function combinedWritesDSRAndMIP(testCase)
            [dataPath, fn, frame] = testCase.writeFrame();
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'skewAngle', testCase.SkewAngle, 'reverse', true);

            expected = deskewRotateFrame3D(frame, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                'reverse', true, 'save16bit', true);
            dsr = readtiff([dataPath, 'DSR/', testCase.Fsn, '.tif']);
            testCase.verifyEqual(dsr, expected);
            testCase.verifyEqual(readtiff([dataPath, 'DSR/MIPs/', testCase.Fsn, '_MIP_z.tif']), max(expected, [], 3));
            testCase.verifyFalse(isfolder([dataPath, 'DS']), 'Combined deskew/rotation should not write DS.');
        end

        function saveZarrMatchesTiff(testCase)
            [dataPath, fn] = testCase.writeFrame();
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'reverse', true);
            tifPath = [dataPath, 'tif/'];
            mkdir(tifPath);
            movefile([dataPath, 'DSR'], tifPath);

            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'reverse', true, 'saveZarr', true);
            testCase.verifyEqual(readzarr([dataPath, 'DSR/', testCase.Fsn, '.zarr']), ...
                readtiff([tifPath, 'DSR/', testCase.Fsn, '.tif']));
            testCase.verifyEqual(readzarr([dataPath, 'DSR/MIPs/', testCase.Fsn, '_MIP_z.zarr']), ...
                readtiff([tifPath, 'DSR/MIPs/', testCase.Fsn, '_MIP_z.tif']));
        end

        function save16bitFalseWritesSingle(testCase)
            [dataPath, fn, frame] = testCase.writeFrame();
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'reverse', true, 'save16bit', false);
            dsr = readtiff([dataPath, 'DSR/', testCase.Fsn, '.tif']);
            testCase.verifyClass(dsr, 'single');
            expected = deskewRotateFrame3D(frame, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, 'reverse', true);
            testCase.verifyEqual(dsr, expected);
        end

        function separatePathWritesDeskewedAndRotated(testCase)
            [dataPath, fn, frame] = testCase.writeFrame();
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'reverse', true, 'DSRCombined', false, 'rotate', true);

            ds = readtiff([dataPath, 'DS/', testCase.Fsn, '.tif']);
            testCase.verifyEqual(ds, uint16(deskewFrame3D(single(frame), testCase.SkewAngle, testCase.Dz, ...
                testCase.XYPixelSize, true, 'Crop', false)));
            testCase.verifyTrue(isfile([dataPath, 'DS/MIPs/', testCase.Fsn, '_MIP_z.tif']));

            % same geometry as the combined transform (beads within 0.15 voxel)
            dsr = readtiff([dataPath, 'DSR/', testCase.Fsn, '.tif']);
            combined = deskewRotateFrame3D(frame, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                'reverse', true, 'save16bit', true);
            testCase.verifyEqual(size(dsr), size(combined));
            [~, expectedPos] = testCase.beadPositions(combined);
            testCase.verifyEqual(pktest.centroids(single(dsr) - 100, expectedPos, 6), ...
                pktest.centroids(single(combined) - 100, expectedPos, 6), 'AbsTol', 0.15);
        end

        function deskewOnlyWritesOnlyDS(testCase)
            [dataPath, fn] = testCase.writeFrame();
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'reverse', true, 'DSRCombined', false, 'rotate', false);
            testCase.verifyTrue(isfile([dataPath, 'DS/', testCase.Fsn, '.tif']));
            testCase.verifyFalse(isfolder([dataPath, 'DSR']));
        end

        function objectiveScanSeparatePathRuns(testCase)
            % regression for 49863e2: theta was undefined for objective
            % scan with DSRCombined=false
            [dataPath, fn] = testCase.writeFrame();
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'objectiveScan', true, 'DSRCombined', false, 'rotate', true);
            theta = testCase.SkewAngle * pi / 180;
            sz = testCase.RawSize;
            q = testCase.Dz / testCase.XYPixelSize;
            outSize = round([sz(1), sz(2) * cos(theta) + sz(3) * q * sin(theta), sz(3) * q * cos(theta) + sz(2) * sin(theta)]);
            testCase.verifyEqual(getImageSize([dataPath, 'DSR/', testCase.Fsn, '.tif']), outSize);
        end

        function largeZSeparatePathWritesDSMIP(testCase)
            % regression for fe83bf6: with more than 1000 slices the
            % separate path deskews in blocks (blockedImage) and the DS MIP
            % code used a misspelled property
            dataPath = testCase.createTempFolder();
            frame = uint16(100 + 1000 * rand(16, 16, 1001));
            fn = [dataPath, testCase.Fsn, '.tif'];
            writetiff(frame, fn);
            dz = 0.01;
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, dz, 'DSRCombined', false, 'rotate', false);

            expected = uint16(deskewFrame3D(single(frame), testCase.SkewAngle, dz, testCase.XYPixelSize, false));
            testCase.verifyEqual(readtiff([dataPath, 'DS/', testCase.Fsn, '.tif']), expected);
            testCase.verifyEqual(readtiff([dataPath, 'DS/MIPs/', testCase.Fsn, '_MIP_z.tif']), max(expected, [], 3));
        end

        function flatFieldCorrectionIsAppliedFirst(testCase)
            [dataPath, fn, frame] = testCase.writeFrame();
            bg = zeros(testCase.RawSize(1 : 2), 'single') + 100;
            ff = bg + single(linspace(500, 1000, testCase.RawSize(2)));
            writetiff(ff, [dataPath, 'ff.tif']);
            writetiff(bg, [dataPath, 'bg.tif']);
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'reverse', true, 'FFCorrection', true, ...
                'FFImage', [dataPath, 'ff.tif'], 'backgroundImage', [dataPath, 'bg.tif'], 'lowerLimit', 0.4);

            corrected = XR_LSFlatFieldCorrection(frame, ff, bg, 'LowerLimit', 0.4);
            expected = deskewRotateFrame3D(corrected, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                'reverse', true, 'save16bit', true);
            testCase.verifyEqual(readtiff([dataPath, 'DSR/', testCase.Fsn, '.tif']), expected);
        end

        function inputBboxAndFlipAreAppliedFirst(testCase)
            [dataPath, fn, frame] = testCase.writeFrame();
            bbox = [9, 5, 3, 56, 90, 72];
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'reverse', true, 'inputBbox', bbox, 'flipZstack', true);

            cropped = flip(frame(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6)), 3);
            expected = deskewRotateFrame3D(cropped, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                'reverse', true, 'save16bit', true);
            testCase.verifyEqual(readtiff([dataPath, 'DSR/', testCase.Fsn, '.tif']), expected);
        end

        function inputAxisOrderXYZIsPermutedToYXZ(testCase)
            [dataPath, ~, frame] = testCase.writeFrame();
            fnXYZ = [dataPath, 'xyz/', testCase.Fsn, '.tif'];
            writetiff(permute(frame, [2, 1, 3]), fnXYZ);
            XR_deskewRotateFrame(fnXYZ, testCase.XYPixelSize, testCase.Dz, 'reverse', true, 'inputAxisOrder', 'xyz');

            expected = deskewRotateFrame3D(frame, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                'reverse', true, 'save16bit', true);
            testCase.verifyEqual(readtiff([dataPath, 'xyz/DSR/', testCase.Fsn, '.tif']), expected);
        end

        function outputAxisOrderZYXPermutesTheResult(testCase)
            % a 3-cycle, which a self-inverse order like 'xyz' can't tell
            % apart from its inverse
            [dataPath, fn, frame] = testCase.writeFrame();
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'reverse', true, 'outputAxisOrder', 'zyx');

            expected = deskewRotateFrame3D(frame, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                'reverse', true, 'save16bit', true);
            testCase.verifyEqual(readtiff([dataPath, 'DSR/', testCase.Fsn, '.tif']), permute(expected, [3, 1, 2]));
        end

        function existingOutputIsNotRecomputed(testCase)
            [dataPath, fn, frame] = testCase.writeFrame();
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'reverse', true);
            dsrFn = [dataPath, 'DSR/', testCase.Fsn, '.tif'];
            first = readtiff(dsrFn);

            writetiff(frame * 2, fn);
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'reverse', true);
            testCase.verifyEqual(readtiff(dsrFn), first);
        end

        function zStageScanSeparateMatchesCombined(testCase)
            % the separate path used to size its output with skewAngle
            % instead of the 90 - skewAngle it rotates by
            [dataPath, fn] = testCase.writeFrame();
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'zStageScan', true);
            combinedSize = getImageSize([dataPath, 'DSR/', testCase.Fsn, '.tif']);
            rmdir([dataPath, 'DSR'], 's');
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'zStageScan', true, 'DSRCombined', false, 'rotate', true);
            testCase.verifyEqual(getImageSize([dataPath, 'DSR/', testCase.Fsn, '.tif']), combinedSize);
        end

        function inputAxisOrderZYXIsPermutedToYXZ(testCase)
            % a 3-cycle, which used to be permuted the wrong way ('zyx' input
            % ended up as 'xzy'; see AxisOrderMappingTest)
            [dataPath, ~, frame] = testCase.writeFrame();
            fnZYX = [dataPath, 'zyx/', testCase.Fsn, '.tif'];
            writetiff(permute(frame, [3, 1, 2]), fnZYX);
            XR_deskewRotateFrame(fnZYX, testCase.XYPixelSize, testCase.Dz, 'reverse', true, 'inputAxisOrder', 'zyx');

            expected = deskewRotateFrame3D(frame, testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                'reverse', true, 'save16bit', true);
            testCase.verifyEqual(readtiff([dataPath, 'zyx/DSR/', testCase.Fsn, '.tif']), expected);
        end

        function privateZarrResultsAreWrittenInPlace(testCase)
            % without a given uuid (the data wrapper) the temporary DS and DSR
            % results are private to the call and their chunks are written
            % without cpp-zarr's temporary files ('uuid' false), except with
            % the short uuids of Windows; with the pipeline's uuid, which
            % copies of a task share, they are kept
            spy = testCase.applyFixture(pktest.FunctionSpyFixture('parallelWriteZarr'));
            runs = {{}, {'DSRCombined', false, 'rotate', false}, ...
                {'DSRCombined', false, 'rotate', false, 'save16bit', false}, ...
                {'uuid', '01234567-89ab-cdef-0123-456789abcdef'}};
            inPlace = [true, true, true, false] & numel(get_uuid()) >= 32;
            for r = 1 : numel(runs)
                [dataPath, fn] = testCase.writeFrame();
                nCalls = numel(spy.getCalls());
                XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'reverse', true, 'saveZarr', true, runs{r}{:});
                calls = spy.getCalls();
                % the result is the write into DS/ or DSR/ (not MIPs/)
                isResult = cellfun(@(c) startsWith(c{1}, {[dataPath, 'DS/'], [dataPath, 'DSR/']}) ...
                    && ~contains(c{1}, 'MIPs'), calls(nCalls + 1 : end));
                testCase.assertEqual(nnz(isResult), 1, sprintf('run %d', r));
                uuidOption = spy.getOption(nCalls + find(isResult), 'uuid');
                if inPlace(r)
                    testCase.verifyEqual(uuidOption, false, sprintf('run %d', r));
                else
                    testCase.verifyEmpty(uuidOption, sprintf('run %d', r));
                end
            end
        end
    end

    methods (Test, TestTags = {'Integration', 'KnownIssue'})
        function backgroundRemovalRuns(testCase)
            % KNOWN ISSUE: BKRemoval=true calls XR_CameraBackgroundRemoval
            % (XR_deskewRotateFrame.m:213), which does not exist in the repo.
            [dataPath, fn] = testCase.writeFrame();
            writetiff(zeros(testCase.RawSize(1 : 2), 'single') + 100, [dataPath, 'bg.tif']);
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'reverse', true, 'BKRemoval', true, ...
                'backgroundImage', [dataPath, 'bg.tif']);
            testCase.verifyTrue(isfile([dataPath, 'DSR/', testCase.Fsn, '.tif']));
        end

        function objectiveScanSeparatePathIsNotDeskewed(testCase)
            % KNOWN ISSUE (may be intended): for objective scan the separate path still
            % deskews (deskewFrame3D is called without an objective scan
            % option in XR_deskewRotateFrame.m:235) before rotating, so beads
            % end up sheared relative to a pure rotation.
            [dataPath, fn, frame] = testCase.writeFrame();
            XR_deskewRotateFrame(fn, testCase.XYPixelSize, testCase.Dz, 'objectiveScan', true, 'DSRCombined', false, 'rotate', true);
            dsr = readtiff([dataPath, 'DSR/', testCase.Fsn, '.tif']);
            % the cubic path of deskewRotateFrame3D is a correct pure rotation
            expected = deskewRotateFrame3D(single(frame), testCase.SkewAngle, testCase.Dz, testCase.XYPixelSize, ...
                'objectiveScan', true, 'interpMethod', 'cubic');
            [~, expectedPos] = testCase.beadPositions(expected);
            testCase.verifyEqual(pktest.centroids(single(dsr) - 100, expectedPos, 6), ...
                pktest.centroids(expected - 100, expectedPos, 6), 'AbsTol', 0.5);
        end
    end

    methods (Access = private)
        function [dataPath, fn, frame] = writeFrame(testCase)
            % uint16 beads on a background of 100
            dataPath = testCase.createTempFolder();
            [xr, zr] = ndgrid([30, 48, 66], [12, 26, 40, 54, 68]);
            raw = [8 + 3.4 * mod((0 : 14)' * 7, 15), xr(:), zr(:)];
            frame = uint16(pktest.beads(testCase.RawSize, raw, 2, 1000, 100));
            fn = [dataPath, testCase.Fsn, '.tif'];
            writetiff(frame, fn);
        end

        function [raw, pos] = beadPositions(~, out)
            % approximate output bead positions from the brightest local
            % maxima of a reference output
            mask = imregionalmax(out) & out > 0.5 * max(out(:));
            [y, x, z] = ind2sub(size(out), find(mask));
            raw = [];
            pos = [y, x, z];
        end
    end
end
