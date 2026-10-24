classdef PSFAnalysisTest < pktest.TestCase
    % XR_psf_analysis_wrapper ("PSF analysis" in demo_useful_tools.m):
    % deskew the PSF, then plot it and its OTF against a Richards & Wolf
    % widefield PSF. Even for a 64 x 64 x 101 PSF this needs ~22 GB of RAM
    % (the plots are simulated on a 1001^3 double grid) and ~80 s, so it is
    % skipped on machines with less than 32 GB (e.g., GitHub runners).

    properties (TestParameter)
        % {PSF size, bead center, sigma [yx, z], z pixel size (nm), PSFsubpix,
        % optional second bead {center, sigma, amplitude}}: centered, padded in
        % z, larger than the 1001^3 grid in x and y after resizing (the
        % grid's data region wraps around when it is centered), and a sharp
        % bead with a broad, slightly dimmer one that ends up outside the
        % grid but has the larger resized maximum (the normalization)
        psfCase = struct( ...
            'centered', {{[64, 64, 101], [32, 32, 51], [1.6, 3], 160, [128, 128, 63]}}, ...
            'paddedZ', {{[64, 60, 21], [30, 33, 11], [1.6, 2], 300, [128, 128, 35]}}, ...
            'largerThanGrid', {{[380, 370, 31], [185, 190, 15], [2, 3], 160, [128, 128, 63]}}, ...
            'maxOutsideGrid', {{[360, 360, 31], [180, 180, 16], [0.5, 0.8], 160, [128, 128, 63], ...
                {[180, 355, 16], [4, 3], 1990}}})
    end

    methods (Test, TestTags = {'Integration', 'Slow'})
        function plotDataMatchWholeGrid(testCase, psfCase)
            % the PSF and OTF on the 1001^3 simulation grid are kept as
            % their data regions (the grid is zero elsewhere); the slices
            % must equal those of the whole grid (the reference needs it)
            testCase.assumeGreaterThanOrEqual(testCase.memoryGB(), 48, 'The whole-grid reference needs ~40 GB of RAM.');
            [sz, center, sigma, zpixsize, PSFsubpix] = psfCase{1 : 5};
            dataPath = testCase.createTempFolder();
            fn = [dataPath, 'psf.tif'];
            im = pktest.beads(sz, center, sigma, 2000, 100) + 5 * rand(sz);
            if numel(psfCase) > 5
                im = im + pktest.beads(sz, psfCase{6}{:}, 0);
            end
            writetiff(uint16(im), fn);
            ref = PSFAnalysisTest.referencePSFSlices(fn, 108, zpixsize, 1.33, 488, PSFsubpix, 1.5);
            names = {'xy_exp_PSF', 'xz_exp_PSF', 'yz_exp_PSF', 'xy_exp_OTF', 'xz_exp_OTF', 'yz_exp_OTF'};
            out = cell(1, numel(names));
            [out{:}] = Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update(fn, 'test', 108, zpixsize, 1.0, 1.33, 488, 515, ...
                PSFsubpix, 0.5, 1.5, false);
            for i = 1 : numel(names)
                testCase.verifyEqual(out{i}, ref.(names{i}), names{i});
            end
        end

        function analysisWritesFiguresAndInfo(testCase)
            testCase.assumeGreaterThanOrEqual(testCase.memoryGB(), 32, 'PSF analysis needs ~22 GB of RAM.');
            dataPath = testCase.createTempFolder();
            psfPath = [dataPath, 'PSF/'];
            mkdir([psfPath, 'RW_PSFs/']);
            fsn = '488_CamA_ch0';
            writetiff(uint16(pktest.beads([64, 64, 101], [32, 32, 51], [1.6, 3], 2000, 100)), [psfPath, fsn, '.tif']);
            rwFn = [psfPath, 'RW_PSFs/rw.tif'];
            writetiff(single(pktest.beads([128, 128, 101], [64, 64, 51], [2, 6], 1, 0)), rwFn);
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);

            XR_psf_analysis_wrapper({psfPath}, 'xyPixelSize', 0.108, 'dz', 0.3, 'skewAngle', 32.45, ...
                'channelPatterns', {'CamA_ch0'}, 'channels', 488, 'deskew', true, 'RWFn', {rwFn}, ...
                'sourceStr', 'test', 'visible', false, 'parseCluster', false, 'masterCompute', true);

            resultPath = [psfPath, 'DS/PSFAnalysis/'];
            % deskewed PSF (the demo's PSF is in skewed space)
            ds = readtiff([psfPath, 'DS/', fsn, '.tif']);
            testCase.verifyEqual(size(ds, [1, 3]), [64, 101]);
            testCase.verifyGreaterThan(size(ds, 2), 64);
            for f = {['wT_', fsn, '.png'], ['wT_', fsn, '.pdf'], ['comp_', fsn, '.png'], ['comp_', fsn, '.pdf'], ...
                    [fsn, '_infos.mat'], 'RW_info.mat'}
                testCase.verifyTrue(exist([resultPath, f{1}], 'file') == 2, sprintf('Missing %s.', f{1}));
            end

            % the job framework counts the task as done, and may cancel its
            % job, once the wT_ PNG exists; the PDF was printed after it
            pngFn = [resultPath, 'wT_', fsn, '.png'];
            testCase.verifyTrue(moves.wasMoved(pngFn), 'The PNG was not moved into place.');
            testCase.verifyTrue(ismember(['wT_', fsn, '.pdf'], moves.filesWhenMoved(pngFn)), 'The PDF was printed after the PNG.');
        end

        function detectionInfoIsMovedIntoPlace(testCase)
            % XR_psf_detection_and_analysis_wrapper's tasks count as done
            % once their info file exists, so it must not be saved in place
            dataPath = testCase.createTempFolder();
            fn = [dataPath, 'beads.tif'];
            writetiff(uint16(pktest.beads([128, 128, 101], [40, 40, 30; 90, 80, 70], [1.6, 3], 2000, 100)), fn);
            resultDir = [dataPath, 'Cropped/'];
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            XR_psf_detection_and_cropping(fn, resultDir, 'cropSize', [32, 32, 21], 'distThresh', [32, 32, 21], ...
                'prefix', 'beads');
            testCase.verifyTrue(moves.wasMoved([resultDir, 'beads.mat']), 'The detection info was saved in place.');
            a = load([resultDir, 'beads.mat']);
            testCase.verifyTrue(isfield(a, 'pcoords'));
        end
    end

    methods (Static, Access = private)
        function [ref] = referencePSFSlices(filenm, xypixsize, zpixsize, index, exc_lambda, PSFsubpix, bgFactor)
            % the PSF and OTF slices of Load_and_Plot_Exp_Overall_xzPSF_xzOTF_update
            % computed on the whole 1001^3 grid, as before only the data
            % regions were kept
            sim_pixsize = 0.1;
            sim_half_FOV = 50;
            PSF3Dexp = cast(readtiff(filenm), 'double');
            A = size(PSF3Dexp);
            backgnd = bgFactor * median(PSF3Dexp(PSF3Dexp > 0));
            PSF3Dexp = max(PSF3Dexp - backgnd, 0);
            size_xy = max([PSFsubpix(1), PSFsubpix(2), A(1), A(2)]);
            hsz = [ceil((size_xy - A(1) - 1) / 2), ceil((size_xy - A(2) - 1) / 2)];
            psize_1 = [hsz(1), hsz(2), 0];
            psize_2 = [size_xy - A(1) - hsz(1), size_xy - A(2) - hsz(2), 0];
            if A(3) < PSFsubpix(3)
                hsz = ceil((PSFsubpix(3) - A(3) - 1) / 2);
                psize_1(3) = hsz;
                psize_2(3) = abs(PSFsubpix(3) - A(3)) - hsz;
            end
            PSF3Dexp = padarray(padarray(PSF3Dexp, psize_1, 0, 'pre'), psize_2, 0, 'post');
            A = size(PSF3Dexp);
            [~, peakInd] = max(PSF3Dexp(:));
            [peaky, peakx, peakz] = ind2sub(A, peakInd);
            PSF3Dexp = circshift(PSF3Dexp, round((A + 1) / 2 - [peaky, peakx, peakz]));

            % PSF on the grid
            sim_xzpixsize = sim_pixsize .* (exc_lambda ./ index);
            magxy = xypixsize ./ sim_xzpixsize;
            magz = zpixsize ./ sim_xzpixsize;
            B = [round(magxy .* A(1)), round(magxy .* A(2)), round(magz .* A(3))];
            PSF3D = imresize3(PSF3Dexp, B, 'linear');
            PSF3D = PSF3D ./ max(PSF3D(:));
            B = size(PSF3D);
            array_size = 2 .* sim_half_FOV ./ sim_pixsize + 1;
            PSF_array = zeros(array_size, array_size, array_size);
            xyi = max(round(sim_half_FOV ./ sim_pixsize - B(1) ./ 2), 1);
            xyf = xyi + min(B(1), array_size) - 1;
            zi = max(round(sim_half_FOV ./ sim_pixsize - B(3) ./ 2), 1);
            zf = zi + min(B(3), array_size) - 1;
            Bc = round((B + 1) / 2);
            lhsz = ceil((min(B, array_size) - 1) / 2);
            rhsz = min(B, array_size) - lhsz - 1;
            PSF_array(xyi : xyf, xyi : xyf, zi : zf) = PSF3D(Bc(1) - lhsz(1) : Bc(1) + rhsz(1), ...
                Bc(2) - lhsz(2) : Bc(2) + rhsz(2), Bc(3) - lhsz(3) : Bc(3) + rhsz(3));
            PSF3D = PSF_array;
            clear PSF_array;
            A = size(PSF3D);
            PSF3D_1 = imgaussfilt3(PSF3D, 2);
            [~, peakInd] = max(PSF3D_1(:));
            clear PSF3D_1;
            [cy, cx, cz] = ind2sub(A, peakInd);
            if abs((A(1) + 1) / 2 - cy) > 0.5 || abs((A(2) + 1) / 2 - cx) > 0.5 || abs((A(3) + 1) / 2 - cz) > 0.5
                PSF3D = circshift(PSF3D, round((A + 1) / 2 - [cy, cx, cz]));
            end
            PSFxy = squeeze(PSF3D(:, :, (A(3) + 1) ./ 2));
            [~, peakypix] = max(sum(PSFxy ./ max(max(PSFxy)), 1));
            maxval = zeros(1, 11);
            for yoffset = -5 : 5
                maxval(yoffset + 6) = max(max(squeeze(PSF3D(:, peakypix + yoffset, :))'));
            end
            [~, maxplane] = max(maxval);
            yz = squeeze(PSF3D(:, peakypix + maxplane - 6, :))';
            hsz = round((size(yz) + 1) / 2);
            ref.yz_exp_PSF = yz(hsz(1) - 50 : hsz(1) + 50, hsz(2) - 50 : hsz(2) + 50);
            [~, peakInd] = max(PSF3D(:));
            [maxyplane, ~, maxzplane] = ind2sub(size(PSF3D), peakInd);
            xy = squeeze(PSF3D(:, :, maxzplane));
            hsz = round((size(xy) + 1) / 2);
            ref.xy_exp_PSF = xy(hsz(1) - 50 : hsz(1) + 50, hsz(2) - 50 : hsz(2) + 50);
            xz = squeeze(PSF3D(maxyplane, :, :))';
            hsz = round((size(xz) + 1) / 2);
            ref.xz_exp_PSF = xz(hsz(1) - 50 : hsz(1) + 50, hsz(2) - 50 : hsz(2) + 50);
            clear PSF3D;

            % OTF on the grid
            C = size(PSF3Dexp);
            pixi = max(round((C - PSFsubpix) ./ 2), 1);
            PSFsubset = PSF3Dexp(pixi(1) : (pixi(1) + PSFsubpix(1) - 1), pixi(2) : (pixi(2) + PSFsubpix(2) - 1), ...
                pixi(3) : (pixi(3) + PSFsubpix(3) - 1));
            OTF3Dexp = abs(fftshift(fftn(PSFsubset)));
            OTF3Dexp = OTF3Dexp ./ max(max(max(OTF3Dexp)));
            C = size(OTF3Dexp);
            sim_half_FOVnm = sim_half_FOV .* (exc_lambda ./ index);
            magxy = (1 ./ (C(1) .* xypixsize ./ 2)) ./ (1 ./ sim_half_FOVnm);
            magz = (1 ./ (C(3) .* zpixsize ./ 2)) ./ (1 ./ sim_half_FOVnm);
            D = [round(magxy .* C(1)), round(magxy .* C(2)), round(magz .* C(3))];
            [x, y, z] = ndgrid(1 : C(1), 1 : C(2), 1 : C(3));
            [xq, yq, zq] = ndgrid(1 : (C(1) ./ D(1)) : C(1), 1 : (C(2) ./ D(2)) : C(2), 1 : (C(3) ./ D(3)) : C(3));
            OTF3D = interpn(x, y, z, OTF3Dexp, xq, yq, zq);
            OTF3D = OTF3D ./ max(max(max(OTF3D)));
            D = size(OTF3D);
            OTF_array = zeros(array_size, array_size, array_size);
            xyi = max(round(sim_half_FOV ./ sim_pixsize - D(1) ./ 2), 1);
            xyf = xyi + min(D(1), array_size) - 1;
            zi = max(round(sim_half_FOV ./ sim_pixsize - D(3) ./ 2), 1);
            zf = zi + min(D(3), array_size) - 1;
            Dc = round((D + 1) / 2);
            lhsz = ceil((min(D, array_size) - 1) / 2);
            rhsz = min(D, array_size) - lhsz - 1;
            OTF_array(xyi : xyf, xyi : xyf, zi : zf) = OTF3D(Dc(1) - lhsz(1) : Dc(1) + rhsz(1), ...
                Dc(2) - lhsz(2) : Dc(2) + rhsz(2), Dc(3) - lhsz(3) : Dc(3) + rhsz(3));
            B = size(OTF_array);
            ref.yz_exp_OTF = squeeze(OTF_array(:, (B(2) + 1) ./ 2, :))';
            ref.xy_exp_OTF = squeeze(OTF_array(:, :, (B(2) + 1) ./ 2));
            ref.xz_exp_OTF = squeeze(OTF_array((B(1) + 1) ./ 2, :, :))';
        end
    end

    methods (Access = private)
        function [gb] = memoryGB(~)
            % memory available to this MATLAB: the Slurm allocation if any,
            % otherwise the machine's memory
            slurmMB = str2double(getenv('SLURM_MEM_PER_NODE'));
            if ~isnan(slurmMB)
                gb = slurmMB / 1024;
            elseif ispc
                m = memory;
                gb = m.MemAvailableAllArrays / 1024 ^ 3;
            elseif exist('/proc/meminfo', 'file')
                kb = regexp(fileread('/proc/meminfo'), 'MemTotal:\s+(\d+)', 'tokens', 'once');
                gb = str2double(kb{1}) / 1024 ^ 2;
            else
                gb = getSystemMemory();
            end
        end
    end
end
