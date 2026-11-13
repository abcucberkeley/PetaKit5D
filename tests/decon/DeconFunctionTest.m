classdef DeconFunctionTest < pktest.TestCase
    % deconvolution building blocks (decon_lucy_function,
    % decon_lucy_omw_function, omw_backprojector_generation, psf_gen_new,
    % decon_psf2otf and RLdecon in memory) on synthetic beads blurred with a
    % known Gaussian PSF.

    properties (Constant)
        ImSize = [48, 48, 40]
        Background = 100
        PSFSigma = [1.5, 2]
        % blurred beads are Gaussian objects of this sigma convolved with the PSF
        ObjectSigma = 0.8
        WienerAlpha = 0.005
    end

    properties (TestParameter)
        method = {'simplified', 'omw'}
        inputClass = {'uint16', 'single'}
        anyInputClass = {'uint16', 'single', 'double'}
        % post-processing options of the decon functions
        deconOptions = struct('uint16Output', {{'save16bit', true}}, 'singleOutput', {{'save16bit', false}}, ...
            'dampErosionOffsetScale', {{'save16bit', false, 'dampFactor', 1.2, 'EdgeErosion', 4, 'deconOffset', 3, 'scaleFactor', 1.5}}, ...
            'bbox', {{'save16bit', true, 'bbox', [3, 4, 2, 44, 45, 36]}})
    end

    properties (SetAccess = private)
        PSF
        BackProjector
    end

    methods (TestClassSetup)
        function createPSFAndBackProjector(testCase)
            testCase.PSF = pktest.deconPSF([31, 31, 31], testCase.PSFSigma);
            [~, testCase.BackProjector] = evalc('omw_backprojector_generation(testCase.PSF, testCase.WienerAlpha, false)');
        end
    end

    methods (Test, TestTags = {'Unit'})
        function deltaPSFReturnsBackgroundSubtractedInput(testCase, method)
            % with a delta PSF the forward model is the identity, so the
            % result is the background subtracted input, also below background
            im = uint16(randi([50, 1000], 24, 24, 16));
            delta = zeros(3, 3, 3, 'single');
            delta(2, 2, 2) = 1;
            out = testCase.runDecon(method, im, delta, delta, 5, 'Background', testCase.Background, 'save16bit', true);
            testCase.verifyEqual(out, uint16(max(single(im) - testCase.Background, 0)));
        end

        function fluxIsConserved(testCase)
            % the RL update J <- J .* H'(I ./ HJ) with a normalized PSF keeps
            % sum(J) == sum(I) (not for OMW, whose back projector is not H')
            im = testCase.blurredBeads();
            out = testCase.runDecon('simplified', im, testCase.PSF, [], 10, 'Background', testCase.Background, 'save16bit', false);
            testCase.verifyEqual(sum(out, 'all', 'double'), sum(max(double(im) - testCase.Background, 0), 'all'), 'RelTol', 1e-4);
        end

        function regionsAtBackgroundStayZero(testCase, method)
            % RL updates are multiplicative, so voxels at or below the
            % background are exactly zero in the result
            im = testCase.blurredBeads();
            out = testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, 5, ...
                'Background', testCase.Background, 'save16bit', false);
            atBackground = im <= testCase.Background;
            testCase.assertGreaterThan(nnz(atBackground), 0.5 * numel(im));
            testCase.verifyEqual(nnz(out(atBackground)), 0);
            testCase.verifyGreaterThanOrEqual(min(out(:)), single(0));
        end

        function outputScalesLinearlyWithInput(testCase, method)
            % a floor of 10 keeps every voxel well above the FFT round-off
            im = testCase.blurredBeads([], [], 10);
            out1 = testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, 5, 'Background', 0, 'save16bit', false);
            out3 = testCase.runDecon(method, 3 * im, testCase.PSF, testCase.BackProjector, 5, 'Background', 0, 'save16bit', false);
            testCase.verifyImagesClose(out3, 3 * out1, 'relTol', 1e-5);
        end

        function pointSourceIsSharpened(testCase, method)
            % deconvolving a blurred bead raises its peak and narrows it in
            % every axis without moving it
            pos = [24.3, 23.6, 20.2];
            im = testCase.blurredBeads(pos, 1000, testCase.Background);
            nIter = 10;
            if strcmp(method, 'omw')
                nIter = 2;
            end
            out = testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, nIter, ...
                'Background', testCase.Background, 'save16bit', false);
            in = max(single(im) - testCase.Background, 0);

            testCase.verifyGreaterThan(max(out(:)), 1.5 * max(in(:)), 'Peak was not raised.');
            [cIn, varIn] = testCase.moments(in, pos);
            [cOut, varOut] = testCase.moments(out, pos);
            testCase.verifyEqual(cOut, cIn, 'AbsTol', 0.05, 'Bead moved.');
            testCase.verifyLessThan(varOut, 0.8 * varIn, sprintf('Bead was not narrowed: variance %s -> %s.', ...
                mat2str(varIn, 3), mat2str(varOut, 3)));
        end

        function edgeErosionZeroesTheBorderShell(testCase, method)
            % for images without zero voxels the outer EdgeErosion voxels
            % on every side are set to zero, and nothing inside
            im = uint16(pktest.texture(testCase.ImSize, 3, [200, 2000]));
            e = 3;
            out = testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, 3, ...
                'Background', 0, 'save16bit', false, 'EdgeErosion', e);

            shell = true(testCase.ImSize);
            shell(e + 1 : end - e, e + 1 : end - e, e + 1 : end - e) = false;
            testCase.verifyEqual(nnz(out(shell)), 0, 'The eroded border is not zero.');
            interior = out(e + 1 : end - e, e + 1 : end - e, e + 1 : end - e);
            testCase.verifyGreaterThan(nnz(interior) / numel(interior), 0.99, 'Erosion removed voxels inside the border.');
        end

        function dampFactorBoundsTheResult(testCase)
            d = 1.2;
            im = testCase.blurredBeads();
            out = testCase.runDecon('simplified', im, testCase.PSF, [], 10, 'Background', testCase.Background, ...
                'save16bit', false, 'dampFactor', d);
            bound = d * max(single(im) - testCase.Background, 0);
            testCase.verifyLessThanOrEqual(max(out - bound, [], 'all'), single(1e-3));
            % the bound is reached, i.e., the damping actually did something
            testCase.verifyGreaterThan(nnz(abs(out - bound) < 1e-3 & bound > 0), 0);
        end

        function dampFactorKeepsVoxelsAtTheBound(testCase, method)
            % a voxel exactly at dampFactor times the background-subtracted
            % input was doubled (both terms of the damping expression took
            % it). With 0 iterations the result is that input, and a factor
            % of 1 + 2^-30 (1 in single precision) puts every voxel at the
            % bound; on the default path (the mex where there is one) and
            % in the MATLAB loop
            im = testCase.blurredBeads();
            expected = max(single(im) - testCase.Background, 0);
            for useMex = {[], false}
                out = testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, 0, 'Background', testCase.Background, ...
                    'save16bit', false, 'dampFactor', 1 + 2 ^ -30, 'useMex', useMex{1});
                testCase.verifyEqual(out, expected, sprintf('useMex %s', mat2str(useMex{1})));
            end
        end

        function iterationCountIsZeroWithoutIterations(testCase)
            % decon_lucy_function's third output, the iterations run, was []
            % for 0 iterations (the loop variable of an empty loop)
            im = testCase.blurredBeads();
            for nIter = [0, 2]
                [~, ~, k] = decon_lucy_function(im, testCase.PSF, nIter, 'Background', testCase.Background, 'useGPU', false);
                testCase.verifyEqual(k, nIter);
            end
        end

        function debugOutputOfFewerThanThreeSteps(testCase, method)
            % RL's debug err_mat errored (istp undefined) with fewer than three
            % saveStep steps of iterations, e.g., 10 iterations with saveStep 5:
            % both methods keep the first floor(NUMIT / saveStep) rows and
            % write the estimates of the steps
            im = testCase.blurredBeads();
            for c = {[10, 5], [2, 5]}
                folder = sprintf('%sdebug_%d_%d/', testCase.createTempFolder(), c{1});
                args = {'Background', testCase.Background, 'useGPU', false, 'save16bit', false, 'debug', true, ...
                    'debug_folder', folder, 'saveStep', c{1}(2)};
                switch method
                    case 'simplified'
                        [~, ~, errMat] = evalc('decon_lucy_function(im, testCase.PSF, c{1}(1), args{:})');
                    case 'omw'
                        [~, ~, errMat] = evalc('decon_lucy_omw_function(im, testCase.PSF, testCase.BackProjector, c{1}(1), args{:})');
                end
                n = floor(c{1}(1) / c{1}(2));
                testCase.verifyEqual(size(errMat), [n, 4], mat2str(c{1}));
                testCase.verifyEqual(errMat(:, [1, 4]), [(1 : n)', ones(n, 1)], mat2str(c{1}));
                files = dir([folder, 'Iter_*.tif']);
                testCase.verifyEqual(reshape(sort({files.name}), 1, []), arrayfun(@(k) sprintf('Iter_%04d.tif', k), ...
                    c{1}(2) * (1 : n), 'unif', 0), mat2str(c{1}));
            end
        end

        function bboxCropsTheResult(testCase, method)
            im = testCase.blurredBeads();
            full = testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, 3, 'Background', testCase.Background, 'save16bit', true);
            bbox = [5, 6, 7, 40, 41, 30];
            out = testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, 3, 'Background', testCase.Background, ...
                'save16bit', true, 'bbox', bbox);
            testCase.verifyEqual(out, full(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6)));
        end

        function gpuRequestWithoutUsableGPUErrors(testCase, method)
            % useGPU true requires a GPU: with GPUs found but none usable
            % (e.g., a B200 in R2024b without forward compatibility, compute
            % capability 10.0) or none found, the decon errors with the
            % counts, CUDA_VISIBLE_DEVICES and MATLAB's reason, the error of
            % selecting GPU 1 (not tried with none found); GPU jobs used to
            % run on the CPU instead, about 60x slower. Also after a call
            % without useGPU cached a CPU OTF. Stubs of gpuDeviceCount and
            % gpuDevice.
            im = testCase.blurredBeads();
            for found = [1, 0]
                folder = testCase.stubGPUs(found);
                clear decon_lucy_function decon_lucy_omw_function
                evalc('testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, 1, ''Background'', testCase.Background, ''useGPU'', [])');
                ME = [];
                try
                    testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, 3, 'Background', testCase.Background, ...
                        'useGPU', true);
                catch ME
                end
                DeconFunctionTest.removeStubs(folder);
                testCase.assertNotEmpty(ME, sprintf('No error with %d GPUs found and none usable.', found));
                testCase.verifyEqual(ME.identifier, 'decon:noUsableGPU');
                testCase.verifySubstring(ME.message, sprintf('No GPU is available (%d found, CUDA_VISIBLE_DEVICES=', found));
                if found
                    testCase.verifySubstring(ME.message, ['Selecting GPU 1 failed: ', ...
                        'The GPU has compute capability 10.0 (test stub), not supported by this MATLAB.']);
                else
                    testCase.verifyFalse(contains(ME.message, 'Selecting GPU'), ME.message);
                end
            end
        end

        function defaultFallsBackToTheCPUWithAWarning(testCase, method)
            % without useGPU (or with []), a GPU if one is usable, else the
            % CPU with a warning with the same reason (stubs as above)
            im = testCase.blurredBeads();
            args = {'Background', testCase.Background, 'save16bit', false};
            outCPU = testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, 3, args{:});
            for found = [1, 0]
                folder = testCase.stubGPUs(found);
                for gpuArg = {{}, {'useGPU', []}}
                    clear decon_lucy_function decon_lucy_omw_function
                    lastwarn('');
                    switch method
                        case 'simplified'
                            out = decon_lucy_function(im, testCase.PSF, 3, args{:}, gpuArg{1}{:});
                        case 'omw'
                            out = decon_lucy_omw_function(im, testCase.PSF, testCase.BackProjector, 3, args{:}, gpuArg{1}{:});
                    end
                    msg = lastwarn;
                    testCase.verifyEqual(out, outCPU);
                    testCase.verifySubstring(msg, sprintf('No GPU is available (%d found', found));
                    testCase.verifySubstring(msg, 'so the deconvolution runs on the CPU.');
                    testCase.verifyEqual(contains(msg, 'Selecting GPU 1 failed: The GPU has compute capability 10.0'), found == 1, msg);
                end
                DeconFunctionTest.removeStubs(folder);
            end
            % the defaults of RLdecon and the in-memory large-file decon too
            if ~strcmp(method, 'simplified')
                return;
            end
            root = testCase.createTempFolder();
            psfFn = [root, 'psf.tif'];
            writetiff(testCase.PSF, psfFn);
            inFn = [root, 'in.tif'];
            writetiff(im, inFn);
            folder = testCase.stubGPUs(0);
            clear decon_lucy_function decon_lucy_omw_function
            lastwarn('');
            [~, out] = evalc(['RLdecon('''', '''', psfFn, 0.108, 0.2, 0.2, ''rawdata'', im, ''RLMethod'', ''simplified'', ', ...
                '''psfGen'', false, ''Background'', 100, ''DeconIter'', 1, ''save16bit'', false)']);
            testCase.verifySubstring(lastwarn, 'so the deconvolution runs on the CPU.', 'RLdecon');
            testCase.verifyNotEmpty(out);
            clear decon_lucy_function decon_lucy_omw_function
            lastwarn('');
            evalc(['RLdecon_large_in_memory(inFn, psfFn, [root, ''out.tif''], 0.108, 0.2, ''Background'', 100, ', ...
                '''dzPSF'', 0.2, ''DeconIter'', 1, ''RLMethod'', ''simplified'', ''psfGen'', false, ', ...
                '''batchSize'', [64, 64, 64], ''parseCluster'', false)']);
            testCase.verifySubstring(lastwarn, 'so the deconvolution runs on the CPU.', 'RLdecon_large_in_memory');
            testCase.verifyTrue(exist([root, 'out.tif'], 'file') == 2);
            DeconFunctionTest.removeStubs(folder);
        end

        function gpuBlockTaskWithoutUsableGPUErrors(testCase)
            % a large-file GPU job's block task (useGPU = GPUJob) fails, so
            % that the framework retries it and the call reports it, instead
            % of running on the CPU; with useGPU false it runs
            root = testCase.createTempFolder();
            im = testCase.blurredBeads();
            inFn = [root, 'in.zarr'];
            outFn = [root, 'out.zarr'];
            writezarr(im, inFn, 'blockSize', [16, 16, 16]);
            createzarr(outFn, dataSize=size(im), blockSize=[16, 16, 16], dtype='uint16');
            psfFn = [root, 'psf.tif'];
            writetiff(testCase.PSF, psfFn);
            mkdir([root, 'flags']);
            flagFn = [root, 'flags/blocks_1_1.mat'];
            bbox = [1, 1, 1, size(im)];
            args = {1, inFn, psfFn, outFn, flagFn, bbox, bbox, 0.108, 0.2, 'Background', testCase.Background, ...
                'DeconIter', 2, 'psfGen', false, 'RLMethod', 'simplified'};
            folder = testCase.stubGPUs(0);
            clear decon_lucy_function
            testCase.verifyError(@() RLdecon_for_zarr_block(args{:}, 'useGPU', true), 'decon:noUsableGPU');
            testCase.verifyFalse(exist(flagFn, 'file') == 2);
            evalc('RLdecon_for_zarr_block(args{:}, ''useGPU'', false)');
            DeconFunctionTest.removeStubs(folder);
            testCase.verifyTrue(exist(flagFn, 'file') == 2);
        end

        function noUsableGPUMessagesGiveMATLABsReason(testCase, method)
            % the GPUs are counted with usable_gpu_num (GPUs newer than this
            % MATLAB supports get forward compatibility), and MATLAB's reason
            % is in the warning of the default (auto: the CPU) and in the
            % error of useGPU true (a stub here: 1 GPU found, none usable)
            im = testCase.blurredBeads();
            outCPU = testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, 3, 'Background', testCase.Background, 'save16bit', false);
            testCase.applyFixture(pktest.FunctionStubFixture('usable_gpu_num', 0, 'Selecting GPU 1 failed: stub reason.', 1));
            clear decon_lucy_function decon_lucy_omw_function
            lastwarn('');
            outAuto = testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, 3, 'Background', testCase.Background, ...
                'save16bit', false, 'useGPU', []);
            testCase.verifyEqual(outAuto, outCPU);
            testCase.verifySubstring(lastwarn, 'No GPU is available (1 found, CUDA_VISIBLE_DEVICES=');
            testCase.verifySubstring(lastwarn, 'so the deconvolution runs on the CPU. Selecting GPU 1 failed: stub reason.');
            clear decon_lucy_function decon_lucy_omw_function
            ME = [];
            try
                testCase.runDecon(method, im, testCase.PSF, testCase.BackProjector, 3, 'Background', testCase.Background, ...
                    'save16bit', false, 'useGPU', true);
            catch ME
            end
            testCase.assertNotEmpty(ME, 'useGPU true ran without a usable GPU.');
            testCase.verifyEqual(ME.identifier, 'decon:noUsableGPU');
            testCase.verifySubstring(ME.message, 'Selecting GPU 1 failed: stub reason.');
        end

        function wienerButterworthBackProjectorWorksWithOMWDecon(testCase)
            % the WB back projector (Guo et al. 2020) used in the decon demo
            [~, psfWB] = evalc('BackProjector(testCase.PSF, ''wiener-butterworth'', 0.0001, 0.001, 10, 1, [], 0)');
            testCase.verifyEqual(size(psfWB), size(testCase.PSF));
            pos = [24, 24, 20];
            im = testCase.blurredBeads(pos, 1000, testCase.Background);
            out = decon_lucy_omw_function(im, testCase.PSF, psfWB, 2, 'Background', testCase.Background, ...
                'useGPU', false, 'save16bit', false);
            in = max(single(im) - testCase.Background, 0);
            testCase.verifyGreaterThan(max(out(:)), 1.5 * max(in(:)));
            testCase.verifyEqual(pktest.centroid(out, pos, 6), pktest.centroid(in, pos, 6), 'AbsTol', 0.05);
        end

        function backProjectorDCGainIsTheWienerFactor(testCase)
            % the apodization window is 1 at DC and the normalized OTF is 1
            % there, so the back projector sums to 1 / (1 + alpha)
            for alpha = [0.005, 0.05]
                [~, b] = evalc('omw_backprojector_generation(testCase.PSF, alpha, false)');
                testCase.verifyEqual(sum(b, 'all', 'double'), 1 / (1 + alpha), 'RelTol', 1e-4);
            end
        end

        function backProjectorOutputsAreConsistent(testCase)
            [~, b, OTF_bp, absOTF, mask] = evalc('omw_backprojector_generation(testCase.PSF, testCase.WienerAlpha, false)');
            sz = size(testCase.PSF);
            testCase.verifyEqual(size(b), sz);
            testCase.verifyEqual(size(OTF_bp), sz);
            c = floor(sz / 2) + 1;
            % fftshifted OTF magnitude peaks at 1 (DC) in the center
            testCase.verifyEqual(double(max(absOTF(:))), 1, 'AbsTol', 1e-6);
            testCase.verifyEqual(double(absOTF(c(1), c(2), c(3))), 1, 'AbsTol', 1e-6);
            testCase.verifyTrue(islogical(mask) && mask(c(1), c(2), c(3)), 'The OTF mask must contain DC.');
            testCase.verifyLessThan(nnz(mask), numel(mask), 'The OTF mask must not cover everything.');
            % the back projector is the inverse FFT of OTF_bp
            testCase.verifyImagesClose(b, fftshift(real(ifftn(OTF_bp))), 'relTol', 1e-5);
        end

        function higherOTFThresholdGivesLargerMask(testCase)
            [~, ~, ~, ~, mask85] = evalc('omw_backprojector_generation(testCase.PSF, testCase.WienerAlpha, false, ''OTFCumThresh'', 0.85)');
            [~, ~, ~, ~, mask95] = evalc('omw_backprojector_generation(testCase.PSF, testCase.WienerAlpha, false, ''OTFCumThresh'', 0.95)');
            testCase.verifyGreaterThan(nnz(mask95), nnz(mask85));
        end

        function backProjectorRejectsInvalidPSF(testCase)
            testCase.verifyError(@() omw_backprojector_generation(zeros(31, 31, 31, 'single'), 0.005, false), ?MException);
            % a delta PSF has a flat OTF, so no OTF support can be segmented
            delta = zeros(31, 31, 31, 'single');
            delta(16, 16, 16) = 1;
            testCase.verifyError(@() evalc('omw_backprojector_generation(delta, 0.005, false)'), ?MException);
        end

        function psfGenCentersOddSizePSF(testCase)
            % without background (zero edge slices) the PSF is only
            % circularly shifted so its peak is at the center
            psf = pktest.beads([31, 31, 31], [13, 18, 14], testCase.PSFSigma, 1e4, 0);
            psf(:, :, [1 : 5, end - 4 : end]) = 0;
            out = psf_gen_new(psf, 0.2, 0.2);
            testCase.verifyEqual(out, circshift(psf, [16, 16, 16] - [13, 18, 14]));
        end

        function psfGenRemovesBackground(testCase)
            % 'masked' method: keep the component around the peak and
            % subtract the background estimated from the border rows
            bg = 100;
            psf = pktest.beads([41, 41, 41], [21, 21, 21], testCase.PSFSigma, 5000, bg) + 3 * randn(41, 41, 41, 'single');
            out = psf_gen_new(psf, 0.2, 0.2);
            ideal = pktest.beads([41, 41, 41], [21, 21, 21], testCase.PSFSigma, 5000, 0);

            testCase.verifyGreaterThanOrEqual(min(out(:)), single(0));
            testCase.verifyEqual(out(1 : 5, :, :), zeros(5, 41, 41, 'single'), 'Background far from the peak remains.');
            testCase.verifyEqual(double(max(out(:))), 5000, 'RelTol', 0.03);
            testCase.verifyEqual(sum(out, 'all', 'double'), sum(ideal, 'all', 'double'), 'RelTol', 0.05);
        end

        function psfGen2DMedianMethodRemovesIsolatedPoints(testCase)
            % 2D 'median' method (it called medfilt, which doesn't exist):
            % subtract medFactor times the median of the positive pixels and
            % keep the component around the peak, without an isolated hot
            % pixel or a line of 1 pixel width from the PSF's edge (which
            % only the median filter cuts off); the peak is moved to the
            % center
            psf = pktest.beads([41, 41, 1], [18, 23, 1], [2, 2], 5000, 100);
            psf(5, 35) = 3000;
            psf(18, 30 : 41) = 3000;
            out = psf_gen_new(psf, 0.2, 0.2, 1.5, 'median');

            core = max(psf - 1.5 * median(psf(psf > 0)), 0);
            expected = circshift(core, [21, 21] - [18, 23]);
            testCase.verifyEqual(size(out), [41, 41]);
            testCase.verifyEqual(out(21, 21), max(core(:)), 'The peak is not at the center.');
            testCase.verifyTrue(all(out == expected | out == 0, 'all'), 'Values other than those of the PSF without background.');
            testCase.verifyEqual(out(8, 33), single(0), 'The isolated hot pixel remains.');
            testCase.verifyEqual(out(21, 31 : 39), zeros(1, 9, 'single'), 'The line from the PSF remains.');
            core(5, 35) = 0;
            core(18, 30 : 41) = 0;
            testCase.verifyEqual(sum(out, 'all', 'double'), sum(core, 'all', 'double'), 'RelTol', 0.01, ...
                'The component around the peak is not kept whole.');
        end

        function psfGenResamplesZToTheDataStep(testCase)
            % a PSF sampled at half the data z step is resampled to the
            % data step by Fourier truncation: half the slices, half the width
            psf = testCase.PSF;
            out = psf_gen_new(psf, 0.1, 0.2);
            testCase.verifyEqual(size(out), [31, 31, 16]);

            zIn = squeeze(sum(double(psf), [1, 2]));
            zOut = squeeze(sum(out, [1, 2]));
            sdIn = sqrt(sum(((1 : numel(zIn))' - sum((1 : numel(zIn))' .* zIn) / sum(zIn)) .^ 2 .* zIn) / sum(zIn));
            sdOut = sqrt(sum(((1 : numel(zOut))' - sum((1 : numel(zOut))' .* zOut) / sum(zOut)) .^ 2 .* zOut) / sum(zOut));
            testCase.verifyEqual(sdOut, sdIn / 2, 'RelTol', 0.1);
            % no resampling if the data step is smaller than the PSF step
            testCase.verifyEqual(psf_gen_new(psf, 0.2, 0.1), psf);
        end

        function psf2otfOfCenteredDeltaIsFlat(testCase)
            sz = [8, 9, 7];
            delta = zeros(sz, 'single');
            c = floor(sz / 2) + 1;
            delta(c(1), c(2), c(3)) = 1;
            otf = decon_psf2otf(delta);
            testCase.verifyEqual(double(abs(otf - 1)), zeros(sz), 'AbsTol', 1e-6);
        end

        function psf2otfPadsToTheOutputSize(testCase)
            psf = rand(5, 5, 5, 'single');
            otf = decon_psf2otf(psf, [16, 16, 16]);
            testCase.verifyEqual(size(otf), [16, 16, 16]);
            testCase.verifyEqual(double(real(otf(1, 1, 1))), sum(psf, 'all', 'double'), 'RelTol', 1e-6);
            % the PSF center moves to (1, 1, 1): undo the shift and compare
            back = circshift(real(ifftn(otf)), [2, 2, 2]);
            testCase.verifyImagesClose(back(1 : 5, 1 : 5, 1 : 5), psf, 'relTol', 1e-5);
            testCase.verifyError(@() decon_psf2otf(rand(9, 9, 9), [8, 8, 8]), ?MException);
        end

        function rldeconInMemoryMatchesDirectFunction(testCase, inputClass)
            % the background is subtracted once, also for float input
            % (regression test for dda6375, which subtracted it twice)
            dataPath = testCase.createTempFolder();
            psfFn = [dataPath, 'psf.tif'];
            writetiff(testCase.PSF, psfFn);
            im = cast(testCase.blurredBeads(), inputClass);
            [~, out] = evalc(['RLdecon('''', '''', psfFn, 0.108, 0.2, 0.2, ''rawdata'', im, ''RLMethod'', ''simplified'', ', ...
                '''psfGen'', false, ''useGPU'', false, ''Background'', 100, ''DeconIter'', 5, ''save16bit'', false)']);
            clear decon_lucy_function
            % RLdecon normalizes the PSF (summed in double) before passing it on
            psf = testCase.PSF ./ sum(testCase.PSF(:), 'double');
            expected = decon_lucy_function(im, psf, 5, 'Background', 100, 'useGPU', false, 'save16bit', false);
            testCase.verifyEqual(out, expected);
        end

        function rldeconWithZeroIterationsReturnsTheInput(testCase)
            dataPath = testCase.createTempFolder();
            psfFn = [dataPath, 'psf.tif'];
            writetiff(testCase.PSF, psfFn);
            im = testCase.blurredBeads();
            [~, out] = evalc(['RLdecon('''', '''', psfFn, 0.108, 0.2, 0.2, ''rawdata'', im, ''RLMethod'', ''simplified'', ', ...
                '''psfGen'', false, ''useGPU'', false, ''Background'', 100, ''DeconIter'', 0)']);
            testCase.verifyEqual(out, im);
        end

        function psfGenKeepsEvenSizePSFAtTheFFTCenter(testCase)
            % psf_gen_new used to center even-size PSFs one voxel off the FFT
            % center, depending on where the peak started
            psf = pktest.deconPSF([32, 32, 32], testCase.PSFSigma);
            out = psf_gen_new(psf, 0.2, 0.2);
            [~, ind] = max(out(:));
            [y, x, z] = ind2sub(size(out), ind);
            testCase.verifyEqual([y, x, z], [17, 17, 17], 'The PSF peak moved away from the FFT center.');
        end
    end

    methods (Test, TestTags = {'Unit'})
        function loopMatchesThePreviousFormulation(testCase, method, anyInputClass, deconOptions)
            % the decon loops run one FFT step per statement (and keep uint16
            % input as uint16 on the GPU) to save GPU memory; the results are
            % bit-identical to the previous loops (tests/decon/private/*_before.m),
            % here on the CPU and in DeconGPUTest on the GPU
            im = cast(testCase.blurredBeads(), anyInputClass);
            args = [{'Background', testCase.Background, 'useGPU', false}, deconOptions];
            clear decon_lucy_function_before decon_lucy_omw_function_before
            for call = 1 : 2
                % the second call uses the cached OTFs
                switch method
                    case 'simplified'
                        expected = decon_lucy_function_before(im, testCase.PSF, 4, args{:});
                        actual = decon_lucy_function(im, testCase.PSF, 4, args{:});
                    case 'omw'
                        expected = decon_lucy_omw_function_before(im, testCase.PSF, testCase.BackProjector, 2, args{:});
                        actual = decon_lucy_omw_function(im, testCase.PSF, testCase.BackProjector, 2, args{:});
                end
                testCase.verifyEqual(actual, expected, sprintf('Call %d.', call));
            end
        end
    end

    methods (Test, TestTags = {'Unit'})
        function changedPSFWithSameEnergyIsNotIgnored(testCase, method)
            % the decon functions cache the OTF in persistent variables; a
            % shifted or mirrored PSF (e.g., a second channel) gets its own
            % OTF, also when its sum of squares (the cache's key before) is
            % the same. For OMW the back projector is the same, so this also
            % checks that a new forward PSF isn't ignored.
            im = testCase.blurredBeads();
            % integer values make sum(psf .^ 2) exact regardless of the
            % summation order, so that the two old keys are identical on
            % every platform
            psf1 = round(testCase.PSF / 1000);
            psf2 = circshift(psf1, [0, 0, 2]);
            testCase.runDecon(method, im, psf1, testCase.BackProjector, 3, 'Background', testCase.Background, 'save16bit', false);
            out = testCase.runDecon(method, im, psf2, testCase.BackProjector, 3, 'Background', testCase.Background, 'save16bit', false);
            clear decon_lucy_function decon_lucy_omw_function
            fresh = testCase.runDecon(method, im, psf2, testCase.BackProjector, 3, 'Background', testCase.Background, 'save16bit', false);
            testCase.verifyEqual(out, fresh, 'The OTF of the previous PSF was reused.');
        end

        function changedBackProjectorWithSameEnergyIsNotIgnored(testCase)
            % OMW with the same forward PSF and a shifted back projector of
            % the same sum of squares (exact for integer values) gets its own
            % back projector OTF
            im = testCase.blurredBeads();
            psf = round(testCase.PSF / 1000);
            psfB2 = circshift(psf, [0, 0, 2]);
            testCase.runDecon('omw', im, psf, psf, 3, 'Background', testCase.Background, 'save16bit', false);
            out = testCase.runDecon('omw', im, psf, psfB2, 3, 'Background', testCase.Background, 'save16bit', false);
            clear decon_lucy_omw_function
            fresh = testCase.runDecon('omw', im, psf, psfB2, 3, 'Background', testCase.Background, 'save16bit', false);
            testCase.verifyEqual(out, fresh, 'The OTF of the previous back projector was reused.');
        end

        function otfIsCachedForTheSamePSFBitsAndThreads(testCase, method)
            % the cached OTF is used for a PSF with the same values (another
            % array), but computed again for one that only differs in the
            % sign of a zero (the key is bit for bit) or with another number
            % of threads (the PSF normalization depends on it); counted
            % with a copy of decon_psf2otf that logs its calls (MATLAB code)
            countFn = testCase.countOTFs();
            nOTF = 1 + strcmp(method, 'omw');
            im = testCase.blurredBeads();
            % a zero whose sign changes below
            psf = testCase.PSF;
            psf(1) = 0;
            args = {3, 'Background', testCase.Background, 'save16bit', false, 'useMex', false};
            expected = testCase.runDecon(method, im, psf, testCase.BackProjector, args{:});
            out = testCase.runDecon(method, im, psf + 0, testCase.BackProjector, args{:});
            testCase.verifyEqual(out, expected);
            testCase.verifyEqual(countFn(), nOTF, 'The OTF of the same PSF was computed again.');

            psfNegZero = psf;
            psfNegZero(1) = -0;
            testCase.runDecon(method, im, psfNegZero, testCase.BackProjector, args{:});
            testCase.verifyEqual(countFn(), 2 * nOTF, 'The OTF of a PSF with another sign of a zero was not computed.');

            nt = maxNumCompThreads;
            testCase.addTeardown(@maxNumCompThreads, nt);
            maxNumCompThreads(1 + (nt == 1));
            testCase.runDecon(method, im, psfNegZero, testCase.BackProjector, args{:});
            testCase.verifyEqual(countFn(), 3 * nOTF, 'The OTF was not computed again with another number of threads.');
        end
    end

    methods (Access = private)
        function [folder] = stubGPUs(testCase, found)
            % gpuDeviceCount and gpuDevice first on the path: found GPUs
            % ("all"), none available, and selecting one fails as for an
            % unsupported GPU (pktest.FunctionStubFixture returns the same
            % values for every call); a teardown takes them off the path
            folder = testCase.createTempFolder();
            fid = fopen([folder, 'gpuDeviceCount.m'], 'w');
            fprintf(fid, '%s\n', 'function [n] = gpuDeviceCount(varargin)', ...
                '% test stub written by DeconFunctionTest: GPUs found, none available', 'n = 0;', ...
                'if nargin > 0 && strcmp(varargin{1}, ''all'')', sprintf('    n = %d;', found), 'end', 'end');
            fclose(fid);
            fid = fopen([folder, 'gpuDevice.m'], 'w');
            fprintf(fid, '%s\n', 'function [varargout] = gpuDevice(varargin)', ...
                '% test stub written by DeconFunctionTest: selecting a GPU fails as for an unsupported one', ...
                'error(''pktest:gpuDevice'', ''The GPU has compute capability 10.0 (test stub),\nnot supported by this MATLAB.'');', 'end');
            fclose(fid);
            warning('off', 'MATLAB:dispatcher:nameConflict');
            addpath(folder, '-begin');
            testCase.addTeardown(@() DeconFunctionTest.removeStubs(folder));
            clear gpuDeviceCount gpuDevice
        end

        function [im] = blurredBeads(testCase, positions, amplitude, background)
            % beads blurred by the test PSF (Gaussian object of sigma
            % ObjectSigma convolved with the Gaussian PSF), as uint16 camera
            % data. Noise-free single data would have tails far below the FFT
            % round-off, where the eps-clamped RL ratio is unstable.
            if nargin < 2 || isempty(positions)
                positions = [12, 14, 10; 24.3, 23.6, 20.2; 36, 30, 28; 15, 35, 25; 33, 12, 15];
            end
            if nargin < 3 || isempty(amplitude)
                amplitude = 1000;
            end
            if nargin < 4 || isempty(background)
                background = testCase.Background;
            end
            sigma = sqrt(testCase.ObjectSigma ^ 2 + testCase.PSFSigma .^ 2);
            im = uint16(pktest.beads(testCase.ImSize, positions, sigma, amplitude, background));
        end

        function [countFn] = countOTFs(testCase)
            % a copy of decon_psf2otf, first on the path while the test
            % runs, that logs its calls; countFn() returns their number
            folder = testCase.createTempFolder();
            logFn = [folder, 'calls.txt'];
            lines = splitlines(fileread(which('decon_psf2otf')));
            lines = [lines(1); {sprintf('fid = fopen(''%s'', ''a''); fprintf(fid, ''1\\n''); fclose(fid);', logFn)}; lines(2 : end)];
            fid = fopen([folder, 'decon_psf2otf.m'], 'w');
            fprintf(fid, '%s\n', lines{:});
            fclose(fid);
            fclose(fopen(logFn, 'w'));
            addpath(folder, '-begin');
            testCase.addTeardown(@() clear('decon_psf2otf'));
            testCase.addTeardown(@rmpath, folder);
            clear decon_psf2otf decon_lucy_function decon_lucy_omw_function
            countFn = @() numel(strfind(fileread(logFn), '1'));
        end

        function [out] = runDecon(~, method, im, psf, psfB, nIter, varargin)
            switch method
                case 'simplified'
                    out = decon_lucy_function(im, psf, nIter, 'useGPU', false, varargin{:});
                case 'omw'
                    out = decon_lucy_omw_function(im, psf, psfB, nIter, 'useGPU', false, varargin{:});
            end
        end

        function [c, v] = moments(~, vol, pos)
            % centroid and per axis variance within 8 voxels of pos
            r = 8;
            lo = max(round(pos) - r, 1);
            hi = min(round(pos) + r, size(vol, 1 : 3));
            win = double(vol(lo(1) : hi(1), lo(2) : hi(2), lo(3) : hi(3)));
            [y, x, z] = ndgrid(lo(1) : hi(1), lo(2) : hi(2), lo(3) : hi(3));
            s = sum(win(:));
            c = [sum(y(:) .* win(:)), sum(x(:) .* win(:)), sum(z(:) .* win(:))] / s;
            v = [sum((y(:) - c(1)) .^ 2 .* win(:)), sum((x(:) - c(2)) .^ 2 .* win(:)), sum((z(:) - c(3)) .^ 2 .* win(:))] / s;
        end
    end

    methods (Static, Access = private)
        function [] = removeStubs(folder)
            % take the stub folder off the path (if still on it; path entries
            % have backslashes on Windows) and the stubs out of memory
            entries = strsplit(path, pathsep);
            ind = find(strcmp(strrep(entries, '\', '/'), regexprep(strrep(folder, '\', '/'), '/$', '')), 1);
            if ~isempty(ind)
                rmpath(entries{ind});
            end
            clear gpuDeviceCount gpuDevice
        end
    end
end
