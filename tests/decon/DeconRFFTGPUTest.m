classdef DeconRFFTGPUTest < pktest.TestCase
    % the GPU half-spectrum FFT decon (decon_lucy_rfft_gpu_mex) against the
    % MATLAB GPU loops of decon_lucy_function and decon_lucy_omw_function
    % (skipped without a GPU or the mex). Unlike the CPU mex, it uses other
    % cuFFT transforms than MATLAB's gpuArray FFTs, so the paths differ by
    % single-precision rounding: up to 3e-5 relative L2 error on these volumes
    % (OMW, whose back projector amplifies high frequencies; RL about 2e-6).

    properties (Constant)
        ImSize = [48, 46, 36]
        Background = 100
        BlurSigma = sqrt(0.8 ^ 2 + [1.5, 2] .^ 2)
        RelTol = 3e-5
    end

    properties (TestParameter)
        method = {'simplified', 'omw'}
        inputClass = {'uint16', 'single'}
        dataSize = {[47, 45, 33], [64, 50], [1, 64, 33]}
        options = {{'dampFactor', 1.2}, {'EdgeErosion', 3}, {'bbox', [3, 4, 5, 40, 41, 30]}, ...
            {'EdgeErosion', 2, 'scaleFactor', 1.5, 'deconOffset', 2, 'bbox', [2, 2, 2, 44, 40, 33], 'save16bit', true}}
    end

    properties (SetAccess = private)
        PSF
        BackProjector
    end

    methods (TestClassSetup)
        function createPSFAndBackProjector(testCase)
            testCase.PSF = pktest.deconPSF([31, 31, 31], [1.5, 2]);
            [~, testCase.BackProjector] = evalc('omw_backprojector_generation(testCase.PSF, 0.005, false)');
        end
    end

    methods (TestMethodSetup)
        function requireGPUAndMex(testCase)
            testCase.assumeGPUAvailable();
            testCase.assumeMexAvailable('decon_lucy_rfft_gpu_mex');
            testCase.addTeardown(@() reset(gpuDevice));
        end
    end

    methods (Test, TestTags = {'Unit', 'GPU'})
        function gpuMexMatchesGpuMatlabLoop(testCase, method, inputClass)
            im = cast(testCase.blurredBeads(testCase.ImSize), inputClass);
            args = {'Background', testCase.Background, 'save16bit', false};
            ref = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', false);
            out = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', true);
            testCase.verifyFalse(isa(out, 'gpuArray'), 'The result must be gathered to the CPU.');
            testCase.verifyClass(out, 'single');
            testCase.verifyImagesClose(out, ref, 'relTol', testCase.RelTol);
            % the default is the MATLAB loop (see decon_use_rfft_mex)
            testCase.verifyEqual(testCase.runDecon(method, im, testCase.PSF, args{:}), ref);
            [~, canUseMex] = decon_use_rfft_mex([], method, im, true, false, testCase.Background);
            testCase.verifyTrue(canUseMex);
        end

        function gpuMexMatchesForOddAnd2DSizes(testCase, method, dataSize)
            sz3 = [dataSize, ones(1, 3 - numel(dataSize))];
            im = testCase.blurredBeads(sz3);
            psf = testCase.cropToSize(testCase.PSF, sz3);
            args = {'Background', testCase.Background, 'save16bit', false};
            ref = testCase.runDecon(method, im, psf, args{:}, 'useMex', false);
            out = testCase.runDecon(method, im, psf, args{:}, 'useMex', true);
            testCase.verifyEqual(size(out), size(ref));
            testCase.verifyImagesClose(out, ref, 'relTol', testCase.RelTol);
        end

        function gpuMexAppliesTheOptions(testCase, method, options)
            im = testCase.blurredBeads(testCase.ImSize);
            args = [{'Background', testCase.Background, 'save16bit', false}, options];
            ref = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', false);
            out = testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', true);
            testCase.verifyClass(out, class(ref));
            testCase.verifyEqual(size(out), size(ref));
            if isa(ref, 'uint16')
                d = abs(double(out) - double(ref));
                testCase.verifyLessThanOrEqual(max(d(:)), 1);
                testCase.verifyLessThan(mean(d(:) > 0), 1e-3);
            else
                testCase.verifyImagesClose(out, ref, 'relTol', testCase.RelTol);
            end
        end

        function gpuBackgroundIsSubtractedAsInTheGpuLoop(testCase, method)
            % uint16 data: saturating subtraction of uint16(background), as the
            % MATLAB GPU code does (unlike the CPU code for non-integer values)
            im = testCase.blurredBeads(testCase.ImSize);
            for bg = [100.5, 99.4]
                out = testCase.runDecon(method, im, testCase.PSF, 'Background', bg, 'save16bit', false, ...
                    'useMex', true, 'NUMIT', 0);
                testCase.verifyEqual(out, single(im - uint16(bg)));
            end
            ims = single(im) + 0.37;
            out = testCase.runDecon(method, ims, testCase.PSF, 'Background', 100.3, 'save16bit', false, ...
                'useMex', true, 'NUMIT', 0);
            testCase.verifyEqual(out, max(ims - 100.3, 0));
        end

        function gpuHalfOTFIsTheFirstHalfOfTheOTF(testCase)
            % the rows up to the middle of the first dimension
            psf = testCase.PSF ./ sum(testCase.PSF(:));
            for sz = {[48, 46, 36], [47, 45, 33]}
                H = decon_lucy_rfft_gpu_mex('otf', gpuArray(psf), sz{1});
                testCase.verifyTrue(isa(H, 'gpuArray'));
                otf = decon_psf2otf(psf, sz{1});
                otf = otf(1 : floor(sz{1}(1) / 2) + 1, :, :);
                testCase.verifyEqual(size(H), size(otf));
                testCase.verifyEqual(double(gather(H)), double(otf), 'AbsTol', 1e-6);
            end
        end

        function outOfMemoryLeavesNoGpuArraysBehind(testCase)
            % an allocation that fails inside the mex (here its third
            % buffer) frees the buffers it made before, so a CPU fallback or
            % the next file has the memory
            g = gpuDevice;
            sz = [256, 256, 128];
            im = gpuArray(uint16(100 + randi(1000, sz)));
            H = decon_lucy_rfft_gpu_mex('otf', gpuArray(testCase.PSF ./ sum(testCase.PSF(:))), sz);
            wait(g);
            % fill the GPU, in arrays of at most 1 GB (gpuArrays have at most
            % 2^31 - 1 elements), up to room for about 2.5 of the mex's 4
            % bytes/voxel buffers
            leave = 2.5 * 4 * prod(sz);
            ballast = {};
            while g.AvailableMemory - leave > 2 ^ 30
                ballast{end + 1} = gpuArray.zeros(2 ^ 28, 1, 'single');
            end
            ballast{end + 1} = gpuArray.zeros(round((g.AvailableMemory - leave) / 4), 1, 'single');
            wait(g);
            before = g.AvailableMemory;
            testCase.verifyError(@() decon_lucy_rfft_gpu_mex('omw', im, H, H, 2, 100), 'parallel:gpu:array:OOM');
            wait(g);
            % a leak is at least one buffer (33 MB); the allocator's own
            % accounting moves by a few MB
            testCase.verifyGreaterThan(g.AvailableMemory, before - 2 * prod(sz), 'GPU memory was left allocated.');
            clear ballast
        end

        function withoutTheGpuMexTheGpuLoopRuns(testCase, method)
            testCase.applyFixture(pktest.FunctionStubFixture('decon_lucy_rfft_gpu_mex'));
            im = testCase.blurredBeads(testCase.ImSize);
            args = {'Background', testCase.Background, 'save16bit', false};
            out = testCase.runDecon(method, im, testCase.PSF, args{:});
            clear decon_lucy_function decon_lucy_omw_function
            testCase.verifyEqual(out, testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', false));
            [~, canUseMex] = decon_use_rfft_mex([], method, im, true, false, testCase.Background);
            testCase.verifyFalse(canUseMex);
            testCase.verifyError(@() testCase.runDecon(method, im, testCase.PSF, args{:}, 'useMex', true), ?MException);
        end
    end

    methods (Access = private)
        function [im] = blurredBeads(testCase, sz)
            n = max(3, round(prod(sz) / 4000));
            pos = (sz - 1) .* rand(n, 3) + 1;
            im = pktest.beads(sz, pos, testCase.BlurSigma, 1000, testCase.Background);
            im = uint16(im + 3 * rand(sz, 'single'));
        end

        function [out] = runDecon(testCase, method, im, psf, varargin)
            % on the GPU; NUMIT can be given as a name-value pair (default 10 for RL, 2 for OMW)
            nIter = 10;
            if strcmp(method, 'omw')
                nIter = 2;
            end
            ind = find(strcmp(varargin(1 : 2 : end), 'NUMIT'));
            if ~isempty(ind)
                nIter = varargin{2 * ind};
                varargin(2 * ind - 1 : 2 * ind) = [];
            end
            switch method
                case 'simplified'
                    out = decon_lucy_function(im, psf, nIter, 'useGPU', true, varargin{:});
                case 'omw'
                    b = testCase.cropToSize(testCase.BackProjector, size(im, 1 : 3));
                    out = decon_lucy_omw_function(im, psf, b, nIter, 'useGPU', true, varargin{:});
            end
        end
    end

    methods (Static, Access = private)
        function [psf] = cropToSize(psf, sz)
            pSize = size(psf, 1 : 3);
            s = max(0, floor((pSize - sz) / 2)) + 1;
            t = s + min(pSize, sz) - 1;
            psf = psf(s(1) : t(1), s(2) : t(2), s(3) : t(3));
        end
    end
end
