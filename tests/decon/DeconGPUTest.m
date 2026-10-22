classdef DeconGPUTest < pktest.TestCase
    % GPU deconvolution must match the CPU result (skipped without a GPU)

    properties (Constant)
        ImSize = [48, 48, 40]
        Background = 100
        PSFSigma = [1.5, 2]
        BlurSigma = sqrt(0.8 ^ 2 + [1.5, 2] .^ 2)
    end

    properties (TestParameter)
        method = {'simplified', 'omw'}
        inputClass = {'uint16', 'single'}
        anyInputClass = {'uint16', 'single', 'double'}
        % an FFT-friendly size and one with prime factors (other cuFFT plans)
        gpuSize = struct('smooth', [48, 48, 40], 'oddPrimes', [61, 53, 37])
    end

    methods (TestMethodSetup)
        function requireGPU(testCase)
            testCase.assumeGPUAvailable();
            testCase.addTeardown(@() reset(gpuDevice));
        end
    end

    methods (Test, TestTags = {'Unit', 'GPU'})
        function gpuMatchesCPU(testCase, method, inputClass)
            psf = pktest.deconPSF([31, 31, 31], testCase.PSFSigma);
            pos = [12, 14, 10; 24.3, 23.6, 20.2; 36, 30, 28; 15, 35, 25; 33, 12, 15];
            im = cast(uint16(pktest.beads(testCase.ImSize, pos, testCase.BlurSigma, 1000, testCase.Background)), inputClass);
            args = {'Background', testCase.Background, 'save16bit', false};
            switch method
                case 'simplified'
                    cpu = decon_lucy_function(im, psf, 5, 'useGPU', false, args{:});
                    gpu = decon_lucy_function(im, psf, 5, 'useGPU', true, args{:});
                case 'omw'
                    [~, b] = evalc('omw_backprojector_generation(psf, 0.005, false)');
                    cpu = decon_lucy_omw_function(im, psf, b, 2, 'useGPU', false, args{:});
                    gpu = decon_lucy_omw_function(im, psf, b, 2, 'useGPU', true, args{:});
            end
            testCase.verifyFalse(isa(gpu, 'gpuArray'), 'The result must be gathered to the CPU.');
            testCase.verifyClass(gpu, class(cpu));
            % different FFT libraries, single precision
            testCase.verifyImagesClose(gpu, cpu, 'relTol', 1e-4);
        end
    end

    methods (Test, TestTags = {'Unit', 'GPU'})
        function gpuLoopMatchesThePreviousFormulation(testCase, method, anyInputClass, gpuSize)
            % the decon loops run one FFT step per statement and keep uint16
            % input as uint16 on the GPU, which saves 6 bytes per voxel of GPU
            % memory; the results are bit-identical to the previous loops
            % (tests/decon/private/*_before.m)
            psf = pktest.deconPSF([31, 31, 31], testCase.PSFSigma);
            pos = [12, 14, 10; 24.3, 23.6, 20.2; 36, 30, 28; 15, 35, 25; 33, 12, 15];
            im = cast(uint16(pktest.beads(gpuSize, pos, testCase.BlurSigma, 1000, testCase.Background)), anyInputClass);
            [~, b] = evalc('omw_backprojector_generation(psf, 0.005, false)');
            optionSets = {{'save16bit', true}, {'save16bit', false}, ...
                {'save16bit', false, 'dampFactor', 1.2, 'EdgeErosion', 4, 'deconOffset', 3, 'scaleFactor', 1.5}, ...
                {'save16bit', true, 'bbox', [3, 4, 2, gpuSize - [4, 3, 2]]}};
            for o = 1 : numel(optionSets)
                args = [{'Background', testCase.Background, 'useGPU', true}, optionSets{o}];
                clear decon_lucy_function decon_lucy_omw_function decon_lucy_function_before decon_lucy_omw_function_before
                for call = 1 : 2
                    % the second call uses the cached OTFs
                    switch method
                        case 'simplified'
                            expected = decon_lucy_function_before(im, psf, 4, args{:});
                            actual = decon_lucy_function(im, psf, 4, args{:});
                        case 'omw'
                            expected = decon_lucy_omw_function_before(im, psf, b, 2, args{:});
                            actual = decon_lucy_omw_function(im, psf, b, 2, args{:});
                    end
                    testCase.verifyEqual(actual, expected, sprintf('Options %d, call %d.', o, call));
                end
            end
        end
    end

    methods (Test, TestTags = {'Integration', 'GPU'})
        function wrapperGPUJobMatchesCPU(testCase)
            root = testCase.createTempFolder();
            dataPath = [root, 'data/'];
            mkdir(dataPath);
            mkdir([root, 'psf']);
            psfFn = [root, 'psf/PSF_CamA.tif'];
            writetiff(uint16(pktest.deconPSF([31, 31, 31], testCase.PSFSigma)), psfFn);
            fsn = 'Scan_Iter_0000_0000_CamA_ch0_CAM1_stack0000_488nm_0000000msec_0000000000msecAbs_000x_000y_000z_0000t';
            pos = [8 + 32 * rand(6, 2), 6 + 28 * rand(6, 1)];
            writetiff(uint16(pktest.beads(testCase.ImSize, pos, testCase.BlurSigma, 1000, testCase.Background)), [dataPath, fsn, '.tif']);

            for gpuJob = [false, true]
                resultDirName = sprintf('decon_gpu_%d', gpuJob);
                evalc(['XR_decon_data_wrapper({dataPath}, ''resultDirName'', resultDirName, ''channelPatterns'', {''CamA_ch0''}, ', ...
                    '''psfFullpaths'', {psfFn}, ''xyPixelSize'', 0.108, ''dz'', 0.2, ''dzPSF'', 0.2, ''background'', 100, ', ...
                    '''RLMethod'', ''omw'', ''deconIter'', 2, ''skewed'', false, ''parseCluster'', false, ''GPUJob'', gpuJob)']);
            end
            % uint16 results: allow rounding differences of 1
            testCase.verifyImagesClose(readtiff([dataPath, 'decon_gpu_1/', fsn, '.tif']), ...
                readtiff([dataPath, 'decon_gpu_0/', fsn, '.tif']), 'relTol', 1e-3, 'absTol', 1);
        end
    end
end
