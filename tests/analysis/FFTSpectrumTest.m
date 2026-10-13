classdef FFTSpectrumTest < pktest.TestCase
    % XR_fftSpectrumComputingWrapper ("FFT analysis" in demo_useful_tools.m).
    % A cosine along x with an integer number of periods has a power
    % spectrum with exactly three nonzero peaks, so the outputs can be
    % checked analytically.

    properties (Constant)
        N = 65
        % periods across the N voxels along x
        K = 8
    end

    methods (Test, TestTags = {'Integration'})
        function spectrumOfCosineHasExpectedPeaks(testCase)
            [dataPath, fsn] = testCase.cosineImage();

            XR_fftSpectrumComputingWrapper({dataPath}, 'resultDirName', 'FFT', 'xyPixelSize', 0.1, 'dz', 0.1, ...
                'outSize', [testCase.N, testCase.N, testCase.N], 'outPixelSize', 0.1, 'save3DStack', false, ...
                'parseCluster', false, 'masterCompute', true);

            fftPath = [dataPath, 'FFT/'];
            % without save3DStack, the stack is an empty placeholder
            dir_info = dir([fftPath, fsn, '.tif']);
            testCase.verifyEqual(dir_info.bytes, 0);

            c = (testCase.N + 1) / 2;
            xy = readtiff([fftPath, 'central_slices/', fsn, '_center_xy.tif']);
            testCase.verifySize(xy, [testCase.N, testCase.N]);
            % amplitude (gamma 0.5) normalized to the DC peak: the cosine
            % has half the DC amplitude at +/-K, i.e., 0.25 in power units
            testCase.verifyEqual(double(xy(c, c)), 1, 'AbsTol', 1e-5);
            testCase.verifyEqual(double(xy(c, c + testCase.K)), sqrt(0.25), 'AbsTol', 1e-4);
            testCase.verifyEqual(double(xy(c, c - testCase.K)), sqrt(0.25), 'AbsTol', 1e-4);
            others = true(size(xy));
            others(c, [c - testCase.K, c, c + testCase.K]) = false;
            % (gamma 0.5 turns single precision round-off of ~1e-6 into ~1e-3)
            testCase.verifyLessThan(max(xy(others)), 0.01, 'Unexpected spectral content.');
            % power spectra of real images are point symmetric about DC
            testCase.verifyImagesClose(xy, rot90(xy, 2), 'relTol', 1e-5);

            for axisName = {'y', 'x', 'z'}
                mip = readtiff(sprintf('%sMIPs/%s_MIP_%s.tif', fftPath, fsn, axisName{1}));
                testCase.verifySize(mip, [testCase.N, testCase.N]);
                testCase.verifyEqual(double(max(mip(:))), 1, 'AbsTol', 1e-5);
            end
            for planeName = {'yz', 'xz'}
                testCase.verifyTrue(exist(sprintf('%scentral_slices/%s_center_%s.tif', fftPath, fsn, planeName{1}), 'file') == 2);
            end
        end

        function save3DStackWritesLogSpectrum(testCase)
            [dataPath, fsn] = testCase.cosineImage();

            XR_fftSpectrumComputingWrapper({dataPath}, 'xyPixelSize', 0.1, 'dz', 0.1, ...
                'outSize', [testCase.N, testCase.N, testCase.N], 'outPixelSize', 0.1, 'save3DStack', true, ...
                'parseCluster', false, 'masterCompute', true);

            spectrum = readtiff([dataPath, 'FFT/', fsn, '.tif']);
            c = (testCase.N + 1) / 2;
            testCase.verifySize(spectrum, [testCase.N, testCase.N, testCase.N]);
            % log10 of the normalized power: 0 at DC, log10(0.25) at +/-K
            testCase.verifyEqual(double(spectrum(c, c, c)), 0, 'AbsTol', 1e-5);
            testCase.verifyEqual(double(spectrum(c, c + testCase.K, c)), log10(0.25), 'AbsTol', 1e-4);
        end

        function resamplesAnisotropicVoxels(testCase)
            % input with dz = 2 * xyPixelSize is resampled to isotropic
            % voxels of outPixelSize before the FFT, so the output size is
            % outSize regardless of the input size
            dataPath = testCase.createTempFolder();
            fsn = 'aniso';
            writetiff(single(pktest.texture([80, 80, 40], 2)), [dataPath, fsn, '.tif']);

            XR_fftSpectrumComputingWrapper({dataPath}, 'xyPixelSize', 0.1, 'dz', 0.2, ...
                'outSize', [testCase.N, testCase.N, testCase.N], 'parseCluster', false, 'masterCompute', true);

            xz = readtiff([dataPath, 'FFT/central_slices/', fsn, '_center_xz.tif']);
            testCase.verifySize(xz, [testCase.N, testCase.N]);
        end
    end

    methods (Access = private)
        function [dataPath, fsn] = cosineImage(testCase)
            dataPath = testCase.createTempFolder();
            fsn = 'cosine';
            [~, x] = ndgrid(1 : testCase.N, 0 : testCase.N - 1, 1 : testCase.N);
            im = single(1000 + 500 * cos(2 * pi * testCase.K * x / testCase.N));
            writetiff(im, [dataPath, fsn, '.tif']);
        end
    end
end
