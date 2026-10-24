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

    properties (TestParameter)
        % {image size, outSize, frame options}: resized image (by [2.16,
        % 2.16, 6]) larger than outSize in all dimensions (cropped), smaller
        % in y, smaller in all (padded), with the 3D stack, background
        % subtraction, cubic interpolation, and an even size equal to it
        spectrumCase = struct( ...
            'crop', {{[31, 30, 12], [61, 61, 61], {}}}, ...
            'padY', {{[28, 30, 12], [61, 61, 61], {}}}, ...
            'padAll', {{[28, 30, 12], [81, 81, 81], {}}}, ...
            'stack', {{[31, 30, 12], [61, 61, 61], {'save3DStack', true}}}, ...
            'background', {{[31, 30, 12], [61, 61, 61], {'background', 1200}}}, ...
            'cubic', {{[31, 30, 12], [61, 61, 61], {'interpMethod', 'cubic'}}}, ...
            'evenEqual', {{[30, 30, 12], [65, 65, 72], {}}})
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

    methods (Test, TestTags = {'Integration'})
        function matchesWholeSpectrumComputation(testCase, spectrumCase)
            % only the region of the resized image that is kept is resized,
            % and the MIPs and central slices come from the unshifted |FFT|
            % with the normalization and gamma applied to them: the same
            % values as from the whole resized image and shifted spectrum
            [sz, outSize, opts] = spectrumCase{:};
            dataPath = testCase.createTempFolder();
            im = uint16(1000 + 300 * pktest.texture(sz, 1.5, [0, 1]) + 30 * randn(sz));
            writetiff(im, [dataPath, 'im.tif']);
            p = struct('background', 0, 'interpMethod', 'linear', 'save3DStack', false);
            for i = 1 : 2 : numel(opts)
                p.(opts{i}) = opts{i + 1};
            end
            [mips, slices, stack] = FFTSpectrumTest.referenceSpectrum(im, 0.108, 0.3, 0.05, outSize, ...
                p.interpMethod, p.background);

            % the resize must not fall back to imresize3 of the whole image
            testCase.applyFixture(pktest.FunctionSpyFixture('imresize3'));
            mkdir([dataPath, 'FFT']);
            XR_fftSpectrumComputingFrame([dataPath, 'im.tif'], [dataPath, 'FFT/im.tif'], 'xyPixelSize', 0.108, ...
                'dz', 0.3, 'outPixelSize', 0.05, 'outSize', outSize, opts{:});
            axisNames = {'y', 'x', 'z'};
            planeNames = {'xy', 'yz', 'xz'};
            for k = 1 : 3
                testCase.verifyEqual(readtiff(sprintf('%sFFT/MIPs/im_MIP_%s.tif', dataPath, axisNames{k})), mips{k}, ...
                    sprintf('MIP along %s', axisNames{k}));
                testCase.verifyEqual(readtiff(sprintf('%sFFT/central_slices/im_center_%s.tif', dataPath, planeNames{k})), ...
                    slices{k}, sprintf('Central %s slice', planeNames{k}));
            end
            if p.save3DStack
                testCase.verifyEqual(readtiff([dataPath, 'FFT/im.tif']), stack);
            end
        end

        function evenSizeThatNeedsCropOrPadErrors(testCase)
            % the center of an even outSize was never defined: cropping or
            % padding the resized image to it errors
            dataPath = testCase.createTempFolder();
            writetiff(uint16(pktest.texture([20, 20, 10])), [dataPath, 'im.tif']);
            mkdir([dataPath, 'FFT']);
            testCase.verifyError(@() XR_fftSpectrumComputingFrame([dataPath, 'im.tif'], [dataPath, 'FFT/im.tif'], ...
                'xyPixelSize', 0.1, 'dz', 0.2, 'outPixelSize', 0.1, 'outSize', [33, 32, 33]), ?MException);
            testCase.verifyFalse(exist([dataPath, 'FFT/im.tif'], 'file') > 0);
        end
    end

    methods (Static, Access = private)
        function [mips, slices, stack] = referenceSpectrum(im, xy, dz, px, N, interpMethod, background)
            % XR_fftSpectrumComputingFrame before the region crop and the
            % unshifted spectrum: resize the whole image, pad/crop to N, and
            % normalize the shifted |FFT|
            if background ~= 0
                im = single(im) - background;
                im = im .* (im > 0);
            end
            im = imresize3(single(im), round(size(im) .* [xy, xy, dz] ./ px), 'method', interpMethod);
            sz = size(im);
            hfN = (N - 1) / 2;
            if any(sz ~= N)
                if any(sz < N)
                    im = padarray(im, max(0, floor((N - sz) / 2)), 0, 'pre');
                    im = padarray(im, max(0, ceil((N - sz) / 2)), 0, 'post');
                end
                c = round((size(im) + 1) / 2);
                im = im(c(1) - hfN(1) : c(1) + hfN(1), c(2) - hfN(2) : c(2) + hfN(2), c(3) - hfN(3) : c(3) + hfN(3));
            end
            im = abs(fftshift(fftn(im)));
            im = im ./ max(im(:));
            stack = log10(single(im));
            im = single(im .^ 0.5);
            mips = {squeeze(max(im, [], 1)), squeeze(max(im, [], 2)), squeeze(max(im, [], 3))};
            [ny, nx, nz] = size(im);
            slices = {im(:, :, round((nz + 1) / 2)), squeeze(im(:, round((nx + 1) / 2), :)), ...
                squeeze(im(round((ny + 1) / 2), :, :))};
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
