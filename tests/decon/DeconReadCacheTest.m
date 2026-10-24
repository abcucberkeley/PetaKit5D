classdef DeconReadCacheTest < pktest.TestCase
    % decon reads the generated PSF and the OMW back projector for each batch
    % of a large file (and each file of a job); decon_readtiff_cached keeps
    % the last ones read and reads a file again once its modification time or
    % size changes

    methods (TestMethodSetup)
        function clearCache(testCase)
            clear decon_readtiff_cached
            testCase.addTeardown(@() clear('decon_readtiff_cached'));
        end
    end

    methods (Test, TestTags = {'Unit'})
        function backProjectorIsReadOnce(testCase)
            % the second call doesn't read the file (readtiff stubbed)
            psfgenPath = testCase.createTempFolder();
            bpFn = sprintf('%s/PSF_back_projector_alpha_%0.6f_otf_cum_thresh_%0.6f.tif', psfgenPath, 0.005, 0.9);
            bp = rand(9, 8, 7, 'single');
            writetiff(bp, bpFn);
            args = {[], psfgenPath, 'PSF', 0.005, 0.9, false, [0.8, 1], 'PSF.tif'};
            [~, first] = evalc('omw_backprojector_files(args{:})');
            testCase.verifyEqual(first, bp);
            testCase.applyFixture(pktest.FunctionStubFixture('readtiff', zeros(2, 'single')));
            [~, second] = evalc('omw_backprojector_files(args{:})');
            testCase.verifyEqual(second, bp);
        end

        function changedFilesAreReadAgain(testCase)
            folder = testCase.createTempFolder();
            fn = [folder, 'psf.tif'];
            a = rand(9, 8, 7, 'single');
            writetiff(a, fn);
            testCase.verifyEqual(decon_readtiff_cached(fn), a);
            % another size
            b = rand(9, 8, 8, 'single');
            writetiff(b, fn);
            testCase.verifyEqual(decon_readtiff_cached(fn), b);
            testCase.verifyEqual(decon_readtiff_cached(fn), b);
            % the same file with another modification time is read again
            % (readtiff stubbed)
            java.io.File(fn).setLastModified(java.lang.System.currentTimeMillis() + 5000);
            testCase.applyFixture(pktest.FunctionStubFixture('readtiff', zeros(2, 'single')));
            testCase.verifyEqual(decon_readtiff_cached(fn), zeros(2, 'single'));
        end

        function severalFilesAreKept(testCase)
            folder = testCase.createTempFolder();
            fns = arrayfun(@(i) sprintf('%sf%d.tif', folder, i), 1 : 5, 'unif', 0);
            ims = cell(1, 5);
            for i = 1 : 5
                ims{i} = rand(4, 5, 3, 'single') + i;
                writetiff(ims{i}, fns{i});
                testCase.verifyEqual(decon_readtiff_cached(fns{i}), ims{i});
            end
            % the last four stay in memory, the first is read again
            testCase.applyFixture(pktest.FunctionStubFixture('readtiff', zeros(2, 'single')));
            for i = 2 : 5
                testCase.verifyEqual(decon_readtiff_cached(fns{i}), ims{i});
            end
            testCase.verifyEqual(decon_readtiff_cached(fns{1}), zeros(2, 'single'));
        end

        function missingFilesError(testCase)
            folder = testCase.createTempFolder();
            testCase.verifyError(@() decon_readtiff_cached([folder, 'missing.tif']), ?MException);
        end
    end
end
