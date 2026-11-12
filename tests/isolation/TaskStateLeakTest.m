classdef TaskStateLeakTest < pktest.TestCase
    % task functions leave the process as they found it, so that the next
    % task in the same process (the job framework's workers, local runs)
    % sees what a new process would: warning states, the global random
    % numbers, and caches that notice changed files (see tests/isolation and
    % speed/reports/isolation.md)

    properties (TestParameter)
        % wrappers that turned the directory exists warning off for good,
        % called on an empty folder (they stop or fail after turning it off)
        mkdirWrapper = struct( ...
            'crop', @(d) XR_crop_dataset(d, [1, 1, 1, 2, 2, 2], 'parseCluster', false), ...
            'resample', @(d) XR_resample_dataset(d, [2, 2, 2], 'parseCluster', false), ...
            'deskewRotate', @(d) XR_deskew_rotate_data_wrapper({d}, 'parseCluster', false), ...
            'decon', @(d) XR_decon_data_wrapper({d}, 'psfFullpaths', {[d, 'psf.tif']}, 'parseCluster', false), ...
            'fft', @(d) XR_fftSpectrumComputingWrapper({d}, 'parseCluster', false), ...
            'resaveZarr', @(d) XR_resave_zarr_wrapper(d, 'parseCluster', false), ...
            'rotate', @(d) XR_RotateFrame3D({[d, 'missing.tif']}, 0.108, 0.5), ...
            'imaris', @(d) XR_imaris_conversion_data_wrapper({d}, 'parseCluster', false), ...
            'pipeline', @(d) XR_microscopeAutomaticProcessing({d}, 'deskew', false, 'rotate', false, ...
                'stitch', false, 'parseCluster', false, 'streaming', false))
    end

    methods (Test, TestTags = {'Unit'})
        function wrapperRestoresTheDirectoryExistsWarning(testCase, mkdirWrapper)
            dataPath = testCase.createTempFolder();
            warning('on', 'MATLAB:MKDIR:DirectoryExists');
            try
                evalc('mkdirWrapper(dataPath)');
            catch
            end
            s = warning('query', 'MATLAB:MKDIR:DirectoryExists');
            testCase.verifyEqual(s.state, 'on', 'The wrapper left the directory exists warning off.');
        end

        function rsyncShufflesWithoutTheGlobalRandomNumbers(testCase)
            % the wrapper shuffled its files with rng(1), which reset the
            % global random numbers of the process; the files are shuffled
            % as before
            testCase.assumeFalse(ispc, 'rsync is not available on Windows.');
            src = [testCase.createTempFolder(), 'src'];
            mkdir(src);
            nF = 40;
            for i = 1 : nF
                fclose(fopen(sprintf('%s/f%04d.txt', src, i), 'w'));
            end
            dest = testCase.createTempFolder();
            testCase.applyFixture(matlab.unittest.fixtures.CurrentFolderFixture(pwd));
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            rng(7, 'twister');
            before = rng;
            evalc('XR_parallel_rsync_wrapper(src, dest, ''parseCluster'', false)');
            testCase.verifyEqual(rng, before, 'The wrapper changed the global random numbers.');

            % batches of one file, in the order of rng(1)'s permutation
            calls = spy.getCalls();
            testCase.assertNotEmpty(calls, 'The job framework was not called.');
            % (the batch script, then pairs of a file and its folder, quoted)
            files = cellfun(@(x) regexp(x, '''([^'']*)''', 'tokens'), calls{1}.funcStrs, 'unif', 0);
            files = cellfun(@(x) x{2}{1}, files, 'unif', 0);
            d = dir([src, '/*.txt']);
            names = sort(strcat(src, '/', {d.name}'));
            rng(1);
            expected = names(randperm(nF));
            testCase.verifyEqual(files(:), expected(:), 'The files are not shuffled as before.');
        end

        function ownRandomStreamGivesThePermutationOfRng1(testCase)
            % XR_parallel_rsync_wrapper's shuffle: a stream of its own seeded
            % with 1 gives the permutation of rng(1) and randperm
            for n = [1, 7, 1000, 123457]
                rng(1);
                expected = randperm(n);
                testCase.verifyEqual(randperm(RandStream('twister', 'Seed', 1), n), expected, sprintf('n = %d', n));
            end
        end

        function uuidWithoutJavaKeepsTheGlobalRandomNumbers(testCase)
            % without Java (compiled and MATLAB jobs), get_uuid took random
            % numbers from the global stream (Windows: after rng('shuffle');
            % elsewhere when uuidgen fails)
            testCase.applyFixture(pktest.FunctionStubFixture('usejava', false));
            testCase.applyFixture(pktest.FunctionStubFixture('system', 1, ''));
            rng(7, 'twister');
            before = rng;
            uuid = get_uuid();
            testCase.verifyEqual(rng, before, 'get_uuid changed the global random numbers.');
            testCase.verifyNotEmpty(uuid);
            testCase.verifyTrue(all(isstrprop(uuid, 'digit')), 'The fallback uuid is not a number.');
        end
    end

    methods (Test, TestTags = {'Integration'})
        function fscAnalysisReadsAChangedFileAgain(testCase)
            % XR_one_image_FSC_analysis kept the image of the last file by
            % its name only, so a file written again under the same name
            % (e.g., a new decon result) was analyzed as the old one
            testCase.assumeTrue(~isempty(ver('curvefit')) && license('test', 'Curve_Fitting_Toolbox'), ...
                'FSC resolution estimation needs the Curve Fitting Toolbox.');
            dataPath = testCase.createTempFolder();
            fn = [dataPath, 'im.tif'];
            args = {'xyPixelSize', 0.1, 'dz', 0.1, 'dr', 2, 'dtheta', pi / 6, 'halfSize', [12, 12, 12]};
            [imA, imB] = testCase.twoImages();
            writetiff(imA, fn, Compression='none');
            fscA = XR_one_image_FSC_analysis(fn, args{:});
            % same size, written later
            pause(1.1);
            writetiff(imB, fn, Compression='none');
            fscB = XR_one_image_FSC_analysis(fn, args{:});
            expected = XR_one_image_FSC_analysis(imB, args{:});
            testCase.verifyEqual(fscB, expected, 'The changed file was not read again.');
            testCase.verifyNotEqual(fscA.fsc, expected.fsc, 'The two images give the same FSC.');
        end

        function fscAnalysisSeesAChangedSizeAtTheSameTime(testCase)
            % a file written again within the same second (modification
            % times have a resolution of 1 s) is told apart by its size
            testCase.assumeTrue(~isempty(ver('curvefit')) && license('test', 'Curve_Fitting_Toolbox'), ...
                'FSC resolution estimation needs the Curve Fitting Toolbox.');
            testCase.assumeFalse(ispc, 'Needs touch -r to set the modification time.');
            dataPath = testCase.createTempFolder();
            fn = [dataPath, 'im.tif'];
            refFn = [dataPath, 'ref.tif'];
            args = {'xyPixelSize', 0.1, 'dz', 0.1, 'dr', 2, 'dtheta', pi / 6, 'halfSize', [12, 12, 12]};
            [imA, imB] = testCase.twoImages();
            writetiff(imA, fn, Compression='none');
            XR_one_image_FSC_analysis(fn, args{:});
            copyfile(fn, refFn);
            writetiff(imB, fn, Compression='lzw');
            status = system(sprintf('touch -r "%s" "%s"', refFn, fn));
            testCase.assumeEqual(status, 0, 'touch failed.');
            a = dir(fn);
            b = dir(refFn);
            testCase.assumeEqual(a.datenum, b.datenum);
            testCase.assumeNotEqual(a.bytes, b.bytes);

            fsc = XR_one_image_FSC_analysis(fn, args{:});
            testCase.verifyEqual(fsc, XR_one_image_FSC_analysis(imB, args{:}), 'The changed file was not read again.');
        end

        function fscFrameLeavesNoImageInTheCache(testCase)
            % the frame task read its image into XR_one_image_FSC_analysis's
            % cache, which kept it (and its memory) after the task: a later
            % call on a file with the same name, size and modification time
            % got it
            testCase.assumeTrue(~isempty(ver('curvefit')) && license('test', 'Curve_Fitting_Toolbox'), ...
                'FSC resolution estimation needs the Curve Fitting Toolbox.');
            testCase.assumeFalse(ispc, 'Needs touch -r to set the modification time.');
            dataPath = testCase.createTempFolder();
            fn = [dataPath, 'im.tif'];
            refFn = [dataPath, 'ref.tif'];
            args = {'xyPixelSize', 0.1, 'dz', 0.1, 'dr', 2, 'dtheta', pi / 6, 'halfSize', [12, 12, 12]};
            [imA, imB] = testCase.twoImages();
            writetiff(imA, fn, Compression='none');
            evalc('XR_one_image_FSC_analysis_frame(fn, [dataPath, ''fsc.mat''], args{:})');
            % the other image under the same name, size and modification time
            copyfile(fn, refFn);
            writetiff(imB, fn, Compression='none');
            status = system(sprintf('touch -r "%s" "%s"', refFn, fn));
            testCase.assumeEqual(status, 0, 'touch failed.');
            a = dir(fn);
            b = dir(refFn);
            testCase.assumeEqual([a.datenum, a.bytes], [b.datenum, b.bytes]);

            fsc = XR_one_image_FSC_analysis(fn, args{:});
            testCase.verifyEqual(fsc, XR_one_image_FSC_analysis(imB, args{:}), ...
                'The frame left the image in the cache.');
        end
    end

    methods (Access = private)
        function [imA, imB] = twoImages(~)
            % two noisy textures of the same size (uncompressed tiffs of
            % the same size on disk)
            rs = RandStream('twister', 'Seed', 3);
            b = imgaussfilt3(randn(rs, 32, 32, 32, 'single'), 1.5);
            signal = 1000 + 300 * b / std(b, 0, 'all');
            imA = single(signal + 20 * randn(rs, 32, 32, 32, 'single'));
            imB = single(signal + 400 * randn(rs, 32, 32, 32, 'single'));
        end
    end
end
