classdef UnmixChannelsTest < pktest.TestCase
    % channel unmixing (demo_useful_tools.m, "Channel Unmixing" sections):
    % out = max(0, sum_c factor(c) * (ch_c - constBackground(c))) + constOffset,
    % computed by unmix_channels_sum (unmix_channels_sum_mex in one pass, or
    % the MATLAB code where the mex is missing)

    properties (Constant)
        Fsn = 'Scan_Iter_0000_0000_%s_ch0_CAM1_stack0000_488nm_0000000msec_0106060251msecAbs_000x_000y_000z_0000t'
    end

    properties (TestParameter)
        outClass = {'uint8', 'uint16', 'single', 'double'}
    end

    methods (Test, TestTags = {'Unit'})
        function mexMatchesTheMatlabCode(testCase, outClass)
            % uint8, uint16, single and double channels (single and double
            % with NaN, Inf, -0, huge values and values that round at .5),
            % mixed, with backgrounds and offsets (also 0, NaN and single),
            % factors and backgrounds that single can't represent, 0 and
            % Inf; 3D, 2D, a vector, sizes around the mex's blocks of 8192,
            % more blocks than threads, and empty
            testCase.assumeMexAvailable('unmix_channels_sum_mex');
            sizes = {[23, 19, 11], [17, 13], [1, 9], [8193, 1], [8191, 2], [0, 4], [1, 1], [700, 600, 3]};
            classes = {{'uint16', 'uint16'}, {'uint8', 'single', 'double'}, {'single'}, {'double', 'uint16'}, {'single', 'uint8'}};
            factors = {[1, -0.1], [0.1, -2.345678901234, 1/3], 1e3, single([-0.25, 1]), [1, 0.5]};
            backgrounds = {[100, 99.7], [0, 99.7, -1/3], [], single([99.5, 0]), [NaN, 1/3]};
            offsets = {1, -3.3, [], single(1.7), 0, NaN, 65535};
            cases = {};
            for i = 1 : numel(sizes)
                for k = 1 : numel(classes)
                    ims = cellfun(@(c) testCase.randomChannel(sizes{i}, c), classes{k}, 'unif', 0);
                    o = offsets{mod(i + k, numel(offsets)) + 1};
                    cases{end + 1} = {ims, factors{k}, backgrounds{k}, o}; %#ok<AGROW>
                end
            end
            outs = cell(size(cases));
            for i = 1 : numel(cases)
                outs{i} = unmix_channels_sum(cases{i}{:}, outClass);
            end
            % the MATLAB code: the mex hidden by a stub that errors
            testCase.applyFixture(pktest.FunctionStubFixture('unmix_channels_sum_mex'));
            testCase.expectMexFallback('unmix_channels_sum_mex');
            for i = 1 : numel(cases)
                expected = unmix_channels_sum(cases{i}{:}, outClass);
                name = sprintf('case %d (%s, %s)', i, mat2str(size(cases{i}{1}{1})), ...
                    strjoin(cellfun(@class, cases{i}{1}, 'unif', 0), ', '));
                testCase.verifyClass(outs{i}, outClass, name);
                testCase.verifySize(outs{i}, size(expected), name);
                if strcmp(outClass, 'single') || strcmp(outClass, 'double')
                    % the same bits (NaN, -0)
                    bits = sprintf('uint%d', 8 * (4 + 4 * strcmp(outClass, 'double')));
                    testCase.verifyEqual(typecast(outs{i}(:), bits), typecast(expected(:), bits), name);
                else
                    testCase.verifyEqual(outs{i}, expected, name);
                end
            end
        end

        function partialSumsGiveTheSameResult(testCase)
            % channels added one at a time to a partial sum (dtype empty),
            % then clipped, offset and cast (no channels and the sum), give
            % the result of one call, bit for bit: uint8, uint16, double and
            % single channels (single ones are added in place by MATLAB),
            % with the mex and with the MATLAB code
            ims = {testCase.randomChannel([23, 19, 11], 'uint16'), testCase.randomChannel([23, 19, 11], 'single'), ...
                testCase.randomChannel([23, 19, 11], 'uint8'), testCase.randomChannel([23, 19, 11], 'double')};
            factors = [1, -0.1, 1/3, -2.345678901234];
            bgs = [100, 99.7, 0, 1/3];
            outClasses = {'uint16', 'single', 'uint8', 'double'};
            for useMex = [true, false]
                if ~useMex
                    testCase.applyFixture(pktest.FunctionStubFixture('unmix_channels_sum_mex'));
                    testCase.expectMexFallback('unmix_channels_sum_mex');
                end
                for k = 1 : numel(outClasses)
                    expected = unmix_channels_sum(ims, factors, bgs, 1.7, outClasses{k});
                    im = [];
                    for c = 1 : numel(ims)
                        im = unmix_channels_sum(ims(c), factors(c), bgs(c), [], '', im);
                    end
                    testCase.verifyClass(im, 'single');
                    out = unmix_channels_sum({}, [], [], 1.7, outClasses{k}, im);
                    name = sprintf('%s, mex %d', outClasses{k}, useMex);
                    testCase.verifyClass(out, outClasses{k}, name);
                    if any(strcmp(outClasses{k}, {'single', 'double'}))
                        bits = sprintf('uint%d', 8 * (4 + 4 * strcmp(outClasses{k}, 'double')));
                        testCase.verifyEqual(typecast(out(:), bits), typecast(expected(:), bits), name);
                    else
                        testCase.verifyEqual(out, expected, name);
                    end
                end
            end
        end

        function sumRoundsAsMatlab(testCase)
            % products with a double factor and a double offset are rounded
            % to single once (not computed in single), sums in single, cast
            % rounds half away from zero (0.49999997 is 0, uint16() gives 1)
            % and saturates, NaN and negative sums are 0
            x = single([3, 0.49999997, 0.5, 1.5, 2.5, 65534.5, 65535.5, Inf, NaN, -5, 1e-3]);
            out = unmix_channels_sum({x}, 1, [], [], 'uint16');
            testCase.verifyEqual(out, uint16([3, 0, 1, 2, 3, 65535, 65535, 65535, 0, 0, 0]));
            out = unmix_channels_sum({x}, 1, [], [], 'uint8');
            testCase.verifyEqual(out, uint8([3, 0, 1, 2, 3, 255, 255, 255, 0, 0, 0]));
            x = single(rand(1, 1000) * 4000);
            out = unmix_channels_sum({x}, 0.1, [], [], 'single');
            testCase.verifyEqual(out, single(double(x) * 0.1));
            testCase.verifyNotEqual(out, x * single(0.1), 'The product was computed in single.');
            out = unmix_channels_sum({x, x}, [1, -0.1], [], [], 'single');
            testCase.verifyEqual(out, single(double(x) + double(single(double(x) * -0.1))));
            x = single(rand(1, 1000) * 1e-3);
            out = unmix_channels_sum({x}, 1, [], 0.1, 'single');
            testCase.verifyEqual(out, single(double(x) + 0.1));
            testCase.verifyNotEqual(out, x + single(0.1), 'The offset was added in single.');
        end

        function mexRejectsInvalidInput(testCase)
            testCase.assumeMexAvailable('unmix_channels_sum_mex');
            ims = {zeros(4, 5, 'uint16'), zeros(4, 5, 'uint16')};
            testCase.verifyError(@() unmix_channels_sum_mex(ims, [1, 2], [0, 0], 0), 'unmix_channels:inputError');
            testCase.verifyError(@() unmix_channels_sum_mex({}, [], [], 0, 'uint16'), 'unmix_channels:inputError');
            testCase.verifyError(@() unmix_channels_sum_mex(ims{1}, 1, 0, 0, 'uint16'), 'unmix_channels:inputError');
            testCase.verifyError(@() unmix_channels_sum_mex(ims, 1, [0, 0], 0, 'uint16'), 'unmix_channels:inputError');
            testCase.verifyError(@() unmix_channels_sum_mex(ims, [1, 2], 0, 0, 'uint16'), 'unmix_channels:inputError');
            testCase.verifyError(@() unmix_channels_sum_mex(ims, single([1, 2]), [0, 0], 0, 'uint16'), 'unmix_channels:inputError');
            testCase.verifyError(@() unmix_channels_sum_mex(ims, [1, 2], [0, 0], [0, 1], 'uint16'), 'unmix_channels:inputError');
            testCase.verifyError(@() unmix_channels_sum_mex({ims{1}, zeros(5, 4, 'uint16')}, [1, 2], [0, 0], 0, 'uint16'), ...
                'unmix_channels:inputError');
            testCase.verifyError(@() unmix_channels_sum_mex({ims{1}, zeros(4, 5, 1, 2, 'uint16')}, [1, 2], [0, 0], 0, 'uint16'), ...
                'unmix_channels:inputError');
            testCase.verifyError(@() unmix_channels_sum_mex({int16(ims{1})}, 1, 0, 0, 'uint16'), 'unmix_channels:dataTypeError');
            testCase.verifyError(@() unmix_channels_sum_mex({complex(single(ims{1}))}, 1, 0, 0, 'uint16'), 'unmix_channels:dataTypeError');
            testCase.verifyError(@() unmix_channels_sum_mex(ims, [1, 2], [0, 0], 0, 'int16'), 'unmix_channels:dataTypeError');
            testCase.verifyError(@() unmix_channels_sum_mex(ims, [1, 2], [0, 0], 0, 1), 'unmix_channels:inputError');
            % partial sums: single, of the channels' size; no channels need one
            testCase.verifyError(@() unmix_channels_sum_mex(ims, [1, 2], [0, 0], 0, '', zeros(4, 5)), 'unmix_channels:dataTypeError');
            testCase.verifyError(@() unmix_channels_sum_mex(ims, [1, 2], [0, 0], 0, '', zeros(5, 4, 'single')), ...
                'unmix_channels:inputError');
            testCase.verifyError(@() unmix_channels_sum_mex({}, [], [], 0, 'uint16', []), 'unmix_channels:inputError');
            testCase.verifyError(@() unmix_channels_sum_mex(ims, [1, 2], [0, 0], 0, 'uint16', [], 1), 'unmix_channels:inputError');
        end

        function otherClassesUseTheMatlabCode(testCase)
            % other channel and result classes are computed in MATLAB, without
            % a mex fallback
            x = int16(randi([-100, 3000], 9, 7, 3));
            y = uint16(randi(3000, 9, 7, 3));
            out1 = unmix_channels_sum({x, y}, [1, -0.1], [10, 5], 2, 'uint16');
            out2 = unmix_channels_sum({y, y}, [1, -0.1], [], 2, 'int32');
            testCase.verifyEmpty(mex_fallback('-log'));
            testCase.verifyEqual(out1, unmix_channels_sum({single(x), y}, [1, -0.1], [10, 5], 2, 'uint16'));
            testCase.verifyEqual(out2, cast(unmix_channels_sum({y, y}, [1, -0.1], [], 2, 'single'), 'int32'));
        end
    end

    methods (Test, TestTags = {'Integration'})
        function linearUnmixingAsInDemo(testCase)
            [dataPath, chA, chB] = testCase.twoChannelData(false);

            XR_unmix_channels_data_wrapper(dataPath, 'unmixFactors', [1, -0.1], 'mode', 'linear', ...
                'channelPatterns', {'CamA', 'CamB'}, 'channelInd', 1, 'constBackground', [100, 100], ...
                'constOffset', 1, 'save16bit', true, 'parseCluster', false, 'masterCompute', true);

            out = readtiff([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.tif']);
            expected = uint16(max(0, (single(chA) - 100) - 0.1 * (single(chB) - 100)) + 1);
            testCase.verifyClass(out, 'uint16');
            testCase.verifyEqual(out, expected);
            % only the unmixed channel is written
            testCase.verifyNumElements(dir([dataPath, 'Unmixed/*.tif']), 1);
        end

        function unmixSecondChannelAsSingle(testCase)
            [dataPath, chA, chB] = testCase.twoChannelData(false);

            XR_unmix_channels_data_wrapper(dataPath, 'unmixFactors', [-0.25, 1], 'channelInd', 2, ...
                'channelPatterns', {'CamA', 'CamB'}, 'save16bit', false, 'parseCluster', false, 'masterCompute', true);

            out = readtiff([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamB'), '.tif']);
            testCase.verifyClass(out, 'single');
            testCase.verifyImagesClose(out, max(0, single(chB) - 0.25 * single(chA)), 'relTol', 1e-6);
        end

        function flatFieldCorrectionBeforeUnmixing(testCase)
            % 1f40f4c added flat-field correction to unmixing; 824f476 fixed
            % flat-field/background images that are not single (as read
            % from uint16 tiffs), which made the mex reject them and fall
            % back to the slow MATLAB path
            [dataPath, chA, chB] = testCase.twoChannelData(false);
            ffPath = [dataPath, 'FF/'];
            mkdir(ffPath);
            [~, x] = ndgrid(1 : 48, 1 : 40);
            bk = uint16(100 * ones(48, 40));
            ffs = {uint16(100 + 200 + 800 * x / 40), uint16(100 + 1000 - 500 * x / 40)};
            ffFns = {[ffPath, 'ff_CamA.tif'], [ffPath, 'ff_CamB.tif']};
            bkFns = {[ffPath, 'bk_CamA.tif'], [ffPath, 'bk_CamB.tif']};
            for c = 1 : 2
                writetiff(ffs{c}, ffFns{c});
                writetiff(bk, bkFns{c});
            end

            consoleOut = evalc(['XR_unmix_channels_data_wrapper(dataPath, ''unmixFactors'', [1, -0.1], ', ...
                '''channelPatterns'', {''CamA'', ''CamB''}, ''FFCorrection'', true, ''lowerLimit'', 0.4, ', ...
                '''FFImagePaths'', ffFns, ''backgroundPaths'', bkFns, ''parseCluster'', false, ''masterCompute'', true)']);

            out = readtiff([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.tif']);
            corrected = cell(1, 2);
            raws = {chA, chB};
            for c = 1 : 2
                ls = single(ffs{c}) - single(bk);
                ls = max(ls / max(ls(:)), 0.4);
                corrected{c} = max((single(raws{c}) - single(bk)) ./ ls, 0);
            end
            expected = max(0, corrected{1} - 0.1 * corrected{2});
            % +/-1 from uint16 rounding
            testCase.verifyImagesClose(out, expected, 'absTol', 1, 'relTol', 1e-3);
            testCase.verifyFalse(contains(consoleOut, 'must be single'), ...
                'Flat-field correction fell back to the MATLAB implementation.');
        end

        function memAllocateIsInGigabytes(testCase)
            % 979fa55 regression: memAllocate (GB) was computed in bytes, so
            % every Slurm job requested the maximum number of CPUs
            dataPath = testCase.twoChannelData(false);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);

            XR_unmix_channels_data_wrapper(dataPath, 'unmixFactors', [1, -0.1], 'channelPatterns', {'CamA', 'CamB'}, ...
                'parseCluster', false, 'masterCompute', true);

            testCase.verifyEqual(spy.getOption(1, 'memAllocate'), prod([48, 40, 24]) * 4 * 2.5 / 1024 ^ 3, 'RelTol', 1e-12);
        end

        function imagesWithoutAllChannelsAreLeftOut(testCase)
            % time points without all channels failed the grouping by
            % channels; their images are left out with a warning, as the
            % wrapper says (and the second folder was grouped with the
            % channels of the first one's files)
            dataPath1 = testCase.twoChannelData(false);
            dataPath2 = testCase.twoChannelData(false);
            for t = 1 : 2
                fn = strrep(sprintf(testCase.Fsn, 'CamA'), 'Iter_0000_', sprintf('Iter_%04d_', t));
                writetiff(zeros(48, 40, 24, 'uint16'), [dataPath2, fn, '.tif']);
            end
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);

            lastwarn('');
            XR_unmix_channels_data_wrapper({dataPath1, dataPath2}, 'unmixFactors', [1, -0.1], ...
                'channelPatterns', {'CamA', 'CamB'}, 'parseCluster', false, 'masterCompute', true);

            testCase.verifySubstring(lastwarn(), 'do not have all channels');
            calls = spy.getCalls();
            testCase.assertNumElements(calls, 1);
            % the first time point of each folder (paths compared with / and
            % single separators: the wrapper's folders are as dir gives
            % them, with \ on Windows)
            unify = @(p) regexprep(strrep(p, '\', '/'), '/+', '/');
            ins = unify(calls{1}.inputFullpaths);
            fnA = ['/', sprintf(testCase.Fsn, 'CamA'), '.tif'];
            testCase.assertNumElements(ins, 2);
            testCase.verifyEqual(ins(:), {[unify(dataPath1), fnA(2 : end)]; [unify(dataPath2), fnA(2 : end)]});
            testCase.verifyNumElements(calls{1}.funcStrs, 2);
            testCase.verifySubstring(unify(calls{1}.funcStrs{2}), [unify(dataPath2), sprintf(testCase.Fsn, 'CamB'), '.tif'], ...
                'Not unmixed with CamB.');
        end

        function zarrInputAndOutput(testCase)
            [dataPath, chA, chB] = testCase.twoChannelData(true);

            XR_unmix_channels_data_wrapper(dataPath, 'unmixFactors', [1, -0.1], 'channelPatterns', {'CamA', 'CamB'}, ...
                'constBackground', 100, 'zarrFile', true, 'saveZarr', true, 'parseCluster', false, 'masterCompute', true);

            out = readzarr([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.zarr']);
            testCase.verifyEqual(out, uint16(max(0, (single(chA) - 100) - 0.1 * (single(chB) - 100))));
        end

        function zarrFunctionMatchesFrameFunction(testCase)
            % XR_unmix_channels_zarr ("Channel Unmixing for Zarr Files")
            % processes in batches; small batches must give the same result
            [dataPath, chA, chB] = testCase.twoChannelData(true);
            zarrFns = cellfun(@(x) [dataPath, sprintf(testCase.Fsn, x), '.zarr'], {'CamA', 'CamB'}, 'unif', 0);

            XR_unmix_channels_zarr(zarrFns, [1, -0.1], 'resultDirName', 'Unmixed', 'mode', 'linear', ...
                'channelInd', 1, 'batchSize', [32, 32, 16], 'blockSize', [16, 16, 16], ...
                'parseCluster', false, 'masterCompute', true);

            out = readzarr([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.zarr']);
            expected = max(0, single(chA) - 0.1 * single(chB));
            testCase.verifyImagesClose(single(out), expected, 'absTol', 1, 'relTol', 1e-3);
        end

        function zarrFunctionKeepsResultOfAnotherRun(testCase)
            % a copy from another run (its own uuid) that finds the result
            % published replaced it with its own, incomplete result
            [dataPath, chA] = testCase.twoChannelData(true);
            zarrFns = cellfun(@(x) [dataPath, sprintf(testCase.Fsn, x), '.zarr'], {'CamA', 'CamB'}, 'unif', 0);
            outFn = [dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.zarr'];
            testCase.applyFixture(pktest.FrameworkSpyFixture(@() writezarr(chA, outFn)));
            XR_unmix_channels_zarr(zarrFns, [1, -0.1], 'resultDirName', 'Unmixed', 'mode', 'linear', ...
                'channelInd', 1, 'batchSize', [32, 32, 16], 'blockSize', [16, 16, 16], ...
                'parseCluster', false, 'uuid', 'test');
            testCase.verifyEqual(readzarr(outFn), chA, 'The published result was replaced.');
            testCase.verifyEmpty(dir([dataPath, 'Unmixed/*_test.zarr']), 'The temporary result was left behind.');
        end

        function largeFileResultsAreTheZarrFilesOfTheTasks(testCase)
            % with largeFile the tasks (XR_unmix_channels_zarr) write zarr,
            % also without saveZarr, but the wrapper gave the job framework
            % .tif results: it waited for files no task writes and ran each
            % task three times (about 6 minutes even in local mode)
            dataPath = testCase.twoChannelData(true);
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture);
            for saveZarr = [false, true]
                XR_unmix_channels_data_wrapper(dataPath, 'unmixFactors', [1, -0.1], 'channelPatterns', {'CamA', 'CamB'}, ...
                    'zarrFile', true, 'largeFile', true, 'saveZarr', saveZarr, 'parseCluster', false, 'masterCompute', true);
            end
            calls = spy.getCalls();
            testCase.assertNumElements(calls, 2);
            unify = @(p) regexprep(strrep(p, '\', '/'), '/+', '/');
            expected = unify([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.zarr']);
            for i = 1 : 2
                testCase.verifyEqual(unify(calls{i}.outputFullpaths), {expected}, sprintf('saveZarr %d', i - 1));
            end
        end

        function largeFileTasksRunOnce(testCase)
            % (as above) the job framework finds the result of a large-file
            % task after its first run
            [dataPath, chA, chB] = testCase.twoChannelData(true);
            consoleOut = evalc(['XR_unmix_channels_data_wrapper(dataPath, ''unmixFactors'', [1, -0.1], ', ...
                '''channelPatterns'', {''CamA'', ''CamB''}, ''zarrFile'', true, ''largeFile'', true, ', ...
                '''batchSize'', [32, 24, 16], ''blockSize'', [16, 8, 8], ''parseCluster'', false, ''masterCompute'', true)']);
            testCase.verifyEqual(count(consoleOut, 'with function XR_unmix_channels_zarr('), 1, 'The task did not run once.');
            out = readzarr([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.zarr']);
            testCase.verifyImagesClose(single(out), max(0, single(chA) - 0.1 * single(chB)), 'absTol', 1, 'relTol', 1e-3);
        end

        function gaussianUnmixingBlursTheOtherChannel(testCase)
            % 'gaussian' mode filters each channel with its sigma (one row
            % per channel, [sy, sx, sz]) before the weighted sum
            [dataPath, chA, chB] = testCase.twoChannelData(true);
            zarrFns = cellfun(@(x) [dataPath, sprintf(testCase.Fsn, x), '.zarr'], {'CamA', 'CamB'}, 'unif', 0);
            sigmas = [0, 0, 0; 1, 1, 1];

            XR_unmix_channels_zarr(zarrFns, [1, -0.1], 'mode', 'gaussian', 'unmixSigmas', sigmas, ...
                'channelInd', 1, 'parseCluster', false, 'masterCompute', true);

            % the output keeps the input data type (uint16)
            out = readzarr([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.zarr']);
            expected = max(0, single(chA) - 0.1 * imgaussfilt3(single(chB), [1, 1, 1]));
            testCase.verifyImagesClose(single(out), expected, 'absTol', 1, 'relTol', 1e-3);
        end

        function gaussianUnmixingThroughDataWrapper(testCase)
            % the wrapper used to format the sigma matrix (one row per channel) as a
            % row vector
            [dataPath, chA, chB] = testCase.twoChannelData(false);
            sigmas = [0, 0, 0; 1, 1, 1];

            XR_unmix_channels_data_wrapper(dataPath, 'unmixFactors', [1, -0.1], 'mode', 'gaussian', ...
                'unmixSigmas', sigmas, 'channelPatterns', {'CamA', 'CamB'}, 'save16bit', false, ...
                'parseCluster', false, 'masterCompute', true);

            out = readtiff([dataPath, 'Unmixed/', sprintf(testCase.Fsn, 'CamA'), '.tif']);
            expected = max(0, single(chA) - 0.1 * imgaussfilt3(single(chB), [1, 1, 1]));
            testCase.verifyImagesClose(out, expected, 'relTol', 1e-5);
        end

        function frameMatchesTheOriginalCode(testCase)
            % the frame function computes the sum in unmix_channels_sum, at
            % once (channels as read that take with the result at most 8
            % bytes per voxel) or one channel at a time: the same values as
            % the per-channel code it replaced, with flat field correction,
            % gaussian filtering, a channel with factor 0 (not read),
            % backgrounds, offsets, and single channels
            [dataPath, chA, chB] = testCase.twoChannelData(false);
            chC = uint16(pktest.texture([48, 40, 24], 3, [0, 900]));
            fsnC = sprintf(testCase.Fsn, 'CamC');
            writetiff(chC, [dataPath, fsnC, '.tif']);
            fns = cellfun(@(x) [dataPath, sprintf(testCase.Fsn, x), '.tif'], {'CamA', 'CamB', 'CamC'}, 'unif', 0);
            [~, x] = ndgrid(1 : 48, 1 : 40);
            ffFns = cell(1, 3);
            bkFns = cell(1, 3);
            for c = 1 : 3
                ffFns{c} = sprintf('%sff_%d.tif', dataPath, c);
                bkFns{c} = sprintf('%sbk_%d.tif', dataPath, c);
                writetiff(uint16(300 + 700 * x / 40 + 100 * c), ffFns{c});
                writetiff(uint16(90 + 5 * c * ones(48, 40)), bkFns{c});
            end
            raws = {chA, chB, chC};
            settings = {
                {'mode', 'linear', 'constBackground', [100, 99.5, 98], 'constOffset', 1, 'save16bit', true}
                {'mode', 'gaussian', 'unmixSigmas', [0, 0, 0; 1, 1, 1; 0.7, 0.7, 2], 'constBackground', [100, 99.5, 98], ...
                    'constOffset', 0.5, 'save16bit', false}
                {'mode', 'linear', 'FFCorrection', true, 'constOffset', -0.3, 'save16bit', true}
                {'mode', 'gaussian', 'unmixSigmas', [0, 0, 0; 1, 1, 1; 0, 0, 0], 'FFCorrection', true, ...
                    'constBackground', [3, 2, 1], 'save16bit', false}
                {'mode', 'linear', 'constBackground', [100, 99.5, 98], 'constOffset', -0.3, 'save16bit', false}
                {'mode', 'linear', 'save16bit', true}};
            factors = {[1, -0.1, 0], [1, -0.1, 0.37], [-0.25, 1, -1/3], [1, 0, -0.1234567], [1, -0.1, 0.37], ...
                [1, -0.1, 0.37]};
            for i = 1 : numel(settings)
                uuid = sprintf('t%d', i);
                resultDirName = sprintf('Unmixed%d', i);
                evalc(['XR_unmix_channels_frame(fns, factors{i}, ''resultDirName'', resultDirName, ''channelInd'', 1, ', ...
                    '''FFImagePaths'', ffFns, ''backgroundPaths'', bkFns, ''uuid'', uuid, settings{i}{:})']);
                out = readtiff([dataPath, resultDirName, '/', sprintf(testCase.Fsn, 'CamA'), '.tif']);
                expected = testCase.originalFrameUnmix(raws, factors{i}, ffFns, bkFns, settings{i}{:});
                testCase.verifyClass(out, class(expected), sprintf('setting %d', i));
                testCase.verifyEqual(out, expected, sprintf('setting %d', i));
            end

            % single channels (e.g., deconvolved data), with their backgrounds
            % subtracted in place
            singlePath = [dataPath, 'single/'];
            mkdir(singlePath);
            sfns = cellfun(@(x) [singlePath, sprintf(testCase.Fsn, x), '.tif'], {'CamA', 'CamB', 'CamC'}, 'unif', 0);
            sraws = cellfun(@(x) single(x) + 0.37, raws, 'unif', 0);
            for c = 1 : 3
                writetiff(sraws{c}, sfns{c});
            end
            for save16bit = [true, false]
                args = {'mode', 'linear', 'constBackground', [100, 99.7, 98], 'constOffset', 1, 'save16bit', save16bit};
                resultDirName = sprintf('UnmixedSingle%d', save16bit);
                evalc(['XR_unmix_channels_frame(sfns, [1, -0.1, 0.37], ''resultDirName'', resultDirName, ', ...
                    '''channelInd'', 1, ''uuid'', ''s'', args{:})']);
                out = readtiff([singlePath, resultDirName, '/', sprintf(testCase.Fsn, 'CamA'), '.tif']);
                testCase.verifyEqual(out, testCase.originalFrameUnmix(sraws, [1, -0.1, 0.37], {}, {}, args{:}), ...
                    sprintf('single channels, save16bit %d', save16bit));
            end
        end

        function zarrBlocksMatchTheOriginalCode(testCase)
            % the linear and gaussian zarr block functions with batches that
            % don't divide the image, as the per-channel code they replaced:
            % uint16 channels summed at once, filtered and single channels one
            % at a time
            [dataPath, chA, chB] = testCase.twoChannelData(true);
            zarrFns = cellfun(@(x) [dataPath, sprintf(testCase.Fsn, x), '.zarr'], {'CamA', 'CamB'}, 'unif', 0);
            factors = [1, -0.1];
            sigmas = [0, 0, 0; 1, 1, 1];
            XR_unmix_channels_zarr(zarrFns, factors, 'resultDirName', 'Linear', 'mode', 'linear', ...
                'channelInd', 1, 'batchSize', [32, 24, 16], 'blockSize', [16, 8, 8], ...
                'parseCluster', false, 'masterCompute', true);
            XR_unmix_channels_zarr(zarrFns, factors, 'resultDirName', 'Gaussian', 'mode', 'gaussian', ...
                'unmixSigmas', sigmas, 'channelInd', 1, 'batchSize', [32, 24, 16], 'blockSize', [16, 8, 8], ...
                'borderSize', [4, 4, 4], 'parseCluster', false, 'masterCompute', true);

            out = readzarr([dataPath, 'Linear/', sprintf(testCase.Fsn, 'CamA'), '.zarr']);
            expected = zeros(size(chA), 'single');
            expected = expected + single(chA) * factors(1);
            expected = expected + single(chB) * factors(2);
            testCase.verifyEqual(out, cast(max(0, expected), 'uint16'));

            % the gaussian blocks filter each batch with its border
            out = readzarr([dataPath, 'Gaussian/', sprintf(testCase.Fsn, 'CamA'), '.zarr']);
            [batchBBoxes, regionBBoxes, localBBoxes] = XR_zarrChunkCoordinatesExtraction(size(chA), ...
                'batchSize', [32, 24, 16], 'blockSize', [16, 8, 8], 'sameBatchSize', false, 'borderSize', [4, 4, 4]);
            expected = zeros(size(chA), 'uint16');
            ims = {chA, chB};
            for i = 1 : size(batchBBoxes, 1)
                b = batchBBoxes(i, :);
                in_batch = zeros(b(4 : 6) - b(1 : 3) + 1, 'single');
                for c = 1 : 2
                    im_c = single(ims{c}(b(1) : b(4), b(2) : b(5), b(3) : b(6)));
                    if ~all(sigmas(c, :) == 0)
                        im_c = imgaussfilt3(im_c, sigmas(c, :));
                    end
                    in_batch = in_batch + im_c * factors(c);
                end
                in_batch = cast(max(0, in_batch), 'uint16');
                l = localBBoxes(i, :);
                r = regionBBoxes(i, :);
                expected(r(1) : r(4), r(2) : r(5), r(3) : r(6)) = in_batch(l(1) : l(4), l(2) : l(5), l(3) : l(6));
            end
            testCase.verifyEqual(out, expected);

            % single channels and result
            singleFns = strrep(zarrFns, '.zarr', '_single.zarr');
            writezarr(single(chA) + 0.37, singleFns{1}, blockSize=[16, 16, 16]);
            writezarr(single(chB) - 0.37, singleFns{2}, blockSize=[16, 16, 16]);
            XR_unmix_channels_zarr(singleFns, factors, 'resultDirName', 'LinearSingle', 'mode', 'linear', ...
                'channelInd', 1, 'batchSize', [32, 24, 16], 'blockSize', [16, 8, 8], ...
                'parseCluster', false, 'masterCompute', true);
            out = readzarr([dataPath, 'LinearSingle/', sprintf(testCase.Fsn, 'CamA'), '_single.zarr']);
            expected = zeros(size(chA), 'single');
            expected = expected + single(single(chA) + 0.37) * factors(1);
            expected = expected + single(single(chB) - 0.37) * factors(2);
            testCase.verifyEqual(out, max(0, expected));
        end
    end

    methods (Access = private)
        function [dataPath, chA, chB] = twoChannelData(testCase, zarrFile)
            % CamB bleeds into CamA
            dataPath = testCase.createTempFolder();
            chB = uint16(pktest.texture([48, 40, 24], 2, [100, 3000]));
            chA = uint16(pktest.texture([48, 40, 24], 2, [100, 1500]) + 0.1 * (single(chB) - 100));
            ims = {chA, chB};
            cams = {'CamA', 'CamB'};
            for c = 1 : 2
                fn = [dataPath, sprintf(testCase.Fsn, cams{c})];
                if zarrFile
                    writezarr(ims{c}, [fn, '.zarr'], blockSize=[16, 16, 16]);
                else
                    writetiff(ims{c}, [fn, '.tif']);
                end
            end
        end
    end

    methods (Static, Access = private)
        function [x] = randomChannel(sz, cls)
            % a channel image of class cls; single and double ones with NaN,
            % Inf, -0, huge values and values that round at .5
            switch cls
                case 'uint8'
                    x = uint8(randi([0, 255], sz));
                case 'uint16'
                    x = uint16(randi([0, 65535], sz));
                otherwise
                    x = (rand(sz) - 0.3) * 70000;
                    edge = [NaN, Inf, -Inf, -0, 0, 0.49999997, 0.5, 65535.5, 254.5, 255.5, 1e-40, 3e38, -3e38, 1e300];
                    n = min(numel(x), 3 * numel(edge));
                    x(randperm(numel(x), n)) = edge(mod(0 : n - 1, numel(edge)) + 1);
                    x = cast(x, cls);
            end
        end

        function [im] = originalFrameUnmix(raws, unmixFactors, ffFns, bkFns, varargin)
            % the per-channel code of XR_unmix_channels_frame before
            % unmix_channels_sum
            ip = inputParser;
            ip.addParameter('mode', 'linear');
            ip.addParameter('unmixSigmas', []);
            ip.addParameter('FFCorrection', false);
            ip.addParameter('constBackground', []);
            ip.addParameter('constOffset', []);
            ip.addParameter('save16bit', true);
            ip.parse(varargin{:});
            pr = ip.Results;

            im = zeros(size(raws{1}), 'single');
            for c = 1 : numel(raws)
                if unmixFactors(c) == 0
                    continue;
                end
                im_c = single(raws{c});
                if pr.FFCorrection
                    im_c = XR_LSFlatFieldCorrection(im_c, readtiff(ffFns{c}), readtiff(bkFns{c}), 'lowerLimit', 0.4, ...
                        'constOffset', 0);
                end
                if ~isempty(pr.constBackground)
                    im_c = im_c - pr.constBackground(c);
                end
                if strcmp(pr.mode, 'gaussian') && ~all(pr.unmixSigmas(c, :) == 0)
                    im_c = imgaussfilt3(im_c, pr.unmixSigmas(c, :));
                end
                im = im + im_c * unmixFactors(c);
            end
            im = max(0, im);
            if ~isempty(pr.constOffset) && pr.constOffset ~= 0
                im = im + pr.constOffset;
            end
            if pr.save16bit
                im = cast(im, 'uint16');
            end
        end
    end
end
