classdef BlockMIPTest < pktest.TestCase
    % MIPs of large-file results from the blocks their tasks write:
    % MIP_zarr_block_mips (with range_max_3d_mex), MIP_zarr_save_block_mips
    % and MIP_zarr_merge_block_mips, which XR_MIP_zarr uses to make the same
    % MIP slabs and MIPs it computes from the zarr, and MIP_zarr_batch_size

    properties (Constant)
        ImSize = [37, 29, 23]
        % MIP batches, and blocks that are not aligned with them
        MIPBatchSize = [8, 6, 5]
        BlockSize = [7, 9, 4]
    end

    properties (TestParameter)
        dtype = {'uint16', 'single'}
    end

    methods (Test, TestTags = {'Unit'})
        function rangeMaxMatchesMaxOfEachRange(testCase, dtype)
            testCase.assumeMexAvailable('range_max_3d_mex');
            for sz = {[13, 11, 9], [17, 9], [1, 12, 7]}
                im = testCase.makeImage(sz{1}, dtype);
                for dim = 1 : 3
                    n = size(im, dim);
                    for starts = {1, unique([1, min(2, n), ceil(n / 2), n])}
                        out = range_max_3d_mex(im, dim, starts{1});
                        testCase.verifySameBits(out, testCase.rangeMax(im, dim, starts{1}), ...
                            sprintf('size %s, dim %d, starts %s', mat2str(sz{1}), dim, mat2str(starts{1})));
                    end
                end
            end
        end

        function rangeMaxRejectsInvalidInput(testCase)
            testCase.assumeMexAvailable('range_max_3d_mex');
            im = uint16(magic(4));
            testCase.verifyError(@() range_max_3d_mex(im, 1, [2, 3]), 'range_max_3d:inputError');
            testCase.verifyError(@() range_max_3d_mex(im, 1, [1, 3, 2]), 'range_max_3d:inputError');
            testCase.verifyError(@() range_max_3d_mex(im, 1, [1, 5]), 'range_max_3d:inputError');
            testCase.verifyError(@() range_max_3d_mex(im, 1, [1, 2.5]), 'range_max_3d:inputError');
            testCase.verifyError(@() range_max_3d_mex(im, 4, 1), 'range_max_3d:inputError');
            testCase.verifyError(@() range_max_3d_mex(int16(im), 1, 1), 'range_max_3d:inputError');
        end

        function mergedSlabsMatchTheImage(testCase, dtype)
            % blocks not aligned with the MIP batches, written by tasks in
            % an order that isn't their order in the image; for single, a
            % region of NaN and of -0 and 0, where max takes the first value
            im = testCase.makeImage(testCase.ImSize, dtype);
            flagFns = testCase.writeFlags(im, testCase.BlockSize, testCase.MIPBatchSize, 4);
            for j = 1 : 3
                slab = MIP_zarr_merge_block_mips(flagFns, j, size(im), testCase.MIPBatchSize, dtype);
                testCase.verifySameBits(slab, testCase.rangeMax(im, j, 1 : testCase.MIPBatchSize(j) : size(im, j)), ...
                    sprintf('slab along axis %d', j));
            end
        end

        function blockMIPsWithoutTheMexMatch(testCase, dtype)
            % platforms without range_max_3d_mex take the max of each part
            im = testCase.makeImage([20, 15, 13], dtype);
            bbox = [3, 2, 4, 22, 16, 16];
            expected = MIP_zarr_block_mips(im, bbox, [8, 6, 5], dtype);
            testCase.applyFixture(pktest.FunctionStubFixture('range_max_3d_mex'));
            evalc('mips = MIP_zarr_block_mips(im, bbox, [8, 6, 5], dtype)');
            for j = 1 : 3
                testCase.verifySameBits(mips{j}, expected{j}, sprintf('axis %d', j));
            end
        end

        function blockMIPsAreCastAsTheZarrStoresThem(testCase)
            % single blocks written to a uint16 zarr (large-file decon with
            % save16bit) are rounded and saturated as written
            folder = testCase.createTempFolder();
            block = single(randn([12, 10, 9]) * 40000 + 30000);
            block(1 : 5, 2, 3) = [0.5, 1.5, -3, NaN, 70000];
            block(:, 4, :) = NaN;
            zarrFn = [folder, 'out.zarr'];
            createzarr(zarrFn, dataSize=[12, 10, 9], blockSize=[4, 4, 4], dtype='uint16');
            writezarr(block, zarrFn, bbox=[1, 1, 1, 12, 10, 9], create=false);
            stored = readzarr(zarrFn);
            mips = MIP_zarr_block_mips(block, [1, 1, 1, 12, 10, 9], [5, 4, 4], 'uint16');
            for j = 1 : 3
                testCase.verifyEqual(mips{j}, testCase.rangeMax(stored, j, 1 : [5, 4, 4] * (j == 1 : 3)' : size(stored, j)), ...
                    sprintf('axis %d', j));
            end
        end

        function flagIsMovedIntoPlace(testCase)
            % the flag marks the task as done, so it must not be seen half
            % written
            folder = testCase.createTempFolder();
            moves = testCase.applyFixture(pktest.MovefileSpyFixture);
            im = uint16(randi(1000, [6, 5, 4]));
            mips = MIP_zarr_block_mips(im, [1, 1, 1, 6, 5, 4], [8, 8, 8], 'uint16');
            MIP_zarr_save_block_mips([folder, 'blocks_1_1.mat'], [1, 1, 1, 6, 5, 4], mips, [8, 8, 8], 1.5);
            testCase.verifyTrue(moves.wasMoved([folder, 'blocks_1_1.mat']), 'The flag was written in place.');
            testCase.verifyEqual(numel(dir([folder, '*.mat'])), 1, 'A temporary flag was left behind.');
            a = load([folder, 'blocks_1_1.mat']);
            testCase.verifyEqual(a.t, 1.5);
            testCase.verifyEqual(a.zMIPs{1}, max(im, [], 3));
        end

        function mergeNeedsMIPsOfAllBlocks(testCase)
            % flags of tasks without MIPs (older code), MIPs for other
            % batches, or blocks that don't cover the image give no slab, so
            % that the caller reads the image
            im = testCase.makeImage(testCase.ImSize, 'uint16');
            sz = size(im);
            flagFns = testCase.writeFlags(im, testCase.BlockSize, testCase.MIPBatchSize, 4);
            testCase.assertNotEmpty(MIP_zarr_merge_block_mips(flagFns, 1, sz, testCase.MIPBatchSize, 'uint16'));

            testCase.verifyEmpty(MIP_zarr_merge_block_mips(flagFns, 1, sz, [16, 6, 5], 'uint16'), 'other batch size');
            % MIPs saved for other batches along z only: their y slab would fit
            a = load(flagFns{1});
            a.mipBatchSize = [8, 6, 6];
            save(flagFns{1}, '-struct', 'a');
            testCase.verifyEmpty(MIP_zarr_merge_block_mips(flagFns, 1, sz, testCase.MIPBatchSize, 'uint16'), ...
                'flag for other batches');
            a.mipBatchSize = testCase.MIPBatchSize;
            save(flagFns{1}, '-struct', 'a');
            % the last flag without the MIP of its last block
            a = load(flagFns{end});
            yMIPs = a.yMIPs;
            a.yMIPs = yMIPs(1 : end - 1);
            save(flagFns{end}, '-struct', 'a');
            testCase.verifyEmpty(MIP_zarr_merge_block_mips(flagFns, 1, sz, testCase.MIPBatchSize, 'uint16'), ...
                'flag with a MIP missing');
            a.yMIPs = yMIPs;
            save(flagFns{end}, '-struct', 'a');
            testCase.verifyEmpty(MIP_zarr_merge_block_mips(flagFns, 1, sz, testCase.MIPBatchSize, 'single'), 'other data type');
            testCase.verifyEmpty(MIP_zarr_merge_block_mips(flagFns(2 : end), 1, sz, testCase.MIPBatchSize, 'uint16'), 'missing flag');
            testCase.verifyEmpty(MIP_zarr_merge_block_mips([flagFns; flagFns(1)], 1, sz, testCase.MIPBatchSize, 'uint16'), ...
                'repeated flag');

            a = load(flagFns{2});
            t = a.t;
            save(flagFns{2}, 't');
            testCase.verifyEmpty(MIP_zarr_merge_block_mips(flagFns, 1, sz, testCase.MIPBatchSize, 'uint16'), 'flag with only t');
            fclose(fopen(flagFns{2}, 'w'));
            testCase.verifyEmpty(MIP_zarr_merge_block_mips(flagFns, 1, sz, testCase.MIPBatchSize, 'uint16'), 'empty flag');
        end

        function mergeNeedsBlocksThatCoverTheImageOnce(testCase)
            % two overlapping blocks and a gap of the same size add up to
            % the image size
            im = uint16(randi(1000, [10, 8, 6]));
            folder = testCase.createTempFolder();
            bboxes = [1, 1, 1, 5, 8, 6; 1, 1, 1, 5, 8, 6];
            mips = [MIP_zarr_block_mips(im(1 : 5, :, :), bboxes(1, :), [4, 4, 4], 'uint16'); ...
                MIP_zarr_block_mips(im(1 : 5, :, :), bboxes(2, :), [4, 4, 4], 'uint16')];
            MIP_zarr_save_block_mips([folder, 'blocks_1_2.mat'], bboxes, mips, [4, 4, 4]);
            for j = 1 : 3
                testCase.verifyEmpty(MIP_zarr_merge_block_mips({[folder, 'blocks_1_2.mat']}, j, size(im), [4, 4, 4], 'uint16'), ...
                    sprintf('axis %d', j));
            end
        end

        function mergeKeepsToTheMemoryLimit(testCase)
            % the y slab and its coverage count (5 x 29 x 23 x 3 bytes) and
            % the blocks' MIPs (blocks overlap 10 batches along y in all:
            % 10 x 29 x 23 x 2 bytes) are 23345 bytes, 2.174e-5 GiB
            im = testCase.makeImage(testCase.ImSize, 'uint16');
            flagFns = testCase.writeFlags(im, testCase.BlockSize, testCase.MIPBatchSize, 4);
            args = {flagFns, 1, size(im), testCase.MIPBatchSize, 'uint16'};
            testCase.verifyEmpty(MIP_zarr_merge_block_mips(args{:}, 2.1e-5), 'over the limit');
            testCase.verifyNotEmpty(MIP_zarr_merge_block_mips(args{:}, 2.2e-5), 'under the limit');
            testCase.verifyNotEmpty(MIP_zarr_merge_block_mips(args{:}, NaN), 'unknown memory');
        end

        function batchSizeIsRoundedToBlocksAndCapped(testCase)
            testCase.verifyEqual(MIP_zarr_batch_size([5000, 3000, 900], [256, 256, 256], 'uint16'), [2048, 2048, 900]);
            testCase.verifyEqual(MIP_zarr_batch_size([5000, 3000, 900], [256, 300, 100], 'uint16'), [2048, 2100, 900]);
            testCase.verifyEqual(MIP_zarr_batch_size([5000, 3000, 900], [256, 256, 256], 'uint16', 'batchSize', [100, 100, 100]), ...
                [256, 256, 256]);
        end

        function batchSizeGrowsForSlabsOver100GB(testCase)
            % the batch grows by a block along the axis of the largest slab
            % until no slab is over 100 GiB (3 steps here), or the batch is
            % (21 steps, to a 104 GiB batch); values of XR_MIP_zarr before
            % the batch size moved to MIP_zarr_batch_size
            testCase.verifyEqual(MIP_zarr_batch_size([100000, 100000, 11000], [256, 256, 256], 'uint16'), [2304, 2304, 2304]);
            testCase.verifyEqual(MIP_zarr_batch_size([200000, 100000, 10000], [256, 256, 256], 'uint16'), [3584, 3584, 4352]);
            testCase.verifyEqual(MIP_zarr_batch_size([200000, 100000, 10000], [256, 256, 256], 'uint16', 'mipSlab', true), ...
                [2048, 2048, 2048]);
        end
    end

    methods (Test, TestTags = {'Integration'})
        function xrMIPZarrMergesTheBlockMIPs(testCase, dtype)
            % the same slab zarrs (data, chunks, data type) and MIPs as from
            % the zarr, without job framework tasks
            [zarrFns, im] = testCase.writeZarrCopies(dtype);
            flagFns = testCase.writeFlags(im, testCase.BlockSize, testCase.MIPBatchSize, 4);
            evalc('XR_MIP_zarr(zarrFns{1}, ''batchSize'', testCase.MIPBatchSize, ''parseCluster'', false)');
            spy = testCase.applyFixture(pktest.FrameworkSpyFixture());
            evalc('XR_MIP_zarr(zarrFns{2}, ''batchSize'', testCase.MIPBatchSize, ''parseCluster'', false, ''blockFlagFullpaths'', flagFns)');
            testCase.verifyEmpty(spy.getCalls(), 'The MIPs were computed by reading the zarr.');
            testCase.verifySameMIPs(zarrFns{1}, zarrFns{2});
        end

        function xrMIPZarrReadsTheZarrIfAFlagHasNoMIPs(testCase)
            % all or nothing: the other flags' MIPs are wrong here, so a
            % partial merge would show
            [zarrFns, im] = testCase.writeZarrCopies('uint16');
            flagFns = testCase.writeFlags(im + 1, testCase.BlockSize, testCase.MIPBatchSize, 4);
            t = 1;
            save(flagFns{end}, 't');
            evalc('XR_MIP_zarr(zarrFns{1}, ''batchSize'', testCase.MIPBatchSize, ''parseCluster'', false)');
            evalc('XR_MIP_zarr(zarrFns{2}, ''batchSize'', testCase.MIPBatchSize, ''parseCluster'', false, ''blockFlagFullpaths'', flagFns)');
            testCase.verifySameMIPs(zarrFns{1}, zarrFns{2});
        end

        function xrMIPZarrReadsTheZarrWithoutMemoryForTheMerge(testCase)
            % the flags' MIPs are wrong here, so a merge would show
            [zarrFns, im] = testCase.writeZarrCopies('uint16');
            flagFns = testCase.writeFlags(im + 1, testCase.BlockSize, testCase.MIPBatchSize, 4);
            evalc('XR_MIP_zarr(zarrFns{1}, ''batchSize'', testCase.MIPBatchSize, ''parseCluster'', false)');
            testCase.applyFixture(pktest.FunctionStubFixture('getSystemMemory', 1e-6));
            evalc('XR_MIP_zarr(zarrFns{2}, ''batchSize'', testCase.MIPBatchSize, ''parseCluster'', false, ''blockFlagFullpaths'', flagFns)');
            testCase.verifySameMIPs(zarrFns{1}, zarrFns{2});
        end

        function xrMIPZarrResumesFromExistingSlabs(testCase)
            % slab zarrs from an earlier run are used for missing MIPs, as
            % before (the flags' MIPs are wrong here)
            [zarrFns, im] = testCase.writeZarrCopies('uint16');
            flagFns = testCase.writeFlags(im + 1, testCase.BlockSize, testCase.MIPBatchSize, 4);
            evalc('XR_MIP_zarr(zarrFns{1}, ''batchSize'', testCase.MIPBatchSize, ''parseCluster'', false)');
            mipFn = [fileparts(zarrFns{1}), '/MIPs/im_MIP_x.tif'];
            expected = readtiff(mipFn);
            delete(mipFn);
            evalc('XR_MIP_zarr(zarrFns{1}, ''batchSize'', testCase.MIPBatchSize, ''parseCluster'', false, ''blockFlagFullpaths'', flagFns)');
            testCase.verifyEqual(readtiff(mipFn), expected);
        end
    end

    methods (Access = private)
        function [flagFns] = writeFlags(testCase, im, blockSize, mipBatchSize, taskSize)
            % flags of tasks that wrote im in blocks of blockSize, taskSize
            % blocks each, as the large-file block functions save them; the
            % blocks are in reverse order, so that the tasks' order isn't
            % the order of their blocks in the image
            folder = testCase.createTempFolder();
            sz = size(im, 1 : 3);
            bSubSz = ceil(sz ./ blockSize);
            [y, x, z] = ndgrid(1 : bSubSz(1), 1 : bSubSz(2), 1 : bSubSz(3));
            bboxes = ([y(:), x(:), z(:)] - 1) .* blockSize + 1;
            bboxes = [bboxes, min(bboxes + blockSize - 1, sz)];
            bboxes = flipud(bboxes);
            flagFns = {};
            for s = 1 : taskSize : size(bboxes, 1)
                inds = s : min(s + taskSize - 1, size(bboxes, 1));
                mips = cell(numel(inds), 3);
                for i = 1 : numel(inds)
                    b = bboxes(inds(i), :);
                    mips(i, :) = MIP_zarr_block_mips(im(b(1) : b(4), b(2) : b(5), b(3) : b(6)), b, mipBatchSize, class(im));
                end
                flagFns{end + 1, 1} = sprintf('%sblocks_%d_%d.mat', folder, inds(1), inds(end)); %#ok<AGROW>
                MIP_zarr_save_block_mips(flagFns{end}, bboxes(inds, :), mips, mipBatchSize, 1);
            end
        end

        function [zarrFns, im] = writeZarrCopies(testCase, dtype)
            % the same zarr in two folders
            im = testCase.makeImage(testCase.ImSize, dtype);
            zarrFns = cell(1, 2);
            for i = 1 : 2
                folder = testCase.createTempFolder();
                zarrFns{i} = [folder, 'im.zarr'];
                writezarr(im, zarrFns{i}, 'blockSize', [4, 3, 5]);
            end
        end

        function [] = verifySameMIPs(testCase, zarrFn1, zarrFn2)
            % slab zarrs and MIPs next to two zarr files
            path1 = [fileparts(zarrFn1), '/MIPs/'];
            path2 = [fileparts(zarrFn2), '/MIPs/'];
            for a = {'y', 'x', 'z'}
                slab1 = [path1, 'im_', a{1}, '.zarr'];
                slab2 = [path2, 'im_', a{1}, '.zarr'];
                testCase.assertTrue(exist(slab2, 'dir') > 0, sprintf('No %s slab.', a{1}));
                testCase.verifyEqual(jsondecode(fileread([slab2, '/.zarray'])), jsondecode(fileread([slab1, '/.zarray'])), ...
                    sprintf('%s slab metadata', a{1}));
                testCase.verifySameBits(readzarr(slab2), readzarr(slab1), sprintf('%s slab', a{1}));
                mipFn = ['im_MIP_', a{1}, '.tif'];
                testCase.verifySameBits(readtiff([path2, mipFn]), readtiff([path1, mipFn]), mipFn);
            end
            flagDirs = dir([path2, 'zarr_flag/']);
            testCase.verifyEmpty(setdiff({flagDirs.name}, {'.', '..'}), 'The MIP flag folder was left behind.');
        end

        function [] = verifySameBits(testCase, actual, expected, name)
            % equal values, and for floats equal bits (NaN payloads, -0)
            testCase.verifyClass(actual, class(expected), name);
            testCase.verifySize(actual, size(expected), name);
            if isa(actual, class(expected)) && isequal(size(actual), size(expected))
                if isfloat(expected)
                    testCase.verifyEqual(typecast(actual(:), 'uint8'), typecast(expected(:), 'uint8'), name);
                else
                    testCase.verifyEqual(actual, expected, name);
                end
            end
        end
    end

    methods (Static, Access = private)
        function [im] = makeImage(sz, dtype)
            % for single: negative values, NaN (with a second payload), and
            % lines of only NaN and of only -0 and 0
            if strcmp(dtype, 'single')
                im = single(randn(sz) * 100);
                im(rand(sz) < 0.05) = NaN;
                im(rand(sz) < 0.05) = -0;
                im(rand(sz) < 0.05) = typecast(uint32(hex2dec('7FC00001')), 'single');
                im(2, :, :) = NaN;
                im(:, 3, :) = -0;
                im(:, 3, 2 : 3 : end) = 0;
                im(4, 5, :) = -0;
                im(4, 5, end) = 0;
            else
                im = cast(randi(60000, sz), dtype);
            end
        end

        function [out] = rangeMax(im, dim, starts)
            % max along dim over the ranges that start at starts
            ends = [starts(2 : end) - 1, size(im, dim)];
            idx = {':', ':', ':'};
            parts = cell(numel(starts), 1);
            for k = 1 : numel(starts)
                idx{dim} = starts(k) : ends(k);
                parts{k} = max(im(idx{:}), [], dim);
            end
            out = cat(dim, parts{:});
        end
    end
end
