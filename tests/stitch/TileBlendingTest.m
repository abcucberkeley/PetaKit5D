classdef TileBlendingTest < pktest.TestCase
    % tile_blending_3d, which blends the tiles of a stitching batch for the
    % blend methods 'none', 'max', 'mean' and 'median' in processStitchBlock:
    % tile_blending_3d_mex in one pass, or the MATLAB code where the mex is
    % missing. A voxel gets the blend of the tiles outside their regions
    % (overlaps other tiles fill), or where that is 0, of all tiles.

    properties (TestParameter)
        blendMethod = {'none', 'max', 'mean', 'median'}
        dtype = {'uint8', 'uint16', 'single'}
    end

    methods (Test, TestTags = {'Unit'})
        function mexMatchesTheMatlabCode(testCase, blendMethod, dtype)
            % tiles covering parts of the batch with zeros, regions (also
            % overlapping and empty ones), for single negative values, NaN,
            % Inf and -0; 3D, 2D and a singleton dimension
            testCase.assumeMexAvailable('tile_blending_3d_mex');
            sizes = {[23, 19, 11], [17, 13], [1, 9, 6], [23, 19, 11]};
            nTiles = [2, 3, 5, 9];
            tims = cell(numel(sizes), 1);
            mboxes = cell(numel(sizes), 1);
            outs = cell(numel(sizes), 1);
            for i = 1 : numel(sizes)
                [tims{i}, mboxes{i}] = testCase.randomBatch(sizes{i}, nTiles(i), dtype);
                outs{i} = tile_blending_3d(tims{i}, mboxes{i}, blendMethod);
            end
            % the MATLAB code: the mex hidden by a stub that errors
            testCase.applyFixture(pktest.FunctionStubFixture('tile_blending_3d_mex'));
            for i = 1 : numel(sizes)
                evalc('expected = tile_blending_3d(tims{i}, mboxes{i}, blendMethod);');
                testCase.verifySameResult(outs{i}, expected, dtype, blendMethod, ...
                    sprintf('size %s, %d tiles', mat2str(sizes{i}), nTiles(i)));
            end
        end

        function voxelsBlendTheTilesOutsideTheirRegions(testCase, blendMethod)
            % 3 tiles of 10, 21 and 40 with zeros, tile 2 counts as 0 in
            % columns 1-2 and tile 3 in column 1, rows 1-2
            tim = zeros(4, 4, 2, 3, 'uint16');
            tim(:, :, :, 1) = 10;
            tim(:, :, :, 2) = 21;
            tim(:, :, :, 3) = 40;
            tim(1, :, :, 1) = 0;
            tim(2 : 3, 3, 1, 2) = 0;
            tim(4, 4, 2, :) = 0;
            mboxes = [2, 1, 1, 1, 4, 2, 2; 3, 1, 1, 1, 2, 1, 2];
            masked = false(size(tim));
            masked(:, 1 : 2, :, 2) = true;
            masked(1 : 2, 1, :, 3) = true;

            out = tile_blending_3d(tim, mboxes, blendMethod);

            expected = zeros(4, 4, 2, 'single');
            for i = 1 : numel(expected)
                [y, x, z] = ind2sub([4, 4, 2], i);
                v = single(squeeze(tim(y, x, z, :)));
                sel = v ~= 0 & ~squeeze(masked(y, x, z, :));
                if ~any(sel)
                    sel = v ~= 0;
                end
                if any(sel)
                    switch blendMethod
                        case 'none'
                            expected(i) = v(find(sel, 1));
                        case 'max'
                            expected(i) = max(v(sel));
                        case 'mean'
                            expected(i) = mean(v(sel));
                        case 'median'
                            expected(i) = median(v(sel));
                    end
                end
            end
            testCase.verifyEqual(single(out), expected);
            % the mex keeps the class of the tiles for 'none' and 'max'
            expectedClass = 'single';
            if any(strcmp(blendMethod, {'none', 'max'})) && exist('tile_blending_3d_mex', 'file') == 3
                expectedClass = 'uint16';
            end
            testCase.verifyClass(out, expectedClass);
        end

        function mexRejectsInvalidInput(testCase)
            testCase.assumeMexAvailable('tile_blending_3d_mex');
            tim = zeros(4, 4, 4, 2, 'uint16');
            testCase.verifyError(@() tile_blending_3d_mex(tim, [3, 1, 1, 1, 2, 2, 2], 'none'), 'tile_blending:inputError');
            testCase.verifyError(@() tile_blending_3d_mex(tim, [1, 0, 1, 1, 2, 2, 2], 'none'), 'tile_blending:inputError');
            testCase.verifyError(@() tile_blending_3d_mex(tim, [1, 1, 1, 1, 5, 2, 2], 'none'), 'tile_blending:inputError');
            testCase.verifyError(@() tile_blending_3d_mex(tim, [1, 1, 1, 1, 2.5, 2, 2], 'none'), 'tile_blending:inputError');
            testCase.verifyError(@() tile_blending_3d_mex(tim, [1, 1, 1, 1, 2, 2], 'none'), 'tile_blending:inputError');
            testCase.verifyError(@() tile_blending_3d_mex(tim, [], 'feather'), 'tile_blending:inputError');
            testCase.verifyError(@() tile_blending_3d_mex(zeros(4, 4, 4, 0, 'uint16'), [], 'none'), 'tile_blending:inputError');
            testCase.verifyError(@() tile_blending_3d_mex(double(tim), [], 'none'), 'tile_blending:dataTypeError');
            testCase.verifyError(@() tile_blending_3d(tim, [], 'None'), ?MException);
        end

        function otherClassesUseTheMatlabCode(testCase)
            % double tiles are blended in single by the MATLAB code
            [tim, mboxes] = testCase.randomBatch([12, 10, 7], 3, 'uint16');
            for method = {'none', 'max', 'mean', 'median'}
                log = evalc('out = tile_blending_3d(double(tim), mboxes, method{1});');
                testCase.verifyEmpty(log, 'The mex was called for double tiles.');
                testCase.verifySameResult(out, tile_blending_3d(single(tim), mboxes, method{1}), 'single', method{1}, method{1});
            end
        end
    end

    methods (Access = private)
        function [] = verifySameResult(testCase, actual, expected, dtype, blendMethod, name)
            % integer 'none' and 'max' keep the class of the tiles (the values
            % of the single result); otherwise the bits of the single result
            if any(strcmp(blendMethod, {'none', 'max'})) && ~strcmp(dtype, 'single')
                testCase.verifyClass(actual, dtype, name);
                actual = single(actual);
            end
            testCase.verifyClass(actual, 'single', name);
            testCase.verifySize(actual, size(expected), name);
            if isa(actual, 'single') && isequal(size(actual), size(expected))
                testCase.verifyEqual(typecast(actual(:), 'uint32'), typecast(expected(:), 'uint32'), name);
            end
        end
    end

    methods (Static, Access = private)
        function [tim, mboxes] = randomBatch(sz, nT, dtype)
            % tiles placed at random parts of the batch (zeros elsewhere),
            % 0 to 3 regions per tile and an empty one
            sz = [sz, ones(1, 3 - numel(sz))];
            tim = zeros([sz, nT], dtype);
            mboxes = zeros(0, 7);
            for j = 1 : nT
                s = arrayfun(@(n) randi(n), sz);
                e = arrayfun(@(a, n) randi([a, n]), s, sz);
                if j == 1
                    s = [1, 1, 1];
                    e = sz;
                end
                tsz = e - s + 1;
                switch dtype
                    case 'uint8'
                        v = randi([0, 255], tsz);
                    case 'uint16'
                        v = randi([0, 65535], tsz);
                    case 'single'
                        v = randn(tsz) .* 10 .^ randi([-2, 4], tsz);
                end
                v(rand(tsz) < 0.3) = 0;
                v = cast(v, dtype);
                if strcmp(dtype, 'single')
                    r = rand(tsz);
                    v(r < 0.05) = NaN;
                    v(r >= 0.05 & r < 0.07) = Inf;
                    v(r >= 0.07 & r < 0.09) = -Inf;
                    v(r >= 0.09 & r < 0.12) = -0;
                end
                tim(s(1) : e(1), s(2) : e(2), s(3) : e(3), j) = v;
                for k = 1 : randi([0, 3])
                    ms = arrayfun(@(a, b) randi([a, b]), s, e);
                    me = arrayfun(@(a, b) randi([a, b]), ms, e);
                    mboxes(end + 1, :) = [j, ms, me];
                end
            end
            mboxes(end + 1, :) = [nT, 2, 1, 1, 1, sz(2), sz(3)];
        end
    end
end
