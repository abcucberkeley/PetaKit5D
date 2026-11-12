classdef FeatherBlendingTest < pktest.TestCase
    % feather_blending_3d, which blends the tiles of a batch of
    % processStitchBlock by their distance maps: the whole batch, or, for a
    % batch that isn't full, the common overlap region of the tiles, written
    % at bbox into the block that holds the tiles. The mex takes uint16 and
    % single tiles with single weights; the MATLAB code blends the others.

    properties (TestParameter)
        % classes of the tiles (and the block) and of the weights that the
        % MATLAB code blends
        classes = struct('uint8Tiles', {{'uint8', 'single'}}, 'doubleTiles', {{'double', 'single'}}, ...
            'uint16TilesDoubleWeights', {{'uint16', 'double'}}, 'singleTilesDoubleWeights', {{'single', 'double'}})
    end

    methods (Test, TestTags = {'Unit'})
        function overlapIsBlendedIntoTheBlock(testCase, classes)
            % the MATLAB code wrote its blend (single) into a copy of a block
            % of another class (indexing4d_mex rejects the mix), which it
            % dropped: the block came back with the tiles' values, so uint8
            % tiles in a batch that isn't full weren't blended
            [tileClass, weightClass] = classes{:};
            if strcmp(weightClass, 'double')
                % the mex rejects the weights as an input error
                testCase.expectMexFallback('feather_blending_3d_with_indexing_mex');
            end
            [f, d, block, bbox] = FeatherBlendingTest.randomBatch(tileClass, weightClass);
            expected = block;
            expected(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6)) = ...
                cast(FeatherBlendingTest.blend(f, d), tileClass);

            [out, mex_compute] = feather_blending_3d(f, d, block, bbox);

            testCase.verifyClass(out, tileClass);
            testCase.verifyEqual(out, expected, 'The overlap is not blended into the block.');
            testCase.verifyFalse(mex_compute);
        end

        function overlapIsBlendedIntoTheBlockWithoutIndexingMex(testCase)
            % without indexing4d_mex (a stub that errors) indexing4d writes
            % into a copy of the block, which must be kept
            testCase.applyFixture(pktest.FunctionStubFixture('indexing4d_mex'));
            testCase.expectMexFallback('indexing4d_mex');
            [f, d, block, bbox] = FeatherBlendingTest.randomBatch('uint8', 'single');
            expected = block;
            expected(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6)) = ...
                cast(FeatherBlendingTest.blend(f, d), 'uint8');

            out = feather_blending_3d(f, d, block, bbox);

            testCase.verifyEqual(out, expected, 'The overlap is not blended into the block.');
        end
    end

    methods (Static, Access = private)
        function [f, d, block, bbox] = randomBatch(tileClass, weightClass)
            % two tiles of the overlap region with values to 200 (tile 2
            % without data in places, tile 1 with data everywhere), their
            % weights, and a block whose region bbox ([y0, x0, z0, y1, x1,
            % z1]) is the overlap
            sz = [6, 7, 5];
            f = cast(randi([1, 200], [sz, 2]), tileClass);
            f(randperm(prod(sz), 40) + prod(sz)) = 0;
            d = cast(rand([sz, 2]) * 0.9 + 0.1, weightClass);
            block = cast(randi([1, 200], [10, 12, 8]), tileClass);
            bbox = [3, 4, 2, 2 + sz(1), 3 + sz(2), 1 + sz(3)];
        end

        function [im] = blend(f, d)
            % weighted mean of the tiles with data
            w = d .* (f ~= 0);
            im = sum(single(f) .* w, 4) ./ sum(w, 4);
        end
    end
end
