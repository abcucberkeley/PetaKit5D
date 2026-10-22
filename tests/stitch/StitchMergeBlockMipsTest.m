classdef StitchMergeBlockMipsTest < pktest.TestCase
    % stitch_merge_block_mips, which merges the z MIPs the blending tasks
    % save in their flag files into the z MIP of the stitched image

    methods (Test, TestTags = {'Unit'})
        function mergedMIPEqualsMIPOfTheImage(testCase)
            im = uint16(randi(60000, [8, 10, 6]));
            flagFns = testCase.writeFlags(im, [4, 5, 3], 3);
            testCase.verifyEqual(stitch_merge_block_mips(flagFns, size(im), 'uint16'), max(im, [], 3));
        end

        function singleKeepsNegativesAndNaN(testCase)
            % max ignores NaN, so a column of NaN stays NaN
            im = single(randn([8, 10, 6]) - 5);
            im(2, 3, :) = NaN;
            im(5, 7, 2) = NaN;
            flagFns = testCase.writeFlags(im, [4, 5, 3], 2);
            testCase.verifyEqual(stitch_merge_block_mips(flagFns, size(im), 'single'), max(im, [], 3));
        end

        function flagWithoutMIPsGivesNoMIP(testCase)
            % e.g., flags of a run whose tasks didn't save MIPs: the caller
            % then reads the image
            im = uint16(randi(60000, [8, 10, 6]));
            flagFns = testCase.writeFlags(im, [4, 5, 3], 3);
            fclose(fopen(flagFns{2}, 'w'));
            testCase.verifyEmpty(stitch_merge_block_mips(flagFns, size(im), 'uint16'));
        end
    end

    methods (Access = private)
        function [flagFns] = writeFlags(testCase, im, batchSize, taskSize)
            % flags as processStitchBlock writes them, with taskSize batches
            % per task
            folder = testCase.createTempFolder();
            bSubSz = ceil(size(im) ./ batchSize);
            numBatch = prod(bSubSz);
            flagFns = {};
            for s = 1 : taskSize : numBatch
                batchInds = s : min(s + taskSize - 1, numBatch);
                mips = cell(numel(batchInds), 1);
                mipBboxes = zeros(numel(batchInds), 4);
                for i = 1 : numel(batchInds)
                    [y, x, z] = ind2sub(bSubSz, batchInds(i));
                    obStart = ([y, x, z] - 1) .* batchSize + 1;
                    obEnd = min(obStart + batchSize - 1, size(im));
                    mips{i} = max(im(obStart(1) : obEnd(1), obStart(2) : obEnd(2), obStart(3) : obEnd(3)), [], 3);
                    mipBboxes(i, :) = [obStart(1 : 2), obEnd(1 : 2)];
                end
                flagFns{end + 1} = sprintf('%sblocks_%d_%d.mat', folder, batchInds(1), batchInds(end));
                save(flagFns{end}, 'batchInds', 'mipBboxes', 'mips');
            end
        end
    end
end
