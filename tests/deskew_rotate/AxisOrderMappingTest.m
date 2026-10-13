classdef AxisOrderMappingTest < pktest.TestCase
    % axis_order_mapping(inputAxisOrder, 'yxz') must give the permutation
    % that brings data stored in inputAxisOrder back to yxz.

    properties (TestParameter)
        % self-inverse permutations of yxz
        swapOrder = {'yxz', 'xyz', 'zxy', 'yzx'}
        % 3-cycles
        cycleOrder = {'zyx', 'xzy'}
    end

    methods (Test, TestTags = {'Unit'})
        function permutesBackToYXZ(testCase, swapOrder)
            testCase.verifyRoundTrip(swapOrder);
        end
    end

    methods (Test, TestTags = {'Unit', 'KnownIssue'})
        function cyclicOrdersPermuteBackToYXZ(testCase, cycleOrder)
            % KNOWN ISSUE: utils/axis_order_mapping.m computes
            % order_mat(i) = strfind(outputAxisOrder, inputAxisOrder(i)),
            % which is the inverse permutation; it only works when the
            % permutation is its own inverse. Fix: swap the two arguments
            % of strfind.
            testCase.verifyRoundTrip(cycleOrder);
        end
    end

    methods (Access = private)
        function [] = verifyRoundTrip(testCase, order)
            vol = reshape(1 : 2 * 3 * 4, [2, 3, 4]);
            % data stored in the given order: dimension i holds axis order(i)
            stored = permute(vol, arrayfun(@(c) strfind('yxz', c), order));
            testCase.verifyEqual(permute(stored, axis_order_mapping(order, 'yxz')), vol, ...
                sprintf('axis_order_mapping(''%s'', ''yxz'') does not restore yxz.', order));
        end
    end
end
