classdef AxisOrderMappingTest < pktest.TestCase
    % axis_order_mapping(axisOrder, 'yxz') is the permutation that stores
    % yxz data in axisOrder (dimension i holds axis axisOrder(i)), which is
    % how XR_deskewRotateFrame applies outputAxisOrder. Bringing data
    % stored in axisOrder back to yxz takes the inverse (ipermute), see the
    % known issue XRDeskewRotateFrameTest/inputAxisOrderZYXIsPermutedToYXZ.

    properties (TestParameter)
        axisOrder = {'yxz', 'xyz', 'zxy', 'yzx', 'zyx', 'xzy'}
    end

    methods (Test, TestTags = {'Unit'})
        function storesYXZDataInTheGivenOrder(testCase, axisOrder)
            vol = reshape(1 : 2 * 3 * 4, [2, 3, 4]);
            stored = permute(vol, axis_order_mapping(axisOrder, 'yxz'));
            % the axes have different lengths, so the size pins the order
            testCase.verifySize(stored, arrayfun(@(c) size(vol, strfind('yxz', c)), axisOrder), ...
                sprintf('axis_order_mapping(''%s'', ''yxz'') does not store yxz data in that order.', axisOrder));
        end
    end
end
