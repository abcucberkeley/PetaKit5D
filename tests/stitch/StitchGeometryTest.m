classdef StitchGeometryTest < pktest.TestCase
    % geometry helpers of the stitching: the axis order check of
    % XR_stitching_frame_zarr_dev_v1

    methods (Test, TestTags = {'Unit'})
        function axisOrderMustBeAPermutationOfXYZ(testCase)
            % the check never raised (~regexpi(...) is false for no match
            % too), so the stitching sorted any order into some permutation;
            % 'z,x,y' was missing from it. It comes before the tiles are
            % read, so the tiles here don't exist: valid orders end in
            % another error later.
            dataPath = testCase.createTempFolder();
            for ao = {'x,x,z', 'a,b,c', 'x,y', 'x,y,z,y', 'x,y;z'}
                testCase.verifyError(@() testCase.frameWithAxisOrder(dataPath, ao{1}), ...
                    'XR_stitching_frame_zarr_dev_v1:axisOrder', ao{1});
            end
            for ao = {'x,y,z', 'y,x,z', 'z,y,x', 'x,z,y', 'z,x,y', 'y,z,x', '-x,y,z', 'x,-y,z', '-z,x,-y', 'yxz', ' y, x, z'}
                id = '';
                try
                    testCase.frameWithAxisOrder(dataPath, ao{1});
                catch ME
                    id = ME.identifier;
                end
                testCase.verifyNotEqual(id, 'XR_stitching_frame_zarr_dev_v1:axisOrder', ao{1});
            end
        end
    end

    methods (Access = private)
        function [] = frameWithAxisOrder(~, dataPath, axisOrder)
            fns = {[dataPath, 't1.tif'], [dataPath, 't2.tif']};
            evalc(['XR_stitching_frame_zarr_dev_v1(fns, [0, 0, 0; 1, 0, 0], ''axisOrder'', axisOrder, ', ...
                '''parseCluster'', false, ''ResultPath'', [dataPath, ''res/out.zarr''])']);
        end
    end
end
