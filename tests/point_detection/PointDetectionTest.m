classdef PointDetectionTest < pktest.TestCase
    % 3D point source detection (demos/demo_point_detection_3d_*.m) on
    % synthetic Gaussian spots with known positions, amplitude and
    % background. The detection mex functions are not built for ARM macOS,
    % so these tests are skipped there (see MexInventoryTest).

    properties (Constant)
        Spots = [16.3, 20.7, 18.2; 16, 48, 32; 48, 24, 40; 40, 40, 20; 32, 32, 48]
        Sigma = [1.5, 1.5]
        Amplitude = 200
        Background = 100
    end

    methods (TestMethodSetup)
        function requireDetectionMex(testCase)
            testCase.assumeMexAvailable('compute_gradient_hessian_mex', 'adtest_mex', 'KDTreeBallQuery');
            % the t-statistic threshold is cached in a persistent variable
            % keyed only on alpha
            clear XR_pointSourceDetection3D
        end
    end

    methods (Test, TestTags = {'Unit'})
        function detectsSpotsAtTheirPositions(testCase)
            vol = testCase.spotVolume();
            [ps, mask] = XR_pointSourceDetection3D(vol, testCase.Sigma, 'BackgroundCorrection', false);

            testCase.verifyClass(mask, 'logical');
            testCase.verifySize(mask, size(vol));
            testCase.verifyNumElements(ps.x, size(testCase.Spots, 1), 'Wrong number of detections.');
            testCase.verifyEqual(testCase.matchDetections([ps.y(:), ps.x(:), ps.z(:)]), testCase.Spots, 'AbsTol', 0.25);
            testCase.verifyEqual(ps.A(:), testCase.Amplitude * ones(numel(ps.A), 1), 'RelTol', 0.1);
            testCase.verifyEqual(ps.c(:), testCase.Background * ones(numel(ps.c), 1), 'AbsTol', 2);
            testCase.verifyClass(ps.isPSF, 'logical');
        end

        function noSpotsInPureNoise(testCase)
            vol = testCase.Background + 5 * randn(64, 64, 64);
            ps = XR_pointSourceDetection3D(vol, testCase.Sigma, 'BackgroundCorrection', false);
            testCase.verifyTrue(isempty(ps) || isempty(ps.x), 'Detected spots in pure noise.');
        end
    end

    methods (Test, TestTags = {'Integration'})
        function runDetectionWritesDetectionsAndMask(testCase)
            % getShortPath (used for messages) only splits paths on filesep,
            % so on Windows it errors for paths with '/' separators, which
            % these tests (and PetaKit5D's other functions) use
            testCase.assumeFalse(ispc, 'getShortPath does not handle ''/'' separators on Windows.');
            [data, dataPath] = testCase.movieData(1);

            XR_runDetection3D(data, 'Sigma', testCase.Sigma, 'DetectionMethod', 'speedup', 'WriteAmira', false, ...
                'Overwrite', true, 'ResultsPath', {[dataPath, 'Analysis/']});

            a = load([dataPath, 'Analysis/Detection3D.mat'], 'frameInfo');
            fi = a.frameInfo(1);
            testCase.verifyNumElements(fi.x, size(testCase.Spots, 1));
            testCase.verifyEqual(testCase.matchDetections([fi.y(:), fi.x(:), fi.z(:)]), testCase.Spots, 'AbsTol', 0.3);
            % coordinates for the tracker, z scaled by zAniso
            testCase.verifyEqual(fi.zCoord(:, 1), fi.z(:) * data.zAniso, 'AbsTol', 1e-10);
            mask = readtiff([dataPath, 'Analysis/Masks/dmask_1.tif']);
            testCase.verifySize(mask, [64, 64, 64]);
            testCase.verifyTrue(any(mask(:)));
        end

        function trackingLinksSpotsAcrossFrames(testCase)
            % detection + tracking (demo step 5.1) on a 3-frame movie with
            % spots moving by a fraction of a voxel per frame
            testCase.assumeFalse(ispc, 'getShortPath does not handle ''/'' separators on Windows.');
            velocity = [0.6, 1.0, 0.4];
            [data, dataPath] = testCase.movieData(3, velocity);
            resultsPath = {[dataPath, 'Analysis/']};

            XR_runDetection3D(data, 'Sigma', testCase.Sigma, 'DetectionMethod', 'speedup', 'WriteAmira', false, ...
                'Overwrite', true, 'ResultsPath', resultsPath);
            XR_runTracking3D(data, loadTrackSettings('Radius', [3, 6], 'MaxGapLength', 2), ...
                'FileName', 'trackedFeatures.mat', 'Overwrite', true, 'ResultsPath', resultsPath);

            a = load([resultsPath{1}, 'trackedFeatures.mat'], 'tracksFinal');
            tracks = a.tracksFinal;
            testCase.verifyNumElements(tracks, size(testCase.Spots, 1), 'Each spot should form one track.');
            for i = 1 : numel(tracks)
                % 8 values per frame (x, y, z, A and their std)
                coords = reshape(tracks(i).tracksCoordAmpCG, 8, [])';
                testCase.verifySize(coords, [3, 8], 'Each track should span all 3 frames.');
                % displacement per frame in [y, x, z]; the tracker works on
                % z scaled by zAniso (all axes in xy pixel units)
                step = mean(diff(coords(:, [2, 1, 3]), 1, 1), 1);
                testCase.verifyEqual(step, velocity .* [1, 1, data.zAniso], 'AbsTol', 0.1);
            end
        end

        function estimateSigmaFromBead(testCase)
            % GU_estimateSigma3D (demo step 2) on a single bead
            dataPath = testCase.createTempFolder();
            vol = pktest.beads([48, 48, 48], [24, 24, 24], [1.6, 2.4], 1000, 100) + 5 * randn(48, 48, 48, 'single');
            writetiff(uint16(vol), [dataPath, 'psf.tif']);

            [sigmaXY, sigmaZ, XYZ] = GU_estimateSigma3D([dataPath, 'psf.tif'], []);

            testCase.verifyEqual(sigmaXY, 1.6, 'AbsTol', 0.1);
            testCase.verifyEqual(sigmaZ, 2.4, 'AbsTol', 0.15);
            testCase.verifyEqual([XYZ.y, XYZ.x, XYZ.z], [24, 24, 24], 'AbsTol', 0.2);
        end
    end

    methods (Test, TestTags = {'Integration', 'KnownIssue'})
        function lowSNRDetectionRuns(testCase)
            % KNOWN ISSUE: DetectionMethod 'lowSNR' (used in both detection
            % demos) calls XR_estimate_overall_background_information,
            % which calls crop_image_by_bounding_box, a function that does
            % not exist in the repository.
            testCase.assumeFalse(ispc, 'getShortPath does not handle ''/'' separators on Windows.');
            [data, dataPath] = testCase.movieData(1);

            XR_runDetection3D(data, 'Sigma', testCase.Sigma, 'DetectionMethod', 'lowSNR', 'WriteAmira', false, ...
                'Overwrite', true, 'ResultsPath', {[dataPath, 'Analysis/']}, 'BackgroundCorrection', true);

            a = load([dataPath, 'Analysis/Detection3D.mat'], 'frameInfo');
            testCase.verifyNumElements(a.frameInfo(1).x, size(testCase.Spots, 1));
        end
    end

    methods (Access = private)
        function [vol] = spotVolume(testCase, shift)
            if nargin < 2
                shift = [0, 0, 0];
            end
            vol = double(pktest.beads([64, 64, 64], testCase.Spots + shift, testCase.Sigma, testCase.Amplitude, ...
                testCase.Background)) + 5 * randn(64, 64, 64);
        end

        function [data, dataPath] = movieData(testCase, nFrames, velocity)
            % minimal movie structure as returned by XR_loadConditionData3D;
            % spots move by velocity [y, x, z] voxels per frame
            if nargin < 3
                velocity = [0, 0, 0];
            end
            dataPath = testCase.createTempFolder();
            chPath = [dataPath, 'ch1/'];
            mkdir(chPath);
            frameFns = cell(1, nFrames);
            for k = 1 : nFrames
                frameFns{k} = sprintf('%sframe_%04d.tif', chPath, k);
                writetiff(uint16(testCase.spotVolume((k - 1) * velocity)), frameFns{k});
            end
            data = struct('source', chPath, 'channels', {{chPath}}, 'movieLength', nFrames, ...
                'framePathsDS', {{frameFns}}, 'zAniso', 1.5);
        end

        function [matched] = matchDetections(testCase, detected)
            % detected [y, x, z] reordered to the nearest true spot
            matched = nan(size(testCase.Spots));
            for i = 1 : size(testCase.Spots, 1)
                [~, j] = min(sum((detected - testCase.Spots(i, :)) .^ 2, 2));
                matched(i, :) = detected(j, :);
            end
        end
    end
end
