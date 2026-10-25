classdef DeconBackProjectorTest < pktest.TestCase
    % omw_backprojector_generation with omw_backprojector_mex (sorted OTF
    % values, connected components, convex hulls of the OTF mask slices,
    % bounds of the window) gives the same outputs, bit for bit, as the code
    % before it (private/omw_backprojector_generation_before), psf_gen_new's
    % median filter in the mex the same PSFs as medfilt3, and the mex modes
    % match the MATLAB functions they replace.

    properties (TestParameter)
        % PSF, alpha, skewed, parameters: Gaussians and lattice-like PSFs
        % (1, 3 or 5 OTF lobes, sheared or not), skewed or not or detected,
        % and several apodization ranges and thresholds
        psfCase = struct( ...
            'gaussian', {{'gaussian', 0.005, false, {}}}, ...
            'gaussianSkewed', {{'gaussian', 0.005, true, {'OTFCumThresh', 0.8}}}, ...
            'elongatedWideTaper', {{'elongated', 0.02, true, {'OTFCumThresh', 0.65, 'hannWinBounds', [0.8, 1.2]}}}, ...
            'elongatedNarrowTaper', {{'elongated', 0.01, false, {'hannWinBounds', [0.95, 1]}}}, ...
            'threeLobes', {{'lobes3', 0.005, true, {}}}, ...
            'threeLobesShearedDetected', {{'lobes3sheared', 0.005, [], {'OTFCumThresh', 0.8}}}, ...
            'fiveLobes', {{'lobes5', 0.005, true, {'OTFCumThresh', 0.85}}}, ...
            'degenerateBounds', {{'gaussian', 0.005, false, {'hannWinBounds', [0.8, 0.8]}}}, ...
            'boundsBeyondTheEdge', {{'elongated', 0.005, false, {'hannWinBounds', [0.2, 2]}}})
    end

    methods (Test, TestTags = {'Unit'})
        function backProjectorMatchesThePreviousCode(testCase, psfCase)
            testCase.assumeMexAvailable('omw_backprojector_mex');
            psf = testCase.makePSF(psfCase{1});
            [expected, actual] = testCase.bothVersions(psf, psfCase{2 : end});
            testCase.verifyOutputsEqual(actual, expected);
        end

        function areaThresholdOfAComponentSizeMatches(testCase)
            % a component of exactly OTFAreaThresh voxels is kept, as
            % bwareaopen keeps it (here the smallest OTF lobe, whose loss
            % changes the labels of the skewed code)
            testCase.assumeMexAvailable('omw_backprojector_mex');
            psf = testCase.makePSF('lobes3');
            psfn = single(psf) ./ sum(single(psf), 'all', 'double');
            a = abs(decon_psf2otf(psfn));
            v = sort(double(a(:)), 'descend');
            t = max(v(find(cumsum(v) > sum(v) * 0.9, 1, 'first')) / a(1), 1e-3);
            M1 = imopen(fftshift(a) > a(1) * t, strel('sphere', 2));
            M2 = imopen(imclose(M1, strel('sphere', 2)), strel('sphere', 2));
            sizes2 = cellfun(@numel, getfield(bwconncomp(M2, 26), 'PixelIdxList'));
            sizes1 = cellfun(@numel, getfield(bwconncomp(M1, 26), 'PixelIdxList'));
            % the threshold only acts in the second area opening
            p = min(sizes2);
            testCase.assumeEqual(numel(sizes2), 3);
            testCase.assumeGreaterThanOrEqual(min(sizes1), p);
            [expected, actual] = testCase.bothVersions(psf, 0.005, true, {'OTFAreaThresh', p});
            testCase.verifyOutputsEqual(actual, expected);
        end

        function backProjectorMatchesWithOneThread(testCase)
            % the slices and components are processed in parallel
            testCase.assumeMexAvailable('omw_backprojector_mex');
            psf = testCase.makePSF('lobes3');
            nt = maxNumCompThreads(1);
            testCase.addTeardown(@maxNumCompThreads, nt);
            [expected, actual] = testCase.bothVersions(psf, 0.005, true, {});
            testCase.verifyOutputsEqual(actual, expected);
        end

        function withoutTheMexTheMatlabCodeRuns(testCase)
            psf = testCase.makePSF('lobes3sheared');
            expected = cell(1, 4);
            actual = cell(1, 4);
            [~, expected{1 : 4}] = evalc('omw_backprojector_generation_before(psf, 0.005, [], ''OTFCumThresh'', 0.8)');
            % a function of that name makes exist() return 2, not 3
            testCase.applyFixture(pktest.FunctionStubFixture('omw_backprojector_mex', []));
            [~, actual{1 : 4}] = evalc('omw_backprojector_generation(psf, 0.005, [], ''OTFCumThresh'', 0.8)');
            testCase.verifyOutputsEqual(actual, expected);
        end

        function theReplacedFunctionsAreNotCalled(testCase)
            % the mex replaces bwconvhull, bwareaopen and bwconncomp (stubs
            % that would change the result if called)
            testCase.assumeMexAvailable('omw_backprojector_mex');
            psf = testCase.makePSF('lobes5');
            expected = cell(1, 4);
            actual = cell(1, 4);
            [~, expected{1 : 4}] = evalc('omw_backprojector_generation_before(psf, 0.005, true, ''OTFCumThresh'', 0.85)');
            testCase.applyFixture(pktest.FunctionStubFixture('bwconvhull', false));
            testCase.applyFixture(pktest.FunctionStubFixture('bwareaopen', false));
            testCase.applyFixture(pktest.FunctionStubFixture('bwconncomp', struct('NumObjects', 0, 'PixelIdxList', {{}})));
            [~, actual{1 : 4}] = evalc('omw_backprojector_generation(psf, 0.005, true, ''OTFCumThresh'', 0.85)');
            testCase.verifyOutputsEqual(actual, expected);
        end

        function convexHullsMatchBwconvhull(testCase)
            testCase.assumeMexAvailable('omw_backprojector_mex');
            for trial = 1 : 60
                sz = [randi([2, 40]), randi([2, 40]), randi([2, 30])];
                switch mod(trial, 4)
                    case 0
                        B = rand(sz) < 0.02 + 0.3 * rand;
                    case 1
                        B = imgaussfilt3(rand(sz), 0.5 + 2 * rand) > 0.5 + 0.03 * randn;
                    case 2
                        B = false(sz);
                        B(randi(prod(sz), randi(20), 1)) = true;
                    case 3
                        B = rand(sz) < 0.6;
                end
                zinds = find(sum(B, [1, 2]))';
                xinds = find(sum(B, [1, 3]));
                yinds = find(sum(B, [2, 3]))';
                expected = B;
                for z = zinds
                    expected(:, :, z) = bwconvhull(expected(:, :, z), "union");
                end
                for y = yinds
                    expected(y, :, :) = bwconvhull(squeeze(expected(y, :, :)), "union");
                end
                for x = xinds
                    expected(:, x, :) = bwconvhull(squeeze(expected(:, x, :)), "union");
                end
                for nt = [1, 4]
                    actual = omw_backprojector_mex('convhull', B, zinds, yinds, xinds, nt);
                    testCase.verifyEqual(actual, expected, sprintf('trial %d, %d threads', trial, nt));
                end
            end
        end

        function sliceHullsMatchBwconvhull(testCase)
            % single slices of many shapes, including lines, points and
            % sloped edges, where the polygon filling rules matter
            testCase.assumeMexAvailable('omw_backprojector_mex');
            for trial = 1 : 600
                S = false(randi([2, 50]), randi([2, 50]));
                switch mod(trial, 6)
                    case 0
                        S = rand(size(S)) < rand;
                    case 1
                        S = imgaussfilt(rand(size(S)), 1 + 3 * rand) > 0.5 + 0.05 * randn;
                    case 2
                        S(randi(numel(S), randi(6), 1)) = true;
                    case 3
                        r0 = randi(size(S, 1));
                        c0 = randi(size(S, 2));
                        S(r0 : min(end, r0 + randi(12)), c0 : min(end, c0 + randi(12))) = true;
                    case 4
                        S(randi(size(S, 1)), :) = rand(1, size(S, 2)) < 0.5;
                    case 5
                        S(:, randi(size(S, 2))) = rand(size(S, 1), 1) < 0.5;
                end
                B = false([size(S), 2]);
                B(:, :, 1) = S;
                actual = omw_backprojector_mex('convhull', B, 1, [], [], 2);
                testCase.verifyEqual(actual(:, :, 1), bwconvhull(S, 'union'), sprintf('trial %d', trial));
            end
        end

        function componentsMatchBwareaopenAndBwconncomp(testCase)
            testCase.assumeMexAvailable('omw_backprojector_mex');
            for trial = 1 : 80
                sz = [randi([2, 30]), randi([2, 30]), randi([2, 25])];
                if mod(trial, 2) == 0
                    B = rand(sz) < 0.05 + 0.4 * rand;
                else
                    B = imgaussfilt3(rand(sz), 0.5 + 2 * rand) > 0.5;
                end
                p = randi([0, 40]);
                CC = bwconncomp(B, 26);
                vols = cellfun(@numel, CC.PixelIdxList);
                for nt = [1, 3]
                    testCase.verifyEqual(omw_backprojector_mex('components', B, 'minsize', p, nt), bwareaopen(B, p));
                    for k = [1, 3]
                        expected = false(sz);
                        if ~isempty(vols)
                            if k == 1
                                [~, mi] = max(vols);
                            else
                                [~, mi] = maxk(vols, k);
                            end
                            expected(cat(1, CC.PixelIdxList{mi})) = true;
                        end
                        testCase.verifyEqual(omw_backprojector_mex('components', B, 'largest', k, nt), expected);
                    end
                end
            end
        end

        function largestComponentTiesGoToTheFirst(testCase)
            % components of equal size: max and maxk pick the first ones in
            % bwconncomp's order
            testCase.assumeMexAvailable('omw_backprojector_mex');
            B = false(20, 20, 20);
            B(2 : 3, 2 : 3, 2) = true;
            B(10 : 11, 2 : 3, 2) = true;
            B(2 : 3, 10 : 11, 15) = true;
            B(15 : 16, 15 : 16, 8) = true;
            CC = bwconncomp(B, 26);
            for k = [1, 3]
                [~, mi] = maxk(cellfun(@numel, CC.PixelIdxList), k);
                expected = false(size(B));
                expected(cat(1, CC.PixelIdxList{mi})) = true;
                testCase.verifyEqual(omw_backprojector_mex('components', B, 'largest', k, 2), expected);
            end
        end

        function sortMatchesSort(testCase)
            testCase.assumeMexAvailable('omw_backprojector_mex');
            for trial = 1 : 20
                a = abs(randn(randi([0, 2e5]), 1, 'single')) .^ (1 + 3 * rand);
                a(rand(size(a)) < 0.1) = 0;
                if mod(trial, 3) == 0
                    a = round(a * 10) / 10;
                end
                testCase.verifyEqual(omw_backprojector_mex('sortdesc', a, 4), sort(double(a(:)), 'descend'));
            end
            % NaN and negative values are left to sort
            testCase.verifyEmpty(omw_backprojector_mex('sortdesc', single([1, NaN, 2]), 2));
            testCase.verifyEmpty(omw_backprojector_mex('sortdesc', single([1, -1]), 2));
        end

        function windowIsDecidedOnlyWhereTheBoundsHold(testCase)
            % the voxels the mex decides have the window value the MATLAB
            % code computes from the nearest edge points
            testCase.assumeMexAvailable('omw_backprojector_mex');
            psz = [40, 36, 30];
            c = (psz + 1) / 2;
            [Y, X, Z] = ndgrid((1 : psz(1)) - c(1), (1 : psz(2)) - c(2), (1 : psz(3)) - c(3));
            % an ellipsoid with a bumpy surface as the OTF mask
            [th, ph, rs] = cart2sph(X, Y, Z);
            mask = rs < 9 ./ sqrt(cos(ph) .^ 2 .* (cos(th) .^ 2 / 1.2 + sin(th) .^ 2 / 0.8) + sin(ph) .^ 2 / 2) .* (1 + 0.1 * cos(5 * th));
            edge = mask & ~imerode(mask, strel('sphere', 1));
            [Xn, Yn, Zn] = sph2cart(th, ph, ones(psz));
            E = [Yn(edge), Xn(edge), Zn(edge)];
            Rn = rs(edge);
            [inds, d] = knnsearch(E, [Yn(:), Xn(:), Zn(:)], 'K', 2);
            % (with a zero distance the code doesn't use the bounds)
            testCase.assumeGreaterThan(min(d(1 : size(E, 1), 2)), 0);
            d1 = d(inds(:, 1), 1);
            d2 = d(inds(:, 2), 2);
            rn = reshape((Rn(inds(:, 1)) .* d2 + Rn(inds(:, 2)) .* d1) ./ (d1 + d2), psz);
            bw = rs ./ rn;
            for bounds = {[0.8, 1], [0.8, 1.2], [0.5, 0.6]}
                l = bounds{1}(1);
                u = bounds{1}(2);
                [w, uInds] = omw_backprojector_mex('window', psz, E, Rn, bounds{1}, 3);
                decided = true(psz);
                decided(uInds) = false;
                testCase.verifyEqual(w(decided), double(bw(decided) < l));
                testCase.verifyTrue(all(bw(decided) < l | bw(decided) >= u));
                testCase.verifyEqual(w(uInds), zeros(numel(uInds), 1));
                % most voxels are decided
                testCase.verifyLessThan(numel(uInds), 0.5 * prod(psz));
            end
        end

        function medianFilterMatchesMedfilt3(testCase)
            testCase.assumeMexAvailable('omw_backprojector_mex');
            shapes = {[1, 1, 1], [1, 5, 7], [5, 1, 7], [5, 7, 1], [2, 2, 2], [17, 16, 15], [16, 33, 2], [40, 3, 9]};
            for s = 1 : numel(shapes)
                sz = shapes{s};
                for kind = 1 : 4
                    switch kind
                        case 1
                            A = single(rand(sz));
                        case 2
                            % ties
                            A = single(randi([0, 4], sz));
                        case 3
                            A = single(randn(sz) * 100);
                            A(rand(sz) < 0.05) = Inf;
                            A(rand(sz) < 0.05) = -Inf;
                        case 4
                            A = zeros(sz, 'single');
                    end
                    for nt = [1, 3]
                        testCase.verifyEqual(omw_backprojector_mex('medfilt3', A, nt), medfilt3(A, [3, 3, 3]), ...
                            sprintf('size %s, kind %d, %d threads', mat2str(sz), kind, nt));
                    end
                end
            end
            % NaN and -0 are left to medfilt3
            testCase.verifyEmpty(omw_backprojector_mex('medfilt3', single(cat(3, [1, NaN; 2, 3], [4, 5; 6, 7])), 2));
            testCase.verifyEmpty(omw_backprojector_mex('medfilt3', single(cat(3, [1, -0; 2, 3], [4, 5; 6, 7])), 2));
        end

        function psfGenerationMatchesMedfilt3(testCase)
            % psf_gen_new with the mex's median filter and with medfilt3
            % (the mex replaced by a function of that name)
            testCase.assumeMexAvailable('omw_backprojector_mex');
            psfs = {testCase.makePSF('lobes3sheared'), testCase.makePSF('elongated')};
            raws = cellfun(@(p) p * 10000 + 100 + single(randi([0, 20], size(p))), psfs, 'unif', 0);
            withMex = cellfun(@(r) psf_gen_new(r, 0.3, 0.3 * sind(32.45), 1.5, 'masked'), raws, 'unif', 0);
            testCase.applyFixture(pktest.FunctionStubFixture('omw_backprojector_mex', []));
            for i = 1 : numel(raws)
                testCase.verifyEqual(withMex{i}, psf_gen_new(raws{i}, 0.3, 0.3 * sind(32.45), 1.5, 'masked'));
            end
        end

        function psfGenerationDoesNotCallMedfilt3(testCase)
            testCase.assumeMexAvailable('omw_backprojector_mex');
            psf = testCase.makePSF('lobes3sheared');
            raw = psf * 10000 + 100 + single(randi([0, 20], size(psf)));
            expected = psf_gen_new(raw, 0.3, 0.3 * sind(32.45), 1.5, 'masked');
            testCase.applyFixture(pktest.FunctionStubFixture('medfilt3', zeros(size(raw), 'single')));
            testCase.verifyEqual(psf_gen_new(raw, 0.3, 0.3 * sind(32.45), 1.5, 'masked'), expected);
        end

        function windowBoundsHoldForSparseEdgePoints(testCase)
            % few edge points with very different radii: a voxel the mex
            % decides is below l times both radii of the 2 edge points
            % nearest to its direction, or at least u times both (then any
            % positive weights give that window value)
            testCase.assumeMexAvailable('omw_backprojector_mex');
            psz = [41, 38, 36];
            c = (psz + 1) / 2;
            [Y, X, Z] = ndgrid((1 : psz(1)) - c(1), (1 : psz(2)) - c(2), (1 : psz(3)) - c(3));
            [th, ph, rs] = cart2sph(X, Y, Z);
            [Xn, Yn, Zn] = sph2cart(th, ph, ones(psz));
            for trial = 1 : 20
                nE = randi([3, 300]);
                E = randn(nE, 3);
                E = E ./ sqrt(sum(E .^ 2, 2));
                R = 2 + 30 * rand(nE, 1);
                bounds = sort(0.2 + 1.5 * rand(1, 2));
                [w, uInds] = omw_backprojector_mex('window', psz, E, R, bounds, 3);
                decided = true(psz);
                decided(uInds) = false;
                inds = knnsearch(E, [Yn(decided), Xn(decided), Zn(decided)], 'K', 2);
                r = rs(decided);
                lo = min(R(inds), [], 2);
                hi = max(R(inds), [], 2);
                one = w(decided) == 1;
                testCase.verifyTrue(all(r(one) < bounds(1) * lo(one)), sprintf('trial %d: window 1', trial));
                testCase.verifyTrue(all(r(~one) >= bounds(2) * hi(~one)), sprintf('trial %d: window 0', trial));
            end
        end

        function invalidMexInputsError(testCase)
            testCase.assumeMexAvailable('omw_backprojector_mex');
            B = rand(8, 8, 8) > 0.5;
            testCase.verifyError(@() omw_backprojector_mex('convhull', double(B), 1, 1, 1, 1), 'omw_backprojector_mex:invalidInput');
            testCase.verifyError(@() omw_backprojector_mex('convhull', B, 9, 1, 1, 1), 'omw_backprojector_mex:invalidInput');
            testCase.verifyError(@() omw_backprojector_mex('convhull', B, [2, 2], 1, 1, 1), 'omw_backprojector_mex:invalidInput');
            testCase.verifyError(@() omw_backprojector_mex('components', B, 'other', 1, 1), 'omw_backprojector_mex:invalidInput');
            testCase.verifyError(@() omw_backprojector_mex('components', B, 'minsize', 1.5, 1), 'omw_backprojector_mex:invalidInput');
            testCase.verifyError(@() omw_backprojector_mex('window', [8, 8, 8], [1, 0, 0], 1, [0.8, 1], 1), 'omw_backprojector_mex:invalidInput');
            testCase.verifyError(@() omw_backprojector_mex('window', [8, 8, 8], [1, 0, 0; 0, 1, 0], [1; 2], [1, 0.8], 1), 'omw_backprojector_mex:invalidInput');
            testCase.verifyError(@() omw_backprojector_mex('sortdesc', double([1, 2]), 1), 'omw_backprojector_mex:invalidInput');
            testCase.verifyError(@() omw_backprojector_mex('medfilt3', rand(4, 4, 4), 1), 'omw_backprojector_mex:invalidInput');
            testCase.verifyError(@() omw_backprojector_mex('other'), 'omw_backprojector_mex:invalidInput');
        end
    end

    methods (Access = private)
        function [expected, actual] = bothVersions(~, psf, alpha, skewed, args)
            expected = cell(1, 4);
            actual = cell(1, 4);
            [~, expected{:}] = evalc('omw_backprojector_generation_before(psf, alpha, skewed, args{:})');
            [~, actual{:}] = evalc('omw_backprojector_generation(psf, alpha, skewed, args{:})');
        end

        function verifyOutputsEqual(testCase, actual, expected)
            names = {'b_omw', 'OTF_bp_omw', 'abs_OTF_c', 'OTF_mask'};
            for i = 1 : 4
                testCase.verifyClass(actual{i}, class(expected{i}), names{i});
                testCase.verifySize(actual{i}, size(expected{i}), names{i});
                testCase.verifyTrue(isequaln(actual{i}, expected{i}), sprintf('%s differs', names{i}));
            end
        end

        function [psf] = makePSF(~, kind)
            % Gaussian envelope, with a modulation along (sheared) z for the
            % OTF side lobes of lattice light sheets
            sz = [48, 48, 48];
            [y, x, z] = ndgrid((1 : sz(1)) - floor(sz(1) / 2) - 1, (1 : sz(2)) - floor(sz(2) / 2) - 1, ...
                (1 : sz(3)) - floor(sz(3) / 2) - 1);
            env = exp(-(y .^ 2 / 2 / 1.5 ^ 2 + x .^ 2 / 2 / 1.5 ^ 2 + z .^ 2 / 2 / 6 ^ 2));
            switch kind
                case 'gaussian'
                    psf = pktest.deconPSF([31, 31, 31], [1.5, 2]);
                    return;
                case 'elongated'
                    psf = exp(-(y .^ 2 / 2 / 1.2 ^ 2 + x .^ 2 / 2 / 1.6 ^ 2 + z .^ 2 / 2 / 4 ^ 2));
                case 'lobes3'
                    psf = env .* (1 + 0.9 * cos(2 * pi * 0.25 * z));
                case 'lobes3sheared'
                    psf = env .* (1 + 0.9 * cos(2 * pi * 0.25 * (z * cosd(30) + x * sind(30))));
                case 'lobes5'
                    psf = env .* (1 + 0.6 * cos(2 * pi * 0.2 * z) + 0.6 * cos(2 * pi * 0.4 * z));
            end
            psf = single(psf);
        end
    end
end
