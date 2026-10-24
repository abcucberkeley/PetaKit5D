classdef DeskewRotateKernelExactTest < pktest.TestCase
    % deskewRotateFrame3D's mex path (sample scan, linear, no resampling)
    % against plain MATLAB code that does the mex's single precision
    % arithmetic in the same order: the skewed space interpolation (dz above
    % the x step threshold; uint16 rows are rounded) and the weighted sums of
    % 4 lines of the warp. On x86 the uint16 line sums add the first
    % floor(n / 16) * 16 values of a line as (a*wa + b*wb) + (c*wc + d*wd)
    % (SIMD) and the rest left to right; single lines are summed left to
    % right. So the outputs must be identical there, which catches any
    % change of the kernel's output. On ARM macOS the mex sums left to right
    % and the compiler may fuse multiply-adds, so only rounding differences
    % are allowed there.

    properties (Constant)
        SkewAngle = 32.45
        XYPixelSize = 0.108
    end

    properties (TestParameter)
        % 0.2: no interpolation (volume_deskew_rotate_warp_mex); 0.3: one
        % interpolated slice per step (the default production path); 0.8:
        % three, so neighboring grid slices are both interpolated
        dz = struct('dz02', 0.2, 'dz03', 0.3, 'dz08', 0.8)
        reverse = struct('reverse', true, 'forward', false)
        % uint16 with save16bit, uint16 with single output, single
        dataType = {'uint16', 'uint16ToSingle', 'single'}
        % y (the line length) with and without a remainder of 16
        rawSize = struct('y61', [61, 40, 30], 'y48', [48, 37, 26])
    end

    methods (Test, TestTags = {'Unit'})
        function wholeVolumeMatchesReference(testCase, dz, reverse, dataType, rawSize)
            [vol, save16bit] = testCase.makeVolume(rawSize, dataType);
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, 'reverse', reverse, ...
                'save16bit', save16bit);

            g = testCase.geometry(rawSize, dz, reverse);
            expected = testCase.referenceKernel(vol, g, [1, 1, 1, g.outSize], save16bit);
            testCase.verifyKernelOutput(out, expected);
        end

        function boundingBoxMatchesReference(testCase, dz, reverse, dataType)
            % an output region, as large-file batches compute them: the line
            % sums start at the region's first y
            rawSize = [61, 40, 30];
            [vol, save16bit] = testCase.makeVolume(rawSize, dataType);
            g = testCase.geometry(rawSize, dz, reverse);
            bbox = [5, 7, 3, 57, g.outSize(2) - 9, g.outSize(3) - 2];
            out = deskewRotateFrame3D(vol, testCase.SkewAngle, dz, testCase.XYPixelSize, 'reverse', reverse, ...
                'save16bit', save16bit, 'bbox', bbox);

            expected = testCase.referenceKernel(vol, g, bbox, save16bit);
            testCase.verifyKernelOutput(out, expected);
        end
    end

    methods (Access = private)
        function [vol, save16bit] = makeVolume(~, rawSize, dataType)
            % full uint16 range, so that rounding differences show
            switch dataType
                case 'uint16'
                    vol = uint16(randi([0, 65535], rawSize));
                    save16bit = true;
                case 'uint16ToSingle'
                    vol = uint16(randi([0, 65535], rawSize));
                    save16bit = false;
                case 'single'
                    vol = single(rand(rawSize) * 1000);
                    save16bit = true;
            end
        end

        function [] = verifyKernelOutput(testCase, out, expected)
            testCase.verifyClass(out, class(expected));
            testCase.verifySize(out, size(expected));
            if strcmp(mexext, 'mexmaca64')
                % left-to-right sums, possibly with fused multiply-adds
                if isa(expected, 'uint16')
                    testCase.verifyLessThanOrEqual(max(abs(double(out(:)) - double(expected(:)))), 2);
                else
                    testCase.verifyEqual(out, expected, 'AbsTol', 1e-4 * max(abs(expected(:))));
                end
                return;
            end
            d = find(out ~= expected);
            msg = '';
            if ~isempty(d)
                [iy, ix, iz] = ind2sub(size(out), d(1));
                msg = sprintf('%d values differ, the first at (%d, %d, %d): %.9g instead of %.9g', numel(d), ...
                    iy, ix, iz, double(out(d(1))), double(expected(d(1))));
            end
            testCase.verifyTrue(isequal(out, expected), msg);
        end

        function [g] = geometry(~, rawSize, dz, reverse)
            % the transform deskewRotateFrame3D passes to the mex (sample
            % scan, no resampling)
            px = DeskewRotateKernelExactTest.XYPixelSize;
            theta = DeskewRotateKernelExactTest.SkewAngle * pi / 180;
            ny = rawSize(1);
            nx = rawSize(2);
            nz = rawSize(3);
            dx = cos(theta) * dz / px;
            zAniso = sin(abs(theta)) * dz / px;
            outSize = round([ny, (nx - 1) * cos(theta) + (nz - 1) * zAniso / sin(abs(theta)), (nx - 1) * sin(abs(theta)) - 4]);

            xStepThresh = 2.0;
            stepsize = 1;
            xstepOrig = abs(dx);
            if abs(dx) > xStepThresh
                if abs(dx) / xStepThresh < 1.5
                    dzout_thresh = px * xStepThresh / cos(theta);
                    dzout = dz / ceil(dz / dzout_thresh);
                    counter = 1;
                    while dzout / dzout_thresh < 0.4
                        dzout = dz / (ceil((dz / dzout_thresh) / (1 / (counter + 1))) * (1 / (counter + 1)));
                        counter = counter + 1;
                    end
                else
                    dzout = dz / ceil(abs(dx) / xStepThresh);
                end
                stepsize = dzout / dz;
                nz = floor(round((nz - 1) / stepsize * 100000) / 100000) + 1;
                dz = dzout;
                dx = cos(theta) * dz / px;
                zAniso = sin(abs(theta)) * dz / px;
            end

            if ~reverse
                xshift = -dx;
                xstep = dx;
            else
                xshift = dx + ceil((nz - 1) * dx);
                xstep = -dx;
                theta = -theta;
            end
            nxDs = ceil((nz - 1) * dx) + nx;
            ds_S = [1 0 0 0; 0 1 0 0; xstep 0 1 0; xshift 0 0 1];
            center = ([ny nxDs nz] + 1) / 2;
            T1 = [1 0 0 0; 0 1 0 0; 0 0 1 0; -center([2 1 3]) 1];
            S = [1 0 0 0; 0 1 0 0; 0 0 zAniso 0; 0 0 0 1];
            R = [cos(theta) 0 -sin(theta) 0; 0 1 0 0; sin(theta) 0 cos(theta) 0; 0 0 0 1];
            T2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; (outSize([2 1 3]) + 1) / 2 1];
            offset = zeros(4, 4);
            offset(4, 1 : 3) = 1;
            ds_S_c = ds_S;
            ds_S_c(4, 1) = ds_S_c(4, 1) + xstep;
            tmat = eye(4) / (ds_S_c * ((T1 + offset) * S * R * (T2 - offset)))';
            tmat = tmat([2, 1, 3, 4], [2, 1, 3, 4]);

            g = struct('outSize', outSize, 'tmat', tmat, 'xstep', xstepOrig, 'stepsize', stepsize, 'reverse', reverse);
        end

        function [out] = referenceKernel(testCase, vol, g, bbox, save16bit)
            % the mex's arithmetic, in its order (C indices are 0-based)
            [ny, nx, nz] = size(vol);
            tmat = g.tmat;
            if ~isa(vol, 'uint16')
                kind = 'single';
            elseif save16bit
                kind = 'uint16';
            else
                kind = 'mixed';
            end
            simd = ~strcmp(mexext, 'mexmaca64');

            % slices of the interpolation grid: original or interpolated
            stepsize = g.stepsize;
            ntotal = floor(round((nz - 1) / stepsize * 100000) / 100000) + 1;
            isOrig = false(1, ntotal);
            idx = zeros(1, ntotal);
            nint = 0;
            for i = 0 : ntotal - 1
                cs = i * stepsize;
                if abs(cs - round(cs)) < 1e-3
                    isOrig(i + 1) = true;
                    idx(i + 1) = round(cs);
                else
                    idx(i + 1) = nint;
                    nint = nint + 1;
                end
            end

            % interpolated slices, every row (uint16 rows are rounded)
            if strcmp(kind, 'uint16')
                im_int = zeros(ny, nx, nint, 'uint16');
            else
                im_int = zeros(ny, nx, nint, 'single');
            end
            xstep = g.xstep;
            xa = ceil(xstep);
            raw = @(zr, r) single(vol(:, r + 1, zr + 1));
            for i = 0 : ntotal - 1
                if isOrig(i + 1)
                    continue;
                end
                zout = stepsize * i;
                zw = zout - floor(zout);
                if g.reverse
                    s = ceil(xstep) - xstep * zw;
                    t = xstep * (1 - zw);
                else
                    s = xstep * (1 - zw);
                    t = ceil(xstep) - xstep * zw;
                end
                sw = s - floor(s);
                s = floor(s);
                tw = t - floor(t);
                t = floor(t);
                swf = double(single(sw));
                twf = double(single(tw));
                zwf = double(single(zw));
                c00 = single((1 - swf) * (1 - zwf));
                c01 = single(swf * (1 - zwf));
                c10 = single((1 - twf) * zwf);
                c11 = single(twf * zwf);
                zorg = floor(zout);
                % uint64 in the mex: no negative limit
                lim = nx - t - 1;
                if lim < 0
                    lim = inf;
                end
                for y = 0 : nx - 1
                    if g.reverse
                        if y < xa - s
                            v = testCase.sum2(raw(zorg + 1, y + t), raw(zorg + 1, y + t + 1), c10, c11, kind);
                        elseif y < lim
                            v = testCase.sum4(raw(zorg, y + s - xa), raw(zorg, y + s - xa + 1), raw(zorg + 1, y + t), ...
                                raw(zorg + 1, y + t + 1), c00, c01, c10, c11, kind, simd);
                        else
                            v = testCase.sum2(raw(zorg, y + s - xa), raw(zorg, y + s - xa + 1), c00, c01, kind);
                        end
                    else
                        if y < xa - s
                            v = testCase.sum2(raw(zorg, y + s), raw(zorg, y + s + 1), c00, c01, kind);
                        elseif y < lim
                            v = testCase.sum4(raw(zorg, y + s), raw(zorg, y + s + 1), raw(zorg + 1, y + t - xa), ...
                                raw(zorg + 1, y + t - xa + 1), c00, c01, c10, c11, kind, simd);
                        else
                            v = testCase.sum2(raw(zorg + 1, y + t - xa), raw(zorg + 1, y + t - xa + 1), c10, c11, kind);
                        end
                    end
                    im_int(:, y + 1, idx(i + 1) + 1) = v;
                end
            end

            % warp: output line (y, z) from rows y0, y0 + 1 of grid slices z0, z0 + 1
            startX = single(bbox(1)) - 1;
            startY = single(bbox(2)) - 1;
            startZ = single(bbox(3)) - 1;
            outShape = floor(double([single(bbox(4)) - startX, single(bbox(5)) - startY, single(bbox(6)) - startZ]));
            x0 = floor(double(startX));
            xr = x0 + 1 : x0 + outShape(1);
            if strcmp(kind, 'uint16')
                out = zeros(outShape, 'uint16');
            else
                out = zeros(outShape, 'single');
            end
            lineOf = @(zq, r) testCase.gridRow(vol, im_int, isOrig, idx, zq, r, xr);
            for z = 0 : outShape(3) - 1
                zf = single(z) + startZ;
                p = tmat(2, 3) * double(zf);
                ty = p + tmat(2, 4);
                if ty < 0 || ty + 1 >= nx
                    if ty + 1 == nx
                        counter = 1;
                        while floor(ty) + 1 == nx
                            ty = ty - 1e-9 * counter;
                            counter = counter * 2;
                        end
                    else
                        continue;
                    end
                end
                y0 = floor(ty);
                wy = single(ty - y0);
                for y = 0 : outShape(2) - 1
                    p1 = tmat(3, 2) * double(single(y) + startY);
                    p2 = tmat(3, 3) * double(zf);
                    tz = (p1 + p2) + tmat(3, 4);
                    if tz < 0 || tz + 1 >= ntotal
                        if tz + 1 == ntotal
                            counter = 1;
                            while floor(tz) + 1 == ntotal
                                tz = tz - 1e-9 * counter;
                                counter = counter * 2;
                            end
                        else
                            continue;
                        end
                    end
                    z0 = floor(tz);
                    wz = single(tz - z0);
                    w00 = (single(1) - wy) * (single(1) - wz);
                    w01 = (single(1) - wy) * wz;
                    w10 = wy * (single(1) - wz);
                    w11 = wy * wz;
                    L00 = lineOf(z0, y0);
                    L01 = lineOf(z0 + 1, y0);
                    L10 = lineOf(z0, y0 + 1);
                    L11 = lineOf(z0 + 1, y0 + 1);
                    switch kind
                        case 'uint16'
                            v = testCase.sum4(L00, L01, L10, L11, w00, w01, w10, w11, 'uint16', simd);
                        case 'single'
                            v = testCase.sum4(L00, L01, L10, L11, w00, w01, w10, w11, 'single', simd);
                        case 'mixed'
                            % the mex's uint16/single row combinations
                            if isOrig(z0 + 1) && ~isOrig(z0 + 2)
                                v = testCase.sum4(L00, L01, L10, L11, w00, w01, w10, w11, 'mixed', simd);
                            elseif ~isOrig(z0 + 1) && isOrig(z0 + 2)
                                v = testCase.sum4(L01, L00, L11, L10, w01, w00, w11, w10, 'mixed', simd);
                            elseif isOrig(z0 + 1)
                                v = testCase.sum4(L00, L01, L10, L11, w00, w01, w10, w11, 'mixed', simd);
                            else
                                v = testCase.sum4(L00, L01, L10, L11, w00, w01, w10, w11, 'single', simd);
                            end
                    end
                    out(:, y + 1, z + 1) = v;
                end
            end
        end

        function [v] = gridRow(~, vol, im_int, isOrig, idx, zq, r, xr)
            % row r (0-based) of grid slice zq (0-based), values xr, single
            if isOrig(zq + 1)
                v = single(vol(xr, r + 1, idx(zq + 1) + 1));
            else
                v = single(im_int(xr, r + 1, idx(zq + 1) + 1));
            end
        end

        function [r] = sum2(~, a, b, wa, wb, kind)
            % compute_weighted_sum_two_vectors*: a*wa + b*wb (+ 0.5, truncated, for uint16)
            pa = a .* wa;
            pb = b .* wb;
            r = pa + pb;
            if strcmp(kind, 'uint16')
                r = r + 0.5;
                r = uint16(floor(r));
            end
        end

        function [r] = sum4(~, a, b, c, d, wa, wb, wc, wd, kind, simd)
            % compute_weighted_sum*: single lines left to right; uint16 lines
            % (and uint16/single mixes) with SIMD pairwise sums for the
            % first floor(n / 16) * 16 values on x86
            n = numel(a);
            nb = 0;
            if simd && ~strcmp(kind, 'single')
                nb = n - mod(n, 16);
            end
            r = zeros(n, 1, 'single');
            if nb > 0
                i = 1 : nb;
                pa = a(i) .* wa;
                pb = b(i) .* wb;
                pc = c(i) .* wc;
                pd = d(i) .* wd;
                ab = pa + pb;
                cd = pc + pd;
                r(i) = ab + cd;
            end
            i = nb + 1 : n;
            pa = wa .* a(i);
            pb = wb .* b(i);
            pc = wc .* c(i);
            pd = wd .* d(i);
            s = pa + pb;
            s = s + pc;
            r(i) = s + pd;
            if strcmp(kind, 'uint16')
                r = r + 0.5;
                r = uint16(floor(r));
            end
        end
    end
end
