function deskewRotateFrame3D_resample_test()
% Tests for deskewRotateFrame3D: matrix-derived mex path, index conjugation, resampling.
%
% Run:  matlab -batch "setup; deskewRotateFrame3D_resample_test"
%
% The warp mex implements one affine class (dsr_mex_accepts); everything else must go
% through imwarp, and whatever goes through the mex must agree with imwarp on the same
% matrix. imwarp calls are detected by shadowing imwarp with a function that throws
% 'dsr_test:imwarpCalled'.
%
% The oracle (dsr_reference / dsr_matrices) re-derives the transform independently of
% deskewRotateFrame3D and applies it with imwarp, for sample and objective scan. The
% duplication is deliberate, in the spirit of demo_geometric_transformation's old-vs-new
% comparison: it checks the transform itself, not just two interpolators fed one matrix.
% Drift between production and this reference is a test failure by design.

angle = 32.45;
dz = 0.35;
xy = 0.098;
rz = sind(angle) * dz / xy;        % z factor that restores the raw sampling density (~1.916)
rs_native = [1, 1, rz];

vol = synthetic_volume(128, 96, 81);

shadow_dir = fullfile(tempname, 'imwarp_shadow');
mkdir(shadow_dir);
fid = fopen(fullfile(shadow_dir, 'imwarp.m'), 'w');
fprintf(fid, 'function varargout = imwarp(varargin)\nerror(''dsr_test:imwarpCalled'', ''imwarp called'');\nend\n');
fclose(fid);
cleanup = onCleanup(@() cleanup_shadow(shadow_dir));

results = {};
results{end+1} = run_case('mex_accepts_classifies', @() case_mex_accepts_classifies(vol, angle, dz, xy, rs_native));
results{end+1} = run_case('native_z_uses_fast_path', @() case_native_z_uses_fast_path(vol, angle, dz, xy, rs_native, shadow_dir));
results{end+1} = run_case('native_z_matches_imwarp_reference', @() case_native_z_matches_imwarp_reference(vol, angle, dz, xy, rs_native, shadow_dir));
results{end+1} = run_case('reverse_false_default_matches_imwarp', @() case_reverse_false_default_matches_imwarp(vol, angle, dz, xy, shadow_dir));
results{end+1} = run_case('x_only_factor_fast_and_matches', @() case_x_only_factor_fast_and_matches(vol, angle, dz, xy, shadow_dir));
results{end+1} = run_case('y_factor_still_imwarp', @() case_y_factor_still_imwarp(vol, angle, dz, xy, shadow_dir));
results{end+1} = run_case('objective_scan_rotation_uses_imwarp_and_matches', @() case_objective_scan_rotation(vol, angle, dz, xy, rs_native, shadow_dir));
results{end+1} = run_case('double_input_uses_imwarp_matches', @() case_double_input_uses_imwarp(vol, angle, dz, xy, rs_native, shadow_dir));
results{end+1} = run_case('native_z_voxel_count', @() case_native_z_voxel_count(vol, angle, dz, xy, rs_native, shadow_dir));
results{end+1} = run_case('identity_factor_is_byte_identical', @() case_identity_factor_is_byte_identical(vol, angle, dz, xy, shadow_dir));
% last: an overrunning kernel can corrupt the heap and take the runner down with it
results{end+1} = run_case('odd_row_single_precision_completes', @() case_odd_row_single_precision(angle, dz, xy, rs_native, shadow_dir));

failed = sum(cellfun(@(r) ~r.ok, results));
fprintf('%d/%d passed\n', numel(results) - failed, numel(results));
if failed > 0
    error('dsr_test:failed', '%d test(s) failed', failed);
end
end


%% ---------------------------------------------------------------- cases

function case_mex_accepts_classifies(vol, angle, dz, xy, rs_native)
sz = size(vol);
T = @(rev, rs, objective) mex_matrix(sz, angle, dz, xy, rev, rs, objective);
assert(dsr_mex_accepts(T(false, [], false)), 'sample scan, reverse=false rejected');
assert(dsr_mex_accepts(T(true, [], false)), 'sample scan, reverse=true rejected');
assert(dsr_mex_accepts(T(true, rs_native, false)), 'z-only factor rejected');
assert(dsr_mex_accepts(T(true, [2, 1, 1], false)), 'x-only factor rejected');
assert(~dsr_mex_accepts(T(true, [1, 2, 1], false)), 'y factor accepted (breaks dim-1 pass-through)');
assert(~dsr_mex_accepts(T(true, [], true)), 'objective-scan rotation accepted (drops T(2,2)=cos)');
Rgen = [cosd(20) 0 -sind(20) 0; 0 cosd(10) 0 0; sind(20) 0 cosd(20) 0; 0 0 0 1];
assert(~dsr_mex_accepts(Rgen), 'generic rotation accepted');
end


function case_native_z_uses_fast_path(vol, angle, dz, xy, rs, shadow_dir)
for reverse = [false, true]
    with_shadow(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
        'reverse', reverse, 'resampleFactor', rs, 'save16bit', true));
end
end


function case_native_z_matches_imwarp_reference(vol, angle, dz, xy, rs, shadow_dir)
for reverse = [false, true]
    out = with_shadow(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
        'reverse', reverse, 'resampleFactor', rs, 'save16bit', true));
    ref = dsr_reference(vol, angle, dz, xy, reverse, rs, false);
    assert_close(out, ref, 2, sprintf('native z, reverse=%d', reverse));
end
end


function case_reverse_false_default_matches_imwarp(vol, angle, dz, xy, shadow_dir)
% Default path, reverse=false: the shear must be applied to 1-based planes like imwarp.
out = with_shadow(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, 'reverse', false, 'save16bit', true));
ref = dsr_reference(vol, angle, dz, xy, false, [], false);
assert_close(out, ref, 2, 'reverse=false default');
end


function case_x_only_factor_fast_and_matches(vol, angle, dz, xy, shadow_dir)
rs = [2, 1, 1];
for reverse = [false, true]
    out = with_shadow(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
        'reverse', reverse, 'resampleFactor', rs, 'save16bit', true));
    ref = dsr_reference(vol, angle, dz, xy, reverse, rs, false);
    assert_close(out, ref, 2, sprintf('x-only factor, reverse=%d', reverse));
end
end


function case_y_factor_still_imwarp(vol, angle, dz, xy, shadow_dir)
assert_calls_imwarp(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
    'resampleFactor', [1, 2, 1], 'save16bit', true), 'resampleFactor=[1 2 1]');
end


function case_objective_scan_rotation(vol, angle, dz, xy, rs_native, shadow_dir)
% Objective scan with rotation: input x depends on output x, which the mex drops.
for rs = {[], rs_native}
    assert_calls_imwarp(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
        'objectiveScan', true, 'resampleFactor', rs{1}, 'save16bit', true), ...
        sprintf('objectiveScan rs=%s', mat2str(rs{1})));
    out = deskewRotateFrame3D(vol, angle, dz, xy, 'objectiveScan', true, 'resampleFactor', rs{1}, 'save16bit', true);
    ref = dsr_reference(vol, angle, dz, xy, false, rs{1}, true);
    assert(isequal(size(out), size(ref)), 'objectiveScan rs=%s: size %s vs %s', mat2str(rs{1}), mat2str(size(out)), mat2str(size(ref)));
    assert_close(out, ref, 2, sprintf('objectiveScan rs=%s', mat2str(rs{1})));
end
end


function case_double_input_uses_imwarp(vol, angle, dz, xy, rs, shadow_dir)
% The mex has no double kernel: double input goes to imwarp with the untouched transform.
vold = double(vol);
assert_calls_imwarp(shadow_dir, @() deskewRotateFrame3D(vold, angle, dz, xy, ...
    'reverse', true, 'resampleFactor', rs, 'save16bit', false), 'double input');
out = deskewRotateFrame3D(vold, angle, dz, xy, 'reverse', true, 'resampleFactor', rs, 'save16bit', false);
ref = dsr_reference(vold, angle, dz, xy, true, rs, false);
assert(isequal(size(out), size(ref)), 'double: size %s vs %s', mat2str(size(out)), mat2str(size(ref)));
assert_close(out, ref, 1e-6 * double(max(ref(:))), 'double input');
end


function case_native_z_voxel_count(vol, angle, dz, xy, rs, shadow_dir)
out_native = with_shadow(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
    'resampleFactor', rs, 'save16bit', true));
out_iso = deskewRotateFrame3D(vol, angle, dz, xy, 'save16bit', true);
n = numel(vol);
assert(numel(out_native) <= 1.2 * n, 'native-z output %d voxels > 1.2x raw %d', numel(out_native), n);
assert(numel(out_iso) >= 1.8 * n, 'isotropic output %d voxels < 1.8x raw %d', numel(out_iso), n);
assert(abs(size(out_iso, 3) / size(out_native, 3) - rs(3)) < 0.1, ...
    'z ratio iso/native %g != rz %g', size(out_iso, 3) / size(out_native, 3), rs(3));
end


function case_identity_factor_is_byte_identical(vol, angle, dz, xy, shadow_dir)
for reverse = [false, true]
    a = deskewRotateFrame3D(vol, angle, dz, xy, 'reverse', reverse, 'save16bit', true);
    b = with_shadow(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
        'reverse', reverse, 'resampleFactor', [1, 1, 1], 'save16bit', true));
    assert(isequal(a, b), 'reverse=%d: rs=[1 1 1] differs from rs=[] in %d voxels', reverse, sum(a(:) ~= b(:)));
end
end


function case_odd_row_single_precision(angle, dz, xy, rs, shadow_dir)
% Row lengths not divisible by 8 must not overrun the AVX kernels (float path), with and
% without a resample factor and with a cropped output width.
for ny = [17, 23]
    v = single(synthetic_volume(ny, 40, 33));   % bounds test: only the row length matters
    for rsi = {[], rs}
        out = with_shadow(shadow_dir, @() deskewRotateFrame3D(v, angle, dz, xy, ...
            'reverse', true, 'resampleFactor', rsi{1}, 'save16bit', false));
        ref = dsr_reference(v, angle, dz, xy, true, rsi{1}, false);
        assert(isequal(size(out), size(ref)), 'ny=%d rs=%s: size %s vs %s', ny, mat2str(rsi{1}), mat2str(size(out)), mat2str(size(ref)));
        assert_close(out, ref, 1e-3 * double(max(ref(:))), sprintf('single ny=%d rs=%s', ny, mat2str(rsi{1})));
    end
    full = dsr_reference(v, angle, dz, xy, true, [], false);
    bbox = [1, 3, 1, ny, 3 + 21, size(full, 3)];  % width 22, not divisible by 8
    out = with_shadow(shadow_dir, @() deskewRotateFrame3D(v, angle, dz, xy, 'reverse', true, 'bbox', bbox, 'save16bit', false));
    ref = full(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6));
    assert_close(out, ref, 1e-3 * double(max(ref(:))), sprintf('single ny=%d bbox', ny));
end
end


%% -------------------------------------------------------------- helpers

function vol = synthetic_volume(ny, nx, nz)
% smooth blobs on a gradient so interpolation differences stay small
rng(7);
[Y, X, Z] = ndgrid(1 : ny, 1 : nx, 1 : nz);
vol = 200 + 0.3 * X + 0.2 * Z;
for k = 1 : 12
    c = [randi([8, max(9, ny - 8)]), randi([8, nx - 8]), randi([6, nz - 6])];
    s = 3 + 4 * rand;
    vol = vol + 3000 * rand * exp(-((Y - c(1)).^2 + (X - c(2)).^2 + (Z - c(3)).^2) / (2 * s^2));
end
vol = uint16(vol);
end


function varargout = with_shadow(shadow_dir, fn)
addpath(shadow_dir, '-begin');
c = onCleanup(@() rmpath(shadow_dir));
[varargout{1 : max(1, nargout)}] = fn();
end


function assert_calls_imwarp(shadow_dir, fn, what)
called = false;
try
    with_shadow(shadow_dir, fn);
catch ME
    if ~strcmp(ME.identifier, 'dsr_test:imwarpCalled')
        rethrow(ME);
    end
    called = true;
end
assert(called, '%s did not call imwarp', what);
end


function assert_close(a, b, tol, what)
% interior only: the mex and imwarp differ in edge handling by design (2-voxel border)
assert(isequal(size(a), size(b)), '%s: size %s vs reference %s', what, mat2str(size(a)), mat2str(size(b)));
a = double(a(3 : end - 2, 3 : end - 2, 3 : end - 2));
b = double(b(3 : end - 2, 3 : end - 2, 3 : end - 2));
d = abs(a(:) - b(:));
assert(max(d) <= tol, '%s: max |diff| %g > %g (rms %.3f)', what, max(d), tol, sqrt(mean(d .^ 2)));
end


function cleanup_shadow(shadow_dir)
if any(strcmp(strsplit(path, pathsep), shadow_dir))
    rmpath(shadow_dir);
end
rmdir(fileparts(shadow_dir), 's');
end


function r = run_case(name, fn)
% prints as it goes: a kernel overrun can kill the process before the summary
r = struct('name', name, 'ok', true, 'msg', '');
try
    fn();
    fprintf('PASS  %s\n', name);
catch ME
    r.ok = false;
    r.msg = ME.message;
    fprintf('FAIL  %s\n      %s\n', name, ME.message);
end
end


function [M, outSize, vol_1] = dsr_matrices(vol, angle, dz, xyPixelSize, Reverse, rs, objectiveScan, materialize)
% The imwarp forward transform (1-based pixel centres) exactly as deskewRotateFrame3D
% composes it: skewed-space interpolation to a finer dz when the per-slice shift exceeds
% xStepThresh (sample scan only), then shear * z-scale * rotate * resample.
xStepThresh = 2.0;
[ny, nx, nz] = size(vol);
theta = angle * pi / 180;
dx = cos(theta) * dz / xyPixelSize;
if objectiveScan
    zAniso = dz / xyPixelSize;
    outSize = round([ny, nx * cos(theta) + nz * zAniso * sin(abs(theta)), nz * zAniso * cos(theta) + nx * sin(abs(theta))]);
else
    zAniso = sin(abs(theta)) * dz / xyPixelSize;
    outSize = round([ny, (nx - 1) * cos(theta) + (nz - 1) * zAniso / sin(abs(theta)), (nx - 1) * sin(abs(theta)) - 4]);
end

vol_1 = vol;
if ~objectiveScan && abs(dx) > xStepThresh
    if abs(dx) / xStepThresh < 1.5
        dzout_thresh = xyPixelSize * xStepThresh / cos(theta);
        dzout = dz / ceil(dz / dzout_thresh);
        counter = 1;
        while dzout / dzout_thresh < 0.4
            dzout = dz / (ceil((dz / dzout_thresh) / (1 / (counter + 1))) * (1 / (counter + 1)));
            counter = counter + 1;
        end
    else
        dzout = dz / ceil(abs(dx) / xStepThresh);
    end
    int_stepsize = dzout / dz;
    if materialize
        if isa(vol, 'single') || isa(vol, 'uint16')
            vol_1 = skewed_space_interp_defined_stepsize_mex(vol, abs(dx), int_stepsize, Reverse, isa(vol, 'uint16'));
        else
            vol_1 = skewed_space_interp_defined_stepsize(vol, abs(dx), int_stepsize, 'Reverse', Reverse);
        end
        nz = size(vol_1, 3);
    else
        nz = floor(round((nz - 1) / int_stepsize * 100000) / 100000) + 1;
    end
    dz = dzout;
    dx = cos(theta) * dz / xyPixelSize;
    zAniso = sin(abs(theta)) * dz / xyPixelSize;
end

if ~Reverse
    xshift = -dx;
    xstep = dx;
else
    xshift = dx + ceil((nz - 1) * dx);
    xstep = -dx;
end
nxDs = ceil((nz - 1) * dx) + nx;
if objectiveScan
    nxDs = nx;
    ds_S = eye(4);
else
    ds_S = [1 0 0 0; 0 1 0 0; xstep 0 1 0; xshift 0 0 1];
end
if Reverse
    theta = -theta;
end
center = ([ny nxDs nz] + 1) / 2;
T1 = [1 0 0 0; 0 1 0 0; 0 0 1 0; -center([2 1 3]) 1];
S = [1 0 0 0; 0 1 0 0; 0 0 zAniso 0; 0 0 0 1];
R = [cos(theta) 0 -sin(theta) 0; 0 1 0 0; sin(theta) 0 cos(theta) 0; 0 0 0 1];
T2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; (outSize([2 1 3]) + 1) / 2 1];
if ~isempty(rs)
    RT1 = [1 0 0 0; 0 1 0 0; 0 0 1 0; -(outSize([2, 1, 3]) + 1) / 2 1];
    RS = [1 / rs(1) 0 0 0; 0 1 / rs(2) 0 0; 0 0 1 / rs(3) 0; 0 0 0 1];
    outSize = round(outSize ./ rs([2, 1, 3]));
    RT2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; (outSize([2, 1, 3]) + 1) / 2 1];
else
    RT1 = eye(4); RS = eye(4); RT2 = eye(4);
end
M = ds_S * (T1 * S * R * T2) * (RT1 * RS * RT2);
end


function T = mex_matrix(sz, angle, dz, xy, Reverse, rs, objectiveScan)
% 0-based backward map in mex axis order: what deskewRotateFrame3D hands the kernel.
M = dsr_matrices(zeros(sz, 'uint16'), angle, dz, xy, Reverse, rs, objectiveScan, false);
P = eye(4); P(4, 1 : 3) = 1;
T = eye(4) / (P * M / P)';
T = T([2, 1, 3, 4], [2, 1, 3, 4]);
end


function volout = dsr_reference(vol, angle, dz, xyPixelSize, Reverse, rs, objectiveScan)
[M, outSize, vol_1] = dsr_matrices(vol, angle, dz, xyPixelSize, Reverse, rs, objectiveScan, true);
RA = imref3d(outSize, 1, 1, 1);
volout = imwarp(vol_1, affine3d(M), 'linear', 'FillValues', 0, 'OutputView', RA);
end
