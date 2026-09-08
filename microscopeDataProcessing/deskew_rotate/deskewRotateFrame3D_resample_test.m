function deskewRotateFrame3D_resample_test()
% Tests for deskewRotateFrame3D with resampleFactor on the fast (mex) path.
%
% Run:  matlab -batch "run setup.m; deskewRotateFrame3D_resample_test"
%
% A z-only (or x-only) resampleFactor keeps the sparsity the DSR mex hard-codes
% (input x depends only on output z; input z on output x,z; y passes through),
% so those factors must run on the fused mex path rather than the materialized
% skewed-space interpolation + imwarp path. A y factor breaks the pass-through
% and must stay on imwarp.
%
% imwarp calls are detected by shadowing imwarp with a function that throws
% 'dsr_test:imwarpCalled' (deskewRotateFrame3D's mex try/catch falls back to
% imwarp, so a broken mex is caught too).

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
results{end+1} = run_case('native_z_uses_fast_path', @() case_native_z_uses_fast_path(vol, angle, dz, xy, rs_native, shadow_dir));
results{end+1} = run_case('native_z_matches_imwarp_reference', @() case_native_z_matches_imwarp_reference(vol, angle, dz, xy, rs_native, shadow_dir));
results{end+1} = run_case('native_z_voxel_count', @() case_native_z_voxel_count(vol, angle, dz, xy, rs_native, shadow_dir));
results{end+1} = run_case('y_factor_still_imwarp', @() case_y_factor_still_imwarp(vol, angle, dz, xy, shadow_dir));
results{end+1} = run_case('objective_scan_unaffected', @() case_objective_scan_unaffected(vol, angle, dz, xy, shadow_dir));
results{end+1} = run_case('identity_factor_is_byte_identical', @() case_identity_factor_is_byte_identical(vol, angle, dz, xy, shadow_dir));

failed = 0;
for i = 1 : numel(results)
    r = results{i};
    if r.ok
        fprintf('PASS  %s\n', r.name);
    else
        failed = failed + 1;
        fprintf('FAIL  %s\n      %s\n', r.name, r.msg);
    end
end
fprintf('%d/%d passed\n', numel(results) - failed, numel(results));
if failed > 0
    error('dsr_test:failed', '%d test(s) failed', failed);
end
end


%% ---------------------------------------------------------------- cases

function case_native_z_uses_fast_path(vol, angle, dz, xy, rs, shadow_dir)
% RED before the fix: any resampleFactor forces the imwarp path.
for reverse = [false, true]
    with_shadow(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
        'reverse', reverse, 'resampleFactor', rs, 'save16bit', true));
end
end


function case_native_z_matches_imwarp_reference(vol, angle, dz, xy, rs, shadow_dir)
% The fused mex with a z-only factor must agree with imwarp on the same affine
% at least as well as the existing (rs = []) fast path does.
% The baseline is measured in-test: with reverse=true the mex and imwarp agree to
% ~1 gray level, so the bound is tight there; with reverse=false the existing
% (rs = []) mex path already differs from imwarp (a pre-existing offset, not
% touched here), so the bound is "no worse than today".
for reverse = [false, true]
    base_fast = deskewRotateFrame3D(vol, angle, dz, xy, 'reverse', reverse, 'save16bit', true);
    base_ref = reference_imwarp(vol, angle, dz, xy, reverse, []);
    [base_max, base_rms] = interior_diff(base_fast, base_ref);

    out_fast = with_shadow(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
        'reverse', reverse, 'resampleFactor', rs, 'save16bit', true));
    out_ref = reference_imwarp(vol, angle, dz, xy, reverse, rs);
    assert(isequal(size(out_fast), size(out_ref)), ...
        'reverse=%d: size %s vs reference %s', reverse, mat2str(size(out_fast)), mat2str(size(out_ref)));
    [err_max, err_rms] = interior_diff(out_fast, out_ref);
    tol_max = max(2, base_max);
    tol_rms = max(1, 1.1 * base_rms);
    assert(err_max <= tol_max, 'reverse=%d: max |diff| %g > %g (rs=[] baseline %g)', reverse, err_max, tol_max, base_max);
    assert(err_rms <= tol_rms, 'reverse=%d: rms diff %.2f > %.2f (rs=[] baseline %.2f)', reverse, err_rms, tol_rms, base_rms);
    fprintf('      reverse=%d: native-z max %g rms %.2f | rs=[] baseline max %g rms %.2f\n', ...
        reverse, err_max, err_rms, base_max, base_rms);
end
end


function case_native_z_voxel_count(vol, angle, dz, xy, rs, shadow_dir)
% Native-z output has ~the raw voxel count (parallelogram box padding only);
% the isotropic default is ~2x.
out_native = with_shadow(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
    'resampleFactor', rs, 'save16bit', true));
out_iso = deskewRotateFrame3D(vol, angle, dz, xy, 'save16bit', true);
n = numel(vol);
assert(numel(out_native) <= 1.2 * n, 'native-z output %d voxels > 1.2x raw %d', numel(out_native), n);
assert(numel(out_iso) >= 1.8 * n, 'isotropic output %d voxels < 1.8x raw %d', numel(out_iso), n);
assert(abs(size(out_iso, 3) / size(out_native, 3) - rs(3)) < 0.1, ...
    'z ratio iso/native %g != rz %g', size(out_iso, 3) / size(out_native, 3), rs(3));
end


function case_y_factor_still_imwarp(vol, angle, dz, xy, shadow_dir)
% Only a z-only factor is on the mex path. An x or y factor (rs(1), rs(2) in this
% function's [x y z] order) must still go through imwarp, as before the change.
for rs = {[1, 2, 1], [2, 1, 1]}
    called = false;
    try
        with_shadow(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
            'resampleFactor', rs{1}, 'save16bit', true));
    catch ME
        if ~strcmp(ME.identifier, 'dsr_test:imwarpCalled')
            rethrow(ME);
        end
        called = true;
    end
    assert(called, 'resampleFactor=%s did not call imwarp', mat2str(rs{1}));
end
end


function case_objective_scan_unaffected(vol, angle, dz, xy, shadow_dir)
out = with_shadow(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
    'objectiveScan', true, 'save16bit', true));
[ny, nx, nz] = size(vol);
theta = angle * pi / 180;
zAniso = dz / xy;
expected = round([ny, nx * cos(theta) + nz * zAniso * sin(abs(theta)), nz * zAniso * cos(theta) + nx * sin(abs(theta))]);
assert(isequal(size(out), expected), 'objective-scan size %s != %s', mat2str(size(out)), mat2str(expected));
end


function case_identity_factor_is_byte_identical(vol, angle, dz, xy, shadow_dir)
% resampleFactor=[1 1 1] goes through the RS/RT1/RT2 block with identity scaling
% and must reproduce the rs=[] output byte for byte: pins the -offset placement
% (any shift from the resample translation would show up here) for both scan
% directions, on the fast path.
for reverse = [false, true]
    a = deskewRotateFrame3D(vol, angle, dz, xy, 'reverse', reverse, 'save16bit', true);
    b = with_shadow(shadow_dir, @() deskewRotateFrame3D(vol, angle, dz, xy, ...
        'reverse', reverse, 'resampleFactor', [1, 1, 1], 'save16bit', true));
    assert(isequal(a, b), 'reverse=%d: rs=[1 1 1] differs from rs=[] in %d voxels', reverse, sum(a(:) ~= b(:)));
end
end


%% -------------------------------------------------------------- helpers

function vol = synthetic_volume(ny, nx, nz)
% smooth blobs on a gradient so interpolation differences stay small
rng(7);
[Y, X, Z] = ndgrid(1 : ny, 1 : nx, 1 : nz);
vol = 200 + 0.3 * X + 0.2 * Z;
for k = 1 : 12
    c = [randi([8, ny - 8]), randi([8, nx - 8]), randi([6, nz - 6])];
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


function cleanup_shadow(shadow_dir)
if any(strcmp(strsplit(path, pathsep), shadow_dir))
    rmpath(shadow_dir);
end
rmdir(fileparts(shadow_dir), 's');
end


function r = run_case(name, fn)
r = struct('name', name, 'ok', true, 'msg', '');
try
    fn();
catch ME
    r.ok = false;
    r.msg = ME.message;
end
end


function [err_max, err_rms] = interior_diff(a, b)
% ignore a 2-voxel border: the mex and imwarp differ in edge handling by design
a = double(a(3 : end - 2, 3 : end - 2, 3 : end - 2));
b = double(b(3 : end - 2, 3 : end - 2, 3 : end - 2));
d = abs(a(:) - b(:));
err_max = max(d);
err_rms = sqrt(mean(d .^ 2));
end


function volout = reference_imwarp(vol, angle, dz, xyPixelSize, Reverse, rs)
% The non-fast branch of deskewRotateFrame3D (skewed-space interpolation to a
% finer dz when the per-slice shift exceeds xStepThresh, then one imwarp with
% shear * z-scale * rotate * resample), written out independently as the oracle.
xStepThresh = 2.0;
[ny, nx, nz] = size(vol);
theta = angle * pi / 180;
dx = cos(theta) * dz / xyPixelSize;
zAniso = sin(abs(theta)) * dz / xyPixelSize;
outSize = round([ny, (nx - 1) * cos(theta) + (nz - 1) * zAniso / sin(abs(theta)), (nx - 1) * sin(abs(theta)) - 4]);

vol_1 = vol;
if abs(dx) > xStepThresh
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
    vol_1 = skewed_space_interp_defined_stepsize_mex(vol, abs(dx), int_stepsize, Reverse, true);
    [ny, nx, nz] = size(vol_1);
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
ds_S = [1 0 0 0; 0 1 0 0; xstep 0 1 0; xshift 0 0 1];
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
RA = imref3d(outSize, 1, 1, 1);
volout = imwarp(vol_1, affine3d(ds_S * (T1 * S * R * T2) * (RT1 * RS * RT2)), 'linear', 'FillValues', 0, 'OutputView', RA);
end
