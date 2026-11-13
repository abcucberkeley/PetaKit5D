function [] = build_variant_mex(outDir, mingwDir, mode)
% for the CI job exit-probe (a diagnostic, not a build of the repository's
% binaries): build the PetaKit5D OpenMP mex files of Windows into outDir from
% copies of their sources changed for mode, with the flags of
% compileAllMexFiles and the MinGW-w64 in mingwDir. The committed mexw64 are
% not touched; a variant's folder goes ahead of them on the path to shadow
% them.
%   mode 'lock':   the mexLock block (macOS only in the committed source) also
%                  applies on Windows, so the mex file is never unloaded
%   mode 'atexit': a mexAtExit handler that ends the OpenMP pool
%                  (omp_pause_resource_all(omp_pause_hard)) on Windows, so the
%                  mex file's OpenMP worker threads end while it is still loaded
%                  (the mex file stays unloadable, unlike with the lock)
%   mode 'both':   both of the above
% The deadlock this probes: a MinGW mex file's static libgomp keeps its OpenMP
% workers alive after the parallel region; when a worker ends through
% FreeLibraryAndExitThread after MATLAB has unloaded the mex file, it waits for
% the Windows loader lock that MATLAB's exiting main thread holds.

if ~ismember(mode, {'lock', 'atexit', 'both'})
    error('build_variant_mex:mode', 'mode must be lock, atexit or both.');
end
setenv('MW_MINGW64_LOC', mingwDir);
evalc('mex(''-setup'', ''C'')');
evalc('mex(''-setup'', ''C++'')');
root = pwd;
srcDir = fullfile(outDir, 'src');
if ~exist(srcDir, 'dir')
    mkdir(srcDir);
end
% folder under microscopeDataProcessing (or a path from the root), source,
% AVX2, extra compiler flags (as in compileAllMexFiles)
srcs = {'deskew_rotate/mex', 'skewed_space_interp_defined_stepsize_mex.c', true, '';
    'deskew_rotate/mex', 'volume_deskew_rotate_warp_mex.cpp', true, '';
    'deskew_rotate/mex', 'skewed_space_interp_volume_deskew_rotate_warp_mex.cpp', true, '';
    'stitch/mex', 'feather_blending_3d_mex.cpp', false, '';
    'stitch/mex', 'feather_blending_3d_with_indexing_mex.cpp', false, '';
    'stitch/mex', 'feather_distance_map_resize_3d_mex.c', false, '';
    'stitch/mex', 'integral_image_3d_mex.cpp', false, '';
    'stitch/mex', 'any_4th_dim_mex.cpp', false, '';
    'stitch/mex', 'tile_blending_3d_mex.cpp', false, ' -ffp-contract=off';
    'tools/MIP/mex', 'max_pooling_3d_mex.cpp', false, '';
    'tools/MIP/mex', 'range_max_3d_mex.cpp', false, '';
    'crop/mex', 'crop3d_mex.c', false, '';
    'crop/mex', 'crop4d_mex.cpp', false, '';
    'crop/mex', 'indexing3d_mex.c', false, '';
    'crop/mex', 'indexing4d_mex.c', false, '';
    'crop/mex', 'indexing4d_crop_mex.cpp', false, '';
    'flatfield_correction/mex', 'flat_field_correction_3d_mex.cpp', false, '';
    'decon/mex', 'decon_mask_edge_erosion_mex.cpp', false, '';
    'decon/mex', 'omw_backprojector_mex.cpp', false, ' -ffp-contract=off';
    'decon/mex', 'decon_lucy_accel_mex.cpp', false, ' -ffp-contract=off';
    'unmix/mex', 'unmix_channels_sum_mex.cpp', false, ' -ffp-contract=off';
    '../pointDetection/mex_functions', 'replace_nan_with_value.c', false, ''};
for i = 1 : size(srcs, 1)
    txt = fileread(fullfile(root, 'microscopeDataProcessing', srcs{i, 1}, srcs{i, 2}));
    txt = variant_mex_source(txt, mode, srcs{i, 2});
    srcFn = fullfile(srcDir, srcs{i, 2});
    fid = fopen(srcFn, 'w');
    fwrite(fid, txt);
    fclose(fid);
    avxFlag = '';
    if srcs{i, 3}
        avxFlag = ' -mavx2';
    end
    if endsWith(srcFn, '.cpp')
        mex('CXXOPTIMFLAGS=-O3 -DNDEBUG', ['CXXFLAGS=$CXXFLAGS -O3 -fopenmp', avxFlag, srcs{i, 4}], ...
            ['LDFLAGS=$LDFLAGS -O3 -fopenmp', avxFlag], '-outdir', outDir, srcFn);
    else
        mex('COPTIMFLAGS=-O3 -DNDEBUG', ['CFLAGS=$CFLAGS -O3 -fopenmp', avxFlag, srcs{i, 4}], ...
            ['LDFLAGS=$LDFLAGS -O3 -fopenmp', avxFlag], '-outdir', outDir, srcFn);
    end
end
fprintf('Built %d %s-variant mex files in %s\n', size(srcs, 1), mode, outDir);

end
