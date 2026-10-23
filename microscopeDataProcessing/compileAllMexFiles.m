function [] = compileAllMexFiles(varargin)
% compile the mex files of PetaKit5D for the current platform
%
% Intended to be run from PetaKit5D/microscopeDataProcessing or pass the
% path to the PetaKit5D base dir. Linux and Windows use the configured mex
% compilers (GCC, MinGW-w64 on Windows, e.g., with MW_MINGW64_LOC set to
% the MinGW folder). ARM macOS uses Homebrew gcc@13 (brew install gcc@13)
% with the GCC runtime libraries linked statically, so that the mex files
% don't need Homebrew to run.

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('PetaKit5DBaseDir', '', @ischar);
ip.parse(varargin{:});

PetaKit5DBaseDir = ip.Results.PetaKit5DBaseDir;

baseDir = pwd;

if(~isempty(PetaKit5DBaseDir))
    cd([PetaKit5DBaseDir '/./microscopeDataProcessing']);
else
    funcFn = which(mfilename);
    path = fileparts(funcFn);
    cd(path);
end

[~,name,~]=fileparts(pwd);
if ~strcmp(name,'microscopeDataProcessing')
    error('Pass the path to the PetaKit5D base dir or be in the PetaKit5D/microscopeDataProcessing dir');
end
srcDir = pwd;

% folder, source file, whether it uses AVX2 on x86, extra mex arguments and
% extra compiler flags
mexFiles = {'deskew_rotate/mex', 'skewed_space_interp_defined_stepsize_mex.c', true, {}, '';
    'deskew_rotate/mex', 'skewed_space_interp_mex.c', false, {}, '';
    'deskew_rotate/mex', 'volume_deskew_rotate_warp_mex.cpp', true, {}, '';
    'deskew_rotate/mex', 'skewed_space_interp_volume_deskew_rotate_warp_mex.cpp', true, {}, '';
    'stitch/mex', 'feather_blending_3d_mex.cpp', false, {}, '';
    'stitch/mex', 'feather_blending_3d_with_indexing_mex.cpp', false, {}, '';
    'stitch/mex', 'feather_distance_map_resize_3d_mex.c', false, {}, '';
    'stitch/mex', 'integral_image_3d_mex.cpp', false, {}, '';
    'stitch/mex', 'any_4th_dim_mex.cpp', false, {}, '';
    'tools/MIP/mex', 'max_pooling_3d_mex.cpp', false, {}, '';
    'tools/MIP/mex', 'range_max_3d_mex.cpp', false, {}, '';
    'crop/mex', 'crop3d_mex.c', false, {}, '';
    'crop/mex', 'crop4d_mex.cpp', false, {}, '';
    'crop/mex', 'indexing3d_mex.c', false, {}, '';
    'crop/mex', 'indexing4d_mex.c', false, {}, '';
    'crop/mex', 'indexing4d_crop_mex.cpp', false, {}, '';
    'flatfield_correction/mex', 'flat_field_correction_3d_mex.cpp', false, {}, ''};

for i = 1 : size(mexFiles, 1)
    cd([srcDir, '/', mexFiles{i, 1}]);
    buildMex(mexFiles{i, 2}, mexFiles{i, 3}, mexFiles{i, 4}, mexFiles{i, 5});
end

% the decon mex files, on Linux only. The CPU mex uses interleaved complex
% arrays, MATLAB's single precision FFTW and no fused multiply-adds, so that
% it computes as MATLAB does; on macOS MATLAB's FFTs differ from it in the
% last bits, and Windows isn't checked. The GPU mex: mexcuda needs MSVC on
% Windows, and macOS has no CUDA.
if isunix && ~ismac
    cd([srcDir, '/decon/mex']);
    buildMex('decon_lucy_rfft_mex.cpp', false, [{'-R2018a'}, matlabFFTWArgs()], ' -ffp-contract=off');
    buildCudaMex('decon_lucy_rfft_gpu_mex.cu');
end

if ismac
    % the mex files must not load libraries from Homebrew
    for i = 1 : size(mexFiles, 1)
        [~, fsn] = fileparts(mexFiles{i, 2});
        mexFn = sprintf('%s/%s/%s.%s', srcDir, mexFiles{i, 1}, fsn, mexext);
        [~, deps] = system(sprintf('otool -L "%s"', mexFn));
        if contains(deps, '/opt/homebrew/')
            warning('%s still loads libraries from Homebrew:\n%s', mexFn, deps);
        end
    end
end

cd(baseDir);

end


function [] = buildMex(srcFn, useAVX, extraArgs, extraFlags)
% build one mex file in the current folder

[~, ~, ext] = fileparts(srcFn);
isCpp = strcmp(ext, '.cpp');
if ~ismac
    avxFlag = '';
    if useAVX
        avxFlag = ' -mavx2';
    end
    if isCpp
        mex('-v', 'CXXOPTIMFLAGS=-O3 -DNDEBUG', ['CXXFLAGS=$CXXFLAGS -O3 -fopenmp', avxFlag, extraFlags], ...
            ['LDFLAGS=$LDFLAGS -O3 -fopenmp', avxFlag], srcFn, extraArgs{:});
    else
        mex('-v', 'COPTIMFLAGS=-O3 -DNDEBUG', ['CFLAGS=$CFLAGS -O3 -fopenmp', avxFlag, extraFlags], ...
            ['LDFLAGS=$LDFLAGS -O3 -fopenmp', avxFlag], srcFn, extraArgs{:});
    end
    return;
end

% ARM macOS with Homebrew gcc@13; libstdc++ is linked statically, and with
% -static-libgcc GCC on macOS also links libgcc and libgomp statically
gccDir = '/opt/homebrew/opt/gcc@13/lib/gcc/13';
macFlags = ['-fno-common -arch arm64 -mmacosx-version-min=11.0 -fexceptions ', ...
    '-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -O3 -fopenmp ', ...
    '-DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD -DMATLAB_MEX_FILE', extraFlags];
if any(strcmp(extraArgs, '-R2018a'))
    macFlags = strrep(macFlags, 'R2017b', 'R2018a');
end
if isCpp
    mex('-v', 'CXX=/opt/homebrew/bin/g++-13', 'CXXOPTIMFLAGS=-O3 -DNDEBUG', ['CXXFLAGS=', macFlags, ' -std=c++11'], ...
        'LDFLAGS=$LDFLAGS -O3 -fopenmp -static-libgcc', [gccDir, '/libstdc++.a'], srcFn, extraArgs{:});
else
    mex('-v', 'CC=/opt/homebrew/bin/gcc-13', 'COPTIMFLAGS=-O3 -DNDEBUG', ['CFLAGS=', macFlags], ...
        'LDFLAGS=$LDFLAGS -O3 -fopenmp -static-libgcc', srcFn, extraArgs{:});
end

end


function [] = buildCudaMex(srcFn)
% build a CUDA mex file in the current folder with the nvcc of the MATLAB CUDA
% toolkit, linked to MATLAB's cuFFT. That nvcc (CUDA 12.2 in R2024b) only
% accepts GCC 12 or older as host compiler, so if g++ is newer, an older one is
% used (the system g++ or g++-12 ... g++-10).

binDir = [matlabroot, '/bin/', computer('arch')];
dir_info = dir([binDir, '/libcufft.so.*']);
if isempty(dir_info)
    error('MATLAB''s cuFFT library (libcufft.so) is not in %s', binDir);
end
fns = {dir_info.name};
[~, ind] = min(cellfun(@numel, fns));

args = {};
hostCompilers = {'g++', '/usr/bin/g++', 'g++-12', 'g++-11', 'g++-10'};
for i = 1 : numel(hostCompilers)
    [status, ver] = system([hostCompilers{i}, ' -dumpversion']);
    if status == 0 && str2double(strtok(strtrim(ver), '.')) <= 12
        if i > 1
            args = {['NVCC_FLAGS=-ccbin ', hostCompilers{i}]};
        end
        break;
    end
end
mexcuda('-v', '-R2018a', args{:}, srcFn, ['-L', binDir], ['-l:', fns{ind}]);

end


function [args] = matlabFFTWArgs()
% mex arguments to link MATLAB's single precision FFTW (libmwfftw3f in the
% MATLAB bin folder, the unversioned file if there is one; fftw3.h is in
% extern/include). MinGW links the dll directly.

binDir = [matlabroot, '/bin/', computer('arch')];
fns = {};
for pattern = {'*mwfftw3f*', '*fftw3f*'}
    dir_info = dir([binDir, '/', pattern{1}]);
    fns = {dir_info.name};
    fns = fns(~contains(fns, {'mpi', 'omp', 'threads'}) & ~endsWith(fns, {'.lib', '.exp', '.manifest', '.a'}));
    if ~isempty(fns)
        break;
    end
end
if isempty(fns)
    error('MATLAB''s single precision FFTW library (libmwfftw3f) is not in %s', binDir);
end
[~, ind] = min(cellfun(@numel, fns));
args = {['-L', binDir], ['-l', regexprep(fns{ind}, '^lib|\.(so|dylib|dll)(\..*)?$', '')]};

end
