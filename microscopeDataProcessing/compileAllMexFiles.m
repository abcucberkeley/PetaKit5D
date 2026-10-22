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

% folder, source file, and whether it uses AVX2 on x86
mexFiles = {'deskew_rotate/mex', 'skewed_space_interp_defined_stepsize_mex.c', true;
    'deskew_rotate/mex', 'skewed_space_interp_mex.c', false;
    'deskew_rotate/mex', 'volume_deskew_rotate_warp_mex.cpp', true;
    'deskew_rotate/mex', 'skewed_space_interp_volume_deskew_rotate_warp_mex.cpp', true;
    'stitch/mex', 'feather_blending_3d_mex.cpp', false;
    'stitch/mex', 'feather_blending_3d_with_indexing_mex.cpp', false;
    'stitch/mex', 'feather_distance_map_resize_3d_mex.c', false;
    'stitch/mex', 'integral_image_3d_mex.cpp', false;
    'stitch/mex', 'any_4th_dim_mex.cpp', false;
    'tools/MIP/mex', 'max_pooling_3d_mex.cpp', false;
    'crop/mex', 'crop3d_mex.c', false;
    'crop/mex', 'crop4d_mex.cpp', false;
    'crop/mex', 'indexing3d_mex.c', false;
    'crop/mex', 'indexing4d_mex.c', false;
    'crop/mex', 'indexing4d_crop_mex.cpp', false;
    'flatfield_correction/mex', 'flat_field_correction_3d_mex.cpp', false};

for i = 1 : size(mexFiles, 1)
    cd([srcDir, '/', mexFiles{i, 1}]);
    buildMex(mexFiles{i, 2}, mexFiles{i, 3});
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


function [] = buildMex(srcFn, useAVX)
% build one mex file in the current folder

[~, ~, ext] = fileparts(srcFn);
isCpp = strcmp(ext, '.cpp');
if ~ismac
    avxFlag = '';
    if useAVX
        avxFlag = ' -mavx2';
    end
    if isCpp
        mex('-v', 'CXXOPTIMFLAGS=-O3 -DNDEBUG', ['CXXFLAGS=$CXXFLAGS -O3 -fopenmp', avxFlag], ...
            ['LDFLAGS=$LDFLAGS -O3 -fopenmp', avxFlag], srcFn);
    else
        mex('-v', 'COPTIMFLAGS=-O3 -DNDEBUG', ['CFLAGS=$CFLAGS -O3 -fopenmp', avxFlag], ...
            ['LDFLAGS=$LDFLAGS -O3 -fopenmp', avxFlag], srcFn);
    end
    return;
end

% ARM macOS with Homebrew gcc@13; libstdc++ is linked statically, and with
% -static-libgcc GCC on macOS also links libgcc and libgomp statically
gccDir = '/opt/homebrew/opt/gcc@13/lib/gcc/13';
macFlags = ['-fno-common -arch arm64 -mmacosx-version-min=11.0 -fexceptions ', ...
    '-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -O3 -fopenmp ', ...
    '-DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD -DMATLAB_MEX_FILE'];
if isCpp
    mex('-v', 'CXX=/opt/homebrew/bin/g++-13', 'CXXOPTIMFLAGS=-O3 -DNDEBUG', ['CXXFLAGS=', macFlags, ' -std=c++11'], ...
        'LDFLAGS=$LDFLAGS -O3 -fopenmp -static-libgcc', [gccDir, '/libstdc++.a'], srcFn);
else
    mex('-v', 'CC=/opt/homebrew/bin/gcc-13', 'COPTIMFLAGS=-O3 -DNDEBUG', ['CFLAGS=', macFlags], ...
        'LDFLAGS=$LDFLAGS -O3 -fopenmp -static-libgcc', srcFn);
end

end
