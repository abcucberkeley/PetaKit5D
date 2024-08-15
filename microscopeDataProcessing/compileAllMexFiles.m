function [] = compileAllMexFiles(varargin)

% Intended to be run from PetaKit5D/microscopeDataProcessing or pass the
% path to the PetaKit5D base dir
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

cd('./deskew_rotate/mex');

% skewed_space_interp_defined_stepsize_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_defined_stepsize_mex.c
else
    mex -v CC="/usr/local/bin/gcc-13" CXX="/usr/local/bin/g++-13" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_defined_stepsize_mex.c
end

% skewed_space_interp_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_mex.c
else
    mex -v CC="/usr/local/bin/gcc-13" CXX="/usr/local/bin/g++-13" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_defined_stepsize_mex.c
end

% volume_deskew_rotate_warp_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' volume_deskew_rotate_warp_mex.cpp
else
    mex -v CC="/usr/local/bin/gcc-13" CXX="/usr/local/bin/g++-13" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' CXXFLAGS='-fno-common -arch arm64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O3 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' LDFLAGS='$LDFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_blending_3d_with_indexing_mex.cpp
end


cd('../../stitch/mex');

% feather_blending_3d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_blending_3d_mex.cpp
else
    mex -v CC="/usr/local/bin/gcc-13" CXX="/usr/local/bin/g++-13" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' CXXFLAGS='-fno-common -arch arm64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O3 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' LDFLAGS='$LDFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_blending_3d_mex.cpp
end

% feather_blending_3d_with_indexing_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_blending_3d_with_indexing_mex.cpp
else
    mex -v CC="/usr/local/bin/gcc-13" CXX="/usr/local/bin/g++-13" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' CXXFLAGS='-fno-common -arch arm64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O3 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' LDFLAGS='$LDFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_blending_3d_with_indexing_mex.cpp
end

% feather_distance_map_resize_3d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_distance_map_resize_3d_mex.c
else
    mex -v CC="/usr/local/bin/gcc-13" CXX="/usr/local/bin/g++-13" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' CXXFLAGS='-fno-common -arch arm64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O3 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_distance_map_resize_3d_mex.c
end

% integral_image_3d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' integral_image_3d_mex.cpp
else
    mex -v CC="/usr/local/bin/gcc-13" CXX="/usr/local/bin/g++-13" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' CXXFLAGS='-fno-common -arch arm64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O3 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' LDFLAGS='$LDFLAGS -O3 -fopenmp' integral_image_3d_mex.cpp
end

% any_4th_dim_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' any_4th_dim_mex.cpp
else
    mex -v CC="/usr/local/bin/gcc-13" CXX="/usr/local/bin/g++-13" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' CXXFLAGS='-fno-common -arch arm64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O3 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' LDFLAGS='$LDFLAGS -O3 -fopenmp' any_4th_dim_mex.cpp
end

cd('../../tools/MIP/mex/');

% max_pooling_3d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' max_pooling_3d_mex.cpp
else
    mex -v CC="/usr/local/bin/gcc-13" CXX="/usr/local/bin/g++-13" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' CXXFLAGS='-fno-common -arch arm64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O3 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' LDFLAGS='$LDFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' max_pooling_3d_mex.cpp
end

cd('../../../crop/mex');

% crop3d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' crop3d_mex.c
else
    mex -v CC="/usr/local/bin/gcc-13" CXX="/usr/local/bin/g++-13" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' crop3d_mex.c
end

% crop4d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' crop4d_mex.cpp
else
    mex -v CXX="/usr/local/bin/g++-13" CXXFLAGS='-fno-common -arch x86_64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O3 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' CXXOPTIMFLAGS="-O3 -DNDEBUG" LDFLAGS='$LDFLAGS -O3 -fopenmp' crop4d_mex.cpp
end

% indexing3d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing3d_mex.c
else
    mex -v CC="/usr/local/bin/gcc-13" CXX="/usr/local/bin/g++-13" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing3d_mex.c
end

% indexing4d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing4d_mex.c
else
    mex -v CC="/usr/local/bin/gcc-13" CXX="/usr/local/bin/g++-13" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing4d_mex.c
end

% indexing4d_crop_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing4d_crop_mex.cpp
else
    mex -v CXX="/usr/local/bin/g++-13" CXXFLAGS='-fno-common -arch x86_64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O3 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' CXXOPTIMFLAGS="-O3 -DNDEBUG" LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing4d_crop_mex.cpp
end

cd('../../flatfield_correction/mex');

% flat_field_correction_3d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' flat_field_correction_3d_mex.cpp
else
    mex -v CXX="/usr/local/bin/g++-13" CXXFLAGS='-fno-common -arch x86_64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O3 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' CXXOPTIMFLAGS="-O3 -DNDEBUG" LDFLAGS='$LDFLAGS -O3 -fopenmp' flat_field_correction_3d_mex.cpp
end

cd(baseDir);

end
