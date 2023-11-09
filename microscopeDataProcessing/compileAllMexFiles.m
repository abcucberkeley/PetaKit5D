function [] = compileAllMexFiles(varargin)

% Intended to be run from LLSM5DTools/microscopeDataProcessing or pass the
% path to the LLSM5DTools base dir
ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('LLSM5DToolsBaseDir', '', @ischar);
ip.parse(varargin{:});

LLSM5DToolsBaseDir = ip.Results.LLSM5DToolsBaseDir;

baseDir = pwd;

if(~isempty(LLSM5DToolsBaseDir))
    cd([LLSM5DToolsBaseDir '/./microscopeDataProcessing']);
end
    
[~,name,~]=fileparts(pwd);
if ~strcmp(name,'microscopeDataProcessing')
    error('Pass the path to the LLSM5DTools base dir or be in the LLSM5DTools/microscopeDataProcessing dir');
end

cd('./deskew_rotate/mex');

% skewed_space_interp_defined_stepsize_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_defined_stepsize_mex.c
else
    mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_defined_stepsize_mex.c
end

% skewed_space_interp_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_mex.c
else
    mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_defined_stepsize_mex.c
end

cd('../../stitch/mex');

% feather_blending_3d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_blending_3d_mex.cpp
else
    mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_blending_3d_mex.cpp
end

% feather_distance_map_resize_3d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_distance_map_resize_3d_mex.c
else
    mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_distance_map_resize_3d_mex.c
end

% max_pooling_3d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' max_pooling_3d_mex.cpp
else
    mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' max_pooling_3d_mex.cpp
end

cd('../../crop/mex');

% crop3d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' crop3d_mex.c
else
    mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' crop3d_mex.c
end

% indexing3d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing3d_mex.c
else
    mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing3d_mex.c
end

% indexing4d_mex
if(~ismac)
    mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing4d_mex.c
else
    mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing4d_mex.c
end

cd(baseDir);

end
