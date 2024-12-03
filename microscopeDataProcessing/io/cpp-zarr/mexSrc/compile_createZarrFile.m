% ADD --disable-new-dtags
if isunix && ~ismac
    releaseFolder = '../linux';
    if ~exist(releaseFolder, 'dir')
        mkdir(releaseFolder);
    end
    mex -outdir ../linux -output createZarrFile.mexa64 -v CXXOPTIMFLAGS="-O2 -DNDEBUG" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O2 -DNDEBUG" CXXFLAGS='$CXXFLAGS -O2 -fopenmp' LDFLAGS='$LDFLAGS -O2 -fopenmp' -I'/clusterfs/fiona/matthewmueller/cppZarrTest/c-zarr/jenkinsBuild/install/include' -L'/clusterfs/fiona/matthewmueller/cppZarrTest/c-zarr/jenkinsBuild/install/lib64' -lcppZarr createzarrfilemex.cpp
    % Need to change the library name because matlab preloads their own version
    % of libstdc++
    % Setting it to libstdc++.so.6.0.30 as of MATLAB R2022b
    system(['patchelf --replace-needed libstdc++.so.6 libstdc++.so.6.0.32 ' releaseFolder '/createZarrFile.mexa64']);
    %system('export LD_LIBRARY_PATH="";patchelf --replace-needed libstdc++.so.6 libstdc++.so.6.0.32 createZarrFile.mexa64');
elseif ismac
    % Might have to do this part in terminal. First change the library
    % linked to libstdc++
    % system('install_name_tool -change @rpath/libgcc_s.1.1.dylib @loader_path/libgcc_s.1.1.0.dylib libstdc++.6.0.32.dylib');
    
    releaseFolder = '../mac';
    if ~exist(releaseFolder, 'dir')
        mkdir(releaseFolder);
    end
    mex -outdir ../mac -output createZarrFile.mexa64 -v CXX="/usr/local/bin/g++-13" CXXOPTIMFLAGS='-O2 -DNDEBUG' LDOPTIMFLAGS='-O2 -DNDEBUG' CXXFLAGS='-fno-common -arch arm64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O2 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' LDFLAGS='$LDFLAGS -O2 -fopenmp' -I'/Users/abcx86mac/c-zarr/jenkinsBuild/install/include/' -L'/Users/abcx86mac/c-zarr/jenkinsBuild/install/lib' /usr/local/opt/gcc/lib/gcc/current/libstdc++.a -lcppZarr createzarrfilemex.cpp

    % We need to change all the current paths to be relative to the mex file
    %system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libstdc++.6.dylib @loader_path/libstdc++.6.0.32.dylib ../mac/createZarrFile.mexmaci64');
    system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libgcc_s.1.1.dylib @loader_path/libgcc_s.1.1.0.dylib ../mac/createZarrFile.mexmaci64');
    system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libgomp.1.dylib @loader_path/libgomp.1.dylib ../mac/createZarrFile.mexmaci64');    
    system('install_name_tool -change @rpath/libcppZarr.dylib @loader_path/libcppZarr.dylib ../mac/createZarrFile.mexmaci64');

    system('chmod 777 ../mac/createZarrFile.mexmaci64');
elseif ispc
    setenv('MW_MINGW64_LOC','C:/mingw64');
    releaseFolder = '../windows';
    if ~exist(releaseFolder, 'dir')
        mkdir(releaseFolder);
    end
    mex -outdir ../windows -output createZarrFile.mexa64 -v CXX="C:/mingw64/bin/g++" CXXOPTIMFLAGS="-DNDEBUG -O2" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O2 -DNDEBUG" CXXFLAGS='$CXXFLAGS -fopenmp -O2' LDFLAGS='$LDFLAGS -fopenmp -O2' -I'C:/Users/matt/Documents/GitHub/c-zarr/jenkinsBuild/install/include' -L'C:/Users/matt/Documents/GitHub/c-zarr/jenkinsBuild/install/lib' C:\mingw64\lib\gcc\x86_64-w64-mingw32\12.2.0\libgcc_eh.a -lcppZarr.dll createzarrfilemex.cpp
end
