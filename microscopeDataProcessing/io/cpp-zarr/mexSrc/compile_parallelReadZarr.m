debug = false;
% ADD --disable-new-dtags
if isunix && ~ismac
    releaseFolder = '../linux';
    if ~exist(releaseFolder, 'dir')
        mkdir(releaseFolder);
    end
    if debug
        fsanitize = false;
        if fsanitize
            mex -outdir ../linux -output parallelReadZarr.mexa64 -v CXXOPTIMFLAGS="" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O0 -g" CXXFLAGS='$CXXFLAGS -fsanitize=address -fopenmp -O0 -g' LDFLAGS='$LDFLAGS -fsanitize=address -fopenmp -O0 -g' -I'/clusterfs/fiona/matthewmueller/cppZarrTest' -I'/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.8.0/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.8.0/lib64' -lblosc2 -lz -luuid parallelreadzarrmex.cpp ../src/zarr.cpp ../src/helperfunctions.cpp ../src/parallelreadzarr.cpp
        else
            mex -outdir ../linux -output parallelReadZarr.mexa64 -v CXXOPTIMFLAGS="" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O0 -g" CXXFLAGS='$CXXFLAGS -fopenmp -O0 -g' LDFLAGS='$LDFLAGS -fopenmp -O0 -g' -I'/clusterfs/fiona/matthewmueller/cppZarrTest' -I'/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.8.0/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.8.0/lib64' -lblosc2 -lz -luuid parallelreadzarrmex.cpp ../src/zarr.cpp ../src/helperfunctions.cpp ../src/parallelreadzarr.cpp
        end
    else
        mex -outdir ../linux -output parallelReadZarr.mexa64 -v CXXOPTIMFLAGS="-DNDEBUG -O2" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O2 -DNDEBUG" CXXFLAGS='$CXXFLAGS -fopenmp -O2' LDFLAGS='$LDFLAGS -fopenmp -O2' -I'/clusterfs/fiona/matthewmueller/cppZarrTest/c-zarr/jenkinsBuild/install/include' -L'/clusterfs/fiona/matthewmueller/cppZarrTest/c-zarr/jenkinsBuild/install/lib64' -lcppZarr parallelreadzarrmex.cpp
    end
    % Need to change the library name because matlab preloads their own version
    % of libstdc++
    % Setting it to libstdc++.so.6.0.32 as of MATLAB R2023a with gcc 13.2
    %system('patchelf --replace-needed libc.so.6 libc.so.6.0.0 parallelreadzarr.mexa64');
    system(['patchelf --replace-needed libstdc++.so.6 libstdc++.so.6.0.32 ' releaseFolder '/parallelReadZarr.mexa64']);
    %system('export LD_LIBRARY_PATH="";patchelf --replace-needed libstdc++.so.6 libstdc++.so.6.0.32 parallelReadZarr.mexa64');
elseif ismac
    % Might have to do this part in terminal. First change the library
    % linked to libstdc++
    % system('install_name_tool -change @rpath/libgcc_s.1.1.dylib @loader_path/libgcc_s.1.1.0.dylib libstdc++.6.0.32.dylib');
    if computer == "MACI64"
        releaseFolder = '../mac';
        if ~exist(releaseFolder, 'dir')
            mkdir(releaseFolder);
        end
        mex -outdir ../mac -output parallelReadZarr.mexmaci64 -v CXX="/usr/local/bin/g++-13" CXXOPTIMFLAGS='-O2 -DNDEBUG' LDOPTIMFLAGS='-O2 -DNDEBUG' CXXFLAGS='-fno-common -arch x86_64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O2 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' LDFLAGS='$LDFLAGS -O2 -fopenmp' -I'/Users/abcx86mac/c-zarr/jenkinsBuild/install/include' -L'/Users/abcx86mac/c-zarr/jenkinsBuild/install/lib' /usr/local/opt/gcc/lib/gcc/current/libstdc++.a -lcppZarr parallelreadzarrmex.cpp
    
        % We need to change all the current paths to be relative to the mex file
        %system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libstdc++.6.dylib @loader_path/libstdc++.6.0.32.dylib ../mac/parallelReadZarr.mexmaci64');
        system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libgcc_s.1.1.dylib @loader_path/libgcc_s.1.1.0.dylib ../mac/parallelReadZarr.mexmaci64');
        system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libgomp.1.dylib @loader_path/libgomp.1.dylib ../mac/parallelReadZarr.mexmaci64');    
        system('install_name_tool -change @rpath/libcppZarr.dylib @loader_path/libcppZarr.dylib ../mac/parallelReadZarr.mexmaci64');
    
        system('chmod 777 ../mac/parallelReadZarr.mexmaci64');
    else
        releaseFolder = '../macArm';
        if ~exist(releaseFolder, 'dir')
            mkdir(releaseFolder);
        end
        mex -outdir ../macArm -output parallelReadZarr.mexmaca64 -v CXX="/opt/homebrew/bin/g++-13" CXXOPTIMFLAGS='-O2 -DNDEBUG' LDOPTIMFLAGS='-O2 -DNDEBUG' CXXFLAGS='-fno-common -arch x86_64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O2 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' LDFLAGS='$LDFLAGS -O2 -fopenmp' -I'/Users/abcarmmac/cpp-zarr/jenkinsBuild/install/include' -L'/Users/abcarmmac/cpp-zarr/jenkinsBuild/install/lib' /opt/homebrew/opt/gcc@13/lib/gcc/13/libstdc++.a -lcppZarr parallelreadzarrmex.cpp
    
        % We need to change all the current paths to be relative to the mex file
        %system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libstdc++.6.dylib @loader_path/libstdc++.6.0.32.dylib ../mac/parallelReadZarr.mexmaci64');
        system('install_name_tool -change /opt/homebrew/opt/gcc@13/lib/gcc/13/libgcc_s.1.1.dylib @loader_path/libgcc_s.1.1.0.dylib ../macArm/parallelReadZarr.mexmaca64');
        system('install_name_tool -change /opt/homebrew/opt/gcc@13/lib/gcc/13/libgomp.1.dylib @loader_path/libgomp.1.dylib ../macArm/parallelReadZarr.mexmaca64');    
        system('install_name_tool -change @rpath/libcppZarr.dylib @loader_path/libcppZarr.dylib ../macArm/parallelReadZarr.mexmaca64');
    
        system('chmod 777 ../macArm/parallelReadZarr.mexmaca64');
    end
elseif ispc
    setenv('MW_MINGW64_LOC','C:/mingw64');
    releaseFolder = '../windows';
    if ~exist(releaseFolder, 'dir')
        mkdir(releaseFolder);
    end
    mex -outdir ../windows -output parallelReadZarr.mexa64 -v CXX="C:/mingw64/bin/g++" CXXOPTIMFLAGS="-DNDEBUG -O2" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O2 -DNDEBUG" CXXFLAGS='$CXXFLAGS -fopenmp -O2' LDFLAGS='$LDFLAGS -fopenmp -O2' -I'C:/Users/matt/Documents/GitHub/c-zarr/jenkinsBuild/install/include' -L'C:/Users/matt/Documents/GitHub/c-zarr/jenkinsBuild/install/lib' C:\mingw64\lib\gcc\x86_64-w64-mingw32\12.2.0\libgcc_eh.a -lcppZarr.dll parallelreadzarrmex.cpp
end
