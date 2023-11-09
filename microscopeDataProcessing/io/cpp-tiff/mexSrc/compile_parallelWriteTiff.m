debug = false;
% ADD --disable-new-dtags
if isunix && ~ismac
    releaseFolder = '../linux';
    if ~exist(releaseFolder, 'dir')
        mkdir(releaseFolder);
    end
    if debug
        %mex -v -g CXXOPTIMFLAGS="" LDOPTIMFLAGS="-g -O0 -Wall -Wextra -Wl',-rpath='''$ORIGIN''''" CXXFLAGS='$CXXFLAGS -Wall -Wextra -g -O0 -fsanitize=address -gdwarf-4 -gstrict-dwarf' LDFLAGS='$LDFLAGS -g -O0 -fopenmp -fsanitize=address' -I'/clusterfs/fiona/matthewmueller/cppZarrTest' -I'/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.8.0/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.8.0/lib64' -lblosc2 -lz -luuid parallelreadzarr.cpp zarr.cpp helperfunctions.cpp parallelreadzarr.cpp
    else
        mex -outdir ../linux -output parallelWriteTiff.mexa64 -v CXXOPTIMFLAGS="-DNDEBUG -O3" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O3 -DNDEBUG" CXXFLAGS='$CXXFLAGS -fopenmp -O3' LDFLAGS='$LDFLAGS -fopenmp -O3' -I'/global/home/groups/software/sl-7.x86_64/modules/libtiff/4.6.0/include' -L'/global/home/groups/software/sl-7.x86_64/modules/libtiff/4.6.0/lib64' -lstdc++ -ltiff parallelwritetiffmex.cpp ../src/lzwencode.cpp ../src/helperfunctions.cpp ../src/parallelwritetiff.cpp
    end

    % Need to change the library name because matlab preloads their own version
    % of libstdc++
    % Setting it to libstdc++.so.6.0.32 as of MATLAB R2022b
    system(['patchelf --replace-needed libstdc++.so.6 libstdc++.so.6.0.32 ' releaseFolder '/parallelWriteTiff.mexa64']);
elseif ismac
    % Might have to do this part in terminal. First change the library
    % linked to libstdc++

    % Can make a custom libtiff dylib where the libaries are linked by
    % loader path
    
    releaseFolder = '../mac';
    if ~exist(releaseFolder, 'dir')
        mkdir(releaseFolder);
    end   
    mex -v -outdir ../mac -output parallelWriteTiff.mexw64 CXX="/usr/local/bin/g++-13" CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3 -DNDEBUG' CXXFLAGS='-fno-common -arch x86_64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O3 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/usr/local/include/' -lstdc++ -ltiff parallelwritetiffmex.cpp ../src/lzwencode.cpp ../src/helperfunctions.cpp ../src/parallelwritetiff.cpp

    % We need to change all the current paths to be relative to the mex file
    system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libstdc++.6.dylib @loader_path/libstdc++.6.0.32.dylib ../mac/parallelWriteTiff.mexmaci64');
    system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libgomp.1.dylib @loader_path/libgomp.1.dylib ../mac/parallelWriteTiff.mexmaci64');
    system('install_name_tool -change @rpath/libtiff.6.dylib @loader_path/libtiff.6.0.1.dylib ../mac/parallelWriteTiff.mexmaci64');
    system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libgcc_s.1.1.dylib @loader_path/libgcc_s.1.1.dylib ../mac/parallelWriteTiff.mexmaci64'); 

    system('chmod 777 ../mac/parallelWriteTiff.mexmaci64');
elseif ispc
    setenv('MW_MINGW64_LOC','C:/mingw64');
    releaseFolder = '../windows';
    if ~exist(releaseFolder, 'dir')
        mkdir(releaseFolder);
    end
    mex  -outdir ../windows -output parallelWriteTiff.mexw64 -v CXX="C:/mingw64/bin/g++" -v CXXOPTIMFLAGS="-DNDEBUG -O3" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O3 -DNDEBUG" CXXFLAGS='$CXXFLAGS -fopenmp -O3' LDFLAGS='$LDFLAGS -fopenmp -O3' -I'C:/Program Files (x86)/tiff/include' -L'C:\Program Files (x86)\tiff\lib' -ltiff.dll parallelwritetiffmex.cpp ../src/lzwencode.cpp ../src/helperfunctions.cpp ../src/parallelwritetiff.cpp
end
