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
            mex -outdir ../linux -output parallelWriteZarr.mexa64 -g -v CXXOPTIMFLAGS="" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O0 -g" CXXFLAGS='$CXXFLAGS -fsanitize=address -fopenmp -O0 -g' LDFLAGS='$LDFLAGS -fsanitize=address -fopenmp -O0 -g' -I'/clusterfs/fiona/matthewmueller/cppZarrTest' -I'/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.8.0/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.8.0/lib64' -I'/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/lib' -lblosc -lblosc2 -lz -luuid parallelwritezarrmex.cpp ../src/parallelwritezarr.cpp ../src/helperfunctions.cpp ../src/zarr.cpp ../src/parallelreadzarr.cpp
        else
            mex -outdir ../linux -output parallelWriteZarr.mexa64 -v CXXOPTIMFLAGS="" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O0 -g" CXXFLAGS='$CXXFLAGS -fopenmp -O0 -g' LDFLAGS='$LDFLAGS -fopenmp -O0 -g' -I'/clusterfs/fiona/matthewmueller/cppZarrTest' -I'/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.8.0/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.8.0/lib64' -I'/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/lib' -lblosc -lblosc2 -lz -luuid parallelwritezarrmex.cpp ../src/parallelwritezarr.cpp ../src/helperfunctions.cpp ../src/zarr.cpp ../src/parallelreadzarr.cpp
        end
    else
        mex -outdir ../linux -output parallelWriteZarr.mexa64 -v CXXOPTIMFLAGS="-DNDEBUG -O2" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O2 -DNDEBUG" CXXFLAGS='$CXXFLAGS -fopenmp -O2' LDFLAGS='$LDFLAGS -fopenmp -O2' -I'/clusterfs/fiona/matthewmueller/cppZarrTest' -I'/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.8.0/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.8.0/lib64' -I'/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/lib' -lblosc -lblosc2 -lz -luuid parallelwritezarrmex.cpp ../src/parallelwritezarr.cpp ../src/helperfunctions.cpp ../src/zarr.cpp ../src/parallelreadzarr.cpp
    end
    % Need to change the library name because matlab preloads their own version
    % of libstdc++
    % Setting it to libstdc++.so.6.0.30 as of MATLAB R2022b
    system(['patchelf --replace-needed libstdc++.so.6 libstdc++.so.6.0.32 ' releaseFolder '/parallelWriteZarr.mexa64']);
    %system('export LD_LIBRARY_PATH="";patchelf --replace-needed libstdc++.so.6 libstdc++.so.6.0.32 parallelWriteZarr.mexa64');
elseif ismac
    % Might have to do this part in terminal. First change the library
    % linked to libstdc++
    % system('install_name_tool -change @rpath/libgcc_s.1.1.dylib @loader_path/libgcc_s.1.1.0.dylib libstdc++.6.0.32.dylib');
    releaseFolder = '../mac';
    if ~exist(releaseFolder, 'dir')
        mkdir(releaseFolder);
    end

    if debug
        %mex -v -g CXX="/usr/local/opt/llvm/bin/clang++" CXXOPTIMFLAGS="" LDOPTIMFLAGS='-g -O0 -Wall -Wextra' CXXFLAGS='$CXXFLAGS -g -O0 -Wall -Wextra -fopenmp' LDFLAGS='$LDFLAGS -g -O0 -fopenmp' '-I/usr/local/include/' -L'/usr/local/lib/' -lblosc -lblosc2 -lz -luuid parallelwritezarr.cpp helperfunctions.cpp zarr.cpp parallelwritezarrread.cpp
        %mex -v -g CXXOPTIMFLAGS='' LDOPTIMFLAGS='-g -O0' CXXFLAGS='$CXXFLAGS -O0 -Xpreprocessor -fopenmp' LDFLAGS='$LDFLAGS -g -O0 -Xpreprocessor -fopenmp' '-I/usr/local/opt/libomp/include' '-I/usr/include' '-I/usr/local/include/' -L'/usr/local/lib/' -L'/usr/local/opt/libomp/lib' -lomp -lblosc -lblosc2 -lz -luuid parallelwritezarr.cpp helperfunctions.cpp zarr.cpp parallelwritezarrread.cpp
    else
        mex -outdir ../mac -output parallelWriteZarr.mexa64 -v CXX="/usr/local/bin/g++-13" CXXOPTIMFLAGS='-O2 -DNDEBUG' LDOPTIMFLAGS='-O2 -DNDEBUG' CXXFLAGS='-fno-common -arch x86_64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O2 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' LDFLAGS='$LDFLAGS -O2 -fopenmp' '-I/usr/local/include/' -lstdc++ -lblosc -lblosc2 -lz -luuid parallelwritezarrmex.cpp ../src/parallelwritezarr.cpp ../src/helperfunctions.cpp ../src/zarr.cpp ../src/parallelreadzarr.cpp
    end
    
    % We need to change all the current paths to be relative to the mex file
    system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libstdc++.6.dylib @loader_path/libstdc++.6.0.32.dylib ../mac/parallelWriteZarr.mexmaci64');
    system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libgcc_s.1.1.dylib @loader_path/libgcc_s.1.1.0.dylib ../mac/parallelWriteZarr.mexmaci64');
    system('install_name_tool -change /usr/local/opt/ossp-uuid/lib/libuuid.16.dylib @loader_path/libuuid.16.22.0.dylib ../mac/parallelWriteZarr.mexmaci64');
    system('install_name_tool -change /usr/local/opt/gcc/lib/gcc/current/libgomp.1.dylib @loader_path/libgomp.1.dylib ../mac/parallelWriteZarr.mexmaci64');    
    system('install_name_tool -change @rpath/libblosc.1.dylib @loader_path/libblosc.1.21.0.0.dylib ../mac/parallelWriteZarr.mexmaci64');
    system('install_name_tool -change @rpath/libblosc2.2.dylib @loader_path/libblosc2.2.8.0.dylib ../mac/parallelWriteZarr.mexmaci64');
    system('install_name_tool -change /usr/lib/libz.1.dylib @loader_path/libz.1.2.13.dylib ../mac/parallelWriteZarr.mexmaci64');    

    system('chmod 777 ../mac/parallelWriteZarr.mexmaci64');
elseif ispc
    setenv('MW_MINGW64_LOC','C:/mingw64');
    releaseFolder = '../windows';
    if ~exist(releaseFolder, 'dir')
        mkdir(releaseFolder);
    end
    mex -outdir ../windows -output parallelWriteZarr.mexa64 -v CXX="C:/mingw64/bin/g++" CXXOPTIMFLAGS="-DNDEBUG -O2" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O2 -DNDEBUG" CXXFLAGS='$CXXFLAGS -fopenmp -O2' LDFLAGS='$LDFLAGS -fopenmp -O2' -I'C:/Program Files (x86)/nlohmann_json/include' -I'C:/Program Files (x86)/blosc2/include' '-LC:/Program Files (x86)/blosc2/lib' -I'C:/Program Files (x86)/blosc/include' '-LC:/Program Files (x86)/blosc/lib' -I'C:/Program Files (x86)/zlib/include' '-LC:/Program Files (x86)/zlib/lib' -lblosc.dll -lblosc2.dll -lzlib.dll parallelwritezarrmex.cpp ../src/parallelwritezarr.cpp ../src/helperfunctions.cpp ../src/zarr.cpp ../src/parallelreadzarr.cpp
end