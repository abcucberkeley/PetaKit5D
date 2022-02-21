sName = 'buildReader.sh';
filepath = fileparts(mfilename('fullpath'));
cd(filepath);
fileattrib([filepath '/' sName],'+x','a');
system([filepath '/' sName]);
mex('-v','COPTIMFLAGS="-O3 -fwrapv -DNDEBUG"', 'CFLAGS=$CFLAGS -O3 -fopenmp', 'LDFLAGS=$LDFLAGS -O3 -fopenmp', ['-I' filepath '/cBlosc2/include/'], ['-L' filepath '/cJSON/lib'] ,'-lblosc2',['-I' filepath '/cJSON/include/'], ['-L' filepath '/cJSON/lib'] ,'-lcjson' ,'parallelReadZarr.c');