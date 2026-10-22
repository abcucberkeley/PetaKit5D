#include <cstdint>
#include <cstring>
#include "mex.h"
#include "../src/helperfunctions.h"
#include "../src/parallelwritezarr.h"
#include "../src/parallelreadzarr.h"
#include "../src/zarr.h"
#include "zarrmexhelpers.h"

//compile
//mex -v COPTIMFLAGS="-DNDEBUG -O3" CFLAGS='$CFLAGS -fopenmp -O3' LDFLAGS='$LDFLAGS -fopenmp -O3' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/include/' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/include/' '-I/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/lib' -lblosc '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/lib64' -lblosc2 '-L/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/lib64' -lcjson -luuid parallelWriteZarr.c helperFunctions.c parallelReadZarr.c

//With zlib
//mex -v COPTIMFLAGS="-DNDEBUG -O3" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O3 -DNDEBUG" CFLAGS='$CFLAGS -fopenmp -O3' LDFLAGS='$LDFLAGS -fopenmp -O3' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/include/' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/include/' '-I/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/zarr/lib' -lblosc '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/lib64' -lblosc2 '-L/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/lib64' -lcjson -luuid -lz parallelWriteZarr.c helperFunctions.c parallelReadZarr.c

//mex -v COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/lib64' -lblosc2 zarrMex.c
//
//Windows
//mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-IC:\Program Files (x86)\bloscZarr\include' '-LC:\Program Files (x86)\bloscZarr\lib' -lblosc '-IC:\Program Files (x86)\cJSON\include\' '-LC:\Program Files (x86)\cJSON\lib' -lcjson '-IC:\Program Files (x86)\blosc\include' '-LC:\Program Files (x86)\blosc\lib' -lblosc2 parallelWriteZarr.c parallelReadZarr.c helperFunctions.c

// TODO: FIX MEMORY LEAKS
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
#if defined(__APPLE__)
    // keep the mex file loaded until MATLAB exits; unloading it crashes MATLAB on exit
    if (!mexIsLocked()) {
        mexLock();
    }
#endif
    if(nrhs < 2) mexErrMsgIdAndTxt("zarr:inputError","This functions requires at least 2 arguments");
    if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("zarr:inputError","The first argument must be a string");
    if(mxIsEmpty(prhs[1])) mexErrMsgIdAndTxt("zarr:inputError","All input data axes must be of at least size 1");

    std::vector<uint64_t> startCoords;
    std::vector<uint64_t> endCoords;
    bool crop = false;
    bool useUuid = true;
    bool sparse = true;
    // Per-axis options are collected first and fit to the array's number of
    // dimensions once it is known
    std::vector<uint64_t> bboxVals, chunksVals, subfoldersVals, chunkShapeVals;
    bool hasChunks = false, hasSubfolders = false, hasChunkShape = false;
    // Applied after the Zarr format, which sets their Zarr v3 defaults
    uint64_t zarrFormat = 0;
    std::string order, dimension_separator;
    // The data type of a new array (by default the data's)
    std::string dtypeOption;

    // Input data dimensions (MATLAB arrays have at least 2)
    const uint64_t nDataDims = (uint64_t)mxGetNumberOfDimensions(prhs[1]);
    const mwSize* dataDimsT = mxGetDimensions(prhs[1]);
    const std::vector<uint64_t> dataDims(dataDimsT, dataDimsT+nDataDims);

    std::string folderName(mxArrayToString(prhs[0]));
    // Handle the tilde character in filenames on Linux/Mac
    #ifndef _WIN32
    folderName = expandTilde(folderName.c_str());
    #endif
    
    // Check if metadata exists that we can use or if we have to create new metadata
    zarr Zarr;
    const bool metadataExists = fileExists(folderName+"/.zarray") || fileExists(folderName+"/zarr.json");
    if(!metadataExists){
        Zarr.set_fileName(folderName);
    }
    else{
        try{
            Zarr = zarr(folderName);
        }
        catch(const std::string &e){
            if(e.find("metadataFileMissing") != std::string::npos){
                mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s for writing. Try checking permissions or the file path.\n",e.substr(e.find(':')+1).c_str());
            }
            else if(e == "metadataIncomplete"){
                mexErrMsgIdAndTxt("zarr:zarrayError","Metadata is incomplete. Check the .zarray or zarr.json file");
            }
            else mexZarrError(e);
        }
    }

    for(int i = 2; i < nrhs; i+=2){
        if(i+1 == nrhs) mexErrMsgIdAndTxt("zarr:inputError","Mismatched argument pair for input number %d\n",i+1);
        if(!mxIsChar(prhs[i])) mexErrMsgIdAndTxt("zarr:inputError","The argument in input location %d is not a string\n",i+1);
        std::string currInput = mxArrayToString(prhs[i]);

        if(currInput == "uuid"){
            useUuid = (bool)*((mxGetPr(prhs[i+1])));
        }
        else if(currInput == "bbox"){
            // Skip bbox if it is empty
            if(!mxGetN(prhs[i+1])) continue;
            // A region [starts ends], or one value per axis for the chunk size
            bboxVals = mexVector(prhs[i+1]);
        }
        else if(currInput == "cname"){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","cname must be a string\n");
            Zarr.set_cname(mxArrayToString(prhs[i+1]));
        }
        // (like cname, used for a new array; a region written into an existing array
        // keeps the array's)
        else if(currInput == "clevel"){
            Zarr.set_clevel(mexClevel(prhs[i+1]));
        }
        else if(currInput == "subfolders"){
            subfoldersVals = mexVector(prhs[i+1]);
            hasSubfolders = true;
        }
        else if(currInput == "chunks"){
            chunksVals = mexVector(prhs[i+1]);
            hasChunks = true;
        }
        else if(currInput == "chunk_shape"){
            chunkShapeVals = mexVector(prhs[i+1]);
            hasChunkShape = true;
        }
        else if(currInput == "sparse"){
            sparse = (bool)*((mxGetPr(prhs[i+1])));
        }
        else if(currInput == "dimension_separator"){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","dimension_separator must be a string\n");
            dimension_separator = mxArrayToString(prhs[i+1]);
            if(dimension_separator != "." && dimension_separator != "/") mexErrMsgIdAndTxt("zarr:inputError","dimension_separator must be a . or /\n");
        }
        else if(currInput == "order"){
            order = mexOrder(prhs[i+1]);
        }
        else if(currInput == "zarr_format"){
            zarrFormat = mexZarrFormat(prhs[i+1]);
        }
        else if(currInput == "dtype"){
            dtypeOption = mexDtype(prhs[i+1]);
        }
        else{
            mexErrMsgIdAndTxt("zarr:inputError","The argument \"%s\" does not match the name of any supported input name.\n \
            Currently Supported Names: uuid, bbox, cname, clevel, subfolders, chunks, chunk_shape, sparse, dimension_separator, order, zarr_format, dtype\n",currInput.c_str());
        }
    }

    // The Zarr format: as given, or the existing array's (2 for a new array). A new
    // format (a new array, or an existing one rewritten in the other format) starts
    // from that format's defaults for order and dimension_separator.
    const uint64_t existingFormat = metadataExists ? Zarr.get_zarr_format() : 0;
    if(zarrFormat && zarrFormat != existingFormat) Zarr.set_zarr_format(zarrFormat);
    if(!order.empty()) Zarr.set_order(order);
    if(!dimension_separator.empty()) Zarr.set_dimension_separator(dimension_separator);

    // Number of dimensions of the array being written: a region written into an
    // existing array uses that array's; otherwise the input data's
    uint64_t nDims = nDataDims;
    std::vector<uint64_t> chunksFromBbox;
    if(!bboxVals.empty()){
        uint64_t nRegion = nDataDims;
        std::vector<uint64_t> fileShape;
        if(metadataExists){
            nRegion = Zarr.get_ndims();
            for(uint64_t d = 0; d < nRegion; d++) fileShape.push_back(Zarr.get_shape(d));
        }
        if(mexParseRegion(bboxVals, nRegion, startCoords, endCoords, fileShape)){
            crop = true;
            nDims = nRegion;
        }
        else if(bboxVals.size() == nDataDims || (bboxVals.size() == 3 && nDataDims < 3)){
            // One value per axis: the chunk size
            chunksFromBbox = mexFitAxes(bboxVals, nDataDims, "bbox");
        }
        else{
            mexErrMsgIdAndTxt("zarr:inputError","bbox must have %d values ([starts ends]) or %d values (chunk size)",
                              (int)(2*nRegion),(int)nDataDims);
        }
    }

    // The data may omit trailing singleton axes (MATLAB drops them)
    std::vector<uint64_t> iDims(nDims,1);
    for(uint64_t d = 0; d < nDataDims; d++){
        if(d < nDims) iDims[d] = dataDims[d];
        else if(dataDims[d] != 1) mexErrMsgIdAndTxt("zarr:inputError","Bounds are invalid for the input data size");
    }
    if(crop){
        for(uint64_t d = 0; d < nDims; d++){
            if(endCoords[d]-startCoords[d] > iDims[d]) mexErrMsgIdAndTxt("zarr:inputError","Bounds are invalid for the input data size");
        }
        // Writing a region into an existing array keeps its format
        if(metadataExists && zarrFormat && zarrFormat != existingFormat){
            mexErrMsgIdAndTxt("zarr:inputError","%s is a Zarr v%d array; writing a region keeps its format%s\n",folderName.c_str(),
                              (int)existingFormat,existingFormat == 2 ? " (convert it with convertZarrToV3 first)" : "");
        }
    }

    if(hasSubfolders) Zarr.set_subfolders(mexFitAxes(subfoldersVals, nDims, "subfolders"));
    if(hasChunks) Zarr.set_chunks(mexFitAxes(chunksVals, nDims, "chunks"));
    if(hasChunkShape){
        Zarr.set_shard(true);
        Zarr.set_chunk_shape(mexFitAxes(chunkShapeVals, nDims, "chunk_shape"));
    }

    void* zarrC = NULL;

    // The data's own type (logical is bool, complex single and double are complex64
    // and complex128), which a new array gets unless the dtype option says otherwise
    const std::string dataDtype = mexClassDtype(prhs[1]);
    if(dataDtype.empty()) mexErrMsgIdAndTxt("zarr:dataTypeError","The input data type is not supported");
    const bool cropExisting = crop && metadataExists;
    if(!cropExisting) Zarr.set_dtype(dtypeOption.empty() ? dataDtype : dtypeOption);

    if(!crop){
        Zarr.set_shape(iDims);
        if(!chunksFromBbox.empty()) Zarr.set_chunks(chunksFromBbox);
        try{
            Zarr.write_zarray();
        }
        catch(const std::string &e){
            if(e == "unsupportedCompressor"){
                mexErrMsgIdAndTxt("zarr:zarrayError","Compressor: \"%s\" is not currently supported\n",Zarr.get_cname().c_str());
            }
            else if(e.find("cannotOpenZarray") != std::string::npos){
                mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s for writing. Try checking permissions and path.\n",e.substr(e.find(':')+1).c_str());
            }
            else mexZarrError(e);
        }
    }
    else{
        Zarr.set_shape(endCoords);

        if(metadataExists){
            for(uint64_t d = 0; d < nDims; d++){
                if(endCoords[d]-startCoords[d] != iDims[d]) mexErrMsgIdAndTxt("zarr:inputError","Bounding box size does not match the size of the input data");
            }
        }
        else {
            try{
                Zarr.write_zarray();
            }
            catch(const std::string &e){
                if(e == "unsupportedCompressor"){
                    mexErrMsgIdAndTxt("zarr:zarrayError","Compressor: \"%s\" is not currently supported\n",Zarr.get_cname().c_str());
                }
                else if(e.find("cannotOpenZarray") != std::string::npos){
                    mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s for writing. Try checking permissions and path.\n",e.substr(e.find(':')+1).c_str());
                }
                else mexZarrError(e);
            }
        }

        try{
            Zarr = zarr(folderName);
        }
        catch(const std::string &e){
            if(e.find("metadataFileMissing") != std::string::npos){
                mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s for writing. Try checking permissions or the file path.\n",e.substr(e.find(':')+1).c_str());
            }
            else if(e == "metadataIncomplete"){
                mexErrMsgIdAndTxt("zarr:zarrayError","Metadata is incomplete. Check the .zarray or zarr.json file");
            }
            else mexZarrError(e);
        }
        
        // (the dtype option names the type of a new array: an existing one keeps its own)
        if(cropExisting && !dtypeOption.empty()){
            const DtypeInfo want = dtypeInfo(dtypeOption), have = dtypeInfo(Zarr.get_dtype());
            if(want.kind != have.kind || want.bytes != have.bytes || want.bfloat16 != have.bfloat16){
                mexErrMsgIdAndTxt("zarr:dataTypeError","%s has data type \"%s\", not \"%s\" (dtype is the data type of a new array)",
                                  folderName.c_str(),Zarr.get_dtype().c_str(),dtypeOption.c_str());
            }
        }
    }

    // The data in the array's data type: converted as MATLAB converts between classes
    // (integers rounded and saturated) when the types differ, and complex data with
    // its real and imaginary parts interleaved, as Zarr stores them
    {
        const DtypeInfo out = dtypeInfo(Zarr.get_dtype()), in = dtypeInfo(dataDtype);
        if(in.kind == 'c' && out.kind != 'c'){
            mexErrMsgIdAndTxt("zarr:dataTypeError","Complex data can only be written into complex64 or complex128 arrays, not \"%s\"",
                              Zarr.get_dtype().c_str());
        }
        bool convert = out.kind != in.kind || out.bytes != in.bytes || out.bfloat16 != in.bfloat16;
#if !(MX_HAS_INTERLEAVED_COMPLEX)
        convert = convert || in.kind == 'c';
#endif
        if(convert && out.kind){
            uint64_t size = 1;
            for(uint64_t d = 0; d < nDims; d++) size *= iDims[d];
            zarrC = mexConvertToDtype(prhs[1], Zarr.get_dtype(), size);
            if(!zarrC) mexErrMsgIdAndTxt("zarr:dataTypeError","Cannot convert the input data to the data type \"%s\"",Zarr.get_dtype().c_str());
        }
    }

    if(!crop){
        startCoords.assign(nDims,0);
        endCoords.assign(nDims,0);
        for(uint64_t d = 0; d < nDims; d++) endCoords[d] = Zarr.get_shape(d);
    }
    for(uint64_t d = 0; d < nDims; d++){
        if(endCoords[d] > Zarr.get_shape(d)) mexErrMsgIdAndTxt("zarr:inputError","Upper bound is invalid");
    }
    std::vector<uint64_t> writeShape(nDims);
    for(uint64_t d = 0; d < nDims; d++) writeShape[d] = endCoords[d]-startCoords[d];

    Zarr.set_chunkInfo(startCoords, endCoords);
    bool err = 0;
    // The write machinery is element-width based, so one generic path covers
    // every supported dtype
    const std::string &dtypeOut = Zarr.get_dtype();
    const uint64_t outBytes = Zarr.dtypeBytes();
    if(!dtypeInfo(dtypeOut).kind){
        free(zarrC);
        mexErrMsgIdAndTxt("zarr:dataTypeError","Data type \"%s\" is not supported",dtypeOut.c_str());
    }
    void* zarrArr = zarrC ? zarrC : mxGetData(prhs[1]);
    err = parallelWriteZarr(Zarr, zarrArr, startCoords, endCoords, writeShape, outBytes*8, useUuid, crop, sparse);

    // zarrC is either a copy for data conversion or NULL
    free(zarrC);

    if(err) mexErrMsgIdAndTxt("zarr:writeError",Zarr.get_errString().c_str());
    
}
