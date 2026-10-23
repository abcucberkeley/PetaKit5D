#include <cstring>
#include "mex.h"
#include "../src/helperfunctions.h"
#include "../src/parallelreadzarr.h"
#include "../src/zarr.h"
#include "zarrmexhelpers.h"

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
    if(!nrhs) mexErrMsgIdAndTxt("zarr:inputError","This functions requires at least 1 argument");
    if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("zarr:inputError","The first argument must be a string");  

    std::vector<uint64_t> startCoords;
    std::vector<uint64_t> endCoords;
    std::vector<uint64_t> bboxVals;
    bool bbox = false;
    bool useCtx = true;
    bool sparse = false;
    int threads = 0;
    std::string folderName(mxArrayToString(prhs[0]));

    for(int i = 1; i < nrhs; i+=2){
        if(i+1 == nrhs) mexErrMsgIdAndTxt("zarr:inputError","Mismatched argument pair for input number %d\n",i+1);
        if(!mxIsChar(prhs[i])) mexErrMsgIdAndTxt("zarr:inputError","The argument in input location %d is not a string\n",i+1);
        std::string currInput = mxArrayToString(prhs[i]);

        if(currInput == "bbox"){
            // Skip bbox if it is empty
            if(!mxGetN(prhs[i+1])) continue;
            // Interpreted once the array's number of dimensions is known
            bbox = true;
            bboxVals = mexVector(prhs[i+1]);
        
        }
        else if(currInput == "sparse"){
            sparse = (bool)*((mxGetPr(prhs[i+1])));
        }
        else if(currInput == "threads"){
            threads = mexThreads(prhs[i+1]);
        }
        else{
            mexErrMsgIdAndTxt("zarr:inputError","The argument \"%s\" does not match the name of any supported input name.\n \
            Currently Supported Names: bbox, sparse, threads\n",currInput.c_str());
        }
    }
    const MexCallThreads callThreadsSet(threads);
    
    // Handle the tilde character in filenames on Linux/Mac
    #ifndef _WIN32
    folderName = expandTilde(folderName.c_str());
    #endif
    zarr Zarr;
    try{
        Zarr = zarr(folderName);
    }
    catch(const std::string &e){
        if(e.find("metadataFileMissing") != std::string::npos){
            mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s for reading. Try checking permissions or the file path.\n",e.substr(e.find(':')+1).c_str());
        }
        else if(e == "metadataIncomplete"){
            mexErrMsgIdAndTxt("zarr:zarrayError","Metadata is incomplete. Check the .zarray or zarr.json file");
        }
        else mexZarrError(e);
    }

    const uint64_t nDims = Zarr.get_ndims();
    if(bbox){
        std::vector<uint64_t> shape(nDims);
        for(uint64_t d = 0; d < nDims; d++) shape[d] = Zarr.get_shape(d);
        if(!mexParseRegion(bboxVals, nDims, startCoords, endCoords, shape)){
            mexErrMsgIdAndTxt("zarr:inputError","bbox must have %d values ([starts ends]) for this %dD array",(int)(2*nDims),(int)nDims);
        }
    }
    else{
        startCoords.assign(nDims,0);
        endCoords.assign(nDims,0);
        for(uint64_t d = 0; d < nDims; d++) endCoords[d] = Zarr.get_shape(d);
    }
    for(uint64_t d = 0; d < nDims; d++){
        if(endCoords[d] > Zarr.get_shape(d)) mexErrMsgIdAndTxt("zarr:inputError","Upper bound is invalid");
    }
    std::vector<uint64_t> readShape(nDims);
    for(uint64_t d = 0; d < nDims; d++) readShape[d] = endCoords[d]-startCoords[d];
    // MATLAB arrays have at least 2 dimensions
    std::vector<mwSize> dim(nDims < 2 ? 2 : nDims, 1);
    for(uint64_t d = 0; d < nDims; d++) dim[d] = readShape[d];
    // TESTING
    /*
    if(Zarr.get_order() == "C"){
        //Zarr.set_order("F");
        dim[0] = readShape[2];
        dim[2] = readShape[0];
        
    }
    */
    
    
    
    
    Zarr.set_chunkInfo(startCoords, endCoords);

    bool err = 0;
    uint64_t readSize = 1;
    for(uint64_t d = 0; d < nDims; d++) readSize *= readShape[d];

    // Map the zarr dtype to the matching MATLAB class. The read machinery is
    // element-width based, so one generic path covers every supported dtype. MATLAB
    // has no 16-bit float class: float16 and bfloat16 arrays read as single. Complex
    // arrays read as complex single and double.
    const std::string dtype = Zarr.get_dtype();
    const DtypeInfo t = dtypeInfo(dtype);
    mxClassID mxClass = mxUNKNOWN_CLASS;
    if(t.kind == 'u'){
        if(t.bytes == 1) mxClass = mxUINT8_CLASS;
        else if(t.bytes == 2) mxClass = mxUINT16_CLASS;
        else if(t.bytes == 4) mxClass = mxUINT32_CLASS;
        else if(t.bytes == 8) mxClass = mxUINT64_CLASS;
    }
    else if(t.kind == 'i'){
        if(t.bytes == 1) mxClass = mxINT8_CLASS;
        else if(t.bytes == 2) mxClass = mxINT16_CLASS;
        else if(t.bytes == 4) mxClass = mxINT32_CLASS;
        else if(t.bytes == 8) mxClass = mxINT64_CLASS;
    }
    else if(t.kind == 'f') mxClass = t.bytes == 8 ? mxDOUBLE_CLASS : mxSINGLE_CLASS;
    else if(t.kind == 'c') mxClass = t.bytes == 16 ? mxDOUBLE_CLASS : mxSINGLE_CLASS;
    else if(t.kind == 'b') mxClass = mxLOGICAL_CLASS;
    if(mxClass == mxUNKNOWN_CLASS) mexErrMsgIdAndTxt("zarr:dataTypeError","Data type \"%s\" is not supported",dtype.c_str());

    // The output starts as the fill value, which is what missing chunks read as
    const uint64_t bytes = Zarr.dtypeBytes();
    uint8_t fillElem[16];
    fillValueElement(Zarr.get_fill_value(), dtype, fillElem);
    bool fill = false;
    for(uint64_t b = 0; b < bytes; b++) fill = fill || fillElem[b];
    const bool cplx = t.kind == 'c';
    // 16-bit floats are read into a buffer and widened; so are complex numbers when
    // MATLAB keeps their real and imaginary parts apart
    bool viaBuffer = t.kind == 'f' && t.bytes == 2;
#if !(MX_HAS_INTERLEAVED_COMPLEX)
    viaBuffer = viaBuffer || cplx;
#endif
    if(mxClass == mxLOGICAL_CLASS) plhs[0] = mxCreateLogicalArray(dim.size(),dim.data());
    else if(fill || viaBuffer) plhs[0] = mxCreateUninitNumericArray(dim.size(),dim.data(),mxClass, cplx ? mxCOMPLEX : mxREAL);
    else plhs[0] = mxCreateNumericArray(dim.size(),dim.data(),mxClass, cplx ? mxCOMPLEX : mxREAL);
    void* zarrArr = viaBuffer ? malloc(readSize*bytes ? readSize*bytes : 1) : mxGetData(plhs[0]);
    if(!zarrArr) mexErrMsgIdAndTxt("zarr:readError","Out of memory");
    if(fill) fillElements(bytes, zarrArr, readSize, fillElem);
    else if(viaBuffer) memset(zarrArr, 0, readSize*bytes);
    err = parallelReadZarr(Zarr, zarrArr,startCoords,endCoords,readShape,bytes*8,useCtx,sparse);
    if(err){
        if(viaBuffer) free(zarrArr);
        mexErrMsgIdAndTxt("zarr:readError",Zarr.get_errString().c_str());
    }
    if(viaBuffer){
        const int64_t n = (int64_t)readSize;
        if(t.kind == 'f'){
            const uint16_t* src = (const uint16_t*)zarrArr;
            float* out = (float*)mxGetData(plhs[0]);
            const bool bf = t.bfloat16;
            #pragma omp parallel for
            for(int64_t i = 0; i < n; i++) out[i] = bf ? bfloat16ToFloat(src[i]) : halfToFloat(src[i]);
        }
#if !(MX_HAS_INTERLEAVED_COMPLEX)
        else if(t.bytes == 8){
            const float* src = (const float*)zarrArr;
            float* re = (float*)mxGetData(plhs[0]);
            float* im = (float*)mxGetImagData(plhs[0]);
            #pragma omp parallel for
            for(int64_t i = 0; i < n; i++){ re[i] = src[2*i]; im[i] = src[2*i+1]; }
        }
        else{
            const double* src = (const double*)zarrArr;
            double* re = (double*)mxGetData(plhs[0]);
            double* im = (double*)mxGetImagData(plhs[0]);
            #pragma omp parallel for
            for(int64_t i = 0; i < n; i++){ re[i] = src[2*i]; im[i] = src[2*i+1]; }
        }
#endif
        free(zarrArr);
    }
}
