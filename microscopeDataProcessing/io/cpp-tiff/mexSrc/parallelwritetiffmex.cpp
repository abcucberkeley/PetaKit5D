#include <cstdint>
#include <cstring>
#include <string>
#include "mex.h"
#include "tiffio.h"
#include "../src/helperfunctions.h"
#include "../src/parallelwritetiff.h"


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if(nrhs < 2) mexErrMsgIdAndTxt("tiff:inputError","This function requires at least 2 arguments");

    // Check if the fileName is a char array or matlab style
    char* fileName = NULL;
    if(!mxIsClass(prhs[0], "string")){
        if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("tiff:inputError","The first argument must be a string");
        fileName = mxArrayToString(prhs[0]);
    }
    else{ 
        mxArray* mString[1];
        mxArray* mCharA[1];

        // Convert string to char array
        mString[0] = mxDuplicateArray(prhs[0]);
        mexCallMATLAB(1, mCharA, 1, mString, "char");
        fileName = mxArrayToString(mCharA[0]);
    }
    if(mxIsEmpty(prhs[1])) mexErrMsgIdAndTxt("tiff:inputError","All input data axes must be of at least size 1");

    std::string compression = "lzw";
    if(nrhs >= 4){
        if(!mxIsClass(prhs[3], "string")){
            if(!mxIsChar(prhs[3])) mexErrMsgIdAndTxt("tiff:inputError","The fourth argument must be a string");
            compression = mxArrayToString(prhs[3]);
        }
        else{
            mxArray* mString[1];
            mxArray* mCharA[1];

            // Convert string to char array
            mString[0] = mxDuplicateArray(prhs[3]);
            mexCallMATLAB(1, mCharA, 1, mString, "char");
            compression = mxArrayToString(mCharA[0]);
        }
    }

    // Optional fifth argument: tile size (0/omitted = striped, the default).
    uint64_t tileSize = 0;
    if(nrhs >= 5){
        if(!mxIsNumeric(prhs[4]) || mxIsEmpty(prhs[4]))
            mexErrMsgIdAndTxt("tiff:inputError","The fifth argument (tile size) must be a number");
        double t = mxGetScalar(prhs[4]);
        if(t < 0) mexErrMsgIdAndTxt("tiff:inputError","Tile size cannot be negative");
        tileSize = (uint64_t)t;
        if(tileSize && tileSize % 16)
            mexErrMsgIdAndTxt("tiff:inputError","Tile size must be a multiple of 16");
    }

    // Optional sixth argument: true writes 3D (y, x, c) data as one RGB (c = 3)
    // or RGBA (c = 4) image instead of c grayscale slices. 4D (y, x, c, z) data
    // is always color.
    bool rgb = false;
    if(nrhs >= 6){
        if(!(mxIsLogical(prhs[5]) || mxIsNumeric(prhs[5])) || mxGetNumberOfElements(prhs[5]) != 1)
            mexErrMsgIdAndTxt("tiff:inputError","The sixth argument (rgb) must be true or false");
        rgb = mxGetScalar(prhs[5]) != 0;
    }

    // Handle the tilde character in filenames on Linux/Mac
    #ifndef _WIN32
    if(strchr(fileName,'~')) fileName = expandTilde(fileName);
    #endif

    // Check if folder exists, if not then make it (recursive if needed)
    char* folderName = strdup(fileName);
    char *lastSlash = NULL;
    #ifdef _WIN32
    lastSlash = strrchr(folderName, '\\');
    #else
    lastSlash = strrchr(folderName, '/');
    #endif
    if(lastSlash){
        *lastSlash = '\0';
        FILE* f = fopen(folderName,"r");
        if(f){
            fclose(f);
        }
        else{
            mkdirRecursive(folderName);
        }
    }
    free(folderName);

    TIFFSetWarningHandler(DummyHandler);
    const char* mode;
    if(nrhs > 2){
         mode = mxArrayToString(prhs[2]);
    }
    else{
        mode = "w";
    }
    int nDims = (int) mxGetNumberOfDimensions(prhs[1]);
    if(nDims < 2 || nDims > 4) mexErrMsgIdAndTxt("tiff:inputError","Data must be 2D, 3D, or 4D (y, x, c, z) color");
    uint64_t* dims = (uint64_t*) mxGetDimensions(prhs[1]);


    uint64_t x = dims[1], y = dims[0], bits = 0, startSlice = 0;

    // mxGetDimensions returns exactly nDims entries, so reading dims[2] for a 2D
    // image is out of bounds -- MATLAB does not pad the array with a trailing 0.
    uint64_t z = (nDims > 2) ? dims[2] : 1;
    uint64_t samples = 1;
    if(nDims == 4 || rgb){
        if(nDims < 3) mexErrMsgIdAndTxt("tiff:inputError","RGB data must be (y, x, c) or (y, x, c, z)");
        samples = dims[2];
        z = (nDims > 3) ? dims[3] : 1;
        if(samples != 3 && samples != 4)
            mexErrMsgIdAndTxt("tiff:inputError","Color data must have 3 (RGB) or 4 (RGBA) samples in its third dimension");
    }
    if(!z){
        z = 1;
    }

    mxClassID mDType = mxGetClassID(prhs[1]);
    uint16_t sampleFormat = SAMPLEFORMAT_UINT;
    if(mDType == mxUINT8_CLASS){        bits = 8;  sampleFormat = SAMPLEFORMAT_UINT; }
    else if(mDType == mxINT8_CLASS){    bits = 8;  sampleFormat = SAMPLEFORMAT_INT; }
    else if(mDType == mxUINT16_CLASS){  bits = 16; sampleFormat = SAMPLEFORMAT_UINT; }
    else if(mDType == mxINT16_CLASS){   bits = 16; sampleFormat = SAMPLEFORMAT_INT; }
    else if(mDType == mxUINT32_CLASS){  bits = 32; sampleFormat = SAMPLEFORMAT_UINT; }
    else if(mDType == mxINT32_CLASS){   bits = 32; sampleFormat = SAMPLEFORMAT_INT; }
    else if(mDType == mxSINGLE_CLASS){  bits = 32; sampleFormat = SAMPLEFORMAT_IEEEFP; }
    else if(mDType == mxUINT64_CLASS){  bits = 64; sampleFormat = SAMPLEFORMAT_UINT; }
    else if(mDType == mxINT64_CLASS){   bits = 64; sampleFormat = SAMPLEFORMAT_INT; }
    else if(mDType == mxDOUBLE_CLASS){  bits = 64; sampleFormat = SAMPLEFORMAT_IEEEFP; }
    else{
        mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
    }

    // 0 lets the writer choose the rows per strip from the image size.
    uint64_t stripSize = 0;
    void* data = mxGetData(prhs[1]);
    uint8_t err = writeTiffParallelWrapper(x,y,z,fileName,data,bits,startSlice,stripSize,mode,true,compression,sampleFormat,tileSize,samples);

    if(err) mexErrMsgIdAndTxt("tiff:tiffError","An Error occured within the write function");
}
