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
    if(nrhs == 3){
        if(!mxIsClass(prhs[2], "string")){
            if(!mxIsChar(prhs[2])) mexErrMsgIdAndTxt("tiff:inputError","The third argument must be a string");
            compression = mxArrayToString(prhs[0]);
        }
        else{
            mxArray* mString[1];
            mxArray* mCharA[1];

            // Convert string to char array
            mString[0] = mxDuplicateArray(prhs[2]);
            mexCallMATLAB(1, mCharA, 1, mString, "char");
            compression = mxArrayToString(mCharA[0]);
        }
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
    if(nDims < 2 || nDims > 3) mexErrMsgIdAndTxt("tiff:inputError","Data must be 2D or 3D");
    uint64_t* dims = (uint64_t*) mxGetDimensions(prhs[1]);


    uint64_t x = dims[1], y = dims[0], z = dims[2], bits = 0, startSlice = 0;

    // For 2D images MATLAB passes in the 3rd dim as 0 so we set it to 1;
    if(!z){
        z = 1;
    }

    mxClassID mDType = mxGetClassID(prhs[1]);
    if(mDType == mxUINT8_CLASS){
        bits = 8;
    }
    else if(mDType == mxUINT16_CLASS){
        bits = 16;
    }
    else if(mDType == mxSINGLE_CLASS){
        bits = 32;
    }
    else if(mDType == mxDOUBLE_CLASS){
        bits = 64;
    }
    else{
        mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
    }

    uint64_t stripSize = 512;
    void* data = (void*)mxGetPr(prhs[1]);
    uint8_t err = writeTiffParallelWrapper(x,y,z,fileName,data,bits,startSlice,stripSize,mode,true,compression);

    if(err) mexErrMsgIdAndTxt("tiff:tiffError","An Error occured within the write function");
}
