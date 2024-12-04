#include <cstdint>
#include <cstring>
#include "mex.h"
#include "tiffio.h"
#include "../src/helperfunctions.h"


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if(nrhs != 1) mexErrMsgIdAndTxt("tiff:inputError","This function requires one argument only");
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

    // Handle the tilde character in filenames on Linux/Mac
    #ifndef _WIN32
    if(strchr(fileName,'~')) fileName = expandTilde(fileName);
    #endif
    
    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) mexErrMsgIdAndTxt("tiff:inputError","File \"%s\" cannot be opened",fileName);
    
    uint64_t x = 1,y = 1,z = 1;    
    if(nrhs == 1){
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &x);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &y);
        z = getImageSizeZ(fileName);
    }
    else{
        mexErrMsgIdAndTxt("tiff:inputError","Function only accepts one input argument");       
    }
  
    TIFFClose(tif);
	if(isImageJIm(fileName)){
		uint64_t tempZ = imageJImGetZ(fileName);
		if(tempZ) z = tempZ;
	}

    plhs[0] = mxCreateNumericMatrix(1,3,mxDOUBLE_CLASS, mxREAL);
    double* dims = (double*)mxGetPr(plhs[0]);
    dims[0] = y;
    dims[1] = x;
    dims[2] = z;
    
}
