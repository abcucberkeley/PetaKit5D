#include <cstdint>
#include <cstring>
#include "mex.h"
#include "tiffio.h"
#include "../src/helperfunctions.h"
#include "../src/parallelreadtiff.h"


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if(nrhs < 1 || nrhs > 2) mexErrMsgIdAndTxt("tiff:inputError","This function takes one or two arguments only");
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

    uint8_t flipXY = 1;
    //uint8_t flipXY = 0;


    //if(nrhs > 2){
    //    flipXY = (uint8_t)*(mxGetPr(prhs[2]));
    //}


    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) mexErrMsgIdAndTxt("tiff:inputError","File \"%s\" cannot be opened",fileName);

    uint64_t x = 1,y = 1,z = 1,bits = 1, startSlice = 0;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &x);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &y);
    
    uint8_t imageJIm = 0;
    if(nrhs == 1){
        z = getImageSizeZ(fileName);
        if(isImageJIm(fileName)){
            imageJIm = 1;
            uint64_t tempZ = imageJImGetZ(fileName);
            if(tempZ) z = tempZ;
        }
    }
    else{
        if(mxGetN(prhs[1]) != 2){
            mexErrMsgIdAndTxt("tiff:inputError","Input range is not 2");
        }
        else{
            startSlice = (uint64_t)*(mxGetPr(prhs[1]))-1;
            z = (uint64_t)*((mxGetPr(prhs[1])+1))-startSlice;
            uint64_t maxSize = 0;
            if(isImageJIm(fileName)){
                imageJIm = 1;
                maxSize = imageJImGetZ(fileName);
            }
            else maxSize = getImageSizeZ(fileName);
            if (startSlice < 0 || startSlice+z > maxSize){
                mexErrMsgIdAndTxt("tiff:rangeOutOfBound","Range is out of bounds");
            }
        }
    }

    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits);
    uint64_t stripSize = 1;
    TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &stripSize);
    TIFFClose(tif);

    uint64_t dim[3];
    dim[0] = y;
    dim[1] = x;
    dim[2] = z;

    // Case for ImageJ
    uint8_t err = 0;
    if(imageJIm){
        if(bits == 8){
            plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxUINT8_CLASS, mxREAL);
            uint8_t* tiff = (uint8_t*)mxGetPr(plhs[0]);
            err = readTiffParallelImageJ(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 16){
            plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxUINT16_CLASS, mxREAL);
            uint16_t* tiff = (uint16_t*)mxGetPr(plhs[0]);
            err = readTiffParallelImageJ(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 32){
            plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxSINGLE_CLASS, mxREAL);
            float* tiff = (float*)mxGetPr(plhs[0]);
            err = readTiffParallelImageJ(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 64){
            plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxDOUBLE_CLASS, mxREAL);
            double* tiff = (double*)mxGetPr(plhs[0]);
            err = readTiffParallelImageJ(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else{
            mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
        }
    }
    // Case for 2D
    else if(z <= 1){
        if(bits == 8){
            plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxUINT8_CLASS, mxREAL);
            uint8_t* tiff = (uint8_t*)mxGetPr(plhs[0]);
            err = readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 16){
            plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxUINT16_CLASS, mxREAL);
            uint16_t* tiff = (uint16_t*)mxGetPr(plhs[0]);
            err = readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 32){
            plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxSINGLE_CLASS, mxREAL);
            float* tiff = (float*)mxGetPr(plhs[0]);
            err = readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 64){
            plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxDOUBLE_CLASS, mxREAL);
            double* tiff = (double*)mxGetPr(plhs[0]);
            err = readTiffParallel2D(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else{
            mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
        }
    }
    // Case for 3D
    else{
        if(bits == 8){
            plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxUINT8_CLASS, mxREAL);
            uint8_t* tiff = (uint8_t*)mxGetPr(plhs[0]);
            err = readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 16){
            plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxUINT16_CLASS, mxREAL);
            uint16_t* tiff = (uint16_t*)mxGetPr(plhs[0]);
            err = readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 32){
            plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxSINGLE_CLASS, mxREAL);
            float* tiff = (float*)mxGetPr(plhs[0]);
            err = readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else if(bits == 64){
            plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxDOUBLE_CLASS, mxREAL);
            double* tiff = (double*)mxGetPr(plhs[0]);
            err = readTiffParallel(x,y,z,fileName, (void*)tiff, bits, startSlice, stripSize, flipXY);
        }
        else{
            mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
        }
    }
    if(err) mexErrMsgIdAndTxt("tiff:tiffError","An Error occured within the read function");
}