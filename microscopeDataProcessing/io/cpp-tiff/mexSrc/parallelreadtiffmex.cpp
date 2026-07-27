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
    uint16_t spp = 1, planar = 1;
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG, &planar);
    TIFFClose(tif);

    uint64_t samples = spp;
    // ImageJ hyperstacks are grayscale (channels live in metadata); RGB ImageJ is out of scope.
    if(imageJIm) samples = 1;
    if(samples > 1 && planar == 2)
        mexErrMsgIdAndTxt("tiff:planarError","Planar (PLANARCONFIG=2) RGB TIFFs are not supported");

    // Channel-last MATLAB layout: 3D (y,x,z) grayscale, or 4D (y,x,c,z) for RGB/RGBA.
    uint64_t dim[4];
    dim[0] = y;
    dim[1] = x;
    mwSize ndim;
    if(samples > 1){ ndim = 4; dim[2] = samples; dim[3] = z; }
    else           { ndim = 3; dim[2] = z; }

    // Pick the MATLAB class from bit depth + sample format (2=signed int, 3=float, else unsigned).
    uint64_t sampleFormat = getSampleFormat(fileName);
    mxClassID classID = mxUINT8_CLASS;
    if(bits == 8)       classID = (sampleFormat == 2) ? mxINT8_CLASS  : mxUINT8_CLASS;
    else if(bits == 16) classID = (sampleFormat == 2) ? mxINT16_CLASS : mxUINT16_CLASS;
    else if(bits == 32) classID = (sampleFormat == 3) ? mxSINGLE_CLASS : (sampleFormat == 2) ? mxINT32_CLASS : mxUINT32_CLASS;
    else if(bits == 64) classID = (sampleFormat == 3) ? mxDOUBLE_CLASS : (sampleFormat == 2) ? mxINT64_CLASS : mxUINT64_CLASS;
    else                mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");

    plhs[0] = mxCreateNumericArray(ndim,(mwSize*)dim,classID, mxREAL);
    void* tiff = mxGetData(plhs[0]);

    // Bytes are moved purely by bit width, so the worker only needs `bits`.
    uint8_t err = 0;
    if(imageJIm)    err = readTiffParallelImageJ(x,y,z,fileName, tiff, bits, startSlice, stripSize, flipXY);
    else if(z <= 1) err = readTiffParallel2D(x,y,z,fileName, tiff, bits, startSlice, stripSize, flipXY);
    else            err = readTiffParallel(x,y,z,fileName, tiff, bits, startSlice, stripSize, flipXY);
    if(err) mexErrMsgIdAndTxt("tiff:tiffError","An Error occured within the read function");
}