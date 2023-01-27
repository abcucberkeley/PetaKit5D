#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing4d_mex.c

void indexing4d_mex(void* orig, void* region, uint64_t startX, uint64_t startY, uint64_t startZ, uint64_t startT, uint64_t endX, uint64_t endY, uint64_t endZ, uint64_t endT, uint64_t origShapeX, uint64_t origShapeY, uint64_t origShapeZ, uint64_t origShapeT, uint64_t shapeX, uint64_t shapeY, uint64_t shapeZ, uint64_t shapeT, uint64_t bits){
    uint64_t bytes = bits/8;

    #pragma omp parallel for collapse(3)
    for(uint64_t t = startT; t < endT; t++){
	    for(uint64_t z = startZ; z < endZ; z++){
		    for(uint64_t y = startY; y < endY; y++){
		        memcpy((uint8_t*)orig+((startX+y*origShapeX+z*origShapeX*origShapeY+t*origShapeX*origShapeY*origShapeZ)*bytes), (uint8_t*)region+(((y-startY)*shapeX+(z-startZ)*shapeX*shapeY+(t-startT)*shapeX*shapeY*shapeZ)*bytes), shapeX*bytes);
		    }
	    }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3) mexErrMsgIdAndTxt("indexing:inputError","Number of input arguments must be 3");
    if (nlhs != 0) mexErrMsgIdAndTxt("indexing:outputError","Number of output arguments must be 0");

    uint64_t startX = 0;
    uint64_t startY = 0;
    uint64_t startZ = 0;
    uint64_t startT = 0;
    uint64_t endX = 0;
    uint64_t endY = 0;
    uint64_t endZ = 0;
    uint64_t endT = 0;
    uint64_t shapeX = 0;
    uint64_t shapeY = 0;
    uint64_t shapeZ = 0;
    uint64_t shapeT = 0;
    uint64_t* dimsA = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t dims[4] = {1, 1, 1, 1};
    uint64_t ndim = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    for(uint64_t t=0; t < ndim; t++) dims[t] = dimsA[t];

    uint64_t origShapeX = dims[0];
    uint64_t origShapeY = dims[1];
    uint64_t origShapeZ = dims[2];
    uint64_t origShapeT = dims[3];

    if(mxGetN(prhs[1]) != 8) mexErrMsgIdAndTxt("indexing:inputError","Input range for bbox is not 8");
    startX = (uint64_t)*(mxGetPr(prhs[1]))-1;
    startY = (uint64_t)*((mxGetPr(prhs[1])+1))-1;
    startZ = (uint64_t)*((mxGetPr(prhs[1])+2))-1;
    startT = (uint64_t)*((mxGetPr(prhs[1])+3))-1;
    endX = (uint64_t)*((mxGetPr(prhs[1])+4));
    endY = (uint64_t)*((mxGetPr(prhs[1])+5));
    endZ = (uint64_t)*((mxGetPr(prhs[1])+6));
    endT = (uint64_t)*((mxGetPr(prhs[1])+7));
        
    uint64_t* rdimsT = (uint64_t*)mxGetDimensions(prhs[2]);
    uint64_t rdims[4] = {1, 1, 1, 1};
    uint64_t nrdim = (uint64_t) mxGetNumberOfDimensions(prhs[2]);
    for(uint64_t t=0; t < nrdim; t++) rdims[t] = rdimsT[t];

    if(startX+1 < 1 || startY+1 < 1 || startZ+1 < 1 || startT+1 < 1) mexErrMsgIdAndTxt("indexing:inputError","Lower bounds must be at least 1");

    if(endX > origShapeX || endY > origShapeY || endZ > origShapeZ || endT > origShapeT) mexErrMsgIdAndTxt("indexing:inputError","Upper bound is invalid");

    uint64_t dim[4];
    shapeX = endX-startX;
    shapeY = endY-startY;
    shapeZ = endZ-startZ;
    shapeT = endT-startT;
    dim[0] = shapeX;
    dim[1] = shapeY;
    dim[2] = shapeZ;
    dim[3] = shapeT;

    // printf("%llu %llu %llu %llu\n", rdims[0], rdims[1], rdims[2], rdims[3]);

    if(rdims[0] != dim[0] || rdims[1] != dim[1] || rdims[2] != dim[2] || rdims[3] != dim[3]) mexErrMsgIdAndTxt("indexing:inputError","Subregion size does not match the bounding box size");
    
    mxClassID mDType = mxGetClassID(prhs[0]);
    mxClassID mDType_region = mxGetClassID(prhs[2]);
    if(mDType != mDType_region) mexErrMsgIdAndTxt("indexing:inputError","The data type of the region does not match that of the data!");
    if(mDType == mxUINT8_CLASS){
        uint64_t bits = 8;
        uint8_t* orig = (uint8_t*)mxGetPr(prhs[0]);
        uint8_t* region = (uint8_t*)mxGetPr(prhs[2]);
        indexing4d_mex((void*)orig, (void*)region, startX,startY,startZ,startT,endX,endY,endZ,endT,origShapeX,origShapeY,origShapeZ,origShapeT,shapeX,shapeY,shapeZ,shapeT,bits);
    }
    else if(mDType == mxUINT16_CLASS){
        uint64_t bits = 16;
        uint16_t* orig = (uint16_t*)mxGetPr(prhs[0]);
        uint16_t* region = (uint16_t*)mxGetPr(prhs[2]);
        indexing4d_mex((void*)orig, (void*)region, startX,startY,startZ,startT,endX,endY,endZ,endT,origShapeX,origShapeY,origShapeZ,origShapeT,shapeX,shapeY,shapeZ,shapeT,bits);
    }
    else if(mDType == mxSINGLE_CLASS){
        uint64_t bits = 32;
        float* orig = (float*)mxGetPr(prhs[0]);
        float* region = (float*)mxGetPr(prhs[2]);
        indexing4d_mex((void*)orig, (void*)region, startX,startY,startZ,startT,endX,endY,endZ,endT,origShapeX,origShapeY,origShapeZ,origShapeT,shapeX,shapeY,shapeZ,shapeT,bits);
    }
    else if(mDType == mxDOUBLE_CLASS){
        uint64_t bits = 64;
        double* orig = (double*)mxGetPr(prhs[0]);
        double* region = (double*)mxGetPr(prhs[2]);
        indexing4d_mex((void*)orig, (void*)region, startX,startY,startZ,startT,endX,endY,endZ,endT,origShapeX,origShapeY,origShapeZ,origShapeT,shapeX,shapeY,shapeZ,shapeT,bits);
    }
    else{
        mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
    }
}
