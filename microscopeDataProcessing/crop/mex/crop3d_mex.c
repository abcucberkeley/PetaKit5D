#include <stdint.h>
#include <string.h>
#include <omp.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' crop3d_mex.c
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' crop3d_mex.c

void crop3d_mex(const void* restrict orig, void* restrict crop, const uint64_t startX, const uint64_t startY, const uint64_t startZ, const uint64_t endX, const uint64_t endY, const uint64_t endZ, const uint64_t origShapeX, const uint64_t origShapeY, const uint64_t origShapeZ, const uint64_t shapeX, const uint64_t shapeY, const uint64_t shapeZ, const uint64_t bits){
    const uint64_t bytes = bits/8;
    const uint64_t cropShapeXY = shapeX*shapeY;
    const uint64_t origShapeXY = origShapeX*origShapeY;
    const int nthread = omp_get_num_threads();
    
    if (nthread > 1){
        #pragma omp parallel for collapse(2)
        for(uint64_t z = startZ; z < endZ; z++){
            for(uint64_t y = startY; y < endY; y++){
                memcpy((uint8_t*)crop+(((y-startY)*shapeX+(z-startZ)*cropShapeXY)*bytes),(uint8_t*)orig+((startX+(y*origShapeX)+(z*origShapeXY))*bytes),shapeX*bytes);
            }
        }
    }
    else{
        for(uint64_t z = startZ; z < endZ; z++){
            for(uint64_t y = startY; y < endY; y++){
                memcpy((uint8_t*)crop+(((y-startY)*shapeX+(z-startZ)*cropShapeXY)*bytes),(uint8_t*)orig+((startX+(y*origShapeX)+(z*origShapeXY))*bytes),shapeX*bytes);
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2) mexErrMsgIdAndTxt("crop:inputError","Number of input arguments must be 2");

    uint64_t startX = 0;
    uint64_t startY = 0;
    uint64_t startZ = 0;
    uint64_t endX = 0;
    uint64_t endY = 0;
    uint64_t endZ = 0;
    uint64_t shapeX = 0;
    uint64_t shapeY = 0;
    uint64_t shapeZ = 0;
    uint64_t* dims = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t numDim = (uint64_t*) mxGetNumberOfDimensions(prhs[0]);
    uint64_t origShapeX = dims[0];
    uint64_t origShapeY = (numDim <= 1) ? 1 : dims[1];    
    uint64_t origShapeZ = (numDim <= 2) ? 1 : dims[2];

    if(mxGetN(prhs[1]) != 6) mexErrMsgIdAndTxt("crop:inputError","Input range for bbox is not 6");
    startX = (uint64_t)*(mxGetPr(prhs[1]))-1;
    startY = (uint64_t)*((mxGetPr(prhs[1])+1))-1;
    startZ = (uint64_t)*((mxGetPr(prhs[1])+2))-1;
    endX = (uint64_t)*((mxGetPr(prhs[1])+3));
    endY = (uint64_t)*((mxGetPr(prhs[1])+4));
    endZ = (uint64_t)*((mxGetPr(prhs[1])+5));
        
    if(startX+1 < 1 || startY+1 < 1 || startZ+1 < 1) mexErrMsgIdAndTxt("crop:inputError","Lower bounds must be at least 1");

    if(endX > origShapeX || endY > origShapeY || endZ > origShapeZ) mexErrMsgIdAndTxt("crop:inputError","Upper bound is invalid");

    uint64_t dim[3];
    shapeX = endX-startX;
    shapeY = endY-startY;
    shapeZ = endZ-startZ;
    dim[0] = shapeX;
    dim[1] = shapeY;
    dim[2] = shapeZ;

    mxClassID mDType = mxGetClassID(prhs[0]);
    if(mDType == mxUINT8_CLASS){
        uint64_t bits = 8;
        plhs[0] = mxCreateNumericArray(3,dim,mxUINT8_CLASS, mxREAL);
        uint8_t* crop = (uint8_t*)mxGetPr(plhs[0]);
        uint8_t* orig = (uint8_t*)mxGetPr(prhs[0]);
        crop3d_mex((void*)orig,(void*)crop,startX,startY,startZ,endX,endY,endZ,origShapeX,origShapeY,origShapeZ,shapeX,shapeY,shapeZ,bits);
    }
    else if(mDType == mxUINT16_CLASS){
        uint64_t bits = 16;
        plhs[0] = mxCreateNumericArray(3,dim,mxUINT16_CLASS, mxREAL);
        uint16_t* crop = (uint16_t*)mxGetPr(plhs[0]);
        uint16_t* orig = (uint16_t*)mxGetPr(prhs[0]);
        uint64_t size = origShapeX*origShapeY*origShapeZ;
        crop3d_mex((void*)orig,(void*)crop,startX,startY,startZ,endX,endY,endZ,origShapeX,origShapeY,origShapeZ,shapeX,shapeY,shapeZ,bits);
    }
    else if(mDType == mxSINGLE_CLASS){
        uint64_t bits = 32;
        plhs[0] = mxCreateNumericArray(3,dim,mxSINGLE_CLASS, mxREAL);
        float* crop = (float*)mxGetPr(plhs[0]);
        float* orig = (float*)mxGetPr(prhs[0]);
        crop3d_mex((void*)orig,(void*)crop,startX,startY,startZ,endX,endY,endZ,origShapeX,origShapeY,origShapeZ,shapeX,shapeY,shapeZ,bits);
    }
    else if(mDType == mxDOUBLE_CLASS){
        uint64_t bits = 64;
        plhs[0] = mxCreateNumericArray(3,dim,mxDOUBLE_CLASS, mxREAL);
        double* crop = (double*)mxGetPr(plhs[0]);
        double* orig = (double*)mxGetPr(prhs[0]);
        crop3d_mex((void*)orig,(void*)crop,startX,startY,startZ,endX,endY,endZ,origShapeX,origShapeY,origShapeZ,shapeX,shapeY,shapeZ,bits);
    }
    else{
        mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
    }
}
