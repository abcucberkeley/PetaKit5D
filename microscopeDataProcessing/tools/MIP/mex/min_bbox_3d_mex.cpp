#include <stdint.h>
#include <string.h>
#include <limits>
#include <omp.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' min_bbox_3d_mex.cpp
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' min_bbox_3d_mex.cpp

template <typename T>
void min_bbox_3d_mex(const T* const &orig, T* const &out, const uint64_t &startX, const uint64_t &startY, const uint64_t &startZ, const uint64_t &endX, const uint64_t &endY, const uint64_t &endZ, const uint64_t &origShapeX, const uint64_t &origShapeY, const uint64_t &origShapeZ){
    const uint64_t origShapeXY = origShapeX * origShapeY;
    const uint16_t nthread = omp_get_max_threads();
    
    T min_value = std::numeric_limits<T>::max();    
    if (nthread > 1){
        // #pragma omp parallel for collapse(3)
        for (uint64_t z = startZ; z < endZ; ++z) {
            for (uint64_t y = startY; y < endY; ++y) {
                for (uint64_t x = startX; x < endX; ++x) {
                    const T* p = orig + (z * origShapeXY + y * origShapeX + x);
                    if (*p < min_value) min_value = *p;
                }
            }
        }
    }
    else{
        for (uint64_t z = startZ; z < endZ; ++z) {
            for (uint64_t y = startY; y < endY; ++y) {
                for (uint64_t x = startX; x < endX; ++x) {
                    const T* p = orig + (z * origShapeXY + y * origShapeX + x);
                    if (*p < min_value) min_value = *p;
                }
            }
        }
    }
    *out = min_value;
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /*  check for proper number of arguments */
    if(nrhs!=2) 
        mexErrMsgTxt("Two input required.");
    if(nlhs!=1) 
        mexErrMsgTxt("One output required.");

    if(mxGetN(prhs[1]) != 6) mexErrMsgIdAndTxt("minbbox:inputError","The length of bbox is not 6");

    /*  get the dimensions of the matrix input x */
    const uint64_t* dims = (uint64_t*)mxGetDimensions(prhs[0]);
    const uint64_t numDim = (uint64_t)mxGetNumberOfDimensions(prhs[0]);
    const uint64_t origShapeX = dims[0];
    const uint64_t origShapeY = (numDim <= 1) ? 1 : dims[1];
    const uint64_t origShapeZ = (numDim <= 2) ? 1 : dims[2];
    
    const uint64_t startX = (uint64_t)* mxGetPr(prhs[1])-1;
    const uint64_t startY = (uint64_t)* (mxGetPr(prhs[1])+1)-1;
    const uint64_t startZ = (uint64_t)* (mxGetPr(prhs[1])+2)-1;
    const uint64_t endX = (uint64_t)* (mxGetPr(prhs[1])+3);
    const uint64_t endY = (uint64_t)* (mxGetPr(prhs[1])+4);
    const uint64_t endZ = (uint64_t)* (mxGetPr(prhs[1])+5);
    if(endX < startX || endY < startY || endZ < startZ || endX > origShapeX || endY > origShapeY || endZ > origShapeZ) mexErrMsgIdAndTxt("minbbox:inputError","Bbox is not valid");

    // printf("%llu %llu %llu %llu %llu %llu %llu %llu %llu\n", poolSizeX, poolSizeY, poolSizeZ, origShapeX, origShapeY, origShapeZ, shapeX, shapeY, shapeZ);
    
    uint64_t dim[1];
    dim[0] = 1;

    mxClassID mDType = mxGetClassID(prhs[0]);
    if(mDType == mxLOGICAL_CLASS){
        plhs[0] = mxCreateLogicalArray(1, (mwSize*)dim);
        bool* out = (bool*)mxGetPr(plhs[0]);
        const bool* orig = (bool*)mxGetPr(prhs[0]);
        min_bbox_3d_mex<bool>(orig, out, startX, startY, startZ, endX, endY, endZ, origShapeX, origShapeY, origShapeZ);
    }
    else if(mDType == mxUINT8_CLASS){
        plhs[0] = mxCreateNumericArray(1,(mwSize*)dim,mxUINT8_CLASS, mxREAL);
        uint8_t* out = (uint8_t*)mxGetPr(plhs[0]);
        const uint8_t* orig = (uint8_t*)mxGetPr(prhs[0]);
        min_bbox_3d_mex<uint8_t>(orig, out, startX, startY, startZ, endX, endY, endZ, origShapeX, origShapeY, origShapeZ);
    }
    else if(mDType == mxUINT16_CLASS){
        plhs[0] = mxCreateNumericArray(1,(mwSize*)dim,mxUINT16_CLASS, mxREAL);
        uint16_t* out = (uint16_t*)mxGetPr(plhs[0]);
        const uint16_t* orig = (uint16_t*)mxGetPr(prhs[0]);
        min_bbox_3d_mex<uint16_t>(orig, out, startX, startY, startZ, endX, endY, endZ, origShapeX, origShapeY, origShapeZ);
    }
    else if(mDType == mxSINGLE_CLASS){
        plhs[0] = mxCreateNumericArray(1,(mwSize*)dim,mxSINGLE_CLASS, mxREAL);
        float* out = (float*)mxGetPr(plhs[0]);
        const float* orig = (float*)mxGetPr(prhs[0]);
        min_bbox_3d_mex<float>(orig, out, startX, startY, startZ, endX, endY, endZ, origShapeX, origShapeY, origShapeZ);
    }
    else if(mDType == mxDOUBLE_CLASS){
        plhs[0] = mxCreateNumericArray(1,(mwSize*)dim,mxDOUBLE_CLASS, mxREAL);
        double* out = (double*)mxGetPr(plhs[0]);
        const double* orig = (double*)mxGetPr(prhs[0]);
        min_bbox_3d_mex<double>(orig, out, startX, startY, startZ, endX, endY, endZ, origShapeX, origShapeY, origShapeZ);
    }
    else{
        mexErrMsgIdAndTxt("mat:dataTypeError","Data type not suppported");
    }
}

