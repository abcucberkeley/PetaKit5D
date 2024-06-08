#include <stdint.h>
#include <string.h>
#include <omp.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' max_pooling_3d_mex.cpp
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' max_pooling_3d_mex.c

template <typename T>
void max_pooling_3d_mex(const T* const &orig, T* const &out, const uint64_t &poolSizeX, const uint64_t &poolSizeY, const uint64_t &poolSizeZ, const uint64_t &origShapeX, const uint64_t &origShapeY, const uint64_t &origShapeZ, const uint64_t &shapeX, const uint64_t &shapeY, const uint64_t &shapeZ){
    const uint64_t origShapeXY = origShapeX * origShapeY;
    const uint64_t outShapeXY = shapeX * shapeY;
    const uint16_t nthread = omp_get_max_threads();
    
    if (nthread > 1){
        #pragma omp parallel for collapse(3)
        for (uint64_t z = 0; z < shapeZ; ++z) {
            for (uint64_t y = 0; y < shapeY; ++y) {
                for (uint64_t x = 0; x < shapeX; ++x) {
                    const uint64_t zp1 = (z + 1) * poolSizeZ;
                    const uint64_t zp1min = zp1 < origShapeZ ? zp1 : origShapeZ;

                    const uint64_t yp1 = (y + 1) * poolSizeY;
                    const uint64_t yp1min = yp1 < origShapeY ? yp1 : origShapeY;
    
                    const uint64_t xp1 = (x + 1) * poolSizeX;
                    const uint64_t xp1min = xp1 < origShapeX ? xp1 : origShapeX;
    
                    T max_value = 0;
                    // disable simd which is very slow for Windows
                    // #pragma omp parallel for simd reduction(max:max_value)
                    for (uint64_t zi = z * poolSizeZ; zi < zp1min; ++zi){
                        for (uint64_t yi = y * poolSizeY; yi < yp1min; yi++) {
                            for (uint64_t xi = x * poolSizeX; xi < xp1min; xi++) {
                                const T* p = orig + (zi * origShapeXY + yi * origShapeX + xi);
                                if (*p > max_value) max_value = *p;
                            }
                        }
                    }
    
                    T* p = out + (z * outShapeXY + y * shapeX + x);
                    *p = max_value;
                }
            }
        }
    }
    else{
        for (uint64_t z = 0; z < shapeZ; ++z) {
            for (uint64_t y = 0; y < shapeY; ++y) {
                for (uint64_t x = 0; x < shapeX; ++x) {
                    const uint64_t zp1 = (z + 1) * poolSizeZ;
                    const uint64_t zp1min = zp1 < origShapeZ ? zp1 : origShapeZ;
    
                    const uint64_t yp1 = (y + 1) * poolSizeY;
                    const uint64_t yp1min = yp1 < origShapeY ? yp1 : origShapeY;
    
                    const uint64_t xp1 = (x + 1) * poolSizeX;
                    const uint64_t xp1min = xp1 < origShapeX ? xp1 : origShapeX;
    
                    T max_value = 0;
                    for (uint64_t zi = z * poolSizeZ; zi < zp1min; ++zi){
                        for (uint64_t yi = y * poolSizeY; yi < yp1min; yi++) {
                            for (uint64_t xi = x * poolSizeX; xi < xp1min; xi++) {
                                const T* p = orig + (zi * origShapeXY + yi * origShapeX + xi);
                                if (*p > max_value) max_value = *p;
                            }
                        }
                    }
    
                    T* p = out + (z * outShapeXY + y * shapeX + x);
                    *p = max_value;
                }
            }
        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /*  check for proper number of arguments */
    if(nrhs!=2) 
        mexErrMsgTxt("Two input required.");
    if(nlhs!=1) 
        mexErrMsgTxt("One output required.");

    if(mxGetN(prhs[1]) != 3) mexErrMsgIdAndTxt("pooling:inputError","Input range for pooling size is not 3");

    /*  get the dimensions of the matrix input x */
    const uint64_t* dims = (uint64_t*)mxGetDimensions(prhs[0]);
    const uint64_t numDim = (uint64_t)mxGetNumberOfDimensions(prhs[0]);
    const uint64_t origShapeX = dims[0];
    const uint64_t origShapeY = (numDim <= 1) ? 1 : dims[1];
    const uint64_t origShapeZ = (numDim <= 2) ? 1 : dims[2];
    
    const uint64_t poolSizeX = (uint64_t)* mxGetPr(prhs[1]);
    const uint64_t poolSizeY = (uint64_t)* (mxGetPr(prhs[1])+1);
    const uint64_t poolSizeZ = (uint64_t)* (mxGetPr(prhs[1])+2);
    if(poolSizeX < 1 || poolSizeY < 1 || poolSizeZ < 1) mexErrMsgIdAndTxt("pooling:inputError","Pooling size must be at least 1");

    const uint64_t shapeX = origShapeX / poolSizeX + (origShapeX % poolSizeX != 0);
    const uint64_t shapeY = origShapeY / poolSizeY + (origShapeY % poolSizeY != 0);
    const uint64_t shapeZ = origShapeZ / poolSizeZ + (origShapeZ % poolSizeZ != 0);
    // printf("%llu %llu %llu %llu %llu %llu %llu %llu %llu\n", poolSizeX, poolSizeY, poolSizeZ, origShapeX, origShapeY, origShapeZ, shapeX, shapeY, shapeZ);
    
    uint64_t dim[3];
    dim[0] = shapeX;
    dim[1] = shapeY;
    dim[2] = shapeZ;

    mxClassID mDType = mxGetClassID(prhs[0]);
    if(mDType == mxLOGICAL_CLASS){
        plhs[0] = mxCreateLogicalArray(3, (mwSize*)dim);
        bool* out = (bool*)mxGetPr(plhs[0]);
        const bool* orig = (bool*)mxGetPr(prhs[0]);
        max_pooling_3d_mex<bool>(orig, out, poolSizeX, poolSizeY, poolSizeZ, origShapeX, origShapeY, origShapeZ, shapeX, shapeY, shapeZ);
    }
    else if(mDType == mxUINT8_CLASS){
        plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxUINT8_CLASS, mxREAL);
        uint8_t* out = (uint8_t*)mxGetPr(plhs[0]);
        const uint8_t* orig = (uint8_t*)mxGetPr(prhs[0]);
        max_pooling_3d_mex<uint8_t>(orig, out, poolSizeX, poolSizeY, poolSizeZ, origShapeX, origShapeY, origShapeZ, shapeX, shapeY, shapeZ);
    }
    else if(mDType == mxUINT16_CLASS){
        plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxUINT16_CLASS, mxREAL);
        uint16_t* out = (uint16_t*)mxGetPr(plhs[0]);
        const uint16_t* orig = (uint16_t*)mxGetPr(prhs[0]);
        max_pooling_3d_mex<uint16_t>(orig, out, poolSizeX, poolSizeY, poolSizeZ, origShapeX, origShapeY, origShapeZ, shapeX, shapeY, shapeZ);
    }
    else if(mDType == mxSINGLE_CLASS){
        plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxSINGLE_CLASS, mxREAL);
        float* out = (float*)mxGetPr(plhs[0]);
        const float* orig = (float*)mxGetPr(prhs[0]);
        max_pooling_3d_mex<float>(orig, out, poolSizeX, poolSizeY, poolSizeZ, origShapeX, origShapeY, origShapeZ, shapeX, shapeY, shapeZ);
    }
    else if(mDType == mxDOUBLE_CLASS){
        plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxDOUBLE_CLASS, mxREAL);
        double* out = (double*)mxGetPr(plhs[0]);
        const double* orig = (double*)mxGetPr(prhs[0]);
        max_pooling_3d_mex<double>(orig, out, poolSizeX, poolSizeY, poolSizeZ, origShapeX, origShapeY, origShapeZ, shapeX, shapeY, shapeZ);
    }
    else{
        mexErrMsgIdAndTxt("mat:dataTypeError","Data type not suppported");
    }
}

