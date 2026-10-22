#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_distance_map_resize_3d_mex.c
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_distance_map_resize_3d_mex.c

// 04/24/2024: change the input distance map in the linear space

float fastPower(float base, int exponent) {
    // Only for integer exponent
    if (exponent < 0) {
        base = 1.0 / base;
        exponent = -exponent;
    }

    float result = 1.0;
    while (exponent > 0) {
        // Check the least significant bit of the exponent
        if (exponent & 1) {
            result *= base;
        }
        // Square the base and right-shift the exponent
        base *= base;
        exponent >>= 1;
    }

    return result;
}


float trilinearInterpolation(const float* dmat, const float x, const float y, const float z, const uint64_t shapeX, const uint64_t shapeY, const uint64_t shapeZ, const uint64_t shapeXY) {
    // const uint64_t shapeXY = shapeX * shapeY;

    uint64_t x0 = (uint64_t)x;
    uint64_t y0 = (uint64_t)y;
    uint64_t z0 = (uint64_t)z;

    uint64_t xi = x0 + 1 < shapeX ? 1 : 0;
    uint64_t yi = y0 + 1 < shapeY ? 1 : 0;
    uint64_t zi = z0 + 1 < shapeZ ? 1 : 0;

    float xd = x - x0;
    float yd = y - y0;
    float zd = z - z0;
    
    uint64_t xyz_ind = x0 + y0 * shapeX + z0 * shapeXY;
    uint64_t ys = yi * shapeX;
    uint64_t zs = zi * shapeXY;

    float c00 = *(dmat + xyz_ind) * (1 - xd) + *(dmat + xyz_ind + xi) * xd;
    float c10 = *(dmat + xyz_ind + ys) * (1 - xd) + *(dmat + xyz_ind + xi + ys) * xd;
    float c01 = *(dmat + xyz_ind + zs) * (1 - xd) + *(dmat + xyz_ind + xi + zs) * xd;
    float c11 = *(dmat + xyz_ind + ys + zs) * (1 - xd) + *(dmat + xyz_ind + xi + ys + zs) * xd;

    float c0 = c00 * (1 - yd) + c10 * yd;
    float c1 = c01 * (1 - yd) + c11 * yd;

    float c = c0 * (1 - zd) + c1 * zd;

    return c;
}


void feather_distance_map_resize_3d_mex(const float* dmat, float* rmat, const float d, const uint64_t shapeX, const uint64_t shapeY, const uint64_t shapeZ, const uint64_t rShapeX, const uint64_t rShapeY, const uint64_t rShapeZ) {
    const uint64_t shapeXY = shapeX * shapeY;
    const uint64_t rShapeXY = rShapeX * rShapeY;
    const int d_int = (int) d;
    const uint16_t nthread = omp_get_max_threads();
    
    const float xfactor = rShapeX > 1 ? ((float)shapeX - 1.0) / ((float)rShapeX - 1.0) : 1.0;
    const float yfactor = rShapeY > 1 ? ((float)shapeY - 1.0) / ((float)rShapeY - 1.0) : 1.0;
    const float zfactor = rShapeZ > 1 ? ((float)shapeZ - 1.0) / ((float)rShapeZ - 1.0) : 1.0;
    
    if (nthread > 1){
        #pragma omp parallel for collapse(3)
        for (uint64_t z = 0; z < rShapeZ; z++) {
            for (uint64_t y = 0; y < rShapeY; y++) {
                for (uint64_t x = 0; x < rShapeX; x++) {
                    float xr = (float)x * xfactor;
                    float yr = (float)y * yfactor;
                    float zr = (float)z * zfactor;
    
                    uint64_t ind_xyz = x + y * rShapeX + z * rShapeXY;
    
                    float dr = trilinearInterpolation(dmat, xr, yr, zr, shapeX, shapeY, shapeZ, shapeXY);
                    *(rmat + ind_xyz) = fastPower(dr, d_int);                
                }
            }
        }
    }
    else{
        for (uint64_t z = 0; z < rShapeZ; z++) {
            for (uint64_t y = 0; y < rShapeY; y++) {
                for (uint64_t x = 0; x < rShapeX; x++) {
                    float xr = (float)x * xfactor;
                    float yr = (float)y * yfactor;
                    float zr = (float)z * zfactor;
    
                    uint64_t ind_xyz = x + y * rShapeX + z * rShapeXY;
    
                    float dr = trilinearInterpolation(dmat, xr, yr, zr, shapeX, shapeY, shapeZ, shapeXY);
                    *(rmat + ind_xyz) = fastPower(dr, d_int);                
                }
            }
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
#ifdef __APPLE__
    // keep the mex file loaded until MATLAB exits: GCC's OpenMP runtime leaves
    // thread exit handlers in it, so unloading it crashes MATLAB on exit
    if (!mexIsLocked()) {
        mexLock();
    }
#endif
    if (nrhs != 3) mexErrMsgIdAndTxt("feather_blending:inputError","Number of input arguments must be 3");
    if (nlhs != 1) mexErrMsgIdAndTxt("feather_blending:outputError","Number of output arguments must be 1");

    uint64_t* dimsA = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t dims[3] = {1, 1, 1};
    uint64_t ndim = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    for(uint64_t t=0; t < ndim; t++){
        dims[t] = dimsA[t];
    }

    uint64_t shapeX = dims[0];
    uint64_t shapeY = dims[1];
    uint64_t shapeZ = dims[2];

    if(mxGetN(prhs[1]) != 3) mexErrMsgIdAndTxt("crop:inputError","Input range for bbox is not 3");

    uint64_t rShapeX = (uint64_t)*(mxGetPr(prhs[1]));
    uint64_t rShapeY = (uint64_t)*(mxGetPr(prhs[1]) + 1);
    uint64_t rShapeZ = (uint64_t)*(mxGetPr(prhs[1]) + 2);

    uint64_t dim[3];
    dim[0] = rShapeX;
    dim[1] = rShapeY;
    dim[2] = rShapeZ;

    if (shapeX > rShapeX || shapeY > rShapeY || shapeZ > rShapeZ) mexErrMsgIdAndTxt("distmap_resize:inputError","The resized dimensions must not be smaller than those of the input distance map!");

    float d = (float) mxGetScalar(prhs[2]);

    mxClassID mDType = mxGetClassID(prhs[0]);
    if(mDType == mxSINGLE_CLASS){
        float* dmat = (float*)mxGetPr(prhs[0]);

        plhs[0] = mxCreateNumericArray(3,(mwSize*)dim,mxSINGLE_CLASS, mxREAL);
        float* rmat = (float*)mxGetPr(plhs[0]);
        
        feather_distance_map_resize_3d_mex(dmat, rmat, d, shapeX, shapeY, shapeZ, rShapeX, rShapeY, rShapeZ);
    }
    else{
        mexErrMsgIdAndTxt("feather_blending:dataTypeError","Data type not suppported");
    }
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
