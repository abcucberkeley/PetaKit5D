#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <type_traits>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' volume_deskew_rotate_warp_mex.cpp
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' volume_deskew_rotate_warp_mex.cpp

// c++ version of simplified geometric transformation for deskew/rotation without resampling, and output to isotropic voxels.

// Function to interpolate the voxel value
float interpolate(const float* const &volume, const uint64_t &shapeX, const uint64_t &shapeY, const uint64_t &shapeZ, const uint64_t &shapeXY, const float &x, const float &y, const float &z) {
    
    int64_t x0 = std::floor(x);
    int64_t y0 = std::floor(y);
    int64_t z0 = std::floor(z);

    float yd = y - y0;
    float zd = z - z0;

    uint64_t ind = x0 + y0 * shapeX + z0 * shapeXY;
    float c000 = *(volume + ind);
    float c001 = *(volume + ind + shapeXY);
    float c010 = *(volume + ind + shapeX);
    float c011 = *(volume + ind + shapeX + shapeXY);

    float c0 = c000 * (1 - yd) + c010 * yd;
    float c1 = c001 * (1 - yd) + c011 * yd;

    return (c0 * (1 - zd) + c1 * zd);
}


uint16_t interpolate_16bit(const uint16_t* const &volume, const uint64_t &shapeX, const uint64_t &shapeY, const uint64_t &shapeZ, const uint64_t &shapeXY, const float &x, const float &y, const float &z) {
    
    int64_t x0 = std::floor(x);
    int64_t y0 = std::floor(y);
    int64_t z0 = std::floor(z);

    float yd = y - y0;
    float zd = z - z0;

    uint64_t ind = x0 + y0 * shapeX + z0 * shapeXY;
    uint16_t c000 = *(volume + ind);
    uint16_t c001 = *(volume + ind + shapeXY);
    uint16_t c010 = *(volume + ind + shapeX);
    uint16_t c011 = *(volume + ind + shapeX + shapeXY);

    float c0 = c000 * (1 - yd) + c010 * yd;
    float c1 = c001 * (1 - yd) + c011 * yd;

    return (uint16_t)(c0 * (1 - zd) + c1 * zd + 0.5);
}


// Perform the transformation on a volume
void transformVolume(const float* const &dmat, const float (&tmat)[4][4], float* &nmat, 
                     const uint64_t &shapeX, const uint64_t &shapeY, const uint64_t &shapeZ, 
                     const float &startX, const float &startY, const float &startZ, 
                     const uint64_t &outShapeX, const uint64_t &outShapeY, const uint64_t &outShapeZ){
    
    uint64_t shapeXY = shapeX * shapeY;
    uint64_t outShapeXY = outShapeX * outShapeY;

    #pragma omp parallel for collapse(2) 
    for (uint64_t z = 0; z < outShapeZ; ++z){
        for (uint64_t y = 0; y < outShapeY; ++y){
            float ty, tz;
            ty = tmat[1][2] * (z + startZ) + tmat[1][3];
	        if (ty < 0 || ty + 1 > shapeY){
		        continue;
	        }

            tz = tmat[2][1] * (y + startY) + tmat[2][2] * (z + startZ) + tmat[2][3];
            if (tz < 0 || tz + 1 > shapeZ){
                continue;
            }

            for (uint64_t x = 0; x < outShapeX; ++x){
                uint64_t ind = x + y * outShapeX + z * outShapeXY;
                *(nmat + ind) = interpolate(dmat, shapeX, shapeY, shapeZ, shapeXY, x+startX, ty, tz);
            }
        }
    }
}


void transformVolume_16bit(const uint16_t* const &dmat, const float (&tmat)[4][4], uint16_t* &nmat, 
                     const uint64_t &shapeX, const uint64_t &shapeY, const uint64_t &shapeZ, 
                     const float &startX, const float &startY, const float &startZ, 
                     const uint64_t &outShapeX, const uint64_t &outShapeY, const uint64_t &outShapeZ){
    
    uint64_t shapeXY = shapeX * shapeY;
    uint64_t outShapeXY = outShapeX * outShapeY;

    #pragma omp parallel for collapse(2) 
    for (uint64_t z = 0; z < outShapeZ; ++z){
        for (uint64_t y = 0; y < outShapeY; ++y){
            float ty, tz;
            ty = tmat[1][2] * (z + startZ) + tmat[1][3];
	        if (ty < 0 || ty + 1 > shapeY){
		        continue;
	        }
            
	        tz = tmat[2][1] * (y + startY) + tmat[2][2] * (z + startZ) + tmat[2][3];
            if (tz < 0 || tz + 1 > shapeZ){
                continue;
            }

            for (uint64_t x = 0; x < outShapeX; ++x){
                uint64_t ind = x + y * outShapeX + z * outShapeXY;
                *(nmat + ind) = interpolate_16bit(dmat, shapeX, shapeY, shapeZ, shapeXY, x+startX, ty, tz);
            }
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if (nrhs != 3) mexErrMsgIdAndTxt("volume_warp:inputError","Number of input arguments must be 3");
    if (nlhs != 1) mexErrMsgIdAndTxt("volume_warp:outputError","Number of output arguments must be 1");
    
    uint64_t* dimsA = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t dims[3] = {1, 1, 1};
    uint64_t ndim = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    for(uint64_t t=0; t < ndim; t++){
        dims[t] = dimsA[t];
    }
    
    uint64_t shapeX = dims[0];
    uint64_t shapeY = dims[1];
    uint64_t shapeZ = dims[2];
    
    if(mxGetNumberOfDimensions(prhs[1]) != 2 || mxGetM(prhs[1]) != 4 || mxGetN(prhs[1]) != 4) mexErrMsgIdAndTxt("volume_warp:inputError","Input transformation matrix must be 4x4!");
    float tmat[4][4];
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            tmat[i][j] = (float)*(mxGetPr(prhs[1])+i+j*4);
        }
    }

    if(mxGetN(prhs[2]) != 6) mexErrMsgIdAndTxt("volume_warp:inputError","Input range for bbox is not 6");
    
    float startX = (float)(*mxGetPr(prhs[2]))-1;
    float startY = (float)(*(mxGetPr(prhs[2]) + 1))-1;
    float startZ = (float)(*(mxGetPr(prhs[2]) + 2))-1;
    float endX = (float)(*(mxGetPr(prhs[2]) + 3));
    float endY = (float)(*(mxGetPr(prhs[2]) + 4));
    float endZ = (float)(*(mxGetPr(prhs[2]) + 5));

    uint64_t outShapeX = (uint64_t) (endX - startX);
    uint64_t outShapeY = (uint64_t) (endY - startY);
    uint64_t outShapeZ = (uint64_t) (endZ - startZ);

    uint64_t outDim[3];
    outDim[0] = outShapeX;
    outDim[1] = outShapeY;
    outDim[2] = outShapeZ;
    
    // printf("%llu %llu %llu %llu %llu %llu\n", shapeX, shapeY, shapeZ, outShapeX, outShapeY, outShapeZ);
    
    mxClassID mDType = mxGetClassID(prhs[0]);
    if(mDType == mxSINGLE_CLASS){
        float* dmat = (float*)mxGetPr(prhs[0]);

        plhs[0] = mxCreateNumericArray(3, (mwSize*)outDim, mxSINGLE_CLASS, mxREAL);
        float* nmat = (float*)mxGetPr(plhs[0]);
        bool isUint16Type = false;

        transformVolume(dmat, tmat, nmat, shapeX, shapeY, shapeZ, startX, startY, startZ, outShapeX, outShapeY, outShapeZ);
    }
    else if(mDType == mxUINT16_CLASS){
        uint16_t* dmat = (uint16_t*)mxGetPr(prhs[0]);

        plhs[0] = mxCreateNumericArray(3, (mwSize*)outDim, mxUINT16_CLASS, mxREAL);
        uint16_t* nmat = (uint16_t*)mxGetPr(plhs[0]);
        bool isUint16Type = true;

        transformVolume_16bit(dmat, tmat, nmat, shapeX, shapeY, shapeZ, startX, startY, startZ, outShapeX, outShapeY, outShapeZ);
    }
    else{
        mexErrMsgIdAndTxt("feather_blending:dataTypeError","Data type not suppported");
    }
}

