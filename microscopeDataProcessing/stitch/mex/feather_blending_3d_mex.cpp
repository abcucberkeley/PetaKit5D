#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_blending_3d_mex.cpp
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' indexing4d_mex.c

template <typename T>
void feather_blending_3d_mex(const T* const &fmat, const float* const &dmat, float* &nmat, const uint64_t &shapeX, const uint64_t &shapeY, const uint64_t &shapeZ, const uint64_t &shapeT) {
    const uint64_t shapeXY = shapeX * shapeY;
    const uint64_t shapeXYZ = shapeX * shapeY * shapeZ;
    const int nthread = omp_get_num_threads();
    
    if (nthread > 1){
        #pragma omp parallel for collapse(3)
        for (uint64_t z = 0; z < shapeZ; z++) {
            for (uint64_t y = 0; y < shapeY; y++) {
                for (uint64_t x = 0; x < shapeX; x++) {
                    float sum_w = 0;
                    float nv = 0;
                    uint64_t ind_xyz = x + y * shapeX + z * shapeXY;
                    
                    // Calculate fmat and dmat pointers once
                    const T* f_ptr = fmat + ind_xyz;
                    const float* d_ptr = dmat + ind_xyz;
    
                    for (uint64_t t = 0; t < shapeT; t++) {
                        if (*f_ptr != 0) {
                            sum_w += *d_ptr;
                            nv += static_cast<float>(*f_ptr) * *d_ptr;
                        }
                        // Move to the next slice in fmat and dmat
                        f_ptr += shapeXYZ;
                        d_ptr += shapeXYZ;
                    }
                    // printf("%f %f\n", nv, sum_w);
                    if (sum_w != 0) {
                        *(nmat + ind_xyz) = nv / sum_w;
                    }
                }
            }
        }
    }
    else{
        for (uint64_t z = 0; z < shapeZ; z++) {
            for (uint64_t y = 0; y < shapeY; y++) {
                for (uint64_t x = 0; x < shapeX; x++) {
                    float sum_w = 0;
                    float nv = 0;
                    uint64_t ind_xyz = x + y * shapeX + z * shapeXY;
                    
                    // Calculate fmat and dmat pointers once
                    const T* f_ptr = fmat + ind_xyz;
                    const float* d_ptr = dmat + ind_xyz;
    
                    for (uint64_t t = 0; t < shapeT; t++) {
                        if (*f_ptr != 0) {
                            sum_w += *d_ptr;
                            nv += static_cast<float>(*f_ptr) * *d_ptr;
                        }
                        // Move to the next slice in fmat and dmat
                        f_ptr += shapeXYZ;
                        d_ptr += shapeXYZ;
                    }
                    // printf("%f %f\n", nv, sum_w);
                    if (sum_w != 0) {
                        *(nmat + ind_xyz) = nv / sum_w;
                    }
                }
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2) mexErrMsgIdAndTxt("feather_blending:inputError","Number of input arguments must be 2");
    if (nlhs != 1) mexErrMsgIdAndTxt("feather_blending:outputError","Number of output arguments must be 1");

    uint64_t* dimsA = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t* dimsB = (uint64_t*)mxGetDimensions(prhs[1]);
    uint64_t dims[4] = {1, 1, 1, 1};
    uint64_t ndim = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    for(uint64_t t=0; t < ndim; t++){
        if(dimsA[t] != dimsB[t]){
            mexErrMsgIdAndTxt("feather_blending:inputError","The size of distance map does not matach that of the data");
        }
        dims[t] = dimsA[t];
    }

    uint64_t shapeX = dims[0];
    uint64_t shapeY = dims[1];
    uint64_t shapeZ = dims[2];
    uint64_t shapeT = dims[3];

    uint64_t dim[3];
    dim[0] = shapeX;
    dim[1] = shapeY;
    dim[2] = shapeZ;

    mxClassID mDType = mxGetClassID(prhs[0]);
    mxClassID mDType_d = mxGetClassID(prhs[1]);
    // if(mDType != mDType_d) mexErrMsgIdAndTxt("feather_blending:inputError","The data type of the distance map does not match that of the data!");
    if(mDType == mxUINT16_CLASS){
        uint16_t* fmat = (uint16_t*)mxGetPr(prhs[0]);
        float* dmat = (float*)mxGetPr(prhs[1]);

        plhs[0] = mxCreateNumericArray(3,dim,mxSINGLE_CLASS, mxREAL);
        float* nmat = (float*)mxGetPr(plhs[0]);
        
        feather_blending_3d_mex<uint16_t>(fmat, dmat, nmat, shapeX, shapeY, shapeZ, shapeT);
    }    
    else if(mDType == mxSINGLE_CLASS){
        float* fmat = (float*)mxGetPr(prhs[0]);
        float* dmat = (float*)mxGetPr(prhs[1]);

        plhs[0] = mxCreateNumericArray(3,dim,mxSINGLE_CLASS, mxREAL);
        float* nmat = (float*)mxGetPr(plhs[0]);
        
        feather_blending_3d_mex<float>(fmat, dmat, nmat, shapeX, shapeY, shapeZ, shapeT);
    }
    else{
        mexErrMsgIdAndTxt("feather_blending:dataTypeError","Data type not suppported");
    }
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
