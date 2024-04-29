#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_blending_3d_with_indexing_mex.cpp
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' feather_blending_3d_with_indexing_mex.c

template <typename T>
void feather_blending_3d_mex(const T* const &fmat, const float* const &dmat, T* &nmat, const uint64_t* const &bbox, const uint64_t* const &dims, const uint64_t* const &odims) {
    const uint64_t shapeX = dims[0];
    const uint64_t shapeY = dims[1];
    const uint64_t shapeZ = dims[2];
    const uint64_t shapeT = dims[3];
    const uint64_t shapeXY = shapeX * shapeY;
    const uint64_t shapeXYZ = shapeXY * shapeZ;
    const uint64_t outShapeX = odims[0];
    const uint64_t outShapeXY = outShapeX * odims[1];    
    const uint64_t bbox_offset = bbox[0] + bbox[1] * outShapeX + bbox[2] * outShapeXY;

    // get the number of thread
    const uint16_t nthread=omp_get_max_threads();

    // printf("%llu %llu %llu %llu %llu %llu %llu %llu %llu %llu\n", shapeX, shapeY, shapeZ, shapeT, shapeXY, shapeXYZ, outShapeX, outShapeXY, bbox_offset, nthread);    

    if (nthread > 1){
        #pragma omp parallel for collapse(3)
        for (uint64_t z = 0; z < shapeZ; z++) {
            for (uint64_t y = 0; y < shapeY; y++) {
                for (uint64_t x = 0; x < shapeX; x++) {
                    float sum_w = 0;
                    float nv = 0;
                    uint64_t ind_xyz = x + y * shapeX + z * shapeXY;
                    uint64_t out_ind_xyz = x + y * outShapeX + z * outShapeXY + bbox_offset;
                    
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
                    if (sum_w != 0) {
                        *(nmat + out_ind_xyz) = static_cast<T>(round(nv / sum_w));
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
                    uint64_t out_ind_xyz = x + y * outShapeX + z * outShapeXY + bbox_offset;
                    
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
                    if (sum_w != 0) {
                        *(nmat + out_ind_xyz) = static_cast<T>(round(nv / sum_w));
                    }
                }
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4) mexErrMsgIdAndTxt("feather_blending:inputError","Number of input arguments must be 4");
    if (nlhs != 0) mexErrMsgIdAndTxt("feather_blending:outputError","Number of output arguments must be 0");

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

    uint64_t* dimsC = (uint64_t*)mxGetDimensions(prhs[2]);
    uint64_t nodim = (uint64_t) mxGetNumberOfDimensions(prhs[2]);
    uint64_t odims[3] = {1, 1, 1};
    for(uint64_t t=0; t < nodim; t++){
        odims[t] = dimsC[t];
    }
    
    uint64_t bbox[6];
    uint64_t nbbox = (uint64_t) mxGetNumberOfElements(prhs[3]);
    if(nbbox != 6) mexErrMsgIdAndTxt("feather_blending:inputError","The size of bounding box must be 1x6.");
    for(uint64_t t=0; t<6; t++){
        if(t < 3) bbox[t] = (uint64_t)*(mxGetPr(prhs[3])+t)-1;
        else  bbox[t] = (uint64_t)*(mxGetPr(prhs[3])+t);
    }

    // check odims and bbox
    for(uint64_t t=0; t < 3; t++){
        if(odims[t] < bbox[t] || dims[t] != bbox[t+3]-bbox[t]) mexErrMsgIdAndTxt("feather_blending:inputError","The bounding box is out of bound.");
    }

    mxClassID mDType = mxGetClassID(prhs[0]);
    mxClassID mDType_d = mxGetClassID(prhs[1]);
    mxClassID mDType_o = mxGetClassID(prhs[2]);
    if(mDType_d != mxSINGLE_CLASS) mexErrMsgIdAndTxt("feather_blending:inputError","The data type of the distance map must be single!");    
    if(mDType != mDType_o) mexErrMsgIdAndTxt("feather_blending:inputError","The data type of the overlap region does not match that of the data!");
    if(mDType == mxUINT16_CLASS){
        uint16_t* fmat = (uint16_t*)mxGetPr(prhs[0]);
        float* dmat = (float*)mxGetPr(prhs[1]);
        uint16_t* nmat = (uint16_t*)mxGetPr(prhs[2]);
        
        feather_blending_3d_mex<uint16_t>(fmat, dmat, nmat, bbox, dims, odims);
    }    
    else if(mDType == mxSINGLE_CLASS){
        float* fmat = (float*)mxGetPr(prhs[0]);
        float* dmat = (float*)mxGetPr(prhs[1]);
        float* nmat = (float*)mxGetPr(prhs[2]);
        
        feather_blending_3d_mex<float>(fmat, dmat, nmat, bbox, dims, odims);
    }
    else{
        mexErrMsgIdAndTxt("feather_blending:dataTypeError","Data type not suppported");
    }
}
