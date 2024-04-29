#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <omp.h>
#include "mex.h"

using namespace std;

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' any_4th_dim_mex.cpp
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' any_4th_dim_mex.c

// change to check if any of the tile contains a valid (nonzero) element in a given voxel in xyz

template <typename T>
void any_4th_dim_bak_mex(const T* const &fmat, const bool* const &amat, bool* &valid_flag, const uint64_t &shapeX, const uint64_t &shapeY, const uint64_t &shapeZ, const uint64_t &shapeT, const uint64_t &numA) {
    const uint64_t shapeXY = shapeX * shapeY;
    const uint64_t shapeXYZ = shapeX * shapeY * shapeZ;
    const int nthread = omp_get_max_threads();
    uint64_t counter = 0;        
    
    if (nthread > 1) {
        for (uint64_t t = 0; t < shapeT; t++){
            if (! *(amat + t)) continue;

            volatile bool flag = false;
            #pragma omp parallel for collapse(3) shared(flag)
            for (uint64_t z = 0; z < shapeZ; z++) {
                for (uint64_t y = 0; y < shapeY; y++) {
                    for (uint64_t x = 0; x < shapeX; x++) {
                        if(flag) continue;
    
                        if (*(fmat + x + y * shapeX + z * shapeXY + t * shapeXYZ) == 0) flag = true;
                    }
                }
            }
            *(valid_flag+counter) = !flag;
            counter += 1;
        }
    }
    else {
        for (uint64_t t = 0; t < shapeT; t++){
            if (! *(amat + t)) continue;

            bool flag = false;        
            for (uint64_t z = 0; z < shapeZ; z++) {
                for (uint64_t y = 0; y < shapeY; y++) {
                    for (uint64_t x = 0; x < shapeX; x++) {
                        if(flag) break;
    
                        if (*(fmat + x + y * shapeX + z * shapeXY + t * shapeXYZ) == 0) flag = true;
                    }
                }
            }
            *(valid_flag+counter) = !flag;
            counter += 1;
        }
    }
}


template <typename T>
void any_4th_dim_mex(const T* const &fmat, const bool* const &amat, bool* &valid_flag, const uint64_t &shapeX, const uint64_t &shapeY, const uint64_t &shapeZ, const uint64_t &shapeT, const uint64_t &numA) {
    const uint64_t shapeXY = shapeX * shapeY;
    const uint64_t shapeXYZ = shapeX * shapeY * shapeZ;
    const int nthread = omp_get_max_threads();
    
    if (nthread > 1) {
        volatile bool flag = false;
        #pragma omp parallel for collapse(3) shared(flag)
        for (uint64_t z = 0; z < shapeZ; z++) {
            for (uint64_t y = 0; y < shapeY; y++) {
                for (uint64_t x = 0; x < shapeX; x++) {
                    if(flag) continue;
                    
                    uint64_t counter = 0;
                    for (uint64_t t = 0; t < shapeT; t++){
                        if (amat[t] && ! *(fmat + x + y * shapeX + z * shapeXY + t * shapeXYZ)) counter ++;
                    }
                    if (counter == numA) flag = true;
                }
            }
        }
        *(valid_flag) = !flag;
    }
    else {
        bool flag = false;
        for (uint64_t z = 0; z < shapeZ; z++) {
            for (uint64_t y = 0; y < shapeY; y++) {
                for (uint64_t x = 0; x < shapeX; x++) {
                    if(flag) break;
                    
                    uint64_t counter = 0;

                    for (uint64_t t = 0; t < shapeT; t++){
                        if (amat[t] && ! *(fmat + x + y * shapeX + z * shapeXY + t * shapeXYZ)) counter ++;
                    }

                    if (counter == numA) flag = true;
                }
            }
        }
        *(valid_flag) = !flag;
    }
}


template <typename T>
void get_4th_dim_unvalid_bbox(const T* const &fmat, const bool* const &amat, uint64_t* &ubbox, const uint64_t &shapeX, const uint64_t &shapeY, const uint64_t &shapeZ, const uint64_t &shapeT, const uint64_t &numA) {
    const uint64_t shapeXY = shapeX * shapeY;
    const uint64_t shapeXYZ = shapeX * shapeY * shapeZ;
    ubbox[0] = shapeX;
    ubbox[1] = shapeY;
    ubbox[2] = shapeZ;
    
    const int nthread = omp_get_max_threads();
    
    if (nthread > 0) {
        #pragma omp parallel for collapse(3) shared(ubbox)
        for (uint64_t z = 0; z < shapeZ; z++) {
            for (uint64_t y = 0; y < shapeY; y++) {
                for (uint64_t x = 0; x < shapeX; x++) {
                    uint64_t counter = 0;
                    for (uint64_t t = 0; t < shapeT; t++){
                        if (amat[t] && ! *(fmat + x + y * shapeX + z * shapeXY + t * shapeXYZ)) counter ++;
                    }
                    if (counter == numA){
                        ubbox[0] = std::min(ubbox[0], x+1);
                        ubbox[1] = std::min(ubbox[1], y+1);
                        ubbox[2] = std::min(ubbox[2], z+1);
                        ubbox[3] = std::max(ubbox[3], x+1);
                        ubbox[4] = std::max(ubbox[4], y+1);
                        ubbox[5] = std::max(ubbox[5], z+1);
                    }
                }
            }
        }
    }
    else {
        for (uint64_t z = 0; z < shapeZ; z++) {
            for (uint64_t y = 0; y < shapeY; y++) {
                for (uint64_t x = 0; x < shapeX; x++) {
                    uint64_t counter = 0;
                    for (uint64_t t = 0; t < shapeT; t++){
                        if (amat[t] && ! *(fmat + x + y * shapeX + z * shapeXY + t * shapeXYZ)) counter ++;
                    }
                    if (counter == numA){
                        ubbox[0] = std::min(ubbox[0], x);
                        ubbox[1] = std::min(ubbox[1], y);
                        ubbox[2] = std::min(ubbox[2], z);
                        ubbox[3] = std::max(ubbox[3], x);
                        ubbox[4] = std::max(ubbox[4], y);
                        ubbox[5] = std::max(ubbox[5], z);
                    }
                }
            }
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2) mexErrMsgIdAndTxt("any_4th_dim_mex:inputError","Number of input arguments must be 2");
    if (nlhs != 2) mexErrMsgIdAndTxt("any_4th_dim_mex:outputError","Number of output arguments must be 2");

    uint64_t* dimsF = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t* dimsA = (uint64_t*)mxGetDimensions(prhs[1]);
    uint64_t dims[4] = {1, 1, 1, 1};
    uint64_t ndim = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    uint64_t shapeA = (uint64_t) mxGetNumberOfElements(prhs[1]);

    for(uint64_t t=0; t < ndim; t++){
        dims[t] = dimsF[t];
    }
    
    uint64_t shapeX = dims[0];
    uint64_t shapeY = dims[1];
    uint64_t shapeZ = dims[2];
    uint64_t shapeT = dims[3];
    
    if(shapeA != shapeT) {
        mexErrMsgIdAndTxt("any_4th_dim_mex:dataTypeError","The second input must be a logical array with same number element of the 4-th dimension of the first input");
    }

    mxClassID mAType = mxGetClassID(prhs[1]);
    if(mAType != mxLOGICAL_CLASS) {
        mexErrMsgIdAndTxt("any_4th_dim_mex:dataTypeError","The second input must be logical type");
    }
    
    uint64_t numA=0;
    bool* amat = (bool*) mxGetPr(prhs[1]);    
    for(uint64_t a=0; a < shapeA; a++){
        if(*(amat+a)){
            numA += 1;
        }
    }
    // printf("%d", numA);

    if(!numA) mexErrMsgIdAndTxt("any_4th_dim_mex:dataRangeError","All elements in the second input are false, no indices to check!");

    mxClassID mDType = mxGetClassID(prhs[0]);
    if(mDType == mxUINT16_CLASS){
        uint16_t* fmat = (uint16_t*)mxGetPr(prhs[0]);
        plhs[0] = mxCreateLogicalMatrix(1, 1);
        bool* valid_flag = (bool*)mxGetPr(plhs[0]);
        uint64_t bSize[2] = {1, 6};
        plhs[1] = mxCreateNumericArray(2, (mwSize*) bSize, mxUINT64_CLASS, mxREAL);
        uint64_t* ubbox = (uint64_t*)mxGetPr(plhs[1]);
        
        any_4th_dim_mex<uint16_t>(fmat, amat, valid_flag, shapeX, shapeY, shapeZ, shapeT, numA);
        if (! *valid_flag) get_4th_dim_unvalid_bbox<uint16_t>(fmat, amat, ubbox, shapeX, shapeY, shapeZ, shapeT, numA);
    }    
    else if(mDType == mxSINGLE_CLASS){
        float* fmat = (float*)mxGetPr(prhs[0]);
        plhs[0] = mxCreateLogicalMatrix(1, 1);
        bool* valid_flag = (bool*)mxGetPr(plhs[0]);
        uint64_t bSize[2] = {1, 6};        
        plhs[1] = mxCreateNumericArray(2, (mwSize*) bSize, mxUINT64_CLASS, mxREAL);
        uint64_t* ubbox = (uint64_t*)mxGetPr(plhs[1]);
        
        any_4th_dim_mex<float>(fmat, amat, valid_flag, shapeX, shapeY, shapeZ, shapeT, numA);
        if (! *valid_flag) get_4th_dim_unvalid_bbox<float>(fmat, amat, ubbox, shapeX, shapeY, shapeZ, shapeT, numA);
        // printf("%llu %llu %llu %B %llu %llu %llu %llu %llu %llu\n", shapeX, shapeY, shapeZ, *valid_flag, ubbox[0], ubbox[1], ubbox[2], ubbox[3], ubbox[4], ubbox[5]);        
    }
    else{
        mexErrMsgIdAndTxt("any_4th_dim_mex:dataTypeError","Data type not suppported");
    }
}       

