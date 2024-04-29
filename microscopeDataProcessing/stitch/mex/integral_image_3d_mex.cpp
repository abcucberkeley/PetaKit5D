#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' integral_image_3d_mex.cpp
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' integral_image_3d_mex.c

template <typename T>
void integral_image_x(const T* const &fmat, float* &nmat, const uint64_t* const &dims, const uint64_t* const &dimsT, const uint64_t* const &dimsI) {
    const uint64_t shapeX = dimsI[0];
    const uint64_t shapeY = dimsI[1];
    const uint64_t shapeZ = dimsI[2];
    const uint64_t shapeXY = shapeX * shapeY;

    const uint64_t AShapeX = dims[0];
    const uint64_t AShapeY = dims[1];
    const uint64_t AShapeZ = dims[2];
    const uint64_t AShapeXY = AShapeX * AShapeY;

    const uint64_t shapeT = dimsT[0];

    const uint16_t nthread=omp_get_max_threads();
    
    if (nthread > 1){
        #pragma omp parallel for collapse(2)
        for (uint64_t z = 0; z < AShapeZ; z++) {
            for (uint64_t y = 0; y < AShapeY; y++) {
                const uint64_t ind_xyz = y * shapeX + z * shapeXY;
                const uint64_t ind_xyz_1 = y * AShapeX + z * AShapeXY;
                for (uint64_t x = 0; x < shapeX; x++) {
                    if (x == 0) *(nmat + ind_xyz) = (float)(*(fmat + ind_xyz_1));
                    else if (x < shapeT) *(nmat + x + ind_xyz) = *(nmat + x + ind_xyz - 1) + (float)(*(fmat + x + ind_xyz_1));
                    else if (x < AShapeX) *(nmat + x + ind_xyz) = *(nmat + x + ind_xyz - 1) + (float)(*(fmat + x + ind_xyz_1)) - (float)(*(fmat + x + ind_xyz_1 - shapeT));
                    else *(nmat + x + ind_xyz) = *(nmat + x + ind_xyz - 1) - (float)(*(fmat + x + ind_xyz_1 - shapeT));
                }
            }
        }
    } else {
        for (uint64_t z = 0; z < AShapeZ; z++) {
            for (uint64_t y = 0; y < AShapeY; y++) {
                const uint64_t ind_xyz = y * shapeX + z * shapeXY;
                const uint64_t ind_xyz_1 = y * AShapeX + z * AShapeXY;
                for (uint64_t x = 0; x < shapeX; x++) {
                    if (x == 0) *(nmat + ind_xyz) = (float)(*(fmat + ind_xyz_1));
                    else if (x < shapeT) *(nmat + x + ind_xyz) = *(nmat + x + ind_xyz - 1) + (float)(*(fmat + x + ind_xyz_1));
                    else if (x < AShapeX) *(nmat + x + ind_xyz) = *(nmat + x + ind_xyz - 1) + (float)(*(fmat + x + ind_xyz_1)) - (float)(*(fmat + x + ind_xyz_1 - shapeT));
                    else *(nmat + x + ind_xyz) = *(nmat + x + ind_xyz - 1) - (float)(*(fmat + x + ind_xyz_1 - shapeT));
                }
            }
        }
    }
}


void integral_image_y(float* &nmat, const uint64_t* const &dims, const uint64_t* const &dimsT, const uint64_t* const &dimsI) {
    const uint64_t shapeX = dimsI[0];
    const uint64_t shapeY = dimsI[1];
    const uint64_t shapeZ = dimsI[2];
    const uint64_t shapeXY = shapeX * shapeY;

    const uint64_t AShapeX = dims[0];
    const uint64_t AShapeY = dims[1];
    const uint64_t AShapeZ = dims[2];
    const uint64_t AShapeXY = AShapeX * AShapeY;

    const uint64_t shapeT = dimsT[1];

    // printf("%llu %llu %llu %llu %llu %llu %llu\n", shapeX, shapeY, shapeZ, AShapeX, AShapeY, AShapeZ, shapeT);

    const uint16_t nthread=omp_get_max_threads();
    
    if (nthread > 1){
        #pragma omp parallel for collapse(2)
        for (uint64_t z = 0; z < AShapeZ; z++) {
            for (uint64_t x = 0; x < shapeX; x++) {
                float* n_ptr = nmat + x + z * shapeXY;
                for (uint64_t y = shapeY - 1; y > 0; y--) {
                    float tmp_post;
                    float tmp = *(n_ptr + y * shapeX);
                    // printf("%llu %llu %llu %llu\n", y, y+1, y-shapeT, y-shapeT);
                    // printf("%f %f\n", *(n_ptr + y * shapeX), tmp);
    
                    if (y == shapeY - 1) *(n_ptr + y * shapeX) = *(n_ptr + (shapeY - shapeT) * shapeX);
                    else if (y >= AShapeY) *(n_ptr + y * shapeX) = *(n_ptr + (y+1) * shapeX) + *(n_ptr + (y - shapeT + 1) * shapeX);
                    else if (y >= shapeT-1) *(n_ptr + y * shapeX)  = *(n_ptr + (y+1) * shapeX) - tmp_post + *(n_ptr + (y - shapeT + 1) * shapeX);                        
                    else *(n_ptr + y * shapeX)  = *(n_ptr + (y+1) * shapeX) - tmp_post;
                    
                    tmp_post = tmp;
    
                    // printf("%f %f\n", *(n_ptr + y * shapeX), tmp);
                }
            }
        }
    } else {
        for (uint64_t z = 0; z < AShapeZ; z++) {
            for (uint64_t x = 0; x < shapeX; x++) {
                float* n_ptr = nmat + x + z * shapeXY;
                for (uint64_t y = shapeY - 1; y > 0; y--) {
                    float tmp_post;
                    float tmp = *(n_ptr + y * shapeX);
                    // printf("%llu %llu %llu %llu\n", y, y+1, y-shapeT, y-shapeT);
                    // printf("%f %f\n", *(n_ptr + y * shapeX), tmp);
    
                    if (y == shapeY - 1) *(n_ptr + y * shapeX) = *(n_ptr + (shapeY - shapeT) * shapeX);
                    else if (y >= AShapeY) *(n_ptr + y * shapeX) = *(n_ptr + (y+1) * shapeX) + *(n_ptr + (y - shapeT + 1) * shapeX);
                    else if (y >= shapeT-1) *(n_ptr + y * shapeX)  = *(n_ptr + (y+1) * shapeX) - tmp_post + *(n_ptr + (y - shapeT + 1) * shapeX);                        
                    else *(n_ptr + y * shapeX)  = *(n_ptr + (y+1) * shapeX) - tmp_post;
                    
                    tmp_post = tmp;
    
                    // printf("%f %f\n", *(n_ptr + y * shapeX), tmp);
                }
            }
        }
    }
}


void integral_image_z(float* &nmat, const uint64_t* const &dims, const uint64_t* const &dimsT, const uint64_t* const &dimsI) {
    const uint64_t shapeX = dimsI[0];
    const uint64_t shapeY = dimsI[1];
    const uint64_t shapeZ = dimsI[2];
    const uint64_t shapeXY = shapeX * shapeY;

    const uint64_t AShapeX = dims[0];
    const uint64_t AShapeY = dims[1];
    const uint64_t AShapeZ = dims[2];
    const uint64_t AShapeXY = AShapeX * AShapeY;

    const uint64_t shapeT = dimsT[2];

    // printf("%llu %llu %llu %llu %llu %llu %llu\n", shapeX, shapeY, shapeZ, AShapeX, AShapeY, AShapeZ, shapeT);

    const uint16_t nthread=omp_get_max_threads();
    
    if (nthread > 1){
        #pragma omp parallel for collapse(2)
        for (uint64_t y = 0; y < shapeY; y++) {
            for (uint64_t x = 0; x < shapeX; x++) {
                float* n_ptr = nmat + x + y * shapeX;
                for (uint64_t z = shapeZ-1; z > 0; z--) {
                    float tmp_post;
                    float tmp = *(n_ptr + z * shapeXY);
                    
                    if (z == shapeZ - 1) *(n_ptr + z * shapeXY) = *(n_ptr + (shapeZ - shapeT) * shapeXY);
                    else if (z >= AShapeZ) *(n_ptr + z * shapeXY) = *(n_ptr + (z+1) * shapeXY) + *(n_ptr + (z - shapeT + 1) * shapeXY);
                    else if (z >= shapeT-1) *(n_ptr + z * shapeXY)  = *(n_ptr + (z+1) * shapeXY) - tmp_post + *(n_ptr + (z - shapeT + 1) * shapeXY);
                    else *(n_ptr + z * shapeXY)  = *(n_ptr + (z+1) * shapeXY) - tmp_post;
                    
                    tmp_post = tmp;
                }
            }
        }
    } else {
        for (uint64_t y = 0; y < shapeY; y++) {
            for (uint64_t x = 0; x < shapeX; x++) {
                float* n_ptr = nmat + x + y * shapeX;
                for (uint64_t z = shapeZ-1; z > 0; z--) {
                    float tmp_post;
                    float tmp = *(n_ptr + z * shapeXY);
                    
                    if (z == shapeZ - 1) *(n_ptr + z * shapeXY) = *(n_ptr + (shapeZ - shapeT) * shapeXY);
                    else if (z >= AShapeZ) *(n_ptr + z * shapeXY) = *(n_ptr + (z+1) * shapeXY) + *(n_ptr + (z - shapeT + 1) * shapeXY);
                    else if (z >= shapeT-1) *(n_ptr + z * shapeXY)  = *(n_ptr + (z+1) * shapeXY) - tmp_post + *(n_ptr + (z - shapeT + 1) * shapeXY);
                    else *(n_ptr + z * shapeXY)  = *(n_ptr + (z+1) * shapeXY) - tmp_post;
                    
                    tmp_post = tmp;
                }
            }
        }
    }
}


template <typename T>
void integral_image_3d_mex(const T* const &fmat, float* &nmat, const uint64_t* const &dims, const uint64_t* const &dimsT, const uint64_t* const &dimsI) {
    
    integral_image_x(fmat, nmat, dims, dimsT, dimsI);
    integral_image_y(nmat, dims, dimsT, dimsI);
    integral_image_z(nmat, dims, dimsT, dimsI);
    
}


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2) mexErrMsgIdAndTxt("integral image 3d:inputError", "Number of input arguments must be 2");
    if (nlhs != 1) mexErrMsgIdAndTxt("integral image 3d:outputError", "Number of output arguments must be 1");

    uint64_t* dimsA = (uint64_t*) mxGetDimensions(prhs[0]);
    uint64_t dims[3] = {1, 1, 1};
    uint64_t ndim = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    for(uint64_t t=0; t < ndim; t++){
        dims[t] = dimsA[t];
    }

    uint64_t nT = (uint64_t) mxGetNumberOfElements(prhs[1]);
    if(nT != 3) mexErrMsgIdAndTxt("integral image 3d:inputError", "The number of elements of size T must be 3");
    uint64_t dimsT[3];
    for (uint16_t t=0; t<3; t++) dimsT[t] = (uint64_t)*(mxGetPr(prhs[1]) + t);
    // printf("%d %d %d\n", dimsT[0], dimsT[1], dimsT[2]);
    for (uint16_t t=0; t<3; t++) if(!dimsT[t]) mexErrMsgIdAndTxt("integral image 3d:inputError", "The elements of size T must be greater than zeros");

    uint64_t dimsI[3];
    dimsI[0] = dims[0] + dimsT[0] - 1;
    dimsI[1] = dims[1] + dimsT[1] - 1;
    dimsI[2] = dims[2] + dimsT[2] - 1;

    mxClassID mDType = mxGetClassID(prhs[0]);
    if(mDType == mxUINT16_CLASS){
        uint16_t* fmat = (uint16_t*) mxGetPr(prhs[0]);

        plhs[0] = mxCreateNumericArray(3, (mwSize*)dimsI, mxSINGLE_CLASS, mxREAL);
        float* nmat = (float*)mxGetPr(plhs[0]);
        
        integral_image_3d_mex<uint16_t>(fmat, nmat, dims, dimsT, dimsI);
    }
    else if(mDType == mxSINGLE_CLASS){
        float* fmat = (float*)mxGetPr(prhs[0]);

        plhs[0] = mxCreateNumericArray(3,(mwSize*)dimsI, mxSINGLE_CLASS, mxREAL);
        float* nmat = (float*)mxGetPr(plhs[0]);

        integral_image_3d_mex<float>(fmat, nmat, dims, dimsT, dimsI);
    }
    else{
        mexErrMsgIdAndTxt("integral image 3d:dataTypeError", "Data type not suppported");
    }
}

