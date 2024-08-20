#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <type_traits>
#include <immintrin.h>
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp -mavx2' LDFLAGS='$LDFLAGS -O3 -fopenmp -mavx2' volume_deskew_rotate_warp_mex.cpp
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' volume_deskew_rotate_warp_mex.c

// c++ version of simplified geometric transformation for deskew/rotation without resampling, and output to isotropic voxels.


#ifdef __APPLE__

void compute_weighted_sum_sequential(const float* a, const float* b, const float* c, const float* d, 
                          float* result, float wa, float wb, float wc, float wd, size_t len) {

    #pragma omp parallel for
    for (size_t i = 0; i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i];
        result[i] = res;
    }
}


void compute_weighted_sum_sequential_16bit(const uint16_t* a, const uint16_t* b, const uint16_t* c, const uint16_t* d,
                                   uint16_t* result, float wa, float wb, float wc, float wd, size_t len) {

    #pragma parallel for
    for (size_t i = 0; i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i] + 0.5f;
        result[i] = (uint16_t)res;
    }
}

#else

void compute_weighted_sum(const float* a, const float* b, const float* c, const float* d, 
                          float* result, float wa, float wb, float wc, float wd, size_t len) {

    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);
    const __m256 wc_vec = _mm256_set1_ps(wc);
    const __m256 wd_vec = _mm256_set1_ps(wd);

    #pragma omp parallel for
    for (size_t i = 0; i < len; i += 8) {
        // Load 8 floats from each vector
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);
        __m256 vc = _mm256_loadu_ps(&c[i]);
        __m256 vd = _mm256_loadu_ps(&d[i]);

        // Multiply each vector by its corresponding weight
        __m256 va_weighted = _mm256_mul_ps(va, wa_vec);
        __m256 vb_weighted = _mm256_mul_ps(vb, wb_vec);
        __m256 vc_weighted = _mm256_mul_ps(vc, wc_vec);
        __m256 vd_weighted = _mm256_mul_ps(vd, wd_vec);

        // Sum the weighted vectors
        __m256 sum = _mm256_add_ps(va_weighted, vb_weighted);
        sum = _mm256_add_ps(sum, vc_weighted);
        sum = _mm256_add_ps(sum, vd_weighted);

        // Store the result
        _mm256_storeu_ps(&result[i], sum);
    }

    // Handle any remaining elements (len not multiple of 8)
    for (size_t i = len - (len % 8); i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i];
        result[i] = res;
    }
}


void compute_weighted_sum_16bit(const uint16_t* a, const uint16_t* b, const uint16_t* c, const uint16_t* d,
                                   uint16_t* result, float wa, float wb, float wc, float wd, size_t len) {
    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);
    const __m256 wc_vec = _mm256_set1_ps(wc);
    const __m256 wd_vec = _mm256_set1_ps(wd);
    const __m256 half_vec = _mm256_set1_ps(0.5f);

    #pragma parallel for
    for (size_t i = 0; i + 16 <= len; i += 16) {
        // Load 16 uint16_t values (two 128-bit chunks)
        __m256i a_chunk = _mm256_loadu_si256((__m256i*)&a[i]);
        __m256i b_chunk = _mm256_loadu_si256((__m256i*)&b[i]);
        __m256i c_chunk = _mm256_loadu_si256((__m256i*)&c[i]);
        __m256i d_chunk = _mm256_loadu_si256((__m256i*)&d[i]);

        // Convert uint16_t to uint32_t and then to float
        __m256i a_low  = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(a_chunk));
        __m256i a_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(a_chunk, 1));
        __m256i b_low  = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(b_chunk));
        __m256i b_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(b_chunk, 1));
        __m256i c_low  = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(c_chunk));
        __m256i c_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(c_chunk, 1));
        __m256i d_low  = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(d_chunk));
        __m256i d_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(d_chunk, 1));

        __m256 a_float_low  = _mm256_cvtepi32_ps(a_low);
        __m256 a_float_high = _mm256_cvtepi32_ps(a_high);
        __m256 b_float_low  = _mm256_cvtepi32_ps(b_low);
        __m256 b_float_high = _mm256_cvtepi32_ps(b_high);
        __m256 c_float_low  = _mm256_cvtepi32_ps(c_low);
        __m256 c_float_high = _mm256_cvtepi32_ps(c_high);
        __m256 d_float_low  = _mm256_cvtepi32_ps(d_low);
        __m256 d_float_high = _mm256_cvtepi32_ps(d_high);

        // Perform weighted sum: wa * a + wb * b + wc * c + wd * d + 0.5
        __m256 result_low  = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(a_float_low, wa_vec), _mm256_mul_ps(b_float_low, wb_vec)),
                                           _mm256_add_ps(_mm256_mul_ps(c_float_low, wc_vec), _mm256_mul_ps(d_float_low, wd_vec)));
        result_low = _mm256_add_ps(result_low, half_vec);

        __m256 result_high = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(a_float_high, wa_vec), _mm256_mul_ps(b_float_high, wb_vec)),
                                           _mm256_add_ps(_mm256_mul_ps(c_float_high, wc_vec), _mm256_mul_ps(d_float_high, wd_vec)));
        result_high = _mm256_add_ps(result_high, half_vec);

        // Convert float back to uint32_t
        __m256i result_int_low  = _mm256_cvtps_epi32(result_low);
        __m256i result_int_high = _mm256_cvtps_epi32(result_high);

        // Pack the results back to uint16_t and store them
        __m256i result_uint16 = _mm256_packus_epi32(result_int_low, result_int_high);
        result_uint16 = _mm256_permute4x64_epi64(result_uint16, _MM_SHUFFLE(3, 1, 2, 0));

        _mm256_storeu_si256((__m256i*)&result[i], result_uint16);
    }

    // Handle remaining elements (if len is not a multiple of 16)
    for (size_t i = len - (len % 16); i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i] + 0.5f;
        result[i] = (uint16_t)res;
    }
}

#endif

// Perform the transformation on a volume
void transformVolume(const float* const &dmat, const float (&tmat)[4][4], float* &nmat, 
                     const uint64_t &shapeX, const uint64_t &shapeY, const uint64_t &shapeZ, 
                     const float &startX, const float &startY, const float &startZ, 
                     const uint64_t &outShapeX, const uint64_t &outShapeY, const uint64_t &outShapeZ){
    
    uint64_t shapeXY = shapeX * shapeY;
    uint64_t outShapeXY = outShapeX * outShapeY;

    int64_t x0 = floor(startX);

    #pragma omp parallel for
    for (uint64_t z = 0; z < outShapeZ; ++z){
        float ty;
        ty = tmat[1][2] * (z + startZ) + tmat[1][3];
        if (ty < 0 || ty + 1 > shapeY){
		    continue;
        }
        int64_t y0 = floor(ty);
        float wy = ty - y0;

        for (uint64_t y = 0; y < outShapeY; ++y){
            float tz;
	        tz = tmat[2][1] * (y + startY) + tmat[2][2] * (z + startZ) + tmat[2][3];
            if (tz < 0 || tz + 1 > shapeZ){
                continue;
            }

            int64_t z0 = floor(tz);
            float wz = tz - z0;

            float w00 = (1 - wy) * (1 - wz);
            float w01 = (1 - wy) * wz;
            float w10 = wy * (1 - wz);
            float w11 = wy * wz;

            uint64_t ind = x0 + y0 * shapeX + z0 * shapeXY;

            #ifndef __APPLE__
            // interpolation with simd vectorization
            compute_weighted_sum(dmat + ind, dmat + ind + shapeXY, dmat + ind + shapeX, dmat + ind + shapeX + shapeXY, 
                                nmat + y * outShapeX + z * outShapeXY, w00, w01, w10, w11, outShapeX);
            #else
            compute_weighted_sum_sequential(dmat + ind, dmat + ind + shapeXY, dmat + ind + shapeX, dmat + ind + shapeX + shapeXY, 
                                nmat + y * outShapeX + z * outShapeXY, w00, w01, w10, w11, outShapeX);
            #endif
        }
    }
}


void transformVolume_16bit(const uint16_t* const &dmat, const float (&tmat)[4][4], uint16_t* &nmat, 
                     const uint64_t &shapeX, const uint64_t &shapeY, const uint64_t &shapeZ, 
                     const float &startX, const float &startY, const float &startZ, 
                     const uint64_t &outShapeX, const uint64_t &outShapeY, const uint64_t &outShapeZ){

    uint64_t shapeXY = shapeX * shapeY;
    uint64_t outShapeXY = outShapeX * outShapeY;
    uint64_t x0 = floor(startX);

    #pragma omp parallel for
    for (uint64_t z = 0; z < outShapeZ; ++z){
        float ty;
        ty = tmat[1][2] * (z + startZ) + tmat[1][3];
        if (ty < 0 || ty + 1 > shapeY){
		    continue;
        }
        int64_t y0 = floor(ty);
        float wy = ty - y0;

        for (uint64_t y = 0; y < outShapeY; ++y){
            float tz;
	        tz = tmat[2][1] * (y + startY) + tmat[2][2] * (z + startZ) + tmat[2][3];
            if (tz < 0 || tz + 1 > shapeZ){
                continue;
            }
            int64_t z0 = floor(tz);
            float wz = tz - z0;

            float w00 = (1 - wy) * (1 - wz);
            float w01 = (1 - wy) * wz;
            float w10 = wy * (1 - wz);
            float w11 = wy * wz;

            uint64_t ind = x0 + y0 * shapeX + z0 * shapeXY;

            #ifndef __APPLE__
                // interpolation with simd vectorization
                compute_weighted_sum_16bit(dmat + ind, dmat + ind + shapeXY, dmat + ind + shapeX, dmat + ind + shapeX + shapeXY, 
                                    nmat + y * outShapeX + z * outShapeXY, w00, w01, w10, w11, outShapeX);
            #else
                compute_weighted_sum_sequential_16bit(dmat + ind, dmat + ind + shapeXY, dmat + ind + shapeX, dmat + ind + shapeX + shapeXY, 
                                    nmat + y * outShapeX + z * outShapeXY, w00, w01, w10, w11, outShapeX);
            #endif
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

    if (startX < 0 || startY < 0 || startZ < 0 || startX > endX || startY > endY || startZ > endZ)
        mexErrMsgIdAndTxt("volume_warp:inputError","Input bbox is invalid, the numbers must be nonnegative and the end ranges must be nonless than the start ranges.");

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
