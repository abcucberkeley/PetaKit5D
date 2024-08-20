#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <inttypes.h>
#include <string.h>
#include <immintrin.h> // For AVX/AVX2
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp -mavx2' LDFLAGS='$LDFLAGS -O3 -fopenmp -mavx2' skewed_space_interp_defined_stepsize_mex.c
// macOS
// mex -v CC="/usr/local/bin/gcc-12" CXX="/usr/local/bin/g++-12" COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_defined_stepsize_mex.c

// use simd for weighted sum calculation


#ifdef __APPLE__

void compute_weighted_sum_two_vectors_sequential(const float* a, const float* b, float* result, 
                                      float wa, float wb, size_t len) {
    #pragma omp parallel for
    for (size_t i = 0; i < len; ++i) {
        float res = wa * a[i] + wb * b[i];
        result[i] = res;
    }
}

void compute_weighted_sum_sequential(const float* a, const float* b, const float* c, const float* d, 
                          float* result, float wa, float wb, float wc, float wd, size_t len) {

    #pragma omp parallel for
    for (size_t i = 0; i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i];
        result[i] = res;
    }
}


void compute_weighted_sum_two_vectors_sequential_16bit(const uint16_t* a, const uint16_t* b, uint16_t* result, float wa, float wb, size_t len) {

    #pragma parallel for
    for (size_t i = 0; i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + 0.5f;
        result[i] = (uint16_t)res;
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

void compute_weighted_sum_two_vectors(const float* a, const float* b, float* result, 
                                      float wa, float wb, size_t len) {

    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);

    #pragma omp parallel for
    for (size_t i = 0; i < len; i += 8) {
        // Load 8 floats from each vector
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);

        // Multiply each vector by its corresponding weight
        __m256 va_weighted = _mm256_mul_ps(va, wa_vec);
        __m256 vb_weighted = _mm256_mul_ps(vb, wb_vec);

        // Sum the weighted vectors
        __m256 sum = _mm256_add_ps(va_weighted, vb_weighted);

        // Store the result
        _mm256_storeu_ps(&result[i], sum);
    }

    // Handle any remaining elements (len not multiple of 8)
    for (size_t i = len - (len % 8); i < len; ++i) {
        float res = wa * a[i] + wb * b[i];
        result[i] = res;
    }
}


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


void compute_weighted_sum_two_vectors_16bit(const uint16_t* a, const uint16_t* b, uint16_t* result, float wa, float wb, size_t len) {

    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);
    const __m256 half_vec = _mm256_set1_ps(0.5f);

    #pragma parallel for
    for (size_t i = 0; i + 16 <= len; i += 16) {
        // Load 16 uint16_t values (two 128-bit chunks)
        __m256i a_chunk = _mm256_loadu_si256((__m256i*)&a[i]);
        __m256i b_chunk = _mm256_loadu_si256((__m256i*)&b[i]);

        // Convert uint16_t to uint32_t and then to float
        __m256i a_low  = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(a_chunk));
        __m256i a_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(a_chunk, 1));
        __m256i b_low  = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(b_chunk));
        __m256i b_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(b_chunk, 1));

        __m256 a_float_low  = _mm256_cvtepi32_ps(a_low);
        __m256 a_float_high = _mm256_cvtepi32_ps(a_high);
        __m256 b_float_low  = _mm256_cvtepi32_ps(b_low);
        __m256 b_float_high = _mm256_cvtepi32_ps(b_high);

        // Perform weighted sum: wa * a + wb * b + wc * c + wd * d + 0.5
        __m256 result_low  = _mm256_add_ps(_mm256_mul_ps(a_float_low, wa_vec), _mm256_mul_ps(b_float_low, wb_vec));
        result_low = _mm256_add_ps(result_low, half_vec);

        __m256 result_high = _mm256_add_ps(_mm256_mul_ps(a_float_high, wa_vec), _mm256_mul_ps(b_float_high, wb_vec));
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
        float res = wa * a[i] + wb * b[i] + 0.5f;
        result[i] = (uint16_t)res;
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


void skewed_space_interp_defined_stepsize(const float* im, float* im_int, const float xstep, const float stepsize, const uint8_t Reverse,
                                          const uint64_t sz[3], const uint64_t outDim[3]){

    const uint64_t nint = outDim[2];
    const uint64_t nxy = sz[0] * sz[1];

    const uint64_t xa = (uint64_t)ceil(xstep);
    double* s_mat = (double*)malloc((nint)*sizeof(double));
    double* t_mat = (double*)malloc((nint)*sizeof(double));
    double* sw_mat = (double*)malloc((nint)*sizeof(double));
    double* tw_mat = (double*)malloc((nint)*sizeof(double));
    double* zw_mat = (double*)malloc((nint)*sizeof(double));

    // copy or interp z indices
    bool* z_copy_mat = (bool*)malloc(nint*sizeof(bool));
    uint64_t* z_ind_mat = (uint64_t*)malloc(nint*sizeof(uint64_t));
    uint64_t n_copy = 0;
    uint64_t n_interp = 0;

    for(uint64_t i = 0; i < nint; i++){
        double s;
        double t;
        double zout;
        double zw;
        // zout = ((double)(sz[2] - 1)) / ((double)(nint - 1)) * ((double)i);
        zout = stepsize * (double)i;
        zw = zout - floor(zout);

        if(Reverse){
            s = ceil(xstep) - xstep * zw;
            t = xstep * (1 - zw);
        }
        else{
            s = xstep * (1 - zw);
            t = ceil(xstep) - xstep * zw;
        }

        //distance to start
        double sw = s - floor(s);
        s = floor(s);
        double tw = t - floor(t);
        t = floor(t);

        if((Reverse && (sw < 1e-3 || (1.0 - sw) < 1e-3 || zw < 1e-3)) || (!Reverse && !tw)){
            z_copy_mat[i] = true;
        } else {
            z_copy_mat[i] = false;
        }

        s_mat[i] = s;
        t_mat[i] = t;
        sw_mat[i] = sw;
        tw_mat[i] = tw;
        zw_mat[i] = zw;

        z_ind_mat[i] = floor(zout);
        // printf("%llu %f %llu %f %f %f %f %f %llu\n", i, zout, z_ind_mat[i], s, t, sw, tw, zw, z_copy_mat[i]);
    }

    #pragma omp parallel for
    for(uint64_t z = 0; z < nint; z++){
        uint64_t zint = z;
        uint64_t zorg = z_ind_mat[zint];
        // printf("%llu %llu % llu\n", z, zint, zorg);

        bool z_copy = z_copy_mat[zint];
        if(z_copy){
            memcpy(im_int+zint*nxy, im+zorg*nxy, nxy*sizeof(float));
            continue;
        }

        uint64_t s = (uint64_t)s_mat[zint];
        uint64_t t = (uint64_t)t_mat[zint];
        float sw = sw_mat[zint];
        float tw = tw_mat[zint];
        float zw = zw_mat[zint];

        float c00_w = (1.0 - sw) * (1.0 - zw);
        float c01_w = sw * (1.0 - zw);
        float c10_w = (1.0 - tw) * zw;
        float c11_w = tw * zw;
        
        #ifndef __APPLE__
        if(Reverse){
            #pragma omp parallel for
            for(uint64_t y = 0; y < sz[1]; y++){
                if(y<xa-s){
                    compute_weighted_sum_two_vectors(im + (y+t) * sz[0] + (zorg+1) * nxy, im + (y+t+1) * sz[0] + (zorg+1) * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c10_w, c11_w, sz[0]);
                }else if(y < sz[1]-t-1){
                    compute_weighted_sum(im + (y+s-xa) * sz[0] + zorg * nxy, im + (y+s-xa+1) * sz[0] + zorg * nxy, 
                                        im + (y+t) * sz[0] + (zorg+1) * nxy, im + (y+t+1) * sz[0] + (zorg+1) * nxy, 
                                        im_int+y*sz[0]+zint*nxy, c00_w, c01_w, c10_w, c11_w, sz[0]);
                } else{
                    compute_weighted_sum_two_vectors(im + (y+s-xa) * sz[0] + zorg * nxy, im + (y+s-xa+1) * sz[0] + zorg * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c00_w, c01_w, sz[0]);
                }
            }
        } else {
            #pragma omp parallel for
            for(uint64_t y = 0; y < sz[1]; y++){
                if(y<xa-s){
                    compute_weighted_sum_two_vectors(im + (y+s) * sz[0] + zorg * nxy, im + (y+s+1) * sz[0] + zorg * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c00_w, c01_w, sz[0]);

                }else if(y < sz[1]-t-1){
                    compute_weighted_sum(im + (y+s) * sz[0] + zorg * nxy, im + (y+s+1) * sz[0] + zorg * nxy, 
                                        im + (y+t-xa) * sz[0] + (zorg+1) * nxy, im + (y+t-xa+1) * sz[0] + (zorg+1) * nxy, 
                                        im_int+y*sz[0]+zint*nxy, c00_w, c01_w, c10_w, c11_w, sz[0]);
                } else{
                    compute_weighted_sum_two_vectors(im + (y+t-xa) * sz[0] + (zorg+1) * nxy, im + (y+t-xa+1) * sz[0] + (zorg+1) * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c10_w, c11_w, sz[0]);
                }
            }
        }
        #else
        if(Reverse){
            #pragma omp parallel for
            for(uint64_t y = 0; y < sz[1]; y++){
                if(y<xa-s){
                    compute_weighted_sum_two_vectors_sequential(im + (y+t) * sz[0] + (zorg+1) * nxy, im + (y+t+1) * sz[0] + (zorg+1) * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c10_w, c11_w, sz[0]);
                }else if(y < sz[1]-t-1){
                    compute_weighted_sum_sequential(im + (y+s-xa) * sz[0] + zorg * nxy, im + (y+s-xa+1) * sz[0] + zorg * nxy, 
                                        im + (y+t) * sz[0] + (zorg+1) * nxy, im + (y+t+1) * sz[0] + (zorg+1) * nxy, 
                                        im_int+y*sz[0]+zint*nxy, c00_w, c01_w, c10_w, c11_w, sz[0]);
                } else{
                    compute_weighted_sum_two_vectors_sequential(im + (y+s-xa) * sz[0] + zorg * nxy, im + (y+s-xa+1) * sz[0] + zorg * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c00_w, c01_w, sz[0]);
                }
            }
        } else {
            #pragma omp parallel for
            for(uint64_t y = 0; y < sz[1]; y++){
                if(y<xa-s){
                    compute_weighted_sum_two_vectors_sequential(im + (y+s) * sz[0] + zorg * nxy, im + (y+s+1) * sz[0] + zorg * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c00_w, c01_w, sz[0]);

                }else if(y < sz[1]-t-1){
                    compute_weighted_sum_sequential(im + (y+s) * sz[0] + zorg * nxy, im + (y+s+1) * sz[0] + zorg * nxy, 
                                        im + (y+t-xa) * sz[0] + (zorg+1) * nxy, im + (y+t-xa+1) * sz[0] + (zorg+1) * nxy, 
                                        im_int+y*sz[0]+zint*nxy, c00_w, c01_w, c10_w, c11_w, sz[0]);
                } else{
                    compute_weighted_sum_two_vectors_sequential(im + (y+t-xa) * sz[0] + (zorg+1) * nxy, im + (y+t-xa+1) * sz[0] + (zorg+1) * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c10_w, c11_w, sz[0]);
                }
            }
        }
        #endif
    }

    free(s_mat);
    free(t_mat);
    free(sw_mat);
    free(tw_mat);
    free(zw_mat);
    free(z_copy_mat);
    free(z_ind_mat);
}


void skewed_space_interp_defined_stepsize_16bit(const uint16_t* im, uint16_t* im_int, const float xstep, const float stepsize, const uint8_t Reverse,
                                          const uint64_t sz[3], const uint64_t outDim[3]){

    const uint64_t nint = outDim[2];
    const uint64_t nxy = sz[0] * sz[1];

    const uint64_t xa = (uint64_t)ceil(xstep);
    double* s_mat = (double*)malloc((nint)*sizeof(double));
    double* t_mat = (double*)malloc((nint)*sizeof(double));
    double* sw_mat = (double*)malloc((nint)*sizeof(double));
    double* tw_mat = (double*)malloc((nint)*sizeof(double));
    double* zw_mat = (double*)malloc((nint)*sizeof(double));

    // copy or interp z indices
    bool* z_copy_mat = (bool*)malloc(nint*sizeof(bool));
    uint64_t* z_ind_mat = (uint64_t*)malloc(nint*sizeof(uint64_t));
    uint64_t n_copy = 0;
    uint64_t n_interp = 0;

    for(uint64_t i = 0; i < nint; i++){
        double s;
        double t;
        double zout;
        double zw;
        // zout = ((double)(sz[2] - 1)) / ((double)(nint - 1)) * ((double)i);
        zout = stepsize * (double)i;
        zw = zout - floor(zout);

        if(Reverse){
            s = ceil(xstep) - xstep * zw;
            t = xstep * (1 - zw);
        }
        else{
            s = xstep * (1 - zw);
            t = ceil(xstep) - xstep * zw;
        }

        //distance to start
        double sw = s - floor(s);
        s = floor(s);
        double tw = t - floor(t);
        t = floor(t);

        if((Reverse && (sw < 1e-3 || (1.0 - sw) < 1e-3 || zw < 1e-3)) || (!Reverse && !tw)){
            z_copy_mat[i] = true;
        } else {
            z_copy_mat[i] = false;
        }

        s_mat[i] = s;
        t_mat[i] = t;
        sw_mat[i] = sw;
        tw_mat[i] = tw;
        zw_mat[i] = zw;

        z_ind_mat[i] = floor(zout);
        // printf("%llu %f %llu %f %f %f %f %f %llu\n", i, zout, z_ind_mat[i], s, t, sw, tw, zw, z_copy_mat[i]);
    }

    #pragma omp parallel for
    for(uint64_t z = 0; z < nint; z++){
        uint64_t zint = z;
        uint64_t zorg = z_ind_mat[zint];
        // printf("%llu %llu % llu\n", z, zint, zorg);

        bool z_copy = z_copy_mat[zint];
        if(z_copy){
            memcpy(im_int+zint*nxy, im+zorg*nxy, nxy*sizeof(uint16_t));
            continue;
        }

        uint64_t s = (uint64_t)s_mat[zint];
        uint64_t t = (uint64_t)t_mat[zint];
        float sw = sw_mat[zint];
        float tw = tw_mat[zint];
        float zw = zw_mat[zint];

        float c00_w = (1.0 - sw) * (1.0 - zw);
        float c01_w = sw * (1.0 - zw);
        float c10_w = (1.0 - tw) * zw;
        float c11_w = tw * zw;

        #ifndef __APPLE__
        if(Reverse){
            #pragma omp parallel for
            for(uint64_t y = 0; y < sz[1]; y++){
                if(y<xa-s){
                    compute_weighted_sum_two_vectors_16bit(im + (y+t) * sz[0] + (zorg+1) * nxy, im + (y+t+1) * sz[0] + (zorg+1) * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c10_w, c11_w, sz[0]);
                }else if(y < sz[1]-t-1){
                    compute_weighted_sum_16bit(im + (y+s-xa) * sz[0] + zorg * nxy, im + (y+s-xa+1) * sz[0] + zorg * nxy, 
                                        im + (y+t) * sz[0] + (zorg+1) * nxy, im + (y+t+1) * sz[0] + (zorg+1) * nxy, 
                                        im_int+y*sz[0]+zint*nxy, c00_w, c01_w, c10_w, c11_w, sz[0]);
                } else{
                    compute_weighted_sum_two_vectors_16bit(im + (y+s-xa) * sz[0] + zorg * nxy, im + (y+s-xa+1) * sz[0] + zorg * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c00_w, c01_w, sz[0]);
                }
            }
        } else {
            #pragma omp parallel for
            for(uint64_t y = 0; y < sz[1]; y++){
                if(y<xa-s){
                    compute_weighted_sum_two_vectors_16bit(im + (y+s) * sz[0] + zorg * nxy, im + (y+s+1) * sz[0] + zorg * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c00_w, c01_w, sz[0]);

                }else if(y < sz[1]-t-1){
                    compute_weighted_sum_16bit(im + (y+s) * sz[0] + zorg * nxy, im + (y+s+1) * sz[0] + zorg * nxy, 
                                        im + (y+t-xa) * sz[0] + (zorg+1) * nxy, im + (y+t-xa+1) * sz[0] + (zorg+1) * nxy, 
                                        im_int+y*sz[0]+zint*nxy, c00_w, c01_w, c10_w, c11_w, sz[0]);
                } else{
                    compute_weighted_sum_two_vectors_16bit(im + (y+t-xa) * sz[0] + (zorg+1) * nxy, im + (y+t-xa+1) * sz[0] + (zorg+1) * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c10_w, c11_w, sz[0]);
                }
            }
        }
        #else
        if(Reverse){
            #pragma omp parallel for
            for(uint64_t y = 0; y < sz[1]; y++){
                if(y<xa-s){
                    compute_weighted_sum_two_vectors_sequential_16bit(im + (y+t) * sz[0] + (zorg+1) * nxy, im + (y+t+1) * sz[0] + (zorg+1) * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c10_w, c11_w, sz[0]);
                }else if(y < sz[1]-t-1){
                    compute_weighted_sum_sequential_16bit(im + (y+s-xa) * sz[0] + zorg * nxy, im + (y+s-xa+1) * sz[0] + zorg * nxy, 
                                        im + (y+t) * sz[0] + (zorg+1) * nxy, im + (y+t+1) * sz[0] + (zorg+1) * nxy, 
                                        im_int+y*sz[0]+zint*nxy, c00_w, c01_w, c10_w, c11_w, sz[0]);
                } else{
                    compute_weighted_sum_two_vectors_sequential_16bit(im + (y+s-xa) * sz[0] + zorg * nxy, im + (y+s-xa+1) * sz[0] + zorg * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c00_w, c01_w, sz[0]);
                }
            }
        } else {
            #pragma omp parallel for
            for(uint64_t y = 0; y < sz[1]; y++){
                if(y<xa-s){
                    compute_weighted_sum_two_vectors_sequential_16bit(im + (y+s) * sz[0] + zorg * nxy, im + (y+s+1) * sz[0] + zorg * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c00_w, c01_w, sz[0]);

                }else if(y < sz[1]-t-1){
                    compute_weighted_sum_sequential_16bit(im + (y+s) * sz[0] + zorg * nxy, im + (y+s+1) * sz[0] + zorg * nxy, 
                                        im + (y+t-xa) * sz[0] + (zorg+1) * nxy, im + (y+t-xa+1) * sz[0] + (zorg+1) * nxy, 
                                        im_int+y*sz[0]+zint*nxy, c00_w, c01_w, c10_w, c11_w, sz[0]);
                } else{
                    compute_weighted_sum_two_vectors_sequential_16bit(im + (y+t-xa) * sz[0] + (zorg+1) * nxy, im + (y+t-xa+1) * sz[0] + (zorg+1) * nxy, 
                                                           im_int+y*sz[0]+zint*nxy, c10_w, c11_w, sz[0]);
                }
            }
        }
        #endif
    }

    free(s_mat);
    free(t_mat);
    free(sw_mat);
    free(tw_mat);
    free(zw_mat);
    free(z_copy_mat);
    free(z_ind_mat);
}


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{

    if (nrhs != 4) mexErrMsgIdAndTxt("volume_warp:inputError","Number of input arguments must be 3");
    if (nlhs != 1) mexErrMsgIdAndTxt("volume_warp:outputError","Number of output arguments must be 1");

    const float xstep = (float)*mxGetPr(prhs[1]);
    const float stepsize = (float)*mxGetPr(prhs[2]);
    const uint8_t Reverse = (uint8_t)mxIsLogicalScalarTrue(prhs[3]);

    uint64_t* dimsA = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t sz[3] = {1, 1, 1};
    uint64_t ndim = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    for(uint64_t t=0; t < ndim; t++){
        sz[t] = dimsA[t];
    }

    uint64_t nint = (uint64_t) (floor(round((float)(sz[2] - 1) / stepsize * 100000) / 100000) + 1);
    uint64_t outDim[3];
    outDim[0] = sz[0];
    outDim[1] = sz[1];
    outDim[2] = nint;

    mxClassID mDType = mxGetClassID(prhs[0]);
    if(mDType == mxSINGLE_CLASS){
        float* im = (float*)mxGetPr(prhs[0]);

        plhs[0] = mxCreateNumericArray(3, (mwSize*)outDim, mxSINGLE_CLASS, mxREAL);
        float* im_int = (float*)mxGetPr(plhs[0]);

        skewed_space_interp_defined_stepsize(im, im_int, xstep, stepsize, Reverse, sz, outDim);
    }
    else if(mDType == mxUINT16_CLASS){
        uint16_t* im = (uint16_t*)mxGetPr(prhs[0]);

        plhs[0] = mxCreateNumericArray(3, (mwSize*)outDim, mxUINT16_CLASS, mxREAL);
        uint16_t* im_int = (uint16_t*)mxGetPr(plhs[0]);

        skewed_space_interp_defined_stepsize_16bit(im, im_int, xstep, stepsize, Reverse, sz, outDim);
    }
    else{
        mexErrMsgIdAndTxt("feather_blending:dataTypeError","Data type not suppported");
    }
}
