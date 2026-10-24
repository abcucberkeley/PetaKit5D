#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <inttypes.h>
#include <string.h>
#ifdef __x86_64__
#include <immintrin.h> // For AVX/AVX2
#endif
#include "mex.h"

// mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp -mavx2' LDFLAGS='$LDFLAGS -O3 -fopenmp -mavx2' skewed_space_interp_volume_deskew_rotate_warp_mex.cpp
// macOS
// mex -v CXX="/usr/local/bin/g++-13" CXXFLAGS='-fno-common -arch x86_64 -mmacosx-version-min=10.15 -fexceptions -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -std=c++11 -O3 -fopenmp -DMATLAB_DEFAULT_RELEASE=R2017b  -DUSE_MEX_CMD   -DMATLAB_MEX_FILE' CXXOPTIMFLAGS="-O3 -DNDEBUG" LDFLAGS='$LDFLAGS -O3 -fopenmp' skewed_space_interp_volume_deskew_rotate_warp_mex.cpp

// combine skewed space interpolation together with volume deskew/rotation, without copying the original data.
// The interpolated slices are not stored: the rows of them each output z slice needs are computed
// on the fly into a small per-thread cache (see transformVolume).


#ifndef __x86_64__

void compute_weighted_sum_two_vectors_sequential(const float* a, const float* b, float* result, 
                                      float wa, float wb, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        float res = wa * a[i] + wb * b[i];
        result[i] = res;
    }
}


void compute_weighted_sum_sequential(const float* a, const float* b, const float* c, const float* d, 
                          float* result, float wa, float wb, float wc, float wd, size_t len) {

    for (size_t i = 0; i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i];
        result[i] = res;
    }
}

// a, c: uint16; b, d: float
void compute_weighted_sum_sequential_16bit_float_to_float(const uint16_t* a, const float* b, const uint16_t* c, const float* d,
                                                          float* result, float wa, float wb, float wc, float wd, size_t len) {

    #pragma parallel for
    for (size_t i = 0; i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i];
        result[i] = res;
    }
}


void compute_weighted_sum_two_vectors_sequential_16bit_to_float(const uint16_t* a, const uint16_t* b, float* result, float wa, float wb, size_t len) {

    #pragma parallel for
    for (size_t i = 0; i < len; ++i) {
        float res = wa * a[i] + wb * b[i];
        result[i] = res;
    }
}


void compute_weighted_sum_sequential_16bit_to_float(const uint16_t* a, const uint16_t* b, const uint16_t* c, const uint16_t* d,
                                   float* result, float wa, float wb, float wc, float wd, size_t len) {

    #pragma parallel for
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

    // full blocks of 8 floats; the loop after this one handles the rest
    const size_t len8 = len - (len % 8);
    for (size_t i = 0; i < len8; i += 8) {
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

    // full blocks of 8 floats; the loop after this one handles the rest
    const size_t len8 = len - (len % 8);
    for (size_t i = 0; i < len8; i += 8) {
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

// a, c: uint16; b, d: float
void compute_weighted_sum_16bit_float_to_float(const uint16_t* a, const float* b, const uint16_t* c, const float* d,
                                               float* result, float wa, float wb, float wc, float wd, size_t len) {
    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);
    const __m256 wc_vec = _mm256_set1_ps(wc);
    const __m256 wd_vec = _mm256_set1_ps(wd);

    #pragma parallel for
    for (size_t i = 0; i + 16 <= len; i += 16) {
        // Load 16 uint16_t values (two 128-bit chunks)
        __m256i a_chunk = _mm256_loadu_si256((__m256i*)&a[i]);
        __m256i c_chunk = _mm256_loadu_si256((__m256i*)&c[i]);

        // Convert uint16_t to uint32_t and then to float
        __m256i a_low  = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(a_chunk));
        __m256i a_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(a_chunk, 1));
        __m256i c_low  = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(c_chunk));
        __m256i c_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(c_chunk, 1));

        __m256 a_float_low  = _mm256_cvtepi32_ps(a_low);
        __m256 a_float_high = _mm256_cvtepi32_ps(a_high);
        __m256 c_float_low  = _mm256_cvtepi32_ps(c_low);
        __m256 c_float_high = _mm256_cvtepi32_ps(c_high);

        __m256 b_float_low  = _mm256_loadu_ps(&b[i]);
        __m256 b_float_high = _mm256_loadu_ps(&b[i+8]);        
        __m256 d_float_low  = _mm256_loadu_ps(&d[i]);
        __m256 d_float_high = _mm256_loadu_ps(&d[i+8]);

        // Perform weighted sum: wa * a + wb * b + wc * c + wd * d
        __m256 result_low  = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(a_float_low, wa_vec), _mm256_mul_ps(b_float_low, wb_vec)),
                                           _mm256_add_ps(_mm256_mul_ps(c_float_low, wc_vec), _mm256_mul_ps(d_float_low, wd_vec)));

        __m256 result_high = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(a_float_high, wa_vec), _mm256_mul_ps(b_float_high, wb_vec)),
                                           _mm256_add_ps(_mm256_mul_ps(c_float_high, wc_vec), _mm256_mul_ps(d_float_high, wd_vec)));

        // Store the result
        _mm256_storeu_ps(&result[i], result_low);
        _mm256_storeu_ps(&result[i + 8], result_high);
    }

    // Handle remaining elements (if len is not a multiple of 16)
    for (size_t i = len - (len % 16); i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i];
        result[i] = res;
    }
}


void compute_weighted_sum_two_vectors_16bit_to_float(const uint16_t* a, const uint16_t* b, float* result, float wa, float wb, size_t len) {

    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);

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

        // Perform weighted sum: wa * a + wb * b
        __m256 result_low  = _mm256_add_ps(_mm256_mul_ps(a_float_low, wa_vec), _mm256_mul_ps(b_float_low, wb_vec));
        __m256 result_high = _mm256_add_ps(_mm256_mul_ps(a_float_high, wa_vec), _mm256_mul_ps(b_float_high, wb_vec));

        // Store the result
        _mm256_storeu_ps(&result[i], result_low);
        _mm256_storeu_ps(&result[i + 8], result_high);
    }

    // Handle remaining elements (if len is not a multiple of 16)
    for (size_t i = len - (len % 16); i < len; ++i) {
        float res = wa * a[i] + wb * b[i];
        result[i] = res;
    }
}


void compute_weighted_sum_16bit_to_float(const uint16_t* a, const uint16_t* b, const uint16_t* c, const uint16_t* d,
                                   float* result, float wa, float wb, float wc, float wd, size_t len) {
    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);
    const __m256 wc_vec = _mm256_set1_ps(wc);
    const __m256 wd_vec = _mm256_set1_ps(wd);

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

        // Perform weighted sum: wa * a + wb * b + wc * c + wd * d
        __m256 result_low  = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(a_float_low, wa_vec), _mm256_mul_ps(b_float_low, wb_vec)),
                                           _mm256_add_ps(_mm256_mul_ps(c_float_low, wc_vec), _mm256_mul_ps(d_float_low, wd_vec)));

        __m256 result_high = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(a_float_high, wa_vec), _mm256_mul_ps(b_float_high, wb_vec)),
                                           _mm256_add_ps(_mm256_mul_ps(c_float_high, wc_vec), _mm256_mul_ps(d_float_high, wd_vec)));

        // Store the result
        _mm256_storeu_ps(&result[i], result_low);
        _mm256_storeu_ps(&result[i + 8], result_high);
    }

    // Handle remaining elements (if len is not a multiple of 16)
    for (size_t i = len - (len % 16); i < len; ++i) {
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

        // Perform weighted sum: wa * a + wb * b + 0.5
        __m256 result_low  = _mm256_add_ps(_mm256_mul_ps(a_float_low, wa_vec), _mm256_mul_ps(b_float_low, wb_vec));
        result_low = _mm256_add_ps(result_low, half_vec);

        __m256 result_high = _mm256_add_ps(_mm256_mul_ps(a_float_high, wa_vec), _mm256_mul_ps(b_float_high, wb_vec));
        result_high = _mm256_add_ps(result_high, half_vec);

        // Convert float back to uint32_t, truncating: with the 0.5 added above, this rounds to nearest
        __m256i result_int_low  = _mm256_cvttps_epi32(result_low);
        __m256i result_int_high = _mm256_cvttps_epi32(result_high);

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

        // Convert float back to uint32_t, truncating: with the 0.5 added above, this rounds to nearest
        __m256i result_int_low  = _mm256_cvttps_epi32(result_low);
        __m256i result_int_high = _mm256_cvttps_epi32(result_high);

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


// line sums by row types: the SIMD versions on x86, the sequential loops on ARM (macOS)

static inline void weighted_sum_2(const float* a, const float* b, float* result, float wa, float wb, size_t len) {
    #ifndef __APPLE__
    compute_weighted_sum_two_vectors(a, b, result, wa, wb, len);
    #else
    compute_weighted_sum_two_vectors_sequential(a, b, result, wa, wb, len);
    #endif
}


static inline void weighted_sum_4(const float* a, const float* b, const float* c, const float* d,
                                  float* result, float wa, float wb, float wc, float wd, size_t len) {
    #ifndef __APPLE__
    compute_weighted_sum(a, b, c, d, result, wa, wb, wc, wd, len);
    #else
    compute_weighted_sum_sequential(a, b, c, d, result, wa, wb, wc, wd, len);
    #endif
}


static inline void weighted_sum_2(const uint16_t* a, const uint16_t* b, float* result, float wa, float wb, size_t len) {
    #ifndef __APPLE__
    compute_weighted_sum_two_vectors_16bit_to_float(a, b, result, wa, wb, len);
    #else
    compute_weighted_sum_two_vectors_sequential_16bit_to_float(a, b, result, wa, wb, len);
    #endif
}


static inline void weighted_sum_4(const uint16_t* a, const uint16_t* b, const uint16_t* c, const uint16_t* d,
                                  float* result, float wa, float wb, float wc, float wd, size_t len) {
    #ifndef __APPLE__
    compute_weighted_sum_16bit_to_float(a, b, c, d, result, wa, wb, wc, wd, len);
    #else
    compute_weighted_sum_sequential_16bit_to_float(a, b, c, d, result, wa, wb, wc, wd, len);
    #endif
}


// a, c: uint16; b, d: float
static inline void weighted_sum_4(const uint16_t* a, const float* b, const uint16_t* c, const float* d,
                                  float* result, float wa, float wb, float wc, float wd, size_t len) {
    #ifndef __APPLE__
    compute_weighted_sum_16bit_float_to_float(a, b, c, d, result, wa, wb, wc, wd, len);
    #else
    compute_weighted_sum_sequential_16bit_float_to_float(a, b, c, d, result, wa, wb, wc, wd, len);
    #endif
}


static inline void weighted_sum_2(const uint16_t* a, const uint16_t* b, uint16_t* result, float wa, float wb, size_t len) {
    #ifndef __APPLE__
    compute_weighted_sum_two_vectors_16bit(a, b, result, wa, wb, len);
    #else
    compute_weighted_sum_two_vectors_sequential_16bit(a, b, result, wa, wb, len);
    #endif
}


static inline void weighted_sum_4(const uint16_t* a, const uint16_t* b, const uint16_t* c, const uint16_t* d,
                                  uint16_t* result, float wa, float wb, float wc, float wd, size_t len) {
    #ifndef __APPLE__
    compute_weighted_sum_16bit(a, b, c, d, result, wa, wb, wc, wd, len);
    #else
    compute_weighted_sum_sequential_16bit(a, b, c, d, result, wa, wb, wc, wd, len);
    #endif
}


// an interpolated slice: between raw slices zorg and zorg + 1, rows shifted by s and t
struct InterpSlice {
    uint64_t zorg;
    uint64_t s;
    uint64_t t;
    float c00_w;
    float c01_w;
    float c10_w;
    float c11_w;
};


// parameters of the interpolated slices (indexed by interp_index)
void skewed_space_interp_slices(const double xstep, const double stepsize, const uint8_t Reverse, const uint64_t* const &interp_index,
                                const bool* const &is_original, const uint64_t ntotal, InterpSlice* const &slices){

    for(uint64_t i = 0; i < ntotal; i++){
        if (is_original[i]) continue;

        double s;
        double t;
        double zout;
        double zw;
        uint64_t zint = interp_index[i];

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

        // single precision distances, weights computed in double
        float swf = sw;
        float twf = tw;
        float zwf = zw;

        InterpSlice &p = slices[zint];
        p.zorg = floor(zout);
        p.s = (uint64_t)s;
        p.t = (uint64_t)t;
        p.c00_w = (1.0 - swf) * (1.0 - zwf);
        p.c01_w = swf * (1.0 - zwf);
        p.c10_w = (1.0 - twf) * zwf;
        p.c11_w = twf * zwf;
    }
}


// row y (all sz[0] values) of an interpolated slice
template <typename T, typename TI>
void skewed_space_interp_row(const T* const &im, TI* const &row, const InterpSlice &p, const uint64_t y, const uint64_t xa,
                             const uint8_t Reverse, const uint64_t sz[3]){

    const uint64_t nxy = sz[0] * sz[1];
    const uint64_t zorg = p.zorg;
    const uint64_t s = p.s;
    const uint64_t t = p.t;

    if(Reverse){
        if(y<xa-s){
            weighted_sum_2(im + (y+t) * sz[0] + (zorg+1) * nxy, im + (y+t+1) * sz[0] + (zorg+1) * nxy,
                           row, p.c10_w, p.c11_w, sz[0]);
        }else if(y < sz[1]-t-1){
            weighted_sum_4(im + (y+s-xa) * sz[0] + zorg * nxy, im + (y+s-xa+1) * sz[0] + zorg * nxy,
                           im + (y+t) * sz[0] + (zorg+1) * nxy, im + (y+t+1) * sz[0] + (zorg+1) * nxy,
                           row, p.c00_w, p.c01_w, p.c10_w, p.c11_w, sz[0]);
        } else{
            weighted_sum_2(im + (y+s-xa) * sz[0] + zorg * nxy, im + (y+s-xa+1) * sz[0] + zorg * nxy,
                           row, p.c00_w, p.c01_w, sz[0]);
        }
    } else {
        if(y<xa-s){
            weighted_sum_2(im + (y+s) * sz[0] + zorg * nxy, im + (y+s+1) * sz[0] + zorg * nxy,
                           row, p.c00_w, p.c01_w, sz[0]);
        }else if(y < sz[1]-t-1){
            weighted_sum_4(im + (y+s) * sz[0] + zorg * nxy, im + (y+s+1) * sz[0] + zorg * nxy,
                           im + (y+t-xa) * sz[0] + (zorg+1) * nxy, im + (y+t-xa+1) * sz[0] + (zorg+1) * nxy,
                           row, p.c00_w, p.c01_w, p.c10_w, p.c11_w, sz[0]);
        } else{
            weighted_sum_2(im + (y+t-xa) * sz[0] + (zorg+1) * nxy, im + (y+t-xa+1) * sz[0] + (zorg+1) * nxy,
                           row, p.c10_w, p.c11_w, sz[0]);
        }
    }
}


// output line from rows y0 and y0 + 1 (shapeX apart) of slices z0 (a) and z0 + 1 (b) of the
// interpolation grid; each slice is original (raw rows) or interpolated (cached rows: ai, bi)
static inline void transform_line(const float* a, const float* ai, const float* b, const float* bi, float* result,
                                  float w00, float w01, float w10, float w11, uint64_t shapeX, uint64_t len) {
    const float* dmat = a ? a : ai;
    const float* dpmat = b ? b : bi;
    weighted_sum_4(dmat, dpmat, dmat + shapeX, dpmat + shapeX, result, w00, w01, w10, w11, len);
}


static inline void transform_line(const uint16_t* a, const uint16_t* ai, const uint16_t* b, const uint16_t* bi, uint16_t* result,
                                  float w00, float w01, float w10, float w11, uint64_t shapeX, uint64_t len) {
    const uint16_t* dmat = a ? a : ai;
    const uint16_t* dpmat = b ? b : bi;
    weighted_sum_4(dmat, dpmat, dmat + shapeX, dpmat + shapeX, result, w00, w01, w10, w11, len);
}


// uint16 input, single output: interpolated rows are single
static inline void transform_line(const uint16_t* a, const float* ai, const uint16_t* b, const float* bi, float* result,
                                  float w00, float w01, float w10, float w11, uint64_t shapeX, uint64_t len) {
    if (a && b){
        weighted_sum_4(a, b, a + shapeX, b + shapeX, result, w00, w01, w10, w11, len);
    } else if (a) {
        weighted_sum_4(a, bi, a + shapeX, bi + shapeX, result, w00, w01, w10, w11, len);
    } else if (b) {
        weighted_sum_4(b, ai, b + shapeX, ai + shapeX, result, w01, w00, w11, w10, len);
    } else {
        weighted_sum_4(ai, bi, ai + shapeX, bi + shapeX, result, w00, w01, w10, w11, len);
    }
}


// Perform the transformation on a volume. Input x depends on output z only, so an output z slice
// reads two rows (y0, y0 + 1) of each slice of the interpolation grid; the two rows of the
// interpolated slices are computed when they are first needed and kept in a small per-thread
// cache (the slice index moves monotonically along output y), instead of interpolating whole
// slices into an intermediate volume first.
#define NCACHE 4

template <typename T, typename TI, typename TO>
void transformVolume(const T* const &im, const double (&tmat)[4][4], TO* const &nmat, const InterpSlice* const &slices,
                     const uint64_t* const &interp_index, const bool* const &is_original, const uint64_t xa, const uint8_t Reverse,
                     const uint64_t sz[3], const uint64_t &shapeZ, const float &startX, const float &startY, const float &startZ,
                     const uint64_t &outShapeX, const uint64_t &outShapeY, const uint64_t &outShapeZ, TI* const &cacheAll){

    const uint64_t shapeX = sz[0];
    const uint64_t shapeY = sz[1];
    const uint64_t shapeXY = shapeX * shapeY;
    const uint64_t outShapeXY = outShapeX * outShapeY;
    const uint64_t x0 = floor(startX);

    #pragma omp parallel
    {
        TI* cache = cacheAll + (uint64_t)omp_get_thread_num() * NCACHE * 2 * shapeX;
        int64_t tag[NCACHE];

        #pragma omp for schedule(dynamic, 1)
        for (uint64_t z = 0; z < outShapeZ; ++z){
            double ty;
            ty = tmat[1][2] * (z + startZ) + tmat[1][3];
            if (ty < 0 || ty + 1 >= shapeY){
                if (ty + 1 == shapeY) {
                    uint64_t counter = 1;
                    while (floor(ty) + 1 == shapeY){
                        ty -= 1e-9 * counter;
                        counter*=2;
                    }
                }
                else continue;
            }
            int64_t y0 = floor(ty);
            float wy = ty - y0;

            for (int k = 0; k < NCACHE; k++) tag[k] = -1;

            for (uint64_t y = 0; y < outShapeY; ++y){
                double tz;
                tz = tmat[2][1] * (y + startY) + tmat[2][2] * (z + startZ) + tmat[2][3];
                if (tz < 0 || tz + 1 >= shapeZ){
                    // handling of edge cases
                    if (tz + 1 == shapeZ) {
                        uint64_t counter = 1;
                        while (floor(tz) + 1 == shapeZ){
                            tz -= 1e-9 * counter;
                            counter*=2;
                        }
                    }
                    else continue;
                }
                uint64_t z0 = floor(tz);
                float wz = tz - z0;

                float w00 = (1 - wy) * (1 - wz);
                float w01 = (1 - wy) * wz;
                float w10 = wy * (1 - wz);
                float w11 = wy * wz;

                // rows y0 and y0 + 1 of slices z0 and z0 + 1, from x0 on
                const T* rows[2] = {NULL, NULL};
                const TI* cached[2] = {NULL, NULL};
                for (int q = 0; q < 2; q++){
                    const uint64_t zq = z0 + q;
                    if (is_original[zq]){
                        rows[q] = im + x0 + y0 * shapeX + interp_index[zq] * shapeXY;
                    } else {
                        const uint64_t zint = interp_index[zq];
                        const int k = zint % NCACHE;
                        TI* c = cache + k * 2 * shapeX;
                        if (tag[k] != (int64_t)zint){
                            skewed_space_interp_row(im, c, slices[zint], y0, xa, Reverse, sz);
                            skewed_space_interp_row(im, c + shapeX, slices[zint], y0 + 1, xa, Reverse, sz);
                            tag[k] = zint;
                        }
                        cached[q] = c + x0;
                    }
                }

                transform_line(rows[0], cached[0], rows[1], cached[1], nmat + y * outShapeX + z * outShapeXY,
                               w00, w01, w10, w11, shapeX, outShapeX);
            }
        }
    }
}


void determine_interpolation_indexing(const uint64_t n, const double stepSize, uint64_t* &interp_index, bool* &is_original, uint64_t &nint){
    uint64_t interp_counter = 0;
    
    for (uint64_t i = 0; i < n; i++){
        double cur_step = (double) i * stepSize;
        double cur_step_round = round(cur_step);
        if (fabs(cur_step - cur_step_round) < 1e-3){
            is_original[i] = true;
            interp_index[i] = (uint64_t) cur_step_round;
        } else {
            is_original[i] = false;
            interp_index[i] = interp_counter;
            interp_counter++;
        }
    }
    nint = interp_counter;    
}




void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]){
#ifdef __APPLE__
    // keep the mex file loaded until MATLAB exits: GCC's OpenMP runtime leaves
    // thread exit handlers in it, so unloading it crashes MATLAB on exit
    if (!mexIsLocked()) {
        mexLock();
    }
#endif

    if (nrhs != 7) mexErrMsgIdAndTxt("volume_warp:inputError","Number of input arguments must be 7");
    if (nlhs != 1) mexErrMsgIdAndTxt("volume_warp:outputError","Number of output arguments must be 1");

    const double xstep = (double)*mxGetPr(prhs[1]);
    const double stepsize = (double)*mxGetPr(prhs[2]);
    const uint8_t Reverse = (uint8_t)mxIsLogicalScalarTrue(prhs[3]);
    const uint8_t save16bit = (uint8_t)mxIsLogicalScalarTrue(prhs[6]);

    uint64_t* dimsA = (uint64_t*)mxGetDimensions(prhs[0]);
    uint64_t sz[3] = {1, 1, 1};
    uint64_t ndim = (uint64_t) mxGetNumberOfDimensions(prhs[0]);
    for(uint64_t t=0; t < ndim; t++){
        sz[t] = dimsA[t];
    }

    uint64_t shapeX = sz[0];

    uint64_t ntotal = (uint64_t) (floor(round((double)(sz[2] - 1) / stepsize * 100000) / 100000) + 1);

    // parameters for deskew/rotation
    if(mxGetNumberOfDimensions(prhs[4]) != 2 || mxGetM(prhs[4]) != 4 || mxGetN(prhs[4]) != 4) mexErrMsgIdAndTxt("volume_warp:inputError","Input transformation matrix must be 4x4!");
    double tmat[4][4];
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            tmat[i][j] = (double)*(mxGetPr(prhs[4])+i+j*4);
        }
    }

    if(mxGetN(prhs[5]) != 6) mexErrMsgIdAndTxt("volume_warp:inputError","Input range for bbox is not 6");

    float startX = (float)(*mxGetPr(prhs[5]))-1;
    float startY = (float)(*(mxGetPr(prhs[5]) + 1))-1;
    float startZ = (float)(*(mxGetPr(prhs[5]) + 2))-1;
    float endX = (float)(*(mxGetPr(prhs[5]) + 3));
    float endY = (float)(*(mxGetPr(prhs[5]) + 4));
    float endZ = (float)(*(mxGetPr(prhs[5]) + 5));

    if (startX < 0 || startY < 0 || startZ < 0 || startX > endX || startY > endY || startZ > endZ)
        mexErrMsgIdAndTxt("volume_warp:inputError","Input bbox is invalid, the numbers must be nonnegative and the end ranges must be nonless than the start ranges.");

    uint64_t outShapeX = (uint64_t) (endX - startX);
    uint64_t outShapeY = (uint64_t) (endY - startY);
    uint64_t outShapeZ = (uint64_t) (endZ - startZ);

    // output lines are read from x0 on in rows of the input (and of the cached interpolated rows)
    if ((uint64_t)floor(startX) + outShapeX > shapeX)
        mexErrMsgIdAndTxt("volume_warp:inputError","Input bbox is invalid, its first range is outside the input.");

    uint64_t outDim[3];
    outDim[0] = outShapeX;
    outDim[1] = outShapeY;
    outDim[2] = outShapeZ;

    mxClassID mDType = mxGetClassID(prhs[0]);
    if(mDType != mxSINGLE_CLASS && mDType != mxUINT16_CLASS){
        mexErrMsgIdAndTxt("feather_blending:dataTypeError","Data type not suppported");
    }

    // determine which slices come from original and which come from interpolated
    uint64_t* interp_index = (uint64_t*)malloc((ntotal)*sizeof(uint64_t));
    bool* is_original = (bool*)malloc((ntotal)*sizeof(bool));
    uint64_t nint = 0;

    determine_interpolation_indexing(ntotal, stepsize, interp_index, is_original, nint);

    InterpSlice* slices = (InterpSlice*)mxMalloc((nint > 0 ? nint : 1) * sizeof(InterpSlice));
    skewed_space_interp_slices(xstep, stepsize, Reverse, interp_index, is_original, ntotal, slices);
    const uint64_t xa = (uint64_t)ceil(xstep);

    // cached rows of interpolated slices, per thread
    const uint64_t cacheSize = (uint64_t)omp_get_max_threads() * NCACHE * 2 * shapeX;

    if(mDType == mxSINGLE_CLASS){
        float* im = (float*)mxGetPr(prhs[0]);

        float* cache = (float*)mxMalloc(cacheSize * sizeof(float));

        plhs[0] = mxCreateNumericArray(3, (mwSize*)outDim, mxSINGLE_CLASS, mxREAL);
        float* nmat = (float*)mxGetPr(plhs[0]);

        transformVolume(im, tmat, nmat, slices, interp_index, is_original, xa, Reverse, sz, ntotal, startX, startY, startZ, outShapeX, outShapeY, outShapeZ, cache);

        mxFree(cache);
    }
    else {
        uint16_t* im = (uint16_t*)mxGetPr(prhs[0]);

        if (save16bit){
            uint16_t* cache = (uint16_t*)mxMalloc(cacheSize * sizeof(uint16_t));

            plhs[0] = mxCreateNumericArray(3, (mwSize*)outDim, mxUINT16_CLASS, mxREAL);
            uint16_t* nmat = (uint16_t*)mxGetPr(plhs[0]);

            transformVolume(im, tmat, nmat, slices, interp_index, is_original, xa, Reverse, sz, ntotal, startX, startY, startZ, outShapeX, outShapeY, outShapeZ, cache);

            mxFree(cache);
        } else {
            float* cache = (float*)mxMalloc(cacheSize * sizeof(float));

            plhs[0] = mxCreateNumericArray(3, (mwSize*)outDim, mxSINGLE_CLASS, mxREAL);
            float* nmat = (float*)mxGetPr(plhs[0]);

            transformVolume(im, tmat, nmat, slices, interp_index, is_original, xa, Reverse, sz, ntotal, startX, startY, startZ, outShapeX, outShapeY, outShapeZ, cache);

            mxFree(cache);
        }
    }

    mxFree(slices);
    free(interp_index);
    free(is_original);
}
