// Weighted sums of 2 or 4 lines, the inner loops of the deskew/rotate kernels.
//
// Moved from the three deskew/rotate mex files (their copies of these functions are
// identical). With AVX2 (x86-64, built with -mavx2 as the mex) the SIMD code of the mex;
// elsewhere scalar code that computes the same operations in the same order, so that every
// platform gives the bits of the x86 mex:
// - float lines: ((wa*a + wb*b) + wc*c) + wd*d, left to right (the SIMD order is the same);
// - uint16 lines: the first floor(len / 16) * 16 values as (a*wa + b*wb) + (c*wc + d*wd), the
//   rest left to right;
// - uint16 output: + 0.5f, then truncated (in the SIMD part a saturating pack to [0, 65535]).
// The kernels are compiled without fused multiply-adds (-ffp-contract=off).
// The mex's nested "#pragma omp parallel for" and "#pragma parallel for" lines are left out:
// these functions only run inside the kernels' parallel loops, where they were serialized.
#pragma once

#include <cstddef>
#include <cstdint>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

namespace pk5d::deskew::ml {

#if defined(__AVX2__)

inline void compute_weighted_sum_two_vectors(const float* a, const float* b, float* result,
                                             float wa, float wb, size_t len) {
    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);

    // full blocks of 8 floats; the loop after this one handles the rest
    const size_t len8 = len - (len % 8);
    for (size_t i = 0; i < len8; i += 8) {
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);

        __m256 va_weighted = _mm256_mul_ps(va, wa_vec);
        __m256 vb_weighted = _mm256_mul_ps(vb, wb_vec);

        __m256 sum = _mm256_add_ps(va_weighted, vb_weighted);

        _mm256_storeu_ps(&result[i], sum);
    }

    // Handle any remaining elements (len not multiple of 8)
    for (size_t i = len - (len % 8); i < len; ++i) {
        float res = wa * a[i] + wb * b[i];
        result[i] = res;
    }
}

inline void compute_weighted_sum(const float* a, const float* b, const float* c, const float* d,
                                 float* result, float wa, float wb, float wc, float wd,
                                 size_t len) {
    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);
    const __m256 wc_vec = _mm256_set1_ps(wc);
    const __m256 wd_vec = _mm256_set1_ps(wd);

    // full blocks of 8 floats; the loop after this one handles the rest
    const size_t len8 = len - (len % 8);
    for (size_t i = 0; i < len8; i += 8) {
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);
        __m256 vc = _mm256_loadu_ps(&c[i]);
        __m256 vd = _mm256_loadu_ps(&d[i]);

        __m256 va_weighted = _mm256_mul_ps(va, wa_vec);
        __m256 vb_weighted = _mm256_mul_ps(vb, wb_vec);
        __m256 vc_weighted = _mm256_mul_ps(vc, wc_vec);
        __m256 vd_weighted = _mm256_mul_ps(vd, wd_vec);

        __m256 sum = _mm256_add_ps(va_weighted, vb_weighted);
        sum = _mm256_add_ps(sum, vc_weighted);
        sum = _mm256_add_ps(sum, vd_weighted);

        _mm256_storeu_ps(&result[i], sum);
    }

    // Handle any remaining elements (len not multiple of 8)
    for (size_t i = len - (len % 8); i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i];
        result[i] = res;
    }
}

inline void copy_uint16_to_float(const uint16_t* src, float* dst, size_t len) {
    for (size_t i = 0; i + 16 <= len; i += 16) {
        // Load 16 uint16_t values (two 128-bit chunks)
        __m256i src_chunk = _mm256_loadu_si256((const __m256i*)&src[i]);

        // Convert uint16_t to uint32_t and then to float
        __m256i src_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(src_chunk));
        __m256i src_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(src_chunk, 1));

        __m256 src_float_low = _mm256_cvtepi32_ps(src_low);
        __m256 src_float_high = _mm256_cvtepi32_ps(src_high);

        _mm256_storeu_ps(&dst[i], src_float_low);
        _mm256_storeu_ps(&dst[i + 8], src_float_high);
    }

    // Handle remaining elements (if len is not a multiple of 16)
    for (size_t i = len - (len % 16); i < len; ++i) {
        dst[i] = (float)src[i];
    }
}

// a, c: uint16; b, d: float
inline void compute_weighted_sum_16bit_float_to_float(const uint16_t* a, const float* b,
                                                      const uint16_t* c, const float* d,
                                                      float* result, float wa, float wb,
                                                      float wc, float wd, size_t len) {
    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);
    const __m256 wc_vec = _mm256_set1_ps(wc);
    const __m256 wd_vec = _mm256_set1_ps(wd);

    for (size_t i = 0; i + 16 <= len; i += 16) {
        // Load 16 uint16_t values (two 128-bit chunks)
        __m256i a_chunk = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i c_chunk = _mm256_loadu_si256((const __m256i*)&c[i]);

        // Convert uint16_t to uint32_t and then to float
        __m256i a_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(a_chunk));
        __m256i a_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(a_chunk, 1));
        __m256i c_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(c_chunk));
        __m256i c_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(c_chunk, 1));

        __m256 a_float_low = _mm256_cvtepi32_ps(a_low);
        __m256 a_float_high = _mm256_cvtepi32_ps(a_high);
        __m256 c_float_low = _mm256_cvtepi32_ps(c_low);
        __m256 c_float_high = _mm256_cvtepi32_ps(c_high);

        __m256 b_float_low = _mm256_loadu_ps(&b[i]);
        __m256 b_float_high = _mm256_loadu_ps(&b[i + 8]);
        __m256 d_float_low = _mm256_loadu_ps(&d[i]);
        __m256 d_float_high = _mm256_loadu_ps(&d[i + 8]);

        // Perform weighted sum: wa * a + wb * b + wc * c + wd * d
        __m256 result_low =
            _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(a_float_low, wa_vec), _mm256_mul_ps(b_float_low, wb_vec)),
                          _mm256_add_ps(_mm256_mul_ps(c_float_low, wc_vec), _mm256_mul_ps(d_float_low, wd_vec)));

        __m256 result_high =
            _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(a_float_high, wa_vec), _mm256_mul_ps(b_float_high, wb_vec)),
                          _mm256_add_ps(_mm256_mul_ps(c_float_high, wc_vec), _mm256_mul_ps(d_float_high, wd_vec)));

        _mm256_storeu_ps(&result[i], result_low);
        _mm256_storeu_ps(&result[i + 8], result_high);
    }

    // Handle remaining elements (if len is not a multiple of 16)
    for (size_t i = len - (len % 16); i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i];
        result[i] = res;
    }
}

inline void compute_weighted_sum_two_vectors_16bit_to_float(const uint16_t* a, const uint16_t* b,
                                                            float* result, float wa, float wb,
                                                            size_t len) {
    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);

    for (size_t i = 0; i + 16 <= len; i += 16) {
        // Load 16 uint16_t values (two 128-bit chunks)
        __m256i a_chunk = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i b_chunk = _mm256_loadu_si256((const __m256i*)&b[i]);

        // Convert uint16_t to uint32_t and then to float
        __m256i a_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(a_chunk));
        __m256i a_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(a_chunk, 1));
        __m256i b_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(b_chunk));
        __m256i b_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(b_chunk, 1));

        __m256 a_float_low = _mm256_cvtepi32_ps(a_low);
        __m256 a_float_high = _mm256_cvtepi32_ps(a_high);
        __m256 b_float_low = _mm256_cvtepi32_ps(b_low);
        __m256 b_float_high = _mm256_cvtepi32_ps(b_high);

        // Perform weighted sum: wa * a + wb * b
        __m256 result_low = _mm256_add_ps(_mm256_mul_ps(a_float_low, wa_vec), _mm256_mul_ps(b_float_low, wb_vec));
        __m256 result_high = _mm256_add_ps(_mm256_mul_ps(a_float_high, wa_vec), _mm256_mul_ps(b_float_high, wb_vec));

        _mm256_storeu_ps(&result[i], result_low);
        _mm256_storeu_ps(&result[i + 8], result_high);
    }

    // Handle remaining elements (if len is not a multiple of 16)
    for (size_t i = len - (len % 16); i < len; ++i) {
        float res = wa * a[i] + wb * b[i];
        result[i] = res;
    }
}

inline void compute_weighted_sum_16bit_to_float(const uint16_t* a, const uint16_t* b,
                                                const uint16_t* c, const uint16_t* d,
                                                float* result, float wa, float wb, float wc,
                                                float wd, size_t len) {
    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);
    const __m256 wc_vec = _mm256_set1_ps(wc);
    const __m256 wd_vec = _mm256_set1_ps(wd);

    for (size_t i = 0; i + 16 <= len; i += 16) {
        // Load 16 uint16_t values (two 128-bit chunks)
        __m256i a_chunk = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i b_chunk = _mm256_loadu_si256((const __m256i*)&b[i]);
        __m256i c_chunk = _mm256_loadu_si256((const __m256i*)&c[i]);
        __m256i d_chunk = _mm256_loadu_si256((const __m256i*)&d[i]);

        // Convert uint16_t to uint32_t and then to float
        __m256i a_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(a_chunk));
        __m256i a_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(a_chunk, 1));
        __m256i b_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(b_chunk));
        __m256i b_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(b_chunk, 1));
        __m256i c_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(c_chunk));
        __m256i c_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(c_chunk, 1));
        __m256i d_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(d_chunk));
        __m256i d_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(d_chunk, 1));

        __m256 a_float_low = _mm256_cvtepi32_ps(a_low);
        __m256 a_float_high = _mm256_cvtepi32_ps(a_high);
        __m256 b_float_low = _mm256_cvtepi32_ps(b_low);
        __m256 b_float_high = _mm256_cvtepi32_ps(b_high);
        __m256 c_float_low = _mm256_cvtepi32_ps(c_low);
        __m256 c_float_high = _mm256_cvtepi32_ps(c_high);
        __m256 d_float_low = _mm256_cvtepi32_ps(d_low);
        __m256 d_float_high = _mm256_cvtepi32_ps(d_high);

        // Perform weighted sum: wa * a + wb * b + wc * c + wd * d
        __m256 result_low =
            _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(a_float_low, wa_vec), _mm256_mul_ps(b_float_low, wb_vec)),
                          _mm256_add_ps(_mm256_mul_ps(c_float_low, wc_vec), _mm256_mul_ps(d_float_low, wd_vec)));

        __m256 result_high =
            _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(a_float_high, wa_vec), _mm256_mul_ps(b_float_high, wb_vec)),
                          _mm256_add_ps(_mm256_mul_ps(c_float_high, wc_vec), _mm256_mul_ps(d_float_high, wd_vec)));

        _mm256_storeu_ps(&result[i], result_low);
        _mm256_storeu_ps(&result[i + 8], result_high);
    }

    // Handle remaining elements (if len is not a multiple of 16)
    for (size_t i = len - (len % 16); i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i];
        result[i] = res;
    }
}

inline void compute_weighted_sum_two_vectors_16bit(const uint16_t* a, const uint16_t* b,
                                                   uint16_t* result, float wa, float wb,
                                                   size_t len) {
    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);
    const __m256 half_vec = _mm256_set1_ps(0.5f);

    for (size_t i = 0; i + 16 <= len; i += 16) {
        // Load 16 uint16_t values (two 128-bit chunks)
        __m256i a_chunk = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i b_chunk = _mm256_loadu_si256((const __m256i*)&b[i]);

        // Convert uint16_t to uint32_t and then to float
        __m256i a_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(a_chunk));
        __m256i a_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(a_chunk, 1));
        __m256i b_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(b_chunk));
        __m256i b_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(b_chunk, 1));

        __m256 a_float_low = _mm256_cvtepi32_ps(a_low);
        __m256 a_float_high = _mm256_cvtepi32_ps(a_high);
        __m256 b_float_low = _mm256_cvtepi32_ps(b_low);
        __m256 b_float_high = _mm256_cvtepi32_ps(b_high);

        // Perform weighted sum: wa * a + wb * b + 0.5
        __m256 result_low = _mm256_add_ps(_mm256_mul_ps(a_float_low, wa_vec), _mm256_mul_ps(b_float_low, wb_vec));
        result_low = _mm256_add_ps(result_low, half_vec);

        __m256 result_high = _mm256_add_ps(_mm256_mul_ps(a_float_high, wa_vec), _mm256_mul_ps(b_float_high, wb_vec));
        result_high = _mm256_add_ps(result_high, half_vec);

        // Convert float back to uint32_t, truncating: with the 0.5 added above, this rounds to nearest
        __m256i result_int_low = _mm256_cvttps_epi32(result_low);
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

inline void compute_weighted_sum_16bit(const uint16_t* a, const uint16_t* b, const uint16_t* c,
                                       const uint16_t* d, uint16_t* result, float wa, float wb,
                                       float wc, float wd, size_t len) {
    // Load weights into SIMD registers
    const __m256 wa_vec = _mm256_set1_ps(wa);
    const __m256 wb_vec = _mm256_set1_ps(wb);
    const __m256 wc_vec = _mm256_set1_ps(wc);
    const __m256 wd_vec = _mm256_set1_ps(wd);
    const __m256 half_vec = _mm256_set1_ps(0.5f);

    for (size_t i = 0; i + 16 <= len; i += 16) {
        // Load 16 uint16_t values (two 128-bit chunks)
        __m256i a_chunk = _mm256_loadu_si256((const __m256i*)&a[i]);
        __m256i b_chunk = _mm256_loadu_si256((const __m256i*)&b[i]);
        __m256i c_chunk = _mm256_loadu_si256((const __m256i*)&c[i]);
        __m256i d_chunk = _mm256_loadu_si256((const __m256i*)&d[i]);

        // Convert uint16_t to uint32_t and then to float
        __m256i a_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(a_chunk));
        __m256i a_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(a_chunk, 1));
        __m256i b_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(b_chunk));
        __m256i b_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(b_chunk, 1));
        __m256i c_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(c_chunk));
        __m256i c_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(c_chunk, 1));
        __m256i d_low = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(d_chunk));
        __m256i d_high = _mm256_cvtepu16_epi32(_mm256_extracti128_si256(d_chunk, 1));

        __m256 a_float_low = _mm256_cvtepi32_ps(a_low);
        __m256 a_float_high = _mm256_cvtepi32_ps(a_high);
        __m256 b_float_low = _mm256_cvtepi32_ps(b_low);
        __m256 b_float_high = _mm256_cvtepi32_ps(b_high);
        __m256 c_float_low = _mm256_cvtepi32_ps(c_low);
        __m256 c_float_high = _mm256_cvtepi32_ps(c_high);
        __m256 d_float_low = _mm256_cvtepi32_ps(d_low);
        __m256 d_float_high = _mm256_cvtepi32_ps(d_high);

        // Perform weighted sum: wa * a + wb * b + wc * c + wd * d + 0.5
        __m256 result_low =
            _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(a_float_low, wa_vec), _mm256_mul_ps(b_float_low, wb_vec)),
                          _mm256_add_ps(_mm256_mul_ps(c_float_low, wc_vec), _mm256_mul_ps(d_float_low, wd_vec)));
        result_low = _mm256_add_ps(result_low, half_vec);

        __m256 result_high =
            _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(a_float_high, wa_vec), _mm256_mul_ps(b_float_high, wb_vec)),
                          _mm256_add_ps(_mm256_mul_ps(c_float_high, wc_vec), _mm256_mul_ps(d_float_high, wd_vec)));
        result_high = _mm256_add_ps(result_high, half_vec);

        // Convert float back to uint32_t, truncating: with the 0.5 added above, this rounds to nearest
        __m256i result_int_low = _mm256_cvttps_epi32(result_low);
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

#else  // no AVX2: the same operations in the same order, one value at a time

// _mm256_cvttps_epi32 then _mm256_packus_epi32: truncation to int32 (out of range and NaN:
// INT32_MIN), saturated to [0, 65535]
inline uint16_t simd_pack_uint16(float v) {
    int32_t i = (v < 2147483648.0f && v >= -2147483648.0f) ? (int32_t)v : INT32_MIN;
    return (uint16_t)(i < 0 ? 0 : (i > 65535 ? 65535 : i));
}

inline void compute_weighted_sum_two_vectors(const float* a, const float* b, float* result,
                                             float wa, float wb, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        float res = wa * a[i] + wb * b[i];
        result[i] = res;
    }
}

inline void compute_weighted_sum(const float* a, const float* b, const float* c, const float* d,
                                 float* result, float wa, float wb, float wc, float wd,
                                 size_t len) {
    for (size_t i = 0; i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i];
        result[i] = res;
    }
}

inline void copy_uint16_to_float(const uint16_t* src, float* dst, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        dst[i] = (float)src[i];
    }
}

// a, c: uint16; b, d: float
inline void compute_weighted_sum_16bit_float_to_float(const uint16_t* a, const float* b,
                                                      const uint16_t* c, const float* d,
                                                      float* result, float wa, float wb,
                                                      float wc, float wd, size_t len) {
    const size_t len16 = len - (len % 16);
    for (size_t i = 0; i < len16; ++i) {
        float ab = (float)a[i] * wa + b[i] * wb;
        float cd = (float)c[i] * wc + d[i] * wd;
        result[i] = ab + cd;
    }
    for (size_t i = len16; i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i];
        result[i] = res;
    }
}

inline void compute_weighted_sum_two_vectors_16bit_to_float(const uint16_t* a, const uint16_t* b,
                                                            float* result, float wa, float wb,
                                                            size_t len) {
    for (size_t i = 0; i < len; ++i) {
        float res = wa * a[i] + wb * b[i];
        result[i] = res;
    }
}

inline void compute_weighted_sum_16bit_to_float(const uint16_t* a, const uint16_t* b,
                                                const uint16_t* c, const uint16_t* d,
                                                float* result, float wa, float wb, float wc,
                                                float wd, size_t len) {
    const size_t len16 = len - (len % 16);
    for (size_t i = 0; i < len16; ++i) {
        float ab = (float)a[i] * wa + (float)b[i] * wb;
        float cd = (float)c[i] * wc + (float)d[i] * wd;
        result[i] = ab + cd;
    }
    for (size_t i = len16; i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i];
        result[i] = res;
    }
}

inline void compute_weighted_sum_two_vectors_16bit(const uint16_t* a, const uint16_t* b,
                                                   uint16_t* result, float wa, float wb,
                                                   size_t len) {
    const size_t len16 = len - (len % 16);
    for (size_t i = 0; i < len16; ++i) {
        float res = (float)a[i] * wa + (float)b[i] * wb;
        result[i] = simd_pack_uint16(res + 0.5f);
    }
    for (size_t i = len16; i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + 0.5f;
        result[i] = (uint16_t)res;
    }
}

inline void compute_weighted_sum_16bit(const uint16_t* a, const uint16_t* b, const uint16_t* c,
                                       const uint16_t* d, uint16_t* result, float wa, float wb,
                                       float wc, float wd, size_t len) {
    const size_t len16 = len - (len % 16);
    for (size_t i = 0; i < len16; ++i) {
        float ab = (float)a[i] * wa + (float)b[i] * wb;
        float cd = (float)c[i] * wc + (float)d[i] * wd;
        result[i] = simd_pack_uint16((ab + cd) + 0.5f);
    }
    for (size_t i = len16; i < len; ++i) {
        float res = wa * a[i] + wb * b[i] + wc * c[i] + wd * d[i] + 0.5f;
        result[i] = (uint16_t)res;
    }
}

#endif

}  // namespace pk5d::deskew::ml
