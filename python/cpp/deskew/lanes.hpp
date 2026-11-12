// The per-voxel arithmetic of the deskew/rotate mex files, on 8 lanes.
//
// The mex sum lines of MATLAB's first dimension (y) with AVX2 and a scalar tail; every voxel's
// value depends only on its operands, its weights and its position in the line. These are the
// rules (the products commute exactly, so only the order of the additions matters):
//   float lines:          ((wa*a + wb*b) + wc*c) + wd*d                       (SEQ)
//   uint16 lines (also with single operands), at positions below floor(len / 16) * 16:
//                         (a*wa + b*wb) + (c*wc + d*wd)                       (PAIR)
//                         and SEQ after
//   uint16 to single, an interpolated then an original slice (the mex swaps the operands):
//                         PAIR, and after it ((wa*a + wb*b) + wd*d) + wc*c    (SEQ_SWAP)
//   two lines:            wa*a + wb*b
//   uint16 output:        + 0.5f, truncated, saturated to [0, 65535] (the SIMD pack; the
//                         scalar tail's (uint16_t) cast gives the same for the values that occur)
// uint16 operands are converted to single exactly, so they are kept as single here. Lanes
// are 8 neighbouring voxels of any axis: the kernels choose the axis that is contiguous in
// their scratch buffers. With AVX2 one __m256 per vector; elsewhere 8 floats, the same
// operations in the same order (the kernels are compiled without fused multiply-adds).
#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

namespace pk5d::deskew {


// float to uint16 as _mm256_cvttps_epi32 then _mm256_packus_epi32 do it: truncated to int32
// (out of range and NaN: INT32_MIN), saturated to [0, 65535]
inline uint16_t sat_u16(float v) {
    const int32_t k = (v < 2147483648.0f && v >= -2147483648.0f) ? (int32_t)v : INT32_MIN;
    return (uint16_t)(k < 0 ? 0 : (k > 65535 ? 65535 : k));
}

// an output value from its single precision value (uint16: truncated and saturated)
template <typename TO>
inline TO to_out(float v);
template <>
inline float to_out<float>(float v) {
    return v;
}
template <>
inline uint16_t to_out<uint16_t>(float v) {
    return sat_u16(v);
}

#if defined(__AVX2__)

struct V8 {
    __m256 v;
};

inline V8 load8(const float* p) { return {_mm256_loadu_ps(p)}; }
inline void store8(float* p, V8 a) { _mm256_storeu_ps(p, a.v); }
inline V8 set8(float w) { return {_mm256_set1_ps(w)}; }
inline V8 zero8() { return {_mm256_setzero_ps()}; }
inline V8 operator+(V8 a, V8 b) { return {_mm256_add_ps(a.v, b.v)}; }
inline V8 operator*(V8 a, V8 b) { return {_mm256_mul_ps(a.v, b.v)}; }
// lanes of a where mask lane i is set (i < n), else b
inline V8 first_lanes(V8 a, V8 b, int n) {
    const __m256i idx = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256 m = _mm256_castsi256_ps(_mm256_cmpgt_epi32(_mm256_set1_epi32(n), idx));
    return {_mm256_blendv_ps(b.v, a.v, m)};
}

// 8 uint16 values as single
inline V8 load8_u16(const uint16_t* p) {
    __m128i x = _mm_loadu_si128((const __m128i*)p);
    return {_mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(x))};
}

// + 0.5f, truncate, saturate (the mex's _mm256_cvttps_epi32 and _mm256_packus_epi32)
inline V8 round_u16(V8 a) {
    __m256i i = _mm256_cvttps_epi32(_mm256_add_ps(a.v, _mm256_set1_ps(0.5f)));
    i = _mm256_max_epi32(_mm256_min_epi32(i, _mm256_set1_epi32(65535)), _mm256_setzero_si256());
    return {_mm256_cvtepi32_ps(i)};
}

// 8 values as uint16: truncated, saturated (sat_u16)
inline void store8_u16(uint16_t* p, V8 a) {
    __m256i i = _mm256_cvttps_epi32(a.v);
    __m128i lo = _mm256_castsi256_si128(i);
    __m128i hi = _mm256_extracti128_si256(i, 1);
    _mm_storeu_si128((__m128i*)p, _mm_packus_epi32(lo, hi));
}

// 8 x 8 transpose of single values: rows r[0..7] in, columns out
inline void transpose8(__m256& r0, __m256& r1, __m256& r2, __m256& r3, __m256& r4, __m256& r5,
                       __m256& r6, __m256& r7) {
    __m256 t0 = _mm256_unpacklo_ps(r0, r1), t1 = _mm256_unpackhi_ps(r0, r1);
    __m256 t2 = _mm256_unpacklo_ps(r2, r3), t3 = _mm256_unpackhi_ps(r2, r3);
    __m256 t4 = _mm256_unpacklo_ps(r4, r5), t5 = _mm256_unpackhi_ps(r4, r5);
    __m256 t6 = _mm256_unpacklo_ps(r6, r7), t7 = _mm256_unpackhi_ps(r6, r7);
    __m256 s0 = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(1, 0, 1, 0));
    __m256 s1 = _mm256_shuffle_ps(t0, t2, _MM_SHUFFLE(3, 2, 3, 2));
    __m256 s2 = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 0, 1, 0));
    __m256 s3 = _mm256_shuffle_ps(t1, t3, _MM_SHUFFLE(3, 2, 3, 2));
    __m256 s4 = _mm256_shuffle_ps(t4, t6, _MM_SHUFFLE(1, 0, 1, 0));
    __m256 s5 = _mm256_shuffle_ps(t4, t6, _MM_SHUFFLE(3, 2, 3, 2));
    __m256 s6 = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(1, 0, 1, 0));
    __m256 s7 = _mm256_shuffle_ps(t5, t7, _MM_SHUFFLE(3, 2, 3, 2));
    r0 = _mm256_permute2f128_ps(s0, s4, 0x20);
    r1 = _mm256_permute2f128_ps(s1, s5, 0x20);
    r2 = _mm256_permute2f128_ps(s2, s6, 0x20);
    r3 = _mm256_permute2f128_ps(s3, s7, 0x20);
    r4 = _mm256_permute2f128_ps(s0, s4, 0x31);
    r5 = _mm256_permute2f128_ps(s1, s5, 0x31);
    r6 = _mm256_permute2f128_ps(s2, s6, 0x31);
    r7 = _mm256_permute2f128_ps(s3, s7, 0x31);
}

// dst[c * ls + j] = row j, value c (c, j < 8): an 8 x 8 block of rows (stride in elements)
// as 8 columns of lanes, ls apart
inline void block_to_lanes(const float* src, size_t stride, float* dst, size_t ls) {
    __m256 r[8];
    for (int j = 0; j < 8; j++) r[j] = _mm256_loadu_ps(src + j * stride);
    transpose8(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
    for (int j = 0; j < 8; j++) _mm256_storeu_ps(dst + ls * j, r[j]);
}

inline void block_to_lanes(const uint16_t* src, size_t stride, float* dst, size_t ls) {
    __m256 r[8];
    for (int j = 0; j < 8; j++) r[j] = load8_u16(src + j * stride).v;
    transpose8(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
    for (int j = 0; j < 8; j++) _mm256_storeu_ps(dst + ls * j, r[j]);
}

// the inverse for the output: src[c * ls + j] (8 columns of lanes) to rows j of 8 values
inline void lanes_to_block(const float* src, size_t ls, float* dst, size_t stride) {
    __m256 r[8];
    for (int j = 0; j < 8; j++) r[j] = _mm256_loadu_ps(src + ls * j);
    transpose8(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
    for (int j = 0; j < 8; j++) _mm256_storeu_ps(dst + j * stride, r[j]);
}

inline void lanes_to_block(const float* src, size_t ls, uint16_t* dst, size_t stride) {
    __m256 r[8];
    for (int j = 0; j < 8; j++) r[j] = _mm256_loadu_ps(src + ls * j);
    transpose8(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]);
    for (int j = 0; j < 8; j++) store8_u16(dst + j * stride, {r[j]});
}

#else

struct V8 {
    float v[8];
};

inline V8 load8(const float* p) {
    V8 r;
    for (int i = 0; i < 8; i++) r.v[i] = p[i];
    return r;
}
inline void store8(float* p, V8 a) {
    for (int i = 0; i < 8; i++) p[i] = a.v[i];
}
inline V8 set8(float w) {
    V8 r;
    for (int i = 0; i < 8; i++) r.v[i] = w;
    return r;
}
inline V8 zero8() { return set8(0.0f); }
inline V8 operator+(V8 a, V8 b) {
    V8 r;
    for (int i = 0; i < 8; i++) r.v[i] = a.v[i] + b.v[i];
    return r;
}
inline V8 operator*(V8 a, V8 b) {
    V8 r;
    for (int i = 0; i < 8; i++) r.v[i] = a.v[i] * b.v[i];
    return r;
}
inline V8 first_lanes(V8 a, V8 b, int n) {
    V8 r;
    for (int i = 0; i < 8; i++) r.v[i] = i < n ? a.v[i] : b.v[i];
    return r;
}
inline V8 load8_u16(const uint16_t* p) {
    V8 r;
    for (int i = 0; i < 8; i++) r.v[i] = (float)p[i];
    return r;
}
inline V8 round_u16(V8 a) {
    V8 r;
    for (int i = 0; i < 8; i++) r.v[i] = (float)sat_u16(a.v[i] + 0.5f);
    return r;
}
inline void store8_u16(uint16_t* p, V8 a) {
    for (int i = 0; i < 8; i++) p[i] = sat_u16(a.v[i]);
}

template <typename T>
inline void block_to_lanes(const T* src, size_t stride, float* dst, size_t ls) {
    for (int j = 0; j < 8; j++)
        for (int c = 0; c < 8; c++) dst[c * ls + j] = (float)src[j * stride + c];
}

inline void lanes_to_block(const float* src, size_t ls, float* dst, size_t stride) {
    for (int j = 0; j < 8; j++)
        for (int c = 0; c < 8; c++) dst[j * stride + c] = src[c * ls + j];
}

inline void lanes_to_block(const float* src, size_t ls, uint16_t* dst, size_t stride) {
    for (int j = 0; j < 8; j++)
        for (int c = 0; c < 8; c++) dst[j * stride + c] = sat_u16(src[c * ls + j]);
}

#endif

// the formulas
inline V8 seq4(V8 a, V8 b, V8 c, V8 d, V8 wa, V8 wb, V8 wc, V8 wd) {
    return ((wa * a + wb * b) + wc * c) + wd * d;
}
inline V8 pair4(V8 a, V8 b, V8 c, V8 d, V8 wa, V8 wb, V8 wc, V8 wd) {
    return (a * wa + b * wb) + (c * wc + d * wd);
}
inline V8 seq4_swap(V8 a, V8 b, V8 c, V8 d, V8 wa, V8 wb, V8 wc, V8 wd) {
    return ((wa * a + wb * b) + wd * d) + wc * c;
}
inline V8 sum2(V8 a, V8 b, V8 wa, V8 wb) { return wa * a + wb * b; }

// + 0.5f: with the truncation of the uint16 stores, the mex's rounding to uint16
inline V8 add_half(V8 a) { return a + set8(0.5f); }

// single values as lanes
inline V8 load8_any(const float* p) { return load8(p); }
inline V8 load8_any(const uint16_t* p) { return load8_u16(p); }

// values of one row in lanes, n <= 8 of them (the others 0)
template <typename T>
inline V8 load_partial(const T* p, int n) {
    float tmp[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < n; i++) tmp[i] = (float)p[i];
    return load8(tmp);
}

template <typename T>
inline void store_partial(T* p, V8 a, int n);

template <>
inline void store_partial<float>(float* p, V8 a, int n) {
    float tmp[8];
    store8(tmp, a);
    for (int i = 0; i < n; i++) p[i] = tmp[i];
}

template <>
inline void store_partial<uint16_t>(uint16_t* p, V8 a, int n) {
    float tmp[8];
    store8(tmp, a);
    for (int i = 0; i < n; i++) p[i] = sat_u16(tmp[i]);
}

inline void store8_any(float* p, V8 a) { store8(p, a); }
inline void store8_any(uint16_t* p, V8 a) { store8_u16(p, a); }

// copy n values to dst, bypassing the caches where the processor can (the output is written
// once); stream_fence() makes such stores visible before the kernel returns
template <typename TO>
inline void stream_rows(const TO* src, TO* dst, uint64_t n) {
#if defined(__AVX2__)
    uint64_t i = 0;
    while (i < n && (reinterpret_cast<uintptr_t>(dst + i) & 31) != 0) {
        dst[i] = src[i];
        ++i;
    }
    constexpr uint64_t step = 32 / sizeof(TO);
    for (; i + step <= n; i += step)
        _mm256_stream_si256((__m256i*)(dst + i), _mm256_loadu_si256((const __m256i*)(src + i)));
    for (; i < n; ++i) dst[i] = src[i];
#else
    memcpy(dst, src, n * sizeof(TO));
#endif
}

inline void stream_fence() {
#if defined(__AVX2__)
    _mm_sfence();
#endif
}

}  // namespace pk5d::deskew
