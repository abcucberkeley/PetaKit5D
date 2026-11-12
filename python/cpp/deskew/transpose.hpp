// The last two axes of a C-ordered volume swapped: (n, r, c) to (n, c, r), in threads. Our
// (z, y, x) volumes become MATLAB's layout ((z, x, y) in C order, y fastest) for the kernels in
// ml/, and back. 8 x 8 blocks (AVX2 shuffles, or plain loops), split among threads by blocks
// of 64 rows of each plane.
#pragma once

#include <algorithm>
#include <cstdint>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

namespace pk5d::deskew {

#if defined(__AVX2__)

// an 8 x 8 block of floats: src rows of stride ss to dst rows of stride ds
inline void transpose_block8(const float* src, uint64_t ss, float* dst, uint64_t ds) {
    __m256 r0 = _mm256_loadu_ps(src), r1 = _mm256_loadu_ps(src + ss);
    __m256 r2 = _mm256_loadu_ps(src + 2 * ss), r3 = _mm256_loadu_ps(src + 3 * ss);
    __m256 r4 = _mm256_loadu_ps(src + 4 * ss), r5 = _mm256_loadu_ps(src + 5 * ss);
    __m256 r6 = _mm256_loadu_ps(src + 6 * ss), r7 = _mm256_loadu_ps(src + 7 * ss);
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
    _mm256_storeu_ps(dst, _mm256_permute2f128_ps(s0, s4, 0x20));
    _mm256_storeu_ps(dst + ds, _mm256_permute2f128_ps(s1, s5, 0x20));
    _mm256_storeu_ps(dst + 2 * ds, _mm256_permute2f128_ps(s2, s6, 0x20));
    _mm256_storeu_ps(dst + 3 * ds, _mm256_permute2f128_ps(s3, s7, 0x20));
    _mm256_storeu_ps(dst + 4 * ds, _mm256_permute2f128_ps(s0, s4, 0x31));
    _mm256_storeu_ps(dst + 5 * ds, _mm256_permute2f128_ps(s1, s5, 0x31));
    _mm256_storeu_ps(dst + 6 * ds, _mm256_permute2f128_ps(s2, s6, 0x31));
    _mm256_storeu_ps(dst + 7 * ds, _mm256_permute2f128_ps(s3, s7, 0x31));
}

// an 8 x 8 block of uint16 values
inline void transpose_block8(const uint16_t* src, uint64_t ss, uint16_t* dst, uint64_t ds) {
    __m128i r[8];
    for (int i = 0; i < 8; i++) r[i] = _mm_loadu_si128((const __m128i*)(src + i * ss));
    __m128i a0 = _mm_unpacklo_epi16(r[0], r[1]), a1 = _mm_unpackhi_epi16(r[0], r[1]);
    __m128i a2 = _mm_unpacklo_epi16(r[2], r[3]), a3 = _mm_unpackhi_epi16(r[2], r[3]);
    __m128i a4 = _mm_unpacklo_epi16(r[4], r[5]), a5 = _mm_unpackhi_epi16(r[4], r[5]);
    __m128i a6 = _mm_unpacklo_epi16(r[6], r[7]), a7 = _mm_unpackhi_epi16(r[6], r[7]);
    __m128i b0 = _mm_unpacklo_epi32(a0, a2), b1 = _mm_unpackhi_epi32(a0, a2);
    __m128i b2 = _mm_unpacklo_epi32(a1, a3), b3 = _mm_unpackhi_epi32(a1, a3);
    __m128i b4 = _mm_unpacklo_epi32(a4, a6), b5 = _mm_unpackhi_epi32(a4, a6);
    __m128i b6 = _mm_unpacklo_epi32(a5, a7), b7 = _mm_unpackhi_epi32(a5, a7);
    _mm_storeu_si128((__m128i*)(dst), _mm_unpacklo_epi64(b0, b4));
    _mm_storeu_si128((__m128i*)(dst + ds), _mm_unpackhi_epi64(b0, b4));
    _mm_storeu_si128((__m128i*)(dst + 2 * ds), _mm_unpacklo_epi64(b1, b5));
    _mm_storeu_si128((__m128i*)(dst + 3 * ds), _mm_unpackhi_epi64(b1, b5));
    _mm_storeu_si128((__m128i*)(dst + 4 * ds), _mm_unpacklo_epi64(b2, b6));
    _mm_storeu_si128((__m128i*)(dst + 5 * ds), _mm_unpackhi_epi64(b2, b6));
    _mm_storeu_si128((__m128i*)(dst + 6 * ds), _mm_unpacklo_epi64(b3, b7));
    _mm_storeu_si128((__m128i*)(dst + 7 * ds), _mm_unpackhi_epi64(b3, b7));
}

#else

template <typename T>
inline void transpose_block8(const T* src, uint64_t ss, T* dst, uint64_t ds) {
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++) dst[j * ds + i] = src[i * ss + j];
}

#endif

// dst (n, c, r) from src (n, r, c), both C-contiguous
template <typename T>
void swap_last_axes(const T* src, T* dst, const uint64_t n, const uint64_t r, const uint64_t c,
                    const int nthreads) {
    constexpr uint64_t RB = 64;  // rows of a task
    constexpr uint64_t CB = 64;  // columns of a tile: dst rows get 64 contiguous values
    const uint64_t nrb = (r + RB - 1) / RB;
    const uint64_t ntasks = n * nrb;
#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
    for (uint64_t task = 0; task < ntasks; ++task) {
        const uint64_t p = task / nrb;
        const uint64_t r0 = (task % nrb) * RB, r1 = std::min(r, r0 + RB);
        const T* s = src + p * r * c;
        T* d = dst + p * r * c;
        const uint64_t r8 = r0 + (r1 - r0) / 8 * 8;
        for (uint64_t j0 = 0; j0 < c; j0 += CB) {
            const uint64_t j1 = std::min(c, j0 + CB);
            const uint64_t j8 = j0 + (j1 - j0) / 8 * 8;
            for (uint64_t j = j0; j < j8; j += 8)
                for (uint64_t i = r0; i < r8; i += 8)
                    transpose_block8(s + i * c + j, c, d + j * r + i, r);
            // the rows and columns past the last whole block
            for (uint64_t j = j8; j < j1; ++j)
                for (uint64_t i = r0; i < r1; ++i) d[j * r + i] = s[i * c + j];
            for (uint64_t j = j0; j < j8; ++j)
                for (uint64_t i = r8; i < r1; ++i) d[j * r + i] = s[i * c + j];
        }
    }
}

}  // namespace pk5d::deskew
