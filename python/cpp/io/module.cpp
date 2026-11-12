// petakit5d._ext.io: the layout change after a tiff read (task IO2).
//
// cpp-tiff's Python package returns a stack as (z, y, x) C order: page by page, each page
// row by row (x fastest). PetaKit5D's volumes are (y, x, z) Fortran order (y fastest), the
// layout the mex fills while it decodes. transpose_pages_inplace turns one into the other in
// the reader's own buffer (what petakit5d.io.read_tiff uses), transpose_pages copies it into
// another array; both page by page, in tiles that fit the first-level cache, split among the
// call's threads. They move bytes only, so the result never depends on the path (SIMD or
// scalar) or the thread count.
#include <pk5d/common.hpp>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

#if defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
#include <emmintrin.h>
#define PK5D_IO_SSE2 1
#endif

namespace nb = nanobind;

namespace {

// tile edge in elements: 64 x 64 x 2 bytes in and out stay in the first-level cache
constexpr int64_t kTile = 64;

// one tile of a page: src is (nx, ny) with x fastest, dst (ny, nx) with y fastest
template <typename T>
inline void tile_scalar(const T *s, T *d, int64_t nx, int64_t ny, int64_t x0, int64_t x1,
                        int64_t y0, int64_t y1) {
    for (int64_t x = x0; x < x1; ++x) {
        T *dc = d + x * ny;
        const T *sc = s + x;
        for (int64_t y = y0; y < y1; ++y) dc[y] = sc[y * nx];
    }
}

#ifdef PK5D_IO_SSE2
// 8 x 8 blocks of 2-byte elements: 8 rows in, 8 columns out
inline void block8_u16(const uint16_t *s, uint16_t *d, int64_t nx, int64_t ny) {
    const __m128i a0 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s));
    const __m128i a1 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + nx));
    const __m128i a2 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + 2 * nx));
    const __m128i a3 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + 3 * nx));
    const __m128i a4 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + 4 * nx));
    const __m128i a5 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + 5 * nx));
    const __m128i a6 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + 6 * nx));
    const __m128i a7 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + 7 * nx));
    const __m128i t0 = _mm_unpacklo_epi16(a0, a1), t1 = _mm_unpackhi_epi16(a0, a1);
    const __m128i t2 = _mm_unpacklo_epi16(a2, a3), t3 = _mm_unpackhi_epi16(a2, a3);
    const __m128i t4 = _mm_unpacklo_epi16(a4, a5), t5 = _mm_unpackhi_epi16(a4, a5);
    const __m128i t6 = _mm_unpacklo_epi16(a6, a7), t7 = _mm_unpackhi_epi16(a6, a7);
    const __m128i u0 = _mm_unpacklo_epi32(t0, t2), u1 = _mm_unpackhi_epi32(t0, t2);
    const __m128i u2 = _mm_unpacklo_epi32(t1, t3), u3 = _mm_unpackhi_epi32(t1, t3);
    const __m128i u4 = _mm_unpacklo_epi32(t4, t6), u5 = _mm_unpackhi_epi32(t4, t6);
    const __m128i u6 = _mm_unpacklo_epi32(t5, t7), u7 = _mm_unpackhi_epi32(t5, t7);
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d), _mm_unpacklo_epi64(u0, u4));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d + ny), _mm_unpackhi_epi64(u0, u4));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d + 2 * ny), _mm_unpacklo_epi64(u1, u5));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d + 3 * ny), _mm_unpackhi_epi64(u1, u5));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d + 4 * ny), _mm_unpacklo_epi64(u2, u6));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d + 5 * ny), _mm_unpackhi_epi64(u2, u6));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d + 6 * ny), _mm_unpacklo_epi64(u3, u7));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d + 7 * ny), _mm_unpackhi_epi64(u3, u7));
}

// 4 x 4 blocks of 4-byte elements
inline void block4_u32(const uint32_t *s, uint32_t *d, int64_t nx, int64_t ny) {
    const __m128i a0 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s));
    const __m128i a1 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + nx));
    const __m128i a2 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + 2 * nx));
    const __m128i a3 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + 3 * nx));
    const __m128i t0 = _mm_unpacklo_epi32(a0, a1), t1 = _mm_unpackhi_epi32(a0, a1);
    const __m128i t2 = _mm_unpacklo_epi32(a2, a3), t3 = _mm_unpackhi_epi32(a2, a3);
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d), _mm_unpacklo_epi64(t0, t2));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d + ny), _mm_unpackhi_epi64(t0, t2));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d + 2 * ny), _mm_unpacklo_epi64(t1, t3));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d + 3 * ny), _mm_unpackhi_epi64(t1, t3));
}

// 2 x 2 blocks of 8-byte elements
inline void block2_u64(const uint64_t *s, uint64_t *d, int64_t nx, int64_t ny) {
    const __m128i a0 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s));
    const __m128i a1 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(s + nx));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d), _mm_unpacklo_epi64(a0, a1));
    _mm_storeu_si128(reinterpret_cast<__m128i *>(d + ny), _mm_unpackhi_epi64(a0, a1));
}

template <typename T>
struct simd_block {
    static constexpr int64_t n = 0;
    static void run(const T *, T *, int64_t, int64_t) {}
};
template <>
struct simd_block<uint16_t> {
    static constexpr int64_t n = 8;
    static void run(const uint16_t *s, uint16_t *d, int64_t nx, int64_t ny) {
        block8_u16(s, d, nx, ny);
    }
};
template <>
struct simd_block<uint32_t> {
    static constexpr int64_t n = 4;
    static void run(const uint32_t *s, uint32_t *d, int64_t nx, int64_t ny) {
        block4_u32(s, d, nx, ny);
    }
};
template <>
struct simd_block<uint64_t> {
    static constexpr int64_t n = 2;
    static void run(const uint64_t *s, uint64_t *d, int64_t nx, int64_t ny) {
        block2_u64(s, d, nx, ny);
    }
};
#endif

template <typename T>
inline void tile(const T *s, T *d, int64_t nx, int64_t ny, int64_t x0, int64_t x1, int64_t y0,
                 int64_t y1) {
#ifdef PK5D_IO_SSE2
    constexpr int64_t b = simd_block<T>::n;
    if constexpr (b > 0) {
        const int64_t xe = x0 + (x1 - x0) / b * b, ye = y0 + (y1 - y0) / b * b;
        for (int64_t x = x0; x < xe; x += b)
            for (int64_t y = y0; y < ye; y += b)
                simd_block<T>::run(s + y * nx + x, d + x * ny + y, nx, ny);
        // the edges: rows below the blocks, then columns right of them
        if (ye < y1) tile_scalar(s, d, nx, ny, x0, xe, ye, y1);
        if (xe < x1) tile_scalar(s, d, nx, ny, xe, x1, y0, y1);
        return;
    }
#endif
    tile_scalar(s, d, nx, ny, x0, x1, y0, y1);
}

template <typename T>
void transpose_pages_t(const T *src, T *dst, int64_t nx, int64_t ny, int64_t nz, int threads) {
    const int64_t nby = (ny + kTile - 1) / kTile, nbx = (nx + kTile - 1) / kTile;
    const int64_t per_page = nby * nbx, ntiles = nz * per_page, page = nx * ny;
    // tiles in the order of the output (y fastest), so that a thread's static share is a
    // contiguous part of it
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int64_t t = 0; t < ntiles; ++t) {
        const int64_t z = t / per_page, r = t % per_page;
        const int64_t bx = r / nby, by = r % nby;
        const int64_t x0 = bx * kTile, y0 = by * kTile;
        const int64_t x1 = x0 + kTile < nx ? x0 + kTile : nx;
        const int64_t y1 = y0 + kTile < ny ? y0 + kTile : ny;
        tile(src + z * page, dst + z * page, nx, ny, x0, x1, y0, y1);
    }
}

// one page through a scratch copy, its tiles in the order of the output
template <typename T>
inline void page_tiles(const T *s, T *d, int64_t nx, int64_t ny, int64_t t0, int64_t t1) {
    const int64_t nby = (ny + kTile - 1) / kTile;
    for (int64_t t = t0; t < t1; ++t) {
        const int64_t bx = t / nby, by = t % nby;
        const int64_t x0 = bx * kTile, y0 = by * kTile;
        tile(s, d, nx, ny, x0, x0 + kTile < nx ? x0 + kTile : nx, y0,
             y0 + kTile < ny ? y0 + kTile : ny);
    }
}

// In place: each page is copied to a scratch buffer that stays in the cache and transposed
// back into its own bytes. No second volume, and no first touch of new memory. With enough
// pages each thread takes whole pages (a scratch page each, if their total stays small);
// otherwise the threads share each page.
template <typename T>
void transpose_pages_inplace_t(T *buf, int64_t nx, int64_t ny, int64_t nz, int threads, int mode) {
    const int64_t page = nx * ny;
    const uint64_t page_bytes = static_cast<uint64_t>(page) * sizeof(T);
    const uint64_t total = page_bytes * static_cast<uint64_t>(nz);
    const uint64_t budget = std::max<uint64_t>(uint64_t(64) << 20, total / 16);
    const bool whole_pages = mode == 1 || (mode == 0 && nz >= threads &&
                                           page_bytes * static_cast<uint64_t>(threads) <= budget);
    if (whole_pages) {
        std::vector<T> scratch(static_cast<size_t>(page) * threads);
#pragma omp parallel num_threads(threads)
        {
#ifdef _OPENMP
            T *sc = scratch.data() + static_cast<size_t>(page) * omp_get_thread_num();
#else
            T *sc = scratch.data();
#endif
            const int64_t ntiles = ((ny + kTile - 1) / kTile) * ((nx + kTile - 1) / kTile);
#pragma omp for schedule(static)
            for (int64_t z = 0; z < nz; ++z) {
                std::memcpy(sc, buf + z * page, page_bytes);
                page_tiles(sc, buf + z * page, nx, ny, 0, ntiles);
            }
        }
        return;
    }
    std::vector<T> scratch(static_cast<size_t>(page));
    T *sc = scratch.data();
    const int64_t ntiles = ((ny + kTile - 1) / kTile) * ((nx + kTile - 1) / kTile);
    // the copy in pieces of 64 KiB
    const int64_t piece = std::max<int64_t>(1, (int64_t(1) << 16) / static_cast<int64_t>(sizeof(T)));
    const int64_t npieces = (page + piece - 1) / piece;
#pragma omp parallel num_threads(threads)
    for (int64_t z = 0; z < nz; ++z) {
        T *p = buf + z * page;
#pragma omp for schedule(static)
        for (int64_t i = 0; i < npieces; ++i) {
            const int64_t a = i * piece, n = std::min(piece, page - a);
            std::memcpy(sc + a, p + a, static_cast<size_t>(n) * sizeof(T));
        }
#pragma omp for schedule(static)
        for (int64_t t = 0; t < ntiles; ++t) page_tiles(sc, p, nx, ny, t, t + 1);
    }
}

using flat_array = nb::ndarray<nb::ndim<1>, nb::c_contig, nb::device::cpu>;

void transpose_pages_inplace(flat_array buf, int64_t nx, int64_t ny, int64_t nz, int num_threads,
                             int mode) {
    pk5d::check_threads(num_threads);
    const size_t bits = buf.dtype().bits * buf.dtype().lanes;
    if (bits != 8 && bits != 16 && bits != 32 && bits != 64)
        throw std::invalid_argument("transpose_pages_inplace: elements of 1, 2, 4 or 8 bytes only");
    if (nx < 0 || ny < 0 || nz < 0 || static_cast<int64_t>(buf.shape(0)) != nx * ny * nz)
        throw std::invalid_argument("transpose_pages_inplace: buf must hold nx * ny * nz elements");
    if (mode < 0 || mode > 2)
        throw std::invalid_argument("transpose_pages_inplace: mode is 0 (auto), 1 or 2");
    if (nx == 0 || ny == 0 || nz == 0) return;
    void *d = buf.data();
    nb::gil_scoped_release release;
    switch (bits) {
        case 8:
            transpose_pages_inplace_t(static_cast<uint8_t *>(d), nx, ny, nz, num_threads, mode);
            break;
        case 16:
            transpose_pages_inplace_t(static_cast<uint16_t *>(d), nx, ny, nz, num_threads, mode);
            break;
        case 32:
            transpose_pages_inplace_t(static_cast<uint32_t *>(d), nx, ny, nz, num_threads, mode);
            break;
        default:
            transpose_pages_inplace_t(static_cast<uint64_t *>(d), nx, ny, nz, num_threads, mode);
            break;
    }
}

using src_array = nb::ndarray<nb::ro, nb::ndim<3>, nb::f_contig, nb::device::cpu>;
using dst_array = nb::ndarray<nb::ndim<3>, nb::f_contig, nb::device::cpu>;

void transpose_pages(src_array src, dst_array dst, int num_threads) {
    pk5d::check_threads(num_threads);
    if (src.dtype() != dst.dtype())
        throw std::invalid_argument("transpose_pages: src and dst have different data types");
    const size_t bits = src.dtype().bits * src.dtype().lanes;
    const int64_t nx = static_cast<int64_t>(src.shape(0)), ny = static_cast<int64_t>(src.shape(1)),
                  nz = static_cast<int64_t>(src.shape(2));
    if (static_cast<int64_t>(dst.shape(0)) != ny || static_cast<int64_t>(dst.shape(1)) != nx ||
        static_cast<int64_t>(dst.shape(2)) != nz)
        throw std::invalid_argument("transpose_pages: dst must have the shape (ny, nx, nz) of src (nx, ny, nz)");
    if (bits != 8 && bits != 16 && bits != 32 && bits != 64)
        throw std::invalid_argument("transpose_pages: elements of 1, 2, 4 or 8 bytes only");
    if (nx == 0 || ny == 0 || nz == 0) return;
    const void *s = src.data();
    void *d = dst.data();
    nb::gil_scoped_release release;
    switch (bits) {
        case 8:
            transpose_pages_t(static_cast<const uint8_t *>(s), static_cast<uint8_t *>(d), nx, ny, nz,
                              num_threads);
            break;
        case 16:
            transpose_pages_t(static_cast<const uint16_t *>(s), static_cast<uint16_t *>(d), nx, ny, nz,
                              num_threads);
            break;
        case 32:
            transpose_pages_t(static_cast<const uint32_t *>(s), static_cast<uint32_t *>(d), nx, ny, nz,
                              num_threads);
            break;
        case 64:
            transpose_pages_t(static_cast<const uint64_t *>(s), static_cast<uint64_t *>(d), nx, ny, nz,
                              num_threads);
            break;
        default:
            break;
    }
}

}  // namespace

NB_MODULE(io, m) {
    m.doc() = "Layout changes of the tiff and zarr readers and writers (petakit5d.io)";
    m.def("transpose_pages", &transpose_pages, nb::arg("src").noconvert(),
          nb::arg("dst").noconvert(), nb::arg("num_threads"),
          "Copy src, (nx, ny, nz) Fortran order (pages of rows, x fastest), into dst,\n"
          "(ny, nx, nz) Fortran order (y fastest), with num_threads threads.");
    m.def("transpose_pages_inplace", &transpose_pages_inplace, nb::arg("buf").noconvert(),
          nb::arg("nx"), nb::arg("ny"), nb::arg("nz"), nb::arg("num_threads"), nb::arg("mode") = 0,
          "Turn buf, nz pages of ny rows of nx elements (x fastest), into the (ny, nx, nz)\n"
          "Fortran-ordered volume (y fastest) in place, with num_threads threads. mode: 0 chooses,\n"
          "1 gives each thread whole pages, 2 has the threads share each page.");
    m.def("simd", []() {
#ifdef PK5D_IO_SSE2
        return true;
#else
        return false;
#endif
    }, "Whether the SSE2 blocks are compiled in.");
}
