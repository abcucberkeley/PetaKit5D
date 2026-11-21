// Order-statistic filters with box windows (medfilt2, medfilt3) on C-contiguous (z, y, x)
// volumes: out(z, y, x) is the k-th smallest (0-based) of the window
// P(z .. z + wz - 1, y .. y + wy - 1, x .. x + wx - 1) of the padded volume P (the caller pads
// as MATLAB does). Equal numbers have equal bits except +0 and -0, so windows without NaN whose
// k-th smallest is not a zero, or whose zeros all have one sign, have one answer whatever the
// algorithm: that answer is MATLAB's. MATLAB's own choice between +0 and -0 depends on the
// position of the values in the window (and is +0 even for windows of -0 only), and its NaN rule
// on the code path (ordfilt2, IPP, mwmedfilt3; measured): here a window with NaN gives its first
// NaN in MATLAB's order (y fastest, then x, then z) and a zero result is +0.
//
// 3 x 3 x 3 single windows without NaN and -0 use the selection network of
// omw_backprojector_mex (Batcher's odd-even merge sort of 32 values, 5 of them +Inf, pruned to
// the 14th smallest), several voxels along x at a time: on x86-64 with GCC or Clang 4 (SSE) or
// 8 (AVX2, chosen at run time), elsewhere 16.
#include <pk5d/common.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

// x86-64 with GCC or Clang: the network's compare-exchanges as the minimum and maximum
// instructions, the instruction set chosen at run time. (GCC makes scalar code of the generic
// form's 16-lane selects where the instruction set has no 64-byte vectors: its vector lowering
// splits their comparisons into 16 scalar ones, where it splits arithmetic into the registers'
// widths. Its 4-lane selects are a compare and three bitwise operations.)
#if (defined(__x86_64__) || defined(_M_X64)) && (defined(__GNUC__) || defined(__clang__))
#define PK5D_MEDIAN_DISPATCH 1
#include <immintrin.h>
#endif

namespace {

namespace nb = nanobind;
using idx_t = std::ptrdiff_t;

template <typename T>
void order_filter(const T* P, T* out, const idx_t py, const idx_t px, const idx_t nz, const idx_t ny,
                  const idx_t nx, const idx_t wz, const idx_t wy, const idx_t wx, const idx_t k, const int nt) {
    const idx_t wn = wz * wy * wx;
    const idx_t nrows = nz * ny;
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
    {
        std::vector<T> buf((size_t)wn);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (idx_t r = 0; r < nrows; ++r) {
            const idx_t z = r / ny, y = r % ny;
            T* o = out + r * nx;
            for (idx_t x = 0; x < nx; ++x) {
                idx_t m = 0;
                bool nan = false;
                // MATLAB's order within the window: y fastest, then x, then z
                for (idx_t dz = 0; dz < wz; ++dz) {
                    for (idx_t dx = 0; dx < wx; ++dx) {
                        const T* s = P + ((z + dz) * py + y) * px + x + dx;
                        for (idx_t dy = 0; dy < wy; ++dy) {
                            const T v = s[dy * px];
                            if constexpr (std::is_floating_point_v<T>) {
                                if (!nan && std::isnan(v)) {
                                    nan = true;
                                    o[x] = v;
                                }
                            }
                            buf[(size_t)m++] = v;
                        }
                    }
                }
                if (nan) continue;
                std::nth_element(buf.begin(), buf.begin() + k, buf.end());
                T v = buf[(size_t)k];
                if constexpr (std::is_floating_point_v<T>) {
                    if (v == 0) v = 0;
                }
                o[x] = v;
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// 3 x 3 x 3 single median (omw_backprojector_mex)

static const int ML = 16;  // the most voxels along x per network pass

// the 14th smallest of v[0 .. 26] in v[13]: CX(a, b) puts the smaller of v[a] and v[b] in v[a],
// the larger in v[b]
#define PK5D_MEDIAN27_NETWORK(CX)                                                              \
    CX(0, 1); CX(2, 3); CX(0, 2); CX(1, 3); CX(1, 2); CX(4, 5);                                \
    CX(6, 7); CX(4, 6); CX(5, 7); CX(5, 6); CX(0, 4); CX(2, 6);                                \
    CX(2, 4); CX(1, 5); CX(3, 7); CX(3, 5); CX(1, 2); CX(3, 4);                                \
    CX(5, 6); CX(8, 9); CX(10, 11); CX(8, 10); CX(9, 11); CX(9, 10);                           \
    CX(12, 13); CX(14, 15); CX(12, 14); CX(13, 15); CX(13, 14); CX(8, 12);                     \
    CX(10, 14); CX(10, 12); CX(9, 13); CX(11, 15); CX(11, 13); CX(9, 10);                      \
    CX(11, 12); CX(13, 14); CX(0, 8); CX(4, 12); CX(4, 8); CX(2, 10);                          \
    CX(6, 14); CX(6, 10); CX(2, 4); CX(6, 8); CX(10, 12); CX(1, 9);                            \
    CX(5, 13); CX(5, 9); CX(3, 11); CX(7, 15); CX(7, 11); CX(3, 5);                            \
    CX(7, 9); CX(11, 13); CX(1, 2); CX(3, 4); CX(5, 6); CX(7, 8);                              \
    CX(9, 10); CX(11, 12); CX(13, 14); CX(16, 17); CX(18, 19); CX(16, 18);                     \
    CX(17, 19); CX(17, 18); CX(20, 21); CX(22, 23); CX(20, 22); CX(21, 23);                    \
    CX(21, 22); CX(16, 20); CX(18, 22); CX(18, 20); CX(17, 21); CX(19, 23);                    \
    CX(19, 21); CX(17, 18); CX(19, 20); CX(21, 22); CX(24, 25); CX(24, 26);                    \
    CX(25, 26); CX(25, 26); CX(16, 24); CX(20, 24); CX(18, 26); CX(22, 26);                    \
    CX(18, 20); CX(22, 24); CX(17, 25); CX(21, 25); CX(19, 21); CX(23, 25);                    \
    CX(17, 18); CX(19, 20); CX(21, 22); CX(23, 24); CX(25, 26); CX(0, 16);                     \
    CX(8, 24); CX(8, 16); CX(4, 20); CX(12, 20); CX(12, 16); CX(2, 18);                        \
    CX(10, 26); CX(10, 18); CX(6, 22); CX(14, 22); CX(14, 18); CX(14, 16);                     \
    CX(1, 17); CX(9, 25); CX(9, 17); CX(5, 21); CX(13, 21); CX(13, 17);                        \
    CX(3, 19); CX(11, 19); CX(7, 23); CX(7, 11); CX(11, 13); CX(13, 14);

#ifdef PK5D_MEDIAN_DISPATCH
// the line pass in SSE (4 voxels a pass) and AVX2 (8)
#define PK5D_MEDIAN_NS isa_sse
#define PK5D_MEDIAN_ISA
#define PK5D_MEDIAN_V __m128
#define PK5D_MEDIAN_OP(f) _mm_##f
#include "median27.inc"
#undef PK5D_MEDIAN_NS
#undef PK5D_MEDIAN_ISA
#undef PK5D_MEDIAN_V
#undef PK5D_MEDIAN_OP
#define PK5D_MEDIAN_NS isa_avx2
#define PK5D_MEDIAN_ISA __attribute__((target("avx2")))
#define PK5D_MEDIAN_V __m256
#define PK5D_MEDIAN_OP(f) _mm256_##f
#include "median27.inc"
#undef PK5D_MEDIAN_NS
#undef PK5D_MEDIAN_ISA
#undef PK5D_MEDIAN_V
#undef PK5D_MEDIAN_OP
#else
#if defined(__GNUC__) || defined(__clang__)
typedef float vfloat __attribute__((vector_size(4 * ML)));
#define PK5D_MEDIAN_CX(a, b) { const vfloat x = v[a], y = v[b]; v[a] = y < x ? y : x; v[b] = y < x ? x : y; }
#else
// compilers without vector extensions (MSVC): the same compare-exchanges lane by lane, the same
// selections, so the same bits
struct vfloat {
    float l[ML];
};
#define PK5D_MEDIAN_CX(a, b)                                                                       \
    {                                                                                              \
        float* p_ = v[a].l;                                                                        \
        float* q_ = v[b].l;                                                                        \
        for (int l_ = 0; l_ < ML; l_++) {                                                          \
            const float x = p_[l_], y = q_[l_];                                                    \
            p_[l_] = y < x ? y : x;                                                                \
            q_[l_] = y < x ? x : y;                                                                \
        }                                                                                          \
    }
#endif

// the 14th smallest of v[0 .. 26] in v[13], for ML voxels at once
inline void median27(vfloat* v) {
    PK5D_MEDIAN27_NETWORK(PK5D_MEDIAN_CX)
}

namespace isa_generic {

// the medians of a line of n0 voxels along x whose windows' lines are rows[0 .. 8] (q = 3 dz + dx),
// ML voxels a pass
void median27_line(const float* const* rows, float* o, const idx_t n0) {
    vfloat v[27];
    for (idx_t i0 = 0; i0 < n0; i0 += ML) {
        for (int q = 0; q < 9; q++) {
            for (int dy = 0; dy < 3; dy++) std::memcpy(&v[q * 3 + dy], rows[q] + i0 + dy, sizeof(vfloat));
        }
        median27(v);
        const idx_t nl = std::min<idx_t>(ML, n0 - i0);
        float med[ML];
        std::memcpy(med, &v[13], sizeof(vfloat));
        std::memcpy(o + i0, med, (size_t)nl * sizeof(float));
    }
}

}  // namespace isa_generic
#endif

// the network's instruction sets: 0 baseline (SSE on x86-64), 1 AVX2
int median27_best_isa() {
#ifdef PK5D_MEDIAN_DISPATCH
    static const int best = [] {
        __builtin_cpu_init();
        return __builtin_cpu_supports("avx2") ? 1 : 0;
    }();
    return best;
#else
    return 0;
#endif
}

// the set a call runs: -1 the best one, else the one asked for (refused if the CPU lacks it)
int pick_isa(int isa) {
    if (isa < 0) return median27_best_isa();
    if (isa > median27_best_isa()) throw std::invalid_argument("this CPU lacks that instruction set");
    return isa;
}

// P padded by 1 on every side (p0 = n0 + 2, p1 = n1 + 2), with ML floats of room after its end
void median27_filter(const float* P, float* out, const idx_t n0, const idx_t n1, const idx_t n2, const int nt,
                     const int isa) {
    const idx_t p0 = n0 + 2, p1 = n1 + 2;
    const idx_t ncols = n1 * n2;
#ifdef PK5D_MEDIAN_DISPATCH
    void (*const line)(const float* const*, float*, idx_t) = isa == 1 ? isa_avx2::median27_line
                                                                       : isa_sse::median27_line;
#else
    (void)isa;
    void (*const line)(const float* const*, float*, idx_t) = isa_generic::median27_line;
#endif
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(dynamic, 4)
#endif
    for (idx_t c = 0; c < ncols; ++c) {
        const idx_t x = c % n1, z = c / n1;
        const float* rows[9];
        int m = 0;
        for (int dz = 0; dz < 3; dz++) {
            for (int dx = 0; dx < 3; dx++) rows[m++] = P + ((z + dz) * p1 + (x + dx)) * p0;
        }
        line(rows, out + c * n0, n0);
    }
}

template <typename T>
using cvol = pk5d::cvolume<T>;

template <typename T>
void order_filter_bind(cvol<const T> P, cvol<T> out, int wz, int wy, int wx, int k, int num_threads) {
    pk5d::check_threads(num_threads);
    const idx_t nz = (idx_t)out.shape(0), ny = (idx_t)out.shape(1), nx = (idx_t)out.shape(2);
    if (wz < 1 || wy < 1 || wx < 1) throw std::invalid_argument("order_filter: window sizes must be positive");
    if ((idx_t)P.shape(0) != nz + wz - 1 || (idx_t)P.shape(1) != ny + wy - 1 || (idx_t)P.shape(2) != nx + wx - 1)
        throw std::invalid_argument("order_filter: the padded volume must be the output size + window - 1");
    if (k < 0 || (idx_t)k >= (idx_t)wz * wy * wx) throw std::invalid_argument("order_filter: order out of range");
    if (nz * ny * nx == 0) return;
    nb::gil_scoped_release release;
    order_filter<T>(P.data(), out.data(), (idx_t)P.shape(1), (idx_t)P.shape(2), nz, ny, nx, wz, wy, wx, k,
                    num_threads);
}

// the 3 x 3 x 3 median of a single volume padded by 1, true when the network applies (no NaN,
// no -0 in P); P has ML floats of room after its end (the caller allocates them). Memory axes 0,
// 1, 2 of median27_filter are x, y, z.
bool median27_bind(cvol<const float> P, cvol<float> out, int num_threads, int isa) {
    pk5d::check_threads(num_threads);
    isa = pick_isa(isa);
    const idx_t nz = (idx_t)out.shape(0), ny = (idx_t)out.shape(1), nx = (idx_t)out.shape(2);
    if ((idx_t)P.shape(0) != nz + 2 || (idx_t)P.shape(1) != ny + 2 || (idx_t)P.shape(2) != nx + 2)
        throw std::invalid_argument("median27: the padded volume must be the output size + 2");
    const float* p = P.data();
    const idx_t np = (idx_t)(P.shape(0) * P.shape(1) * P.shape(2));
    nb::gil_scoped_release release;
    int bad = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) schedule(static) reduction(| : bad)
#endif
    for (idx_t i = 0; i < np; i++) bad |= (p[i] != p[i]) || (p[i] == 0 && std::signbit(p[i]));
    if (bad) return false;
    if (nz * ny * nx == 0) return true;
    median27_filter(p, out.data(), nx, ny, nz, num_threads, isa);
    return true;
}

// medfilt3(A) of a single volume (3 x 3 x 3, symmetric padding: the border voxel repeated):
// the padded copy made here in parallel (with room for the network's last pass); false (out
// untouched) when A has NaN or -0
bool median27_symmetric_bind(cvol<const float> A, cvol<float> out, int num_threads, int isa) {
    pk5d::check_threads(num_threads);
    isa = pick_isa(isa);
    const idx_t nz = (idx_t)A.shape(0), ny = (idx_t)A.shape(1), nx = (idx_t)A.shape(2);
    if ((idx_t)out.shape(0) != nz || (idx_t)out.shape(1) != ny || (idx_t)out.shape(2) != nx)
        throw std::invalid_argument("median27_symmetric: A and out differ in size");
    const float* a = A.data();
    const idx_t n = nz * ny * nx;
    nb::gil_scoped_release release;
    int bad = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) schedule(static) reduction(| : bad)
#endif
    for (idx_t i = 0; i < n; i++) bad |= (a[i] != a[i]) || (a[i] == 0 && std::signbit(a[i]));
    if (bad) return false;
    if (n == 0) return true;
    const idx_t pz = nz + 2, py = ny + 2, px = nx + 2;
    std::unique_ptr<float[]> P(new float[(size_t)(pz * py * px + ML)]);
    float* p = P.get();
    std::memset(p + pz * py * px, 0, sizeof(float) * ML);
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads) schedule(static)
#endif
    for (idx_t r = 0; r < pz * py; r++) {
        const idx_t z = std::min(std::max<idx_t>(r / py - 1, 0), nz - 1);
        const idx_t y = std::min(std::max<idx_t>(r % py - 1, 0), ny - 1);
        const float* s = a + (z * ny + y) * nx;
        float* d = p + r * px;
        d[0] = s[0];
        std::memcpy(d + 1, s, sizeof(float) * (size_t)nx);
        d[nx + 1] = s[nx - 1];
    }
    median27_filter(p, out.data(), nx, ny, nz, num_threads, isa);
    return true;
}

}  // namespace

PK5D_BINDINGS(bind_median) {
    m.def("median27_symmetric_single", &median27_symmetric_bind, nb::arg("a").noconvert(),
          nb::arg("out").noconvert(), nb::arg("num_threads"), nb::arg("isa") = -1,
          "medfilt3(A) of a C-contiguous (z, y, x) single volume (3 x 3 x 3, symmetric padding) "
          "with a selection network; False (out untouched) when A has NaN or -0. isa: the "
          "network's instruction set (median27_best_isa; -1 the best one).");
    const char* doc =
        "out(z, y, x) = the k-th smallest (0-based) of the wz x wy x wx window of the padded "
        "C-contiguous volume P at (z, y, x) (P: out's size + window - 1); a window with NaN gives "
        "its first NaN in MATLAB's order, a zero result is +0.";
#define PK5D_ORDER(T, name)                                                                         \
    m.def(name, &order_filter_bind<T>, nb::arg("p").noconvert(), nb::arg("out").noconvert(), \
          nb::arg("wz"), nb::arg("wy"), nb::arg("wx"), nb::arg("k"), nb::arg("num_threads"), doc);
    PK5D_ORDER(bool, "order_filter_bool")
    PK5D_ORDER(uint8_t, "order_filter_uint8")
    PK5D_ORDER(int8_t, "order_filter_int8")
    PK5D_ORDER(uint16_t, "order_filter_uint16")
    PK5D_ORDER(int16_t, "order_filter_int16")
    PK5D_ORDER(uint32_t, "order_filter_uint32")
    PK5D_ORDER(int32_t, "order_filter_int32")
    PK5D_ORDER(float, "order_filter_single")
    PK5D_ORDER(double, "order_filter_double")
#undef PK5D_ORDER
    m.def("median27_single", &median27_bind, nb::arg("p").noconvert(), nb::arg("out").noconvert(),
          nb::arg("num_threads"), nb::arg("isa") = -1,
          "The 3 x 3 x 3 median of a single volume padded by 1 (P: out's size + 2, with 16 floats "
          "of room after its data) with a selection network; False (out untouched) when P has NaN "
          "or -0. isa: the network's instruction set (median27_best_isa; -1 the best one).");
    m.def("median27_best_isa", &median27_best_isa,
          "The instruction set the median network runs: 0 the baseline (SSE on x86-64), 1 AVX2.");
}
