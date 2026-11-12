// The Euclidean distance transform of bwdist: for each voxel of a C-contiguous (z, y, x) bool
// volume the distance to the nearest true voxel, as single: sqrt(single(d2)) of the exact
// squared distance d2 (Inf without true voxels). Separable and exact in integers: the squared
// distance to the nearest feature along x (the contiguous axis), then the lower envelope of the
// parabolas (u - i)^2 + f(i) along y and along z (Meijster, Roerdink and Hesselink, 2000, with
// integer floor division, so that no rounding enters). The squared distance is unique, so the
// result doesn't depend on the order of the passes or the thread count.
#include <pk5d/common.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

namespace nb = nanobind;
using idx_t = std::ptrdiff_t;

inline int64_t floordiv(const int64_t a, const int64_t b) {  // b > 0
    int64_t q = a / b;
    if ((a % b) != 0 && a < 0) --q;
    return q;
}

// d[x] = min over i with f[i] < INF of (x - i)^2 + f[i], for x = 0 .. n - 1 (INF if no such i).
// s, t: scratch of n entries each.
template <typename T>
void envelope(const T* f, T* d, const idx_t n, const T INF, idx_t* s, idx_t* t) {
    idx_t q = -1;
    auto F = [&](const idx_t x, const idx_t i) -> int64_t {
        const int64_t dx = (int64_t)(x - i);
        return dx * dx + (int64_t)f[i];
    };
    for (idx_t u = 0; u < n; ++u) {
        if (f[u] == INF) continue;
        while (q >= 0 && F(t[q], s[q]) > F(t[q], u)) --q;
        if (q < 0) {
            q = 0;
            s[0] = u;
            t[0] = 0;
        } else {
            const int64_t i = (int64_t)s[q], uu = (int64_t)u;
            const int64_t w = 1 + floordiv(uu * uu - i * i + (int64_t)f[u] - (int64_t)f[i], 2 * (uu - i));
            if (w < (int64_t)n) {
                ++q;
                s[q] = u;
                t[q] = (idx_t)w;
            }
        }
    }
    if (q < 0) {
        for (idx_t x = 0; x < n; ++x) d[x] = INF;
        return;
    }
    for (idx_t x = n - 1; x >= 0; --x) {
        d[x] = (T)F(x, s[q]);
        if (x == t[q]) --q;
    }
}

// squared distance along the contiguous axis (x) within each line
template <typename T>
void pass_contig(const uint8_t* A, T* W, const idx_t n1, const idx_t ncol, const T INF, const int nt) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static)
#endif
    for (idx_t c = 0; c < ncol; ++c) {
        const uint8_t* a = A + c * n1;
        T* w = W + c * n1;
        idx_t last = -1;
        for (idx_t y = 0; y < n1; ++y) {
            if (a[y]) last = y;
            w[y] = last < 0 ? INF : (T)(y - last);
        }
        last = -1;
        for (idx_t y = n1 - 1; y >= 0; --y) {
            if (a[y]) last = y;
            if (last >= 0 && (w[y] == INF || (T)(last - y) < w[y])) w[y] = (T)(last - y);
        }
        for (idx_t y = 0; y < n1; ++y) {
            if (w[y] != INF) w[y] = w[y] * w[y];
        }
    }
}

template <typename T>
inline float distance_of(const T d2, const T INF) {
    return d2 == INF ? std::numeric_limits<float>::infinity() : std::sqrt((float)d2);
}

// the envelope along an axis of stride `stride` and length n, for the lines that start at
// o * outer_stride + y (o < nouter, y < n1); lines are taken in tiles of L consecutive starts
// (consecutive in memory), so that the gather reads whole cache lines (L smaller when there are
// few tiles for the threads). With out, the last pass writes the distances (sqrt(single(d2)))
// there instead of the squared distances into W; W may then be out's own memory (a tile reads
// all its values before it writes any, and both go through memcpy, which may alias).
template <typename T>
void pass_axis(T* W, const idx_t n, const idx_t stride, const idx_t n1, const idx_t nouter,
               const idx_t outer_stride, const T INF, const int nt, float* out = nullptr) {
    idx_t L = 16;
    while (L > 4 && ((n1 + L - 1) / L) * nouter < 8 * (idx_t)nt) L /= 2;
    const idx_t ntiles_y = (n1 + L - 1) / L;
    const idx_t ntiles = ntiles_y * nouter;
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
    {
        std::vector<T> f((size_t)(L * n)), d((size_t)n);
        std::vector<idx_t> s((size_t)n), t((size_t)n);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (idx_t tile = 0; tile < ntiles; ++tile) {
            const idx_t y0 = (tile % ntiles_y) * L, o = tile / ntiles_y;
            const idx_t ly = std::min(L, n1 - y0);
            const idx_t off = o * outer_stride + y0;
            for (idx_t i = 0; i < n; ++i) {
                const unsigned char* src = reinterpret_cast<const unsigned char*>(W + off + i * stride);
                for (idx_t k = 0; k < ly; ++k) std::memcpy(&f[(size_t)(k * n + i)], src + k * sizeof(T), sizeof(T));
            }
            for (idx_t k = 0; k < ly; ++k) {
                envelope<T>(f.data() + k * n, d.data(), n, INF, s.data(), t.data());
                std::memcpy(f.data() + k * n, d.data(), sizeof(T) * (size_t)n);
            }
            if (out != nullptr) {
                for (idx_t i = 0; i < n; ++i) {
                    unsigned char* dst = reinterpret_cast<unsigned char*>(out + off + i * stride);
                    for (idx_t k = 0; k < ly; ++k) {
                        const float v = distance_of<T>(f[(size_t)(k * n + i)], INF);
                        std::memcpy(dst + k * sizeof(float), &v, sizeof(float));
                    }
                }
            } else {
                for (idx_t i = 0; i < n; ++i) {
                    T* dst = W + off + i * stride;
                    for (idx_t k = 0; k < ly; ++k) dst[k] = f[(size_t)(k * n + i)];
                }
            }
        }
    }
}

template <typename T>
void edt(const uint8_t* A, float* out, const idx_t n1, const idx_t n2, const idx_t n3, const int nt) {
    const T INF = std::numeric_limits<T>::max();
    const idx_t nv = n1 * n2 * n3;
    // squared distances of 4 bytes are computed in out's memory (no allocation, no page faults
    // beyond out's own); wider ones in a buffer of their own
    std::unique_ptr<T[]> own;
    T* W;
    if (sizeof(T) == sizeof(float)) {
        W = reinterpret_cast<T*>(out);
    } else {
        own.reset(new T[(size_t)nv]);  // every element is written by the first pass
        W = own.get();
    }
    pass_contig<T>(A, W, n1, n2 * n3, INF, nt);
    if (n3 > 1) {
        if (n2 > 1) pass_axis<T>(W, n2, n1, n1, n3, n1 * n2, INF, nt);
        pass_axis<T>(W, n3, n1 * n2, n1, n2, n1, INF, nt, out);
    } else if (n2 > 1) {
        pass_axis<T>(W, n2, n1, n1, n3, n1 * n2, INF, nt, out);
    } else {
        for (idx_t i = 0; i < nv; ++i) {
            T v;
            std::memcpy(&v, reinterpret_cast<const unsigned char*>(W + i), sizeof(T));
            const float d = distance_of<T>(v, INF);
            std::memcpy(reinterpret_cast<unsigned char*>(out + i), &d, sizeof(float));
        }
    }
}

using cbvol = pk5d::cvolume<const bool>;
using fvol = pk5d::cvolume<float>;

void bwdist_edt(cbvol a, fvol out, int num_threads) {
    pk5d::check_threads(num_threads);
    for (int d = 0; d < 3; ++d) {
        if (a.shape(d) != out.shape(d)) throw std::invalid_argument("bwdist_edt: a and out differ in size");
    }
    const idx_t n1 = (idx_t)a.shape(2), n2 = (idx_t)a.shape(1), n3 = (idx_t)a.shape(0);  // x, y, z
    if (n1 * n2 * n3 == 0) return;
    const uint8_t* A = reinterpret_cast<const uint8_t*>(a.data());
    float* O = out.data();
    // the largest squared distance, for the storage type of the squared distances
    const double dmax = (double)(n1 - 1) * (n1 - 1) + (double)(n2 - 1) * (n2 - 1) + (double)(n3 - 1) * (n3 - 1);
    nb::gil_scoped_release release;
    if (dmax < 4294967295.0) {
        edt<uint32_t>(A, O, n1, n2, n3, num_threads);
    } else {
        edt<uint64_t>(A, O, n1, n2, n3, num_threads);
    }
}

}  // namespace

PK5D_BINDINGS(bind_edt) {
    m.def("bwdist_edt", &bwdist_edt, nb::arg("a").noconvert(), nb::arg("out").noconvert(),
          nb::arg("num_threads"),
          "Euclidean distance of every voxel of a C-contiguous (z, y, x) bool volume to the nearest true voxel, "
          "sqrt(single(d2)) of the exact squared distance (Inf without true voxels), into out "
          "(single).");
}
