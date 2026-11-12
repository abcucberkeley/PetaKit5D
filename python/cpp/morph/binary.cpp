// Binary erosion and dilation by a sequence of neighborhoods, as MATLAB's morphop applies a
// (decomposed) structuring element: the image is padded by pad_before / pad_after with the
// outside value (true for erosion, false for dilation), eroded or dilated by each neighborhood
// in turn over the whole padded domain, with pixels outside it counting as the outside value,
// and cropped back. Erosion: B(p) = AND of A(p + o) over the offsets o; dilation:
// B(p) = OR of A(p - o) (the reflected neighborhood). Volumes are C-contiguous (z, y, x);
// offsets are rows (dz, dy, dx) from the neighborhood's center (floor((size + 1) / 2) in
// MATLAB's terms), pads (z, y, x). Each output voxel is an AND or OR of a fixed set of voxels,
// so the visiting order doesn't matter.
//
// Inside, memory axis 0 is x (contiguous), 1 is y, 2 is z. Neighborhoods without dz offsets
// and without padding along z act on each z plane alone (planes per thread, or the columns of a
// plane split among the threads when there are few planes); others on the whole volume. The
// result doesn't depend on the thread count.
#include <pk5d/common.hpp>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

namespace nb = nanobind;
using idx_t = std::ptrdiff_t;

struct Plan {
    idx_t n[3];             // image size along memory axes 0, 1, 2 (x, y, z)
    idx_t pb[3], pa[3];     // padding before / after
    idx_t P[3];             // padded domain
    idx_t G[3];             // guard band (outside value, never written)
    idx_t B[3];             // buffer size P + 2 G
    idx_t s1, s2;           // buffer strides of memory axes 1 and 2
    bool erode;
    bool slices;            // each z slice alone (buffers of one slice)
    std::vector<idx_t> start;  // step k: deltas[start[k] .. start[k + 1])
    std::vector<idx_t> delta;  // linear buffer offsets of the neighbors to read
};

#if defined(__GNUC__) || defined(__clang__)
// 16 bytes (SSE2 or NEON registers) combined over all the shifted lines, then stored once
typedef uint8_t v16 __attribute__((vector_size(16)));
#define PK5D_MORPH_V16 1
#endif

// one line (padded domain) of a step: out = AND / OR of in shifted by the deltas
template <bool ERODE>
inline void column_step(const uint8_t* in, uint8_t* out, const idx_t len, const idx_t* d, const idx_t nd) {
    if (nd == 0) {
        std::memset(out, ERODE ? 1 : 0, (size_t)len);
        return;
    }
    idx_t i = 0;
#ifdef PK5D_MORPH_V16
    for (; i + 16 <= len; i += 16) {
        v16 v;
        std::memcpy(&v, in + d[0] + i, 16);
        for (idx_t k = 1; k < nd; ++k) {
            v16 t;
            std::memcpy(&t, in + d[k] + i, 16);
            v = ERODE ? (v & t) : (v | t);
        }
        std::memcpy(out + i, &v, 16);
    }
#else
    // without vector extensions: the first line copied, the others combined in place
    if (len >= 16) {
        std::memcpy(out, in + d[0], (size_t)len);
        for (idx_t k = 1; k < nd; ++k) {
            const uint8_t* t = in + d[k];
            for (idx_t j = 0; j < len; ++j) out[j] = ERODE ? (uint8_t)(out[j] & t[j]) : (uint8_t)(out[j] | t[j]);
        }
        return;
    }
#endif
    for (; i < len; ++i) {
        uint8_t v = in[d[0] + i];
        for (idx_t k = 1; k < nd; ++k) v = ERODE ? (uint8_t)(v & in[d[k] + i]) : (uint8_t)(v | in[d[k] + i]);
        out[i] = v;
    }
}

inline void column_step(const Plan& p, const int k, const uint8_t* in, uint8_t* out, const idx_t base) {
    const idx_t* d = p.delta.data() + p.start[k];
    const idx_t nd = p.start[k + 1] - p.start[k];
    if (p.erode) {
        column_step<true>(in + base, out + base, p.P[0], d, nd);
    } else {
        column_step<false>(in + base, out + base, p.P[0], d, nd);
    }
}

// buffer index of padded-domain voxel (0, x, z) of memory axes (0, 1, 2) (z = 0 in plane buffers)
inline idx_t col_base(const Plan& p, const idx_t x, const idx_t z) {
    return (z + p.G[2]) * p.s2 + (x + p.G[1]) * p.s1 + p.G[0];
}

// padded-domain line (x, z) (memory axes 1 and 2) of the buffer from the image (z: padded-domain plane; zi: the
// image plane it holds, or -1 for padding)
inline void load_col(const Plan& p, const uint8_t* A, uint8_t* buf, const idx_t x, const idx_t zb, const idx_t zi) {
    uint8_t* col = buf + col_base(p, x, zb);
    const uint8_t v0 = p.erode ? 1 : 0;
    const idx_t xi = x - p.pb[1];
    if (zi < 0 || xi < 0 || xi >= p.n[1]) {
        std::memset(col, v0, (size_t)p.P[0]);
        return;
    }
    std::memset(col, v0, (size_t)p.pb[0]);
    std::memset(col + p.pb[0] + p.n[0], v0, (size_t)p.pa[0]);
    // NumPy's bool arrays hold 0 and 1 only
    std::memcpy(col + p.pb[0], A + (zi * p.n[1] + xi) * p.n[0], (size_t)p.n[0]);
}

inline void store_col(const Plan& p, const uint8_t* buf, uint8_t* Bout, const idx_t xi, const idx_t zb, const idx_t zi) {
    std::memcpy(Bout + (zi * p.n[1] + xi) * p.n[0], buf + col_base(p, xi + p.pb[1], zb) + p.pb[0], (size_t)p.n[0]);
}

// all steps of slice z, serially, with the buffers a and b
void run_slice(const Plan& p, const uint8_t* A, uint8_t* Bout, const idx_t z, uint8_t* a, uint8_t* b) {
    for (idx_t x = 0; x < p.P[1]; ++x) load_col(p, A, a, x, 0, z);
    const int ns = (int)p.start.size() - 1;
    for (int k = 0; k < ns; ++k) {
        for (idx_t x = 0; x < p.P[1]; ++x) column_step(p, k, a, b, col_base(p, x, 0));
        std::swap(a, b);
    }
    for (idx_t x = 0; x < p.n[1]; ++x) store_col(p, a, Bout, x, 0, z);
}

void run(const Plan& p, const uint8_t* A, uint8_t* Bout, const int nt) {
    const int ns = (int)p.start.size() - 1;
    const uint8_t v0 = p.erode ? 1 : 0;
    if (p.slices) {
        const size_t bsz = (size_t)(p.B[0] * p.B[1]);
        const bool per_thread = nt == 1 || p.n[2] >= 2 * (idx_t)nt;
        if (per_thread) {
            // each thread makes its own plane buffers (in parallel, in its memory)
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
            {
                std::unique_ptr<uint8_t[]> a(new uint8_t[bsz]), b(new uint8_t[bsz]);
                std::memset(a.get(), v0, bsz);
                std::memset(b.get(), v0, bsz);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
                for (idx_t z = 0; z < p.n[2]; ++z) run_slice(p, A, Bout, z, a.get(), b.get());
            }
            return;
        }
        std::vector<uint8_t> abuf(bsz, v0), bbuf(bsz, v0);
        for (idx_t z = 0; z < p.n[2]; ++z) {
            uint8_t* a = abuf.data();
            uint8_t* b = bbuf.data();
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
            {
                uint8_t* ta = a;
                uint8_t* tb = b;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
                for (idx_t x = 0; x < p.P[1]; ++x) load_col(p, A, ta, x, 0, z);
                for (int k = 0; k < ns; ++k) {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
                    for (idx_t x = 0; x < p.P[1]; ++x) column_step(p, k, ta, tb, col_base(p, x, 0));
                    std::swap(ta, tb);
                }
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
                for (idx_t x = 0; x < p.n[1]; ++x) store_col(p, ta, Bout, x, 0, z);
            }
        }
        return;
    }
    // whole volume
    const size_t bsz = (size_t)(p.B[0] * p.B[1] * p.B[2]);
    std::unique_ptr<uint8_t[]> abuf(new uint8_t[bsz]), bbuf(new uint8_t[bsz]);
    uint8_t* a = abuf.get();
    uint8_t* b = bbuf.get();
    const idx_t ncols = p.P[1] * p.P[2];
    const idx_t nchunks = 64 * (idx_t)nt;
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
    {
        // the outside value everywhere (guard bands included), filled in parallel
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (idx_t c = 0; c < nchunks; ++c) {
            const size_t c0 = bsz * (size_t)c / (size_t)nchunks, c1 = bsz * (size_t)(c + 1) / (size_t)nchunks;
            std::memset(a + c0, v0, c1 - c0);
            std::memset(b + c0, v0, c1 - c0);
        }
        uint8_t* ta = a;
        uint8_t* tb = b;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (idx_t c = 0; c < ncols; ++c) {
            const idx_t x = c % p.P[1], z = c / p.P[1];
            const idx_t zi = z - p.pb[2];
            load_col(p, A, ta, x, z, zi >= 0 && zi < p.n[2] ? zi : -1);
        }
        for (int k = 0; k < ns; ++k) {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (idx_t c = 0; c < ncols; ++c) {
                column_step(p, k, ta, tb, col_base(p, c % p.P[1], c / p.P[1]));
            }
            std::swap(ta, tb);
        }
        const idx_t nimg = p.n[1] * p.n[2];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (idx_t c = 0; c < nimg; ++c) {
            const idx_t x = c % p.n[1], z = c / p.n[1];
            store_col(p, ta, Bout, x, z + p.pb[2], z);
        }
    }
}

using bvol = pk5d::cvolume<bool>;
using cbvol = pk5d::cvolume<const bool>;
using offsets_t = nb::ndarray<const int64_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using ivec = nb::ndarray<const int64_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

void binary_morph(cbvol a, bvol out, offsets_t offsets, ivec counts, bool erode, ivec pad_before,
                  ivec pad_after, int num_threads) {
    pk5d::check_threads(num_threads);
    for (int d = 0; d < 3; ++d) {
        if (a.shape(d) != out.shape(d)) throw std::invalid_argument("binary_morph: a and out differ in size");
    }
    if (offsets.shape(1) != 3) throw std::invalid_argument("binary_morph: offsets must have 3 columns");
    if (pad_before.shape(0) != 3 || pad_after.shape(0) != 3)
        throw std::invalid_argument("binary_morph: pads must have 3 elements");
    Plan p;
    p.erode = erode;
    for (int d = 0; d < 3; ++d) {
        p.n[d] = (idx_t)a.shape(2 - d);
        p.pb[d] = (idx_t)pad_before(2 - d);
        p.pa[d] = (idx_t)pad_after(2 - d);
        if (p.pb[d] < 0 || p.pa[d] < 0) throw std::invalid_argument("binary_morph: negative padding");
        p.P[d] = p.n[d] + p.pb[d] + p.pa[d];
        p.G[d] = 0;
    }
    const idx_t nsteps = (idx_t)counts.shape(0);
    idx_t total = 0;
    for (idx_t k = 0; k < nsteps; ++k) {
        if (counts(k) < 0) throw std::invalid_argument("binary_morph: negative count");
        total += (idx_t)counts(k);
    }
    if (total != (idx_t)offsets.shape(0)) throw std::invalid_argument("binary_morph: counts don't add up to the offsets");
    bool has_dz = false;
    for (idx_t i = 0; i < total; ++i) {
        for (int d = 0; d < 3; ++d) {
            const idx_t o = (idx_t)offsets(i, 2 - d);
            p.G[d] = std::max(p.G[d], o < 0 ? -o : o);
        }
        has_dz = has_dz || offsets(i, 0) != 0;
    }
    if (p.n[0] * p.n[1] * p.n[2] == 0) return;
    p.slices = !has_dz && p.pb[2] == 0 && p.pa[2] == 0;
    if (p.slices) p.G[2] = 0;
    for (int d = 0; d < 3; ++d) p.B[d] = p.P[d] + 2 * p.G[d];
    if (p.slices) p.B[2] = 1;
    p.s1 = p.B[0];
    p.s2 = p.B[0] * p.B[1];
    p.start.push_back(0);
    idx_t i = 0;
    const idx_t sign = erode ? 1 : -1;
    for (idx_t k = 0; k < nsteps; ++k) {
        for (idx_t j = 0; j < (idx_t)counts(k); ++j, ++i) {
            p.delta.push_back(sign * ((idx_t)offsets(i, 2) + (idx_t)offsets(i, 1) * p.s1 + (idx_t)offsets(i, 0) * p.s2));
        }
        p.start.push_back((idx_t)p.delta.size());
    }
    const uint8_t* A = reinterpret_cast<const uint8_t*>(a.data());
    uint8_t* B = reinterpret_cast<uint8_t*>(out.data());
    nb::gil_scoped_release release;
    run(p, A, B, num_threads);
}

}  // namespace

PK5D_BINDINGS(bind_binary) {
    m.def("binary_morph", &binary_morph, nb::arg("a").noconvert(), nb::arg("out").noconvert(),
          nb::arg("offsets").noconvert(), nb::arg("counts").noconvert(), nb::arg("erode"),
          nb::arg("pad_before").noconvert(), nb::arg("pad_after").noconvert(), nb::arg("num_threads"),
          "Binary erosion (AND, outside true) or dilation (OR of the reflected offsets, outside "
          "false) of a C-contiguous (z, y, x) bool volume by a sequence of neighborhoods "
          "(offsets: rows (dz, dy, dx), counts: rows per step), on the volume padded by "
          "pad_before / pad_after ((z, y, x)) with the outside value and cropped back, as "
          "MATLAB's morphop; writes out.");
}
