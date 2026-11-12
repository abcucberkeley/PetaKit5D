// Connected components of C-contiguous (z, y, x) bool volumes, numbered as MATLAB numbers them.
// MATLAB's linear order is y fastest, then x, then z: voxel (z, y, x) has MATLAB index
// t = y + ny * (x + nx * z).
//
// label_first: components in the order of their first voxel in MATLAB's order (bwconncomp,
//   bwlabel). Union-find over C indices with the smaller index as the root, in slabs of z planes
//   per thread joined across the slabs; one pass in C order numbers the roots; then the
//   components are renumbered by their first voxel in MATLAB's order (a component's first
//   voxel is in its first z plane, so both orders agree on the plane). Doesn't depend on the
//   thread count.
// label_bwlabeln: bwlabeln's numbering (measured against MATLAB R2024b): voxels in MATLAB's order;
//   a voxel takes the provisional label of its first labeled earlier neighbor (the neighbors in
//   the column-major order of MATLAB's connectivity matrix) and the sets of its other labeled
//   earlier neighbors are united with its set in that order, by size in provisional labels
//   (ties: the voxel's set keeps its root); a voxel without labeled earlier neighbors gets a new
//   provisional label; the roots are numbered in increasing provisional order. Serial (the
//   result depends on the union order); it runs on a copy of the mask in MATLAB's layout
//   ((z, x, y) planes, made by a blocked transpose per plane) and writes the labels back the
//   same way.
#include <pk5d/common.hpp>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {

namespace nb = nanobind;
using idx_t = std::ptrdiff_t;

// a neighbor offset along memory axes 0, 1, 2 (fastest first)
struct Nb {
    int d0, d1, d2;
};

template <typename I>
inline I find_root(I* parent, I i) {
    while (parent[i] != i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}

// union by the smaller index: the root of a component is its first voxel in memory order (and a
// parent always comes before its child in memory order)
template <typename I>
inline void unite_min(I* parent, const I a, const I b) {
    const I ra = find_root(parent, a), rb = find_root(parent, b);
    if (ra < rb) {
        parent[rb] = ra;
    } else if (rb < ra) {
        parent[ra] = rb;
    }
}

// labels in memory order (1-based by first voxel in memory order; 0 background); n0, n1, n2: the
// sizes along memory axes 0, 1, 2; nb: the earlier neighbors in memory order
template <typename I>
int64_t label_memory_order(const uint8_t* A, I* lab, const idx_t n0, const idx_t n1, const idx_t n2,
                           const std::vector<Nb>& nb, const int nthreads) {
    const idx_t np = n0 * n1;
    const int nt = (int)std::max<idx_t>(1, std::min<idx_t>(nthreads, n2));
    std::vector<idx_t> zs(nt + 1);
    for (int t = 0; t <= nt; t++) zs[t] = n2 * t / nt;
    I* parent = lab;
    auto link = [&](const idx_t k2, const idx_t k1, const idx_t k0, const idx_t kmin) {
        const I i = (I)(k2 * np + k1 * n0 + k0);
        for (const Nb& o : nb) {
            const idx_t a0 = k0 + o.d0, a1 = k1 + o.d1, a2 = k2 + o.d2;
            if (a0 < 0 || a0 >= n0 || a1 < 0 || a1 >= n1 || a2 < kmin || a2 >= n2) continue;
            const I j = (I)(a2 * np + a1 * n0 + a0);
            if (A[j]) unite_min(parent, i, j);
        }
    };
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static, 1)
#endif
    for (int t = 0; t < nt; t++) {
        for (idx_t k2 = zs[t]; k2 < zs[t + 1]; k2++) {
            for (idx_t k1 = 0; k1 < n1; k1++) {
                for (idx_t k0 = 0; k0 < n0; k0++) {
                    const idx_t i = k2 * np + k1 * n0 + k0;
                    parent[i] = A[i] ? (I)i : (I)-1;
                    if (A[i]) link(k2, k1, k0, zs[t]);
                }
            }
        }
    }
    // the first plane of each slab with its earlier neighbors in the slab before
    for (int t = 1; t < nt; t++) {
        const idx_t k2 = zs[t];
        for (idx_t k1 = 0; k1 < n1; k1++) {
            for (idx_t k0 = 0; k0 < n0; k0++) {
                if (A[k2 * np + k1 * n0 + k0]) link(k2, k1, k0, k2 - 1);
            }
        }
    }
    // in memory order a voxel's parent comes first and is already resolved: its entry holds
    // -(label) of its root
    const idx_t n = np * n2;
    I nl = 0;
    for (idx_t i = 0; i < n; i++) {
        const I p = parent[i];
        if (p < 0) {
            lab[i] = 0;
            continue;
        }
        const I l = (p == (I)i) ? ++nl : -lab[p];
        lab[i] = -l;
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
    for (idx_t i = 0; i < n; i++) lab[i] = -lab[i];
    return (int64_t)nl;
}

// renumber the labels of a C (nz, ny, nx) label volume by first voxel in MATLAB's order
template <typename I>
void renumber_matlab_order(I* lab, const idx_t nz, const idx_t ny, const idx_t nx, const int64_t num,
                           const int nt) {
    if (num <= 1) return;
    std::vector<I> perm((size_t)num + 1, 0);
    I next = 0;
    for (idx_t z = 0; z < nz && next < (I)num; z++) {
        const I* p = lab + z * ny * nx;
        for (idx_t x = 0; x < nx; x++) {
            for (idx_t y = 0; y < ny; y++) {
                const I l = p[y * nx + x];
                if (l != 0 && perm[(size_t)l] == 0) perm[(size_t)l] = ++next;
            }
        }
    }
    const idx_t n = nz * ny * nx;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static)
#endif
    for (idx_t i = 0; i < n; i++) lab[i] = perm[(size_t)lab[i]];
}

// bwlabeln on a mask in MATLAB's layout (memory axes y, x, z), provisional and then final labels
// into P (same layout); returns the number of components and the final label of each
// provisional label in fin
template <typename I>
int64_t bwlabeln_scan(const uint8_t* A, I* P, const idx_t n0, const idx_t n1, const idx_t n2,
                      const std::vector<Nb>& nb, std::vector<I>& fin) {
    const idx_t np = n0 * n1;
    std::vector<I> parent(1, 0), size(1, 0);
    auto root = [&](I a) {
        while (parent[a] != a) {
            parent[a] = parent[parent[a]];
            a = parent[a];
        }
        return a;
    };
    for (idx_t k2 = 0; k2 < n2; k2++) {
        for (idx_t k1 = 0; k1 < n1; k1++) {
            for (idx_t k0 = 0; k0 < n0; k0++) {
                const idx_t i = k2 * np + k1 * n0 + k0;
                if (!A[i]) {
                    P[i] = 0;
                    continue;
                }
                I l = 0;
                for (const Nb& o : nb) {
                    const idx_t a0 = k0 + o.d0, a1 = k1 + o.d1, a2 = k2 + o.d2;
                    if (a0 < 0 || a0 >= n0 || a1 < 0 || a1 >= n1 || a2 < 0 || a2 >= n2) continue;
                    const I m = P[a2 * np + a1 * n0 + a0];
                    if (m == 0) continue;
                    if (l == 0) {
                        l = m;
                        continue;
                    }
                    const I ra = root(l), rb = root(m);
                    if (ra == rb) continue;
                    if (size[ra] >= size[rb]) {
                        parent[rb] = ra;
                        size[ra] += size[rb];
                    } else {
                        parent[ra] = rb;
                        size[rb] += size[ra];
                    }
                }
                if (l == 0) {
                    l = (I)parent.size();
                    parent.push_back(l);
                    size.push_back(1);
                }
                P[i] = l;
            }
        }
    }
    // roots numbered in increasing provisional order
    const I nprov = (I)parent.size();
    fin.assign((size_t)nprov, 0);
    I nl = 0;
    for (I k = 1; k < nprov; k++) {
        if (root(k) == k) fin[(size_t)k] = ++nl;
    }
    for (I k = 1; k < nprov; k++) fin[(size_t)k] = fin[(size_t)root(k)];
    return (int64_t)nl;
}

// blocked transpose of each z plane: src (nz, ny, nx) C to dst (nz, nx, ny) C, or back
template <typename S, typename D, typename F>
void transpose_planes(const S* src, D* dst, const idx_t nz, const idx_t ny, const idx_t nx, const int nt, F f) {
    const idx_t B = 32;
    const idx_t nby = (ny + B - 1) / B, nbx = (nx + B - 1) / B;
    const idx_t ntiles = nz * nby * nbx;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static)
#endif
    for (idx_t t = 0; t < ntiles; t++) {
        const idx_t z = t / (nby * nbx), r = t % (nby * nbx);
        const idx_t y0 = (r / nbx) * B, x0 = (r % nbx) * B;
        const idx_t y1 = std::min(ny, y0 + B), x1 = std::min(nx, x0 + B);
        const S* s = src + z * ny * nx;
        D* d = dst + z * ny * nx;
        for (idx_t y = y0; y < y1; y++) {
            for (idx_t x = x0; x < x1; x++) d[x * ny + y] = f(s[y * nx + x]);
        }
    }
}

template <typename D, typename I>
void untranspose_labels(const I* P, D* out, const std::vector<I>& fin, const idx_t nz, const idx_t ny,
                        const idx_t nx, const int nt) {
    const idx_t B = 32;
    const idx_t nby = (ny + B - 1) / B, nbx = (nx + B - 1) / B;
    const idx_t ntiles = nz * nby * nbx;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static)
#endif
    for (idx_t t = 0; t < ntiles; t++) {
        const idx_t z = t / (nby * nbx), r = t % (nby * nbx);
        const idx_t y0 = (r / nbx) * B, x0 = (r % nbx) * B;
        const idx_t y1 = std::min(ny, y0 + B), x1 = std::min(nx, x0 + B);
        const I* s = P + z * ny * nx;
        D* d = out + z * ny * nx;
        for (idx_t y = y0; y < y1; y++) {
            for (idx_t x = x0; x < x1; x++) d[y * nx + x] = (D)fin[(size_t)s[x * ny + y]];
        }
    }
}

using cbvol = pk5d::cvolume<const bool>;
using offsets_t = nb::ndarray<const int64_t, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
template <typename T>
using cvol = pk5d::cvolume<T>;

// rows of 3 offsets, in the order given; reverse: rows (a, b, c) are axes 2, 1, 0 of memory
std::vector<Nb> parse_offsets(offsets_t offsets, const bool reverse) {
    if (offsets.shape(1) != 3) throw std::invalid_argument("offsets must have 3 columns");
    std::vector<Nb> nb((size_t)offsets.shape(0));
    for (size_t k = 0; k < nb.size(); ++k) {
        const int a = (int)offsets(k, 0), b = (int)offsets(k, 1), c = (int)offsets(k, 2);
        if (a < -1 || a > 1 || b < -1 || b > 1 || c < -1 || c > 1)
            throw std::invalid_argument("offsets must be within -1 .. 1");
        nb[k] = reverse ? Nb{c, b, a} : Nb{a, b, c};
    }
    return nb;
}

template <typename I>
void check_sizes(cbvol a, const idx_t* shape, const char* name) {
    for (int d = 0; d < 3; ++d) {
        if ((idx_t)a.shape(d) != shape[d]) throw std::invalid_argument(std::string(name) + ": sizes differ");
    }
    const double n = (double)a.shape(0) * (double)a.shape(1) * (double)a.shape(2);
    if (n >= (double)std::numeric_limits<I>::max()) throw std::invalid_argument(std::string(name) + ": too many voxels");
}

// first-voxel numbering (MATLAB's order); offsets: rows (dz, dy, dx), the earlier neighbors in C order
template <typename I>
int64_t label_first_bind(cbvol a, cvol<I> lab, offsets_t offsets, int num_threads) {
    pk5d::check_threads(num_threads);
    const idx_t shape[3] = {(idx_t)lab.shape(0), (idx_t)lab.shape(1), (idx_t)lab.shape(2)};
    check_sizes<I>(a, shape, "label_first");
    const std::vector<Nb> nb = parse_offsets(offsets, true);
    const idx_t nz = shape[0], ny = shape[1], nx = shape[2];
    if (nz * ny * nx == 0) return 0;
    const uint8_t* A = reinterpret_cast<const uint8_t*>(a.data());
    I* L = lab.data();
    nb::gil_scoped_release release;
    const int64_t num = label_memory_order<I>(A, L, nx, ny, nz, nb, num_threads);
    renumber_matlab_order<I>(L, nz, ny, nx, num, num_threads);
    return num;
}

// labels in memory (C) order, any numbering (sizes, holes); offsets as label_first_bind
template <typename I>
int64_t label_any_bind(cbvol a, cvol<I> lab, offsets_t offsets, int num_threads) {
    pk5d::check_threads(num_threads);
    const idx_t shape[3] = {(idx_t)lab.shape(0), (idx_t)lab.shape(1), (idx_t)lab.shape(2)};
    check_sizes<I>(a, shape, "label_any");
    const std::vector<Nb> nb = parse_offsets(offsets, true);
    if (shape[0] * shape[1] * shape[2] == 0) return 0;
    const uint8_t* A = reinterpret_cast<const uint8_t*>(a.data());
    nb::gil_scoped_release release;
    return label_memory_order<I>(A, lab.data(), shape[2], shape[1], shape[0], nb, num_threads);
}

// bwlabeln's numbering into a double label volume; offsets: rows (dy, dx, dz), MATLAB's earlier
// neighbors in the column-major order of its connectivity matrix
template <typename I>
int64_t bwlabeln_bind(cbvol a, cvol<double> out, offsets_t offsets, int num_threads) {
    pk5d::check_threads(num_threads);
    const idx_t shape[3] = {(idx_t)out.shape(0), (idx_t)out.shape(1), (idx_t)out.shape(2)};
    check_sizes<I>(a, shape, "bwlabeln");
    const std::vector<Nb> nb = parse_offsets(offsets, false);
    const idx_t nz = shape[0], ny = shape[1], nx = shape[2];
    if (nz * ny * nx == 0) return 0;
    const uint8_t* A = reinterpret_cast<const uint8_t*>(a.data());
    double* O = out.data();
    nb::gil_scoped_release release;
    const size_t n = (size_t)(nz * ny * nx);
    std::vector<uint8_t> M(n);
    transpose_planes(A, M.data(), nz, ny, nx, num_threads, [](const uint8_t v) { return (uint8_t)(v != 0); });
    std::vector<I> P(n);
    std::vector<I> fin;
    const int64_t num = bwlabeln_scan<I>(M.data(), P.data(), ny, nx, nz, nb, fin);
    untranspose_labels<double, I>(P.data(), O, fin, nz, ny, nx, num_threads);
    return num;
}

// the C flat indices of the voxels of each label in MATLAB's order, grouped by label: label k = 1
// .. num at idx[starts[k - 1] .. starts[k])
template <typename I>
nb::tuple group_by_label(nb::ndarray<const I, nb::ndim<3>, nb::c_contig, nb::device::cpu> lab, int64_t num) {
    const idx_t nz = (idx_t)lab.shape(0), ny = (idx_t)lab.shape(1), nx = (idx_t)lab.shape(2);
    const idx_t n = nz * ny * nx;
    if (num < 0) throw std::invalid_argument("group_by_label: negative num");
    const I* L = lab.data();
    for (idx_t i = 0; i < n; ++i) {
        if (L[i] < 0 || (int64_t)L[i] > num) throw std::invalid_argument("group_by_label: label out of range");
    }
    int64_t* starts = new int64_t[(size_t)num + 1];
    std::fill(starts, starts + num + 1, 0);
    int64_t total = 0;
    int64_t* idx = nullptr;
    {
        nb::gil_scoped_release release;
        for (idx_t i = 0; i < n; ++i) {
            if (L[i] > 0) starts[L[i]]++;
        }
        for (int64_t k = 1; k <= num; ++k) starts[k] += starts[k - 1];
        total = starts[num];
        idx = new int64_t[(size_t)std::max<int64_t>(total, 1)];
        std::vector<int64_t> pos(starts, starts + num);
        for (idx_t z = 0; z < nz; z++) {
            for (idx_t x = 0; x < nx; x++) {
                for (idx_t y = 0; y < ny; y++) {
                    const idx_t i = (z * ny + y) * nx + x;
                    const I l = L[i];
                    if (l > 0) idx[pos[(size_t)(l - 1)]++] = (int64_t)i;
                }
            }
        }
    }
    nb::capsule own_s(starts, [](void* p) noexcept { delete[] static_cast<int64_t*>(p); });
    nb::capsule own_i(idx, [](void* p) noexcept { delete[] static_cast<int64_t*>(p); });
    auto a_idx = nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(idx, {(size_t)total}, own_i);
    auto a_st = nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(starts, {(size_t)num + 1}, own_s);
    return nb::make_tuple(a_idx, a_st);
}

}  // namespace

PK5D_BINDINGS(bind_label) {
    const char* first =
        "Label the components of a C-contiguous (z, y, x) bool volume into lab (0 background), "
        "numbered by first voxel in MATLAB's linear order (y fastest, then x, then z); offsets: "
        "rows (dz, dy, dx) of the connectivity's neighbors before the center in C order. Returns "
        "the number of components.";
    const char* any = "As label_first_*, numbered in C order (for sizes and holes).";
    const char* bwl =
        "bwlabeln's labels (double) of a C-contiguous (z, y, x) bool volume into out; offsets: rows "
        "(dy, dx, dz) of MATLAB's connectivity matrix before its center in MATLAB's order, in the "
        "matrix's column-major order. Returns the number of components.";
    const char* grp =
        "(idx, starts): the C flat indices of the voxels of labels 1 .. num in MATLAB's linear "
        "order, label k at idx[starts[k - 1]:starts[k]].";
    m.def("label_first_int32", &label_first_bind<int32_t>, nb::arg("a").noconvert(), nb::arg("lab").noconvert(),
          nb::arg("offsets").noconvert(), nb::arg("num_threads"), first);
    m.def("label_first_int64", &label_first_bind<int64_t>, nb::arg("a").noconvert(), nb::arg("lab").noconvert(),
          nb::arg("offsets").noconvert(), nb::arg("num_threads"), first);
    m.def("label_any_int32", &label_any_bind<int32_t>, nb::arg("a").noconvert(), nb::arg("lab").noconvert(),
          nb::arg("offsets").noconvert(), nb::arg("num_threads"), any);
    m.def("label_any_int64", &label_any_bind<int64_t>, nb::arg("a").noconvert(), nb::arg("lab").noconvert(),
          nb::arg("offsets").noconvert(), nb::arg("num_threads"), any);
    m.def("bwlabeln_int32", &bwlabeln_bind<int32_t>, nb::arg("a").noconvert(), nb::arg("out").noconvert(),
          nb::arg("offsets").noconvert(), nb::arg("num_threads"), bwl);
    m.def("bwlabeln_int64", &bwlabeln_bind<int64_t>, nb::arg("a").noconvert(), nb::arg("out").noconvert(),
          nb::arg("offsets").noconvert(), nb::arg("num_threads"), bwl);
    m.def("group_by_label_int32", &group_by_label<int32_t>, nb::arg("lab").noconvert(), nb::arg("num"), grp);
    m.def("group_by_label_int64", &group_by_label<int64_t>, nb::arg("lab").noconvert(), nb::arg("num"), grp);
}
