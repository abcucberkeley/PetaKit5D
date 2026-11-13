// Steps of the OMW back projector generation (omw_backprojector_mex's 'window', 'components' and
// 'sortdesc'), with the mex's results bit for bit.
//
// window: for the voxels of a psz grid with center (psz + 1) / 2, decides where the apodization window of
//   omw_backprojector_generation is 1 or 0 whatever the nearest edge points are: the window of a voxel at
//   radius r is (r / rn < l) + taper (l <= r / rn < u), where rn is a convex combination of the radii R of
//   the 2 edge points nearest to the voxel's direction (E: their unit vectors). The unit vectors are binned
//   in cells of 1/64 of [-1, 1]^3; for each cell, the 2 nearest edge points of any direction in the cell are
//   within the distance of the 2nd nearest edge point from the cell center plus twice the cell's half
//   diagonal, so the min and max of R over the edge points within that distance bound rn. mask is 1 for
//   voxels with r < l * min R (with a margin of 1e-9 for rounding), 0 elsewhere; inds are the MATLAB linear
//   indices of the voxels in between (r < u * max R), whose window has to be computed. The vectors keep the
//   mex's component order (y, x, z), so every distance is the mex's.
// components: the 26-connected components of a 3D mask in MATLAB's layout that bwareaopen(A, n) keeps
//   ('minsize': at least n voxels) or that are the n largest (ties to the component first in bwconncomp's
//   order, that of their first voxels).
// sortdesc: sort(double(a(:)), 'descend') for single values without NaN or negative values (radix sort;
//   the sorted values don't depend on the algorithm).
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

#include "decon.hpp"

namespace pk5d::decon {
namespace {

// --- window bounds ---------------------------------------------------------------------------------

constexpr int G = 128;            // cells per axis of [-1, 1]
constexpr double CELL = 2.0 / G;  // cell size
constexpr int MAXSHELL = 8;       // cells searched around a cell before all edge points are

inline int cell_of(double v) {
    const int c = static_cast<int>(std::floor((v + 1.0) * (G / 2.0)));
    return c < 0 ? 0 : (c >= G ? G - 1 : c);
}

inline idx_t cell_index(int a, int b, int c) { return (static_cast<idx_t>(c) * G + b) * G + a; }

// calls f(cell) for the cells at Chebyshev distance sh from cell (a, b, c)
template <class F>
void visit_shell(int a, int b, int c, int sh, F f) {
    for (int kc = std::max(0, c - sh); kc <= std::min(G - 1, c + sh); kc++) {
        const bool zFace = std::abs(kc - c) == sh;
        for (int kb = std::max(0, b - sh); kb <= std::min(G - 1, b + sh); kb++) {
            if (zFace || std::abs(kb - b) == sh) {
                for (int ka = std::max(0, a - sh); ka <= std::min(G - 1, a + sh); ka++) f(cell_index(ka, kb, kc));
            } else {
                if (a - sh >= 0) f(cell_index(a - sh, kb, kc));
                if (a + sh <= G - 1) f(cell_index(a + sh, kb, kc));
            }
        }
    }
}

// --- connected components (26-connected), MATLAB's layout ---------------------------------------------

inline int32_t find_root(int32_t* parent, int32_t i) {
    while (parent[i] != i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}

// union by the smaller index, so that the root of a component is its first voxel in linear order (and a
// parent is always before its child)
inline void unite(int32_t* parent, int32_t a, int32_t b) {
    const int32_t ra = find_root(parent, a), rb = find_root(parent, b);
    if (ra < rb) {
        parent[rb] = ra;
    } else if (rb < ra) {
        parent[ra] = rb;
    }
}

// labels of the components of A in lab (1-based, in the order of their first voxels; 0: background) and
// their sizes (sizes[label])
void label_components(const uint8_t* A, idx_t n1, idx_t n2, idx_t n3, int32_t* lab, std::vector<idx_t>& sizes,
                      int nThreads) {
    const idx_t np = n1 * n2;
    // offsets (x, y, z; x along dimension 1) of the 13 neighbors before a voxel in linear order
    int nb[13][3];
    int m = 0;
    for (int dz = -1; dz <= 0; dz++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                if (dz < 0 || (dz == 0 && dy < 0) || (dz == 0 && dy == 0 && dx < 0)) {
                    nb[m][0] = dx;
                    nb[m][1] = dy;
                    nb[m][2] = dz;
                    m++;
                }
            }
        }
    }
    // components within slabs of planes in parallel, then joined across the slabs
    const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(nThreads, n3)));
    std::vector<idx_t> zs(static_cast<size_t>(nt + 1));
    for (int t = 0; t <= nt; t++) zs[t] = n3 * t / nt;
    int32_t* parent = lab;
    auto link = [&](idx_t z, idx_t y, idx_t x, idx_t zmin) {
        const int32_t i = static_cast<int32_t>(z * np + y * n1 + x);
        for (int k = 0; k < 13; k++) {
            const idx_t xx = x + nb[k][0], yy = y + nb[k][1], zz = z + nb[k][2];
            if (xx < 0 || xx >= n1 || yy < 0 || yy >= n2 || zz < zmin) continue;
            const int32_t j = static_cast<int32_t>(zz * np + yy * n1 + xx);
            if (A[j]) unite(parent, i, j);
        }
    };
#pragma omp parallel for num_threads(nt) schedule(static, 1)
    for (int t = 0; t < nt; t++) {
        for (idx_t z = zs[t]; z < zs[t + 1]; z++) {
            for (idx_t y = 0; y < n2; y++) {
                for (idx_t x = 0; x < n1; x++) {
                    const idx_t i = z * np + y * n1 + x;
                    parent[i] = A[i] ? static_cast<int32_t>(i) : -1;
                    if (A[i]) link(z, y, x, zs[t]);
                }
            }
        }
    }
    for (int t = 1; t < nt; t++) {
        const idx_t z = zs[t];
        for (idx_t y = 0; y < n2; y++) {
            for (idx_t x = 0; x < n1; x++) {
                if (A[z * np + y * n1 + x]) link(z, y, x, z - 1);
            }
        }
    }
    // in linear order a voxel's parent is final (its root) before the voxel, so one step finds the root;
    // a root gets the next label, and the other voxels the label already in their root's entry
    sizes.assign(1, 0);
    const idx_t n = np * n3;
    int32_t nl = 0;
    for (idx_t i = 0; i < n; i++) {
        const int32_t p = parent[i];
        if (p < 0) {
            lab[i] = 0;
            continue;
        }
        int32_t l;
        if (p == static_cast<int32_t>(i)) {
            l = ++nl;
            sizes.push_back(0);
        } else {
            // the parent was visited: its entry holds -(label) of the root
            l = -lab[p];
        }
        sizes[l]++;
        lab[i] = -l;
    }
#pragma omp parallel for num_threads(nThreads) schedule(static)
    for (idx_t i = 0; i < n; i++) lab[i] = -lab[i];
}

}  // namespace

void omw_window(const idx_t* psz, const double* E, const double* R, idx_t nE, double l, double u, double* mask,
                std::vector<int64_t>& uinds, int nThreads) {
    const idx_t nCells = static_cast<idx_t>(G) * G * G;
    // edge points by cell (CSR)
    std::vector<int32_t> start(static_cast<size_t>(nCells + 1), 0);
    std::vector<idx_t> ecell(static_cast<size_t>(nE));
    for (idx_t j = 0; j < nE; j++) {
        ecell[j] = cell_index(cell_of(E[j]), cell_of(E[j + nE]), cell_of(E[j + 2 * nE]));
        start[ecell[j] + 1]++;
    }
    for (idx_t i = 0; i < nCells; i++) start[i + 1] += start[i];
    std::vector<int32_t> fill(start.begin(), start.end() - 1);
    std::vector<int32_t> pts(static_cast<size_t>(nE));
    for (idx_t j = 0; j < nE; j++) pts[fill[ecell[j]]++] = static_cast<int32_t>(j);

    // cells that meet the unit sphere (the directions of the voxels are unit vectors)
    std::vector<int32_t> slot(static_cast<size_t>(nCells), -1);
    std::vector<idx_t> sphereCells;
    for (int c = 0; c < G; c++) {
        for (int b = 0; b < G; b++) {
            for (int a = 0; a < G; a++) {
                const int ijk[3] = {a, b, c};
                double dmin = 0, dmax = 0;
                for (int t = 0; t < 3; t++) {
                    const double lo = -1.0 + ijk[t] * CELL, hi = lo + CELL;
                    const double nearv = (lo > 0) ? lo : ((hi < 0) ? hi : 0.0);
                    const double farv = std::max(std::fabs(lo), std::fabs(hi));
                    dmin += nearv * nearv;
                    dmax += farv * farv;
                }
                if (dmin <= (1 + 1e-6) * (1 + 1e-6) && dmax >= (1 - 1e-6) * (1 - 1e-6)) {
                    slot[cell_index(a, b, c)] = static_cast<int32_t>(sphereCells.size());
                    sphereCells.push_back(cell_index(a, b, c));
                }
            }
        }
    }
    const idx_t nS = static_cast<idx_t>(sphereCells.size());
    std::vector<double> Rlo(static_cast<size_t>(nS)), Rhi(static_cast<size_t>(nS));
    // half diagonal of a cell, with a margin for the directions of the voxels as MATLAB computes them
    const double rho = std::sqrt(3.0) * (CELL / 2 + 1e-9);

#pragma omp parallel for num_threads(nThreads) schedule(dynamic, 64)
    for (idx_t s = 0; s < nS; s++) {
        const idx_t ci = sphereCells[s];
        const int a = static_cast<int>(ci % G), b = static_cast<int>((ci / G) % G), c = static_cast<int>(ci / (static_cast<idx_t>(G) * G));
        const double ctr[3] = {-1.0 + (a + 0.5) * CELL, -1.0 + (b + 0.5) * CELL, -1.0 + (c + 0.5) * CELL};
        // distance to the 2nd nearest edge point, over shells of cells (points in shells beyond sh are
        // farther than (sh + 0.5) cells), or over all edge points if they are far
        double d1 = INFINITY, d2 = INFINITY;
        auto nearest2 = [&](idx_t j) {
            const double e0 = E[j] - ctr[0], e1 = E[j + nE] - ctr[1], e2 = E[j + 2 * nE] - ctr[2];
            const double d = std::sqrt(e0 * e0 + e1 * e1 + e2 * e2);
            if (d < d1) {
                d2 = d1;
                d1 = d;
            } else if (d < d2) {
                d2 = d;
            }
        };
        bool found = false;
        for (int sh = 0; sh <= MAXSHELL && !found; sh++) {
            visit_shell(a, b, c, sh, [&](idx_t k) {
                for (int32_t q = start[k]; q < start[k + 1]; q++) nearest2(pts[q]);
            });
            found = d2 <= (sh + 0.5) * CELL;
        }
        if (!found) {
            d1 = d2 = INFINITY;
            for (idx_t j = 0; j < nE; j++) nearest2(j);
        }
        // the 2 nearest edge points of any direction in the cell are within rad of the cell center
        const double rad = d2 + 2 * rho + 1e-12;
        const int S = static_cast<int>(std::ceil(rad / CELL + 0.5));
        double lo = INFINITY, hi = -INFINITY;
        auto bound = [&](idx_t j) {
            const double e0 = E[j] - ctr[0], e1 = E[j + nE] - ctr[1], e2 = E[j + 2 * nE] - ctr[2];
            if (std::sqrt(e0 * e0 + e1 * e1 + e2 * e2) <= rad) {
                lo = std::min(lo, R[j]);
                hi = std::max(hi, R[j]);
            }
        };
        if (S <= MAXSHELL) {
            for (int sh = 0; sh <= S; sh++) {
                visit_shell(a, b, c, sh, [&](idx_t k) {
                    for (int32_t q = start[k]; q < start[k + 1]; q++) bound(pts[q]);
                });
            }
        } else {
            for (idx_t j = 0; j < nE; j++) bound(j);
        }
        Rlo[s] = lo;
        Rhi[s] = hi;
    }

    // voxels: 1 where the window is 1, 0 where it is 0, 2 where it has to be computed; first with the
    // bounds of all edge radii (on the squared radius, with a margin for the square root). The mask is
    // written in the port's layout (plane k: rows i, columns j), the indices in MATLAB's order.
    const idx_t n1 = psz[0], n2 = psz[1], n3 = psz[2];
    const double cy = (n1 + 1) / 2.0, cx = (n2 + 1) / 2.0, cz = (n3 + 1) / 2.0;
    double Rmin = INFINITY, Rmax = -INFINITY;
    for (idx_t j = 0; j < nE; j++) {
        Rmin = std::min(Rmin, R[j]);
        Rmax = std::max(Rmax, R[j]);
    }
    const double in2 = (l * Rmin * (1 - 1e-9)) * (l * Rmin * (1 - 1e-9)) * (1 - 1e-12);
    const double out2 = (u * Rmax * (1 + 1e-9)) * (u * Rmax * (1 + 1e-9)) * (1 + 1e-12);
    std::vector<std::vector<int64_t>> planeInds(static_cast<size_t>(n3));
#pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1)
    for (idx_t k = 0; k < n3; k++) {
        const double Z = (k + 1) - cz;
        std::vector<uint8_t> cls(static_cast<size_t>(n1 * n2));
        double* mk = mask + k * n1 * n2;
        for (idx_t j = 0; j < n2; j++) {
            const double X = (j + 1) - cx;
            for (idx_t i = 0; i < n1; i++) {
                const double Y = (i + 1) - cy;
                const double r2 = Y * Y + X * X + Z * Z;
                uint8_t v = 2;
                if (r2 > 0 && r2 < in2) {
                    v = 1;
                } else if (r2 >= out2) {
                    v = 0;
                } else if (r2 > 0) {
                    const double r = std::sqrt(r2);
                    const double inv = 1.0 / r;
                    const int32_t sl = slot[cell_index(cell_of(Y * inv), cell_of(X * inv), cell_of(Z * inv))];
                    if (sl >= 0 && Rlo[sl] <= Rhi[sl]) {
                        if (r < l * Rlo[sl] * (1 - 1e-9)) {
                            v = 1;
                        } else if (r >= u * Rhi[sl] * (1 + 1e-9)) {
                            v = 0;
                        }
                    }
                }
                cls[j * n1 + i] = v;
                mk[i * n2 + j] = (v == 1) ? 1.0 : 0.0;
            }
        }
        std::vector<int64_t>& pi = planeInds[k];
        const idx_t base = k * n1 * n2;
        for (idx_t i = 0; i < n1 * n2; i++) {
            if (cls[i] == 2) pi.push_back(static_cast<int64_t>(base + i));
        }
    }
    uinds.clear();
    for (const auto& pi : planeInds) uinds.insert(uinds.end(), pi.begin(), pi.end());
}

void omw_components(const uint8_t* A, idx_t n1, idx_t n2, idx_t n3, bool largest, idx_t n, uint8_t* B,
                    int nThreads) {
    const idx_t nv = n1 * n2 * n3;
    std::vector<int32_t> lab(static_cast<size_t>(nv > 0 ? nv : 1));
    std::vector<idx_t> sizes;
    label_components(A, n1, n2, n3, lab.data(), sizes, nThreads);
    const idx_t nl = static_cast<idx_t>(sizes.size()) - 1;
    std::vector<uint8_t> keep(static_cast<size_t>(nl + 1), 0);
    if (largest) {
        // the n largest, ties to the earlier component
        std::vector<idx_t> order(static_cast<size_t>(nl));
        for (idx_t k = 0; k < nl; k++) order[k] = k + 1;
        std::stable_sort(order.begin(), order.end(), [&](idx_t a, idx_t b) { return sizes[a] > sizes[b]; });
        for (idx_t k = 0; k < std::min(n, nl); k++) keep[order[k]] = 1;
    } else {
        for (idx_t k = 1; k <= nl; k++) keep[k] = sizes[k] >= n;
    }
#pragma omp parallel for num_threads(nThreads) schedule(static)
    for (idx_t i = 0; i < nv; i++) B[i] = keep[lab[i]];
}

bool omw_sortdesc(const float* a, idx_t n, double* out, int nThreads) {
    // nonnegative floats sort as their bits
    std::vector<uint32_t> k0(static_cast<size_t>(n)), k1(static_cast<size_t>(n));
    int bad = 0;
#pragma omp parallel for num_threads(nThreads) schedule(static) reduction(| : bad)
    for (idx_t i = 0; i < n; i++) {
        uint32_t b;
        std::memcpy(&b, a + i, 4);
        // NaN and negative values (-0 included) are left to MATLAB's sort
        bad |= (b > 0x7F800000u);
        k0[i] = b;
    }
    if (bad) return false;
    // LSD radix sort, 4 passes of 8 bits, ascending
    const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(nThreads, n / 65536 + 1)));
    std::vector<idx_t> hist(static_cast<size_t>(nt) * 256);
    uint32_t* src = k0.data();
    uint32_t* dst = k1.data();
    for (int pass = 0; pass < 4; pass++) {
        const int sh = 8 * pass;
        std::fill(hist.begin(), hist.end(), 0);
#pragma omp parallel num_threads(nt)
        {
            const int t = omp_get_thread_num();
            const idx_t i0 = n * t / nt, i1 = n * (t + 1) / nt;
            idx_t* h = hist.data() + static_cast<size_t>(t) * 256;
            for (idx_t i = i0; i < i1; i++) h[(src[i] >> sh) & 255]++;
#pragma omp barrier
#pragma omp single
            {
                idx_t sum = 0;
                for (int b = 0; b < 256; b++) {
                    for (int tt = 0; tt < nt; tt++) {
                        const idx_t c = hist[static_cast<size_t>(tt) * 256 + b];
                        hist[static_cast<size_t>(tt) * 256 + b] = sum;
                        sum += c;
                    }
                }
            }
            for (idx_t i = i0; i < i1; i++) dst[h[(src[i] >> sh) & 255]++] = src[i];
        }
        std::swap(src, dst);
    }
#pragma omp parallel for num_threads(nThreads) schedule(static)
    for (idx_t i = 0; i < n; i++) {
        float f;
        std::memcpy(&f, src + (n - 1 - i), 4);
        out[i] = static_cast<double>(f);
    }
    return true;
}

}  // namespace pk5d::decon
