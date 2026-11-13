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
// window values (R02, omw_backprojector_generation's omw_window_with_mex): window's mask with the
//   window of the voxels it leaves computed as the MATLAB code computes it: the directions and radii by
//   MATLAB's cart2sph and sph2cart (fdlibm and msun in C), the 2 nearest edge points by the squared
//   distance summed over y, x, z (ties to the lower index; found on a grid of cells), rn weighted with
//   the distances of the first nE voxels (MATLAB's indexing), the cos^2 taper.
// abs (R02): MATLAB's abs of complex single values, msun's hypotf (measured), MATLAB's NaN for a NaN part.
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

// --- the window values and |OTF| of omw_backprojector_generation (R02) ------------------------------
//
// MATLAB's elementary functions as petakit5d.matlabcompat.elementary writes them out: hypot is FreeBSD
// msun's e_hypot.c, atan2, atan, sin and cos fdlibm 5.3's (e_atan2.c, s_atan.c, e_rem_pio2.c, k_sin.c,
// k_cos.c), single hypotf msun's e_hypotf.c, with MATLAB's rules for NaN and two zeros; one IEEE
// operation per C operation (the module is built without fused multiply-adds). Their notice:
//   Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//   Developed at SunSoft, a Sun Microsystems, Inc. business.
//   Permission to use, copy, modify, and distribute this software is freely granted, provided
//   that this notice is preserved.
namespace {

inline int32_t hi_of(double x) {
    uint64_t b;
    std::memcpy(&b, &x, 8);
    return static_cast<int32_t>(b >> 32);
}
inline uint32_t lo_of(double x) {
    uint64_t b;
    std::memcpy(&b, &x, 8);
    return static_cast<uint32_t>(b & 0xFFFFFFFFu);
}
inline double from_hi(int32_t hw) {
    const uint64_t b = static_cast<uint64_t>(static_cast<uint32_t>(hw)) << 32;
    double x;
    std::memcpy(&x, &b, 8);
    return x;
}
inline double with_hi(double x, int32_t hw) {
    const uint64_t b = (static_cast<uint64_t>(static_cast<uint32_t>(hw)) << 32) | lo_of(x);
    double y;
    std::memcpy(&y, &b, 8);
    return y;
}
inline double from_bits(uint64_t b) {
    double x;
    std::memcpy(&x, &b, 8);
    return x;
}

const double kPosNaN = from_bits(0x7FF8000000000000ull);
const double kMatlabNaN = from_bits(0xFFF8000000000000ull);

// msun hypot, +Inf for an infinite argument, the positive NaN for a NaN one (MATLAB's)
double m_hypot(double x, double y) {
    if (std::isnan(x) || std::isnan(y)) return kPosNaN;
    if (std::isinf(x) || std::isinf(y)) return INFINITY;
    int32_t ha = hi_of(x) & 0x7FFFFFFF, hb = hi_of(y) & 0x7FFFFFFF;
    double a = x, b = y;
    if (hb > ha) {
        a = y;
        b = x;
        std::swap(ha, hb);
    }
    a = std::fabs(a);
    b = std::fabs(b);
    if ((ha - hb) > 0x3C00000) return a + b;
    int k = 0;
    if (ha > 0x5F300000) {
        ha -= 0x25800000;
        hb -= 0x25800000;
        k += 600;
        a = with_hi(a, ha);
        b = with_hi(b, hb);
    }
    if (hb < 0x23D00000) {
        if (hb <= 0x000FFFFF) {
            if ((hb | static_cast<int32_t>(lo_of(b))) == 0) return a;
            const double t1 = from_hi(0x7FD00000);  // 2^1022
            b *= t1;
            a *= t1;
            k -= 1022;
            ha = hi_of(a);
            hb = hi_of(b);
            if (hb > ha) {
                std::swap(a, b);
                std::swap(ha, hb);
            }
        } else {
            ha += 0x25800000;
            hb += 0x25800000;
            k -= 600;
            a = with_hi(a, ha);
            b = with_hi(b, hb);
        }
    }
    double w = a - b;
    if (w > b) {
        const double t1 = from_hi(ha);
        const double t2 = a - t1;
        w = std::sqrt(t1 * t1 - (b * (-b) - t2 * (a + t1)));
    } else {
        a = a + a;
        const double y1 = from_hi(hb);
        const double y2 = b - y1;
        const double t1 = from_hi(ha + 0x00100000);
        const double t2 = a - t1;
        w = std::sqrt(t1 * y1 - (w * (-w) - (t1 * y2 + t2 * b)));
    }
    if (k != 0) return from_hi((1023 + k) << 20) * w;
    return w;
}

const double kAtanHi[4] = {4.63647609000806093515e-01, 7.85398163397448278999e-01, 9.82793723247329054082e-01,
                           1.57079632679489655800e00};
const double kAtanLo[4] = {2.26987774529616870924e-17, 3.06161699786838301793e-17, 1.39033110312309984516e-17,
                           6.12323399573676603587e-17};
const double kAT[11] = {3.33333333333329318027e-01,  -1.99999999998764832476e-01, 1.42857142725034663711e-01,
                        -1.11111104054623557880e-01, 9.09088713343650656196e-02,  -7.69187620504482999495e-02,
                        6.66107313738753120669e-02,  -5.83357013379057348645e-02, 4.97687799461593236017e-02,
                        -3.65315727442169155270e-02, 1.62858201153657823623e-02};
constexpr double kTiny = 1.0e-300;
constexpr double kPiO4 = 7.8539816339744827900e-01;
constexpr double kPiO2 = 1.5707963267948965580e00;
constexpr double kPi = 3.1415926535897931160e00;
constexpr double kPiLo = 1.2246467991473531772e-16;

double fd_atan(double x) {
    const int32_t hx = hi_of(x), ix = hx & 0x7FFFFFFF;
    if (std::isnan(x)) return x + x;
    if (ix >= 0x44100000) return hx > 0 ? kAtanHi[3] + kAtanLo[3] : -kAtanHi[3] - kAtanLo[3];
    int id;
    if (ix < 0x3FDC0000) {
        if (ix < 0x3E200000) return x;
        id = -1;
    } else {
        x = std::fabs(x);
        if (ix < 0x3FF30000) {
            if (ix < 0x3FE60000) {
                id = 0;
                x = (2.0 * x - 1.0) / (2.0 + x);
            } else {
                id = 1;
                x = (x - 1.0) / (x + 1.0);
            }
        } else if (ix < 0x40038000) {
            id = 2;
            x = (x - 1.5) / (1.0 + 1.5 * x);
        } else {
            id = 3;
            x = -1.0 / x;
        }
    }
    const double z = x * x, w = z * z;
    const double s1 = z * (kAT[0] + w * (kAT[2] + w * (kAT[4] + w * (kAT[6] + w * (kAT[8] + w * kAT[10])))));
    const double s2 = w * (kAT[1] + w * (kAT[3] + w * (kAT[5] + w * (kAT[7] + w * kAT[9]))));
    if (id < 0) return x - x * (s1 + s2);
    const double zz = kAtanHi[id] - ((x * (s1 + s2) - kAtanLo[id]) - x);
    return hx < 0 ? -zz : zz;
}

// fdlibm atan2, +0 for two zeros and MATLAB's NaN for a NaN (MATLAB's)
double m_atan2(double y, double x) {
    if (std::isnan(x) || std::isnan(y)) return kMatlabNaN;
    const int32_t hx = hi_of(x), hy = hi_of(y), ix = hx & 0x7FFFFFFF, iy = hy & 0x7FFFFFFF;
    const uint32_t lx = lo_of(x), ly = lo_of(y);
    const bool xzero = (ix | static_cast<int32_t>(lx)) == 0, yzero = (iy | static_cast<int32_t>(ly)) == 0;
    if (xzero && yzero) return 0.0;
    if (x == 1.0) return fd_atan(y);
    const int m = ((hy >> 31) & 1) | ((hx >> 30) & 2);
    if (yzero) {
        if (m < 2) return y;
        return m == 2 ? kPi + kTiny : -kPi - kTiny;
    }
    if (xzero) return hy < 0 ? -kPiO2 - kTiny : kPiO2 + kTiny;
    if (ix == 0x7FF00000) {
        if (iy == 0x7FF00000) {
            const double v[4] = {kPiO4 + kTiny, -kPiO4 - kTiny, 3.0 * kPiO4 + kTiny, -3.0 * kPiO4 - kTiny};
            return v[m];
        }
        const double v[4] = {0.0, -0.0, kPi + kTiny, -kPi - kTiny};
        return v[m];
    }
    if (iy == 0x7FF00000) return hy < 0 ? -kPiO2 - kTiny : kPiO2 + kTiny;
    const int k = (iy - ix) >> 20;
    double z;
    if (k > 60) {
        z = kPiO2 + 0.5 * kPiLo;
    } else if (hx < 0 && k < -60) {
        z = 0.0;
    } else {
        z = fd_atan(std::fabs(y / x));
    }
    switch (m) {
        case 0:
            return z;
        case 1:
            return -z;
        case 2:
            return kPi - (z - kPiLo);
        default:
            return (z - kPiLo) - kPi;
    }
}

const double kS[6] = {-1.66666666666666324348e-01, 8.33333333332248946124e-03, -1.98412698298579493134e-04,
                      2.75573137070700676789e-06,  -2.50507602534068634195e-08, 1.58969099521155010221e-10};
const double kC[6] = {4.16666666666666019037e-02,  -1.38888888888741095749e-03, 2.48015872894767294178e-05,
                      -2.75573143513906633035e-07, 2.08757232129817482790e-09,  -1.13596475577881948265e-11};

double k_sin(double x, double y, int iy) {
    if ((hi_of(x) & 0x7FFFFFFF) < 0x3E400000) return x;
    const double z = x * x, v = z * x;
    const double r = kS[1] + z * (kS[2] + z * (kS[3] + z * (kS[4] + z * kS[5])));
    if (iy == 0) return x + v * (kS[0] + z * r);
    return x - ((z * (0.5 * y - v * r) - y) - v * kS[0]);
}

double k_cos(double x, double y) {
    const int32_t ix = hi_of(x) & 0x7FFFFFFF;
    if (ix < 0x3E400000) return 1.0;
    const double z = x * x;
    const double r = z * (kC[0] + z * (kC[1] + z * (kC[2] + z * (kC[3] + z * (kC[4] + z * kC[5])))));
    if (ix < 0x3FD33333) return 1.0 - (0.5 * z - (z * r - x * y));
    const double qx = ix > 0x3FE90000 ? 0.28125 : from_hi(ix - 0x00200000);
    const double hz = 0.5 * z - qx, a = 1.0 - qx;
    return a - (hz - (z * r - x * y));
}

constexpr double kInvPio2 = 6.36619772367581382433e-01;
constexpr double kPio2_1 = 1.57079632673412561417e00, kPio2_1t = 6.07710050650619224932e-11;
constexpr double kPio2_2 = 6.07710050630396597660e-11, kPio2_2t = 2.02226624879595063154e-21;
constexpr double kPio2_3 = 2.02226624871116645580e-21, kPio2_3t = 8.47842766036889956997e-32;
const int32_t kNpio2Hw[32] = {0x3FF921FB, 0x400921FB, 0x4012D97C, 0x401921FB, 0x401F6A7A, 0x4022D97C,
                              0x4025FDBB, 0x402921FB, 0x402C463A, 0x402F6A7A, 0x4031475C, 0x4032D97C,
                              0x40346B9C, 0x4035FDBB, 0x40378FDB, 0x403921FB, 0x403AB41B, 0x403C463A,
                              0x403DD85A, 0x403F6A7A, 0x40407E4C, 0x4041475C, 0x4042106C, 0x4042D97C,
                              0x4043A28C, 0x40446B9C, 0x404534AC, 0x4045FDBB, 0x4046C6CB, 0x40478FDB,
                              0x404858EB, 0x404921FB};

// fdlibm __ieee754_rem_pio2 for pi/4 < |x| <= 2^20 pi/2
int rem_pio2(double x, double* y) {
    const int32_t hx = hi_of(x), ix = hx & 0x7FFFFFFF;
    if (ix < 0x4002D97C) {
        const double sg = hx > 0 ? 1.0 : -1.0;
        double z = x - sg * kPio2_1;
        if (ix != 0x3FF921FB) {
            y[0] = z - sg * kPio2_1t;
            y[1] = (z - y[0]) - sg * kPio2_1t;
        } else {
            z = z - sg * kPio2_2;
            y[0] = z - sg * kPio2_2t;
            y[1] = (z - y[0]) - sg * kPio2_2t;
        }
        return hx > 0 ? 1 : -1;
    }
    const double t = std::fabs(x);
    const int n = static_cast<int>(t * kInvPio2 + 0.5);
    const double fn = static_cast<double>(n);
    double r = t - fn * kPio2_1, w = fn * kPio2_1t;
    if (n < 32 && ix != kNpio2Hw[n - 1]) {
        y[0] = r - w;
    } else {
        const int32_t j = ix >> 20;
        y[0] = r - w;
        if (j - ((hi_of(y[0]) >> 20) & 0x7FF) > 16) {
            double tt = r;
            w = fn * kPio2_2;
            r = tt - w;
            w = fn * kPio2_2t - ((tt - r) - w);
            y[0] = r - w;
            if (j - ((hi_of(y[0]) >> 20) & 0x7FF) > 49) {
                tt = r;
                w = fn * kPio2_3;
                r = tt - w;
                w = fn * kPio2_3t - ((tt - r) - w);
                y[0] = r - w;
            }
        }
    }
    y[1] = (r - y[0]) - w;
    if (hx < 0) {
        y[0] = -y[0];
        y[1] = -y[1];
        return -n;
    }
    return n;
}

// fdlibm sin and cos (x - x for Inf and NaN; above 2^20 pi/2 the C library's, as elementary)
double fd_sincos(double x, bool sine) {
    const int32_t ix = hi_of(x) & 0x7FFFFFFF;
    if (ix <= 0x3FE921FB) return sine ? k_sin(x, 0.0, 0) : k_cos(x, 0.0);
    if (ix >= 0x7FF00000) return x - x;
    if (ix > 0x413921FB) return sine ? std::sin(x) : std::cos(x);
    double y[2];
    const int q = rem_pio2(x, y) & 3;
    const double s = k_sin(y[0], y[1], 1), c = k_cos(y[0], y[1]);
    if (sine) return q == 0 ? s : (q == 1 ? c : (q == 2 ? -s : -c));
    return q == 0 ? c : (q == 1 ? -s : (q == 2 ? -c : s));
}

inline float f_of(int32_t w) {
    float f;
    std::memcpy(&f, &w, 4);
    return f;
}

// msun hypotf of |parts| with the high words ha >= hb (signs cleared), MATLAB's NaN for a NaN part
float m_hypotf(int32_t ha, int32_t hb) {
    constexpr int32_t HI = static_cast<int32_t>(0xFFFFF000u);
    if ((ha - hb) > 0x0F000000) return f_of(ha) + f_of(hb);
    if (ha >= 0x7F800000) return ha > 0x7F800000 ? f_of(0x7FC00000) : INFINITY;
    float a = f_of(ha), b = f_of(hb);
    int k = 0;
    if (ha > 0x58800000) {
        ha -= 0x22000000;
        hb -= 0x22000000;
        k = 68;
        a = f_of(ha);
        b = f_of(hb);
    }
    if (hb < 0x26800000) {
        if (hb == 0) return f_of(ha);
        if (hb <= 0x007FFFFF) {
            // the high words are not read again after this scaling
            a *= 0x1p126f;
            b *= 0x1p126f;
            k = -126;
        } else {
            ha += 0x22000000;
            hb += 0x22000000;
            k = -68;
            a = f_of(ha);
            b = f_of(hb);
        }
    }
    float w = a - b, h;
    if (w > b) {
        const float t1 = f_of(ha & HI);
        h = std::sqrt(t1 * t1 - (b * (-b) - (a - t1) * (a + t1)));
    } else {
        const float y1 = f_of(hb & HI), t1 = f_of((ha + 0x00800000) & HI);
        h = std::sqrt(t1 * y1 - (w * (-w) - (t1 * (b - y1) + ((a + a) - t1) * b)));
    }
    if (k != 0) h = f_of(0x3F800000 + (k << 23)) * h;
    return h;
}

// the unit direction (y, x, z) and the radius of the offset (dx, dy, dz): MATLAB's cart2sph and
// sph2cart with radius 1, given the plane's hypot(dx, dy) and cos and sin of atan2(dy, dx)
inline void direction(double h, double caz, double saz, double dz, double* v, double* r) {
    const double el = m_atan2(dz, h);
    *r = m_hypot(h, dz);
    const double ce = fd_sincos(el, false);
    v[0] = ce * saz;
    v[1] = ce * caz;
    v[2] = fd_sincos(el, true);
}

}  // namespace

void omw_abs(const float* c, idx_t n, float* out, int nThreads) {
#pragma omp parallel for num_threads(nThreads) schedule(static)
    for (idx_t i = 0; i < n; i++) {
        int32_t hx, hy;
        std::memcpy(&hx, c + 2 * i, 4);
        std::memcpy(&hy, c + 2 * i + 1, 4);
        hx &= 0x7FFFFFFF;
        hy &= 0x7FFFFFFF;
        out[i] = m_hypotf(std::max(hx, hy), std::min(hx, hy));
    }
}

bool omw_window_values(const idx_t* psz, const int64_t* et, idx_t nE, double l, double u, double* mask,
                       int nThreads) {
    const idx_t n1 = psz[0], n2 = psz[1], n3 = psz[2], np = n1 * n2;
    const double cy = (n1 + 1) / 2.0, cx = (n2 + 1) / 2.0, cz = (n3 + 1) / 2.0;
    // the plane's hypot(dx, dy) and the cosine and sine of atan2(dy, dx), MATLAB's linear order
    std::vector<double> H(static_cast<size_t>(np)), CA(static_cast<size_t>(np)), SA(static_cast<size_t>(np));
#pragma omp parallel for num_threads(nThreads) schedule(static)
    for (idx_t p = 0; p < np; p++) {
        const double dy = static_cast<double>(p % n1 + 1) - cy, dx = static_cast<double>(p / n1 + 1) - cx;
        H[p] = m_hypot(dx, dy);
        const double az = m_atan2(dy, dx);
        CA[p] = fd_sincos(az, false);
        SA[p] = fd_sincos(az, true);
    }
    auto dir = [&](int64_t t, double* v, double* r) {
        const idx_t p = static_cast<idx_t>(t % np);
        direction(H[p], CA[p], SA[p], static_cast<double>(t / np + 1) - cz, v, r);
    };
    // the edge points' directions (columns y, x, z) and radii
    std::vector<double> E(static_cast<size_t>(3 * nE)), R(static_cast<size_t>(nE));
#pragma omp parallel for num_threads(nThreads) schedule(static)
    for (idx_t j = 0; j < nE; j++) {
        double v[3];
        dir(et[j], v, &R[j]);
        E[j] = v[0];
        E[j + nE] = v[1];
        E[j + 2 * nE] = v[2];
    }
    std::vector<int64_t> uinds;
    omw_window(psz, E.data(), R.data(), nE, l, u, mask, uinds, nThreads);

    // the 2 nearest edge points (squared distances summed over y, x, z in order, ties to the lower
    // index) of the first nE voxels and of the voxels to compute, on a grid of cells of about 2 points
    const idx_t nQ = nE + static_cast<idx_t>(uinds.size());
    const int g = static_cast<int>(std::max(8.0, std::min(160.0, std::sqrt(nE / 6.2832))));
    const double cell = 2.0 / g;
    auto cof = [&](double v) {
        const int c = static_cast<int>(std::floor((v + 1.0) * (g / 2.0)));
        return c < 0 ? 0 : (c >= g ? g - 1 : c);
    };
    const idx_t nc = static_cast<idx_t>(g) * g * g;
    std::vector<int32_t> cstart(static_cast<size_t>(nc + 1), 0), cpts(static_cast<size_t>(nE));
    std::vector<idx_t> ec(static_cast<size_t>(nE));
    for (idx_t j = 0; j < nE; j++) {
        ec[j] = (static_cast<idx_t>(cof(E[j + 2 * nE])) * g + cof(E[j + nE])) * g + cof(E[j]);
        cstart[ec[j] + 1]++;
    }
    for (idx_t c = 0; c < nc; c++) cstart[c + 1] += cstart[c];
    {
        std::vector<int32_t> fill(cstart.begin(), cstart.end() - 1);
        for (idx_t j = 0; j < nE; j++) cpts[fill[ec[j]]++] = static_cast<int32_t>(j);
    }
    std::vector<int32_t> nn(static_cast<size_t>(2 * nQ));
    std::vector<double> dd(static_cast<size_t>(2 * nE));  // the distances of the first nE voxels
    std::vector<double> rq(static_cast<size_t>(nQ - nE));
#pragma omp parallel for num_threads(nThreads) schedule(dynamic, 1024)
    for (idx_t q = 0; q < nQ; q++) {
        const int64_t t = q < nE ? static_cast<int64_t>(q) : uinds[static_cast<size_t>(q - nE)];
        double v[3], r;
        dir(t, v, &r);
        if (q >= nE) rq[q - nE] = r;
        double b1 = INFINITY, b2 = INFINITY;
        int32_t i1 = -1, i2 = -1;
        auto consider = [&](int32_t j) {
            const double e0 = v[0] - E[j], e1 = v[1] - E[j + nE], e2 = v[2] - E[j + 2 * nE];
            const double d2 = (e0 * e0 + e1 * e1) + e2 * e2;
            if (d2 < b1 || (d2 == b1 && j < i1)) {
                b2 = b1;
                i2 = i1;
                b1 = d2;
                i1 = j;
            } else if (d2 < b2 || (d2 == b2 && j < i2)) {
                b2 = d2;
                i2 = j;
            }
        };
        const int a = cof(v[0]), b = cof(v[1]), c = cof(v[2]);
        bool done = false;
        for (int sh = 0; sh < g && !done; sh++) {
            for (int kc = std::max(0, c - sh); kc <= std::min(g - 1, c + sh); kc++) {
                const bool zf = std::abs(kc - c) == sh;
                for (int kb = std::max(0, b - sh); kb <= std::min(g - 1, b + sh); kb++) {
                    const bool face = zf || std::abs(kb - b) == sh;
                    for (int ka = std::max(0, a - sh); ka <= std::min(g - 1, a + sh); ka++) {
                        if (!face && std::abs(ka - a) != sh) continue;
                        const idx_t ci = (static_cast<idx_t>(kc) * g + kb) * g + ka;
                        for (int32_t s = cstart[ci]; s < cstart[ci + 1]; s++) consider(cpts[s]);
                    }
                }
            }
            // points not visited are at least sh cells away along an axis
            const double gap = sh * cell;
            done = i2 >= 0 && b2 < gap * gap * (1 - 1e-9);
        }
        nn[2 * q] = i1;
        nn[2 * q + 1] = i2;
        if (q < nE) {
            dd[2 * q] = std::sqrt(b1);
            dd[2 * q + 1] = std::sqrt(b2);
        }
    }
    // the bounds assume positive weights and radii; otherwise MATLAB's code for every voxel
    for (idx_t j = 0; j < nE; j++) {
        if (!(dd[2 * j + 1] > 0) || !(R[j] > 0)) return false;
    }
    // rn from the 2 nearest edge points, weighted with the distances of the first voxels (MATLAB's
    // indexing), and the window (x < l) + cos(pi (x - l) / 2 / (u - l))^2 (l <= x < u)
#pragma omp parallel for num_threads(nThreads) schedule(static)
    for (idx_t q = nE; q < nQ; q++) {
        const int32_t j1 = nn[2 * q], j2 = nn[2 * q + 1];
        const double d1 = dd[2 * j1], d2 = dd[2 * j2 + 1];
        const double rn = (R[j1] * d2 + R[j2] * d1) / (d1 + d2);
        const double x = rq[q - nE] / rn;
        double w = x < l ? 1.0 : 0.0;
        if (x >= l && x < u) {
            const double cs = fd_sincos(((kPi * (x - l)) / 2.0) / (u - l), false);
            w = cs * cs + w;
        }
        const int64_t t = uinds[static_cast<size_t>(q - nE)];
        const idx_t z = static_cast<idx_t>(t / np), p = static_cast<idx_t>(t % np);
        mask[z * np + (p % n1) * n2 + p / n1] = w;
    }
    return true;
}

}  // namespace pk5d::decon
