// Skewed space interpolation to a defined step size (skewed_space_interp_defined_stepsize_mex)
// on C-ordered (z, y, x) volumes: single to single, uint16 to single, uint16 to uint16. Grid
// slice z' at zout = stepsize * z' is a copy of a raw slice (within 1e-3) or interpolates raw
// slices zorg and zorg + 1 shifted in x. The mex sums y lines (one x column at a time); here
// rows of x are contiguous, so the lanes run along x and every voxel keeps its arithmetic: the
// formula of its x column's range, and PAIR for rows y below floor(ny / 16) * 16 (lanes.hpp).
// xstep is single precision, as in the mex. Grid slices are split among threads statically.
#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

#include "lanes.hpp"
#include "skewed_warp.hpp"

namespace pk5d::deskew {

// parameters of the grid slices (the mex's interp_slice_parameters)
inline void interp_slice_parameters(const float xstep, const double stepsize, const uint8_t Reverse,
                                    const uint64_t nint, double* s_mat, double* t_mat,
                                    double* sw_mat, double* tw_mat, double* zw_mat,
                                    bool* z_copy_mat, uint64_t* z_ind_mat) {
    // ceil of the single xstep in double, as C's ceil in the mex
    const double xstep_ceil = ceil((double)xstep);
    for (uint64_t i = 0; i < nint; i++) {
        double s;
        double t;
        double zout;
        double zw;
        zout = stepsize * (double)i;
        zw = zout - floor(zout);

        if (fabs(zout - round(zout)) < 1e-3) {
            z_copy_mat[i] = true;
            z_ind_mat[i] = (uint64_t)round(zout);
            continue;
        }

        if (Reverse) {
            s = xstep_ceil - xstep * zw;
            t = xstep * (1 - zw);
        } else {
            s = xstep * zw;
            t = xstep_ceil - xstep * (1 - zw);
        }

        // distance to start
        double sw = s - floor(s);
        s = floor(s);
        double tw = t - floor(t);
        t = floor(t);

        z_copy_mat[i] = false;
        s_mat[i] = s;
        t_mat[i] = t;
        sw_mat[i] = sw;
        tw_mat[i] = tw;
        zw_mat[i] = zw;
        z_ind_mat[i] = floor(zout);
    }
}

// a raw slice copied into the grid, by type
inline void copy_slice(const float* src, float* dst, uint64_t n) { memcpy(dst, src, n * sizeof(float)); }
inline void copy_slice(const uint16_t* src, uint16_t* dst, uint64_t n) {
    memcpy(dst, src, n * sizeof(uint16_t));
}
inline void copy_slice(const uint16_t* src, float* dst, uint64_t n) {
    for (uint64_t i = 0; i < n; ++i) dst[i] = (float)src[i];
}

// out[i] = wa * a[i] + wb * b[i] (rounded for uint16 out), i < n
template <typename T, typename TO>
inline void row_sum2(const T* a, const T* b, TO* out, float wa, float wb, uint64_t n) {
    const V8 va = set8(wa), vb = set8(wb);
    uint64_t i = 0;
    for (; i + 8 <= n; i += 8) {
        V8 v = sum2(load8_any(a + i), load8_any(b + i), va, vb);
        if (sizeof(TO) == 2) v = add_half(v);
        store8_any(out + i, v);
    }
    if (i < n) {
        const int m = (int)(n - i);
        V8 v = sum2(load_partial(a + i, m), load_partial(b + i, m), va, vb);
        if (sizeof(TO) == 2) v = add_half(v);
        store_partial(out + i, v, m);
    }
}

// the 4-term sums of a row (PAIR or SEQ)
template <bool PAIR, typename T, typename TO>
inline void row_sum4_mode(const T* a, const T* b, const T* c, const T* d, TO* out, float wa,
                          float wb, float wc, float wd, uint64_t n) {
    const V8 va = set8(wa), vb = set8(wb), vc = set8(wc), vd = set8(wd);
    uint64_t i = 0;
    for (; i < n; i += 8) {
        const int m = (int)std::min<uint64_t>(8, n - i);
        V8 A, B, C, D;
        if (m == 8) {
            A = load8_any(a + i), B = load8_any(b + i), C = load8_any(c + i), D = load8_any(d + i);
        } else {
            A = load_partial(a + i, m), B = load_partial(b + i, m), C = load_partial(c + i, m);
            D = load_partial(d + i, m);
        }
        V8 v = PAIR ? pair4(A, B, C, D, va, vb, vc, vd) : seq4(A, B, C, D, va, vb, vc, vd);
        if (sizeof(TO) == 2) v = add_half(v);
        if (m == 8)
            store8_any(out + i, v);
        else
            store_partial(out + i, v, m);
    }
}

template <typename T, typename TO>
inline void row_sum4(const T* a, const T* b, const T* c, const T* d, TO* out, float wa, float wb,
                     float wc, float wd, uint64_t n, bool pair) {
    if (pair)
        row_sum4_mode<true>(a, b, c, d, out, wa, wb, wc, wd, n);
    else
        row_sum4_mode<false>(a, b, c, d, out, wa, wb, wc, wd, n);
}

template <typename T, typename TO>
void skewed_space_interp_defined_stepsize(const T* im, TO* im_int, const float xstep,
                                          const double stepsize, const uint8_t Reverse,
                                          const uint64_t nz, const uint64_t ny, const uint64_t nx,
                                          const uint64_t nint, const int nthreads) {
    const uint64_t nxy = ny * nx;
    const uint64_t xa = (uint64_t)ceil((double)xstep);
    const uint64_t ny16 = ny - (ny % 16);
    std::vector<double> s_mat(nint), t_mat(nint), sw_mat(nint), tw_mat(nint), zw_mat(nint);
    std::unique_ptr<bool[]> z_copy_mat(new bool[nint > 0 ? nint : 1]);
    std::vector<uint64_t> z_ind_mat(nint);
    (void)nz;

    interp_slice_parameters(xstep, stepsize, Reverse, nint, s_mat.data(), t_mat.data(),
                            sw_mat.data(), tw_mat.data(), zw_mat.data(), z_copy_mat.get(),
                            z_ind_mat.data());

#pragma omp parallel for num_threads(nthreads)
    for (uint64_t z = 0; z < nint; z++) {
        const uint64_t zint = z;
        const uint64_t zorg = z_ind_mat[zint];
        if (z_copy_mat[zint]) {
            copy_slice(im + zorg * nxy, im_int + zint * nxy, nxy);
            continue;
        }

        const uint64_t s = (uint64_t)s_mat[zint];
        const uint64_t t = (uint64_t)t_mat[zint];
        const float sw = sw_mat[zint];
        const float tw = tw_mat[zint];
        const float zw = zw_mat[zint];

        const float c00_w = (1.0 - sw) * (1.0 - zw);
        const float c01_w = sw * (1.0 - zw);
        const float c10_w = (1.0 - tw) * zw;
        const float c11_w = tw * zw;

        // x ranges of the mex's three branches (its rows are x columns here)
        const uint64_t e1 = std::min<uint64_t>(xa - s, nx);
        const uint64_t e2 = std::max<uint64_t>(e1, nx - t - 1);
        const T* s0 = im + zorg * nxy;
        const T* s1 = im + (zorg + 1) * nxy;
        for (uint64_t y = 0; y < ny; y++) {
            const T* r0 = s0 + y * nx;
            const T* r1 = s1 + y * nx;
            TO* o = im_int + zint * nxy + y * nx;
            // the 4-term sums: pairwise in the rows below floor(ny / 16) * 16, except for
            // single input (always left to right)
            const bool pair = sizeof(T) == 2 && y < ny16;
            if (Reverse) {
                row_sum2(r1 + t, r1 + t + 1, o, c10_w, c11_w, e1);
                row_sum4(r0 + e1 + s - xa, r0 + e1 + s - xa + 1, r1 + e1 + t, r1 + e1 + t + 1, o + e1,
                         c00_w, c01_w, c10_w, c11_w, e2 - e1, pair);
                row_sum2(r0 + e2 + s - xa, r0 + e2 + s - xa + 1, o + e2, c00_w, c01_w, nx - e2);
            } else {
                row_sum2(r0 + s, r0 + s + 1, o, c00_w, c01_w, e1);
                row_sum4(r0 + e1 + s, r0 + e1 + s + 1, r1 + e1 + t - xa, r1 + e1 + t - xa + 1, o + e1,
                         c00_w, c01_w, c10_w, c11_w, e2 - e1, pair);
                row_sum2(r1 + e2 + t - xa, r1 + e2 + t - xa + 1, o + e2, c10_w, c11_w, nx - e2);
            }
        }
    }
}

}  // namespace pk5d::deskew
