// Skewed space interpolation to a defined step size: skewed_space_interp_defined_stepsize_mex.c
// (whole interpolated slices; its three type versions as one template: single to single,
// uint16 to single, uint16 to uint16). The loops, the split of the slices among threads
// (static) and the arithmetic are the mex's; xstep is single precision, as there. The mex's
// nested "#pragma omp parallel for" on the row loop is left out (it was serialized).
#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

#include "line_sums.hpp"
#include "skewed_warp.hpp"

namespace pk5d::deskew::ml {

// parameters of the slices of the interpolation grid at zout = stepsize * i, as in the combined
// mex (double precision step): a slice within 1e-3 of a raw slice is a copy of it; the others
// interpolate raw slices zorg and zorg + 1, read at x - xstep * zw (reverse) or x + xstep * zw
// (forward) and at x + xstep * (1 - zw) or x - xstep * (1 - zw), as in
// skewed_space_interp_defined_stepsize.m (s and t are ceil(xstep) - xstep < 1 apart before
// flooring either way, which the row ranges rely on)
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
inline void copy_slice(const float* src, float* dst, uint64_t n) {
    memcpy(dst, src, n * sizeof(float));
}

inline void copy_slice(const uint16_t* src, uint16_t* dst, uint64_t n) {
    memcpy(dst, src, n * sizeof(uint16_t));
}

inline void copy_slice(const uint16_t* src, float* dst, uint64_t n) {
    copy_uint16_to_float(src, dst, n);
}

template <typename T, typename TO>
void skewed_space_interp_defined_stepsize(const T* im, TO* im_int, const float xstep,
                                          const double stepsize, const uint8_t Reverse,
                                          const uint64_t sz[3], const uint64_t outDim[3],
                                          const int nthreads) {
    const uint64_t nint = outDim[2];
    const uint64_t nxy = sz[0] * sz[1];

    const uint64_t xa = (uint64_t)ceil((double)xstep);
    std::vector<double> s_mat(nint), t_mat(nint), sw_mat(nint), tw_mat(nint), zw_mat(nint);

    // copy or interp z indices
    std::unique_ptr<bool[]> z_copy_mat(new bool[nint > 0 ? nint : 1]);
    std::vector<uint64_t> z_ind_mat(nint);

    interp_slice_parameters(xstep, stepsize, Reverse, nint, s_mat.data(), t_mat.data(),
                            sw_mat.data(), tw_mat.data(), zw_mat.data(), z_copy_mat.get(),
                            z_ind_mat.data());

#pragma omp parallel for num_threads(nthreads)
    for (uint64_t z = 0; z < nint; z++) {
        uint64_t zint = z;
        uint64_t zorg = z_ind_mat[zint];

        bool z_copy = z_copy_mat[zint];
        if (z_copy) {
            copy_slice(im + zorg * nxy, im_int + zint * nxy, nxy);
            continue;
        }

        uint64_t s = (uint64_t)s_mat[zint];
        uint64_t t = (uint64_t)t_mat[zint];
        float sw = sw_mat[zint];
        float tw = tw_mat[zint];
        float zw = zw_mat[zint];

        float c00_w = (1.0 - sw) * (1.0 - zw);
        float c01_w = sw * (1.0 - zw);
        float c10_w = (1.0 - tw) * zw;
        float c11_w = tw * zw;

        if (Reverse) {
            for (uint64_t y = 0; y < sz[1]; y++) {
                if (y < xa - s) {
                    weighted_sum_2(im + (y + t) * sz[0] + (zorg + 1) * nxy,
                                   im + (y + t + 1) * sz[0] + (zorg + 1) * nxy,
                                   im_int + y * sz[0] + zint * nxy, c10_w, c11_w, sz[0]);
                } else if (y < sz[1] - t - 1) {
                    weighted_sum_4(im + (y + s - xa) * sz[0] + zorg * nxy,
                                   im + (y + s - xa + 1) * sz[0] + zorg * nxy,
                                   im + (y + t) * sz[0] + (zorg + 1) * nxy,
                                   im + (y + t + 1) * sz[0] + (zorg + 1) * nxy,
                                   im_int + y * sz[0] + zint * nxy, c00_w, c01_w, c10_w, c11_w,
                                   sz[0]);
                } else {
                    weighted_sum_2(im + (y + s - xa) * sz[0] + zorg * nxy,
                                   im + (y + s - xa + 1) * sz[0] + zorg * nxy,
                                   im_int + y * sz[0] + zint * nxy, c00_w, c01_w, sz[0]);
                }
            }
        } else {
            for (uint64_t y = 0; y < sz[1]; y++) {
                if (y < xa - s) {
                    weighted_sum_2(im + (y + s) * sz[0] + zorg * nxy,
                                   im + (y + s + 1) * sz[0] + zorg * nxy,
                                   im_int + y * sz[0] + zint * nxy, c00_w, c01_w, sz[0]);

                } else if (y < sz[1] - t - 1) {
                    weighted_sum_4(im + (y + s) * sz[0] + zorg * nxy,
                                   im + (y + s + 1) * sz[0] + zorg * nxy,
                                   im + (y + t - xa) * sz[0] + (zorg + 1) * nxy,
                                   im + (y + t - xa + 1) * sz[0] + (zorg + 1) * nxy,
                                   im_int + y * sz[0] + zint * nxy, c00_w, c01_w, c10_w, c11_w,
                                   sz[0]);
                } else {
                    weighted_sum_2(im + (y + t - xa) * sz[0] + (zorg + 1) * nxy,
                                   im + (y + t - xa + 1) * sz[0] + (zorg + 1) * nxy,
                                   im_int + y * sz[0] + zint * nxy, c10_w, c11_w, sz[0]);
                }
            }
        }
    }
}

}  // namespace pk5d::deskew::ml
