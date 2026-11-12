// Skewed space interpolation together with the deskew/rotate warp, without an interpolated
// volume: skewed_space_interp_volume_deskew_rotate_warp_mex.cpp (the fused kernel). The rows
// of the interpolated slices each output z slice needs are computed on the fly into a small
// per-thread cache. Loops, schedule, cache and arithmetic are the mex's.
#pragma once

#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "line_sums.hpp"

namespace pk5d::deskew::ml {

// line sums by row types (the overloads of the mex)

inline void weighted_sum_2(const float* a, const float* b, float* result, float wa, float wb,
                           size_t len) {
    compute_weighted_sum_two_vectors(a, b, result, wa, wb, len);
}

inline void weighted_sum_4(const float* a, const float* b, const float* c, const float* d,
                           float* result, float wa, float wb, float wc, float wd, size_t len) {
    compute_weighted_sum(a, b, c, d, result, wa, wb, wc, wd, len);
}

inline void weighted_sum_2(const uint16_t* a, const uint16_t* b, float* result, float wa,
                           float wb, size_t len) {
    compute_weighted_sum_two_vectors_16bit_to_float(a, b, result, wa, wb, len);
}

inline void weighted_sum_4(const uint16_t* a, const uint16_t* b, const uint16_t* c,
                           const uint16_t* d, float* result, float wa, float wb, float wc,
                           float wd, size_t len) {
    compute_weighted_sum_16bit_to_float(a, b, c, d, result, wa, wb, wc, wd, len);
}

// a, c: uint16; b, d: float
inline void weighted_sum_4(const uint16_t* a, const float* b, const uint16_t* c, const float* d,
                           float* result, float wa, float wb, float wc, float wd, size_t len) {
    compute_weighted_sum_16bit_float_to_float(a, b, c, d, result, wa, wb, wc, wd, len);
}

inline void weighted_sum_2(const uint16_t* a, const uint16_t* b, uint16_t* result, float wa,
                           float wb, size_t len) {
    compute_weighted_sum_two_vectors_16bit(a, b, result, wa, wb, len);
}

inline void weighted_sum_4(const uint16_t* a, const uint16_t* b, const uint16_t* c,
                           const uint16_t* d, uint16_t* result, float wa, float wb, float wc,
                           float wd, size_t len) {
    compute_weighted_sum_16bit(a, b, c, d, result, wa, wb, wc, wd, len);
}

// an interpolated slice: between raw slices zorg and zorg + 1, rows shifted by s and t
struct InterpSlice {
    uint64_t zorg;
    uint64_t s;
    uint64_t t;
    float c00_w;
    float c01_w;
    float c10_w;
    float c11_w;
};

// number of slices of the interpolation grid (none for an empty input), in double
inline double interp_grid_size(const uint64_t nz, const double stepsize) {
    return nz == 0 ? 0 : floor(round((double)(nz - 1) / stepsize * 100000) / 100000) + 1;
}

// which slices of the grid are raw slices (within 1e-3) and which are interpolated
inline void determine_interpolation_indexing(const uint64_t n, const double stepSize,
                                             uint64_t* interp_index, bool* is_original,
                                             uint64_t& nint) {
    uint64_t interp_counter = 0;

    for (uint64_t i = 0; i < n; i++) {
        double cur_step = (double)i * stepSize;
        double cur_step_round = round(cur_step);
        if (fabs(cur_step - cur_step_round) < 1e-3) {
            is_original[i] = true;
            interp_index[i] = (uint64_t)cur_step_round;
        } else {
            is_original[i] = false;
            interp_index[i] = interp_counter;
            interp_counter++;
        }
    }
    nint = interp_counter;
}

// parameters of the interpolated slices (indexed by interp_index)
inline void skewed_space_interp_slices(const double xstep, const double stepsize,
                                       const uint8_t Reverse, const uint64_t* interp_index,
                                       const bool* is_original, const uint64_t ntotal,
                                       InterpSlice* slices) {
    for (uint64_t i = 0; i < ntotal; i++) {
        if (is_original[i]) continue;

        double s;
        double t;
        double zout;
        double zw;
        uint64_t zint = interp_index[i];

        zout = stepsize * (double)i;
        zw = zout - floor(zout);

        // slice zorg is read at x - xstep * zw (reverse) or x + xstep * zw (forward), slice
        // zorg + 1 at x + xstep * (1 - zw) or x - xstep * (1 - zw), as in
        // skewed_space_interp_defined_stepsize.m (s and t are ceil(xstep) - xstep < 1 apart
        // before flooring either way, which the row ranges of skewed_space_interp_row rely on)
        if (Reverse) {
            s = ceil(xstep) - xstep * zw;
            t = xstep * (1 - zw);
        } else {
            s = xstep * zw;
            t = ceil(xstep) - xstep * (1 - zw);
        }

        // distance to start
        double sw = s - floor(s);
        s = floor(s);
        double tw = t - floor(t);
        t = floor(t);

        // single precision distances, weights computed in double
        float swf = sw;
        float twf = tw;
        float zwf = zw;

        InterpSlice& p = slices[zint];
        p.zorg = floor(zout);
        p.s = (uint64_t)s;
        p.t = (uint64_t)t;
        p.c00_w = (1.0 - swf) * (1.0 - zwf);
        p.c01_w = swf * (1.0 - zwf);
        p.c10_w = (1.0 - twf) * zwf;
        p.c11_w = twf * zwf;
    }
}

// row y (all sz[0] values) of an interpolated slice
template <typename T, typename TI>
void skewed_space_interp_row(const T* const im, TI* const row, const InterpSlice& p,
                             const uint64_t y, const uint64_t xa, const uint8_t Reverse,
                             const uint64_t sz[3]) {
    const uint64_t nxy = sz[0] * sz[1];
    const uint64_t zorg = p.zorg;
    const uint64_t s = p.s;
    const uint64_t t = p.t;

    if (Reverse) {
        if (y < xa - s) {
            weighted_sum_2(im + (y + t) * sz[0] + (zorg + 1) * nxy,
                           im + (y + t + 1) * sz[0] + (zorg + 1) * nxy, row, p.c10_w, p.c11_w,
                           sz[0]);
        } else if (y < sz[1] - t - 1) {
            weighted_sum_4(im + (y + s - xa) * sz[0] + zorg * nxy,
                           im + (y + s - xa + 1) * sz[0] + zorg * nxy,
                           im + (y + t) * sz[0] + (zorg + 1) * nxy,
                           im + (y + t + 1) * sz[0] + (zorg + 1) * nxy, row, p.c00_w, p.c01_w,
                           p.c10_w, p.c11_w, sz[0]);
        } else {
            weighted_sum_2(im + (y + s - xa) * sz[0] + zorg * nxy,
                           im + (y + s - xa + 1) * sz[0] + zorg * nxy, row, p.c00_w, p.c01_w,
                           sz[0]);
        }
    } else {
        if (y < xa - s) {
            weighted_sum_2(im + (y + s) * sz[0] + zorg * nxy, im + (y + s + 1) * sz[0] + zorg * nxy,
                           row, p.c00_w, p.c01_w, sz[0]);
        } else if (y < sz[1] - t - 1) {
            weighted_sum_4(im + (y + s) * sz[0] + zorg * nxy, im + (y + s + 1) * sz[0] + zorg * nxy,
                           im + (y + t - xa) * sz[0] + (zorg + 1) * nxy,
                           im + (y + t - xa + 1) * sz[0] + (zorg + 1) * nxy, row, p.c00_w,
                           p.c01_w, p.c10_w, p.c11_w, sz[0]);
        } else {
            weighted_sum_2(im + (y + t - xa) * sz[0] + (zorg + 1) * nxy,
                           im + (y + t - xa + 1) * sz[0] + (zorg + 1) * nxy, row, p.c10_w,
                           p.c11_w, sz[0]);
        }
    }
}

// output line from rows y0 and y0 + 1 (shapeX apart) of slices z0 (a) and z0 + 1 (b) of the
// interpolation grid; each slice is original (raw rows) or interpolated (cached rows: ai, bi)
inline void transform_line(const float* a, const float* ai, const float* b, const float* bi,
                           float* result, float w00, float w01, float w10, float w11,
                           uint64_t shapeX, uint64_t len) {
    const float* dmat = a ? a : ai;
    const float* dpmat = b ? b : bi;
    weighted_sum_4(dmat, dpmat, dmat + shapeX, dpmat + shapeX, result, w00, w01, w10, w11, len);
}

inline void transform_line(const uint16_t* a, const uint16_t* ai, const uint16_t* b,
                           const uint16_t* bi, uint16_t* result, float w00, float w01, float w10,
                           float w11, uint64_t shapeX, uint64_t len) {
    const uint16_t* dmat = a ? a : ai;
    const uint16_t* dpmat = b ? b : bi;
    weighted_sum_4(dmat, dpmat, dmat + shapeX, dpmat + shapeX, result, w00, w01, w10, w11, len);
}

// uint16 input, single output: interpolated rows are single
inline void transform_line(const uint16_t* a, const float* ai, const uint16_t* b, const float* bi,
                           float* result, float w00, float w01, float w10, float w11,
                           uint64_t shapeX, uint64_t len) {
    if (a && b) {
        weighted_sum_4(a, b, a + shapeX, b + shapeX, result, w00, w01, w10, w11, len);
    } else if (a) {
        weighted_sum_4(a, bi, a + shapeX, bi + shapeX, result, w00, w01, w10, w11, len);
    } else if (b) {
        weighted_sum_4(b, ai, b + shapeX, ai + shapeX, result, w01, w00, w11, w10, len);
    } else {
        weighted_sum_4(ai, bi, ai + shapeX, bi + shapeX, result, w00, w01, w10, w11, len);
    }
}

// Perform the transformation on a volume. Input x depends on output z only, so an output z
// slice reads two rows (y0, y0 + 1) of each slice of the interpolation grid; the two rows of
// the interpolated slices are computed when they are first needed and kept in a small
// per-thread cache (the slice index moves monotonically along output y), instead of
// interpolating whole slices into an intermediate volume first.
constexpr int NCACHE = 4;

template <typename T, typename TI, typename TO>
void transformVolume(const T* const im, const double (&tmat)[4][4], TO* const nmat,
                     const InterpSlice* const slices, const uint64_t* const interp_index,
                     const bool* const is_original, const uint64_t xa, const uint8_t Reverse,
                     const uint64_t sz[3], const uint64_t shapeZ, const float startX,
                     const float startY, const float startZ, const uint64_t outShapeX,
                     const uint64_t outShapeY, const uint64_t outShapeZ, const uint64_t lineLen,
                     TI* const cacheAll, const int nthreads) {
    // a bbox past the input's rows: nothing to compute
    if (lineLen == 0) return;

    const uint64_t shapeX = sz[0];
    const uint64_t shapeY = sz[1];
    const uint64_t shapeXY = shapeX * shapeY;
    const uint64_t outShapeXY = outShapeX * outShapeY;
    const uint64_t x0 = floor(startX);

#pragma omp parallel num_threads(nthreads)
    {
#ifdef _OPENMP
        const uint64_t thread = (uint64_t)omp_get_thread_num();
#else
        const uint64_t thread = 0;
#endif
        TI* cache = cacheAll + thread * NCACHE * 2 * shapeX;
        int64_t tag[NCACHE];

#pragma omp for schedule(dynamic, 1)
        for (uint64_t z = 0; z < outShapeZ; ++z) {
            double ty;
            ty = tmat[1][2] * (z + startZ) + tmat[1][3];
            if (ty < 0 || ty + 1 >= shapeY) {
                if (ty + 1 == shapeY && shapeY > 1) {  // on the last row (a single row: outside)
                    uint64_t counter = 1;
                    while (floor(ty) + 1 == shapeY) {
                        ty -= 1e-9 * counter;
                        counter *= 2;
                    }
                } else
                    continue;
            }
            int64_t y0 = floor(ty);
            float wy = ty - y0;

            for (int k = 0; k < NCACHE; k++) tag[k] = -1;

            for (uint64_t y = 0; y < outShapeY; ++y) {
                double tz;
                tz = tmat[2][1] * (y + startY) + tmat[2][2] * (z + startZ) + tmat[2][3];
                if (tz < 0 || tz + 1 >= shapeZ) {
                    // handling of edge cases
                    if (tz + 1 == shapeZ && shapeZ > 1) {  // on the last slice (a single slice: outside)
                        uint64_t counter = 1;
                        while (floor(tz) + 1 == shapeZ) {
                            tz -= 1e-9 * counter;
                            counter *= 2;
                        }
                    } else
                        continue;
                }
                uint64_t z0 = floor(tz);
                float wz = tz - z0;

                float w00 = (1 - wy) * (1 - wz);
                float w01 = (1 - wy) * wz;
                float w10 = wy * (1 - wz);
                float w11 = wy * wz;

                // rows y0 and y0 + 1 of slices z0 and z0 + 1, from x0 on
                const T* rows[2] = {nullptr, nullptr};
                const TI* cached[2] = {nullptr, nullptr};
                for (int q = 0; q < 2; q++) {
                    const uint64_t zq = z0 + q;
                    if (is_original[zq]) {
                        rows[q] = im + x0 + y0 * shapeX + interp_index[zq] * shapeXY;
                    } else {
                        const uint64_t zint = interp_index[zq];
                        const int k = zint % NCACHE;
                        TI* c = cache + k * 2 * shapeX;
                        if (tag[k] != (int64_t)zint) {
                            skewed_space_interp_row(im, c, slices[zint], y0, xa, Reverse, sz);
                            skewed_space_interp_row(im, c + shapeX, slices[zint], y0 + 1, xa,
                                                    Reverse, sz);
                            tag[k] = zint;
                        }
                        cached[q] = c + x0;
                    }
                }

                transform_line(rows[0], cached[0], rows[1], cached[1],
                               nmat + y * outShapeX + z * outShapeXY, w00, w01, w10, w11, shapeX,
                               lineLen);
            }
        }
    }
}

// the whole fused kernel: grid indexing, slice parameters, per-thread cache, warp
template <typename T, typename TI, typename TO>
void skewed_space_interp_warp(const T* const im, const uint64_t sz[3], const double xstep,
                              const double stepsize, const uint8_t Reverse,
                              const double (&tmat)[4][4], TO* const nmat, const float startX,
                              const float startY, const float startZ, const uint64_t outShapeX,
                              const uint64_t outShapeY, const uint64_t outShapeZ,
                              const uint64_t lineLen, const int nthreads) {
    const uint64_t ntotal = (uint64_t)interp_grid_size(sz[2], stepsize);

    // determine which slices come from original and which come from interpolated
    std::vector<uint64_t> interp_index(ntotal);
    std::unique_ptr<bool[]> is_original(new bool[ntotal > 0 ? ntotal : 1]);
    uint64_t nint = 0;
    determine_interpolation_indexing(ntotal, stepsize, interp_index.data(), is_original.get(),
                                     nint);

    std::vector<InterpSlice> slices(nint > 0 ? nint : 1);
    skewed_space_interp_slices(xstep, stepsize, Reverse, interp_index.data(), is_original.get(),
                               ntotal, slices.data());
    const uint64_t xa = (uint64_t)ceil(xstep);

    // cached rows of interpolated slices, per thread
    std::vector<TI> cache((uint64_t)nthreads * NCACHE * 2 * sz[0]);

    transformVolume<T, TI, TO>(im, tmat, nmat, slices.data(), interp_index.data(),
                               is_original.get(), xa, Reverse, sz, ntotal, startX, startY, startZ,
                               outShapeX, outShapeY, outShapeZ, lineLen, cache.data(), nthreads);
}

}  // namespace pk5d::deskew::ml
