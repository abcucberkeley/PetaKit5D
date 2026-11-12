// The deskew/rotate warp without interpolation in skewed space: transformVolume of
// volume_deskew_rotate_warp_mex.cpp (its three type versions as one template: single to
// single, uint16 to single, uint16 to uint16; the loops, the schedule and the arithmetic are
// the mex's). Output line (y, z) reads input rows y0, y0 + 1 of slices z0, z0 + 1 from row
// x0 = floor(startX) on, lineLen values.
#pragma once

#include <cmath>
#include <cstdint>

#include "line_sums.hpp"

namespace pk5d::deskew::ml {

// line sums by type: single, uint16 to single, uint16 to uint16
inline void warp_line(const float* a, const float* b, const float* c, const float* d, float* r,
                      float wa, float wb, float wc, float wd, uint64_t len) {
    compute_weighted_sum(a, b, c, d, r, wa, wb, wc, wd, len);
}

inline void warp_line(const uint16_t* a, const uint16_t* b, const uint16_t* c, const uint16_t* d,
                      float* r, float wa, float wb, float wc, float wd, uint64_t len) {
    compute_weighted_sum_16bit_to_float(a, b, c, d, r, wa, wb, wc, wd, len);
}

inline void warp_line(const uint16_t* a, const uint16_t* b, const uint16_t* c, const uint16_t* d,
                      uint16_t* r, float wa, float wb, float wc, float wd, uint64_t len) {
    compute_weighted_sum_16bit(a, b, c, d, r, wa, wb, wc, wd, len);
}

// Perform the transformation on a volume
template <typename T, typename TO>
void transformVolume(const T* const dmat, const double (&tmat)[4][4], TO* const nmat,
                     const uint64_t shapeX, const uint64_t shapeY, const uint64_t shapeZ,
                     const float startX, const float startY, const float startZ,
                     const uint64_t outShapeX, const uint64_t outShapeY, const uint64_t outShapeZ,
                     const uint64_t lineLen, const int nthreads) {
    // a bbox past the input's rows: nothing to compute
    if (lineLen == 0) return;

    uint64_t shapeXY = shapeX * shapeY;
    uint64_t outShapeXY = outShapeX * outShapeY;

    int64_t x0 = floor(startX);

#pragma omp parallel for schedule(dynamic, 1) num_threads(nthreads)
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

            int64_t z0 = floor(tz);
            float wz = tz - z0;

            float w00 = (1 - wy) * (1 - wz);
            float w01 = (1 - wy) * wz;
            float w10 = wy * (1 - wz);
            float w11 = wy * wz;

            uint64_t ind = x0 + y0 * shapeX + z0 * shapeXY;

            // interpolation with simd vectorization
            warp_line(dmat + ind, dmat + ind + shapeXY, dmat + ind + shapeX,
                      dmat + ind + shapeX + shapeXY, nmat + y * outShapeX + z * outShapeXY, w00,
                      w01, w10, w11, lineLen);
        }
    }
}

}  // namespace pk5d::deskew::ml
