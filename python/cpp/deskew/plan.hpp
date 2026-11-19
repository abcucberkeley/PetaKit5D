// Output geometry of the deskew/rotate warps.
//
// Names follow the mex, whose arrays are MATLAB's (y, x, z): the mex's "x" is MATLAB's y (the
// lines it sums), its "y" MATLAB's x, its "z" z. Here volumes are C-ordered (z, y, x), so the
// mex's lines are y columns of stride nx. The geometry is the mex's: output slice zo reads input
// columns x = col, col + 1 (col from tmat row 2, the nudge at the last column included); output
// line (zo, xo) reads slices z0, z0 + 1 of the (interpolation) grid (tmat row 3, the nudge at
// the last slice included), with the weights w00 = (1 - wy)(1 - wz), w01 = (1 - wy) wz,
// w10 = wy (1 - wz), w11 = wy wz in single precision. tmat keeps MATLAB's convention.
#pragma once

#include <cmath>
#include <cstdint>

namespace pk5d::deskew {

// The mex's geometry. tmat rows 2 and 3 give, for output slice z (start startZ), the input
// column ty and, for its line y (start startY), the grid slice tz; the nudges keep the last
// column and the last slice inside.
struct Geometry {
    double tmat[4][4];
    uint64_t shapeY, shapeZ;  // input columns, grid slices
    float startY, startZ;

    // column of slice z and its weight wy; false: the slice is outside
    bool slice(const uint64_t z, int64_t& y0, float& wy) const {
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
                return false;
        }
        y0 = floor(ty);
        wy = ty - y0;
        return true;
    }

    // grid slice of line y of slice z and the weights w00, w01, w10, w11; false: outside
    bool line(const uint64_t z, const uint64_t y, const float wy, int64_t& z0, float (&w)[4]) const {
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
                return false;
        }
        z0 = floor(tz);
        float wz = tz - z0;
        w[0] = (1 - wy) * (1 - wz);
        w[1] = (1 - wy) * wz;
        w[2] = wy * (1 - wz);
        w[3] = wy * wz;
        return true;
    }
};

}  // namespace pk5d::deskew
