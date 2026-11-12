// Resize of a feather distance map to a block's size by trilinear interpolation, then raised
// to the integer weight degree: feather_distance_map_resize_3d_mex.c, and its older variant
// feather_distance_map_resize_3d_w_space_mex.c (the map is first taken back to the linear
// space with powf(d, 1 / wd)). The single-precision arithmetic is the mex's: the source
// coordinate of each axis is (float)i * factor with factor = (n - 1) / (rn - 1) computed in
// double, and the interpolation goes along y (MATLAB's first dimension), then x, then z.
//
// dmat: (nz, ny, nx), out: (rz, ry, rx), C order; the work is split by output (z, y) lines.
#include "stitch.hpp"

#include <cmath>
#include <memory>

namespace {

using namespace pk5d_stitch;

// (int) of a float as the x86-64 mex computes it (cvttss2si: INT32_MIN for NaN and values
// outside int32)
inline int x86_float_to_int(float v) {
    return (v > -2147483904.0f && v < 2147483648.0f) ? (int)v : INT32_MIN;
}

inline float fastPower(float base, int exponent) {
    // Only for integer exponent
    if (exponent < 0) {
        base = 1.0 / base;
        // -INT32_MIN wraps to INT32_MIN on x86 (the loop below then does nothing)
        exponent = (int)(0u - (unsigned)exponent);
    }

    float result = 1.0;
    while (exponent > 0) {
        if (exponent & 1) {
            result *= base;
        }
        base *= base;
        exponent >>= 1;
    }

    return result;
}

// the mex's factor of one axis (n: the map's size, rn: the output's)
inline float axis_factor(uint64_t n, uint64_t rn) {
    return rn > 1 ? ((float)n - 1.0) / ((float)rn - 1.0) : 1.0;
}

void resize_loop(const float *dmat, float *rmat, const int d_int, const uint64_t *dims,
                 const uint64_t *rdims, int nthread) {
    const uint64_t nZ = dims[0], nY = dims[1], nX = dims[2];
    const uint64_t rZ = rdims[0], rY = rdims[1], rX = rdims[2];
    const uint64_t nYX = nY * nX;
    // the mex's xfactor, yfactor, zfactor: MATLAB's y, x and z
    const float yfactor = axis_factor(nY, rY);
    const float xfactor = axis_factor(nX, rX);
    const float zfactor = axis_factor(nZ, rZ);
    const int64_t nLines = (int64_t)(rZ * rY);

#pragma omp parallel for schedule(static) num_threads(nthread) if (nthread > 1)
    for (int64_t l = 0; l < nLines; l++) {
        const uint64_t z = (uint64_t)l / rY, y = (uint64_t)l % rY;
        const float zr = (float)z * zfactor;
        const float yr = (float)y * yfactor;
        const uint64_t z0 = (uint64_t)zr, y0 = (uint64_t)yr;
        const uint64_t zs = z0 + 1 < nZ ? nYX : 0;
        const uint64_t ys = y0 + 1 < nY ? nX : 0;
        const float zd = zr - z0;
        const float yd = yr - y0;
        float *o = rmat + (z * rY + y) * rX;
        for (uint64_t x = 0; x < rX; x++) {
            const float xr = (float)x * xfactor;
            const uint64_t x0 = (uint64_t)xr;
            const uint64_t xs = x0 + 1 < nX ? 1 : 0;
            const float xd = xr - x0;
            const float *p = dmat + z0 * nYX + y0 * nX + x0;

            // the mex's c00, c10, c01, c11: along y first, then x, then z
            float c00 = *(p) * (1 - yd) + *(p + ys) * yd;
            float c10 = *(p + xs) * (1 - yd) + *(p + ys + xs) * yd;
            float c01 = *(p + zs) * (1 - yd) + *(p + ys + zs) * yd;
            float c11 = *(p + xs + zs) * (1 - yd) + *(p + ys + xs + zs) * yd;

            float c0 = c00 * (1 - xd) + c10 * xd;
            float c1 = c01 * (1 - xd) + c11 * xd;

            float c = c0 * (1 - zd) + c1 * zd;
            o[x] = fastPower(c, d_int);
        }
    }
}

// feather_distance_map_resize_3d_mex
void feather_distance_map_resize_3d(const float *dmat, float *rmat, const float d,
                                    const uint64_t *dims, const uint64_t *rdims, int nthread) {
    resize_loop(dmat, rmat, x86_float_to_int(d), dims, rdims, nthread);
}

// feather_distance_map_resize_3d_w_space_mex
void feather_distance_map_resize_3d_w_space(const float *dmat, float *rmat, const float d,
                                            const uint64_t *dims, const uint64_t *rdims,
                                            int nthread) {
    const int64_t n = (int64_t)(dims[0] * dims[1] * dims[2]);
    const float d_inv = 1. / (float)d;

    // restore to the linear space
    std::unique_ptr<float[]> dmat_1(new float[n > 0 ? n : 1]);
    float *d1 = dmat_1.get();
#pragma omp parallel for schedule(static) num_threads(nthread) if (nthread > 1)
    for (int64_t i = 0; i < n; i++) {
        d1[i] = powf(dmat[i], d_inv);
    }
    resize_loop(d1, rmat, x86_float_to_int(d), dims, rdims, nthread);
}

template <void (*K)(const float *, float *, const float, const uint64_t *, const uint64_t *, int)>
auto resize_binding() {
    return [](in1d<float> dmat, out1d<float> out, std::vector<uint64_t> dims,
              std::vector<uint64_t> rdims, float wd, int num_threads) {
        pk5d::check_threads(num_threads);
        need_dims(dims, 3, "dims");
        need_dims(rdims, 3, "rdims");
        need_size(dmat, prod(dims), "distance map");
        need_size(out, prod(rdims), "output");
        for (int i = 0; i < 3; i++) need(dims[i] <= rdims[i], "the output must not be smaller");
        need(prod(dims) > 0 || prod(rdims) == 0, "the distance map is empty");
        nb::gil_scoped_release release;
        K(dmat.data(), out.data(), wd, dims.data(), rdims.data(), num_threads);
    };
}

}  // namespace

PK5D_BINDINGS(bind_distance_map_resize) {
    m.def("feather_distance_map_resize_3d", resize_binding<feather_distance_map_resize_3d>(),
          nb::arg("dmat").noconvert(), nb::arg("out").noconvert(), nb::arg("dims"),
          nb::arg("rdims"), nb::arg("wd"), nb::arg("num_threads"),
          "Trilinear resize of dmat ((z, y, x) sizes dims) into out (sizes rdims), raised to "
          "int(wd).");
    m.def("feather_distance_map_resize_3d_w_space",
          resize_binding<feather_distance_map_resize_3d_w_space>(), nb::arg("dmat").noconvert(),
          nb::arg("out").noconvert(), nb::arg("dims"), nb::arg("rdims"), nb::arg("wd"),
          nb::arg("num_threads"),
          "As feather_distance_map_resize_3d, with dmat ** (1 / wd) interpolated.");
}
