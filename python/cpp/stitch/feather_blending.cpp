// Feather blending of the tiles of a stitching batch: feather_blending_3d_mex.cpp and
// feather_blending_3d_with_indexing_mex.cpp. Per voxel, the tiles with a nonzero value are
// averaged with their distance-map weights, in tile order and in single precision:
//   for each tile t: if f != 0 { sum_w += d; nv += (float)f * d; }   then nv / sum_w
// with integer results rounded by + 0.5 in double and a truncating conversion. Voxels without a
// weight keep the output's value (0 in a new output).
//
// tiles, weights: (t, z, y, x), C order; out: (oz, oy, ox), the region written at offset
// (z0, y0, x0). Each voxel is computed as in the mex (the tile loop inside, in tile order); the
// work is split by (z, y) lines, statically, each thread running a plain function on its share
// (an OpenMP region around the loop kept GCC from optimizing it: 1.24x the mex at 1 thread), and
// the tile loop has a fixed trip count for 1 to 8 tiles (unrolled; any count runs the same
// operations).
#include "stitch.hpp"

#include <type_traits>

namespace {

using namespace pk5d_stitch;

// (uint16_t) of a double as the x86-64 mex computes it (cvttsd2si: truncation to int32,
// INT32_MIN for NaN and values outside int32, then the low 16 bits), written out so that other
// processors give the same values for weights that leave [0, 65535]
inline uint16_t x86_double_to_uint16(double v) {
    const int32_t i = (v > -2147483649.0 && v < 2147483648.0) ? (int32_t)v : INT32_MIN;
    return (uint16_t)(uint32_t)i;
}

// the mex's round_output: integer classes get + 0.5 (in double) and a truncating conversion
template <typename T>
inline T round_output(float value) {
    if constexpr (std::is_same<T, uint16_t>::value) {
        return x86_double_to_uint16(value + 0.5);
    } else {
        return static_cast<T>(value);
    }
}

// lines [l0, l1) of the region; NT > 0: that many tiles, fixed at compile time
template <typename T, int NT>
NB_NOINLINE void feather_lines(const T *fmat, const float *dmat, T *nmat, const uint64_t nT_rt,
                               const uint64_t nY, const uint64_t nX, const uint64_t nZYX,
                               const uint64_t oY, const uint64_t oX, const uint64_t oz,
                               const uint64_t oy, const uint64_t ox, const int64_t l0,
                               const int64_t l1) {
    const uint64_t nT = NT > 0 ? (uint64_t)NT : nT_rt;
    for (int64_t l = l0; l < l1; l++) {
        const uint64_t z = (uint64_t)l / nY, y = (uint64_t)l % nY;
        const uint64_t line = (uint64_t)l * nX;
        T *o = nmat + ((z + oz) * oY + (y + oy)) * oX + ox;
        for (uint64_t x = 0; x < nX; x++) {
            float sum_w = 0;
            float nv = 0;
            const T *f_ptr = fmat + line + x;
            const float *d_ptr = dmat + line + x;
            for (uint64_t t = 0; t < nT; t++) {
                if (*f_ptr != 0) {
                    // Inf * 0, Inf / Inf: MATLAB's NaN on every processor (fpexact.hpp)
                    const float fv = static_cast<float>(*f_ptr), d = *d_ptr;
                    sum_w = PK5D_X86_NAN(sum_w + d, sum_w, d);
                    const float p = PK5D_X86_NAN(fv * d, fv, d);
                    nv = PK5D_X86_NAN(nv + p, nv, p);
                }
                f_ptr += nZYX;
                d_ptr += nZYX;
            }
            if (sum_w != 0) {
                o[x] = round_output<T>(PK5D_X86_NAN(nv / sum_w, nv, sum_w));
            }
        }
    }
}

template <typename T>
void feather_range(const T *fmat, const float *dmat, T *nmat, const uint64_t *dims,
                   const uint64_t *odims, const uint64_t *offset, const int64_t l0,
                   const int64_t l1) {
    const uint64_t nT = dims[0], nY = dims[2], nX = dims[3];
    const uint64_t nZYX = dims[1] * nY * nX;
    const uint64_t oY = odims[1], oX = odims[2];
#define PK5D_FEATHER(N)                                                                      \
    feather_lines<T, N>(fmat, dmat, nmat, nT, nY, nX, nZYX, oY, oX, offset[0], offset[1], \
                        offset[2], l0, l1)
    switch (nT) {
        case 1: PK5D_FEATHER(1); break;
        case 2: PK5D_FEATHER(2); break;
        case 3: PK5D_FEATHER(3); break;
        case 4: PK5D_FEATHER(4); break;
        case 5: PK5D_FEATHER(5); break;
        case 6: PK5D_FEATHER(6); break;
        case 7: PK5D_FEATHER(7); break;
        case 8: PK5D_FEATHER(8); break;
        default: PK5D_FEATHER(0); break;
    }
#undef PK5D_FEATHER
}

template <typename T>
void feather_blending_3d(const T *fmat, const float *dmat, T *nmat, const uint64_t *dims,
                         const uint64_t *odims, const uint64_t *offset, int nthread) {
    const int64_t nLines = (int64_t)(dims[1] * dims[2]);
    if (nthread <= 1) {
        feather_range<T>(fmat, dmat, nmat, dims, odims, offset, 0, nLines);
        return;
    }
#pragma omp parallel num_threads(nthread)
    {
        // OpenMP's static split of the lines among the team
        int64_t k = 0, m = 1;
#ifdef _OPENMP
        k = omp_get_thread_num();
        m = omp_get_num_threads();
#endif
        const int64_t q = nLines / m, r = nLines % m;
        const int64_t l0 = k * q + (k < r ? k : r), l1 = l0 + q + (k < r ? 1 : 0);
        feather_range<T>(fmat, dmat, nmat, dims, odims, offset, l0, l1);
    }
}

// dims: the tiles' (t, z, y, x) sizes; odims: the output's (z, y, x) sizes; offset: 0-based
// start of the region in the output
template <typename T>
void bind_feather(nb::module_ &m) {
    m.def(
        "feather_blending_3d",
        [](in1d<T> f, in1d<float> d, out1d<T> out, std::vector<uint64_t> dims,
           std::vector<uint64_t> odims, std::vector<uint64_t> offset, int num_threads) {
            pk5d::check_threads(num_threads);
            need_dims(dims, 4, "dims");
            need_dims(odims, 3, "odims");
            need_dims(offset, 3, "offset");
            need_size(f, prod(dims), "tiles");
            need_size(d, prod(dims), "distance map");
            need_size(out, prod(odims), "output");
            for (int i = 0; i < 3; i++)
                need(offset[i] + dims[i + 1] <= odims[i] || prod(dims) == 0,
                     "the region does not fit in the output");
            nb::gil_scoped_release release;
            feather_blending_3d<T>(f.data(), d.data(), out.data(), dims.data(), odims.data(),
                                   offset.data(), num_threads);
        },
        nb::arg("f").noconvert(), nb::arg("d").noconvert(), nb::arg("out").noconvert(),
        nb::arg("dims"), nb::arg("odims"), nb::arg("offset"), nb::arg("num_threads"),
        "Feather blending of the tiles f ((t, z, y, x) sizes dims) with weights d into out "
        "((z, y, x) sizes odims) at offset (0-based).");
}

}  // namespace

PK5D_BINDINGS(bind_feather_blending) {
    bind_feather<uint16_t>(m);
    bind_feather<float>(m);
}
