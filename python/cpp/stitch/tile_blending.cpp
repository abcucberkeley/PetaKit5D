// Blending of the tiles of a stitching batch for the blend methods 'none', 'max', 'mean' and
// 'median' in one pass: tile_blending_3d_mex.cpp, with its per-voxel code and single-precision
// arithmetic unchanged. tiles: (t, z, y, x), C order, 0 = no data, uint8, uint16 or single.
// boxes: regions (t, z0, y0, x0, z1, y1, x1), 0-based, end-exclusive, where a tile counts as 0.
// A voxel gets the blend of the tiles outside their regions, or, where that is 0, the blend of
// all tiles. 'none' and 'max' give the class of the tiles, 'mean' and 'median' single. The work
// is split by (z, y) lines along x (the mex's lines run along MATLAB's first dimension, y); a
// line's masks come from the regions that cross it.
//
// The single arithmetic follows MATLAB's: values are summed in tile order, the mean of the two
// middle values of an even count is a + (b - a) / 2, or (a + b) / 2 if their signs differ or one
// is infinite; no fused multiply-adds (-ffp-contract=off, as the mex is built).
#include "stitch.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <type_traits>

namespace {

using namespace pk5d_stitch;

enum Method { NONE, MAX, MEAN, MEDIAN };

// a region where a tile counts as 0 (0-based, end-exclusive, not empty)
struct Box {
    uint64_t t, z0, y0, x0, z1, y1, x1;
};

// value of tile j at the voxel: f and m point to the voxel's value and mask in tile 0, nXYZ
// and nX are the strides to the next tile (the tile's volume, the line's length)
#define TVAL(j) ((MASK && m[(j) * nX]) ? (T)0 : f[(j) * nXYZ])

// MATLAB's sign
static inline int sgn(const float v) { return (v > 0) - (v < 0); }

// MATLAB's median of n > 0 values without NaN
static inline float medianOf(float *buf, const uint64_t n) {
    std::sort(buf, buf + n);
    if (n % 2 == 1) {
        return buf[n / 2];
    }
    const float a = buf[n / 2 - 1];
    const float b = buf[n / 2];
    if (sgn(a) != sgn(b) || std::isinf(a) || std::isinf(b)) {
        return (a + b) / 2.0f;
    }
    return a + (b - a) / 2.0f;
}

// the last step of the MATLAB code: nv, or nvf where nv is 0, as nv .* (1 - z) + nvf .* z
static inline float fillZeros(const float nv, const float nvf) {
    const float z = (nv == 0.0f && nvf != 0.0f) ? 1.0f : 0.0f;
    return nv * (1.0f - z) + nvf * z;
}

// Integer tiles: their values are exact in single, so 'none' is the first nonzero tile value
// and 'max' the largest, outside the regions or else of all tiles; 'mean' and 'median' are
// those of the nonzero values.
template <typename T, bool MASK>
static inline T noneInt(const T *f, const uint8_t *m, const uint64_t nX, const uint64_t nXYZ,
                        const uint64_t nT) {
    for (uint64_t j = 0; j < nT; j++) {
        const T v = TVAL(j);
        if (v != 0) return v;
    }
    for (uint64_t j = 0; j < nT; j++) {
        const T v = f[j * nXYZ];
        if (v != 0) return v;
    }
    return 0;
}

template <typename T, bool MASK>
static inline T maxInt(const T *f, const uint8_t *m, const uint64_t nX, const uint64_t nXYZ,
                       const uint64_t nT) {
    T mx = 0;
    for (uint64_t j = 0; j < nT; j++) {
        const T v = TVAL(j);
        if (v > mx) mx = v;
    }
    if (mx == 0) {
        for (uint64_t j = 0; j < nT; j++) {
            const T v = f[j * nXYZ];
            if (v > mx) mx = v;
        }
    }
    return mx;
}

template <typename T, bool MASK>
static inline float meanInt(const T *f, const uint8_t *m, const uint64_t nX, const uint64_t nXYZ,
                            const uint64_t nT) {
    float s = 0.0f;
    float n = 0.0f;
    for (uint64_t j = 0; j < nT; j++) {
        const T v = TVAL(j);
        if (v != 0) {
            s = s + (float)v;
            n = n + 1.0f;
        }
    }
    if (n == 0.0f) {
        for (uint64_t j = 0; j < nT; j++) {
            const T v = f[j * nXYZ];
            if (v != 0) {
                s = s + (float)v;
                n = n + 1.0f;
            }
        }
    }
    return n == 0.0f ? 0.0f : s / n;
}

template <typename T, bool MASK>
static inline float medianInt(const T *f, const uint8_t *m, const uint64_t nX,
                              const uint64_t nXYZ, const uint64_t nT, float *buf) {
    uint64_t n = 0;
    for (uint64_t j = 0; j < nT; j++) {
        const T v = TVAL(j);
        if (v != 0) buf[n++] = (float)v;
    }
    if (n == 0) {
        for (uint64_t j = 0; j < nT; j++) {
            const T v = f[j * nXYZ];
            if (v != 0) buf[n++] = (float)v;
        }
    }
    return n == 0 ? 0.0f : medianOf(buf, n);
}

// Single tiles: the operations of the MATLAB code (NaN, Inf and -0 included).
// 'none': nv = t1, then nv = nv + tj .* (nv == 0) for the other tiles.
template <typename T, bool MASK>
static inline float noneSingle(const T *f, const uint8_t *m, const uint64_t nX,
                               const uint64_t nXYZ, const uint64_t nT) {
    float nv = TVAL(0);
    float nvf = f[0];
    for (uint64_t j = 1; j < nT; j++) {
        nv = nv + TVAL(j) * (nv == 0.0f ? 1.0f : 0.0f);
        nvf = nvf + f[j * nXYZ] * (nvf == 0.0f ? 1.0f : 0.0f);
    }
    if (std::isnan(nvf)) nvf = 0.0f;
    return fillZeros(nv, nvf);
}

// 'max', 'mean' and 'median' leave out 0 and NaN; NaN results become 0
template <typename T, bool MASK>
static inline float maxSingle(const T *f, const uint8_t *m, const uint64_t nX,
                              const uint64_t nXYZ, const uint64_t nT) {
    float nv = NAN;
    float nvf = NAN;
    for (uint64_t j = 0; j < nT; j++) {
        const float v = TVAL(j);
        if (v != 0.0f && !std::isnan(v) && (std::isnan(nv) || v > nv)) nv = v;
        const float vf = f[j * nXYZ];
        if (vf != 0.0f && !std::isnan(vf) && (std::isnan(nvf) || vf > nvf)) nvf = vf;
    }
    if (std::isnan(nv)) nv = 0.0f;
    if (std::isnan(nvf)) nvf = 0.0f;
    return fillZeros(nv, nvf);
}

template <typename T, bool MASK>
static inline float meanSingle(const T *f, const uint8_t *m, const uint64_t nX,
                               const uint64_t nXYZ, const uint64_t nT) {
    float s = 0.0f, n = 0.0f, sf = 0.0f, nf = 0.0f;
    for (uint64_t j = 0; j < nT; j++) {
        const float v = TVAL(j);
        if (v != 0.0f && !std::isnan(v)) {
            s = s + v;
            n = n + 1.0f;
        }
        const float vf = f[j * nXYZ];
        if (vf != 0.0f && !std::isnan(vf)) {
            sf = sf + vf;
            nf = nf + 1.0f;
        }
    }
    float nv = s / n;
    float nvf = sf / nf;
    if (std::isnan(nv)) nv = 0.0f;
    if (std::isnan(nvf)) nvf = 0.0f;
    return fillZeros(nv, nvf);
}

template <typename T, bool MASK>
static inline float medianSingle(const T *f, const uint8_t *m, const uint64_t nX,
                                 const uint64_t nXYZ, const uint64_t nT, float *buf) {
    uint64_t n = 0;
    for (uint64_t j = 0; j < nT; j++) {
        const float v = TVAL(j);
        if (v != 0.0f && !std::isnan(v)) buf[n++] = v;
    }
    float nv = n == 0 ? NAN : medianOf(buf, n);
    n = 0;
    for (uint64_t j = 0; j < nT; j++) {
        const float v = f[j * nXYZ];
        if (v != 0.0f && !std::isnan(v)) buf[n++] = v;
    }
    float nvf = n == 0 ? NAN : medianOf(buf, n);
    if (std::isnan(nv)) nv = 0.0f;
    if (std::isnan(nvf)) nvf = 0.0f;
    return fillZeros(nv, nvf);
}

// one line of voxels along x; f and o point to its first voxel, mask to the masks of the line
// (nX per tile)
template <typename T, typename O, Method M, bool MASK>
static void blendLine(const T *f, const uint8_t *mask, O *o, const uint64_t nX,
                      const uint64_t nXYZ, const uint64_t nT, float *buf) {
    const bool isInt = !std::is_floating_point<T>::value;
    for (uint64_t x = 0; x < nX; x++) {
        const T *fx = f + x;
        const uint8_t *mx = MASK ? mask + x : NULL;
        switch (M) {
            case NONE:
                o[x] = isInt ? (O)noneInt<T, MASK>(fx, mx, nX, nXYZ, nT)
                             : (O)noneSingle<T, MASK>(fx, mx, nX, nXYZ, nT);
                break;
            case MAX:
                o[x] = isInt ? (O)maxInt<T, MASK>(fx, mx, nX, nXYZ, nT)
                             : (O)maxSingle<T, MASK>(fx, mx, nX, nXYZ, nT);
                break;
            case MEAN:
                o[x] = isInt ? (O)meanInt<T, MASK>(fx, mx, nX, nXYZ, nT)
                             : (O)meanSingle<T, MASK>(fx, mx, nX, nXYZ, nT);
                break;
            case MEDIAN:
                o[x] = isInt ? (O)medianInt<T, MASK>(fx, mx, nX, nXYZ, nT, buf)
                             : (O)medianSingle<T, MASK>(fx, mx, nX, nXYZ, nT, buf);
                break;
        }
    }
}

template <typename T, typename O, Method M>
static void blendTiles(const T *f, O *o, const uint64_t nT, const uint64_t nZ, const uint64_t nY,
                       const uint64_t nX, const std::vector<Box> &boxes, int nthread) {
    const uint64_t nXYZ = nZ * nY * nX;
    const int64_t nLines = (int64_t)(nZ * nY);
    const uint64_t nBoxes = boxes.size();

#pragma omp parallel num_threads(nthread)
    {
        std::vector<uint8_t> mask(nT * nX);
        std::vector<float> buf(nT);

#pragma omp for schedule(static)
        for (int64_t l = 0; l < nLines; l++) {
            const uint64_t z = (uint64_t)l / nY;
            const uint64_t y = (uint64_t)l % nY;
            // masks of the regions that cross the line
            bool masked = false;
            for (uint64_t k = 0; k < nBoxes; k++) {
                const Box &b = boxes[k];
                if (z < b.z0 || z >= b.z1 || y < b.y0 || y >= b.y1) continue;
                if (!masked) {
                    memset(mask.data(), 0, nT * nX);
                    masked = true;
                }
                memset(mask.data() + b.t * nX + b.x0, 1, b.x1 - b.x0);
            }
            const uint64_t offset = (z * nY + y) * nX;
            if (masked) {
                blendLine<T, O, M, true>(f + offset, mask.data(), o + offset, nX, nXYZ, nT,
                                         buf.data());
            } else {
                blendLine<T, O, M, false>(f + offset, NULL, o + offset, nX, nXYZ, nT,
                                          buf.data());
            }
        }
    }
}

template <typename T, typename O, Method M>
void bind_one(nb::module_ &m, const char *name) {
    m.def(
        name,
        [](in1d<T> f, out1d<O> out, std::vector<uint64_t> dims, std::vector<int64_t> boxes,
           int num_threads) {
            pk5d::check_threads(num_threads);
            need_dims(dims, 4, "dims");
            need(dims[0] > 0, "there are no tiles");
            need_size(f, prod(dims), "tiles");
            need_size(out, dims[1] * dims[2] * dims[3], "output");
            need(boxes.size() % 7 == 0, "boxes must have 7 values each");
            std::vector<Box> bx;
            for (size_t k = 0; k < boxes.size(); k += 7) {
                const int64_t *v = boxes.data() + k;
                need(v[0] >= 0 && (uint64_t)v[0] < dims[0] && v[1] >= 0 && v[2] >= 0 &&
                         v[3] >= 0 && v[1] < v[4] && v[2] < v[5] && v[3] < v[6] &&
                         (uint64_t)v[4] <= dims[1] && (uint64_t)v[5] <= dims[2] &&
                         (uint64_t)v[6] <= dims[3],
                     "a region is empty or outside the tiles");
                bx.push_back(Box{(uint64_t)v[0], (uint64_t)v[1], (uint64_t)v[2], (uint64_t)v[3],
                                 (uint64_t)v[4], (uint64_t)v[5], (uint64_t)v[6]});
            }
            nb::gil_scoped_release release;
            blendTiles<T, O, M>(f.data(), out.data(), dims[0], dims[1], dims[2], dims[3], bx,
                                num_threads);
        },
        nb::arg("f").noconvert(), nb::arg("out").noconvert(), nb::arg("dims"), nb::arg("boxes"),
        nb::arg("num_threads"),
        "Blend the tiles f ((t, z, y, x) sizes dims) into out; boxes: flat (t, z0, y0, x0, z1, "
        "y1, x1) rows, 0-based, end-exclusive, not empty and inside the tiles.");
}

template <typename T>
void bind_class(nb::module_ &m) {
    bind_one<T, T, NONE>(m, "tile_blending_3d_none");
    bind_one<T, T, MAX>(m, "tile_blending_3d_max");
    bind_one<T, float, MEAN>(m, "tile_blending_3d_mean");
    bind_one<T, float, MEDIAN>(m, "tile_blending_3d_median");
}

}  // namespace

PK5D_BINDINGS(bind_tile_blending) {
    bind_class<uint8_t>(m);
    bind_class<uint16_t>(m);
    bind_class<float>(m);
}
