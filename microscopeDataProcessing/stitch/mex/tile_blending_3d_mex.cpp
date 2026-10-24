#include <stdint.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <omp.h>
#include "mex.h"

// nv_block = tile_blending_3d_mex(tim_f_block, mboxes, BlendMethod)
//
// Blends the tiles of a stitching batch for the blend methods 'none', 'max',
// 'mean' and 'median' in one pass, with the results of the MATLAB code in
// processStitchBlock (tile_blending_3d.m). tim_f_block: tiles in the 4th
// dimension (0 = no data), uint8, uint16 or single. mboxes: one row
// [tile, y1, x1, z1, y2, x2, z2] (1-based, inclusive) per region where the
// tile counts as 0 (overlap regions other tiles fill), can be empty. A voxel
// gets the blend of the tiles outside their regions, or, where that is 0, the
// blend of all tiles. The result is of the class of the tiles for 'none' and
// 'max' (it is one of their values), single for 'mean' and 'median'.
//
// The single arithmetic follows MATLAB's: values are summed in tile order,
// the mean of the two middle values of an even count is a + (b - a) / 2, or
// (a + b) / 2 if their signs differ or one is infinite, and no fused
// multiply-adds (-ffp-contract=off).
//
// mex -v CXXOPTIMFLAGS="-O3 -DNDEBUG" CXXFLAGS='$CXXFLAGS -O3 -fopenmp -ffp-contract=off' LDFLAGS='$LDFLAGS -O3 -fopenmp' tile_blending_3d_mex.cpp

enum Method {NONE, MAX, MEAN, MEDIAN};

// a region where a tile counts as 0 (0-based, inclusive)
struct Box {
    uint64_t t, x1, y1, z1, x2, y2, z2;
};

// value of tile j at the voxel: f and m point to the voxel's value and mask
// in tile 0, nXYZ and nX are the strides to the next tile
#define TVAL(j) ((MASK && m[(j) * nX]) ? (T)0 : f[(j) * nXYZ])

// MATLAB's sign
static inline int sgn(const float v) {
    return (v > 0) - (v < 0);
}

// MATLAB's median of n > 0 values without NaN
static inline float medianOf(float* buf, const uint64_t n) {
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

// the last step of the MATLAB code: nv, or nvf where nv is 0, as
// nv .* (1 - z) + nvf .* z
static inline float fillZeros(const float nv, const float nvf) {
    const float z = (nv == 0.0f && nvf != 0.0f) ? 1.0f : 0.0f;
    return nv * (1.0f - z) + nvf * z;
}

// Integer tiles: their values are exact in single, so 'none' is the first
// nonzero tile value and 'max' the largest, outside the regions or else of
// all tiles; 'mean' and 'median' are those of the nonzero values.
template <typename T, bool MASK>
static inline T noneInt(const T* f, const uint8_t* m, const uint64_t nX, const uint64_t nXYZ, const uint64_t nT) {
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
static inline T maxInt(const T* f, const uint8_t* m, const uint64_t nX, const uint64_t nXYZ, const uint64_t nT) {
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
static inline float meanInt(const T* f, const uint8_t* m, const uint64_t nX, const uint64_t nXYZ, const uint64_t nT) {
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
static inline float medianInt(const T* f, const uint8_t* m, const uint64_t nX, const uint64_t nXYZ, const uint64_t nT, float* buf) {
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
static inline float noneSingle(const T* f, const uint8_t* m, const uint64_t nX, const uint64_t nXYZ, const uint64_t nT) {
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
static inline float maxSingle(const T* f, const uint8_t* m, const uint64_t nX, const uint64_t nXYZ, const uint64_t nT) {
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
static inline float meanSingle(const T* f, const uint8_t* m, const uint64_t nX, const uint64_t nXYZ, const uint64_t nT) {
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
static inline float medianSingle(const T* f, const uint8_t* m, const uint64_t nX, const uint64_t nXYZ, const uint64_t nT, float* buf) {
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

// one line of voxels along the first dimension; f and o point to its first
// voxel, mask to the masks of the line (nX per tile)
template <typename T, typename O, Method M, bool MASK>
static void blendLine(const T* f, const uint8_t* mask, O* o, const uint64_t nX, const uint64_t nXYZ, const uint64_t nT, float* buf) {
    const bool isInt = !std::is_floating_point<T>::value;
    for (uint64_t x = 0; x < nX; x++) {
        const T* fx = f + x;
        const uint8_t* mx = MASK ? mask + x : NULL;
        switch (M) {
            case NONE:
                o[x] = isInt ? (O)noneInt<T, MASK>(fx, mx, nX, nXYZ, nT) : (O)noneSingle<T, MASK>(fx, mx, nX, nXYZ, nT);
                break;
            case MAX:
                o[x] = isInt ? (O)maxInt<T, MASK>(fx, mx, nX, nXYZ, nT) : (O)maxSingle<T, MASK>(fx, mx, nX, nXYZ, nT);
                break;
            case MEAN:
                o[x] = isInt ? (O)meanInt<T, MASK>(fx, mx, nX, nXYZ, nT) : (O)meanSingle<T, MASK>(fx, mx, nX, nXYZ, nT);
                break;
            case MEDIAN:
                o[x] = isInt ? (O)medianInt<T, MASK>(fx, mx, nX, nXYZ, nT, buf) : (O)medianSingle<T, MASK>(fx, mx, nX, nXYZ, nT, buf);
                break;
        }
    }
}

template <typename T, typename O, Method M>
static void blendTiles(const T* f, O* o, const uint64_t nX, const uint64_t nY, const uint64_t nZ, const uint64_t nT, const std::vector<Box> &boxes) {
    const uint64_t nXY = nX * nY;
    const uint64_t nXYZ = nXY * nZ;
    const int64_t nLines = (int64_t)(nY * nZ);
    const uint64_t nBoxes = boxes.size();

    #pragma omp parallel
    {
        std::vector<uint8_t> mask(nT * nX);
        std::vector<float> buf(nT);

        #pragma omp for schedule(static)
        for (int64_t l = 0; l < nLines; l++) {
            const uint64_t y = (uint64_t)l % nY;
            const uint64_t z = (uint64_t)l / nY;
            // masks of the regions that cross the line
            bool masked = false;
            for (uint64_t k = 0; k < nBoxes; k++) {
                const Box &b = boxes[k];
                if (y < b.y1 || y > b.y2 || z < b.z1 || z > b.z2) continue;
                if (!masked) {
                    memset(mask.data(), 0, nT * nX);
                    masked = true;
                }
                memset(mask.data() + b.t * nX + b.x1, 1, b.x2 - b.x1 + 1);
            }
            const uint64_t offset = y * nX + z * nXY;
            if (masked) {
                blendLine<T, O, M, true>(f + offset, mask.data(), o + offset, nX, nXYZ, nT, buf.data());
            }
            else {
                blendLine<T, O, M, false>(f + offset, NULL, o + offset, nX, nXYZ, nT, buf.data());
            }
        }
    }
}

template <typename T>
static void blendTilesMethod(const mxArray* fArr, mxArray* oArr, const Method method, const uint64_t* dims, const std::vector<Box> &boxes) {
    const T* f = (const T*)mxGetData(fArr);
    void* o = mxGetData(oArr);
    switch (method) {
        case NONE:
            blendTiles<T, T, NONE>(f, (T*)o, dims[0], dims[1], dims[2], dims[3], boxes);
            break;
        case MAX:
            blendTiles<T, T, MAX>(f, (T*)o, dims[0], dims[1], dims[2], dims[3], boxes);
            break;
        case MEAN:
            blendTiles<T, float, MEAN>(f, (float*)o, dims[0], dims[1], dims[2], dims[3], boxes);
            break;
        case MEDIAN:
            blendTiles<T, float, MEDIAN>(f, (float*)o, dims[0], dims[1], dims[2], dims[3], boxes);
            break;
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
#ifdef __APPLE__
    // keep the mex file loaded until MATLAB exits: GCC's OpenMP runtime leaves
    // thread exit handlers in it, so unloading it crashes MATLAB on exit
    if (!mexIsLocked()) {
        mexLock();
    }
#endif
    if (nrhs != 3) mexErrMsgIdAndTxt("tile_blending:inputError", "Number of input arguments must be 3");
    if (nlhs > 1) mexErrMsgIdAndTxt("tile_blending:outputError", "Number of output arguments must be 1");

    const mxArray* fArr = prhs[0];
    const mxClassID mDType = mxGetClassID(fArr);
    if ((mDType != mxUINT8_CLASS && mDType != mxUINT16_CLASS && mDType != mxSINGLE_CLASS) || mxIsComplex(fArr) || mxIsSparse(fArr)) {
        mexErrMsgIdAndTxt("tile_blending:dataTypeError", "Data type not supported (uint8, uint16 or single)");
    }
    const uint64_t ndim = (uint64_t)mxGetNumberOfDimensions(fArr);
    if (ndim > 4) mexErrMsgIdAndTxt("tile_blending:inputError", "The tiles must be a 4D array (tiles in the 4th dimension)");
    const mwSize* dimsA = mxGetDimensions(fArr);
    uint64_t dims[4] = {1, 1, 1, 1};
    for (uint64_t i = 0; i < ndim; i++) dims[i] = (uint64_t)dimsA[i];
    if (dims[3] == 0) mexErrMsgIdAndTxt("tile_blending:inputError", "There are no tiles");

    if (!mxIsChar(prhs[2])) mexErrMsgIdAndTxt("tile_blending:inputError", "The blend method must be a char array");
    char methodStr[16];
    Method method = NONE;
    if (mxGetString(prhs[2], methodStr, sizeof(methodStr)) != 0) {
        mexErrMsgIdAndTxt("tile_blending:inputError", "Unsupported blending method");
    }
    if (strcmp(methodStr, "none") == 0) method = NONE;
    else if (strcmp(methodStr, "max") == 0) method = MAX;
    else if (strcmp(methodStr, "mean") == 0) method = MEAN;
    else if (strcmp(methodStr, "median") == 0) method = MEDIAN;
    else mexErrMsgIdAndTxt("tile_blending:inputError", "Unsupported blending method: %s", methodStr);

    // regions where tiles count as 0; empty ones are skipped
    const mxArray* bArr = prhs[1];
    std::vector<Box> boxes;
    if (!mxIsEmpty(bArr)) {
        if (!mxIsDouble(bArr) || mxIsComplex(bArr) || mxGetNumberOfDimensions(bArr) != 2 || mxGetN(bArr) != 7) {
            mexErrMsgIdAndTxt("tile_blending:inputError", "The regions must be a K x 7 double array [tile, y1, x1, z1, y2, x2, z2]");
        }
        const uint64_t K = (uint64_t)mxGetM(bArr);
        const double* b = mxGetPr(bArr);
        for (uint64_t k = 0; k < K; k++) {
            double v[7];
            for (uint64_t c = 0; c < 7; c++) {
                v[c] = b[k + c * K];
                if (!(v[c] == std::floor(v[c])) || std::fabs(v[c]) > 9.0e15) {
                    mexErrMsgIdAndTxt("tile_blending:inputError", "The regions must have integer coordinates");
                }
            }
            if (v[4] < v[1] || v[5] < v[2] || v[6] < v[3]) continue;
            if (v[0] < 1 || v[0] > (double)dims[3] || v[1] < 1 || v[2] < 1 || v[3] < 1 ||
                    v[4] > (double)dims[0] || v[5] > (double)dims[1] || v[6] > (double)dims[2]) {
                mexErrMsgIdAndTxt("tile_blending:inputError", "Region %d is outside the tiles", (int)(k + 1));
            }
            Box box;
            box.t = (uint64_t)v[0] - 1;
            box.x1 = (uint64_t)v[1] - 1;
            box.y1 = (uint64_t)v[2] - 1;
            box.z1 = (uint64_t)v[3] - 1;
            box.x2 = (uint64_t)v[4] - 1;
            box.y2 = (uint64_t)v[5] - 1;
            box.z2 = (uint64_t)v[6] - 1;
            boxes.push_back(box);
        }
    }

    const mxClassID oDType = (method == NONE || method == MAX) ? mDType : mxSINGLE_CLASS;
    mwSize odims[3] = {(mwSize)dims[0], (mwSize)dims[1], (mwSize)dims[2]};
    plhs[0] = mxCreateUninitNumericArray(3, odims, oDType, mxREAL);

    if (mDType == mxUINT8_CLASS) {
        blendTilesMethod<uint8_t>(fArr, plhs[0], method, dims, boxes);
    }
    else if (mDType == mxUINT16_CLASS) {
        blendTilesMethod<uint16_t>(fArr, plhs[0], method, dims, boxes);
    }
    else {
        blendTilesMethod<float>(fArr, plhs[0], method, dims, boxes);
    }
}
