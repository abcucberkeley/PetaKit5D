#include <stdint.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "mex.h"

// im = unmix_channels_sum_mex(ims, unmixFactors, constBackground, constOffset, dtype)
// im = unmix_channels_sum_mex(ims, unmixFactors, constBackground, constOffset, dtype, im0)
//
// The unmixed image of the channel images ims (a cell array of real arrays of
// the same size: uint8, uint16, single or double) in one pass, with the
// results of the MATLAB code in unmix_channels_sum.m:
//   im = im0 (single, or 0)
//   for each channel c: im = im + (single(ims{c}) - constBackground(c)) * unmixFactors(c)
//   im = max(0, im), then im + constOffset, then cast(im, dtype)
// unmixFactors and constBackground: one value per channel, constOffset: a
// scalar (doubles; a background or offset of 0 is not applied, which gives
// the same values). dtype: 'uint8', 'uint16', 'single' or 'double', or ''
// for the partial sum (single, not clipped, without the offset). im0: a
// partial sum to add the channels to (single, of their size, or [] for 0);
// ims can be empty if it is given.
//
// The arithmetic follows MATLAB's: an operation of a single and a double is
// computed in double and rounded to single, two singles are added in single,
// max(0, NaN) is 0, and cast rounds half away from zero, saturates and maps
// NaN to 0. No fused multiply-adds (-ffp-contract=off).
//
// mex -v CXXOPTIMFLAGS="-O3 -DNDEBUG" CXXFLAGS='$CXXFLAGS -O3 -fopenmp -ffp-contract=off' LDFLAGS='$LDFLAGS -O3 -fopenmp' unmix_channels_sum_mex.cpp

// elements per block of work (the sum of a block stays in the cache)
static const uint64_t BLOCK = 8192;

struct Channel {
    const void* data;
    mxClassID cls;
    double factor;
    double background;
};

// acc = acc + (single(x) - b) * f
template <typename T>
static void addChannel(float* acc, const T* x, const uint64_t n, const double b, const double f) {
    if (b != 0.0) {
        for (uint64_t i = 0; i < n; i++) {
            const float v = (float)((double)(float)x[i] - b);
            acc[i] = acc[i] + (float)((double)v * f);
        }
    }
    else {
        for (uint64_t i = 0; i < n; i++) {
            acc[i] = acc[i] + (float)((double)(float)x[i] * f);
        }
    }
}

// MATLAB's cast of a single value
template <typename O>
static inline O castTo(const float v);

template <>
inline uint16_t castTo<uint16_t>(const float v) {
    // NaN, negative values and values below 0.5 are 0
    if (!(v >= 0.5f)) return 0;
    if (v >= 65535.0f) return 65535;
    return (uint16_t)((double)v + 0.5);
}

template <>
inline uint8_t castTo<uint8_t>(const float v) {
    if (!(v >= 0.5f)) return 0;
    if (v >= 255.0f) return 255;
    return (uint8_t)((double)v + 0.5);
}

template <>
inline float castTo<float>(const float v) {
    return v;
}

template <>
inline double castTo<double>(const float v) {
    return (double)v;
}

// the sum of the channels (from im0, or 0), then clipped, offset and cast
// (finish), or as it is
template <typename O>
static void unmixChannels(const std::vector<Channel> &channels, const float* im0, O* out, const uint64_t n,
        const double offset, const bool finish) {
    const int64_t nBlocks = (int64_t)((n + BLOCK - 1) / BLOCK);
    #pragma omp parallel
    {
        std::vector<float> acc(BLOCK);
        #pragma omp for schedule(static)
        for (int64_t k = 0; k < nBlocks; k++) {
            const uint64_t s = (uint64_t)k * BLOCK;
            const uint64_t m = std::min(BLOCK, n - s);
            float* a = acc.data();
            if (im0) {
                std::copy(im0 + s, im0 + s + m, a);
            }
            else {
                std::fill(a, a + m, 0.0f);
            }
            for (const Channel &c : channels) {
                switch (c.cls) {
                    case mxUINT8_CLASS:
                        addChannel<uint8_t>(a, (const uint8_t*)c.data + s, m, c.background, c.factor);
                        break;
                    case mxUINT16_CLASS:
                        addChannel<uint16_t>(a, (const uint16_t*)c.data + s, m, c.background, c.factor);
                        break;
                    case mxSINGLE_CLASS:
                        addChannel<float>(a, (const float*)c.data + s, m, c.background, c.factor);
                        break;
                    default:
                        addChannel<double>(a, (const double*)c.data + s, m, c.background, c.factor);
                        break;
                }
            }
            O* o = out + s;
            if (!finish) {
                for (uint64_t i = 0; i < m; i++) {
                    o[i] = castTo<O>(a[i]);
                }
            }
            else if (offset != 0.0) {
                for (uint64_t i = 0; i < m; i++) {
                    // max(0, v): NaN and -0 become 0
                    const float v = a[i] > 0.0f ? a[i] : 0.0f;
                    o[i] = castTo<O>((float)((double)v + offset));
                }
            }
            else {
                for (uint64_t i = 0; i < m; i++) {
                    o[i] = castTo<O>(a[i] > 0.0f ? a[i] : 0.0f);
                }
            }
        }
    }
}

// the values of a real double vector with n elements
static const double* doubleVector(const mxArray* a, const uint64_t n, const char* name) {
    if (!mxIsDouble(a) || mxIsComplex(a) || mxIsSparse(a) || (uint64_t)mxGetNumberOfElements(a) != n) {
        mexErrMsgIdAndTxt("unmix_channels:inputError", "%s must be a real double array with %d elements", name, (int)n);
    }
    return mxGetPr(a);
}

// whether two arrays have the same size
static bool sameSize(const mxArray* a, const mxArray* b) {
    return mxGetNumberOfDimensions(a) == mxGetNumberOfDimensions(b) &&
        memcmp(mxGetDimensions(a), mxGetDimensions(b), mxGetNumberOfDimensions(a) * sizeof(mwSize)) == 0;
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
    if (nrhs != 5 && nrhs != 6) mexErrMsgIdAndTxt("unmix_channels:inputError", "Number of input arguments must be 5 or 6");
    if (nlhs > 1) mexErrMsgIdAndTxt("unmix_channels:outputError", "Number of output arguments must be 1");

    // the partial sum to start from
    const mxArray* im0Arr = (nrhs == 6 && !mxIsEmpty(prhs[5])) ? prhs[5] : NULL;
    if (im0Arr && (!mxIsSingle(im0Arr) || mxIsComplex(im0Arr) || mxIsSparse(im0Arr))) {
        mexErrMsgIdAndTxt("unmix_channels:dataTypeError", "The partial sum must be a real single array");
    }

    const mxArray* imsArr = prhs[0];
    if (!mxIsCell(imsArr) || (mxIsEmpty(imsArr) && !im0Arr)) {
        mexErrMsgIdAndTxt("unmix_channels:inputError", "The channel images must be a non-empty cell array (or empty with a partial sum)");
    }
    const uint64_t nc = (uint64_t)mxGetNumberOfElements(imsArr);
    const double* factors = doubleVector(prhs[1], nc, "The unmixing factors");
    const double* backgrounds = doubleVector(prhs[2], nc, "The backgrounds");
    const double offset = *doubleVector(prhs[3], 1, "The offset");

    // the result's class, or the partial sum (single)
    if (!mxIsChar(prhs[4])) mexErrMsgIdAndTxt("unmix_channels:inputError", "The output data type must be a char array");
    char dtypeStr[16];
    if (mxGetString(prhs[4], dtypeStr, sizeof(dtypeStr)) != 0) {
        mexErrMsgIdAndTxt("unmix_channels:dataTypeError", "Unsupported output data type");
    }
    const bool finish = dtypeStr[0] != '\0';
    mxClassID oDType = mxSINGLE_CLASS;
    if (!finish || strcmp(dtypeStr, "single") == 0) oDType = mxSINGLE_CLASS;
    else if (strcmp(dtypeStr, "uint8") == 0) oDType = mxUINT8_CLASS;
    else if (strcmp(dtypeStr, "uint16") == 0) oDType = mxUINT16_CLASS;
    else if (strcmp(dtypeStr, "double") == 0) oDType = mxDOUBLE_CLASS;
    else mexErrMsgIdAndTxt("unmix_channels:dataTypeError", "Unsupported output data type: %s (uint8, uint16, single or double)", dtypeStr);

    // the channels and the partial sum, all of the same size
    const mxArray* first = nc > 0 ? mxGetCell(imsArr, 0) : im0Arr;
    if (first == NULL) mexErrMsgIdAndTxt("unmix_channels:inputError", "Channel 1 has no image");
    if (im0Arr && !sameSize(im0Arr, first)) {
        mexErrMsgIdAndTxt("unmix_channels:inputError", "The partial sum has a different size than channel 1");
    }
    std::vector<Channel> channels(nc);
    for (uint64_t c = 0; c < nc; c++) {
        const mxArray* im = mxGetCell(imsArr, c);
        if (im == NULL) mexErrMsgIdAndTxt("unmix_channels:inputError", "Channel %d has no image", (int)(c + 1));
        const mxClassID cls = mxGetClassID(im);
        if ((cls != mxUINT8_CLASS && cls != mxUINT16_CLASS && cls != mxSINGLE_CLASS && cls != mxDOUBLE_CLASS) ||
                mxIsComplex(im) || mxIsSparse(im)) {
            mexErrMsgIdAndTxt("unmix_channels:dataTypeError", "Channel %d: data type not supported (uint8, uint16, single or double)", (int)(c + 1));
        }
        if (!sameSize(im, first)) {
            mexErrMsgIdAndTxt("unmix_channels:inputError", "Channel %d has a different size than channel 1", (int)(c + 1));
        }
        channels[c].data = mxGetData(im);
        channels[c].cls = cls;
        channels[c].factor = factors[c];
        channels[c].background = backgrounds[c];
    }

    const uint64_t n = (uint64_t)mxGetNumberOfElements(first);
    const mwSize* dimsA = mxGetDimensions(first);
    std::vector<mwSize> dims(dimsA, dimsA + mxGetNumberOfDimensions(first));
    plhs[0] = mxCreateUninitNumericArray(dims.size(), dims.data(), oDType, mxREAL);
    if (n == 0) return;

    const float* im0 = im0Arr ? (const float*)mxGetData(im0Arr) : NULL;
    switch (oDType) {
        case mxUINT8_CLASS:
            unmixChannels<uint8_t>(channels, im0, (uint8_t*)mxGetData(plhs[0]), n, offset, finish);
            break;
        case mxUINT16_CLASS:
            unmixChannels<uint16_t>(channels, im0, (uint16_t*)mxGetData(plhs[0]), n, offset, finish);
            break;
        case mxSINGLE_CLASS:
            unmixChannels<float>(channels, im0, (float*)mxGetData(plhs[0]), n, offset, finish);
            break;
        default:
            unmixChannels<double>(channels, im0, (double*)mxGetData(plhs[0]), n, offset, finish);
            break;
    }
}
