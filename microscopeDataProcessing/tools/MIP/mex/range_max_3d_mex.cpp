#include <stdint.h>
#include <math.h>
#include <vector>
#include <omp.h>
#include "mex.h"

// out = range_max_3d_mex(im, dim, starts)
//
// Max of a 3D array along dimension dim (1, 2 or 3) over consecutive ranges
// that start at starts (1-based and increasing, starting with 1; a range ends
// where the next one starts, the last one at the end of the array). out has the
// size of im with numel(starts) along dim, and range r is
// max(im(..., starts(r) : starts(r + 1) - 1, ...), [], dim). As MATLAB's max,
// NaN are ignored (a range of NaN gives its first NaN), and of equal values
// (e.g. 0 and -0) the first one is kept. Types: uint8, uint16, single, double.
//
// mex -v CXXOPTIMFLAGS="-O3 -DNDEBUG" CXXFLAGS='$CXXFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' range_max_3d_mex.cpp

template <typename T>
static inline T maxUpdate(const T m, const T v){
    return v > m ? v : m;
}

// a NaN maximum is replaced by the first number
static inline float maxUpdate(const float m, const float v){
    return (v > m || (m != m && v == v)) ? v : m;
}

static inline double maxUpdate(const double m, const double v){
    return (v > m || (m != m && v == v)) ? v : m;
}

template <typename T>
static void rangeMax3D(const T* const im, T* const out, const uint64_t ny, const uint64_t nx, const uint64_t nz,
                       const int dim, const std::vector<uint64_t> &starts, const std::vector<uint64_t> &ends){
    const int64_t nr = (int64_t)starts.size();
    if(dim == 1){
        // each column (x, z) is contiguous
        const int64_t ncol = (int64_t)(nx * nz);
        #pragma omp parallel for schedule(static)
        for(int64_t c = 0; c < ncol; c++){
            const T* col = im + c * ny;
            T* o = out + c * nr;
            for(int64_t r = 0; r < nr; r++){
                T m = col[starts[r]];
                for(uint64_t i = starts[r] + 1; i < ends[r]; i++){
                    m = maxUpdate(m, col[i]);
                }
                o[r] = m;
            }
        }
    }
    else if(dim == 2){
        // per z slice and range, the columns of the range merged into one, in
        // parts of rows
        const uint64_t chunk = 4096;
        const uint64_t nyc = (ny + chunk - 1) / chunk;
        const int64_t ntask = (int64_t)(nz * nr * nyc);
        #pragma omp parallel for schedule(static)
        for(int64_t t = 0; t < ntask; t++){
            const uint64_t yc = t % nyc;
            const uint64_t r = (t / nyc) % nr;
            const uint64_t z = t / (nyc * nr);
            const uint64_t y0 = yc * chunk;
            const uint64_t y1 = y0 + chunk < ny ? y0 + chunk : ny;
            T* o = out + (z * nr + r) * ny;
            const T* p = im + (z * nx + starts[r]) * ny;
            for(uint64_t y = y0; y < y1; y++){
                o[y] = p[y];
            }
            for(uint64_t x = starts[r] + 1; x < ends[r]; x++){
                p = im + (z * nx + x) * ny;
                for(uint64_t y = y0; y < y1; y++){
                    o[y] = maxUpdate(o[y], p[y]);
                }
            }
        }
    }
    else{
        // per range, the z slices of the range merged into one, in parts
        const uint64_t nxy = ny * nx;
        const uint64_t chunk = 16384;
        const uint64_t nc = (nxy + chunk - 1) / chunk;
        const int64_t ntask = (int64_t)(nr * nc);
        #pragma omp parallel for schedule(static)
        for(int64_t t = 0; t < ntask; t++){
            const uint64_t ci = t % nc;
            const uint64_t r = t / nc;
            const uint64_t i0 = ci * chunk;
            const uint64_t i1 = i0 + chunk < nxy ? i0 + chunk : nxy;
            T* o = out + r * nxy;
            const T* p = im + starts[r] * nxy;
            for(uint64_t i = i0; i < i1; i++){
                o[i] = p[i];
            }
            for(uint64_t z = starts[r] + 1; z < ends[r]; z++){
                p = im + z * nxy;
                for(uint64_t i = i0; i < i1; i++){
                    o[i] = maxUpdate(o[i], p[i]);
                }
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#ifdef __APPLE__
    // keep the mex file loaded until MATLAB exits: GCC's OpenMP runtime leaves
    // thread exit handlers in it, so unloading it crashes MATLAB on exit
    if (!mexIsLocked()) {
        mexLock();
    }
#endif
    if(nrhs != 3) mexErrMsgIdAndTxt("range_max_3d:inputError", "Three inputs required: im, dim, starts.");
    if(nlhs > 1) mexErrMsgIdAndTxt("range_max_3d:outputError", "One output only.");

    const mxClassID classID = mxGetClassID(prhs[0]);
    if(mxIsComplex(prhs[0]) || (classID != mxUINT8_CLASS && classID != mxUINT16_CLASS
       && classID != mxSINGLE_CLASS && classID != mxDOUBLE_CLASS)){
        mexErrMsgIdAndTxt("range_max_3d:inputError", "im must be a real uint8, uint16, single or double array.");
    }
    const uint64_t numDim = (uint64_t)mxGetNumberOfDimensions(prhs[0]);
    if(numDim > 3) mexErrMsgIdAndTxt("range_max_3d:inputError", "im must have at most 3 dimensions.");
    const mwSize* dims = mxGetDimensions(prhs[0]);
    uint64_t sz[3] = {(uint64_t)dims[0], (uint64_t)dims[1], numDim > 2 ? (uint64_t)dims[2] : 1};
    if(sz[0] * sz[1] * sz[2] == 0) mexErrMsgIdAndTxt("range_max_3d:inputError", "im must not be empty.");

    if(!mxIsNumeric(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1) mexErrMsgIdAndTxt("range_max_3d:inputError", "dim must be 1, 2 or 3.");
    const double dimVal = mxGetScalar(prhs[1]);
    if(dimVal != 1 && dimVal != 2 && dimVal != 3) mexErrMsgIdAndTxt("range_max_3d:inputError", "dim must be 1, 2 or 3.");
    const int dim = (int)dimVal;
    const uint64_t n = sz[dim - 1];

    if(!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxGetNumberOfElements(prhs[2]) == 0){
        mexErrMsgIdAndTxt("range_max_3d:inputError", "starts must be a nonempty double vector.");
    }
    const uint64_t nr = (uint64_t)mxGetNumberOfElements(prhs[2]);
    const double* s = mxGetPr(prhs[2]);
    std::vector<uint64_t> starts(nr), ends(nr);
    for(uint64_t r = 0; r < nr; r++){
        const double prev = r == 0 ? 0 : s[r - 1];
        if(s[r] <= prev || s[r] > (double)n || s[r] != floor(s[r]) || (r == 0 && s[r] != 1)){
            mexErrMsgIdAndTxt("range_max_3d:inputError", "starts must be increasing integers from 1 to at most the size of im along dim.");
        }
        starts[r] = (uint64_t)s[r] - 1;
    }
    for(uint64_t r = 0; r < nr; r++){
        ends[r] = r + 1 < nr ? starts[r + 1] : n;
    }

    mwSize outDims[3] = {(mwSize)sz[0], (mwSize)sz[1], (mwSize)sz[2]};
    outDims[dim - 1] = (mwSize)nr;
    plhs[0] = mxCreateUninitNumericArray(3, outDims, classID, mxREAL);

    switch(classID){
        case mxUINT8_CLASS:
            rangeMax3D<uint8_t>((const uint8_t*)mxGetData(prhs[0]), (uint8_t*)mxGetData(plhs[0]), sz[0], sz[1], sz[2], dim, starts, ends);
            break;
        case mxUINT16_CLASS:
            rangeMax3D<uint16_t>((const uint16_t*)mxGetData(prhs[0]), (uint16_t*)mxGetData(plhs[0]), sz[0], sz[1], sz[2], dim, starts, ends);
            break;
        case mxSINGLE_CLASS:
            rangeMax3D<float>((const float*)mxGetData(prhs[0]), (float*)mxGetData(plhs[0]), sz[0], sz[1], sz[2], dim, starts, ends);
            break;
        default:
            rangeMax3D<double>((const double*)mxGetData(prhs[0]), (double*)mxGetData(plhs[0]), sz[0], sz[1], sz[2], dim, starts, ends);
            break;
    }
}
