// The acceleration step of decon_lucy_function's Richardson-Lucy iterations on
// the CPU in two passes over the data, and sums in double, with the sums in an
// order that doesn't depend on the number of threads:
//
// [lambda, J_4, s] = decon_lucy_accel_mex('lambda', J_2, Y, J_4, nThreads)
//   s      : [sum((J_2(:) - Y(:)) .* J_4, 'double'), sum(J_4 .* J_4, 'double')]
//   lambda : max(min(s(1) / (s(2) + eps), 1), 0)
//   J_4    : the new J_4, J_2(:) - Y(:) (a column)
//
// Y = decon_lucy_accel_mex('predict', J_2, J_3, lambda, nThreads)
//   max(J_2 + lambda * (J_2 - J_3), 0), the size of J_2
//
// s = decon_lucy_accel_mex('sum', x, nThreads)
//   sum(x(:), 'double')
//
//   J_2, Y, J_3, J_4 : real single or double arrays of one class with the same number of elements
//   x                : real single or double array
//   lambda           : real double scalar
//   nThreads         : threads
//
// The order of the sums is MATLAB's at one thread for arrays of 89000 or more
// elements: blocks of 4096 elements in index order (the last one shorter),
// each summed in four lanes (element t of a block in lane mod(t, 4), in index
// order), the lanes added as ((l1 + l2) + l3) + l4 and the block sums in index
// order. Each block is summed by one thread, so the sums don't depend on the
// number of threads (MATLAB splits the blocks among its threads, which changes
// the last bits). The accumulators start at -0, so a sum of -0s is -0 as in
// MATLAB; an empty sum is 0. (Below 89000 elements MATLAB sums in one block of
// four lanes, below 44 in index order, on one thread at any thread count;
// decon_sum_double and decon_lucy_accel use MATLAB's sum there.)
//
// The element-wise arithmetic is MATLAB's: differences and products in the
// class of the data (single products converted to double for the sums), a
// double scalar times a single array computed in double and rounded to single,
// max(v, 0) is 0 for NaN and keeps -0 (MATLAB returns the first argument on a
// tie), and min(NaN, 1) is 1.
//
// mex -v CXXOPTIMFLAGS="-O3 -DNDEBUG" CXXFLAGS='$CXXFLAGS -O3 -fopenmp -ffp-contract=off' LDFLAGS='$LDFLAGS -O3 -fopenmp' decon_lucy_accel_mex.cpp
// (no fused multiply-adds, so that products are rounded as MATLAB's)

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <float.h>
#include <cmath>
#include <vector>
#include <omp.h>
#include "mex.h"

typedef ptrdiff_t idx_t;

// elements per block of the sums
static const idx_t kBlock = 4096;

// the lanes of a block: ((l1 + l2) + l3) + l4
static inline double add_lanes(const double* a) {
    return ((a[0] + a[1]) + a[2]) + a[3];
}

// the block sums in index order
static double add_blocks(const std::vector<double>& b) {
    if (b.empty()) {
        return 0.0;
    }
    double s = -0.0;
    for (size_t k = 0; k < b.size(); k++) {
        s = s + b[k];
    }
    return s;
}

static inline idx_t n_blocks(const idx_t n) {
    return (n + kBlock - 1) / kBlock;
}

// sum(x(:), 'double')
template <typename T>
static double sum_blocks(const T* x, const idx_t n, const int nt) {
    const idx_t nb = n_blocks(n);
    std::vector<double> bs(nb);
    #pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const T* p = x + b * kBlock;
        const idx_t m = b < nb - 1 ? kBlock : n - b * kBlock;
        double a[4] = {-0.0, -0.0, -0.0, -0.0};
        idx_t t = 0;
        for (; t + 4 <= m; t += 4) {
            a[0] = a[0] + (double)p[t];
            a[1] = a[1] + (double)p[t + 1];
            a[2] = a[2] + (double)p[t + 2];
            a[3] = a[3] + (double)p[t + 3];
        }
        for (; t < m; t++) {
            a[t & 3] = a[t & 3] + (double)p[t];
        }
        bs[b] = add_lanes(a);
    }
    return add_blocks(bs);
}

// the sums of (J_2 - Y) .* J_4 and J_4 .* J_4 (products in T, summed in
// double), and the new J_4 = J_2 - Y in out
template <typename T>
static void lambda_sums(const T* __restrict__ j2, const T* __restrict__ y, const T* __restrict__ j4, T* __restrict__ out,
        const idx_t n, const int nt, double* s) {
    const idx_t nb = n_blocks(n);
    std::vector<double> b1(nb), b2(nb);
    #pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t b = 0; b < nb; b++) {
        const idx_t i0 = b * kBlock;
        const idx_t m = b < nb - 1 ? kBlock : n - i0;
        const T* __restrict__ pj2 = j2 + i0;
        const T* __restrict__ py = y + i0;
        const T* __restrict__ pj4 = j4 + i0;
        T* __restrict__ po = out + i0;
        double a1[4] = {-0.0, -0.0, -0.0, -0.0};
        double a2[4] = {-0.0, -0.0, -0.0, -0.0};
        idx_t t = 0;
        for (; t + 4 <= m; t += 4) {
            for (int l = 0; l < 4; l++) {
                const T d = pj2[t + l] - py[t + l];
                const T v = pj4[t + l];
                const T p1 = d * v;
                const T p2 = v * v;
                a1[l] = a1[l] + (double)p1;
                a2[l] = a2[l] + (double)p2;
                po[t + l] = d;
            }
        }
        for (; t < m; t++) {
            const T d = pj2[t] - py[t];
            const T v = pj4[t];
            const T p1 = d * v;
            const T p2 = v * v;
            a1[t & 3] = a1[t & 3] + (double)p1;
            a2[t & 3] = a2[t & 3] + (double)p2;
            po[t] = d;
        }
        b1[b] = add_lanes(a1);
        b2[b] = add_lanes(a2);
    }
    s[0] = add_blocks(b1);
    s[1] = add_blocks(b2);
}

// max(J_2 + lambda * (J_2 - J_3), 0): the product in double, rounded to T
template <typename T>
static void predict(const T* __restrict__ j2, const T* __restrict__ j3, const double lambda, T* __restrict__ y, const idx_t n,
        const int nt) {
    #pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t i = 0; i < n; i++) {
        const T d = j2[i] - j3[i];
        const T t = (T)(lambda * (double)d);
        const T u = j2[i] + t;
        // max(u, 0): NaN is 0, -0 is kept
        y[i] = u >= (T)0 ? u : (T)0;
    }
}

// lambda from the sums: max(min(s(1) / (s(2) + eps), 1), 0), with MATLAB's
// min (NaN gives 1) and max (-0 is kept)
static double clamped_lambda(const double* s) {
    const double v = s[0] / (s[1] + DBL_EPSILON);
    const double m = std::isnan(v) ? 1.0 : (v < 1.0 ? v : 1.0);
    return m >= 0.0 ? m : 0.0;
}

static void check_data(const mxArray* a, const char* name) {
    if (!(mxIsSingle(a) || mxIsDouble(a)) || mxIsComplex(a) || mxIsSparse(a)) {
        mexErrMsgIdAndTxt("decon_lucy_accel_mex:inputError", "%s must be a real single or double array.", name);
    }
}

static void check_like(const mxArray* a, const mxArray* ref, const char* name) {
    check_data(a, name);
    if (mxGetClassID(a) != mxGetClassID(ref) || mxGetNumberOfElements(a) != mxGetNumberOfElements(ref)) {
        mexErrMsgIdAndTxt("decon_lucy_accel_mex:inputError", "%s must have the class and number of elements of J_2.", name);
    }
}

static int threads_arg(const mxArray* a) {
    if (!mxIsNumeric(a) || mxGetNumberOfElements(a) != 1 || mxIsComplex(a)) {
        mexErrMsgIdAndTxt("decon_lucy_accel_mex:inputError", "nThreads must be a real numeric scalar.");
    }
    const double v = mxGetScalar(a);
    if (!(v >= 1)) {
        return 1;
    }
    return v > 4096 ? 4096 : (int)v;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
#ifdef __APPLE__
    // keep the mex file loaded until MATLAB exits: GCC's OpenMP runtime leaves
    // thread exit handlers in it, so unloading it crashes MATLAB on exit
    if (!mexIsLocked()) {
        mexLock();
    }
#endif
    if (nrhs < 1 || !mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("decon_lucy_accel_mex:inputError", "The first input must be 'lambda', 'predict' or 'sum'.");
    }
    char method[16];
    mxGetString(prhs[0], method, sizeof(method));

    if (strcmp(method, "lambda") == 0) {
        if (nrhs != 5 || nlhs > 3) {
            mexErrMsgIdAndTxt("decon_lucy_accel_mex:inputError", "Usage: [lambda, J_4, s] = decon_lucy_accel_mex('lambda', J_2, Y, J_4, nThreads)");
        }
        const mxArray* j2 = prhs[1];
        check_data(j2, "J_2");
        check_like(prhs[2], j2, "Y");
        check_like(prhs[3], j2, "J_4");
        const int nt = threads_arg(prhs[4]);
        const idx_t n = (idx_t)mxGetNumberOfElements(j2);
        const bool isSingle = mxIsSingle(j2);
        mxArray* out = mxCreateUninitNumericMatrix((mwSize)n, 1, mxGetClassID(j2), mxREAL);
        double s[2];
        if (isSingle) {
            lambda_sums((const float*)mxGetData(j2), (const float*)mxGetData(prhs[2]), (const float*)mxGetData(prhs[3]),
                (float*)mxGetData(out), n, nt, s);
        } else {
            lambda_sums((const double*)mxGetData(j2), (const double*)mxGetData(prhs[2]), (const double*)mxGetData(prhs[3]),
                (double*)mxGetData(out), n, nt, s);
        }
        plhs[0] = mxCreateDoubleScalar(clamped_lambda(s));
        if (nlhs > 1) {
            plhs[1] = out;
        } else {
            mxDestroyArray(out);
        }
        if (nlhs > 2) {
            plhs[2] = mxCreateDoubleMatrix(1, 2, mxREAL);
            double* so = (double*)mxGetData(plhs[2]);
            so[0] = s[0];
            so[1] = s[1];
        }
    } else if (strcmp(method, "predict") == 0) {
        if (nrhs != 5 || nlhs > 1) {
            mexErrMsgIdAndTxt("decon_lucy_accel_mex:inputError", "Usage: Y = decon_lucy_accel_mex('predict', J_2, J_3, lambda, nThreads)");
        }
        const mxArray* j2 = prhs[1];
        check_data(j2, "J_2");
        check_like(prhs[2], j2, "J_3");
        const mxArray* la = prhs[3];
        if (!mxIsDouble(la) || mxGetNumberOfElements(la) != 1 || mxIsComplex(la) || mxIsSparse(la)) {
            mexErrMsgIdAndTxt("decon_lucy_accel_mex:inputError", "lambda must be a real double scalar.");
        }
        const double lambda = mxGetScalar(la);
        const int nt = threads_arg(prhs[4]);
        const idx_t n = (idx_t)mxGetNumberOfElements(j2);
        const mwSize* dimsA = mxGetDimensions(j2);
        std::vector<mwSize> dims(dimsA, dimsA + mxGetNumberOfDimensions(j2));
        plhs[0] = mxCreateUninitNumericArray(dims.size(), dims.data(), mxGetClassID(j2), mxREAL);
        if (mxIsSingle(j2)) {
            predict((const float*)mxGetData(j2), (const float*)mxGetData(prhs[2]), lambda, (float*)mxGetData(plhs[0]), n, nt);
        } else {
            predict((const double*)mxGetData(j2), (const double*)mxGetData(prhs[2]), lambda, (double*)mxGetData(plhs[0]), n, nt);
        }
    } else if (strcmp(method, "sum") == 0) {
        if (nrhs != 3 || nlhs > 1) {
            mexErrMsgIdAndTxt("decon_lucy_accel_mex:inputError", "Usage: s = decon_lucy_accel_mex('sum', x, nThreads)");
        }
        const mxArray* x = prhs[1];
        check_data(x, "x");
        const int nt = threads_arg(prhs[2]);
        const idx_t n = (idx_t)mxGetNumberOfElements(x);
        double s;
        if (mxIsSingle(x)) {
            s = sum_blocks((const float*)mxGetData(x), n, nt);
        } else {
            s = sum_blocks((const double*)mxGetData(x), n, nt);
        }
        plhs[0] = mxCreateDoubleScalar(s);
    } else {
        mexErrMsgIdAndTxt("decon_lucy_accel_mex:inputError", "Unknown method %s, use 'lambda', 'predict' or 'sum'.", method);
    }
}
