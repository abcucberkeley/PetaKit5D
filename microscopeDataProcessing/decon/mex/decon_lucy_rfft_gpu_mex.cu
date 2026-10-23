// GPU version of decon_lucy_rfft_mex: the Richardson-Lucy iterations of
// decon_lucy_omw_function ('omw') and the FFT steps of decon_lucy_function
// ('rlstep') with cuFFT real-to-complex FFTs (half spectra, half OTFs) and
// fused kernels for the element-wise steps, on gpuArrays. The results differ
// from MATLAB's gpuArray FFTs in the last bits (other cuFFT plans).
//
// H = decon_lucy_rfft_gpu_mex('otf', psf, outSize)
//   psf      : single, normalized PSF (gpuArray or host), at most outSize in each dimension
//   outSize  : data size, 2 or 3 elements
//   H        : complex single gpuArray half OTF, size [floor(outSize(1) / 2) + 1, outSize(2), outSize(3)]
//
// [J, t] = decon_lucy_rfft_gpu_mex('omw', I, H_f, H_b, nIter, background)
//   I          : uint16 or single raw data (gpuArray or host, 2D or 3D); the iterations use
//                single(I - uint16(background)) for uint16 and max(I - background, 0) for single,
//                as the GPU code of the decon functions
//   H_f, H_b   : half OTF gpuArrays from 'otf'
//   J          : single gpuArray result after nIter iterations
//
// [J, t] = decon_lucy_rfft_gpu_mex('rlstep', Y, I, H)
//   one update of decon_lucy_function: max(Y .* real(ifftn(conj(H) .* fftn(I ./ max(real(ifftn(H .* fftn(Y))), eps)))), 0)
//   Y          : single gpuArray
//   I          : background-subtracted data, uint16 or single gpuArray
//
//   t          : seconds in [planning, r2c FFTs, c2r FFTs, element-wise kernels] (with t, every
//                step is synchronized to time it)
//
// Out of GPU memory raises parallel:gpu:array:OOM (arrays), parallel:gpu:array:OOMForOperation or
// parallel:gpu:array:FFTPlanInternalError (cuFFT plans), as gpuArray operations do.
//
// Build (compileAllMexFiles builds it on Linux; nvcc of R2024b needs GCC 12 or older, e.g.,
// NVCC_FLAGS='-ccbin /usr/bin/g++'):
// mexcuda -R2018a decon_lucy_rfft_gpu_mex.cu -L<matlabroot>/bin/glnxa64 -l:libcufft.so.11

#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>
#include <cufft.h>

typedef long long idx_t;

static const int kThreads = 256;

static int grid_size(const idx_t n) {
    const idx_t b = (n + kThreads - 1) / kThreads;
    return (int)(b < 65535 * 32 ? (b > 0 ? b : 1) : 65535 * 32);
}

// background-subtracted input as in the GPU code: single(I - uint16(background)) or max(I - background, 0)
struct InputU16 {
    const uint16_t* raw;
    uint16_t bg;
    __device__ __forceinline__ float operator()(const idx_t i) const {
        const uint16_t v = raw[i];
        return v > bg ? (float)(uint16_t)(v - bg) : 0.0f;
    }
};

struct InputSingle {
    const float* raw;
    double bg;
    __device__ __forceinline__ float operator()(const idx_t i) const {
        const float v = (float)((double)raw[i] - bg);
        return v > 0.0f ? v : 0.0f;
    }
};

template <typename In>
__global__ void init_kernel(float* j, float* y, const In in, const idx_t n) {
    for (idx_t i = blockIdx.x * (idx_t)blockDim.x + threadIdx.x; i < n; i += (idx_t)blockDim.x * gridDim.x) {
        const float v = in(i);
        j[i] = v;
        if (y != NULL) {
            y[i] = v;
        }
    }
}

__global__ void multiply_kernel(cufftComplex* s, const cufftComplex* h, const idx_t nh, const bool conjugate) {
    for (idx_t i = blockIdx.x * (idx_t)blockDim.x + threadIdx.x; i < nh; i += (idx_t)blockDim.x * gridDim.x) {
        const float a = s[i].x, b = s[i].y, c = h[i].x, e = h[i].y;
        if (conjugate) {
            s[i].x = __fadd_rn(__fmul_rn(a, c), __fmul_rn(b, e));
            s[i].y = __fsub_rn(__fmul_rn(b, c), __fmul_rn(a, e));
        } else {
            s[i].x = __fsub_rn(__fmul_rn(a, c), __fmul_rn(b, e));
            s[i].y = __fadd_rn(__fmul_rn(a, e), __fmul_rn(b, c));
        }
    }
}

// R = I ./ max(real(ifftn(.)), eps), with R holding the unscaled inverse FFT
template <typename In>
__global__ void ratio_kernel(float* r, const In in, const float invN, const idx_t n) {
    const float epsf = (float)DBL_EPSILON;
    for (idx_t i = blockIdx.x * (idx_t)blockDim.x + threadIdx.x; i < n; i += (idx_t)blockDim.x * gridDim.x) {
        float cx = __fmul_rn(r[i], invN);
        cx = cx > epsf ? cx : epsf;
        r[i] = __fdiv_rn(in(i), cx);
    }
}

// out = max(real(ifftn(.)) .* J, 0), with R holding the unscaled inverse FFT (out may be r or j)
__global__ void update_kernel(float* out, const float* r, const float* j, const float invN, const idx_t n) {
    for (idx_t i = blockIdx.x * (idx_t)blockDim.x + threadIdx.x; i < n; i += (idx_t)blockDim.x * gridDim.x) {
        const float v = __fmul_rn(__fmul_rn(r[i], invN), j[i]);
        out[i] = v > 0.0f ? v : 0.0f;
    }
}

// zero-padded PSF circularly shifted by -floor(psfSize / 2)
__global__ void scatter_psf_kernel(float* r, const float* psf, const int ps0, const int ps1, const int ps2,
        const int d0, const int d1, const int d2) {
    const idx_t np = (idx_t)ps0 * ps1 * ps2;
    for (idx_t i = blockIdx.x * (idx_t)blockDim.x + threadIdx.x; i < np; i += (idx_t)blockDim.x * gridDim.x) {
        const int y = (int)(i % ps0);
        const int x = (int)((i / ps0) % ps1);
        const int z = (int)(i / ((idx_t)ps0 * ps1));
        const idx_t o = (idx_t)((y - ps0 / 2 + d0) % d0) + (idx_t)((x - ps1 / 2 + d1) % d1) * d0
            + (idx_t)((z - ps2 / 2 + d2) % d2) * d0 * d1;
        r[o] = psf[i];
    }
}

// GPU arrays and cuFFT plans of one call, freed at the end and on errors, also
// when MATLAB raises an error in an allocation (the destructor)
struct Resources {
    mxGPUArray* a[8];
    const mxGPUArray* in[4];
    cufftHandle p[2];
    bool hasPlan[2];
    cudaEvent_t e0, e1;
    bool timing;
    double t[4];
    Resources(const bool timing_) : timing(timing_) {
        for (int i = 0; i < 8; i++) a[i] = NULL;
        for (int i = 0; i < 4; i++) in[i] = NULL;
        for (int i = 0; i < 2; i++) hasPlan[i] = false;
        for (int i = 0; i < 4; i++) t[i] = 0;
        if (timing) {
            cudaEventCreate(&e0);
            cudaEventCreate(&e1);
        }
    }
    ~Resources() {
        release();
    }
    void release() {
        for (int i = 0; i < 2; i++) {
            if (hasPlan[i]) cufftDestroy(p[i]);
            hasPlan[i] = false;
        }
        for (int i = 0; i < 8; i++) {
            if (a[i] != NULL) mxGPUDestroyGPUArray(a[i]);
            a[i] = NULL;
        }
        for (int i = 0; i < 4; i++) {
            if (in[i] != NULL) mxGPUDestroyGPUArray(in[i]);
            in[i] = NULL;
        }
        if (timing) {
            cudaEventDestroy(e0);
            cudaEventDestroy(e1);
            timing = false;
        }
    }
    void start() {
        if (timing) cudaEventRecord(e0);
    }
    void stop(const int k) {
        if (timing) {
            float ms = 0;
            cudaEventRecord(e1);
            cudaEventSynchronize(e1);
            cudaEventElapsedTime(&ms, e0, e1);
            t[k] += ms / 1000.0;
        }
    }
};

static void fail(Resources& res, const char* id, const char* msg) {
    res.release();
    mexErrMsgIdAndTxt(id, msg);
}

static void check_cuda(Resources& res, const char* what) {
    const cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        res.release();
        mexErrMsgIdAndTxt("decon_lucy_rfft_gpu_mex:cudaError", "%s: %s", what, cudaGetErrorString(err));
    }
}

static void check_cufft(Resources& res, const cufftResult r, const char* what) {
    if (r == CUFFT_SUCCESS) {
        return;
    }
    res.release();
    if (r == CUFFT_ALLOC_FAILED) {
        mexErrMsgIdAndTxt("parallel:gpu:array:OOMForOperation", "Out of memory on device for the cuFFT plans (%s).", what);
    }
    if (r == CUFFT_INTERNAL_ERROR && strcmp(what, "cufftMakePlanMany") == 0) {
        // as MATLAB's FFTs report a cuFFT plan that fails, e.g., when the GPU memory is too small
        mexErrMsgIdAndTxt("parallel:gpu:array:FFTPlanInternalError", "cuFFT failed to make a plan (%s).", what);
    }
    mexErrMsgIdAndTxt("decon_lucy_rfft_gpu_mex:cufftError", "cuFFT error %d in %s.", (int)r, what);
}

static float* new_real(Resources& res, const int slot, const mwSize nd, const mwSize* dims) {
    res.a[slot] = mxGPUCreateGPUArray(nd, dims, mxSINGLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
    return (float*)mxGPUGetData(res.a[slot]);
}

// FFT dimensions in row-major order: MATLAB's [ny, nx, nz] are cuFFT's (nz, nx, ny), so the halved
// (last) cuFFT dimension is MATLAB's first; size-1 dimensions other than the last are left out
static int fft_dims(const int* d, int* n) {
    int rank = 0;
    if (d[2] > 1) n[rank++] = d[2];
    if (d[1] > 1) n[rank++] = d[1];
    n[rank++] = d[0];
    return rank;
}

// plans for r2c (slot 0) and c2r (slot 1) with one shared work area
static void make_plans(Resources& res, const int* d, const bool both) {
    int n[3];
    const int rank = fft_dims(d, n);
    size_t ws[2] = {0, 0};
    const cufftType types[2] = {CUFFT_R2C, CUFFT_C2R};
    for (int i = 0; i < (both ? 2 : 1); i++) {
        check_cufft(res, cufftCreate(&res.p[i]), "cufftCreate");
        res.hasPlan[i] = true;
        check_cufft(res, cufftSetAutoAllocation(res.p[i], 0), "cufftSetAutoAllocation");
        check_cufft(res, cufftMakePlanMany(res.p[i], rank, n, NULL, 1, 0, NULL, 1, 0, types[i], 1, &ws[i]), "cufftMakePlanMany");
    }
    const size_t w = ws[0] > ws[1] ? ws[0] : ws[1];
    if (w > 0) {
        // in doubles: gpuArrays have at most 2^31 - 1 elements, and the work area can be larger than 2 GB
        const mwSize wd[1] = {(mwSize)((w + 7) / 8)};
        res.a[7] = mxGPUCreateGPUArray(1, wd, mxDOUBLE_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
        void* wp = mxGPUGetData(res.a[7]);
        for (int i = 0; i < (both ? 2 : 1); i++) {
            check_cufft(res, cufftSetWorkArea(res.p[i], wp), "cufftSetWorkArea");
        }
    }
}

template <typename In>
static mxGPUArray* omw_iterations(Resources& res, const In in, const cufftComplex* hf, const cufftComplex* hb,
        const int nIter, const int* d, const mwSize nd, const mwSize* dims, const mwSize* dimsH) {
    const idx_t n = (idx_t)d[0] * d[1] * d[2];
    const idx_t nh = (idx_t)(d[0] / 2 + 1) * d[1] * d[2];
    const float invN = (float)(1.0 / (double)n);
    const int g = grid_size(n), gh = grid_size(nh);

    float* j = new_real(res, 0, nd, dims);
    init_kernel<<<g, kThreads>>>(j, (float*)NULL, in, n);
    if (nIter == 0) {
        check_cuda(res, "init");
        mxGPUArray* out = res.a[0];
        res.a[0] = NULL;
        return out;
    }
    float* r = new_real(res, 1, nd, dims);
    res.a[2] = mxGPUCreateGPUArray(3, dimsH, mxSINGLE_CLASS, mxCOMPLEX, MX_GPU_DO_NOT_INITIALIZE);
    cufftComplex* s = (cufftComplex*)mxGPUGetData(res.a[2]);
    res.start();
    make_plans(res, d, true);
    res.stop(0);

    for (int k = 1; k <= nIter; k++) {
        res.start();
        check_cufft(res, cufftExecR2C(res.p[0], j, s), "cufftExecR2C");
        res.stop(1);
        res.start();
        multiply_kernel<<<gh, kThreads>>>(s, hf, nh, false);
        res.stop(3);
        res.start();
        check_cufft(res, cufftExecC2R(res.p[1], s, r), "cufftExecC2R");
        res.stop(2);
        res.start();
        ratio_kernel<<<g, kThreads>>>(r, in, invN, n);
        res.stop(3);
        res.start();
        check_cufft(res, cufftExecR2C(res.p[0], r, s), "cufftExecR2C");
        res.stop(1);
        res.start();
        multiply_kernel<<<gh, kThreads>>>(s, hb, nh, false);
        res.stop(3);
        res.start();
        check_cufft(res, cufftExecC2R(res.p[1], s, r), "cufftExecC2R");
        res.stop(2);
        res.start();
        update_kernel<<<g, kThreads>>>(j, r, j, invN, n);
        res.stop(3);
    }
    check_cuda(res, "omw iterations");
    mxGPUArray* out = res.a[0];
    res.a[0] = NULL;
    return out;
}

// one simplified RL update: max(Y .* real(ifftn(conj(H) .* fftn(I ./ max(real(ifftn(H .* fftn(Y))), eps)))), 0)
template <typename In>
static mxGPUArray* rl_step(Resources& res, const float* y, const In in, const cufftComplex* h, const int* d, const mwSize nd,
        const mwSize* dims, const mwSize* dimsH) {
    const idx_t n = (idx_t)d[0] * d[1] * d[2];
    const idx_t nh = (idx_t)(d[0] / 2 + 1) * d[1] * d[2];
    const float invN = (float)(1.0 / (double)n);
    const int g = grid_size(n), gh = grid_size(nh);

    float* r = new_real(res, 0, nd, dims);
    res.a[2] = mxGPUCreateGPUArray(3, dimsH, mxSINGLE_CLASS, mxCOMPLEX, MX_GPU_DO_NOT_INITIALIZE);
    cufftComplex* s = (cufftComplex*)mxGPUGetData(res.a[2]);
    res.start();
    make_plans(res, d, true);
    res.stop(0);

    // an out-of-place real-to-complex cuFFT leaves its input (Y) as it is
    res.start();
    check_cufft(res, cufftExecR2C(res.p[0], (float*)y, s), "cufftExecR2C");
    res.stop(1);
    res.start();
    multiply_kernel<<<gh, kThreads>>>(s, h, nh, false);
    res.stop(3);
    res.start();
    check_cufft(res, cufftExecC2R(res.p[1], s, r), "cufftExecC2R");
    res.stop(2);
    res.start();
    ratio_kernel<<<g, kThreads>>>(r, in, invN, n);
    res.stop(3);
    res.start();
    check_cufft(res, cufftExecR2C(res.p[0], r, s), "cufftExecR2C");
    res.stop(1);
    res.start();
    multiply_kernel<<<gh, kThreads>>>(s, h, nh, true);
    res.stop(3);
    res.start();
    check_cufft(res, cufftExecC2R(res.p[1], s, r), "cufftExecC2R");
    res.stop(2);
    res.start();
    update_kernel<<<g, kThreads>>>(r, r, y, invN, n);
    res.stop(3);
    check_cuda(res, "rl step");
    mxGPUArray* out = res.a[0];
    res.a[0] = NULL;
    return out;
}

// size of a gpuArray as 3 int dimensions
static void get_size(Resources& res, const mxGPUArray* a, int* d) {
    const mwSize nd = mxGPUGetNumberOfDimensions(a);
    const mwSize* dims = mxGPUGetDimensions(a);
    bool valid = nd <= 3;
    for (int i = 0; i < 3; i++) {
        const mwSize v = i < (int)nd ? dims[i] : 1;
        valid = valid && v <= 2147483647;
        d[i] = (int)v;
    }
    mxFree((void*)dims);
    if (!valid) {
        fail(res, "decon_lucy_rfft_gpu_mex:inputError", "Only 2D and 3D data with dimensions below 2^31 are supported.");
    }
}

static const mxGPUArray* check_otf(Resources& res, const int slot, const mxArray* a, const int* d) {
    if (!mxIsGPUArray(a)) {
        fail(res, "decon_lucy_rfft_gpu_mex:inputError", "The OTFs must be complex single half OTF gpuArrays.");
    }
    res.in[slot] = mxGPUCreateFromMxArray(a);
    const mxGPUArray* h = res.in[slot];
    if (mxGPUGetClassID(h) != mxSINGLE_CLASS || mxGPUGetComplexity(h) != mxCOMPLEX) {
        fail(res, "decon_lucy_rfft_gpu_mex:inputError", "The OTFs must be complex single half OTF gpuArrays.");
    }
    int dh[3];
    get_size(res, h, dh);
    if (dh[0] != d[0] / 2 + 1 || dh[1] != d[1] || dh[2] != d[2]) {
        fail(res, "decon_lucy_rfft_gpu_mex:inputError", "The OTF size must be [floor(ny / 2) + 1, nx, nz] of the data.");
    }
    return h;
}

// count from a numeric scalar as a loop 1 : v runs it (0 for v < 1 or NaN)
static int scalar_int(const mxArray* a, const char* name) {
    if (!mxIsNumeric(a) || mxGetNumberOfElements(a) != 1 || mxIsComplex(a) || mxIsGPUArray(a)) {
        mexErrMsgIdAndTxt("decon_lucy_rfft_gpu_mex:inputError", "%s must be a real numeric scalar.", name);
    }
    const double v = mxGetScalar(a);
    if (!(v >= 1)) {
        return 0;
    }
    if (v > 2147483647.0) {
        mexErrMsgIdAndTxt("decon_lucy_rfft_gpu_mex:inputError", "%s is too large.", name);
    }
    return (int)v;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    mxInitGPU();
    if (nrhs < 1 || !mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("decon_lucy_rfft_gpu_mex:inputError", "The first input must be 'otf', 'omw' or 'rlstep'.");
    }
    char method[16];
    mxGetString(prhs[0], method, sizeof(method));
    Resources res(nlhs > 1);

    if (strcmp(method, "otf") == 0) {
        if (nrhs != 3) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_gpu_mex:inputError", "Usage: H = decon_lucy_rfft_gpu_mex('otf', psf, outSize)");
        }
        const mxArray* os = prhs[2];
        if (!mxIsDouble(os) || mxIsComplex(os) || mxIsGPUArray(os) || mxGetNumberOfElements(os) < 1 || mxGetNumberOfElements(os) > 3) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_gpu_mex:inputError", "outSize must be a double vector of 1 to 3 elements.");
        }
        res.in[0] = mxGPUCreateFromMxArray(prhs[1]);
        const mxGPUArray* psf = res.in[0];
        if (mxGPUGetClassID(psf) != mxSINGLE_CLASS || mxGPUGetComplexity(psf) != mxREAL || mxGPUGetNumberOfElements(psf) == 0) {
            fail(res, "decon_lucy_rfft_gpu_mex:inputError", "The PSF must be a real, non-empty single array.");
        }
        int ps[3], d[3];
        get_size(res, psf, ps);
        const double* osv = mxGetDoubles(os);
        for (int i = 0; i < 3; i++) {
            const double v = i < (int)mxGetNumberOfElements(os) ? osv[i] : 1;
            if (!(v >= 1) || v > 2147483647.0 || v != (double)(int)v) {
                fail(res, "decon_lucy_rfft_gpu_mex:inputError", "outSize must contain positive integers.");
            }
            d[i] = (int)v;
            if (ps[i] > d[i]) {
                fail(res, "decon_lucy_rfft_gpu_mex:outSizeIsSmallerThanPsfSize", "outSizeIsSmallerThanPsfSize");
            }
        }
        const mwSize dims[3] = {(mwSize)d[0], (mwSize)d[1], (mwSize)d[2]};
        mwSize dh[3] = {(mwSize)(d[0] / 2 + 1), (mwSize)d[1], (mwSize)d[2]};
        const mwSize ndh = (mxGetNumberOfElements(os) >= 3 && d[2] > 1) ? 3 : 2;
        res.a[0] = mxGPUCreateGPUArray(3, dims, mxSINGLE_CLASS, mxREAL, MX_GPU_INITIALIZE_VALUES);
        float* r = (float*)mxGPUGetData(res.a[0]);
        res.a[1] = mxGPUCreateGPUArray(ndh, dh, mxSINGLE_CLASS, mxCOMPLEX, MX_GPU_DO_NOT_INITIALIZE);
        cufftComplex* s = (cufftComplex*)mxGPUGetData(res.a[1]);
        res.start();
        make_plans(res, d, false);
        res.stop(0);
        res.start();
        const idx_t np = (idx_t)ps[0] * ps[1] * ps[2];
        scatter_psf_kernel<<<grid_size(np), kThreads>>>(r, (const float*)mxGPUGetDataReadOnly(psf), ps[0], ps[1], ps[2], d[0], d[1], d[2]);
        res.stop(3);
        res.start();
        check_cufft(res, cufftExecR2C(res.p[0], r, s), "cufftExecR2C");
        res.stop(1);
        check_cuda(res, "otf");
        plhs[0] = mxGPUCreateMxArrayOnGPU(res.a[1]);
    } else if (strcmp(method, "omw") == 0) {
        if (nrhs != 6) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_gpu_mex:inputError", "Usage: J = decon_lucy_rfft_gpu_mex('omw', I, H_f, H_b, nIter, background)");
        }
        const int nIter = scalar_int(prhs[4], "nIter");
        const mxArray* bga = prhs[5];
        if (!mxIsNumeric(bga) || mxGetNumberOfElements(bga) != 1 || mxIsComplex(bga) || mxIsGPUArray(bga)) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_gpu_mex:inputError", "background must be a real numeric scalar.");
        }
        const double bg = mxGetScalar(bga);
        res.in[0] = mxGPUCreateFromMxArray(prhs[1]);
        const mxGPUArray* im = res.in[0];
        const mxClassID cls = mxGPUGetClassID(im);
        if (!(cls == mxUINT16_CLASS || cls == mxSINGLE_CLASS) || mxGPUGetComplexity(im) != mxREAL || mxGPUGetNumberOfElements(im) == 0) {
            fail(res, "decon_lucy_rfft_gpu_mex:inputError", "The data must be a real, non-empty uint16 or single array.");
        }
        int d[3];
        get_size(res, im, d);
        const mwSize nd = mxGPUGetNumberOfDimensions(im);
        const mwSize* dims = mxGPUGetDimensions(im);
        const mwSize dimsH[3] = {(mwSize)(d[0] / 2 + 1), (mwSize)d[1], (mwSize)d[2]};
        const cufftComplex* h1 = (const cufftComplex*)mxGPUGetDataReadOnly(check_otf(res, 1, prhs[2], d));
        const cufftComplex* h2 = (const cufftComplex*)mxGPUGetDataReadOnly(check_otf(res, 2, prhs[3], d));

        mxGPUArray* out;
        if (cls == mxUINT16_CLASS) {
            InputU16 in;
            in.raw = (const uint16_t*)mxGPUGetDataReadOnly(im);
            // uint16(background): rounded half away from zero, saturated, NaN is 0
            in.bg = (uint16_t)(bg != bg || bg <= 0 ? 0 : (bg >= 65535 ? 65535 : floor(bg + 0.5)));
            out = omw_iterations(res, in, h1, h2, nIter, d, nd, dims, dimsH);
        } else {
            InputSingle in;
            in.raw = (const float*)mxGPUGetDataReadOnly(im);
            in.bg = bg;
            out = omw_iterations(res, in, h1, h2, nIter, d, nd, dims, dimsH);
        }
        mxFree((void*)dims);
        plhs[0] = mxGPUCreateMxArrayOnGPU(out);
        mxGPUDestroyGPUArray(out);
    } else if (strcmp(method, "rlstep") == 0) {
        if (nrhs != 4) {
            mexErrMsgIdAndTxt("decon_lucy_rfft_gpu_mex:inputError", "Usage: J = decon_lucy_rfft_gpu_mex('rlstep', Y, I, H)");
        }
        res.in[0] = mxGPUCreateFromMxArray(prhs[1]);
        res.in[3] = mxGPUCreateFromMxArray(prhs[2]);
        const mxGPUArray* ya = res.in[0];
        const mxGPUArray* im = res.in[3];
        const mxClassID cls = mxGPUGetClassID(im);
        int d[3], di[3];
        if (mxGPUGetClassID(ya) != mxSINGLE_CLASS || mxGPUGetComplexity(ya) != mxREAL || mxGPUGetNumberOfElements(ya) == 0
                || !(cls == mxUINT16_CLASS || cls == mxSINGLE_CLASS) || mxGPUGetComplexity(im) != mxREAL) {
            fail(res, "decon_lucy_rfft_gpu_mex:inputError", "Y must be a real single and I a real uint16 or single array.");
        }
        get_size(res, ya, d);
        get_size(res, im, di);
        if (d[0] != di[0] || d[1] != di[1] || d[2] != di[2]) {
            fail(res, "decon_lucy_rfft_gpu_mex:inputError", "Y and I must have the same size.");
        }
        const mwSize nd = mxGPUGetNumberOfDimensions(ya);
        const mwSize* dims = mxGPUGetDimensions(ya);
        const mwSize dimsH[3] = {(mwSize)(d[0] / 2 + 1), (mwSize)d[1], (mwSize)d[2]};
        const cufftComplex* h = (const cufftComplex*)mxGPUGetDataReadOnly(check_otf(res, 1, prhs[3], d));
        const float* y = (const float*)mxGPUGetDataReadOnly(ya);

        // I is background subtracted already (no background here)
        mxGPUArray* out;
        if (cls == mxUINT16_CLASS) {
            InputU16 in;
            in.raw = (const uint16_t*)mxGPUGetDataReadOnly(im);
            in.bg = 0;
            out = rl_step(res, y, in, h, d, nd, dims, dimsH);
        } else {
            InputSingle in;
            in.raw = (const float*)mxGPUGetDataReadOnly(im);
            in.bg = 0;
            out = rl_step(res, y, in, h, d, nd, dims, dimsH);
        }
        mxFree((void*)dims);
        plhs[0] = mxGPUCreateMxArrayOnGPU(out);
        mxGPUDestroyGPUArray(out);
    } else {
        mexErrMsgIdAndTxt("decon_lucy_rfft_gpu_mex:inputError", "Unknown method %s, use 'otf', 'omw' or 'rlstep'.", method);
    }

    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(1, 4, mxREAL);
        double* t = mxGetDoubles(plhs[1]);
        for (int i = 0; i < 4; i++) {
            t[i] = res.t[i];
        }
    }
    res.release();
}
