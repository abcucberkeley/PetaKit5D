// GPU version of decon_mask_edge_erosion_mex for gpuArray masks (the GPU decon): the same erosions with the
// same result as imerode (on the CPU and on gpuArrays): each z slice padded with true by padUL and padLR,
// eroded by each 2D neighborhood in turn (pixels outside the padded slice count as true) and cropped back, then
// the erosion along z by a line of 2 * rz + 1 voxels (voxels outside the volume count as true).
//
// B = decon_mask_edge_erosion_gpu_mex(A, nhoods, padUL, padLR, rz, clearBorder)
//   A           : logical gpuArray (or host array), 2D or 3D
//   nhoods      : cell array of 2D neighborhoods (logical or numeric, each with a nonzero element), in order
//   padUL       : [rows, columns] of true padding before each slice (0 for none)
//   padLR       : [rows, columns] of true padding after each slice
//   rz          : radius of the erosion along z (0 for none; ignored for 2D)
//   clearBorder : true to set the first and last rows and columns (3D: and slices) to false first
//   B           : logical gpuArray
//
// Build (compileAllMexFiles builds it on Linux; nvcc of R2024b needs GCC 12 or older, e.g.,
// NVCC_FLAGS='-ccbin /usr/bin/g++'):
// mexcuda decon_mask_edge_erosion_gpu_mex.cu

#include <stdint.h>
#include <string.h>
#include <vector>
#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cuda_runtime.h>

typedef long long idx_t;

static const int kMaxOffsets = 4096;

// offsets (row, column) of the nonzero elements of all neighborhoods
__constant__ int2 c_off[kMaxOffsets];

struct Dims {
    idx_t h, w, d;      // mask
    idx_t pu0, pu1;     // padding before (rows, columns)
    idx_t hp, wp;       // padded slice
    int clearBorder;
};

// blocks of 64 rows x 4 columns of a slice, slices along z (grid-stride loops)
static const dim3 kBlock(64, 4, 1);

static dim3 grid_3d(const idx_t rows, const idx_t cols, const idx_t slices) {
    const idx_t gx = (rows + kBlock.x - 1) / kBlock.x, gy = (cols + kBlock.y - 1) / kBlock.y;
    return dim3((unsigned)(gx < 65535 ? gx : 65535), (unsigned)(gy < 65535 ? gy : 65535),
        (unsigned)(slices < 65535 ? slices : 65535));
}

// the padded volume: true outside the slice, the (normalized) mask inside, with the border cleared
__global__ void load_kernel(uint8_t* p, const uint8_t* a, const Dims g) {
    for (idx_t z = blockIdx.z; z < g.d; z += gridDim.z) {
        const bool clearSlice = g.clearBorder && g.d > 1 && (z == 0 || z == g.d - 1);
        for (idx_t c = blockIdx.y * (idx_t)blockDim.y + threadIdx.y; c < g.wp; c += (idx_t)gridDim.y * blockDim.y) {
            const idx_t x = c - g.pu1;
            for (idx_t r = blockIdx.x * (idx_t)blockDim.x + threadIdx.x; r < g.hp; r += (idx_t)gridDim.x * blockDim.x) {
                const idx_t y = r - g.pu0;
                uint8_t v = 1;
                if (y >= 0 && y < g.h && x >= 0 && x < g.w) {
                    const bool clear = clearSlice || (g.clearBorder && (y == 0 || y == g.h - 1 || x == 0 || x == g.w - 1));
                    v = clear ? 0 : (a[y + g.h * (x + g.w * z)] != 0);
                }
                p[r + g.hp * (c + g.wp * z)] = v;
            }
        }
    }
}

// erosion by the neighborhood with offsets c_off[i0 .. i0 + n - 1] of the padded volume, from a to b
__global__ void erode_kernel(uint8_t* b, const uint8_t* a, const int i0, const int n, const Dims g) {
    for (idx_t z = blockIdx.z; z < g.d; z += gridDim.z) {
        const uint8_t* as = a + g.hp * g.wp * z;
        uint8_t* bs = b + g.hp * g.wp * z;
        for (idx_t c = blockIdx.y * (idx_t)blockDim.y + threadIdx.y; c < g.wp; c += (idx_t)gridDim.y * blockDim.y) {
            for (idx_t r = blockIdx.x * (idx_t)blockDim.x + threadIdx.x; r < g.hp; r += (idx_t)gridDim.x * blockDim.x) {
                uint8_t v = 1;
                for (int k = 0; k < n; k++) {
                    const int2 o = c_off[i0 + k];
                    const idx_t rr = r + o.x, cc = c + o.y;
                    if (rr >= 0 && rr < g.hp && cc >= 0 && cc < g.wp) {
                        v &= as[rr + g.hp * cc];
                    }
                }
                bs[r + g.hp * c] = v;
            }
        }
    }
}

// the cropped slices eroded along z by 2 * rz + 1 voxels: one thread per voxel column, counting the false
// voxels in its window
__global__ void z_kernel(uint8_t* out, const uint8_t* p, const idx_t rz, const Dims g) {
    const idx_t n2 = g.h * g.w, np2 = g.hp * g.wp;
    for (idx_t x = blockIdx.y * (idx_t)blockDim.y + threadIdx.y; x < g.w; x += (idx_t)gridDim.y * blockDim.y) {
        for (idx_t y = blockIdx.x * (idx_t)blockDim.x + threadIdx.x; y < g.h; y += (idx_t)gridDim.x * blockDim.x) {
            const uint8_t* col = p + (y + g.pu0) + g.hp * (x + g.pu1);
            uint8_t* o = out + y + g.h * x;
            idx_t cnt = 0;
            const idx_t zi = rz < g.d - 1 ? rz : g.d - 1;
            for (idx_t z = 0; z <= zi; z++) cnt += 1 - col[np2 * z];
            for (idx_t z = 0; z < g.d; z++) {
                o[n2 * z] = cnt == 0;
                if (z + rz + 1 < g.d) cnt += 1 - col[np2 * (z + rz + 1)];
                if (z - rz >= 0) cnt -= 1 - col[np2 * (z - rz)];
            }
        }
    }
}

// GPU arrays of one call, freed at the end and before errors
struct Resources {
    const mxGPUArray* in;
    mxGPUArray* a[3];
    Resources() : in(NULL) {
        for (int i = 0; i < 3; i++) a[i] = NULL;
    }
    ~Resources() {
        release();
    }
    void release() {
        if (in != NULL) mxGPUDestroyGPUArray(in);
        in = NULL;
        for (int i = 0; i < 3; i++) {
            if (a[i] != NULL) mxGPUDestroyGPUArray(a[i]);
            a[i] = NULL;
        }
    }
};

static void fail(Resources& res, const char* msg) {
    res.release();
    mexErrMsgIdAndTxt("decon_mask_edge_erosion_gpu_mex:inputError", msg);
}

static void check_cuda(Resources& res, const char* what) {
    const cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        res.release();
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_gpu_mex:cudaError", "%s: %s", what, cudaGetErrorString(err));
    }
}

static double scalar_arg(const mxArray* a, const char* name) {
    if (!(mxIsNumeric(a) || mxIsLogical(a)) || mxIsComplex(a) || mxIsGPUArray(a) || mxGetNumberOfElements(a) != 1) {
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_gpu_mex:inputError", "%s must be a real numeric scalar.", name);
    }
    return mxGetScalar(a);
}

static void pad_arg(const mxArray* a, const char* name, idx_t* v) {
    if (!mxIsDouble(a) || mxIsComplex(a) || mxIsGPUArray(a) || mxGetNumberOfElements(a) != 2) {
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_gpu_mex:inputError", "%s must be a double vector of 2 elements.", name);
    }
    const double* pv = mxGetPr(a);
    for (int i = 0; i < 2; i++) {
        if (!(pv[i] >= 0) || pv[i] > 1e9 || pv[i] != (double)(idx_t)pv[i]) {
            mexErrMsgIdAndTxt("decon_mask_edge_erosion_gpu_mex:inputError", "%s must contain nonnegative integers.", name);
        }
        v[i] = (idx_t)pv[i];
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    mxInitGPU();
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_gpu_mex:inputError",
            "Usage: B = decon_mask_edge_erosion_gpu_mex(A, nhoods, padUL, padLR, rz, clearBorder)");
    }
    const mxArray* nh = prhs[1];
    if (!mxIsCell(nh)) {
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_gpu_mex:inputError", "nhoods must be a cell array.");
    }
    idx_t pu[2], pl[2];
    pad_arg(prhs[2], "padUL", pu);
    pad_arg(prhs[3], "padLR", pl);
    const double rzv = scalar_arg(prhs[4], "rz");
    if (!(rzv >= 0) || rzv > 1e9 || rzv != (double)(idx_t)rzv) {
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_gpu_mex:inputError", "rz must be a nonnegative integer.");
    }
    const bool clearBorder = scalar_arg(prhs[5], "clearBorder") != 0;

    // offsets of the nonzero elements of the neighborhoods (center floor((size + 1) / 2))
    const int ns = (int)mxGetNumberOfElements(nh);
    std::vector<int2> off;
    std::vector<int> start(1, 0);
    for (int k = 0; k < ns; k++) {
        const mxArray* e = mxGetCell(nh, k);
        if (e == NULL || !(mxIsLogical(e) || mxIsNumeric(e)) || mxIsComplex(e) || mxIsGPUArray(e) || mxIsSparse(e)
                || mxGetNumberOfDimensions(e) > 2 || mxIsEmpty(e)) {
            mexErrMsgIdAndTxt("decon_mask_edge_erosion_gpu_mex:inputError", "Each neighborhood must be a real 2D array.");
        }
        mxArray* tmp = NULL;
        const mxArray* lg = e;
        if (!mxIsLogical(e)) {
            mxArray* in = (mxArray*)e;
            if (mexCallMATLAB(1, &tmp, 1, &in, "logical") != 0) {
                mexErrMsgIdAndTxt("decon_mask_edge_erosion_gpu_mex:inputError", "A neighborhood can't be converted to logical.");
            }
            lg = tmp;
        }
        const idx_t m = (idx_t)mxGetM(lg), n = (idx_t)mxGetN(lg);
        const mxLogical* v = mxGetLogicals(lg);
        const idx_t cy = (m + 1) / 2 - 1, cx = (n + 1) / 2 - 1;
        int cnt = 0;
        for (idx_t j = 0; j < n; j++) {
            for (idx_t i = 0; i < m; i++) {
                if (v[i + j * m]) {
                    off.push_back(make_int2((int)(i - cy), (int)(j - cx)));
                    cnt++;
                }
            }
        }
        if (tmp != NULL) mxDestroyArray(tmp);
        if (cnt == 0) {
            mexErrMsgIdAndTxt("decon_mask_edge_erosion_gpu_mex:inputError", "Each neighborhood must have a nonzero element.");
        }
        start.push_back(start.back() + cnt);
    }
    if ((int)off.size() > kMaxOffsets) {
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_gpu_mex:inputError", "The neighborhoods have more than %d elements.", kMaxOffsets);
    }

    Resources res;
    res.in = mxGPUCreateFromMxArray(prhs[0]);
    const mwSize nd = mxGPUGetNumberOfDimensions(res.in);
    if (mxGPUGetClassID(res.in) != mxLOGICAL_CLASS || mxGPUIsSparse(res.in) || nd > 3) {
        fail(res, "The mask must be a full 2D or 3D logical array.");
    }
    const mwSize* dims = mxGPUGetDimensions(res.in);
    Dims g;
    g.h = (idx_t)dims[0];
    g.w = (idx_t)dims[1];
    g.d = nd > 2 ? (idx_t)dims[2] : 1;
    g.pu0 = pu[0];
    g.pu1 = pu[1];
    g.hp = g.h + pu[0] + pl[0];
    g.wp = g.w + pu[1] + pl[1];
    g.clearBorder = clearBorder ? 1 : 0;
    res.a[0] = mxGPUCreateGPUArray(nd, dims, mxLOGICAL_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
    mxFree((void*)dims);
    const idx_t n = g.h * g.w * g.d;
    if (n > 0) {
        if (!off.empty()) {
            cudaMemcpyToSymbol(c_off, off.data(), sizeof(int2) * off.size());
            check_cuda(res, "offsets");
        }
        const mwSize pd[3] = {(mwSize)g.hp, (mwSize)g.wp, (mwSize)g.d};
        res.a[1] = mxGPUCreateGPUArray(3, pd, mxUINT8_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
        res.a[2] = mxGPUCreateGPUArray(3, pd, mxUINT8_CLASS, mxREAL, MX_GPU_DO_NOT_INITIALIZE);
        uint8_t* p = (uint8_t*)mxGPUGetData(res.a[1]);
        uint8_t* q = (uint8_t*)mxGPUGetData(res.a[2]);
        load_kernel<<<grid_3d(g.hp, g.wp, g.d), kBlock>>>(p, (const uint8_t*)mxGPUGetDataReadOnly(res.in), g);
        for (int k = 0; k < ns; k++) {
            erode_kernel<<<grid_3d(g.hp, g.wp, g.d), kBlock>>>(q, p, start[k], start[k + 1] - start[k], g);
            uint8_t* t = p;
            p = q;
            q = t;
        }
        const idx_t rz = g.d > 1 ? (idx_t)rzv : 0;
        z_kernel<<<grid_3d(g.h, g.w, 1), kBlock>>>((uint8_t*)mxGPUGetData(res.a[0]), p, rz, g);
        check_cuda(res, "erosion");
    }
    plhs[0] = mxGPUCreateMxArrayOnGPU(res.a[0]);
}
