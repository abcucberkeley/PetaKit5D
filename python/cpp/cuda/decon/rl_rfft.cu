// The half-spectrum Richardson-Lucy steps of decon_lucy_rfft_gpu_mex
// (microscopeDataProcessing/decon/mex/decon_lucy_rfft_gpu_mex.cu) for petakit5d.kernels.decon_gpu:
// the mex's element-wise kernels with its arithmetic, and the same steps for the port's (z, y, x)
// volumes next to arrays in MATLAB's layout.
//
// Layouts: the mex's arrays are MATLAB's (y, x, z) Fortran volumes, which are the C arrays
// (z, x, y) ("internal" arrays); the port's volumes are (z, y, x) C arrays. The linear kernels (as
// the mex's) run on arrays of one layout. The *_t kernels read a (z, y, x) volume next to an
// internal array, or write one from it, through 32 x 32 tiles of each z plane in shared memory, so
// that both are read and written coalesced. Every value is computed from the same operands by the
// same operations as in the mex, so neither the layout nor the grid changes a bit. No atomics, no
// reductions: every kernel is element-wise.
//
// Compiled by NVRTC with --fmad=false (petakit5d.gpu.kernels); the rounding intrinsics of the mex
// are kept.

typedef long long idx_t;

#define TILE 32
#define ROWS 8

// grid-stride loop over 0 .. n - 1, as the mex's kernels
#define FOR_EACH(i, n) \
    for (idx_t i = blockIdx.x * (idx_t)blockDim.x + threadIdx.x; i < (n); i += (idx_t)blockDim.x * gridDim.x)

// the background-subtracted input as the mex's InputU16 and InputSingle: single(I - uint16(bg))
// (saturating; bg is the rounded uint16 value) and max(single(double(I) - bg), 0)
template <typename T> struct Input;

template <> struct Input<unsigned short> {
    const unsigned short* raw;
    unsigned short bg;
    __device__ Input(const unsigned short* r, const double b) : raw(r), bg((unsigned short)b) {}
    __device__ __forceinline__ float operator()(const idx_t i) const {
        const unsigned short v = raw[i];
        return v > bg ? (float)(unsigned short)(v - bg) : 0.0f;
    }
};

template <> struct Input<float> {
    const float* raw;
    double bg;
    __device__ Input(const float* r, const double b) : raw(r), bg(b) {}
    __device__ __forceinline__ float operator()(const idx_t i) const {
        const float v = (float)((double)raw[i] - bg);
        return v > 0.0f ? v : 0.0f;
    }
};

// MATLAB's eps (double), as single: 2^-52
#define EPS_F 0x1p-52f

__device__ __forceinline__ float ratio(const float in, const float r, const float inv_n) {
    float cx = __fmul_rn(r, inv_n);
    cx = cx > EPS_F ? cx : EPS_F;
    return __fdiv_rn(in, cx);
}

__device__ __forceinline__ float update(const float r, const float j, const float inv_n) {
    const float v = __fmul_rn(__fmul_rn(r, inv_n), j);
    return v > 0.0f ? v : 0.0f;
}

// ---- the mex's kernels (one layout) ------------------------------------------------------------

// j = the background-subtracted input (init_kernel)
template <typename T>
__global__ void init_kernel(float* __restrict__ j, const T* __restrict__ raw, const double bg, const idx_t n) {
    const Input<T> in(raw, bg);
    FOR_EACH(i, n) {
        j[i] = in(i);
    }
}

// s = s .* h or s .* conj(h), MATLAB's complex product without fused multiply-adds (multiply_kernel)
extern "C" __global__ void multiply_kernel(float2* s, const float2* __restrict__ h, const idx_t nh,
                                           const int conjugate) {
    FOR_EACH(i, nh) {
        const float2 u = s[i], w = h[i];
        const float a = u.x, b = u.y, c = w.x, e = w.y;
        if (conjugate) {
            s[i] = make_float2(__fadd_rn(__fmul_rn(a, c), __fmul_rn(b, e)), __fsub_rn(__fmul_rn(b, c), __fmul_rn(a, e)));
        } else {
            s[i] = make_float2(__fsub_rn(__fmul_rn(a, c), __fmul_rn(b, e)), __fadd_rn(__fmul_rn(a, e), __fmul_rn(b, c)));
        }
    }
}

// r = I ./ max(real(ifftn(.)), eps), r holding the unscaled inverse FFT (ratio_kernel)
template <typename T>
__global__ void ratio_kernel(float* r, const T* __restrict__ raw, const double bg, const float inv_n, const idx_t n) {
    const Input<T> in(raw, bg);
    FOR_EACH(i, n) {
        r[i] = ratio(in(i), r[i], inv_n);
    }
}

// out = max(real(ifftn(.)) .* j, 0), r holding the unscaled inverse FFT; out may be r or j (update_kernel)
extern "C" __global__ void update_kernel(float* out, const float* r, const float* j, const float inv_n, const idx_t n) {
    FOR_EACH(i, n) {
        out[i] = update(r[i], j[i], inv_n);
    }
}

// the zero-padded PSF circularly shifted by -floor(psf size / 2) (scatter_psf_kernel); the PSF is a
// (z, y, x) C array, the output has the strides (sz, sy, sx): MATLAB's layout or the port's
extern "C" __global__ void scatter_psf_kernel(float* __restrict__ r, const float* __restrict__ psf, const int pz,
                                              const int py, const int px, const int nz, const int ny, const int nx,
                                              const idx_t sz, const idx_t sy, const idx_t sx) {
    const idx_t np = (idx_t)pz * py * px;
    FOR_EACH(i, np) {
        const int x = (int)(i % px);
        const int y = (int)((i / px) % py);
        const int z = (int)(i / ((idx_t)px * py));
        r[(idx_t)((y - py / 2 + ny) % ny) * sy + (idx_t)((x - px / 2 + nx) % nx) * sx
          + (idx_t)((z - pz / 2 + nz) % nz) * sz] = psf[i];
    }
}

// ---- (z, y, x) volumes next to internal arrays (z, x, y), by tiles of each z plane ---------------
// Block (32, 8), grid (ceil(nx / 32), ceil(ny / 32), planes); z planes beyond the grid in a loop.
// Phase "zyx": thread (tx, ty) visits x = x0 + tx, y = y0 + ty + k (coalesced in a (z, y, x)
// plane); phase "int": y = y0 + tx, x = x0 + ty + k (coalesced in an internal plane).

#define PLANES for (idx_t z = blockIdx.z; z < nz; z += gridDim.z)

// dst (internal) = src ((z, y, x)): a transpose of each plane
template <typename T>
__global__ void zyx_to_int(T* __restrict__ dst, const T* __restrict__ src, const int ny, const int nx, const int nz) {
    __shared__ T tile[TILE][TILE + 1];
    const int x0 = blockIdx.x * TILE, y0 = blockIdx.y * TILE;
    PLANES {
        const idx_t plane = z * ny * nx;
        for (int k = 0; k < TILE; k += ROWS) {
            const int x = x0 + threadIdx.x, y = y0 + threadIdx.y + k;
            if (x < nx && y < ny) tile[threadIdx.y + k][threadIdx.x] = src[plane + (idx_t)y * nx + x];
        }
        __syncthreads();
        for (int k = 0; k < TILE; k += ROWS) {
            const int y = y0 + threadIdx.x, x = x0 + threadIdx.y + k;
            if (x < nx && y < ny) dst[plane + (idx_t)x * ny + y] = tile[threadIdx.x][threadIdx.y + k];
        }
        __syncthreads();
    }
}

// dst ((z, y, x)) = src (internal)
template <typename T>
__global__ void int_to_zyx(T* __restrict__ dst, const T* __restrict__ src, const int ny, const int nx, const int nz) {
    __shared__ T tile[TILE][TILE + 1];
    const int x0 = blockIdx.x * TILE, y0 = blockIdx.y * TILE;
    PLANES {
        const idx_t plane = z * ny * nx;
        for (int k = 0; k < TILE; k += ROWS) {
            const int y = y0 + threadIdx.x, x = x0 + threadIdx.y + k;
            if (x < nx && y < ny) tile[threadIdx.y + k][threadIdx.x] = src[plane + (idx_t)x * ny + y];
        }
        __syncthreads();
        for (int k = 0; k < TILE; k += ROWS) {
            const int x = x0 + threadIdx.x, y = y0 + threadIdx.y + k;
            if (x < nx && y < ny) dst[plane + (idx_t)y * nx + x] = tile[threadIdx.x][threadIdx.y + k];
        }
        __syncthreads();
    }
}

// j (internal) = the background-subtracted input of a (z, y, x) volume
template <typename T>
__global__ void init_t(float* __restrict__ j, const T* __restrict__ raw, const double bg, const int ny, const int nx,
                       const int nz) {
    __shared__ float tile[TILE][TILE + 1];
    const Input<T> in(raw, bg);
    const int x0 = blockIdx.x * TILE, y0 = blockIdx.y * TILE;
    PLANES {
        const idx_t plane = z * ny * nx;
        for (int k = 0; k < TILE; k += ROWS) {
            const int x = x0 + threadIdx.x, y = y0 + threadIdx.y + k;
            if (x < nx && y < ny) tile[threadIdx.y + k][threadIdx.x] = in(plane + (idx_t)y * nx + x);
        }
        __syncthreads();
        for (int k = 0; k < TILE; k += ROWS) {
            const int y = y0 + threadIdx.x, x = x0 + threadIdx.y + k;
            if (x < nx && y < ny) j[plane + (idx_t)x * ny + y] = tile[threadIdx.x][threadIdx.y + k];
        }
        __syncthreads();
    }
}

// r (internal) = I ./ max(r * inv_n, eps), I the background-subtracted input of a (z, y, x) volume
template <typename T>
__global__ void ratio_t(float* r, const T* __restrict__ raw, const double bg, const float inv_n, const int ny,
                        const int nx, const int nz) {
    __shared__ float tile[TILE][TILE + 1];
    const Input<T> in(raw, bg);
    const int x0 = blockIdx.x * TILE, y0 = blockIdx.y * TILE;
    PLANES {
        const idx_t plane = z * ny * nx;
        for (int k = 0; k < TILE; k += ROWS) {
            const int x = x0 + threadIdx.x, y = y0 + threadIdx.y + k;
            if (x < nx && y < ny) tile[threadIdx.y + k][threadIdx.x] = in(plane + (idx_t)y * nx + x);
        }
        __syncthreads();
        for (int k = 0; k < TILE; k += ROWS) {
            const int y = y0 + threadIdx.x, x = x0 + threadIdx.y + k;
            if (x < nx && y < ny) {
                const idx_t o = plane + (idx_t)x * ny + y;
                r[o] = ratio(tile[threadIdx.x][threadIdx.y + k], r[o], inv_n);
            }
        }
        __syncthreads();
    }
}

// out ((z, y, x)) = max(r * inv_n .* j, 0) with r and j internal (OMW's last update)
extern "C" __global__ void update_out_t(float* __restrict__ out, const float* __restrict__ r,
                                        const float* __restrict__ j, const float inv_n, const int ny, const int nx,
                                        const int nz) {
    __shared__ float tile[TILE][TILE + 1];
    const int x0 = blockIdx.x * TILE, y0 = blockIdx.y * TILE;
    PLANES {
        const idx_t plane = z * ny * nx;
        for (int k = 0; k < TILE; k += ROWS) {
            const int y = y0 + threadIdx.x, x = x0 + threadIdx.y + k;
            if (x < nx && y < ny) {
                const idx_t o = plane + (idx_t)x * ny + y;
                tile[threadIdx.y + k][threadIdx.x] = update(r[o], j[o], inv_n);
            }
        }
        __syncthreads();
        for (int k = 0; k < TILE; k += ROWS) {
            const int x = x0 + threadIdx.x, y = y0 + threadIdx.y + k;
            if (x < nx && y < ny) out[plane + (idx_t)y * nx + x] = tile[threadIdx.x][threadIdx.y + k];
        }
        __syncthreads();
    }
}

// out ((z, y, x)) = max(r * inv_n .* y, 0) with r internal and y a (z, y, x) volume (RL's update)
extern "C" __global__ void update_y_t(float* __restrict__ out, const float* __restrict__ r,
                                      const float* __restrict__ yv, const float inv_n, const int ny, const int nx,
                                      const int nz) {
    __shared__ float tile[TILE][TILE + 1];
    const int x0 = blockIdx.x * TILE, y0 = blockIdx.y * TILE;
    PLANES {
        const idx_t plane = z * ny * nx;
        for (int k = 0; k < TILE; k += ROWS) {
            const int y = y0 + threadIdx.x, x = x0 + threadIdx.y + k;
            if (x < nx && y < ny) tile[threadIdx.y + k][threadIdx.x] = r[plane + (idx_t)x * ny + y];
        }
        __syncthreads();
        for (int k = 0; k < TILE; k += ROWS) {
            const int x = x0 + threadIdx.x, y = y0 + threadIdx.y + k;
            if (x < nx && y < ny) {
                const idx_t o = plane + (idx_t)y * nx + x;
                out[o] = update(tile[threadIdx.x][threadIdx.y + k], yv[o], inv_n);
            }
        }
        __syncthreads();
    }
}
