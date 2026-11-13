"""Richardson-Lucy deconvolution on the GPU: the gpuArray branches of decon_lucy_function.m
(``decon_lucy``) and decon_lucy_omw_function.m (``decon_lucy_omw``).

``petakit5d.decon.rl`` hands a call here when ``use_gpu`` is True or None and a usable GPU was
selected; the keywords are the registered functions'. Called directly, ``use_gpu`` False runs on
the CPU (``rl.lucy`` / ``rl.lucy_omw``), None runs on the CPU with a warning when no GPU is usable,
True raises ``NoUsableGPUError`` (identifier ``decon:noUsableGPU``) then.

**Paths, as MATLAB chooses them on the GPU** (decon_use_rfft_mex with a GPU):

- the default (``use_mex`` None or False) is MATLAB's own gpuArray loop on full complex spectra:
  ``fftn``/``ifftn`` with MATLAB's cuFFT plans (``petakit5d.gpu.fft``: C2C plans of MATLAB's
  ``(y, x, z)`` layout, the spectra kept in that layout) and MATLAB's element-wise steps in run-time
  compiled kernels. MATLAB uses the GPU mex by default for no method.
- ``use_mex=True``: the half spectra of decon_lucy_rfft_gpu_mex (``kernels.decon_gpu.LucyRFFTGPU``):
  OMW's iterations in the mex, RL's ``rlstep`` in the mex with the acceleration in MATLAB's GPU
  code; for uint16 or single data, a scalar background and no debug output, else MATLAB's error.
  ``layout='native'`` (Python only) runs the mex's steps on half spectra of the (z, y, x) arrays.

**Bits.** MATLAB's gpuArray arithmetic as R04's probes measured it (A100, B200, Titan V): every
operation rounded once (no fused multiply-adds, also in complex products), single with double
computed in double and rounded to single (``I - Background``, ``J * scaleFactor + deconOffset``,
``dampFactor * I``, double data with single spectra), IEEE division, subnormals kept,
``max(-0, 0)`` and ``max(NaN, 0)`` +0, ``uint16`` of single rounded half away from zero
(0.49999997 gives 0, where MATLAB's CPU ``uint16`` gives 1), NaN results as CUDA's. The FFTs are
MATLAB's with cuFFT 11.0.8 on the same GPU model (F19). RL's acceleration factor and the PSF
normalization are MATLAB's gpuArray sums (``matlab_gpu_sum``). uint16 data stay uint16 on the
GPU (``I - uint16(Background)`` with MATLAB's rounding of the background; the mex path rounds it as
the mex does).

**Memory.** The data are put in MATLAB's layout once on the GPU (tiled transposes, no host
copies), the loops run in place on one spectrum (OMW: the two OTFs, the data, J and a spectrum; RL:
the OTF, the data, J_2, J_3, J_4, Y and a spectrum) and the FFT plan's work area. ``memory_need``
gives a call's bytes; a call that doesn't fit raises before it allocates anything. The result is
post-processed, converted to uint16 and cropped to ``bbox`` on the GPU, then copied to the host.

**Overlap.** The iterations are queued without waiting for the GPU (lambda stays on the GPU), so
``overlap_fcn`` runs on the CPU while the GPU computes; then the post-processing and the copy.

**Errors.** Out of GPU memory raises ``GPUOutOfMemoryError`` (a ``MemoryError``) with MATLAB's
identifiers (``parallel:gpu:array:OOM``, ``parallel:gpu:array:OOMForOperation``,
``parallel:gpu:array:FFTPlanInternalError``), the ones XR_RLdeconFrame3D catches to deconvolve
the file on the CPU instead; the arrays of the failed call are released first. Boxes outside the
result raise MATLAB's gpuArray indexing errors; integer data other than uint16 MATLAB's
``MATLAB:mixedClasses`` of its GPU loop.
"""

from __future__ import annotations

import math
import threading
import traceback
from collections.abc import Callable, Sequence
from typing import Any, NamedTuple

import numpy as np

__all__ = [
    "MemoryNeed",
    "clear_otf_cache",
    "decon_lucy",
    "decon_lucy_omw",
    "matlab_gpu_sum",
    "memory_need",
    "sum_plan",
    "use_gpu_kernels",
]

_OOM = "parallel:gpu:array:OOM"
_OOM_PLAN = "parallel:gpu:array:OOMForOperation"
_PLAN_FAILED = "parallel:gpu:array:FFTPlanInternalError"
_EXCEEDS = "parallel:array:SubsrefSubscriptExceedsDim"
_BAD_SUBSCRIPT = "parallel:array:SubsrefInvalidSubscript"
_MIXED = "MATLAB:mixedClasses"
_THREADS = 256
_TILE, _ROWS = 32, 8
_MAX_GRID_Z = 65535
_SLACK = 64 << 20  # cuFFT's own buffers and the allocator's rounding

# ---------------------------------------------------------------------------------------------
# the CUDA kernels: MATLAB's gpuArray arithmetic (compiled with --fmad=false, no flush to zero)
# ---------------------------------------------------------------------------------------------

_SOURCE = r"""
typedef long long idx_t;
#define TILE 32
#define ROWS 8
#define EPS_F 2.220446049250313e-16f
#define LOOP(i, n) for (idx_t i = blockIdx.x * (idx_t)blockDim.x + threadIdx.x; i < (n); \
                        i += (idx_t)blockDim.x * gridDim.x)

// max(v, 0) of a gpuArray: -0 and NaN give +0
__device__ __forceinline__ float pos(const float v) { return v > 0.0f ? v : 0.0f; }
__device__ __forceinline__ double posd(const double v) { return v > 0.0 ? v : 0.0; }
// real(ifftn(s)): the real part divided by N in double, rounded to single
__device__ __forceinline__ float re_n(const float2 s,
                                      const double n) { return (float)((double)s.x / n); }
// cast(I, 'like', J) for single J
__device__ __forceinline__ float as_f(const float v) { return v; }
__device__ __forceinline__ float as_f(const unsigned short v) { return (float)v; }
__device__ __forceinline__ float as_f(const double v) { return (float)v; }
// MATLAB's min: a NaN gives the other value, a tie the first
__device__ __forceinline__ float mmin(const float a, const float b) {
    if (a != a) return b;
    if (b != b) return a;
    return b < a ? b : a;
}
// uint16 of a single gpuArray: half away from zero, saturated, NaN 0
__device__ __forceinline__ unsigned short to_u16(const float v) {
    if (!(v > 0.0f)) return 0;
    if (v >= 65535.0f) return 65535;
    return (unsigned short)roundf(v);
}
__device__ __forceinline__ void put(float* o, const float v) { *o = v; }
__device__ __forceinline__ void put(unsigned short* o, const float v) { *o = to_u16(v); }

// ---- the data in MATLAB's layout (z, x, y), from (z, y, x) through 32 x 32 tiles of each plane,
// background subtracted: uint16 I - uint16(bg) saturating, single max(single(double(I) - bg), 0),
// double max(I - bg, 0)
__device__ __forceinline__ unsigned short prep_v(const unsigned short v, const double bg,
                                                 const unsigned short b) {
    return v > b ? (unsigned short)(v - b) : (unsigned short)0;
}
__device__ __forceinline__ float prep_v(const float v, const double bg, const unsigned short b) {
    return pos((float)((double)v - bg));
}
__device__ __forceinline__ double prep_v(const double v, const double bg, const unsigned short b) {
    return posd(v - bg);
}
template <typename T>
__global__ void prep_t(T* __restrict__ out, const T* __restrict__ in, const int ny, const int nx,
                       const int nz,
                       const double bg, const unsigned short b16, const int sub) {
    __shared__ T tile[TILE][TILE + 1];
    const int bx = blockIdx.x * TILE, by = blockIdx.y * TILE;
    for (int z = blockIdx.z; z < nz; z += gridDim.z) {
        const idx_t plane = (idx_t)z * ny * nx;
        for (int r = threadIdx.y; r < TILE; r += ROWS) {
            const int y = by + r, x = bx + threadIdx.x;
            if (y < ny && x < nx) {
                const T v = in[plane + (idx_t)y * nx + x];
                tile[r][threadIdx.x] = sub ? prep_v(v, bg, b16) : v;
            }
        }
        __syncthreads();
        for (int r = threadIdx.y; r < TILE; r += ROWS) {
            const int x = bx + r, y = by + threadIdx.x;
            if (y < ny && x < nx) out[plane + (idx_t)x * ny + y] = tile[threadIdx.x][r];
        }
        __syncthreads();
    }
}
template <typename T>
__global__ void prep_l(T* __restrict__ a, const idx_t n, const double bg,
                       const unsigned short b16) {
    LOOP(i, n) a[i] = prep_v(a[i], bg, b16);
}

// single(I) of uint16 data, a copy of single data, double data rounded
template <typename T>
__global__ void to_single(float* __restrict__ out, const T* __restrict__ in, const idx_t n) {
    LOOP(i, n) out[i] = as_f(in[i]);
}

// PSF ./ sum(PSF(:), 'double'): single(double(PSF) / s)
extern "C" __global__ void normalize(float* __restrict__ out, const float* __restrict__ in,
                                     const idx_t n,
                                     const double* __restrict__ s) {
    const double d = s[0];
    LOOP(i, n) out[i] = (float)((double)in[i] / d);
}

// decon_psf2otf's PSF, padded at the end and shifted by -floor(size / 2), as complex(psf, 0) in
// MATLAB's layout (the spectrum is zero elsewhere)
extern "C" __global__ void scatter_psf(float2* __restrict__ out, const float* __restrict__ psf,
                                       const int pz,
                                       const int py, const int px, const int nz, const int ny,
                                               const int nx) {
    const idx_t np = (idx_t)pz * py * px;
    LOOP(i, np) {
        const int x = (int)(i % px), y = (int)((i / px) % py), z = (int)(i / ((idx_t)px * py));
        const int zo = ((z - pz / 2) % nz + nz) % nz, yo = ((y - py / 2) % ny + ny) % ny,
                  xo = ((x - px / 2) % nx + nx) % nx;
        out[((idx_t)zo * nx + xo) * ny + yo] = make_float2(psf[i], 0.0f);
    }
}

// fftn's input of real data: complex(v, 0)
extern "C" __global__ void spec_real(float2* __restrict__ s, const float* __restrict__ v,
                                     const idx_t n) {
    LOOP(i, n) s[i] = make_float2(v[i], 0.0f);
}

// s = h .* s (conj: conj(h) .* s), MATLAB's complex product: each product and sum rounded
extern "C" __global__ void cmul(float2* __restrict__ s, const float2* __restrict__ h,
                                const idx_t n,
                                const int conj) {
    LOOP(i, n) {
        const float2 a = s[i], b = h[i];
        const float bi = conj ? -b.y : b.y;
        s[i] = make_float2(__fsub_rn(__fmul_rn(b.x, a.x), __fmul_rn(bi, a.y)),
                           __fadd_rn(__fmul_rn(b.x, a.y), __fmul_rn(bi, a.x)));
    }
}
// a single OTF times a double spectrum: in double, rounded to single (into s)
extern "C" __global__ void cmul_d(float2* __restrict__ s, const double2* __restrict__ sd,
                                  const float2* __restrict__ h, const idx_t n) {
    LOOP(i, n) {
        const double2 a = sd[i];
        const double br = (double)h[i].x, bi = (double)h[i].y;
        s[i] = make_float2((float)__dsub_rn(__dmul_rn(br, a.x), __dmul_rn(bi, a.y)),
                           (float)__dadd_rn(__dmul_rn(br, a.y), __dmul_rn(bi, a.x)));
    }
}

// complex(I ./ max(real(ifftn(s)), eps), 0) in place: the ratio, the next fftn's input; cast(I,
// 'like', J) single (ratio), or double I with a double J (ratio_d: in double, rounded)
template <typename T>
__global__ void ratio(float2* __restrict__ s, const T* __restrict__ im, const idx_t n,
                      const double nn) {
    LOOP(i, n) {
        float r = re_n(s[i], nn);
        r = r > EPS_F ? r : EPS_F;
        s[i] = make_float2(__fdiv_rn(as_f(im[i]), r), 0.0f);
    }
}
extern "C" __global__ void ratio_d(float2* __restrict__ s, const double* __restrict__ im,
                                   const idx_t n,
                                   const double nn) {
    LOOP(i, n) {
        float r = re_n(s[i], nn);
        r = r > EPS_F ? r : EPS_F;
        s[i] = make_float2((float)__ddiv_rn(im[i], (double)r), 0.0f);
    }
}

// OMW: J = max(real(ifftn(s)) .* J, 0) (omw_update_d: double J, the new J single)
extern "C" __global__ void omw_update(float* __restrict__ j, const float2* __restrict__ s,
                                      const idx_t n,
                                      const double nn) {
    LOOP(i, n) j[i] = pos(__fmul_rn(re_n(s[i], nn), j[i]));
}
extern "C" __global__ void omw_update_d(float* __restrict__ j, const double* __restrict__ jd,
                                        const float2* __restrict__ s, const idx_t n,
                                                const double nn) {
    LOOP(i, n) j[i] = pos((float)__dmul_rn((double)re_n(s[i], nn), jd[i]));
}
// RL: J_2 = max(Y .* real(ifftn(s)), 0) (rl_update_d: double Y)
extern "C" __global__ void rl_update(float* __restrict__ jn, const float* __restrict__ y,
                                     const float2* __restrict__ s, const idx_t n,
                                             const double nn) {
    LOOP(i, n) jn[i] = pos(__fmul_rn(y[i], re_n(s[i], nn)));
}
extern "C" __global__ void rl_update_d(float* __restrict__ jn, const double* __restrict__ y,
                                       const float2* __restrict__ s, const idx_t n,
                                               const double nn) {
    LOOP(i, n) jn[i] = pos((float)__dmul_rn(y[i], (double)re_n(s[i], nn)));
}

// RL's prediction Y = max(J_2 + lambda * (J_2 - J_3), 0) (spec: also complex(Y, 0) into s);
// predict_fd: single J_2 and double J_3 (double data's second iteration; lambda is 0 there)
extern "C" __global__ void predict(float* __restrict__ y, float2* __restrict__ s,
                                   const float* __restrict__ j2,
                                   const float* __restrict__ j3, const float* __restrict__ lam,
                                           const idx_t n,
                                   const int spec) {
    const float l = *lam;
    LOOP(i, n) {
        const float v = pos(__fadd_rn(j2[i], __fmul_rn(l, __fsub_rn(j2[i], j3[i]))));
        y[i] = v;
        if (spec) s[i] = make_float2(v, 0.0f);
    }
}
extern "C" __global__ void predict_fd(float* __restrict__ y, float2* __restrict__ s,
                                      const float* __restrict__ j2,
                                      const double* __restrict__ j3, const idx_t n) {
    LOOP(i, n) {
        const float d = (float)__dsub_rn((double)j2[i], j3[i]);
        const float v = pos(__fadd_rn(j2[i], (float)__dmul_rn(0.0, (double)d)));
        y[i] = v;
        s[i] = make_float2(v, 0.0f);
    }
}
// J_4 = J_2(:) - Y(:) (diff_fd: double Y, in double and rounded)
extern "C" __global__ void diff(float* __restrict__ j4, const float* __restrict__ j2,
                                const float* __restrict__ y,
                                const idx_t n) {
    LOOP(i, n) j4[i] = __fsub_rn(j2[i], y[i]);
}
extern "C" __global__ void diff_fd(float* __restrict__ j4, const float* __restrict__ j2,
                                   const double* __restrict__ y, const idx_t n) {
    LOOP(i, n) j4[i] = (float)__dsub_rn((double)j2[i], y[i]);
}
// lambda = max(min(s0 / single(double(s1) + eps), 1), 0) (min of NaN and 1 is 1)
extern "C" __global__ void lambda_from_sums(float* __restrict__ lam,
                                            const float* __restrict__ sums) {
    const float d = (float)__dadd_rn((double)sums[1], 2.220446049250313e-16);
    float l = __fdiv_rn(sums[0], d);
    l = l < 1.0f ? l : 1.0f;
    *lam = l > 0.0f ? l : 0.0f;
}

// ---- MATLAB's gpuArray sum (see matlab_gpu_sum). The element k of a sum: K 0 x[k] (in the
// accumulator's class), K 1 (a[k] - b[k]) .* c[k] and K 2 c[k] .* c[k] (single products, lambda's)
template <typename T> __device__ __forceinline__ T add(const T a, const T b);
template <> __device__ __forceinline__ float add(const float a,
                                                 const float b) { return __fadd_rn(a, b); }
template <> __device__ __forceinline__ double add(const double a,
                                                  const double b) { return __dadd_rn(a, b); }
template <typename T, typename TI, int K>
__device__ __forceinline__ T elem(const TI* a, const float* b, const float* c, const idx_t k) {
    if (K == 1) return (T)__fmul_rn(__fsub_rn((float)a[k], b[k]), c[k]);
    if (K == 2) return (T)__fmul_rn(c[k], c[k]);
    return (T)a[k];
}
// n <= 128: one block of 128 threads, one element each, a halving tree in shared memory
template <typename T, typename TI, int K>
__global__ void small_sum(T* __restrict__ out, const TI* __restrict__ a,
                          const float* __restrict__ b,
                          const float* __restrict__ c, const idx_t n) {
    __shared__ T v[128];
    const int t = threadIdx.x;
    v[t] = t < n ? elem<T, TI, K>(a, b, c, t) : (T)0;
    __syncthreads();
    for (int s = 64; s > 0; s >>= 1) {
        if (t < s) v[t] = add(v[t], v[t + s]);
        __syncthreads();
    }
    if (t == 0) out[0] = v[0];
}
// chunks of `chunk` elements, one block of 256 threads each: thread t adds the chunk's elements
// t, t + 256, ... in order; then each warp adds its lanes by shuffles at offsets 1, 2, 4, 8, 16
// and thread 0 the warp sums in order (CUB's block reduction); out[chunk]
template <typename T>
__device__ __forceinline__ T warp_sum(T v, const int last_lane) {
    const int lane = threadIdx.x & 31;
#pragma unroll
    for (int off = 1; off < 32; off <<= 1) {
        const T o = __shfl_down_sync(0xffffffffu, v, off);
        if (lane + off <= last_lane) v = add(v, o);
    }
    return v;
}
template <typename T, typename TI, int K>
__global__ void chunk_sum(T* __restrict__ out, const TI* __restrict__ a,
                          const float* __restrict__ b, const float* __restrict__ c,
                          const idx_t n, const idx_t chunk) {
    __shared__ T wsum[8];
    const int t = threadIdx.x;
    const idx_t beg = (idx_t)blockIdx.x * chunk;
    const idx_t end = beg + chunk < n ? beg + chunk : n;
    idx_t k = beg + t;
    T agg = (T)0;
    if (k < end) {
        agg = elem<T, TI, K>(a, b, c, k);
        for (k += 256; k < end; k += 256) agg = add(agg, elem<T, TI, K>(a, b, c, k));
    }
    const idx_t left = end - beg;
    const int num_valid = left >= 256 ? 256 : (int)left;
    const int warp = t >> 5, lane = t & 31;
    const T v = warp_sum(agg, (warp * 32 + 32 <= num_valid ? 32 : num_valid - warp * 32) - 1);
    if (lane == 0) wsum[warp] = v;
    __syncthreads();
    if (t == 0) {
        T r = wsum[0];
        for (int w = 1; w < 8; ++w)
            if (w * 32 < num_valid) r = add(r, wsum[w]);
        out[blockIdx.x] = r;
    }
}

// the debug error metric's elements (double(J) - double(I)) .^ 2
template <typename T>
__global__ void sqdiff(double* __restrict__ out, const float* __restrict__ j,
                       const T* __restrict__ im,
                       const idx_t n) {
    LOOP(i, n) {
        const double d = __dsub_rn((double)j[i], (double)im[i]);
        out[i] = __dmul_rn(d, d);
    }
}

// ---- post-processing, uint16 and the box: dampFactor, scaleFactor and deconOffset (single with
// double: in double, rounded), the eroded mask, only the box written. J and I in MATLAB's layout
// (post_t: tiles of each z plane) or (z, y, x) (post_l); the mask and the output (z, y, x)
#define F_DAMP 1
#define F_SCALE 2
#define F_OFFSET 4
#define F_MASK 8
template <typename TI>
__device__ __forceinline__ float post_v(float v, const TI* im, const idx_t k, const double damp,
                                        const double scale, const double offset, const int flags) {
    if (flags & F_DAMP) {
        const float di = (float)__dmul_rn(damp, (double)as_f(im[k]));
        const float lo = __fmul_rn(v, v < di ? 1.0f : 0.0f);
        const float hi = __fmul_rn(di, v >= di ? 1.0f : 0.0f);
        v = __fadd_rn(lo, mmin(v, hi));
    }
    if (flags & F_SCALE) {
        v = (float)__dmul_rn((double)v, scale);
        if (flags & F_OFFSET) v = (float)__dadd_rn((double)v, offset);
    }
    return v;
}
template <typename TI, typename TO>
__global__ void post_t(TO* __restrict__ out, const float* __restrict__ j,
                       const TI* __restrict__ im,
                       const bool* __restrict__ mask, const int ny, const int nx, const int z0,
                               const int y0,
                       const int x0, const int oz, const int oy, const int ox, const double damp,
                       const double scale, const double offset, const int flags) {
    __shared__ float tile[TILE][TILE + 1];
    const int bx = blockIdx.x * TILE, by = blockIdx.y * TILE;
    for (int zz = blockIdx.z; zz < oz; zz += gridDim.z) {
        const idx_t plane = (idx_t)(z0 + zz) * nx * ny;
        for (int r = threadIdx.y; r < TILE; r += ROWS) {
            const int y = by + threadIdx.x, x = bx + r;
            if (y < oy && x < ox) {
                const idx_t k = plane + (idx_t)(x0 + x) * ny + (y0 + y);
                tile[r][threadIdx.x] = post_v(j[k], im, k, damp, scale, offset, flags);
            }
        }
        __syncthreads();
        for (int r = threadIdx.y; r < TILE; r += ROWS) {
            const int y = by + r, x = bx + threadIdx.x;
            if (y < oy && x < ox) {
                float v = tile[threadIdx.x][r];
                if (flags & F_MASK) v = __fmul_rn(v,
                    mask[plane + (idx_t)(y0 + y) * nx + (x0 + x)] ? 1.0f : 0.0f);
                put(out + ((idx_t)zz * oy + y) * ox + x, v);
            }
        }
        __syncthreads();
    }
}
template <typename TI, typename TO>
__global__ void post_l(TO* __restrict__ out, const float* __restrict__ j,
                       const TI* __restrict__ im,
                       const bool* __restrict__ mask, const int ny, const int nx, const int z0,
                               const int y0,
                       const int x0, const int oz, const int oy, const int ox, const double damp,
                       const double scale, const double offset, const int flags) {
    const idx_t no = (idx_t)oz * oy * ox;
    LOOP(i, no) {
        const int x = (int)(i % ox), y = (int)((i / ox) % oy), z = (int)(i / ((idx_t)ox * oy));
        const idx_t k = ((idx_t)(z0 + z) * ny + (y0 + y)) * nx + (x0 + x);
        float v = post_v(j[k], im, k, damp, scale, offset, flags);
        if (flags & F_MASK) v = __fmul_rn(v, mask[k] ? 1.0f : 0.0f);
        put(out + i, v);
    }
}
"""

_U16, _F32, _F64 = "unsigned short", "float", "double"
_CT = {np.dtype(np.uint16): _U16, np.dtype(np.float32): _F32, np.dtype(np.float64): _F64}
_SUMS = (
    (_F32, _F32, 0),
    (_F32, _F32, 1),
    (_F32, _F32, 2),
    (_F64, _F32, 0),
    (_F64, _F64, 0),
)
_NAMES = (
    *(f"prep_t<{t}>" for t in (_U16, _F32, _F64)),
    *(f"prep_l<{t}>" for t in (_U16, _F32, _F64)),
    *(f"to_single<{t}>" for t in (_U16, _F32, _F64)),
    *(f"ratio<{t}>" for t in (_U16, _F32, _F64)),
    *(f"small_sum<{a}, {b}, {k}>" for a, b, k in _SUMS),
    *(f"chunk_sum<{a}, {b}, {k}>" for a, b, k in _SUMS),
    *(f"sqdiff<{t}>" for t in (_U16, _F32, _F64)),
    *(f"post_{v}<{a}, {b}>" for v in ("t", "l") for a in (_U16, _F32, _F64) for b in (_F32, _U16)),
)

_modules: dict[int, Any] = {}
_modules_lock = threading.Lock()


def _cp() -> Any:
    from petakit5d.gpu._cupy import import_cupy

    return import_cupy()


def _mod() -> Any:
    """the compiled kernels on the current GPU"""
    cp = _cp()
    dev = cp.cuda.Device().id
    with _modules_lock:
        m = _modules.get(dev)
        if m is None:
            from petakit5d.gpu import kernels

            m = kernels.load_source(_SOURCE, name="rl_gpu.cu", name_expressions=_NAMES)
            _modules[dev] = m
        return m


def _grid(n: int) -> int:
    """blocks of 256 threads for grid-stride loops over n elements"""
    return max(1, min(-(-int(n) // _THREADS), 1 << 20))


def _tiles(ny: int, nx: int, nz: int) -> tuple[int, int, int]:
    return (-(-nx // _TILE), -(-ny // _TILE), max(1, min(nz, _MAX_GRID_Z)))


_SMS: dict[int, int] = {}


def _sms() -> int:
    """the current GPU's multiprocessors"""
    cp = _cp()
    dev = cp.cuda.Device().id
    if dev not in _SMS:
        _SMS[dev] = int(cp.cuda.runtime.deviceGetAttribute(16, dev))  # MultiProcessorCount
    return _SMS[dev]


# ---------------------------------------------------------------------------------------------
# MATLAB's gpuArray sum
# ---------------------------------------------------------------------------------------------


def sum_plan(n: int, sms: int) -> tuple[int, int]:
    """(chunks, chunk length) of MATLAB's gpuArray sum of n > 128 elements on a GPU with ``sms``
    multiprocessors: one chunk up to 8192 elements, else chunks of 1408 * k elements with k the
    smallest that gives at most SMs * ceil(1024 / SMs) chunks (R04's probes: exact on 80, 108 and
    148 multiprocessors, Titan V, A100 and B200; other counts are the model's prediction)"""
    if n <= 8192:
        return 1, int(n)
    sms = max(1, int(sms))
    gmax = sms * -(-1024 // sms)
    c = 1408 * -(-n // (1408 * gmax))
    return -(-n // c), c


def _reduce(kind: str, n: int, a: Any, b: Any, c: Any, out: Any, *, double: bool = False) -> None:
    """MATLAB's gpuArray sum of n elements into out[0] (on the GPU, no wait): kind 'x' a(:) (with
    double: sum(a(:), 'double')), 'num' sum((a(:) - b(:)) .* c), 'sq' sum(c .* c)"""
    cp = _cp()
    m = _mod()
    acc = _F64 if double else _F32
    src = _CT[a.dtype]
    k = {"x": 0, "num": 1, "sq": 2}[kind]
    b = a if b is None else b
    c = a if c is None else c
    # (b and c are float pointers; they are only read for kinds 1 and 2)
    if b.dtype != cp.float32:
        b = c = cp.empty(1, dtype=cp.float32)
    if n <= 128:
        m[f"small_sum<{acc}, {src}, {k}>"](1, 128, (out, a, b, c, np.int64(n)))
        return
    g, length = sum_plan(n, _sms())
    fn = m[f"chunk_sum<{acc}, {src}, {k}>"]
    if g == 1:
        fn(1, 256, (out, a, b, c, np.int64(n), np.int64(length)))
        return
    parts = cp.empty(g, dtype=cp.float64 if double else cp.float32)
    fn(g, 256, (parts, a, b, c, np.int64(n), np.int64(length)))
    _reduce("x", g, parts, None, None, out, double=double)


def matlab_gpu_sum(x: Any, *, double: bool = False) -> Any:
    """MATLAB's gpuArray ``sum(x(:))`` of a single CuPy array (``double``: ``sum(x(:),
    'double')``, each element converted to double first), as a 1-element CuPy array, in MATLAB's
    order (R04's probes of MATLAB R2024b on Titan V, A100 and B200, every addition one IEEE add):

    - up to 128 elements: a halving tree over 128 lanes (lane t holds element t);
    - up to 8192: one block of 256 threads, thread t adding elements t, t + 256, ... in order; the
      lanes of each warp added by shuffles at offsets 1, 2, 4, 8, 16 (lane l takes lane l + off
      while that is a lane with elements), the warp sums in order (CUB's block reduction);
    - more: contiguous chunks (``sum_plan``), each summed as one block, the chunk sums by these
      rules again.
    """
    cp = _cp()
    a = cp.ascontiguousarray(x).reshape(-1)
    if a.dtype != cp.float32:
        raise ValueError("matlab_gpu_sum: a single array")
    out = cp.zeros(1, dtype=cp.float64 if double else cp.float32)
    if a.size:
        _reduce("x", a.size, a, None, None, out, double=double)
    return out


# ---------------------------------------------------------------------------------------------
# path, shapes, arguments
# ---------------------------------------------------------------------------------------------


def use_gpu_kernels(
    use_mex: bool | None, method: str, data: np.ndarray, debug: bool, background: Any
) -> tuple[bool, bool]:
    """(use, can): decon_use_rfft_mex on the GPU. ``can``: uint16 or single, real, non-empty 2D or
    3D data, a real scalar background, no debug output. None runs MATLAB's loop (MATLAB's default
    on the GPU), True the half-spectrum kernels (else MATLAB's error), False the loop."""
    if method not in ("simplified", "omw"):
        raise ValueError(f"unknown method {method!r}")
    scalar_bg = background is None or (
        np.size(background) == 1 and not np.iscomplexobj(np.asarray(background))
    )
    can = (
        data.dtype in (np.uint16, np.float32)
        and data.size > 0
        and data.ndim in (2, 3)
        and scalar_bg
        and not debug
    )
    if use_mex is None:
        return False, can
    if use_mex and not can:
        from petakit5d.core.arrays import matlab_class, size_to_matlab
        from petakit5d.decon.rl import DeconError

        msize = "[" + " ".join(str(v) for v in size_to_matlab(data.shape)) + "]"
        raise DeconError(
            "",
            f"decon_lucy_rfft_gpu_mex does not support the data (class {matlab_class(data)}, size "
            f"{msize}) or the options (debug on the GPU or a non-scalar background).",
        )
    return bool(use_mex), can


class _Geom(NamedTuple):
    """data of shape (nz, ny, nx) or (ny, nx); MATLAB's layout is the C array mshape (z, x, y)"""

    shape: tuple[int, ...]
    zyx: tuple[int, int, int]
    mshape: tuple[int, ...]

    @property
    def n(self) -> int:
        return math.prod(self.zyx)


def _geom(shape: Sequence[int]) -> _Geom:
    s = tuple(int(v) for v in shape)
    nz, ny, nx = (1, *s) if len(s) == 2 else s
    return _Geom(s, (nz, ny, nx), (nz, nx, ny) if len(s) == 3 else (nx, ny))


def _scalar(v: Any) -> float:
    return float(np.asarray(v).reshape(-1)[0])


def _box(bbox: Any, g: _Geom) -> tuple[int, int, int, int, int, int]:
    """(z0, y0, x0, oz, oy, ox) of ``J(bbox(1) : bbox(4), bbox(2) : bbox(5), bbox(3) : bbox(6))``
    on the GPU: a reversed range is empty, a nonempty range starting before the first element or
    ending past the last raises MATLAB's gpuArray errors (a 2D result's third range must be
    1 : 1)"""
    nz, ny, nx = g.zyx
    if bbox is None or len(bbox) == 0:
        return 0, 0, 0, nz, ny, nx
    from petakit5d.decon.rl import DeconError

    b = list(bbox)
    if len(b) == 4:
        b = [0, b[0], b[1], 1, b[2], b[3]]
    if len(b) != 6:
        raise ValueError(f"bbox must have 6 values (z0, y0, x0, z1, y1, x1), not {len(b)}")
    starts, ends = [], []
    for s, e, size in zip(b[:3], b[3:], (nz, ny, nx), strict=True):
        starts.append(int(s))
        ends.append(size if e is None else int(e))
    # MATLAB's dimensions in its order: y, x, z
    for k in (1, 2, 0):
        s, e, size = starts[k], ends[k], (nz, ny, nx)[k]
        if e <= s:
            continue
        if s < 0:
            raise DeconError(
                _BAD_SUBSCRIPT,
                "Subscript indices must either be real positive integers or logicals.",
            )
        if e > size:
            raise DeconError(_EXCEEDS, "Index exceeds matrix dimension.")
    ext = [max(0, e - s) for s, e in zip(starts, ends, strict=True)]
    return starts[0], starts[1], starts[2], ext[0], ext[1], ext[2]


# ---------------------------------------------------------------------------------------------
# the OTF cache (MATLAB's persistent variables, on the GPU)
# ---------------------------------------------------------------------------------------------


class _Entry(NamedTuple):
    key: tuple  # (data shape, mex or loop, layout, GPU)
    psfs: tuple[np.ndarray, ...]
    otfs: tuple[Any, ...]


_CACHE: dict[str, _Entry] = {}
_CACHE_LOCK = threading.Lock()


def clear_otf_cache() -> None:
    """Forget the OTFs cached on the GPU (MATLAB's ``clear decon_lucy_function
    decon_lucy_omw_function``)."""
    with _CACHE_LOCK:
        _CACHE.clear()


def _cached_otf_bytes(method: str) -> int:
    with _CACHE_LOCK:
        e = _CACHE.get(method)
        return sum(int(o.nbytes) for o in e.otfs) if e is not None else 0


def _normalized(p: np.ndarray) -> Any:
    """``PSF ./ sum(PSF(:), 'double')`` on the GPU"""
    cp = _cp()
    d = cp.asarray(p)
    s = matlab_gpu_sum(d, double=True)
    out = cp.empty_like(d)
    _mod()["normalize"](_grid(d.size), _THREADS, (out, d, np.int64(d.size), s))
    return out


def _full_otf(psf: Any, g: _Geom, plan: Any) -> Any:
    """decon_psf2otf on the GPU: MATLAB's full OTF, complex64 in MATLAB's layout (an all-zero PSF:
    zeros)"""
    cp = _cp()
    from petakit5d.decon.rl import DeconError

    ps = (1,) * (3 - psf.ndim) + tuple(int(v) for v in psf.shape)
    if any(p > d for p, d in zip(ps, g.zyx, strict=True)):
        raise DeconError("", "outSizeIsSmallerThanPsfSize")
    h = cp.zeros(g.mshape, dtype=cp.complex64)
    if not np.any(psf.get() != 0):  # (the PSF is small: checked on the host)
        return h
    nz, ny, nx = g.zyx
    _mod()["scatter_psf"](
        _grid(psf.size),
        _THREADS,
        (h, psf, *(np.int32(v) for v in ps), *(np.int32(v) for v in (nz, ny, nx))),
    )
    plan.forward(h, h)
    return h


def _otfs(method, psfs, g, use_mex, layout, keep_plans):
    """the OTFs on the GPU, cached for the same PSF bits, data size, path, layout and GPU"""
    cp = _cp()
    key = (g.shape, use_mex, layout if use_mex else "matlab", cp.cuda.Device().id)
    with _CACHE_LOCK:
        e = _CACHE.get(method)
        if e is not None and e.key == key and len(e.psfs) == len(psfs):
            if all(
                a.shape == b.shape and a.tobytes() == b.tobytes()
                for a, b in zip(e.psfs, psfs, strict=True)
            ):
                return e.otfs
        _CACHE.pop(method, None)
    norm = [_normalized(p) for p in psfs]
    if use_mex:
        from petakit5d.kernels.decon_gpu import LucyRFFTGPU

        k = LucyRFFTGPU(g.shape, layout=layout, keep_plans=keep_plans)
        otfs = tuple(k.otf(p) for p in norm)
    else:
        plan = _complex_plan(g, keep_plans)
        try:
            otfs = tuple(_full_otf(p, g, plan) for p in norm)
        finally:
            if not keep_plans:
                plan.close()
    with _CACHE_LOCK:
        _CACHE[method] = _Entry(key, tuple(p.copy() for p in psfs), otfs)
    return otfs


def _complex_plan(g: _Geom, keep: bool) -> Any:
    from petakit5d.gpu import fft

    return fft.get_complex_plan(g.mshape) if keep else fft.ComplexFFTPlan(g.mshape)


# ---------------------------------------------------------------------------------------------
# memory
# ---------------------------------------------------------------------------------------------


class MemoryNeed(NamedTuple):
    """GPU bytes of a decon call at its peak (``total``), of which ``otfs`` (kept between calls of
    the same PSFs and size) and ``plan`` (the FFT work area, kept between calls with keep_plans);
    ``per_voxel``: total / voxels."""

    total: int
    otfs: int
    plan: int
    per_voxel: float


def memory_need(
    method: str,
    shape: Sequence[int],
    dtype: Any = np.uint16,
    *,
    use_mex: bool = False,
    mask: bool = False,
    plan_bytes: int | None = None,
) -> MemoryNeed:
    """The GPU memory of a decon call (``method`` 'omw' or 'simplified') of data of ``shape`` and
    ``dtype`` at its peak, from R04's buffers (measured against the peaks, see the R04 report):
    the loop on full spectra (default) or the half spectra (``use_mex``); ``mask`` for deconOffset
    or EdgeErosion (one byte per voxel). ``plan_bytes``: the FFT plan's work area when known (else
    one spectrum: cuFFT's work area of C2C plans at FFT-friendly sizes; sizes with large prime
    factors take up to twice that)."""
    if method not in ("simplified", "omw"):
        raise ValueError(f"unknown method {method!r}")
    g = _geom(shape)
    n = g.n
    dt = np.dtype(dtype)
    im = 2 if dt == np.uint16 else (8 if dt == np.float64 else 4)
    if use_mex:
        nh = math.prod((*g.mshape[:-1], g.mshape[-1] // 2 + 1)) * 8
        work = nh if plan_bytes is None else int(plan_bytes)
        otfs = nh * (2 if method == "omw" else 1)
        # OMW: the data, J, a real scratch and a half spectrum; RL: I, J_2, J_3, J_4, Y and
        # rlstep's result and half spectrum
        peak = (n * im + 8 * n + nh) if method == "omw" else (n * im + 20 * n + nh)
    else:
        work = 8 * n if plan_bytes is None else int(plan_bytes)
        otfs = 8 * n * (2 if method == "omw" else 1)
        dbl = dt == np.float64
        spec = 8 * n + (16 * n if dbl else 0)  # double data: a double spectrum in the first FFT
        # OMW: the data, J and a spectrum; RL: the data, J_2, J_3, J_4, Y (double data: J_2, J_3
        # and Y double in the first iteration) and a spectrum
        peak = (
            (n * im + 4 * n + spec)
            if method == "omw"
            else (n * im + (16 + (8 if dbl else 0)) * n + spec)
        )
    peak += n if mask else 0
    total = peak + otfs + work
    return MemoryNeed(int(total), int(otfs), int(work), total / n)


def _check_memory(need: MemoryNeed, cached: int) -> None:
    """GPUOutOfMemoryError (parallel:gpu:array:OOM) before anything is allocated when the GPU
    lacks the call's memory (cached: bytes of the call already on the GPU)"""
    from petakit5d.gpu.memory import available_bytes
    from petakit5d.kernels.decon_gpu import GPUOutOfMemoryError

    avail = available_bytes()
    more = need.total - cached
    if more + _SLACK > avail:
        raise GPUOutOfMemoryError(
            _OOM,
            f"Out of memory on device: the deconvolution needs {more / 2**30:.2f} GiB more on the "
            f"GPU, {avail / 2**30:.2f} GiB are free.",
        )


def _run(fn: Callable[..., Any], *args: Any, **kw: Any) -> Any:
    """fn(*args, **kw), with out-of-memory errors as MATLAB's GPU code gives them; the failed call's
    arrays are released before the error reaches the caller"""
    cp = _cp()
    from petakit5d.gpu.fft import CuFFTMemoryError
    from petakit5d.kernels.decon_gpu import GPUOutOfMemoryError

    err: Exception | None = None
    try:
        return fn(*args, **kw)
    except cp.cuda.memory.OutOfMemoryError as e:
        err = GPUOutOfMemoryError(_OOM, f"Out of memory on device. {e}")
        traceback.clear_frames(e.__traceback__)
    except CuFFTMemoryError as e:
        err = (
            GPUOutOfMemoryError(_PLAN_FAILED, f"cuFFT failed to make a plan ({e}).")
            if e.code == 5
            else GPUOutOfMemoryError(_OOM_PLAN, f"Out of memory on device for the FFT plans ({e}).")
        )
        traceback.clear_frames(e.__traceback__)
    except GPUOutOfMemoryError as e:
        traceback.clear_frames(e.__traceback__)
        err = e
    raise err


# ---------------------------------------------------------------------------------------------
# the data on the GPU
# ---------------------------------------------------------------------------------------------

_CHUNK = 64 << 20  # bytes per pinned staging buffer (two, reused by every call of the process)
_staging: list[Any] = []
_staging_lock = threading.Lock()
_pool: Any = None
_workers = 1


def _buffers() -> tuple[list[np.ndarray], list[int], Any]:
    """the two pinned host buffers (as byte arrays and their addresses) and the copy threads"""
    global _pool, _workers
    cp = _cp()
    with _staging_lock:
        if not _staging:
            for _ in range(2):
                mem = cp.cuda.alloc_pinned_memory(_CHUNK)
                _staging.append((mem, np.frombuffer(mem, np.uint8, _CHUNK)))
        if _pool is None:
            from concurrent.futures import ThreadPoolExecutor

            from petakit5d.core.threads import resolve_threads

            _workers = min(4, resolve_threads(None))
            _pool = ThreadPoolExecutor(_workers)
    return [b for _, b in _staging], [m.ptr for m, _ in _staging], _pool


def _copy(dst: np.ndarray, src: np.ndarray, pool: Any) -> None:
    """dst[:] = src for byte arrays, in pieces on the copy threads"""
    n = dst.size
    k = max(1, min(_workers, n >> 22))
    step = -(-n // k)
    list(pool.map(lambda a: np.copyto(dst[a : a + step], src[a : a + step]), range(0, n, step)))


def _to_gpu(a: np.ndarray) -> Any:
    """a C-contiguous host array on the GPU: through two pinned buffers, the copy of each piece
    into one overlapping the transfer of the other (a large pageable copy is about half as fast)"""
    cp = _cp()
    a = np.ascontiguousarray(a)
    if a.nbytes <= _CHUNK:
        return cp.asarray(a)
    out = cp.empty(a.shape, dtype=a.dtype)
    bufs, ptrs, pool = _buffers()
    src = a.reshape(-1).view(np.uint8)
    stream = cp.cuda.get_current_stream()
    done = [None, None]
    for i, off in enumerate(range(0, a.nbytes, _CHUNK)):
        b, n = i % 2, min(_CHUNK, a.nbytes - off)
        if done[b] is not None:
            done[b].synchronize()
        _copy(bufs[b][:n], src[off : off + n], pool)
        cp.cuda.runtime.memcpyAsync(out.data.ptr + off, ptrs[b], n, 1, stream.ptr)
        done[b] = cp.cuda.Event()
        done[b].record(stream)
    stream.synchronize()
    return out


def _to_host(a: Any) -> np.ndarray:
    """a C-contiguous GPU array on the host, through the pinned buffers as _to_gpu"""
    cp = _cp()
    if a.nbytes <= _CHUNK:
        return a.get()
    out = np.empty(a.shape, dtype=a.dtype)
    dst = out.reshape(-1).view(np.uint8)
    bufs, ptrs, pool = _buffers()
    stream = cp.cuda.get_current_stream()
    pending = None  # (buffer, offset, bytes, event) of the piece in flight
    for i, off in enumerate(range(0, a.nbytes, _CHUNK)):
        b, n = i % 2, min(_CHUNK, a.nbytes - off)
        cp.cuda.runtime.memcpyAsync(ptrs[b], a.data.ptr + off, n, 2, stream.ptr)
        ev = cp.cuda.Event()
        ev.record(stream)
        if pending is not None:
            pb, po, pn, pe = pending
            pe.synchronize()
            _copy(dst[po : po + pn], bufs[pb][:pn], pool)
        pending = (b, off, n, ev)
    pb, po, pn, pe = pending
    pe.synchronize()
    _copy(dst[po : po + pn], bufs[pb][:pn], pool)
    return out


def _prep(g: _Geom, raw: Any, bg: Any, scalar: bool, matlab_layout: bool) -> Any:
    """the background-subtracted data on the GPU (MATLAB's layout, or (z, y, x)): uint16 ``I -
    uint16(Background)`` (MATLAB's rounding of the background, saturating), else ``max(I -
    Background, 0)``; ``raw`` is the (z, y, x) data on the GPU (changed in place where it is not
    put in MATLAB's layout)"""
    cp = _cp()
    from petakit5d.core import dtypes

    t = _CT[raw.dtype]
    m = _mod()
    nz, ny, nx = g.zyx
    if scalar:
        b16 = (
            np.uint16(np.asarray(dtypes.uint16(float(bg))).reshape(-1)[0])
            if raw.dtype == cp.uint16
            else np.uint16(0)
        )
        bgv, sub = np.float64(bg), np.int32(1)
    else:
        # a background per voxel (MATLAB's implicit expansion): subtracted here, as MATLAB does
        bgv, b16, sub = np.float64(0), np.uint16(0), np.int32(0)
        b = np.asarray(bg)
        if raw.dtype == cp.uint16:
            bb = cp.asarray(np.broadcast_to(np.asarray(dtypes.uint16(b)), g.shape).astype(np.int32))
            raw = cp.maximum(raw.astype(cp.int32) - bb, 0).astype(cp.uint16)
        else:
            bd = cp.asarray(np.broadcast_to(b.astype(np.float64), g.shape))
            v = raw.astype(cp.float64) - bd
            v = cp.where(v > 0, v, 0)
            raw = v if raw.dtype == cp.float64 else v.astype(cp.float32)
            del bd, v
    if not matlab_layout:
        if sub:
            m[f"prep_l<{t}>"](_grid(raw.size), _THREADS, (raw, np.int64(raw.size), bgv, b16))
        return raw
    out = cp.empty(g.mshape, dtype=raw.dtype)
    m[f"prep_t<{t}>"](
        _tiles(ny, nx, nz),
        (_TILE, _ROWS),
        (out, raw, np.int32(ny), np.int32(nx), np.int32(nz), bgv, b16, sub),
    )
    return out


def _single(a: Any) -> Any:
    """single(a): a new single array"""
    cp = _cp()
    out = cp.empty(a.shape, dtype=cp.float32)
    _mod()[f"to_single<{_CT[a.dtype]}>"](_grid(a.size), _THREADS, (out, a, np.int64(a.size)))
    return out


def _host_zyx(a: Any, matlab_layout: bool) -> np.ndarray:
    """an array on the GPU as a (z, y, x) host array"""
    h = a.get()
    return np.ascontiguousarray(np.swapaxes(h, -1, -2)) if matlab_layout else h


class _DebugGPU:
    """the debug outputs of MATLAB's GPU loops: err_mat with the error metric's gpuArray mean, the
    estimates as uint16 of the gathered single array"""

    def __init__(self, n_iter: int, save_step: Any, folder: str, n: int, matlab_layout: bool):
        from petakit5d.core.threads import resolve_threads
        from petakit5d.decon.rl import _Debug

        self.dbg = _Debug(n_iter, save_step, folder, resolve_threads(None))
        self.scratch = _cp().empty(n, dtype=np.float64)
        self.ml = matlab_layout

    def __call__(self, k: int, j: Any, im: Any) -> None:
        cp = _cp()
        n = j.size
        _mod()[f"sqdiff<{_CT[im.dtype]}>"](_grid(n), _THREADS, (self.scratch, j, im, np.int64(n)))
        out = cp.zeros(1, dtype=cp.float64)
        _reduce("x", n, self.scratch, None, None, out, double=True)
        self.dbg.bookkeeping(k, float(out.get()[0]) / n)
        self.dbg.files(k, lambda: _host_zyx(j, self.ml))


# ---------------------------------------------------------------------------------------------
# MATLAB's gpuArray loops on full spectra
# ---------------------------------------------------------------------------------------------


def _omw_loop(g: _Geom, im: Any, hf: Any, hb: Any, n_iter: int, plan: Any, dbg: Any) -> Any:
    """decon_lucy_omw_function's loop: J (single, MATLAB's layout) after n_iter iterations from J
    = I (single(I) for uint16); I in MATLAB's layout (uint16, single, or double: its first
    iteration in double)"""
    cp = _cp()
    m = _mod()
    n = g.n
    nn, gn, n64, zero = np.float64(n), _grid(n), np.int64(n), np.int32(0)
    if n_iter == 0:
        return im if im.dtype == cp.float64 else _single(im)
    s = cp.empty(g.mshape, dtype=cp.complex64)
    ratio = m[f"ratio<{_CT[im.dtype]}>"]
    first = 1
    if im.dtype == cp.float64:
        # J = I double: fftn(J) double, times the single OTF in double; I ./ CX in double; the new J
        # single (real(CX) .* J in double)
        sd = cp.fft.fftn(im.astype(cp.complex128), axes=tuple(range(im.ndim)))
        m["cmul_d"](gn, _THREADS, (s, sd, hf, n64))
        del sd
        plan.inverse(s, s)
        m["ratio_d"](gn, _THREADS, (s, im, n64, nn))
        plan.forward(s, s)
        m["cmul"](gn, _THREADS, (s, hb, n64, zero))
        plan.inverse(s, s)
        j = cp.empty(g.mshape, dtype=cp.float32)
        m["omw_update_d"](gn, _THREADS, (j, im, s, n64, nn))
        if dbg is not None:
            dbg(1, j, im)
        first = 2
    else:
        j = _single(im)
    for k in range(first, n_iter + 1):
        m["spec_real"](gn, _THREADS, (s, j, n64))
        plan.forward(s, s)
        m["cmul"](gn, _THREADS, (s, hf, n64, zero))
        plan.inverse(s, s)
        ratio(gn, _THREADS, (s, im, n64, nn))
        plan.forward(s, s)
        m["cmul"](gn, _THREADS, (s, hb, n64, zero))
        plan.inverse(s, s)
        m["omw_update"](gn, _THREADS, (j, s, n64, nn))
        if dbg is not None:
            dbg(k, j, im)
    return j


def _lambda(n: int, j2: Any, y: Any, j4: Any, sums: Any, lam: Any) -> None:
    """lambda = max(min(sum((J_2(:) - Y(:)) .* J_4) / (sum(J_4 .* J_4) + eps), 1), 0) on the GPU
    (MATLAB's gpuArray sums), then J_4 = J_2(:) - Y(:)"""
    m = _mod()
    _reduce("num", n, j2, y, j4, sums[0:1])
    _reduce("sq", n, j2, y, j4, sums[1:2])
    m["lambda_from_sums"](1, 1, (lam, sums))
    m["diff"](_grid(n), _THREADS, (j4, j2, y, np.int64(n)))


def _rl_loop(g: _Geom, im: Any, h: Any, n_iter: int, plan: Any, dbg: Any, rlstep: Any = None):
    """decon_lucy_function's loop: J_2 (single; MATLAB's layout, or (z, y, x) with the native
    half spectra) after n_iter iterations; the update on full spectra (plan) or rlstep(Y, I, H)
    (the mex's half spectra); I uint16, single, or double (its first iterations mixed)"""
    cp = _cp()
    m = _mod()
    n = g.n
    nn, gn, n64 = np.float64(n), _grid(n), np.int64(n)
    shape = im.shape
    lam = cp.zeros(1, dtype=cp.float32)
    sums = cp.zeros(2, dtype=cp.float32)
    if n_iter == 0:
        return im if im.dtype == cp.float64 else _single(im)
    s = cp.empty(shape, dtype=cp.complex64) if rlstep is None else None
    first = 1
    if im.dtype == cp.float64:
        # k = 1: Y = J_2 = I double (+0 for -0); fftn(Y) double; H .* (.) in double; I ./ (.) in
        # double; the new J_2 = max(Y .* real(.), 0) in double, rounded
        y64 = cp.where(im > 0, im, 0.0)
        sd = cp.fft.fftn(y64.astype(cp.complex128), axes=tuple(range(im.ndim)))
        m["cmul_d"](gn, _THREADS, (s, sd, h, n64))
        del sd
        plan.inverse(s, s)
        m["ratio_d"](gn, _THREADS, (s, im, n64, nn))
        plan.forward(s, s)
        m["cmul"](gn, _THREADS, (s, h, n64, np.int32(1)))
        plan.inverse(s, s)
        j2 = cp.empty(shape, dtype=cp.float32)
        m["rl_update_d"](gn, _THREADS, (j2, y64, s, n64, nn))
        if dbg is not None:
            dbg(1, j2, im)
        if n_iter == 1:
            return j2
        # k = 2: J_4 = J_2 - Y (Y double), Y = max(J_2 + 0 * (J_2 - J_3), 0) with J_3 = I double
        j4 = cp.empty(shape, dtype=cp.float32)
        m["diff_fd"](gn, _THREADS, (j4, j2, y64, n64))
        y = cp.empty(shape, dtype=cp.float32)
        m["predict_fd"](gn, _THREADS, (y, s, j2, y64, n64))
        del y64
        plan.forward(s, s)
        m["cmul"](gn, _THREADS, (s, h, n64, np.int32(0)))
        plan.inverse(s, s)
        m["ratio<double>"](gn, _THREADS, (s, im, n64, nn))
        plan.forward(s, s)
        m["cmul"](gn, _THREADS, (s, h, n64, np.int32(1)))
        plan.inverse(s, s)
        j3 = j2
        j2 = cp.empty(shape, dtype=cp.float32)
        m["rl_update"](gn, _THREADS, (j2, y, s, n64, nn))
        if dbg is not None:
            dbg(2, j2, im)
        first = 3
    else:
        j2 = _single(im)
        j3 = cp.zeros(shape, dtype=cp.float32)
        j4 = cp.zeros(shape, dtype=cp.float32)
        y = cp.empty(shape, dtype=cp.float32)
    ratio = m[f"ratio<{_CT[im.dtype]}>"]
    for k in range(first, n_iter + 1):
        if k > 2:
            _lambda(n, j2, y, j4, sums, lam)
        elif k == 2:
            m["diff"](gn, _THREADS, (j4, j2, y, n64))
        if rlstep is not None:
            m["predict"](gn, _THREADS, (y, y, j2, j3, lam, n64, np.int32(0)))
            j3 = j2
            j2 = rlstep(y, im, h)
        else:
            m["predict"](gn, _THREADS, (y, s, j2, j3, lam, n64, np.int32(1)))
            plan.forward(s, s)
            m["cmul"](gn, _THREADS, (s, h, n64, np.int32(0)))
            plan.inverse(s, s)
            ratio(gn, _THREADS, (s, im, n64, nn))
            plan.forward(s, s)
            m["cmul"](gn, _THREADS, (s, h, n64, np.int32(1)))
            plan.inverse(s, s)
            j3, j2 = j2, j3
            m["rl_update"](gn, _THREADS, (j2, y, s, n64, nn))
        if dbg is not None:
            dbg(k, j2, im)
    return j2


# ---------------------------------------------------------------------------------------------
# post-processing on the GPU
# ---------------------------------------------------------------------------------------------


def _post(
    g: _Geom,
    j: Any,
    matlab_layout: bool,
    damp_in: Any,
    mask: Any,
    *,
    damp_factor: Any,
    scale_factor: Any,
    decon_offset: Any,
    edge_erosion: Any,
    save_16bit: bool,
    bbox: Any,
) -> np.ndarray:
    """dampFactor, EdgeErosion, scaleFactor, deconOffset, save16bit and bbox as MATLAB's GPU code
    applies them, on the GPU, and the result copied to the host"""
    cp = _cp()
    damp, scale, offset, erosion = (
        _scalar(v) for v in (damp_factor, scale_factor, decon_offset, edge_erosion)
    )
    flags = 0
    if damp > 1:
        flags |= 1
    if erosion > 0:
        from petakit5d.kernels.decon_gpu import decon_mask_edge_erosion_gpu

        mask = decon_mask_edge_erosion_gpu(mask, edge_erosion)
    if scale != 1 or offset != 0 or erosion > 0:
        flags |= 2
        if offset != 0 or erosion > 0:
            flags |= 4 | 8
    z0, y0, x0, oz, oy, ox = _box(bbox, g)
    out = cp.empty((oz, oy, ox), dtype=cp.uint16 if save_16bit else cp.float32)
    if out.size:
        nz, ny, nx = g.zyx
        im = damp_in if damp_in is not None else j
        msk = mask if mask is not None else cp.empty(1, dtype=cp.bool_)
        name = f"post_{'t' if matlab_layout else 'l'}<{_CT[im.dtype]}, {_CT[out.dtype]}>"
        args = (
            out,
            j,
            im,
            msk,
            np.int32(ny),
            np.int32(nx),
            np.int32(z0),
            np.int32(y0),
            np.int32(x0),
            np.int32(oz),
            np.int32(oy),
            np.int32(ox),
            np.float64(damp),
            np.float64(scale),
            np.float64(offset),
            np.int32(flags),
        )
        if matlab_layout:
            _mod()[name](_tiles(oy, ox, oz), (_TILE, _ROWS), args)
        else:
            _mod()[name](_grid(out.size), _THREADS, args)
    h = _to_host(out)
    return h.reshape(h.shape[1:]) if len(g.shape) == 2 else h


# ---------------------------------------------------------------------------------------------
# the decon functions
# ---------------------------------------------------------------------------------------------


def _decon(method: str, i, psfs, numit, kw: dict) -> Any:
    from petakit5d.decon import rl

    use_gpu = kw["use_gpu"]
    cpu = rl.lucy if method == "simplified" else rl.lucy_omw
    reg = {k: v for k, v in kw.items() if k not in ("layout", "keep_plans")}
    if use_gpu is False or rl._gpu_device(use_gpu) is None:
        reg["use_gpu"] = False
        return cpu(i, *psfs, numit, **reg)
    layout = kw.get("layout", "matlab")
    if layout not in ("matlab", "native"):
        raise ValueError(f"layout must be 'matlab' or 'native', not {layout!r}")
    shape_in = np.shape(i)
    data = rl._data(i)
    scalar_bg, bg = rl._background(kw["background"])
    ps = tuple(rl._psf(p) for p in psfs)
    debug = bool(kw["debug"])
    use, _ = use_gpu_kernels(kw["use_mex"], method, data, debug, None if scalar_bg else bg)
    out, err_mat, n_iter, overlap = _run(
        _gpu_body, method, data, ps, numit, scalar_bg, bg, use, layout, kw, debug
    )
    if len(shape_in) == 3 and out.ndim == 2:
        out = out.reshape((1, *out.shape))
    if method == "simplified":
        return rl.RLResult(out, err_mat, n_iter, overlap)
    return rl.OMWResult(out, err_mat, overlap)


def _gpu_body(method, data, ps, numit, scalar_bg, bg, use, layout, kw, debug):
    cp = _cp()
    from petakit5d.decon import rl
    from petakit5d.decon.rl import DeconError
    from petakit5d.gpu import fft

    g = _geom(data.shape)
    keep = bool(kw.get("keep_plans", True))
    native = use and layout == "native"
    offset, erosion = _scalar(kw["decon_offset"]), _scalar(kw["edge_erosion"])
    need_mask = offset != 0 or erosion > 0
    # the memory check before anything is allocated (with the cached plan's own work area). The
    # order of the allocations is the one with the lowest device peak (measured, R04 report): the
    # OTFs first, then the data; for the GPU mex's OMW the data first, so that it doesn't split the
    # OTFs' freed scratch, a block CuPy's pool then can't give to the iterations' larger arrays
    plan_bytes = None
    if keep:
        if use:
            from petakit5d.kernels.decon_gpu import geometry

            plan_bytes = fft.get_plan(geometry(g.shape, layout).fft_shape).work_bytes
        else:
            plan_bytes = _complex_plan(g, True).work_bytes
    need = memory_need(
        method, g.shape, data.dtype, use_mex=use, mask=need_mask, plan_bytes=plan_bytes
    )
    _check_memory(need, _cached_otf_bytes(method) + (fft.cached_plan_bytes() if keep else 0))
    data_first = use and method == "omw"
    if data_first:
        raw = _to_gpu(data)
        mask = (raw != 0) if need_mask else None
    otfs = _otfs(method, ps, g, use, layout, keep)
    supported = data.dtype in (np.uint16, np.float32, np.float64)
    im = None
    if not data_first:
        raw = _to_gpu(data)
        mask = (raw != 0) if need_mask else None
        if supported:
            im = _prep(g, raw, bg, scalar_bg, not native)
            raw = None
    n_iter = rl._iterations(numit)
    if not supported:
        # MATLAB's GPU loop keeps other classes as they are; its first iteration mixes them with
        # single spectra
        if n_iter > 0:
            raise DeconError(
                _MIXED,
                "Integers can only be combined with integers of the same class, or scalar doubles.",
            )
        raise ValueError(f"{data.dtype} data without iterations are not supported on the GPU")
    damp = _scalar(kw["damp_factor"]) > 1
    overlap_fcn = kw["overlap_fcn"]
    if use and method == "omw":
        from petakit5d.kernels.decon_gpu import LucyRFFTGPU

        k = LucyRFFTGPU(g.shape, layout=layout, keep_plans=keep)
        j = k.omw(raw, otfs[0], otfs[1], numit, bg)
        damp_in = k.omw(raw, otfs[0], otfs[1], 0, bg) if damp else None
        del raw
        err_mat = np.zeros((rl._iterations(numit), 4))
        overlap = overlap_fcn() if overlap_fcn is not None else None
        res = _post(g, j, False, damp_in, mask, **_post_kw(kw))
        return res, err_mat, n_iter, overlap
    dbg = _DebugGPU(n_iter, kw["save_step"], kw["debug_folder"], g.n, not native) if debug else None
    if use:
        from petakit5d.kernels.decon_gpu import LucyRFFTGPU

        k = LucyRFFTGPU(g.shape, layout=layout, keep_plans=keep)
        if native:

            def rlstep(y, i, h):
                return k.rlstep(y, i, h, internal=True)
        else:

            def rlstep(y, i, h):
                return k.rlstep(y.T, i.T, h, internal=True).T

        j = _rl_loop(g, im, otfs[0], n_iter, None, dbg, rlstep)
    else:
        plan = _complex_plan(g, keep)
        try:
            if method == "omw":
                j = _omw_loop(g, im, otfs[0], otfs[1], n_iter, plan, dbg)
            else:
                j = _rl_loop(g, im, otfs[0], n_iter, plan, dbg)
        finally:
            if not keep:
                plan.close()
    if debug:
        err_mat = dbg.dbg.rows_kept(n_iter)
    else:
        err_mat = np.zeros((n_iter, 4))
    overlap = overlap_fcn() if overlap_fcn is not None else None
    if j.dtype == cp.float64:
        # double data and no iterations: J = I stays double (MATLAB's post-processing of a double
        # array, the same on the CPU)
        jh = _host_zyx(j, not native)
        mh = mask.get() if mask is not None else None
        res = rl._post(
            jh, lambda: jh, mh, **_post_kw(kw), shape_in=data.shape, threads=rl._threads(None)
        )
        return res, err_mat, n_iter, overlap
    res = _post(g, j, not native, im if damp else None, mask, **_post_kw(kw))
    return res, err_mat, n_iter, overlap


def _post_kw(kw: dict) -> dict:
    return {
        k: kw[k]
        for k in (
            "damp_factor",
            "scale_factor",
            "decon_offset",
            "edge_erosion",
            "save_16bit",
            "bbox",
        )
    }


def decon_lucy(
    i: Any,
    psf: Any,
    numit: Any,
    *,
    background: Any = None,
    use_gpu: bool | None = None,
    save_16bit: bool = True,
    damp_factor: Any = 1,
    scale_factor: Any = 1,
    decon_offset: Any = 0,
    edge_erosion: Any = 0,
    bbox: Sequence[Any] | None = None,
    debug: bool = False,
    debug_folder: str = "./debug/",
    save_step: Any = 1,
    use_mex: bool | None = None,
    overlap_fcn: Callable[[], Any] | None = None,
    layout: str = "matlab",
    keep_plans: bool = True,
):
    """decon_lucy_function on the GPU (``petakit5d.decon.rl.decon_lucy``'s keywords; see the
    module): an ``RLResult``. ``layout`` (the half spectra's, with ``use_mex``) and ``keep_plans``
    (keep the FFT plan and its work area between calls) are Python-only."""
    kw = dict(locals())
    del kw["i"], kw["psf"], kw["numit"]
    return _decon("simplified", i, (psf,), numit, kw)


def decon_lucy_omw(
    i: Any,
    psf_f: Any,
    psf_b: Any,
    numit: Any,
    *,
    background: Any = None,
    use_gpu: bool | None = None,
    save_16bit: bool = True,
    damp_factor: Any = 1,
    scale_factor: Any = 1,
    decon_offset: Any = 0,
    edge_erosion: Any = 0,
    bbox: Sequence[Any] | None = None,
    debug: bool = False,
    debug_folder: str = "./debug/",
    save_step: Any = 1,
    use_mex: bool | None = None,
    overlap_fcn: Callable[[], Any] | None = None,
    layout: str = "matlab",
    keep_plans: bool = True,
):
    """decon_lucy_omw_function on the GPU (``petakit5d.decon.rl.decon_lucy_omw``'s keywords; see
    the module): an ``OMWResult``."""
    kw = dict(locals())
    del kw["i"], kw["psf_f"], kw["psf_b"], kw["numit"]
    return _decon("omw", i, (psf_f, psf_b), numit, kw)
