// The deconvolution kernels of petakit5d._ext.decon (moved from the mex files of
// microscopeDataProcessing/decon/mex): the half-spectrum Richardson-Lucy steps of
// decon_lucy_rfft_mex in MATLAB's layout (rfft_matlab.cpp) and on the port's (z, y, x) layout
// (rfft_native.cpp), the acceleration step and fixed-order sums of decon_lucy_accel_mex
// (accel.cpp), the mask erosion of decon_mask_edge_erosion_mex (erosion.cpp), and the window
// bounds, components and sort of omw_backprojector_mex (omw.cpp).
//
// MATLAB's layout of a (z, y, x) volume is its (y, x, z) Fortran array: the C array (z, x, y),
// y fastest. Element (z, y, x) has MATLAB's linear index y + ny * (x + nx * z).
#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <new>
#include <vector>

namespace pk5d::decon {

using idx_t = std::ptrdiff_t;
using cfloat = std::complex<float>;

// --- buffers ---------------------------------------------------------------------------------------

// malloc'd storage, freed unless released (the bindings hand released storage to NumPy)
template <typename T>
struct Owned {
    T* p = nullptr;
    Owned() = default;
    explicit Owned(size_t n) { alloc(n); }
    Owned(const Owned&) = delete;
    Owned& operator=(const Owned&) = delete;
    ~Owned() { std::free(p); }
    T* alloc(size_t n) {
        std::free(p);
        p = static_cast<T*>(std::malloc(sizeof(T) * (n > 0 ? n : 1)));
        if (p == nullptr) throw std::bad_alloc();
        return p;
    }
    T* release() {
        T* q = p;
        p = nullptr;
        return q;
    }
};

// --- layouts ---------------------------------------------------------------------------------------

// dst[z][c][r] = src[z][r][c] for nz planes of rows x cols elements of `bytes` bytes (1, 2, 4 or 8):
// the port's (z, y, x) to MATLAB's layout with rows = ny, cols = nx, and back with rows = nx,
// cols = ny. Blocked, planes and blocks split among the threads.
void transpose_planes(const void* src, void* dst, idx_t nz, idx_t rows, idx_t cols, int bytes, int nt);

// --- decon_lucy_rfft_mex in MATLAB's layout, on a backend's FFT table (rfft_matlab.cpp) ------------

// seconds in planning, forward transforms, inverse transforms, other
struct Times {
    double t[4] = {0, 0, 0, 0};
};

// the call's sizes: d = {ny, nx, nz} (nz = 1 for 2D), threads, MATLAB's split grain, FFT backend
struct Call {
    int d[3] = {1, 1, 1};
    int nt = 1;
    idx_t grain = 0;
    int backend = 0;
};

// The FFT plans, the lines MATLAB transforms alone and the scratch of one call: made (and destroyed) with
// the interpreter lock held, since the planner of the FFT library is shared with petakit5d._ext.fft, whose
// calls hold it; the steps below run without it.
struct RfftPlan;
struct RfftPlanDeleter {
    void operator()(RfftPlan* p) const;
};
using RfftPlanPtr = std::unique_ptr<RfftPlan, RfftPlanDeleter>;
RfftPlanPtr matlab_plan(const Call& c, Times& tm);

// the half OTF [ny, nx, nz / 2 + 1] (2D: [ny, nx / 2 + 1]) in MATLAB's layout, of the (z, y, x) PSF
// of size ps = {pz, py, px}
void matlab_otf(RfftPlan& c, const float* psf, const int* ps, cfloat* h, Times& tm);
// OMW iterations; raw and j0 (or null) and the result in MATLAB's layout, or all in the port's
// (z, y, x) layout (native_io); the result: n floats from malloc
template <typename T>
float* matlab_omw(RfftPlan& c, const T* raw, double bg, const cfloat* hf, const cfloat* hb, int n_iter,
                  const float* j0, bool native_io, Times& tm);
// one RL step (the same layouts)
float* matlab_rlstep(RfftPlan& c, const float* y, const float* im, const cfloat* h, bool native_io, Times& tm);
// rows 1 .. n / 2 + 1 of fft(X, [], 1) of the n x m single matrix X (columns contiguous) into y
// (columns of n / 2 + 1 complex values); the plan's size is {n, m, 1}
void matlab_fft1(RfftPlan& c, const float* x, cfloat* y);

// --- the same steps on the port's layout with the native half spectra (rfft_native.cpp) ------------

// shape: (nz, ny, nx) or (ny, nx); the half spectra (..., nx / 2 + 1)
void native_otf(const std::vector<std::int64_t>& shape, const float* psf, const std::vector<std::int64_t>& ps,
                cfloat* h, int nt, Times& tm);
template <typename T>
float* native_omw(const std::vector<std::int64_t>& shape, const T* raw, double bg, const cfloat* hf,
                  const cfloat* hb, int n_iter, const float* j0, int nt, Times& tm);
float* native_rlstep(const std::vector<std::int64_t>& shape, const float* y, const float* im, const cfloat* h,
                     int nt, Times& tm);

// --- decon_lucy_accel_mex (accel.cpp) ----------------------------------------------------------------

// the sums of (J_2 - Y) .* J_4 and J_4 .* J_4 in double (s), and J_4 = J_2 - Y into out; zyx: the
// arrays are (planes, ny, nx) C arrays summed in MATLAB's linear order, else flat in memory order
template <typename T>
void accel_lambda(const T* j2, const T* y, const T* j4, T* out, idx_t n, bool zyx, idx_t ny, idx_t nx, int nt,
                  double* s);
// max(J_2 + lambda * (J_2 - J_3), 0)
template <typename T>
void accel_predict(const T* j2, const T* j3, double lambda, T* out, idx_t n, int nt);
// sum(x(:), 'double') in blocks of 4096 (any size, the mex's order)
template <typename T>
double accel_sum(const T* x, idx_t n, bool zyx, idx_t ny, idx_t nx, int nt);
// lambda from the sums, clamped as MATLAB's max(min(., 1), 0)
double clamped_lambda(const double* s);

// --- decon_mask_edge_erosion_mex on (z, y, x) masks (erosion.cpp) -------------------------------------

struct ErosionPlan {
    std::vector<int> start, dr, dc;  // neighborhood k: offsets [start[k], start[k + 1]) (row = y, column = x)
    idx_t pu[2] = {0, 0}, pl[2] = {0, 0};  // padding before / after: (rows, columns)
    idx_t rz = 0;
    bool clear_border = false;
};
void mask_edge_erosion(const uint8_t* a, uint8_t* b, idx_t nz, idx_t ny, idx_t nx, const ErosionPlan& p, int nt);

// --- omw_backprojector_mex (omw.cpp) -----------------------------------------------------------------

// 'window' on a grid of psz = {n1, n2, n3} (MATLAB's [ny, nx, nz]) with edge directions E (nE x 3, columns
// y, x, z, row-major) and radii R: mask (C layout (n3, n1, n2)), the MATLAB linear indices (0-based) of the
// voxels to compute, in increasing order
void omw_window(const idx_t* psz, const double* e, const double* r, idx_t ne, double l, double u, double* mask,
                std::vector<int64_t>& inds, int nt);
// 'components' of a mask in MATLAB's layout (n1 fastest): keep those with at least n voxels, or the n
// largest (ties to the first in MATLAB's order)
void omw_components(const uint8_t* a, idx_t n1, idx_t n2, idx_t n3, bool largest, idx_t n, uint8_t* b, int nt);
// 'sortdesc': false for NaN or negative values (-0 included)
bool omw_sortdesc(const float* a, idx_t n, double* out, int nt);
// omw_backprojector_generation's window (R02): 'window''s mask for the edge voxels et (MATLAB linear
// indices, increasing) with the window of the voxels it leaves computed as MATLAB computes it (cart2sph,
// sph2cart, the 2 nearest edge points, the radii weighted with the first voxels' distances); false where
// MATLAB computes every voxel (a zero weight or radius)
bool omw_window_values(const idx_t* psz, const int64_t* et, idx_t ne, double l, double u, double* mask, int nt);
// MATLAB's abs of n complex singles (interleaved): msun's hypotf, MATLAB's NaN for a NaN part
void omw_abs(const float* c, idx_t n, float* out, int nt);

}  // namespace pk5d::decon
