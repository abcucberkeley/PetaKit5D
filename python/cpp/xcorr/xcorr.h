// The registration's normalized cross-correlation (normxcorr3_fast and normxcorr3_window_mex of PetaKit5D) in
// a window of its 'full' output.
//
// The kernels work in MATLAB's layout (y fastest, then x, then z: "MATLAB layout" below), in which MATLAB's FFT
// stages, its groups of lines and its cumulative sums are defined; ncc.cpp moves the port's (z, y, x) C arrays
// into it and the results back. With the exact FFT backend the results are MATLAB's bit for bit.
//
//   window.cpp  normxcorr3_window_mex's core, moved unchanged from the mex: the correlation
//                 corr = real(ifftn(fftn(A0, szOut) .* fftn(T0(end : -1 : 1, end : -1 : 1, end : -1 : 1), szOut)))
//               with MATLAB's FFTs (on the fft group's backends), and the window sums of normxcorr3_fast in
//               double (box_sums: cumsum along y, then x, then z, P(hi) - P(lo - 1))
//   ncc.cpp     the rest of normxcorr3_fast: the means, the template's sum and standard deviation and the
//               image's variance as MATLAB's sum, mean, std and var compute them, the NCC in double, the
//               first maximum, and the layout changes
#pragma once

#include <cstddef>
#include <cstdint>

namespace pk5d::xcorr {

using idx_t = std::ptrdiff_t;

// the window of the full output in MATLAB's order (y, x, z): w0 .. w1 (0-based, inclusive), wn elements
struct Win {
    int w0[3], w1[3], wn[3];
};

// --- window.cpp (MATLAB layout) ---------------------------------------------------------------------------

// normxcorr3_window_mex('ncc', T0, A0, bbox, nt, grain, fast): the correlation of T0 (size ts) and A0 (size a)
// in the window W (corr: window size, single) on FFT backend `backend` (pk5d::fft::kPermissive, kExact);
// tm: seconds of [planning, forward T0, forward A0, inverse] (or null). The caller checks: 1 <= ts <= a, the
// full size ts + a - 1 >= 2 along y and x, W inside it. Throws std::bad_alloc, std::runtime_error.
void window_corr(const float* T0, const int* ts, const float* A0, const int* a, const Win& W, int nt, idx_t grain,
    bool fast, int backend, float* corr, double* tm);

// The window sums of normxcorr3_fast in double, for the window positions of W: of A0 and A0 .* A0 over the
// image's elements that a template at each position overlaps (zero outside A0), and of T0 over the
// template's elements that overlap the image. Null outputs are skipped.
template <typename T>
void window_sums(const T* T0, const int* ts, const T* A0, const int* a, const Win& W, double* sA, double* sA2,
    double* sT, int nt);

// --- ncc.cpp (MATLAB layout unless named otherwise) --------------------------------------------------------

// MATLAB's sums over x(:) (the vector rule, as MATLAB R2024b adds on one thread): in index order up to 43
// elements; up to 88,999 in four lanes (element t in lane t mod 4, each lane in index order, added as
// ((l0 + l1) + l2) + l3); from 89,000 in blocks of 4,096 elements, each in four lanes, the block sums in
// index order. The element function f gives each addend in the class of the sum.
struct Stats {
    double cT, cA;    // cast(mean(T(:), 'double'), 'like', T), cast(mean(A(:), 'double'), 'like', A)
    double sumT0;     // sum(T0(:), 'double')
    double denomT;    // std(double(T0(:)))
    double varA1;     // var(A(:), 1), in A's class
};

// The statistics of T and A (MATLAB layout, n elements each) and T0 = T - cT, A0 = A - cA (in the class) into
// T0, A0 (may be T, A).
template <typename T>
Stats prepare(const T* Tm, idx_t nT, const T* Am, idx_t nA, T* T0, T* A0, int nt);

// the NCC of normxcorr3_fast's last lines in double, in the window (MATLAB layout), into C (window size, A's
// class): 0 where the window's variance is at most eps(class) * var(A(:), 1) or the template is flat
template <typename T, typename C>
void ncc_combine(const C* corr, const double* sA, const double* sA2, const double* sT, const int* ts,
    const int* a, const Win& W, const Stats& st, T* out, int nt);

// [m, i] = max(C(:)): the first largest value in linear order, NaN ignored (all NaN: index 0)
template <typename T>
idx_t first_max(const T* C, idx_t n, int nt);

// (z, y, x) C array (nz, ny, nx) <-> MATLAB layout, with nt threads
template <typename T>
void c_to_matlab(const T* src, idx_t nz, idx_t ny, idx_t nx, T* dst, int nt);
template <typename T>
void matlab_to_c(const T* src, idx_t nz, idx_t ny, idx_t nx, T* dst, int nt);

}  // namespace pk5d::xcorr
