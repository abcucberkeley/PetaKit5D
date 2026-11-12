// MATLAB R2024b's fft, ifft, fftn and ifftn, and the half-spectrum transforms of real data that FFT-based
// kernels convolve with, on strided arrays (the port's (z, y, x) C arrays: MATLAB's dimension 1 is y, of
// stride nx; 2 is x, of stride 1; 3 is z; ...). On the exact backend (FFTW 3.3.8 with MATLAB's codelets)
// the results are MATLAB's bit for bit (Linux x86-64, the C library's cos and sin as where the references
// were made); on the permissive one (pocketfft) they differ in the last bits:
//   - 1D plans per dimension (rank 1, FFTW_ESTIMATE, contiguous lines in scratch, out of place), one for
//     four lines per execute and one for a single line;
//   - per stage (the lines along one dimension) MATLAB runs lines in groups of four and some alone: when it
//     runs its FFTs in parallel and n * N * (bits of n) > 98304 (N lines of length n) it halves the line
//     range until pieces have at most 32 lines; each piece runs groups of four from its start and its last
//     (length % 4) lines alone; along a dimension d > 1 groups also restart every prod(size(1 : d - 1))
//     lines. Line j is numbered in MATLAB's order (its linear index among the stage's lines). FFTW computes
//     each lane of a group as it would alone in another lane, so only which lines run alone matters: the
//     driver visits lines in the memory order of the C array and packs the grouped ones in fours;
//   - real input: real-to-complex transforms along the first dimension that is transformed, conjugate fill of
//     each line (imaginary part negated), complex transforms along the others; fftn then makes the result
//     exactly Hermitian (an element whose mirror has the lower MATLAB linear index becomes the mirror's
//     conjugate, its imaginary part 0 - im; self-conjugate elements get imaginary part 0);
//   - inverse: complex transforms, or for Hermitian input (exact test, or 'symmetric') complex transforms
//     along the transformed dimensions but MATLAB's last one (the last of its size, 2 dimensions at least:
//     for a column vector the second, of length 1) and complex-to-real along that; ifftn's 'symmetric'
//     first applies fftn's fill to the padded data; every inverse stage divides by its length;
//   - padding and truncation to the requested size happen before the transforms;
//   - a complex result whose imaginary parts are all zero is returned as real.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace pk5d::fft {

using idx_t = std::int64_t;

// most lines per piece of MATLAB's split of a stage among its threads when it runs FFTs in parallel
constexpr idx_t kGrainParallel = 32;

// an array: data, MATLAB's dimensions (y, x, z, ...; at least 2), element strides per dimension (a complex
// element counts once), class
struct Array {
    void* data = nullptr;
    std::vector<idx_t> dims;
    std::vector<idx_t> strides;
    bool is_double = false;
    bool is_complex = false;
    idx_t numel() const {
        idx_t n = 1;
        for (idx_t d : dims) n *= d;
        return n;
    }
};

// the strides of the port's C layout for MATLAB's dimensions D: x (dimension 2) fastest, then y, then 3, 4, ...
std::vector<idx_t> native_strides(const std::vector<idx_t>& D);

// Result buffers are allocated by the driver through this callback (the binding hands them to NumPy): storage
// for numel elements of the given class, laid out with native_strides(dims).
using Allocator = void* (*)(void* ctx, const std::vector<idx_t>& dims, bool is_double, bool is_complex);

struct Options {
    int threads = 1;
    idx_t grain = kGrainParallel;  // 0: MATLAB's single-thread mode (no split)
    int backend = 0;               // kPermissive or kExact (backend.h)
};

// Y = fft(X, n, dim): n < 0 for size(X, dim); dim 0-based (MATLAB's dimension - 1; beyond X's: length 1)
Array fft(const Array& x, idx_t n, int dim, const Options& o, Allocator alloc, void* ctx);
// X = ifft(Y, n, dim[, 'symmetric'])
Array ifft(const Array& x, idx_t n, int dim, bool symmetric, const Options& o, Allocator alloc, void* ctx);
// Y = fftn(X, siz): siz (MATLAB's order) empty for size(X); otherwise its length is at least ndims(X)
Array fftn(const Array& x, const std::vector<idx_t>& siz, const Options& o, Allocator alloc, void* ctx);
// X = ifftn(Y, siz[, 'symmetric'])
Array ifftn(const Array& x, const std::vector<idx_t>& siz, bool symmetric, const Options& o, Allocator alloc, void* ctx);

// MATLAB's half spectrum of real data of dimensions D, through MATLAB's driver on the backend of o: the
// elements of fftn(x) with index <= n / 2 along MATLAB's last transformed dimension (the part ifftn's
// complex-to-real stage reads); MATLAB's bits on the exact backend. half_dims gives its dimensions.
std::vector<idx_t> half_dims(const std::vector<idx_t>& D);
Array half_forward(const Array& x, const Options& o, Allocator alloc, void* ctx);
// real(ifftn(S)) of the conjugate-symmetric spectrum whose half is h (dims half_dims(D)), D the real size;
// h is overwritten
Array half_inverse(Array& h, const std::vector<idx_t>& D, const Options& o, Allocator alloc, void* ctx);
// h .* g (conjugate: conj(g) .* h), MATLAB's complex product (no fused multiply-add), in place in h
void half_multiply(Array& h, const Array& g, bool conjugate, int threads);

// whether MATLAB transforms line j (0-based) of a stage of N lines of length n alone (seg: restart period of
// the groups, 0 for none)
bool line_alone(idx_t j, idx_t N, idx_t n, idx_t seg, idx_t grain);

// the plan the driver uses on a backend: kind "r2c", "c2c-1" (forward), "c2c+1" (backward), "c2r"; howmany 1
// or 4; its description (FFTW's sprint_plan)
std::string plan_string(int backend, bool is_double, const std::string& kind, idx_t n, int howmany);

// drops the cached plans of a backend (before it changes)
void clear_plans(int backend);

}  // namespace pk5d::fft
