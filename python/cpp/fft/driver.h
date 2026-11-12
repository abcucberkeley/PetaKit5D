// MATLAB R2024b's fft, ifft, fftn and ifftn on FFTW, bit for bit (Linux x86-64):
//   - 1D FFTW plans per dimension (rank 1, FFTW_ESTIMATE, contiguous lines in scratch, out of place), one
//     for four lines per execute and one for a single line;
//   - per stage (the lines along one dimension) MATLAB runs lines in groups of four and some alone: when it
//     runs its FFTs in parallel and n * N * (bits of n) > 98304 (N lines of length n) it halves the line
//     range until pieces have at most 32 lines; each piece runs groups of four from its start and its last
//     (length % 4) lines alone; along a dimension d > 1 groups also restart every prod(size(1 : d - 1))
//     lines. The lanes of a group are independent, so only which lines run alone matters;
//   - real input: real-to-complex transforms along the first dimension that is transformed, conjugate fill of
//     each line (imaginary part negated), complex transforms along the others; fftn then makes the result
//     exactly Hermitian (an element whose mirror has the lower linear index becomes the mirror's conjugate,
//     its imaginary part 0 - im; self-conjugate elements get imaginary part 0);
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

// an array in Fortran order: data, dimensions, class
struct Array {
    void* data = nullptr;
    std::vector<idx_t> dims;
    bool is_double = false;
    bool is_complex = false;
    idx_t numel() const {
        idx_t n = 1;
        for (idx_t d : dims) n *= d;
        return n;
    }
};

// Result buffers are allocated by the driver through this callback (the binding hands them to NumPy):
// returns storage for numel elements of the given class.
using Allocator = void* (*)(void* ctx, const std::vector<idx_t>& dims, bool is_double, bool is_complex);

struct Options {
    int threads = 1;
    idx_t grain = kGrainParallel;  // 0: MATLAB's single-thread mode (no split)
};

// Y = fft(X, n, dim): n < 0 for size(X, dim); dim 0-based (dims beyond X's are length 1)
Array fft(const Array& x, idx_t n, int dim, const Options& o, Allocator alloc, void* ctx);
// X = ifft(Y, n, dim[, 'symmetric'])
Array ifft(const Array& x, idx_t n, int dim, bool symmetric, const Options& o, Allocator alloc, void* ctx);
// Y = fftn(X, siz): siz empty for size(X); otherwise its length is at least ndims(X)
Array fftn(const Array& x, const std::vector<idx_t>& siz, const Options& o, Allocator alloc, void* ctx);
// X = ifftn(Y, siz[, 'symmetric'])
Array ifftn(const Array& x, const std::vector<idx_t>& siz, bool symmetric, const Options& o, Allocator alloc, void* ctx);

// whether MATLAB transforms line j (0-based) of a stage of N lines of length n alone (seg: restart period of
// the groups, 0 for none)
bool line_alone(idx_t j, idx_t N, idx_t n, idx_t seg, idx_t grain);

// sprint_plan of the plan the driver uses: kind "r2c", "c2c-1" (forward), "c2c+1" (backward), "c2r";
// howmany 1 or 4
std::string plan_string(bool is_double, const std::string& kind, idx_t n, int howmany);

// drops the cached plans (before the backend changes)
void clear_plans();

}  // namespace pk5d::fft
