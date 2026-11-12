// The permissive half spectrum: plain real-to-complex transforms of (..., z, y, x) C arrays on pocketfft's own
// multi-dimensional code, the spectrum halved along x (the last axis: (..., z, y, x / 2 + 1)), and the
// complex-to-real inverse divided by the number of elements. The results do not depend on the thread count.
#pragma once

#include <cstdint>
#include <vector>

namespace pk5d::fft {

// the half spectrum's shape of a real array of shape `shape` (C order, the last axis halved)
std::vector<std::int64_t> native_half_shape(const std::vector<std::int64_t>& shape);

// x (real, C-contiguous, shape) -> h (complex, C-contiguous, native_half_shape(shape))
void native_forward(const float* x, const std::vector<std::int64_t>& shape, float* h, int threads);
void native_forward(const double* x, const std::vector<std::int64_t>& shape, double* h, int threads);

// h (complex, native_half_shape(shape); overwritten) -> y (real, shape): real(ifftn) of the conjugate-
// symmetric spectrum whose half is h
void native_inverse(float* h, const std::vector<std::int64_t>& shape, float* y, int threads);
void native_inverse(double* h, const std::vector<std::int64_t>& shape, double* y, int threads);

}  // namespace pk5d::fft
