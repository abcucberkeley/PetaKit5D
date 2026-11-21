// The permissive half spectrum on pocketfft's multi-dimensional transforms (its own threads; a namespace of its
// own, so that this copy of the header, with threads, stays apart from the driver's backend copy without them).
#include "halfspec.h"

#include <complex>
#include <stdexcept>

#include "pk5d/x87.hpp"

#define POCKETFFT_NAMESPACE pk5d_pocketfft_mt
#define POCKETFFT_CACHE_SIZE 16
#include "pocketfft_hdronly.h"

#if defined(__MINGW32__)
#include <mutex>

// MinGW-w64: pocketfft's thread pool is a static object of each module that links this file. At the
// process's exit Windows ends every other thread first, then runs the modules' static destructors, and
// destroying the pool's condition variables (winpthreads) waits forever for the threads that waited on
// them: the process never exits. So the pool's threads are ended at the interpreter's exit instead
// (Py_AtExit, registered at the first multi-threaded transform), while they still run.
extern "C" int Py_AtExit(void (*func)(void));
#endif

namespace pk5d::fft {
namespace {

#if defined(__MINGW32__)
std::mutex& pool_mutex() {
    static std::mutex m;
    return m;
}
bool pool_running = false;
bool pool_ended = false;

void end_pool() {
    std::lock_guard<std::mutex> g(pool_mutex());
    if (pool_running) {
        pk5d_pocketfft_mt::detail::threading::get_pool().shutdown();
        pool_running = false;
        pool_ended = true;
    }
}

// before a transform on nt threads: the pool runs, and is ended at the interpreter's exit
void use_pool(size_t nt) {
    if (nt <= 1) return;
    std::lock_guard<std::mutex> g(pool_mutex());
    if (pool_running) return;
    static bool registered = false;
    if (!registered) registered = Py_AtExit(end_pool) == 0;
    if (pool_ended) {
        pk5d_pocketfft_mt::detail::threading::get_pool().restart();
    } else {
        pk5d_pocketfft_mt::detail::threading::get_pool();
    }
    pool_running = true;
}
#else
inline void use_pool(size_t) {}
#endif

using shape_t = pk5d_pocketfft_mt::shape_t;
using stride_t = pk5d_pocketfft_mt::stride_t;

shape_t to_shape(const std::vector<std::int64_t>& s) {
    shape_t r;
    for (std::int64_t v : s) {
        if (v < 1) throw std::invalid_argument("a half spectrum needs every length >= 1");
        r.push_back(static_cast<size_t>(v));
    }
    if (r.empty()) throw std::invalid_argument("a half spectrum needs at least one axis");
    return r;
}

// C strides in bytes of a shape
stride_t c_strides(const shape_t& s, size_t elem) {
    stride_t st(s.size());
    std::ptrdiff_t m = static_cast<std::ptrdiff_t>(elem);
    for (size_t k = s.size(); k-- > 0;) {
        st[k] = m;
        m *= static_cast<std::ptrdiff_t>(s[k]);
    }
    return st;
}

shape_t all_axes(size_t nd) {
    shape_t a(nd);
    for (size_t k = 0; k < nd; k++) a[k] = k;
    return a;
}

template <typename T>
void forward_t(const T* x, const std::vector<std::int64_t>& shape, T* h, int threads) {
    use_pool(static_cast<size_t>(threads > 0 ? threads : 1));
    pk5d::x87_extended x87;  // pocketfft's plans (pk5d/x87.hpp), made on this thread
    const shape_t s = to_shape(shape);
    shape_t sh = s;
    sh.back() = s.back() / 2 + 1;
    pk5d_pocketfft_mt::r2c(s, c_strides(s, sizeof(T)), c_strides(sh, sizeof(std::complex<T>)), all_axes(s.size()),
        pk5d_pocketfft_mt::FORWARD, x, reinterpret_cast<std::complex<T>*>(h), T(1),
        static_cast<size_t>(threads > 0 ? threads : 1));
}

// pocketfft's multi-dimensional c2r: complex-to-complex along the other axes, then complex-to-real along the
// last, the factor applied there; here the first part runs in place on h, without pocketfft's temporary copy
template <typename T>
void inverse_t(T* h, const std::vector<std::int64_t>& shape, T* y, int threads) {
    use_pool(static_cast<size_t>(threads > 0 ? threads : 1));
    pk5d::x87_extended x87;
    const shape_t s = to_shape(shape);
    shape_t sh = s;
    sh.back() = s.back() / 2 + 1;
    double n = 1;
    for (size_t v : s) n *= static_cast<double>(v);
    const size_t nt = static_cast<size_t>(threads > 0 ? threads : 1);
    const stride_t hst = c_strides(sh, sizeof(std::complex<T>));
    auto* hc = reinterpret_cast<std::complex<T>*>(h);
    if (s.size() > 1) {
        const shape_t first(all_axes(s.size() - 1));
        pk5d_pocketfft_mt::c2c(sh, hst, hst, first, pk5d_pocketfft_mt::BACKWARD, hc, hc, T(1), nt);
    }
    pk5d_pocketfft_mt::c2r(s, hst, c_strides(s, sizeof(T)), s.size() - 1, pk5d_pocketfft_mt::BACKWARD, hc, y,
        static_cast<T>(1.0 / n), nt);
}

}  // namespace

std::vector<std::int64_t> native_half_shape(const std::vector<std::int64_t>& shape) {
    std::vector<std::int64_t> r(shape);
    if (r.empty()) throw std::invalid_argument("a half spectrum needs at least one axis");
    r.back() = r.back() / 2 + 1;
    return r;
}

void native_forward(const float* x, const std::vector<std::int64_t>& shape, float* h, int threads) {
    forward_t(x, shape, h, threads);
}
void native_forward(const double* x, const std::vector<std::int64_t>& shape, double* h, int threads) {
    forward_t(x, shape, h, threads);
}
void native_inverse(float* h, const std::vector<std::int64_t>& shape, float* y, int threads) {
    inverse_t(h, shape, y, threads);
}
void native_inverse(double* h, const std::vector<std::int64_t>& shape, double* y, int threads) {
    inverse_t(h, shape, y, threads);
}

}  // namespace pk5d::fft
