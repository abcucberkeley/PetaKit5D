// The permissive backend: the part of FFTW's guru64 API that MATLAB's fft driver calls (rank-1 plans of
// several lines, out of place, new-array execute), on pocketfft (BSD-3, third_party/pocketfft).
//
// Lines are transformed VLEN at a time with pocketfft's vector types (GCC and Clang vector extensions; MSVC
// has none, so there every line runs alone), the rest one at a time. pocketfft computes every lane with the
// operations of a single line, so the results do not depend on the vector width, on which lines share an
// execute, or on the thread count: the SSE2, AVX and scalar builds give the same bits.
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>

#include "backend.h"

#define POCKETFFT_NAMESPACE pk5d_pocketfft
#define POCKETFFT_NO_MULTITHREADING
#define POCKETFFT_CACHE_SIZE 0
#include "pocketfft_hdronly.h"

namespace pk5d::fft {
namespace {

using pk5d_pocketfft::detail::arr;
using pk5d_pocketfft::detail::cmplx;
using pk5d_pocketfft::detail::pocketfft_c;
using pk5d_pocketfft::detail::pocketfft_r;
using pk5d_pocketfft::detail::VLEN;

enum Kind { kC2C = 0, kR2C = 1, kC2R = 2 };

template <typename T>
struct Plan {
    int kind = kC2C;
    bool fwd = true;
    size_t n = 0;
    std::ptrdiff_t is = 1, os = 1;  // strides within a line (FFTW's units: complex or real elements)
    size_t hm = 1;
    std::ptrdiff_t idist = 0, odist = 0;
    std::unique_ptr<pocketfft_c<T>> c;
    std::unique_ptr<pocketfft_r<T>> r;
};

template <typename T>
Plan<T>* make_plan(int kind, int rank, const Iodim* dims, int hrank, const Iodim* hdims, int sign) {
    if (rank != 1 || hrank > 1 || dims[0].n < 1) return nullptr;
    auto* p = new Plan<T>;
    p->kind = kind;
    p->fwd = kind == kR2C || (kind == kC2C && sign < 0);
    p->n = static_cast<size_t>(dims[0].n);
    p->is = dims[0].is;
    p->os = dims[0].os;
    if (hrank == 1) {
        p->hm = static_cast<size_t>(hdims[0].n);
        p->idist = hdims[0].is;
        p->odist = hdims[0].os;
    }
    try {
        if (kind == kC2C) {
            p->c.reset(new pocketfft_c<T>(p->n));
        } else {
            p->r.reset(new pocketfft_r<T>(p->n));
        }
    } catch (...) {
        delete p;
        return nullptr;
    }
    return p;
}

template <typename V, size_t W, typename T>
inline T lane(const V& v, size_t l) {
    if constexpr (W == 1) {
        (void)l;
        return v;
    } else {
        return v[l];
    }
}

template <typename V, size_t W, typename T>
inline void set_lane(V& v, size_t l, T x) {
    if constexpr (W == 1) {
        (void)l;
        v = x;
    } else {
        v[l] = x;
    }
}

// W lines from line j0 (W = 1: scalar; else pocketfft's vector type, lane l = line j0 + l)
template <typename T, typename V, size_t W>
void run_c2c(const Plan<T>& p, const T* in, T* out, size_t j0) {
    thread_local arr<cmplx<V>> buf;
    if (buf.size() < p.n) buf.resize(p.n);
    cmplx<V>* b = buf.data();
    for (size_t l = 0; l < W; l++) {
        const T* s = in + 2 * (static_cast<std::ptrdiff_t>(j0 + l) * p.idist);
        for (size_t k = 0; k < p.n; k++) {
            set_lane<V, W, T>(b[k].r, l, s[2 * k * p.is]);
            set_lane<V, W, T>(b[k].i, l, s[2 * k * p.is + 1]);
        }
    }
    p.c->exec(b, T(1), p.fwd);
    for (size_t l = 0; l < W; l++) {
        T* t = out + 2 * (static_cast<std::ptrdiff_t>(j0 + l) * p.odist);
        for (size_t k = 0; k < p.n; k++) {
            t[2 * k * p.os] = lane<V, W, T>(b[k].r, l);
            t[2 * k * p.os + 1] = lane<V, W, T>(b[k].i, l);
        }
    }
}

// real to complex: pocketfft's half-complex order (r0, r1, i1, r2, i2, ...) into FFTW's n / 2 + 1 values
template <typename T, typename V, size_t W>
void run_r2c(const Plan<T>& p, const T* in, T* out, size_t j0) {
    thread_local arr<V> buf;
    if (buf.size() < p.n) buf.resize(p.n);
    V* b = buf.data();
    const size_t n = p.n;
    for (size_t l = 0; l < W; l++) {
        const T* s = in + static_cast<std::ptrdiff_t>(j0 + l) * p.idist;
        for (size_t k = 0; k < n; k++) set_lane<V, W, T>(b[k], l, s[k * p.is]);
    }
    p.r->exec(b, T(1), true);
    for (size_t l = 0; l < W; l++) {
        T* t = out + 2 * (static_cast<std::ptrdiff_t>(j0 + l) * p.odist);
        t[0] = lane<V, W, T>(b[0], l);
        t[1] = T(0);
        size_t i = 1, ii = 1;
        for (; i + 1 < n; i += 2, ++ii) {
            t[2 * ii * p.os] = lane<V, W, T>(b[i], l);
            t[2 * ii * p.os + 1] = lane<V, W, T>(b[i + 1], l);
        }
        if (i < n) {
            t[2 * ii * p.os] = lane<V, W, T>(b[i], l);
            t[2 * ii * p.os + 1] = T(0);
        }
    }
}

// complex to real: FFTW's n / 2 + 1 values (the imaginary parts of the first and, n even, the last ignored)
template <typename T, typename V, size_t W>
void run_c2r(const Plan<T>& p, const T* in, T* out, size_t j0) {
    thread_local arr<V> buf;
    if (buf.size() < p.n) buf.resize(p.n);
    V* b = buf.data();
    const size_t n = p.n;
    for (size_t l = 0; l < W; l++) {
        const T* s = in + 2 * (static_cast<std::ptrdiff_t>(j0 + l) * p.idist);
        set_lane<V, W, T>(b[0], l, s[0]);
        size_t i = 1, ii = 1;
        for (; i + 1 < n; i += 2, ++ii) {
            set_lane<V, W, T>(b[i], l, s[2 * ii * p.is]);
            set_lane<V, W, T>(b[i + 1], l, s[2 * ii * p.is + 1]);
        }
        if (i < n) set_lane<V, W, T>(b[i], l, s[2 * ii * p.is]);
    }
    p.r->exec(b, T(1), false);
    for (size_t l = 0; l < W; l++) {
        T* t = out + static_cast<std::ptrdiff_t>(j0 + l) * p.odist;
        for (size_t k = 0; k < n; k++) t[k * p.os] = lane<V, W, T>(b[k], l);
    }
}

template <typename T>
void execute(void* plan, void* in, void* out) {
    const Plan<T>* p = static_cast<const Plan<T>*>(plan);
    const T* src = static_cast<const T*>(in);
    T* dst = static_cast<T*>(out);
    size_t j = 0;
#ifndef POCKETFFT_NO_VECTORS
    constexpr size_t W = VLEN<T>::val;
    using V = pk5d_pocketfft::detail::vtype_t<T>;
    for (; j + W <= p->hm; j += W) {
        if (p->kind == kC2C) {
            run_c2c<T, V, W>(*p, src, dst, j);
        } else if (p->kind == kR2C) {
            run_r2c<T, V, W>(*p, src, dst, j);
        } else {
            run_c2r<T, V, W>(*p, src, dst, j);
        }
    }
#endif
    for (; j < p->hm; j++) {
        if (p->kind == kC2C) {
            run_c2c<T, T, 1>(*p, src, dst, j);
        } else if (p->kind == kR2C) {
            run_r2c<T, T, 1>(*p, src, dst, j);
        } else {
            run_c2r<T, T, 1>(*p, src, dst, j);
        }
    }
}

template <typename T>
char* sprint(void* plan) {
    const Plan<T>* p = static_cast<const Plan<T>*>(plan);
    char s[160];
    std::snprintf(s, sizeof(s), "(pocketfft-%s-%s %zu x%zu)",
                  p->kind == kC2C ? "c2c" : (p->kind == kR2C ? "r2c" : "c2r"), p->fwd ? "fwd" : "bwd", p->n, p->hm);
    char* r = static_cast<char*>(std::malloc(std::strlen(s) + 1));
    if (r != nullptr) std::strcpy(r, s);
    return r;
}

template <typename T>
void* plan_dft(int r, const Iodim* d, int hr, const Iodim* h, void*, void*, int sign, unsigned) {
    return make_plan<T>(kC2C, r, d, hr, h, sign);
}
template <typename T>
void* plan_r2c(int r, const Iodim* d, int hr, const Iodim* h, void*, void*, unsigned) {
    return make_plan<T>(kR2C, r, d, hr, h, -1);
}
template <typename T>
void* plan_c2r(int r, const Iodim* d, int hr, const Iodim* h, void*, void*, unsigned) {
    return make_plan<T>(kC2R, r, d, hr, h, 1);
}
template <typename T>
void destroy(void* p) {
    delete static_cast<Plan<T>*>(p);
}
void free_string(void* s) { std::free(s); }

template <typename T>
FftwApi table(const char* version) {
    FftwApi a;
    a.plan_dft = plan_dft<T>;
    a.plan_r2c = plan_r2c<T>;
    a.plan_c2r = plan_c2r<T>;
    a.exec_dft = execute<T>;
    a.exec_r2c = execute<T>;
    a.exec_c2r = execute<T>;
    a.destroy = destroy<T>;
    a.sprint = sprint<T>;
    a.free_string = free_string;
    a.version = version;
    return a;
}

}  // namespace

FftwApi pocketfft_api(int prec) {
    return prec == 0 ? table<float>("pocketfft-cpp-4d865cba") : table<double>("pocketfft-cpp-4d865cba");
}

}  // namespace pk5d::fft
