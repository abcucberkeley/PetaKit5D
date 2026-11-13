// The Richardson-Lucy steps of decon_lucy_rfft_mex on the port's (z, y, x) arrays with the permissive
// backend's half spectra (cpp/fft/halfspec.h: pocketfft's real-to-complex transforms, halved along x): the
// same element-wise arithmetic per voxel (the ratio, the update, MATLAB's complex product), other FFTs (the
// last bits differ from MATLAB's; docs/decisions/0006-fft-backends.md). The element-wise steps are passes
// of their own; nothing depends on the thread count.
#include <omp.h>

#include <cfloat>
#include <stdexcept>

#include "decon.hpp"
#include "halfspec.h"

namespace pk5d::decon {
namespace {

using shape_t = std::vector<std::int64_t>;

struct Timer {
    Times& t;
    double t0 = 0;
    explicit Timer(Times& tt) : t(tt) {}
    void start() { t0 = omp_get_wtime(); }
    void stop(int k) { t.t[k] += omp_get_wtime() - t0; }
};

idx_t numel(const shape_t& s) {
    idx_t n = 1;
    for (auto v : s) n *= static_cast<idx_t>(v);
    return n;
}

// s .* h (conjugate: conj(h) .* s) in place, MATLAB's complex product (no fused multiply-add)
void multiply(cfloat* s, const cfloat* h, idx_t n, bool conjugate, int nt) {
    float* p = reinterpret_cast<float*>(s);
    const float* q = reinterpret_cast<const float*>(h);
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t i = 0; i < n; i++) {
        const float a = p[2 * i], b = p[2 * i + 1], cr = q[2 * i], ci = q[2 * i + 1];
        if (conjugate) {
            p[2 * i] = a * cr + b * ci;
            p[2 * i + 1] = b * cr - a * ci;
        } else {
            p[2 * i] = a * cr - b * ci;
            p[2 * i + 1] = a * ci + b * cr;
        }
    }
}

template <typename T>
inline float input_value(const T* raw, idx_t i, double bg) {
    const float v = static_cast<float>(static_cast<double>(raw[i]) - bg);
    return v > 0.0f ? v : 0.0f;
}

// r = max(I - bg, 0) ./ max(r, eps)
template <typename T>
void ratio(float* r, const T* raw, double bg, idx_t n, int nt) {
    const float epsf = static_cast<float>(DBL_EPSILON);
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t i = 0; i < n; i++) {
        const float v = r[i];
        const float cx = v > epsf ? v : epsf;
        r[i] = input_value(raw, i, bg) / cx;
    }
}

// out = max(v .* J, 0)
void update(float* out, const float* v, const float* j, idx_t n, int nt) {
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t i = 0; i < n; i++) {
        const float w = v[i] * j[i];
        out[i] = w > 0.0f ? w : 0.0f;
    }
}

// real(ifftn(fftn(x) .* h)) (conjugate: conj(h)) into out (may be x); s: the half spectrum's buffer
void convolve(const shape_t& shape, const float* x, const cfloat* h, bool conjugate, cfloat* s, float* out, idx_t nh,
              int nt, Timer& tm) {
    tm.start();
    pk5d::fft::native_forward(x, shape, reinterpret_cast<float*>(s), nt);
    tm.stop(1);
    tm.start();
    multiply(s, h, nh, conjugate, nt);
    pk5d::fft::native_inverse(reinterpret_cast<float*>(s), shape, out, nt);
    tm.stop(2);
}

}  // namespace

void native_otf(const shape_t& shape, const float* psf, const shape_t& ps, cfloat* h, int nt, Times& times) {
    Timer tm(times);
    if (shape.size() != ps.size() || shape.empty()) throw std::invalid_argument("native_otf: sizes");
    // as (z, y, x) with z = 1 for an image
    const idx_t nz = shape.size() == 3 ? shape[0] : 1, ny = shape[shape.size() - 2], nx = shape.back();
    const idx_t pz = ps.size() == 3 ? ps[0] : 1, py = ps[ps.size() - 2], px = ps.back();
    const idx_t n = nz * ny * nx;
    tm.start();
    Owned<float> pb(static_cast<size_t>(n));
    float* p = pb.p;
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t i = 0; i < n; i++) p[i] = 0.0f;
    // the PSF zero-padded to the data size, shifted by -floor(psfSize / 2) along every axis
    const idx_t cz = pz / 2, cy = py / 2, cx = px / 2;
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t z = 0; z < pz; z++) {
        const idx_t oz = (z - cz + nz) % nz;
        for (idx_t y = 0; y < py; y++) {
            float* row = p + (oz * ny + (y - cy + ny) % ny) * nx;
            const float* src = psf + (z * py + y) * px;
            for (idx_t x = 0; x < px; x++) row[(x - cx + nx) % nx] = src[x];
        }
    }
    tm.stop(3);
    tm.start();
    pk5d::fft::native_forward(p, shape, reinterpret_cast<float*>(h), nt);
    tm.stop(1);
}

template <typename T>
float* native_omw(const shape_t& shape, const T* raw, double bg, const cfloat* hf, const cfloat* hb, int n_iter,
                  const float* j0, int nt, Times& times) {
    Timer tm(times);
    const idx_t n = numel(shape);
    const idx_t nh = numel(pk5d::fft::native_half_shape(shape));
    Owned<float> rb(static_cast<size_t>(n));
    float* r = rb.p;
    if (n_iter == 0) {
        tm.start();
#pragma omp parallel for num_threads(nt) schedule(static)
        for (idx_t i = 0; i < n; i++) r[i] = j0 != nullptr ? j0[i] : input_value(raw, i, bg);
        tm.stop(3);
        return rb.release();
    }
    tm.start();
    Owned<float> jb(static_cast<size_t>(n));
    float* j = jb.p;
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t i = 0; i < n; i++) j[i] = j0 != nullptr ? j0[i] : input_value(raw, i, bg);
    Owned<cfloat> sb(static_cast<size_t>(nh));
    tm.stop(3);
    for (int k = 1; k <= n_iter; k++) {
        convolve(shape, j, hf, false, sb.p, r, nh, nt, tm);
        tm.start();
        ratio(r, raw, bg, n, nt);
        tm.stop(2);
        convolve(shape, r, hb, false, sb.p, r, nh, nt, tm);
        tm.start();
        update(k == n_iter ? r : j, r, j, n, nt);
        tm.stop(2);
    }
    return rb.release();
}

float* native_rlstep(const shape_t& shape, const float* y, const float* im, const cfloat* h, int nt, Times& times) {
    Timer tm(times);
    const idx_t n = numel(shape);
    const idx_t nh = numel(pk5d::fft::native_half_shape(shape));
    Owned<float> rb(static_cast<size_t>(n));
    Owned<cfloat> sb(static_cast<size_t>(nh));
    float* r = rb.p;
    convolve(shape, y, h, false, sb.p, r, nh, nt, tm);
    tm.start();
    ratio(r, im, 0.0, n, nt);
    tm.stop(2);
    convolve(shape, r, h, true, sb.p, r, nh, nt, tm);
    tm.start();
    update(r, r, y, n, nt);
    tm.stop(2);
    return rb.release();
}

template float* native_omw<uint16_t>(const shape_t&, const uint16_t*, double, const cfloat*, const cfloat*, int,
                                     const float*, int, Times&);
template float* native_omw<float>(const shape_t&, const float*, double, const cfloat*, const cfloat*, int,
                                  const float*, int, Times&);

}  // namespace pk5d::decon
