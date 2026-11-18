// MATLAB's abs of complex single values (R02's kernel of the decon extension, moved here for
// the OTF threshold of omw_backprojector_generation and the FFT spectra): msun's hypotf
// (e_hypotf.c in single precision, measured), MATLAB's NaN for a NaN part.
//
// ops_abs_single(c, out, num_threads): c the n complex values as 2n interleaved singles (re, im),
//   out n singles.
#include "ops_common.hpp"

#include <algorithm>
#include <cstring>

namespace pk5d::ops {
namespace {

inline float f_of(int32_t w) {
    float f;
    std::memcpy(&f, &w, 4);
    return f;
}

// msun hypotf of |parts| with the high words ha >= hb (signs cleared), MATLAB's NaN for a NaN part
float m_hypotf(int32_t ha, int32_t hb) {
    constexpr int32_t HI = static_cast<int32_t>(0xFFFFF000u);
    if ((ha - hb) > 0x0F000000) return f_of(ha) + f_of(hb);
    if (ha >= 0x7F800000) return ha > 0x7F800000 ? f_of(0x7FC00000) : INFINITY;
    float a = f_of(ha), b = f_of(hb);
    int k = 0;
    if (ha > 0x58800000) {
        ha -= 0x22000000;
        hb -= 0x22000000;
        k = 68;
        a = f_of(ha);
        b = f_of(hb);
    }
    if (hb < 0x26800000) {
        if (hb == 0) return f_of(ha);
        if (hb <= 0x007FFFFF) {
            // the high words are not read again after this scaling
            a *= 0x1p126f;
            b *= 0x1p126f;
            k = -126;
        } else {
            ha += 0x22000000;
            hb += 0x22000000;
            k = -68;
            a = f_of(ha);
            b = f_of(hb);
        }
    }
    float w = a - b, h;
    if (w > b) {
        const float t1 = f_of(ha & HI);
        h = std::sqrt(t1 * t1 - (b * (-b) - (a - t1) * (a + t1)));
    } else {
        const float y1 = f_of(hb & HI), t1 = f_of((ha + 0x00800000) & HI);
        h = std::sqrt(t1 * y1 - (w * (-w) - (t1 * (b - y1) + ((a + a) - t1) * b)));
    }
    if (k != 0) h = f_of(0x3F800000 + (k << 23)) * h;
    return h;
}

void abs_single(const flat_in &c, flat_out &out, int num_threads, bool bind) {
    if (class_of(c.dtype()) != mclass::single || class_of(out.dtype()) != mclass::single)
        throw std::invalid_argument("ops_abs_single: c and out must be single");
    const idx_t n = static_cast<idx_t>(out.shape(0));
    if (static_cast<idx_t>(c.shape(0)) != 2 * n)
        throw std::invalid_argument("ops_abs_single: c must hold 2 values per output");
    const int nt = team(n, num_threads);
    const float *pc = static_cast<const float *>(c.data());
    float *po = static_cast<float *>(out.data());
    nb::gil_scoped_release release;
    parallel(nt, bind, [&](int k, int nk) {
        idx_t lo, hi;
        share(n, k, nk, lo, hi);
        for (idx_t i = lo; i < hi; i++) {
            int32_t hx, hy;
            std::memcpy(&hx, pc + 2 * i, 4);
            std::memcpy(&hy, pc + 2 * i + 1, 4);
            hx &= 0x7FFFFFFF;
            hy &= 0x7FFFFFFF;
            po[i] = m_hypotf(std::max(hx, hy), std::min(hx, hy));
        }
    });
}

}  // namespace
}  // namespace pk5d::ops

PK5D_BINDINGS(bind_ops_abs) {
    namespace nb = nanobind;
    m.def("ops_abs_single", &pk5d::ops::abs_single, nb::arg("c").noconvert(),
          nb::arg("out").noconvert(), nb::arg("num_threads"), nb::arg("bind") = true,
          "out = abs of the complex singles c (interleaved re, im): msun's hypotf, MATLAB's NaN "
          "for a NaN part.");
}
