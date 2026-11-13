// MATLAB's slice rule of sums (K08; N05's docs/decisions/0005-sums-and-math-library.md): the sum
// along a dimension that a non-singleton one precedes (sum(x, 2), sum(x, 3)).
//
// ops_slice_sum(x, out, chunked, center, num_threads):
//   x (O, R, I) C-contiguous, out (O, I) C-contiguous single or double, center None or (O, I) of
//   x's class. out[o, i] = x[o, 0, i] + x[o, 1, i] + ... + x[o, R-1, i] in index order, in the
//   class of out (single sums of single data only); chunked: in chunks of 1024 slices, each in
//   order, the chunk sums added in order ((c0 + c1) + c2) + ...; with a center the values added
//   are abs(x - center[o, i]).^2 in the class of x (var's sum; out of the same class). An empty
//   sum (R = 0) is 0. Each output is summed by one thread, so threads never change a result;
//   the payload of a NaN result is left to the caller.
//
// Inner sizes of 16 and more: tiles of 4096 outputs along i, every slice added to the whole tile
// (the additions across i are independent and vectorize). Smaller inner sizes: 16 outputs at a
// time, each a chain along r; for single data with an inner size of 1 (the port's sum(x, 2))
// four rows per vector through 4 x 4 transposes.
#include "ops_sum.hpp"

namespace pk5d::ops {
namespace {

using x3_in = nb::ndarray<nb::ro, nb::ndim<3>, nb::c_contig, nb::device::cpu>;
using o2_out = nb::ndarray<nb::ndim<2>, nb::c_contig, nb::device::cpu>;
using c2_in = nb::ndarray<nb::ro, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

constexpr idx_t kSliceChunk = 1024;
constexpr idx_t kTile = 4096;  // outputs per task along i (16 KB of single; 1024-16384 measured)
constexpr idx_t kGroup = 16;   // outputs per task for small inner sizes

template <typename T, typename Acc, bool Center>
inline Acc slice_value(T v, T c) {
    if constexpr (Center) {
        const T d = std::fabs(v - c);
        return d * d;
    } else {
        (void)c;
        return static_cast<Acc>(v);
    }
}

// outputs [i0, i0 + w) of row o
template <typename T, typename Acc, bool Center>
void slice_tile(const T *x, idx_t R, idx_t I, idx_t o, idx_t i0, idx_t w, Acc *out,
                const T *center, idx_t C, Acc *part, Acc *total) {
    const T *c = Center ? center + o * I + i0 : nullptr;
    auto row = [&](idx_t r) { return x + (o * R + r) * I + i0; };
    auto val = [&](const T *p, idx_t i) {
        return slice_value<T, Acc, Center>(p[i], Center ? c[i] : T(0));
    };
    bool have_total = false;
    for (idx_t r0 = 0; r0 < R; r0 += C) {
        // one chunk, slices [r0, r1): its first slice, then 4 slices per pass (the partial sum
        // stays in a register for the 4 additions, in order)
        const idx_t r1 = r0 + C < R ? r0 + C : R;
        const T *p0 = row(r0);
        for (idx_t i = 0; i < w; i++) part[i] = val(p0, i);
        idx_t r = r0 + 1;
        for (; r + 4 <= r1; r += 4) {
            const T *a = row(r), *b = row(r + 1), *d = row(r + 2), *e = row(r + 3);
            for (idx_t i = 0; i < w; i++) {
                Acc v = part[i];
                v = v + val(a, i);
                v = v + val(b, i);
                v = v + val(d, i);
                v = v + val(e, i);
                part[i] = v;
            }
        }
        for (; r < r1; r++) {
            const T *a = row(r);
            for (idx_t i = 0; i < w; i++) part[i] = part[i] + val(a, i);
        }
        if (r1 < R) {
            if (have_total) {
                for (idx_t i = 0; i < w; i++) total[i] = total[i] + part[i];
            } else {
                for (idx_t i = 0; i < w; i++) total[i] = part[i];
                have_total = true;
            }
        }
    }
    Acc *dst = out + o * I + i0;
    if (have_total) {
        for (idx_t i = 0; i < w; i++) dst[i] = total[i] + part[i];
    } else {
        for (idx_t i = 0; i < w; i++) dst[i] = part[i];
    }
}

// outputs q0 .. q0 + g - 1 (q = o * I + i), each a chain along r
template <typename T, typename Acc, bool Center>
void slice_group(const T *x, idx_t R, idx_t I, idx_t q0, idx_t g, Acc *out, const T *center,
                 idx_t C) {
    const T *base[kGroup];
    T cv[kGroup];
    Acc part[kGroup], total[kGroup];
    for (idx_t k = 0; k < g; k++) {
        const idx_t q = q0 + k, o = q / I, i = q % I;
        base[k] = x + o * R * I + i;
        cv[k] = Center ? center[q] : T(0);
    }
    bool have_total = false;
    for (idx_t r = 0; r < R; r++) {
        const idx_t off = r * I;
        if (r % C == 0) {
            if (r > 0) {
                for (idx_t k = 0; k < g; k++)
                    total[k] = have_total ? total[k] + part[k] : part[k];
                have_total = true;
            }
            for (idx_t k = 0; k < g; k++)
                part[k] = slice_value<T, Acc, Center>(base[k][off], cv[k]);
        } else {
            for (idx_t k = 0; k < g; k++)
                part[k] = part[k] + slice_value<T, Acc, Center>(base[k][off], cv[k]);
        }
    }
    for (idx_t k = 0; k < g; k++) out[q0 + k] = have_total ? total[k] + part[k] : part[k];
}

#ifdef PK5D_OPS_SSE2
// rows o0 .. o0 + 4 * ng - 1 of single data with an inner size of 1 (sum(x, 2) of the port's
// arrays: each row a chain along r), ng groups of 4 rows at a time: 4 x 4 blocks are transposed
// in registers, so that one vector addition adds the next slice of 4 rows, in order; the ng
// chains are independent (latency hidden). A chunk's sum starts at -0 (-0 + v is v).
template <bool Center>
void slice_rows_f(const float *x, idx_t R, idx_t o0, idx_t ng, float *out, const float *center,
                  idx_t C) {
    constexpr int kMax = 4;
    const __m128 sign = _mm_set1_ps(-0.0f), mz = _mm_set1_ps(-0.0f);
    __m128 part[kMax], total[kMax], cq[kMax];
    const float *row[kMax][4];
    for (idx_t g = 0; g < ng; g++) {
        for (int k = 0; k < 4; k++) row[g][k] = x + (o0 + 4 * g + k) * R;
        part[g] = mz;
        cq[g] = Center ? _mm_loadu_ps(center + o0 + 4 * g) : _mm_setzero_ps();
    }
    bool have_total = false;
    auto close_chunk = [&]() {
        for (idx_t g = 0; g < ng; g++) {
            total[g] = have_total ? _mm_add_ps(total[g], part[g]) : part[g];
            part[g] = mz;
        }
        have_total = true;
    };
    auto value = [&](__m128 v, idx_t g) {
        if constexpr (Center) {
            const __m128 d = _mm_andnot_ps(sign, _mm_sub_ps(v, cq[g]));
            return _mm_mul_ps(d, d);
        } else {
            (void)g;
            return v;
        }
    };
    idx_t r = 0;
    for (; r + 4 <= R; r += 4) {
        // chunk boundaries fall on multiples of 4 (C is 1024, or R without chunks)
        if (r > 0 && r % C == 0) close_chunk();
        for (idx_t g = 0; g < ng; g++) {
            __m128 a = _mm_loadu_ps(row[g][0] + r), b = _mm_loadu_ps(row[g][1] + r);
            __m128 c = _mm_loadu_ps(row[g][2] + r), d = _mm_loadu_ps(row[g][3] + r);
            _MM_TRANSPOSE4_PS(a, b, c, d);
            part[g] = _mm_add_ps(part[g], value(a, g));
            part[g] = _mm_add_ps(part[g], value(b, g));
            part[g] = _mm_add_ps(part[g], value(c, g));
            part[g] = _mm_add_ps(part[g], value(d, g));
        }
    }
    for (; r < R; r++) {
        if (r > 0 && r % C == 0) close_chunk();
        for (idx_t g = 0; g < ng; g++) {
            const __m128 v = _mm_setr_ps(row[g][0][r], row[g][1][r], row[g][2][r], row[g][3][r]);
            part[g] = _mm_add_ps(part[g], value(v, g));
        }
    }
    for (idx_t g = 0; g < ng; g++)
        _mm_storeu_ps(out + o0 + 4 * g, have_total ? _mm_add_ps(total[g], part[g]) : part[g]);
}
#endif

template <typename T, typename Acc, bool Center>
void slice_sum_t(const T *x, idx_t O, idx_t R, idx_t I, Acc *out, bool chunked, const T *center,
                 int num_threads, bool bind) {
    const idx_t nout = O * I;
    if (nout == 0) return;
    if (R == 0) {
        for (idx_t q = 0; q < nout; q++) out[q] = Acc(0);
        return;
    }
    const idx_t C = chunked ? kSliceChunk : R;
    const bool wide = I >= kGroup;
    const idx_t tiles = wide ? (I + kTile - 1) / kTile : 0;
    const idx_t tasks = wide ? O * tiles : (nout + kGroup - 1) / kGroup;
    const int nt = team(nout * R, num_threads);
    std::atomic<idx_t> next{0};
    parallel(nt < tasks ? nt : static_cast<int>(tasks), bind, [&](int, int) {
        std::vector<Acc> part(wide ? static_cast<size_t>(kTile) : 0);
        std::vector<Acc> total(wide ? static_cast<size_t>(kTile) : 0);
        for (;;) {
            const idx_t t = next.fetch_add(1, std::memory_order_relaxed);
            if (t >= tasks) break;
            if (wide) {
                const idx_t o = t / tiles, i0 = (t % tiles) * kTile;
                const idx_t w = I - i0 < kTile ? I - i0 : kTile;
                slice_tile<T, Acc, Center>(x, R, I, o, i0, w, out, center, C, part.data(),
                                           total.data());
            } else {
                const idx_t q0 = t * kGroup;
                const idx_t g = nout - q0 < kGroup ? nout - q0 : kGroup;
#ifdef PK5D_OPS_SSE2
                if constexpr (std::is_same_v<T, float> && std::is_same_v<Acc, float>) {
                    if (I == 1 && g >= 4) {
                        const idx_t ng = g / 4;
                        slice_rows_f<Center>(x, R, q0, ng, out, center, C);
                        if (g > 4 * ng)
                            slice_group<T, Acc, Center>(x, R, I, q0 + 4 * ng, g - 4 * ng, out,
                                                        center, C);
                        continue;
                    }
                }
#endif
                slice_group<T, Acc, Center>(x, R, I, q0, g, out, center, C);
            }
        }
    });
}

void slice_sum(const x3_in &x, o2_out &out, bool chunked, const c2_in &center, int num_threads,
               bool bind) {
    check_threads(num_threads);
    const idx_t O = static_cast<idx_t>(x.shape(0)), R = static_cast<idx_t>(x.shape(1)),
                I = static_cast<idx_t>(x.shape(2));
    if (static_cast<idx_t>(out.shape(0)) != O || static_cast<idx_t>(out.shape(1)) != I)
        throw std::invalid_argument("ops_slice_sum: out must have the shape (O, I) of x");
    const bool has_center = center.is_valid();
    if (has_center) {
        if (static_cast<idx_t>(center.shape(0)) != O || static_cast<idx_t>(center.shape(1)) != I)
            throw std::invalid_argument("ops_slice_sum: center must have the shape (O, I)");
        if (center.dtype() != x.dtype())
            throw std::invalid_argument("ops_slice_sum: center must have the class of x");
    }
    const mclass co = class_of(out.dtype());
    if (co != mclass::single && co != mclass::dbl)
        throw std::invalid_argument("ops_slice_sum: out must be single or double");
    const void *px = x.data();
    const void *pc = has_center ? center.data() : nullptr;
    void *po = out.data();
    nb::gil_scoped_release release;
    with_class(class_of(x.dtype()), "ops_slice_sum input", [&](auto ts) {
        using T = typename decltype(ts)::type;
        if constexpr (std::is_same_v<T, int64_t> || std::is_same_v<T, uint64_t>) {
            bad_class("ops_slice_sum input (64-bit integers)");
        } else {
            const T *xs = static_cast<const T *>(px);
            const T *cs = static_cast<const T *>(pc);
            if (co == mclass::single) {
                if constexpr (std::is_same_v<T, float>) {
                    float *o = static_cast<float *>(po);
                    if (has_center)
                        slice_sum_t<float, float, true>(xs, O, R, I, o, chunked, cs, num_threads,
                                                        bind);
                    else
                        slice_sum_t<float, float, false>(xs, O, R, I, o, chunked, nullptr,
                                                         num_threads, bind);
                } else {
                    throw std::invalid_argument("ops_slice_sum: single sums need single data");
                }
            } else {
                double *o = static_cast<double *>(po);
                if (has_center) {
                    if constexpr (std::is_same_v<T, double>)
                        slice_sum_t<double, double, true>(xs, O, R, I, o, chunked, cs,
                                                          num_threads, bind);
                    else
                        throw std::invalid_argument(
                            "ops_slice_sum: abs(x - center).^2 is summed in the class of x");
                } else {
                    slice_sum_t<T, double, false>(xs, O, R, I, o, chunked, nullptr, num_threads,
                                                  bind);
                }
            }
        }
    });
}

}  // namespace
}  // namespace pk5d::ops

PK5D_BINDINGS(bind_ops_slice) {
    namespace nb = nanobind;
    m.def("ops_slice_sum", &pk5d::ops::slice_sum, nb::arg("x").noconvert(),
          nb::arg("out").noconvert(), nb::arg("chunked"), nb::arg("center").noconvert().none(),
          nb::arg("num_threads"), nb::arg("bind") = true,
          "MATLAB's slice rule: out[o, i] = sum over r of x[o, r, i] in index order (chunks of "
          "1024 when chunked), of the values or of abs(x - center).^2.");
}
