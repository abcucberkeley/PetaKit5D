// The resize passes, compiled once per instruction set (module.cpp includes this file with
// PK5D_RS_NS, the namespace, and PK5D_RS_ISA, the target attribute). Each output value is one
// sum in double, in tap order, from +0: acc = ((0 + w0 x0) + w1 x1) + ..., every tap
// multiplied (a weight of 0 included, so Inf * 0 gives NaN as in MATLAB), then rounded to the
// class (single: to nearest; integers: floor(acc + 0.5), or ceil(acc - 0.5) below 0, saturated,
// NaN to 0). No fused multiply-adds (-ffp-contract=off): the vector width changes no bits.
// A sum that comes out NaN is computed again by nan_sum, which gives MATLAB's NaN payload
// whatever the compiler's operand order and the processor's NaN rules.
namespace PK5D_RS_NS {

// MATLAB's NaN of a NaN sum (measured on R2024b): the first NaN in tap order stays, as the
// sum's first operand keeps its NaN: a NaN input's payload (quieted), MATLAB's NaN
// (0xFFF8000000000000) for Inf * 0 and for Inf - Inf. s points at the line's first value,
// stride apart.
template <typename T>
PK5D_RS_ISA double nan_sum(const T *s, i64 stride, const double *wi, const i64 *ii, i64 P) {
    double acc = 0.0;
    for (i64 k = 0; k < P; k++) {
        const double x = static_cast<double>(s[ii[k] * stride]);
        if (x != x) return quiet_nan(x);
        if (wi[k] == 0 && std::isinf(x)) return matlab_nan();
        const double r = acc + wi[k] * x;
        if (r != r) return matlab_nan();
        acc = r;
    }
    return acc;
}

template <typename T>
PK5D_RS_ISA inline T to_class(double a) {
    if constexpr (std::is_same_v<T, double>) {
        return a;
    } else if constexpr (std::is_same_v<T, float>) {
        return static_cast<float>(a);
    } else {
        constexpr double lo = static_cast<double>(std::numeric_limits<T>::min());
        constexpr double hi = static_cast<double>(std::numeric_limits<T>::max());
        if (!(a == a)) return T(0);
        const double r = a >= 0 ? std::floor(a + 0.5) : std::ceil(a - 0.5);
        if (r <= lo) return std::numeric_limits<T>::min();
        if (r >= hi) return std::numeric_limits<T>::max();
        return static_cast<T>(r);
    }
}

// along an axis of an array viewed as (outer, n_in, inner) -> (outer, n_out, inner); w and ind
// are n_out x P (ind 0-based along the axis)
template <typename T>
PK5D_RS_ISA void resize_lines(const T *src, T *dst, i64 outer, i64 n_in, i64 inner, i64 n_out,
                              const double *w, const i64 *ind, i64 P, int nt) {
    if (outer == 0 || n_out == 0 || inner == 0) return;
    if (inner == 1) {
        // along the contiguous axis: a gather per output value
#pragma omp parallel for schedule(static) num_threads(nt)
        for (i64 o = 0; o < outer; o++) {
            const T *s = src + o * n_in;
            T *d = dst + o * n_out;
            for (i64 i = 0; i < n_out; i++) {
                const double *wi = w + i * P;
                const i64 *ii = ind + i * P;
                double acc = 0.0;
                for (i64 k = 0; k < P; k++) acc = acc + wi[k] * static_cast<double>(s[ii[k]]);
                if (acc != acc) acc = nan_sum(s, 1, wi, ii, P);
                d[i] = to_class<T>(acc);
            }
        }
        return;
    }
    // along a strided axis: whole rows of `inner` values, in blocks that fit the cache
    constexpr i64 B = 1024;
    const i64 nb = (inner + B - 1) / B;
    const i64 total = outer * n_out * nb;
#pragma omp parallel num_threads(nt)
    {
        alignas(64) double acc[B];
#pragma omp for schedule(static)
        for (i64 t = 0; t < total; t++) {
            const i64 b = t % nb;
            const i64 oi = t / nb;
            const i64 i = oi % n_out;
            const i64 o = oi / n_out;
            const i64 j0 = b * B;
            const i64 len = std::min(B, inner - j0);
            for (i64 j = 0; j < len; j++) acc[j] = 0.0;
            const double *wi = w + i * P;
            const i64 *ii = ind + i * P;
            for (i64 k = 0; k < P; k++) {
                const double wk = wi[k];
                const T *s = src + (o * n_in + ii[k]) * inner + j0;
                for (i64 j = 0; j < len; j++) acc[j] = acc[j] + wk * static_cast<double>(s[j]);
            }
            bool nan = false;
            for (i64 j = 0; j < len; j++) nan |= acc[j] != acc[j];
            if (nan) {
                const T *line = src + o * n_in * inner + j0;
                for (i64 j = 0; j < len; j++)
                    if (acc[j] != acc[j]) acc[j] = nan_sum(line + j, inner, wi, ii, P);
            }
            T *d = dst + (o * n_out + i) * inner + j0;
            for (i64 j = 0; j < len; j++) d[j] = to_class<T>(acc[j]);
        }
    }
}

// index copies along an axis (the passes whose weights are one-hot): out(i) = in(idx[i])
template <typename T>
PK5D_RS_ISA void gather_lines(const T *src, T *dst, i64 outer, i64 n_in, i64 inner, i64 n_out,
                              const i64 *idx, int nt) {
    if (outer == 0 || n_out == 0 || inner == 0) return;
    if (inner == 1) {
#pragma omp parallel for schedule(static) num_threads(nt)
        for (i64 o = 0; o < outer; o++) {
            const T *s = src + o * n_in;
            T *d = dst + o * n_out;
            for (i64 i = 0; i < n_out; i++) d[i] = s[idx[i]];
        }
        return;
    }
    const i64 total = outer * n_out;
#pragma omp parallel for schedule(static) num_threads(nt)
    for (i64 t = 0; t < total; t++) {
        const i64 i = t % n_out;
        const i64 o = t / n_out;
        std::memcpy(dst + (o * n_out + i) * inner, src + (o * n_in + idx[i]) * inner,
                    static_cast<size_t>(inner) * sizeof(T));
    }
}

}  // namespace PK5D_RS_NS
