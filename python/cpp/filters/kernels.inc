// The kernels of petakit5d._ext.filters, included by module.cpp once per instruction set: in
// namespace PK5D_FILT_NS, each kernel with the attribute PK5D_FILT_ISA (a target attribute with
// FMA, or nothing for the baseline). A fused multiply-add is std::fma, which rounds once on every
// processor and in every instruction set (with the FMA target it is the instruction, without it
// the C library's correctly rounded function), and a plain one is a product and a sum, each
// rounded (-ffp-contract=off): so every set gives the same bits.
namespace PK5D_FILT_NS {

// acc = fma(w, s, acc) or acc + w * s over n values
template <typename T>
PK5D_FILT_ISA void mac_row(T *acc, const T *s, T w, i64 n, bool fused) {
    if (fused) {
        for (i64 x = 0; x < n; x++) acc[x] = std::fma(w, s[x], acc[x]);
    } else {
        for (i64 x = 0; x < n; x++) {
            const T p = w * s[x];
            acc[x] = acc[x] + p;
        }
    }
}

template <typename T>
PK5D_FILT_ISA void mac_const(T *acc, T v, T w, i64 n, bool fused) {
    if (fused) {
        for (i64 x = 0; x < n; x++) acc[x] = std::fma(w, v, acc[x]);
    } else {
        const T p = w * v;
        for (i64 x = 0; x < n; x++) acc[x] = acc[x] + p;
    }
}

// One output row of a tap program over rows of n values: rows(q) gives the input row of padded
// index q (nullptr: the constant value). The chains of the program are added in their order; the
// taps after a -3 continue from the sum of the chains before it.
template <typename T, typename Rows>
PK5D_FILT_ISA void run_program(const int32_t *pr, i64 plen, i64 i, i64 pre, i64 n_in, const T *w,
                               T value, Rows rows, T *acc, T *total, i64 n) {
    bool first = true;
    bool open = false;
    for (i64 x = 0; x < n; x++) acc[x] = T(0);
    for (i64 k = 0; k <= plen; k++) {
        const int32_t op = k < plen ? pr[k] : int32_t(-2);
        if (op == -1 || op == -2 || op == -3) {
            if (open) {
                if (first) {
                    for (i64 x = 0; x < n; x++) total[x] = acc[x];
                    first = false;
                } else {
                    for (i64 x = 0; x < n; x++) total[x] = total[x] + acc[x];
                }
                for (i64 x = 0; x < n; x++) acc[x] = T(0);
                open = false;
            }
            if (op == -2) break;
            if (op == -3 && !first) {
                for (i64 x = 0; x < n; x++) acc[x] = total[x];
                first = true;
                open = true;
            }
            continue;
        }
        // a tap (fused or plain), or a fused tap of the padding just past one end with weight
        // w[k] (the edge value with replicate padding, the constant with a constant: the padded
        // taps at that end as one product)
        bool fused = op >= 0;
        i64 t, q;
        if (op <= -kEdgeHigh) {
            t = -(i64)op - kEdgeHigh;
            q = n_in;
            fused = true;
        } else if (op <= -kEdgeLow) {
            t = -(i64)op - kEdgeLow;
            q = -1;
            fused = true;
        } else {
            t = fused ? op : -(i64)op - 10;
            q = i + t - pre;
        }
        const T *s = rows(q);
        if (s != nullptr)
            mac_row(acc, s, w[t], n, fused);
        else
            mac_const(acc, value, w[t], n, fused);
        open = true;
    }
    if (first)
        for (i64 x = 0; x < n; x++) total[x] = T(0);
}

// the number of taps of a program that is one fused chain of kernel taps, else 0
inline i64 simple_taps(const int32_t *pr, i64 plen, i64 span) {
    i64 k = 0;
    for (; k < plen && pr[k] != -2; k++)
        if (pr[k] < 0 || pr[k] >= span) return 0;
    return k;
}

// nr (1 to 4) output rows i0.. of one fused chain over the taps pr[0 .. ntaps), each output
// row's values in the order of the chain; blocks of the row kept in registers
template <typename T, typename Rows>
PK5D_FILT_ISA void chain_rows(const int32_t *pr, i64 ntaps, i64 i0, i64 nr, i64 pre, const T *w,
                              T value, Rows rows, T *d, i64 dstride, i64 nx) {
    constexpr i64 XB = 64 / (i64)sizeof(T);
    std::vector<const T *> ptr((size_t)(4 * ntaps));
    for (i64 r = 0; r < nr; r++)
        for (i64 k = 0; k < ntaps; k++) ptr[r * ntaps + k] = rows(i0 + r + pr[k] - pre);
    for (i64 x0 = 0; x0 < nx; x0 += XB) {
        const i64 xb = std::min(XB, nx - x0);
        T acc[4][XB];
        for (i64 r = 0; r < 4; r++)
            for (i64 x = 0; x < XB; x++) acc[r][x] = T(0);
        if (xb == XB) {
            for (i64 k = 0; k < ntaps; k++) {
                const T wk = w[pr[k]];
                for (i64 r = 0; r < nr; r++) {
                    const T *s = ptr[r * ntaps + k];
                    if (s != nullptr) {
                        s += x0;
                        for (i64 x = 0; x < XB; x++) acc[r][x] = std::fma(wk, s[x], acc[r][x]);
                    } else {
                        for (i64 x = 0; x < XB; x++) acc[r][x] = std::fma(wk, value, acc[r][x]);
                    }
                }
            }
        } else {
            for (i64 k = 0; k < ntaps; k++) {
                const T wk = w[pr[k]];
                for (i64 r = 0; r < nr; r++) {
                    const T *s = ptr[r * ntaps + k];
                    for (i64 x = 0; x < xb; x++)
                        acc[r][x] = std::fma(wk, s != nullptr ? s[x0 + x] : value, acc[r][x]);
                }
            }
        }
        for (i64 r = 0; r < nr; r++)
            for (i64 x = 0; x < xb; x++) d[r * dstride + x0 + x] = acc[r][x];
    }
}

// n outputs along a padded line of one fused chain: out[x] = chain over the taps of
// w[t] * b[x + t], blocks of outputs kept in registers
template <typename T>
PK5D_FILT_ISA void chain_line(const int32_t *pr, i64 ntaps, const T *b, const T *w, T *d, i64 n) {
    constexpr i64 XB = 128 / (i64)sizeof(T);
    for (i64 x0 = 0; x0 < n; x0 += XB) {
        const i64 xb = std::min(XB, n - x0);
        T acc[XB];
        for (i64 x = 0; x < XB; x++) acc[x] = T(0);
        if (xb == XB) {
            for (i64 k = 0; k < ntaps; k++) {
                const T wk = w[pr[k]];
                const T *s = b + x0 + pr[k];
                for (i64 x = 0; x < XB; x++) acc[x] = std::fma(wk, s[x], acc[x]);
            }
        } else {
            for (i64 k = 0; k < ntaps; k++) {
                const T wk = w[pr[k]];
                const T *s = b + x0 + pr[k];
                for (i64 x = 0; x < xb; x++) acc[x] = std::fma(wk, s[x], acc[x]);
            }
        }
        for (i64 x = 0; x < xb; x++) d[x0 + x] = acc[x];
    }
}

// Tap programs along the middle axis of src (outer, n_in, inner) into dst (outer, n_out, inner).
template <typename T>
PK5D_FILT_ISA void line_program(const T *src, T *dst, i64 outer, i64 n_in, i64 inner, i64 n_out,
                                const T *w, const int32_t *prog, i64 plen, const int32_t *rowmap,
                                i64 pre, i64 span, int mode, T value, int num_threads) {
    if (outer == 0 || n_out == 0 || inner == 0) return;
    if (inner > 1) {
        const i64 tw = std::min<i64>(inner, 512);
        const i64 ntile = (inner + tw - 1) / tw;
        const i64 rc = 64;  // output rows per task (a multiple of the blocks of 4)
        const i64 nrc = (n_out + rc - 1) / rc;
        const i64 ntask = outer * ntile * nrc;
#pragma omp parallel num_threads(num_threads)
        {
            std::vector<T> acc((size_t)tw);
#pragma omp for schedule(static)
            for (i64 task = 0; task < ntask; task++) {
                const i64 o = task / (ntile * nrc);
                const i64 rem = task % (ntile * nrc);
                const i64 x0 = (rem / nrc) * tw;
                const i64 r0 = (rem % nrc) * rc;
                const i64 r1 = std::min(r0 + rc, n_out);
                const i64 nx = std::min(tw, inner - x0);
                const T *base = src + o * n_in * inner + x0;
                auto rows = [&](i64 q) -> const T * {
                    const i64 r = pad_index(q, n_in, mode);
                    return r < 0 ? nullptr : base + r * inner;
                };
                i64 i = r0;
                while (i < r1) {
                    const int32_t *pr = prog + (i64)rowmap[i] * plen;
                    const i64 nt_ = simple_taps(pr, plen, span);
                    if (nt_ > 0) {
                        // rows of one fused chain, four at a time with the sums in registers
                        i64 j = i + 1;
                        while (j < r1 && j - i < 4 && rowmap[j] == rowmap[i]) j++;
                        chain_rows(pr, nt_, i, j - i, pre, w, value, rows,
                                   dst + (o * n_out + i) * inner + x0, inner, nx);
                        i = j;
                        continue;
                    }
                    T *d = dst + (o * n_out + i) * inner + x0;
                    run_program(pr, plen, i, pre, n_in, w, value, rows, acc.data(), d, nx);
                    i++;
                }
            }
        }
        return;
    }
    // lines along the contiguous axis: each line padded into a buffer; a run of rows with one
    // fused chain is computed tap by tap across the run (vectorized), other rows one by one
    const i64 np = n_out + span - 1;  // padded indices i - pre .. i - pre + span - 1
#pragma omp parallel num_threads(num_threads)
    {
        std::vector<T> buf((size_t)np), acc((size_t)std::max<i64>(n_out, 1));
#pragma omp for schedule(static)
        for (i64 l = 0; l < outer; l++) {
            const T *line = src + l * n_in;
            for (i64 p = 0; p < np; p++) {
                const i64 r = pad_index(p - pre, n_in, mode);
                buf[p] = r < 0 ? value : line[r];
            }
            auto rows = [&](i64 q) -> const T * {
                const i64 p = q + pre;
                return p < 0 || p >= np ? nullptr : buf.data() + p;
            };
            T *d = dst + l * n_out;
            i64 i = 0;
            while (i < n_out) {
                const int32_t prow = rowmap[i];
                i64 j = i + 1;
                while (j < n_out && rowmap[j] == prow) j++;
                const int32_t *pr = prog + (i64)prow * plen;
                bool simple = true;
                for (i64 k = 0; k < plen && pr[k] != -2; k++)
                    if (pr[k] < 0 || pr[k] >= span) simple = false;
                if (simple) {
                    chain_line(pr, simple_taps(pr, plen, span), buf.data() + i, w, d + i, j - i);
                } else {
                    for (i64 r = i; r < j; r++) {
                        T tmp;
                        run_program(pr, plen, r, pre, n_in, w, value, rows, acc.data(), &tmp, 1);
                        d[r] = tmp;
                    }
                }
                i = j;
            }
        }
    }
}

// out[z, y, x] = the chain over the taps k, in their order, of w[k] * P[z + dz, y + dy, x + dx]
// (P padded by the caller), fused or plain.
template <typename T>
PK5D_FILT_ISA void taps_nd(const T *P, i64 py, i64 px, T *out, i64 oz, i64 oy, i64 ox,
                           const i64 *offs, const T *w, i64 K, bool fused, int num_threads) {
    if (oz == 0 || oy == 0 || ox == 0) return;
    const i64 tw = std::min<i64>(ox, 1024);
    const i64 ntile = (ox + tw - 1) / tw;
    const i64 ntask = oz * oy * ntile;
#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (i64 task = 0; task < ntask; task++) {
        const i64 row = task / ntile;
        const i64 x0 = (task % ntile) * tw;
        const i64 nx = std::min(tw, ox - x0);
        const i64 z = row / oy, y = row % oy;
        T *o = out + row * ox + x0;
        for (i64 x = 0; x < nx; x++) o[x] = T(0);
        for (i64 k = 0; k < K; k++) {
            const T *s = P + ((z + offs[3 * k]) * py + (y + offs[3 * k + 1])) * px + offs[3 * k + 2] + x0;
            mac_row(o, s, w[k], nx, fused);
        }
    }
}

}  // namespace PK5D_FILT_NS
