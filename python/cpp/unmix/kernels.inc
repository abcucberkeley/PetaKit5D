// The loops of petakit5d._ext.unmix, included by module.cpp once per instruction set: in
// namespace PK5D_UX_NS, each function with the attribute PK5D_UX_ISA (a target attribute, or
// nothing for the baseline). Every step is an IEEE operation on one element (conversions between
// integers, float and double, subtraction, multiplication, addition, comparisons and
// selections), so each set gives the same bits.
namespace PK5D_UX_NS {

// acc = acc + (single(x) - b) * f, each single-double step in double, rounded to single once;
// a NaN made by a step (Inf - Inf, 0 * Inf) is x86's (PK5D_X86_NAN)
template <typename T>
PK5D_UX_ISA void add_channel(float *__restrict acc, const T *__restrict x, const u64 n,
                             const double b, const double f) {
    if (b != 0.0) {
        for (u64 i = 0; i < n; i++) {
            const double xd = (double)(float)x[i];
            const float v = (float)PK5D_X86_NAN(xd - b, xd, b);
            const double vd = (double)v;
            const float p = (float)PK5D_X86_NAN(vd * f, vd, f);
            acc[i] = PK5D_X86_NAN(acc[i] + p, acc[i], p);
        }
    } else {
        for (u64 i = 0; i < n; i++) {
            const double xd = (double)(float)x[i];
            const float p = (float)PK5D_X86_NAN(xd * f, xd, f);
            acc[i] = PK5D_X86_NAN(acc[i] + p, acc[i], p);
        }
    }
}

// MATLAB's cast of a single value: round half away from zero (the value plus 0.5 in double,
// truncated), saturate, NaN and values below 0.5 to 0; written as selections so the loops
// vectorize (the same values as the mex's branches)
template <typename O>
PK5D_UX_ISA inline O cast_to(const float v) {
    if constexpr (std::is_same_v<O, uint16_t> || std::is_same_v<O, uint8_t>) {
        constexpr double hi = std::is_same_v<O, uint16_t> ? 65535.0 : 255.0;
        double d = (double)v + 0.5;
        d = d >= 1.0 ? d : 0.0;  // NaN, negative values and values below 0.5
        d = d <= hi ? d : hi;
        return (O)(int32_t)d;
    } else {
        return (O)v;
    }
}

// elements [s, s + m): the sum of the channels from im0 (or 0), then finished or as it is
template <typename O>
PK5D_UX_ISA void block(const channel *chans, const u64 nc, const float *im0, O *out, float *a,
                       const u64 s, const u64 m, const double offset, const bool finish) {
    if (im0) {
        std::copy(im0 + s, im0 + s + m, a);
    } else {
        std::fill(a, a + m, 0.0f);
    }
    for (u64 k = 0; k < nc; k++) {
        const channel &c = chans[k];
        switch (c.c) {
            case cls::u8:
                add_channel(a, (const uint8_t *)c.data + s, m, c.background, c.factor);
                break;
            case cls::u16:
                add_channel(a, (const uint16_t *)c.data + s, m, c.background, c.factor);
                break;
            case cls::f32:
                add_channel(a, (const float *)c.data + s, m, c.background, c.factor);
                break;
            default:
                add_channel(a, (const double *)c.data + s, m, c.background, c.factor);
                break;
        }
    }
    O *o = out + s;
    if (!finish) {
        for (u64 i = 0; i < m; i++) o[i] = cast_to<O>(a[i]);
    } else if (offset != 0.0) {
        for (u64 i = 0; i < m; i++) {
            // max(0, v): NaN and -0 become +0
            const double vd = (double)(a[i] > 0.0f ? a[i] : 0.0f);
            o[i] = cast_to<O>((float)PK5D_X86_NAN(vd + offset, vd, offset));
        }
    } else {
        for (u64 i = 0; i < m; i++) o[i] = cast_to<O>(a[i] > 0.0f ? a[i] : 0.0f);
    }
}

// the blocks [b0, b1) of a call, with one scratch sum of a block
template <typename O>
PK5D_UX_ISA void blocks(const channel *chans, const u64 nc, const float *im0, O *out, const u64 n,
                        const u64 b0, const u64 b1, const double offset, const bool finish) {
    std::vector<float> acc(kBlock);
    for (u64 b = b0; b < b1; b++) {
        const u64 s = b * kBlock;
        block<O>(chans, nc, im0, out, acc.data(), s, std::min(kBlock, n - s), offset, finish);
    }
}

// out = x - b for single x (MATLAB's single minus double), elements [s, e); o may be x
PK5D_UX_ISA inline void subtract(const float *x, float *o, const double b,
                                 const u64 s, const u64 e) {
    for (u64 i = s; i < e; i++) {
        const double xd = (double)x[i];
        o[i] = (float)PK5D_X86_NAN(xd - b, xd, b);
    }
}

}  // namespace PK5D_UX_NS
