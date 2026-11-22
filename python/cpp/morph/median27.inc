// The 3 x 3 x 3 median network's line pass on x86-64 with GCC or Clang, included by median.cpp
// once per instruction set: in namespace PK5D_MEDIAN_NS, with the attribute PK5D_MEDIAN_ISA, on
// one register PK5D_MEDIAN_V of floats with the intrinsics PK5D_MEDIAN_OP(loadu_ps), (storeu_ps),
// (min_ps) and (max_ps). min_ps(y, x) is y < x ? y : x and max_ps(x, y) is x > y ? x : y, lane by
// lane and for every input: the selections of the generic form's compare-exchanges.
namespace PK5D_MEDIAN_NS {

// the medians of a line of n0 voxels along x whose windows' lines are rows[0 .. 8] (q = 3 dz + dx),
// read at i .. i + W + 1 by the pass at i
PK5D_MEDIAN_ISA void median27_line(const float* const* rows, float* o, const idx_t n0) {
    constexpr int W = (int)(sizeof(PK5D_MEDIAN_V) / sizeof(float));
    for (idx_t i0 = 0; i0 < n0; i0 += W) {
        PK5D_MEDIAN_V v[27];
        for (int q = 0; q < 9; q++) {
            for (int dy = 0; dy < 3; dy++) v[q * 3 + dy] = PK5D_MEDIAN_OP(loadu_ps)(rows[q] + i0 + dy);
        }
#define CX(a, b)                                    \
    {                                               \
        const PK5D_MEDIAN_V x = v[a], y = v[b];     \
        v[a] = PK5D_MEDIAN_OP(min_ps)(y, x);        \
        v[b] = PK5D_MEDIAN_OP(max_ps)(x, y);        \
    }
        PK5D_MEDIAN27_NETWORK(CX)
#undef CX
        if (n0 - i0 >= W) {
            PK5D_MEDIAN_OP(storeu_ps)(o + i0, v[13]);
        } else {
            float med[W];
            PK5D_MEDIAN_OP(storeu_ps)(med, v[13]);
            std::memcpy(o + i0, med, (size_t)(n0 - i0) * sizeof(float));
        }
    }
}

}  // namespace PK5D_MEDIAN_NS
