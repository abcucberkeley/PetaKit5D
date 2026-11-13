// The loops of petakit5d._ext.flatfield, included by module.cpp once per instruction set: in
// namespace PK5D_FF_NS, each function with the attribute PK5D_FF_ISA (a target attribute, or
// nothing for the baseline). Every step is an IEEE operation on one voxel (subtraction,
// division, addition, conversions between float, double and int32, comparisons and
// selections), so each set gives the same bits.
namespace PK5D_FF_NS {

// one run of n voxels along x: out = round_output(max((in - b) / f, 0) + offset or b)
template <typename T, typename U, bool Off>
PK5D_FF_ISA void ff_row(const T *__restrict in, const float *__restrict f,
                        const float *__restrict b, U *__restrict out, u64 n, float off) {
    for (u64 i = 0; i < n; i++) {
        float v = (static_cast<float>(in[i]) - b[i]) / f[i];
        v = PK5D_GE_SELECT(v, 0.0f);  // NaN gives +0, -0 stays
        out[i] = round_output<U>(Off ? PK5D_X86_NAN(v + off, v, off) : PK5D_X86_NAN(v + b[i], v, b[i]));
    }
}

// tasks [t0, t1) of a job (module.cpp): task t is the block t % nyb of rows of the z chunk
// t / nyb; the block's rows of the flat field and the background stay in cache over the chunk
template <typename T, typename U, bool Off>
PK5D_FF_ISA void ff_tasks_off(const ff_job<T, U> &j, u64 t0, u64 t1) {
    for (u64 t = t0; t < t1; t++) {
        const u64 yb = t % j.nyb, zc = t / j.nyb;
        const u64 y0 = yb * j.rows, y1 = std::min(y0 + j.rows, j.ny);
        const u64 z0 = zc * j.zc, z1 = std::min(z0 + j.zc, j.nz);
        const u64 n = (y1 - y0) * j.nx;
        const float *f = j.flat + y0 * j.nx;
        const float *b = j.bg + y0 * j.nx;
        for (u64 z = z0; z < z1; z++) {
            const u64 o = (z * j.ny + y0) * j.nx;
            ff_row<T, U, Off>(j.in + o, f, b, j.out + o, n, j.offset);
        }
    }
}

template <typename T, typename U>
PK5D_FF_ISA void ff_tasks(const ff_job<T, U> &j, u64 t0, u64 t1) {
    if (j.has_offset)
        ff_tasks_off<T, U, true>(j, t0, t1);
    else
        ff_tasks_off<T, U, false>(j, t0, t1);
}

}  // namespace PK5D_FF_NS
