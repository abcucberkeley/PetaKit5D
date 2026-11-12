// The kernels of petakit5d._ext.mip, included by module.cpp once per instruction set: in
// namespace PK5D_MIP_NS, each kernel with the attribute PK5D_MIP_ISA (a target attribute, or
// nothing for the baseline). The loops (and the OpenMP regions GCC and Clang outline from
// them) are compiled for that set; the helpers of module.cpp are inlined into them. Every
// operation is a comparison and a selection, so each set gives the same bits.
namespace PK5D_MIP_NS {

template <typename T>
PK5D_MIP_ISA void max_pooling_3d(const T *in, T *out, u64 nz, u64 ny, u64 nx, u64 pz, u64 py, u64 px,
                    u64 mz, u64 my, u64 mx, int num_threads) {
    // output voxels along x per task: about 16 K input voxels per row
    const u64 cx = std::max<u64>(1, 16384 / px);
    const u64 nxc = (mx + cx - 1) / cx;
    const int64_t ntask = (int64_t)(nxc * my * mz);
    const u64 nbuf = px == 1 ? 0 : std::min(cx * px, nx);
#pragma omp parallel num_threads(num_threads)
    {
        std::unique_ptr<T[]> buffer(nbuf ? new T[nbuf] : nullptr);
#pragma omp for schedule(static)
        for (int64_t t = 0; t < ntask; t++) {
            const u64 xc = (u64)t % nxc;
            const u64 row = (u64)t / nxc;  // output row index y + my * z
            const u64 yo = row % my;
            const u64 zo = row / my;
            const u64 xo0 = xc * cx;
            const u64 xo1 = std::min(xo0 + cx, mx);
            const u64 x0 = xo0 * px;
            const u64 x1 = std::min(xo1 * px, nx);
            T *o = out + row * mx;
            // the row maximum goes into the output itself for px == 1
            T *r = px == 1 ? o + xo0 : buffer.get();
            const u64 n = x1 - x0;
            for (u64 i = 0; i < n; i++) r[i] = T(0);
            const u64 z1 = std::min((zo + 1) * pz, nz);
            const u64 y1 = std::min((yo + 1) * py, ny);
            for (u64 zi = zo * pz; zi < z1; zi++) {
                for (u64 yi = yo * py; yi < y1; yi++) {
                    const T *p = in + (zi * ny + yi) * nx + x0;
                    for (u64 i = 0; i < n; i++) r[i] = take_greater(r[i], p[i]);
                }
            }
            if (px == 1) continue;
            for (u64 xo = xo0; xo < xo1; xo++) {
                const u64 a = xo * px - x0;
                const u64 b = std::min((xo + 1) * px, nx) - x0;
                T m = T(0);
                for (u64 i = a; i < b; i++) m = take_greater(m, r[i]);
                o[xo] = m;
            }
        }
    }
}

template <typename T>
PK5D_MIP_ISA void range_max_3d(const T *im, T *out, u64 nz, u64 ny, u64 nx, int axis,
                  const std::vector<u64> &starts, const std::vector<u64> &ends,
                  int num_threads) {
    const int64_t nr = (int64_t)starts.size();
    if (axis == 2) {
        // along x: each row (z, y) is contiguous
        const int64_t nrow = (int64_t)(nz * ny);
#pragma omp parallel for schedule(static) num_threads(num_threads)
        for (int64_t c = 0; c < nrow; c++) {
            const T *row = im + c * nx;
            T *o = out + c * nr;
            for (int64_t r = 0; r < nr; r++) o[r] = fold_max(row + starts[r], ends[r] - starts[r]);
        }
    } else if (axis == 1) {
        // along y: per z plane and range, the rows of the range merged into one, in parts
        const u64 chunk = task_chunk(nx, nz * nr, num_threads, 256, 4096);
        const u64 nxc = (nx + chunk - 1) / chunk;
        const int64_t ntask = (int64_t)(nz * nr * nxc);
#pragma omp parallel for schedule(static) num_threads(num_threads)
        for (int64_t t = 0; t < ntask; t++) {
            const u64 xc = t % nxc;
            const u64 r = (t / nxc) % nr;
            const u64 z = t / (nxc * nr);
            const u64 x0 = xc * chunk;
            const u64 x1 = x0 + chunk < nx ? x0 + chunk : nx;
            T *o = out + (z * nr + r) * nx;
            const T *p = im + (z * ny + starts[r]) * nx;
            for (u64 x = x0; x < x1; x++) o[x] = p[x];
            for (u64 y = starts[r] + 1; y < ends[r]; y++) {
                p = im + (z * ny + y) * nx;
                for (u64 x = x0; x < x1; x++) o[x] = max_update(o[x], p[x]);
            }
        }
    } else {
        // along z: per range, the planes of the range merged into one, in parts (the mex's)
        const u64 nyx = ny * nx;
        const u64 chunk = 16384;
        const u64 nc = (nyx + chunk - 1) / chunk;
        const int64_t ntask = (int64_t)(nr * nc);
#pragma omp parallel for schedule(static) num_threads(num_threads)
        for (int64_t t = 0; t < ntask; t++) {
            const u64 ci = t % nc;
            const u64 r = t / nc;
            const u64 i0 = ci * chunk;
            const u64 i1 = i0 + chunk < nyx ? i0 + chunk : nyx;
            T *o = out + r * nyx;
            const T *p = im + starts[r] * nyx;
            for (u64 i = i0; i < i1; i++) o[i] = p[i];
            for (u64 z = starts[r] + 1; z < ends[r]; z++) {
                p = im + z * nyx;
                for (u64 i = i0; i < i1; i++) o[i] = max_update(o[i], p[i]);
            }
        }
    }
}

template <typename T>
PK5D_MIP_ISA T min_box_3d(const T *im, u64 ny, u64 nx, u64 z0, u64 y0, u64 x0, u64 z1, u64 y1, u64 x1,
             int num_threads) {
    constexpr int lanes = 32;
    const T start = std::numeric_limits<T>::max();
    const u64 by = y1 - y0, bz = z1 - z0;
    const u64 nrow = (x1 > x0) ? by * bz : 0;
    if (nrow == 0) return start;
    const u64 nparts = std::min<u64>((u64)num_threads, nrow);
    // (not std::vector: vector<bool> packs bits, which threads cannot write apart)
    std::unique_ptr<T[]> part(new T[nparts]);
    const u64 n = x1 - x0;
    const u64 nv = n / lanes * lanes;
#pragma omp parallel for schedule(static, 1) num_threads((int)nparts)
    for (int64_t q = 0; q < (int64_t)nparts; q++) {
        const u64 c0 = nrow * (u64)q / nparts;
        const u64 c1 = nrow * ((u64)q + 1) / nparts;
        T acc[lanes];
        for (int k = 0; k < lanes; k++) acc[k] = start;
        for (u64 c = c0; c < c1; c++) {
            const u64 y = y0 + c % by;
            const u64 z = z0 + c / by;
            const T *p = im + (z * ny + y) * nx + x0;
            for (u64 i = 0; i < nv; i += lanes) {
                for (int k = 0; k < lanes; k++) acc[k] = take_smaller(acc[k], p[i + k]);
            }
            for (u64 i = nv; i < n; i++) acc[0] = take_smaller(acc[0], p[i]);
        }
        T m = start;
        for (int k = 0; k < lanes; k++) m = take_smaller(m, acc[k]);
        part[q] = m;
    }
    T m = start;
    for (u64 q = 0; q < nparts; q++) m = take_smaller(m, part[q]);
    if constexpr (std::numeric_limits<T>::is_iec559) {
        if (m == T(0)) {
            // the first zero in MATLAB's order (strided; only for a zero minimum)
            for (u64 z = z0; z < z1; z++) {
                for (u64 x = x0; x < x1; x++) {
                    for (u64 y = y0; y < y1; y++) {
                        const T v = im[(z * ny + y) * nx + x];
                        if (v == T(0)) return v;
                    }
                }
            }
        }
    }
    return m;
}

}  // namespace PK5D_MIP_NS
