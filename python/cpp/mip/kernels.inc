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

// the range maxima along all three axes in one pass (range_max_3d along each axis): out0
// (nr0, n1, n2), out1 (n0, nr1, n2), out2 (n0, n1, nr2), every value range_max_3d's. The
// threads take runs of whole blocks of rows (axis 1); a block is visited plane after plane, so
// its axis-0 maxima (one plane's rows of the block) stay in the cache, and every value folds
// its range in increasing index: along axis 2 in fold_max, along axes 0 and 1 by max_update in
// the order of the planes and rows. A thread whose first row lies inside a range of axis 1
// that began in an earlier thread's rows folds its rows of that range into a part of its own;
// the parts are merged into the range's maximum in the order of the threads (of the rows),
// which gives the fold's value (max_update of the two folds is the fold of both runs).
template <typename T>
PK5D_MIP_ISA void range_max_axes(const T *im, T *out0, T *out1, T *out2, u64 n0, u64 n1, u64 n2,
                                 const axis_ranges &a0, const axis_ranges &a1,
                                 const axis_ranges &a2, int num_threads, u64 block_rows) {
    const u64 nr1 = a1.starts.size(), nr2 = a2.starts.size();
    // rows per block: about 128 KB of a plane (or the caller's)
    const u64 rb = std::min<u64>(
        n1, block_rows ? block_rows : std::max<u64>(1, 131072 / (n2 * sizeof(T))));
    const u64 nblocks = (n1 + rb - 1) / rb;
    const int nt = (int)std::min<u64>((u64)num_threads, nblocks);
    std::vector<u64> row0(nt + 1);
    for (int t = 0; t <= nt; t++) row0[t] = std::min(n1, nblocks * (u64)t / (u64)nt * rb);
    // the part of thread t (t >= 1): its rows of the range that holds its first row, if that
    // range began before them
    std::vector<char> has_part(nt, 0);
    for (int t = 1; t < nt; t++) has_part[t] = a1.starts[a1.of[row0[t]]] < row0[t];
    std::unique_ptr<T[]> parts(nt > 1 ? new T[(u64)(nt - 1) * n0 * n2] : nullptr);
    // (a run of rows per iteration: every run is done whatever the team's size)
#pragma omp parallel for schedule(static, 1) num_threads(nt)
    for (int t = 0; t < nt; t++) {
        const u64 y0 = row0[t], y1 = row0[t + 1];
        const u64 rp = a1.of[y0];
        const u64 yp = has_part[t] ? std::min(a1.ends[rp], y1) : y0;  // rows [y0, yp): the part
        T *part = has_part[t] ? parts.get() + (u64)(t - 1) * n0 * n2 : nullptr;
        for (u64 ya = y0; ya < y1; ya += rb) {
            const u64 yb = std::min(ya + rb, y1);
            for (u64 z = 0; z < n0; z++) {
                const u64 q = a0.of[z];
                const bool first0 = z == a0.starts[q];
                const T *src = im + (z * n1 + ya) * n2;
                T *acc = out0 + (q * n1 + ya) * n2;
                for (u64 y = ya; y < yb; y++) {
                    const T *row = src + (y - ya) * n2;
                    // the row of the next plane on its way (the block's rows of a plane are
                    // one run; the next plane's run lies a plane further)
                    if (z + 1 < n0) {
                        const char *next = reinterpret_cast<const char *>(row + n1 * n2);
                        for (u64 b = 0; b < n2 * sizeof(T); b += 64) prefetch_line(next + b);
                    }
                    T *o2 = out2 + (z * n1 + y) * nr2;
                    for (u64 r = 0; r < nr2; r++)
                        o2[r] = fold_max(row + a2.starts[r], a2.ends[r] - a2.starts[r]);
                    const u64 r = a1.of[y];
                    T *o1 = y < yp ? part + z * n2 : out1 + (z * nr1 + r) * n2;
                    const bool first1 = y == std::max(a1.starts[r], y0);
                    T *o0 = acc + (y - ya) * n2;
                    if (first0 && first1) {
                        for (u64 k = 0; k < n2; k++) o0[k] = o1[k] = row[k];
                    } else if (first0) {
                        for (u64 k = 0; k < n2; k++) {
                            o0[k] = row[k];
                            o1[k] = max_update(o1[k], row[k]);
                        }
                    } else if (first1) {
                        for (u64 k = 0; k < n2; k++) {
                            o0[k] = max_update(o0[k], row[k]);
                            o1[k] = row[k];
                        }
                    } else {
                        for (u64 k = 0; k < n2; k++) {
                            o0[k] = max_update(o0[k], row[k]);
                            o1[k] = max_update(o1[k], row[k]);
                        }
                    }
                }
            }
        }
    }
    if (nt < 2) return;
    // the parts, in the order of the threads
    const int64_t nz = (int64_t)n0;
#pragma omp parallel for schedule(static) num_threads(nt)
    for (int64_t z = 0; z < nz; z++) {
        for (int t = 1; t < nt; t++) {
            if (!has_part[t]) continue;
            const T *p = parts.get() + ((u64)(t - 1) * n0 + (u64)z) * n2;
            T *o = out1 + ((u64)z * nr1 + a1.of[row0[t]]) * n2;
            for (u64 k = 0; k < n2; k++) o[k] = max_update(o[k], p[k]);
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
