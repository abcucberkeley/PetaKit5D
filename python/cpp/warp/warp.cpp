// petakit5d._ext.warp: MATLAB R2024b's imwarp for 3D affine transforms (linear, nearest, cubic),
// value for value.
//
// Arrays are (z, y, x), C-contiguous. tc is imwarp's composite matrix in MATLAB's convention:
// row vectors [x y z 1] of output subscripts times tc give input coordinates, tc[r][c] is
// MATLAB's tComp(r + 1, c + 1). The output voxel (z, y, x) has MATLAB's subscripts
// i = y + 1 (row), j = x + 1 (column), k = z + 1 (plane); its value depends only on them, so
// the order in which voxels are visited is free.
//
// linear and nearest (MATLAB's warp3d): the input coordinate along axis a (0: x, 1: y, 2: z)
// is 0-based, q = i*tc[1][a] + (j*tc[0][a] + (k*tc[2][a] + (tc[3][a] - 1))). Inside the
// volume when 0 <= q <= n - 1 on every axis, else the fill value. linear interpolates in
// double along x, then y, then z, each step a*(1 - w) + b*w with a and b at floor(q) and
// ceil(q); integer classes take the truncation of v - 0.5 below 0, else of v + 0.5, saturated
// (measured: not the floor that MATLAB's code-generation version of warp3d has). nearest takes
// the sample at floor(q + 0.5).
// cubic (MATLAB's interp3 'cubic', cubic convolution): the coordinate is 1-based,
// p = ((tc[0][a]*j + tc[1][a]*i) + tc[2][a]*k) + tc[3][a], inside when 1 <= p <= n. With
// m = floor(p) (n - 1 at p = n) and t = p - m the weights of samples m - 1 .. m + 2 are
// ((-t/2 + 1)t - 1/2)t, (3t/2 - 5/2)t^2 + 1, ((-3t/2 + 2)t + 1/2)t and (t/2 - 1/2)t^2, each sum
// of four taken left to right, along x, then y, then z. Samples one step outside the volume are
// 3*(f1 - f2) + f3 from the nearest three, made along x first, then along y from those values,
// then along z. The result is cast as MATLAB's cast: round half away from zero and saturate.
#include <pk5d/common.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

namespace nb = nanobind;

namespace {

using idx = std::int64_t;

using pk5d::cvolume;
using mat44 = nb::ndarray<const double, nb::shape<4, 4>, nb::c_contig, nb::device::cpu>;

struct Geometry {
    idx pz, py, px;  // input planes, rows, columns
    idx nz, ny, nx;  // output
};

// products of the subscripts with the matrix, as MATLAB forms them, per axis a
struct Terms {
    std::vector<double> I[3], J[3], K[3];
    Terms(const double (&tc)[4][4], const Geometry& g) {
        for (int a = 0; a < 3; ++a) {
            I[a].resize(g.ny);
            J[a].resize(g.nx);
            K[a].resize(g.nz);
            for (idx i = 0; i < g.ny; ++i) I[a][i] = static_cast<double>(i + 1) * tc[1][a];
            for (idx j = 0; j < g.nx; ++j) J[a][j] = static_cast<double>(j + 1) * tc[0][a];
            for (idx k = 0; k < g.nz; ++k) K[a][k] = static_cast<double>(k + 1) * tc[2][a];
        }
    }
};

// floor(x) for 0 <= x < 2^62: the truncation (std::floor is a library call on baseline x86-64)
inline idx floor_nn(double x) { return static_cast<idx>(x); }

template <typename T>
inline T round_trunc_half(double v) {
    // warp3d's conversion: the truncation of v + 0.5 (of v - 0.5 below 0), saturated
    if constexpr (std::is_floating_point_v<T>) {
        return static_cast<T>(v);
    } else {
        const double r = v < 0 ? v - 0.5 : v + 0.5;
        constexpr double lo = static_cast<double>(std::numeric_limits<T>::min());
        constexpr double hi = static_cast<double>(std::numeric_limits<T>::max());
        if (r <= lo - 1) return std::numeric_limits<T>::min();
        if (r >= hi + 1) return std::numeric_limits<T>::max();
        return static_cast<T>(static_cast<idx>(r));
    }
}

template <typename T>
inline T cast_like(double v) {
    // MATLAB's cast of a double: round half away from zero, saturate, NaN to 0
    if constexpr (std::is_floating_point_v<T>) {
        return static_cast<T>(v);
    } else {
        constexpr double lo = static_cast<double>(std::numeric_limits<T>::min());
        constexpr double hi = static_cast<double>(std::numeric_limits<T>::max());
        if (std::isnan(v)) return T(0);
        if (v <= lo) return std::numeric_limits<T>::min();
        if (v >= hi) return std::numeric_limits<T>::max();
        // round half away from zero, exactly (v - t is exact)
        idx t = static_cast<idx>(v);
        const double d = v - static_cast<double>(t);
        if (d >= 0.5)
            ++t;
        else if (d <= -0.5)
            --t;
        return static_cast<T>(t);
    }
}

// The order in which the output is visited (each voxel's value does not depend on it). In C
// order an input plane lies megabytes from the next, and for the transforms that rotate or
// resample x and z, consecutive voxels along x read other planes. So the voxels are visited by
// bands of input planes: an item is BK output planes of BI output rows; with one cursor per
// (row, plane) it visits the voxels whose input plane qt(k, i, j) falls in a band of BAND planes,
// band after band. qt is monotone in j (a sum of a multiple of j and constants), so a band is a
// run of j for each row and plane; voxels outside the input in z go with the first or last band.
// The few planes of a band and the rows of an item stay in cache, and their rows are read along
// x. dir is qt's step along j (its sign gives the order of the bands).
constexpr idx BI = 4, BK = 128, BAND = 8;

template <typename QT, typename F>
void for_voxels(const Geometry& g, int threads, double dir, QT&& qt, F&& fn) {
    const idx nkb = (g.nz + BK - 1) / BK, nib = (g.ny + BI - 1) / BI;
    const idx items = nkb * nib;  // row blocks fastest
    const idx nband = (g.pz + BAND - 1) / BAND;  // bands 0 .. nband - 1, and -1 / nband outside
#ifdef _OPENMP
#pragma omp parallel num_threads(threads)
#endif
    {
        std::vector<idx> cur(static_cast<size_t>(BI * BK));
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (idx it = 0; it < items; ++it) {
            const idx i0 = (it % nib) * BI, k0 = (it / nib) * BK;
            const idx i1 = i0 + BI < g.ny ? i0 + BI : g.ny, k1 = k0 + BK < g.nz ? k0 + BK : g.nz;
            const idx nk = k1 - k0;
            for (idx c = 0; c < (i1 - i0) * nk; ++c) cur[c] = 0;
            for (idx s = 0; s <= nband + 1; ++s) {
                // ascending qt: band b takes qt < (b + 1) * BAND (b = -1: below 0), the last
                // pass (b = nband) the rest; descending: band b takes qt >= b * BAND (b = nband
                // first: above the volume), the last pass (b = -1) the rest
                const idx b = dir >= 0 ? s - 1 : nband - s;
                const double lim = static_cast<double>(dir >= 0 ? (b + 1) * BAND : b * BAND);
                for (idx i = i0; i < i1; ++i)
                    for (idx k = k0; k < k1; ++k) {
                        idx& jc = cur[(i - i0) * nk + (k - k0)];
                        const idx j0 = jc;
                        if (j0 >= g.nx) continue;
                        // the run's end, the first j not in the band: from the linear estimate,
                        // then made exact by evaluating qt around it. The last pass takes all
                        // that is left, so that every voxel is visited once whatever qt is
                        // (+-inf or NaN where a nearly singular T overflows: those are fill)
                        auto in_band = [&](idx j) {
                            const double q = qt(k, i, j);
                            return dir >= 0 ? q < lim : q >= lim;
                        };
                        idx je;
                        if (s == nband + 1) {
                            je = g.nx;
                        } else if (dir == 0) {
                            je = in_band(j0) ? g.nx : j0;
                        } else {
                            const double q0 = qt(k, i, j0);
                            double est = static_cast<double>(j0) + (lim - q0) / dir;
                            // below j0, or NaN (from a NaN qt): j0
                            if (!(est >= static_cast<double>(j0))) est = static_cast<double>(j0);
                            je = est > static_cast<double>(g.nx) ? g.nx : static_cast<idx>(est);
                            while (je > j0 && !in_band(je - 1)) --je;
                            while (je < g.nx && in_band(je)) ++je;
                        }
                        for (idx j = j0; j < je; ++j) fn(k, i, j);
                        jc = je;
                    }
            }
        }
    }
}

template <typename T, bool Nearest>
void warp3d(const T* in, T* out, const double (&tc)[4][4], const Geometry& g, T fill,
            int threads) {
    const Terms tm(tc, g);
    std::vector<double> P[3];
    for (int a = 0; a < 3; ++a) {
        P[a].resize(g.nz);
        const double c = tc[3][a] - 1.0;
        for (idx k = 0; k < g.nz; ++k) P[a][k] = tm.K[a][k] + c;
    }
    const double nf = static_cast<double>(g.py - 1), ns = static_cast<double>(g.px - 1),
                 nt = static_cast<double>(g.pz - 1);
    const double *If = tm.I[1].data(), *Is = tm.I[0].data(), *It = tm.I[2].data();
    const double *Jf = tm.J[1].data(), *Js = tm.J[0].data(), *Jt = tm.J[2].data();
    const double *Pf = P[1].data(), *Ps = P[0].data(), *Pt = P[2].data();
    const idx plane = g.py * g.px;
    auto qt_of = [&](idx k, idx i, idx j) { return It[i] + (Jt[j] + Pt[k]); };
    auto voxel = [&](idx k, idx i, idx j) {
        T* o = out + (k * g.ny + i) * g.nx;
        {
            const double qf = If[i] + (Jf[j] + Pf[k]);
            const double qs = Is[i] + (Js[j] + Ps[k]);
            const double qt = qt_of(k, i, j);
            if (!(qf >= 0 && qf <= nf && qs >= 0 && qs <= ns && qt >= 0 && qt <= nt)) {
                o[j] = fill;
                return;
            }
            if constexpr (Nearest) {
                const idx f = floor_nn(qf + 0.5);
                const idx s = floor_nn(qs + 0.5);
                const idx t = floor_nn(qt + 0.5);
                o[j] = in[t * plane + f * g.px + s];
            } else {
                // q >= 0 here: floor is the truncation, ceil is floor + 1 unless q is whole
                const idx F0 = floor_nn(qf), S0 = floor_nn(qs), T0 = floor_nn(qt);
                const double fd = qf - static_cast<double>(F0), sd = qs - static_cast<double>(S0),
                             td = qt - static_cast<double>(T0);
                // the sample (f0, s0, t0) and the steps to f1, s1, t1 (one stride or none)
                const T* p00 = in + (T0 * plane + F0 * g.px + S0);
                const idx ds = sd > 0, df = fd > 0 ? g.px : 0, dt = td > 0 ? plane : 0;
                const T* p10 = p00 + df;  // (f1, t0)
                const T* p01 = p00 + dt;  // (f0, t1)
                const T* p11 = p01 + df;  // (f1, t1)
                const double oms = 1 - sd;
                const double c00 = static_cast<double>(p00[0]) * oms + static_cast<double>(p00[ds]) * sd;
                const double c10 = static_cast<double>(p10[0]) * oms + static_cast<double>(p10[ds]) * sd;
                const double c01 = static_cast<double>(p01[0]) * oms + static_cast<double>(p01[ds]) * sd;
                const double c11 = static_cast<double>(p11[0]) * oms + static_cast<double>(p11[ds]) * sd;
                const double omf = 1 - fd;
                const double c0 = c00 * omf + c10 * fd;
                const double c1 = c01 * omf + c11 * fd;
                o[j] = round_trunc_half<T>(c0 * (1 - td) + c1 * td);
            }
        }
    };
    for_voxels(g, threads, tc[0][2], qt_of, voxel);
}

inline void cubic_weights(double t, double (&w)[4]) {
    const double sq = t * t;
    w[0] = ((-0.5 * t + 1) * t - 0.5) * t;
    w[1] = (1.5 * t - 2.5) * sq + 1;
    w[2] = ((-1.5 * t + 2) * t + 0.5) * t;
    w[3] = (0.5 * t - 0.5) * sq;
}

inline double sum4(const double (&w)[4], const double (&v)[4]) {
    return ((w[0] * v[0] + w[1] * v[1]) + w[2] * v[2]) + w[3] * v[3];
}

// a sample of the volume extended by one step on every side (1-based x, y, z in 0 .. n + 1):
// the end values are made along x, then y, then z, each from the values of the step before
template <typename T>
struct Extended {
    const T* in;
    const Geometry& g;
    double v(idx x, idx y, idx z) const {  // all inside
        return static_cast<double>(in[((z - 1) * g.py + (y - 1)) * g.px + (x - 1)]);
    }
    double e1(idx x, idx y, idx z) const {
        if (x == 0) return 3 * (v(1, y, z) - v(2, y, z)) + v(3, y, z);
        if (x == g.px + 1) return 3 * (v(g.px, y, z) - v(g.px - 1, y, z)) + v(g.px - 2, y, z);
        return v(x, y, z);
    }
    double e2(idx x, idx y, idx z) const {
        if (y == 0) return 3 * (e1(x, 1, z) - e1(x, 2, z)) + e1(x, 3, z);
        if (y == g.py + 1) return 3 * (e1(x, g.py, z) - e1(x, g.py - 1, z)) + e1(x, g.py - 2, z);
        return e1(x, y, z);
    }
    double e3(idx x, idx y, idx z) const {
        if (z == 0) return 3 * (e2(x, y, 1) - e2(x, y, 2)) + e2(x, y, 3);
        if (z == g.pz + 1) return 3 * (e2(x, y, g.pz) - e2(x, y, g.pz - 1)) + e2(x, y, g.pz - 2);
        return e2(x, y, z);
    }
};

template <typename T>
void cubic3d(const T* in, T* out, const double (&tc)[4][4], const Geometry& g, double fill,
             int threads) {
    const Terms tm(tc, g);
    const Extended<T> ext{in, g};
    const double n[3] = {static_cast<double>(g.px), static_cast<double>(g.py),
                         static_cast<double>(g.pz)};
    const double* Jv[3] = {tm.J[0].data(), tm.J[1].data(), tm.J[2].data()};
    const double* Iv[3] = {tm.I[0].data(), tm.I[1].data(), tm.I[2].data()};
    const double* Kv[3] = {tm.K[0].data(), tm.K[1].data(), tm.K[2].data()};
    auto voxel = [&](idx k, idx i, idx j) {
        T* o = out + (k * g.ny + i) * g.nx;
        {
            double p[3];
            for (int a = 0; a < 3; ++a) p[a] = ((Jv[a][j] + Iv[a][i]) + Kv[a][k]) + tc[3][a];
            if (!(p[0] >= 1 && p[0] <= n[0] && p[1] >= 1 && p[1] <= n[1] && p[2] >= 1 &&
                  p[2] <= n[2])) {
                o[j] = cast_like<T>(fill);
                return;
            }
            idx m[3];
            double w[3][4];
            for (int a = 0; a < 3; ++a) {
                double f = static_cast<double>(floor_nn(p[a]));  // p >= 1 here
                if (f >= n[a]) f = n[a] - 1;
                m[a] = static_cast<idx>(f);
                cubic_weights(p[a] - f, w[a]);
            }
            const bool inner = m[0] >= 2 && m[0] + 2 <= g.px && m[1] >= 2 && m[1] + 2 <= g.py &&
                               m[2] >= 2 && m[2] + 2 <= g.pz;
            double lz[4];
            for (int c = 0; c < 4; ++c) {
                double ly[4];
                for (int b = 0; b < 4; ++b) {
                    double lx[4];
                    const idx y = m[1] - 1 + b, z = m[2] - 1 + c;
                    if (inner) {
                        const T* row = in + ((z - 1) * g.py + (y - 1)) * g.px + (m[0] - 2);
                        for (int a = 0; a < 4; ++a) lx[a] = static_cast<double>(row[a]);
                    } else {
                        for (int a = 0; a < 4; ++a) lx[a] = ext.e3(m[0] - 1 + a, y, z);
                    }
                    ly[b] = sum4(w[0], lx);
                }
                lz[c] = sum4(w[1], ly);
            }
            o[j] = cast_like<T>(sum4(w[2], lz));
        }
    };
    auto pt = [&](idx k, idx i, idx j) { return ((Jv[2][j] + Iv[2][i]) + Kv[2][k]) + tc[3][2]; };
    for_voxels(g, threads, tc[0][2], pt, voxel);
}

// A fresh output's pages are faulted in (and zeroed by the system) on their first write. The
// kernels write neighbouring rows from different threads, which then wait for each other's
// page faults; touching one byte per page first, each thread its own contiguous range, makes
// the faults parallel. The kernel writes every voxel afterwards, so values are not affected.
inline void prefault(void* p, size_t bytes, int threads) {
    unsigned char* b = static_cast<unsigned char*>(p);
    const idx page = 4096;
    const idx n = static_cast<idx>((bytes + page - 1) / page);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (idx q = 0; q < n; ++q) b[q * page] = 0;
}

template <typename T>
void warp_affine(const cvolume<const T>& vol, const mat44& tcm, cvolume<T>& out, int method, double fill,
                 int threads) {
    pk5d::check_threads(threads);
    Geometry g{static_cast<idx>(vol.shape(0)), static_cast<idx>(vol.shape(1)),
               static_cast<idx>(vol.shape(2)), static_cast<idx>(out.shape(0)),
               static_cast<idx>(out.shape(1)), static_cast<idx>(out.shape(2))};
    if (g.pz < 1 || g.py < 1 || g.px < 1) throw nb::value_error("empty input volume");
    if (method == 2 && (g.pz < 3 || g.py < 3 || g.px < 3))
        throw nb::value_error("cubic needs at least 3 samples along every axis");
    double tc[4][4];
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c) tc[r][c] = tcm(r, c);
    const T* in = vol.data();
    T* o = out.data();
    nb::gil_scoped_release release;
    if (g.nz * g.ny * g.nx == 0) return;
    prefault(o, static_cast<size_t>(g.nz * g.ny * g.nx) * sizeof(T), threads);
    switch (method) {
        case 0: warp3d<T, false>(in, o, tc, g, static_cast<T>(fill), threads); break;
        case 1: warp3d<T, true>(in, o, tc, g, static_cast<T>(fill), threads); break;
        case 2: cubic3d<T>(in, o, tc, g, fill, threads); break;
        default: throw std::invalid_argument("method: 0 linear, 1 nearest, 2 cubic");
    }
}

template <typename T>
void bind_type(nb::module_& m) {
    m.def("warp_affine", &warp_affine<T>, nb::arg("vol").noconvert(), nb::arg("tcomp").noconvert(),
          nb::arg("out").noconvert(), nb::arg("method"), nb::arg("fill"), nb::arg("num_threads"));
}

}  // namespace

NB_MODULE(warp, m) {
    m.doc() = "MATLAB R2024b's imwarp for 3D affine transforms, voxel for voxel";
    bind_type<std::uint8_t>(m);
    bind_type<std::int16_t>(m);
    bind_type<std::uint16_t>(m);
    bind_type<float>(m);
    bind_type<double>(m);
}
