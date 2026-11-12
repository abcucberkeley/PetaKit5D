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

// Visit the output in tiles of BK planes x BJ columns of one row: the input voxels a tile reads
// (two rows, a patch of x and z for the transforms that keep y) stay in cache. Each voxel's value
// does not depend on the order. fn(k, i, j0, j1) computes out[k, i, j0:j1].
constexpr idx BK = 16, BJ = 64;

template <typename F>
void for_tiles(const Geometry& g, int threads, F&& fn) {
    const idx nkb = (g.nz + BK - 1) / BK, njb = (g.nx + BJ - 1) / BJ;
    const idx items = nkb * g.ny;  // (k block, row), rows fastest
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4) num_threads(threads)
#endif
    for (idx it = 0; it < items; ++it) {
        const idx i = it % g.ny, k0 = (it / g.ny) * BK;
        const idx k1 = k0 + BK < g.nz ? k0 + BK : g.nz;
        for (idx jb = 0; jb < njb; ++jb) {
            const idx j0 = jb * BJ, j1 = j0 + BJ < g.nx ? j0 + BJ : g.nx;
            for (idx k = k0; k < k1; ++k) fn(k, i, j0, j1);
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
    for_tiles(g, threads, [&](idx k, idx i, idx j0, idx j1) {
        const double If = tm.I[1][i], Is = tm.I[0][i], It = tm.I[2][i];
        const double Pf = P[1][k], Ps = P[0][k], Pt = P[2][k];
        const double* Jf = tm.J[1].data();
        const double* Js = tm.J[0].data();
        const double* Jt = tm.J[2].data();
        T* o = out + (k * g.ny + i) * g.nx;
        for (idx j = j0; j < j1; ++j) {
            const double qf = If + (Jf[j] + Pf);
            const double qs = Is + (Js[j] + Ps);
            const double qt = It + (Jt[j] + Pt);
            if (!(qf >= 0 && qf <= nf && qs >= 0 && qs <= ns && qt >= 0 && qt <= nt)) {
                o[j] = fill;
                continue;
            }
            if constexpr (Nearest) {
                const idx f = floor_nn(qf + 0.5);
                const idx s = floor_nn(qs + 0.5);
                const idx t = floor_nn(qt + 0.5);
                o[j] = in[(t * g.py + f) * g.px + s];
            } else {
                // q >= 0 here: floor is the truncation, ceil is floor + 1 unless q is whole
                const idx F0 = floor_nn(qf), S0 = floor_nn(qs), T0 = floor_nn(qt);
                const double fd = qf - static_cast<double>(F0), sd = qs - static_cast<double>(S0),
                             td = qt - static_cast<double>(T0);
                const idx F1 = F0 + (fd > 0), S1 = S0 + (sd > 0), T1 = T0 + (td > 0);
                const T* p00 = in + (T0 * g.py + F0) * g.px;  // (f0, t0)
                const T* p10 = in + (T0 * g.py + F1) * g.px;  // (f1, t0)
                const T* p01 = in + (T1 * g.py + F0) * g.px;  // (f0, t1)
                const T* p11 = in + (T1 * g.py + F1) * g.px;  // (f1, t1)
                const double oms = 1 - sd;
                const double c00 = static_cast<double>(p00[S0]) * oms + static_cast<double>(p00[S1]) * sd;
                const double c10 = static_cast<double>(p10[S0]) * oms + static_cast<double>(p10[S1]) * sd;
                const double c01 = static_cast<double>(p01[S0]) * oms + static_cast<double>(p01[S1]) * sd;
                const double c11 = static_cast<double>(p11[S0]) * oms + static_cast<double>(p11[S1]) * sd;
                const double omf = 1 - fd;
                const double c0 = c00 * omf + c10 * fd;
                const double c1 = c01 * omf + c11 * fd;
                o[j] = round_trunc_half<T>(c0 * (1 - td) + c1 * td);
            }
        }
    });
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
    for_tiles(g, threads, [&](idx k, idx i, idx j0, idx j1) {
        T* o = out + (k * g.ny + i) * g.nx;
        for (idx j = j0; j < j1; ++j) {
            double p[3];
            for (int a = 0; a < 3; ++a) p[a] = ((tm.J[a][j] + tm.I[a][i]) + tm.K[a][k]) + tc[3][a];
            if (!(p[0] >= 1 && p[0] <= n[0] && p[1] >= 1 && p[1] <= n[1] && p[2] >= 1 &&
                  p[2] <= n[2])) {
                o[j] = cast_like<T>(fill);
                continue;
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
    });
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
