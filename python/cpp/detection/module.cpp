// petakit5d._ext.detection: the Gaussian fits of the 3-D point-source detection (task C02),
// written from the published method (Aguet et al. 2013) and MATLAB R2024b's behaviour measured as
// a black box (python/docs/detection_scope.md; no GPL or MathWorks source):
// - fit_cd: XR_fitGaussian3D's 'CD_Newton' fit of A exp(-((x-x0)^2 + (y-y0)^2)/(2 s^2)
//   - (z-z0)^2/(2 sz^2)) + c to windows (sigmas fixed): Newton steps on the position (the
//   gradient and Hessian of the RSS as compute_gradient_hessian_mex sums them, in voxel order),
//   a damped step where the Hessian is not positive definite, the step cut to the window's box,
//   an Armijo line search (1/4 of the predicted decrease, at most 27 halvings), then A by least
//   squares at the new position and c from the new A; restarts from the window's LoG maxima
//   (python/src/petakit5d/detection/fit.py); the result's covariance from the Hessian in
//   (x, y, z, A, c) and the residual tests.
// - ad_stat: the Anderson-Darling statistic of a sample in its given order (adtest_mex's).
// Sums over a window's voxels are MATLAB's (index order up to 43 values, four lanes above), the
// exponential MATLAB's (fdlibm 5.3, e_exp.c: Sun Microsystems' permissive notice). Every fit runs
// on one thread in a fixed order, so results do not depend on the thread count.
#include <pk5d/common.hpp>

#include "ad_libm.hpp"
#include "mexfit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nb = nanobind;
using i64 = int64_t;

namespace {

// -- fdlibm 5.3's __ieee754_exp (MATLAB's exp; exp(1) is e rounded, as MATLAB's) -------------
// Its notice, e_exp.c's:
//   Copyright (C) 1993, 2004 by Sun Microsystems, Inc. All rights reserved.
//   Developed at SunSoft, a Sun Microsystems, Inc. business.
//   Permission to use, copy, modify, and distribute this software is freely granted, provided
//   that this notice is preserved.

inline uint32_t hi_word(double x) {
    uint64_t u;
    std::memcpy(&u, &x, 8);
    return (uint32_t)(u >> 32);
}
inline double with_hi(double x, uint32_t hi) {
    uint64_t u;
    std::memcpy(&u, &x, 8);
    u = (u & 0xffffffffull) | ((uint64_t)hi << 32);
    std::memcpy(&x, &u, 8);
    return x;
}

double fdlibm_exp(double x) {
    static const double halF[2] = {0.5, -0.5}, o_threshold = 7.09782712893383973096e+02,
                        u_threshold = -7.45133219101941108420e+02,
                        ln2HI[2] = {6.93147180369123816490e-01, -6.93147180369123816490e-01},
                        ln2LO[2] = {1.90821492927058770002e-10, -1.90821492927058770002e-10},
                        invln2 = 1.44269504088896338700e+00, P1 = 1.66666666666666019037e-01,
                        P2 = -2.77777777770155933842e-03, P3 = 6.61375632143793436117e-05,
                        P4 = -1.65339022054652515390e-06, P5 = 4.13813679705723846039e-08,
                        twom1000 = 9.33263618503218878990e-302;
    if (x == 1.0) return 2.718281828459045;
    double y, hi = 0.0, lo = 0.0, c, t;
    int k = 0;
    uint32_t hx = hi_word(x);
    const int xsb = (hx >> 31) & 1;
    hx &= 0x7fffffff;
    if (hx >= 0x40862E42) {  // |x| >= 709.78...
        if (hx >= 0x7ff00000) {
            if (std::isnan(x)) return x + x;
            return xsb == 0 ? x : 0.0;
        }
        if (x > o_threshold) return std::numeric_limits<double>::infinity();
        if (x < u_threshold) return 0.0;
    }
    if (hx > 0x3fd62e42) {  // |x| > 0.5 ln2
        if (hx < 0x3FF0A2B2) {  // and < 1.5 ln2
            hi = x - ln2HI[xsb];
            lo = ln2LO[xsb];
            k = 1 - xsb - xsb;
        } else {
            k = (int)(invln2 * x + halF[xsb]);
            t = k;
            hi = x - t * ln2HI[0];
            lo = t * ln2LO[0];
        }
        x = hi - lo;
    } else if (hx < 0x3e300000) {  // |x| < 2^-28
        return 1.0 + x;
    } else {
        k = 0;
    }
    t = x * x;
    c = x - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
    if (k == 0) return 1.0 - ((x * c) / (c - 2.0) - x);
    y = 1.0 - ((lo - (x * c) / (2.0 - c)) - hi);
    const uint32_t hy = hi_word(y);
    if (k >= -1021) return with_hi(y, hy + ((uint32_t)k << 20));
    return with_hi(y, hy + ((uint32_t)(k + 1000) << 20)) * twom1000;
}

// -- MATLAB's sum of a vector: index order up to 43 values, else four lanes ((l1+l2)+l3)+l4 ---

template <typename F>
inline double msum(int n, F f) {
    if (n <= 43) {
        double s = 0.0;
        for (int i = 0; i < n; i++) s += f(i);
        return s;
    }
    double l[4] = {f(0), f(1), f(2), f(3)};
    for (int i = 4; i < n; i++) l[i & 3] += f(i);
    return ((l[0] + l[1]) + l[2]) + l[3];
}

// -- small dense linear algebra ----------------------------------------------------------------

// the inverse of an n x n matrix (Gauss-Jordan with partial pivoting); NaN if singular
void inverse(int n, const double *A, double *inv) {
    std::vector<double> M(A, A + n * n);
    for (int i = 0; i < n * n; i++) inv[i] = 0.0;
    for (int i = 0; i < n; i++) inv[i * n + i] = 1.0;
    for (int k = 0; k < n; k++) {
        int p = k;
        for (int i = k + 1; i < n; i++)
            if (std::fabs(M[i * n + k]) > std::fabs(M[p * n + k])) p = i;
        if (M[p * n + k] == 0.0) {
            for (int i = 0; i < n * n; i++) inv[i] = std::numeric_limits<double>::quiet_NaN();
            return;
        }
        if (p != k)
            for (int j = 0; j < n; j++) {
                std::swap(M[k * n + j], M[p * n + j]);
                std::swap(inv[k * n + j], inv[p * n + j]);
            }
        const double d = M[k * n + k];
        for (int j = 0; j < n; j++) {
            M[k * n + j] /= d;
            inv[k * n + j] /= d;
        }
        for (int i = 0; i < n; i++) {
            if (i == k) continue;
            const double f = M[i * n + k];
            if (f == 0.0) continue;
            for (int j = 0; j < n; j++) {
                M[i * n + j] -= f * M[k * n + j];
                inv[i * n + j] -= f * inv[k * n + j];
            }
        }
    }
}

// the eigenvalues of a symmetric 3 x 3 matrix (cyclic Jacobi)
void eig_sym3(const double A[9], double w[3]) {
    double a[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) a[i][j] = A[i * 3 + j];
    for (int sweep = 0; sweep < 50; sweep++) {
        const double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        const double dia = a[0][0] * a[0][0] + a[1][1] * a[1][1] + a[2][2] * a[2][2];
        if (off <= 1e-300 || off <= dia * 1e-34) break;
        for (int p = 0; p < 2; p++)
            for (int q = p + 1; q < 3; q++) {
                if (a[p][q] == 0.0) continue;
                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < 3; k++) {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; k++) {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
    }
    w[0] = a[0][0];
    w[1] = a[1][1];
    w[2] = a[2][2];
}

// -- the fit -------------------------------------------------------------------------------------

struct Window {
    int n;  // valid voxels
    const double *d, *X, *Y, *Z;
    double lo[3], hi[3];  // the box: 0 (the window's first voxel) to the last valid voxel per axis
};

struct State {
    double x[3], A, c;
};

// the model, the residuals and the RSS at a position (amplitude and background given)
struct Eval {
    std::vector<double> g, r;
    double rss;
};

inline void model(const Window &w, const double p[3], double s, double sz, std::vector<double> &g) {
    const double ds = 2 * (s * s), dz = 2 * (sz * sz);
    for (int i = 0; i < w.n; i++) {
        const double ex = w.X[i] - p[0], ey = w.Y[i] - p[1], ez = w.Z[i] - p[2];
        g[i] = fdlibm_exp(-(ex * ex + ey * ey) / ds - (ez * ez) / dz);
    }
}

inline double rss_of(const Window &w, const std::vector<double> &g, double A, double c) {
    return msum(w.n, [&](int i) {
        const double r = (w.d[i] - A * g[i]) - c;
        return r * r;
    });
}

// gradient and Hessian of the RSS in (x, y, z): compute_gradient_hessian_mex's sums (sequential)
void grad_hess(const Window &w, const double p[3], double s, double sz, const std::vector<double> &g,
               double A, double c, double G[3], double H[9]) {
    double gx = 0, gy = 0, gz = 0, hxx = 0, hyy = 0, hzz = 0, hxy = 0, hxz = 0, hyz = 0;
    const double s2 = s * s, sz2 = sz * sz;
    for (int i = 0; i < w.n; i++) {
        const double Ag = A * g[i];
        const double r = (w.d[i] - Ag) - c;
        const double dx = (w.X[i] - p[0]) / s2, dy = (w.Y[i] - p[1]) / s2, dz = (w.Z[i] - p[2]) / sz2;
        const double rA = r * Ag;
        gx += rA * dx;
        gy += rA * dy;
        gz += rA * dz;
        const double ww = Ag * Ag - rA;
        hxx += ww * dx * dx + rA / s2;
        hyy += ww * dy * dy + rA / s2;
        hzz += ww * dz * dz + rA / sz2;
        hxy += ww * dx * dy;
        hxz += ww * dx * dz;
        hyz += ww * dy * dz;
    }
    G[0] = -2 * gx;
    G[1] = -2 * gy;
    G[2] = -2 * gz;
    H[0] = 2 * hxx, H[1] = 2 * hxy, H[2] = 2 * hxz;
    H[3] = 2 * hxy, H[4] = 2 * hyy, H[5] = 2 * hyz;
    H[6] = 2 * hxz, H[7] = 2 * hyz, H[8] = 2 * hzz;
}

// x = M \ b for a 3 x 3 M: LU with partial pivoting (the first largest pivot), divisions,
// no fused multiply-adds, column-oriented substitutions (the order MATLAB's results show)
void lu3(const double M[9], const double b[3], double x[3]) {
    double a[3][3] = {{M[0], M[1], M[2]}, {M[3], M[4], M[5]}, {M[6], M[7], M[8]}};
    double v[3] = {b[0], b[1], b[2]};
    for (int k = 0; k < 3; k++) {
        int p = k;
        for (int i = k + 1; i < 3; i++)
            if (std::fabs(a[i][k]) > std::fabs(a[p][k])) p = i;
        if (p != k) {
            for (int j = 0; j < 3; j++) std::swap(a[k][j], a[p][j]);
            std::swap(v[k], v[p]);
        }
        const double piv = a[k][k];
        for (int i = k + 1; i < 3; i++) a[i][k] = a[i][k] / piv;
        for (int i = k + 1; i < 3; i++)
            for (int j = k + 1; j < 3; j++) a[i][j] = a[i][j] - a[i][k] * a[k][j];
    }
    for (int j = 0; j < 3; j++)
        for (int i = j + 1; i < 3; i++) v[i] = v[i] - a[i][j] * v[j];
    for (int j = 2; j >= 0; j--) {
        v[j] = v[j] / a[j][j];
        for (int i = 0; i < j; i++) v[i] = v[i] - a[i][j] * v[j];
    }
    x[0] = v[0], x[1] = v[1], x[2] = v[2];
}

// the step direction: Newton when H's leading minors are positive (Sylvester); when every
// eigenvalue is negative -(-H + 10 I) \ G; otherwise -(H + lam I) \ G with
// lam = |min eig| + max(min |eig| / 2, 10)
void direction(const double G[3], const double H[9], double d[3]) {
    const double m2 = H[0] * H[4] - H[1] * H[3];
    const double m3 = H[0] * (H[4] * H[8] - H[5] * H[7]) - H[1] * (H[3] * H[8] - H[5] * H[6]) +
                      H[2] * (H[3] * H[7] - H[4] * H[6]);
    double M[9];
    std::memcpy(M, H, sizeof(M));
    if (!(H[0] > 0 && m2 > 0 && m3 > 0)) {
        double w[3];
        eig_sym3(H, w);
        if (w[0] < 0 && w[1] < 0 && w[2] < 0) {
            for (int k = 0; k < 9; k++) M[k] = -H[k];
            M[0] += 10.0, M[4] += 10.0, M[8] += 10.0;
        } else {
            double mn = w[0], mabs = std::fabs(w[0]);
            for (int k = 1; k < 3; k++) {
                mn = std::min(mn, w[k]);
                mabs = std::min(mabs, std::fabs(w[k]));
            }
            const double lam = std::fabs(mn) + std::max(0.5 * mabs, 10.0);
            M[0] += lam, M[4] += lam, M[8] += lam;
        }
    }
    double x[3];
    lu3(M, G, x);
    d[0] = -x[0], d[1] = -x[1], d[2] = -x[2];
}

// -- the restarts: starts from the window's LoG maxima and bright components ------------------

// the 3 x 3 x 3 Laplacian of a Gaussian (fspecial3('log', [3 3 3], [s s sz])), [y][x][z]
void log_kernel(double s, double sz, double K[27]) {
    double h[27], sum = 0.0;
    for (int yi = 0; yi < 3; yi++)
        for (int xi = 0; xi < 3; xi++)
            for (int zi = 0; zi < 3; zi++) {
                const double X = xi - 1.0, Y = yi - 1.0, Z = zi - 1.0;
                h[(yi * 3 + xi) * 3 + zi] =
                    fdlibm_exp(-(X * X / (2 * s * s) + Y * Y / (2 * s * s) + Z * Z / (2 * sz * sz)));
            }
    for (int q = 0; q < 27; q++) sum += h[q];
    double tot = 0.0;
    for (int yi = 0; yi < 3; yi++)
        for (int xi = 0; xi < 3; xi++)
            for (int zi = 0; zi < 3; zi++) {
                const double X = xi - 1.0, Y = yi - 1.0, Z = zi - 1.0;
                const int q = (yi * 3 + xi) * 3 + zi;
                K[q] = h[q] / sum *
                       (X * X / (s * s * s * s) + Y * Y / (s * s * s * s) + Z * Z / (sz * sz * sz * sz) -
                        (2 / (s * s) + 1 / (sz * sz)));
                tot += K[q];
            }
    for (int q = 0; q < 27; q++) K[q] -= tot / 27;
}

struct FullWindow {
    const double *v;  // MATLAB's column-major order, NaN for voxels left out
    int ny, nx, nz;
    inline double at(int y, int x, int z) const { return v[y + ny * (x + nx * z)]; }
};

// the starts in window coordinates [x y z], at most 16: the local maxima (5^3, ties kept) of
// -LoG(window) on bright voxels (above the floor(0.8 n)-th smallest value) by decreasing value,
// topped up to 3 with positive ones below, without p0; then the value-weighted centroids of the
// bright 18-connected components of 2 voxels or more, larger ones first
std::vector<std::array<double, 3>> starts(const FullWindow &W, const double p0[3], double s, double sz) {
    const int ny = W.ny, nx = W.nx, nz = W.nz, N = ny * nx * nz;
    double K[27];
    log_kernel(s, sz, K);
    std::vector<double> L(N);
    auto cl = [](int a, int n) { return a < 0 ? 0 : (a >= n ? n - 1 : a); };
    for (int z = 0; z < nz; z++)
        for (int x = 0; x < nx; x++)
            for (int y = 0; y < ny; y++) {
                double acc = 0.0;
                for (int dy = -1; dy <= 1; dy++)
                    for (int dx = -1; dx <= 1; dx++)
                        for (int dz = -1; dz <= 1; dz++)
                            acc += K[((1 - dy) * 3 + (1 - dx)) * 3 + (1 - dz)] *
                                   W.at(cl(y + dy, ny), cl(x + dx, nx), cl(z + dz, nz));
                L[y + ny * (x + nx * z)] = -acc;
            }
    std::vector<int> mx;
    std::vector<double> mv;
    for (int q = 0; q < N; q++) {
        const double lq = std::isnan(L[q]) ? -1e8 : L[q];
        if (!(lq > -1e8) || lq == 0.0) continue;  // MATLAB finds the nonzero maxima
        const int y = q % ny, x = (q / ny) % nx, z = q / (ny * nx);
        bool ok = true;
        for (int dz = -2; dz <= 2 && ok; dz++)
            for (int dx = -2; dx <= 2 && ok; dx++)
                for (int dy = -2; dy <= 2 && ok; dy++) {
                    const int yy = y + dy, xx = x + dx, zz = z + dz;
                    if (yy < 0 || yy >= ny || xx < 0 || xx >= nx || zz < 0 || zz >= nz) continue;
                    const double o = L[yy + ny * (xx + nx * zz)];
                    if ((std::isnan(o) ? -1e8 : o) > lq) ok = false;
                }
        if (ok) {
            mx.push_back(q);
            mv.push_back(lq);
        }
    }
    std::vector<int> ord(mx.size());
    for (size_t i = 0; i < ord.size(); i++) ord[i] = (int)i;
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return mv[a] > mv[b]; });
    std::vector<double> vals;
    for (int q = 0; q < N; q++)
        if (!std::isnan(W.v[q])) vals.push_back(W.v[q]);
    std::sort(vals.begin(), vals.end());
    const int kth = (int)std::floor(0.8 * (double)vals.size()) - 1;
    const double T = vals[std::max(kth, 0)];
    std::vector<char> bright(N, 0);
    for (int q = 0; q < N; q++) bright[q] = !std::isnan(W.v[q]) && W.v[q] > T;
    std::vector<int> sel, below;
    for (int i : ord) {
        if (bright[mx[i]])
            sel.push_back(mx[i]);
        else if (mv[i] > 0)
            below.push_back(mx[i]);
    }
    for (int q : below) {
        if (sel.size() >= 3) break;
        sel.push_back(q);
    }
    std::vector<std::array<double, 3>> out;
    for (int q : sel) {
        const double y = q % ny, x = (q / ny) % nx, z = q / (ny * nx);
        if (x == p0[0] && y == p0[1] && z == p0[2]) continue;
        out.push_back({x, y, z});
    }
    std::vector<int> lab(N, 0);
    std::vector<std::vector<int>> comps;
    for (int q = 0; q < N; q++) {
        if (!bright[q] || lab[q]) continue;
        comps.emplace_back();
        const int id = (int)comps.size();
        std::vector<int> stack{q};
        lab[q] = id;
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            comps.back().push_back(u);
            const int y = u % ny, x = (u / ny) % nx, z = u / (ny * nx);
            for (int dz = -1; dz <= 1; dz++)
                for (int dx = -1; dx <= 1; dx++)
                    for (int dy = -1; dy <= 1; dy++) {
                        const int nb = std::abs(dx) + std::abs(dy) + std::abs(dz);
                        if (nb == 0 || nb == 3) continue;  // 18-connected
                        const int yy = y + dy, xx = x + dx, zz = z + dz;
                        if (yy < 0 || yy >= ny || xx < 0 || xx >= nx || zz < 0 || zz >= nz) continue;
                        const int v = yy + ny * (xx + nx * zz);
                        if (bright[v] && !lab[v]) {
                            lab[v] = id;
                            stack.push_back(v);
                        }
                    }
        }
    }
    std::vector<int> cord;
    for (size_t j = 0; j < comps.size(); j++)
        if (comps[j].size() >= 2) cord.push_back((int)j);
    std::stable_sort(cord.begin(), cord.end(),
                     [&](int a, int b) { return comps[a].size() > comps[b].size(); });
    for (int j : cord) {
        std::vector<int> &cv = comps[j];
        std::sort(cv.begin(), cv.end());  // MATLAB's column-major order
        const int m = (int)cv.size();
        const double sw = msum(m, [&](int i) { return W.v[cv[i]]; });
        const double cx = msum(m, [&](int i) { return W.v[cv[i]] * (double)((cv[i] / ny) % nx); });
        const double cy = msum(m, [&](int i) { return W.v[cv[i]] * (double)(cv[i] % ny); });
        const double cz = msum(m, [&](int i) { return W.v[cv[i]] * (double)(cv[i] / (ny * nx)); });
        out.push_back({cx / sw, cy / sw, cz / sw});
    }
    if (out.size() > 16) out.resize(16);
    return out;
}

struct FitResult {
    State st;
    double rss_rep;  // the RSS the result reports (the best record's)
    int iters;
};

// the fit from p0: Newton steps on the position cut to the box, an Armijo line search (1/4 of
// the predicted decrease, at most 27 halvings), A and c by least squares; a switch to the next
// start when A < 0 or A exceeds the data's range, or when the box cut the step to nothing; the
// best record (A > 0 and a lower RSS, from the start's state) is the result; stops: the RSS drop
// below 1e-7 of the first RSS with steps below 1e-4 (held on the box, the direction it was held
// from), a drop below 1e-10, 31 iterations without a better record, max_iter iterations
FitResult iterate(const Window &w, const FullWindow &W, State st, double s, double sz, int max_iter,
                  std::vector<double> &g, std::vector<double> &gt) {
    const int n = w.n;
    double dmin = w.d[0], dmax = w.d[0];
    for (int i = 1; i < n; i++) dmin = std::min(dmin, w.d[i]), dmax = std::max(dmax, w.d[i]);
    const double range = dmax - dmin;
    const double p0[3] = {st.x[0], st.x[1], st.x[2]};
    model(w, st.x, s, sz, g);
    double rss = rss_of(w, g, st.A, st.c);
    const double rss0 = rss;
    double best_rss = rss;
    bool have_best = st.A > 0;
    State best = st;
    std::vector<std::array<double, 3>> ms;
    bool ms_done = false;
    size_t si = 0;
    int it_best = 0, it = 0;
    for (; it < max_iter; it++) {
        model(w, st.x, s, sz, g);
        rss = rss_of(w, g, st.A, st.c);
        double G[3], H[9], dv[3];
        grad_hess(w, st.x, s, sz, g, st.A, st.c, G, H);
        direction(G, H, dv);
        double alpha = 1.0;
        for (int j = 0; j < 3; j++) {
            if (st.x[j] + dv[j] > w.hi[j])
                alpha = std::min(alpha, (w.hi[j] - st.x[j]) / dv[j]);
            else if (st.x[j] + dv[j] < w.lo[j])
                alpha = std::min(alpha, (w.lo[j] - st.x[j]) / dv[j]);
        }
        const double step[3] = {alpha * dv[0], alpha * dv[1], alpha * dv[2]};
        const bool abad = st.A < 0 || st.A > range;
        bool switched = false;
        const double pstart[3] = {st.x[0], st.x[1], st.x[2]};
        double pn[3];
        if (abad || alpha == 0) {
            if (!ms_done) {
                ms = starts(W, p0, s, sz);
                ms_done = true;
            }
            if (si < ms.size()) {
                pn[0] = ms[si][0], pn[1] = ms[si][1], pn[2] = ms[si][2];
                si++;
                switched = true;
            }
        }
        if (switched) {
            model(w, pn, s, sz, gt);
        } else {
            const double gd = (G[0] * step[0] + G[1] * step[1]) + G[2] * step[2];
            double lam = 1.0, trial[3];
            for (int nh = 0;; nh++) {
                for (int j = 0; j < 3; j++) trial[j] = st.x[j] + lam * step[j];
                model(w, trial, s, sz, gt);
                const double rt = rss_of(w, gt, st.A, st.c);
                if (rt <= rss + 0.25 * lam * gd || nh >= 27) break;
                lam /= 2;
            }
            // positions within 1e-12 of the box are on it (where MATLAB's directions land; MATLAB
            // lands an ulp inside at times, and takes one more step there: not reproduced. On
            // 23,969 traced fits, within 1e-12 of MATLAB: 23,510 so, 23,376 with positions put
            // on the box only from beyond it, 23,273 without either; task C02d)
            bool snapped = false;
            for (int j = 0; j < 3; j++) {
                pn[j] = trial[j];
                if (std::fabs(trial[j] - w.lo[j]) <= 1e-12)
                    pn[j] = w.lo[j], snapped = snapped || trial[j] != w.lo[j];
                else if (std::fabs(trial[j] - w.hi[j]) <= 1e-12)
                    pn[j] = w.hi[j], snapped = snapped || trial[j] != w.hi[j];
            }
            if (snapped) model(w, pn, s, sz, gt);
        }
        const double A0 = st.A;
        const double num = msum(n, [&](int i) { return ((w.d[i] - A0 * gt[i]) - st.c) * gt[i]; });
        const double den = msum(n, [&](int i) { return gt[i] * gt[i]; });
        st.A = st.A + num / den;
        const double A1 = st.A;
        st.c = st.c + msum(n, [&](int i) { return (w.d[i] - A1 * gt[i]) - st.c; }) / n;
        const double rss_new = rss_of(w, gt, st.A, st.c);
        // the step of the stop test: the position's move, or, held on the box (the step cut to
        // nothing), the direction it was held from (a fit on the box goes on while A and c
        // still move; task C02d)
        double dp = 0.0;
        for (int j = 0; j < 3; j++)
            dp = std::max(dp, std::fabs(!switched && alpha == 0 ? dv[j] : pn[j] - pstart[j]));
        std::memcpy(st.x, pn, sizeof(pn));
        g.swap(gt);
        if (st.A > 0 && rss_new < best_rss) {
            best_rss = rss_new;
            best = st;
            have_best = true;
            it_best = it;
        }
        if (it - it_best >= 31) {
            it++;
            break;
        }
        const double drop = (rss - rss_new) / rss0;
        if ((drop < 1e-7 && dp < 1e-4) || std::fabs(drop) < 1e-10) {
            it++;
            break;
        }
    }
    return {have_best ? best : st, best_rss, it};
}

// the Anderson-Darling statistic of x in its given order against N(mu, sigma^2)
double ad_stat(const double *x, int n, double mu, double sigma) {
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        const double zi = (x[i] - mu) / sigma, zj = (x[n - 1 - i] - mu) / sigma;
        const double Fi = 0.5 * pk5d_ad::erfc(-zi / std::sqrt(2.0));
        const double Fj = 0.5 * pk5d_ad::erfc(-zj / std::sqrt(2.0));
        s += (2.0 * (i + 1) - 1.0) * (pk5d_ad::log(Fi) + pk5d_ad::log(1.0 - Fj));
    }
    return -n - s / n;
}

// the outputs at the result: the Hessian in (x, y, z, A, c) and its inverse, the parameters'
// standard deviations from the reported RSS, the residuals' statistics at the result
void finish(const Window &w, const State &st, double rss_rep, double s, double sz, double *prm, double *pstd,
            double *C, double *stats, std::vector<double> &g, double ad_crit) {
    model(w, st.x, s, sz, g);
    const int n = w.n;
    std::vector<double> r(n);
    for (int i = 0; i < n; i++) r[i] = (w.d[i] - st.A * g[i]) - st.c;
    double G[3], H[9];
    grad_hess(w, st.x, s, sz, g, st.A, st.c, G, H);
    const double s2 = s * s, sz2 = sz * sz;
    double H5[25];
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++) H5[a * 5 + b] = H[a * 3 + b];
    H5[3 * 5 + 3] = 2 * msum(n, [&](int i) { return g[i] * g[i]; });
    H5[4 * 5 + 4] = 2.0 * n;
    H5[3 * 5 + 4] = H5[4 * 5 + 3] = 2 * msum(n, [&](int i) { return g[i]; });
    const double sc[3] = {s2, s2, sz2};
    const double *P[3] = {w.X, w.Y, w.Z};
    for (int a = 0; a < 3; a++) {
        auto dk = [&](int i) { return (P[a][i] - st.x[a]) / sc[a]; };
        // the position-amplitude entries as MATLAB's result has them: the Gauss-Newton part
        // without its factor A, minus the residual part
        H5[a * 5 + 3] = H5[3 * 5 + a] = 2 * msum(n, [&](int i) { return g[i] * g[i] * dk(i); }) -
                                        2 * msum(n, [&](int i) { return r[i] * g[i] * dk(i); });
        H5[a * 5 + 4] = H5[4 * 5 + a] = 2 * msum(n, [&](int i) { return st.A * g[i] * dk(i); });
    }
    inverse(5, H5, C);
    for (int k = 0; k < 7; k++) pstd[k] = 0.0;
    const int idx[5] = {0, 1, 2, 3, 6};
    // a negative diagonal entry (C of a fit that ends where the RSS is not convex, as on the
    // box at a border plane) gives the magnitude of its complex root, as MATLAB's results
    // have it (task C02d: the pstd of such fits, from MATLAB's own parameters, within 1e-9)
    for (int k = 0; k < 5; k++) {
        const double v = C[k * 5 + k] * rss_rep / (n - 6);
        pstd[idx[k]] = std::sqrt(std::fabs(v));
    }
    prm[0] = st.x[0], prm[1] = st.x[1], prm[2] = st.x[2], prm[3] = st.A, prm[4] = s, prm[5] = sz,
    prm[6] = st.c;
    const double sd = std::sqrt(msum(n, [&](int i) { return r[i] * r[i]; }) / (n - 1));
    stats[0] = rss_rep;
    stats[1] = msum(n, [&](int i) { return r[i]; }) / n;
    stats[2] = sd;
    std::sort(r.begin(), r.end());
    const double a2 = ad_stat(r.data(), n, 0.0, sd);
    stats[3] = a2 > ad_crit ? 1.0 : 0.0;
    stats[4] = a2;
}

}  // namespace

// -- bindings ------------------------------------------------------------------------------------

using f64_in = nb::ndarray<const double, nb::c_contig, nb::device::cpu>;
using f64_out = nb::ndarray<double, nb::c_contig, nb::device::cpu>;
using i64_in = nb::ndarray<const int64_t, nb::c_contig, nb::device::cpu>;

// windows: (k, m) values in MATLAB's column-major window order (NaN: excluded), dims (k, 3) the
// windows' [ny nx nz], p0 (k, 7) [x y z A s sz c] in window coordinates (0-based);
// outputs prm (k, 7), pstd (k, 7), C (k, 25), stats (k, 6): RSS, mean, std, hAD, A2, iterations
void bind_fit_cd(f64_in vals, i64_in dims, f64_in p0, f64_out prm, f64_out pstd, f64_out C, f64_out stats,
                 int max_iter, double ad_crit, int num_threads) {
    pk5d::check_threads(num_threads);
    if (vals.ndim() != 2 || dims.ndim() != 2 || p0.ndim() != 2) throw std::invalid_argument("2-D inputs");
    const i64 k = (i64)vals.shape(0), m = (i64)vals.shape(1);
    if ((i64)dims.shape(0) != k || dims.shape(1) != 3) throw std::invalid_argument("dims must be (k, 3)");
    if ((i64)p0.shape(0) != k || p0.shape(1) != 7) throw std::invalid_argument("p0 must be (k, 7)");
    if ((i64)prm.size() != 7 * k || (i64)pstd.size() != 7 * k || (i64)C.size() != 25 * k ||
        (i64)stats.size() != 6 * k)
        throw std::invalid_argument("output sizes");
    for (i64 i = 0; i < k; i++) {
        const int64_t *dd = dims.data() + 3 * i;
        if (dd[0] < 1 || dd[1] < 1 || dd[2] < 1 || dd[0] * dd[1] * dd[2] > m)
            throw std::invalid_argument("a window's dims do not fit its row");
    }
    const double *V = vals.data(), *P0 = p0.data();
    const int64_t *D = dims.data();
    double *PR = prm.data(), *PS = pstd.data(), *CC = C.data(), *ST = stats.data();
    {
        nb::gil_scoped_release release;
#pragma omp parallel for schedule(dynamic, 16) num_threads(num_threads)
        for (i64 i = 0; i < k; i++) {
            const i64 ny = D[3 * i], nx = D[3 * i + 1], nz = D[3 * i + 2];
            const i64 nv = ny * nx * nz;
            std::vector<double> d, X, Y, Z;
            d.reserve(nv), X.reserve(nv), Y.reserve(nv), Z.reserve(nv);
            for (i64 q = 0; q < nv; q++) {
                const double v = V[i * m + q];
                if (std::isnan(v)) continue;
                d.push_back(v);
                Y.push_back((double)(q % ny));
                X.push_back((double)((q / ny) % nx));
                Z.push_back((double)(q / (ny * nx)));
            }
            // the box: from the window's first voxel to its last valid one per axis (a window
            // whose last planes are all NaN, as at the deskewed frame's empty wedges, keeps its
            // fits off them, its first planes not; task C02d, MATLAB's traced fits)
            Window w{(int)d.size(), d.data(), X.data(), Y.data(), Z.data(), {0, 0, 0}, {0, 0, 0}};
            if (!d.empty()) {
                const double *P3[3] = {X.data(), Y.data(), Z.data()};
                for (int a = 0; a < 3; a++) w.hi[a] = *std::max_element(P3[a], P3[a] + d.size());
            }
            const double *p = P0 + 7 * i;
            double *pr = PR + 7 * i, *ps = PS + 7 * i, *c5 = CC + 25 * i, *st = ST + 6 * i;
            if (w.n < 7) {
                for (int j = 0; j < 7; j++) pr[j] = ps[j] = std::numeric_limits<double>::quiet_NaN();
                for (int j = 0; j < 25; j++) c5[j] = std::numeric_limits<double>::quiet_NaN();
                for (int j = 0; j < 6; j++) st[j] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }
            std::vector<double> g(w.n), gt(w.n);
            State s0{{p[0], p[1], p[2]}, p[3], p[6]};
            const FullWindow W{V + i * m, (int)ny, (int)nx, (int)nz};
            FitResult fr = iterate(w, W, s0, p[4], p[5], max_iter, g, gt);
            finish(w, fr.st, fr.rss_rep, p[4], p[5], pr, ps, c5, st, g, ad_crit);
            st[5] = fr.iters;
        }
    }
}

double bind_ad_stat(f64_in x, double mu, double sigma) {
    return ad_stat(x.data(), (int)x.size(), mu, sigma);
}

// the fitGaussian3D mex (mexfit.hpp) on k windows: vals (k, m) in MATLAB's column-major window
// order (NaN: left out), dims (k, 3) [ny nx nz], p0 (k, 7) [x y z A s sz c], free (np,) the
// free parameters' indices in that order; outputs prm (k, 7), pstd (k, 7) (NaN for fixed
// ones), C (k, np, np), stats (k, 8): RSS, mean, std, hAD, A2, iterations, evaluations, status;
// res (k, m) the residuals (model - data, NaN where left out) when res.size() > 0
void bind_fit_mex(f64_in vals, i64_in dims, f64_in p0, i64_in free, f64_out prm, f64_out pstd, f64_out C,
                  f64_out stats, f64_out res, double epsabs, double epsrel, int maxiter,
                  double ad_crit, int num_threads) {
    pk5d::check_threads(num_threads);
    if (vals.ndim() != 2 || dims.ndim() != 2 || p0.ndim() != 2) throw std::invalid_argument("2-D inputs");
    const i64 k = (i64)vals.shape(0), m = (i64)vals.shape(1), np = (i64)free.size();
    if ((i64)dims.shape(0) != k || dims.shape(1) != 3) throw std::invalid_argument("dims must be (k, 3)");
    if ((i64)p0.shape(0) != k || p0.shape(1) != 7) throw std::invalid_argument("p0 must be (k, 7)");
    if (np < 1 || np > 7) throw std::invalid_argument("1 to 7 free parameters");
    if ((i64)prm.size() != 7 * k || (i64)pstd.size() != 7 * k || (i64)C.size() != np * np * k ||
        (i64)stats.size() != 8 * k || (res.size() != 0 && (i64)res.size() != k * m))
        throw std::invalid_argument("output sizes");
    std::vector<int> fr(np);
    for (i64 j = 0; j < np; j++) {
        fr[j] = (int)free.data()[j];
        if (fr[j] < 0 || fr[j] > 6 || (j > 0 && fr[j] <= fr[j - 1]))
            throw std::invalid_argument("free parameters: increasing indices 0 to 6");
    }
    const i64 *D = dims.data();
    for (i64 i = 0; i < k; i++) {
        if (D[3 * i] < 1 || D[3 * i + 1] < 1 || D[3 * i + 2] < 1 || D[3 * i] * D[3 * i + 1] * D[3 * i + 2] > m)
            throw std::invalid_argument("a window's dims do not fit its row");
    }
    const double *V = vals.data(), *P0 = p0.data();
    double *PR = prm.data(), *PS = pstd.data(), *CC = C.data(), *ST = stats.data();
    double *RS = res.size() ? res.data() : nullptr;
    mexfit::Options o;
    o.epsabs = epsabs, o.epsrel = epsrel, o.maxiter = maxiter;
    {
        nb::gil_scoped_release release;
#pragma omp parallel for schedule(dynamic, 8) num_threads(num_threads)
        for (i64 i = 0; i < k; i++) {
            mexfit::Win w;
            w.ny = (int)D[3 * i], w.nx = (int)D[3 * i + 1], w.nz = (int)D[3 * i + 2];
            const i64 nv = (i64)w.ny * w.nx * w.nz;
            const double *v = V + i * m;
            for (i64 q = 0; q < nv; q++) {
                if (std::isnan(v[q])) continue;
                w.d.push_back(v[q]);
                w.yi.push_back((int)(q % w.ny));
                w.xi.push_back((int)((q / w.ny) % w.nx));
                w.zi.push_back((int)(q / ((i64)w.ny * w.nx)));
                w.pos.push_back((int)q);
            }
            w.n = (int)w.d.size();
            double *st = ST + 8 * i;
            const mexfit::Result r = mexfit::fit(w, P0 + 7 * i, fr.data(), (int)np, o);
            const int n = w.n;
            double rss = 0.0, sm = 0.0;
            for (int q = 0; q < n; q++) rss += r.f[q] * r.f[q];
            for (int q = 0; q < n; q++) sm += r.f[q];
            const double mean = sm / n;
            double ss = 0.0;
            for (int q = 0; q < n; q++) ss += (r.f[q] - mean) * (r.f[q] - mean);
            const double sd = std::sqrt(ss / (n - 1));
            for (int j = 0; j < 7; j++) {
                PR[7 * i + j] = r.prm[j];
                PS[7 * i + j] = std::numeric_limits<double>::quiet_NaN();
            }
            const double scale = rss / (double)(n - np - 1);
            for (i64 j = 0; j < np; j++) PS[7 * i + fr[j]] = std::sqrt(r.C[j * np + j] * scale);
            for (i64 q = 0; q < np * np; q++) CC[i * np * np + q] = r.C[q];
            // the Anderson-Darling test of the sorted residuals over their standard deviation
            std::vector<double> sr(r.f);
            std::sort(sr.begin(), sr.end());
            const double den = sd * std::sqrt(2.0);
            std::vector<double> F(n);
            for (int q = 0; q < n; q++) F[q] = 0.5 * (1.0 + pk5d_ad::erf(sr[q] / den));
            double S = 0.0;
            for (int q = 0; q < n; q++)
                S += (2.0 * (q + 1) - 1.0) * (pk5d_ad::log(F[q]) + pk5d_ad::log(1.0 - F[n - 1 - q]));
            const double a2 = -n - S / n;
            st[0] = rss, st[1] = mean, st[2] = sd, st[3] = a2 > ad_crit ? 1.0 : 0.0, st[4] = a2;
            st[5] = r.iterations, st[6] = r.nfev, st[7] = r.status;
            if (RS) {
                double *o_ = RS + i * m;
                for (i64 q = 0; q < m; q++) o_[q] = std::numeric_limits<double>::quiet_NaN();
                for (int q = 0; q < n; q++) o_[w.pos[q]] = r.f[q];
            }
        }
    }
}

double bind_cr_exp(double x) { return mexfit::cr_exp(x); }

// the lowSNR and 'speedup' detectors' voxel-wise model statistics in one threaded pass, each
// voxel's arithmetic in the order of the expressions in detect.py (so the results are those of
// the NumPy expressions): A, c, the residuals' sigma, the t statistic and its degrees of freedom,
// from the window sums fg (of the Gaussian times the volume), fu (of the volume), fu2 (of its
// square). gfac multiplies fg first ('speedup': the normalized kernel's sum; 1 for lowSNR, which
// skips the product). df's fourth powers are (x x)(x x).
void bind_prefit_stats(f64_in fg, f64_in fu, f64_in fu2, double gsum, double g2sum, double n, double C11,
                       double k_level, double sq, double gfac, f64_out A, f64_out c, f64_out sres, f64_out T,
                       f64_out df, int num_threads) {
    pk5d::check_threads(num_threads);
    const i64 N = (i64)fg.size();
    if ((i64)fu.size() != N || (i64)fu2.size() != N || (i64)A.size() != N || (i64)c.size() != N ||
        (i64)sres.size() != N || (i64)T.size() != N || (i64)df.size() != N)
        throw std::invalid_argument("sizes");
    const double *FG = fg.data(), *FU = fu.data(), *FU2 = fu2.data();
    double *PA = A.data(), *PC = c.data(), *PS = sres.data(), *PT = T.data(), *PD = df.data();
    const bool scaled = gfac != 1.0;
    const double denA = g2sum - gsum * gsum / n;
    nb::gil_scoped_release release;
#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (i64 i = 0; i < N; i++) {
        const double fgi = scaled ? FG[i] * gfac : FG[i];
        // A = (fg - gsum fu / n) / (g2sum - gsum^2 / n) ('speedup': fg gsum - gsum fu / n)
        const double a = scaled ? (FG[i] * gfac - (gsum * FU[i]) / n) / denA : (FG[i] - (gsum * FU[i]) / n) / denA;
        const double cc = (FU[i] - a * gsum) / n;
        // RSS = A A g2sum - 2 A (fg - c gsum) + fu2 + n c c - 2 c fu, clipped at 0
        double t = ((a * a) * g2sum) - ((2.0 * a) * (fgi - cc * gsum));
        t = t + FU2[i];
        t = t + ((n * cc) * cc);
        t = t - ((2.0 * cc) * FU[i]);
        const double rss = t < 0 ? 0.0 : t;
        const double sr = std::sqrt(rss / (n - 1));
        const double sA = std::sqrt((rss / (n - 3)) * C11);
        const double se = (sr / sq) * k_level;
        const double u = sA * sA + se * se;
        const double sc = std::sqrt(u / n);
        PA[i] = a;
        PC[i] = cc;
        PS[i] = sr;
        PT[i] = (a - sr * k_level) / sc;
        const double sA2 = sA * sA, se2 = se * se;
        PD[i] = ((n - 1) * (u * u)) / (sA2 * sA2 + se2 * se2);
    }
}

// XR_locmax3d as the detectors call it (detect.locmax3d): with NaN as -1e8, a voxel is a
// maximum when every other voxel of its window (clipped at the faces; half sizes hz, hy, hx) is
// strictly smaller; out = the replaced value times 1 at maxima and 0 elsewhere (so -0 at NaN and
// negative non-maxima); (nz, ny, nx)
void bind_locmax3d(f64_in vol, int hz, int hy, int hx, f64_out out, int num_threads) {
    pk5d::check_threads(num_threads);
    if (vol.ndim() != 3) throw std::invalid_argument("vol must be (nz, ny, nx)");
    const i64 nz = (i64)vol.shape(0), ny = (i64)vol.shape(1), nx = (i64)vol.shape(2);
    if ((i64)out.size() != nz * ny * nx) throw std::invalid_argument("out size");
    const double *V = vol.data();
    double *O = out.data();
    auto val = [&](i64 i) {
        const double v = V[i];
        return std::isnan(v) ? -1e8 : v;
    };
    nb::gil_scoped_release release;
#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (i64 row = 0; row < nz * ny; row++) {
        const i64 z = row / ny, y = row % ny;
        const i64 z0 = std::max<i64>(z - hz, 0), z1 = std::min<i64>(z + hz, nz - 1);
        const i64 y0 = std::max<i64>(y - hy, 0), y1 = std::min<i64>(y + hy, ny - 1);
        for (i64 x = 0; x < nx; x++) {
            const i64 i = row * nx + x;
            const double v = val(i);
            const i64 x0 = std::max<i64>(x - hx, 0), x1 = std::min<i64>(x + hx, nx - 1);
            bool is_max = true;
            for (i64 zz = z0; zz <= z1 && is_max; zz++)
                for (i64 yy = y0; yy <= y1 && is_max; yy++) {
                    const i64 base = (zz * ny + yy) * nx;
                    for (i64 xx = x0; xx <= x1; xx++) {
                        if (zz == z && yy == y && xx == x) continue;
                        if (!(val(base + xx) < v)) {
                            is_max = false;
                            break;
                        }
                    }
                }
            O[i] = v * (is_max ? 1.0 : 0.0);
        }
    }
}

// the lowSNR detector's mask: T > -t or the chi-square test (fu2 - 2 c fu + n c^2) / sres^2 > chi,
// without the voxels within border of a face; (nz, ny, nx) arrays
void bind_lowsnr_mask(f64_in T, f64_in c, f64_in fu, f64_in fu2, f64_in sres, double n, double tthr, double chi,
                      int border, nb::ndarray<uint8_t, nb::c_contig, nb::device::cpu> mask, int num_threads) {
    pk5d::check_threads(num_threads);
    if (T.ndim() != 3) throw std::invalid_argument("(nz, ny, nx) arrays");
    const i64 nz = (i64)T.shape(0), ny = (i64)T.shape(1), nx = (i64)T.shape(2);
    const double *PT = T.data(), *PC = c.data(), *FU = fu.data(), *FU2 = fu2.data(), *PS = sres.data();
    uint8_t *M = mask.data();
    nb::gil_scoped_release release;
#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (i64 row = 0; row < nz * ny; row++) {
        const i64 z = row / ny, y = row % ny;
        const bool edge = z < border || z >= nz - border || y < border || y >= ny - border;
        for (i64 x = 0; x < nx; x++) {
            const i64 i = row * nx + x;
            bool m = false;
            if (!edge && x >= border && x < nx - border) {
                const double cc = PC[i];
                const double S = ((FU2[i] - (2.0 * cc) * FU[i]) + (n * cc) * cc);
                m = (PT[i] > -tthr) || (S / (PS[i] * PS[i]) > chi);
            }
            M[i] = m ? 1 : 0;
        }
    }
}


// conv3fast's separable filter (estGaussianAmplitude3D's), measured bit for bit: per axis the
// symmetric kernel given by its half h[0..w], out(i) = h0 v(i) + h1 (v(i-1) + v(i+1)) + ...
// summed from the centre out, the borders mirrored without repeating the edge voxel
// (v(-k) = v(k)); the passes along x, then y, then z. vol and out are (nz, ny, nx); work is
// a buffer of the same size.
inline i64 mirror(i64 i, i64 n) {
    if (n == 1) return 0;
    const i64 p = 2 * (n - 1);
    i64 j = i % p;
    if (j < 0) j += p;
    return j < n ? j : p - j;
}

void bind_conv3fast(f64_in vol, f64_in hy, f64_in hx, f64_in hz, f64_out out, f64_out work, int num_threads) {
    pk5d::check_threads(num_threads);
    if (vol.ndim() != 3) throw std::invalid_argument("vol must be (nz, ny, nx)");
    const i64 nz = (i64)vol.shape(0), ny = (i64)vol.shape(1), nx = (i64)vol.shape(2);
    if ((i64)out.size() != nz * ny * nx || (i64)work.size() != nz * ny * nx)
        throw std::invalid_argument("out and work must be the volume's size");
    const i64 wy = (i64)hy.size() - 1, wx = (i64)hx.size() - 1, wz = (i64)hz.size() - 1;
    if (wy < 0 || wx < 0 || wz < 0) throw std::invalid_argument("empty kernel");
    const double *V = vol.data(), *Hy = hy.data(), *Hx = hx.data(), *Hz = hz.data();
    double *O = out.data(), *T = work.data();
    const i64 plane = ny * nx;
    nb::gil_scoped_release release;
    // x: vol -> out
#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (i64 row = 0; row < nz * ny; row++) {
        const double *v = V + row * nx;
        double *o = O + row * nx;
        for (i64 x = 0; x < nx; x++) {
            double acc = Hx[0] * v[x];
            for (i64 k = 1; k <= wx; k++) acc = acc + Hx[k] * (v[mirror(x - k, nx)] + v[mirror(x + k, nx)]);
            o[x] = acc;
        }
    }
    // y: out -> work
#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (i64 row = 0; row < nz * ny; row++) {
        const i64 z = row / ny, y = row % ny;
        const double *c = O + z * plane + y * nx;
        double *o = T + z * plane + y * nx;
        for (i64 x = 0; x < nx; x++) o[x] = Hy[0] * c[x];
        for (i64 k = 1; k <= wy; k++) {
            const double *a = O + z * plane + mirror(y - k, ny) * nx, *b = O + z * plane + mirror(y + k, ny) * nx;
            const double h = Hy[k];
            for (i64 x = 0; x < nx; x++) o[x] = o[x] + h * (a[x] + b[x]);
        }
    }
    // z: work -> out
#pragma omp parallel for schedule(static) num_threads(num_threads)
    for (i64 row = 0; row < nz * ny; row++) {
        const i64 z = row / ny, y = row % ny;
        const double *c = T + z * plane + y * nx;
        double *o = O + z * plane + y * nx;
        for (i64 x = 0; x < nx; x++) o[x] = Hz[0] * c[x];
        for (i64 k = 1; k <= wz; k++) {
            const double *a = T + mirror(z - k, nz) * plane + y * nx, *b = T + mirror(z + k, nz) * plane + y * nx;
            const double h = Hz[k];
            for (i64 x = 0; x < nx; x++) o[x] = o[x] + h * (a[x] + b[x]);
        }
    }
}

double bind_exp(double x) { return fdlibm_exp(x); }

NB_MODULE(PK5D_MODULE_NAME, m) {
    m.doc() = "the Gaussian fits of the 3-D point-source detection (CD_Newton) and their tests";
    m.def("fit_cd", &bind_fit_cd, nb::arg("vals").noconvert(), nb::arg("dims").noconvert(),
          nb::arg("p0").noconvert(), nb::arg("prm").noconvert(), nb::arg("pstd").noconvert(),
          nb::arg("C").noconvert(), nb::arg("stats").noconvert(), nb::arg("max_iter"), nb::arg("ad_crit"),
          nb::arg("num_threads"));
    m.def("ad_stat", &bind_ad_stat, nb::arg("x").noconvert(), nb::arg("mu"), nb::arg("sigma"));
    m.def("exp", &bind_exp, nb::arg("x"));
    m.def("cr_exp", &bind_cr_exp, nb::arg("x"));
    m.def("prefit_stats", &bind_prefit_stats, nb::arg("fg").noconvert(), nb::arg("fu").noconvert(),
          nb::arg("fu2").noconvert(), nb::arg("gsum"), nb::arg("g2sum"), nb::arg("n"), nb::arg("C11"),
          nb::arg("k_level"), nb::arg("sq"), nb::arg("gfac"), nb::arg("A").noconvert(), nb::arg("c").noconvert(),
          nb::arg("sres").noconvert(), nb::arg("T").noconvert(), nb::arg("df").noconvert(), nb::arg("num_threads"));
    m.def("locmax3d", &bind_locmax3d, nb::arg("vol").noconvert(), nb::arg("hz"), nb::arg("hy"), nb::arg("hx"),
          nb::arg("out").noconvert(), nb::arg("num_threads"));
    m.def("lowsnr_mask", &bind_lowsnr_mask, nb::arg("T").noconvert(), nb::arg("c").noconvert(),
          nb::arg("fu").noconvert(), nb::arg("fu2").noconvert(), nb::arg("sres").noconvert(), nb::arg("n"),
          nb::arg("tthr"), nb::arg("chi"), nb::arg("border"), nb::arg("mask").noconvert(), nb::arg("num_threads"));
    m.def("conv3fast", &bind_conv3fast, nb::arg("vol").noconvert(), nb::arg("hy").noconvert(),
          nb::arg("hx").noconvert(), nb::arg("hz").noconvert(), nb::arg("out").noconvert(),
          nb::arg("work").noconvert(), nb::arg("num_threads"));
    m.def("fit_mex", &bind_fit_mex, nb::arg("vals").noconvert(), nb::arg("dims").noconvert(),
          nb::arg("p0").noconvert(), nb::arg("free").noconvert(), nb::arg("prm").noconvert(),
          nb::arg("pstd").noconvert(), nb::arg("C").noconvert(), nb::arg("stats").noconvert(),
          nb::arg("res").noconvert(), nb::arg("epsabs"), nb::arg("epsrel"), nb::arg("maxiter"),
          nb::arg("ad_crit"), nb::arg("num_threads"));
}
