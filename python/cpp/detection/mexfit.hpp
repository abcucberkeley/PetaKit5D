// The fitGaussian3D mex's fits (task C02b), reimplemented from the method and from MATLAB
// R2024b's measured behaviour (black-box probes: inputs and outputs, and the mex's own calls of
// exp, erf and log observed with an interposed math library; no source of the mex, which has
// none in the repository, nor of any GPL library):
// - the model A exp(-(X-x)^2/(2 s^2)) exp(-(Y-y)^2/(2 s^2)) exp(-(Z-z)^2/(2 sz^2)) + c on the
//   window's valid voxels (0-based coordinates from the window's first voxel), evaluated per
//   axis and multiplied amplitude first, ((A ex) ey) ez + c - d (the measured residuals);
//   at most 500 iterations, an iteration ending after ten rejected trials;
// - the parameters of the mode's letters (x y z A s r c, r = sz) fitted by the Levenberg-
//   Marquardt method in Moré's trust-region form (J. J. Moré, "The Levenberg-Marquardt
//   algorithm: implementation and theory", 1978; the MINPACK-1 routines lmder, lmpar, qrfac,
//   qrsolv and enorm as published by Moré, Garbow and Hillstrom, 1980): the Jacobian at the
//   start and at every accepted point, column scaling from the Jacobian's norms, steps until
//   the change of every parameter is below 1e-8 + 1e-8 |p| or the step can no longer reduce the
//   residuals (the measured evaluation sequence);
// - the covariance inv(J'J) at the result, prmStd = sqrt(diag(C) RSS / (n - p - 1)), and the
//   residuals' RSS, mean, standard deviation and Anderson-Darling test (sorted residuals over
//   their standard deviation, mean 0).
// The exponential is correctly rounded (double-double), as the platform-independent stand-in
// for the C library's: MATLAB's mex calls glibc's, which is correctly rounded for 99.85 % of the
// arguments the probes recorded.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace mexfit {

// -- double-double arithmetic ------------------------------------------------------------------

struct dd {
    double hi, lo;
};

inline dd two_sum(double a, double b) {
    const double s = a + b, bb = s - a;
    return {s, (a - (s - bb)) + (b - bb)};
}
inline dd fast_two_sum(double a, double b) {  // |a| >= |b|
    const double s = a + b;
    return {s, b - (s - a)};
}
inline void split(double a, double &h, double &l) {  // Veltkamp's split (|a| < 2^996)
    const double c = 134217729.0 * a;  // 2^27 + 1
    h = c - (c - a);
    l = a - h;
}
inline dd two_prod(double a, double b) {  // Dekker's product, no fused multiply-add
    const double p = a * b;
    double ah, al, bh, bl;
    split(a, ah, al);
    split(b, bh, bl);
    return {p, ((ah * bh - p) + ah * bl + al * bh) + al * bl};
}
inline dd dd_add(dd a, dd b) {
    dd s = two_sum(a.hi, b.hi);
    dd t = two_sum(a.lo, b.lo);
    s.lo += t.hi;
    s = fast_two_sum(s.hi, s.lo);
    s.lo += t.lo;
    return fast_two_sum(s.hi, s.lo);
}
inline dd dd_mul(dd a, dd b) {
    dd p = two_prod(a.hi, b.hi);
    p.lo += a.hi * b.lo + a.lo * b.hi;
    return fast_two_sum(p.hi, p.lo);
}
inline dd dd_mul_d(dd a, double b) {
    dd p = two_prod(a.hi, b);
    p.lo += a.lo * b;
    return fast_two_sum(p.hi, p.lo);
}

// exp(x) correctly rounded (to nearest) but in cases closer to a rounding boundary than about
// 2^-100 relative: x = k ln2 + r, |r| <= ln2/2, r in double-double from a three-part ln2;
// exp(r / 256) by its Taylor series, squared eight times; scaled by 2^k
inline double cr_exp(double x) {
    if (std::isnan(x)) return x + x;
    if (x > 709.782712893384) return std::numeric_limits<double>::infinity();
    if (x < -745.1332191019412) return 0.0;
    if (x == 0.0) return 1.0;
    static const double INV_LN2 = 1.4426950408889634;
    // ln2 = L1 + L2 + L3 to about 2^-117, L1 and L2 of 32 significant bits (k L1, k L2 exact)
    static const double L1 = 0x1.62e42fee00000p-1, L2 = 0x1.a39ef35600000p-33, L3 = 0x1.93c7673007e5fp-65;
    const double k = std::nearbyint(x * INV_LN2);
    // r = x - k L1 - k L2 - k L3 (the first two exact)
    dd r = two_sum(x - k * L1, -(k * L2));
    r = dd_add(r, two_prod(-k, L3));
    // exp(t), t = r / 256: 1 + t + t^2/2! + ... + t^11/11! (|t| < 1.4e-3: the next term < 1e-40)
    const dd t = {std::ldexp(r.hi, -8), std::ldexp(r.lo, -8)};
    static const double INVF_HI[12] = {1.0,
                                       1.0,
                                       0.5,
                                       0.16666666666666666,
                                       0.041666666666666664,
                                       0.008333333333333333,
                                       0.001388888888888889,
                                       0.0001984126984126984,
                                       2.48015873015873e-05,
                                       2.7557319223985893e-06,
                                       2.755731922398589e-07,
                                       2.505210838544172e-08};
    static const double INVF_LO[12] = {0.0,
                                       0.0,
                                       0.0,
                                       9.25185853854297e-18,
                                       2.3129646346357427e-18,
                                       1.1564823173178714e-19,
                                       -5.300543954373577e-20,
                                       1.7209558293420705e-22,
                                       2.1511947866775882e-23,
                                       -1.858393274046472e-22,
                                       2.3767714622250297e-23,
                                       -1.448814070935912e-24};
    dd p = {INVF_HI[11], INVF_LO[11]};
    for (int i = 10; i >= 0; i--) p = dd_add(dd_mul(p, t), dd{INVF_HI[i], INVF_LO[i]});
    for (int i = 0; i < 8; i++) p = dd_mul(p, p);
    const int ki = (int)k;
    if (ki >= -1021) {
        const double y = p.hi + p.lo;  // the double nearest to the double-double
        return std::ldexp(y, ki);
    }
    // subnormal results: round once at the subnormal's precision
    const double y = std::ldexp(p.hi, ki) + std::ldexp(p.lo, ki);
    return y;
}

// -- the model and its Jacobian --------------------------------------------------------------

enum { PX = 0, PY, PZ, PA, PS, PR, PC };

struct Win {
    int ny, nx, nz, n;
    std::vector<double> d;        // valid values in MATLAB's column-major order
    std::vector<int> xi, yi, zi;  // their voxel indices
    std::vector<int> pos;         // their offsets in the window (column-major)
};

// per-axis exponentials at parameters p
inline void axes_exp(const Win &w, const double *p, std::vector<double> &ex, std::vector<double> &ey,
                     std::vector<double> &ez) {
    const double s2 = 2 * p[PS] * p[PS], r2 = 2 * p[PR] * p[PR];
    ex.resize(w.nx);
    ey.resize(w.ny);
    ez.resize(w.nz);
    for (int k = 0; k < w.nx; k++) {
        const double t = k - p[PX];
        ex[k] = cr_exp(-(t * t) / s2);
    }
    for (int k = 0; k < w.ny; k++) {
        const double t = k - p[PY];
        ey[k] = cr_exp(-(t * t) / s2);
    }
    for (int k = 0; k < w.nz; k++) {
        const double t = k - p[PZ];
        ez[k] = cr_exp(-(t * t) / r2);
    }
}

// residuals f = ((A ex) ey) ez + c - d
inline void residuals(const Win &w, const double *p, std::vector<double> &f, std::vector<double> &ex,
                      std::vector<double> &ey, std::vector<double> &ez) {
    axes_exp(w, p, ex, ey, ez);
    f.resize(w.n);
    const double A = p[PA], c = p[PC];
    for (int i = 0; i < w.n; i++) f[i] = ((A * ex[w.xi[i]]) * ey[w.yi[i]]) * ez[w.zi[i]] + c - w.d[i];
}

// Jacobian columns of the free parameters (column-major, n rows), at p
inline void jacobian(const Win &w, const double *p, const int *free, int np, std::vector<double> &J,
                     std::vector<double> &ex, std::vector<double> &ey, std::vector<double> &ez) {
    axes_exp(w, p, ex, ey, ez);
    J.resize((size_t)w.n * np);
    const double A = p[PA], s = p[PS], r = p[PR];
    const double s2 = s * s, r2 = r * r, s3 = s2 * s, r3 = r2 * r;
    for (int i = 0; i < w.n; i++) {
        const double g = (ex[w.xi[i]] * ey[w.yi[i]]) * ez[w.zi[i]];
        const double Ag = ((A * ex[w.xi[i]]) * ey[w.yi[i]]) * ez[w.zi[i]];
        const double dx = w.xi[i] - p[PX], dy = w.yi[i] - p[PY], dz = w.zi[i] - p[PZ];
        for (int j = 0; j < np; j++) {
            double v;
            switch (free[j]) {
                case PX: v = Ag * dx / s2; break;
                case PY: v = Ag * dy / s2; break;
                case PZ: v = Ag * dz / r2; break;
                case PA: v = g; break;
                case PS: v = Ag * (dx * dx + dy * dy) / s3; break;
                case PR: v = Ag * (dz * dz) / r3; break;
                default: v = 1.0; break;
            }
            J[(size_t)j * w.n + i] = v;
        }
    }
}

// -- MINPACK's building blocks (Moré, Garbow, Hillstrom 1980) ----------------------------------
// enorm, qrfac, qrsolv and lmpar below, and fit's iteration, follow MINPACK-1's lmder and its
// routines (one IEEE operation per C operation; the module is built without fused multiply-adds).
// Their notice:
//   Minpack Copyright Notice (1999) University of Chicago. All rights reserved
//   Redistribution and use in source and binary forms, with or without modification, are
//   permitted provided that the following conditions are met:
//   1. Redistributions of source code must retain the above copyright notice, this list of
//   conditions and the following disclaimer.
//   2. Redistributions in binary form must reproduce the above copyright notice, this list of
//   conditions and the following disclaimer in the documentation and/or other materials provided
//   with the distribution.
//   3. The end-user documentation included with the redistribution, if any, must include the
//   following acknowledgment: "This product includes software developed by the University of
//   Chicago, as Operator of Argonne National Laboratory." Alternately, this acknowledgment may
//   appear in the software itself, if and wherever such third-party acknowledgments normally
//   appear.
//   4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF ANY KIND. THE
//   COPYRIGHT HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND THEIR
//   EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO ANY
//   IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE OR
//   NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
//   COMPLETENESS, OR USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF THE SOFTWARE
//   WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4) DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
//   UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL BE CORRECTED.
//   5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT HOLDER, THE UNITED STATES, THE
//   UNITED STATES DEPARTMENT OF ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
//   INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF ANY KIND OR NATURE, INCLUDING BUT
//   NOT LIMITED TO LOSS OF PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER SUCH
//   LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT (INCLUDING NEGLIGENCE OR STRICT
//   LIABILITY), OR OTHERWISE, EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE POSSIBILITY OF
//   SUCH LOSS OR DAMAGES.

// the Euclidean norm (enorm: sums of squares in three ranges against overflow and underflow)
inline double enorm(int n, const double *x) {
    const double rdwarf = 3.834e-20, rgiant = 1.304e19;
    double s1 = 0, s2 = 0, s3 = 0, x1max = 0, x3max = 0;
    const double agiant = rgiant / n;
    for (int i = 0; i < n; i++) {
        const double xabs = std::fabs(x[i]);
        if (xabs > rdwarf && xabs < agiant) {
            s2 += xabs * xabs;
        } else if (xabs <= rdwarf) {
            if (xabs > x3max) {
                const double t = x3max / xabs;
                s3 = 1 + s3 * t * t;
                x3max = xabs;
            } else if (xabs != 0) {
                const double t = xabs / x3max;
                s3 += t * t;
            }
        } else {
            if (xabs > x1max) {
                const double t = x1max / xabs;
                s1 = 1 + s1 * t * t;
                x1max = xabs;
            } else {
                const double t = xabs / x1max;
                s1 += t * t;
            }
        }
    }
    if (s1 != 0) return x1max * std::sqrt(s1 + (s2 / x1max) / x1max);
    if (s2 != 0) {
        if (s2 >= x3max) return std::sqrt(s2 * (1 + (x3max / s2) * (x3max * s3)));
        return std::sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
    }
    return x3max * std::sqrt(s3);
}

// QR with column pivoting by Householder reflections (qrfac): a is m x n column-major
// (overwritten: R above the diagonal, the reflections' vectors below), ipvt the pivots,
// rdiag the diagonal of R, acnorm the columns' norms
inline void qrfac(int m, int n, double *a, int *ipvt, double *rdiag, double *acnorm, double *wa) {
    const double epsmch = std::numeric_limits<double>::epsilon();
    for (int j = 0; j < n; j++) {
        acnorm[j] = enorm(m, a + (size_t)j * m);
        rdiag[j] = acnorm[j];
        wa[j] = rdiag[j];
        ipvt[j] = j;
    }
    const int minmn = std::min(m, n);
    for (int j = 0; j < minmn; j++) {
        int kmax = j;
        for (int k = j; k < n; k++)
            if (rdiag[k] > rdiag[kmax]) kmax = k;
        if (kmax != j) {
            for (int i = 0; i < m; i++) std::swap(a[(size_t)j * m + i], a[(size_t)kmax * m + i]);
            rdiag[kmax] = rdiag[j];
            wa[kmax] = wa[j];
            std::swap(ipvt[j], ipvt[kmax]);
        }
        double *aj = a + (size_t)j * m;
        double ajnorm = enorm(m - j, aj + j);
        if (ajnorm != 0) {
            if (aj[j] < 0) ajnorm = -ajnorm;
            for (int i = j; i < m; i++) aj[i] /= ajnorm;
            aj[j] += 1.0;
            for (int k = j + 1; k < n; k++) {
                double *ak = a + (size_t)k * m;
                double sum = 0;
                for (int i = j; i < m; i++) sum += aj[i] * ak[i];
                const double temp = sum / aj[j];
                for (int i = j; i < m; i++) ak[i] -= temp * aj[i];
                if (rdiag[k] != 0) {
                    double t = ak[j] / rdiag[k];
                    rdiag[k] *= std::sqrt(std::max(0.0, 1 - t * t));
                    t = rdiag[k] / wa[k];
                    if (0.05 * t * t <= epsmch) {
                        rdiag[k] = enorm(m - j - 1, ak + j + 1);
                        wa[k] = rdiag[k];
                    }
                }
            }
        }
        rdiag[j] = -ajnorm;
    }
}

// solve the least-squares system [R; diag(d)] z = [qtb; 0] given R's QR (qrsolv); r is n x n
// column-major with ldr = n, its upper triangle R, the lower triangle overwritten
inline void qrsolv(int n, double *r, const int *ipvt, const double *diag, const double *qtb, double *x,
                   double *sdiag, double *wa) {
    for (int j = 0; j < n; j++) {
        for (int i = j; i < n; i++) r[(size_t)j * n + i] = r[(size_t)i * n + j];
        x[j] = r[(size_t)j * n + j];
        wa[j] = qtb[j];
    }
    for (int j = 0; j < n; j++) {
        const int l = ipvt[j];
        if (diag[l] != 0) {
            for (int k = j; k < n; k++) sdiag[k] = 0;
            sdiag[j] = diag[l];
            double qtbpj = 0;
            for (int k = j; k < n; k++) {
                if (sdiag[k] == 0) continue;
                double sin_, cos_;
                const double rkk = r[(size_t)k * n + k];
                if (std::fabs(rkk) < std::fabs(sdiag[k])) {
                    const double cotan = rkk / sdiag[k];
                    sin_ = 0.5 / std::sqrt(0.25 + 0.25 * cotan * cotan);
                    cos_ = sin_ * cotan;
                } else {
                    const double tan_ = sdiag[k] / rkk;
                    cos_ = 0.5 / std::sqrt(0.25 + 0.25 * tan_ * tan_);
                    sin_ = cos_ * tan_;
                }
                r[(size_t)k * n + k] = cos_ * rkk + sin_ * sdiag[k];
                const double temp = cos_ * wa[k] + sin_ * qtbpj;
                qtbpj = -sin_ * wa[k] + cos_ * qtbpj;
                wa[k] = temp;
                for (int i = k + 1; i < n; i++) {
                    const double t2 = cos_ * r[(size_t)k * n + i] + sin_ * sdiag[i];
                    sdiag[i] = -sin_ * r[(size_t)k * n + i] + cos_ * sdiag[i];
                    r[(size_t)k * n + i] = t2;
                }
            }
        }
        sdiag[j] = r[(size_t)j * n + j];
        r[(size_t)j * n + j] = x[j];
    }
    int nsing = n;
    for (int j = 0; j < n; j++) {
        if (sdiag[j] == 0 && nsing == n) nsing = j;
        if (nsing < n) wa[j] = 0;
    }
    for (int k = 0; k < nsing; k++) {
        const int j = nsing - 1 - k;
        double sum = 0;
        for (int i = j + 1; i < nsing; i++) sum += r[(size_t)j * n + i] * wa[i];
        wa[j] = (wa[j] - sum) / sdiag[j];
    }
    for (int j = 0; j < n; j++) x[ipvt[j]] = wa[j];
}

// the Levenberg-Marquardt parameter for the trust region delta (lmpar); x is the Gauss-Newton
// or damped direction (its negative is the step), r as in qrsolv
inline void lmpar(int n, double *r, const int *ipvt, const double *diag, const double *qtb, double delta,
                  double &par, double *x, double *sdiag, double *wa1, double *wa2) {
    const double dwarf = std::numeric_limits<double>::min();
    int nsing = n;
    for (int j = 0; j < n; j++) {
        wa1[j] = qtb[j];
        if (r[(size_t)j * n + j] == 0 && nsing == n) nsing = j;
        if (nsing < n) wa1[j] = 0;
    }
    for (int k = 0; k < nsing; k++) {
        const int j = nsing - 1 - k;
        wa1[j] /= r[(size_t)j * n + j];
        const double temp = wa1[j];
        for (int i = 0; i < j; i++) wa1[i] -= r[(size_t)j * n + i] * temp;
    }
    for (int j = 0; j < n; j++) x[ipvt[j]] = wa1[j];
    int iter = 0;
    for (int j = 0; j < n; j++) wa2[j] = diag[j] * x[j];
    double dxnorm = enorm(n, wa2);
    double fp = dxnorm - delta;
    if (fp <= 0.1 * delta) {
        par = 0;
        return;
    }
    double parl = 0;
    if (nsing >= n) {
        for (int j = 0; j < n; j++) {
            const int l = ipvt[j];
            wa1[j] = diag[l] * (wa2[l] / dxnorm);
        }
        for (int j = 0; j < n; j++) {
            double sum = 0;
            for (int i = 0; i < j; i++) sum += r[(size_t)j * n + i] * wa1[i];
            wa1[j] = (wa1[j] - sum) / r[(size_t)j * n + j];
        }
        const double temp = enorm(n, wa1);
        parl = ((fp / delta) / temp) / temp;
    }
    for (int j = 0; j < n; j++) {
        double sum = 0;
        for (int i = 0; i <= j; i++) sum += r[(size_t)j * n + i] * qtb[i];
        wa1[j] = sum / diag[ipvt[j]];
    }
    const double gnorm = enorm(n, wa1);
    double paru = gnorm / delta;
    if (paru == 0) paru = dwarf / std::min(delta, 0.1);
    par = std::max(par, parl);
    par = std::min(par, paru);
    if (par == 0) par = gnorm / dxnorm;
    for (;;) {
        iter++;
        if (par == 0) par = std::max(dwarf, 0.001 * paru);
        const double temp = std::sqrt(par);
        for (int j = 0; j < n; j++) wa1[j] = temp * diag[j];
        qrsolv(n, r, ipvt, wa1, qtb, x, sdiag, wa2);
        for (int j = 0; j < n; j++) wa2[j] = diag[j] * x[j];
        dxnorm = enorm(n, wa2);
        const double temp0 = fp;
        fp = dxnorm - delta;
        if (std::fabs(fp) <= 0.1 * delta || (parl == 0 && fp <= temp0 && temp0 < 0) || iter == 10) break;
        for (int j = 0; j < n; j++) {
            const int l = ipvt[j];
            wa1[j] = diag[l] * (wa2[l] / dxnorm);
        }
        for (int j = 0; j < n; j++) {
            wa1[j] /= sdiag[j];
            const double t = wa1[j];
            for (int i = j + 1; i < n; i++) wa1[i] -= r[(size_t)j * n + i] * t;
        }
        const double t2 = enorm(n, wa1);
        const double parc = ((fp / delta) / t2) / t2;
        if (fp > 0) parl = std::max(parl, par);
        if (fp < 0) paru = std::min(paru, par);
        par = std::max(parl, par + parc);
    }
    if (iter == 0) par = 0;
}

// -- the fit -----------------------------------------------------------------------------------

struct Options {
    double epsabs = 1e-8, epsrel = 1e-8, factor = 100.0;
    int maxiter = 500;
};

struct Result {
    double prm[7];
    std::vector<double> C;  // np x np
    std::vector<double> f;  // residuals at the result
    int iterations = 0, nfev = 0, status = 0;
};

// QR of J (copied), then (R'R)^-1 with the pivots: the covariance inv(J'J)
inline void covariance(int m, int n, const std::vector<double> &J, std::vector<double> &C) {
    std::vector<double> a(J);
    std::vector<int> ipvt(n);
    std::vector<double> rdiag(n), acnorm(n), wa(n);
    qrfac(m, n, a.data(), ipvt.data(), rdiag.data(), acnorm.data(), wa.data());
    // R (n x n) upper triangular: R(i, j) = a[j*m + i] for i < j, R(j, j) = rdiag[j]
    std::vector<double> R((size_t)n * n, 0.0), Ri((size_t)n * n, 0.0);
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < j; i++) R[(size_t)j * n + i] = a[(size_t)j * m + i];
        R[(size_t)j * n + j] = rdiag[j];
    }
    // Ri = R^-1 (upper triangular), column by column
    for (int j = 0; j < n; j++) {
        Ri[(size_t)j * n + j] = 1.0 / R[(size_t)j * n + j];
        for (int i = j - 1; i >= 0; i--) {
            double sum = 0;
            for (int k = i + 1; k <= j; k++) sum += R[(size_t)k * n + i] * Ri[(size_t)j * n + k];
            Ri[(size_t)j * n + i] = -sum / R[(size_t)i * n + i];
        }
    }
    // (R'R)^-1 = Ri Ri', then undo the pivots: C(ipvt i, ipvt j)
    C.assign((size_t)n * n, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0;
            for (int k = std::max(i, j); k < n; k++) sum += Ri[(size_t)k * n + i] * Ri[(size_t)k * n + j];
            C[(size_t)ipvt[j] * n + ipvt[i]] = sum;
        }
}

// one fit: Moré's trust-region Levenberg-Marquardt on the free parameters
inline Result fit(const Win &w, const double *p0, const int *free, int np, const Options &o) {
    Result res;
    std::memcpy(res.prm, p0, sizeof(res.prm));
    const int m = w.n, n = np;
    const double eps = std::numeric_limits<double>::epsilon();
    std::vector<double> ex, ey, ez, f, ft, J, a, qtf(n), diag(n), wa1(n), wa2(n), wa3(n), wa4(m), rdiag(n),
        acnorm(n), sdiag(n), rr((size_t)n * n);
    std::vector<int> ipvt(n);
    double p[7], pt[7];
    std::memcpy(p, p0, sizeof(p));
    residuals(w, p, f, ex, ey, ez);
    jacobian(w, p, free, n, J, ex, ey, ez);
    res.nfev = 1;
    double fnorm = enorm(m, f.data());
    double par = 0, delta = 0, xnorm = 0;
    int iter = 1;
    bool first = true;
    for (int outer = 0; outer < o.maxiter; outer++) {
        // QR of J, diag and delta on the first iteration, qtf = Q' f
        a = J;
        qrfac(m, n, a.data(), ipvt.data(), rdiag.data(), acnorm.data(), wa2.data());
        if (first) {
            for (int j = 0; j < n; j++) {
                diag[j] = acnorm[j] == 0 ? 1.0 : acnorm[j];
                wa3[j] = diag[j] * p[free[j]];
            }
            xnorm = enorm(n, wa3.data());
            delta = o.factor * xnorm;
            if (delta == 0) delta = o.factor;
            first = false;
        }
        std::memcpy(wa4.data(), f.data(), sizeof(double) * m);
        for (int j = 0; j < n; j++) {
            double *aj = a.data() + (size_t)j * m;
            if (aj[j] != 0) {
                double sum = 0;
                for (int i = j; i < m; i++) sum += aj[i] * wa4[i];
                const double temp = -sum / aj[j];
                for (int i = j; i < m; i++) wa4[i] += aj[i] * temp;
            }
            aj[j] = rdiag[j];
            qtf[j] = wa4[j];
        }
        double gnorm = 0;
        if (fnorm != 0) {
            for (int j = 0; j < n; j++) {
                const int l = ipvt[j];
                if (acnorm[l] != 0) {
                    double sum = 0;
                    for (int i = 0; i <= j; i++) sum += a[(size_t)j * m + i] * (qtf[i] / fnorm);
                    gnorm = std::max(gnorm, std::fabs(sum / acnorm[l]));
                }
            }
        }
        for (int j = 0; j < n; j++) diag[j] = std::max(diag[j], acnorm[j]);
        // R as an n x n column-major matrix for lmpar
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) rr[(size_t)j * n + i] = a[(size_t)j * m + i];
        bool accepted = false;
        int status = 0, rejected = 0;
        for (;;) {
            lmpar(n, rr.data(), ipvt.data(), diag.data(), qtf.data(), delta, par, wa1.data(), sdiag.data(),
                  wa2.data(), wa3.data());
            // restore R's upper triangle overwritten by qrsolv (it keeps the strict upper part)
            for (int j = 0; j < n; j++) rr[(size_t)j * n + j] = a[(size_t)j * m + j];
            std::memcpy(pt, p, sizeof(pt));
            for (int j = 0; j < n; j++) {
                wa1[j] = -wa1[j];
                pt[free[j]] = p[free[j]] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            const double pnorm = enorm(n, wa3.data());
            if (iter == 1) delta = std::min(delta, pnorm);
            residuals(w, pt, ft, ex, ey, ez);
            res.nfev++;
            const double fnorm1 = enorm(m, ft.data());
            double actred = -1;
            if (0.1 * fnorm1 < fnorm) {
                const double t = fnorm1 / fnorm;
                actred = 1 - t * t;
            }
            for (int j = 0; j < n; j++) {
                wa3[j] = 0;
                const double temp = wa1[ipvt[j]];
                for (int i = 0; i <= j; i++) wa3[i] += a[(size_t)j * m + i] * temp;
            }
            const double temp1 = enorm(n, wa3.data()) / fnorm;
            const double temp2 = (std::sqrt(par) * pnorm) / fnorm;
            const double prered = temp1 * temp1 + temp2 * temp2 / 0.5;
            const double dirder = -(temp1 * temp1 + temp2 * temp2);
            const double ratio = prered != 0 ? actred / prered : 0;
            if (ratio <= 0.25) {
                double temp = actred >= 0 ? 0.5 : 0.5 * dirder / (dirder + 0.5 * actred);
                if (0.1 * fnorm1 >= fnorm || temp < 0.1) temp = 0.1;
                delta = temp * std::min(delta, pnorm / 0.1);
                par /= temp;
            } else if (par == 0 || ratio >= 0.75) {
                delta = pnorm / 0.5;
                par = 0.5 * par;
            }
            if (ratio >= 1e-4) {
                std::memcpy(p, pt, sizeof(p));
                f.swap(ft);
                for (int j = 0; j < n; j++) wa2[j] = diag[j] * p[free[j]];
                xnorm = enorm(n, wa2.data());
                fnorm = fnorm1;
                iter++;
                jacobian(w, p, free, n, J, ex, ey, ez);
                res.nfev++;
                accepted = true;
                // the step taken, for the convergence test
                break;
            }
            if (std::fabs(actred) <= eps && prered <= eps && 0.5 * ratio <= 1) {
                status = 1;  // the residuals cannot be reduced further
                break;
            }
            if (delta <= eps * xnorm) {
                status = 2;
                break;
            }
            if (gnorm <= eps) {
                status = 3;
                break;
            }
            if (++rejected >= 10) {  // ten rejected trials in a row end the fit (measured)
                status = 4;
                break;
            }
        }
        res.iterations = outer + 1;
        if (!accepted) {
            res.status = status;
            break;
        }
        bool conv = true;
        for (int j = 0; j < n; j++)
            if (!(std::fabs(wa1[j]) < o.epsabs + o.epsrel * std::fabs(p[free[j]]))) conv = false;
        if (conv) break;
    }
    std::memcpy(res.prm, p, sizeof(p));
    res.f = f;
    covariance(m, n, J, res.C);
    return res;
}

}  // namespace mexfit
