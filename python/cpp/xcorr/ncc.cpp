// The rest of normxcorr3_fast (MATLAB 4a3bef3) around the mex's window: the means of the template and the
// image, the template's sum and standard deviation and the image's variance as MATLAB's fixed_order_stat
// computes them over x(:), the NCC in double, the first maximum, and the moves between the port's (z, y, x)
// C arrays and MATLAB's layout. normxcorr3_fast is the lab's adaptation of Daniel Eaton's normxcorr3 (MATLAB
// File Exchange 73946, BSD license: THIRD_PARTY_NOTICES.md), after J. P. Lewis, "Fast Normalized
// Cross-Correlation", Vision Interface 1995.
//
// MATLAB's reductions over x(:) at one thread follow the vector rule (docs/decisions/0005): in index order up
// to 43 elements; up to 88,999 in four lanes (element t in lane t mod 4, each lane in index order, added as
// ((l0 + l1) + l2) + l3); from 89,000 in blocks of 4,096 elements, each in four lanes, the block sums in index
// order. Every addition is one IEEE operation in the class of the sum; the blocks may be summed by any thread.
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "xcorr.h"

namespace pk5d::xcorr {
namespace {

constexpr idx_t kSequentialUpto = 43;
constexpr idx_t kBlocksFrom = 89000;
constexpr idx_t kBlock = 4096;

// elements b .. e - 1 (e > b) in four lanes (lane of element t: (t - b) mod 4; b is a multiple of 4), each
// lane from its first element (-0 for an empty lane: x + -0 is x), added as ((l0 + l1) + l2) + l3
template <typename A, typename F>
A four_lanes(const idx_t b, const idx_t e, const F& f) {
    A l[4];
    for (int k = 0; k < 4; k++) l[k] = b + k < e ? f(b + k) : A(-0.0);
    for (idx_t i = b + 4; i < e; i += 4) {
        const idx_t m = e - i < 4 ? e - i : 4;
        for (idx_t k = 0; k < m; k++) l[k] = l[k] + f(i + k);
    }
    return ((l[0] + l[1]) + l[2]) + l[3];
}

// MATLAB's sum of the n addends f(0) .. f(n - 1) in the class A (the vector rule)
template <typename A, typename F>
A vector_sum(const idx_t n, const F& f, const int nt) {
    if (n <= 0) return A(0);
    if (n <= kSequentialUpto) {
        A s = f(0);
        for (idx_t i = 1; i < n; i++) s = s + f(i);
        return s;
    }
    if (n < kBlocksFrom) return four_lanes<A>(0, n, f);
    const idx_t nb = (n + kBlock - 1) / kBlock;
    std::vector<A> bs(nb);
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t k = 0; k < nb; k++) {
        const idx_t e = (k + 1) * kBlock < n ? (k + 1) * kBlock : n;
        bs[k] = four_lanes<A>(k * kBlock, e, f);
    }
    A s = bs[0];
    for (idx_t k = 1; k < nb; k++) s = s + bs[k];
    return s;
}

template <typename T>
inline T eps_of() {
    return std::numeric_limits<T>::epsilon();
}

}  // namespace

template <typename T>
Stats prepare(const T* Tm, const idx_t nT, const T* Am, const idx_t nA, T* T0, T* A0, const int nt) {
    Stats st;
    // cT = cast(mean(T(:), 'double'), 'like', T); cA likewise
    const double sT = vector_sum<double>(nT, [&](idx_t i) { return (double)Tm[i]; }, nt);
    const double sA = vector_sum<double>(nA, [&](idx_t i) { return (double)Am[i]; }, nt);
    const T cT = (T)(sT / (double)nT);
    const T cA = (T)(sA / (double)nA);
    st.cT = cT;
    st.cA = cA;
    // var(A(:), 1) in A's class: mean in the class (the sum's double / n, rounded), then the sum of the squared
    // absolute deviations in the class, divided by n in double and rounded
    const T sAc = vector_sum<T>(nA, [&](idx_t i) { return Am[i]; }, nt);
    const T q = (T)((double)sAc / (double)nA);
    const T s2 = vector_sum<T>(nA, [&](idx_t i) {
        const T d = std::abs((T)(Am[i] - q));
        return (T)(d * d);
    }, nt);
    st.varA1 = (T)((double)s2 / (double)nA);
    // T0 = T - cT, A0 = A - cA (in the class)
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t i = 0; i < nA; i++) A0[i] = (T)(Am[i] - cA);
#pragma omp parallel for num_threads(nt) schedule(static) if (nT > 65536)
    for (idx_t i = 0; i < nT; i++) T0[i] = (T)(Tm[i] - cT);
    // sum(T0(:), 'double'), std(double(T0(:))) (the mean s / n, the squared absolute deviations, / (n - 1); one
    // element: 0)
    st.sumT0 = vector_sum<double>(nT, [&](idx_t i) { return (double)T0[i]; }, nt);
    const double mu = st.sumT0 / (double)nT;
    const double ss = vector_sum<double>(nT, [&](idx_t i) {
        const double d = std::abs((double)T0[i] - mu);
        return d * d;
    }, nt);
    st.denomT = std::sqrt(ss / (nT > 1 ? (double)(nT - 1) : 1.0));
    return st;
}

template <typename T, typename C>
void ncc_combine(const C* corr, const double* sA, const double* sA2, const double* sT, const int* ts, const int* a,
    const Win& W, const Stats& st, T* out, const int nt) {
    const double pSzT = (double)ts[0] * (double)ts[1] * (double)ts[2];
    // the image's elements a template at each window position overlaps, per dimension
    std::vector<double> cnt[3];
    for (int d = 0; d < 3; d++) {
        cnt[d].resize(W.wn[d]);
        for (int j = 0; j < W.wn[d]; j++) {
            const int k = W.w0[d] + j + 1;
            const int lo = k - ts[d] + 1 > 1 ? k - ts[d] + 1 : 1;
            const int hi = k < a[d] ? k : a[d];
            cnt[d][j] = (double)(hi - lo + 1);
        }
    }
    const double cAd = st.cA;
    const double cA2 = cAd * cAd;
    const double sumT0 = st.sumT0;
    const double mT = sumT0 / pSzT;
    const double pm1 = pSzT - 1;
    const double denomT = st.denomT;
    // eps(class(A)) * var(A(:), 1), in A's class
    const double thr = (double)(T)(eps_of<T>() * (T)st.varA1);
    const idx_t wy = W.wn[0], wx = W.wn[1], wz = W.wn[2];
#pragma omp parallel for num_threads(nt) schedule(static) collapse(2)
    for (idx_t z = 0; z < wz; z++) {
        for (idx_t x = 0; x < wx; x++) {
            const double c12 = cnt[1][x];
            const double c3 = cnt[2][z];
            const idx_t L0 = wy * (x + wx * z);
            for (idx_t y = 0; y < wy; y++) {
                const idx_t L = L0 + y;
                // nPad = pSzT - cnt{1}(:) .* cnt{2}(:)' .* reshape(cnt{3}, 1, 1, [])
                const double nPad = pSzT - (cnt[0][y] * c12) * c3;
                const double sB = sA[L] - nPad * cAd;
                const double sB2 = sA2[L] + nPad * cA2;
                const double num = (((double)corr[L] - cAd * (sumT0 - sT[L])) - sB * mT) / pm1;
                const double q = sB2 - (sB * sB) / pSzT;
                // max(q, 0): NaN gives 0
                const double varA = (q > 0 ? q : 0.0) / pm1;
                out[L] = (varA > thr && denomT > 0) ? (T)(num / (denomT * std::sqrt(varA))) : T(0);
            }
        }
    }
}

template <typename T>
idx_t first_max(const T* C, const idx_t n, const int nt) {
    if (n <= 0) return 0;
    // pieces in order; the first largest of each piece, then the first largest of the pieces
    const idx_t np = n < ((idx_t)1 << 16) ? 1 : (nt > 1 ? 4 * (idx_t)nt : 1);
    const idx_t len = (n + np - 1) / np;
    std::vector<idx_t> best(np, -1);
#pragma omp parallel for num_threads(nt) schedule(static) if (np > 1)
    for (idx_t p = 0; p < np; p++) {
        const idx_t e = (p + 1) * len < n ? (p + 1) * len : n;
        idx_t bi = -1;
        T bv = 0;
        for (idx_t i = p * len; i < e; i++) {
            const T v = C[i];
            if (v != v) continue;
            if (bi < 0 || v > bv) {
                bi = i;
                bv = v;
            }
        }
        best[p] = bi;
    }
    idx_t bi = -1;
    for (idx_t p = 0; p < np; p++) {
        const idx_t i = best[p];
        if (i >= 0 && (bi < 0 || C[i] > C[bi])) bi = i;
    }
    return bi < 0 ? 0 : bi;
}

// tiles of 32 x 32 per plane: dst (MATLAB: y fastest) from src (C: x fastest), or back
template <typename T>
void c_to_matlab(const T* src, const idx_t nz, const idx_t ny, const idx_t nx, T* dst, const int nt) {
    constexpr idx_t B = 32;
    const idx_t ty = (ny + B - 1) / B, tx = (nx + B - 1) / B;
    const idx_t plane = ny * nx;
#pragma omp parallel for num_threads(nt) schedule(static) collapse(2) if (nz * plane > 65536)
    for (idx_t z = 0; z < nz; z++) {
        for (idx_t t = 0; t < ty * tx; t++) {
            const idx_t y0 = (t / tx) * B, x0 = (t % tx) * B;
            const idx_t y1 = y0 + B < ny ? y0 + B : ny, x1 = x0 + B < nx ? x0 + B : nx;
            const T* s = src + plane * z;
            T* d = dst + plane * z;
            for (idx_t x = x0; x < x1; x++) {
                for (idx_t y = y0; y < y1; y++) d[y + ny * x] = s[x + nx * y];
            }
        }
    }
}

template <typename T>
void matlab_to_c(const T* src, const idx_t nz, const idx_t ny, const idx_t nx, T* dst, const int nt) {
    constexpr idx_t B = 32;
    const idx_t ty = (ny + B - 1) / B, tx = (nx + B - 1) / B;
    const idx_t plane = ny * nx;
#pragma omp parallel for num_threads(nt) schedule(static) collapse(2) if (nz * plane > 65536)
    for (idx_t z = 0; z < nz; z++) {
        for (idx_t t = 0; t < ty * tx; t++) {
            const idx_t y0 = (t / tx) * B, x0 = (t % tx) * B;
            const idx_t y1 = y0 + B < ny ? y0 + B : ny, x1 = x0 + B < nx ? x0 + B : nx;
            const T* s = src + plane * z;
            T* d = dst + plane * z;
            for (idx_t y = y0; y < y1; y++) {
                for (idx_t x = x0; x < x1; x++) d[x + nx * y] = s[y + ny * x];
            }
        }
    }
}

template Stats prepare<float>(const float*, idx_t, const float*, idx_t, float*, float*, int);
template Stats prepare<double>(const double*, idx_t, const double*, idx_t, double*, double*, int);
template void ncc_combine<float, float>(const float*, const double*, const double*, const double*, const int*,
    const int*, const Win&, const Stats&, float*, int);
template void ncc_combine<double, double>(const double*, const double*, const double*, const double*, const int*,
    const int*, const Win&, const Stats&, double*, int);
template idx_t first_max<float>(const float*, idx_t, int);
template idx_t first_max<double>(const double*, idx_t, int);
template void c_to_matlab<float>(const float*, idx_t, idx_t, idx_t, float*, int);
template void c_to_matlab<double>(const double*, idx_t, idx_t, idx_t, double*, int);
template void matlab_to_c<float>(const float*, idx_t, idx_t, idx_t, float*, int);
template void matlab_to_c<double>(const double*, idx_t, idx_t, idx_t, double*, int);

}  // namespace pk5d::xcorr
