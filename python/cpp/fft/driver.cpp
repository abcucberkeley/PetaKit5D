// MATLAB R2024b's fft drivers on FFTW (see driver.h for the rules this reproduces).
#include "driver.h"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <map>
#include <stdexcept>
#include <tuple>

#include "backend.h"

namespace pk5d::fft {
namespace {

enum Kind { kR2C = 0, kC2CF = 1, kC2CB = 2, kC2R = 3 };

constexpr unsigned kEstimate = 1U << 6;  // FFTW_ESTIMATE, MATLAB's planner flags
constexpr int kB = 4;                    // lines per execute of MATLAB's grouped plans
constexpr int kInOff = 16;               // MATLAB's scratch offsets (bytes modulo 64)
constexpr int kOutOff = 32;

// --- buffers: 64-byte aligned base plus an offset (FFTW sees the alignment modulo 16) ---------------------

struct Buf {
    void* base = nullptr;
    char* p = nullptr;
    Buf() = default;
    Buf(const Buf&) = delete;
    Buf& operator=(const Buf&) = delete;
    ~Buf() { std::free(base); }
    void alloc(size_t bytes, int off) {
        std::free(base);
        base = std::malloc(bytes + 128);
        if (base == nullptr) throw std::bad_alloc();
        p = reinterpret_cast<char*>(((reinterpret_cast<uintptr_t>(base) + 63) & ~static_cast<uintptr_t>(63)) + off);
    }
};

// --- plans: rank 1, contiguous lines, FFTW_ESTIMATE, out of place, as MATLAB makes them --------------------

using PlanKey = std::tuple<int, int, idx_t, int>;  // precision, kind, n, lines per execute

std::map<PlanKey, void*>& plan_cache() {
    static std::map<PlanKey, void*> m;
    return m;
}

void dists(int kind, idx_t n, idx_t& idist, idx_t& odist) {
    const idx_t nh = n / 2 + 1;
    idist = kind == kC2R ? nh : n;
    odist = kind == kR2C ? nh : n;
}

void* get_plan(int prec, int kind, idx_t n, int hm) {
    std::lock_guard<std::mutex> g(planner_lock());
    const PlanKey key(prec, kind, n, hm);
    auto it = plan_cache().find(key);
    if (it != plan_cache().end()) return it->second;
    const FftwApi& a = backend().api[prec];
    idx_t idist, odist;
    dists(kind, n, idist, odist);
    const size_t es = prec ? sizeof(double) : sizeof(float);
    Buf bi, bo;
    bi.alloc(static_cast<size_t>(hm) * idist * (kind == kR2C ? 1 : 2) * es, kInOff);
    bo.alloc(static_cast<size_t>(hm) * odist * (kind == kC2R ? 1 : 2) * es, kOutOff);
    const Iodim d = {n, 1, 1};
    const Iodim h = {hm, idist, odist};
    void* p;
    if (kind == kR2C) {
        p = a.plan_r2c(1, &d, 1, &h, bi.p, bo.p, kEstimate);
    } else if (kind == kC2R) {
        p = a.plan_c2r(1, &d, 1, &h, bi.p, bo.p, kEstimate);
    } else {
        p = a.plan_dft(1, &d, 1, &h, bi.p, bo.p, kind == kC2CF ? -1 : 1, kEstimate);
    }
    if (p == nullptr) throw std::runtime_error("FFTW planning failed");
    plan_cache()[key] = p;
    return p;
}

inline void execute(int prec, int kind, void* plan, void* in, void* out) {
    const FftwApi& a = backend().api[prec];
    if (kind == kR2C) {
        a.exec_r2c(plan, in, out);
    } else if (kind == kC2R) {
        a.exec_c2r(plan, in, out);
    } else {
        a.exec_dft(plan, in, out);
    }
}

// --- MATLAB's split of a stage ---------------------------------------------------------------------------

bool stage_split(idx_t grain, idx_t n, idx_t N) {
    int bits = 0;
    for (idx_t v = n; v > 0; v >>= 1) bits++;
    return grain > 0 && static_cast<double>(n) * static_cast<double>(N) * bits > 98304.0;
}

// line j of N alone, given the effective grain (0: no split) and the segment period seg (0: none)
inline bool alone_g(idx_t j, idx_t N, idx_t g, idx_t seg) {
    idx_t a = 0, b = N;
    if (g > 0) {
        while (b - a > g) {
            const idx_t mid = a + (b - a) / 2;
            if (j < mid) {
                b = mid;
            } else {
                a = mid;
            }
        }
    }
    idx_t s = a, e = b;
    if (seg > 0) {
        const idx_t s0 = (j / seg) * seg;
        s = std::max(a, s0);
        e = std::min(b, s0 + seg);
    }
    return j >= e - (e - s) % kB;
}

// --- one stage: the lines along dimension d ----------------------------------------------------------------

template <typename R>
struct Stage {
    int prec;
    int kind;
    // source: real or complex, dims S (Fortran order); padded with zeros or truncated to the destination's
    // extent in every dimension
    const R* src;
    bool src_complex;
    std::vector<idx_t> S;
    // destination: dims D; D[d] is the transform length n (C2C, R2C: n complex values per line with the
    // conjugate fill; C2R: n real values)
    R* dst;
    std::vector<idx_t> D;
    int d;
    idx_t n;
    // MATLAB's line count of the stage is L * Ufull; only lines with i_hi < D's U are computed (the
    // destination may hold a prefix of the full array along its slowest dimensions)
    idx_t Ufull;
    bool scale;       // divide by n (inverse stages)
    idx_t kwrite = -1;  // outputs written per line (the first ones; -1: all)
};

// per-thread scratch of a stage
template <typename R>
struct Scratch {
    Buf tin, tout, in1, out1;
    std::vector<idx_t> sbase, slot;
    std::vector<int> alone;
    std::vector<idx_t> coord;
};

// the geometry of a stage, shared by its tiles
struct Geom {
    int nd, d;
    idx_t n, nh, L, U, Nfull, g, seg, idist, odist, need, avail, sstep, T, tilesPerU, ntiles, kw;
    bool same;
    std::vector<idx_t> sstride;
};

// one tile of adjacent lines: gather (zero-padded), MATLAB's groups of four and lines alone, scatter.
// KIND and the source's complexity are compile-time, so the copy loops have no branches.
template <typename R, int KIND, bool SRCC>
void run_tile(const Stage<R>& st, const Geom& G, Scratch<R>& w, void* p4, void* p1, const idx_t q) {
    constexpr int ci = KIND == kR2C ? 1 : 2;   // R values per scratch input element
    constexpr int co = KIND == kC2R ? 1 : 2;   // per scratch output element
    constexpr int cs = SRCC ? 2 : 1;           // per source element
    constexpr int cd = KIND == kC2R ? 1 : 2;   // per destination element
    const idx_t L = G.L, n = G.n, nh = G.nh, idist = G.idist, odist = G.odist;
    idx_t ihi0, ilo0, Tq;
    if (L > 1) {
        ihi0 = q / G.tilesPerU;
        ilo0 = (q % G.tilesPerU) * G.T;
        Tq = std::min(G.T, L - ilo0);
    } else {
        ihi0 = q * G.T;
        ilo0 = 0;
        Tq = std::min(G.T, G.U - ihi0);
    }
    // source line offsets (-1: zero line), MATLAB's alone flags, slots (grouped lines first)
    idx_t ng = 0;
    for (idx_t t = 0; t < Tq; t++) {
        const idx_t ilo = L > 1 ? ilo0 + t : 0;
        const idx_t ihi = L > 1 ? ihi0 : ihi0 + t;
        w.alone[t] = alone_g(ilo + L * ihi, G.Nfull, G.g, G.seg) ? 1 : 0;
        if (!w.alone[t]) ng++;
        if (G.same) {
            w.sbase[t] = ilo + L * st.S[G.d] * ihi;
        } else {
            idx_t r = ilo;
            for (int k = 0; k < G.d; k++) {
                w.coord[k] = r % st.D[k];
                r /= st.D[k];
            }
            r = ihi;
            for (int k = G.d + 1; k < G.nd; k++) {
                w.coord[k] = r % st.D[k];
                r /= st.D[k];
            }
            idx_t off = 0;
            bool in = true;
            for (int k = 0; k < G.nd; k++) {
                if (k == G.d) continue;
                if (w.coord[k] >= st.S[k]) in = false;
                off += w.coord[k] * G.sstride[k];
            }
            w.sbase[t] = in ? off : -1;
        }
    }
    const idx_t g4 = (ng + kB - 1) / kB * kB;
    {
        idx_t sg = 0, sa = g4;
        for (idx_t t = 0; t < Tq; t++) w.slot[t] = w.alone[t] ? sa++ : sg++;
    }
    R* const tI = reinterpret_cast<R*>(w.tin.p);
    R* const tO = reinterpret_cast<R*>(w.tout.p);
    const idx_t avail = G.avail, need = G.need, sstep = G.sstep;
    // gather into tI (slot s at s * idist elements)
    if (L > 1) {
        bool zero = false;
        for (idx_t t = 0; t < Tq; t++) zero = zero || w.sbase[t] < 0;
        if (!zero && G.same) {
            // adjacent lines: element k of line t at src[(sbase[0] + t + L * k) * cs]
            const R* s0 = st.src + w.sbase[0] * cs;
            for (idx_t k = 0; k < avail; k++) {
                const R* s = s0 + L * k * cs;
                for (idx_t t = 0; t < Tq; t++) {
                    R* e = tI + (w.slot[t] * idist + k) * ci;
                    e[0] = s[t * cs];
                    if constexpr (ci == 2) e[1] = SRCC ? s[t * cs + 1] : R(0);
                }
            }
        } else {
            for (idx_t k = 0; k < avail; k++) {
                for (idx_t t = 0; t < Tq; t++) {
                    R* e = tI + (w.slot[t] * idist + k) * ci;
                    if (w.sbase[t] < 0) {
                        e[0] = 0;
                        if constexpr (ci == 2) e[1] = 0;
                        continue;
                    }
                    const R* s = st.src + (w.sbase[t] + sstep * k) * cs;
                    e[0] = s[0];
                    if constexpr (ci == 2) e[1] = SRCC ? s[1] : R(0);
                }
            }
        }
    } else {
        for (idx_t t = 0; t < Tq; t++) {
            R* e = tI + w.slot[t] * idist * ci;
            if (w.sbase[t] < 0) {
                std::memset(e, 0, sizeof(R) * avail * ci);
                continue;
            }
            const R* s = st.src + w.sbase[t] * cs;
            if (ci == cs && sstep == 1) {
                std::memcpy(e, s, sizeof(R) * avail * ci);
            } else {
                for (idx_t k = 0; k < avail; k++) {
                    e[k * ci] = s[sstep * k * cs];
                    if constexpr (ci == 2) e[k * ci + 1] = SRCC ? s[sstep * k * cs + 1] : R(0);
                }
            }
        }
    }
    if (avail < need) {
        for (idx_t t = 0; t < Tq; t++) {
            std::memset(tI + (w.slot[t] * idist + avail) * ci, 0, sizeof(R) * (need - avail) * ci);
        }
    }
    // a group of fewer than four lines is filled up with copies of its last line
    for (idx_t s = ng; s < g4; s++) {
        std::memcpy(tI + s * idist * ci, tI + (ng - 1) * idist * ci, sizeof(R) * idist * ci);
    }
    for (idx_t s = 0; s < g4; s += kB) {
        execute(st.prec, KIND, p4, tI + s * idist * ci, tO + s * odist * co);
    }
    for (idx_t s = g4; s < g4 + (Tq - ng); s++) {
        std::memcpy(w.in1.p, tI + s * idist * ci, sizeof(R) * idist * ci);
        execute(st.prec, KIND, p1, w.in1.p, w.out1.p);
        std::memcpy(tO + s * odist * co, w.out1.p, sizeof(R) * odist * co);
    }
    // scatter: element k of line t at dst[(i_lo + L * (k + n * i_hi)) * cd]; the first kw elements of a line
    const idx_t kw = G.kw;
    const R fn = static_cast<R>(n);
    const bool scale = st.scale;
    if (L > 1) {
        R* const db = st.dst + (ilo0 + L * n * ihi0) * cd;
        for (idx_t k = 0; k < kw; k++) {
            R* drow = db + L * k * cd;
            if constexpr (KIND == kR2C) {
                const idx_t kk = k < nh ? k : n - k;
                if (k < nh) {
                    for (idx_t t = 0; t < Tq; t++) {
                        const R* lo = tO + w.slot[t] * odist * co;
                        drow[2 * t] = lo[2 * kk];
                        drow[2 * t + 1] = lo[2 * kk + 1];
                    }
                } else {
                    for (idx_t t = 0; t < Tq; t++) {
                        const R* lo = tO + w.slot[t] * odist * co;
                        drow[2 * t] = lo[2 * kk];
                        drow[2 * t + 1] = -lo[2 * kk + 1];
                    }
                }
            } else if constexpr (KIND == kC2R) {
                for (idx_t t = 0; t < Tq; t++) {
                    const R v = tO[w.slot[t] * odist + k];
                    drow[t] = scale ? v / fn : v;
                }
            } else {
                if (scale) {
                    for (idx_t t = 0; t < Tq; t++) {
                        const R* lo = tO + (w.slot[t] * odist + k) * 2;
                        drow[2 * t] = lo[0] / fn;
                        drow[2 * t + 1] = lo[1] / fn;
                    }
                } else {
                    for (idx_t t = 0; t < Tq; t++) {
                        const R* lo = tO + (w.slot[t] * odist + k) * 2;
                        drow[2 * t] = lo[0];
                        drow[2 * t + 1] = lo[1];
                    }
                }
            }
        }
    } else {
        for (idx_t t = 0; t < Tq; t++) {
            const R* lo = tO + w.slot[t] * odist * co;
            R* e = st.dst + (n * (ihi0 + t)) * cd;
            if constexpr (KIND == kR2C) {
                std::memcpy(e, lo, sizeof(R) * 2 * std::min(nh, kw));
                for (idx_t k = nh; k < kw; k++) {
                    e[2 * k] = lo[2 * (n - k)];
                    e[2 * k + 1] = -lo[2 * (n - k) + 1];
                }
            } else {
                if (!scale) {
                    std::memcpy(e, lo, sizeof(R) * kw * cd);
                } else {
                    for (idx_t k = 0; k < kw * cd; k++) e[k] = lo[k] / fn;
                }
            }
        }
    }
}

template <typename R, int KIND, bool SRCC>
void run_stage_k(const Stage<R>& st, const Geom& G, const Options& o) {
    constexpr int ci = KIND == kR2C ? 1 : 2;
    constexpr int co = KIND == kC2R ? 1 : 2;
    void* p4 = get_plan(st.prec, KIND, G.n, kB);
    void* p1 = get_plan(st.prec, KIND, G.n, 1);
    const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(o.threads, G.ntiles)));
    std::atomic<bool> failed(false);
#pragma omp parallel num_threads(nt)
    {
        Scratch<R> w;
        try {
            const size_t es = sizeof(R);
            w.tin.alloc((2 * G.T + kB) * G.idist * ci * es, kInOff);
            w.tout.alloc((2 * G.T + kB) * G.odist * co * es, kOutOff);
            w.in1.alloc(G.idist * ci * es, kInOff);
            w.out1.alloc(G.odist * co * es, kOutOff);
            w.sbase.resize(static_cast<size_t>(G.T));
            w.slot.resize(static_cast<size_t>(G.T));
            w.alone.resize(static_cast<size_t>(G.T));
            w.coord.resize(static_cast<size_t>(G.nd));
        } catch (...) {
            failed = true;
        }
#pragma omp for schedule(dynamic, 1)
        for (idx_t q = 0; q < G.ntiles; q++) {
            if (failed) continue;
            run_tile<R, KIND, SRCC>(st, G, w, p4, p1, q);
        }
    }
    if (failed) throw std::bad_alloc();
}

template <typename R>
void run_stage(const Stage<R>& st, const Options& o) {
    Geom G;
    G.nd = static_cast<int>(st.D.size());
    G.d = st.d;
    G.n = st.n;
    G.nh = st.n / 2 + 1;
    G.L = 1;
    G.U = 1;
    for (int k = 0; k < G.d; k++) G.L *= st.D[k];
    for (int k = G.d + 1; k < G.nd; k++) G.U *= st.D[k];
    if (G.L * G.U == 0 || G.n == 0) return;
    G.Nfull = G.L * st.Ufull;
    G.g = stage_split(o.grain, G.n, G.Nfull) ? o.grain : 0;
    G.seg = G.L > 1 ? G.L : 0;
    dists(st.kind, G.n, G.idist, G.odist);
    // elements of the source line that are read (the rest is zero)
    G.need = st.kind == kC2R ? G.nh : G.n;
    G.avail = std::min(st.S[G.d], G.need);
    G.sstride.assign(G.nd + 1, 1);
    for (int k = 0; k < G.nd; k++) G.sstride[k + 1] = G.sstride[k] * st.S[k];
    G.same = true;  // the source's layout equals the destination's apart from dimension d
    for (int k = 0; k < G.nd; k++) {
        if (k != G.d && st.S[k] != st.D[k]) G.same = false;
    }
    G.sstep = G.sstride[G.d];
    G.kw = st.kwrite < 0 ? G.n : std::min(st.kwrite, G.n);
    // lines per tile: adjacent lines (consecutive i_lo of one i_hi; for L == 1 consecutive i_hi)
    G.T = std::max<idx_t>(kB, std::min<idx_t>(64, ((16384 / std::max<idx_t>(G.n, 1)) + 3) / 4 * 4));
    if (G.L > 1) G.T = std::min(G.T, G.L);
    if (G.L == 1) G.T = std::min(G.T, G.U);
    G.tilesPerU = G.L > 1 ? (G.L + G.T - 1) / G.T : 1;
    G.ntiles = G.L > 1 ? G.tilesPerU * G.U : (G.U + G.T - 1) / G.T;
    switch (st.kind) {
        case kR2C:
            run_stage_k<R, kR2C, false>(st, G, o);
            break;
        case kC2R:
            if (st.src_complex) {
                run_stage_k<R, kC2R, true>(st, G, o);
            } else {
                run_stage_k<R, kC2R, false>(st, G, o);
            }
            break;
        case kC2CF:
            if (st.src_complex) {
                run_stage_k<R, kC2CF, true>(st, G, o);
            } else {
                run_stage_k<R, kC2CF, false>(st, G, o);
            }
            break;
        default:
            if (st.src_complex) {
                run_stage_k<R, kC2CB, true>(st, G, o);
            } else {
                run_stage_k<R, kC2CB, false>(st, G, o);
            }
            break;
    }
}

// --- MATLAB's checks and fills ----------------------------------------------------------------------------

// linear index of the N-D mirror (-k mod n per dimension) of the column (dims 1 .. nd - 1) with index col
idx_t mirror_column(idx_t col, const std::vector<idx_t>& D) {
    idx_t m = 0, stride = 1, r = col;
    for (size_t k = 1; k < D.size(); k++) {
        const idx_t c = r % D[k];
        r /= D[k];
        m += ((D[k] - c) % D[k]) * stride;
        stride *= D[k];
    }
    return m;
}

// fftn of real data: an element whose mirror has the lower linear index becomes the mirror's conjugate, a
// self-conjugate element gets imaginary part 0. MATLAB's conjugate here is 0 - im: a zero imaginary part
// becomes +0 (the conjugate fill of a line in fft negates it instead).
template <typename R>
void hermitian_fill(R* x, const std::vector<idx_t>& D, int threads) {
    const idx_t n0 = D[0];
    idx_t ncol = 1;
    for (size_t k = 1; k < D.size(); k++) ncol *= D[k];
    const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(threads, ncol)));
    const R zero = 0;
#pragma omp parallel for num_threads(nt) schedule(dynamic, 16)
    for (idx_t col = 0; col < ncol; col++) {
        const idx_t mc = mirror_column(col, D);
        if (mc > col) continue;
        R* c = x + 2 * n0 * col;
        if (mc < col) {
            // every element: the mirror column, reversed (element 0 maps to element 0)
            const R* m = x + 2 * n0 * mc;
            c[0] = m[0];
            c[1] = zero - m[1];
            for (idx_t i = 1; i < n0; i++) {
                c[2 * i] = m[2 * (n0 - i)];
                c[2 * i + 1] = zero - m[2 * (n0 - i) + 1];
            }
            continue;
        }
        // the column is its own mirror: element i mirrors n0 - i
        c[1] = 0;
        for (idx_t i = 1; i < n0; i++) {
            const idx_t mi = n0 - i;
            if (mi < i) {
                c[2 * i] = c[2 * mi];
                c[2 * i + 1] = zero - c[2 * mi + 1];
            } else if (mi == i) {
                c[2 * i + 1] = 0;
            }
        }
    }
}

// element c (coordinates in the padded or truncated view D of the source S) or zero outside the source
template <typename R>
inline void view_value(const R* src, bool cplx, const std::vector<idx_t>& S, const std::vector<idx_t>& sstride,
        const std::vector<idx_t>& c, R& re, R& im) {
    idx_t off = 0;
    for (size_t k = 0; k < S.size(); k++) {
        if (c[k] >= S[k]) {
            re = 0;
            im = 0;
            return;
        }
        off += c[k] * sstride[k];
    }
    re = cplx ? src[2 * off] : src[off];
    im = cplx ? src[2 * off + 1] : R(0);
}

// x == conj(y), IEEE comparisons (NaN fails, -0 equals +0); real data: x == y
template <typename R, bool CPLX>
inline bool conj_equal(const R* x, const R* y) {
    if constexpr (CPLX) {
        return x[0] == y[0] && x[1] == -y[1];
    } else {
        return x[0] == y[0];
    }
}

// the test of is_hermitian when the view is the source itself, each pair of mirrors once
template <typename R, bool CPLX>
bool is_hermitian_same(const R* src, const std::vector<idx_t>& D, int d, int threads) {
    constexpr int c = CPLX ? 2 : 1;
    const int nd = static_cast<int>(D.size());
    std::atomic<bool> ok(true);
    if (d < 0) {
        const idx_t n0 = D[0];
        idx_t ncol = 1;
        for (int k = 1; k < nd; k++) ncol *= D[k];
        const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(threads, ncol)));
#pragma omp parallel for num_threads(nt) schedule(dynamic, 16)
        for (idx_t col = 0; col < ncol; col++) {
            if (!ok) continue;
            const idx_t mc = mirror_column(col, D);
            if (mc < col) continue;
            const R* a = src + c * n0 * col;
            const R* b = src + c * n0 * mc;
            bool good = conj_equal<R, CPLX>(a, b);
            if (mc > col) {
                for (idx_t i = 1; i < n0 && good; i++) good = conj_equal<R, CPLX>(a + c * i, b + c * (n0 - i));
            } else {
                for (idx_t i = 1; i <= n0 - i && good; i++) good = conj_equal<R, CPLX>(a + c * i, a + c * (n0 - i));
            }
            if (!good) ok = false;
        }
        return ok;
    }
    // lines along d: the block of L elements at index k against the one at n - k, for k <= n / 2
    idx_t L = 1, U = 1;
    for (int k = 0; k < d; k++) L *= D[k];
    for (int k = d + 1; k < nd; k++) U *= D[k];
    const idx_t n = D[d];
    const idx_t nk = n / 2 + 1;
    const idx_t ntask = U * nk;
    const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(threads, ntask)));
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t t = 0; t < ntask; t++) {
        if (!ok) continue;
        const idx_t ihi = t / nk, k = t % nk;
        const idx_t km = k == 0 ? 0 : n - k;
        const R* a = src + c * L * (k + n * ihi);
        const R* b = src + c * L * (km + n * ihi);
        bool good = true;
        for (idx_t i = 0; i < L && good; i++) good = conj_equal<R, CPLX>(a + c * i, b + c * i);
        if (!good) ok = false;
    }
    return ok;
}

// MATLAB's exact test for conjugate symmetry, of the source padded or truncated to D: along dimension d of
// every line (d >= 0), or N-D (d < 0): x(k) == conj(x(mirror(k))) for every element (IEEE comparisons)
template <typename R>
bool is_hermitian(const R* src, bool cplx, const std::vector<idx_t>& S, const std::vector<idx_t>& D, int d, int threads) {
    const int nd = static_cast<int>(D.size());
    idx_t total = 1;
    for (idx_t v : D) total *= v;
    if (total == 0) return true;
    if (S == D) {
        return cplx ? is_hermitian_same<R, true>(src, D, d, threads) : is_hermitian_same<R, false>(src, D, d, threads);
    }
    std::vector<idx_t> sstride(nd + 1, 1);
    for (int k = 0; k < nd; k++) sstride[k + 1] = sstride[k] * S[k];
    bool same = S == D;
    std::atomic<bool> ok(true);
    const idx_t n0 = D[0];
    const idx_t ncol = total / n0;
    const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(threads, ncol)));
#pragma omp parallel num_threads(nt)
    {
        std::vector<idx_t> c(nd), cm(nd);
#pragma omp for schedule(static)
        for (idx_t col = 0; col < ncol; col++) {
            if (!ok) continue;
            // coordinates of the column, and of the mirror column
            idx_t r = col;
            for (int k = 1; k < nd; k++) {
                c[k] = r % D[k];
                r /= D[k];
                cm[k] = (d < 0 || d == k) ? (D[k] - c[k]) % D[k] : c[k];
            }
            for (idx_t i = 0; i < n0 && ok; i++) {
                c[0] = i;
                cm[0] = (d < 0 || d == 0) ? (n0 - i) % n0 : i;
                R a, b, am, bm;
                if (same) {
                    idx_t off = 0, offm = 0;
                    for (int k = 0; k < nd; k++) {
                        off += c[k] * sstride[k];
                        offm += cm[k] * sstride[k];
                    }
                    a = cplx ? src[2 * off] : src[off];
                    b = cplx ? src[2 * off + 1] : R(0);
                    am = cplx ? src[2 * offm] : src[offm];
                    bm = cplx ? src[2 * offm + 1] : R(0);
                } else {
                    view_value(src, cplx, S, sstride, c, a, b);
                    view_value(src, cplx, S, sstride, cm, am, bm);
                }
                if (!(a == am && b == -bm)) ok = false;
            }
        }
    }
    return ok;
}

// a complex result whose imaginary parts are all zero is returned as real
template <typename R>
bool all_imag_zero(const R* x, idx_t numel, int threads) {
    std::atomic<bool> ok(true);
    const idx_t chunk = 1 << 16;
    const idx_t nchunks = (numel + chunk - 1) / chunk;
    const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(threads, nchunks)));
#pragma omp parallel for num_threads(nt) schedule(dynamic, 1)
    for (idx_t q = 0; q < nchunks; q++) {
        if (!ok) continue;
        const idx_t e = std::min(numel, (q + 1) * chunk);
        for (idx_t i = q * chunk; i < e; i++) {
            if (x[2 * i + 1] != 0) {
                ok = false;
                break;
            }
        }
    }
    return ok;
}

// --- arrays -----------------------------------------------------------------------------------------------

struct Ctx {
    Allocator alloc;
    void* ctx;
    int threads;
};

// Touches each page of a new buffer, each thread a contiguous part, so that the kernel maps its pages without
// the threads contending for the same page tables, as they do when a stage first writes a buffer in its
// scattered order (on a node whose memory has no huge pages left that costs several times the transform).
void prefault(void* p, size_t bytes, int threads) {
    constexpr size_t kPage = 4096;
    if (p == nullptr || threads < 2 || bytes < (size_t(8) << 20)) return;
    volatile char* c = static_cast<volatile char*>(p);
    const idx_t npages = static_cast<idx_t>((bytes + kPage - 1) / kPage);
#pragma omp parallel for num_threads(threads) schedule(static)
    for (idx_t i = 0; i < npages; i++) c[static_cast<size_t>(i) * kPage] = 0;
}

Array make(const Ctx& c, const std::vector<idx_t>& dims, bool is_double, bool is_complex) {
    Array a;
    a.dims = dims;
    a.is_double = is_double;
    a.is_complex = is_complex;
    a.data = c.alloc(c.ctx, dims, is_double, is_complex);
    if (a.data == nullptr && a.numel() > 0) throw std::bad_alloc();
    prefault(a.data, static_cast<size_t>(a.numel()) * (is_double ? 8 : 4) * (is_complex ? 2 : 1), c.threads);
    return a;
}

template <typename R>
Array realify(const Ctx& cx, const Array& y, int threads) {
    if (!y.is_complex || !all_imag_zero(static_cast<const R*>(y.data), y.numel(), threads)) return y;
    Array r = make(cx, y.dims, y.is_double, false);
    const R* s = static_cast<const R*>(y.data);
    R* t = static_cast<R*>(r.data);
    const idx_t n = y.numel();
#pragma omp parallel for num_threads(std::max(1, threads)) schedule(static)
    for (idx_t i = 0; i < n; i++) t[i] = s[2 * i];
    return r;
}

// x padded or truncated to D (zeros outside), same class
template <typename R>
Array resized(const Ctx& cx, const Array& x, const std::vector<idx_t>& D) {
    Array y = make(cx, D, x.is_double, x.is_complex);
    const int nd = static_cast<int>(D.size());
    const int c = x.is_complex ? 2 : 1;
    // the source's dimensions, those beyond its own of length 1
    std::vector<idx_t> S(x.dims);
    S.resize(D.size(), 1);
    std::vector<idx_t> sstride(nd + 1, 1);
    for (int k = 0; k < nd; k++) sstride[k + 1] = sstride[k] * S[k];
    std::vector<idx_t> co(nd, 0);
    const R* s = static_cast<const R*>(x.data);
    R* t = static_cast<R*>(y.data);
    const idx_t n = y.numel();
    for (idx_t i = 0; i < n; i++) {
        idx_t r = i, off = 0;
        bool in = true;
        for (int k = 0; k < nd; k++) {
            co[k] = r % D[k];
            r /= D[k];
            if (co[k] >= S[k]) in = false;
            off += co[k] * sstride[k];
        }
        for (int j = 0; j < c; j++) t[c * i + j] = in ? s[c * off + j] : R(0);
    }
    return y;
}

// the source (dims S, real or complex) padded with zeros or truncated to dims H, as complex values in dst
template <typename R>
void copy_view(const R* src, bool cplx, const std::vector<idx_t>& S, R* dst, const std::vector<idx_t>& H, int threads) {
    const int nd = static_cast<int>(H.size());
    std::vector<idx_t> sstride(nd + 1, 1);
    for (int k = 0; k < nd; k++) sstride[k + 1] = sstride[k] * S[k];
    const idx_t n0 = H[0];
    const idx_t m0 = std::min(H[0], S[0]);
    idx_t ncol = 1;
    for (int k = 1; k < nd; k++) ncol *= H[k];
    const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(threads, ncol)));
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t col = 0; col < ncol; col++) {
        R* d = dst + 2 * n0 * col;
        idx_t r = col, off = 0;
        bool in = true;
        for (int k = 1; k < nd; k++) {
            const idx_t c = r % H[k];
            r /= H[k];
            if (c >= S[k]) in = false;
            off += c * sstride[k];
        }
        idx_t i = 0;
        if (in) {
            for (; i < m0; i++) {
                d[2 * i] = cplx ? src[2 * (off + i)] : src[off + i];
                d[2 * i + 1] = cplx ? src[2 * (off + i) + 1] : R(0);
            }
        }
        for (; i < n0; i++) {
            d[2 * i] = 0;
            d[2 * i + 1] = 0;
        }
    }
}

std::vector<idx_t> dims_to(const std::vector<idx_t>& d, size_t nd) {
    std::vector<idx_t> r(d);
    while (r.size() < nd) r.push_back(1);
    return r;
}

template <typename R>
Array fft_t(const Array& x, idx_t n, int dim, const Options& o, const Ctx& cx) {
    const size_t nd = std::max<size_t>(x.dims.size(), static_cast<size_t>(dim) + 1);
    const std::vector<idx_t> S = dims_to(x.dims, nd);
    std::vector<idx_t> D = S;
    if (n >= 0) D[dim] = n;
    const idx_t len = D[dim];
    // MATLAB returns the data unchanged along a dimension of length 1 (no transform)
    if (S[dim] == 1 && len == 1) return resized<R>(cx, x, D);
    Array y = make(cx, D, x.is_double, true);
    if (y.numel() == 0) return realify<R>(cx, y, o.threads);
    Stage<R> st{x.is_double ? 1 : 0, x.is_complex ? kC2CF : kR2C, static_cast<const R*>(x.data), x.is_complex, S,
        static_cast<R*>(y.data), D, dim, len, 0, false};
    idx_t U = 1;
    for (size_t k = dim + 1; k < nd; k++) U *= D[k];
    st.Ufull = U;
    run_stage(st, o);
    return realify<R>(cx, y, o.threads);
}

template <typename R>
Array ifft_t(const Array& x, idx_t n, int dim, bool symmetric, const Options& o, const Ctx& cx) {
    const size_t nd = std::max<size_t>(x.dims.size(), static_cast<size_t>(dim) + 1);
    const std::vector<idx_t> S = dims_to(x.dims, nd);
    std::vector<idx_t> D = S;
    if (n >= 0) D[dim] = n;
    const idx_t len = D[dim];
    idx_t U = 1;
    for (size_t k = dim + 1; k < nd; k++) U *= D[k];
    if (S[dim] == 1 && len == 1 && !symmetric) return resized<R>(cx, x, D);
    const bool herm = symmetric
        || is_hermitian(static_cast<const R*>(x.data), x.is_complex, S, D, dim, o.threads);
    if (herm) {
        Array y = make(cx, D, x.is_double, false);
        if (y.numel() == 0) return y;
        Stage<R> st{x.is_double ? 1 : 0, kC2R, static_cast<const R*>(x.data), x.is_complex, S, static_cast<R*>(y.data),
            D, dim, len, U, true};
        run_stage(st, o);
        return y;
    }
    Array y = make(cx, D, x.is_double, true);
    if (y.numel() == 0) return realify<R>(cx, y, o.threads);
    Stage<R> st{x.is_double ? 1 : 0, kC2CB, static_cast<const R*>(x.data), x.is_complex, S, static_cast<R*>(y.data), D,
        dim, len, U, true};
    run_stage(st, o);
    return realify<R>(cx, y, o.threads);
}

template <typename R>
Array fftn_t(const Array& x, const std::vector<idx_t>& siz, const Options& o, const Ctx& cx) {
    const size_t nd = std::max(x.dims.size(), siz.size());
    const std::vector<idx_t> S = dims_to(x.dims, nd);
    const std::vector<idx_t> D = siz.empty() ? S : dims_to(siz, nd);
    std::vector<int> tdims;
    for (size_t k = 0; k < nd; k++) {
        if (D[k] > 1) tdims.push_back(static_cast<int>(k));
    }
    idx_t total = 1;
    for (idx_t v : D) total *= v;
    if (tdims.empty() || total == 0) {
        Array r = resized<R>(cx, x, D);
        return realify<R>(cx, r, o.threads);
    }
    Array y = make(cx, D, x.is_double, true);
    const int prec = x.is_double ? 1 : 0;
    for (size_t i = 0; i < tdims.size(); i++) {
        const int d = tdims[i];
        idx_t U = 1;
        for (size_t k = d + 1; k < nd; k++) U *= D[k];
        Stage<R> st;
        st.prec = prec;
        st.d = d;
        st.n = D[d];
        st.dst = static_cast<R*>(y.data);
        st.D = D;
        st.Ufull = U;
        st.scale = false;
        // real input: the Hermitian fill below writes the elements past n / 2 along the last dimension
        if (!x.is_complex && i + 1 == tdims.size()) st.kwrite = D[d] / 2 + 1;
        if (i == 0) {
            st.kind = x.is_complex ? kC2CF : kR2C;
            st.src = static_cast<const R*>(x.data);
            st.src_complex = x.is_complex;
            st.S = S;
        } else {
            st.kind = kC2CF;
            st.src = static_cast<const R*>(y.data);
            st.src_complex = true;
            st.S = D;
        }
        run_stage(st, o);
    }
    if (!x.is_complex) hermitian_fill(static_cast<R*>(y.data), D, o.threads);
    return realify<R>(cx, y, o.threads);
}

template <typename R>
Array ifftn_t(const Array& x, const std::vector<idx_t>& siz, bool symmetric, const Options& o, const Ctx& cx) {
    const size_t nd = std::max({x.dims.size(), siz.size(), size_t(2)});
    const std::vector<idx_t> S = dims_to(x.dims, nd);
    const std::vector<idx_t> D = siz.empty() ? S : dims_to(siz, nd);
    std::vector<int> tdims;
    for (size_t k = 0; k < nd; k++) {
        if (D[k] > 1) tdims.push_back(static_cast<int>(k));
    }
    idx_t total = 1;
    for (idx_t v : D) total *= v;
    const int prec = x.is_double ? 1 : 0;
    if (tdims.empty() || total == 0) {
        // nothing to transform; 'symmetric' keeps the real part
        Array r = resized<R>(cx, x, D);
        if (symmetric && r.is_complex) {
            Array t = make(cx, D, x.is_double, false);
            const R* s = static_cast<const R*>(r.data);
            R* q = static_cast<R*>(t.data);
            for (idx_t i = 0; i < total; i++) q[i] = s[2 * i];
            return t;
        }
        return realify<R>(cx, r, o.threads);
    }
    const bool herm = symmetric || is_hermitian(static_cast<const R*>(x.data), x.is_complex, S, D, -1, o.threads);
    if (herm) {
        // complex stages along the transformed dimensions but MATLAB's last one, on the part of the array the
        // last stage reads (the first n / 2 + 1 elements along it), then complex-to-real along it. MATLAB's
        // last dimension is that of its size (2 dimensions at least): for a column vector the second, of
        // length 1, whose complex-to-real transform keeps the real part.
        int last = 1;
        for (int d : tdims) last = std::max(last, d);
        const idx_t nl = D[last];
        std::vector<int> cdims;
        for (int d : tdims) {
            if (d != last) cdims.push_back(d);
        }
        Array y = make(cx, D, x.is_double, false);
        if (cdims.empty() && !symmetric) {
            Stage<R> st{prec, kC2R, static_cast<const R*>(x.data), x.is_complex, S, static_cast<R*>(y.data), D, last,
                nl, 1, true};
            run_stage(st, o);
            return y;
        }
        std::vector<idx_t> H = D;
        H[last] = nl / 2 + 1;
        Array h = make(cx, H, x.is_double, true);
        // 'symmetric': MATLAB first makes the padded or truncated data conjugate symmetric with fftn's fill; in
        // the part read here that changes only the planes k = 0 and n / 2 along the last dimension
        const R* first = static_cast<const R*>(x.data);
        bool first_complex = x.is_complex;
        std::vector<idx_t> firstS = S;
        if (symmetric) {
            copy_view(static_cast<const R*>(x.data), x.is_complex, S, static_cast<R*>(h.data), H, o.threads);
            idx_t P = 1;
            for (int k = 0; k < last; k++) P *= D[k];
            const std::vector<idx_t> Dp(D.begin(), D.begin() + last);
            const int nplanes = (nl % 2 == 0 && nl / 2 > 0) ? 2 : 1;
            for (int q = 0; q < nplanes; q++) {
                R* plane = static_cast<R*>(h.data) + 2 * P * (q == 0 ? 0 : nl / 2);
                hermitian_fill(plane, Dp, o.threads);
            }
            first = static_cast<const R*>(h.data);
            first_complex = true;
            firstS = H;
        }
        for (size_t i = 0; i < cdims.size(); i++) {
            const int d = cdims[i];
            idx_t U = 1;
            for (size_t k = d + 1; k < nd; k++) U *= D[k];
            Stage<R> st;
            st.prec = prec;
            st.kind = kC2CB;
            st.d = d;
            st.n = D[d];
            st.dst = static_cast<R*>(h.data);
            st.D = H;
            st.Ufull = U;
            st.scale = true;
            if (i == 0) {
                st.src = first;
                st.src_complex = first_complex;
                st.S = firstS;
            } else {
                st.src = static_cast<const R*>(h.data);
                st.src_complex = true;
                st.S = H;
            }
            run_stage(st, o);
        }
        if (nl > 1) {
            Stage<R> st{prec, kC2R, static_cast<const R*>(h.data), true, H, static_cast<R*>(y.data), D, last, nl, 1,
                true};
            run_stage(st, o);
        } else {
            const R* q = static_cast<const R*>(h.data);
            R* r = static_cast<R*>(y.data);
#pragma omp parallel for num_threads(std::max(1, o.threads)) schedule(static)
            for (idx_t i = 0; i < total; i++) r[i] = q[2 * i];
        }
        return y;
    }
    Array y = make(cx, D, x.is_double, true);
    for (size_t i = 0; i < tdims.size(); i++) {
        const int d = tdims[i];
        idx_t U = 1;
        for (size_t k = d + 1; k < nd; k++) U *= D[k];
        Stage<R> st;
        st.prec = prec;
        st.kind = kC2CB;
        st.d = d;
        st.n = D[d];
        st.dst = static_cast<R*>(y.data);
        st.D = D;
        st.Ufull = U;
        st.scale = true;
        if (i == 0) {
            st.src = static_cast<const R*>(x.data);
            st.src_complex = x.is_complex;
            st.S = S;
        } else {
            st.src = static_cast<const R*>(y.data);
            st.src_complex = true;
            st.S = D;
        }
        run_stage(st, o);
    }
    return realify<R>(cx, y, o.threads);
}

}  // namespace

bool line_alone(idx_t j, idx_t N, idx_t n, idx_t seg, idx_t grain) {
    return alone_g(j, N, stage_split(grain, n, N) ? grain : 0, seg);
}

Array fft(const Array& x, idx_t n, int dim, const Options& o, Allocator alloc, void* ctx) {
    const Ctx c{alloc, ctx, o.threads};
    return x.is_double ? fft_t<double>(x, n, dim, o, c) : fft_t<float>(x, n, dim, o, c);
}

Array ifft(const Array& x, idx_t n, int dim, bool symmetric, const Options& o, Allocator alloc, void* ctx) {
    const Ctx c{alloc, ctx, o.threads};
    return x.is_double ? ifft_t<double>(x, n, dim, symmetric, o, c) : ifft_t<float>(x, n, dim, symmetric, o, c);
}

Array fftn(const Array& x, const std::vector<idx_t>& siz, const Options& o, Allocator alloc, void* ctx) {
    const Ctx c{alloc, ctx, o.threads};
    return x.is_double ? fftn_t<double>(x, siz, o, c) : fftn_t<float>(x, siz, o, c);
}

Array ifftn(const Array& x, const std::vector<idx_t>& siz, bool symmetric, const Options& o, Allocator alloc, void* ctx) {
    const Ctx c{alloc, ctx, o.threads};
    return x.is_double ? ifftn_t<double>(x, siz, symmetric, o, c) : ifftn_t<float>(x, siz, symmetric, o, c);
}

std::string plan_string(bool is_double, const std::string& kind, idx_t n, int howmany) {
    int k;
    if (kind == "r2c") {
        k = kR2C;
    } else if (kind == "c2c-1") {
        k = kC2CF;
    } else if (kind == "c2c+1") {
        k = kC2CB;
    } else if (kind == "c2r") {
        k = kC2R;
    } else {
        throw std::invalid_argument("kind must be r2c, c2c-1, c2c+1 or c2r");
    }
    const int prec = is_double ? 1 : 0;
    void* p = get_plan(prec, k, n, howmany);
    std::lock_guard<std::mutex> g(planner_lock());
    const FftwApi& a = backend().api[prec];
    char* s = a.sprint(p);
    std::string r = s != nullptr ? s : "";
    if (s != nullptr) a.free_string(s);
    for (char& ch : r) {
        if (ch == '\n') ch = ' ';
    }
    return r;
}

void clear_plans() {
    std::lock_guard<std::mutex> g(planner_lock());
    for (auto& kv : plan_cache()) {
        backend().api[std::get<0>(kv.first)].destroy(kv.second);
    }
    plan_cache().clear();
}

}  // namespace pk5d::fft
