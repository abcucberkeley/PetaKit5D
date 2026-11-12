// MATLAB R2024b's fft drivers on strided arrays (see driver.h for the rules this reproduces).
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

using dims_t = std::vector<idx_t>;

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

using PlanKey = std::tuple<int, int, int, idx_t, int>;  // backend, precision, kind, n, lines per execute

std::map<PlanKey, void*>& plan_cache() {
    static std::map<PlanKey, void*> m;
    return m;
}

void dists(int kind, idx_t n, idx_t& idist, idx_t& odist) {
    const idx_t nh = n / 2 + 1;
    idist = kind == kC2R ? nh : n;
    odist = kind == kR2C ? nh : n;
}

void* get_plan(int be, int prec, int kind, idx_t n, int hm) {
    std::lock_guard<std::mutex> g(planner_lock());
    const PlanKey key(be, prec, kind, n, hm);
    auto it = plan_cache().find(key);
    if (it != plan_cache().end()) return it->second;
    const FftwApi& a = backend(be).api[prec];
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
    if (p == nullptr) throw std::runtime_error("FFT planning failed");
    plan_cache()[key] = p;
    return p;
}

inline void execute(const FftwApi& a, int kind, void* plan, void* in, void* out) {
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
    int be = 0;
    int prec = 0;
    int kind = kC2CF;
    // source: real or complex, MATLAB's dimensions S and element strides ss; padded with zeros or truncated to
    // the destination's extent in every dimension (a dimension beyond the source's own has length 1)
    const R* src = nullptr;
    bool src_complex = false;
    dims_t S, ss;
    // destination: dimensions D and strides ds; along d, kw <= D[d] elements of the n of a line are written
    R* dst = nullptr;
    dims_t D, ds;
    // MATLAB's dimensions for its numbering of the lines (D, or more along the dimensions after d when the
    // destination holds a prefix of MATLAB's array there)
    dims_t Dm;
    int d = 0;
    idx_t n = 0;
    bool scale = false;  // divide by n (inverse stages)
    idx_t kwrite = -1;   // outputs written per line (the first ones; -1: n)
};

// per-thread scratch of a stage
template <typename R>
struct Scratch {
    Buf tin, tout, in1, out1;
    std::vector<idx_t> soff, doff, slot;
    std::vector<int> alone;
    std::vector<R*> ein;         // a line's input in the scratch
    std::vector<const R*> eout;  // its output
};

// the geometry of a stage, shared by its tiles: the lines' other dimensions in the destination's memory
// order (the fastest one, f, inside a tile), MATLAB's numbering of the lines
struct Geom {
    int nd = 0, d = 0, f = -1;
    idx_t n = 0, nh = 0, Nfull = 0, g = 0, seg = 0, idist = 0, odist = 0, need = 0, avail = 0, kw = 0;
    idx_t T = 1, nlines = 1, ntiles = 1;
    std::vector<int> outer;  // the other dimensions but f, fastest first
    dims_t mstride;          // MATLAB's line index stride of each dimension (0 for d)
    bool gather_by_line = true, scatter_by_line = true;
};

template <typename R>
inline void zero_line(R* e, idx_t count) {
    std::memset(e, 0, sizeof(R) * count);
}

// one tile: T consecutive lines in the destination's memory order (f fastest, then the outer dimensions):
// gather (zero-padded), MATLAB's groups of four and lines alone, scatter. KIND and the source's complexity
// are compile-time, so the copy loops have no branches.
template <typename R, int KIND, bool SRCC>
void run_tile(const Stage<R>& st, const Geom& G, Scratch<R>& w, const FftwApi& api, void* p4, void* p1, idx_t q) {
    constexpr int ci = KIND == kR2C ? 1 : 2;   // R values per scratch input element
    constexpr int co = KIND == kC2R ? 1 : 2;   // per scratch output element
    constexpr int cs = SRCC ? 2 : 1;           // per source element
    constexpr int cd = KIND == kC2R ? 1 : 2;   // per destination element
    const idx_t n = G.n, nh = G.nh, idist = G.idist, odist = G.odist;
    const int d = G.d, f = G.f;
    // the tile's lines: offsets in the source and the destination, MATLAB's line index, alone or grouped
    const idx_t l0 = q * G.T;
    const idx_t Tq = std::min(G.T, G.nlines - l0);
    idx_t ng = 0;
    for (idx_t t = 0; t < Tq; t++) {
        idx_t r = l0 + t, so = 0, dof = 0, j = 0;
        bool in = true;
        if (f >= 0) {
            const idx_t c = r % st.D[f];
            r /= st.D[f];
            if (c >= st.S[f]) in = false;
            so += c * st.ss[f];
            dof += c * st.ds[f];
            j += c * G.mstride[f];
        }
        for (int k : G.outer) {
            const idx_t c = r % st.D[k];
            r /= st.D[k];
            if (c >= st.S[k]) in = false;
            so += c * st.ss[k];
            dof += c * st.ds[k];
            j += c * G.mstride[k];
        }
        w.soff[t] = in ? so : -1;
        w.doff[t] = dof;
        w.alone[t] = alone_g(j, G.Nfull, G.g, G.seg) ? 1 : 0;
        if (!w.alone[t]) ng++;
    }
    const idx_t g4 = (ng + kB - 1) / kB * kB;
    R* const tI = reinterpret_cast<R*>(w.tin.p);
    R* const tO = reinterpret_cast<R*>(w.tout.p);
    {
        idx_t sg = 0, sa = g4;
        for (idx_t t = 0; t < Tq; t++) {
            w.slot[t] = w.alone[t] ? sa++ : sg++;
            w.ein[t] = tI + w.slot[t] * idist * ci;
            w.eout[t] = tO + w.slot[t] * odist * co;
        }
    }
    // the tile's lines are neighbours along f in the source and in the destination (the usual case: a tile
    // that does not wrap to the next row of f and has no padding lines)
    bool adj = f >= 0 && w.soff[0] >= 0;
    if (adj) {
        const idx_t sf = st.ss[f], df = st.ds[f];
        for (idx_t t = 1; t < Tq && adj; t++) {
            adj = w.soff[t] == w.soff[0] + t * sf && w.doff[t] == w.doff[0] + t * df;
        }
    }
    R* const* const ein = w.ein.data();
    const R* const* const eout = w.eout.data();
    const idx_t avail = G.avail, need = G.need, sd = st.ss[d];
    // gather into the scratch
    if (G.gather_by_line) {
        for (idx_t t = 0; t < Tq; t++) {
            R* e = ein[t];
            if (w.soff[t] < 0) {
                zero_line(e, avail * ci);
                continue;
            }
            const R* s = st.src + w.soff[t] * cs;
            if (ci == cs && sd == 1) {
                std::memcpy(e, s, sizeof(R) * avail * ci);
            } else {
                for (idx_t k = 0; k < avail; k++) {
                    e[k * ci] = s[sd * k * cs];
                    if constexpr (ci == 2) e[k * ci + 1] = SRCC ? s[sd * k * cs + 1] : R(0);
                }
            }
        }
    } else if (adj) {
        // element k of line t at s0[(sd * k + sf * t) * cs]: each row of the tile read in one piece
        const R* const s0 = st.src + w.soff[0] * cs;
        const idx_t sk = sd * cs, sf = st.ss[f] * cs;
        for (idx_t k = 0; k < avail; k++) {
            const R* s = s0 + sk * k;
            if constexpr (ci == 2) {
                for (idx_t t = 0; t < Tq; t++) {
                    R* e = ein[t] + 2 * k;
                    e[0] = s[sf * t];
                    e[1] = SRCC ? s[sf * t + 1] : R(0);
                }
            } else {
                for (idx_t t = 0; t < Tq; t++) ein[t][k] = s[sf * t];
            }
        }
    } else {
        // lines adjacent in memory: element k of line t at src[soff[t] + sd * k]
        for (idx_t t = 0; t < Tq; t++) {
            if (w.soff[t] < 0) zero_line(ein[t], avail * ci);
        }
        for (idx_t k = 0; k < avail; k++) {
            for (idx_t t = 0; t < Tq; t++) {
                if (w.soff[t] < 0) continue;
                const R* s = st.src + (w.soff[t] + sd * k) * cs;
                R* e = ein[t] + k * ci;
                e[0] = s[0];
                if constexpr (ci == 2) e[1] = SRCC ? s[1] : R(0);
            }
        }
    }
    if (avail < need) {
        for (idx_t t = 0; t < Tq; t++) zero_line(tI + (w.slot[t] * idist + avail) * ci, (need - avail) * ci);
    }
    // a group of fewer than four lines is filled up with copies of its last line
    for (idx_t s = ng; s < g4; s++) {
        std::memcpy(tI + s * idist * ci, tI + (ng - 1) * idist * ci, sizeof(R) * idist * ci);
    }
    for (idx_t s = 0; s < g4; s += kB) {
        execute(api, KIND, p4, tI + s * idist * ci, tO + s * odist * co);
    }
    for (idx_t s = g4; s < g4 + (Tq - ng); s++) {
        std::memcpy(w.in1.p, tI + s * idist * ci, sizeof(R) * idist * ci);
        execute(api, KIND, p1, w.in1.p, w.out1.p);
        std::memcpy(tO + s * odist * co, w.out1.p, sizeof(R) * odist * co);
    }
    // scatter: element k of line t at dst[doff[t] + ds[d] * k], the first kw elements of a line. A real-to-
    // complex line is filled to its full length by conjugation (the imaginary part negated); an inverse stage
    // divides by n.
    const idx_t kw = G.kw;
    const idx_t dd = st.ds[d];
    const R fn = static_cast<R>(n);
    const bool scale = st.scale;
    auto put = [&](idx_t t, idx_t k) {
        R* e = st.dst + (w.doff[t] + dd * k) * cd;
        if constexpr (KIND == kR2C) {
            const R* lo = eout[t];
            if (k < nh) {
                e[0] = lo[2 * k];
                e[1] = lo[2 * k + 1];
            } else {
                e[0] = lo[2 * (n - k)];
                e[1] = -lo[2 * (n - k) + 1];
            }
        } else if constexpr (KIND == kC2R) {
            const R v = eout[t][k];
            e[0] = scale ? v / fn : v;
        } else {
            const R* lo = eout[t] + 2 * k;
            if (scale) {
                e[0] = lo[0] / fn;
                e[1] = lo[1] / fn;
            } else {
                e[0] = lo[0];
                e[1] = lo[1];
            }
        }
    };
    if (G.scatter_by_line && dd == 1) {
        // contiguous destination lines
        for (idx_t t = 0; t < Tq; t++) {
            const R* lo = eout[t];
            R* e = st.dst + w.doff[t] * cd;
            if constexpr (KIND == kR2C) {
                std::memcpy(e, lo, sizeof(R) * 2 * std::min(nh, kw));
                for (idx_t k = nh; k < kw; k++) {
                    e[2 * k] = lo[2 * (n - k)];
                    e[2 * k + 1] = -lo[2 * (n - k) + 1];
                }
            } else if (scale) {
                const idx_t m = kw * cd;
                for (idx_t j = 0; j < m; j++) e[j] = lo[j] / fn;
            } else {
                std::memcpy(e, lo, sizeof(R) * kw * cd);
            }
        }
    } else if (G.scatter_by_line) {
        for (idx_t t = 0; t < Tq; t++) {
            for (idx_t k = 0; k < kw; k++) put(t, k);
        }
    } else if (adj) {
        // element k of line t at d0[(dd * k + df * t) * cd]: each row of the tile written in one piece
        R* const d0 = st.dst + w.doff[0] * cd;
        const idx_t dk = dd * cd, df = st.ds[f] * cd;
        for (idx_t k = 0; k < kw; k++) {
            R* row = d0 + dk * k;
            if constexpr (KIND == kR2C) {
                if (k < nh) {
                    for (idx_t t = 0; t < Tq; t++) {
                        row[df * t] = eout[t][2 * k];
                        row[df * t + 1] = eout[t][2 * k + 1];
                    }
                } else {
                    const idx_t kk = n - k;
                    for (idx_t t = 0; t < Tq; t++) {
                        row[df * t] = eout[t][2 * kk];
                        row[df * t + 1] = -eout[t][2 * kk + 1];
                    }
                }
            } else if constexpr (KIND == kC2R) {
                if (scale) {
                    for (idx_t t = 0; t < Tq; t++) row[df * t] = eout[t][k] / fn;
                } else {
                    for (idx_t t = 0; t < Tq; t++) row[df * t] = eout[t][k];
                }
            } else {
                if (scale) {
                    for (idx_t t = 0; t < Tq; t++) {
                        row[df * t] = eout[t][2 * k] / fn;
                        row[df * t + 1] = eout[t][2 * k + 1] / fn;
                    }
                } else {
                    for (idx_t t = 0; t < Tq; t++) {
                        row[df * t] = eout[t][2 * k];
                        row[df * t + 1] = eout[t][2 * k + 1];
                    }
                }
            }
        }
    } else {
        for (idx_t k = 0; k < kw; k++) {
            for (idx_t t = 0; t < Tq; t++) put(t, k);
        }
    }
}

template <typename R, int KIND, bool SRCC>
void run_stage_k(const Stage<R>& st, const Geom& G, const Options& o) {
    constexpr int ci = KIND == kR2C ? 1 : 2;
    constexpr int co = KIND == kC2R ? 1 : 2;
    void* p4 = get_plan(st.be, st.prec, KIND, G.n, kB);
    void* p1 = get_plan(st.be, st.prec, KIND, G.n, 1);
    const FftwApi& api = backend(st.be).api[st.prec];
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
            w.soff.resize(static_cast<size_t>(G.T));
            w.doff.resize(static_cast<size_t>(G.T));
            w.slot.resize(static_cast<size_t>(G.T));
            w.alone.resize(static_cast<size_t>(G.T));
            w.ein.resize(static_cast<size_t>(G.T));
            w.eout.resize(static_cast<size_t>(G.T));
        } catch (...) {
            failed = true;
        }
#pragma omp for schedule(dynamic, 1)
        for (idx_t q = 0; q < G.ntiles; q++) {
            if (failed) continue;
            run_tile<R, KIND, SRCC>(st, G, w, api, p4, p1, q);
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
    idx_t lines = 1;
    for (int k = 0; k < G.nd; k++) {
        if (k != G.d) lines *= st.D[k];
    }
    if (lines == 0 || G.n == 0) return;
    // MATLAB's numbering: line index j = sum of c[k] * mstride[k] over the other dimensions, in MATLAB's
    // order (those before d first), with MATLAB's sizes Dm
    G.mstride.assign(G.nd, 0);
    idx_t L = 1;
    for (int k = 0; k < G.d; k++) {
        G.mstride[k] = L;
        L *= st.Dm[k];
    }
    idx_t m = L;
    for (int k = G.d + 1; k < G.nd; k++) {
        G.mstride[k] = m;
        m *= st.Dm[k];
    }
    G.Nfull = m;
    G.g = stage_split(o.grain, G.n, G.Nfull) ? o.grain : 0;
    G.seg = L > 1 ? L : 0;
    dists(st.kind, G.n, G.idist, G.odist);
    // elements of the source line that are read (the rest is zero)
    G.need = st.kind == kC2R ? G.nh : G.n;
    G.avail = std::min(st.S[G.d], G.need);
    G.kw = st.kwrite < 0 ? G.n : std::min(st.kwrite, G.n);
    // the other dimensions of length > 1 in the destination's memory order; f, the fastest, inside a tile
    std::vector<int> od;
    for (int k = 0; k < G.nd; k++) {
        if (k != G.d && st.D[k] > 1) od.push_back(k);
    }
    std::stable_sort(od.begin(), od.end(), [&](int a, int b) { return st.ds[a] < st.ds[b]; });
    G.f = od.empty() ? -1 : od[0];
    G.outer.assign(od.size() > 1 ? od.begin() + 1 : od.end(), od.end());
    G.nlines = lines;
    // copy loops: line by line when a line's elements are closer in memory than neighbouring lines
    G.gather_by_line = G.f < 0 || st.ss[G.d] <= st.ss[G.f];
    G.scatter_by_line = G.f < 0 || st.ds[G.d] <= st.ds[G.f];
    // lines per tile (the tiles do not change which lines MATLAB groups); a copy across lines takes them 16 at a
    // time or more, a cache line or more of each row, while the tile's scratch stays in the core's cache
    G.T = std::max<idx_t>(kB, std::min<idx_t>(64, ((16384 / std::max<idx_t>(G.n, 1)) + 3) / 4 * 4));
    if ((!G.gather_by_line || !G.scatter_by_line) && G.n <= 4096) G.T = std::min<idx_t>(64, (G.T + 15) / 16 * 16);
    G.T = std::min(G.T, G.nlines);
    G.ntiles = (G.nlines + G.T - 1) / G.T;
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

// --- element access of strided arrays ------------------------------------------------------------------

// planes: MATLAB's dimensions 3, 4, ... (index p in MATLAB's order); within a plane rows along dimension 1
// (y) and columns along dimension 2 (x)
struct Planes {
    int nd = 0;
    idx_t n0 = 1, n1 = 1, P = 1;
    idx_t s0 = 0, s1 = 0;
    dims_t D, st;
    explicit Planes(const dims_t& dims, const dims_t& strides) : D(dims), st(strides) {
        nd = static_cast<int>(D.size());
        n0 = nd > 0 ? D[0] : 1;
        n1 = nd > 1 ? D[1] : 1;
        s0 = nd > 0 ? st[0] : 0;
        s1 = nd > 1 ? st[1] : 0;
        for (int k = 2; k < nd; k++) P *= D[k];
    }
    // offset of plane p and, mirror true, of its N-D mirror plane (index pm)
    idx_t offset(idx_t p, bool mirror, idx_t* pm) const {
        idx_t off = 0, mp = 0, mul = 1;
        for (int k = 2; k < nd; k++) {
            const idx_t c = p % D[k];
            p /= D[k];
            const idx_t cm = mirror ? (D[k] - c) % D[k] : c;
            off += cm * st[k];
            mp += cm * mul;
            mul *= D[k];
        }
        if (pm != nullptr) *pm = mp;
        return off;
    }
};

// fftn of real data: an element whose mirror has the lower MATLAB linear index becomes the mirror's
// conjugate, a self-conjugate element gets imaginary part 0. MATLAB's conjugate here is 0 - im: a zero
// imaginary part becomes +0 (the conjugate fill of a line in fft negates it instead).
// The elements it reads are never written (an element whose mirror has the higher index keeps its value),
// so the rows (dimension 2, x, in a row of the native layout) are filled in any order, in parallel.
template <typename R>
void hermitian_fill(R* x, const dims_t& D, const dims_t& st, int threads) {
    const Planes pl(D, st);
    const idx_t n0 = pl.n0, n1 = pl.n1, s0 = pl.s0, s1 = pl.s1;
    const idx_t rows = pl.P * n0;
    const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(threads, rows * n1 / 16384 + 1)));
    const R zero = 0;
#pragma omp parallel for num_threads(nt) schedule(dynamic, 8)
    for (idx_t r = 0; r < rows; r++) {
        const idx_t p = r / n0, y = r % n0;
        idx_t pm = 0;
        const idx_t mo = pl.offset(p, true, &pm);
        if (pm > p) continue;
        const idx_t my = (n0 - y) % n0;
        R* c = x + 2 * (pl.offset(p, false, nullptr) + y * s0);
        const R* m = x + 2 * (mo + my * s0);
        // element i of the row takes the conjugate of element mi of the mirror row (i = 0 mirrors 0, i >= 1
        // mirrors n1 - i)
        auto take = [&](idx_t i, idx_t mi) {
            c[2 * i * s1] = m[2 * mi * s1];
            c[2 * i * s1 + 1] = zero - m[2 * mi * s1 + 1];
        };
        if (pm < p) {
            // the mirror plane comes first in MATLAB's order: every element is its mirror's conjugate
            take(0, 0);
            for (idx_t i = 1; i < n1; i++) take(i, n1 - i);
            continue;
        }
        // the plane is its own mirror: within it MATLAB's index compares x first, then y. Elements past n1 / 2
        // take their mirror's conjugate; those that are their own mirror along x (0, and n1 / 2 for even n1)
        // take the mirror row's when it comes first, or are self-conjugate (imaginary part 0)
        auto self = [&](idx_t i) {
            if (my < y) {
                take(i, i);
            } else if (my == y) {
                c[2 * i * s1 + 1] = 0;
            }
        };
        self(0);
        if (n1 % 2 == 0 && n1 > 1) self(n1 / 2);
        for (idx_t i = n1 / 2 + 1; i < n1; i++) take(i, n1 - i);
    }
}

// MATLAB's exact test for conjugate symmetry of the source (dims S, strides ss) padded or truncated to D:
// along dimension d of every line (d >= 0), or N-D (d < 0): x(k) == conj(x(mirror(k))) for every element
// (IEEE comparisons: NaN fails, -0 equals +0); real data: x(k) == x(mirror(k))
template <typename R>
bool is_hermitian(const R* src, bool cplx, const dims_t& S, const dims_t& ss, const dims_t& D, int d, int threads) {
    const int nd = static_cast<int>(D.size());
    idx_t total = 1;
    for (idx_t v : D) total *= v;
    if (total == 0) return true;
    const Planes pl(D, dims_t(nd, 0));
    const idx_t n0 = pl.n0, n1 = pl.n1;
    const bool m0 = d < 0 || d == 0, m1 = d < 0 || d == 1;
    bool same = true;  // no padding or truncation: every element is the source's
    for (int k = 0; k < nd; k++) same = same && S[k] == D[k];
    std::atomic<bool> ok(true);
    const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(threads, pl.P * n0 * n1 / 16384 + 1)));
#pragma omp parallel num_threads(nt)
    {
        dims_t c(nd, 0), cm(nd, 0);
#pragma omp for schedule(dynamic, 16)
        for (idx_t t = 0; t < pl.P * n0; t++) {
            if (!ok) continue;
            const idx_t p = t / n0, y = t % n0;
            idx_t r = p;
            bool plane_in = true, mplane_in = true;
            idx_t po = 0, pmo = 0;
            for (int k = 2; k < nd; k++) {
                c[k] = r % D[k];
                r /= D[k];
                cm[k] = (d < 0 || d == k) ? (D[k] - c[k]) % D[k] : c[k];
                if (c[k] >= S[k]) plane_in = false;
                if (cm[k] >= S[k]) mplane_in = false;
                po += c[k] * ss[k];
                pmo += cm[k] * ss[k];
            }
            const idx_t my = m0 ? (n0 - y) % n0 : y;
            if (same) {
                // a row against its mirror row, read backwards along x; each pair of rows once (the test is
                // symmetric: a == am and b == -bm hold for the mirror as well), half of a row that is its own
                // mirror
                idx_t pm = 0;
                for (int k = nd - 1; k >= 2; k--) pm = pm * D[k] + cm[k];
                if (pm < p || (pm == p && my < y)) continue;
                const idx_t s1 = nd > 1 ? ss[1] : 0;
                const idx_t iend = (pm == p && my == y && m1) ? n1 / 2 + 1 : n1;
                bool good = true;
                if (cplx) {
                    const R* a0 = src + 2 * (po + y * ss[0]);
                    const R* b0 = src + 2 * (pmo + my * ss[0]);
                    for (idx_t i = 0; i < iend; i++) {
                        const idx_t mi = m1 ? (i == 0 ? 0 : n1 - i) : i;
                        good &= (a0[2 * i * s1] == b0[2 * mi * s1]) & (a0[2 * i * s1 + 1] == -b0[2 * mi * s1 + 1]);
                    }
                } else {
                    const R* a0 = src + po + y * ss[0];
                    const R* b0 = src + pmo + my * ss[0];
                    for (idx_t i = 0; i < iend; i++) {
                        const idx_t mi = m1 ? (i == 0 ? 0 : n1 - i) : i;
                        good &= a0[i * s1] == b0[mi * s1];
                    }
                }
                if (!good) ok = false;
                continue;
            }
            const bool yin = plane_in && y < S[0], myin = mplane_in && my < S[0];
            for (idx_t i = 0; i < n1; i++) {
                const idx_t mi = m1 ? (n1 - i) % n1 : i;
                R a = 0, b = 0, am = 0, bm = 0;
                if (yin && (nd < 2 || i < S[1])) {
                    const idx_t off = po + y * ss[0] + (nd > 1 ? i * ss[1] : 0);
                    a = cplx ? src[2 * off] : src[off];
                    b = cplx ? src[2 * off + 1] : R(0);
                }
                if (myin && (nd < 2 || mi < S[1])) {
                    const idx_t off = pmo + my * ss[0] + (nd > 1 ? mi * ss[1] : 0);
                    am = cplx ? src[2 * off] : src[off];
                    bm = cplx ? src[2 * off + 1] : R(0);
                }
                if (!(a == am && b == -bm)) {
                    ok = false;
                    break;
                }
            }
        }
    }
    return ok;
}

// a complex result whose imaginary parts are all zero is returned as real (contiguous arrays)
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

Array make(const Ctx& c, const dims_t& dims, bool is_double, bool is_complex) {
    Array a;
    a.dims = dims;
    a.strides = native_strides(dims);
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

dims_t dims_to(const dims_t& d, size_t nd) {
    dims_t r(d);
    while (r.size() < nd) r.push_back(1);
    return r;
}

// the strides of x for nd dimensions (those beyond its own: 0, their coordinate is always 0)
dims_t strides_to(const Array& x, size_t nd) {
    dims_t r(x.strides);
    while (r.size() < nd) r.push_back(0);
    return r;
}

// the source (dims S, strides ss; real or complex) padded with zeros or truncated to the native array dst of
// dims D; complex_out: dst is complex (a real source gets imaginary parts 0)
template <typename R>
void copy_padded(const R* src, bool cplx, const dims_t& S, const dims_t& ss, R* dst, bool complex_out, const dims_t& D,
        int threads) {
    const dims_t ds = native_strides(D);
    const Planes pl(D, ds);
    const int nd = static_cast<int>(D.size());
    const int cs = cplx ? 2 : 1, cd = complex_out ? 2 : 1;
    const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(threads, pl.P * pl.n0)));
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t t = 0; t < pl.P * pl.n0; t++) {
        const idx_t p = t / pl.n0, y = t % pl.n0;
        idx_t r = p, so = 0;
        bool in = y < S[0];
        for (int k = 2; k < nd; k++) {
            const idx_t c = r % D[k];
            r /= D[k];
            if (c >= S[k]) in = false;
            so += c * ss[k];
        }
        R* e = dst + cd * (pl.offset(p, false, nullptr) + y * ds[0]);
        const idx_t m1 = in ? std::min(pl.n1, nd > 1 ? S[1] : 1) : 0;
        const R* s = src + cs * (so + y * ss[0]);
        const idx_t s1 = nd > 1 ? ss[1] : 0;
        idx_t i = 0;
        for (; i < m1; i++) {
            e[cd * i] = s[cs * i * s1];
            if (complex_out) e[cd * i + 1] = cplx ? s[cs * i * s1 + 1] : R(0);
        }
        for (; i < pl.n1; i++) {
            e[cd * i] = 0;
            if (complex_out) e[cd * i + 1] = 0;
        }
    }
}

// x padded or truncated to D (zeros outside), same class, native layout
template <typename R>
Array resized(const Ctx& cx, const Array& x, const dims_t& D) {
    Array y = make(cx, D, x.is_double, x.is_complex);
    copy_padded(static_cast<const R*>(x.data), x.is_complex, dims_to(x.dims, D.size()), strides_to(x, D.size()),
        static_cast<R*>(y.data), x.is_complex, D, cx.threads);
    return y;
}

template <typename R>
Stage<R> stage(int be, int prec, int kind, int d, idx_t n) {
    Stage<R> st;
    st.be = be;
    st.prec = prec;
    st.kind = kind;
    st.d = d;
    st.n = n;
    return st;
}

template <typename R>
Array fft_t(const Array& x, idx_t n, int dim, const Options& o, const Ctx& cx) {
    const size_t nd = std::max<size_t>({x.dims.size(), static_cast<size_t>(dim) + 1, size_t(2)});
    const dims_t S = dims_to(x.dims, nd);
    dims_t D = S;
    if (n >= 0) D[dim] = n;
    const idx_t len = D[dim];
    // MATLAB returns the data unchanged along a dimension of length 1 (no transform)
    if (S[dim] == 1 && len == 1) return resized<R>(cx, x, D);
    Array y = make(cx, D, x.is_double, true);
    if (y.numel() == 0) return realify<R>(cx, y, o.threads);
    Stage<R> st = stage<R>(o.backend, x.is_double ? 1 : 0, x.is_complex ? kC2CF : kR2C, dim, len);
    st.src = static_cast<const R*>(x.data);
    st.src_complex = x.is_complex;
    st.S = S;
    st.ss = strides_to(x, nd);
    st.dst = static_cast<R*>(y.data);
    st.D = D;
    st.ds = y.strides;
    st.Dm = D;
    run_stage(st, o);
    return realify<R>(cx, y, o.threads);
}

template <typename R>
Array ifft_t(const Array& x, idx_t n, int dim, bool symmetric, const Options& o, const Ctx& cx) {
    const size_t nd = std::max<size_t>({x.dims.size(), static_cast<size_t>(dim) + 1, size_t(2)});
    const dims_t S = dims_to(x.dims, nd);
    const dims_t ss = strides_to(x, nd);
    dims_t D = S;
    if (n >= 0) D[dim] = n;
    const idx_t len = D[dim];
    if (S[dim] == 1 && len == 1 && !symmetric) return resized<R>(cx, x, D);
    const bool herm = symmetric || is_hermitian(static_cast<const R*>(x.data), x.is_complex, S, ss, D, dim, o.threads);
    Array y = make(cx, D, x.is_double, !herm);
    if (y.numel() == 0) return herm ? y : realify<R>(cx, y, o.threads);
    Stage<R> st = stage<R>(o.backend, x.is_double ? 1 : 0, herm ? kC2R : kC2CB, dim, len);
    st.src = static_cast<const R*>(x.data);
    st.src_complex = x.is_complex;
    st.S = S;
    st.ss = ss;
    st.dst = static_cast<R*>(y.data);
    st.D = D;
    st.ds = y.strides;
    st.Dm = D;
    st.scale = true;
    run_stage(st, o);
    return herm ? y : realify<R>(cx, y, o.threads);
}

// the dimensions fftn and ifftn transform (length > 1)
std::vector<int> transformed(const dims_t& D) {
    std::vector<int> t;
    for (size_t k = 0; k < D.size(); k++) {
        if (D[k] > 1) t.push_back(static_cast<int>(k));
    }
    return t;
}

// the forward stages of fftn on x into y (dims D); kw_last >= 0: the last stage writes its first kw_last
// elements into dst instead (dims Dd, a prefix of D along the last transformed dimension)
template <typename R>
void forward_stages(const Array& x, const dims_t& S, const dims_t& ss, Array& y, const dims_t& D,
        const std::vector<int>& tdims, R* last_dst, const dims_t& last_D, idx_t kw_last, const Options& o) {
    const int prec = x.is_double ? 1 : 0;
    for (size_t i = 0; i < tdims.size(); i++) {
        const int d = tdims[i];
        Stage<R> st = stage<R>(o.backend, prec, i == 0 ? (x.is_complex ? kC2CF : kR2C) : kC2CF, d, D[d]);
        if (i == 0) {
            st.src = static_cast<const R*>(x.data);
            st.src_complex = x.is_complex;
            st.S = S;
            st.ss = ss;
        } else {
            st.src = static_cast<const R*>(y.data);
            st.src_complex = true;
            st.S = D;
            st.ss = y.strides;
        }
        st.dst = static_cast<R*>(y.data);
        st.D = D;
        st.ds = y.strides;
        st.Dm = D;
        if (i + 1 == tdims.size() && kw_last >= 0) {
            st.kwrite = kw_last;
            if (last_dst != nullptr) {
                st.dst = last_dst;
                st.D = last_D;
                st.ds = native_strides(last_D);
            }
        }
        run_stage(st, o);
    }
}

template <typename R>
Array fftn_t(const Array& x, const dims_t& siz, const Options& o, const Ctx& cx) {
    const size_t nd = std::max({x.dims.size(), siz.size(), size_t(2)});
    const dims_t S = dims_to(x.dims, nd);
    const dims_t D = siz.empty() ? S : dims_to(siz, nd);
    const std::vector<int> tdims = transformed(D);
    idx_t total = 1;
    for (idx_t v : D) total *= v;
    if (tdims.empty() || total == 0) {
        Array r = resized<R>(cx, x, D);
        return realify<R>(cx, r, o.threads);
    }
    Array y = make(cx, D, x.is_double, true);
    // real input: the Hermitian fill below writes the elements past n / 2 along the last dimension
    const idx_t kw = x.is_complex ? -1 : D[tdims.back()] / 2 + 1;
    forward_stages<R>(x, S, strides_to(x, nd), y, D, tdims, nullptr, D, kw, o);
    if (!x.is_complex) hermitian_fill(static_cast<R*>(y.data), D, y.strides, o.threads);
    return realify<R>(cx, y, o.threads);
}

// MATLAB's last dimension of ifftn's complex-to-real stage: the last transformed one, 1 (x) at least
int last_dim(const std::vector<int>& tdims) {
    int last = 1;
    for (int d : tdims) last = std::max(last, d);
    return last;
}

// the self-conjugate planes (index 0 and n / 2) along dimension last of the native array h (dims H), as
// fftn's fill sees them
template <typename R>
void fill_planes(R* h, const dims_t& H, idx_t nl, int last, int threads) {
    const dims_t hs = native_strides(H);
    const dims_t Dp(H.begin(), H.begin() + last);
    const dims_t sp(hs.begin(), hs.begin() + last);
    const int nplanes = (nl % 2 == 0 && nl / 2 > 0) ? 2 : 1;
    for (int q = 0; q < nplanes; q++) {
        R* plane = h + 2 * hs[last] * (q == 0 ? 0 : nl / 2);
        hermitian_fill(plane, Dp, sp, threads);
    }
}

// complex stages of an inverse along cdims into h (dims H, a prefix of D along last), the first one from the
// source (src, dims S, strides ss: padded or truncated) or from h itself (src null), then complex-to-real along
// last into y (dims D)
template <typename R>
void inverse_stages(const R* src, bool src_complex, const dims_t& S, const dims_t& ss, Array& h, const dims_t& H,
        Array& y, const dims_t& D, const std::vector<int>& cdims, int last, const Options& o) {
    const int prec = h.is_double ? 1 : 0;
    for (size_t i = 0; i < cdims.size(); i++) {
        const int d = cdims[i];
        Stage<R> st = stage<R>(o.backend, prec, kC2CB, d, D[d]);
        if (i == 0 && src != nullptr) {
            st.src = src;
            st.src_complex = src_complex;
            st.S = S;
            st.ss = ss;
        } else {
            st.src = static_cast<const R*>(h.data);
            st.src_complex = true;
            st.S = H;
            st.ss = h.strides;
        }
        st.dst = static_cast<R*>(h.data);
        st.D = H;
        st.ds = h.strides;
        st.Dm = D;
        st.scale = true;
        run_stage(st, o);
    }
    const idx_t nl = D[last];
    if (nl > 1) {
        Stage<R> st = stage<R>(o.backend, prec, kC2R, last, nl);
        st.src = static_cast<const R*>(h.data);
        st.src_complex = true;
        st.S = H;
        st.ss = h.strides;
        st.dst = static_cast<R*>(y.data);
        st.D = D;
        st.ds = y.strides;
        st.Dm = D;
        st.scale = true;
        run_stage(st, o);
    } else {
        // length 1 along the last dimension: the real part
        const R* q = static_cast<const R*>(h.data);
        R* r = static_cast<R*>(y.data);
        const idx_t total = y.numel();
#pragma omp parallel for num_threads(std::max(1, o.threads)) schedule(static)
        for (idx_t i = 0; i < total; i++) r[i] = q[2 * i];
    }
}

template <typename R>
Array ifftn_t(const Array& x, const dims_t& siz, bool symmetric, const Options& o, const Ctx& cx) {
    const size_t nd = std::max({x.dims.size(), siz.size(), size_t(2)});
    const dims_t S = dims_to(x.dims, nd);
    const dims_t ss = strides_to(x, nd);
    const dims_t D = siz.empty() ? S : dims_to(siz, nd);
    const std::vector<int> tdims = transformed(D);
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
    const bool herm = symmetric || is_hermitian(static_cast<const R*>(x.data), x.is_complex, S, ss, D, -1, o.threads);
    if (herm) {
        // complex stages along the transformed dimensions but MATLAB's last one, on the part of the array the
        // last stage reads (the first n / 2 + 1 elements along it), then complex-to-real along it. MATLAB's
        // last dimension is that of its size (2 dimensions at least): for a column vector the second, of
        // length 1, whose complex-to-real transform keeps the real part.
        const int last = last_dim(tdims);
        const idx_t nl = D[last];
        std::vector<int> cdims;
        for (int d : tdims) {
            if (d != last) cdims.push_back(d);
        }
        Array y = make(cx, D, x.is_double, false);
        if (cdims.empty() && !symmetric && nl > 1) {
            Stage<R> st = stage<R>(o.backend, prec, kC2R, last, nl);
            st.src = static_cast<const R*>(x.data);
            st.src_complex = x.is_complex;
            st.S = S;
            st.ss = ss;
            st.dst = static_cast<R*>(y.data);
            st.D = D;
            st.ds = y.strides;
            st.Dm = D;
            st.scale = true;
            run_stage(st, o);
            return y;
        }
        dims_t H = D;
        H[last] = nl / 2 + 1;
        Array h = make(cx, H, x.is_double, true);
        if (symmetric || cdims.empty()) {
            // the source padded or truncated into h; 'symmetric': MATLAB first makes the padded or truncated
            // data conjugate symmetric with fftn's fill, which in this part changes only the planes 0 and n / 2
            // along the last dimension
            copy_padded(static_cast<const R*>(x.data), x.is_complex, S, ss, static_cast<R*>(h.data), true, H,
                o.threads);
            if (symmetric) fill_planes(static_cast<R*>(h.data), H, nl, last, o.threads);
            inverse_stages<R>(nullptr, true, H, h.strides, h, H, y, D, cdims, last, o);
        } else {
            inverse_stages<R>(static_cast<const R*>(x.data), x.is_complex, S, ss, h, H, y, D, cdims, last, o);
        }
        return y;
    }
    Array y = make(cx, D, x.is_double, true);
    for (size_t i = 0; i < tdims.size(); i++) {
        const int d = tdims[i];
        Stage<R> st = stage<R>(o.backend, prec, kC2CB, d, D[d]);
        if (i == 0) {
            st.src = static_cast<const R*>(x.data);
            st.src_complex = x.is_complex;
            st.S = S;
            st.ss = ss;
        } else {
            st.src = static_cast<const R*>(y.data);
            st.src_complex = true;
            st.S = D;
            st.ss = y.strides;
        }
        st.dst = static_cast<R*>(y.data);
        st.D = D;
        st.ds = y.strides;
        st.Dm = D;
        st.scale = true;
        run_stage(st, o);
    }
    return realify<R>(cx, y, o.threads);
}

// --- MATLAB's half spectrum of real data ----------------------------------------------------------------

template <typename R>
Array half_forward_t(const Array& x, const Options& o, const Ctx& cx) {
    const size_t nd = std::max<size_t>(x.dims.size(), 2);
    const dims_t D = dims_to(x.dims, nd);
    const dims_t ss = strides_to(x, nd);
    const dims_t H = half_dims(D);
    const std::vector<int> tdims = transformed(D);
    Array h = make(cx, H, x.is_double, true);
    if (tdims.empty()) {
        copy_padded(static_cast<const R*>(x.data), false, D, ss, static_cast<R*>(h.data), true, H, o.threads);
        return h;
    }
    const int last = last_dim(tdims);
    if (tdims.back() != last) {
        // no transform along MATLAB's last dimension (length 1): the half is the whole spectrum
        forward_stages<R>(x, D, ss, h, D, tdims, nullptr, D, -1, o);
    } else if (tdims.size() == 1) {
        // one transformed dimension: its real-to-complex stage writes the half directly
        Stage<R> st = stage<R>(o.backend, x.is_double ? 1 : 0, kR2C, last, D[last]);
        st.src = static_cast<const R*>(x.data);
        st.S = D;
        st.ss = ss;
        st.dst = static_cast<R*>(h.data);
        st.D = H;
        st.ds = h.strides;
        st.Dm = D;
        st.kwrite = H[last];
        run_stage(st, o);
    } else {
        // the full complex array up to the last stage, whose first n / 2 + 1 outputs go into h
        Array y = make(cx, D, x.is_double, true);
        forward_stages<R>(x, D, ss, y, D, tdims, static_cast<R*>(h.data), H, H[last], o);
    }
    // fftn's fill changes in this part only the self-conjugate planes 0 and n / 2 along the last dimension
    fill_planes(static_cast<R*>(h.data), H, D[last], last, o.threads);
    return h;
}

template <typename R>
Array half_inverse_t(Array& h, const dims_t& D0, const Options& o, const Ctx& cx) {
    const dims_t D = dims_to(D0, std::max<size_t>(D0.size(), 2));
    const dims_t H = half_dims(D);
    const std::vector<int> tdims = transformed(D);
    Array y = make(cx, D, h.is_double, false);
    const int last = last_dim(tdims);
    std::vector<int> cdims;
    for (int d : tdims) {
        if (d != last) cdims.push_back(d);
    }
    inverse_stages<R>(nullptr, true, H, h.strides, h, H, y, D, cdims, last, o);
    return y;
}

template <typename R>
void half_multiply_t(Array& h, const Array& g, bool conjugate, int threads) {
    R* a = static_cast<R*>(h.data);
    const R* b = static_cast<const R*>(g.data);
    const idx_t n = h.numel();
    const int nt = static_cast<int>(std::max<idx_t>(1, std::min<idx_t>(threads, (n >> 14) + 1)));
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t i = 0; i < n; i++) {
        const R pr = a[2 * i], pi = a[2 * i + 1], cr = b[2 * i], ci = b[2 * i + 1];
        if (conjugate) {
            a[2 * i] = pr * cr + pi * ci;
            a[2 * i + 1] = pi * cr - pr * ci;
        } else {
            a[2 * i] = pr * cr - pi * ci;
            a[2 * i + 1] = pr * ci + pi * cr;
        }
    }
}

}  // namespace

std::vector<idx_t> native_strides(const std::vector<idx_t>& D) {
    const size_t nd = D.size();
    std::vector<idx_t> s(nd, 1);
    if (nd == 0) return s;
    if (nd == 1) return s;
    s[1] = 1;
    s[0] = D[1];
    idx_t m = D[0] * D[1];
    for (size_t k = 2; k < nd; k++) {
        s[k] = m;
        m *= D[k];
    }
    return s;
}

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

std::vector<idx_t> half_dims(const std::vector<idx_t>& D0) {
    std::vector<idx_t> D = dims_to(D0, std::max<size_t>(D0.size(), 2));
    const std::vector<int> tdims = transformed(D);
    const int last = last_dim(tdims);
    D[last] = D[last] / 2 + 1;
    return D;
}

Array half_forward(const Array& x, const Options& o, Allocator alloc, void* ctx) {
    if (x.is_complex) throw std::invalid_argument("half_forward takes real data");
    const Ctx c{alloc, ctx, o.threads};
    return x.is_double ? half_forward_t<double>(x, o, c) : half_forward_t<float>(x, o, c);
}

Array half_inverse(Array& h, const std::vector<idx_t>& D, const Options& o, Allocator alloc, void* ctx) {
    if (!h.is_complex) throw std::invalid_argument("half_inverse takes a complex half spectrum");
    if (h.dims != half_dims(D)) throw std::invalid_argument("the half spectrum's size does not match");
    const Ctx c{alloc, ctx, o.threads};
    return h.is_double ? half_inverse_t<double>(h, D, o, c) : half_inverse_t<float>(h, D, o, c);
}

void half_multiply(Array& h, const Array& g, bool conjugate, int threads) {
    if (h.dims != g.dims || h.is_double != g.is_double || !h.is_complex || !g.is_complex) {
        throw std::invalid_argument("the spectra differ in size or class");
    }
    if (h.is_double) {
        half_multiply_t<double>(h, g, conjugate, threads);
    } else {
        half_multiply_t<float>(h, g, conjugate, threads);
    }
}

std::string plan_string(int be, bool is_double, const std::string& kind, idx_t n, int howmany) {
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
    void* p = get_plan(be, prec, k, n, howmany);
    std::lock_guard<std::mutex> g(planner_lock());
    const FftwApi& a = backend(be).api[prec];
    char* s = a.sprint(p);
    std::string r = s != nullptr ? s : "";
    if (s != nullptr) a.free_string(s);
    for (char& ch : r) {
        if (ch == '\n') ch = ' ';
    }
    return r;
}

void clear_plans(int be) {
    std::lock_guard<std::mutex> g(planner_lock());
    for (auto it = plan_cache().begin(); it != plan_cache().end();) {
        if (std::get<0>(it->first) == be) {
            backend(be).api[std::get<1>(it->first)].destroy(it->second);
            it = plan_cache().erase(it);
        } else {
            ++it;
        }
    }
}

}  // namespace pk5d::fft
