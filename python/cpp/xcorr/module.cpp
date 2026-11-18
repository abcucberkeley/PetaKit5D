// petakit5d._ext.xcorr: normxcorr3_fast's NCC in a window of its full output for the port's (z, y, x) C arrays
// (xcorr.h). The arrays move into MATLAB's layout (y fastest) at entry and the results back at exit; windows are
// given in MATLAB's order and 1-based ([s1, s2, s3, t1, t2, t3] of y, x, z), as Python's wrapper converts them.
// The module uses the process's one FFT state (backends, planner lock, plan cache) of petakit5d._ext.fft.
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <new>
#include <string>
#include <tuple>
#include <vector>

#include "backend.h"
#include "driver.h"
#include "shared.h"
#include "xcorr.h"

namespace nb = nanobind;
using namespace pk5d::xcorr;

namespace {

template <typename T>
using InArray = nb::ndarray<const T, nb::c_contig, nb::device::cpu>;
template <typename T>
using OutArray = nb::ndarray<nb::numpy, T>;

// uninitialized storage of n elements, owned until handed to NumPy
template <typename T>
struct Owned {
    T* p = nullptr;
    explicit Owned(size_t n) {
        p = static_cast<T*>(std::malloc(sizeof(T) * (n > 0 ? n : 1)));
        if (p == nullptr) throw std::bad_alloc();
    }
    ~Owned() { std::free(p); }
    Owned(const Owned&) = delete;
    Owned& operator=(const Owned&) = delete;
    // a C-contiguous NumPy array (nz, ny, nx) that owns the storage
    OutArray<T> release(size_t nz, size_t ny, size_t nx) {
        T* q = p;
        p = nullptr;
        nb::capsule owner(q, [](void* r) noexcept { std::free(r); });
        const size_t shape[3] = {nz, ny, nx};
        return OutArray<T>(q, 3, shape, owner);
    }
};

// MATLAB's size (y, x, z) of a (z, y, x) or (y, x) array, each below 2^30
template <typename T>
void msize(const InArray<T>& v, int* d, const char* name) {
    if (v.ndim() != 2 && v.ndim() != 3) {
        throw nb::value_error((std::string(name) + " must be a (z, y, x) volume or a (y, x) image").c_str());
    }
    const size_t nz = v.ndim() == 3 ? v.shape(0) : 1;
    const size_t ny = v.shape(v.ndim() - 2), nx = v.shape(v.ndim() - 1);
    const size_t s[3] = {ny, nx, nz};
    for (int i = 0; i < 3; i++) {
        if (s[i] < 1) throw nb::value_error((std::string(name) + " must not be empty").c_str());
        if (s[i] > 1073741823) {
            throw nb::value_error((std::string("the dimensions of ") + name + " must be smaller than 2^30").c_str());
        }
        d[i] = static_cast<int>(s[i]);
    }
}

// the sizes of a pair and of its full output; checks the template fits the image
template <typename T>
void pair_sizes(const InArray<T>& T_, const InArray<T>& A, int* ts, int* a, int* d) {
    msize(T_, ts, "T");
    msize(A, a, "A");
    for (int i = 0; i < 3; i++) {
        if (ts[i] > a[i]) throw nb::value_error("the template must not be larger than the image");
        d[i] = ts[i] + a[i] - 1;
    }
}

// the window from [s1, s2, s3, t1, t2, t3] (MATLAB's order, 1-based, inclusive) in an output of size d
Win window(const std::vector<int64_t>& w, const int* d) {
    if (w.size() != 6) throw nb::value_error("the window must have 6 elements");
    Win W;
    for (int i = 0; i < 3; i++) {
        const int64_t s = w[i], e = w[i + 3];
        if (!(s >= 1) || !(e <= d[i]) || !(s <= e)) {
            throw nb::value_error("the window must hold integers with 1 <= s <= t <= the size of the full output");
        }
        W.w0[i] = static_cast<int>(s) - 1;
        W.w1[i] = static_cast<int>(e) - 1;
        W.wn[i] = W.w1[i] - W.w0[i] + 1;
    }
    return W;
}

int threads_arg(int nt) { return nt < 1 ? 1 : (nt > 4096 ? 4096 : nt); }
idx_t grain_arg(int64_t g) { return g < 1 ? 0 : static_cast<idx_t>(g); }
idx_t numel(const int* d) { return (idx_t)d[0] * d[1] * d[2]; }

// a (z, y, x) array (MATLAB size d) into MATLAB's layout
template <typename T>
void to_m(const T* src, const int* d, T* dst, int nt) {
    c_to_matlab(src, d[2], d[0], d[1], dst, nt);
}

template <typename T>
OutArray<T> window_out(const T* m, const Win& W, int nt) {
    Owned<T> o((size_t)W.wn[0] * W.wn[1] * W.wn[2]);
    matlab_to_c(m, W.wn[2], W.wn[0], W.wn[1], o.p, nt);
    return o.release(W.wn[2], W.wn[0], W.wn[1]);
}

std::tuple<double, double, double, double, double> stats_tuple(const Stats& s) {
    return {s.cT, s.cA, s.sumT0, s.denomT, s.varA1};
}

Stats stats_from(const std::tuple<double, double, double, double, double>& t) {
    Stats s;
    s.cT = std::get<0>(t);
    s.cA = std::get<1>(t);
    s.sumT0 = std::get<2>(t);
    s.denomT = std::get<3>(t);
    s.varA1 = std::get<4>(t);
    return s;
}

// normxcorr3_fast(T, A, 'full', bbox) of single data with the mex: C in the window (z, y, x), the 0-based
// linear index of its first maximum in MATLAB's order, the statistics, and seconds of [moves and statistics,
// planning, forward T0, forward A0, inverse, window sums, NCC]
std::tuple<OutArray<float>, int64_t, std::tuple<double, double, double, double, double>, std::vector<double>> py_ncc(
    const InArray<float>& T_, const InArray<float>& A, const std::vector<int64_t>& win, int threads, int64_t grain,
    bool fast, int backend) {
    int ts[3], a[3], d[3];
    pair_sizes(T_, A, ts, a, d);
    if (d[0] < 2 || d[1] < 2) throw nb::value_error("the full output must have at least 2 elements along y and x");
    const Win W = window(win, d);
    const int nt = threads_arg(threads);
    pk5d::fft::backend(backend);  // throws if it is not loaded
    const idx_t nw = (idx_t)W.wn[0] * W.wn[1] * W.wn[2];
    Owned<float> Tm(numel(ts)), Am(numel(a)), corr(nw), C(nw);
    Owned<double> sA(nw), sA2(nw), sT(nw);
    std::vector<double> tm(7, 0.0);
    Stats st;
    idx_t imax = 0;
    {
        nb::gil_scoped_release g;
        const auto clock = [] {
            return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
        };
        double t = clock();
        to_m(T_.data(), ts, Tm.p, nt);
        to_m(A.data(), a, Am.p, nt);
        st = prepare(Tm.p, numel(ts), Am.p, numel(a), Tm.p, Am.p, nt);
        tm[0] = clock() - t;
        double tc[4];
        window_corr(Tm.p, ts, Am.p, a, W, nt, grain_arg(grain), fast, backend, corr.p, tc);
        for (int i = 0; i < 4; i++) tm[1 + i] = tc[i];
        t = clock();
        window_sums(Tm.p, ts, Am.p, a, W, sA.p, sA2.p, sT.p, nt);
        tm[5] = clock() - t;
        t = clock();
        ncc_combine(corr.p, sA.p, sA2.p, sT.p, ts, a, W, st, C.p, nt);
        imax = first_max(C.p, nw, nt);
        tm[6] = clock() - t;
    }
    OutArray<float> out = window_out(C.p, W, nt);
    return {out, (int64_t)imax, stats_tuple(st), tm};
}

// normxcorr3_window_mex('ncc', T0, A0, bbox, nt, grain, fast): (corr, sA, sA2, sT) in the window, (z, y, x)
std::tuple<OutArray<float>, OutArray<double>, OutArray<double>, OutArray<double>> py_window(const InArray<float>& T0,
    const InArray<float>& A0, const std::vector<int64_t>& win, int threads, int64_t grain, bool fast, int backend) {
    int ts[3], a[3], d[3];
    pair_sizes(T0, A0, ts, a, d);
    if (d[0] < 2 || d[1] < 2) throw nb::value_error("the full output must have at least 2 elements along y and x");
    const Win W = window(win, d);
    const int nt = threads_arg(threads);
    pk5d::fft::backend(backend);
    const idx_t nw = (idx_t)W.wn[0] * W.wn[1] * W.wn[2];
    Owned<float> Tm(numel(ts)), Am(numel(a)), corr(nw);
    Owned<double> sA(nw), sA2(nw), sT(nw);
    {
        nb::gil_scoped_release g;
        to_m(T0.data(), ts, Tm.p, nt);
        to_m(A0.data(), a, Am.p, nt);
        window_corr(Tm.p, ts, Am.p, a, W, nt, grain_arg(grain), fast, backend, corr.p, nullptr);
        window_sums(Tm.p, ts, Am.p, a, W, sA.p, sA2.p, sT.p, nt);
    }
    return {window_out(corr.p, W, nt), window_out(sA.p, W, nt), window_out(sA2.p, W, nt), window_out(sT.p, W, nt)};
}

// the statistics of normxcorr3_fast and T0 = T - cT, A0 = A - cA ((z, y, x), the class of the data)
template <typename T>
std::tuple<OutArray<T>, OutArray<T>, std::tuple<double, double, double, double, double>> py_prepare(
    const InArray<T>& T_, const InArray<T>& A, int threads) {
    int ts[3], a[3], d[3];
    pair_sizes(T_, A, ts, a, d);
    const int nt = threads_arg(threads);
    Owned<T> Tm(numel(ts)), Am(numel(a)), Tc(numel(ts)), Ac(numel(a));
    Stats st;
    {
        nb::gil_scoped_release g;
        to_m(T_.data(), ts, Tm.p, nt);
        to_m(A.data(), a, Am.p, nt);
        st = prepare(Tm.p, numel(ts), Am.p, numel(a), Tm.p, Am.p, nt);
        matlab_to_c(Tm.p, ts[2], ts[0], ts[1], Tc.p, nt);
        matlab_to_c(Am.p, a[2], a[0], a[1], Ac.p, nt);
    }
    return {Tc.release(ts[2], ts[0], ts[1]), Ac.release(a[2], a[0], a[1]), stats_tuple(st)};
}

// the window sums and the NCC from the correlation in the window (corr: (z, y, x), the class of the data, made
// by MATLAB's code without the mex) and T0, A0: (C, index of the first maximum in MATLAB's order)
template <typename T>
std::tuple<OutArray<T>, int64_t> py_finish(const InArray<T>& corr, const InArray<T>& T0, const InArray<T>& A0,
    const std::vector<int64_t>& win, const std::tuple<double, double, double, double, double>& stats, int threads) {
    int ts[3], a[3], d[3];
    pair_sizes(T0, A0, ts, a, d);
    const Win W = window(win, d);
    int cw[3];
    msize(corr, cw, "corr");
    for (int i = 0; i < 3; i++) {
        if (cw[i] != W.wn[i]) throw nb::value_error("corr must have the window's size");
    }
    const int nt = threads_arg(threads);
    const idx_t nw = (idx_t)W.wn[0] * W.wn[1] * W.wn[2];
    const Stats st = stats_from(stats);
    Owned<T> Tm(numel(ts)), Am(numel(a)), cm(nw), C(nw);
    Owned<double> sA(nw), sA2(nw), sT(nw);
    idx_t imax = 0;
    {
        nb::gil_scoped_release g;
        to_m(T0.data(), ts, Tm.p, nt);
        to_m(A0.data(), a, Am.p, nt);
        to_m(corr.data(), cw, cm.p, nt);
        window_sums(Tm.p, ts, Am.p, a, W, sA.p, sA2.p, sT.p, nt);
        ncc_combine(cm.p, sA.p, sA2.p, sT.p, ts, a, W, st, C.p, nt);
        imax = first_max(C.p, nw, nt);
    }
    return {window_out(C.p, W, nt), (int64_t)imax};
}

}  // namespace

NB_MODULE(xcorr, m) {
    // the process's one FFT state (backends, planner lock, plan cache): petakit5d._ext.fft's
    pk5d::fft::attach_to_fft_module();
    pk5d::fft::def_shared_state_hooks(m);
    m.doc() = "normxcorr3_fast's NCC in a window of its full output, for (z, y, x) C arrays";
    m.def("ncc", &py_ncc, nb::arg("T").noconvert(), nb::arg("A").noconvert(), nb::arg("win"), nb::arg("threads"),
        nb::arg("grain"), nb::arg("fast"), nb::arg("backend"),
        "(C, first maximum's MATLAB index, statistics, seconds) of single data in the window");
    m.def("window", &py_window, nb::arg("T0").noconvert(), nb::arg("A0").noconvert(), nb::arg("win"),
        nb::arg("threads"), nb::arg("grain"), nb::arg("fast"), nb::arg("backend"),
        "normxcorr3_window_mex('ncc', ...): (corr, sA, sA2, sT) in the window");
    m.def("prepare", &py_prepare<float>, nb::arg("T").noconvert(), nb::arg("A").noconvert(), nb::arg("threads"));
    m.def("prepare", &py_prepare<double>, nb::arg("T").noconvert(), nb::arg("A").noconvert(), nb::arg("threads"),
        "(T0, A0, (cT, cA, sumT0, denomT, varA1))");
    m.def("finish", &py_finish<float>, nb::arg("corr").noconvert(), nb::arg("T0").noconvert(),
        nb::arg("A0").noconvert(), nb::arg("win"), nb::arg("stats"), nb::arg("threads"));
    m.def("finish", &py_finish<double>, nb::arg("corr").noconvert(), nb::arg("T0").noconvert(),
        nb::arg("A0").noconvert(), nb::arg("win"), nb::arg("stats"), nb::arg("threads"),
        "(C, first maximum's MATLAB index) from the window's correlation");
    m.attr("PERMISSIVE") = static_cast<int>(pk5d::fft::kPermissive);
    m.attr("EXACT") = static_cast<int>(pk5d::fft::kExact);
    m.attr("GRAIN_PARALLEL") = pk5d::fft::kGrainParallel;
}
