// petakit5d._ext.decon: the bindings of the deconvolution kernels (decon.hpp). The Python side
// (petakit5d.kernels.decon) checks the arguments with the mex's rules; here only sizes are checked. Arrays
// are C-contiguous; arrays in MATLAB's layout come as flat views. The interpreter lock is released around
// the loops; FFT plans are made and destroyed with it held (decon.hpp, RfftPlan).
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <complex>
#include <cstdlib>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "backend.h"
#include "decon.hpp"

namespace nb = nanobind;
using namespace pk5d::decon;

namespace {

template <typename T>
using flat_in = nb::ndarray<const T, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
template <typename T>
using nd_in = nb::ndarray<const T, nb::c_contig, nb::device::cpu>;
using any_in = nb::ndarray<nb::ro, nb::c_contig, nb::device::cpu>;
using any_out = nb::ndarray<nb::c_contig, nb::device::cpu>;
using times_t = std::tuple<double, double, double, double>;

template <typename T>
nb::ndarray<nb::numpy> wrap(T* p, std::vector<size_t> shape) {
    nb::capsule owner(p, [](void* q) noexcept { std::free(q); });
    return nb::ndarray<nb::numpy>(p, shape.size(), shape.data(), owner, nullptr, nb::dtype<T>());
}

nb::ndarray<nb::numpy> wrap_bool(uint8_t* p, std::vector<size_t> shape) {
    nb::capsule owner(p, [](void* q) noexcept { std::free(q); });
    return nb::ndarray<nb::numpy>(p, shape.size(), shape.data(), owner, nullptr, nb::dtype<bool>());
}

times_t times_of(const Times& t) { return {t.t[0], t.t[1], t.t[2], t.t[3]}; }

template <typename A>
size_t numel(const A& a) {
    size_t n = 1;
    for (size_t k = 0; k < a.ndim(); k++) n *= a.shape(k);
    return n;
}

template <typename A>
std::vector<std::int64_t> shape_of(const A& a) {
    std::vector<std::int64_t> s;
    for (size_t k = 0; k < a.ndim(); k++) s.push_back(static_cast<std::int64_t>(a.shape(k)));
    return s;
}

std::vector<size_t> sizes(const std::vector<std::int64_t>& s) { return std::vector<size_t>(s.begin(), s.end()); }

Call make_call(const std::vector<int>& d, int threads, std::int64_t grain, int be) {
    if (d.size() != 3 || d[0] < 1 || d[1] < 1 || d[2] < 1) throw std::invalid_argument("decon: the size is (ny, nx, nz)");
    Call c;
    c.d[0] = d[0];
    c.d[1] = d[1];
    c.d[2] = d[2];
    c.nt = threads > 0 ? threads : 1;
    c.grain = grain > 0 ? grain : 0;
    c.backend = be;
    return c;
}

size_t half_count(const Call& c) {
    const size_t ny = c.d[0], nx = c.d[1], nz = c.d[2];
    return nz > 1 ? ny * nx * (nz / 2 + 1) : ny * (nx / 2 + 1);
}

size_t voxels(const Call& c) { return static_cast<size_t>(c.d[0]) * c.d[1] * c.d[2]; }

void check_size(size_t have, size_t want, const char* what) {
    if (have != want) throw std::invalid_argument(std::string("decon: ") + what + " has the wrong number of elements");
}

// --- MATLAB's layout ------------------------------------------------------------------------------

std::pair<nb::ndarray<nb::numpy>, times_t> py_matlab_otf(const nd_in<float>& psf, const std::vector<int>& d,
                                                          int threads, std::int64_t grain, int be) {
    const Call call = make_call(d, threads, grain, be);
    if (psf.ndim() != 3) throw std::invalid_argument("decon: the PSF is (pz, py, px)");
    const int ps[3] = {static_cast<int>(psf.shape(0)), static_cast<int>(psf.shape(1)), static_cast<int>(psf.shape(2))};
    if (ps[0] > call.d[2] || ps[1] > call.d[0] || ps[2] > call.d[1]) throw std::invalid_argument("decon: the PSF is larger than the data");
    Times tm;
    RfftPlanPtr plan = matlab_plan(call, tm);
    Owned<cfloat> h(half_count(call));
    {
        nb::gil_scoped_release g;
        matlab_otf(*plan, psf.data(), ps, h.p, tm);
    }
    const size_t nh = half_count(call);
    return {wrap(h.release(), {nh}), times_of(tm)};
}

template <typename T>
std::pair<nb::ndarray<nb::numpy>, times_t> py_matlab_omw(const flat_in<T>& raw, double bg, const flat_in<cfloat>& hf,
                                                          const flat_in<cfloat>& hb, int n_iter,
                                                          const std::optional<flat_in<float>>& j0,
                                                          const std::vector<int>& d, int threads, std::int64_t grain,
                                                          int be, bool native_io) {
    const Call call = make_call(d, threads, grain, be);
    check_size(raw.shape(0), voxels(call), "the data");
    check_size(hf.shape(0), half_count(call), "an OTF");
    check_size(hb.shape(0), half_count(call), "an OTF");
    if (j0) check_size(j0->shape(0), voxels(call), "J0");
    Times tm;
    RfftPlanPtr plan = matlab_plan(call, tm);
    float* r;
    {
        nb::gil_scoped_release g;
        r = matlab_omw(*plan, raw.data(), bg, hf.data(), hb.data(), n_iter, j0 ? j0->data() : nullptr, native_io, tm);
    }
    return {wrap(r, {voxels(call)}), times_of(tm)};
}

std::pair<nb::ndarray<nb::numpy>, times_t> py_matlab_rlstep(const flat_in<float>& y, const flat_in<float>& im,
                                                             const flat_in<cfloat>& h, const std::vector<int>& d,
                                                             int threads, std::int64_t grain, int be, bool native_io) {
    const Call call = make_call(d, threads, grain, be);
    check_size(y.shape(0), voxels(call), "Y");
    check_size(im.shape(0), voxels(call), "I");
    check_size(h.shape(0), half_count(call), "the OTF");
    Times tm;
    RfftPlanPtr plan = matlab_plan(call, tm);
    float* r;
    {
        nb::gil_scoped_release g;
        r = matlab_rlstep(*plan, y.data(), im.data(), h.data(), native_io, tm);
    }
    return {wrap(r, {voxels(call)}), times_of(tm)};
}

nb::ndarray<nb::numpy> py_matlab_fft1(const flat_in<float>& x, int n, int m, int threads, std::int64_t grain, int be) {
    const Call call = make_call({n, m, 1}, threads, grain, be);
    check_size(x.shape(0), static_cast<size_t>(n) * m, "X");
    Times tm;
    RfftPlanPtr plan = matlab_plan(call, tm);
    const size_t nh = static_cast<size_t>(n / 2 + 1) * m;
    Owned<cfloat> y(nh);
    {
        nb::gil_scoped_release g;
        matlab_fft1(*plan, x.data(), y.p);
    }
    return wrap(y.release(), {nh});
}

// --- the port's layout, native half spectra ------------------------------------------------------------

std::vector<std::int64_t> half_shape(std::vector<std::int64_t> s) {
    s.back() = s.back() / 2 + 1;
    return s;
}

std::pair<nb::ndarray<nb::numpy>, times_t> py_native_otf(const nd_in<float>& psf, const std::vector<std::int64_t>& shape,
                                                          int threads) {
    if (psf.ndim() != shape.size() || (shape.size() != 2 && shape.size() != 3)) {
        throw std::invalid_argument("decon: the PSF and the data have 2 or 3 axes");
    }
    const std::vector<std::int64_t> ps = shape_of(psf);
    for (size_t k = 0; k < shape.size(); k++) {
        if (ps[k] > shape[k] || shape[k] < 1) throw std::invalid_argument("decon: the PSF is larger than the data");
    }
    const std::vector<std::int64_t> hs = half_shape(shape);
    size_t nh = 1;
    for (auto v : hs) nh *= static_cast<size_t>(v);
    Owned<cfloat> h(nh);
    Times tm;
    {
        nb::gil_scoped_release g;
        native_otf(shape, psf.data(), ps, h.p, threads > 0 ? threads : 1, tm);
    }
    return {wrap(h.release(), sizes(hs)), times_of(tm)};
}

template <typename T>
std::pair<nb::ndarray<nb::numpy>, times_t> py_native_omw(const nd_in<T>& raw, double bg, const nd_in<cfloat>& hf,
                                                          const nd_in<cfloat>& hb, int n_iter,
                                                          const std::optional<nd_in<float>>& j0, int threads) {
    const std::vector<std::int64_t> shape = shape_of(raw);
    if (shape.size() != 2 && shape.size() != 3) throw std::invalid_argument("decon: 2D or 3D data");
    const std::vector<std::int64_t> hs = half_shape(shape);
    if (shape_of(hf) != hs || shape_of(hb) != hs) throw std::invalid_argument("decon: an OTF's shape is not the data's half spectrum");
    if (j0 && shape_of(*j0) != shape) throw std::invalid_argument("decon: J0 has another shape");
    Times tm;
    float* r;
    {
        nb::gil_scoped_release g;
        r = native_omw(shape, raw.data(), bg, hf.data(), hb.data(), n_iter, j0 ? j0->data() : nullptr,
                       threads > 0 ? threads : 1, tm);
    }
    return {wrap(r, sizes(shape)), times_of(tm)};
}

std::pair<nb::ndarray<nb::numpy>, times_t> py_native_rlstep(const nd_in<float>& y, const nd_in<float>& im,
                                                             const nd_in<cfloat>& h, int threads) {
    const std::vector<std::int64_t> shape = shape_of(y);
    if (shape.size() != 2 && shape.size() != 3) throw std::invalid_argument("decon: 2D or 3D data");
    if (shape_of(im) != shape) throw std::invalid_argument("decon: Y and I differ in shape");
    if (shape_of(h) != half_shape(shape)) throw std::invalid_argument("decon: the OTF's shape is not the data's half spectrum");
    Times tm;
    float* r;
    {
        nb::gil_scoped_release g;
        r = native_rlstep(shape, y.data(), im.data(), h.data(), threads > 0 ? threads : 1, tm);
    }
    return {wrap(r, sizes(shape)), times_of(tm)};
}

// --- acceleration ------------------------------------------------------------------------------------

template <typename T>
using flat_out = nb::ndarray<T, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

// the result array: out (returned as None) or a new one
template <typename T>
T* result_buffer(std::optional<flat_out<T>>& out, Owned<T>& own, size_t n) {
    if (out) {
        check_size(out->shape(0), n, "out");
        return out->data();
    }
    return own.alloc(n);
}

template <typename T>
nb::object result_object(std::optional<flat_out<T>>& out, Owned<T>& own, size_t n) {
    if (out) return nb::none();
    return nb::cast(wrap(own.release(), {n}));
}

template <typename T>
std::tuple<double, nb::object, double, double> py_accel_lambda(const flat_in<T>& j2, const flat_in<T>& y,
                                                               const flat_in<T>& j4, bool zyx, std::int64_t ny,
                                                               std::int64_t nx, int threads,
                                                               std::optional<flat_out<T>> out) {
    const idx_t n = static_cast<idx_t>(j2.shape(0));
    check_size(y.shape(0), static_cast<size_t>(n), "Y");
    check_size(j4.shape(0), static_cast<size_t>(n), "J_4");
    if (zyx && (ny < 1 || nx < 1 || n % (ny * nx) != 0)) throw std::invalid_argument("decon: (planes, ny, nx) sizes");
    Owned<T> own;
    // (out may be J_4 itself: each element is read before it is written)
    T* o = result_buffer(out, own, static_cast<size_t>(n));
    double s[2];
    {
        nb::gil_scoped_release g;
        accel_lambda(j2.data(), y.data(), j4.data(), o, n, zyx, ny, nx, threads > 0 ? threads : 1, s);
    }
    return {clamped_lambda(s), result_object(out, own, static_cast<size_t>(n)), s[0], s[1]};
}

template <typename T>
nb::object py_accel_predict(const flat_in<T>& j2, const flat_in<T>& j3, double lambda, int threads,
                            std::optional<flat_out<T>> out) {
    const idx_t n = static_cast<idx_t>(j2.shape(0));
    check_size(j3.shape(0), static_cast<size_t>(n), "J_3");
    Owned<T> own;
    T* o = result_buffer(out, own, static_cast<size_t>(n));
    {
        nb::gil_scoped_release g;
        accel_predict(j2.data(), j3.data(), lambda, o, n, threads > 0 ? threads : 1);
    }
    return result_object(out, own, static_cast<size_t>(n));
}

template <typename T>
double py_accel_sum(const flat_in<T>& x, bool zyx, std::int64_t ny, std::int64_t nx, int threads) {
    const idx_t n = static_cast<idx_t>(x.shape(0));
    if (zyx && n > 0 && (ny < 1 || nx < 1 || n % (ny * nx) != 0)) throw std::invalid_argument("decon: (planes, ny, nx) sizes");
    nb::gil_scoped_release g;
    return accel_sum(x.data(), n, zyx, ny, nx, threads > 0 ? threads : 1);
}

// --- erosion, back projector -----------------------------------------------------------------------------

nb::ndarray<nb::numpy> py_mask_edge_erosion(const nb::ndarray<const uint8_t, nb::ndim<3>, nb::c_contig, nb::device::cpu>& a,
                                            const std::vector<int>& start, const std::vector<int>& dr,
                                            const std::vector<int>& dc, const std::vector<std::int64_t>& pu,
                                            const std::vector<std::int64_t>& pl, std::int64_t rz, bool clear_border,
                                            int threads) {
    ErosionPlan p;
    if (start.empty() || start.front() != 0 || static_cast<size_t>(start.back()) != dr.size() || dr.size() != dc.size() ||
        pu.size() != 2 || pl.size() != 2 || rz < 0) {
        throw std::invalid_argument("decon: an erosion plan");
    }
    p.start = start;
    p.dr = dr;
    p.dc = dc;
    p.pu[0] = pu[0];
    p.pu[1] = pu[1];
    p.pl[0] = pl[0];
    p.pl[1] = pl[1];
    p.rz = rz;
    p.clear_border = clear_border;
    const idx_t nz = a.shape(0), ny = a.shape(1), nx = a.shape(2);
    Owned<uint8_t> b(static_cast<size_t>(nz * ny * nx));
    {
        nb::gil_scoped_release g;
        mask_edge_erosion(a.data(), b.p, nz, ny, nx, p, threads > 0 ? threads : 1);
    }
    return wrap_bool(b.release(), {static_cast<size_t>(nz), static_cast<size_t>(ny), static_cast<size_t>(nx)});
}

std::pair<nb::ndarray<nb::numpy>, nb::ndarray<nb::numpy>> py_omw_window(const std::vector<std::int64_t>& psz,
                                                                        const flat_in<double>& e,
                                                                        const flat_in<double>& r, double l, double u,
                                                                        int threads) {
    if (psz.size() != 3) throw std::invalid_argument("decon: psz is (n1, n2, n3)");
    const idx_t ne = static_cast<idx_t>(r.shape(0));
    check_size(e.shape(0), static_cast<size_t>(3 * ne), "E");
    const idx_t ps[3] = {psz[0], psz[1], psz[2]};
    const size_t nv = static_cast<size_t>(ps[0] * ps[1] * ps[2]);
    Owned<double> mask(nv);
    std::vector<int64_t> inds;
    {
        nb::gil_scoped_release g;
        omw_window(ps, e.data(), r.data(), ne, l, u, mask.p, inds, threads > 0 ? threads : 1);
    }
    Owned<int64_t> io(inds.size());
    if (!inds.empty()) std::memcpy(io.p, inds.data(), sizeof(int64_t) * inds.size());
    return {wrap(mask.release(), {static_cast<size_t>(ps[2]), static_cast<size_t>(ps[0]), static_cast<size_t>(ps[1])}),
            wrap(io.release(), {inds.size()})};
}

nb::ndarray<nb::numpy> py_omw_components(const nb::ndarray<const uint8_t, nb::ndim<3>, nb::c_contig, nb::device::cpu>& a,
                                         bool largest, std::int64_t n, int threads) {
    // (z, y, x) in, MATLAB's layout inside, (z, y, x) out
    const idx_t nz = a.shape(0), ny = a.shape(1), nx = a.shape(2);
    const size_t nv = static_cast<size_t>(nz * ny * nx);
    Owned<uint8_t> am(nv), bm(nv), b(nv);
    const int nt = threads > 0 ? threads : 1;
    {
        nb::gil_scoped_release g;
        transpose_planes(a.data(), am.p, nz, ny, nx, 1, nt);
        omw_components(am.p, ny, nx, nz, largest, n, bm.p, nt);
        transpose_planes(bm.p, b.p, nz, nx, ny, 1, nt);
    }
    return wrap_bool(b.release(), {static_cast<size_t>(nz), static_cast<size_t>(ny), static_cast<size_t>(nx)});
}

nb::object py_omw_sortdesc(const flat_in<float>& a, int threads) {
    const idx_t n = static_cast<idx_t>(a.shape(0));
    Owned<double> out(static_cast<size_t>(n));
    bool ok;
    {
        nb::gil_scoped_release g;
        ok = omw_sortdesc(a.data(), n, out.p, threads > 0 ? threads : 1);
    }
    if (!ok) return nb::none();
    return nb::cast(wrap(out.release(), {static_cast<size_t>(n)}));
}

// --- layouts, backends ------------------------------------------------------------------------------------

void py_transpose(const any_in& src, any_out& dst, std::int64_t rows, std::int64_t cols, int threads) {
    const size_t n = numel(src);
    if (numel(dst) != n || src.itemsize() != dst.itemsize() || rows < 0 || cols < 0 ||
        (rows * cols > 0 && n % static_cast<size_t>(rows * cols) != 0)) {
        throw std::invalid_argument("decon: transpose sizes");
    }
    const idx_t nz = rows * cols > 0 ? static_cast<idx_t>(n) / (rows * cols) : 0;
    nb::gil_scoped_release g;
    transpose_planes(src.data(), dst.data(), nz, rows, cols, static_cast<int>(src.itemsize()), threads > 0 ? threads : 1);
}

void py_load_exact(const std::string& libf, const std::string& libd) { pk5d::fft::load_exact(libf, libd); }

}  // namespace

NB_MODULE(decon, m) {
    m.doc() = "the deconvolution kernels: half-spectrum Richardson-Lucy steps, acceleration, mask erosion, back projector steps";
    m.def("matlab_otf", &py_matlab_otf, nb::arg("psf").noconvert(), nb::arg("d"), nb::arg("threads"), nb::arg("grain"),
          nb::arg("backend"));
    m.def("matlab_omw_uint16", &py_matlab_omw<uint16_t>, nb::arg("raw").noconvert(), nb::arg("bg"),
          nb::arg("hf").noconvert(), nb::arg("hb").noconvert(), nb::arg("n_iter"), nb::arg("j0").noconvert().none(),
          nb::arg("d"), nb::arg("threads"), nb::arg("grain"), nb::arg("backend"), nb::arg("native_io"));
    m.def("matlab_omw_single", &py_matlab_omw<float>, nb::arg("raw").noconvert(), nb::arg("bg"),
          nb::arg("hf").noconvert(), nb::arg("hb").noconvert(), nb::arg("n_iter"), nb::arg("j0").noconvert().none(),
          nb::arg("d"), nb::arg("threads"), nb::arg("grain"), nb::arg("backend"), nb::arg("native_io"));
    m.def("matlab_rlstep", &py_matlab_rlstep, nb::arg("y").noconvert(), nb::arg("im").noconvert(),
          nb::arg("h").noconvert(), nb::arg("d"), nb::arg("threads"), nb::arg("grain"), nb::arg("backend"),
          nb::arg("native_io"));
    m.def("matlab_fft1", &py_matlab_fft1, nb::arg("x").noconvert(), nb::arg("n"), nb::arg("m"), nb::arg("threads"),
          nb::arg("grain"), nb::arg("backend"));
    m.def("native_otf", &py_native_otf, nb::arg("psf").noconvert(), nb::arg("shape"), nb::arg("threads"));
    m.def("native_omw_uint16", &py_native_omw<uint16_t>, nb::arg("raw").noconvert(), nb::arg("bg"),
          nb::arg("hf").noconvert(), nb::arg("hb").noconvert(), nb::arg("n_iter"), nb::arg("j0").noconvert().none(),
          nb::arg("threads"));
    m.def("native_omw_single", &py_native_omw<float>, nb::arg("raw").noconvert(), nb::arg("bg"),
          nb::arg("hf").noconvert(), nb::arg("hb").noconvert(), nb::arg("n_iter"), nb::arg("j0").noconvert().none(),
          nb::arg("threads"));
    m.def("native_rlstep", &py_native_rlstep, nb::arg("y").noconvert(), nb::arg("im").noconvert(),
          nb::arg("h").noconvert(), nb::arg("threads"));
    m.def("accel_lambda_single", &py_accel_lambda<float>, nb::arg("j2").noconvert(), nb::arg("y").noconvert(),
          nb::arg("j4").noconvert(), nb::arg("zyx"), nb::arg("ny"), nb::arg("nx"), nb::arg("threads"),
          nb::arg("out").noconvert().none());
    m.def("accel_lambda_double", &py_accel_lambda<double>, nb::arg("j2").noconvert(), nb::arg("y").noconvert(),
          nb::arg("j4").noconvert(), nb::arg("zyx"), nb::arg("ny"), nb::arg("nx"), nb::arg("threads"),
          nb::arg("out").noconvert().none());
    m.def("accel_predict_single", &py_accel_predict<float>, nb::arg("j2").noconvert(), nb::arg("j3").noconvert(),
          nb::arg("lam"), nb::arg("threads"), nb::arg("out").noconvert().none());
    m.def("accel_predict_double", &py_accel_predict<double>, nb::arg("j2").noconvert(), nb::arg("j3").noconvert(),
          nb::arg("lam"), nb::arg("threads"), nb::arg("out").noconvert().none());
    m.def("accel_sum_single", &py_accel_sum<float>, nb::arg("x").noconvert(), nb::arg("zyx"), nb::arg("ny"),
          nb::arg("nx"), nb::arg("threads"));
    m.def("accel_sum_double", &py_accel_sum<double>, nb::arg("x").noconvert(), nb::arg("zyx"), nb::arg("ny"),
          nb::arg("nx"), nb::arg("threads"));
    m.def("mask_edge_erosion", &py_mask_edge_erosion, nb::arg("a").noconvert(), nb::arg("start"), nb::arg("dr"),
          nb::arg("dc"), nb::arg("pu"), nb::arg("pl"), nb::arg("rz"), nb::arg("clear_border"), nb::arg("threads"));
    m.def("omw_window", &py_omw_window, nb::arg("psz"), nb::arg("e").noconvert(), nb::arg("r").noconvert(),
          nb::arg("l"), nb::arg("u"), nb::arg("threads"));
    m.def("omw_components", &py_omw_components, nb::arg("a").noconvert(), nb::arg("largest"), nb::arg("n"),
          nb::arg("threads"));
    m.def("omw_sortdesc", &py_omw_sortdesc, nb::arg("a").noconvert(), nb::arg("threads"));
    m.def("transpose_planes", &py_transpose, nb::arg("src").noconvert(), nb::arg("dst").noconvert(), nb::arg("rows"),
          nb::arg("cols"), nb::arg("threads"),
          "dst[z][c][r] = src[z][r][c] for planes of rows x cols elements");
    m.def("load_exact", &py_load_exact, nb::arg("libf"), nb::arg("libd"),
          "loads the exact backend's single and double precision libraries into this module's backend table");
    m.def("exact_loaded", &pk5d::fft::exact_loaded);
    m.attr("PERMISSIVE") = static_cast<int>(pk5d::fft::kPermissive);
    m.attr("EXACT") = static_cast<int>(pk5d::fft::kExact);
    m.attr("GRAIN_PARALLEL") = 32;
}
