// petakit5d._ext.fft: MATLAB's fft, ifft, fftn and ifftn and the half-spectrum transforms (driver.h,
// halfspec.h) for the port's C-ordered (..., z, y, x) NumPy arrays, on the permissive or the exact backend.
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <complex>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "backend.h"
#include "driver.h"
#include "halfspec.h"
#include "shared.h"

namespace nb = nanobind;
using namespace pk5d::fft;

namespace {

using InArray = nb::ndarray<nb::ro, nb::c_contig, nb::device::cpu>;

// buffers allocated by the driver during one call; the result's is handed to NumPy, the others freed
struct Pool {
    std::vector<void*> bufs;
    ~Pool() {
        for (void* p : bufs) std::free(p);
    }
    void* release(void* p) {
        for (void*& q : bufs) {
            if (q == p) q = nullptr;
        }
        return p;
    }
};

void* pool_alloc(void* ctx, const std::vector<idx_t>& dims, bool is_double, bool is_complex) {
    size_t n = 1;
    for (idx_t d : dims) n *= static_cast<size_t>(d);
    const size_t bytes = n * (is_double ? 8 : 4) * (is_complex ? 2 : 1);
    void* p = std::malloc(bytes > 0 ? bytes : 16);
    if (p != nullptr) static_cast<Pool*>(ctx)->bufs.push_back(p);
    return p;
}

template <typename A>
void dtype_of(const A& a, bool& is_double, bool& is_complex) {
    const nb::dlpack::dtype t = a.dtype();
    if (t == nb::dtype<float>()) {
        is_double = false;
        is_complex = false;
    } else if (t == nb::dtype<double>()) {
        is_double = true;
        is_complex = false;
    } else if (t == nb::dtype<std::complex<float>>()) {
        is_double = false;
        is_complex = true;
    } else if (t == nb::dtype<std::complex<double>>()) {
        is_double = true;
        is_complex = true;
    } else {
        throw nb::type_error("the array must be float32, float64, complex64 or complex128");
    }
}

using OutArray = nb::ndarray<nb::c_contig, nb::device::cpu>;

// MATLAB's dimensions (y, x, z, d4, ...) of a C array's shape (..., d4, z, y, x): a 1-D array is a row (1, n),
// a 0-D one 1 x 1
std::vector<idx_t> matlab_dims(const std::vector<idx_t>& shape) {
    const size_t nd = shape.size();
    if (nd == 0) return {1, 1};
    if (nd == 1) return {1, shape[0]};
    std::vector<idx_t> D = {shape[nd - 2], shape[nd - 1]};
    for (size_t k = nd - 2; k-- > 0;) D.push_back(shape[k]);
    return D;
}

template <typename A>
Array view(const A& a) {
    Array x;
    x.data = const_cast<void*>(static_cast<const void*>(a.data()));
    dtype_of(a, x.is_double, x.is_complex);
    std::vector<idx_t> shape;
    for (size_t k = 0; k < a.ndim(); k++) shape.push_back(static_cast<idx_t>(a.shape(k)));
    x.dims = matlab_dims(shape);
    x.strides = native_strides(x.dims);
    return x;
}

// MATLAB's dimensions D as the C array's shape (..., D[3], D[2], D[0], D[1])
std::vector<size_t> py_shape(const std::vector<idx_t>& D) {
    std::vector<size_t> s;
    for (size_t k = D.size(); k-- > 2;) s.push_back(static_cast<size_t>(D[k]));
    s.push_back(static_cast<size_t>(D[0]));
    s.push_back(static_cast<size_t>(D.size() > 1 ? D[1] : 1));
    return s;
}

nb::ndarray<nb::numpy> wrap(void* p, const std::vector<size_t>& shape, bool is_double, bool is_complex) {
    nb::capsule owner(p, [](void* q) noexcept { std::free(q); });
    nb::dlpack::dtype t;
    if (is_complex) {
        t = is_double ? nb::dtype<std::complex<double>>() : nb::dtype<std::complex<float>>();
    } else {
        t = is_double ? nb::dtype<double>() : nb::dtype<float>();
    }
    return nb::ndarray<nb::numpy>(p, shape.size(), shape.data(), owner, nullptr, t);
}

nb::ndarray<nb::numpy> to_numpy(Pool& pool, const Array& y) {
    return wrap(pool.release(y.data), py_shape(y.dims), y.is_double, y.is_complex);
}

Options options(int threads, idx_t grain, int be) {
    Options o;
    o.threads = threads > 0 ? threads : 1;
    o.grain = grain > 0 ? grain : 0;
    o.backend = be;
    backend(be);  // throws when the exact backend is not loaded
    return o;
}

nb::ndarray<nb::numpy> py_fft(const InArray& a, idx_t n, int dim, int threads, idx_t grain, int be) {
    const Array x = view(a);
    if (dim < 0) throw nb::value_error("dim must be at least 0");
    const Options o = options(threads, grain, be);
    Pool pool;
    Array y;
    {
        nb::gil_scoped_release g;
        y = fft(x, n, dim, o, pool_alloc, &pool);
    }
    return to_numpy(pool, y);
}

nb::ndarray<nb::numpy> py_ifft(const InArray& a, idx_t n, int dim, bool symmetric, int threads, idx_t grain, int be) {
    const Array x = view(a);
    if (dim < 0) throw nb::value_error("dim must be at least 0");
    const Options o = options(threads, grain, be);
    Pool pool;
    Array y;
    {
        nb::gil_scoped_release g;
        y = ifft(x, n, dim, symmetric, o, pool_alloc, &pool);
    }
    return to_numpy(pool, y);
}

nb::ndarray<nb::numpy> py_fftn(const InArray& a, const std::vector<idx_t>& siz, int threads, idx_t grain, int be) {
    const Array x = view(a);
    const Options o = options(threads, grain, be);
    Pool pool;
    Array y;
    {
        nb::gil_scoped_release g;
        y = fftn(x, siz, o, pool_alloc, &pool);
    }
    return to_numpy(pool, y);
}

nb::ndarray<nb::numpy> py_ifftn(const InArray& a, const std::vector<idx_t>& siz, bool symmetric, int threads,
        idx_t grain, int be) {
    const Array x = view(a);
    const Options o = options(threads, grain, be);
    Pool pool;
    Array y;
    {
        nb::gil_scoped_release g;
        y = ifftn(x, siz, symmetric, o, pool_alloc, &pool);
    }
    return to_numpy(pool, y);
}

// MATLAB's half spectrum (the exact form): its C shape, forward, inverse, product
std::vector<size_t> py_half_shape(const std::vector<idx_t>& shape) { return py_shape(half_dims(matlab_dims(shape))); }

nb::ndarray<nb::numpy> py_half_forward(const InArray& a, int threads, idx_t grain, int be) {
    const Array x = view(a);
    if (x.is_complex) throw nb::type_error("the half spectrum takes real data");
    const Options o = options(threads, grain, be);
    Pool pool;
    Array h;
    {
        nb::gil_scoped_release g;
        h = half_forward(x, o, pool_alloc, &pool);
    }
    return to_numpy(pool, h);
}

nb::ndarray<nb::numpy> py_half_inverse(const OutArray& s, const std::vector<idx_t>& shape, int threads, idx_t grain,
        int be) {
    // the real result's MATLAB dimensions, and the spectrum as MATLAB's half of them
    Array h = view(s);
    if (!h.is_complex) throw nb::type_error("the half spectrum must be complex");
    const std::vector<idx_t> D = matlab_dims(shape);
    const std::vector<idx_t> H = half_dims(D);
    idx_t want = 1;
    for (idx_t v : H) want *= v;
    if (want != h.numel()) throw nb::value_error("the half spectrum's size does not match the shape");
    h.dims = H;
    h.strides = native_strides(H);
    const Options o = options(threads, grain, be);
    Pool pool;
    Array y;
    {
        nb::gil_scoped_release g;
        y = half_inverse(h, D, o, pool_alloc, &pool);
    }
    return to_numpy(pool, y);
}

void py_half_multiply(const OutArray& s, const InArray& t, bool conjugate, int threads) {
    Array h = view(s);
    const Array g = view(t);
    if (h.dims != g.dims) throw nb::value_error("the spectra differ in shape");
    nb::gil_scoped_release r;
    half_multiply(h, g, conjugate, threads > 0 ? threads : 1);
}

// the permissive (native) half spectrum: halved along x
nb::ndarray<nb::numpy> py_native_forward(const InArray& a, int threads) {
    bool dbl, cplx;
    dtype_of(a, dbl, cplx);
    if (cplx) throw nb::type_error("the half spectrum takes real data");
    std::vector<std::int64_t> shape;
    for (size_t k = 0; k < a.ndim(); k++) shape.push_back(static_cast<std::int64_t>(a.shape(k)));
    if (shape.empty()) shape.push_back(1);
    const std::vector<std::int64_t> hs = native_half_shape(shape);
    size_t n = 1;
    for (std::int64_t v : hs) n *= static_cast<size_t>(v);
    void* p = std::malloc(std::max<size_t>(16, n * (dbl ? 16 : 8)));
    if (p == nullptr) throw std::bad_alloc();
    try {
        nb::gil_scoped_release g;
        if (dbl) {
            native_forward(static_cast<const double*>(a.data()), shape, static_cast<double*>(p), threads);
        } else {
            native_forward(static_cast<const float*>(a.data()), shape, static_cast<float*>(p), threads);
        }
    } catch (...) {
        std::free(p);
        throw;
    }
    std::vector<size_t> s(hs.begin(), hs.end());
    return wrap(p, s, dbl, true);
}

nb::ndarray<nb::numpy> py_native_inverse(const OutArray& a, const std::vector<std::int64_t>& shape, int threads) {
    bool dbl, cplx;
    dtype_of(a, dbl, cplx);
    if (!cplx) throw nb::type_error("the half spectrum must be complex");
    std::vector<std::int64_t> sh(shape);
    if (sh.empty()) sh.push_back(1);
    const std::vector<std::int64_t> hs = native_half_shape(sh);
    size_t nh = 1, n = 1;
    for (std::int64_t v : hs) nh *= static_cast<size_t>(v);
    for (std::int64_t v : sh) n *= static_cast<size_t>(v);
    size_t have = 1;
    for (size_t k = 0; k < a.ndim(); k++) have *= a.shape(k);
    if (have != nh) throw nb::value_error("the half spectrum's size does not match the shape");
    void* p = std::malloc(std::max<size_t>(16, n * (dbl ? 8 : 4)));
    if (p == nullptr) throw std::bad_alloc();
    try {
        nb::gil_scoped_release g;
        if (dbl) {
            native_inverse(static_cast<double*>(a.data()), sh, static_cast<double*>(p), threads);
        } else {
            native_inverse(static_cast<float*>(a.data()), sh, static_cast<float*>(p), threads);
        }
    } catch (...) {
        std::free(p);
        throw;
    }
    std::vector<size_t> s(sh.begin(), sh.end());
    return wrap(p, s, dbl, false);
}

void py_load_exact(const std::string& libf, const std::string& libd) {
    load_exact(libf, libd);
}

void py_unload_exact() {
    unload_exact();
}

std::tuple<std::string, std::string, std::string, std::string> py_backend_info(int be) {
    const Backend& b = backend(be);
    return {b.name, b.where, b.api[0].version, b.api[1].version};
}

}  // namespace

NB_MODULE(fft, m) {
    m.doc() = "MATLAB R2024b's fft, ifft, fftn and ifftn and half-spectrum transforms for (..., z, y, x) C arrays";
    m.def("fft", &py_fft, nb::arg("x").noconvert(), nb::arg("n"), nb::arg("dim"), nb::arg("threads"), nb::arg("grain"),
        nb::arg("backend"), "fft(x, n, dim): n < 0 for the size along MATLAB's dimension dim + 1");
    m.def("ifft", &py_ifft, nb::arg("x").noconvert(), nb::arg("n"), nb::arg("dim"), nb::arg("symmetric"),
        nb::arg("threads"), nb::arg("grain"), nb::arg("backend"));
    m.def("fftn", &py_fftn, nb::arg("x").noconvert(), nb::arg("siz"), nb::arg("threads"), nb::arg("grain"),
        nb::arg("backend"), "siz in MATLAB's order (y, x, z, ...)");
    m.def("ifftn", &py_ifftn, nb::arg("x").noconvert(), nb::arg("siz"), nb::arg("symmetric"), nb::arg("threads"),
        nb::arg("grain"), nb::arg("backend"));
    m.def("half_shape", &py_half_shape, nb::arg("shape"), "the C shape of MATLAB's half spectrum of a real array");
    m.def("half_forward", &py_half_forward, nb::arg("x").noconvert(), nb::arg("threads"), nb::arg("grain"),
        nb::arg("backend"));
    m.def("half_inverse", &py_half_inverse, nb::arg("h").noconvert(), nb::arg("shape"), nb::arg("threads"),
        nb::arg("grain"), nb::arg("backend"), "real(ifftn) of the spectrum whose half is h (overwritten)");
    m.def("half_multiply", &py_half_multiply, nb::arg("h").noconvert(), nb::arg("g").noconvert(),
        nb::arg("conjugate"), nb::arg("threads"), "h .* g (conj(g) .* h) in place, MATLAB's complex product");
    m.def("native_forward", &py_native_forward, nb::arg("x").noconvert(), nb::arg("threads"));
    m.def("native_inverse", &py_native_inverse, nb::arg("h").noconvert(), nb::arg("shape"), nb::arg("threads"));
    m.def("load_exact", &py_load_exact, nb::arg("libf"), nb::arg("libd"),
        "load the exact backend: FFTW libraries for single (libf) and double (libd) precision");
    m.def("unload_exact", &py_unload_exact);
    m.def("exact_loaded", &exact_loaded);
    m.def("backend_info", &py_backend_info, nb::arg("backend"), "(name, where, single version, double version)");
    m.def("plan_string", &plan_string, nb::arg("backend"), nb::arg("is_double"), nb::arg("kind"), nb::arg("n"),
        nb::arg("howmany"), "the description of the plan used for kind r2c, c2c-1, c2c+1 or c2r");
    m.def("line_alone", &line_alone, nb::arg("j"), nb::arg("N"), nb::arg("n"), nb::arg("seg"), nb::arg("grain"),
        "whether MATLAB transforms line j of a stage of N lines of length n alone");
    m.attr("GRAIN_PARALLEL") = kGrainParallel;
    m.attr("PERMISSIVE") = static_cast<int>(kPermissive);
    m.attr("EXACT") = static_cast<int>(kExact);
    // the process's FFT state (backends, planner lock, plan cache) for the other extension modules (shared.h)
    m.attr("_shared_state") = nb::capsule(&own_state(), kSharedCapsule);
    def_shared_state_hooks(m);
}
