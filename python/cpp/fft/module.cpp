// petakit5d._ext.fft: MATLAB's fft, ifft, fftn and ifftn (driver.h) for Fortran-ordered NumPy arrays.
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <complex>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "backend.h"
#include "driver.h"

namespace nb = nanobind;
using namespace pk5d::fft;

namespace {

using InArray = nb::ndarray<nb::ro, nb::f_contig, nb::device::cpu>;

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

Array view(const InArray& a) {
    Array x;
    x.data = const_cast<void*>(a.data());
    const nb::dlpack::dtype t = a.dtype();
    if (t == nb::dtype<float>()) {
        x.is_double = false;
        x.is_complex = false;
    } else if (t == nb::dtype<double>()) {
        x.is_double = true;
        x.is_complex = false;
    } else if (t == nb::dtype<std::complex<float>>()) {
        x.is_double = false;
        x.is_complex = true;
    } else if (t == nb::dtype<std::complex<double>>()) {
        x.is_double = true;
        x.is_complex = true;
    } else {
        throw nb::type_error("the array must be float32, float64, complex64 or complex128");
    }
    for (size_t k = 0; k < a.ndim(); k++) x.dims.push_back(static_cast<idx_t>(a.shape(k)));
    if (x.dims.empty()) x.dims.push_back(1);
    return x;
}

nb::ndarray<nb::numpy> to_numpy(Pool& pool, const Array& y) {
    void* p = pool.release(y.data);
    nb::capsule owner(p, [](void* q) noexcept { std::free(q); });
    std::vector<size_t> shape;
    std::vector<int64_t> strides;
    int64_t s = 1;
    for (idx_t d : y.dims) {
        shape.push_back(static_cast<size_t>(d));
        strides.push_back(s);
        s *= d;
    }
    nb::dlpack::dtype t;
    if (y.is_complex) {
        t = y.is_double ? nb::dtype<std::complex<double>>() : nb::dtype<std::complex<float>>();
    } else {
        t = y.is_double ? nb::dtype<double>() : nb::dtype<float>();
    }
    return nb::ndarray<nb::numpy>(p, shape.size(), shape.data(), owner, strides.data(), t);
}

Options options(int threads, idx_t grain) {
    Options o;
    o.threads = threads > 0 ? threads : 1;
    o.grain = grain > 0 ? grain : 0;
    return o;
}

nb::ndarray<nb::numpy> py_fft(const InArray& a, idx_t n, int dim, int threads, idx_t grain) {
    const Array x = view(a);
    if (dim < 0) throw nb::value_error("dim must be at least 0");
    Pool pool;
    Array y;
    {
        nb::gil_scoped_release g;
        y = fft(x, n, dim, options(threads, grain), pool_alloc, &pool);
    }
    return to_numpy(pool, y);
}

nb::ndarray<nb::numpy> py_ifft(const InArray& a, idx_t n, int dim, bool symmetric, int threads, idx_t grain) {
    const Array x = view(a);
    if (dim < 0) throw nb::value_error("dim must be at least 0");
    Pool pool;
    Array y;
    {
        nb::gil_scoped_release g;
        y = ifft(x, n, dim, symmetric, options(threads, grain), pool_alloc, &pool);
    }
    return to_numpy(pool, y);
}

nb::ndarray<nb::numpy> py_fftn(const InArray& a, const std::vector<idx_t>& siz, int threads, idx_t grain) {
    const Array x = view(a);
    Pool pool;
    Array y;
    {
        nb::gil_scoped_release g;
        y = fftn(x, siz, options(threads, grain), pool_alloc, &pool);
    }
    return to_numpy(pool, y);
}

nb::ndarray<nb::numpy> py_ifftn(const InArray& a, const std::vector<idx_t>& siz, bool symmetric, int threads, idx_t grain) {
    const Array x = view(a);
    Pool pool;
    Array y;
    {
        nb::gil_scoped_release g;
        y = ifftn(x, siz, symmetric, options(threads, grain), pool_alloc, &pool);
    }
    return to_numpy(pool, y);
}

void py_set_backend(const std::string& libf, const std::string& libd) {
    clear_plans();
    set_backend(libf, libd);
}

std::tuple<std::string, std::string, std::string> py_backend_info() {
    const Backend& b = backend();
    return {b.name, b.api[0].version, b.api[1].version};
}

}  // namespace

NB_MODULE(fft, m) {
    m.doc() = "MATLAB R2024b's fft, ifft, fftn and ifftn on FFTW, bit for bit (Fortran-ordered arrays)";
    m.def("fft", &py_fft, nb::arg("x").noconvert(), nb::arg("n"), nb::arg("dim"), nb::arg("threads"), nb::arg("grain"),
        "fft(x, n, dim): n < 0 for the size along dim (0-based)");
    m.def("ifft", &py_ifft, nb::arg("x").noconvert(), nb::arg("n"), nb::arg("dim"), nb::arg("symmetric"), nb::arg("threads"),
        nb::arg("grain"));
    m.def("fftn", &py_fftn, nb::arg("x").noconvert(), nb::arg("siz"), nb::arg("threads"), nb::arg("grain"));
    m.def("ifftn", &py_ifftn, nb::arg("x").noconvert(), nb::arg("siz"), nb::arg("symmetric"), nb::arg("threads"), nb::arg("grain"));
    m.def("set_backend", &py_set_backend, nb::arg("libf"), nb::arg("libd"),
        "use the FFTW libraries at libf (single) and libd (double); empty strings: the built-in FFTW");
    m.def("backend_info", &py_backend_info, "(name, single version, double version) of the FFTW in use");
    m.def("plan_string", &plan_string, nb::arg("is_double"), nb::arg("kind"), nb::arg("n"), nb::arg("howmany"),
        "FFTW's description of the plan used for kind r2c, c2c-1, c2c+1 or c2r");
    m.def("line_alone", &line_alone, nb::arg("j"), nb::arg("N"), nb::arg("n"), nb::arg("seg"), nb::arg("grain"),
        "whether MATLAB transforms line j of a stage of N lines of length n alone");
    m.attr("GRAIN_PARALLEL") = kGrainParallel;
}
