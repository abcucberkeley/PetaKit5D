// petakit5d._ext.core: what the build did (compiler, flags, OpenMP) and checks that the
// kernels compute as intended (no fused multiply-adds, subnormals kept, thread counts).
// The element-wise primitives of ops/ (K08) register into this module (PK5D_BINDINGS).
#include <pk5d/common.hpp>

#include <nanobind/stl/string.h>

#include <cfloat>
#include <string>

namespace nb = nanobind;

#ifndef PK5D_BUILD_FLAGS
#define PK5D_BUILD_FLAGS ""
#endif
#ifndef PK5D_BUILD_CONFIG
#define PK5D_BUILD_CONFIG ""
#endif
#ifndef PK5D_OPENMP_DESC
#define PK5D_OPENMP_DESC "none"
#endif
#ifndef PK5D_OPENMP_SHARED
#define PK5D_OPENMP_SHARED ""
#endif

namespace {

std::string compiler_name() {
#if defined(__clang__)
    return std::string("clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("gcc ") + __VERSION__;
#elif defined(_MSC_VER)
    return "msvc " + std::to_string(_MSC_FULL_VER);
#else
    return "unknown";
#endif
}

// a * b + c in a loop the compiler may vectorize, with inputs it cannot fold: a fused
// multiply-add keeps the exact product (1 - 2^-54, 1 - 2^-26 for single), a separate
// multiply rounds it to 1, so the sum is 0 only without contraction
template <typename T>
NB_NOINLINE bool muladd_fused(T a0, T b0, T c0) {
    constexpr int n = 64;
    volatile T va = a0, vb = b0, vc = c0;
    T a[n], b[n], c[n], r[n];
    for (int i = 0; i < n; ++i) {
        a[i] = va;
        b[i] = vb;
        c[i] = vc;
    }
    for (int i = 0; i < n; ++i) r[i] = a[i] * b[i] + c[i];
    T x = va, y = vb, z = vc;
    volatile T s0 = x * y + z;
    bool fused = s0 != T(0);
    for (int i = 0; i < n; ++i) fused = fused || r[i] != T(0);
    return fused;
}

bool fma_contracted() {
    return muladd_fused<double>(1.0 + 0x1p-27, 1.0 - 0x1p-27, -1.0) ||
           muladd_fused<float>(1.0f + 0x1p-13f, 1.0f - 0x1p-13f, -1.0f);
}

// flush-to-zero or denormals-are-zero (set by fast-math start-up code of any library)
bool subnormals_flushed() {
    volatile float tiny = FLT_MIN;
    volatile float half = 0.5f;
    volatile float sub = tiny * half;
    volatile double dtiny = DBL_MIN;
    volatile double dsub = dtiny * 0.5;
    return sub == 0.0f || dsub == 0.0 || !(sub * 2.0f == tiny);
}

int omp_max_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

int omp_num_procs() {
#ifdef _OPENMP
    return omp_get_num_procs();
#else
    return pk5d::affinity_cpu_count();
#endif
}

// threads that ran a parallel region asked for num_threads
int parallel_threads(int num_threads) {
    pk5d::check_threads(num_threads);
    int seen = 1;
    {
        nb::gil_scoped_release release;
#ifdef _OPENMP
#pragma omp parallel num_threads(num_threads)
        {
#pragma omp single
            seen = omp_get_num_threads();
        }
#endif
    }
    return seen;
}

}  // namespace

NB_MODULE(core, m) {
    m.doc() = "Build information and run-time checks of the petakit5d kernels";

    m.def("build_info", []() {
        nb::dict d;
        d["compiler"] = compiler_name();
        d["flags"] = std::string(PK5D_BUILD_FLAGS);
        d["config"] = std::string(PK5D_BUILD_CONFIG);
        d["openmp"] = std::string(PK5D_OPENMP_DESC);
        // Linux: the I/O packages' libgomp the modules share (its SONAME), "off", or empty
        d["openmp_shared"] = std::string(PK5D_OPENMP_SHARED);
#ifdef _OPENMP
        d["openmp_version"] = _OPENMP;
#else
        d["openmp_version"] = 0;
#endif
#if defined(__AVX2__)
        d["avx2"] = true;
#else
        d["avx2"] = false;
#endif
        d["nanobind"] = std::to_string(NB_VERSION_MAJOR) + "." + std::to_string(NB_VERSION_MINOR) +
                        "." + std::to_string(NB_VERSION_PATCH);
        return d;
    }, "Compiler, flags and OpenMP runtime the extension modules were built with.");
    m.def("affinity_cpu_count", &pk5d::affinity_cpu_count,
          "CPUs this process may run on (affinity mask; all online CPUs on macOS).");
    m.def("omp_max_threads", &omp_max_threads,
          "OpenMP's default team size (OMP_NUM_THREADS, else the CPUs it sees).");
    m.def("omp_num_procs", &omp_num_procs, "CPUs the OpenMP runtime sees.");
    m.def("parallel_threads", &parallel_threads, nb::arg("num_threads"),
          "Threads that ran a parallel region asked for num_threads.");
    m.def("fma_contracted", &fma_contracted,
          "True if a * b + c was compiled as a fused multiply-add (it must not be).");
    m.def("subnormals_flushed", &subnormals_flushed,
          "True if subnormal numbers are flushed to zero in this thread (they must not be).");

    pk5d::run_bindings(m);
}
