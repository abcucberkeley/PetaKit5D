// Helpers shared by the kernel groups: array types, thread counts, binding registration.
//
// Kernels take (z, y, x) C-contiguous arrays (cvolume, cimage: x fastest, as in a tiff page;
// see python/docs/conventions.md), check class and layout at the binding,
// release the interpreter lock around their loops and take the thread count as an argument
// (never omp_set_num_threads: no global state). Array arguments are bound with
// nb::arg("name").noconvert(): without it nanobind silently copies an array of another
// layout or class into the one asked for, a hidden copy of a whole volume.
#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#elif defined(__linux__)
#include <sched.h>
#include <unistd.h>
#else
#include <unistd.h>
#endif

namespace pk5d {

namespace nb = nanobind;

// C-contiguous CPU arrays, the convention: a volume (nz, ny, nx), an image (ny, nx), a vector.
// Element (z, y, x) of a cvolume v is v.data()[x + nx * (y + ny * z)] (v.shape(0) is nz);
// MATLAB's linear index of the same voxel is y + ny * (x + nx * z).
template <typename T>
using cvolume = nb::ndarray<T, nb::ndim<3>, nb::c_contig, nb::device::cpu>;
template <typename T>
using cimage = nb::ndarray<T, nb::ndim<2>, nb::c_contig, nb::device::cpu>;
template <typename T>
using vector1d = nb::ndarray<T, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

// Fortran-ordered arrays in MATLAB's layout (y, x, z): kept only until the last module has moved
// to cvolume / cimage; new code does not use them.
template <typename T>
using fvolume = nb::ndarray<T, nb::ndim<3>, nb::f_contig, nb::device::cpu>;
template <typename T>
using fimage = nb::ndarray<T, nb::ndim<2>, nb::f_contig, nb::device::cpu>;

// CPUs this process may run on: the affinity mask (Slurm's -c, taskset, cgroup cpusets).
// cgroup CPU quotas are applied on the Python side (petakit5d.core).
inline int affinity_cpu_count() {
#if defined(__linux__)
    for (int n = 1024; n <= (1 << 20); n *= 2) {
        cpu_set_t *set = CPU_ALLOC(n);
        if (set == nullptr) break;
        size_t size = CPU_ALLOC_SIZE(n);
        CPU_ZERO_S(size, set);
        if (sched_getaffinity(0, size, set) == 0) {
            int count = CPU_COUNT_S(size, set);
            CPU_FREE(set);
            return count > 0 ? count : 1;
        }
        CPU_FREE(set);
    }
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? static_cast<int>(n) : 1;
#elif defined(_WIN32)
    DWORD_PTR process_mask = 0, system_mask = 0;
    if (GetProcessAffinityMask(GetCurrentProcess(), &process_mask, &system_mask) && process_mask) {
        int count = 0;
        for (; process_mask; process_mask &= process_mask - 1) ++count;
        return count;
    }
    DWORD n = GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
    return n > 0 ? static_cast<int>(n) : 1;
#else
    // macOS has no affinity masks
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? static_cast<int>(n) : 1;
#endif
}

// A kernel's thread count: the caller's, which must be positive
inline int check_threads(int num_threads) {
    if (num_threads < 1)
        throw std::invalid_argument("num_threads must be at least 1, got " +
                                    std::to_string(num_threads));
    return num_threads;
}

// Bindings in several source files of one group: each file writes
//   PK5D_BINDINGS(bind_something) { m.def(...); }
// and the file with NB_MODULE calls pk5d::run_bindings(m).
using binder_fn = void (*)(nb::module_ &);

inline std::vector<binder_fn> &binders() {
    static std::vector<binder_fn> list;
    return list;
}

struct register_binder {
    explicit register_binder(binder_fn fn) { binders().push_back(fn); }
};

inline void run_bindings(nb::module_ &m) {
    for (binder_fn fn : binders()) fn(m);
}

}  // namespace pk5d

#define PK5D_BINDINGS(name)                                      \
    static void name(::nanobind::module_ &);                     \
    static ::pk5d::register_binder name##_registered_(name);     \
    static void name(::nanobind::module_ &m)
