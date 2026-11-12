// Thread binding for one parallel region of the K08 kernels, never for the process.
//
// libgomp ignores proc_bind(...) when no place list exists (OMP_PLACES / OMP_PROC_BIND unset),
// and setting those would bind every thread of the process (and the Python thread to one CPU).
// Instead each thread of a team pins itself, for the region only, to one CPU of the calling
// thread's affinity mask and restores its own mask when the region ends: thread k of nk takes
// CPU k * m / nk of the mask's m CPUs, listed alternately from each processor package, so that
// small teams use every socket. Nothing is pinned when the user set OpenMP's binding variables
// (theirs wins), when the mask has fewer CPUs than the team, for one thread, or outside Linux.
#pragma once

#include <pk5d/common.hpp>

#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>

#if defined(__linux__)
#include <pthread.h>
#include <sched.h>
#endif

namespace pk5d::ops {

#if defined(__linux__)

// the processor package of a CPU (sysfs, read once; -1 if unknown)
inline int cpu_package(int cpu) {
    static const std::vector<int> table = [] {
        std::vector<int> t(CPU_SETSIZE, -1);
        for (int c = 0; c < CPU_SETSIZE; c++) {
            char path[96];
            std::snprintf(path, sizeof(path),
                          "/sys/devices/system/cpu/cpu%d/topology/physical_package_id", c);
            FILE *f = std::fopen(path, "r");
            if (f == nullptr) continue;
            int p = -1;
            if (std::fscanf(f, "%d", &p) == 1) t[static_cast<size_t>(c)] = p;
            std::fclose(f);
        }
        return t;
    }();
    return cpu >= 0 && cpu < CPU_SETSIZE ? table[static_cast<size_t>(cpu)] : -1;
}

// whether the user asked OpenMP to bind threads (then we leave binding to it)
inline bool user_binding() {
    static const bool set = std::getenv("OMP_PROC_BIND") != nullptr ||
                            std::getenv("OMP_PLACES") != nullptr ||
                            std::getenv("GOMP_CPU_AFFINITY") != nullptr;
    return set;
}

// the CPUs of the calling thread's mask, alternating between packages
inline std::vector<int> spread_cpus() {
    cpu_set_t s;
    CPU_ZERO(&s);
    if (pthread_getaffinity_np(pthread_self(), sizeof(s), &s) != 0) return {};
    std::map<int, std::vector<int>> by_package;
    for (int c = 0; c < CPU_SETSIZE; c++)
        if (CPU_ISSET(c, &s)) by_package[cpu_package(c)].push_back(c);
    std::vector<int> out;
    for (size_t i = 0;; i++) {
        bool any = false;
        for (auto &kv : by_package) {
            if (i < kv.second.size()) {
                out.push_back(kv.second[i]);
                any = true;
            }
        }
        if (!any) break;
    }
    return out;
}

// the CPUs a team of nt threads pins to, or none
struct bind_plan {
    std::vector<int> cpus;
    bool on = false;
    bind_plan(int nt, bool bind) {
        if (!bind || nt < 2 || user_binding()) return;
        cpus = spread_cpus();
        on = static_cast<int>(cpus.size()) >= nt;
    }
};

// pins the current thread (thread k of nk) for its lifetime, then restores its mask
class bind_guard {
public:
    bind_guard(const bind_plan &p, int k, int nk) {
        if (!p.on || nk < 2) return;
        CPU_ZERO(&saved_);
        if (pthread_getaffinity_np(pthread_self(), sizeof(saved_), &saved_) != 0) return;
        const size_t m = p.cpus.size();
        cpu_set_t one;
        CPU_ZERO(&one);
        CPU_SET(p.cpus[static_cast<size_t>(k) * m / static_cast<size_t>(nk)], &one);
        active_ = pthread_setaffinity_np(pthread_self(), sizeof(one), &one) == 0;
    }
    ~bind_guard() {
        if (active_) pthread_setaffinity_np(pthread_self(), sizeof(saved_), &saved_);
    }
    bind_guard(const bind_guard &) = delete;
    bind_guard &operator=(const bind_guard &) = delete;

private:
    cpu_set_t saved_;
    bool active_ = false;
};

#else

struct bind_plan {
    bind_plan(int, bool) {}
};

class bind_guard {
public:
    bind_guard(const bind_plan &, int, int) {}
};

#endif

// fn(k, nk) on a team of nt threads (thread k of nk), each pinned for the region when bind
template <typename Fn>
void parallel(int nt, bool bind, Fn &&fn) {
#ifdef _OPENMP
    if (nt < 2) {
        fn(0, 1);
        return;
    }
    const bind_plan plan(nt, bind);
#pragma omp parallel num_threads(nt)
    {
        const int k = omp_get_thread_num(), nk = omp_get_num_threads();
        const bind_guard g(plan, k, nk);
        fn(k, nk);
    }
#else
    (void)nt;
    (void)bind;
    fn(0, 1);
#endif
}

// the part [lo, hi) of n items for thread k of nk (static, contiguous)
inline void share(std::ptrdiff_t n, int k, int nk, std::ptrdiff_t &lo, std::ptrdiff_t &hi) {
    lo = n * k / nk;
    hi = n * (k + 1) / nk;
}

}  // namespace pk5d::ops
