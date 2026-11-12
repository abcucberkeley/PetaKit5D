// Checks of the per-region thread binding (K08's ops_bind.hpp), for the tests.
//
// ops_team_cpus(num_threads, bind): for each thread of a team (by thread number), the CPUs of
//   its affinity mask inside the region (one CPU when bound); empty lists outside Linux.
#include "ops_common.hpp"

#include <nanobind/stl/vector.h>

#include <vector>

namespace pk5d::ops {
namespace {

std::vector<int> my_cpus() {
    std::vector<int> out;
#if defined(__linux__)
    cpu_set_t s;
    CPU_ZERO(&s);
    if (pthread_getaffinity_np(pthread_self(), sizeof(s), &s) == 0)
        for (int c = 0; c < CPU_SETSIZE; c++)
            if (CPU_ISSET(c, &s)) out.push_back(c);
#endif
    return out;
}

std::vector<std::vector<int>> team_cpus(int num_threads, bool bind) {
    check_threads(num_threads);
    std::vector<std::vector<int>> out(static_cast<size_t>(num_threads));
    nb::gil_scoped_release release;
    parallel(num_threads, bind, [&](int k, int) { out[static_cast<size_t>(k)] = my_cpus(); });
    return out;
}

}  // namespace
}  // namespace pk5d::ops

PK5D_BINDINGS(bind_ops_threads) {
    namespace nb = nanobind;
    m.def("ops_team_cpus", &pk5d::ops::team_cpus, nb::arg("num_threads"), nb::arg("bind"),
          "The CPUs of each team thread's affinity mask inside a region of the ops kernels.");
}
