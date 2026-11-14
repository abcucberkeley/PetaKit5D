// Linked against a candidate OpenMP runtime before the modules share it (cmake/Pk5dOpenMP.cmake):
// the constructs and calls the kernels use, so that the executable links only if the runtime
// defines every entry point the compiler emits for them.
#include <omp.h>

#include <cstddef>
#include <cstdio>

static int cell[4][256];

int main(int argc, char **) {
    const int n = 64 + argc;
    double top = -1.0, stamp = 0.0;
    float low = 1e30f;
    bool any = false;
#pragma omp parallel num_threads(2)
    {
        const int t = omp_get_thread_num() + omp_get_num_threads();
#pragma omp for schedule(static) nowait
        for (int i = 0; i < n; i++)
            cell[0][i] = t;
#pragma omp for schedule(dynamic)
        for (int i = 0; i < n; i++)
            cell[1][i] = t;
#pragma omp for schedule(dynamic) collapse(2)
        for (std::size_t j = 2; j < 4; j++)
            for (std::size_t i = 0; i < static_cast<std::size_t>(n); i++)
                cell[j][i] = t;
#pragma omp single
        stamp = omp_get_wtime();
#pragma omp barrier
    }
#pragma omp parallel for num_threads(2) schedule(dynamic, 4) reduction(max : top) \
    reduction(min : low) reduction(| : any)
    for (std::size_t i = 0; i < static_cast<std::size_t>(n); i++) {
        top = top > cell[0][i] ? top : cell[0][i];
        low = low < cell[1][i] ? low : static_cast<float>(cell[1][i]);
        any = any | (cell[2][i] == cell[3][i]);
    }
    std::printf("%d %d %g %g %g %d\n", omp_get_max_threads(), omp_get_num_procs(), stamp, top,
                static_cast<double>(low), static_cast<int>(any));
    return 0;
}
