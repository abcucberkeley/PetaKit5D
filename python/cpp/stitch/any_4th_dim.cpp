// Whether the major tiles of a batch have data (a nonzero value; NaN counts) at every voxel, and
// the box of the voxels where none has: any_4th_dim_mex.cpp. The answer and the box (minima and
// maxima over the voxels) do not depend on the order the voxels are visited in: each thread
// keeps its own bounds, combined at the end (a min/max reduction), and the scan stops early
// once a voxel without data is found.
//
// tiles: (t, z, y, x), C order; split by (z, y) lines.
#include "stitch.hpp"

#include <algorithm>
#include <utility>

#include <nanobind/stl/pair.h>

namespace {

using namespace pk5d_stitch;

// voxel i (offset in a tile) has no data in any major tile
template <typename T>
inline bool no_data(const T *f, const bool *major, uint64_t nT, uint64_t nZYX, uint64_t i) {
    for (uint64_t t = 0; t < nT; t++) {
        if (major[t] && f[t * nZYX + i]) return false;
    }
    return true;
}

template <typename T>
bool all_have_data(const T *f, const bool *major, const uint64_t *dims, int nthread) {
    const uint64_t nT = dims[0], nZ = dims[1], nY = dims[2], nX = dims[3];
    const uint64_t nZYX = nZ * nY * nX;
    const int64_t nLines = (int64_t)(nZ * nY);
    volatile bool found = false;
#pragma omp parallel for schedule(static) num_threads(nthread) if (nthread > 1) shared(found)
    for (int64_t l = 0; l < nLines; l++) {
        if (found) continue;
        const uint64_t line = (uint64_t)l * nX;
        for (uint64_t x = 0; x < nX; x++) {
            if (no_data(f, major, nT, nZYX, line + x)) {
                found = true;
                break;
            }
        }
    }
    return !found;
}

// 0-based (z0, y0, x0, z1, y1, x1), end-exclusive, of the voxels without data
template <typename T>
std::vector<uint64_t> no_data_box(const T *f, const bool *major, const uint64_t *dims,
                                  int nthread) {
    const uint64_t nT = dims[0], nZ = dims[1], nY = dims[2], nX = dims[3];
    const uint64_t nZYX = nZ * nY * nX;
    const int64_t nLines = (int64_t)(nZ * nY);
    uint64_t z0 = nZ, y0 = nY, x0 = nX, z1 = 0, y1 = 0, x1 = 0;
#pragma omp parallel for schedule(static) num_threads(nthread) if (nthread > 1) \
    reduction(min : z0, y0, x0) reduction(max : z1, y1, x1)
    for (int64_t l = 0; l < nLines; l++) {
        const uint64_t z = (uint64_t)l / nY, y = (uint64_t)l % nY;
        const uint64_t line = (uint64_t)l * nX;
        for (uint64_t x = 0; x < nX; x++) {
            if (no_data(f, major, nT, nZYX, line + x)) {
                z0 = std::min(z0, z);
                y0 = std::min(y0, y);
                x0 = std::min(x0, x);
                z1 = std::max(z1, z + 1);
                y1 = std::max(y1, y + 1);
                x1 = std::max(x1, x + 1);
            }
        }
    }
    return {z0, y0, x0, z1, y1, x1};
}

template <typename T>
void bind_any(nb::module_ &m) {
    m.def(
        "any_4th_dim",
        [](in1d<T> f, in1d<bool> major, std::vector<uint64_t> dims, int num_threads) {
            pk5d::check_threads(num_threads);
            need_dims(dims, 4, "dims");
            need_size(f, prod(dims), "tiles");
            need_size(major, dims[0], "major");
            bool any_major = false;
            for (uint64_t a = 0; a < dims[0]; a++) any_major = any_major || major.data()[a];
            need(any_major, "no major tile");
            std::vector<uint64_t> box(6, 0);
            bool valid;
            {
                nb::gil_scoped_release release;
                valid = all_have_data<T>(f.data(), major.data(), dims.data(), num_threads);
                if (!valid) box = no_data_box<T>(f.data(), major.data(), dims.data(), num_threads);
            }
            return std::make_pair(valid, box);
        },
        nb::arg("f").noconvert(), nb::arg("major").noconvert(), nb::arg("dims"),
        nb::arg("num_threads"),
        "(valid, box): whether every voxel has data in a major tile, else the 0-based, "
        "end-exclusive (z0, y0, x0, z1, y1, x1) of the voxels without (zeros when valid).");
}

}  // namespace

PK5D_BINDINGS(bind_any_4th_dim) {
    bind_any<uint16_t>(m);
    bind_any<float>(m);
}
