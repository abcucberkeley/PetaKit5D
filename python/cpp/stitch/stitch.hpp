// Shared by the stitching kernels: array types and checks of the bindings.
//
// Arrays reach the kernels as flat C-contiguous views of (z, y, x) volumes and (t, z, y, x)
// tile stacks (x fastest) with their sizes as separate arguments; the Python side
// (petakit5d.kernels.stitch) checks classes, sizes and layout first, as the mex gateways did.
// The checks here only keep the kernels inside their arrays.
//
// The mex files index MATLAB's (y, x, z) arrays with "x" for MATLAB's first dimension (y). The
// kernels here keep each voxel's arithmetic and every order that decides bits (tiles in order,
// the integral image's passes along y, then x, then z, the trilinear interpolation along y,
// then x, then z); only the order in which voxels are visited follows the C layout.
#pragma once

#include <pk5d/common.hpp>

#include <nanobind/stl/vector.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace pk5d_stitch {

namespace nb = nanobind;

template <typename T>
using in1d = nb::ndarray<const T, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
template <typename T>
using out1d = nb::ndarray<T, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

inline void need(bool ok, const std::string &what) {
    if (!ok) throw std::invalid_argument(what);
}

inline uint64_t prod(const std::vector<uint64_t> &d) {
    uint64_t n = 1;
    for (uint64_t v : d) n *= v;
    return n;
}

inline void need_dims(const std::vector<uint64_t> &d, size_t n, const char *name) {
    need(d.size() == n, std::string(name) + " must have " + std::to_string(n) + " sizes");
}

template <typename A>
inline void need_size(const A &a, uint64_t n, const char *name) {
    need((uint64_t)a.shape(0) == n, std::string(name) + " has " + std::to_string(a.shape(0)) +
                                        " elements, expected " + std::to_string(n));
}

}  // namespace pk5d_stitch
