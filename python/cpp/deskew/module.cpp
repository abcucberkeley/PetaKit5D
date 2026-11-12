// petakit5d._ext.deskew: the deskew/rotate kernels of the three mex files, for C-ordered
// (z, y, x) volumes, and (the *_ml functions, the mex's loops in ml/) for MATLAB's layout, the
// C-ordered (z, x, y) arrays; swap_axes turns one layout into the other
//   volume_deskew_rotate_warp_mex                       -> warp
//   skewed_space_interp_volume_deskew_rotate_warp_mex   -> skewed_warp
//   skewed_space_interp_defined_stepsize_mex            -> interp_defined_stepsize
// The argument checks of the mex gateways are in petakit5d.kernels.deskew (with MATLAB's error
// ids); these bindings check class and layout, take the output array (zeros) from the caller
// and release the interpreter lock around the loops. Also built without AVX2 as
// petakit5d._ext.deskew_scalar on x86-64, for the tests.
#include <pk5d/common.hpp>

#include <nanobind/stl/string.h>

#include <cstdint>
#include <stdexcept>
#include <string>

#include "interp.hpp"
#include "ml/interp.hpp"
#include "ml/skewed_warp.hpp"
#include "ml/warp.hpp"
#include "skewed_warp.hpp"
#include "transpose.hpp"
#include "warp.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace {

// C-contiguous (z, y, x) volumes: element (z, y, x) is data[x + nx * (y + ny * z)]
template <typename T>
using in_volume = pk5d::cvolume<const T>;
template <typename T>
using out_volume = pk5d::cvolume<T>;
using matrix4 = nb::ndarray<const double, nb::shape<4, 4>, nb::c_contig, nb::device::cpu>;

// tmat[i][j] = M(i + 1, j + 1) of MATLAB's 4 x 4 matrix (the mex reads it column-major)
void read_tmat(const matrix4& m, double (&tmat)[4][4]) {
    const double* p = m.data();
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) tmat[i][j] = p[i * 4 + j];
}

template <typename A>
void shape3(const A& a, uint64_t (&sz)[3]) {
    for (int i = 0; i < 3; i++) sz[i] = (uint64_t)a.shape(i);
}

// the output lines read input rows y = floor(start_y) + i, i < line_len
void check_line_len(uint64_t line_len, float start_y, uint64_t ny, uint64_t oy) {
    const uint64_t y0 = (uint64_t)floor(start_y);
    if (line_len > oy || (line_len > 0 && (y0 >= ny || line_len > ny - y0)))
        throw std::invalid_argument("line_len does not fit the input and output rows");
}

template <typename T, typename TO>
void warp(const in_volume<T>& vol, const matrix4& m, out_volume<TO>& out, float start_y,
          float start_x, float start_z, uint64_t line_len, int num_threads) {
    pk5d::check_threads(num_threads);
    double tmat[4][4];
    read_tmat(m, tmat);
    uint64_t sz[3], osz[3];
    shape3(vol, sz);
    shape3(out, osz);
    check_line_len(line_len, start_y, sz[1], osz[1]);
    nb::gil_scoped_release release;
    pk5d::deskew::warp_volume<T, TO>(vol.data(), sz[0], sz[1], sz[2], tmat, out.data(), osz[0],
                                     osz[1], osz[2], start_y, start_x, start_z, line_len,
                                     num_threads);
}

template <typename T, typename TO>
void skewed_warp(const in_volume<T>& vol, double xstep, double stepsize, bool reverse,
                 const matrix4& m, out_volume<TO>& out, float start_y, float start_x, float start_z,
                 uint64_t line_len, int num_threads) {
    pk5d::check_threads(num_threads);
    double tmat[4][4];
    read_tmat(m, tmat);
    uint64_t sz[3], osz[3];
    shape3(vol, sz);
    shape3(out, osz);
    if (!(xstep > 0 && xstep < 9.0e15) || !(stepsize > 0 && stepsize < 9.0e15))
        throw std::invalid_argument("xstep and the step size must be positive and finite");
    if (sz[2] < (uint64_t)ceil(xstep) + 2)
        throw std::invalid_argument("nx must be at least ceil(xstep) + 2");
    if (!(pk5d::deskew::interp_grid_size(sz[0], stepsize) < 9.0e15))
        throw std::invalid_argument("the step size is too small");
    check_line_len(line_len, start_y, sz[1], osz[1]);
    nb::gil_scoped_release release;
    pk5d::deskew::skewed_space_interp_warp<T, TO>(
        vol.data(), sz[0], sz[1], sz[2], xstep, stepsize, reverse ? 1 : 0, tmat, out.data(),
        osz[0], osz[1], osz[2], start_y, start_x, start_z, line_len, num_threads);
}

template <typename T, typename TO>
void interp_defined_stepsize(const in_volume<T>& vol, float xstep, double stepsize, bool reverse,
                             out_volume<TO>& out, int num_threads) {
    pk5d::check_threads(num_threads);
    uint64_t sz[3], osz[3];
    shape3(vol, sz);
    shape3(out, osz);
    if (!(xstep > 0 && xstep < 9.0e15) || !(stepsize > 0 && stepsize < 9.0e15))
        throw std::invalid_argument("xstep and the step size must be positive and finite");
    if (sz[2] < (uint64_t)ceil((double)xstep) + 2)
        throw std::invalid_argument("nx must be at least ceil(xstep) + 2");
    const double nint = pk5d::deskew::interp_grid_size(sz[0], stepsize);
    if (!(nint < 9.0e15) || osz[1] != sz[1] || osz[2] != sz[2] || (double)osz[0] != nint)
        throw std::invalid_argument("the output must be (the grid size) x ny x nx");
    nb::gil_scoped_release release;
    pk5d::deskew::skewed_space_interp_defined_stepsize<T, TO>(
        vol.data(), out.data(), xstep, stepsize, reverse ? 1 : 0, sz[0], sz[1], sz[2], osz[0],
        num_threads);
}

// The kernels on MATLAB's layout (ml/: the mex's loops): volumes in C order with axes (z, x, y),
// y fastest, as MATLAB's (y, x, z) arrays are in memory. The mex's sizes are (y, x, z).
template <typename A>
void mshape3(const A& a, uint64_t (&sz)[3]) {
    sz[0] = (uint64_t)a.shape(2);
    sz[1] = (uint64_t)a.shape(1);
    sz[2] = (uint64_t)a.shape(0);
}

template <typename T, typename TO>
void warp_ml(const in_volume<T>& vol, const matrix4& m, out_volume<TO>& out, float start_y,
             float start_x, float start_z, uint64_t line_len, int num_threads) {
    pk5d::check_threads(num_threads);
    double tmat[4][4];
    read_tmat(m, tmat);
    uint64_t sz[3], osz[3];
    mshape3(vol, sz);
    mshape3(out, osz);
    check_line_len(line_len, start_y, sz[0], osz[0]);
    nb::gil_scoped_release release;
    pk5d::deskew::ml::transformVolume<T, TO>(vol.data(), tmat, out.data(), sz[0], sz[1], sz[2],
                                             start_y, start_x, start_z, osz[0], osz[1], osz[2],
                                             line_len, num_threads);
}

template <typename T, typename TI, typename TO>
void skewed_warp_ml(const in_volume<T>& vol, double xstep, double stepsize, bool reverse,
                    const matrix4& m, out_volume<TO>& out, float start_y, float start_x,
                    float start_z, uint64_t line_len, int num_threads) {
    pk5d::check_threads(num_threads);
    double tmat[4][4];
    read_tmat(m, tmat);
    uint64_t sz[3], osz[3];
    mshape3(vol, sz);
    mshape3(out, osz);
    if (!(xstep > 0 && xstep < 9.0e15) || !(stepsize > 0 && stepsize < 9.0e15))
        throw std::invalid_argument("xstep and the step size must be positive and finite");
    if (sz[1] < (uint64_t)ceil(xstep) + 2)
        throw std::invalid_argument("nx must be at least ceil(xstep) + 2");
    if (!(pk5d::deskew::ml::interp_grid_size(sz[2], stepsize) < 9.0e15))
        throw std::invalid_argument("the step size is too small");
    check_line_len(line_len, start_y, sz[0], osz[0]);
    nb::gil_scoped_release release;
    pk5d::deskew::ml::skewed_space_interp_warp<T, TI, TO>(
        vol.data(), sz, xstep, stepsize, reverse ? 1 : 0, tmat, out.data(), start_y, start_x,
        start_z, osz[0], osz[1], osz[2], line_len, num_threads);
}

template <typename T, typename TO>
void interp_ml(const in_volume<T>& vol, float xstep, double stepsize, bool reverse,
               out_volume<TO>& out, int num_threads) {
    pk5d::check_threads(num_threads);
    uint64_t sz[3], osz[3];
    mshape3(vol, sz);
    mshape3(out, osz);
    if (!(xstep > 0 && xstep < 9.0e15) || !(stepsize > 0 && stepsize < 9.0e15))
        throw std::invalid_argument("xstep and the step size must be positive and finite");
    if (sz[1] < (uint64_t)ceil((double)xstep) + 2)
        throw std::invalid_argument("nx must be at least ceil(xstep) + 2");
    const double nint = pk5d::deskew::ml::interp_grid_size(sz[2], stepsize);
    if (!(nint < 9.0e15) || osz[0] != sz[0] || osz[1] != sz[1] || (double)osz[2] != nint)
        throw std::invalid_argument("the output must be (the grid size) x nx x ny");
    nb::gil_scoped_release release;
    pk5d::deskew::ml::skewed_space_interp_defined_stepsize<T, TO>(
        vol.data(), out.data(), xstep, stepsize, reverse ? 1 : 0, sz, osz, num_threads);
}

// dst (n, c, r) = src (n, r, c) with the last two axes swapped
template <typename T>
void swap_axes(const in_volume<T>& src, out_volume<T>& dst, int num_threads) {
    pk5d::check_threads(num_threads);
    const uint64_t n = src.shape(0), r = src.shape(1), c = src.shape(2);
    if (dst.shape(0) != n || dst.shape(1) != c || dst.shape(2) != r)
        throw std::invalid_argument("dst must be the shape of src with the last two axes swapped");
    nb::gil_scoped_release release;
    pk5d::deskew::swap_last_axes(src.data(), dst.data(), n, r, c, num_threads);
}

const char* WARP_DOC =
    "warp(vol, tmat, out, start_y, start_x, start_z, line_len, num_threads)\n\n"
    "volume_deskew_rotate_warp_mex's kernel: fills out (zeros, (z, y, x) C order) from vol "
    "(single, or uint16 to single or uint16 out). tmat: MATLAB's 4 x 4 matrix (C order, "
    "tmat[i, j] = M(i+1, j+1)); start_*: the box's starts in single precision (the mex's "
    "startX, startY, startZ are y, x, z); line_len: rows y computed per output line.";

const char* SKEWED_DOC =
    "skewed_warp(vol, xstep, stepsize, reverse, tmat, out, start_y, start_x, start_z, "
    "line_len, num_threads)\n\n"
    "skewed_space_interp_volume_deskew_rotate_warp_mex's fused kernel (xstep in double).";

const char* INTERP_DOC =
    "interp_defined_stepsize(vol, xstep, stepsize, reverse, out, num_threads)\n\n"
    "skewed_space_interp_defined_stepsize_mex's kernel (xstep in single); out is "
    "interp_grid_size(nz, stepsize) x ny x nx.";

}  // namespace

NB_MODULE(PK5D_MODULE_NAME, m) {
    m.doc() = "Deskew/rotate kernels (the deskew/rotate mex files of PetaKit5D), (z, y, x) C order";

    m.def("warp", &warp<float, float>, "vol"_a.noconvert(), "tmat"_a.noconvert(),
          "out"_a.noconvert(), "start_y"_a, "start_x"_a, "start_z"_a, "line_len"_a,
          "num_threads"_a, WARP_DOC);
    m.def("warp", &warp<uint16_t, float>, "vol"_a.noconvert(), "tmat"_a.noconvert(),
          "out"_a.noconvert(), "start_y"_a, "start_x"_a, "start_z"_a, "line_len"_a,
          "num_threads"_a);
    m.def("warp", &warp<uint16_t, uint16_t>, "vol"_a.noconvert(), "tmat"_a.noconvert(),
          "out"_a.noconvert(), "start_y"_a, "start_x"_a, "start_z"_a, "line_len"_a,
          "num_threads"_a);

    m.def("skewed_warp", &skewed_warp<float, float>, "vol"_a.noconvert(), "xstep"_a,
          "stepsize"_a, "reverse"_a, "tmat"_a.noconvert(), "out"_a.noconvert(), "start_y"_a,
          "start_x"_a, "start_z"_a, "line_len"_a, "num_threads"_a, SKEWED_DOC);
    m.def("skewed_warp", &skewed_warp<uint16_t, float>, "vol"_a.noconvert(), "xstep"_a,
          "stepsize"_a, "reverse"_a, "tmat"_a.noconvert(), "out"_a.noconvert(), "start_y"_a,
          "start_x"_a, "start_z"_a, "line_len"_a, "num_threads"_a);
    m.def("skewed_warp", &skewed_warp<uint16_t, uint16_t>, "vol"_a.noconvert(), "xstep"_a,
          "stepsize"_a, "reverse"_a, "tmat"_a.noconvert(), "out"_a.noconvert(), "start_y"_a,
          "start_x"_a, "start_z"_a, "line_len"_a, "num_threads"_a);

    m.def("interp_defined_stepsize", &interp_defined_stepsize<float, float>, "vol"_a.noconvert(),
          "xstep"_a, "stepsize"_a, "reverse"_a, "out"_a.noconvert(), "num_threads"_a, INTERP_DOC);
    m.def("interp_defined_stepsize", &interp_defined_stepsize<uint16_t, float>,
          "vol"_a.noconvert(), "xstep"_a, "stepsize"_a, "reverse"_a, "out"_a.noconvert(),
          "num_threads"_a);
    m.def("interp_defined_stepsize", &interp_defined_stepsize<uint16_t, uint16_t>,
          "vol"_a.noconvert(), "xstep"_a, "stepsize"_a, "reverse"_a, "out"_a.noconvert(),
          "num_threads"_a);

    // MATLAB's layout: (z, x, y) C-ordered arrays
    m.def("warp_ml", &warp_ml<float, float>, "vol"_a.noconvert(), "tmat"_a.noconvert(),
          "out"_a.noconvert(), "start_y"_a, "start_x"_a, "start_z"_a, "line_len"_a,
          "num_threads"_a, "warp on MATLAB's layout: vol and out (z, x, y), y fastest");
    m.def("warp_ml", &warp_ml<uint16_t, float>, "vol"_a.noconvert(), "tmat"_a.noconvert(),
          "out"_a.noconvert(), "start_y"_a, "start_x"_a, "start_z"_a, "line_len"_a,
          "num_threads"_a);
    m.def("warp_ml", &warp_ml<uint16_t, uint16_t>, "vol"_a.noconvert(), "tmat"_a.noconvert(),
          "out"_a.noconvert(), "start_y"_a, "start_x"_a, "start_z"_a, "line_len"_a,
          "num_threads"_a);
    m.def("skewed_warp_ml", &skewed_warp_ml<float, float, float>, "vol"_a.noconvert(), "xstep"_a,
          "stepsize"_a, "reverse"_a, "tmat"_a.noconvert(), "out"_a.noconvert(), "start_y"_a,
          "start_x"_a, "start_z"_a, "line_len"_a, "num_threads"_a,
          "skewed_warp on MATLAB's layout: vol and out (z, x, y), y fastest");
    m.def("skewed_warp_ml", &skewed_warp_ml<uint16_t, float, float>, "vol"_a.noconvert(),
          "xstep"_a, "stepsize"_a, "reverse"_a, "tmat"_a.noconvert(), "out"_a.noconvert(),
          "start_y"_a, "start_x"_a, "start_z"_a, "line_len"_a, "num_threads"_a);
    m.def("skewed_warp_ml", &skewed_warp_ml<uint16_t, uint16_t, uint16_t>, "vol"_a.noconvert(),
          "xstep"_a, "stepsize"_a, "reverse"_a, "tmat"_a.noconvert(), "out"_a.noconvert(),
          "start_y"_a, "start_x"_a, "start_z"_a, "line_len"_a, "num_threads"_a);
    m.def("interp_defined_stepsize_ml", &interp_ml<float, float>, "vol"_a.noconvert(), "xstep"_a, "stepsize"_a,
          "reverse"_a, "out"_a.noconvert(), "num_threads"_a,
          "interp_defined_stepsize on MATLAB's layout: vol and out (z, x, y), y fastest");
    m.def("interp_defined_stepsize_ml", &interp_ml<uint16_t, float>, "vol"_a.noconvert(), "xstep"_a, "stepsize"_a,
          "reverse"_a, "out"_a.noconvert(), "num_threads"_a);
    m.def("interp_defined_stepsize_ml", &interp_ml<uint16_t, uint16_t>, "vol"_a.noconvert(), "xstep"_a,
          "stepsize"_a, "reverse"_a, "out"_a.noconvert(), "num_threads"_a);
    m.def("swap_axes", &swap_axes<float>, "src"_a.noconvert(), "dst"_a.noconvert(),
          "num_threads"_a, "dst (n, c, r) = src (n, r, c): the last two axes swapped, in threads");
    m.def("swap_axes", &swap_axes<uint16_t>, "src"_a.noconvert(), "dst"_a.noconvert(),
          "num_threads"_a);

    m.def("interp_grid_size", &pk5d::deskew::interp_grid_size, "nz"_a, "stepsize"_a,
          "floor(round((nz - 1) / stepsize * 100000) / 100000) + 1 (0 for nz = 0), in double: "
          "the slices of the interpolation grid");

    m.def("simd", []() {
#if defined(__AVX2__)
        return std::string("avx2");
#else
        return std::string("scalar");
#endif
    }, "the lanes compiled: 'avx2' or 'scalar' (the same operations in the same order)");
}
