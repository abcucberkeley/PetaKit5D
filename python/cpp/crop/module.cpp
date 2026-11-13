// petakit5d._ext.crop: the copies of crop3d_mex, crop4d_mex, indexing3d_mex, indexing4d_mex and
// indexing4d_crop_mex (microscopeDataProcessing/crop/mex) on C-contiguous arrays.
//
// All five mex functions copy a box of one array into another with memcpy, run by run along
// the fastest axis: a crop copies the box of the data into a new array, the indexing functions
// copy a region (or a box of it) into the box of the data, in place. One kernel does both: it
// copies the box of extent (t, z, y, x) at src_origin of src into dst at dst_origin. Arrays are
// seen as 4-D (t, z, y, x) C-contiguous byte buffers with an element size; argument checks and
// MATLAB's error identifiers are on the Python side (petakit5d.kernels.crop).
//
// As the mex's continuity branches, axes along which the box spans both arrays whole are merged
// into one run (a whole x row gives runs of planes, whole planes runs of volumes). The bytes of
// all runs are cut into one contiguous share per thread, so a box of a few long runs is split
// too (the mex gives each thread whole rows or volumes), and each memcpy takes at most `piece`
// bytes (the caller's choice: one huge call makes the C library use non-temporal stores, which
// are slower into a fresh array whose pages the copy faults in). memcpy is the C library's, which picks
// its instructions for the CPU at run time; a copy has no arithmetic, so the result is the same
// bits on every CPU and for every thread count.
#include <pk5d/common.hpp>

#include "ops_bind.hpp"

#include <nanobind/stl/vector.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace nb = nanobind;
using u64 = uint64_t;

namespace {

struct copy_job {
    const uint8_t *src;
    uint8_t *dst;
    u64 sstr[4], dstr[4];  // byte strides of (t, z, y, x)
    u64 ext[4];            // extent of the box
    int outer;             // axes 0 .. outer - 1 step over runs; the rest is one run
    u64 run;               // bytes of a run
    u64 nrun;              // runs
    u64 piece;             // bytes of one memcpy at most
};

// the bytes [lo, hi) of the runs, in run order
void copy_part(const copy_job &j, u64 lo, u64 hi) {
    if (lo >= hi) return;
    u64 r = lo / j.run, off = lo % j.run;
    // the run's index along the outer axes (mixed radix, last outer axis fastest)
    u64 idx[4] = {0, 0, 0, 0};
    for (int a = j.outer - 1; a >= 0; a--) {
        idx[a] = r % j.ext[a];
        r /= j.ext[a];
    }
    u64 pos = lo;
    while (pos < hi) {
        u64 so = 0, d = 0;
        for (int a = 0; a < j.outer; a++) {
            so += idx[a] * j.sstr[a];
            d += idx[a] * j.dstr[a];
        }
        const u64 len = std::min(std::min(j.run - off, hi - pos), j.piece);
        std::memcpy(j.dst + d + off, j.src + so + off, len);
        pos += len;
        off += len;
        if (off < j.run) continue;
        off = 0;
        for (int a = j.outer - 1; a >= 0; a--) {
            if (++idx[a] < j.ext[a]) break;
            idx[a] = 0;
        }
    }
}

std::vector<u64> four(const std::vector<u64> &v, const char *what) {
    if (v.size() != 4) throw std::invalid_argument(std::string(what) + " must have 4 values");
    return v;
}

void bind_copy_box(nb::ndarray<const uint8_t, nb::ndim<1>, nb::c_contig, nb::device::cpu> src,
                   std::vector<u64> src_shape, std::vector<u64> src_origin,
                   nb::ndarray<uint8_t, nb::ndim<1>, nb::c_contig, nb::device::cpu> dst,
                   std::vector<u64> dst_shape, std::vector<u64> dst_origin,
                   std::vector<u64> extent, u64 itemsize, int num_threads, bool bind,
                   u64 piece) {
    pk5d::check_threads(num_threads);
    const auto ss = four(src_shape, "src_shape"), so = four(src_origin, "src_origin");
    const auto ds = four(dst_shape, "dst_shape"), dor = four(dst_origin, "dst_origin");
    const auto e = four(extent, "extent");
    if (itemsize == 0) throw std::invalid_argument("itemsize must be positive");
    u64 sn = itemsize, dn = itemsize;
    for (int a = 0; a < 4; a++) {
        sn *= ss[a];
        dn *= ds[a];
        if (so[a] > ss[a] || e[a] > ss[a] - so[a])
            throw std::invalid_argument("the box is outside src");
        if (dor[a] > ds[a] || e[a] > ds[a] - dor[a])
            throw std::invalid_argument("the box is outside dst");
    }
    if (src.shape(0) != sn) throw std::invalid_argument("src has another size than src_shape");
    if (dst.shape(0) != dn) throw std::invalid_argument("dst has another size than dst_shape");
    if (e[0] == 0 || e[1] == 0 || e[2] == 0 || e[3] == 0) return;

    copy_job j{};
    j.src = src.data();
    j.dst = dst.data();
    j.sstr[3] = j.dstr[3] = itemsize;
    for (int a = 2; a >= 0; a--) {
        j.sstr[a] = j.sstr[a + 1] * ss[a + 1];
        j.dstr[a] = j.dstr[a + 1] * ds[a + 1];
    }
    u64 sbase = 0, dbase = 0;
    for (int a = 0; a < 4; a++) {
        j.ext[a] = e[a];
        sbase += so[a] * j.sstr[a];
        dbase += dor[a] * j.dstr[a];
    }
    j.src += sbase;
    j.dst += dbase;
    // merge the axes the box spans whole in both arrays into the run
    int m = 3;
    j.run = e[3] * itemsize;
    while (m > 0 && e[m] == ss[m] && e[m] == ds[m]) {
        m--;
        j.run *= e[m];
    }
    j.outer = m;
    j.piece = piece == 0 ? ~u64(0) : piece;
    j.nrun = 1;
    for (int a = 0; a < m; a++) j.nrun *= e[a];

    const u64 total = j.nrun * j.run;
    const int nt = (int)std::min<u64>((u64)num_threads, total);
    nb::gil_scoped_release release;
    pk5d::ops::parallel(nt, bind, [&](int k, int nk) {
        std::ptrdiff_t lo, hi;
        pk5d::ops::share((std::ptrdiff_t)total, k, nk, lo, hi);
        copy_part(j, (u64)lo, (u64)hi);
    });
}

}  // namespace

NB_MODULE(crop, m) {
    using namespace nb::literals;
    m.doc() = "Box copies of the crop and indexing mex functions on C-contiguous arrays";
    m.def("copy_box", &bind_copy_box, "src"_a.noconvert(), "src_shape"_a, "src_origin"_a,
          "dst"_a.noconvert(), "dst_shape"_a, "dst_origin"_a, "extent"_a, "itemsize"_a,
          "num_threads"_a, "bind"_a = true, "piece"_a = 0,
          "Copy the box of extent (t, z, y, x) at src_origin of src into dst at dst_origin. src "
          "and dst are the bytes (1-D uint8, C order) of C-contiguous 4-D arrays of the given "
          "shapes and element size; they must not overlap.");
}
