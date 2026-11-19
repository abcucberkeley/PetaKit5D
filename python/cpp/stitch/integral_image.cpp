// Integral image of a 3D image over a window of size T, in single precision:
// integral_image_3d_mex.cpp. The mex's three passes keep their axis, direction, order and
// arithmetic: running window sums along y (MATLAB's first dimension) from the start, then the
// differences along x and then along z from the end. Each element goes through the mex's
// operations in the mex's order; the passes visit the other axes in the C layout's order (and
// the y and z passes handle several x at once).
//
// a: (nz, ny, nx), out: (nz + tz - 1, ny + ty - 1, nx + tx - 1), C order, zeroed by the caller
// (the mex's mxCreateNumericArray).
//
// NaN with single input: the mex gives the first NaN operand of each addition and subtraction
// as written, else x86's NaN for Inf - Inf (measured). An input with NaN is computed again with
// each operation's NaN made so (pk5d::x86_nans), since compilers swap the operands of a + b.
// Where the input holds no NaN, every NaN is an invalid operation's, and which operand's it is
// changes no bits: x86's NaN, except on ARM, where an output with NaN is computed again too.
// The sums of uint16 input are finite.
#include "stitch.hpp"

#include <algorithm>
#include <type_traits>

namespace {

using namespace pk5d_stitch;

constexpr uint64_t XBLOCK = 512;  // x values per work item of the y and z passes

// a + b and a - b; X: with x86's NaN for the operands as written
template <bool X>
inline float add(float a, float b) {
    const float r = a + b;
    if constexpr (X) return pk5d::x86_nans(r, a, b);
    else return r;
}

template <bool X>
inline float sub(float a, float b) {
    const float r = a - b;
    if constexpr (X) return pk5d::x86_nans(r, a, b);
    else return r;
}

#ifdef PK5D_FP_ARM
// whether any of the n values is NaN
bool any_nan(const float *v, uint64_t n, int nthread) {
    int nan = 0;
#pragma omp parallel for schedule(static) num_threads(nthread) if (nthread > 1) reduction(| : nan)
    for (int64_t i = 0; i < (int64_t)n; i++) nan |= v[i] != v[i];
    return nan != 0;
}
#endif

// the mex's integral_image_x: along MATLAB's first dimension (y), forward, for the input's
// z and x. Returns whether a single input holds NaN (each value is read once as f or fy; not
// asked with X).
template <typename T, bool X>
bool pass_y(const T *a, float *n, const uint64_t *dims, const uint64_t *dimsT,
            const uint64_t *dimsI, int nthread) {
    const uint64_t nZ = dims[0], nY = dims[1], nX = dims[2];
    const uint64_t iY = dimsI[1], iX = dimsI[2];
    const uint64_t tY = dimsT[1];
    const uint64_t nB = (nX + XBLOCK - 1) / XBLOCK;
    constexpr bool check = std::is_same_v<T, float> && !X;
    int has_nan = 0;

#pragma omp parallel for collapse(2) schedule(static) num_threads(nthread) if (nthread > 1) \
    reduction(| : has_nan)
    for (uint64_t z = 0; z < nZ; z++) {
        for (uint64_t b = 0; b < nB; b++) {
            const uint64_t x0 = b * XBLOCK, x1 = std::min(nX, x0 + XBLOCK);
            const T *f = a + z * nY * nX;
            float *out = n + z * iY * iX;
            for (uint64_t x = x0; x < x1; x++) {
                out[x] = (float)(f[x]);
                if constexpr (check) has_nan |= f[x] != f[x];
            }
            for (uint64_t y = 1; y < iY; y++) {
                float *o = out + y * iX;
                const float *op = out + (y - 1) * iX;
                if (y < tY) {
                    const T *fy = f + y * nX;
                    for (uint64_t x = x0; x < x1; x++) {
                        o[x] = add<X>(op[x], (float)(fy[x]));
                        if constexpr (check) has_nan |= fy[x] != fy[x];
                    }
                } else if (y < nY) {
                    const T *fy = f + y * nX, *ft = f + (y - tY) * nX;
                    for (uint64_t x = x0; x < x1; x++) {
                        o[x] = sub<X>(add<X>(op[x], (float)(fy[x])), (float)(ft[x]));
                        if constexpr (check) has_nan |= fy[x] != fy[x];
                    }
                } else {
                    const T *ft = f + (y - tY) * nX;
                    for (uint64_t x = x0; x < x1; x++) o[x] = sub<X>(op[x], (float)(ft[x]));
                }
            }
        }
    }
    return has_nan != 0;
}

// the mex's integral_image_y: along MATLAB's second dimension (x), backward, for the input's z
// and every y of the output; tmp_post is the value of x + 1 before it was overwritten
template <bool X>
void pass_x(float *n, const uint64_t *dims, const uint64_t *dimsT, const uint64_t *dimsI,
            int nthread) {
    const uint64_t nZ = dims[0], nX = dims[2];
    const uint64_t iY = dimsI[1], iX = dimsI[2];
    const uint64_t tX = dimsT[2];
    const int64_t nRows = (int64_t)(nZ * iY);

#pragma omp parallel for schedule(static) num_threads(nthread) if (nthread > 1)
    for (int64_t r = 0; r < nRows; r++) {
        float *o = n + (uint64_t)r * iX;
        float tmp_post = 0;
        for (uint64_t x = iX - 1; x > 0; x--) {
            float tmp = o[x];

            if (x == iX - 1) o[x] = o[iX - tX];
            else if (x >= nX) o[x] = add<X>(o[x + 1], o[x - tX + 1]);
            else if (x >= tX - 1) o[x] = add<X>(sub<X>(o[x + 1], tmp_post), o[x - tX + 1]);
            else o[x] = sub<X>(o[x + 1], tmp_post);

            tmp_post = tmp;
        }
    }
}

// the mex's integral_image_z: along z, backward, for every y and x of the output
template <bool X>
void pass_z(float *n, const uint64_t *dims, const uint64_t *dimsT, const uint64_t *dimsI,
            int nthread) {
    const uint64_t nZ = dims[0];
    const uint64_t iZ = dimsI[0], iY = dimsI[1], iX = dimsI[2];
    const uint64_t tZ = dimsT[0];
    const uint64_t plane = iY * iX;
    const uint64_t nB = (iX + XBLOCK - 1) / XBLOCK;

#pragma omp parallel num_threads(nthread) if (nthread > 1)
    {
        std::vector<float> post(XBLOCK), cur(XBLOCK);
#pragma omp for collapse(2) schedule(static)
        for (uint64_t y = 0; y < iY; y++) {
            for (uint64_t b = 0; b < nB; b++) {
                const uint64_t x0 = b * XBLOCK, w = std::min(iX, x0 + XBLOCK) - x0;
                float *col = n + y * iX + x0;
                float *tmp_post = post.data();
                float *tmp = cur.data();
                for (uint64_t x = 0; x < w; x++) tmp_post[x] = 0;
                for (uint64_t z = iZ - 1; z > 0; z--) {
                    float *o = col + z * plane;
                    for (uint64_t x = 0; x < w; x++) tmp[x] = o[x];
                    if (z == iZ - 1) {
                        const float *s = col + (iZ - tZ) * plane;
                        for (uint64_t x = 0; x < w; x++) o[x] = s[x];
                    } else if (z >= nZ) {
                        const float *s1 = o + plane, *s2 = col + (z - tZ + 1) * plane;
                        for (uint64_t x = 0; x < w; x++) o[x] = add<X>(s1[x], s2[x]);
                    } else if (z >= tZ - 1) {
                        const float *s1 = o + plane, *s2 = col + (z - tZ + 1) * plane;
                        for (uint64_t x = 0; x < w; x++)
                            o[x] = add<X>(sub<X>(s1[x], tmp_post[x]), s2[x]);
                    } else {
                        const float *s1 = o + plane;
                        for (uint64_t x = 0; x < w; x++) o[x] = sub<X>(s1[x], tmp_post[x]);
                    }
                    std::swap(tmp, tmp_post);
                }
            }
        }
    }
}

template <typename T>
void bind_integral(nb::module_ &m) {
    m.def(
        "integral_image_3d",
        [](in1d<T> a, out1d<float> out, std::vector<uint64_t> dims, std::vector<uint64_t> dimsT,
           int num_threads) {
            pk5d::check_threads(num_threads);
            need_dims(dims, 3, "dims");
            need_dims(dimsT, 3, "dimsT");
            for (int i = 0; i < 3; i++)
                need(dimsT[i] >= 1 && dimsT[i] <= dims[i], "T must be in 1 .. the input's size");
            std::vector<uint64_t> dimsI(3);
            for (int i = 0; i < 3; i++) dimsI[i] = dims[i] + dimsT[i] - 1;
            need_size(a, prod(dims), "input");
            need_size(out, prod(dimsI), "output");
            nb::gil_scoped_release release;
            const uint64_t *d = dims.data(), *dt = dimsT.data(), *di = dimsI.data();
            const bool has_nan = pass_y<T, false>(a.data(), out.data(), d, dt, di, num_threads);
            if constexpr (std::is_same_v<T, float>) {
                if (has_nan) {
                    // again, with each NaN as written: pass_y writes what it wrote, the rest
                    // is still zero
                    pass_y<T, true>(a.data(), out.data(), d, dt, di, num_threads);
                    pass_x<true>(out.data(), d, dt, di, num_threads);
                    pass_z<true>(out.data(), d, dt, di, num_threads);
                    return;
                }
            }
            pass_x<false>(out.data(), d, dt, di, num_threads);
            pass_z<false>(out.data(), d, dt, di, num_threads);
#ifdef PK5D_FP_ARM
            // ARM: the NaN of Inf - Inf has the sign bit clear; again from a zeroed output
            if constexpr (std::is_same_v<T, float>) {
                if (any_nan(out.data(), prod(dimsI), num_threads)) {
                    std::fill(out.data(), out.data() + prod(dimsI), 0.0f);
                    pass_y<T, true>(a.data(), out.data(), d, dt, di, num_threads);
                    pass_x<true>(out.data(), d, dt, di, num_threads);
                    pass_z<true>(out.data(), d, dt, di, num_threads);
                }
            }
#endif
        },
        nb::arg("a").noconvert(), nb::arg("out").noconvert(), nb::arg("dims"),
        nb::arg("dims_t"), nb::arg("num_threads"),
        "Integral image of a ((z, y, x) sizes dims) over windows dims_t into the zeroed out "
        "(sizes dims + dims_t - 1).");
}

}  // namespace

PK5D_BINDINGS(bind_integral_image) {
    bind_integral<uint16_t>(m);
    bind_integral<float>(m);
}
