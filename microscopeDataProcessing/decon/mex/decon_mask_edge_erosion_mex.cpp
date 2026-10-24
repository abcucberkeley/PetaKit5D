// Edge erosion of the decon masks (decon_mask_edge_erosion) with the same result as MATLAB's imerode, bit for
// bit: each z slice of a logical mask is eroded by a sequence of 2D neighborhoods as imerode applies a
// decomposed structuring element (e.g., the periodic lines of strel('disk', r)): the slice is padded with true
// by padUL and padLR (morphop's padding for decompositions), eroded by each neighborhood in turn with pixels
// outside the padded slice counting as true, and cropped back. Then the mask is eroded along z by a line of
// 2 * rz + 1 voxels (voxels outside the volume count as true), as imerode with true(2 * rz + 1, 1) along z.
// The erosion by a neighborhood is B(y, x) = all of A(y + i - cy, x + j - cx) over its nonzero elements (i, j),
// with the center (cy, cx) = floor((size + 1) / 2), as imerode's.
//
// B = decon_mask_edge_erosion_mex(A, nhoods, padUL, padLR, rz, clearBorder, nThreads)
//   A           : logical 2D or 3D mask
//   nhoods      : cell array of 2D neighborhoods (logical or numeric, each with a nonzero element), in order
//   padUL       : [rows, columns] of true padding before each slice (0 for none)
//   padLR       : [rows, columns] of true padding after each slice
//   rz          : radius of the erosion along z (0 for none; ignored for 2D)
//   clearBorder : true to set the first and last rows and columns (3D: and slices) to false first
//   nThreads    : threads
//
// Build: mex -v CXXOPTIMFLAGS="-O3 -DNDEBUG" CXXFLAGS='$CXXFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' decon_mask_edge_erosion_mex.cpp
// (compileAllMexFiles builds it on all platforms)

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <vector>
#include <omp.h>
#include "mex.h"

typedef ptrdiff_t idx_t;

// the neighborhoods as offsets (nhood k: [start[k], start[k + 1])), padding and the guard band of the buffers
struct Plan {
    int ns;
    std::vector<int> start, dr, dc;
    idx_t h, w, d;         // mask size
    idx_t pu[2], pl[2];    // padding before / after (rows, columns)
    idx_t gr, gc;          // guard rows / columns around the padded slice (true, never written)
    idx_t hp, wp, hg, wg;  // padded slice and buffer size
    bool clearBorder;
};

// buffer element of row r, column c of the padded slice
static inline idx_t bidx(const Plan& p, const idx_t r, const idx_t c) {
    return (c + p.gc) * p.hg + r + p.gr;
}

// columns c0 .. c1 - 1 of the padded slice from slice z of A (clear: the whole slice false)
static void load_cols(const Plan& p, const uint8_t* src, const bool clear, uint8_t* a, const idx_t c0, const idx_t c1) {
    for (idx_t c = c0; c < c1; c++) {
        uint8_t* col = a + bidx(p, 0, c);
        const idx_t x = c - p.pu[1];
        if (x < 0 || x >= p.w) {
            memset(col, 1, p.hp);
            continue;
        }
        memset(col, 1, p.pu[0]);
        memset(col + p.pu[0] + p.h, 1, p.pl[0]);
        uint8_t* cc = col + p.pu[0];
        if (clear || (p.clearBorder && (x == 0 || x == p.w - 1))) {
            memset(cc, 0, p.h);
            continue;
        }
        const uint8_t* s = src + x * p.h;
        for (idx_t y = 0; y < p.h; y++) {
            cc[y] = s[y] != 0;
        }
        if (p.clearBorder) {
            cc[0] = 0;
            cc[p.h - 1] = 0;
        }
    }
}

// erosion by neighborhood k of columns c0 .. c1 - 1 of the padded slice, from a to b, 16 rows at a time (GCC
// vector extension: SSE2 / NEON) with the source columns of the offsets in srcs (scratch)
typedef uint8_t vbytes __attribute__((vector_size(16)));

static void erode_cols(const Plan& p, const int k, const uint8_t* a, uint8_t* b, const idx_t c0, const idx_t c1,
        const uint8_t** srcs) {
    const int i0 = p.start[k], n = p.start[k + 1] - i0;
    const idx_t hp = p.hp, nv = hp / 16 * 16;
    for (idx_t c = c0; c < c1; c++) {
        uint8_t* o = b + bidx(p, 0, c);
        for (int i = 0; i < n; i++) {
            srcs[i] = a + bidx(p, p.dr[i0 + i], c + p.dc[i0 + i]);
        }
        idx_t y = 0;
        for (; y < nv; y += 16) {
            vbytes v, t;
            memcpy(&v, srcs[0] + y, 16);
            for (int i = 1; i < n; i++) {
                memcpy(&t, srcs[i] + y, 16);
                v &= t;
            }
            memcpy(o + y, &v, 16);
        }
        for (; y < hp; y++) {
            uint8_t v = srcs[0][y];
            for (int i = 1; i < n; i++) {
                v &= srcs[i][y];
            }
            o[y] = v;
        }
    }
}

// columns x0 .. x1 - 1 of the cropped slice from the padded slice a to dst
static void store_cols(const Plan& p, const uint8_t* a, uint8_t* dst, const idx_t x0, const idx_t x1) {
    for (idx_t x = x0; x < x1; x++) {
        memcpy(dst + x * p.h, a + bidx(p, p.pu[0], x + p.pu[1]), p.h);
    }
}

// the 2D erosions of slice z from A to B with the buffers a and b (guards true); returns nothing, the result is
// in B. Serial over the columns.
static void erode_slice(const Plan& p, const uint8_t* A, uint8_t* B, const idx_t z, uint8_t* a, uint8_t* b) {
    const idx_t n2 = p.h * p.w;
    const bool clear = p.clearBorder && p.d > 1 && (z == 0 || z == p.d - 1);
    std::vector<const uint8_t*> srcs(p.dr.size());
    load_cols(p, A + z * n2, clear, a, 0, p.wp);
    for (int k = 0; k < p.ns; k++) {
        erode_cols(p, k, a, b, 0, p.wp, srcs.data());
        uint8_t* t = a;
        a = b;
        b = t;
    }
    store_cols(p, a, B + z * n2, 0, p.w);
}

// the same with the columns split among the threads (for few slices)
static void erode_slice_parallel(const Plan& p, const uint8_t* A, uint8_t* B, const idx_t z, uint8_t* a, uint8_t* b,
        const int nt) {
    const idx_t n2 = p.h * p.w;
    const bool clear = p.clearBorder && p.d > 1 && (z == 0 || z == p.d - 1);
    #pragma omp parallel num_threads(nt)
    {
        uint8_t* ta = a;
        uint8_t* tb = b;
        std::vector<const uint8_t*> srcs(p.dr.size());
        #pragma omp for schedule(static)
        for (idx_t c = 0; c < p.wp; c++) {
            load_cols(p, A + z * n2, clear, ta, c, c + 1);
        }
        for (int k = 0; k < p.ns; k++) {
            #pragma omp for schedule(static)
            for (idx_t c = 0; c < p.wp; c++) {
                erode_cols(p, k, ta, tb, c, c + 1, srcs.data());
            }
            uint8_t* t = ta;
            ta = tb;
            tb = t;
        }
        #pragma omp for schedule(static)
        for (idx_t x = 0; x < p.w; x++) {
            store_cols(p, ta, B + z * n2, x, x + 1);
        }
    }
}

// erosion along z by 2 * rz + 1 voxels, in place, with per-thread copies of blocks of voxel columns; T counts
// the false voxels in the window of each voxel column
template <typename T>
static void erode_z(const Plan& p, uint8_t* B, const idx_t rz, const int nt) {
    const idx_t n2 = p.h * p.w, d = p.d;
    const idx_t bs = 1024;
    const idx_t nb = (n2 + bs - 1) / bs;
    std::vector<uint8_t*> blk(nt);
    std::vector<T*> cnt(nt);
    for (int t = 0; t < nt; t++) {
        blk[t] = (uint8_t*)mxMalloc((size_t)(bs * d));
        cnt[t] = (T*)mxMalloc(sizeof(T) * (size_t)bs);
    }
    #pragma omp parallel for num_threads(nt) schedule(dynamic, 4)
    for (idx_t ib = 0; ib < nb; ib++) {
        const int t = omp_get_thread_num();
        const idx_t p0 = ib * bs;
        const idx_t L = n2 - p0 < bs ? n2 - p0 : bs;
        uint8_t* __restrict v = blk[t];
        T* __restrict c = cnt[t];
        for (idx_t z = 0; z < d; z++) {
            memcpy(v + z * L, B + z * n2 + p0, L);
        }
        for (idx_t i = 0; i < L; i++) c[i] = 0;
        const idx_t zi = rz < d - 1 ? rz : d - 1;
        for (idx_t z = 0; z <= zi; z++) {
            const uint8_t* __restrict s = v + z * L;
            for (idx_t i = 0; i < L; i++) c[i] += (T)(1 - s[i]);
        }
        for (idx_t z = 0; z < d; z++) {
            uint8_t* __restrict o = B + z * n2 + p0;
            for (idx_t i = 0; i < L; i++) o[i] = c[i] == 0;
            if (z + rz + 1 < d) {
                const uint8_t* __restrict s = v + (z + rz + 1) * L;
                for (idx_t i = 0; i < L; i++) c[i] += (T)(1 - s[i]);
            }
            if (z - rz >= 0) {
                const uint8_t* __restrict s = v + (z - rz) * L;
                for (idx_t i = 0; i < L; i++) c[i] -= (T)(1 - s[i]);
            }
        }
    }
}

static double scalar_arg(const mxArray* a, const char* name) {
    if (!(mxIsNumeric(a) || mxIsLogical(a)) || mxIsComplex(a) || mxGetNumberOfElements(a) != 1) {
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_mex:inputError", "%s must be a real numeric scalar.", name);
    }
    return mxGetScalar(a);
}

static void pad_arg(const mxArray* a, const char* name, idx_t* v) {
    if (!mxIsDouble(a) || mxIsComplex(a) || mxGetNumberOfElements(a) != 2) {
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_mex:inputError", "%s must be a double vector of 2 elements.", name);
    }
    const double* pv = mxGetPr(a);
    for (int i = 0; i < 2; i++) {
        if (!(pv[i] >= 0) || pv[i] > 1e9 || pv[i] != (double)(idx_t)pv[i]) {
            mexErrMsgIdAndTxt("decon_mask_edge_erosion_mex:inputError", "%s must contain nonnegative integers.", name);
        }
        v[i] = (idx_t)pv[i];
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
#ifdef __APPLE__
    // keep the mex file loaded until MATLAB exits: GCC's OpenMP runtime leaves
    // thread exit handlers in it, so unloading it crashes MATLAB on exit
    if (!mexIsLocked()) {
        mexLock();
    }
#endif
    if (nrhs != 7) {
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_mex:inputError",
            "Usage: B = decon_mask_edge_erosion_mex(A, nhoods, padUL, padLR, rz, clearBorder, nThreads)");
    }
    const mxArray* ma = prhs[0];
    if (!mxIsLogical(ma) || mxIsSparse(ma) || mxGetNumberOfDimensions(ma) > 3) {
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_mex:inputError", "The mask must be a full 2D or 3D logical array.");
    }
    const mxArray* nh = prhs[1];
    if (!mxIsCell(nh)) {
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_mex:inputError", "nhoods must be a cell array.");
    }
    Plan p;
    const mwSize nd = mxGetNumberOfDimensions(ma);
    const mwSize* dims = mxGetDimensions(ma);
    p.h = (idx_t)dims[0];
    p.w = (idx_t)dims[1];
    p.d = nd > 2 ? (idx_t)dims[2] : 1;
    pad_arg(prhs[2], "padUL", p.pu);
    pad_arg(prhs[3], "padLR", p.pl);
    const double rzv = scalar_arg(prhs[4], "rz");
    if (!(rzv >= 0) || rzv > 1e9 || rzv != (double)(idx_t)rzv) {
        mexErrMsgIdAndTxt("decon_mask_edge_erosion_mex:inputError", "rz must be a nonnegative integer.");
    }
    const idx_t rz = (idx_t)rzv;
    p.clearBorder = scalar_arg(prhs[5], "clearBorder") != 0;
    const double ntv = scalar_arg(prhs[6], "nThreads");
    const int nt = ntv >= 1 && ntv <= 4096 ? (int)ntv : 1;

    // offsets of the nonzero elements of the neighborhoods
    p.ns = (int)mxGetNumberOfElements(nh);
    p.start.push_back(0);
    p.gr = 0;
    p.gc = 0;
    for (int k = 0; k < p.ns; k++) {
        const mxArray* e = mxGetCell(nh, k);
        if (e == NULL || !(mxIsLogical(e) || mxIsNumeric(e)) || mxIsComplex(e) || mxIsSparse(e)
                || mxGetNumberOfDimensions(e) > 2 || mxIsEmpty(e)) {
            mexErrMsgIdAndTxt("decon_mask_edge_erosion_mex:inputError", "Each neighborhood must be a real 2D array.");
        }
        const mxArray* lg = e;
        mxArray* tmp = NULL;
        if (!mxIsLogical(e)) {
            mxArray* in = (mxArray*)e;
            if (mexCallMATLAB(1, &tmp, 1, &in, "logical") != 0) {
                mexErrMsgIdAndTxt("decon_mask_edge_erosion_mex:inputError", "A neighborhood can't be converted to logical.");
            }
            lg = tmp;
        }
        const idx_t m = (idx_t)mxGetM(lg), n = (idx_t)mxGetN(lg);
        const mxLogical* v = mxGetLogicals(lg);
        const idx_t cy = (m + 1) / 2 - 1, cx = (n + 1) / 2 - 1;
        int cntk = 0;
        for (idx_t j = 0; j < n; j++) {
            for (idx_t i = 0; i < m; i++) {
                if (v[i + j * m]) {
                    const int r = (int)(i - cy), c = (int)(j - cx);
                    p.dr.push_back(r);
                    p.dc.push_back(c);
                    p.gr = p.gr > (r < 0 ? -r : r) ? p.gr : (r < 0 ? -r : r);
                    p.gc = p.gc > (c < 0 ? -c : c) ? p.gc : (c < 0 ? -c : c);
                    cntk++;
                }
            }
        }
        if (tmp != NULL) mxDestroyArray(tmp);
        if (cntk == 0) {
            mexErrMsgIdAndTxt("decon_mask_edge_erosion_mex:inputError", "Each neighborhood must have a nonzero element.");
        }
        p.start.push_back(p.start.back() + cntk);
    }

    plhs[0] = mxCreateLogicalArray(nd, dims);
    const idx_t n2 = p.h * p.w;
    if (n2 * p.d == 0) {
        return;
    }
    const uint8_t* A = (const uint8_t*)mxGetLogicals(ma);
    uint8_t* B = (uint8_t*)mxGetLogicals(plhs[0]);
    p.hp = p.h + p.pu[0] + p.pl[0];
    p.wp = p.w + p.pu[1] + p.pl[1];
    p.hg = p.hp + 2 * p.gr;
    p.wg = p.wp + 2 * p.gc;
    const size_t bsz = (size_t)(p.hg * p.wg);

    // slices: per thread with whole slices when there are enough of them, otherwise the columns split
    const bool perThread = nt == 1 || (p.d >= 2 * (idx_t)nt && 2 * bsz * (size_t)nt <= ((size_t)1 << 32));
    if (perThread) {
        std::vector<uint8_t*> bufs(2 * nt);
        for (int t = 0; t < 2 * nt; t++) {
            bufs[t] = (uint8_t*)mxMalloc(bsz);
            memset(bufs[t], 1, bsz);
        }
        #pragma omp parallel for num_threads(nt) schedule(dynamic, 1)
        for (idx_t z = 0; z < p.d; z++) {
            const int t = omp_get_thread_num();
            erode_slice(p, A, B, z, bufs[2 * t], bufs[2 * t + 1]);
        }
    } else {
        uint8_t* a = (uint8_t*)mxMalloc(bsz);
        uint8_t* b = (uint8_t*)mxMalloc(bsz);
        memset(a, 1, bsz);
        memset(b, 1, bsz);
        for (idx_t z = 0; z < p.d; z++) {
            erode_slice_parallel(p, A, B, z, a, b, nt);
        }
    }

    if (p.d > 1 && rz > 0) {
        const idx_t wmax = 2 * rz + 1 < p.d ? 2 * rz + 1 : p.d;
        if (wmax <= 255) {
            erode_z<uint8_t>(p, B, rz, nt);
        } else if (wmax <= 65535) {
            erode_z<uint16_t>(p, B, rz, nt);
        } else {
            erode_z<uint32_t>(p, B, rz, nt);
        }
    }
}
