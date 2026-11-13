// The mask erosion of decon_mask_edge_erosion_gpu_mex
// (microscopeDataProcessing/decon/mex/decon_mask_edge_erosion_gpu_mex.cu) for
// petakit5d.kernels.decon_gpu, on the port's (z, y, x) C arrays: each z slice padded with true by
// padUL and padLR, eroded by each 2D neighborhood in turn (pixels outside the padded slice count as
// true), cropped back, then eroded along z by a line of 2 * rz + 1 voxels (voxels outside the
// volume count as true). Logical AND and counts only, so the result doesn't depend on the order in
// which voxels are visited; the loops are the mex's, with columns (x) on threadIdx.x for the C
// layout. No atomics.
//
// Masks are bytes 0 or 1. offsets: (row, column) pairs of the nonzero elements of all
// neighborhoods (the mex kept them in constant memory).

typedef long long idx_t;

// rows (y) and columns (x) of a slice in grid-stride loops, slices along z
#define SLICES for (idx_t z = blockIdx.z; z < d; z += gridDim.z)
#define ROWS_OF(r, n) for (idx_t r = blockIdx.y * (idx_t)blockDim.y + threadIdx.y; r < (n); r += (idx_t)gridDim.y * blockDim.y)
#define COLS_OF(c, n) for (idx_t c = blockIdx.x * (idx_t)blockDim.x + threadIdx.x; c < (n); c += (idx_t)gridDim.x * blockDim.x)

// the padded volume p (d, hp, wp): true outside the slice, the (normalized) mask a (d, h, w) inside,
// with the border cleared (load_kernel)
extern "C" __global__ void load_kernel(unsigned char* __restrict__ p, const unsigned char* __restrict__ a,
                                       const idx_t h, const idx_t w, const idx_t d, const idx_t pu0,
                                       const idx_t pu1, const idx_t hp, const idx_t wp, const int clear_border) {
    SLICES {
        const bool clear_slice = clear_border && d > 1 && (z == 0 || z == d - 1);
        ROWS_OF(r, hp) {
            const idx_t y = r - pu0;
            COLS_OF(c, wp) {
                const idx_t x = c - pu1;
                unsigned char v = 1;
                if (y >= 0 && y < h && x >= 0 && x < w) {
                    const bool clear = clear_slice || (clear_border && (y == 0 || y == h - 1 || x == 0 || x == w - 1));
                    v = clear ? 0 : (a[x + w * (y + h * z)] != 0);
                }
                p[c + wp * (r + hp * z)] = v;
            }
        }
    }
}

// b = a eroded by the neighborhood with the offsets off[i0 .. i0 + n - 1], in the padded volume
// (erode_kernel)
extern "C" __global__ void erode_kernel(unsigned char* __restrict__ b, const unsigned char* __restrict__ a,
                                        const int2* __restrict__ off, const int i0, const int n, const idx_t d,
                                        const idx_t hp, const idx_t wp) {
    SLICES {
        const unsigned char* as = a + hp * wp * z;
        unsigned char* bs = b + hp * wp * z;
        ROWS_OF(r, hp) {
            COLS_OF(c, wp) {
                unsigned char v = 1;
                for (int k = 0; k < n; k++) {
                    const int2 o = off[i0 + k];
                    const idx_t rr = r + o.x, cc = c + o.y;
                    if (rr >= 0 && rr < hp && cc >= 0 && cc < wp) {
                        v &= as[cc + wp * rr];
                    }
                }
                bs[c + wp * r] = v;
            }
        }
    }
}

// out (d, h, w) = the cropped slices of p eroded along z by 2 * rz + 1 voxels: one thread per voxel
// column, counting the false voxels in its window (z_kernel)
extern "C" __global__ void z_kernel(unsigned char* __restrict__ out, const unsigned char* __restrict__ p,
                                    const idx_t rz, const idx_t h, const idx_t w, const idx_t d, const idx_t pu0,
                                    const idx_t pu1, const idx_t hp, const idx_t wp) {
    const idx_t n2 = h * w, np2 = hp * wp;
    ROWS_OF(y, h) {
        COLS_OF(x, w) {
            const unsigned char* col = p + (x + pu1) + wp * (y + pu0);
            unsigned char* o = out + x + w * y;
            idx_t cnt = 0;
            const idx_t zi = rz < d - 1 ? rz : d - 1;
            for (idx_t z = 0; z <= zi; z++) cnt += 1 - col[np2 * z];
            for (idx_t z = 0; z < d; z++) {
                o[n2 * z] = cnt == 0;
                if (z + rz + 1 < d) cnt += 1 - col[np2 * (z + rz + 1)];
                if (z - rz >= 0) cnt -= 1 - col[np2 * (z - rz)];
            }
        }
    }
}
