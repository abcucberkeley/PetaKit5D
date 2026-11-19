// The deskew/rotate warp without interpolation in skewed space (volume_deskew_rotate_warp_mex)
// on C-ordered (z, y, x) volumes: single to single, uint16 to single, uint16 to uint16, through
// the blocked kernel (blocked.hpp) with a grid of raw slices only.
#pragma once

#include <cstdint>

#include "blocked.hpp"

namespace pk5d::deskew {

// the grid of the warp alone: grid slice q is raw slice q
inline void identity_grid(const uint64_t nz, Grid& G) {
    G.ntotal = nz;
    G.nint = 0;
    G.xa = 0;
    G.interp_index.resize(nz);
    G.is_original.reset(new bool[nz > 0 ? nz : 1]);
    for (uint64_t q = 0; q < nz; ++q) {
        G.interp_index[q] = q;
        G.is_original[q] = true;
    }
    G.slices.assign(1, InterpSlice{});
}

// block_rows and window 0: chosen by the kernel
template <typename T, typename TO>
void warp_volume(const T* const vol, const uint64_t nz, const uint64_t ny, const uint64_t nx,
                 const double (&tmat)[4][4], TO* const out, const uint64_t oz, const uint64_t oy,
                 const uint64_t ox, const float start_y, const float start_x, const float start_z,
                 const uint64_t line_len, const bool fill, const int nthreads,
                 const uint64_t block_rows, const uint64_t window) {
    Grid G;
    identity_grid(nz, G);
    blocked_warp<T, TO>(vol, nz, ny, nx, G, false, 0, tmat, out, oz, oy, ox, start_y, start_x,
                        start_z, line_len, fill, nthreads, block_rows, window);
}

}  // namespace pk5d::deskew
