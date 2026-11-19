// The grid of the fused kernel (skewed_space_interp_volume_deskew_rotate_warp_mex): which
// slices of the interpolation grid are raw slices, and the parameters of the interpolated ones
// (the raw slices zorg, zorg + 1, their shifts s, t in x and the weights), as the mex computes
// them; the kernels compute the interpolated values they need from them (blocked.hpp, ml/).
#pragma once

#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

namespace pk5d::deskew {

// an interpolated slice: between raw slices zorg and zorg + 1, rows shifted by s and t
struct InterpSlice {
    uint64_t zorg;
    uint64_t s;
    uint64_t t;
    float c00_w;
    float c01_w;
    float c10_w;
    float c11_w;
};

// number of slices of the interpolation grid (none for an empty input), in double
inline double interp_grid_size(const uint64_t nz, const double stepsize) {
    return nz == 0 ? 0 : floor(round((double)(nz - 1) / stepsize * 100000) / 100000) + 1;
}

// which slices of the grid are raw slices (within 1e-3) and which are interpolated
inline void determine_interpolation_indexing(const uint64_t n, const double stepSize,
                                             uint64_t* interp_index, bool* is_original,
                                             uint64_t& nint) {
    uint64_t interp_counter = 0;

    for (uint64_t i = 0; i < n; i++) {
        double cur_step = (double)i * stepSize;
        double cur_step_round = round(cur_step);
        if (fabs(cur_step - cur_step_round) < 1e-3) {
            is_original[i] = true;
            interp_index[i] = (uint64_t)cur_step_round;
        } else {
            is_original[i] = false;
            interp_index[i] = interp_counter;
            interp_counter++;
        }
    }
    nint = interp_counter;
}

// parameters of the interpolated slices (indexed by interp_index)
inline void skewed_space_interp_slices(const double xstep, const double stepsize,
                                       const uint8_t Reverse, const uint64_t* interp_index,
                                       const bool* is_original, const uint64_t ntotal,
                                       InterpSlice* slices) {
    for (uint64_t i = 0; i < ntotal; i++) {
        if (is_original[i]) continue;

        double s;
        double t;
        double zout;
        double zw;
        uint64_t zint = interp_index[i];

        zout = stepsize * (double)i;
        zw = zout - floor(zout);

        // slice zorg is read at x - xstep * zw (reverse) or x + xstep * zw (forward), slice
        // zorg + 1 at x + xstep * (1 - zw) or x - xstep * (1 - zw), as in
        // skewed_space_interp_defined_stepsize.m (s and t are ceil(xstep) - xstep < 1 apart
        // before flooring either way, which the column ranges of interp_column_lanes rely on)
        if (Reverse) {
            s = ceil(xstep) - xstep * zw;
            t = xstep * (1 - zw);
        } else {
            s = xstep * zw;
            t = ceil(xstep) - xstep * (1 - zw);
        }

        // distance to start
        double sw = s - floor(s);
        s = floor(s);
        double tw = t - floor(t);
        t = floor(t);

        // single precision distances, weights computed in double
        float swf = sw;
        float twf = tw;
        float zwf = zw;

        InterpSlice& p = slices[zint];
        p.zorg = floor(zout);
        p.s = (uint64_t)s;
        p.t = (uint64_t)t;
        p.c00_w = (1.0 - swf) * (1.0 - zwf);
        p.c01_w = swf * (1.0 - zwf);
        p.c10_w = (1.0 - twf) * zwf;
        p.c11_w = twf * zwf;
    }
}

// the 4-column mode of 8 rows from row y (n rows below floor(ny / 16) * 16 add pairwise)
inline int pair_mode(const uint64_t y, const uint64_t ny16) {
    if (y + 8 <= ny16) return 1;
    if (y >= ny16) return 0;
    return (int)(ny16 - y) + 2;
}

// the grid of the fused kernel: indexing and the interpolated slices' parameters
struct Grid {
    uint64_t ntotal = 0, nint = 0, xa = 0;
    std::vector<uint64_t> interp_index;
    std::unique_ptr<bool[]> is_original;
    std::vector<InterpSlice> slices;
};

inline void make_grid(const uint64_t nz, const double xstep, const double stepsize,
                      const uint8_t Reverse, Grid& G) {
    G.ntotal = (uint64_t)interp_grid_size(nz, stepsize);
    G.interp_index.assign(G.ntotal, 0);
    G.is_original.reset(new bool[G.ntotal > 0 ? G.ntotal : 1]);
    determine_interpolation_indexing(G.ntotal, stepsize, G.interp_index.data(), G.is_original.get(),
                                     G.nint);
    G.slices.assign(G.nint > 0 ? G.nint : 1, InterpSlice{});
    skewed_space_interp_slices(xstep, stepsize, Reverse, G.interp_index.data(), G.is_original.get(),
                               G.ntotal, G.slices.data());
    G.xa = (uint64_t)ceil(xstep);
}

}  // namespace pk5d::deskew
