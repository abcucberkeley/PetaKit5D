// Skewed space interpolation together with the deskew/rotate warp
// (skewed_space_interp_volume_deskew_rotate_warp_mex, the fused kernel) on C-ordered (z, y, x)
// volumes: no interpolated volume; a unit computes the values of the interpolated grid slices
// its lines read (columns col, col + 1 of its slices, 8 rows) from the raw values it copied.
// Each voxel, interpolated or output, has the mex's operands, weights and arithmetic.
#pragma once

#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

#include "lanes.hpp"
#include "plan.hpp"
#include "warp.hpp"

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
        // before flooring either way, which the column ranges of interp_column rely on)
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

// The values of column r of an interpolated slice for 8 lanes (rows): raw(z, c) gives the
// lanes of raw slice z, column c. Two raw columns at the edges, four inside (the mex's row
// ranges; "rows" there are columns here). mode4: the 4-column sum, 0 SEQ, 1 PAIR, or n + 2
// for PAIR in the first n lanes and SEQ after (the rows below floor(ny / 16) * 16 add
// pairwise); round16: uint16 values (rounded).
// cw: the slice's weights c00, c01, c10, c11 in lanes
template <typename Raw>
inline V8 interp_column(const Raw& raw, const InterpSlice& p, const V8 (&cw)[4], const uint64_t r,
                        const uint64_t xa, const uint8_t Reverse, const uint64_t nx, const int mode4,
                        const bool round16) {
    const uint64_t zorg = p.zorg, s = p.s, t = p.t;
    const V8 c00 = cw[0], c01 = cw[1], c10 = cw[2], c11 = cw[3];
    V8 v;
    bool four = false;
    V8 a, b, c, d;
    if (Reverse) {
        if (r < xa - s) {
            v = sum2(raw(zorg + 1, r + t), raw(zorg + 1, r + t + 1), c10, c11);
        } else if (r < nx - t - 1) {
            a = raw(zorg, r + s - xa);
            b = raw(zorg, r + s - xa + 1);
            c = raw(zorg + 1, r + t);
            d = raw(zorg + 1, r + t + 1);
            four = true;
        } else {
            v = sum2(raw(zorg, r + s - xa), raw(zorg, r + s - xa + 1), c00, c01);
        }
    } else {
        if (r < xa - s) {
            v = sum2(raw(zorg, r + s), raw(zorg, r + s + 1), c00, c01);
        } else if (r < nx - t - 1) {
            a = raw(zorg, r + s);
            b = raw(zorg, r + s + 1);
            c = raw(zorg + 1, r + t - xa);
            d = raw(zorg + 1, r + t - xa + 1);
            four = true;
        } else {
            v = sum2(raw(zorg + 1, r + t - xa), raw(zorg + 1, r + t - xa + 1), c10, c11);
        }
    }
    if (four) {
        if (mode4 == 0)
            v = seq4(a, b, c, d, c00, c01, c10, c11);
        else if (mode4 == 1)
            v = pair4(a, b, c, d, c00, c01, c10, c11);
        else
            v = first_lanes(pair4(a, b, c, d, c00, c01, c10, c11),
                            seq4(a, b, c, d, c00, c01, c10, c11), mode4 - 2);
    }
    return round16 ? round_u16(v) : v;
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

// the fused kernel. kind: 0 single, 1 uint16 to single, 2 uint16 to uint16 (interpolated
// values rounded to uint16, as the mex's uint16 cache)
template <typename T, typename TO>
void skewed_space_interp_warp(const T* const vol, const uint64_t nz, const uint64_t ny,
                              const uint64_t nx, const double xstep, const double stepsize,
                              const uint8_t Reverse, const double (&tmat)[4][4], TO* const out,
                              const uint64_t oz, const uint64_t oy, const uint64_t ox,
                              const float start_y, const float start_x, const float start_z,
                              const uint64_t line_len, const int nthreads) {
    if (line_len == 0) return;
    Grid G;
    make_grid(nz, xstep, stepsize, Reverse, G);
    Geometry GE{{}, nx, G.ntotal, start_x, start_z};
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) GE.tmat[i][j] = tmat[i][j];
    Slices SL;
    compute_slices(GE, ox, oz, SL, nthreads);
    const std::vector<Group> groups = make_groups(SL, GROUP_SPAN, GROUP_SLICES);
    if (groups.empty()) return;

    constexpr int kind = warp_kind<T, TO>();
    const int64_t xa = (int64_t)G.xa;
    const uint64_t y0 = (uint64_t)floor(start_y);
    const uint64_t len16 = line_len - (line_len % 16);
    const uint64_t ny16 = ny - (ny % 16);
    const uint64_t nblocks = (line_len + LANES - 1) / LANES;
    const uint64_t nchunks = (nblocks + UNIT_BLOCKS - 1) / UNIT_BLOCKS;
    const uint64_t nunits = groups.size() * nchunks;
    const uint64_t nyx = ny * nx;
    const uint64_t* idx = G.interp_index.data();
    const bool* orig = G.is_original.get();

    // per group: the raw slices and the interpolated slices its grid slices need
    struct Need {
        int64_t rlo, rhi, ilo, ihi;
    };
    std::vector<Need> need(groups.size());
    uint64_t sraw = 0, sint = 0;
    for (size_t gi = 0; gi < groups.size(); ++gi) {
        const Group& g = groups[gi];
        Need n{INT64_MAX, -1, INT64_MAX, -1};
        for (int64_t q = g.glo; q <= g.ghi; ++q) {
            if (orig[q]) {
                n.rlo = std::min<int64_t>(n.rlo, idx[q]);
                n.rhi = std::max<int64_t>(n.rhi, idx[q]);
            } else {
                const InterpSlice& p = G.slices[idx[q]];
                n.rlo = std::min<int64_t>(n.rlo, p.zorg);
                n.rhi = std::max<int64_t>(n.rhi, p.zorg + 1);
                n.ilo = std::min<int64_t>(n.ilo, idx[q]);
                n.ihi = std::max<int64_t>(n.ihi, idx[q]);
            }
        }
        need[gi] = n;
        int64_t clo, chi;
        group_columns(g, xa + 1, nx, clo, chi);
        sraw = std::max(sraw, (uint64_t)(n.rhi - n.rlo + 1) * (chi - clo + 1) * LANES);
        if (n.ihi >= 0)
            sint = std::max(sint, (uint64_t)(n.ihi - n.ilo + 1) * (g.cmax - g.cmin + 2) * LANES);
    }

#pragma omp parallel num_threads(nthreads)
    {
        std::vector<float> S(sraw), I(sint > 0 ? sint : LANES), Tl(ox * LANES);
        std::vector<TO> R(ox * LANES);
        GroupLines GL;

#pragma omp for schedule(dynamic, 1)
        for (uint64_t u = 0; u < nunits; ++u) {
            struct Fence {
                ~Fence() { stream_fence(); }
            } fence;
            const size_t gi = u / nchunks;
            const Group& g = groups[gi];
            const Need& n = need[gi];
            group_lines(GE, SL, g, ox, GL);
            int64_t clo, chi;
            group_columns(g, xa + 1, nx, clo, chi);
            const int64_t w = chi - clo + 1;
            const int64_t wi = g.cmax + 2 - g.cmin;  // columns cmin .. cmax + 1 of interpolated slices
            auto grid = [&](uint64_t q, int64_t c) -> const float* {
                if (orig[q]) return S.data() + (((int64_t)idx[q] - n.rlo) * w + (c - clo)) * LANES;
                return I.data() + (((int64_t)idx[q] - n.ilo) * wi + (c - g.cmin)) * LANES;
            };
            const uint64_t b1 = std::min(nblocks, (u % nchunks + 1) * UNIT_BLOCKS);
            for (uint64_t b = (u % nchunks) * UNIT_BLOCKS; b < b1; ++b) {
                const uint64_t i0 = b * LANES;
                const int nl = (int)std::min<uint64_t>(LANES, line_len - i0);
                const uint64_t yrow = y0 + i0;
                const bool pair = i0 + LANES <= len16;
                // the interpolated slices' 4-column mode of each half (absolute rows y)
                const int mode4[2] = {kind == 0 ? 0 : pair_mode(yrow, ny16),
                                      kind == 0 ? 0 : pair_mode(yrow + 8, ny16)};
                copy_group(vol, nyx, nx, yrow, nl, n.rlo, n.rhi, clo, chi, S.data());
                if (n.ihi >= 0) {
                    for (int h = 0; h < LANES; h += 8) {
                        auto raw = [&](uint64_t z, uint64_t c) {
                            return load8(S.data() +
                                         (((int64_t)z - n.rlo) * w + ((int64_t)c - clo)) * LANES + h);
                        };
                        for (int64_t ii = n.ilo; ii <= n.ihi; ++ii) {
                            const InterpSlice& p = G.slices[ii];
                            const V8 cw[4] = {set8(p.c00_w), set8(p.c01_w), set8(p.c10_w),
                                              set8(p.c11_w)};
                            for (int64_t r = g.cmin; r <= g.cmax + 1; ++r)
                                store8(I.data() + ((ii - n.ilo) * wi + (r - g.cmin)) * LANES + h,
                                       interp_column(raw, p, cw, (uint64_t)r, G.xa, Reverse, nx,
                                                     mode4[h / 8], kind == 2));
                        }
                    }
                }
                for (uint64_t zo = g.zo_begin; zo < g.zo_end; ++zo) {
                    if (SL.col[zo] < 0) continue;
                    const int64_t c = SL.col[zo];
                    const uint64_t k0 = (zo - g.zo_begin) * ox;
                    for (uint64_t xo = 0; xo < ox; ++xo) {
                        const uint64_t k = k0 + xo;
                        const int64_t z0 = GL.z0[k];
                        float* t = Tl.data() + xo * LANES;
                        if (z0 < 0) {
                            store8(t, zero8());
                            store8(t + 8, zero8());
                            continue;
                        }
                        const float* pa = grid(z0, c);
                        const float* pb = grid(z0 + 1, c);
                        const V8 w00 = set8(GL.w[4 * k]), w01 = set8(GL.w[4 * k + 1]);
                        const V8 w10 = set8(GL.w[4 * k + 2]), w11 = set8(GL.w[4 * k + 3]);
                        // the line sum the mex picks for the type and, for uint16 to single,
                        // the slice types (original or interpolated)
                        int f = 0;  // SEQ
                        if (kind == 2) {
                            f = pair ? 1 : 0;
                        } else if (kind == 1) {
                            const bool oa = orig[z0], ob = orig[z0 + 1];
                            if (oa || ob) f = pair ? 1 : (oa ? 0 : 2);
                        }
                        line16(pa, pb, t, [&](V8 a, V8 b, V8 cc, V8 d) {
                            V8 v = f == 1   ? pair4(a, b, cc, d, w00, w01, w10, w11)
                                   : f == 0 ? seq4(a, b, cc, d, w00, w01, w10, w11)
                                            : seq4_swap(a, b, cc, d, w00, w01, w10, w11);
                            return kind == 2 ? add_half(v) : v;
                        });
                    }
                    write_lanes(Tl.data(), ox, nl, out + zo * oy * ox + i0 * ox, R.data());
                }
            }
        }
    }
}

}  // namespace pk5d::deskew
