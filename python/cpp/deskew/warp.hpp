// The deskew/rotate warp without interpolation in skewed space (volume_deskew_rotate_warp_mex)
// on C-ordered (z, y, x) volumes: single to single, uint16 to single, uint16 to uint16.
// Each output voxel has the mex's operands, weights and arithmetic (plan.hpp, lanes.hpp).
#pragma once

#include <cstdint>
#include <vector>

#include "lanes.hpp"
#include "plan.hpp"

namespace pk5d::deskew {

// the arithmetic of the line sums: 0 single (SEQ), 1 uint16 to single, 2 uint16 to uint16
template <typename T, typename TO>
constexpr int warp_kind() {
    return sizeof(T) == 4 ? 0 : (sizeof(TO) == 4 ? 1 : 2);
}

// group span (columns) and slices: a group's scratch is (its grid slices) x (span + 2) x 16
constexpr int64_t GROUP_SPAN = 30;
constexpr uint64_t GROUP_SLICES = 64;

// the scratch of a group's columns [clo, chi] for its slices z in [zlo, zhi]: lanes of (z, c)
// at S + ((z - zlo) * w + (c - clo)) * LANES
template <typename T>
inline void copy_group(const T* vol, uint64_t nyx, uint64_t nx, uint64_t yrow, int nl, int64_t zlo,
                       int64_t zhi, int64_t clo, int64_t chi, float* S) {
    const int64_t w = chi - clo + 1;
    for (int64_t z = zlo; z <= zhi; ++z)
        copy_lanes(vol + z * nyx, nx, yrow, nl, clo, chi, S + (z - zlo) * w * LANES);
}

// 16 voxels of one output line from the lanes of a (z0, col), b (z0 + 1, col), c (z0, col + 1),
// d (z0 + 1, col + 1); f(a, b, c, d) is the line's formula
template <typename F>
inline void line16(const float* pa, const float* pb, float* out, F f) {
    for (int h = 0; h < LANES; h += 8)
        store8(out + h, f(load8(pa + h), load8(pb + h), load8(pa + LANES + h), load8(pb + LANES + h)));
}

// vol: nz x ny x nx; out: oz x oy x ox (zeros). Line i of output slice zo, line xo is output
// voxel (zo, i, xo); it reads input rows y = y0 + i, i < line_len (y0 = floor(start_y), the
// mex's startX)
template <typename T, typename TO>
void warp_volume(const T* const vol, const uint64_t nz, const uint64_t ny, const uint64_t nx,
                 const double (&tmat)[4][4], TO* const out, const uint64_t oz, const uint64_t oy,
                 const uint64_t ox, const float start_y, const float start_x, const float start_z,
                 const uint64_t line_len, const int nthreads) {
    // a bbox past the input's rows: nothing to compute
    if (line_len == 0) return;
    Geometry G{{}, nx, nz, start_x, start_z};
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) G.tmat[i][j] = tmat[i][j];
    Slices SL;
    compute_slices(G, ox, oz, SL, nthreads);
    const std::vector<Group> groups = make_groups(SL, GROUP_SPAN, GROUP_SLICES);
    if (groups.empty()) return;

    constexpr int kind = warp_kind<T, TO>();
    const uint64_t y0 = (uint64_t)floor(start_y);
    const uint64_t len16 = line_len - (line_len % 16);
    const uint64_t nblocks = (line_len + LANES - 1) / LANES;
    const uint64_t nchunks = (nblocks + UNIT_BLOCKS - 1) / UNIT_BLOCKS;
    const uint64_t nunits = groups.size() * nchunks;
    const uint64_t nyx = ny * nx;
    uint64_t scratch = 0;
    for (const Group& g : groups) {
        int64_t clo, chi;
        group_columns(g, 0, nx, clo, chi);
        scratch = std::max(scratch, (uint64_t)(g.ghi - g.glo + 1) * (chi - clo + 1) * LANES);
    }

#pragma omp parallel num_threads(nthreads)
    {
        std::vector<float> S(scratch), Tl(ox * LANES);
        std::vector<TO> R(ox * LANES);
        GroupLines GL;

#pragma omp for schedule(dynamic, 1)
        for (uint64_t u = 0; u < nunits; ++u) {
            struct Fence {
                ~Fence() { stream_fence(); }
            } fence;
            const Group& g = groups[u / nchunks];
            group_lines(G, SL, g, ox, GL);
            int64_t clo, chi;
            group_columns(g, 0, nx, clo, chi);
            const int64_t w = chi - clo + 1;
            const uint64_t b1 = std::min(nblocks, (u % nchunks + 1) * UNIT_BLOCKS);
            for (uint64_t b = (u % nchunks) * UNIT_BLOCKS; b < b1; ++b) {
                const uint64_t i0 = b * LANES;
                const int nl = (int)std::min<uint64_t>(LANES, line_len - i0);
                // floor(len / 16) * 16: a block is all PAIR or all SEQ
                const bool pair = kind != 0 && i0 + LANES <= len16;
                copy_group(vol, nyx, nx, y0 + i0, nl, g.glo, g.ghi, clo, chi, S.data());
                for (uint64_t zo = g.zo_begin; zo < g.zo_end; ++zo) {
                    if (SL.col[zo] < 0) continue;
                    const int64_t c = SL.col[zo] - clo;
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
                        const float* pa = S.data() + ((z0 - g.glo) * w + c) * LANES;
                        const float* pb = pa + w * LANES;
                        const V8 w00 = set8(GL.w[4 * k]), w01 = set8(GL.w[4 * k + 1]);
                        const V8 w10 = set8(GL.w[4 * k + 2]), w11 = set8(GL.w[4 * k + 3]);
                        if (pair) {
                            line16(pa, pb, t, [&](V8 a, V8 b, V8 cc, V8 d) {
                                V8 v = pair4(a, b, cc, d, w00, w01, w10, w11);
                                return kind == 2 ? add_half(v) : v;
                            });
                        } else {
                            line16(pa, pb, t, [&](V8 a, V8 b, V8 cc, V8 d) {
                                V8 v = seq4(a, b, cc, d, w00, w01, w10, w11);
                                return kind == 2 ? add_half(v) : v;
                            });
                        }
                    }
                    write_lanes(Tl.data(), ox, nl, out + zo * oy * ox + i0 * ox, R.data());
                }
            }
        }
    }
}

}  // namespace pk5d::deskew
