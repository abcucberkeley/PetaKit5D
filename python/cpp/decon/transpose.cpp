// Plane transposes between the port's (z, y, x) layout and MATLAB's (z, x, y): copies, so the values
// don't change; blocked so that both the reads and the writes stay in cache lines.
#include <cstring>
#include <stdexcept>

#include "decon.hpp"

namespace pk5d::decon {
namespace {

constexpr idx_t kTile = 32;

template <typename T>
void transpose_t(const T* src, T* dst, idx_t nz, idx_t rows, idx_t cols, int nt) {
    const idx_t tr = (rows + kTile - 1) / kTile, tc = (cols + kTile - 1) / kTile;
    const idx_t ntask = nz * tr;
    const idx_t plane = rows * cols;
#pragma omp parallel for num_threads(nt) schedule(static)
    for (idx_t task = 0; task < ntask; task++) {
        const idx_t z = task / tr;
        const idx_t r0 = (task % tr) * kTile;
        const idx_t r1 = r0 + kTile < rows ? r0 + kTile : rows;
        const T* s = src + z * plane;
        T* d = dst + z * plane;
        for (idx_t b = 0; b < tc; b++) {
            const idx_t c0 = b * kTile;
            const idx_t c1 = c0 + kTile < cols ? c0 + kTile : cols;
            for (idx_t c = c0; c < c1; c++) {
                T* dc = d + c * rows;
                for (idx_t r = r0; r < r1; r++) dc[r] = s[r * cols + c];
            }
        }
    }
}

}  // namespace

void transpose_planes(const void* src, void* dst, idx_t nz, idx_t rows, idx_t cols, int bytes, int nt) {
    if (nz * rows * cols == 0) return;
    if (rows == 1 || cols == 1) {
        std::memcpy(dst, src, static_cast<size_t>(nz * rows * cols) * static_cast<size_t>(bytes));
        return;
    }
    switch (bytes) {
        case 1:
            transpose_t(static_cast<const uint8_t*>(src), static_cast<uint8_t*>(dst), nz, rows, cols, nt);
            break;
        case 2:
            transpose_t(static_cast<const uint16_t*>(src), static_cast<uint16_t*>(dst), nz, rows, cols, nt);
            break;
        case 4:
            transpose_t(static_cast<const uint32_t*>(src), static_cast<uint32_t*>(dst), nz, rows, cols, nt);
            break;
        case 8:
            transpose_t(static_cast<const uint64_t*>(src), static_cast<uint64_t*>(dst), nz, rows, cols, nt);
            break;
        default:
            throw std::invalid_argument("transpose_planes: elements of 1, 2, 4 or 8 bytes");
    }
}

}  // namespace pk5d::decon
