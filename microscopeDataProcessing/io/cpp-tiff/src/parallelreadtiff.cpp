#include <cstdint>
#include <cmath>
#include <cstring>
#include <vector>
#include <atomic>
#include <memory>
#include <algorithm>
#include <limits.h>
#include <omp.h>
#include "tiffio.h"
#include "../src/helperfunctions.h"

// RAII for libtiff handles: TIFFClose runs on every exit path, so the worker
// functions don't carry a naked TIFF* with a manual close. (C++11: no make_unique.)
struct TiffCloser { void operator()(TIFF* t) const noexcept { if(t) TIFFClose(t); } };
using TiffPtr = std::unique_ptr<TIFF, TiffCloser>;


// Cache-blocked x/y transpose of one decoded strip into the column-major
// (MATLAB) output. The decoded strip `buf` is row-major (buf[j + k*x]); the
// output is column-major (out[j*y + (k+rowOff) + sliceOff]). A plain nested
// loop has a large stride (x or y) on whichever array is innermost, so one
// array is streamed through cache one element per line. Tiling by TILE keeps a
// TILE x TILE block of BOTH arrays resident, turning the strided accesses into
// cache hits. j = output column (input x), k = output row within strip (input y).
template <typename T>
static inline void flipTransposeStrip(T* out, const T* buf, int64_t x, int64_t y,
                                      int64_t stripRows, int64_t rowOff, int64_t sliceOff) {
    constexpr int64_t TILE = 128;
    for (int64_t jt = 0; jt < x; jt += TILE) {
        const int64_t jEnd = (jt + TILE < x) ? jt + TILE : x;
        for (int64_t kt = 0; kt < stripRows; kt += TILE) {
            const int64_t kEnd = (kt + TILE < stripRows) ? kt + TILE : stripRows;
            for (int64_t j = jt; j < jEnd; j++) {
                T* o = out + (j * y) + rowOff + sliceOff;   // contiguous in k
                const T* b = buf + j;                       // strided in k (step x)
                for (int64_t k = kt; k < kEnd; k++) o[k] = b[k * x];
            }
        }
    }
}

// 8- and 16-bit RGB: for each block of JB pixels, KB rows at a time are copied
// from the chunky rows (contiguous reads) into a small local block, which is
// then written out as 128 bytes of consecutive rows per output column. Reading
// one sample at a time straight down a column instead stepped 6 KB per read for
// 2048-pixel 8-bit RGB rows, and how fast that ran depended on where the loop
// landed in the library (5-7% between builds); blocked, those reads are
// 1.03-1.25x faster than before. For 32- and 64-bit samples it was no faster,
// and 8-bit RGBA reads with every core busy were 20% slower.
template <typename T, int S>
static void flipTransposeStripColor(T* out, const T* buf, int64_t x, int64_t y,
                                    int64_t stripRows, int64_t rowOff, int64_t sliceOff) {
    constexpr int64_t KB = 128 / (int64_t)sizeof(T), JB = 16;
    T block[S][JB][KB];
    const int64_t plane = x * y, inRow = x * S;
    for (int64_t jt = 0; jt < x; jt += JB) {
        const int64_t jn = (jt + JB < x) ? JB : x - jt;
        for (int64_t kt = 0; kt < stripRows; kt += KB) {
            const int64_t kn = (kt + KB < stripRows) ? KB : stripRows - kt;
            for (int64_t k = 0; k < kn; k++) {
                const T* row = buf + (kt + k) * inRow + jt * S;
                for (int64_t j = 0; j < jn; j++)
                    for (int s = 0; s < S; s++) block[s][j][k] = row[j * S + s];
            }
            for (int s = 0; s < S; s++)
                for (int64_t j = 0; j < jn; j++)
                    memcpy(out + sliceOff + s * plane + (jt + j) * y + rowOff + kt, block[s][j], (size_t)kn * sizeof(T));
        }
    }
}

// RGB/multi-sample variant: de-interleaves a chunky strip (input
// buf[(k*x + j)*samples + c]) into the column-major (y, x, c, z) MATLAB output
// while transposing x/y. Only used when samples>1; the grayscale path keeps
// using flipTransposeStrip above, so its hot loop is unchanged.
template <typename T>
static inline void flipTransposeStripRGB(T* out, const T* buf, int64_t x, int64_t y,
                                         int64_t stripRows, int64_t rowOff, int64_t sliceOff,
                                         int64_t samples) {
    if (sizeof(T) <= 2 && samples == 3) { flipTransposeStripColor<T, 3>(out, buf, x, y, stripRows, rowOff, sliceOff); return; }
    constexpr int64_t TILE = 128;
    const int64_t chanStride = x * y;     // output channel-plane stride
    const int64_t inRow = x * samples;    // input (chunky) row stride, in elements
    for (int64_t c = 0; c < samples; c++) {
        T* outc = out + c * chanStride + sliceOff;
        const T* bufc = buf + c;
        for (int64_t jt = 0; jt < x; jt += TILE) {
            const int64_t jEnd = (jt + TILE < x) ? jt + TILE : x;
            for (int64_t kt = 0; kt < stripRows; kt += TILE) {
                const int64_t kEnd = (kt + TILE < stripRows) ? kt + TILE : stripRows;
                for (int64_t j = jt; j < jEnd; j++) {
                    T* o = outc + (j * y) + rowOff;   // contiguous in k
                    const T* b = bufc + j * samples;  // strided in k (step inRow)
                    for (int64_t k = kt; k < kEnd; k++) o[k] = b[k * inRow];
                }
            }
        }
    }
}

// Number of consecutive strips to decode before each flip transpose. One
// transpose writes a few elements into every output column, so its cost per
// column (a cache line, and a page once columns are 4 KB apart) is paid again
// for every strip. Strips of only a row or two are common: libtiff splits a
// single large uncompressed strip into one-row strips, and many writers use
// ~8 KB strips. Transposing runs of about 16 rows made flip reads of one- and
// two-row strips 2-12x faster (Xeon and Apple M2, 1-48 threads); longer runs
// gained little and were slower with 8-row strips. Strips of 16+ rows are
// transposed one at a time as before. The byte cap bounds the per-thread
// buffer for very wide rows.
static int64_t stripsPerTranspose(uint64_t stripSize, uint64_t rowBytes){
    const uint64_t targetRows = 16;
    const uint64_t maxBytes = (uint64_t)4 << 20;
    const uint64_t rows = std::min<uint64_t>(targetRows, std::max<uint64_t>(1, maxBytes / std::max<uint64_t>(1, rowBytes)));
    return (int64_t)std::max<uint64_t>(1, rows / std::max<uint64_t>(1, stripSize));
}

// Whether the strips of the current directory can be copied straight from the
// file: no compression, normal bit order, and libtiff's strip size equal to
// rows x row bytes (so no subsampling or padding to undo).
static bool rawStripsOK(TIFF* tif, uint64_t stripSize, uint64_t rowBytes){
    uint16_t compression = COMPRESSION_NONE, fillOrder = FILLORDER_MSB2LSB;
    TIFFGetFieldDefaulted(tif, TIFFTAG_COMPRESSION, &compression);
    TIFFGetFieldDefaulted(tif, TIFFTAG_FILLORDER, &fillOrder);
    return compression == COMPRESSION_NONE && fillOrder == FILLORDER_MSB2LSB && !TIFFIsTiled(tif)
           && (uint64_t)TIFFStripSize(tif) == stripSize*rowBytes;
}

// Bytes per read when a flip read copies uncompressed strips into its buffer.
static const uint64_t rawFlipBlockBytes = (uint64_t)1 << 20;

// Read strips [first, last) of the current directory with one call when they
// are stored back to back with exactly the expected sizes, swapping bytes the
// way libtiff would. libtiff splits a large uncompressed strip into ~8 KB
// strips and reads each one separately; on NFS every read call has a fixed
// cost, which made stacks of small uncompressed slices up to 1.9x slower than
// v1.6.2's memory-mapped reads. base is the first strip of the plane being
// read (0 unless the file is planar). Returns false if the layout does not
// match or the read comes up short, so the caller falls back to libtiff strip
// by strip, which also reports real errors.
static bool readRawStrips(TIFF* tif, int64_t first, int64_t last, int64_t base, uint64_t stripSize, uint64_t y,
                          uint64_t rowBytes, uint64_t bits, uint8_t* dst){
    uint64_t *offsets = NULL, *counts = NULL;
    if(!TIFFGetField(tif, TIFFTAG_STRIPOFFSETS, &offsets) || !TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &counts)
       || !offsets || !counts || (uint64_t)last > TIFFNumberOfStrips(tif)) return false;
    uint64_t total = 0;
    for(int64_t i = first; i < last; i++){
        const uint64_t rows = std::min<uint64_t>(stripSize, y - (uint64_t)(i - base)*stripSize);
        if(counts[i] != rows*rowBytes || offsets[i] != offsets[first] + total) return false;
        total += counts[i];
    }
    thandle_t h = TIFFClientdata(tif);
    if(TIFFGetSeekProc(tif)(h, (toff_t)offsets[first], SEEK_SET) != (toff_t)offsets[first]) return false;
    if((uint64_t)TIFFGetReadProc(tif)(h, dst, (tmsize_t)total) != total) return false;
    if(TIFFIsByteSwapped(tif)){
        switch(bits){
            case 16: TIFFSwabArrayOfShort(reinterpret_cast<uint16_t*>(dst), (tmsize_t)(total/2)); break;
            case 32: TIFFSwabArrayOfLong(reinterpret_cast<uint32_t*>(dst), (tmsize_t)(total/4)); break;
            case 64: TIFFSwabArrayOfLong8(reinterpret_cast<uint64_t*>(dst), (tmsize_t)(total/8)); break;
        }
    }
    return true;
}

// ---- Tiled TIFF support -------------------------------------------------
// A tile is a sub-rectangle rather than a full-width strip, so its rows are not
// contiguous in the output (the no-flip path cannot decode straight into the
// output like it does for strips) and libtiff pads edge tiles out to a full
// tile. validW/validH bound every copy to the part inside the image.

// Copy one decoded tile into the row-major (no-flip) output, row by row.
static inline void copyTileNoFlip(uint8_t* out, const uint8_t* buf, int64_t x,
                                  int64_t tileX0, int64_t tileY0,
                                  int64_t validW, int64_t validH,
                                  int64_t tw, int64_t sliceOff,
                                  int64_t samples, int64_t bytes) {
    const int64_t px = samples * bytes;      // bytes per pixel
    const int64_t inRow = tw * px;           // padded tile row stride
    const int64_t outRow = x * px;           // image row stride
    const int64_t copyBytes = validW * px;
    uint8_t* o = out + sliceOff + ((tileY0 * x) + tileX0) * px;
    const uint8_t* b = buf;
    for (int64_t r = 0; r < validH; r++) {
        memcpy(o, b, (size_t)copyBytes);
        o += outRow;
        b += inRow;
    }
}

// Tile equivalent of flipTransposeStrip: transpose one tile into the
// column-major (y, x, z) MATLAB output. j = output column (input x),
// k = output row within the tile (input y).
template <typename T>
static inline void flipTransposeTile(T* out, const T* buf, int64_t y,
                                     int64_t tileX0, int64_t tileY0,
                                     int64_t validW, int64_t validH,
                                     int64_t tw, int64_t sliceOff) {
    constexpr int64_t TILE = 128;
    for (int64_t jt = 0; jt < validW; jt += TILE) {
        const int64_t jEnd = (jt + TILE < validW) ? jt + TILE : validW;
        for (int64_t kt = 0; kt < validH; kt += TILE) {
            const int64_t kEnd = (kt + TILE < validH) ? kt + TILE : validH;
            for (int64_t j = jt; j < jEnd; j++) {
                T* o = out + ((tileX0 + j) * y) + tileY0 + sliceOff;  // contiguous in k
                const T* b = buf + j;                                 // strided in k (step tw)
                for (int64_t k = kt; k < kEnd; k++) o[k] = b[k * tw];
            }
        }
    }
}

// Tile equivalent of flipTransposeStripRGB (chunky multi-sample tiles).
template <typename T>
static inline void flipTransposeTileRGB(T* out, const T* buf, int64_t x, int64_t y,
                                        int64_t tileX0, int64_t tileY0,
                                        int64_t validW, int64_t validH,
                                        int64_t tw, int64_t sliceOff, int64_t samples) {
    constexpr int64_t TILE = 128;
    const int64_t chanStride = x * y;        // output channel-plane stride
    const int64_t inRow = tw * samples;      // input (chunky) tile row stride
    for (int64_t c = 0; c < samples; c++) {
        T* outc = out + c * chanStride + sliceOff;
        const T* bufc = buf + c;
        for (int64_t jt = 0; jt < validW; jt += TILE) {
            const int64_t jEnd = (jt + TILE < validW) ? jt + TILE : validW;
            for (int64_t kt = 0; kt < validH; kt += TILE) {
                const int64_t kEnd = (kt + TILE < validH) ? kt + TILE : validH;
                for (int64_t j = jt; j < jEnd; j++) {
                    T* o = outc + ((tileX0 + j) * y) + tileY0;
                    const T* b = bufc + j * samples;
                    for (int64_t k = kt; k < kEnd; k++) o[k] = b[k * inRow];
                }
            }
        }
    }
}

// ---- Planar color -------------------------------------------------------
// A planar (PLANARCONFIG=2) file stores each color as its own image: every
// strip or tile of the first color, then of the next, each holding one sample
// per pixel. The output layouts are the same as for chunky files. The flip
// (MATLAB) layout already keeps each color in its own (y, x) plane, so a
// planar strip or tile is transposed into it like a grayscale one. The no-flip
// layout keeps each pixel's colors together, so the planes are interleaved.

// Flip: transpose `rows` decoded rows of color p, starting at image row
// rowOff, into that color's plane of output slice sliceIdx.
static void flipPlanarRows(void* out, const uint8_t* buf, int64_t x, int64_t y, int64_t bits, int64_t samples,
                           int64_t rowOff, int64_t rows, int64_t sliceIdx, int64_t p){
    const int64_t planeOff = (sliceIdx*samples + p)*x*y;
    switch(bits){
        case 8:  flipTransposeStrip<uint8_t>(static_cast<uint8_t*>(out), buf, x, y, rows, rowOff, planeOff); break;
        case 16: flipTransposeStrip<uint16_t>(static_cast<uint16_t*>(out), reinterpret_cast<const uint16_t*>(buf), x, y, rows, rowOff, planeOff); break;
        case 32: flipTransposeStrip<uint32_t>(static_cast<uint32_t*>(out), reinterpret_cast<const uint32_t*>(buf), x, y, rows, rowOff, planeOff); break;
        case 64: flipTransposeStrip<uint64_t>(static_cast<uint64_t*>(out), reinterpret_cast<const uint64_t*>(buf), x, y, rows, rowOff, planeOff); break;
    }
}

// With the number of colors fixed at compile time the compiler writes whole
// pixels at once: RGBA and 16- and 32-bit RGB interleave 1.7-10x faster than
// one color at a time.
template <typename T, int S>
static void interleaveFixed(T* out, const T* buf, uint64_t x, uint64_t rows, uint64_t planeStride){
    for(uint64_t r = 0; r < rows; r++){
        T* o = out + r*x*S;
        const T* b = buf + r*x;
        for(uint64_t j = 0; j < x; j++)
            for(int p = 0; p < S; p++) o[j*S + p] = b[p*planeStride + j];
    }
}

// 8-bit RGB, which the compiler does not vectorize: four pixels at a time as
// three 32-bit words, 2.7x faster than a byte at a time.
static void interleaveRGB8(uint8_t* out, const uint8_t* buf, uint64_t x, uint64_t rows, uint64_t planeStride){
    for(uint64_t r = 0; r < rows; r++){
        uint8_t* o = out + r*x*3;
        const uint8_t *R = buf + r*x, *G = R + planeStride, *B = G + planeStride;
        uint64_t j = 0;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
        for(; j + 4 <= x; j += 4){
            uint32_t r4, g4, b4;
            memcpy(&r4, R + j, 4); memcpy(&g4, G + j, 4); memcpy(&b4, B + j, 4);
            const uint32_t w0 = (r4 & 0xFF) | (g4 & 0xFF) << 8 | (b4 & 0xFF) << 16 | (r4 >> 8 & 0xFF) << 24;
            const uint32_t w1 = (g4 >> 8 & 0xFF) | (b4 >> 8 & 0xFF) << 8 | (r4 >> 16 & 0xFF) << 16 | (g4 >> 16 & 0xFF) << 24;
            const uint32_t w2 = (b4 >> 16 & 0xFF) | (r4 >> 24) << 8 | (g4 >> 24) << 16 | (b4 >> 24) << 24;
            memcpy(o + j*3, &w0, 4); memcpy(o + j*3 + 4, &w1, 4); memcpy(o + j*3 + 8, &w2, 4);
        }
#endif
        for(; j < x; j++){ o[j*3] = R[j]; o[j*3 + 1] = G[j]; o[j*3 + 2] = B[j]; }
    }
}

template <typename T>
static void interleavePlanes(T* out, const T* buf, uint64_t x, uint64_t rows, uint64_t samples, uint64_t planeStride){
    if(samples == 3){ interleaveFixed<T, 3>(out, buf, x, rows, planeStride); return; }
    if(samples == 4){ interleaveFixed<T, 4>(out, buf, x, rows, planeStride); return; }
    for(uint64_t r = 0; r < rows; r++){
        T* o = out + r*x*samples;
        for(uint64_t p = 0; p < samples; p++){
            const T* b = buf + p*planeStride + r*x;
            for(uint64_t j = 0; j < x; j++) o[j*samples + p] = b[j];
        }
    }
}

// No-flip: interleave `rows` decoded rows of every color (color p's rows at
// buf + p*planeBytes), starting at image row rowOff, into output slice sliceIdx.
static void interleavePlanarRows(void* out, const uint8_t* buf, uint64_t x, uint64_t y, uint64_t bits, uint64_t samples,
                                 uint64_t rowOff, uint64_t rows, uint64_t sliceIdx, uint64_t planeBytes){
    const uint64_t first = (sliceIdx*y + rowOff)*x*samples;   // in samples
    switch(bits){
        case 8:
            if(samples == 3) interleaveRGB8(static_cast<uint8_t*>(out) + first, buf, x, rows, planeBytes);
            else interleavePlanes<uint8_t>(static_cast<uint8_t*>(out) + first, buf, x, rows, samples, planeBytes);
            break;
        case 16: interleavePlanes<uint16_t>(static_cast<uint16_t*>(out) + first, reinterpret_cast<const uint16_t*>(buf), x, rows, samples, planeBytes/2); break;
        case 32: interleavePlanes<uint32_t>(static_cast<uint32_t*>(out) + first, reinterpret_cast<const uint32_t*>(buf), x, rows, samples, planeBytes/4); break;
        case 64: interleavePlanes<uint64_t>(static_cast<uint64_t*>(out) + first, reinterpret_cast<const uint64_t*>(buf), x, rows, samples, planeBytes/8); break;
    }
}

// Read strips [sBegin, sEnd) of every color of the current directory into
// output slice sliceIdx, `step` strips at a time. buffer holds one block of
// blockBytes (flip), or one per color, blockBytes apart (no-flip). raw reads
// back-to-back uncompressed strips with one call (see readRawStrips). Returns
// the strip that could not be read, or -1.
static int64_t readPlanarStrips(TIFF* tif, uint8_t* buffer, uint64_t blockBytes, void* out,
                                uint64_t x, uint64_t y, uint64_t bits, uint64_t samples, uint64_t stripSize,
                                int64_t sBegin, int64_t sEnd, int64_t step, bool raw, int64_t sliceIdx, uint8_t flipXY){
    const uint64_t rowBytes = x*(bits/8);
    const int64_t stripsPerPlane = (int64_t)((y + stripSize - 1)/stripSize);
    for(int64_t first = sBegin; first < sEnd; first += step){
        const int64_t last = std::min<int64_t>(first + step, sEnd);
        const int64_t rowOff = first*(int64_t)stripSize;
        const int64_t rows = std::min<int64_t>((last - first)*(int64_t)stripSize, (int64_t)y - rowOff);
        for(uint64_t p = 0; p < samples; p++){
            const int64_t base = (int64_t)p*stripsPerPlane;
            uint8_t* dst = buffer + (flipXY ? 0 : p*blockBytes);
            if(!(raw && readRawStrips(tif, base + first, base + last, base, stripSize, y, rowBytes, bits, dst)))
                for(int64_t i = first; i < last; i++)
                    if(TIFFReadEncodedStrip(tif, (uint32_t)(base + i), dst + (i - first)*stripSize*rowBytes, stripSize*rowBytes) < 0)
                        return base + i;
            if(flipXY) flipPlanarRows(out, dst, (int64_t)x, (int64_t)y, (int64_t)bits, (int64_t)samples, rowOff, rows, sliceIdx, (int64_t)p);
        }
        if(!flipXY) interleavePlanarRows(out, buffer, x, y, bits, samples, (uint64_t)rowOff, (uint64_t)rows, (uint64_t)sliceIdx, blockBytes);
    }
    return -1;
}

// Place one decoded tile of color p: transposed into that color's plane
// (flip), or spread over each pixel's samples (no-flip).
template <typename T>
static void placePlanarTileT(T* out, const T* buf, int64_t x, int64_t y, int64_t tileX0, int64_t tileY0,
                             int64_t validW, int64_t validH, int64_t tw, int64_t sliceIdx, int64_t samples,
                             int64_t p, uint8_t flipXY){
    if(flipXY){
        flipTransposeTile<T>(out, buf, y, tileX0, tileY0, validW, validH, tw, (sliceIdx*samples + p)*x*y);
        return;
    }
    for(int64_t r = 0; r < validH; r++){
        T* o = out + ((sliceIdx*y + tileY0 + r)*x + tileX0)*samples + p;
        const T* b = buf + r*tw;
        for(int64_t c = 0; c < validW; c++) o[c*samples] = b[c];
    }
}

static void placePlanarTile(void* out, const uint8_t* buf, int64_t x, int64_t y, int64_t bits, int64_t tileX0, int64_t tileY0,
                            int64_t validW, int64_t validH, int64_t tw, int64_t sliceIdx, int64_t samples, int64_t p, uint8_t flipXY){
    switch(bits){
        case 8:  placePlanarTileT<uint8_t>(static_cast<uint8_t*>(out), buf, x, y, tileX0, tileY0, validW, validH, tw, sliceIdx, samples, p, flipXY); break;
        case 16: placePlanarTileT<uint16_t>(static_cast<uint16_t*>(out), reinterpret_cast<const uint16_t*>(buf), x, y, tileX0, tileY0, validW, validH, tw, sliceIdx, samples, p, flipXY); break;
        case 32: placePlanarTileT<uint32_t>(static_cast<uint32_t*>(out), reinterpret_cast<const uint32_t*>(buf), x, y, tileX0, tileY0, validW, validH, tw, sliceIdx, samples, p, flipXY); break;
        case 64: placePlanarTileT<uint64_t>(static_cast<uint64_t*>(out), reinterpret_cast<const uint64_t*>(buf), x, y, tileX0, tileY0, validW, validH, tw, sliceIdx, samples, p, flipXY); break;
    }
}

// Decode tiles [tileStart, tileEnd) of the current directory into the output.
// Shared by the 2D and 3D readers (they differ only in how tiles are split
// across workers). For a planar file, tile t of every color is read. Returns 1
// and sets failedTile if a tile cannot be read.
static uint8_t readTilesInto(TIFF* tif, uint8_t* buffer, int64_t tileBufBytes,
                             void* out, int64_t x, int64_t y, int64_t bits,
                             int64_t samples, bool planar, int64_t tw, int64_t tl,
                             int64_t tilesAcross, int64_t tileStart, int64_t tileEnd,
                             int64_t sliceIdx, uint8_t flipXY, int64_t* failedTile) {
    const int64_t bytes = bits / 8;
    const int64_t tilesPerPlane = tilesAcross * ((y + tl - 1) / tl);
    for (int64_t t = tileStart; t < tileEnd; t++) {
        const int64_t tileY0 = (t / tilesAcross) * tl;
        const int64_t tileX0 = (t % tilesAcross) * tw;
        int64_t validH = y - tileY0; if (validH > tl) validH = tl;
        int64_t validW = x - tileX0; if (validW > tw) validW = tw;
        if (validW <= 0 || validH <= 0) continue;

        if (planar) {
            for (int64_t p = 0; p < samples; p++) {
                const int64_t tile = p * tilesPerPlane + t;
                if (TIFFReadEncodedTile(tif, (uint32_t)tile, buffer, tileBufBytes) < 0) {
                    *failedTile = tile;
                    return 1;
                }
                placePlanarTile(out, buffer, x, y, bits, tileX0, tileY0, validW, validH, tw, sliceIdx, samples, p, flipXY);
            }
            continue;
        }

        if (TIFFReadEncodedTile(tif, (uint32_t)t, buffer, tileBufBytes) < 0) {
            *failedTile = t;
            return 1;
        }
        if (!flipXY) {
            copyTileNoFlip((uint8_t*)out, buffer, x, tileX0, tileY0, validW, validH, tw,
                           sliceIdx * y * x * samples * bytes, samples, bytes);
            continue;
        }
        if (samples == 1) {
            const int64_t sliceOff = sliceIdx * x * y;
            switch (bits) {
                case 8:
                    flipTransposeTile<uint8_t>((uint8_t*)out, (const uint8_t*)buffer, y, tileX0, tileY0, validW, validH, tw, sliceOff);
                    break;
                case 16:
                    flipTransposeTile<uint16_t>((uint16_t*)out, (const uint16_t*)buffer, y, tileX0, tileY0, validW, validH, tw, sliceOff);
                    break;
                case 32:
                    flipTransposeTile<float>((float*)out, (const float*)buffer, y, tileX0, tileY0, validW, validH, tw, sliceOff);
                    break;
                case 64:
                    flipTransposeTile<double>((double*)out, (const double*)buffer, y, tileX0, tileY0, validW, validH, tw, sliceOff);
                    break;
            }
        }
        else {
            const int64_t sliceOff = sliceIdx * x * y * samples;
            switch (bits) {
                case 8:
                    flipTransposeTileRGB<uint8_t>((uint8_t*)out, (const uint8_t*)buffer, x, y, tileX0, tileY0, validW, validH, tw, sliceOff, samples);
                    break;
                case 16:
                    flipTransposeTileRGB<uint16_t>((uint16_t*)out, (const uint16_t*)buffer, x, y, tileX0, tileY0, validW, validH, tw, sliceOff, samples);
                    break;
                case 32:
                    flipTransposeTileRGB<float>((float*)out, (const float*)buffer, x, y, tileX0, tileY0, validW, validH, tw, sliceOff, samples);
                    break;
                case 64:
                    flipTransposeTileRGB<double>((double*)out, (const double*)buffer, x, y, tileX0, tileY0, validW, validH, tw, sliceOff, samples);
                    break;
            }
        }
    }
    return 0;
}

// Open the file for one reader thread. On Linux this skips libtiff's memory
// map. With every thread mapping the whole file, a read ended with each thread
// unmapping it; the kernel serializes those unmaps and they stall the other
// threads' page faults, which cost up to ~0.4 s at the end of a 2 s read.
// Reading through read() only adds a copy of the compressed strips. macOS
// measured no difference and Windows is unmeasured, so they keep the mapping.
static TIFF* openReaderTiff(const char* fileName){
#ifdef __linux__
    return TIFFOpen(fileName, "rm");
#else
    return TIFFOpen(fileName, "r");
#endif
}

// Place one decoded scanline (row `row` of output slice `slice`) in the output.
// No-flip keeps the chunky row-major (z, y, x, c) layout the strip path writes;
// flip writes the column-major MATLAB (y, x, c, z) layout, matching
// flipTransposeStrip and flipTransposeStripRGB. Bit-exact copies, so the
// sample format does not matter.
template <typename T>
static void placeScanlineFlipped(T* out, const T* line, uint64_t x, uint64_t y,
                                 uint64_t samples, uint64_t row, uint64_t slice){
    T* base = out + slice*x*y*samples + row;
    for(uint64_t j = 0; j < x; j++)
        for(uint64_t c = 0; c < samples; c++)
            base[j*y + c*x*y] = line[j*samples + c];
}

static void placeScanline(void* out, const void* line, uint64_t x, uint64_t y, uint64_t samples,
                          uint64_t bits, uint64_t row, uint64_t slice, uint8_t flipXY){
    const uint64_t rowBytes = x*samples*(bits/8);
    if(!flipXY){
        memcpy(static_cast<uint8_t*>(out) + (slice*y + row)*rowBytes, line, rowBytes);
        return;
    }
    switch(bits){
        case 8:  placeScanlineFlipped(static_cast<uint8_t*>(out),  static_cast<const uint8_t*>(line),  x, y, samples, row, slice); break;
        case 16: placeScanlineFlipped(static_cast<uint16_t*>(out), static_cast<const uint16_t*>(line), x, y, samples, row, slice); break;
        case 32: placeScanlineFlipped(static_cast<uint32_t*>(out), static_cast<const uint32_t*>(line), x, y, samples, row, slice); break;
        case 64: placeScanlineFlipped(static_cast<uint64_t*>(out), static_cast<const uint64_t*>(line), x, y, samples, row, slice); break;
    }
}

// Place one decoded scanline of color p of a planar file (one sample per pixel).
template <typename T>
static void placePlanarScanlineT(T* out, const T* line, uint64_t x, uint64_t y, uint64_t samples, uint64_t p,
                                 uint64_t row, uint64_t slice, uint8_t flipXY){
    if(flipXY){
        T* base = out + (slice*samples + p)*x*y + row;
        for(uint64_t j = 0; j < x; j++) base[j*y] = line[j];
        return;
    }
    T* base = out + (slice*y + row)*x*samples + p;
    for(uint64_t j = 0; j < x; j++) base[j*samples] = line[j];
}

static void placePlanarScanline(void* out, const void* line, uint64_t x, uint64_t y, uint64_t samples, uint64_t p,
                                uint64_t bits, uint64_t row, uint64_t slice, uint8_t flipXY){
    switch(bits){
        case 8:  placePlanarScanlineT(static_cast<uint8_t*>(out),  static_cast<const uint8_t*>(line),  x, y, samples, p, row, slice, flipXY); break;
        case 16: placePlanarScanlineT(static_cast<uint16_t*>(out), static_cast<const uint16_t*>(line), x, y, samples, p, row, slice, flipXY); break;
        case 32: placePlanarScanlineT(static_cast<uint32_t*>(out), static_cast<const uint32_t*>(line), x, y, samples, p, row, slice, flipXY); break;
        case 64: placePlanarScanlineT(static_cast<uint64_t*>(out), static_cast<const uint64_t*>(line), x, y, samples, p, row, slice, flipXY); break;
    }
}

// Colors stored as separate planes, each read on its own (1 if chunky).
static uint16_t planesOf(TIFF* tif){
    uint16_t spp = 1, planarConfig = PLANARCONFIG_CONTIG;
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG, &planarConfig);
    return (spp > 1 && planarConfig == PLANARCONFIG_SEPARATE) ? spp : 1;
}

static bool setDirectoryWithRetry(TIFF* tif, uint64_t dir){
    for(int attempt = 0; attempt < 3; attempt++)
        if(TIFFSetDirectory(tif, (tdir_t)dir)) return true;
    return false;
}

// Backup method in case there are errors reading strips: re-read every slice a
// scanline at a time, each worker taking a block of whole slices. Reports an
// error if a row still cannot be read rather than leaving it unwritten.
uint8_t readTiffParallelBak(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint8_t flipXY){
    int32_t numWorkers = omp_get_max_threads();
    int32_t batchSize = (z-1)/numWorkers+1;
    uint64_t bytes = bits/8;

    int32_t w;
    std::atomic<uint8_t> err{0};
    char errString[10000];
    #pragma omp parallel for
    for(w = 0; w < numWorkers; w++){
        if(err) continue;
        const int64_t first = (int64_t)startSlice + (int64_t)w*batchSize;
        const int64_t last = std::min<int64_t>(first + batchSize, (int64_t)(startSlice + z));
        if(first >= last) continue;

        TiffPtr tif(openReaderTiff(fileName));
        if(!tif){
            #pragma omp critical
            {
                err = 1;
                sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
            }
            continue;
        }
        uint16_t spp = 1;
        TIFFGetFieldDefaulted(tif.get(), TIFFTAG_SAMPLESPERPIXEL, &spp);
        const uint16_t planes = planesOf(tif.get());
        std::unique_ptr<uint8_t[]> line(new uint8_t[x*spp*bytes]);
        for(int64_t dir = first; dir < last && !err; dir++){
            if(!setDirectoryWithRetry(tif.get(), (uint64_t)dir)){
                #pragma omp critical
                {
                    err = 1;
                    sprintf(errString,"Thread %d: TIFFSetDirectory(%ld) failed\n",w,(long)dir);
                }
                break;
            }
            for(uint16_t p = 0; p < planes && !err; p++)
            for(uint64_t i = 0; i < y; i++){
                if(TIFFReadScanline(tif.get(), line.get(), (uint32_t)i, p) < 0){
                    #pragma omp critical
                    {
                        err = 1;
                        sprintf(errString,"Thread %d: Row %lu of directory %ld cannot be read\n",w,(unsigned long)i,(long)dir);
                    }
                    break;
                }
                if(planes > 1) placePlanarScanline(tiff, line.get(), x, y, spp, p, bits, i, (uint64_t)(dir - (int64_t)startSlice), flipXY);
                else placeScanline(tiff, line.get(), x, y, spp, bits, i, (uint64_t)(dir - (int64_t)startSlice), flipXY);
            }
        }
    }
    if(err){
        printf("%s", errString);
    }
    return err;
}

uint8_t readTiffParallel(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint8_t flipXY, uint16_t compression){
    uint64_t bytes = bits/8;
    // More rows per strip than the image has (2^32 - 1, the TIFF default, is
    // common) is one strip. The buffers below were sized from it, so a flip
    // read of such a compressed image asked for terabytes and crashed.
    stripSize = std::max<uint64_t>(1, std::min(stripSize, y));
    // Pieces per slice. Each slice was read by one thread, so a stack of a few
    // more slices than threads left most of them idle at the end: on 48
    // threads, 49 slices of 2048 x 2048 LZW took 1.7x as long per slice as 48,
    // and 24 slices left half the threads idle. Such compressed stacks have
    // their slices cut into pieces of whole strips (or tiles), handed out like
    // slices: 49 slices then read 1.4-1.5x faster and 24 slices 1.3-1.6x.
    // A piece costs its thread a directory read, and in the flip layout the
    // pieces of a slice all write to its pages: 24 slices on 24 threads, cut
    // up for no gain in balance, read 10% slower. So slices are only cut when
    // that is predicted to read the stack at least 1.2x faster. Uncompressed
    // stacks keep their schedule below.
    const uint64_t threads = (uint64_t)std::max(1, omp_get_max_threads());
    const uint64_t kPiecesPerThread = 8, kMinSliceBytes = (uint64_t)1 << 20;
    uint64_t parts = 1;
    if(z > 1 && z < kPiecesPerThread*threads && x*y*bytes >= kMinSliceBytes && compression != COMPRESSION_NONE){
        uint64_t units = (y + stripSize - 1)/stripSize;
        bool tiledFile = false;
        // Callers that have the file open pass its compression; otherwise, and
        // for tiles, look at the first directory.
        TiffPtr probe;
        if(!compression || units == 1){
            probe.reset(openReaderTiff(fileName));
            if(probe){
                TIFFGetFieldDefaulted(probe.get(), TIFFTAG_COMPRESSION, &compression);
                tiledFile = TIFFIsTiled(probe.get());
            }
        }
        if(tiledFile){
            uint32_t tw = 0, tl = 0;
            TIFFGetField(probe.get(), TIFFTAG_TILEWIDTH, &tw);
            TIFFGetField(probe.get(), TIFFTAG_TILELENGTH, &tl);
            units = (tw && tl) ? ((x + tw - 1)/tw)*((y + tl - 1)/tl) : 1;
        }
        if(compression && compression != COMPRESSION_NONE){
            // Rounds of `threads` pieces with every slice cut into p.
            auto rounds = [&](uint64_t p){ return (double)((z*p + threads - 1)/threads); };
            const uint64_t pmax = std::min<uint64_t>(units, (kPiecesPerThread*threads + z - 1)/z);
            // Equal pieces (p divides the strips, or each piece has 8 or more,
            // so the remainder hardly matters) take rounds(p)/p slices' time.
            uint64_t best = 1;
            for(uint64_t p = 2; p <= pmax; p++)
                if((units % p == 0 || units >= 8*p) && rounds(p)/p < rounds(best)/best - 1e-9) best = p;
            if(rounds(1) + 1e-9 < 1.2*rounds(best)/best) best = 1;
            // Unequal pieces take as long as the largest. With equal pieces
            // only, slices of a prime number of strips stayed whole: cpp-tiff
            // writes 25 slices of 2048 x 2048 with 11 strips each, which 24
            // threads read 1.45-1.65x faster in 6 pieces. Where unequal pieces
            // only tied with equal ones, reads took up to 20% longer, so they
            // are used when predicted 1.2x faster.
            auto span = [&](uint64_t p){ return rounds(p)*(double)((units + p - 1)/p)/(double)units; };
            uint64_t uneven = best;
            for(uint64_t p = 2; p <= pmax; p++)
                if(span(p) < span(uneven) - 1e-9) uneven = p;
            if(span(best) + 1e-9 >= 1.2*span(uneven)) best = uneven;
            parts = best;
        }
    }
    // At most one thread per slice (or piece) reads, so a short stack does not
    // open the file once for every core. The team is still every thread, the
    // others only waiting at the end: libgomp ends the threads a smaller team
    // leaves out and starts them again for the next larger one, which made a
    // loop reading 4- and 12-slice stacks in turn 1.2-1.3x as slow on 24 and 48
    // cores.
    const int32_t team = (int32_t)std::max<uint64_t>(1, std::min<uint64_t>(threads, z*parts));
    const int32_t pool = (int32_t)std::max<uint64_t>(threads, (uint64_t)team);

    uint8_t errBak = 0;
    std::atomic<uint8_t> err{0};
    char errString[10000];
#ifdef __linux__
    // Blocks for the work-stealing schedule below: block b covers slices
    // [front, back) relative to startSlice, packed as front << 32 | back so one
    // compare-and-swap updates both ends.
    std::unique_ptr<std::atomic<uint64_t>[]> blocks(new std::atomic<uint64_t>[team]);
    for(int32_t b = 0; b < team; b++)
        blocks[b].store(((z*(uint64_t)b/team) << 32) | (z*(uint64_t)(b + 1)/team));
    bool uncompressedStack = false;   // set by one thread below, read by all
#endif
    // Each thread opens the file once and reads whichever slices it is handed.
    #pragma omp parallel num_threads(pool)
    {
        const int32_t w = omp_get_thread_num();
        const bool reader = w < team;
        uint8_t outCounter = 0;
        TiffPtr tif(reader ? openReaderTiff(fileName) : nullptr);   // RAII: TIFFClose on every exit path
        while(reader && !tif){
            tif.reset(openReaderTiff(fileName));
            if(outCounter == 3){
                #pragma omp critical
                {
                    err = 1;
                    sprintf(errString,"Thread %d: TIFFOpen failed for \"%s\"\n",w,fileName);
                }
                break;
            }
            outCounter++;
        }

        // Samples per pixel (1=grayscale, 3=RGB, 4=RGBA). xs is the interleaved
        // row width in samples; xs==x when grayscale, so those paths are unchanged.
        uint16_t spp = 1;
        if(tif) TIFFGetFieldDefaulted(tif.get(), TIFFTAG_SAMPLESPERPIXEL, &spp);
        const uint64_t samples = spp;
        const uint64_t xs = x * samples;
        // Planar color stores each color as its own image (see readPlanarStrips),
        // so a strip row holds one sample per pixel.
        const bool planar = tif && planesOf(tif.get()) > 1;
        const uint64_t stripRowBytes = (planar ? x : xs)*bytes;

        // Tiled TIFFs store sub-rectangles instead of full-width strips.
        const bool tiled = tif && TIFFIsTiled(tif.get());
        uint32_t tw = 0, tl = 0;
        uint64_t tileBufBytes = 0, tilesAcross = 0, tilesDown = 0;
        if(tiled){
            TIFFGetField(tif.get(), TIFFTAG_TILEWIDTH, &tw);
            TIFFGetField(tif.get(), TIFFTAG_TILELENGTH, &tl);
            const int64_t tileBytes = (int64_t)TIFFTileSize(tif.get());
            // Bail out rather than divide by zero below or read nothing at all.
            if(!tw || !tl || tileBytes <= 0){
                #pragma omp critical
                {
                    err = 1;
                    sprintf(errString,"Thread %d: tiled image has invalid tile geometry\n",w);
                }
            }
            else{
                tileBufBytes = (uint64_t)tileBytes;
                tilesAcross = (x + tw - 1) / tw;
                tilesDown = (y + tl - 1) / tl;
            }
        }

        // One scratch buffer per thread, reused for every slice it reads. Tiles
        // always decode into it; strips only on the flip path or when planar,
        // since the no-flip path decodes each chunky strip straight into the
        // output. The flip path decodes runs of strips into it and transposes
        // each run at once (see stripsPerTranspose); a planar no-flip read keeps
        // a run of every color in it to interleave. new[] (not vector or
        // make_unique) keeps it RAII-owned without zero-filling before the read.
        const bool buffered = flipXY || planar;
        const uint64_t bufferPlanes = (planar && !flipXY) ? samples : 1;
        const int64_t nStrips = (int64_t)((y + stripSize - 1) / stripSize);   // per color when planar
        const int64_t run = (!tiled && buffered) ? std::max<int64_t>(1, std::min<int64_t>(stripsPerTranspose(stripSize, stripRowBytes), nStrips)) : 1;
        // Uncompressed slices are read with one call per slice (no-flip) or per
        // ~1 MB block (flip, planar); see readRawStrips. A planar no-flip block
        // holds every color in that 1 MB: with 1 MB per color the buffer left
        // the L2 cache and those reads took 1.7x as long as chunky RGB. Sized
        // from the first directory.
        const bool rawCandidate = tif && !tiled && rawStripsOK(tif.get(), stripSize, stripRowBytes);
        const int64_t rawRun = buffered ? std::max<int64_t>(run, std::min<int64_t>(nStrips, (int64_t)std::max<uint64_t>(1, rawFlipBlockBytes / (stripSize*stripRowBytes*bufferPlanes))))
                                        : nStrips;
        const uint64_t blockBytes = (uint64_t)(rawCandidate ? rawRun : run)*stripSize*stripRowBytes;
        std::unique_ptr<uint8_t[]> buffer;
        if(reader && tiled && tileBufBytes) buffer.reset(new uint8_t[tileBufBytes]);
        else if(reader && !tiled && buffered) buffer.reset(new uint8_t[bufferPlanes*blockBytes]);

        // Strips (per color when planar) or tiles in a slice; readPart takes a
        // range of them.
        const int64_t units = tiled ? (int64_t)(tilesAcross*tilesDown) : nStrips;

        // Read strips (or tiles) [uBegin, uEnd) of one slice into their place
        // in the output.
        auto readPart = [&](int64_t dir, int64_t uBegin, int64_t uEnd){
            if(err) return;
            // (only readers are handed slices; one whose open failed set err)
            if(!tif){
                #pragma omp critical
                {
                    err = 1;
                    sprintf(errString,"Thread %d: handed a slice without a file handle\n",w);
                }
                return;
            }

            // A fresh handle is already on directory 0, and TIFFSetDirectory
            // would read it again.
            uint8_t counter = 0;
            while(TIFFCurrentDirectory(tif.get()) != (tdir_t)dir && !TIFFSetDirectory(tif.get(), (uint64_t)dir)){
                if(++counter == 3){
                    #pragma omp critical
                    {
                        err = 1;
                        sprintf(errString,"Thread %d: TIFFSetDirectory(%ld) failed\n",w,(long)dir);
                    }
                    return;
                }
            }

            if(tiled){
                int64_t failedTile = -1;
                if(readTilesInto(tif.get(), buffer.get(), (int64_t)tileBufBytes, tiff,
                                 (int64_t)x, (int64_t)y, (int64_t)bits, (int64_t)samples, planar,
                                 (int64_t)tw, (int64_t)tl, (int64_t)tilesAcross,
                                 uBegin, uEnd,
                                 dir - (int64_t)startSlice, flipXY, &failedTile)){
                    #pragma omp critical
                    {
                        // No errBak: the scanline fallback cannot read tiled images.
                        err = 1;
                        sprintf(errString,"Thread %d: Tile %ld cannot be read\n",w,(long)failedTile);
                    }
                }
                return;
            }

            const bool raw = rawCandidate && rawStripsOK(tif.get(), stripSize, stripRowBytes);
            const int64_t step = raw ? rawRun : run;
            if(planar){
                const int64_t failedStrip = readPlanarStrips(tif.get(), buffer.get(), blockBytes, tiff, x, y, bits, samples, stripSize,
                                                             uBegin, uEnd, step, raw, dir - (int64_t)startSlice, flipXY);
                if(failedStrip >= 0){
                    #pragma omp critical
                    {
                        errBak = 1;
                        err = 1;
                        sprintf(errString,"Thread %d: Strip %ld cannot be read\n",w,(long)failedStrip);
                    }
                }
                return;
            }
            for(int64_t first = uBegin; first < uEnd; first += step)
            {
                const int64_t last = std::min<int64_t>(first + step, uEnd);
                uint8_t* blockDst = flipXY ? buffer.get()
                    : static_cast<uint8_t*>(tiff) + (((first*stripSize*xs)+((dir-startSlice)*(y*xs)))*bytes);
                if(!(raw && readRawStrips(tif.get(), first, last, 0, stripSize, y, xs*bytes, bits, blockDst)))
                for(int64_t i = first; i < last; i++){
                    // No-flip: decode the strip straight into the output array (no
                    // bounce buffer, no full-image memcpy). Flip: decode into the
                    // buffer and transpose below.
                    void* stripDst = flipXY
                        ? static_cast<void*>(buffer.get() + (i - first)*stripSize*xs*bytes)
                        : static_cast<void*>(static_cast<uint8_t*>(tiff) + (((i*stripSize*xs)+((dir-startSlice)*(y*xs)))*bytes));
                    int64_t cBytes = TIFFReadEncodedStrip(tif.get(), i, stripDst, stripSize*xs*bytes);
                    if(cBytes < 0){
                        #pragma omp critical
                        {
                            errBak = 1;
                            err = 1;
                            sprintf(errString,"Thread %d: Strip %ld cannot be read\n",w,(long)i);
                        }
                        return;
                    }
                }
                if(!flipXY){
                    continue;
                }
                // Flip x/y into the column-major MATLAB output with a cache-
                // blocked transpose (see flipTransposeStrip). stripRows handles
                // a partial final strip of the image.
                const int64_t rowOff = first*stripSize;
                const int64_t stripRows = std::min<int64_t>((last - first)*(int64_t)stripSize, (int64_t)y - rowOff);
                int64_t sliceOff = (dir-startSlice)*(x*y);
                if(samples == 1){
                    switch(bits){
                        case 8:
                            flipTransposeStrip<uint8_t>(static_cast<uint8_t*>(tiff), buffer.get(), x,y,stripRows,rowOff,sliceOff);
                            break;
                        case 16:
                            flipTransposeStrip<uint16_t>(static_cast<uint16_t*>(tiff), reinterpret_cast<uint16_t*>(buffer.get()), x,y,stripRows,rowOff,sliceOff);
                            break;
                        case 32:
                            flipTransposeStrip<float>(static_cast<float*>(tiff), reinterpret_cast<float*>(buffer.get()), x,y,stripRows,rowOff,sliceOff);
                            break;
                        case 64:
                            flipTransposeStrip<double>(static_cast<double*>(tiff), reinterpret_cast<double*>(buffer.get()), x,y,stripRows,rowOff,sliceOff);
                            break;
                    }
                } else {
                    int64_t sliceOffC = (dir-startSlice)*(x*y*samples);
                    switch(bits){
                        case 8:
                            flipTransposeStripRGB<uint8_t>(static_cast<uint8_t*>(tiff), buffer.get(), x,y,stripRows,rowOff,sliceOffC,samples);
                            break;
                        case 16:
                            flipTransposeStripRGB<uint16_t>(static_cast<uint16_t*>(tiff), reinterpret_cast<uint16_t*>(buffer.get()), x,y,stripRows,rowOff,sliceOffC,samples);
                            break;
                        case 32:
                            flipTransposeStripRGB<float>(static_cast<float*>(tiff), reinterpret_cast<float*>(buffer.get()), x,y,stripRows,rowOff,sliceOffC,samples);
                            break;
                        case 64:
                            flipTransposeStripRGB<double>(static_cast<double*>(tiff), reinterpret_cast<double*>(buffer.get()), x,y,stripRows,rowOff,sliceOffC,samples);
                            break;
                    }
                }
            }
        };
        auto readOne = [&](int64_t dir){ readPart(dir, 0, units); };

        // Pieces of slices (see parts above), handed out one at a time on Linux
        // like whole compressed slices below, and in contiguous blocks
        // elsewhere: the first piece of every slice, then the second, and so
        // on. In the flip layout each piece of a slice writes to all of the
        // slice's output pages, so this way the pieces read at the same time
        // are of different slices rather than all faulting in the same pages.
        if(parts > 1){
            const int64_t items = (int64_t)(z*parts), perSlice = (int64_t)parts, slices = (int64_t)z;
            auto readItem = [&](int64_t k){
                const int64_t piece = k / slices;
                readPart((int64_t)startSlice + k % slices, units*piece/perSlice, units*(piece + 1)/perSlice);
            };
            // Fewer pieces than threads: one each, to the readers.
            if(team < pool){
                #pragma omp for schedule(static, 1)
                for(int64_t k = 0; k < items; k++) readItem(k);
            }
            else{
#ifdef __linux__
                #pragma omp for schedule(dynamic, 1)
#else
                #pragma omp for schedule(static)
#endif
                for(int64_t k = 0; k < items; k++) readItem(k);
            }
        }
        else{
            // On Linux, once there are more slices than threads:
            // - Compressed stacks hand slices out one at a time, so a thread that
            //   decodes slowly does not hold up the whole read. With the unmapped
            //   reads (openReaderTiff) that cut read time by 21-34% on our
            //   benchmarks.
            // - Uncompressed slices take about the same time each, so each thread
            //   reads its own contiguous block front to back and then takes slices
            //   from the back of other blocks. One at a time made every thread walk
            //   the whole IFD chain (libtiff can only step forward to a directory it
            //   has not seen); on NFS, where every step is a read call, that alone
            //   kept 1000-slice uncompressed stacks 6-10% behind v1.6.2, while
            //   contiguous blocks matched or beat it.
            // With no more slices than threads each thread simply takes one.
            // Elsewhere each thread keeps one contiguous block, as before: with a
            // memory-mapped file the balancing was no faster on macOS.
#ifdef __linux__
            #pragma omp single
            {
                uint16_t comp = COMPRESSION_NONE;
                if(tif) TIFFGetFieldDefaulted(tif.get(), TIFFTAG_COMPRESSION, &comp);
                uncompressedStack = tif && comp == COMPRESSION_NONE;
            }
            if(z > (uint64_t)team && uncompressedStack){
                // Take one slice from the front (owner) or the back (thief) of block b.
                auto take = [&](int32_t b, bool front, uint64_t& slice){
                    uint64_t cur = blocks[b].load();
                    for(;;){
                        const uint64_t f = cur >> 32, e = cur & 0xFFFFFFFFu;
                        if(f >= e) return false;
                        const uint64_t next = front ? ((f + 1) << 32) | e : (f << 32) | (e - 1);
                        if(blocks[b].compare_exchange_weak(cur, next)){ slice = front ? f : e - 1; return true; }
                    }
                };
                // libtiff remembers the offset of every IFD a handle has stepped
                // over, so jumping back to a slice at or below the furthest one this
                // thread has read is free, while going past it walks the chain.
                // Thieves therefore pick the block whose last slice is cheapest to
                // reach, and the fullest block among equally cheap ones.
                uint64_t slice, reached = 0;
                while(!err && take(w, true, slice)){ readOne((int64_t)(startSlice + slice)); reached = std::max(reached, slice); }
                while(!err){
                    int32_t victim = -1; uint64_t bestCost = UINT64_MAX, bestLeft = 0;
                    for(int32_t b = 0; b < team; b++){
                        const uint64_t cur = blocks[b].load(), f = cur >> 32, e = cur & 0xFFFFFFFFu;
                        if(e <= f) continue;
                        const uint64_t cost = e - 1 <= reached ? 0 : e - 1 - reached;
                        if(cost < bestCost || (cost == bestCost && e - f > bestLeft)){ bestCost = cost; bestLeft = e - f; victim = b; }
                    }
                    if(victim < 0) break;
                    if(take(victim, false, slice)){ readOne((int64_t)(startSlice + slice)); reached = std::max(reached, slice); }
                }
            }
            else if(z > (uint64_t)team){
                #pragma omp for schedule(dynamic, 1)
                for(int64_t dir = (int64_t)startSlice; dir < (int64_t)(startSlice+z); dir++) readOne(dir);
            }
            else
#endif
            {
                // (contiguous blocks; with no more slices than threads, slice k
                // goes to thread k, a reader)
                #pragma omp for schedule(static)
                for(int64_t dir = (int64_t)startSlice; dir < (int64_t)(startSlice+z); dir++) readOne(dir);
            }
        }
    }
    if(err){
        if(errBak) return readTiffParallelBak(x, y, z, fileName, tiff, bits, startSlice, flipXY);
        else {
			printf("%s", errString);
		}
    }
	return err;
}

// Backup method in case there are errors reading strips of a single image:
// re-read it a scanline at a time. Compressed strips can only be decoded from
// their first row, so workers split the image by whole strips and read each
// one in row order. Reports an error if a row still cannot be read rather than
// leaving it unwritten.
uint8_t readTiffParallel2DBak(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint8_t flipXY){
    (void)z;   // a single image
    uint64_t bytes = bits/8;
    uint32_t rowsPerStrip = 0;
    {
        TiffPtr probe(openReaderTiff(fileName));
        if(probe && setDirectoryWithRetry(probe.get(), startSlice))
            TIFFGetFieldDefaulted(probe.get(), TIFFTAG_ROWSPERSTRIP, &rowsPerStrip);
    }
    if(!rowsPerStrip || rowsPerStrip > y) rowsPerStrip = (uint32_t)y;
    const uint64_t numStrips = (y + rowsPerStrip - 1) / rowsPerStrip;
    const int32_t numWorkers = (int32_t)std::max<uint64_t>(1, std::min<uint64_t>((uint64_t)omp_get_max_threads(), numStrips));
    const uint64_t stripsPerWorker = (numStrips + numWorkers - 1) / numWorkers;

    int32_t w;
    std::atomic<uint8_t> err{0};
    char errString[10000];
    // Every thread of the calling thread's team, one worker each (see
    // readTiffParallel2D).
    #pragma omp parallel for num_threads(omp_get_max_threads()) schedule(static, 1)
    for(w = 0; w < numWorkers; w++){
        if(err) continue;
        const uint64_t rowBegin = std::min<uint64_t>((uint64_t)w*stripsPerWorker*rowsPerStrip, y);
        const uint64_t rowEnd = std::min<uint64_t>(((uint64_t)w + 1)*stripsPerWorker*rowsPerStrip, y);
        if(rowBegin >= rowEnd) continue;

        TiffPtr tif(openReaderTiff(fileName));
        if(!tif){
            #pragma omp critical
            {
                err = 1;
                sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
            }
            continue;
        }
        if(!setDirectoryWithRetry(tif.get(), startSlice)){
            #pragma omp critical
            {
                err = 1;
                sprintf(errString,"Thread %d: TIFFSetDirectory(%lu) failed\n",w,(unsigned long)startSlice);
            }
            continue;
        }
        uint16_t spp = 1;
        TIFFGetFieldDefaulted(tif.get(), TIFFTAG_SAMPLESPERPIXEL, &spp);
        const uint16_t planes = planesOf(tif.get());
        std::unique_ptr<uint8_t[]> line(new uint8_t[x*spp*bytes]);
        for(uint16_t p = 0; p < planes && !err; p++)
        for(uint64_t i = rowBegin; i < rowEnd; i++){
            if(TIFFReadScanline(tif.get(), line.get(), (uint32_t)i, p) < 0){
                #pragma omp critical
                {
                    err = 1;
                    sprintf(errString,"Thread %d: Row %lu cannot be read\n",w,(unsigned long)i);
                }
                break;
            }
            if(planes > 1) placePlanarScanline(tiff, line.get(), x, y, spp, p, bits, i, 0, flipXY);
            else placeScanline(tiff, line.get(), x, y, spp, bits, i, 0, flipXY);
        }
    }
    if(err){
        printf("%s", errString);
    }
    return err;
}


uint8_t readTiffParallel2D(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint8_t flipXY){
    // More rows per strip than the image has is one strip (see readTiffParallel).
    stripSize = std::max<uint64_t>(1, std::min(stripSize, y));
    int32_t numWorkers = omp_get_max_threads();
    uint64_t stripsPerDir = (uint64_t)ceil((double)y/(double)stripSize);
    int32_t batchSize = (stripsPerDir-1)/numWorkers+1;
    uint64_t bytes = bits/8;

    int32_t w;
    std::atomic<uint8_t> err{0};
    uint8_t errBak = 0;
    char errString[10000];
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif){
        printf("File \"%s\" cannot be opened\n",fileName);
        return 1;
    }
    uint16_t spp = 1;
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    const uint64_t samples = spp;          // 1=grayscale, 3=RGB, 4=RGBA
    const uint64_t xs = x * samples;       // interleaved row width (xs==x when grayscale)
    // Planar color: each color is its own image, one sample per pixel per row,
    // read through the buffer in both layouts (see readPlanarStrips).
    const bool planar = planesOf(tif) > 1;
    const uint64_t stripRowBytes = (planar ? x : xs)*bytes;
    const bool buffered = flipXY || planar;

    // Tiled TIFFs store sub-rectangles instead of full-width strips. A single 2D
    // image is one directory, so here the workers are split across tiles.
    const bool tiled = TIFFIsTiled(tif);
    uint32_t tw = 0, tl = 0;
    uint64_t tileBufBytes = 0, tilesAcross = 0, tilesDown = 0, numTiles = 0;
    if(tiled){
        TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tw);
        TIFFGetField(tif, TIFFTAG_TILELENGTH, &tl);
        const int64_t tileBytes = (int64_t)TIFFTileSize(tif);
        // Bail out rather than divide by zero below or read nothing at all.
        if(!tw || !tl || tileBytes <= 0){
            printf("Tiled image \"%s\" has invalid tile geometry\n",fileName);
            TIFFClose(tif);
            return 1;
        }
        tileBufBytes = (uint64_t)tileBytes;
        tilesAcross = (x + tw - 1) / tw;
        tilesDown = (y + tl - 1) / tl;
        numTiles = tilesAcross * tilesDown;
    }
    const int32_t tileBatch = tiled ? (int32_t)((numTiles-1)/numWorkers+1) : 0;

    // Only as many workers as there are strips or tiles to read: a small image
    // has just a few, and every extra worker opened the file and selected its
    // directory for nothing. The team is still every thread (one worker each,
    // the rest wait at the end): libgomp ends the threads a smaller team leaves
    // out and starts them again for the next larger one, which made a loop
    // reading two image sizes in turn twice as slow on 48 cores.
    const uint64_t workItems = tiled ? numTiles : stripsPerDir;
    const uint64_t perWorker = tiled ? (uint64_t)tileBatch : (uint64_t)batchSize;
    const int32_t activeWorkers = (int32_t)std::max<uint64_t>(1, std::min<uint64_t>((uint64_t)numWorkers, (workItems + perWorker - 1) / perWorker));
    #pragma omp parallel for num_threads(numWorkers) schedule(static, 1)
    for(w = 0; w < activeWorkers; w++){

        uint8_t outCounter = 0;
        TIFF* tif = openReaderTiff(fileName);
        while(!tif){
            tif = openReaderTiff(fileName);
            if(outCounter == 3){
                #pragma omp critical
                {
                    err = 1;
                    sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
                }
                break;
            }
            outCounter++;
        }

        if(!tif) continue;   // err is already set; never hand libtiff a null handle
        // This worker's strips. The flip path decodes runs of them into the
        // buffer and transposes each run at once (see stripsPerTranspose).
        const int64_t sBegin = (int64_t)w*batchSize;
        const int64_t sEnd = std::min<int64_t>(sBegin + batchSize, (int64_t)stripsPerDir);
        const int64_t run = (!tiled && buffered) ? std::max<int64_t>(1, std::min<int64_t>(stripsPerTranspose(stripSize, stripRowBytes), sEnd - sBegin)) : 1;
        // A fresh handle is already on directory 0, and TIFFSetDirectory would
        // read it again.
        uint8_t counter = 0;
        while(TIFFCurrentDirectory(tif) != (tdir_t)startSlice && !TIFFSetDirectory(tif, startSlice) && counter<3){
            printf("Thread %d: File \"%s\" Directory \"%d\" failed to open. Try %d\n",w,fileName,0,counter+1);
            counter++;
            if(counter == 3){
                #pragma omp critical
                {
                    err = 1;
                    sprintf(errString,"Thread %d: File \"%s\" cannot be opened\n",w,fileName);
                }
            }
        }
        // An uncompressed image reads this worker's strips with one call
        // (no-flip) or per ~1 MB block (flip, planar); see readRawStrips.
        const bool raw = !tiled && rawStripsOK(tif, stripSize, stripRowBytes);
        const uint64_t bufferPlanes = (planar && !flipXY) ? samples : 1;
        const int64_t step = !raw ? run
            : buffered ? std::max<int64_t>(run, std::min<int64_t>(sEnd - sBegin, (int64_t)std::max<uint64_t>(1, rawFlipBlockBytes / (stripSize*stripRowBytes*bufferPlanes))))
                       : std::max<int64_t>(1, sEnd - sBegin);
        // Tiles always decode into the buffer; chunky strips only on the flip
        // path, since the no-flip path decodes each strip straight into the
        // output. A planar no-flip read keeps a block of every color.
        const uint64_t blockBytes = (uint64_t)step*stripSize*stripRowBytes;
        void* buffer = (tiled || buffered) ? malloc(tiled ? tileBufBytes : bufferPlanes*blockBytes) : NULL;
        if(tiled){
            int64_t tStart = (int64_t)w*tileBatch;
            int64_t tEnd = tStart + tileBatch;
            if(tEnd > (int64_t)numTiles) tEnd = (int64_t)numTiles;
            int64_t failedTile = -1;
            if(!err && tStart < tEnd &&
               readTilesInto(tif, (uint8_t*)buffer, (int64_t)tileBufBytes, tiff,
                             (int64_t)x, (int64_t)y, (int64_t)bits, (int64_t)samples, planar,
                             (int64_t)tw, (int64_t)tl, (int64_t)tilesAcross,
                             tStart, tEnd, 0, flipXY, &failedTile)){
                #pragma omp critical
                {
                    // No errBak: the scanline fallback cannot read tiled images.
                    err = 1;
                    sprintf(errString,"Thread %d: Tile %ld cannot be read\n",w,(long)failedTile);
                }
            }
        }
        else if(planar){
            const int64_t failedStrip = (!err && sBegin < sEnd)
                ? readPlanarStrips(tif, static_cast<uint8_t*>(buffer), blockBytes, tiff, x, y, bits, samples, stripSize, sBegin, sEnd, step, raw, 0, flipXY)
                : -1;
            if(failedStrip >= 0){
                #pragma omp critical
                {
                    errBak = 1;
                    err = 1;
                    sprintf(errString,"Thread %d: Strip %ld cannot be read\n",w,(long)failedStrip);
                }
            }
        }
        else for(int64_t first = sBegin; first < sEnd && !err; first += step){
            const int64_t last = std::min<int64_t>(first + step, sEnd);
            uint8_t* blockDst = flipXY ? static_cast<uint8_t*>(buffer) : static_cast<uint8_t*>(tiff) + first*stripSize*xs*bytes;
            if(!(raw && readRawStrips(tif, first, last, 0, stripSize, y, xs*bytes, bits, blockDst)))
            for(int64_t i = first; i < last; i++){
                // No-flip: decode the strip straight into the output, as the 3D
                // path does. Flip: decode into the buffer and transpose below.
                void* stripDst = flipXY
                    ? static_cast<void*>(static_cast<uint8_t*>(buffer) + (i - first)*stripSize*xs*bytes)
                    : static_cast<void*>(static_cast<uint8_t*>(tiff) + i*stripSize*xs*bytes);
                int64_t cBytes = TIFFReadEncodedStrip(tif, i, stripDst, stripSize*xs*bytes);
                if(cBytes < 0){
                    #pragma omp critical
                    {
                        errBak = 1;
                        err = 1;
                        sprintf(errString,"Thread %d: Strip %ld cannot be read\n",w,(long)i);
                    }
                    break;
                }
            }
            if(!flipXY || err) continue;
            // Same cache-blocked transpose as the 3D path, single image so
            // the slice offset is 0. samples==1 keeps the grayscale path.
            const int64_t rowOff = first*stripSize;
            const int64_t stripRows = std::min<int64_t>((last - first)*(int64_t)stripSize, (int64_t)y - rowOff);
            if(samples == 1){
                switch(bits){
                    case 8:
                        flipTransposeStrip<uint8_t>((uint8_t*)tiff,(uint8_t*)buffer,x,y,stripRows,rowOff,0);
                        break;
                    case 16:
                        flipTransposeStrip<uint16_t>((uint16_t*)tiff,(uint16_t*)buffer,x,y,stripRows,rowOff,0);
                        break;
                    case 32:
                        flipTransposeStrip<float>((float*)tiff,(float*)buffer,x,y,stripRows,rowOff,0);
                        break;
                    case 64:
                        flipTransposeStrip<double>((double*)tiff,(double*)buffer,x,y,stripRows,rowOff,0);
                        break;
                }
            } else {
                switch(bits){
                    case 8:
                        flipTransposeStripRGB<uint8_t>((uint8_t*)tiff,(uint8_t*)buffer,x,y,stripRows,rowOff,0,samples);
                        break;
                    case 16:
                        flipTransposeStripRGB<uint16_t>((uint16_t*)tiff,(uint16_t*)buffer,x,y,stripRows,rowOff,0,samples);
                        break;
                    case 32:
                        flipTransposeStripRGB<float>((float*)tiff,(float*)buffer,x,y,stripRows,rowOff,0,samples);
                        break;
                    case 64:
                        flipTransposeStripRGB<double>((double*)tiff,(double*)buffer,x,y,stripRows,rowOff,0,samples);
                        break;
                }
            }
        }
        free(buffer);
        TIFFClose(tif);
    }
    // Close the outer handle as well. It was only needed for the tag reads
    // above (each worker opens its own), and leaving it open leaked a handle
    // per 2D read -- harmless on Linux but it keeps the file locked on Windows.
    TIFFClose(tif);

    if(err) {
        if(errBak) return readTiffParallel2DBak(x, y, z, fileName, tiff, bits, startSlice, flipXY);
        else printf("%s", errString);
    }
	return err;
}


// Seek to an absolute byte offset; 0 on success. fseek takes a long, which is
// 32 bits on Windows, so offsets past 2 GB need the 64-bit variants.
static int seekFile(FILE* fp, uint64_t offset){
#ifdef _WIN32
    return _fseeki64(fp, (long long)offset, SEEK_SET);
#else
    return fseeko(fp, (off_t)offset, SEEK_SET);
#endif
}

// Reading images saved by ImageJ
uint8_t readTiffParallelImageJ(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint8_t flipXY){
    std::atomic<uint8_t> err{0};
    FILE *fp = fopen(fileName, "rb");
    if(!fp){ 
		printf("File \"%s\" cannot be opened from Disk\n",fileName);
		err = 1;
		return err;
	}
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif){
		printf("File \"%s\" cannot be opened\n",fileName);
		fclose(fp);
		err = 1;
		return err;
	}
    uint64_t* offsets = NULL;
    TIFFGetField(tif, TIFFTAG_STRIPOFFSETS, &offsets);
    if(!offsets){
        printf("Could not get offsets from the tiff file\n");
        TIFFClose(tif);
        fclose(fp);
        err = 1;
        return err;
    }
    const uint64_t offset = offsets[0];

    // ImageJ RGB stacks are chunky interleaved (spp=3), so every size below
    // scales by samples. Grayscale (samples==1) is unchanged.
    uint16_t spp = 1;
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    const uint64_t samples = spp;

    // The raw fread below bypasses libtiff, so swap bytes only when the file's
    // byte order differs from the host's (ImageJ itself writes big-endian, but
    // other programs write ImageJ-format files little-endian).
    const int needSwap = TIFFIsByteSwapped(tif);

    TIFFClose(tif);

    uint64_t bytes = bits/8;
    // The data is one contiguous block, so a z-range is just an offset into it.
    if(seekFile(fp, offset + startSlice*x*y*samples*bytes)){
        printf("File \"%s\": cannot seek to the image data\n",fileName);
        fclose(fp);
        err = 1;
        return err;
    }
    //#pragma omp parallel for
    /*
    for(uint64_t i = 0; i < z; i++){
    uint64_t cOffset = x*y*bytes*i;
    //pread(fd,tiff+cOffset,x*y*bytes,offset+cOffset);
    read(fd,tiff+cOffset,x*y*bytes);
    }*/
    // fread at most INT_MAX bytes at a time. A short read means the file ends
    // before its image data does (e.g. an interrupted copy), so report that
    // instead of retrying forever or returning the unread part of the buffer.
    const uint64_t tBytes = x*y*z*samples*bytes;
    uint64_t chunk = 0;
    while(chunk < tBytes){
        const uint64_t want = std::min<uint64_t>(tBytes - chunk, INT_MAX);
        const uint64_t got = fread(static_cast<uint8_t*>(tiff)+chunk,1,want,fp);
        chunk += got;
        if(got < want) break;
    }
    fclose(fp);
    if(chunk < tBytes){
        printf("File \"%s\" ends before its image data does\n",fileName);
        err = 1;
        return err;
    }
    // Swap endianess for types greater than 8 bits when the file's byte order
    // differs from the host's.
    if(bits > 8 && needSwap){
        #pragma omp parallel for
        for(uint64_t i = 0; i < x*y*z*samples; i++){
            switch(bits){
                case 16:
                    //((uint16_t*)tiff)[i] = ((((uint16_t*)tiff)[i] & 0xff) >> 8) | (((uint16_t*)tiff)[i] << 8);
                    //((uint16_t*)tiff)[i] = bswap_16(((uint16_t*)tiff)[i]);
                    ((uint16_t*)tiff)[i] = ((((uint16_t*)tiff)[i] << 8) & 0xff00) | ((((uint16_t*)tiff)[i] >> 8) & 0x00ff);
                    break;
                case 32:
                    //((num & 0xff000000) >> 24) | ((num & 0x00ff0000) >> 8) | ((num & 0x0000ff00) << 8) | (num << 24)
                    //((float*)tiff)[i] = bswap_32(((float*)tiff)[i]);
                    ((uint32_t*)tiff)[i] = ((((uint32_t*)tiff)[i] << 24) & 0xff000000 ) |
                        ((((uint32_t*)tiff)[i] <<  8) & 0x00ff0000 ) |
                        ((((uint32_t*)tiff)[i] >>  8) & 0x0000ff00 ) |
                        ((((uint32_t*)tiff)[i] >> 24) & 0x000000ff );
                    break;
                case 64:
                    //((double*)tiff)[i] = bswap_64(((double*)tiff)[i]);
                    ((uint64_t*)tiff)[i] = ( (((uint64_t*)tiff)[i] << 56) & 0xff00000000000000UL ) |
                        ( (((uint64_t*)tiff)[i] << 40) & 0x00ff000000000000UL ) |
                        ( (((uint64_t*)tiff)[i] << 24) & 0x0000ff0000000000UL ) |
                        ( (((uint64_t*)tiff)[i] <<  8) & 0x000000ff00000000UL ) |
                        ( (((uint64_t*)tiff)[i] >>  8) & 0x00000000ff000000UL ) |
                        ( (((uint64_t*)tiff)[i] >> 24) & 0x0000000000ff0000UL ) |
                        ( (((uint64_t*)tiff)[i] >> 40) & 0x000000000000ff00UL ) |
                        ( (((uint64_t*)tiff)[i] >> 56) & 0x00000000000000ffUL );
                    break;
            }

        }
    }
    // Find a way to do this in-place without making a copy
    if(flipXY){
        uint64_t size = x*y*z*samples*(bits/8);
        void* tiffC = malloc(size);
        memcpy(tiffC,tiff,size);
        if(samples == 1){
        #pragma omp parallel for
        for(uint64_t k = 0; k < z; k++){
            for(uint64_t j = 0; j < x; j++){
                for(uint64_t i = 0; i < y; i++){
                    switch(bits){
                        case 8:
                            ((uint8_t*)tiff)[i+(j*y)+(k*x*y)] = ((uint8_t*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                        case 16:
                            ((uint16_t*)tiff)[i+(j*y)+(k*x*y)] = ((uint16_t*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                        case 32:
                            ((float*)tiff)[i+(j*y)+(k*x*y)] = ((float*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                        case 64:
                            ((double*)tiff)[i+(j*y)+(k*x*y)] = ((double*)tiffC)[j+(i*x)+(k*x*y)];
                            break;
                    }
                }
            }
        }
        }
        else{
            // De-interleave each chunky slice into the column-major (y,x,c,z)
            // MATLAB layout, same mapping as the strip readers use.
            #pragma omp parallel for
            for(uint64_t k = 0; k < z; k++){
                const uint64_t off = k*x*y*samples;
                switch(bits){
                    case 8:
                        flipTransposeStripRGB<uint8_t>((uint8_t*)tiff, (uint8_t*)tiffC + off, x, y, y, 0, off, samples);
                        break;
                    case 16:
                        flipTransposeStripRGB<uint16_t>((uint16_t*)tiff, (uint16_t*)tiffC + off, x, y, y, 0, off, samples);
                        break;
                    case 32:
                        flipTransposeStripRGB<float>((float*)tiff, (float*)tiffC + off, x, y, y, 0, off, samples);
                        break;
                    case 64:
                        flipTransposeStripRGB<double>((double*)tiff, (double*)tiffC + off, x, y, y, 0, off, samples);
                        break;
                }
            }
        }
        free(tiffC);
    }
	return err;
}


// tiff pointer guaranteed to be NULL or the correct size array for the tiff file
void* readTiffParallelInfo(const char* fileName, const TiffInfo &info, void* tiff, uint8_t flipXY, const std::vector<uint64_t> &zRange)
{
	TIFFSetWarningHandler(DummyHandler);
	uint64_t x = info.x, y = info.y, z = info.z, bits = info.bits, startSlice = 0;
	// libtiff's own value, which is the TIFF default (2^32 - 1, one strip) when
	// the tag is missing; the readers limit it to the image height.
	const uint64_t stripSize = info.stripSize;
	const uint64_t samples = info.samples;   // 1=grayscale, 3=RGB, 4=RGBA (chunky or planar, both read channel-last)
	const bool tiled = info.tiled;
	const uint16_t compression = info.compression;
	// ImageJ images with ImageJ metadata: the ImageJ reader pulls the pixels off
	// disk as one contiguous block, which is not how tiled or planar data is
	// laid out, so those go through the normal readers instead.
	const uint8_t imageJIm = info.imageJ;

    if(zRange.size()){
        if(zRange.size() == 2){
            startSlice = zRange[0];
            z = zRange[1] - zRange[0];   // number of slices to read
        }
        else{
            startSlice = zRange[0];
            z = 1;                        // single slice
        }
    }


	// Bytes are moved purely by bit width, so one allocation covers every
	// sample format.
	if(bits != 8 && bits != 16 && bits != 32 && bits != 64) return NULL;
	const bool allocated = !tiff;
	if(allocated){
		tiff = malloc(x*y*z*samples*(bits/8));
		if(!tiff) return NULL;
	}
	// A stack of a few large slices with many strips each (two 4096 x 4096
	// planes, say) left most threads idle in the 3D reader, which gives each
	// slice one thread: LZW took 108 ms against 28 ms read slice by slice with
	// every thread. Each slice read that way has a fixed cost, since every
	// thread opens the file and walks to the slice, so this is only for large
	// slices (4 MB, or 16 MB uncompressed, where reading is only a copy) near
	// the start of the file.
	const uint64_t sliceStrips = (y + stripSize - 1) / std::max<uint64_t>(1, stripSize);
	const uint64_t sliceBytesAll = x*y*samples*(bits/8);
	const bool sliceBySlice = !imageJIm && !tiled && z > 1 && startSlice + z <= 64
	                          && sliceBytesAll >= (compression == COMPRESSION_NONE ? (uint64_t)16 << 20 : (uint64_t)4 << 20)
	                          && 2*z <= std::min<uint64_t>((uint64_t)omp_get_max_threads(), sliceStrips);
	uint8_t err = 0;
	if(imageJIm)    err = readTiffParallelImageJ(x,y,z,fileName, tiff, bits, startSlice, stripSize, flipXY);
	else if(z <= 1) err = readTiffParallel2D(x,y,z,fileName, tiff, bits, startSlice, stripSize, flipXY);
	else if(sliceBySlice){
		for(uint64_t s = 0; s < z && !err; s++)
			err = readTiffParallel2D(x,y,1,fileName, static_cast<uint8_t*>(tiff) + s*sliceBytesAll, bits, startSlice + s, stripSize, flipXY);
	}
	else            err = readTiffParallel(x,y,z,fileName, tiff, bits, startSlice, stripSize, flipXY, compression);
	// A failed read leaves part of the buffer unwritten, so report it instead of
	// returning partial data. The reader has already printed why.
	if(err){
		if(allocated) free(tiff);
		return NULL;
	}
	return tiff;
}

// One open and one walk of the directories for everything the read needs (it
// was three opens and a walk here).
void* readTiffParallelWrapperHelper(const char* fileName, void* tiff, uint8_t flipXY, const std::vector<uint64_t> &zRange = {})
{
	TiffInfo info;
	if(getTiffInfo(fileName, info)) return NULL;
	return readTiffParallelInfo(fileName, info, tiff, flipXY, zRange);
}

void* readTiffParallelWrapper(const char* fileName, const std::vector<uint64_t> &zRange)
{
	return readTiffParallelWrapperHelper(fileName,NULL,1,zRange);
}

void* readTiffParallelWrapperNoXYFlip(const char* fileName, const std::vector<uint64_t> &zRange)
{
	return readTiffParallelWrapperHelper(fileName,NULL,0,zRange);
}

// tTiff doesn't matter as tiff is set in the function
void readTiffParallelWrapperSet(const char* fileName, void* tiff){
	void* tTiff = readTiffParallelWrapperHelper(fileName,tiff,0);
}
