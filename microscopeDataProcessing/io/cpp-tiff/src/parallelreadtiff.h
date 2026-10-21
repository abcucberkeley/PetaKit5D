#ifndef READTIFFPARALLEL_H
#define READTIFFPARALLEL_H

#include <cstdint>
#include <vector>

// compression is the file's TIFF compression, or 0 if the caller does not know it.
uint8_t readTiffParallel(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint8_t flipXY, uint16_t compression = 0);

uint8_t readTiffParallel2D(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint8_t flipXY);

uint8_t readTiffParallelImageJ(uint64_t x, uint64_t y, uint64_t z, const char* fileName, void* tiff, uint64_t bits, uint64_t startSlice, uint64_t stripSize, uint8_t flipXY);

struct TiffInfo;

// Reads slices of a file whose getTiffInfo is info (all of them, or zRange:
// {first} or {first, end}) into tiff, or into a new buffer when tiff is NULL,
// row-major (z, y, x[, c]) when flipXY is 0. Returns the buffer, or NULL if
// the read fails.
void* readTiffParallelInfo(const char* fileName, const TiffInfo &info, void* tiff, uint8_t flipXY, const std::vector<uint64_t> &zRange = {});

void* readTiffParallelWrapper(const char* fileName, const std::vector<uint64_t> &zRange = {});

void* readTiffParallelWrapperNoXYFlip(const char* fileName, const std::vector<uint64_t> &zRange = {});

void readTiffParallelWrapperSet(const char* fileName, void* tiff);

#endif
