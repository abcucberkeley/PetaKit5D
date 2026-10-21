#ifndef WRITETIFFPARALLEL_H
#define WRITETIFFPARALLEL_H

#include <cstdint>
#include <string>

// samples is the samples per pixel: 1 = grayscale, 3 = RGB, 4 = RGBA. Color
// data keeps each pixel's samples together: row-major (z, y, x, c), or MATLAB
// column-major (y, x, c, z) when transpose is true.
uint8_t writeTiffSingle(const uint64_t x, const uint64_t y, const uint64_t z, const char* fileName, const void* tiff, const void* tiffOld, const uint64_t bits, const uint64_t startSlice, const uint64_t stripSize, const char* mode, const bool transpose, const std::string &compression, uint16_t sampleFormat = 0, uint64_t tileSize = 0, uint64_t samples = 1);

uint8_t writeTiffParallel(const uint64_t x, const uint64_t y, const uint64_t z, const char* fileName, const void* tiff, const void* tiffOld, const uint64_t bits, const uint64_t startSlice, const uint64_t stripSize, const char* mode, const bool transpose, const std::string &compression, uint16_t sampleFormat = 0, uint64_t tileSize = 0, uint64_t samples = 1);

// stripSize is rows per strip; 0 chooses them from the image size, as the
// Python and MATLAB writers do.
uint8_t writeTiffParallelWrapper(const uint64_t x, const uint64_t y, const uint64_t z, const char* fileName, const void* data, const uint64_t bits, const uint64_t startSlice, const uint64_t stripSize, const char* mode, const bool transpose, const std::string &compression, uint16_t sampleFormat = 0, uint64_t tileSize = 0, uint64_t samples = 1);

// Returns 0 on success.
uint8_t writeTiffParallelHelper(const char* fileName, const void* tiffOld, uint64_t bits, const char* mode, uint64_t x, uint64_t y, uint64_t z, uint64_t startSlice, const bool transpose, const std::string &compression, uint16_t sampleFormat = 0, uint64_t tileSize = 0, uint64_t samples = 1);

#endif
