#ifndef PARALLELWRITEZARR_H
#define PARALLELWRITEZARR_H
#include <cstdint>
#include "zarr.h"

uint8_t parallelWriteZarr(zarr &Zarr, void* zarrArr,
                          const std::vector<uint64_t> &startCoords,
                          const std::vector<uint64_t> &endCoords,
                          const std::vector<uint64_t> &writeShape,
                          const uint64_t bits, const bool useUuid,
                          const bool crop, const bool sparse=true);

// Same as above, for an input array in any memory layout: inStrides gives the
// distance in elements between neighboring elements along each axis (the version
// above takes an F-order array, where the first axis is contiguous)
uint8_t parallelWriteZarr(zarr &Zarr, void* zarrArr,
                          const std::vector<uint64_t> &startCoords,
                          const std::vector<uint64_t> &endCoords,
                          const std::vector<uint64_t> &writeShape,
                          const std::vector<uint64_t> &inStrides,
                          const uint64_t bits, const bool useUuid,
                          const bool crop, const bool sparse=true);
#endif
