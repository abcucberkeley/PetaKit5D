#ifndef PARALLELREADZARR_H
#define PARALLELREADZARR_H
#include <cstdint>
#include "zarr.h"
uint8_t parallelReadZarr(zarr &Zarr, void* zarrArr,
                         const std::vector<uint64_t> &startCoords, 
                         const std::vector<uint64_t> &endCoords,
                         const std::vector<uint64_t> &readShape,
                         const uint64_t bits,
                         const bool useCtx=false,
                         const bool sparse=false);

void* parallelReadZarrWriteWrapper(zarr Zarr, const bool &crop,
                              std::vector<uint64_t> startCoords, 
                              std::vector<uint64_t> endCoords);

void* readZarrParallelHelper(const char* folderName, 
							 uint64_t startX, uint64_t startY, uint64_t startZ,
							 uint64_t endX, uint64_t endY, uint64_t endZ,
							 uint8_t imageJIm);
#endif
