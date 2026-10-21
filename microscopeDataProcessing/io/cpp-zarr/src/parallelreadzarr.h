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

// Same as above, with the output in C order (last axis contiguous) when cOrder
// is true instead of F order (first axis contiguous)
uint8_t parallelReadZarr(zarr &Zarr, void* zarrArr,
                         const std::vector<uint64_t> &startCoords,
                         const std::vector<uint64_t> &endCoords,
                         const std::vector<uint64_t> &readShape,
                         const uint64_t bits,
                         const bool useCtx,
                         const bool sparse,
                         const bool cOrder);

// Same as above, for an output that does not already hold the fill value (memory
// the caller reads into, as it is): the chunks that are missing (and the inner chunks
// missing from their shards) get the element at fillMissing (as fillValueElement
// makes it), and the others their data. With fillMissing NULL, the output is taken
// to hold the fill value already, as above.
uint8_t parallelReadZarr(zarr &Zarr, void* zarrArr,
                         const std::vector<uint64_t> &startCoords,
                         const std::vector<uint64_t> &endCoords,
                         const std::vector<uint64_t> &readShape,
                         const uint64_t bits,
                         const bool useCtx,
                         const bool sparse,
                         const bool cOrder,
                         const void* fillMissing);

void* parallelReadZarrWriteWrapper(zarr Zarr, const bool &crop,
                              std::vector<uint64_t> startCoords,
                              std::vector<uint64_t> endCoords);

// Same as above, returning the array in C order when cOrder is true
void* parallelReadZarrWriteWrapper(zarr Zarr, const bool &crop,
                              std::vector<uint64_t> startCoords,
                              std::vector<uint64_t> endCoords,
                              const bool cOrder);

// Decode one stored chunk into dst (exactly dstLen bytes): check and drop the
// crc32c checksum after it when chunkChecksum, then decompress with compressor
// ("none", "blosc", "gzip" for gzip or zlib streams, or "zstd"). Returns 0, or 1
// with the reason in why.
uint8_t decodeChunk(const std::string &compressor, const bool chunkChecksum,
                    const void* src, uint64_t srcLen, void* dst, const uint64_t dstLen,
                    std::string &why);

void* readZarrParallelHelper(const char* folderName, 
							 uint64_t startX, uint64_t startY, uint64_t startZ,
							 uint64_t endX, uint64_t endY, uint64_t endZ,
							 uint8_t imageJIm);
#endif
