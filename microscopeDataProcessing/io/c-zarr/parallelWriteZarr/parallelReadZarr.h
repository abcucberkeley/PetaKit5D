#ifndef PARALLELREADZARR_H
#define PARALLELREADZARR_H
#include <stdint.h>

uint64_t dTypeToBits(char* dtype);

void* parallelReadZarrWrapper(char* folderName,uint8_t crop, uint64_t startX, uint64_t startY, uint64_t startZ, uint64_t endX, uint64_t endY,uint64_t endZ);

#endif // PARALLELREADZARR_H
