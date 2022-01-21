#ifndef PARALLELREADTIFF_H
#define PARALLELREADTIFF_H
#include "parallelReadTiff.c"
void DummyHandler(const char* module, const char* fmt, va_list ap);

void* mallocDynamic(uint64_t x, uint64_t bits);

void readTiffParallel(uint64_t x, uint64_t y, uint64_t z, char* fileName, void* tiff, uint64_t bits, uint64_t startSlice);

void* readTiffParallelWrapper(char* fileName);

#endif // PARALLELREADTIFF_H
