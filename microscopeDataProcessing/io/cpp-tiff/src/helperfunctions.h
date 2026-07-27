#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H
#include <cstdint>
#include <stdarg.h>

#ifndef _WIN32
char* expandTilde(char* path);
#endif

void mkdirRecursive(const char *dir);

void DummyHandler(const char* module, const char* fmt, va_list ap);

uint8_t isImageJIm(const char* fileName);

uint64_t imageJImGetZ(const char* fileName);

uint32_t getImageSizeZ(const char* fileName);

uint64_t* getImageSize(const char* fileName);

uint64_t getDataType(const char* fileName);

// Samples per pixel: 1 = grayscale, 3 = RGB, 4 = RGBA. Defaults to 1.
uint64_t getSamplesPerPixel(const char* fileName);

// Sample format: 1 = unsigned int, 2 = signed int, 3 = IEEE float. Defaults to 1.
uint64_t getSampleFormat(const char* fileName);

#endif
