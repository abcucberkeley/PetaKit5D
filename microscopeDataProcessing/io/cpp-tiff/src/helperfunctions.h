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

uint64_t* getImageSize(const char* fileName);

uint64_t getDataType(const char* fileName);

#endif
