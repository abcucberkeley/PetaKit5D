#ifndef LZWENCODE_H
#define LZWENCODE_H
#include <stdint.h>
#include "tiffio.h"
uint64_t lzwEncode(uint8_t* unCompr, uint8_t* compr, tmsize_t len);
#endif // LZWENCODE_H
