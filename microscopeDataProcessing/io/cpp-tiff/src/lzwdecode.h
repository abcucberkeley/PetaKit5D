#ifndef LZWDECODE_H
#define LZWDECODE_H
#include <cstdint>

// Decodes one LZW strip or tile (TIFF's codes, MSB first) into exactly occ
// bytes of out, the bytes libtiff would decode. Returns false, having written
// any part of out, when libtiff would report the data as bad or decode it
// differently (old-style codes, a code not in the table, too little data):
// the caller then decodes it with libtiff.
bool lzwDecode(const uint8_t* in, uint64_t inBytes, uint8_t* out, uint64_t occ);

#endif // LZWDECODE_H
