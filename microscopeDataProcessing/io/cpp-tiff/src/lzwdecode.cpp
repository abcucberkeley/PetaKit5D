#include "lzwdecode.h"
#include <cstring>

namespace {

enum { CODE_CLEAR = 256, CODE_EOI = 257, CODE_FIRST = 258, BITS_MIN = 9, BITS_MAX = 12,
       // libtiff's table size: 4095 codes plus room for streams whose clear code is late
       CSIZE = 4095 + 1024 };

inline uint64_t load64BigEndian(const uint8_t* p){
    uint64_t v;
    memcpy(&v, p, 8);
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return v;
#else
    return __builtin_bswap64(v);
#endif
}

}

// libtiff decodes a code by walking its string's chain of table entries back
// to front. Here every string is copied instead from where it was already
// decoded: a new code's string is the previous code's string plus one byte,
// and in the output that byte follows it, so each code is just an offset and a
// length into the output.
bool lzwDecode(const uint8_t* in, uint64_t inBytes, uint8_t* out, uint64_t occ){
    // Old-style (bit-reversed) codes, which libtiff decodes separately
    if(inBytes < 2 || (in[0] == 0 && (in[1] & 1)) || occ == 0 || occ > UINT32_MAX) return false;
    uint32_t start[CSIZE];
    uint16_t len[CSIZE];
    uint8_t* const end = out + occ;
    uint8_t* op = out;
    const uint8_t* bp = in;
    const uint8_t* const bend = in + inBytes;
    uint64_t bitbuf = 0;
    int bits = 0;
    int nbits = BITS_MIN;
    uint32_t mask = (1u << BITS_MIN) - 1, freeEnt = 0, maxcode = mask - 1;
    uint32_t oldStart = 0, oldLen = 0;   // the previous code's string, in the output

    // The next code, MSB first. A stream that ends without a code to finish
    // the output is libtiff's error to report.
#define NEXT_CODE(code) {                                                              \
        if(bits < nbits){                                                              \
            if(bend - bp >= 8){                                                        \
                const int take = (63 - bits) >> 3;                                     \
                bitbuf = (bitbuf << (take*8)) | (load64BigEndian(bp) >> (64 - take*8)); \
                bp += take; bits += take*8;                                            \
            }                                                                          \
            else{                                                                      \
                while(bits < nbits && bp < bend){ bitbuf = (bitbuf << 8) | *bp++; bits += 8; } \
                if(bits < nbits) return false;                                         \
            }                                                                          \
        }                                                                              \
        code = (uint32_t)(bitbuf >> (bits - nbits)) & mask;                            \
        bits -= nbits;                                                                 \
    }
    // A new table entry: the previous string plus one byte. The code width grows
    // one code early, as libtiff's (and Aldus's) decoder does.
#define ADD_ENTRY {                                                                    \
        start[freeEnt] = oldStart; len[freeEnt] = (uint16_t)(oldLen + 1);              \
        if(++freeEnt > maxcode){                                                       \
            if(++nbits > BITS_MAX) nbits = BITS_MAX;                                   \
            mask = (1u << nbits) - 1; maxcode = mask - 1;                              \
        }                                                                              \
    }

    uint32_t code;
    NEXT_CODE(code)
    // libtiff requires a clear code first
    if(code != CODE_CLEAR) return false;
    for(;;){
        if(code == CODE_CLEAR){
            freeEnt = CODE_FIRST; nbits = BITS_MIN; mask = (1u << BITS_MIN) - 1; maxcode = mask - 1;
            do { NEXT_CODE(code) } while(code == CODE_CLEAR);
            if(code >= CODE_EOI) return false;   // the end before any data, or a code not in the table
            oldStart = (uint32_t)(op - out); oldLen = 1;
            *op++ = (uint8_t)code;
            if(op == end) return true;
            NEXT_CODE(code)
            continue;
        }
        // Not enough data for the output (libtiff reports it), or a code not
        // yet in the table, or a full table that was not cleared
        if(code == CODE_EOI || freeEnt >= CSIZE || (code > freeEnt && code >= CODE_FIRST)) return false;
        if(code < 256){
            ADD_ENTRY
            oldStart = (uint32_t)(op - out); oldLen = 1;
            *op++ = (uint8_t)code;
            if(op == end) return true;
            NEXT_CODE(code)
            continue;
        }
        uint32_t s, l;
        if(code == freeEnt){ s = oldStart; l = oldLen + 1; }   // the string being added: the previous one plus its first byte
        else { s = start[code]; l = len[code]; }
        ADD_ENTRY
        const uint32_t pos = (uint32_t)(op - out);
        oldStart = pos; oldLen = l;
        const uint8_t* src = out + s;
        const uint64_t room = (uint64_t)(end - op);
        if(room >= (uint64_t)l + 16 && pos - s >= 16){
            // 16 bytes at a time; the bytes copied past the string are
            // overwritten by the codes that follow
            uint8_t* d = op;
            int64_t n = l;
            do { memcpy(d, src, 16); d += 16; src += 16; n -= 16; } while(n > 0);
            op += l;
        }
        else if(room <= l){
            // The string fills the rest of the output: libtiff stops there too
            for(uint64_t i = 0; i < room; i++) op[i] = src[i];
            return true;
        }
        else{
            // Overlapping (a string repeating itself) or near the end
            for(uint32_t i = 0; i < l; i++) op[i] = src[i];
            op += l;
        }
        NEXT_CODE(code)
    }
#undef NEXT_CODE
#undef ADD_ENTRY
}
