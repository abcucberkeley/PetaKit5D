#include "lzwencode.h"
#include <cstdint>
#include <cstdlib>
#include <assert.h>
#include <atomic>
#include <mutex>
#include <vector>
#include "tiffio.h"
typedef size_t WordType;

/*
 * NB: The 5.0 spec describes a different algorithm than Aldus
 *     implements.  Specifically, Aldus does code length transitions
 *     one code earlier than should be done (for real LZW).
 *     Earlier versions of this library implemented the correct
 *     LZW algorithm, but emitted codes in a bit order opposite
 *     to the TIFF spec.  Thus, to maintain compatibility w/ Aldus
 *     we interpret MSB-LSB ordered codes to be images written w/
 *     old versions of this library, but otherwise adhere to the
 *     Aldus "off by one" algorithm.
 *
 * Future revisions to the TIFF spec are expected to "clarify this issue".
 */

#define MAXCODE(n)	((1L<<(n))-1)
/*
 * The TIFF spec specifies that encoded bit
 * strings range from 9 to 12 bits.
 */
#define BITS_MIN        9               /* start with 9 bits */
#define BITS_MAX        12              /* max of 12 bit strings */
/* predefined codes */
#define CODE_CLEAR      256             /* code to clear string table */
#define CODE_EOI        257             /* end-of-information code */
#define CODE_FIRST      258             /* first free code entry */
#define CODE_MAX        MAXCODE(BITS_MAX)

#define CHECK_GAP	10000		/* enc_ratio check interval */

#define	CALCRATIO(sp, rat) {					\
    if (incount > 0x007fffff) { /* NB: shift will overflow */\
    rat = outcount >> 8;				\
    rat = (rat == 0 ? 0x7fffffff : incount/rat);	\
    } else							\
    rat = (incount<<8) / outcount;			\
    }

/* Explicit 0xff masking to make icc -check=conversions happy */
#define	PutNextCode(op, c) {					\
    nextdata = (nextdata << nbits) | c;			\
    nextbits += nbits;					\
    *op++ = (unsigned char)((nextdata >> (nextbits-8))&0xff);		\
    nextbits -= 8;						\
    if (nextbits >= 8) {					\
    *op++ = (unsigned char)((nextdata >> (nextbits-8))&0xff);	\
    nextbits -= 8;					\
    }							\
    outcount += nbits;					\
    }

#define HSIZE           9001L           /* 91% occupancy */
#define HSHIFT          (13-8)

typedef uint16_t hcode_t;			/* codes fit in 16 bits */
typedef struct {
    long	hash;
    hcode_t	code;
} hash_t;

/*
 * Reset encoding hash table.
 */
static void
cl_hash(hash_t* enc_hashtab)
{
    hash_t *hp = &enc_hashtab[HSIZE-1];
    long i = HSIZE-8;

    do {
        i -= 8;
        hp[-7].hash = -1;
        hp[-6].hash = -1;
        hp[-5].hash = -1;
        hp[-4].hash = -1;
        hp[-3].hash = -1;
        hp[-2].hash = -1;
        hp[-1].hash = -1;
        hp[ 0].hash = -1;
        hp -= 8;
    } while (i >= 0);
    for (i += 8; i > 0; i--, hp--)
        hp->hash = -1;
}

// libtiff's encoder, with its 9001-entry hash table (144 KB) as the string
// table. Used for short strips, and until a 2 MB direct table pays for itself
// (see lzwEncode). Neither encoder is inlined: inlined into lzwEncode, how
// GCC 13 compiled their loops depended on the code around them (6-10% either
// way).
__attribute__((noinline)) static uint64_t encodeWithHash(uint8_t* unCompr, uint8_t* compr, tmsize_t cc){
    hash_t* enc_hashtab = (hash_t*)_TIFFmalloc(HSIZE*sizeof (hash_t));
    if(!enc_hashtab) return 0;
    uint8_t* bp = unCompr;
    long fcode;
    hash_t *hp;
    int h, c;
    hcode_t ent;
    long disp;
    tmsize_t incount = 0, outcount = 0, checkpoint = CHECK_GAP, enc_ratio = 0;
    WordType nextdata = 0;
    long nextbits = 0;
    int free_ent = CODE_FIRST, maxcode = MAXCODE(BITS_MIN), nbits = BITS_MIN;
    uint8_t* op = compr;

    cl_hash(enc_hashtab);
    ent = (hcode_t) -1;
    if (cc > 0) {
        PutNextCode(op, CODE_CLEAR);
        ent = *bp++; cc--; incount++;
    }
    while (cc > 0) {
        c = *bp++; cc--; incount++;
        fcode = ((long)c << BITS_MAX) + ent;
        h = (c << HSHIFT) ^ ent;	/* xor hashing */
        #ifdef _WINDOWS
        /*
        * Check hash index for an overflow.
        */
        if (h >= HSIZE)
            h -= HSIZE;
        #endif
        hp = &enc_hashtab[h];
        if (hp->hash == fcode) {
            ent = hp->code;
            continue;
        }
        if (hp->hash >= 0) {
            /*
            * Primary hash failed, check secondary hash.
            */
            disp = HSIZE - h;
            if (h == 0)
                disp = 1;
            do {
                /*
                 * Avoid pointer arithmetic because of
                 * wraparound problems with segments.
                 */
                if ((h -= disp) < 0)
                    h += HSIZE;
                hp = &enc_hashtab[h];
                if (hp->hash == fcode) {
                    ent = hp->code;
                    goto hit;
                }
            } while (hp->hash >= 0);
        }
        /*
         * New entry, emit code and add to table.
         */
        PutNextCode(op, ent);
        ent = (hcode_t)c;
        hp->code = (hcode_t)(free_ent++);
        hp->hash = fcode;
        if (free_ent == CODE_MAX-1) {
            /* table is full, emit clear code and reset */
            cl_hash(enc_hashtab);
            enc_ratio = 0;
            incount = 0;
            outcount = 0;
            free_ent = CODE_FIRST;
            PutNextCode(op, CODE_CLEAR);
            nbits = BITS_MIN;
            maxcode = MAXCODE(BITS_MIN);
        } else {
            /*
             * If the next entry is going to be too big for
             * the code size, then increase it, if possible.
             */
            if (free_ent > maxcode) {
                nbits++;
                assert(nbits <= BITS_MAX);
                maxcode = (int) MAXCODE(nbits);
            } else if (incount >= checkpoint) {
                tmsize_t rat;
                /*
                 * Check compression ratio and, if things seem
                 * to be slipping, clear the hash table and
                 * reset state.  The compression ratio is a
                 * 24+8-bit fractional number.
                 */
                checkpoint = incount+CHECK_GAP;
                CALCRATIO(sp, rat);
                if (rat <= enc_ratio) {
                    cl_hash(enc_hashtab);
                    enc_ratio = 0;
                    incount = 0;
                    outcount = 0;
                    free_ent = CODE_FIRST;
                    PutNextCode(op, CODE_CLEAR);
                    nbits = BITS_MIN;
                    maxcode = MAXCODE(BITS_MIN);
                } else
                    enc_ratio = rat;
            }
        }
hit:
        ;
    }
    if (ent != (hcode_t) -1) {
        PutNextCode(op, ent);
        free_ent++;
        if (free_ent == CODE_MAX-1) {
            /* table is full, emit clear code and reset */
            PutNextCode(op, CODE_CLEAR);
            nbits = BITS_MIN;
        } else if (free_ent > maxcode) {
            nbits++;
            assert(nbits <= BITS_MAX);
        }
    }
    PutNextCode(op, CODE_EOI);
    /* Explicit 0xff masking to make icc -check=conversions happy */
    if (nextbits > 0)
        *op++ = (unsigned char)((nextdata << (8-nextbits))&0xff);
    (void)outcount;
    _TIFFfree(enc_hashtab);
    return (uint64_t)(op - compr);
}

namespace {

// The string table. child[c << 12 | ent] is the code of string ent followed by
// byte c, or 0 if that string is not in the table (codes below 258 are never
// added, so 0 is free). libtiff finds strings through a 9001-entry hash table
// instead; indexing directly finds the same strings, so the codes and the
// output are the same, but encodes 1.3-3x faster. The table is laid out by
// byte first because image data repeats few byte values in some positions
// (the high bytes of 16-bit samples), which keeps those lookups close together.
// A reset empties only the entries added, so the 2 MB table is never cleared
// as a whole.
struct LzwTable {
    uint16_t child[256 << BITS_MAX];
    uint32_t added[1 << BITS_MAX];
};

// Tables kept between encodes, so a strip does not allocate and fault in 2 MB.
// There are as many as the most strips ever encoded at once with them.
struct TablePool {
    std::mutex m;
    std::vector<LzwTable*> tables;
    ~TablePool(){
        for(LzwTable* t : tables) free(t);
    }
};
TablePool pool;

// A new table costs its first touch: faulting in its 2 MB took about a
// millisecond, and threads doing it at once wait on each other (8 tables 5 ms,
// 47 tables 18 ms). Made for every thread at once, that made the first LZW
// writes of a process slower than with libtiff's hash table (2048x2048 images
// 1.4-1.7x, 553 MB stacks with 48 threads 1.1x). So tables are earned, a new
// one only after the process has encoded another 4 MB with the hash table
// (about 40 ms of work), and made one at a time while the other threads go on
// with the hash table.
const uint64_t kEarnBytes = (uint64_t)4 << 20;
const tmsize_t kTableMinStrip = (tmsize_t)128 << 10;
std::atomic<uint64_t> hashEncoded(0), tablesMade(0);
std::atomic<bool> making(false);

// A kept table, or a new one if the process has earned it and no other thread
// is making one, or NULL.
LzwTable* takeTable(){
    {
        std::lock_guard<std::mutex> g(pool.m);
        if(!pool.tables.empty()){
            LzwTable* t = pool.tables.back();
            pool.tables.pop_back();
            return t;
        }
    }
    uint64_t made = tablesMade.load(std::memory_order_relaxed);
    bool idle = false;
    if(made >= hashEncoded.load(std::memory_order_relaxed)/kEarnBytes || !making.compare_exchange_strong(idle, true)) return NULL;
    LzwTable* t = NULL;
    if(tablesMade.compare_exchange_strong(made, made + 1)){
        t = (LzwTable*)calloc(1, sizeof(LzwTable));
        // Fault it all in now, not during the encode, so the next table waits
        if(t) for(size_t i = 0; i < sizeof(LzwTable); i += 4096) reinterpret_cast<volatile uint8_t*>(t)[i] = 0;
    }
    making.store(false);
    return t;
}

void returnTable(LzwTable* t){
    std::lock_guard<std::mutex> g(pool.m);
    pool.tables.push_back(t);
}


}

// libtiff's encoder with a direct table instead of its hash table; the
// table is left empty.
__attribute__((noinline)) static uint64_t encodeWithTable(LzwTable* table, uint8_t* unCompr, uint8_t* compr, tmsize_t cc){
    uint16_t* child = table->child;
    uint32_t* added = table->added;
    int nAdded = 0;
#define RESET_TABLE { for(int i_ = 0; i_ < nAdded; i_++) child[added[i_]] = 0; nAdded = 0; }

    uint8_t* bp = unCompr;
    uint8_t* op = compr;
    tmsize_t incount = 0, outcount = 0, checkpoint = CHECK_GAP, enc_ratio = 0;
    WordType nextdata = 0;
    long nextbits = 0;
    int free_ent = CODE_FIRST, maxcode = MAXCODE(BITS_MIN), nbits = BITS_MIN;
    int ent = -1;

    if (cc > 0) {
        PutNextCode(op, CODE_CLEAR);
        ent = *bp++; cc--; incount++;
    }
    while (cc > 0) {
        const uint32_t c = *bp++; cc--; incount++;
        const uint32_t idx = (c << BITS_MAX) | (uint32_t)ent;
        const uint16_t code = child[idx];
        if (code) {
            ent = code;
            continue;
        }
        /*
         * New entry, emit code and add to table.
         */
        PutNextCode(op, ent);
        ent = (int)c;
        child[idx] = (uint16_t)(free_ent++);
        added[nAdded++] = idx;
        if (free_ent == CODE_MAX-1) {
            /* table is full, emit clear code and reset */
            RESET_TABLE
            enc_ratio = 0;
            incount = 0;
            outcount = 0;
            free_ent = CODE_FIRST;
            PutNextCode(op, CODE_CLEAR);
            nbits = BITS_MIN;
            maxcode = MAXCODE(BITS_MIN);
        } else {
            /*
             * If the next entry is going to be too big for
             * the code size, then increase it, if possible.
             */
            if (free_ent > maxcode) {
                nbits++;
                assert(nbits <= BITS_MAX);
                maxcode = (int) MAXCODE(nbits);
            } else if (incount >= checkpoint) {
                tmsize_t rat;
                /*
                 * Check compression ratio and, if things seem
                 * to be slipping, clear the hash table and
                 * reset state.  The compression ratio is a
                 * 24+8-bit fractional number.
                 */
                checkpoint = incount+CHECK_GAP;
                CALCRATIO(sp, rat);
                if (rat <= enc_ratio) {
                    RESET_TABLE
                    enc_ratio = 0;
                    incount = 0;
                    outcount = 0;
                    free_ent = CODE_FIRST;
                    PutNextCode(op, CODE_CLEAR);
                    nbits = BITS_MIN;
                    maxcode = MAXCODE(BITS_MIN);
                } else
                    enc_ratio = rat;
            }
        }
    }
    RESET_TABLE
#undef RESET_TABLE

    if (ent != -1) {
        // Flush the final pending code.
        PutNextCode(op, ent);
        ent = -1;
        free_ent ++;

        if (free_ent == CODE_MAX-1) {
            /* table is full, emit clear code and reset */
            PutNextCode(op, CODE_CLEAR);
            nbits = BITS_MIN;
        } else {
            /*
             * If the next entry is going to be too big for
             * the code size, then increase it, if possible.
             */
            if (free_ent > maxcode) {
                nbits++;
                assert(nbits <= BITS_MAX);
            }
        }
    }
    PutNextCode(op, CODE_EOI);
    /* Explicit 0xff masking to make icc -check=conversions happy */
    if (nextbits > 0){
        *op++ = (unsigned char)((nextdata << (8-nextbits))&0xff);
        outcount += nbits;
    }
    (void)outcount;
    // The bytes actually written. Rounding the bit count up instead could count
    // a byte or two of the uninitialized buffer past the last code, so the same
    // data did not always produce the same file.
    return (uint64_t)(op - compr);
}

// Encodes len bytes as one LZW strip (libtiff's encoder: the same codes, code
// widths and clear codes, so the same bytes, whichever string table finds
// them). Returns the bytes written to compr, or 0 if no table could be had.
uint64_t lzwEncode(uint8_t* unCompr, uint8_t* compr, tmsize_t cc){
    // A short strip starts with a cold table: on a 24-core Skylake node, many
    // small files written with every core took 1.1x as long with the tables
    // (23 of them, 46 MB, crowding the caches) as with the hash table.
    if(cc < kTableMinStrip) return encodeWithHash(unCompr, compr, cc);
    LzwTable* table = takeTable();
    if(!table){
        hashEncoded.fetch_add((uint64_t)cc, std::memory_order_relaxed);
        return encodeWithHash(unCompr, compr, cc);
    }
    const uint64_t n = encodeWithTable(table, unCompr, compr, cc);
    returnTable(table);
    return n;
}
