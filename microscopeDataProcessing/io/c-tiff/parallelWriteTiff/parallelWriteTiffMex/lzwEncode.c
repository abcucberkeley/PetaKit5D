#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include "tiffio.h"
#include "lzwEncode.h"
#define SIZEOF_WORDTYPE SIZEOF_SIZE_T
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
#define LZW_COMPAT              /* include backwards compatibility code */

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
#define HSIZE           9001L           /* 91% occupancy */
#define HSHIFT          (13-8)
#ifdef LZW_COMPAT
/* NB: +1024 is for compatibility with old files */
#define CSIZE           (MAXCODE(BITS_MAX)+1024L)
#else
#define CSIZE           (MAXCODE(BITS_MAX)+1L)
#endif


/*
 * Encoding-specific state.
 */
typedef uint16_t hcode_t;			/* codes fit in 16 bits */
typedef struct {
    long	hash;
    hcode_t	code;
} hash_t;

#define CHECK_GAP	10000		/* enc_ratio check interval */

/*
 * The 4 here insures there is space for 2 max-sized
 * codes in LZWEncode and LZWPostDecode.
 */
//uint8_t* enc_rawlimit = tif->tif_rawdata + tif->tif_rawdatasize-1 - 4;


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

/*
 * Reset encoding hash table.
 */
static void
cl_hash(hash_t* enc_hashtab)
{
    register hash_t *hp = &enc_hashtab[HSIZE-1];
    register long i = HSIZE-8;

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


uint64_t lzwEncode(uint8_t* unCompr, uint8_t* compr, tmsize_t cc){
    hash_t* enc_hashtab = (hash_t*)_TIFFmalloc(HSIZE*sizeof (hash_t));
    //compr = (uint8_t*)malloc(cc);
    uint8_t* bp = unCompr;
    register long fcode;
    register hash_t *hp;
    register int h, c;
    hcode_t ent;
    long disp;
    tmsize_t incount, outcount, checkpoint;
    uint64_t totalOut = 0;
    WordType nextdata;
    long nextbits;
    int free_ent, maxcode, nbits;
    uint8_t* op;
    uint8_t* limit;

    //(void) s;
    unsigned short lzw_nbits = BITS_MIN;
    unsigned short lzw_maxcode = MAXCODE(BITS_MIN);
    unsigned short lzw_free_ent = CODE_FIRST;
    long lzw_nextbits = 0;
    WordType lzw_nextdata = 0;

    int enc_oldcode = (hcode_t) -1;	/* generates CODE_CLEAR in LZWEncode */
    tmsize_t enc_checkpoint = CHECK_GAP;
    tmsize_t enc_ratio = 0;
    tmsize_t enc_incount = 0;
    tmsize_t enc_outcount = 0;
    cl_hash(enc_hashtab);

    ent = (hcode_t)enc_oldcode;

    incount = enc_incount;
    outcount = enc_outcount;
    checkpoint = enc_checkpoint;
    nextdata = lzw_nextdata;
    nextbits = lzw_nextbits;
    free_ent = lzw_free_ent;
    maxcode = lzw_maxcode;
    nbits = lzw_nbits;
    //op = tif->tif_rawcp;
    op = compr;
    //limit = enc_rawlimit;


    if (ent == (hcode_t) -1 && cc > 0) {
        /*
             * NB: This is safe because it can only happen
             *     at the start of a strip where we know there
             *     is space in the data buffer.
             */
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
        /*
                     * Verify there is space in the buffer for the code
                     * and any potential Clear code that might be emitted
                     * below.  The value of limit is setup so that there
                     * are at least 4 bytes free--room for 2 codes.
                     */
        /*
                    if (op > limit) {
                        tif->tif_rawcc = (tmsize_t)(op - tif->tif_rawdata);
                        if( !TIFFFlushData1(tif) )
                                        return 0;
                        op = tif->tif_rawdata;
                    }*/
        PutNextCode(op, ent);
        ent = (hcode_t)c;
        hp->code = (hcode_t)(free_ent++);
        hp->hash = fcode;
        if (free_ent == CODE_MAX-1) {
            /* table is full, emit clear code and reset */
            cl_hash(enc_hashtab);
            enc_ratio = 0;
            incount = 0;
            totalOut += outcount;
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
                    totalOut += outcount;
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



    // POST ENCODE
    /*
    if (op > enc_rawlimit) {
        tif->tif_rawcc = (tmsize_t)(op - tif->tif_rawdata);
        if( !TIFFFlushData1(tif) )
                    return 0;
        op = tif->tif_rawdata;
    }*/
    if (enc_oldcode != (hcode_t) -1) {
        free_ent = lzw_free_ent;

        PutNextCode(op, enc_oldcode);
        enc_oldcode = (hcode_t) -1;
        free_ent ++;

        if (free_ent == CODE_MAX-1) {
            /* table is full, emit clear code and reset */
            totalOut += outcount;
            outcount = 0;
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
    //tif->tif_rawcc = (tmsize_t)(op - tif->tif_rawdata);

    //(void)outcount;

    //return outcount*8;
    totalOut += outcount;
    //return (totalOut/8);
    _TIFFfree(enc_hashtab);
    return (uint64_t)ceil((double)totalOut/8.0);//+4;
}



