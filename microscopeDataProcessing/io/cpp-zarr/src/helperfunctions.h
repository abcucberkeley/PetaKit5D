#ifndef CPPZARR_HELPERFUNCTIONS_H
#define CPPZARR_HELPERFUNCTIONS_H
#include <cstdint>
#include <string>
#include <vector>

#ifndef _WIN32
const char* expandTilde(const char* path);
#endif
    
#ifdef _WIN32
char* strndup (const char *s, size_t n);

int _vscprintf_so(const char * format, va_list pargs);

int vasprintf(char **strp, const char *fmt, va_list ap);

int asprintf(char *strp[], const char *fmt, ...);
#endif

std::string generateUUID();

void mkdirRecursive(const char *dir);

bool fileExists(const std::string &fileName);

// Move src onto dst, replacing dst if it already exists (the POSIX rename()
// behavior; on Windows rename() fails instead when dst exists). Returns true on
// success.
bool renameReplace(const std::string &src, const std::string &dst);

void makeDimensionFolders(const std::string &fileName);

bool isLittleEndian();

// The number of threads the calling thread's reads and writes use (OpenMP's setting
// for that thread: OMP_NUM_THREADS, or the CPUs it may run on), and setting it for
// the calling thread only. (The library's own OpenMP runtime: the Python and MATLAB
// bindings need not link one.)
int callThreads();
void setCallThreads(int n);

// A data type as cpp-zarr keeps it: the dtype string of a Zarr v2 .zarray ("<u2",
// "|b1", "<f2", "<c16", ...), or "bfloat16", the name Zarr metadata gives that type
// (">bfloat16" for a Zarr v3 array stored big-endian)
struct DtypeInfo {
    // 'b' bool, 'i' signed and 'u' unsigned integers, 'f' floats, 'c' complex
    // numbers; 0 for a data type cpp-zarr does not support
    char kind = 0;
    // Bytes per element: 1 (bool, 8-bit integers), 2, 4 or 8; complex numbers are
    // 8 or 16 (both parts)
    uint64_t bytes = 0;
    // Stored most significant byte first (each part of a complex number)
    bool bigEndian = false;
    // A 2-byte float that is bfloat16 instead of IEEE half precision (float16)
    bool bfloat16 = false;
    // Bytes of each value whose byte order is swapped: a complex number's parts
    uint64_t swapBytes() const { return kind == 'c' ? bytes/2 : bytes; }
};
DtypeInfo dtypeInfo(const std::string &dtype);

// Whether data of the dtype is stored in the other byte order than this machine's
bool oppositeEndianness(const std::string &dtype);

// IEEE half precision (float16) and bfloat16 values, as their bits, from and to
// float or double, rounding to the nearest value (ties to even). A double becomes
// float16 directly, as NumPy converts, and bfloat16 through float, as ml_dtypes
// does.
uint16_t floatToHalf(float v);
uint16_t doubleToHalf(double v);
float halfToFloat(uint16_t h);
uint16_t floatToBfloat16(float v);
float bfloat16ToFloat(uint16_t b);

// Parse a zarr fill_value string for memset-style fills. Returns 0 for anything
// unparsable or outside int range (e.g. the int64 min/max strings that
// Infinity/-Infinity fill values normalize to), instead of throwing like stoi.
// Note: memset-based fills are only exact for 0 (and -1) on multi-byte dtypes;
// any other nonzero fill repeats its low byte across the element.
int fillValueToInt(const std::string &fillValue);

void swapArrayEndianness(void* array, const size_t elementSize, const size_t numElements);

// The same for n elements of a dtype: each part of a complex number is swapped
void swapDtypeEndianness(void* array, const DtypeInfo &dtype, const size_t n);

// Copy an n0 x n1 x n2 box of elements between two strided layouts (strides in
// elements, element size `bytes` = 1, 2, 4, 8 or 16):
//     dst[i*d0 + j*d1 + k*d2] = src[i*s0 + j*s1 + k*s2]
// Used to convert between F-order (x fastest) and C-order (z fastest) data: i and
// k are processed in small tiles so each tile's source and destination cache lines
// stay in L1, instead of striding through memory one element at a time.
void copyBoxTransposed(const uint64_t bytes, const void* src, void* dst,
                       const uint64_t n0, const uint64_t n1, const uint64_t n2,
                       const uint64_t s0, const uint64_t s1, const uint64_t s2,
                       const uint64_t d0, const uint64_t d1, const uint64_t d2);

// Set n elements of size `bytes` to the element whose bytes are at `elem`.
void fillElements(const uint64_t bytes, void* dst, const uint64_t n, const void* elem);

// CRC32C (Castagnoli) checksum, as Zarr v3 uses for shard indexes.
uint32_t crc32c(const uint8_t* data, size_t length);

// The fill_value string as one element of the dtype ("<u2", ">f4", "|b1", "<c8",
// "bfloat16", ...), in the dtype's own byte order: a number, "NaN", "Infinity",
// "-Infinity" or a "0x..." bit pattern (Zarr v3), or for complex numbers a JSON pair
// of those, [real, imaginary]. Writes 1 to 16 bytes; unparsable values give 0.
void fillValueElement(const std::string &fillValue, const std::string &dtype, uint8_t elem[16]);

// True when the fill_value string is a numeric zero ("0", "0.0", ...): chunks of
// zeros can then be left unwritten, since every reader fills a missing chunk with
// zeros. ("NaN", null and nonzero fill values are not.)
bool fillValueIsZero(const std::string &fillValue);

// True when the fill value is an element whose bytes are all zero in the dtype
// (not -0.0, nor a complex number with a nonzero part)
bool fillValueIsZero(const std::string &fillValue, const std::string &dtype);

// Copy an N-dimensional box of elements between two strided layouts (strides in
// elements, element size `bytes` = 1, 2, 4, 8 or 16):
//     dst[sum_d i_d*dstStrides[d]] = src[sum_d i_d*srcStrides[d]],  0 <= i_d < extents[d]
// When both sides are contiguous along the same axis (e.g. F to F) the box is
// copied as contiguous runs; when they differ (F <-> C order) it uses the tiled
// transpose of copyBoxTransposed, one plane at a time. The runs are visited with
// the smallest destination stride innermost, or with spreadWrites the largest
// destination stride innermost (spreads consecutive writes across pages, which
// is faster when the destination is freshly allocated memory).
void copyBoxND(const uint64_t bytes, const void* src, void* dst,
               const std::vector<uint64_t> &extents,
               const std::vector<uint64_t> &srcStrides,
               const std::vector<uint64_t> &dstStrides,
               const bool spreadWrites = false);

// Set every element of an N-dimensional box (extents, strides in elements) of dst
// to the element whose bytes are at `elem`
void fillBoxND(const uint64_t bytes, void* dst,
               const std::vector<uint64_t> &extents,
               const std::vector<uint64_t> &strides,
               const void* elem);

// Set every element of a packed N-dimensional buffer (full extents `full`,
// layout `strides`, F or C order) that lies outside the box [lo, hi) to the
// element whose bytes are at `elem`.
void fillOutsideBoxND(const uint64_t bytes, void* buf,
                      const std::vector<uint64_t> &full,
                      const std::vector<uint64_t> &strides,
                      const std::vector<uint64_t> &lo,
                      const std::vector<uint64_t> &hi,
                      const void* elem);

#endif
