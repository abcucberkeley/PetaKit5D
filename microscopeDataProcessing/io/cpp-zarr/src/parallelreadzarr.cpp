#include <algorithm>
#include <chrono>
#include <memory>
#include <limits>
#include <atomic>
#include <fstream>
#include <functional>
#include <cstdint>
#include <cstdlib>
#include <omp.h>
#include "parallelreadzarr.h"
#include "blosc2.h"
#include "zarr.h"
#include "helperfunctions.h"
#include "parallelwork.h"
#include "zlib.h"
#include "zstd.h"
#ifndef _WIN32
#include <cerrno>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#ifdef __linux__
#include <dirent.h>
#include <mutex>
#include <sys/syscall.h>
#include <sys/sysmacros.h>
#include <sys/uio.h>
// preadv2 (Linux 4.6), which the page-cache probe below uses, where the system headers
// predate it (as the 3.10 headers of manylinux2014, which Python wheels can be built on);
// a kernel without it fails the call, and the probe then has no answer
#ifndef SYS_preadv2
#if defined(__x86_64__)
#define SYS_preadv2 327
#elif defined(__aarch64__)
#define SYS_preadv2 286
#endif
#endif
#endif

#if defined(__GNUC__)
#define CPPZARR_NOINLINE __attribute__((noinline))
#else
#define CPPZARR_NOINLINE
#endif

namespace {
// A shard file a thread keeps open while it reads from it: read at offsets with pread, one
// system call a read (a seek and a read with a stream), and its descriptor can be given
// read-ahead advice; a stream on Windows, which has no pread
class ShardFile {
#ifdef _WIN32
    std::ifstream f;
#else
    int fd = -1;
#endif
public:
    ShardFile(){}
    ShardFile(const ShardFile&) = delete;
    ShardFile& operator=(const ShardFile&) = delete;
    ~ShardFile(){ close(); }
    bool open(const std::string &path){
        close();
    #ifdef _WIN32
        f.open(path, std::ios::binary);
        return f.is_open();
    #else
        fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
        return fd >= 0;
    #endif
    }
    void close(){
    #ifdef _WIN32
        if(f.is_open()) f.close();
        f.clear();
    #else
        if(fd >= 0) ::close(fd);
        fd = -1;
    #endif
    }
    // Its size in bytes, or -1
    int64_t size(){
    #ifdef _WIN32
        f.clear();
        f.seekg(0, std::ios::end);
        return f ? (int64_t)f.tellg() : -1;
    #else
        struct stat st;
        return fstat(fd, &st) == 0 ? (int64_t)st.st_size : -1;
    #endif
    }
    // Read n bytes at offset off into buf; false if they could not all be read
    bool readAt(uint64_t off, void* buf, uint64_t n){
    #ifdef _WIN32
        f.clear();
        f.seekg((std::streamoff)off, std::ios::beg);
        f.read(static_cast<char*>(buf), (std::streamsize)n);
        return (bool)f;
    #else
        uint8_t* p = static_cast<uint8_t*>(buf);
        while(n > 0){
            const ssize_t r = pread(fd, p, n, (off_t)off);
            if(r < 0 && errno == EINTR) continue;
            if(r <= 0) return false;
            p += r;
            off += (uint64_t)r;
            n -= (uint64_t)r;
        }
        return true;
    #endif
    }
    // Its descriptor (-1 on Windows)
    int handle() const {
    #ifdef _WIN32
        return -1;
    #else
        return fd;
    #endif
    }
};

#ifdef __linux__
// What reading ahead needs to know of a file system's block device, found once per device:
// whether it is an NVMe drive, or a device-mapper or md device (LVM, RAID) on NVMe drives
// only (its sysfs path, or those of the devices under it, name an NVMe controller); and the
// most one read-ahead call reads, the larger of its largest request and its read-ahead
// (max_sectors_kb and read_ahead_kb, a partition's those of its disk), 0 if not known.
// (128 KB on many NVMe drives: asked in 512 KB steps, such a drive read a quarter.)
struct BlockDevice {
    bool nvme = false;
    uint64_t readAheadMax = 0;
};
BlockDevice blockDevice(const dev_t dev){
    static std::mutex m;
    static std::vector<std::pair<dev_t, BlockDevice>> known;
    {
        std::lock_guard<std::mutex> g(m);
        for(const auto &k : known) if(k.first == dev) return k.second;
    }
    std::function<bool(const std::string&, int)> nvme = [&](const std::string &path, const int depth) -> bool {
        char* real = realpath(path.c_str(), nullptr);
        if(!real) return false;
        const std::string r(real);
        free(real);
        if(r.find("/nvme") != std::string::npos) return true;
        DIR* d = depth < 4 ? opendir((r+"/slaves").c_str()) : nullptr;
        if(!d) return false;
        bool any = false, all = true;
        while(const struct dirent* e = readdir(d)){
            if(e->d_name[0] == '.') continue;
            any = true;
            all = all && nvme("/sys/class/block/"+std::string(e->d_name), depth+1);
        }
        closedir(d);
        return any && all;
    };
    const std::string sys = "/sys/dev/block/"+std::to_string(major(dev))+":"+std::to_string(minor(dev));
    BlockDevice v;
    v.nvme = nvme(sys, 0);
    auto kb = [](const std::string &path) -> uint64_t {
        std::ifstream f(path);
        uint64_t n = 0;
        return (f >> n) ? n : 0;
    };
    char* real = realpath(sys.c_str(), nullptr);
    if(real){
        const std::string r(real);
        free(real);
        for(const std::string &q : {r+"/queue/", r+"/../queue/"}){
            const uint64_t most = std::max(kb(q+"max_sectors_kb"), kb(q+"read_ahead_kb"));
            if(most){
                v.readAheadMax = most << 10;
                break;
            }
        }
    }
    std::lock_guard<std::mutex> g(m);
    known.push_back({dev, v});
    return v;
}
#endif

// A reader thread's buffers and decoders, kept in a pool between calls (see
// ScratchPool). Buffers over 1 MB (and the decoder used with them) are not kept.
struct ReadScratch {
    std::unique_ptr<uint8_t[]> stored, decoded;
    uint64_t storedCap = 0, decodedCap = 0;
    std::vector<uint64_t> shardIndex;
    // The blocks of a blosc chunk to leave out when decoding part of it, and a part of a
    // chunk decoded on its own (up to 1 MB)
    std::unique_ptr<bool[]> skip;
    uint64_t skipCap = 0;
    std::unique_ptr<uint8_t[]> part;
    uint64_t partCap = 0;
    // The largest blocks dctx has decoded parts of chunks in (its buffers grow with them)
    uint64_t largestBlock = 0;
    // A one-thread blosc decoder and a zstd decoder
    blosc2_context* dctx = nullptr;
    ZSTD_DCtx* zctx = nullptr;
    ReadScratch(){}
    ReadScratch(const ReadScratch&) = delete;
    ReadScratch& operator=(const ReadScratch&) = delete;
    ~ReadScratch(){
        if(dctx) blosc2_free_ctx(dctx);
        if(zctx) ZSTD_freeDCtx(zctx);
    }
    // Room for at least n bytes; when it has to grow it takes least at once (as much as
    // any chunk of the array can need): a buffer that grows from call to call is mapped
    // and unmapped again each time, which stalls every thread of the process
    uint8_t* storedBuffer(const uint64_t n, const uint64_t least){
        if(storedCap < n){
            const uint64_t cap = std::max(n, least);
            stored.reset();
            stored.reset(new uint8_t[cap]);
            storedCap = cap;
        }
        return stored.get();
    }
    uint8_t* decodedBuffer(const uint64_t n){
        if(decodedCap < n){
            decoded.reset();
            decoded.reset(new uint8_t[n]);
            decodedCap = n;
        }
        return decoded.get();
    }
    uint8_t* partBuffer(const uint64_t n){
        if(partCap < n){
            part.reset();
            part.reset(new uint8_t[n]);
            partCap = n;
        }
        return part.get();
    }
    void shrink(){
        const uint64_t keep = 1 << 20;
        if(storedCap > keep){
            stored.reset();
            storedCap = 0;
        }
        // (blosc's buffers grow with the blocks of the chunks it decodes)
        if(decodedCap > keep || largestBlock > keep/4){
            decoded.reset();
            decodedCap = 0;
            if(dctx) blosc2_free_ctx(dctx);
            dctx = nullptr;
        }
        largestBlock = 0;
        if(shardIndex.capacity()*sizeof(uint64_t) > keep) std::vector<uint64_t>().swap(shardIndex);
    }
};
ScratchPool<ReadScratch> readScratchPool;

// Chunks of this many bytes or more are decoded in pieces or a part at a time (a test
// build sets it lower, to take those ways with small arrays)
#ifndef CPPZARR_BIG_CHUNK
#define CPPZARR_BIG_CHUNK (1 << 20)
#endif
const uint64_t bigChunk = CPPZARR_BIG_CHUNK;
}

// zarrArr should be initialized to all zeros if you have empty chunks
uint8_t parallelReadZarr(zarr &Zarr, void* zarrArr,
                         const std::vector<uint64_t> &startCoords, 
                         const std::vector<uint64_t> &endCoords,
                         const std::vector<uint64_t> &readShape,
                         const uint64_t bits,
                         const bool useCtx,
                         const bool sparse)
{
    return parallelReadZarr(Zarr, zarrArr, startCoords, endCoords, readShape, bits, useCtx, sparse, false);
}

uint8_t parallelReadZarr(zarr &Zarr, void* zarrArr,
                         const std::vector<uint64_t> &startCoords,
                         const std::vector<uint64_t> &endCoords,
                         const std::vector<uint64_t> &readShape,
                         const uint64_t bits,
                         const bool useCtx,
                         const bool sparse,
                         const bool cOrder)
{
    return parallelReadZarr(Zarr, zarrArr, startCoords, endCoords, readShape, bits, useCtx, sparse, cOrder, nullptr);
}

uint8_t parallelReadZarr(zarr &Zarr, void* zarrArr,
                         const std::vector<uint64_t> &startCoords,
                         const std::vector<uint64_t> &endCoords,
                         const std::vector<uint64_t> &readShape,
                         const uint64_t bits,
                         const bool useCtx,
                         const bool sparse,
                         const bool cOrder,
                         const void* fillMissing)
{
    // Linux's automatic NUMA balancing is kept off the threads that read, for the call (see
    // LocalMemoryPolicy): the calling thread for all of it, and every thread of each
    // parallel region, also those with nothing to do (any thread can mark the process's
    // pages for it). Without huge pages, its hinting faults made reads of big arrays up to
    // 5% slower; kept off only the threads with work, reads of one chunk took up to 1.5
    // times as long.
    const LocalMemoryPolicy memoryPolicy;
    // A 0-dimensional array is a single element: sharding does not apply
    if(Zarr.get_ndims() == 0 && Zarr.get_shard()){
        Zarr.set_errString("Sharding is not supported for 0-dimensional arrays\n");
        return 1;
    }
    const uint64_t bytes = (bits/8);
    const int32_t numWorkers = omp_get_max_threads();

    // Without contexts, blosc decompresses with its global state and as many threads
    if(!useCtx){
        blosc2_init();
        blosc2_set_nthreads(numWorkers);
    }

    // A sharded array is read with its inner chunk shape (the shard shape is
    // restored on return). Each shard's index is at the start or, by default, the
    // end of its file, optionally followed by a CRC32C of it.
    const uint64_t nDims = Zarr.get_ndims();
    std::vector<uint64_t> shardChunks;
    if(Zarr.get_shard()){
        std::vector<uint64_t> innerChunks(nDims);
        for(uint64_t d = 0; d < nDims; d++){
            shardChunks.push_back(Zarr.get_chunks(d));
            innerChunks[d] = Zarr.get_chunk_shape(d);
        }
        Zarr.set_chunks(innerChunks);
    }
    struct RestoreChunks {
        zarr &z;
        const std::vector<uint64_t> chunks;
        ~RestoreChunks(){ if(!chunks.empty()) z.set_chunks(chunks); }
    } restoreChunks{Zarr, shardChunks};
    const uint64_t nInner = Zarr.get_shard() ? Zarr.get_numChunksPerShard() : 0;
    const bool indexAtStart = Zarr.get_shardIndexAtStart();
    const bool indexChecksum = Zarr.get_shardIndexChecksum();
    const uint64_t emptyEntry = std::numeric_limits<uint64_t>::max();
    // How each stored chunk is encoded: a compressor (v2's, or a Zarr v3 codec),
    // and a CRC32C checksum after it (Zarr v3 crc32c codec)
    const std::string compressor = Zarr.get_compressor();
    const bool chunkChecksum = Zarr.get_chunkChecksum();
    if(compressor != "none" && compressor != "blosc" && compressor != "gzip" && compressor != "zstd"){
        Zarr.set_errString("The \""+compressor+"\" compressor is not supported\n");
        return 1;
    }

    uint64_t s = 1;
    for(uint64_t d = 0; d < nDims; d++) s *= Zarr.get_chunks(d);
    const uint64_t sB = s*bytes;

    // Element strides of a decompressed chunk (its stored axes, slowest first: C or
    // F order, or a Zarr v3 transpose) and of the output (F or C order)
    std::vector<uint64_t> chunkStrides(nDims), outStrides(nDims);
    {
        const std::vector<uint64_t> axes = Zarr.get_chunkAxisOrder();
        uint64_t acc = 1;
        for(int64_t i = (int64_t)nDims-1; i >= 0; i--){ chunkStrides[axes[i]] = acc; acc *= Zarr.get_chunks(axes[i]); }
        acc = 1;
        if(cOrder){
            for(int64_t d = (int64_t)nDims-1; d >= 0; d--){ outStrides[d] = acc; acc *= readShape[d]; }
        }
        else{
            for(uint64_t d = 0; d < nDims; d++){ outStrides[d] = acc; acc *= readShape[d]; }
        }
    }

    // Sparse reads skip all-zero chunks, which the output already holds when the
    // fill value is zero
    void* zeroChunkUnc = NULL;
    if(sparse && !fillMissing && fillValueIsZero(Zarr.get_fill_value(), Zarr.get_dtype())){
        zeroChunkUnc = calloc(s,bytes);
    }

    std::atomic<int> err{0};
    std::string errString;
    auto fail = [&](const std::string &message){
        #pragma omp critical
        {
            err = 1;
            errString = message;
        }
    };

    // Chunks are handed out in the output's memory order (first axis fastest for an
    // F-order output), so each thread's block of chunks fills a part of the output of
    // its own. The chunk list runs last axis fastest: for an F-order output the chunks
    // that share each page of it were far apart in the list, and with some thread
    // counts several threads started on the same rows and raced to fault in the same
    // new pages (up to 3x slower, with most huge pages lost). The list is a grid of
    // chunks, or of shards that each hold a grid of inner chunks (shard by shard),
    // both last axis fastest.
    const int64_t listed = (int64_t)Zarr.get_numChunks();
    std::vector<uint64_t> outerFirst(nDims), outerCount(nDims), outerStride(nDims), innerCount(nDims, 1), innerStride(nDims, 1);
    uint64_t nOuter = 1, nPer = 1;
    for(uint64_t d = 0; d < nDims; d++){
        if(Zarr.get_shard()) innerCount[d] = (shardChunks[d]+Zarr.get_chunks(d)-1)/Zarr.get_chunks(d);
        const uint64_t extent = innerCount[d]*Zarr.get_chunks(d);
        const uint64_t first = startCoords[d]/extent, last = endCoords[d]/extent+(endCoords[d]%extent ? 1 : 0);
        outerFirst[d] = first;
        outerCount[d] = last > first ? last-first : 0;
        nOuter *= outerCount[d];
        nPer *= innerCount[d];
    }
    for(int64_t d = (int64_t)nDims-1; d >= 0; d--){
        outerStride[d] = d+1 < (int64_t)nDims ? outerStride[d+1]*outerCount[d+1] : 1;
        innerStride[d] = d+1 < (int64_t)nDims ? innerStride[d+1]*innerCount[d+1] : 1;
    }
    const bool regular = nOuter*nPer == (uint64_t)listed;
    bool reorder = !cOrder && nDims > 1 && regular;
    // The list position of the chunk at position q of the output's order
    auto listIndex = [&](uint64_t q) -> int64_t {
        if(!reorder) return (int64_t)q;
        uint64_t o = q/nPer, c = q%nPer, f = 0, g = 0;
        for(uint64_t d = 0; d < nDims; d++){
            f += (o%outerCount[d])*outerStride[d];
            o /= outerCount[d];
            g += (c%innerCount[d])*innerStride[d];
            c /= innerCount[d];
        }
        return (int64_t)(f*nPer+g);
    };

    // The chunks to read: the whole list, or for a sharded array (whose list has every
    // inner chunk of the shards the region touches) the inner chunks the region
    // touches. (Inner chunks just past a region that ends on their boundary used to be
    // read and decoded too.) An inner chunk's indices follow from its place in the list,
    // without parsing its name; a parallel pass over the names made small reads wait
    // for any thread the system had paused.
    std::vector<int64_t> items;
    if(Zarr.get_shard()){
        for(int64_t q = 0; q < listed && !err; q++){
            const int64_t f = listIndex(q);
            bool in = true;
            if(regular){
                uint64_t o = (uint64_t)f/nPer, c = (uint64_t)f%nPer;
                for(int64_t d = (int64_t)nDims-1; d >= 0 && in; d--){
                    const uint64_t idx = (outerFirst[d]+o%outerCount[d])*innerCount[d]+c%innerCount[d];
                    o /= outerCount[d];
                    c /= innerCount[d];
                    const uint64_t cs = Zarr.get_chunk_shape(d);
                    in = idx*cs < endCoords[d] && (idx+1)*cs > startCoords[d];
                }
            }
            else{
                try{
                    const std::vector<uint64_t> cAV = Zarr.get_chunkAxisVals(Zarr.get_chunkNames(f));
                    in = cAV.size() == nDims;
                    for(uint64_t d = 0; d < nDims && in; d++){
                        const uint64_t cs = Zarr.get_chunk_shape(d);
                        in = cAV[d]*cs < endCoords[d] && (cAV[d]+1)*cs > startCoords[d];
                    }
                }
                catch(const std::exception &e){
                    fail(std::string("Read error: ")+e.what()+"\n");
                }
            }
            if(in) items.push_back(f);
        }
    }
    const int64_t nItems = Zarr.get_shard() ? (int64_t)items.size() : listed;

    // A chunk's slowest stored slices (along its slowest axis more than one element
    // thick): a blosc chunk is decoded only for the slices inside the region, a range
    // of its bytes
    uint64_t sliceAxis = 0, sliceBytes = sB;
    {
        const std::vector<uint64_t> axes = Zarr.get_chunkAxisOrder();
        for(uint64_t i = 0; i < nDims; i++){
            if(Zarr.get_chunks(axes[i]) > 1){
                sliceAxis = axes[i];
                sliceBytes = chunkStrides[sliceAxis]*bytes;
                break;
            }
        }
    }
    // The bytes [from, to) of chunk cAV's decoded form that hold its slices inside the
    // region
    auto bytesInRegion = [&](const std::vector<uint64_t> &cAV, uint64_t &from, uint64_t &to){
        from = 0;
        to = sB;
        if(sliceBytes == sB) return;
        const uint64_t c = Zarr.get_chunks(sliceAxis), b = cAV[sliceAxis]*c;
        from = (std::max(b, startCoords[sliceAxis])-b)*sliceBytes;
        to = (std::min(b+c, endCoords[sliceAxis])-b)*sliceBytes;
    };

    // Reading ahead (Linux): the kernel is asked to read parts of the files into the page
    // cache in the background, so that more reads are in flight and they are bigger. On a
    // local disk (a file system on a block device, whose device number has a major number;
    // NFS, tmpfs and FUSE have 0) a shard's parts are; and on NVMe drives also each thread's
    // next chunk file (on a SATA SSD that made cold reads up to 2.5% slower). NFS clients
    // read ahead well themselves: reading ahead whole chunk files or big parts of shards
    // there made cold reads up to 1.5 times as long. (Found out when a read first needs to
    // know: stat can ask an NFS server, and statfs always does, which made small reads up
    // to 65% slower.) 0: not a local disk, 1: a local disk, 2: an NVMe drive
    // (Reading ahead is asked for in steps of the most a call reads, the larger of the
    // device's largest request and its read-ahead, up to 2 MB: 1.25 MB on the cluster's
    // disks, where 512 KB steps made cold reads of a whole sharded array on a SATA SSD 3%
    // slower, and 128 KB on many NVMe drives. On NFS, where a call reads 1 MB or more, 512 KB
    // steps were faster than 1 MB ones.)
    std::atomic<int> diskKind{-1};
    std::atomic<uint64_t> readAheadStep{512 << 10};
    auto disk = [&]() -> int {
        int v = diskKind.load(std::memory_order_acquire);
        if(v < 0){
            v = 0;
            #ifdef __linux__
            struct stat st;
            if(stat(Zarr.get_fileName().c_str(), &st) == 0 && major(st.st_dev) != 0){
                const BlockDevice b = blockDevice(st.st_dev);
                v = b.nvme ? 2 : 1;
                if(b.readAheadMax) readAheadStep.store(std::max<uint64_t>(4096, std::min<uint64_t>(2 << 20, b.readAheadMax)), std::memory_order_relaxed);
            }
            #endif
            diskKind.store(v, std::memory_order_release);
        }
        return v;
    };
    // Ask for the byte ranges [first, second) of an open file
    auto readAhead = [&](const int fd, const std::vector<std::pair<uint64_t, uint64_t>> &ranges){
    #ifdef __linux__
        const uint64_t step = readAheadStep.load(std::memory_order_relaxed);
        for(const auto &r : ranges){
            for(uint64_t at = r.first; at < r.second; at += step){
                posix_fadvise(fd, (off_t)at, (off_t)std::min(step, r.second-at), POSIX_FADV_WILLNEED);
            }
        }
    #else
        (void)fd;
        (void)ranges;
    #endif
    };
    // Ask for a whole file
    auto readAheadFile = [&](const std::string &fileName){
    #ifdef __linux__
        const int fd = open(fileName.c_str(), O_RDONLY | O_CLOEXEC);
        if(fd < 0) return;
        struct stat st;
        if(fstat(fd, &st) == 0 && st.st_size > 0) readAhead(fd, {{0, (uint64_t)st.st_size}});
        close(fd);
    #else
        (void)fileName;
    #endif
    };
    // Whether the page cache holds the start of chunk files across the list (four of them),
    // asked with one-byte reads that may not wait for the disk; true where that cannot be
    // asked
    auto filesCached = [&]() -> bool {
    #if defined(__linux__) && defined(SYS_preadv2)
        try{
            for(int64_t k = 0; k < 4; k++){
                const std::string &name = Zarr.get_chunkNames((listed-1)*k/3);
                const std::string fileName = Zarr.get_fileName()+"/"+Zarr.get_subfoldersString(Zarr.get_chunkAxisVals(name))+"/"+
                                             Zarr.chunkKey(name);
                const int fd = open(fileName.c_str(), O_RDONLY | O_CLOEXEC);
                if(fd < 0) continue;
                uint8_t b = 0;
                struct iovec v = {&b, 1};
                const long n = syscall(SYS_preadv2, fd, &v, 1, 0L, 0L, 8 /* RWF_NOWAIT */);
                const int e = errno;
                close(fd);
                if(n < 0 && e == EAGAIN) return false;
            }
        }
        catch(const std::exception&){}
    #endif
        return true;
    };
    // On a local disk other than an NVMe drive, chunk files the page cache does not hold are
    // handed out in the list's order instead (the order they were written in, so mostly the
    // disk's), so that each thread reads files that lie together: in the output's order the
    // threads jumped across the disk, and cold F-order reads of the 16,384 files of an array
    // of 64^3 chunks from a SATA SSD took 1.2-1.4 times as long (as long as in v1.6.0).
    // Only reads of at least 8 chunk files a thread, which count as cold if any of four
    // files across the list is not in the page cache.
    if(reorder && !Zarr.get_shard() && listed >= 8*(int64_t)numWorkers && disk() == 1 && !filesCached()) reorder = false;
    // Shards (by their place in the region's grid of shards) a thread has read ahead in this call
    std::unique_ptr<std::atomic<uint8_t>[]> shardOpened;
    if(Zarr.get_shard() && regular && nOuter > 0){
        shardOpened.reset(new std::atomic<uint8_t>[nOuter]);
        for(uint64_t o = 0; o < nOuter; o++) shardOpened[o].store(0, std::memory_order_relaxed);
    }

    // A thread's shard being read: kept open, with its index loaded
    struct ShardState {
        ShardFile file;
        std::string path;
        bool open = false, indexRead = false;
        uint64_t* index = nullptr;
    };
    // A thread that opens the shard asks to read ahead the inner chunks of it the
    // region needs (once for each shard), neighbors (less than 256 KB apart)
    // together: each read their own as threads came to them, the inner chunks of a
    // shard read whole from a cold SATA SSD 2.1 times as long as with this, and a
    // timepoint of a few inner chunks of each of four shards up to 1.4 times as
    // long from NFS. Shards that need at most 8 MB only when the thread's index
    // read took long, so came from the disk or the server (from the page cache it
    // takes under 25 us, from VAST over 700 us and from an NVMe drive about 100):
    // asked for data already in the page cache, NFS made such reads up to 19%
    // slower. (A thread whose index came from the page cache, which another thread
    // may just have read it into, leaves them to one whose index did not.) Bigger
    // ones only on local disks (see above).
    // (Not inlined: inlined in readChunk, it made one-chunk reads at 8 threads 6% slower, though
    // they never came here.)
    auto readAheadShard = [&](const ShardFile &file, const uint64_t* index, const std::vector<uint64_t> &cAV,
                              const double indexTook, const uint64_t indexBytes) CPPZARR_NOINLINE {
        uint64_t o = 0;
        for(uint64_t d = 0; d < nDims; d++) o += (cAV[d]/innerCount[d]-outerFirst[d])*outerStride[d];
        const int kind = o < nOuter && shardOpened[o].load(std::memory_order_relaxed) == 0 ? disk() : -1;
        const bool slow = indexTook > (kind > 0 ? 40e-6 : 300e-6)+indexBytes/2e9;
        if(kind >= 0 && (slow || kind > 0)){
            // (the shard's inner chunks inside the region, last axis fastest)
            std::vector<uint64_t> lo(nDims), hi(nDims), at(nDims);
            bool any = true;
            for(uint64_t d = 0; d < nDims; d++){
                const uint64_t cs = Zarr.get_chunk_shape(d), first = (cAV[d]/innerCount[d])*innerCount[d];
                lo[d] = std::max(first, startCoords[d]/cs);
                hi[d] = std::min(first+innerCount[d], (endCoords[d]+cs-1)/cs);
                any = any && hi[d] > lo[d];
                at[d] = lo[d];
            }
            std::vector<std::pair<uint64_t, uint64_t>> ranges;
            uint64_t total = 0;
            while(any){
                const uint64_t p = Zarr.get_chunkShardPosition(at);
                const uint64_t off = index[2*p], n = index[(2*p)+1];
                if(!(off == emptyEntry && n == emptyEntry) && n > 0){
                    ranges.push_back({off, off+n});
                    total += n;
                }
                int64_t d = (int64_t)nDims-1;
                for(; d >= 0; d--){
                    if(++at[d] < hi[d]) break;
                    at[d] = lo[d];
                }
                any = d >= 0;
            }
            uint8_t opened = 0;
            if(ranges.size() > 1 && (total <= (8 << 20) ? slow : kind > 0) && shardOpened[o].compare_exchange_strong(opened, 1)){
                std::sort(ranges.begin(), ranges.end());
                std::vector<std::pair<uint64_t, uint64_t>> runs;
                for(const auto &r : ranges){
                    if(!runs.empty() && r.first <= runs.back().second+(256 << 10)) runs.back().second = std::max(runs.back().second, r.second);
                    else runs.push_back(r);
                }
                readAhead(file.handle(), runs);
            }
        }
    };
    // Read a chunk's stored bytes (the chunk named name, at indices cAV) into the buffer
    // alloc(n) gives, and check the checksum after them if it has one. Returns their
    // length (without the checksum), -1 for a chunk that was never written (it reads as
    // the fill value, which the output already holds), or -2 after an error.
    auto readChunk = [&](ShardState &ss, const std::string &name, const std::vector<uint64_t> &cAV,
                         const std::string &subfolderName, const std::function<uint8_t*(uint64_t)> &alloc,
                         uint8_t* &data) -> int64_t {
        std::string fileName;
        int64_t fileLen = 0;
        if(!Zarr.get_shard()){
            fileName = Zarr.get_fileName()+"/"+subfolderName+"/"+Zarr.chunkKey(name);
            std::ifstream file(fileName, std::ios::binary);
            if(!file.is_open()) return -1;
            file.seekg(0, std::ios::end);
            fileLen = (int64_t)file.tellg();
            data = alloc((uint64_t)fileLen);
            file.seekg(0, std::ios::beg);
            file.read(reinterpret_cast<char*>(data), fileLen);
        }
        else{
            fileName = Zarr.get_fileName()+"/"+subfolderName+"/"+Zarr.chunkKey(Zarr.chunkNameToShardName(name));
            if(fileName != ss.path){
                ss.open = ss.file.open(fileName);
                ss.path = fileName;
                ss.indexRead = false;
            }
            // A shard that was never written reads as the fill value
            if(!ss.open) return -1;
            // Its index, checked: its checksum, and every entry inside the file
            if(!ss.indexRead){
                const auto indexStart = std::chrono::steady_clock::now();
                const int64_t size = ss.file.size();
                const uint64_t shardLen = size > 0 ? (uint64_t)size : 0;
                const uint64_t entryBytes = 2*nInner*sizeof(uint64_t);
                const uint64_t indexBytes = entryBytes+(indexChecksum ? 4 : 0);
                std::string problem;
                if(shardLen < indexBytes) problem = "is smaller than its index";
                else{
                    // (the entries and their checksum in one read: the buffer has room for both)
                    uint32_t crc = 0;
                    const bool got = ss.file.readAt(indexAtStart ? 0 : shardLen-indexBytes, ss.index, indexBytes);
                    if(got && indexChecksum) memcpy(&crc, reinterpret_cast<const uint8_t*>(ss.index)+entryBytes, sizeof(crc));
                    if(!got) problem = "could not be read";
                    else if(indexChecksum && crc != crc32c(reinterpret_cast<const uint8_t*>(ss.index), entryBytes)){
                        problem = "has an index that does not match its checksum";
                    }
                    else{
                        // chunk data lies after the index when it is at the start, else before it
                        const uint64_t lo = indexAtStart ? indexBytes : 0;
                        const uint64_t hi = indexAtStart ? shardLen : shardLen-indexBytes;
                        for(uint64_t i = 0; i < nInner && problem.empty(); i++){
                            const uint64_t off = ss.index[2*i], n = ss.index[(2*i)+1];
                            if(off == emptyEntry && n == emptyEntry) continue;
                            if(off < lo || off > hi || n > hi-off) problem = "has an index entry outside the file";
                        }
                    }
                }
                if(!problem.empty()){
                    ss.path.clear();
                    fail("The shard "+fileName+" "+problem+" (it may be damaged, or not match the array's metadata)\n");
                    return -2;
                }
                ss.indexRead = true;
                const double indexTook = std::chrono::duration<double>(std::chrono::steady_clock::now()-indexStart).count();
                if(shardOpened) readAheadShard(ss.file, ss.index, cAV, indexTook, indexBytes);
            }
            const uint64_t pos = Zarr.get_chunkShardPosition(cAV);
            const uint64_t offset = ss.index[2*pos], nbytes = ss.index[(2*pos)+1];
            // An inner chunk that was never written reads as the fill value
            if(offset == emptyEntry && nbytes == emptyEntry) return -1;
            fileLen = (int64_t)nbytes;
            data = alloc(nbytes);
            if(!ss.file.readAt(offset, data, nbytes)){
                ss.path.clear();
                fail("The shard "+fileName+" could not be read (it may be damaged, or not match the array's metadata)\n");
                return -2;
            }
        }
        // A crc32c codec's checksum follows the chunk's bytes
        if(chunkChecksum){
            uint32_t storedCrc = 0;
            if(fileLen >= 4) memcpy(&storedCrc, data+fileLen-4, 4);
            if(fileLen < 4 || storedCrc != crc32c(data, (size_t)fileLen-4)){
                fail("The chunk "+fileName+" does not match its checksum (it may be damaged)\n");
                return -2;
            }
            fileLen -= 4;
        }
        if(compressor == "none" && fileLen != (int64_t)sB){
            fail("Uncompressed chunk is "+std::to_string(fileLen)+" bytes instead of "+std::to_string(sB)+
                 ". ChunkName: "+Zarr.get_fileName()+"/"+subfolderName+"/"+name+"\n");
            return -2;
        }
        return fileLen;
    };
    // Decode the bytes [from, to) of a chunk (len stored bytes at src) into the same
    // place of dst, which has room for all of it: blosc decodes only the blocks they
    // are in (with contexts; straight into dst, with a mask of the blocks to leave out),
    // the other codecs all of the chunk. Reports an error and returns false if it cannot
    // be decoded.
    auto decodeChunkRange = [&](ReadScratch &sc, blosc2_context* dctx, const uint8_t* src, const int64_t len,
                                uint8_t* dst, const uint64_t from, const uint64_t to,
                                const std::string &subfolderName, const std::string &name) -> bool {
        bool decodedOk = true;
        int uncErr = 0;
        if(compressor == "zstd"){
            const size_t n = ZSTD_decompressDCtx(sc.zctx, dst, sB, src, (size_t)len);
            decodedOk = !ZSTD_isError(n) && n == sB;
        }
        else if(compressor == "blosc"){
            bool done = false;
            if(useCtx && (from > 0 || to < sB) && len >= BLOSC_MIN_HEADER_LENGTH){
                int32_t nbytes = 0, blocksize = 0;
                if(blosc2_cbuffer_sizes(src, &nbytes, nullptr, &blocksize) >= 0 && (uint64_t)nbytes == sB && blocksize > 0){
                    const uint64_t bs = (uint64_t)blocksize, nBlocks = (sB+bs-1)/bs;
                    const uint64_t first = from/bs, last = (to+bs-1)/bs;
                    if(last-first < nBlocks){
                        if(sc.skipCap < nBlocks){
                            sc.skip.reset(new bool[nBlocks]);
                            sc.skipCap = nBlocks;
                        }
                        for(uint64_t k = 0; k < nBlocks; k++) sc.skip[k] = k < first || k >= last;
                        // (a chunk the mask does not fit is decoded whole below)
                        if(blosc2_set_maskout(dctx, sc.skip.get(), (int)nBlocks) == 0){
                            done = blosc2_decompress_ctx(dctx, src, (int32_t)len, dst, (int32_t)sB) == (int)sB;
                        }
                    }
                }
            }
            if(!done){
                const int n = useCtx ? blosc2_decompress_ctx(dctx, src, (int32_t)len, dst, (int32_t)sB)
                                     : blosc2_decompress(src, (int32_t)len, dst, (int32_t)sB);
                decodedOk = n >= 0;
            }
        }
        else if(compressor == "gzip"){
            // gzip or zlib streams (header detected), possibly several gzip members
            uint64_t used = 0, made = 0;
            decodedOk = len > 0;
            while(decodedOk && used < (uint64_t)len){
                z_stream stream;
                memset(&stream, 0, sizeof(stream));
                uncErr = inflateInit2(&stream, 32);
                if(uncErr != Z_OK){
                    decodedOk = false;
                    break;
                }
                stream.next_in = const_cast<uint8_t*>(src)+used;
                stream.avail_in = (uInt)(len-used);
                stream.next_out = dst+made;
                stream.avail_out = (uInt)(sB-made);
                uncErr = inflate(&stream, Z_NO_FLUSH);
                used = (uint64_t)len-stream.avail_in;
                made = sB-stream.avail_out;
                inflateEnd(&stream);
                decodedOk = uncErr == Z_STREAM_END;
            }
        }
        if(!decodedOk){
            fail("Decompression error. Error code: "+std::to_string(uncErr)+" ChunkName: "+
                 Zarr.get_fileName()+"/"+subfolderName+"/"+name+"\n");
        }
        return decodedOk;
    };
    // Copy the part of a decoded chunk inside the read region (and inside [lo, hi) of
    // the output's axis sAxis, when given) into the output: contiguous runs for chunks
    // in the output's order, a tiled transpose otherwise. chunkData holds the chunk's
    // decoded bytes from byte offset on.
    auto copyChunkAt = [&](const uint8_t* chunkData, const uint64_t offset, const std::vector<uint64_t> &cAV,
                           std::vector<uint64_t> &boxExt, const int64_t sAxis, const uint64_t sLo, const uint64_t sHi){
        uint64_t srcOff = 0, dstOff = 0;
        for(uint64_t d = 0; d < nDims; d++){
            const uint64_t b = cAV[d]*Zarr.get_chunks(d);
            uint64_t lo = std::max(b, startCoords[d]);
            uint64_t hi = std::min(b+Zarr.get_chunks(d), endCoords[d]);
            if((int64_t)d == sAxis){
                lo = std::max(lo, startCoords[d]+sLo);
                hi = std::min(hi, startCoords[d]+sHi);
            }
            if(hi <= lo) return;
            boxExt[d] = hi-lo;
            srcOff += (lo-b)*chunkStrides[d];
            dstOff += (lo-startCoords[d])*outStrides[d];
        }
        copyBoxND(bytes, chunkData+(srcOff*bytes-offset), (uint8_t*)zarrArr+dstOff*bytes,
                  boxExt, chunkStrides, outStrides, true);
    };
    auto copyChunk = [&](const uint8_t* chunkData, const std::vector<uint64_t> &cAV, std::vector<uint64_t> &boxExt,
                         const int64_t sAxis, const uint64_t sLo, const uint64_t sHi){
        copyChunkAt(chunkData, 0, cAV, boxExt, sAxis, sLo, sHi);
    };
    // When the output does not hold the fill value already (fillMissing), the chunks
    // found missing are noted, and their part of the output set to it after the read
    std::mutex missingMutex;
    std::vector<std::vector<uint64_t>> missing;
    auto noteMissing = [&](const std::vector<uint64_t> &cAV){
        if(!fillMissing) return;
        std::lock_guard<std::mutex> g(missingMutex);
        missing.push_back(cAV);
    };
    auto fillChunk = [&](const std::vector<uint64_t> &cAV, std::vector<uint64_t> &boxExt){
        uint64_t dstOff = 0;
        for(uint64_t d = 0; d < nDims; d++){
            const uint64_t b = cAV[d]*Zarr.get_chunks(d);
            const uint64_t lo = std::max(b, startCoords[d]);
            const uint64_t hi = std::min(b+Zarr.get_chunks(d), endCoords[d]);
            if(hi <= lo) return;
            boxExt[d] = hi-lo;
            dstOff += (lo-startCoords[d])*outStrides[d];
        }
        fillBoxND(bytes, (uint8_t*)zarrArr+dstOff*bytes, boxExt, outStrides, fillMissing);
    };
    // A chunk's stored bytes: at most its decoded size and a codec's overhead (when that
    // is small enough for the pool to keep)
    const uint64_t storedRoom = sB+sB/128+4096 <= (1 << 20) ? sB+sB/128+4096 : 0;
    // A thread's buffers and decoders from the pool: a one-thread blosc decoder (with
    // contexts) and a zstd decoder as the array needs
    auto takeScratch = [&]() -> std::unique_ptr<ReadScratch> {
        std::unique_ptr<ReadScratch> sc = readScratchPool.take(omp_get_thread_num());
        sc->shardIndex.resize(2*nInner+1);
        if(compressor == "blosc" && useCtx && !sc->dctx) sc->dctx = blosc2_create_dctx(BLOSC2_DPARAMS_DEFAULTS);
        if(compressor == "zstd" && !sc->zctx) sc->zctx = ZSTD_createDCtx();
        return sc;
    };

    // The output's slowest axis (more than one element thick)
    int64_t outSlow = 0;
    if(cOrder){
        for(int64_t d = 0; d < (int64_t)nDims; d++) if(readShape[d] > 1){ outSlow = d; break; }
    }
    else{
        for(int64_t d = (int64_t)nDims-1; d >= 0; d--) if(readShape[d] > 1){ outSlow = d; break; }
    }
    // Big blosc chunks whose slices run along the output's slowest axis are decoded a part
    // (a few blocks, whole slices) at a time into a small buffer, and each part copied out
    // while it is in the cache. Decoded whole, they went through memory twice more, into
    // buffers of their size that every call faulted in again, which made reads take up to
    // 1.9 times as long (reads of a few chunks up to 3 times). (A chunk whose slices run
    // along a faster axis of the output is copied with a tiled transpose, which a part at
    // a time would spread over every page of its box.)
    const bool fusable = compressor == "blosc" && useCtx && sB >= bigChunk && nDims > 0 &&
                         sliceBytes < sB && (int64_t)sliceAxis == outSlow;
    // The unit a fusable chunk (len stored bytes at src) is decoded in: whole blocks that
    // hold whole slices, at most 1 MB, and blosc's item size, which ranges are counted in.
    // Returns 0 for a chunk that is to be decoded whole.
    auto fusedUnit = [&](const uint8_t* src, const int64_t len, uint64_t &typesize) -> uint64_t {
        int32_t nbytes = 0, blocksize = 0;
        size_t ts = 0;
        int flags = 0;
        if(len < BLOSC_MIN_HEADER_LENGTH || blosc2_cbuffer_sizes(src, &nbytes, nullptr, &blocksize) < 0 ||
           (uint64_t)nbytes != sB || blocksize <= 0) return 0;
        blosc1_cbuffer_metainfo(src, &ts, &flags);
        const uint64_t bs = (uint64_t)blocksize;
        uint64_t a = bs, b = sliceBytes;
        while(b){
            const uint64_t t = a%b;
            a = b;
            b = t;
        }
        const uint64_t unit = bs/a*sliceBytes;
        if(ts == 0 || bs%ts != 0 || sB%ts != 0 || unit > (1 << 20)) return 0;
        typesize = ts;
        return unit;
    };
    // Decode a fusable chunk's slices [s0, s1) a part at a time (about 256 KB of whole
    // units) and copy each part into the output; a sparse read leaves out parts that are
    // all zeros. The parts start on blocks: blosc (2.23 to at least 3.3.6) returns the
    // next block for a range that starts inside a block and is one block long. With spread,
    // the first part depends on it, so that threads reading chunks that share output pages
    // fault in different pages at a time (from the same part on, they raced to fault them
    // in, which made reads take up to twice as long).
    auto copyFused = [&](ReadScratch &sc, const uint8_t* src, const int64_t len, const uint64_t unit, const uint64_t typesize,
                         const std::vector<uint64_t> &cAV, std::vector<uint64_t> &boxExt, const uint64_t s0, const uint64_t s1,
                         const uint64_t spread, const std::string &subfolderName, const std::string &name) -> bool {
        const uint64_t partBytes = std::max<uint64_t>(1, (256 << 10)/unit)*unit;
        uint8_t* buf = sc.partBuffer(partBytes);
        int32_t nbytes = 0, blocksize = 0;
        blosc2_cbuffer_sizes(src, &nbytes, nullptr, &blocksize);
        sc.largestBlock = std::max<uint64_t>(sc.largestBlock, (uint64_t)blocksize);
        const uint64_t ra = s0*sliceBytes, rb = s1*sliceBytes;
        const uint64_t pa = (ra/unit)*unit, pb = std::min(sB, ((rb+unit-1)/unit)*unit);
        const uint64_t nParts = (pb-pa+partBytes-1)/partBytes;
        const uint64_t first = (((spread*0x9E3779B97F4A7C15ULL) >> 32)*nParts) >> 32;
        const uint64_t b = cAV[sliceAxis]*Zarr.get_chunks(sliceAxis);
        for(uint64_t q = 0; q < nParts; q++){
            const uint64_t ga = pa+((first+q)%nParts)*partBytes, ge = std::min(pb, ga+partBytes);
            const int n = blosc2_getitem_ctx(sc.dctx, src, (int32_t)len, (int)(ga/typesize), (int)((ge-ga)/typesize), buf, (int32_t)(ge-ga));
            if(n != (int)(ge-ga)){
                fail("Decompression error. Error code: "+std::to_string(n)+" ChunkName: "+
                     Zarr.get_fileName()+"/"+subfolderName+"/"+name+"\n");
                return false;
            }
            const uint64_t za = std::max(ga, ra), zb = std::min(ge, rb);
            if(zeroChunkUnc && !memcmp(zeroChunkUnc, buf+(za-ga), zb-za)) continue;
            copyChunkAt(buf, ga, cAV, boxExt, sliceAxis, b+za/sliceBytes-startCoords[sliceAxis], b+zb/sliceBytes-startCoords[sliceAxis]);
        }
        return true;
    };

    // Big zstd chunks (the zstd codec itself) whose slices run along the output's slowest
    // axis are decoded a part (whole slices, about 256 KB) at a time with zstd's streaming
    // decoder, up to the last slice the region needs, and each part copied out while it is
    // in the cache. Decoded whole, they went through memory twice more, into buffers of
    // their size that every call faulted in again, which made region reads up to 1.2 times
    // as long, and the buffers held them to fewer threads. Only a chunk of one frame (all of
    // its bytes) that says it decodes to a chunk's bytes, as cpp-zarr and other Zarr writers
    // write them; one that does not say is decoded to its end, and others whole
    // (which reports any that do not decode to a chunk's bytes, as before). A chunk decoded
    // to its end is checked as one decoded whole is: its size, and its checksum if it has
    // one; past a region's last slice it is not decoded.
    const bool zstdFusable = compressor == "zstd" && sB >= bigChunk && nDims > 0 && sliceBytes < sB && (int64_t)sliceAxis == outSlow;
    // How a chunk (len stored bytes at src) is decoded: 0 whole, 1 a part at a time up to the
    // region's last slice, 2 a part at a time to its end
    auto zstdStreamed = [&](const uint8_t* src, const int64_t len) -> int {
        if(ZSTD_findFrameCompressedSize(src, (size_t)len) != (size_t)len) return 0;
        const unsigned long long n = ZSTD_getFrameContentSize(src, (size_t)len);
        return n == sB ? 1 : n == ZSTD_CONTENTSIZE_UNKNOWN ? 2 : 0;
    };
    auto copyZstdStream = [&](ReadScratch &sc, const uint8_t* src, const int64_t len, const bool toEnd,
                              const std::vector<uint64_t> &cAV, std::vector<uint64_t> &boxExt, const uint64_t s0, const uint64_t s1,
                              const std::string &subfolderName, const std::string &name) -> bool {
        auto bad = [&](const std::string &why){
            fail("Decompression error ("+why+"). ChunkName: "+Zarr.get_fileName()+"/"+subfolderName+"/"+name+"\n");
            return false;
        };
        const uint64_t partBytes = std::max<uint64_t>(1, (256 << 10)/sliceBytes)*sliceBytes;
        uint8_t* buf = sc.partBuffer(partBytes);
        ZSTD_DCtx_reset(sc.zctx, ZSTD_reset_session_only);
        ZSTD_inBuffer in = {src, (size_t)len, 0};
        const uint64_t ra = s0*sliceBytes, rb = s1*sliceBytes, until = toEnd ? sB : rb;
        const uint64_t b = cAV[sliceAxis]*Zarr.get_chunks(sliceAxis);
        size_t r = 1;
        for(uint64_t ga = 0; ga < until; ga += partBytes){
            const uint64_t ge = std::min(sB, ga+partBytes);
            ZSTD_outBuffer out = {buf, (size_t)(ge-ga), 0};
            while(out.pos < out.size){
                const size_t inBefore = in.pos, outBefore = out.pos;
                r = ZSTD_decompressStream(sc.zctx, &out, &in);
                if(ZSTD_isError(r)) return bad(ZSTD_getErrorName(r));
                if(in.pos == inBefore && out.pos == outBefore) return bad("the chunk ends early");
            }
            const uint64_t za = std::max(ga, ra), zb = std::min(ge, rb);
            if(zb <= za) continue;
            if(zeroChunkUnc && !memcmp(zeroChunkUnc, buf+(za-ga), zb-za)) continue;
            copyChunkAt(buf, ga, cAV, boxExt, sliceAxis, b+za/sliceBytes-startCoords[sliceAxis], b+zb/sliceBytes-startCoords[sliceAxis]);
        }
        // Decoded to its end, the frame must end there (and its checksum match)
        if(until == sB){
            uint8_t more = 0;
            while(r != 0){
                ZSTD_outBuffer out = {&more, 1, 0};
                const size_t inBefore = in.pos;
                r = ZSTD_decompressStream(sc.zctx, &out, &in);
                if(ZSTD_isError(r)) return bad(ZSTD_getErrorName(r));
                if(out.pos > 0) return bad("the chunk decodes to more than a chunk's bytes");
                if(r != 0 && in.pos == inBefore) return bad("the chunk ends early");
            }
        }
        return true;
    };

    // Chunks of 1 MB or more that do not spread evenly over the threads (with contexts):
    // fewer chunks than threads, or so few more that handing out whole chunks would leave
    // threads idle for more than a fifth of the call (say 65 chunks for 48 threads, where
    // a few threads read two while the rest wait). Then every thread takes part in all of
    // the work, a group of chunks at a time. Each chunk is read by one thread, and the
    // output is copied a slab at a time (along its slowest axis), each slab from every
    // chunk in it, so that no two threads write the same pages. Fusable chunks are decoded
    // by the copies, each the parts that hold its slab; others are first decoded in pieces
    // (a blosc chunk a range of its blocks at a time). (One blosc decoder of several
    // threads per chunk started its threads in every call, which made reads of a few
    // chunks up to 27% slower, and threads that copied whole chunks raced to fault in the
    // output pages that neighboring chunks share.) Chunks that only decode whole (other
    // codecs than blosc), one thread each, take this way only when at most a quarter of
    // the threads would have one: they share out just their copies, which with more of
    // them took longer than they saved (16 zstd chunks on 48 threads 13% slower without
    // huge pages). Otherwise a chunk at a time per thread overlaps reading, decoding and
    // copying across the threads: with more chunks, and with small ones, whose uneven
    // share costs less than the threads' waiting for each other's chunks.
    const int64_t perThread = nItems > 0 ? (nItems+numWorkers-1)/numWorkers : 0;
    if(useCtx && nItems > 0 && sB >= bigChunk &&
       (compressor == "blosc" ? nItems < numWorkers || 4*perThread*(int64_t)numWorkers > 5*nItems : 4*nItems <= numWorkers)){
        // The output's slabs: whole numbers of 2 MB (huge pages) of it where possible,
        // about two per thread
        const int64_t sAxis = outSlow;
        const uint64_t extent = nDims ? readShape[sAxis] : 1;
        const uint64_t planeBytes = (nDims ? outStrides[sAxis] : 1)*bytes;
        const uint64_t planesPerPage = std::max<uint64_t>(1, ((2 << 20)+planeBytes-1)/planeBytes);
        uint64_t slab = std::max<uint64_t>(planesPerPage, (extent+2*numWorkers-1)/(2*numWorkers));
        slab = ((slab+planesPerPage-1)/planesPerPage)*planesPerPage;
        const uint64_t nSlabs = (extent+slab-1)/slab;

        // Each chunk's place in the list and its indices, in the output's order
        std::vector<int64_t> chunkF(nItems);
        std::vector<std::vector<uint64_t>> chunkAV(nItems);
        for(int64_t i = 0; i < nItems && !err; i++){
            try{
                chunkF[i] = Zarr.get_shard() ? items[i] : listIndex(i);
                chunkAV[i] = Zarr.get_chunkAxisVals(Zarr.get_chunkNames(chunkF[i]));
            }
            catch(const std::exception &e){
                fail(std::string("Read error: ")+e.what()+"\n");
            }
        }

        // A group's chunks are in slots: buffers from the pool (small chunks' buffers are
        // kept between calls), as many as the chunks of the fewest rounds of one chunk per
        // thread (33 for 65 chunks on 48 threads): more would take as long, and fault in
        // more fresh memory, which made reads up to 30% slower
        const int64_t group = (nItems+perThread-1)/perThread;
        const int64_t nGroups = (nItems+group-1)/group;
        std::vector<std::unique_ptr<ReadScratch>> slots(group);
        try{
            for(auto &sl : slots) sl = readScratchPool.take(-1);
        }
        catch(const std::exception &e){
            fail(std::string("Read error: ")+e.what()+"\n");
        }
        std::vector<uint8_t*> storedOf(nItems, nullptr), decodedOf(nItems, nullptr);
        std::vector<int64_t> storedLen(nItems, -1);
        // A chunk's pieces: its bytes in the region, [pieceFrom, pieceTo), in runs of
        // pieceStep bytes (whole blosc blocks) from pieceBase, the start of the block
        // pieceFrom is in; or all of them in one piece, for a pieceStep of 0. About two
        // pieces per thread in all.
        std::vector<uint64_t> pieceFrom(nItems, 0), pieceTo(nItems, 0), pieceBase(nItems, 0), pieceStep(nItems, 0);
        std::vector<int64_t> perChunk(nGroups);
        // Each group's reads and copies (every slab, each from the chunks in it; a copy from
        // a chunk that was never written does nothing) are handed out in contiguous blocks,
        // its pieces from a shared counter. A thread waits only for the chunks it needs:
        // read before decoding a piece of it, decoded before copying from it, and copied out
        // before reading the next group's chunk into its slot. (Waiting for all of the threads
        // between the steps also waited for any of them the system had paused, which made
        // reads under other load up to 20% slower.)
        std::vector<std::unique_ptr<WorkRanges>> readWork(nGroups), copyWork(nGroups);
        std::vector<std::vector<std::pair<uint64_t, int64_t>>> copies(nGroups);
        std::unique_ptr<std::atomic<int64_t>[]> nextPiece(new std::atomic<int64_t>[nGroups]);
        std::unique_ptr<std::atomic<int>[]> readDone(new std::atomic<int>[nItems]);
        std::unique_ptr<std::atomic<int64_t>[]> piecesLeft(new std::atomic<int64_t>[nItems]);
        std::unique_ptr<std::atomic<int64_t>[]> copiesLeft(new std::atomic<int64_t>[nItems]);
        // (in a sparse read, the pieces of a chunk that are not all zeros)
        std::unique_ptr<std::atomic<int64_t>[]> nonZero(new std::atomic<int64_t>[nItems]);
        try{
            for(int64_t g = 0; g < nGroups; g++){
                const int64_t g0 = g*group, n = std::min(group, nItems-g0);
                perChunk[g] = std::max<int64_t>(1, (2*(int64_t)numWorkers+n-1)/n);
                for(int64_t i = g0; i < g0+n; i++){
                    readDone[i].store(0);
                    piecesLeft[i].store(perChunk[g]);
                    copiesLeft[i].store(0);
                    nonZero[i].store(0);
                }
                for(uint64_t sl = 0; sl < nSlabs; sl++){
                    const uint64_t lo = startCoords.empty() ? 0 : startCoords[sAxis]+sl*slab;
                    const uint64_t hi = startCoords.empty() ? 1 : std::min(endCoords[sAxis], lo+slab);
                    for(int64_t i = g0; i < g0+n; i++){
                        if(nDims){
                            const uint64_t c = Zarr.get_chunks(sAxis), b = chunkAV[i][sAxis]*c;
                            if(b >= hi || b+c <= lo) continue;
                        }
                        copies[g].push_back({sl, i});
                        copiesLeft[i].fetch_add(1);
                    }
                }
                readWork[g].reset(new WorkRanges(numWorkers, (err || !slots.back()) ? 0 : n));
                copyWork[g].reset(new WorkRanges(numWorkers, (err || !slots.back()) ? 0 : (int64_t)copies[g].size()));
                nextPiece[g].store((err || !slots.back()) ? n*perChunk[g] : 0);
            }
        }
        catch(const std::exception &e){
            fail(std::string("Read error: ")+e.what()+"\n");
        }
        auto failed = [&]{ return err != 0; };
        // Fusable chunks: each copy decodes the parts of its chunk that hold its slab (and
        // waits only for the chunk to be read). A chunk that cannot be decoded in parts
        // (unusual blocks) is decoded whole, once, by the first copy that needs it.
        if(fusable){
            std::vector<uint64_t> fUnit(nItems, 0), fType(nItems, 0);
            std::vector<uint8_t> wholeZero(nItems, 0);
            // (0: not decoded, 1: being decoded, 2: decoded, 3: failed)
            std::unique_ptr<std::atomic<int>[]> wholeState(new std::atomic<int>[nItems]);
            for(int64_t i = 0; i < nItems; i++) wholeState[i].store(0);
            #pragma omp parallel num_threads(numWorkers)
            {
                const LocalMemoryPolicy memoryPolicy;
                const int w = omp_get_thread_num();
                std::unique_ptr<ReadScratch> scratch;
                ShardState ss;
                std::vector<uint64_t> boxExt(nDims);
                try{
                    scratch = takeScratch();
                    ss.index = scratch->shardIndex.data();
                }
                catch(const std::exception &e){
                    fail(std::string("Read error: ")+e.what()+"\n");
                }
                for(int64_t g = 0; g < nGroups && scratch && !err; g++){
                    const int64_t g0 = g*group;
                    // 1: read the group's chunks, each into its slot once the chunk before it
                    // there has been copied out
                    for(int64_t k = readWork[g]->next(w); k >= 0 && !err; k = readWork[g]->next(w)){
                        const int64_t i = g0+k;
                        try{
                            if(g > 0 && !waitFor([&]{ return copiesLeft[i-group].load(std::memory_order_acquire) == 0; }, failed)) break;
                            const std::string &name = Zarr.get_chunkNames(chunkF[i]);
                            const std::string subfolderName = Zarr.get_subfoldersString(Zarr.get_shard() ? Zarr.chunkToShard(chunkAV[i]) : chunkAV[i]);
                            uint8_t* data = nullptr;
                            const int64_t len = readChunk(ss, name, chunkAV[i], subfolderName, [&](const uint64_t m){
                                return slots[k]->storedBuffer(m ? m : 1, storedRoom);
                            }, data);
                            if(len == -2) break;
                            if(len == -1) noteMissing(chunkAV[i]);
                            storedOf[i] = data;
                            storedLen[i] = len;
                            fUnit[i] = len >= 0 ? fusedUnit(data, len, fType[i]) : 0;
                            if(len >= 0 && !fUnit[i]) decodedOf[i] = slots[k]->decodedBuffer(sB);
                            readDone[i].store(1, std::memory_order_release);
                        }
                        catch(const std::exception &e){
                            fail(std::string("Read error: ")+e.what()+"\n");
                            break;
                        }
                    }
                    // 2: copy the output, a slab at a time, decoding what each copy needs
                    for(int64_t k = err ? -1 : copyWork[g]->next(w); k >= 0 && !err; k = err ? -1 : copyWork[g]->next(w)){
                        const uint64_t sl = copies[g][k].first;
                        const int64_t i = copies[g][k].second;
                        try{
                            if(!waitFor([&]{ return readDone[i].load(std::memory_order_acquire) != 0; }, failed)) break;
                            if(storedLen[i] >= 0 && fUnit[i]){
                                // (the slab's slices of the chunk)
                                const uint64_t c = Zarr.get_chunks(sAxis), b = chunkAV[i][sAxis]*c;
                                const uint64_t lo = std::max(b, startCoords[sAxis]+sl*slab);
                                const uint64_t hi = std::min(std::min(b+c, endCoords[sAxis]), startCoords[sAxis]+(sl+1)*slab);
                                if(hi > lo && !copyFused(*scratch, storedOf[i], storedLen[i], fUnit[i], fType[i], chunkAV[i], boxExt, lo-b, hi-b, 0,
                                                         Zarr.get_subfoldersString(chunkAV[i]), Zarr.get_chunkNames(chunkF[i]))) break;
                            }
                            else if(storedLen[i] >= 0){
                                int whole = 0;
                                if(wholeState[i].compare_exchange_strong(whole, 1)){
                                    uint64_t from = 0, to = sB;
                                    bytesInRegion(chunkAV[i], from, to);
                                    const bool ok = decodeChunkRange(*scratch, scratch->dctx, storedOf[i], storedLen[i], decodedOf[i], from, to,
                                                                     Zarr.get_subfoldersString(chunkAV[i]), Zarr.get_chunkNames(chunkF[i]));
                                    wholeZero[i] = ok && zeroChunkUnc && !memcmp(zeroChunkUnc, decodedOf[i]+from, to-from);
                                    wholeState[i].store(ok ? 2 : 3, std::memory_order_release);
                                }
                                else if(!waitFor([&]{ return wholeState[i].load(std::memory_order_acquire) >= 2; }, failed)) break;
                                if(wholeState[i].load(std::memory_order_acquire) != 2) break;
                                if(!wholeZero[i]) copyChunk(decodedOf[i], chunkAV[i], boxExt, sAxis, sl*slab, std::min(extent, (sl+1)*slab));
                            }
                            copiesLeft[i].fetch_sub(1, std::memory_order_release);
                        }
                        catch(const std::exception &e){
                            fail(std::string("Read error: ")+e.what()+"\n");
                            break;
                        }
                    }
                }
                try{
                    readScratchPool.give(omp_get_thread_num(), std::move(scratch));
                }
                catch(const std::exception&){}
            }
        }
        else{
            #pragma omp parallel num_threads(numWorkers)
            {
                const LocalMemoryPolicy memoryPolicy;
                const int w = omp_get_thread_num();
                std::unique_ptr<ReadScratch> scratch;
                ShardState ss;
                std::vector<uint64_t> boxExt(nDims);
                try{
                    scratch = takeScratch();
                    ss.index = scratch->shardIndex.data();
                }
                catch(const std::exception &e){
                    fail(std::string("Read error: ")+e.what()+"\n");
                }
                for(int64_t g = 0; g < nGroups && scratch && !err; g++){
                    const int64_t g0 = g*group, n = std::min(group, nItems-g0);
                    // 1: read the group's chunks, each into its slot once the chunk before it
                    // there has been copied out
                    for(int64_t k = readWork[g]->next(w); k >= 0 && !err; k = readWork[g]->next(w)){
                        const int64_t i = g0+k;
                        try{
                            if(g > 0 && !waitFor([&]{ return piecesLeft[i-group].load(std::memory_order_acquire) == 0 &&
                                                             copiesLeft[i-group].load(std::memory_order_acquire) == 0; }, failed)) break;
                            const std::string &name = Zarr.get_chunkNames(chunkF[i]);
                            const std::string subfolderName = Zarr.get_subfoldersString(Zarr.get_shard() ? Zarr.chunkToShard(chunkAV[i]) : chunkAV[i]);
                            uint8_t* data = nullptr;
                            const int64_t len = readChunk(ss, name, chunkAV[i], subfolderName, [&](const uint64_t m){
                                return slots[k]->storedBuffer(m ? m : 1, storedRoom);
                            }, data);
                            if(len == -2) break;
                            if(len == -1) noteMissing(chunkAV[i]);
                            storedOf[i] = data;
                            storedLen[i] = len;
                            if(len >= 0){
                                if(compressor != "none") decodedOf[i] = slots[k]->decodedBuffer(sB);
                                // Its pieces: in runs of whole blocks of a blosc chunk (as
                                // decodeChunkRange decodes them)
                                uint64_t from = 0, to = sB;
                                bytesInRegion(chunkAV[i], from, to);
                                pieceFrom[i] = from;
                                pieceTo[i] = to;
                                int32_t nbytes = 0, blocksize = 0;
                                if(compressor == "blosc" && len >= BLOSC_MIN_HEADER_LENGTH &&
                                   blosc2_cbuffer_sizes(data, &nbytes, nullptr, &blocksize) >= 0 && (uint64_t)nbytes == sB && blocksize > 0){
                                    const uint64_t bs = (uint64_t)blocksize;
                                    const uint64_t blocks = (to+bs-1)/bs-from/bs;
                                    pieceBase[i] = (from/bs)*bs;
                                    pieceStep[i] = std::max<uint64_t>(1, (blocks+perChunk[g]-1)/perChunk[g])*bs;
                                }
                            }
                            readDone[i].store(1, std::memory_order_release);
                        }
                        catch(const std::exception &e){
                            fail(std::string("Read error: ")+e.what()+"\n");
                            break;
                        }
                    }
                    // 2: decode them, a piece at a time (each chunk's pieces in turn)
                    const int64_t nPieces = n*perChunk[g];
                    for(int64_t p = err ? nPieces : nextPiece[g].fetch_add(1); p < nPieces && !err; p = nextPiece[g].fetch_add(1)){
                        const int64_t i = g0+p/perChunk[g], j = p%perChunk[g];
                        try{
                            if(!waitFor([&]{ return readDone[i].load(std::memory_order_acquire) != 0; }, failed)) break;
                            uint64_t a = 0, e = 0;
                            if(storedLen[i] >= 0){
                                if(pieceStep[i]){
                                    a = std::max(pieceFrom[i], pieceBase[i]+j*pieceStep[i]);
                                    e = std::min(pieceTo[i], pieceBase[i]+(j+1)*pieceStep[i]);
                                }
                                else if(j == 0){
                                    a = pieceFrom[i];
                                    e = pieceTo[i];
                                }
                            }
                            if(e > a){
                                const uint8_t* chunkData = storedOf[i];
                                if(compressor != "none"){
                                    if(!decodeChunkRange(*scratch, scratch->dctx, storedOf[i], storedLen[i], decodedOf[i], a, e,
                                                         Zarr.get_subfoldersString(chunkAV[i]), Zarr.get_chunkNames(chunkF[i]))) break;
                                    chunkData = decodedOf[i];
                                }
                                if(zeroChunkUnc && memcmp(zeroChunkUnc, chunkData+a, e-a)) nonZero[i].fetch_add(1, std::memory_order_relaxed);
                            }
                            piecesLeft[i].fetch_sub(1, std::memory_order_acq_rel);
                        }
                        catch(const std::exception &e){
                            fail(std::string("Read error: ")+e.what()+"\n");
                            break;
                        }
                    }
                    // 3: copy the output, a slab at a time, from each chunk once it is decoded
                    for(int64_t k = err ? -1 : copyWork[g]->next(w); k >= 0 && !err; k = err ? -1 : copyWork[g]->next(w)){
                        const uint64_t sl = copies[g][k].first;
                        const int64_t i = copies[g][k].second;
                        try{
                            if(!waitFor([&]{ return piecesLeft[i].load(std::memory_order_acquire) == 0; }, failed)) break;
                            // (a sparse read leaves out chunks that are all zeros)
                            if(storedLen[i] >= 0 && !(zeroChunkUnc && nonZero[i].load(std::memory_order_relaxed) == 0)){
                                const uint8_t* chunkData = compressor == "none" ? storedOf[i] : decodedOf[i];
                                copyChunk(chunkData, chunkAV[i], boxExt, sAxis, sl*slab, std::min(extent, (sl+1)*slab));
                            }
                            copiesLeft[i].fetch_sub(1, std::memory_order_release);
                        }
                        catch(const std::exception &e){
                            fail(std::string("Read error: ")+e.what()+"\n");
                            break;
                        }
                    }
                }
                try{
                    readScratchPool.give(omp_get_thread_num(), std::move(scratch));
                }
                catch(const std::exception&){}
            }
        }
        try{
            for(auto &sl : slots) readScratchPool.give(-1, std::move(sl));
        }
        catch(const std::exception&){}
    }
    // Otherwise a chunk at a time per thread: read, decoded and copied. Chunks of 1 MB or
    // more that are decoded whole get as many threads as the fewest rounds of one chunk per
    // thread need (22 for 65 chunks on 24 threads): more would take as long, and fault in
    // more fresh buffers of their size (which the pool does not keep), which made reads
    // without huge pages up to 11% slower.
    else{
        const int32_t team = (int32_t)std::max<int64_t>(1, sB >= bigChunk && !fusable && !zstdFusable && perThread > 0 ?
                                                           (nItems+perThread-1)/perThread : std::min<int64_t>(numWorkers, nItems));
        WorkRanges work(team, nItems);
        // (in a region of all of the threads, those past the team doing nothing: OpenMP ends
        // the threads a smaller region leaves out, and starts them again for the next bigger
        // one, which made reads of small regions of changing sizes up to 8 times slower)
        #pragma omp parallel num_threads(numWorkers)
        {
            const LocalMemoryPolicy memoryPolicy;
            if(omp_get_thread_num() < team)
            {
                const int w = omp_get_thread_num();
                // This thread's buffers (a chunk's stored bytes and the chunk decoded) and
                // decoders, from the pool, and the shard it is reading
                std::unique_ptr<ReadScratch> scratch;
                ShardState ss;
                uint8_t* decoded = nullptr;
                std::vector<uint64_t> boxExt(nDims);
                try{
                    scratch = takeScratch();
                    ss.index = scratch->shardIndex.data();
                }
                catch(const std::exception &e){
                    fail(std::string("Read error: ")+e.what()+"\n");
                }
                // An exception cannot leave the parallel region: report it as an error
                for(int64_t it = (err || !scratch) ? -1 : work.next(w); it >= 0; it = err ? -1 : work.next(w)){
                    const int64_t f = Zarr.get_shard() ? items[it] : listIndex(it);
                    try{
                        const std::string &name = Zarr.get_chunkNames(f);
                        const std::vector<uint64_t> cAV = Zarr.get_chunkAxisVals(name);
                        // (a shard's subfolder comes from the shard's own indices)
                        const std::string subfolderName = Zarr.get_subfoldersString(Zarr.get_shard() ? Zarr.chunkToShard(cAV) : cAV);
                        uint8_t* stored = nullptr;
                        const int64_t fileLen = readChunk(ss, name, cAV, subfolderName, [&](const uint64_t n){
                            return scratch->storedBuffer(n, storedRoom);
                        }, stored);
                        if(fileLen == -2) break;
                        // On an NVMe drive, the thread's next chunk file is read ahead while this one
                        // is decoded, which keeps more reads in flight: cold whole-array reads took
                        // up to 6% longer without it
                        if(!Zarr.get_shard()){
                            const int64_t nx = work.peek(w);
                            if(nx >= 0 && disk() == 2){
                                const std::string &nextName = Zarr.get_chunkNames(listIndex(nx));
                                readAheadFile(Zarr.get_fileName()+"/"+Zarr.get_subfoldersString(Zarr.get_chunkAxisVals(nextName))+"/"+Zarr.chunkKey(nextName));
                            }
                        }
                        if(fileLen == -1){
                            noteMissing(cAV);
                            continue;
                        }

                        // A fusable chunk a part at a time, from a part that depends on its place
                        // in the list
                        uint64_t from = 0, to = sB, typesize = 0;
                        bytesInRegion(cAV, from, to);
                        const uint64_t unit = fusable ? fusedUnit(stored, fileLen, typesize) : 0;
                        if(unit){
                            if(!copyFused(*scratch, stored, fileLen, unit, typesize, cAV, boxExt, from/sliceBytes, to/sliceBytes, (uint64_t)it,
                                          subfolderName, name)) break;
                            continue;
                        }
                        const int streamed = zstdFusable ? zstdStreamed(stored, fileLen) : 0;
                        if(streamed){
                            if(!copyZstdStream(*scratch, stored, fileLen, streamed == 2, cAV, boxExt, from/sliceBytes, to/sliceBytes,
                                               subfolderName, name)) break;
                            continue;
                        }
                        // Otherwise decompressed (a blosc chunk only for its slices inside the
                        // region); uncompressed chunks are used straight from the file buffer
                        const uint8_t* chunkData = stored;
                        if(compressor != "none"){
                            if(!decoded) decoded = scratch->decodedBuffer(sB);
                            if(!decodeChunkRange(*scratch, scratch->dctx, stored, fileLen, decoded, from, to, subfolderName, name)) break;
                            chunkData = decoded;
                        }

                        // If the chunk is all zeros (in the region's slices), skip it
                        if(zeroChunkUnc && !memcmp(zeroChunkUnc, chunkData+from, to-from)) continue;

                        copyChunk(chunkData, cAV, boxExt, -1, 0, 0);
                    }
                    catch(const std::exception &e){
                        fail(std::string("Read error: ")+e.what()+"\n");
                        break;
                    }
                }
                try{
                    readScratchPool.give(omp_get_thread_num(), std::move(scratch));
                }
                catch(const std::exception&){}
            }
        }
    }
    // The missing chunks' part of an output that did not hold the fill value
    if(fillMissing && !err && !missing.empty()){
        const int64_t nMissing = (int64_t)missing.size();
        #pragma omp parallel num_threads(numWorkers)
        {
            const LocalMemoryPolicy memoryPolicy;
            std::vector<uint64_t> boxExt(nDims);
            #pragma omp for schedule(dynamic)
            for(int64_t i = 0; i < nMissing; i++) fillChunk(missing[i], boxExt);
        }
    }
    if(!useCtx){
        blosc2_destroy();
    }
    free(zeroChunkUnc);

    if(err){
        Zarr.set_errString(errString);
        return 1;
    }
    uint64_t readSize = 1;
    for(uint64_t d = 0; d < nDims; d++) readSize *= readShape[d];
    // (each part of a complex number on its own)
    if(oppositeEndianness(Zarr.get_dtype())) swapDtypeEndianness(zarrArr, dtypeInfo(Zarr.get_dtype()), readSize);
    return 0;
}

uint8_t decodeChunk(const std::string &compressor, const bool chunkChecksum,
                    const void* src, uint64_t srcLen, void* dst, const uint64_t dstLen,
                    std::string &why){
    const uint8_t* in = static_cast<const uint8_t*>(src);
    if(chunkChecksum){
        uint32_t stored = 0;
        if(srcLen >= 4) memcpy(&stored, in+srcLen-4, 4);
        if(srcLen < 4 || stored != crc32c(in, (size_t)srcLen-4)){
            why = "does not match its checksum (it may be damaged)";
            return 1;
        }
        srcLen -= 4;
    }
    if(compressor == "none"){
        if(srcLen != dstLen){
            why = "is "+std::to_string(srcLen)+" bytes instead of "+std::to_string(dstLen);
            return 1;
        }
        memcpy(dst, in, dstLen);
        return 0;
    }
    if(compressor == "zstd"){
        const size_t n = ZSTD_decompress(dst, dstLen, in, srcLen);
        if(ZSTD_isError(n) || n != dstLen){
            why = "could not be decompressed (zstd)";
            return 1;
        }
        return 0;
    }
    if(compressor == "blosc"){
        blosc2_context* dctx = blosc2_create_dctx(BLOSC2_DPARAMS_DEFAULTS);
        const int n = blosc2_decompress_ctx(dctx, in, (int32_t)srcLen, dst, (int32_t)dstLen);
        blosc2_free_ctx(dctx);
        if(n != (int)dstLen){
            why = "could not be decompressed (blosc error "+std::to_string(n)+")";
            return 1;
        }
        return 0;
    }
    if(compressor == "gzip"){
        // gzip or zlib streams (header detected), possibly several gzip members
        z_stream stream;
        memset(&stream, 0, sizeof(stream));
        uint64_t used = 0, made = 0;
        while(used < srcLen){
            if(inflateInit2(&stream, 32) != Z_OK){ why = "could not be decompressed (zlib)"; return 1; }
            stream.next_in = const_cast<uint8_t*>(in)+used;
            stream.avail_in = (uInt)(srcLen-used);
            stream.next_out = static_cast<uint8_t*>(dst)+made;
            stream.avail_out = (uInt)(dstLen-made);
            const int r = inflate(&stream, Z_FINISH);
            used = srcLen-stream.avail_in;
            made = dstLen-stream.avail_out;
            inflateEnd(&stream);
            if(r != Z_STREAM_END){ why = "could not be decompressed (zlib error "+std::to_string(r)+")"; return 1; }
        }
        if(made != dstLen){ why = "decompressed to "+std::to_string(made)+" bytes instead of "+std::to_string(dstLen); return 1; }
        return 0;
    }
    why = "uses the \""+compressor+"\" compressor, which is not supported";
    return 1;
}

// TODO: FIX MEMORY LEAKS
// Wrapper used by parallelWriteZarr
void* parallelReadZarrWriteWrapper(zarr Zarr, const bool &crop,
                              std::vector<uint64_t> startCoords, 
                              std::vector<uint64_t> endCoords){
    return parallelReadZarrWriteWrapper(Zarr, crop, startCoords, endCoords, false);
}

void* parallelReadZarrWriteWrapper(zarr Zarr, const bool &crop,
                              std::vector<uint64_t> startCoords,
                              std::vector<uint64_t> endCoords,
                              const bool cOrder){
   
    const uint64_t nDims = Zarr.get_ndims();
    if(!crop){
        startCoords.assign(nDims, 0);
        endCoords.assign(nDims, 0);
        for(uint64_t d = 0; d < nDims; d++) endCoords[d] = Zarr.get_shape(d);
    }

    std::vector<uint64_t> readShape(nDims);
    uint64_t readSize = 1;
    for(uint64_t d = 0; d < nDims; d++){
        readShape[d] = endCoords[d]-startCoords[d];
        readSize *= readShape[d];
    }

    Zarr.set_chunkInfo(startCoords, endCoords);

    // The chunk read/copy machinery is element-width based, so one generic path
    // covers every supported dtype: booleans, signed/unsigned 8/16/32/64-bit
    // integers, 16/32/64-bit floats and bfloat16, and complex numbers
    const std::string &dtype = Zarr.get_dtype();
    const uint64_t bytes = Zarr.dtypeBytes();
    if(!dtypeInfo(dtype).kind) return NULL;

    // The output starts as the fill value, which is what missing chunks read as
    uint8_t fillElem[16];
    fillValueElement(Zarr.get_fill_value(), dtype, fillElem);
    void* zarrArr = nullptr;
    if(std::any_of(fillElem, fillElem+bytes, [](uint8_t b){ return b != 0; })){
        zarrArr = malloc(readSize*bytes);
        if(zarrArr) fillElements(bytes, zarrArr, readSize, fillElem);
    }
    else zarrArr = calloc(readSize,bytes);
    if(!zarrArr) return NULL;
    uint8_t err = parallelReadZarr(Zarr, zarrArr,startCoords,endCoords,readShape,bytes*8,true,false,cOrder);
    if(err){
        free(zarrArr);
        return NULL;
    }
    return zarrArr;
}

void* readZarrParallelHelper(const char* folderName, uint64_t startX, uint64_t startY, uint64_t startZ, uint64_t endX, uint64_t endY, uint64_t endZ, uint8_t imageJIm){
    // Arrays that cannot be read (missing or unsupported metadata) give NULL too
    zarr Zarr;
    try{
        Zarr = zarr(folderName);
    }
    catch(const std::string &){
        return NULL;
    }
    // This helper's interface is 3D
    if(Zarr.get_ndims() != 3) return NULL;
    void* zarrArr = parallelReadZarrWriteWrapper(Zarr, true,
                              {startX, startY, startZ},
                              {endX, endY, endZ});
    // Unsupported dtype or read error: return NULL rather than crashing below
    if(!zarrArr) return NULL;
    // May need to add a check for if the data is f order or c order for ImageJ
    // For the c order images I have tested, we also have to do this flip for now
    if(imageJIm /*&& (order == 'F' || order == 'f')*/){
        // The buffer holds the read window, not the whole array, so allocate and
        // index with the window's dimensions
        const uint64_t rx = endX-startX;
        const uint64_t ry = endY-startY;
        const uint64_t rz = endZ-startZ;
        const uint64_t es = Zarr.dtypeBytes();
        void* zarrArrC = malloc(rx*ry*rz*es);
		#pragma omp parallel for
        for(uint64_t k = 0; k < rz; k++){
            for(uint64_t j = 0; j < ry; j++){
                for(uint64_t i = 0; i < rx; i++){
                    switch(es){
                        case 1:
                            ((uint8_t*)zarrArrC)[j+(i*ry)+(k*ry*rx)] = ((uint8_t*)zarrArr)[i+(j*rx)+(k*rx*ry)];
                            break;
                        case 2:
                            ((uint16_t*)zarrArrC)[j+(i*ry)+(k*ry*rx)] = ((uint16_t*)zarrArr)[i+(j*rx)+(k*rx*ry)];
                            break;
                        case 4:
                            ((float*)zarrArrC)[j+(i*ry)+(k*ry*rx)] = ((float*)zarrArr)[i+(j*rx)+(k*rx*ry)];
                            break;
                        case 8:
                            ((double*)zarrArrC)[j+(i*ry)+(k*ry*rx)] = ((double*)zarrArr)[i+(j*rx)+(k*rx*ry)];
                            break;
                        default:
                            memcpy((uint8_t*)zarrArrC+(j+(i*ry)+(k*ry*rx))*es, (uint8_t*)zarrArr+(i+(j*rx)+(k*rx*ry))*es, es);
                    }
                }
            }
        }
		free(zarrArr);
        return zarrArrC;
    }
    return zarrArr;
}

