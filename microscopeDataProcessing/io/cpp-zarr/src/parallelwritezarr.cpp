#include <algorithm>
#include <thread>
#include <mutex>
#include <memory>
#include <limits>
#include <chrono>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <omp.h>
#ifdef _WIN32
#include <sys/time.h>
#endif
#include <fstream>
#include "blosc.h"
#include "parallelreadzarr.h"
#include "parallelwritezarr.h"
#include "helperfunctions.h"
#include "parallelwork.h"
#include "zarr.h"
#include "zlib.h"
#include "zstd.h"

// Open a chunk or shard file for writing, making its folder first if it is missing
// (an array whose metadata another program wrote may not have its folders yet)
static bool openForWrite(std::ofstream &file, const std::string &path, const std::ios::openmode mode){
    file.open(path, std::ios::binary | mode);
    if(file.is_open()) return true;
    const size_t slash = path.find_last_of("/\\");
    if(slash == std::string::npos) return false;
    mkdirRecursive(path.substr(0, slash).c_str());
    file.clear();
    file.open(path, std::ios::binary | mode);
    return file.is_open();
}

uint8_t parallelWriteZarr(zarr &Zarr, void* zarrArr,
                          const std::vector<uint64_t> &startCoords,
                          const std::vector<uint64_t> &endCoords,
                          const std::vector<uint64_t> &writeShape,
                          const uint64_t bits, const bool useUuid,
                          const bool crop, const bool sparse){
    // The input is in F order: the first axis is contiguous
    std::vector<uint64_t> inStrides(writeShape.size());
    uint64_t acc = 1;
    for(uint64_t d = 0; d < writeShape.size(); d++){ inStrides[d] = acc; acc *= writeShape[d]; }
    return parallelWriteZarr(Zarr, zarrArr, startCoords, endCoords, writeShape, inStrides, bits, useUuid, crop, sparse);
}

// What one chunk becomes when it is stored: nothing (a sharded array's inner chunk
// past the array's edge, an all-zero chunk of a sparse write, or an inner chunk a
// crop write leaves alone that its shard does not hold), or its bytes, with the
// CRC32C to put after them when the array has one
struct StoredChunk {
    bool empty = true;
    bool zero = false;          // all zero in a sparse write (any file it had is removed)
    const uint8_t* data = nullptr;
    uint64_t n = 0;
    bool hasCrc = false;
    uint32_t crc = 0;
};

// A shard being written. Its inner chunks are prepared by any of the threads, in any
// order, and added to its file in index order (one that is ready early waits in
// memory), so the shard's bytes are the same as when one thread wrote all of it. One
// thread at a time adds them (the appender, without holding m while it writes); a
// thread that hands over a chunk meanwhile leaves it for the appender. The inner
// chunks a crop write does not touch are not prepared at all: the appender copies
// their bytes from the shard's current file, a run of them at a time.
struct ShardOut {
    std::mutex m;
    bool started = false, opened = false, done = false, failed = false, appending = false;
    std::string path, writePath;
    std::ofstream file;
    std::vector<uint64_t> index;
    uint64_t bytes = 0;
    std::atomic<int64_t> next{0};
    // 0: not yet, 1: nothing stored, 2: bytes waiting, 3: its current bytes, unchanged
    std::vector<uint8_t> ready;
    std::vector<std::vector<uint8_t>> waiting;
    // The number of inner chunks to prepare before each position (and in all)
    std::vector<int64_t> preparedBefore;
    // Its inner chunks to prepare: toPrepare[prepNext, prepEnd), taken in order
    std::atomic<int64_t> prepNext{0};
    int64_t prepEnd = 0;
    // Crop writes: the shard's current file (one handle, read under oldM) and index
    std::mutex oldM;
    std::ifstream old;
    std::vector<uint64_t> oldIndex;
};

namespace {
// A writer thread's buffers and zstd encoder, kept in a pool between calls (see
// ScratchPool). Buffers for chunks over 1 MB, and an encoder over 8 MB, are not kept.
struct WriteScratch {
    std::unique_ptr<uint8_t[]> unc, comp;
    uint64_t uncCap = 0, compCap = 0;
    std::vector<uint8_t> stored, copyBuf;
    ZSTD_CCtx* zctx = nullptr;
    WriteScratch(){}
    WriteScratch(const WriteScratch&) = delete;
    WriteScratch& operator=(const WriteScratch&) = delete;
    ~WriteScratch(){
        if(zctx) ZSTD_freeCCtx(zctx);
    }
    // A chunk's uncompressed (n) and compressed (nc) buffers
    void buffers(const uint64_t n, const uint64_t nc){
        if(uncCap < n){
            unc.reset();
            unc.reset(new uint8_t[n]);
            uncCap = n;
        }
        if(compCap < nc){
            comp.reset();
            comp.reset(new uint8_t[nc]);
            compCap = nc;
        }
    }
    void shrink(){
        const uint64_t keep = 1 << 20;
        if(uncCap > keep){
            unc.reset();
            comp.reset();
            uncCap = compCap = 0;
        }
        if(stored.capacity() > keep) std::vector<uint8_t>().swap(stored);
        if(copyBuf.capacity() > keep) std::vector<uint8_t>().swap(copyBuf);
        if(zctx && ZSTD_sizeof_CCtx(zctx) > (8 << 20)){
            ZSTD_freeCCtx(zctx);
            zctx = nullptr;
        }
    }
};
ScratchPool<WriteScratch> writeScratchPool;
}

uint8_t parallelWriteZarr(zarr &Zarr, void* zarrArr,
                          const std::vector<uint64_t> &startCoords,
                          const std::vector<uint64_t> &endCoords,
                          const std::vector<uint64_t> &writeShape,
                          const std::vector<uint64_t> &inStrides,
                          const uint64_t bits, const bool useUuid,
                          const bool crop, const bool sparse){
    // A 0-dimensional array is a single element: sharding does not apply
    if(Zarr.get_ndims() == 0 && Zarr.get_shard()){
        Zarr.set_errString("Sharding is not supported for 0-dimensional arrays\n");
        return 1;
    }
    // Chunks are written exactly as the array's metadata (v2 or v3) says: its
    // compressor (blosc, gzip, numcodecs zlib, zstd or none) and a crc32c checksum
    // after each chunk if it has one. Arrays in the other byte order are not
    // written, which would end up with mixed chunks.
    const std::string compressor = Zarr.get_compressor();
    if(compressor != "none" && compressor != "blosc" && compressor != "gzip" && compressor != "zstd"){
        Zarr.set_errString("Writing arrays with the \""+Zarr.get_cname()+"\" compressor is not supported\n");
        return 1;
    }
    // (a numcodecs zlib stream has a zlib header instead of gzip's)
    const bool zlibStream = compressor == "gzip" && Zarr.get_cname() != "gzip";
    const bool chunkChecksum = Zarr.get_chunkChecksum();
    if(oppositeEndianness(Zarr.get_dtype())){
        Zarr.set_errString("Writing arrays of data type \""+Zarr.get_dtype()+"\" (opposite byte order) is not supported\n");
        return 1;
    }
    const uint64_t bytes = (bits/8);
    const int32_t numWorkers = omp_get_max_threads();

    // A sharded array is written with its inner chunk shape. The shard shape is
    // restored on return.
    std::vector<uint64_t> shardChunks;
    if(Zarr.get_shard()){
        std::vector<uint64_t> innerChunks(Zarr.get_ndims());
        for(uint64_t d = 0; d < innerChunks.size(); d++){
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

    const uint64_t nDims = Zarr.get_ndims();
    std::vector<uint64_t> chunkDims(nDims);
    uint64_t s = 1;
    for(uint64_t d = 0; d < nDims; d++){
        chunkDims[d] = Zarr.get_chunks(d);
        s *= chunkDims[d];
    }
    const uint64_t sB = s*bytes;

    if(writeShape.size() != nDims || inStrides.size() != nDims){
        Zarr.set_errString("The input has "+std::to_string(writeShape.size())+" dimensions and "+
                           std::to_string(inStrides.size())+" strides but the array has "+
                           std::to_string(nDims)+" dimensions\n");
        return 1;
    }

    // Element strides of an uncompressed chunk (its stored axes, slowest first: C
    // or F order, or a Zarr v3 transpose)
    std::vector<uint64_t> chunkStrides(nDims);
    {
        const std::vector<uint64_t> axes = Zarr.get_chunkAxisOrder();
        uint64_t acc = 1;
        for(int64_t i = (int64_t)nDims-1; i >= 0; i--){ chunkStrides[axes[i]] = acc; acc *= chunkDims[axes[i]]; }
    }

    // The fill value as one element, in every order and data type: what the
    // elements of a new chunk outside the written region get (they read back as
    // the fill value, as missing chunks do), what inner chunks a crop write
    // completes get when they are missing from their shard, and what v2 chunks a
    // crop write completes get when they do not exist yet
    uint8_t fillElem[16];
    fillValueElement(Zarr.get_fill_value(), Zarr.get_dtype(), fillElem);

    // blosc (both blosc1 and blosc2) stores each chunk as a single frame with
    // 32-bit size fields, so one chunk cannot exceed BLOSC_MAX_BUFFERSIZE
    // (~2 GB) uncompressed. Exceeding it makes blosc_compress_ctx fail and, left
    // unchecked, wrote empty chunk files -- silently losing data. Reject it up
    // front with an actionable message. (gzip uses zlib's streaming API and is
    // not subject to this limit, so it is allowed through.)
    if(compressor == "blosc" && sB > (uint64_t)BLOSC_MAX_BUFFERSIZE){
        std::string dimsStr, productStr;
        for(uint64_t d = 0; d < nDims; d++){
            dimsStr += (d ? "x" : "")+std::to_string(chunkDims[d]);
            productStr += (d ? "*chunk[" : "chunk[")+std::to_string(d)+"]";
        }
        Zarr.set_errString("Chunk is too large for the \""+Zarr.get_cname()+
            "\" compressor: "+dimsStr+
            " x "+std::to_string(bytes)+" bytes/element = "+std::to_string(sB)+
            " bytes exceeds the blosc limit of "+std::to_string((uint64_t)BLOSC_MAX_BUFFERSIZE)+
            " bytes (~2 GB). Reduce the chunk size so that "+
            productStr+"*dtypeBytes < "+
            std::to_string((uint64_t)BLOSC_MAX_BUFFERSIZE+1)+", or use the gzip compressor.\n");
        return 1;
    }

    const std::string uuid(generateUUID());
    // A shard goes to a temporary file that is then moved into place when asked
    // (useUuid), and always when cropping: the shard's current contents are read
    // while its new ones are written
    const bool shardTmpFile = useUuid || crop;
    // Where a shard's index goes (the start or, by default, the end of its file)
    // and whether a CRC32C of it follows, as the array's metadata says
    const bool indexAtStart = Zarr.get_shardIndexAtStart();
    const bool indexChecksum = Zarr.get_shardIndexChecksum();
    const bool sharded = Zarr.get_shard();
    const int64_t perShard = sharded ? (int64_t)Zarr.get_numChunksPerShard() : 1;
    const uint64_t emptyEntry = std::numeric_limits<uint64_t>::max();
    const uint64_t indexBytes = 2*perShard*sizeof(uint64_t)+(indexChecksum ? sizeof(uint32_t) : 0);

    const uint8_t* const typedFill = fillElem;
    // Blosc's shuffle: as a Zarr v3 array's codec says, and always byte shuffle in
    // v2 arrays, as cpp-zarr has always written them
    const uint64_t shuffle = Zarr.get_shuffle();
    const int bloscShuffle = Zarr.get_zarr_format() != 3 || shuffle > 2 ? BLOSC_SHUFFLE : (int)shuffle;

    // Sparse writes leave all-zero chunks unwritten, when that is how a missing
    // chunk reads back (a fill value whose bytes are all zero)
    void* zeroChunkUnc = NULL;
    if(sparse && fillValueIsZero(Zarr.get_fill_value(), Zarr.get_dtype())){
        zeroChunkUnc = calloc(s,bytes);
    }

    std::atomic<int> err{0};
    std::string errString;
    auto reportError = [&](const std::string &message){
        #pragma omp critical
        {
            err = 1;
            errString = message;
        }
    };

    // Sharded arrays: the inner chunks to prepare are those the write touches. The rest
    // of each shard's inner chunks are past the array's edge (not stored), or keep the
    // bytes they have (crop writes).
    const int64_t numChunks = (int64_t)Zarr.get_numChunks();
    // (0: untouched, 1: touched, 2: past the array's edge)
    std::vector<uint8_t> touched;
    std::vector<int64_t> toPrepare;
    if(sharded){
        touched.assign(numChunks, 0);
        // An inner chunk's indices follow from its place in the list (shard by shard, each
        // a grid of inner chunks, both last axis fastest), without parsing its name: a
        // parallel pass over the names was another wait for every thread, which with other
        // load on the machine made small writes take up to 1.6 times as long
        std::vector<uint64_t> outerFirst(nDims), outerCount(nDims), innerCount(nDims), idx(nDims);
        uint64_t nOuter = 1, nPer = 1;
        for(uint64_t d = 0; d < nDims; d++){
            innerCount[d] = (shardChunks[d]+chunkDims[d]-1)/chunkDims[d];
            const uint64_t extent = innerCount[d]*chunkDims[d];
            const uint64_t first = startCoords[d]/extent, last = endCoords[d]/extent+(endCoords[d]%extent ? 1 : 0);
            outerFirst[d] = first;
            outerCount[d] = last > first ? last-first : 0;
            nOuter *= outerCount[d];
            nPer *= innerCount[d];
        }
        const bool regular = nOuter*nPer == (uint64_t)numChunks;
        for(int64_t f = 0; f < numChunks && !err; f++){
            try{
                if(regular){
                    uint64_t o = (uint64_t)f/nPer, c = (uint64_t)f%nPer;
                    for(int64_t d = (int64_t)nDims-1; d >= 0; d--){
                        idx[d] = (outerFirst[d]+o%outerCount[d])*innerCount[d]+c%innerCount[d];
                        o /= outerCount[d];
                        c /= innerCount[d];
                    }
                }
                else idx = Zarr.get_chunkAxisVals(Zarr.get_chunkNames(f));
                uint8_t kind = 1;
                for(uint64_t d = 0; d < nDims && d < idx.size(); d++){
                    const uint64_t lo = idx[d]*chunkDims[d];
                    if(lo >= Zarr.get_shape(d)){
                        kind = 2;
                        break;
                    }
                    if(lo >= endCoords[d] || lo+chunkDims[d] <= startCoords[d]) kind = 0;
                }
                touched[f] = kind;
            }
            catch(const std::exception &e){
                reportError(std::string("Write error: ")+e.what()+"\n");
            }
        }
        for(int64_t f = 0; f < numChunks; f++) if(touched[f] == 1) toPrepare.push_back(f);
    }
    const int64_t numWork = sharded ? (int64_t)toPrepare.size() : numChunks;

    // No more threads than chunks to prepare. With fewer chunks than threads, blosc
    // compresses each chunk's blocks with the threads left over (for a sharded array,
    // counting every inner chunk of its shards, as before: blosc's threads store the
    // blocks in the order they finish, so this keeps the bytes the same)
    const int32_t team = (int32_t)std::max<int64_t>(1, std::min<int64_t>(numWorkers, numWork));
    int32_t nBloscThreads = 1;
    if(numChunks > 0 && numWorkers > numChunks) nBloscThreads = (int32_t)((numWorkers+numChunks-1)/numChunks);

    // gzip/deflate can expand incompressible data beyond blosc's
    // BLOSC_MAX_OVERHEAD guarantee, so size the compressed-chunk buffer to a
    // bound big enough for every codec (compressBound also covers gzip).
    const uint64_t chunkCCap = std::max<uint64_t>(compressBound(sB), ZSTD_compressBound(sB)) + 64;

    // The shards of the write: one thread to a shard while there are shards nobody has
    // started (it adds each inner chunk to the file itself, as soon as it is ready),
    // then the threads that run out help with the shards that have the most inner
    // chunks left, so that every thread has work even with fewer shards than threads.
    // (Handing out every shard's inner chunks to all the threads made the one thread
    // adding to a shard's file the bottleneck: up to 8x slower with many small inner
    // chunks.) A thread stays within a window of prepared inner chunks past the ones
    // its shard has added to its file, which bounds the memory of the ones that wait:
    // about two per thread, or 64 MB of small ones (a window of 16 small inner chunks
    // kept threads waiting while a shard's unchanged inner chunks were copied).
    const int64_t numShards = sharded ? numChunks/perShard : 0;
    std::unique_ptr<ShardOut[]> shards(sharded ? new ShardOut[numShards] : nullptr);
    for(int64_t k = 0; k < (int64_t)toPrepare.size(); k++){
        ShardOut &sh = shards[toPrepare[k]/perShard];
        if(sh.prepEnd == 0) sh.prepNext.store(k);
        sh.prepEnd = k+1;
    }
    std::atomic<int64_t> nextShard{0};
    const int64_t window = std::max<int64_t>({16, 2*(int64_t)team, (int64_t)((64ULL << 20)/std::max<uint64_t>(sB, 1))});

    // Open a shard (once, by whichever thread first has one of its inner chunks): its
    // path, an empty index, and for a crop write its current file and index, checked
    // as the reader checks it (a shard that does not exist yet has only empty
    // entries). Its inner chunks that are not prepared are ready from the start.
    auto startShard = [&](ShardOut &sh, const int64_t f) -> bool {
        std::lock_guard<std::mutex> g(sh.m);
        if(sh.started) return !sh.failed;
        sh.started = true;
        const std::vector<uint64_t> sAV = Zarr.chunkToShard(Zarr.get_chunkAxisVals(Zarr.get_chunkNames(f)));
        sh.path = Zarr.get_fileName()+"/"+Zarr.get_subfoldersString(sAV)+"/"+
                  Zarr.chunkKey(Zarr.chunkNameToShardName(Zarr.get_chunkNames(f)));
        sh.writePath = shardTmpFile ? sh.path+uuid : sh.path;
        sh.index.assign(2*perShard, emptyEntry);
        sh.bytes = indexAtStart ? indexBytes : 0;
        sh.ready.assign(perShard, 0);
        sh.waiting.resize(perShard);
        sh.oldIndex.assign(2*perShard, emptyEntry);
        if(crop) sh.old.open(sh.path, std::ios::binary);
        if(sh.old.is_open()){
            sh.old.seekg(0, std::ios::end);
            const uint64_t len = (uint64_t)sh.old.tellg();
            const uint64_t entryBytes = sh.oldIndex.size()*sizeof(uint64_t);
            std::string problem;
            if(len < indexBytes) problem = "is smaller than its index";
            else{
                sh.old.seekg(indexAtStart ? 0 : (std::streamoff)(len-indexBytes), std::ios::beg);
                sh.old.read(reinterpret_cast<char*>(sh.oldIndex.data()), entryBytes);
                uint32_t crc = 0;
                if(indexChecksum) sh.old.read(reinterpret_cast<char*>(&crc), sizeof(crc));
                if(!sh.old) problem = "could not be read";
                else if(indexChecksum && crc != crc32c(reinterpret_cast<const uint8_t*>(sh.oldIndex.data()), entryBytes)){
                    problem = "has an index that does not match its checksum";
                }
                else{
                    const uint64_t lo = indexAtStart ? indexBytes : 0, hi = indexAtStart ? len : len-indexBytes;
                    for(int64_t i = 0; i < perShard && problem.empty(); i++){
                        const uint64_t off = sh.oldIndex[2*i], nb = sh.oldIndex[(2*i)+1];
                        if(off == emptyEntry && nb == emptyEntry) continue;
                        if(off < lo || off > hi || nb > hi-off) problem = "has an index entry outside the file";
                    }
                }
            }
            if(!problem.empty()){
                sh.failed = true;
                reportError("The shard "+sh.path+" "+problem+" (it may be damaged, or not match the array's metadata)\n");
                return false;
            }
        }
        const int64_t first = (f/perShard)*perShard;
        sh.preparedBefore.assign(perShard+1, 0);
        for(int64_t pos = 0; pos < perShard; pos++){
            sh.preparedBefore[pos+1] = sh.preparedBefore[pos]+(touched[first+pos] == 1);
            if(touched[first+pos] == 1) continue;
            const bool isStored = touched[first+pos] == 0 &&
                                  !(sh.oldIndex[2*pos] == emptyEntry && sh.oldIndex[(2*pos)+1] == emptyEntry);
            sh.ready[pos] = isStored ? 3 : 1;
        }
        return true;
    };
    // The stored bytes of a shard's current inner chunk at position pos
    auto readOldInner = [&](ShardOut &sh, const int64_t pos, std::vector<uint8_t> &stored) -> bool {
        std::lock_guard<std::mutex> g(sh.oldM);
        stored.resize(sh.oldIndex[(2*pos)+1]);
        sh.old.clear();
        sh.old.seekg((std::streamoff)sh.oldIndex[2*pos], std::ios::beg);
        sh.old.read(reinterpret_cast<char*>(stored.data()), (std::streamsize)stored.size());
        if(!sh.old){
            reportError("Could not read the shard "+sh.path+"\n");
            return false;
        }
        return true;
    };
    // The shard's new file, opened by the first inner chunk added to it (with room for
    // an index at the start, written when the shard is finished)
    auto openShardFile = [&](ShardOut &sh) -> bool {
        if(sh.opened) return true;
        if(!openForWrite(sh.file, sh.writePath, std::ios::out | std::ios::trunc)){
            reportError("Check permissions or filepath. Cannot write to path: "+sh.writePath+"\n");
            return false;
        }
        sh.opened = true;
        if(indexAtStart){
            const std::vector<char> room(indexBytes, 0);
            sh.file.write(room.data(), room.size());
        }
        return true;
    };
    // Add an inner chunk to its shard's file (as its appender): its bytes and the
    // checksum after them, and where they are in the index
    auto appendToShard = [&](ShardOut &sh, const int64_t pos, const uint8_t* data, const uint64_t n, const uint32_t* crc) -> bool {
        if(!openShardFile(sh)) return false;
        sh.file.write(reinterpret_cast<const char*>(data), n);
        if(crc) sh.file.write(reinterpret_cast<const char*>(crc), sizeof(*crc));
        if(!sh.file){
            reportError("Could not write all of "+sh.writePath+" (is the disk full?)\n");
            return false;
        }
        const uint64_t total = n+(crc ? sizeof(*crc) : 0);
        sh.index[pos*2] = sh.bytes;
        sh.index[(pos*2)+1] = total;
        sh.bytes += total;
        return true;
    };
    // Add the current bytes of a shard's inner chunks [p, q), which lie one after
    // another in its current file, to its new file (as its appender), through buf
    auto copyOldInner = [&](ShardOut &sh, const int64_t p, const int64_t q, std::vector<uint8_t> &buf) -> bool {
        if(!openShardFile(sh)) return false;
        const uint64_t from = sh.oldIndex[2*p];
        const uint64_t total = sh.oldIndex[2*(q-1)]+sh.oldIndex[(2*(q-1))+1]-from;
        const uint64_t piece = std::min<uint64_t>(total, 8 << 20);
        if(buf.size() < piece) buf.resize(piece);
        for(uint64_t done = 0; done < total;){
            const uint64_t n = std::min(piece, total-done);
            {
                std::lock_guard<std::mutex> g(sh.oldM);
                sh.old.clear();
                sh.old.seekg((std::streamoff)(from+done), std::ios::beg);
                sh.old.read(reinterpret_cast<char*>(buf.data()), (std::streamsize)n);
                if(!sh.old){
                    reportError("Could not read the shard "+sh.path+"\n");
                    return false;
                }
            }
            sh.file.write(reinterpret_cast<const char*>(buf.data()), n);
            if(!sh.file){
                reportError("Could not write all of "+sh.writePath+" (is the disk full?)\n");
                return false;
            }
            done += n;
        }
        for(int64_t r = p; r < q; r++){
            sh.index[2*r] = sh.bytes+(sh.oldIndex[2*r]-from);
            sh.index[(2*r)+1] = sh.oldIndex[(2*r)+1];
        }
        sh.bytes += total;
        return true;
    };
    // Write a shard's index (and its CRC32C) after its inner chunks, or in the room
    // left for it at the start, then move the file into place (the shard's current
    // file is closed first: Windows cannot replace an open file)
    auto finishShard = [&](ShardOut &sh) -> bool {
        {
            std::lock_guard<std::mutex> g(sh.oldM);
            if(sh.old.is_open()) sh.old.close();
        }
        const uint32_t crc = crc32c(reinterpret_cast<const uint8_t*>(sh.index.data()), sh.index.size()*sizeof(uint64_t));
        if(!sh.opened){
            if(!openForWrite(sh.file, sh.writePath, std::ios::out | std::ios::trunc)){
                reportError("Check permissions or filepath. Cannot write to path: "+sh.writePath+"\n");
                return false;
            }
            sh.opened = true;
        }
        else if(indexAtStart) sh.file.seekp(0, std::ios::beg);
        sh.file.write(reinterpret_cast<const char*>(sh.index.data()), sh.index.size()*sizeof(uint64_t));
        if(indexChecksum) sh.file.write(reinterpret_cast<const char*>(&crc), sizeof(crc));
        sh.file.close();
        if(!sh.file){
            if(shardTmpFile) remove(sh.writePath.c_str());
            reportError("Could not write all of "+sh.writePath+" (is the disk full?)\n");
            return false;
        }
        if(shardTmpFile && !renameReplace(sh.writePath, sh.path)){
            remove(sh.writePath.c_str());
            reportError("Cannot move the temporary file into place: "+sh.path+"\n");
            return false;
        }
        sh.done = true;
        return true;
    };
    // As the shard's appender: add its ready inner chunks from next on to its file
    // (m is not held while writing), and finish the shard after its last one. Stops at
    // one that is not ready yet, which its thread then adds; after an error,
    // unfinished shards are not moved into place.
    auto appendReady = [&](ShardOut &sh, std::vector<uint8_t> &buf){
        while(true){
            int64_t first, last;
            std::vector<std::vector<uint8_t>> taken;
            {
                std::lock_guard<std::mutex> g(sh.m);
                first = sh.next.load();
                if(first == perShard){
                    if(!err && !sh.failed && !finishShard(sh)) sh.failed = true;
                    std::vector<uint8_t>().swap(sh.ready);
                    std::vector<std::vector<uint8_t>>().swap(sh.waiting);
                    sh.appending = false;
                    return;
                }
                last = first;
                if(!err && !sh.failed) while(last < perShard && sh.ready[last]) last++;
                if(last == first){
                    sh.appending = false;
                    return;
                }
                for(int64_t p = first; p < last; p++) if(sh.ready[p] == 2) taken.push_back(std::move(sh.waiting[p]));
            }
            bool ok = true;
            size_t k = 0;
            for(int64_t p = first; p < last && ok;){
                if(sh.ready[p] == 1) p++;
                else if(sh.ready[p] == 2){
                    // (the checksum, if any, is already at the end of the waiting bytes)
                    ok = appendToShard(sh, p, taken[k].data(), taken[k].size(), nullptr);
                    std::vector<uint8_t>().swap(taken[k]);
                    k++;
                    p++;
                }
                else{
                    // a run of inner chunks kept as they are, one after another in the
                    // current file, copied together
                    int64_t q = p+1;
                    while(q < last && sh.ready[q] == 3 &&
                          sh.oldIndex[2*q] == sh.oldIndex[2*(q-1)]+sh.oldIndex[(2*(q-1))+1]) q++;
                    ok = copyOldInner(sh, p, q, buf);
                    p = q;
                }
            }
            std::lock_guard<std::mutex> g(sh.m);
            if(!ok){
                sh.failed = true;
                sh.appending = false;
                return;
            }
            sh.next.store(last);
        }
    };
    // Hand a prepared inner chunk to its shard: added to its file right away when it is
    // next and no other thread is adding to it, else left for the appender (copied,
    // unless it is empty). Then, unless another thread is at it, this thread adds what
    // is ready after it.
    auto depositInShard = [&](ShardOut &sh, const int64_t pos, const StoredChunk &c, std::vector<uint8_t> &buf){
        bool direct = false;
        {
            std::lock_guard<std::mutex> g(sh.m);
            if(err || sh.failed) return;
            if(!sh.appending && pos == sh.next.load()) direct = true;
            else{
                sh.ready[pos] = c.empty ? 1 : 2;
                if(!c.empty){
                    std::vector<uint8_t> &w = sh.waiting[pos];
                    w.assign(c.data, c.data+c.n);
                    if(c.hasCrc){
                        const uint8_t* p = reinterpret_cast<const uint8_t*>(&c.crc);
                        w.insert(w.end(), p, p+sizeof(c.crc));
                    }
                }
                if(sh.appending || !sh.ready[sh.next.load()]) return;
            }
            sh.appending = true;
        }
        if(direct){
            const bool ok = c.empty || appendToShard(sh, pos, c.data, c.n, c.hasCrc ? &c.crc : nullptr);
            std::lock_guard<std::mutex> g(sh.m);
            if(!ok){
                sh.failed = true;
                sh.appending = false;
                return;
            }
            sh.next.store(pos+1);
        }
        appendReady(sh, buf);
    };
    // After a shard is started: add its first inner chunks if they need no preparing
    auto kickShard = [&](ShardOut &sh, std::vector<uint8_t> &buf){
        {
            std::lock_guard<std::mutex> g(sh.m);
            const int64_t next = sh.next.load();
            if(sh.appending || err || sh.failed || next >= perShard || !sh.ready[next]) return;
            sh.appending = true;
        }
        appendReady(sh, buf);
    };

    WorkRanges work(team, sharded ? 0 : numChunks);
    // (in a region of all of the threads, those past the team doing nothing: OpenMP ends
    // the threads a smaller region leaves out, and starts them again for the next bigger
    // one, every call here, after the region of all of them above)
    #pragma omp parallel num_threads(numWorkers)
    if(omp_get_thread_num() < team)
    {
        const LocalMemoryPolicy memoryPolicy;
        const int w = omp_get_thread_num();
        // This thread's buffers and zstd encoder, from the pool
        std::unique_ptr<WriteScratch> scratch;
        uint8_t* chunkUnC = nullptr;
        uint8_t* chunkC = nullptr;
        ZSTD_CCtx* zctx = nullptr;
        std::vector<uint8_t> stored, copyBuf;
        try{
            scratch = writeScratchPool.take(w);
            scratch->buffers(sB, chunkCCap);
            chunkUnC = scratch->unc.get();
            chunkC = scratch->comp.get();
            if(compressor == "zstd"){
                if(!scratch->zctx) scratch->zctx = ZSTD_createCCtx();
                zctx = scratch->zctx;
            }
            stored.swap(scratch->stored);
            copyBuf.swap(scratch->copyBuf);
        }
        catch(const std::exception &e){
            reportError("Not enough memory for a chunk of "+std::to_string(sB)+" bytes\n");
        }
        std::vector<uint64_t> cAV;
        std::vector<uint64_t> boxLo(nDims), boxHi(nDims), boxExt(nDims), inArray(nDims);
        const std::vector<uint64_t> zeros(nDims, 0);

        // Make chunk f ready to store: assembled (in its stored axis order) from the
        // input, and for a crop write from the chunk's current contents where the
        // write only partly covers it; then compressed. sh and pos: its shard and
        // place in it, for a sharded array.
        auto prepare = [&](const int64_t f, ShardOut* sh, const int64_t pos, StoredChunk &out) -> bool {
            const std::string &name = Zarr.get_chunkNames(f);
            cAV = Zarr.get_chunkAxisVals(name);
            // When cropping into an existing array, a chunk the written region only
            // partly covers (inside the array) needs its current contents
            bool partial = false;
            for(uint64_t d = 0; d < nDims; d++){
                partial = partial || cAV[d]*chunkDims[d] < startCoords[d] ||
                          ((cAV[d]+1)*chunkDims[d] > endCoords[d] && endCoords[d] < Zarr.get_shape(d));
            }
            const std::string subfolderName = Zarr.get_subfoldersString(cAV);
            // (chunkUnC already holds the chunk's current contents)
            bool existing = false;
            if(crop && partial && sh){
                // In a shard: an inner chunk the write partly covers is decoded from its
                // stored bytes (a missing one is the fill value)
                const bool isStored = !(sh->oldIndex[2*pos] == emptyEntry && sh->oldIndex[(2*pos)+1] == emptyEntry);
                if(isStored){
                    std::string why;
                    if(!readOldInner(*sh, pos, stored)) return false;
                    if(decodeChunk(compressor, chunkChecksum, stored.data(), stored.size(), chunkUnC, sB, why)){
                        reportError("An inner chunk of the shard "+sh->path+" "+why+"\n");
                        return false;
                    }
                }
                else fillElements(bytes, chunkUnC, s, typedFill);
                existing = true;
            }
            else if(crop && partial){
                // The chunk's current contents, decoded (the fill value if it has none),
                // and the fill value past the array's edge
                const std::string path = Zarr.get_fileName()+"/"+subfolderName+"/"+Zarr.chunkKey(name);
                std::ifstream file(path, std::ios::binary);
                if(file.is_open()){
                    file.seekg(0, std::ios::end);
                    const std::streamoff len = file.tellg();
                    stored.resize((size_t)len);
                    file.seekg(0, std::ios::beg);
                    file.read(reinterpret_cast<char*>(stored.data()), len);
                    file.close();
                    std::string why;
                    if(!file || decodeChunk(compressor, chunkChecksum, stored.data(), stored.size(), chunkUnC, sB, why)){
                        reportError("Error in Writer Read. Chunk: "+name+(why.empty() ? std::string() : " ("+why+")")+"\n");
                        return false;
                    }
                }
                else fillElements(bytes, chunkUnC, s, typedFill);
                for(uint64_t d = 0; d < nDims; d++){
                    const uint64_t b = cAV[d]*chunkDims[d];
                    inArray[d] = Zarr.get_shape(d) > b ? std::min(chunkDims[d], Zarr.get_shape(d)-b) : 0;
                }
                fillOutsideBoxND(bytes, chunkUnC, chunkDims, chunkStrides, zeros, inArray, fillElem);
                existing = true;
            }
            // Elements inside the written region come from the input; outside it, from
            // the chunk's existing contents when cropping into an existing array, and
            // otherwise the fill value
            uint64_t srcOff = 0, dstOff = 0;
            bool fullChunk = true;
            for(uint64_t d = 0; d < nDims; d++){
                const uint64_t b = cAV[d]*chunkDims[d];
                const uint64_t lo = std::max(b, startCoords[d]);
                const uint64_t hi = std::min(b+chunkDims[d], endCoords[d]);
                boxLo[d] = lo-b;
                boxHi[d] = hi > lo ? hi-b : lo-b;
                boxExt[d] = boxHi[d]-boxLo[d];
                fullChunk = fullChunk && boxLo[d] == 0 && boxHi[d] == chunkDims[d];
                srcOff += (lo-startCoords[d])*inStrides[d];
                dstOff += boxLo[d]*chunkStrides[d];
            }
            if(!fullChunk && !existing) fillOutsideBoxND(bytes, chunkUnC, chunkDims, chunkStrides, boxLo, boxHi, fillElem);
            copyBoxND(bytes, (const uint8_t*)zarrArr+srcOff*bytes, chunkUnC+dstOff*bytes, boxExt, inStrides, chunkStrides);

            // An all-zero chunk is left unwritten (it reads back as zeros)
            if(zeroChunkUnc && !memcmp(zeroChunkUnc,chunkUnC,sB)){
                out.zero = true;
                return true;
            }

            // Compress as the array's metadata says
            int64_t csize = 0;
            // Uncompressed arrays store the chunk as is
            const uint8_t* chunkOut = chunkC;
            if(compressor == "none"){
                chunkOut = chunkUnC;
                csize = sB;
            }
            else if(compressor == "blosc"){
                csize = blosc_compress_ctx(Zarr.get_clevel(), bloscShuffle, bytes, sB, chunkUnC, chunkC, sB+BLOSC_MAX_OVERHEAD,Zarr.get_cname().c_str(),0,nBloscThreads);
                // A non-positive return means blosc could not compress the chunk
                // (e.g. it exceeds BLOSC_MAX_BUFFERSIZE). Fail loudly instead of
                // writing an empty/garbage chunk. (The chunk-size guard above
                // should already have caught the oversize case.)
                if(csize <= 0){
                    reportError("Compression error (blosc returned "+std::to_string(csize)+
                        "). ChunkName: "+Zarr.get_fileName()+"/"+subfolderName+"/"+name+"\n");
                    return false;
                }
            }
            else if(compressor == "zstd"){
                const size_t n = ZSTD_compressCCtx(zctx, chunkC, chunkCCap, chunkUnC, sB, (int)Zarr.get_clevel());
                if(ZSTD_isError(n)){
                    reportError("Compression error (zstd: "+std::string(ZSTD_getErrorName(n))+
                        "). ChunkName: "+Zarr.get_fileName()+"/"+subfolderName+"/"+name+"\n");
                    return false;
                }
                csize = (int64_t)n;
            }
            else{
                // gzip, or a zlib stream (numcodecs zlib)
                z_stream stream;
                memset(&stream, 0, sizeof(stream));
                stream.next_in = chunkUnC;
                stream.next_out = chunkC;
                stream.avail_in = sB;
                stream.avail_out = chunkCCap;
                int cErr = deflateInit2(&stream, Zarr.get_clevel(), Z_DEFLATED, zlibStream ? MAX_WBITS : MAX_WBITS + 16, MAX_MEM_LEVEL, Z_DEFAULT_STRATEGY);
                bool deflated = cErr == Z_OK;
                if(deflated){
                    cErr = deflate(&stream, Z_FINISH);
                    deflated = cErr == Z_STREAM_END;
                    const int endErr = deflateEnd(&stream);
                    if(deflated && endErr != Z_OK){
                        deflated = false;
                        cErr = endErr;
                    }
                }
                if(!deflated){
                    reportError("Compression error. Error code: "+std::to_string(cErr)+" ChunkName: "+
                                Zarr.get_fileName()+"/"+subfolderName+"/"+name+"\n");
                    return false;
                }
                csize = chunkCCap-stream.avail_out;
            }
            out.empty = false;
            out.data = chunkOut;
            out.n = (uint64_t)csize;
            // A crc32c codec's checksum follows the compressed chunk
            if(chunkChecksum){
                out.hasCrc = true;
                out.crc = crc32c(chunkOut, (size_t)csize);
            }
            return true;
        };

        try{
            if(!sharded){
                // Chunk files: contiguous blocks of chunks with work stealing
                for(int64_t f = (err || !chunkUnC) ? -1 : work.next(w); f >= 0; f = err ? -1 : work.next(w)){
                    StoredChunk c;
                    if(!prepare(f, nullptr, 0, c)) break;
                    const std::string &name = Zarr.get_chunkNames(f);
                    const std::string fileName(Zarr.get_fileName()+"/"+Zarr.get_subfoldersString(cAV)+"/"+Zarr.chunkKey(name));
                    // An all-zero chunk is not written, and any earlier copy of it removed
                    if(c.zero){
                        remove(fileName.c_str());
                        continue;
                    }
                    const std::string tmpName = useUuid ? fileName+uuid : fileName;
                    std::ofstream file;
                    if(!openForWrite(file, tmpName, std::ios::trunc)){
                        reportError("Check permissions or filepath. Cannot write to path: "+tmpName+"\n");
                        break;
                    }
                    file.write(reinterpret_cast<const char*>(c.data), c.n);
                    if(c.hasCrc) file.write(reinterpret_cast<const char*>(&c.crc), sizeof(c.crc));
                    file.close();
                    if(!file){
                        if(useUuid) remove(tmpName.c_str());
                        reportError("Could not write all of "+tmpName+" (is the disk full?)\n");
                        break;
                    }
                    if(useUuid && !renameReplace(tmpName, fileName)){
                        remove(tmpName.c_str());
                        reportError("Cannot move the temporary file into place: "+fileName+"\n");
                        break;
                    }
                }
            }
            else if(chunkUnC){
                // A shard nobody has started, else the one with the most inner chunks
                // left to take (of those with as many, the first from this thread's
                // place on, so helpers spread out); then its inner chunks, one at a time
                while(!err){
                    int64_t s = nextShard.fetch_add(1);
                    if(s >= numShards){
                        s = -1;
                        int64_t most = 0;
                        for(int64_t j = 0; j < numShards; j++){
                            const int64_t i = (w+j)%numShards;
                            const int64_t left = shards[i].prepEnd-shards[i].prepNext.load();
                            if(left > most){
                                most = left;
                                s = i;
                            }
                        }
                        if(s < 0) break;
                    }
                    ShardOut &sh = shards[s];
                    while(!err){
                        const int64_t k = sh.prepNext.fetch_add(1);
                        if(k >= sh.prepEnd) break;
                        const int64_t f = toPrepare[k];
                        const int64_t pos = f%perShard;
                        if(!startShard(sh, f)) break;
                        kickShard(sh, copyBuf);
                        // (not too far ahead of what the shard has added to its file)
                        for(int us = 50; !err && sh.preparedBefore[pos]-sh.preparedBefore[sh.next.load()] > window; us = std::min(2*us, 1000)){
                            std::this_thread::sleep_for(std::chrono::microseconds(us));
                        }
                        StoredChunk c;
                        if(!prepare(f, &sh, pos, c)) break;
                        if(c.zero) c = StoredChunk();
                        depositInShard(sh, pos, c, copyBuf);
                    }
                }
            }
        }
        catch(const std::exception &e){
            reportError(std::string("Write error: ")+e.what()+"\n");
        }
        try{
            if(scratch){
                stored.swap(scratch->stored);
                copyBuf.swap(scratch->copyBuf);
            }
            writeScratchPool.give(w, std::move(scratch));
        }
        catch(const std::exception&){}
    }
    free(zeroChunkUnc);

    // After an error, a shard that was not finished is not moved into place (its
    // temporary file is removed)
    for(int64_t i = 0; i < numShards; i++){
        ShardOut &sh = shards[i];
        if(!sh.started || sh.done) continue;
        if(sh.old.is_open()) sh.old.close();
        if(sh.opened){
            sh.file.close();
            if(shardTmpFile) remove(sh.writePath.c_str());
        }
    }

    if(err) {
        Zarr.set_errString(errString);
        return 1;
    }
    else return 0;
}
