#ifndef ZARR_H
#define ZARR_H
#include <cstdint>
#include <vector>
#include <string>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

class zarr
{
public:
    zarr();
    zarr(const std::string &fileName);
    zarr(const std::string &fileName, const std::vector<uint64_t> &chunks,
         uint64_t blocksize, uint64_t clevel, const std::string &cname,
         const std::string &id, uint64_t shuffle, const std::string &dimension_separator, const std::string &dtype,
         const std::string &fill_value, const std::vector<std::string> &filters,
         const std::string &order, const std::vector<uint64_t> &shape,
         uint64_t zarr_format, const std::vector<uint64_t> &subfolders,
		 const bool shard, const std::vector<uint64_t> &chunk_shape);
    ~zarr();
    void write_zarray();
    const std::string &get_fileName() const;
    void set_fileName(const std::string &fileName);
    const uint64_t &get_chunks(const uint64_t &index) const;
    void set_chunks(const std::vector<uint64_t> &chunks);
    const uint64_t &get_clevel() const;
    void set_clevel(const uint64_t &clevel);
    const std::string &get_cname() const;
    void set_cname(const std::string &cname);
    const std::string &get_dimension_separator() const;
    void set_dimension_separator(const std::string &dimension_separator);
    const std::string &get_dtype() const;
    void set_dtype(const std::string &dtype);
    const std::string &get_fill_value() const;
    void set_fill_value(const std::string &fill_value);
    void set_fill_value(const int64_t &fill_value);
    const std::string &get_order() const;
    void set_order(const std::string &order);
    const uint64_t &get_shape(const uint64_t &index) const;
    void set_shape(const std::vector<uint64_t> &shape);
    const std::string get_subfoldersString(const std::vector<uint64_t> &cAV) const;
    void set_subfolders(const std::vector<uint64_t> &subfolders);

    const bool &get_shard() const;
    void set_shard(const bool shard);
    const uint64_t &get_chunk_shape(const uint64_t &index) const;
    void set_chunk_shape(const std::vector<uint64_t> &chunk_shape);
    const uint64_t &get_numShards() const;
    const uint64_t &get_numChunksPerShard() const;
    const std::string chunkNameToShardName(const std::string &chunkName) const;

    const std::vector<uint64_t> chunkToShard(const std::vector<uint64_t> &cAV) const;
    // Shard index layout (Zarr v3 sharding_indexed index_location and
    // index_codecs): at the start of each shard instead of the end (the
    // default), and whether a CRC32C checksum follows it (the default)
    bool get_shardIndexAtStart() const;
    bool get_shardIndexChecksum() const;
    const uint64_t get_ShardPosition(const std::vector<uint64_t> &cAV) const;

    const uint64_t get_chunkShardPosition(const std::vector<uint64_t> &cAV) const;

    // Zarr format of the array's metadata: 2 (.zarray) or 3 (zarr.json)
    uint64_t get_zarr_format() const;
    // Write Zarr v3 metadata (zarr.json) instead of v2 (.zarray) from
    // write_zarray(). Also sets v3's defaults, C order and "/" in chunk names
    // ("c/0/0/0"); set_order or set_dimension_separator after it changes them.
    void set_zarr_format(const uint64_t zarr_format);
    // Blosc shuffle: 0 none, 1 byte, 2 bit
    uint64_t get_shuffle() const;
    // Turn the v2 array this object was opened from into a Zarr v3 array by
    // writing a zarr.json that describes its existing chunk files (no data is
    // rewritten) and removing its .zarray. Throws "v3Unsupported:why" for what
    // v3 cannot describe (subfolders, zlib, shard shapes that are not a multiple
    // of the inner chunk shape).
    void convert_to_v3();
    // Path of a chunk (or shard) file inside the array's folder: its name, after
    // the "c" prefix of Zarr v3's default chunk key encoding
    const std::string chunkKey(const std::string &chunkName) const;
    // Axes of a stored chunk from slowest to fastest: 0..n-1 for C order,
    // n-1..0 for F order, any permutation for a Zarr v3 transpose codec
    const std::vector<uint64_t> get_chunkAxisOrder() const;
    // How each stored chunk is compressed: "none", "blosc", "gzip" or "zstd"
    // (Zarr v3 zstd codec, or numcodecs Zstd in v2)
    const std::string get_compressor() const;
    // Whether a CRC32C checksum follows each stored chunk (Zarr v3 crc32c codec)
    bool get_chunkChecksum() const;

    const std::vector<uint64_t> get_chunkAxisVals(const std::string &fileName) const;
    void set_chunkInfo(const std::vector<uint64_t> &startCoords,
                             const std::vector<uint64_t> &endCoords);
    const std::string &get_chunkNames(const uint64_t &index) const;
    const uint64_t &get_numChunks() const;
    const std::string &get_errString() const;
    void set_errString(const std::string &errString);
    const uint64_t dtypeBytes() const;
    // Number of dimensions of the array (length of shape)
    uint64_t get_ndims() const;
    // Fit chunks, subfolders and chunk_shape to the number of dimensions: extra
    // trailing values (e.g. the 3D defaults on a 2D array) are dropped, missing
    // ones get chunks 256 on the first three axes and 1 after, subfolders 0 and
    // chunk_shape 1.
    void normalizeDims();
private:
    void parseZarrJson();
    json v3Metadata(const bool existingChunks) const;
    void set_jsonValues();
    void write_jsonValues();
    void set_shardData();
    uint64_t fastCeilDiv(uint64_t num, uint64_t denom);
    void createSubfolders();
    std::string fileName;
    json zarray;
    std::vector<uint64_t> chunks;
    uint64_t blocksize;
    uint64_t clevel;
    std::string cname;
    std::string id;
    uint64_t shuffle;
    std::string dimension_separator;
    std::string dtype;
    std::string fill_value;
    std::vector<std::string> filters;
    std::string order;
    std::vector<uint64_t> shape;
    uint64_t zarr_format;
    std::vector<uint64_t> subfolders;
    
    // Sharding
    bool shard;
    std::vector<uint64_t> chunk_shape;
    std::vector<uint64_t> chunksPerShard;
    std::vector<uint64_t> shards;
    uint64_t numShards;
    uint64_t numChunksPerShard;

    std::vector<std::string> chunkNames;
    uint64_t numChunks;

    std::string errString;
};
#endif
