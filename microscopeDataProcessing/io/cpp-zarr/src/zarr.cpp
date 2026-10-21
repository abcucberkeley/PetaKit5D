#include <cmath>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <omp.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <stdarg.h>
#include <sys/time.h>
#endif
#include <fstream>
#include "zarr.h"
#include "helperfunctions.h"

// A metadata fill_value as the zarr class keeps it: integers in decimal, other
// numbers in full precision, "NaN", "Infinity", "-Infinity" and Zarr v3's "0x..."
// bit patterns as they are, a complex number's [real, imaginary] pair as its JSON
// text, and "0" for no fill value (null)
static std::string fillValueString(const json &v){
    if(v.is_null()) return "0";
    if(v.is_boolean()) return v.get<bool>() ? "1" : "0";
    if(v.is_number_unsigned()) return std::to_string(v.get<uint64_t>());
    if(v.is_number_integer()) return std::to_string(v.get<int64_t>());
    if(v.is_number_float()){
        const double d = v.get<double>();
        if(std::isfinite(d) && d == std::floor(d) && std::fabs(d) < 9e18) return std::to_string((int64_t)d);
        char buf[32];
        snprintf(buf, sizeof(buf), "%.17g", d);
        return buf;
    }
    if(v.is_string()){
        const std::string s = v.get<std::string>();
        return s.empty() || s == "null" ? "0" : s;
    }
    if(v.is_array() && v.size() == 2) return v.dump();
    throw std::string("metadataIncomplete");
}

// One float of a fill value as JSON writes it: a number, or "NaN", "Infinity",
// "-Infinity" as strings. A Zarr v2 .zarray has no "0x..." bit patterns: those
// become the value they stand for (bytes: 2, 4 or 8; bfloat16 for 2-byte bfloat16).
static json floatFillJson(const std::string &f, const uint64_t bytes, const bool bfloat16, const bool bitsAllowed){
    if(f == "NaN" || f == "Infinity" || f == "-Infinity") return f;
    const bool bits = f.size() > 2 && f[0] == '0' && (f[1] == 'x' || f[1] == 'X');
    if(bits && bitsAllowed) return f;
    double v = 0;
    if(bits){
        uint8_t elem[16];
        const std::string dtype = bytes == 2 ? (bfloat16 ? "bfloat16" : "<f2") : bytes == 4 ? "<f4" : "<f8";
        fillValueElement(f, isLittleEndian() ? dtype : (bfloat16 ? ">bfloat16" : ">"+dtype.substr(1)), elem);
        if(bytes == 8) memcpy(&v, elem, 8);
        else if(bytes == 4){ float x; memcpy(&x, elem, 4); v = x; }
        else{ uint16_t h; memcpy(&h, elem, 2); v = bfloat16 ? bfloat16ToFloat(h) : halfToFloat(h); }
        if(std::isnan(v)) return "NaN";
        if(std::isinf(v)) return v > 0 ? "Infinity" : "-Infinity";
    }
    else v = strtod(f.c_str(), NULL);
    // (whole numbers as integers, as cpp-zarr has always written them)
    if(std::isfinite(v) && v == std::floor(v) && std::fabs(v) < 9e15 && !(v == 0 && std::signbit(v))) return (int64_t)v;
    return v;
}

// A complex number's fill value as JSON writes it: [real, imaginary]
static json complexFillJson(const std::string &f, const uint64_t bytes, const bool bitsAllowed){
    json parts = json::array();
    try{
        const json v = json::parse(f);
        if(v.is_array() && v.size() == 2){
            for(const json &p : v) parts.push_back(p.is_string() ? p.get<std::string>() : fillValueString(p));
        }
    }
    catch(...){}
    if(parts.empty()) parts = json::array({f, "0"});
    json out = json::array();
    for(const json &p : parts) out.push_back(floatFillJson(p.get<std::string>(), bytes/2, false, bitsAllowed));
    return out;
}

// ---- Zarr v3 metadata (zarr.json) ---------------------------------------------

namespace {

[[noreturn]] void unsupportedV3(const std::string &why){
    throw std::string("metadataUnsupported:"+why);
}

// A Zarr v3 extension point (codec, chunk grid, chunk key encoding): an object
// {"name": ..., "configuration": {...}}, or just its name as a string
std::string extName(const json &e){
    if(e.is_string()) return e.get<std::string>();
    if(e.is_object()){
        const auto n = e.find("name");
        if(n != e.end() && n->is_string()) return n->get<std::string>();
    }
    throw std::string("metadataIncomplete");
}

json extConfig(const json &e){
    if(e.is_object()){
        const auto c = e.find("configuration");
        if(c != e.end() && c->is_object()) return *c;
    }
    return json::object();
}

// What a list of Zarr v3 codecs does to a chunk: array -> array codecs
// (transpose), one array -> bytes codec (bytes, or sharding_indexed), then
// bytes -> bytes codecs (one compressor, then a crc32c checksum)
struct V3Codecs {
    std::vector<uint64_t> axisOrder;   // transpose: stored axes, slowest first (empty: C order)
    bool bigEndian = false;           // bytes
    bool sharded = false;             // sharding_indexed instead of bytes...
    json sharding;                    // ...and its configuration
    std::string compressor = "none";  // blosc, gzip or zstd
    json compressorConfig;
    bool checksum = false;            // crc32c
};

V3Codecs parseV3Codecs(const json &codecs, const uint64_t nDims, const bool allowSharding){
    if(!codecs.is_array() || codecs.empty()) unsupportedV3("its list of codecs is empty");
    V3Codecs c;
    size_t i = 0;
    for(; i < codecs.size() && extName(codecs[i]) == "transpose"; i++){
        if(!c.axisOrder.empty()) unsupportedV3("it has more than one transpose codec");
        const json order = extConfig(codecs[i]).value("order", json());
        std::vector<uint64_t> axes(nDims);
        if(order.is_string() && (order == "C" || order == "F")){
            for(uint64_t d = 0; d < nDims; d++) axes[d] = order == "C" ? d : nDims-1-d;
        }
        else if(order.is_array() && order.size() == nDims){
            std::vector<bool> seen(nDims, false);
            for(uint64_t d = 0; d < nDims; d++){
                if(!order[d].is_number_integer() || order[d].get<int64_t>() < 0 ||
                   order[d].get<uint64_t>() >= nDims || seen[order[d].get<uint64_t>()]){
                    unsupportedV3("its transpose order is not a permutation of the axes");
                }
                axes[d] = order[d].get<uint64_t>();
                seen[axes[d]] = true;
            }
        }
        else unsupportedV3("its transpose order is not a permutation of the axes");
        c.axisOrder = axes;
    }
    if(i == codecs.size()) unsupportedV3("it has no bytes codec");
    const std::string arrayToBytes = extName(codecs[i]);
    if(arrayToBytes == "bytes"){
        const std::string endian = extConfig(codecs[i]).value("endian", "little");
        if(endian != "little" && endian != "big") unsupportedV3("its bytes codec has endian \""+endian+"\"");
        c.bigEndian = endian == "big";
    }
    else if(arrayToBytes == "sharding_indexed"){
        if(!allowSharding) unsupportedV3("it has shards inside shards");
        if(!c.axisOrder.empty()) unsupportedV3("it has a transpose codec outside its shards");
        c.sharded = true;
        c.sharding = extConfig(codecs[i]);
    }
    else unsupportedV3("the \""+arrayToBytes+"\" codec is not supported");
    for(i++; i < codecs.size(); i++){
        const std::string name = extName(codecs[i]);
        if(c.sharded) unsupportedV3("it has codecs after sharding_indexed");
        if(c.checksum) unsupportedV3("it has codecs after crc32c");
        if(name == "crc32c") c.checksum = true;
        else if(name == "blosc" || name == "gzip" || name == "zstd"){
            if(c.compressor != "none") unsupportedV3("it has more than one compressor");
            c.compressor = name;
            c.compressorConfig = extConfig(codecs[i]);
        }
        else unsupportedV3("the \""+name+"\" codec is not supported");
    }
    return c;
}

// Whether zarray holds Zarr v3 metadata (an array opened from, or written as, a
// zarr.json): a new array set to Zarr v3 has v2-style metadata until written
bool isV3Json(const json &z){
    return z.is_object() && z.contains("codecs") && z.contains("chunk_key_encoding") &&
           z.value("zarr_format", json(0)) == 3;
}

// The codecs of each stored chunk: a sharded array's inner codecs
json chunkCodecsV3(const json &meta){
    const json &codecs = meta.at("codecs");
    if(codecs.is_array() && !codecs.empty() && extName(codecs[0]) == "sharding_indexed"){
        return extConfig(codecs[0]).at("codecs");
    }
    return codecs;
}

// The fields of Zarr v3 array metadata; others must be understood, unless they
// say otherwise (must_understand: false)
const char* const v3Fields[] = {"zarr_format", "node_type", "shape", "data_type", "chunk_grid",
                                "chunk_key_encoding", "fill_value", "codecs", "attributes",
                                "dimension_names", "storage_transformers"};

} // namespace

// Zarr v3 array metadata (zarr.json, already in zarray) into the fields a v2
// .zarray fills. Anything cpp-zarr cannot read throws "metadataUnsupported:why".
void zarr::parseZarrJson(){
    try{
        const json &m = zarray;
        if(!m.is_object() || !m.contains("zarr_format") || m.at("zarr_format") != 3) throw std::string("metadataIncomplete");
        const std::string nodeType = m.value("node_type", "");
        if(nodeType == "group") unsupportedV3("it is a Zarr v3 group, not an array");
        if(nodeType != "array") throw std::string("metadataIncomplete");
        for(auto it = m.begin(); it != m.end(); ++it){
            bool known = false;
            for(const char* f : v3Fields) known = known || it.key() == f;
            if(!known && !(it->is_object() && it->value("must_understand", true) == false)){
                unsupportedV3("its metadata has a field cpp-zarr does not know (\""+it.key()+"\")");
            }
        }
        const auto transformers = m.find("storage_transformers");
        if(transformers != m.end() && !(transformers->is_array() && transformers->empty())){
            unsupportedV3("it has storage transformers");
        }

        shape = m.at("shape").get<std::vector<uint64_t>>();
        const uint64_t n = shape.size();
        const json &grid = m.at("chunk_grid");
        if(extName(grid) != "regular") unsupportedV3("its chunk grid is not regular");
        chunks = extConfig(grid).at("chunk_shape").get<std::vector<uint64_t>>();
        if(chunks.size() != n) throw std::string("metadataIncomplete");
        for(uint64_t c : chunks) if(!c) throw std::string("metadataIncomplete");

        // Chunk keys: "c/0/0/0" (default) or "0.0.0" (v2), with either separator
        const json &keys = m.at("chunk_key_encoding");
        const std::string encoding = extName(keys);
        if(encoding != "default" && encoding != "v2") unsupportedV3("its chunk key encoding \""+encoding+"\" is not supported");
        dimension_separator = extConfig(keys).value("separator", encoding == "default" ? "/" : ".");
        if(dimension_separator != "/" && dimension_separator != ".") unsupportedV3("its chunk key separator is not \"/\" or \".\"");

        // The data type (before the codecs: a type cpp-zarr does not read, such as
        // strings, often comes with codecs it does not know either)
        {
            const json &dataType = m.at("data_type");
            const char* const known[] = {"bool", "int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64",
                                         "float16", "float32", "float64", "complex64", "complex128", "bfloat16"};
            bool ok = false;
            if(dataType.is_string()) for(const char* k : known) ok = ok || dataType.get<std::string>() == k;
            if(!ok) unsupportedV3("its data type "+(dataType.is_string() ? "\""+dataType.get<std::string>()+"\" " : std::string())+
                                  "is not supported");
        }

        // Codecs, and for a sharded array its inner chunks and shard index
        const V3Codecs outer = parseV3Codecs(m.at("codecs"), n, true);
        V3Codecs inner = outer;
        shard = outer.sharded;
        chunk_shape.clear();
        if(shard){
            if(n == 0) unsupportedV3("sharding is not supported for 0-dimensional arrays");
            chunk_shape = outer.sharding.at("chunk_shape").get<std::vector<uint64_t>>();
            if(chunk_shape.size() != n) throw std::string("metadataIncomplete");
            for(uint64_t d = 0; d < n; d++){
                if(!chunk_shape[d] || chunks[d]%chunk_shape[d]) unsupportedV3("its inner chunk shape does not divide its shard shape");
            }
            inner = parseV3Codecs(outer.sharding.at("codecs"), n, false);
            const json indexCodecs = outer.sharding.value("index_codecs", json::array({{{"name", "bytes"}, {"configuration", {{"endian", "little"}}}}, {{"name", "crc32c"}}}));
            bool indexOk = indexCodecs.is_array() && !indexCodecs.empty() && indexCodecs.size() <= 2 &&
                           extName(indexCodecs[0]) == "bytes" && extConfig(indexCodecs[0]).value("endian", "little") == "little";
            if(indexOk && indexCodecs.size() == 2) indexOk = extName(indexCodecs[1]) == "crc32c";
            if(!indexOk) unsupportedV3("its shard index codecs are not little-endian bytes and an optional crc32c");
            const std::string location = outer.sharding.value("index_location", "end");
            if(location != "start" && location != "end") unsupportedV3("its shard index location is not \"start\" or \"end\"");
        }

        // Data type, and the byte order of the stored chunks
        const json &dataType = m.at("data_type");
        if(!dataType.is_string()) unsupportedV3("its data type is not supported");
        const std::string t = dataType.get<std::string>();
        char kind = 0;
        int size = 0;
        if(t == "bool"){ kind = 'b'; size = 1; }
        else if(t == "int8" || t == "int16" || t == "int32" || t == "int64"){ kind = 'i'; size = std::stoi(t.substr(3))/8; }
        else if(t == "uint8" || t == "uint16" || t == "uint32" || t == "uint64"){ kind = 'u'; size = std::stoi(t.substr(4))/8; }
        else if(t == "float16" || t == "float32" || t == "float64"){ kind = 'f'; size = std::stoi(t.substr(5))/8; }
        else if(t == "complex64" || t == "complex128"){ kind = 'c'; size = std::stoi(t.substr(7))/8; }
        // (bfloat16: a Zarr extension data type)
        else if(t == "bfloat16"){ kind = 'f'; size = 2; }
        else unsupportedV3("its data type \""+t+"\" is not supported");
        if(t == "bfloat16") dtype = std::string(inner.bigEndian ? ">" : "")+"bfloat16";
        else dtype = std::string(size == 1 ? "|" : inner.bigEndian ? ">" : "<")+kind+std::to_string(size);

        // Each chunk's storage order: C, F, or another transpose
        order = "C";
        if(!inner.axisOrder.empty() && n > 1){
            bool reversed = true;
            for(uint64_t d = 0; d < n; d++) reversed = reversed && inner.axisOrder[d] == n-1-d;
            if(reversed) order = "F";
        }

        // Compressor, as the v2 fields name it
        const json &cc = inner.compressorConfig;
        cname = "none"; id = ""; clevel = 0; shuffle = 0; blocksize = 0;
        if(inner.compressor == "blosc"){
            cname = cc.at("cname").get<std::string>();
            id = "blosc";
            clevel = cc.value("clevel", 5);
            blocksize = cc.value("blocksize", 0);
            const json sh = cc.value("shuffle", json("noshuffle"));
            if(sh.is_string()) shuffle = sh == "shuffle" ? 1 : sh == "bitshuffle" ? 2 : 0;
            else shuffle = sh.get<uint64_t>();
        }
        else if(inner.compressor == "gzip"){
            cname = "gzip";
            id = "gzip";
            clevel = cc.value("level", 1);
        }
        else if(inner.compressor == "zstd"){
            cname = "zstd";
            id = "zstd";
            // (zstd levels can be negative: kept as two's complement, as a v2
            // .zarray's numcodecs Zstd level is)
            clevel = (uint64_t)cc.value("level", (int64_t)0);
        }

        // Fill value: a number, true/false, for floats "NaN", "Infinity",
        // "-Infinity" or a "0x..." bit pattern, and for complex numbers a pair of
        // those, [real, imaginary]
        const json &fv = m.at("fill_value");
        // (one float of the fill value, of partBytes bytes)
        auto floatFill = [](const json &v, const int partBytes){
            if(v.is_number()) return true;
            if(!v.is_string()) return false;
            const std::string f = v.get<std::string>();
            const bool bits = f.size() == 2+2*(size_t)partBytes && f[0] == '0' && (f[1] == 'x' || f[1] == 'X');
            return f == "NaN" || f == "Infinity" || f == "-Infinity" || bits;
        };
        if(kind == 'b'){
            if(!fv.is_boolean()) throw std::string("metadataIncomplete");
            fill_value = fv.get<bool>() ? "1" : "0";
        }
        else if(kind == 'c'){
            if(!fv.is_array() || fv.size() != 2 || !floatFill(fv[0], size/2) || !floatFill(fv[1], size/2)){
                throw std::string("metadataIncomplete");
            }
            fill_value = fv.dump();
        }
        else if(fv.is_string()){
            if(kind != 'f' || !floatFill(fv, size)) throw std::string("metadataIncomplete");
            fill_value = fv.get<std::string>();
        }
        else if(fv.is_number()) fill_value = fillValueString(fv);
        else throw std::string("metadataIncomplete");

        zarr_format = 3;
        subfolders.assign(n, 0);
        normalizeDims();
    }
    catch(const std::string &){
        throw;
    }
    catch(...){
        throw std::string("metadataIncomplete");
    }
}

// Create a blank zarr object with default values. New arrays use '/' chunk keys
// (nested folders): on network file systems, creating or replacing files in one
// folder happens one at a time per node, so arrays with many chunks in a single
// folder wrote up to 37x slower on NFS. Existing arrays keep their own separator,
// and a .zarray without one means '.'.
zarr::zarr() :
fileName(""), chunks({256,256,256}), blocksize(0),
clevel(5), cname("lz4"), id("blosc"), shuffle(1), dtype("<u2"),
dimension_separator("/"), fill_value("0"), filters({}), order("F"),
shape({0,0,0}), zarr_format(2), subfolders({0,0,0}), shard(false),
chunk_shape({1,1,1})
{
    set_jsonValues();
}

zarr::zarr(const std::string &fileName) :
fileName(fileName), chunks({256,256,256}), blocksize(0),
clevel(5), cname("lz4"), id("blosc"), shuffle(1), dtype("<u2"),
fill_value("0"), filters({}), order("F"), shape({0,0,0}),
zarr_format(2), subfolders({0,0,0}), shard(false), chunk_shape({1,1,1})
{
    // A v2 array (.zarray) is read as before; without one, a Zarr v3 array (zarr.json)
    if(!fileExists(fileName+"/.zarray")){
        if(fileExists(fileName+"/zarr.json")){
            try{
                std::ifstream f(fileName+"/zarr.json");
                zarray = json::parse(f);
            }
            catch(...){
                throw std::string("metadataIncomplete");
            }
            parseZarrJson();
            return;
        }
        throw std::string("metadataFileMissing:"+fileName);
        //mexErrMsgIdAndTxt("zarr:zarrayError","Metadata file in \"%s\" is missing. Does the file exist?",fileName.c_str());
    }
    try{
        std::ifstream f(fileName+"/.zarray");
        zarray = json::parse(f);
    }
    catch(...){
        throw std::string("metadataIncomplete");
    }

    try{
        // Check for Sharding
        // Sharding should be false by default
        try{
            // Assuming the shard info is at the beginning for now
            if(zarray.at("codecs").at(0).at("name") == "sharding_indexed"){
                chunk_shape = zarray.at("codecs").at(0).at("configuration").at("chunk_shape").get<std::vector<uint64_t>>();
                shard = true;
            }
        }
        catch(...){

        }
        chunks = zarray.at("chunks").get<std::vector<uint64_t>>();
        // A null compressor means the chunks are stored uncompressed (zarr-python
        // always does this for 0-dimensional arrays)
        if(zarray.at("compressor").is_null()){
            cname = "none";
            clevel = 0;
            blocksize = 0;
            id = "";
            shuffle = 0;
        }
        else try{
            // Try blosc compression types
            cname = zarray.at("compressor").at("cname");
            clevel = zarray.at("compressor").at("clevel");
            blocksize = zarray.at("compressor").at("blocksize");
            id = zarray.at("compressor").at("id");
            shuffle = zarray.at("compressor").at("shuffle");
        }
        catch(...){
            // Other compressors with a level: gzip, and the numcodecs zlib and zstd
            // codecs (read only)
            clevel = zarray.at("compressor").at("level");
            cname = zarray.at("compressor").at("id");
            blocksize = 0;
            id = cname;
            shuffle = 0;
        }
        // If dimension_separator does not exist then assume it is "."
        try{
            dimension_separator = zarray.at("dimension_separator");
            if(dimension_separator != "/" && dimension_separator != "."){
                throw std::string("metadataIncomplete:dimension_separatorIncorrect");
            }
        }
        catch(...){
            dimension_separator = ".";
        }

        dtype = zarray.at("dtype");
        if(!dtypeInfo(dtype).kind) throw std::string("metadataUnsupported:its data type \""+dtype+"\" is not supported");
        fill_value = fillValueString(zarray.at("fill_value"));
        // Filters (numcodecs Delta, Shuffle, ...) change the stored values, and
        // are not applied: such arrays cannot be read
        const auto filtersIt = zarray.find("filters");
        if(filtersIt != zarray.end() && !filtersIt->is_null() && !(filtersIt->is_array() && filtersIt->empty())){
            throw std::string("metadataUnsupported:it has filters, which are not supported");
        }
        //filters = "";

        order = zarray.at("order");
        shape = zarray.at("shape").get<std::vector<uint64_t>>();
        zarr_format = zarray.at("zarr_format");
    }
    catch(const std::string &e){
        // (a reason the array cannot be read)
        if(e.rfind("metadataUnsupported:", 0) == 0) throw;
        throw std::string("metadataIncomplete");
    }
    catch(...){
        throw std::string("metadataIncomplete");
        //mexErrMsgIdAndTxt("zarr:zarrayError","Metadata is incomplete. Check the .zarray file");
    }
    try{
        subfolders = zarray.at("subfolders").get<std::vector<uint64_t>>();
    }
    catch(...){
        subfolders = {0,0,0};
    }
}

// Create a new zarr file with .zarray metadata file (no actual data)
zarr::zarr(const std::string &fileName, const std::vector<uint64_t> &chunks,
           uint64_t blocksize, uint64_t clevel, const std::string &cname,
           const std::string &id, uint64_t shuffle, const std::string &dimension_separator, const std::string &dtype,
           const std::string &fill_value, const std::vector<std::string> &filters,
           const std::string &order, const std::vector<uint64_t> &shape,
           uint64_t zarr_format, const std::vector<uint64_t> &subfolders,
		   const bool shard, const std::vector<uint64_t> &chunk_shape) :
fileName(fileName), chunks(chunks), blocksize(blocksize), clevel(clevel),
cname(cname), id(id), shuffle(shuffle), dimension_separator(dimension_separator),
dtype(dtype), fill_value(fill_value), filters(filters), order(order), shape(shape),
zarr_format(zarr_format), subfolders(subfolders), shard(shard), chunk_shape(chunk_shape)
{
    // Handle the tilde character in filenames on Linux/Mac
    #ifndef _WIN32
    this->fileName = expandTilde(this->fileName.c_str());
    #endif
    set_jsonValues();
}

zarr::~zarr(){

}

// Write a metadata file (fileName/name) through a temporary file, making the
// array's folder if needed
static void writeJsonFile(const std::string &folder, const std::string &name, const json &meta){
    const std::string final(folder+"/"+name);
    if(!fileExists(final)) mkdirRecursive(folder.c_str());
    const std::string tmp(final+generateUUID());
    std::ofstream o(tmp);
    if(!o.good()) throw std::string("cannotOpenZarray:"+tmp);
    o << std::setw(4) << meta << std::endl;
    o.close();
    if(!o || !renameReplace(tmp, final)){
        remove(tmp.c_str());
        throw std::string("cannotOpenZarray:"+final);
    }
}

// Whether a file holds a JSON object (a .zattrs whose attributes can be moved)
static bool isJsonObjectFile(const std::string &path){
    try{
        std::ifstream f(path);
        return json::parse(f).is_object();
    }
    catch(...){
        return false;
    }
}

// Once a v2 array's folder has its zarr.json: remove its .zarray, and its .zattrs
// when its attributes were moved into the zarr.json
static void removeV2Metadata(const std::string &folder, const bool movedZattrs){
    remove((folder+"/.zarray").c_str());
    if(movedZattrs) remove((folder+"/.zattrs").c_str());
}

// Write the current Metadata to the .zarray and create subfolders if needed
void zarr::write_zarray(){
    // Zarr v3: zarr.json, replacing a v2 array's .zarray at the same path (no
    // dimension folders are made ahead: the writer makes them as needed)
    if(zarr_format == 3){
        normalizeDims();
        // (v3Metadata takes a .zattrs's attributes when there is no zarr.json)
        const bool hadZarray = fileExists(fileName+"/.zarray");
        const bool movedZattrs = !fileExists(fileName+"/zarr.json") && isJsonObjectFile(fileName+"/.zattrs");
        const json meta = v3Metadata(false);
        writeJsonFile(fileName, "zarr.json", meta);
        if(hadZarray) removeV2Metadata(fileName, movedZattrs);
        zarray = meta;
        return;
    }
    // A Zarr v3 array (zarr.json) is never turned into a v2 one in place
    if(!fileExists(fileName+"/.zarray") && fileExists(fileName+"/zarr.json")){
        throw std::string("zarrV3NotWritable:"+fileName);
    }
    normalizeDims();
    createSubfolders();
    set_jsonValues();
    write_jsonValues();
}

void zarr::set_zarr_format(const uint64_t zarr_format){
    if(zarr_format != 2 && zarr_format != 3) throw std::string("zarrFormatUnsupported:"+std::to_string(zarr_format));
    this->zarr_format = zarr_format;
    if(zarr_format == 3){
        order = "C";
        dimension_separator = "/";
    }
}

uint64_t zarr::get_shuffle() const{
    return shuffle;
}

void zarr::convert_to_v3(){
    if(zarr_format != 2 || !fileExists(fileName+"/.zarray")){
        throw std::string("v3Unsupported:"+fileName+" is not a Zarr v2 array");
    }
    const bool movedZattrs = isJsonObjectFile(fileName+"/.zattrs");
    const json meta = v3Metadata(true);
    writeJsonFile(fileName, "zarr.json", meta);
    removeV2Metadata(fileName, movedZattrs);
    zarray = meta;
    zarr_format = 3;
}

// Zarr v3 metadata (zarr.json) for this array: a new one (chunk files named
// "c/0/0/0"), or one describing an existing v2 array's chunk files (named as v2
// names them) to convert it without rewriting data. Throws "v3Unsupported:why"
// for what Zarr v3 cannot describe.
json zarr::v3Metadata(const bool existingChunks) const{
    auto unsupported = [](const std::string &why){ return std::string("v3Unsupported:")+why; };
    const uint64_t n = shape.size();
    const DtypeInfo info = dtypeInfo(dtype);
    const uint64_t bytes = info.bytes;
    const char kind = info.kind;
    std::string dataType;
    if(kind == 'b') dataType = "bool";
    else if(kind == 'u' || kind == 'i') dataType = std::string(kind == 'u' ? "uint" : "int")+std::to_string(8*bytes);
    else if(kind == 'f') dataType = info.bfloat16 ? "bfloat16" : "float"+std::to_string(8*bytes);
    else if(kind == 'c') dataType = "complex"+std::to_string(8*bytes);
    else throw unsupported("its data type \""+dtype+"\" is not supported");
    if(std::any_of(subfolders.begin(), subfolders.end(), [](uint64_t v){ return v != 0; })){
        throw unsupported("subfolders cannot be described in Zarr v3");
    }

    // Rewriting a v3 array's metadata keeps what the fields do not hold: its chunk
    // key encoding, chunk checksums, and an axis order other than C or F
    const bool rewrite = !existingChunks && zarr_format == 3 && isV3Json(zarray);

    // Each chunk's codecs: a transpose for F order (or another axis order), the
    // bytes, then the compressor
    json chunkCodecs = json::array();
    std::vector<uint64_t> axes;
    if(order != "C" && n > 1){
        axes.resize(n);
        for(uint64_t d = 0; d < n; d++) axes[d] = n-1-d;
    }
    else if(rewrite && n > 1){
        // (C and F order are what the order field says)
        const std::vector<uint64_t> old = get_chunkAxisOrder();
        bool cOrder = true, fOrder = true;
        for(uint64_t d = 0; d < n; d++){
            cOrder = cOrder && old[d] == d;
            fOrder = fOrder && old[d] == n-1-d;
        }
        if(!cOrder && !fOrder) axes = old;
    }
    if(!axes.empty()) chunkCodecs.push_back({{"name", "transpose"}, {"configuration", {{"order", axes}}}});
    if(bytes > 1) chunkCodecs.push_back({{"name", "bytes"}, {"configuration", {{"endian", info.bigEndian ? "big" : "little"}}}});
    else chunkCodecs.push_back({{"name", "bytes"}});
    const std::string compressor = get_compressor();
    if(compressor == "blosc"){
        if(cname != "lz4" && cname != "lz4hc" && cname != "blosclz" && cname != "zstd" && cname != "zlib"){
            if(!existingChunks) throw std::string("unsupportedCompressor");
            throw unsupported("the blosc compressor \""+cname+"\" is not supported");
        }
        // (numcodecs' automatic shuffle, -1, is bit shuffle for 1-byte data)
        const std::string shuffleName = shuffle == 0 ? "noshuffle" : shuffle == 1 ? "shuffle" :
                                        shuffle == 2 || bytes == 1 ? "bitshuffle" : "shuffle";
        chunkCodecs.push_back({{"name", "blosc"}, {"configuration", {{"cname", cname}, {"clevel", clevel},
            {"shuffle", shuffleName}, {"typesize", bytes}, {"blocksize", blocksize}}}});
    }
    else if(compressor == "gzip" && cname == "gzip") chunkCodecs.push_back({{"name", "gzip"}, {"configuration", {{"level", clevel}}}});
    else if(compressor == "zstd") chunkCodecs.push_back({{"name", "zstd"}, {"configuration", {{"level", (int64_t)clevel}, {"checksum", false}}}});
    else if(compressor != "none"){
        if(!existingChunks) throw std::string("unsupportedCompressor");
        throw unsupported("the \""+cname+"\" compressor has no Zarr v3 codec");
    }
    if(rewrite && get_chunkChecksum()) chunkCodecs.push_back({{"name", "crc32c"}});

    // Sharding: the inner chunks inside a shard, with its index
    json codecs = chunkCodecs;
    if(shard){
        if(n == 0) throw unsupported("sharding is not supported for 0-dimensional arrays");
        for(uint64_t d = 0; d < n; d++){
            if(!chunk_shape[d] || chunks[d]%chunk_shape[d]) throw unsupported("its shard shape is not a multiple of its inner chunk shape");
        }
        json indexCodecs = json::array({{{"name", "bytes"}, {"configuration", {{"endian", "little"}}}}});
        if(get_shardIndexChecksum()) indexCodecs.push_back({{"name", "crc32c"}});
        codecs = json::array({{{"name", "sharding_indexed"}, {"configuration", {
            {"chunk_shape", chunk_shape}, {"codecs", chunkCodecs}, {"index_codecs", indexCodecs},
            {"index_location", get_shardIndexAtStart() ? "start" : "end"}}}}});
    }

    // The fill value, as Zarr v3 writes it for the data type
    json fill;
    if(kind == 'b') fill = fill_value == "1" || fill_value == "true";
    else if(kind == 'f'){
        const bool bits = fill_value.size() > 2 && fill_value[0] == '0' && (fill_value[1] == 'x' || fill_value[1] == 'X');
        if(fill_value == "NaN" || fill_value == "Infinity" || fill_value == "-Infinity" || bits) fill = fill_value;
        else fill = strtod(fill_value.c_str(), NULL);
    }
    else if(kind == 'c') fill = complexFillJson(fill_value, bytes, true);
    else if(kind == 'u') fill = (uint64_t)strtoull(fill_value.c_str(), NULL, 10);
    else fill = (int64_t)strtoll(fill_value.c_str(), NULL, 10);

    // Attributes are kept, as a v2 array's .zattrs is when its .zarray is
    // rewritten: those of the zarr.json being replaced (and its dimension names),
    // or the .zattrs of a v2 array becoming a v3 one
    json attributes = json::object(), dimensionNames;
    auto readJson = [](const std::string &path){
        try{
            std::ifstream f(path);
            return json::parse(f);
        }
        catch(...){
            return json();
        }
    };
    if(!existingChunks && fileExists(fileName+"/zarr.json")){
        const json old = readJson(fileName+"/zarr.json");
        if(old.is_object()){
            const auto a = old.find("attributes");
            if(a != old.end() && a->is_object()) attributes = *a;
            const auto names = old.find("dimension_names");
            if(names != old.end() && names->is_array() && names->size() == n) dimensionNames = *names;
        }
    }
    else if(fileExists(fileName+"/.zattrs")){
        const json a = readJson(fileName+"/.zattrs");
        if(a.is_object()) attributes = a;
    }
    std::string keyEncoding = existingChunks ? "v2" : "default";
    if(rewrite) keyEncoding = extName(zarray.at("chunk_key_encoding"));

    json meta;
    meta["zarr_format"] = 3;
    meta["node_type"] = "array";
    meta["shape"] = shape;
    meta["data_type"] = dataType;
    meta["chunk_grid"] = {{"name", "regular"}, {"configuration", {{"chunk_shape", chunks}}}};
    meta["chunk_key_encoding"] = {{"name", keyEncoding}, {"configuration", {{"separator", dimension_separator}}}};
    meta["fill_value"] = fill;
    meta["codecs"] = codecs;
    meta["attributes"] = attributes;
    if(!dimensionNames.is_null()) meta["dimension_names"] = dimensionNames;
    return meta;
}

const std::string &zarr::get_fileName() const{
    return fileName;
}

void zarr::set_fileName(const std::string &fileName){
    // Handle the tilde character in filenames on Linux/Mac
    this->fileName = fileName;
    #ifndef _WIN32
    this->fileName = expandTilde(this->fileName.c_str());
    #endif
}

const uint64_t &zarr::get_chunks(const uint64_t &index) const{
    return chunks[index];
}

void zarr::set_chunks(const std::vector<uint64_t> &chunks){
    this->chunks = chunks;
}

const uint64_t &zarr::get_clevel() const{
    return clevel;
}

void zarr::set_clevel(const uint64_t &clevel){
    this->clevel = clevel;
}

const std::string &zarr::get_cname() const{
    return cname;
}

void zarr::set_cname(const std::string &cname){
    this->cname = cname;
}

const std::string &zarr::get_dimension_separator() const{
    return dimension_separator;
}

void zarr::set_dimension_separator(const std::string &dimension_separator){
    this->dimension_separator = dimension_separator;
}

const std::string &zarr::get_dtype() const{
    return dtype;
}

void zarr::set_dtype(const std::string &dtype){
    // 1-byte types have no byte order: "|u1", "|i1" and "|b1", as NumPy names them
    // (zarr-python 3 does not open v2 arrays whose dtype is "<u1")
    const DtypeInfo t = dtypeInfo(dtype);
    if(t.kind && t.bytes == 1 && dtype[0] != '|') this->dtype = std::string("|")+dtype.substr(1);
    else this->dtype = dtype;
}

const std::string &zarr::get_fill_value() const{
    return fill_value;
}

void zarr::set_fill_value(const std::string &fill_value){
    this->fill_value = fill_value;
}

void zarr::set_fill_value(const int64_t &fill_value){
    this->fill_value = std::to_string(fill_value);
}

const std::string &zarr::get_order() const{
    return order;
}

void zarr::set_order(const std::string &order){
    this->order = order;
}

const uint64_t &zarr::get_shape(const uint64_t &index) const{
    return shape[index];
}

void zarr::set_shape(const std::vector<uint64_t> &shape){
    this->shape = shape;
}

const uint64_t zarr::dtypeBytes() const{
    return dtypeInfo(dtype).bytes;
}

// Set the values of the JSON file to the current member values
void zarr::set_jsonValues(){
    zarray.clear();
    zarray["chunks"] = chunks;

    if(cname == "lz4" || cname == "blosclz" || cname == "lz4hc" || cname == "zlib" || cname == "zstd"){
        zarray["compressor"]["blocksize"] = blocksize;
        zarray["compressor"]["clevel"] = clevel;
        zarray["compressor"]["cname"] = cname;
        zarray["compressor"]["id"] = id;
        zarray["compressor"]["shuffle"] = shuffle;
    }
    else if(cname == "gzip"){
        zarray["compressor"]["id"] = cname;
        zarray["compressor"]["level"] = clevel;
    }
    // Uncompressed chunks
    else if(cname == "none") zarray["compressor"] = nullptr;
    else throw std::string("unsupportedCompressor"); 
    
    // dimension_separator only if dimension_separator is "/"
    if(dimension_separator == "/") zarray["dimension_separator"] = dimension_separator;

    zarray["dtype"] = dtype;
    // The fill value: floats as numbers or "NaN", "Infinity" and "-Infinity",
    // complex numbers as [real, imaginary], booleans as true or false, integers as
    // integers (and the Infinity of older metadata as the int64 limits)
    const DtypeInfo info = dtypeInfo(dtype);
    if(fill_value == "null") zarray["fill_value"] = nullptr;
    else if(info.kind == 'f') zarray["fill_value"] = floatFillJson(fill_value, info.bytes, info.bfloat16, false);
    else if(info.kind == 'c') zarray["fill_value"] = complexFillJson(fill_value, info.bytes, false);
    else if(info.kind == 'b') zarray["fill_value"] = fill_value == "1" || fill_value == "true";
    else if(fill_value == "NaN") zarray["fill_value"] = fill_value;
    else if(fill_value == "Infinity") zarray["fill_value"] = std::numeric_limits<int64_t>::max();
    else if(fill_value == "-Infinity") zarray["fill_value"] = std::numeric_limits<int64_t>::min();
    else if(info.kind == 'u' && !fill_value.empty() && fill_value[0] != '-') zarray["fill_value"] = (uint64_t)std::stoull(fill_value);
    else zarray["fill_value"] = std::stoll(fill_value);
    zarray["filters"] = nullptr;
    zarray["order"] = order;    
    zarray["shape"] = shape;

    // zarr_format just 2 for now
    zarray["zarr_format"] = 2;
    
    // Only add the subfolder parameter if subfolders is not all zeros
    if(!std::all_of(subfolders.begin(),
               subfolders.end(),
               [](int i){return !i;})){
        zarray["subfolders"] = subfolders;
    }
    

    // Sharding
    if(shard){
        // Create the inner JSON objects
        json innerConfig;
        json innerCodec;
        json innerCodecConfig;

        // Set the inner JSON values
        // blosc
        if(cname == "lz4" || cname == "blosclz" || cname == "lz4hc" || cname == "zlib" || cname == "zstd"){
            innerCodec["name"] = id;
            innerCodecConfig["cname"] = cname;
            innerCodecConfig["clevel"] = clevel;
            innerCodecConfig["shuffle"] = "shuffle";
            innerCodecConfig["typesize"] = dtypeBytes();
            innerCodecConfig["blocksize"] = blocksize;
        }
        // gzip
        else if(cname == "gzip"){
            innerCodec["name"] = cname;
            innerCodecConfig["level"] = clevel;
        }
        
        innerCodec["configuration"] = innerCodecConfig;

        innerConfig["chunk_shape"] = chunk_shape;
        innerConfig["codecs"] = {innerCodec};
        
        // codecs is an array of json objects
        zarray["codecs"] = {{{"name", "sharding_indexed"}, 
                            {"configuration", innerConfig}}};
    }
}

// Write the current JSON to disk
void zarr::write_jsonValues(){
    // If the .zarray file does not exist then build the zarr fileName path recursively
    const std::string fileNameFinal(fileName+"/.zarray");

    if(!fileExists(fileNameFinal)){
        mkdirRecursive(fileName.c_str());
    }

    const std::string uuid = generateUUID();
    const std::string fnFull(fileName+"/.zarray"+uuid);

    std::ofstream o(fnFull);
    if(!o.good()) throw std::string("cannotOpenZarray:"+fnFull);
    o << std::setw(4) << zarray << std::endl;
    o.close();

    if(!renameReplace(fnFull, fileNameFinal)){
        remove(fnFull.c_str());
        throw std::string("cannotOpenZarray:"+fileNameFinal);
    }
}

const std::string zarr::get_subfoldersString(const std::vector<uint64_t> &cAV) const{
    if(std::all_of(subfolders.begin(), subfolders.end(), [](uint64_t i){return !i;})) return "";

    std::string name;
    for(uint64_t d = 0; d < cAV.size(); d++){
        const uint64_t v = (d < subfolders.size() && subfolders[d] > 0) ? cAV[d]/subfolders[d] : 0;
        if(d) name += "_";
        name += std::to_string(v);
    }
    return name;
}

void zarr::set_subfolders(const std::vector<uint64_t> &subfolders){
    this->subfolders = subfolders;
    zarray["subfolders"] = subfolders;
}

void zarr::set_shardData(){
    normalizeDims();
    const uint64_t n = shape.size();
    shards = std::vector<uint64_t>(n);
    chunksPerShard = std::vector<uint64_t>(n);
    numShards = 1;
    numChunksPerShard = 1;
    for(uint64_t i = 0; i < n; i++){
        chunksPerShard[i] = fastCeilDiv(chunks[i],chunk_shape[i]);
        // A shard holds whole inner chunks, so it spans its shape rounded up to a
        // multiple of the inner chunk shape
        const uint64_t extent = chunksPerShard[i]*chunk_shape[i];
        shards[i] = shape[i]/extent + (shape[i]%extent ? 1 : 0);
        numShards *= shards[i];
        numChunksPerShard *= chunksPerShard[i];
    }
}

const bool &zarr::get_shard() const{
    return shard;
}

void zarr::set_shard(const bool shard){
    this->shard = shard;
}

const uint64_t &zarr::get_chunk_shape(const uint64_t &index) const{
    return chunk_shape[index];
}

void zarr::set_chunk_shape(const std::vector<uint64_t> &chunk_shape){
    this->chunk_shape = chunk_shape;
    // The shard layout for the current shape, keeping the per-axis values as given:
    // they are fit to the number of dimensions once the shape is known (set_shape
    // may come after this; fit to the default 3D shape, the axes after the third
    // were lost)
    const std::vector<uint64_t> givenChunks(chunks), givenSubfolders(subfolders), givenChunkShape(chunk_shape);
    set_shardData();
    chunks = givenChunks;
    subfolders = givenSubfolders;
    this->chunk_shape = givenChunkShape;
}

const uint64_t &zarr::get_numShards() const{
    return numShards;
}

const uint64_t &zarr::get_numChunksPerShard() const{
    return numChunksPerShard;
}

// Fast ceiling for the subfolder function
uint64_t zarr::fastCeilDiv(uint64_t num, uint64_t denom){
    return 1 + ((num - 1) / denom);
}

// Name of the chunk (or subfolder) at flat index i of a grid of `count` indices
// starting at `first`, last axis fastest, joined by `sep`
static std::string gridIndexName(uint64_t i, const std::vector<uint64_t> &first,
                                 const std::vector<uint64_t> &count, const std::string &sep){
    const uint64_t n = count.size();
    // Stack storage for the usual numbers of dimensions, heap beyond that
    uint64_t local[64];
    std::vector<uint64_t> heap;
    uint64_t* idx = local;
    if(n > 64){ heap.resize(n); idx = heap.data(); }
    for(int64_t d = (int64_t)n-1; d >= 0; d--){
        idx[d] = first[d] + i % count[d];
        i /= count[d];
    }
    std::string name;
    for(uint64_t d = 0; d < n; d++){
        if(d) name += sep;
        name += std::to_string(idx[d]);
    }
    return name;
}

// Create subfolder "chunks"
void zarr::createSubfolders(){
    normalizeDims();
    const uint64_t n = shape.size();

    // (The folders of '/' chunk keys are made when chunks are first written in
    // them, as for Zarr v3, not all when the array is created: making every folder
    // of a big array ahead took seconds on NFS, even when only a few chunks of it
    // were written.)

    // If all elements are zero then we don't make subfolders
    if(std::all_of(subfolders.begin(),
                   subfolders.end(),
                   [](uint64_t i){return !i;}))
    {
        return;
    }

    // Use shards instead of chunks when sharding
    if(shard) set_shardData();
    std::vector<uint64_t> nSubfolders(n,1);
    uint64_t total = 1;
    for(uint64_t d = 0; d < n; d++){
        const uint64_t nChunks = shard ? shards[d] : fastCeilDiv(shape[d],chunks[d]);
        if(subfolders[d] > 0) nSubfolders[d] = fastCeilDiv(nChunks,subfolders[d]);
        total *= nSubfolders[d];
    }

    // Create subfolders
    const std::vector<uint64_t> zeros(n,0);
    #pragma omp parallel for
    for(uint64_t i = 0; i < total; i++){
        const std::string currName(fileName+"/"+gridIndexName(i,zeros,nSubfolders,"_"));
        mkdirRecursive(currName.c_str());
    }
}

const std::string zarr::chunkNameToShardName(const std::string &chunkName) const{
    const std::vector<uint64_t> cAV = get_chunkAxisVals(chunkName);
    if(cAV.empty()) return "0";
    std::string name;
    for(uint64_t d = 0; d < cAV.size() && d < chunksPerShard.size(); d++){
        if(d) name += dimension_separator;
        name += std::to_string(cAV[d]/chunksPerShard[d]);
    }
    return name;
}

// Indices of the shard that holds the inner chunk with indices cAV
const std::vector<uint64_t> zarr::chunkToShard(const std::vector<uint64_t> &cAV) const{
    std::vector<uint64_t> s(cAV.size());
    for(uint64_t d = 0; d < cAV.size() && d < chunksPerShard.size(); d++) s[d] = cAV[d]/chunksPerShard[d];
    return s;
}

// The sharding codec's configuration (Zarr v3 sharding_indexed), if any
static const json* shardingConfig(const json &zarray){
    const auto codecs = zarray.find("codecs");
    if(codecs == zarray.end() || !codecs->is_array() || codecs->empty()) return nullptr;
    const json &c = codecs->at(0);
    if(!c.is_object() || c.value("name", "") != "sharding_indexed") return nullptr;
    const auto config = c.find("configuration");
    return config != c.end() && config->is_object() ? &*config : nullptr;
}

bool zarr::get_shardIndexAtStart() const{
    const json* config = shardingConfig(zarray);
    return config && config->value("index_location", "end") == "start";
}

bool zarr::get_shardIndexChecksum() const{
    const json* config = shardingConfig(zarray);
    if(!config) return true;
    const auto indexCodecs = config->find("index_codecs");
    if(indexCodecs == config->end() || !indexCodecs->is_array()) return true;
    for(const json &c : *indexCodecs){
        if(c.is_object() && c.value("name", "") == "crc32c") return true;
    }
    return false;
}

uint64_t zarr::get_zarr_format() const{
    return zarr_format;
}

const std::string zarr::chunkKey(const std::string &chunkName) const{
    if(zarr_format != 3 || !isV3Json(zarray)) return chunkName;
    const auto keys = zarray.find("chunk_key_encoding");
    if(extName(*keys) != "default") return chunkName;
    // A 0-dimensional array's single chunk is "c"
    return shape.empty() ? std::string("c") : "c"+dimension_separator+chunkName;
}

const std::vector<uint64_t> zarr::get_chunkAxisOrder() const{
    const uint64_t n = shape.size();
    std::vector<uint64_t> axes(n);
    for(uint64_t d = 0; d < n; d++) axes[d] = order == "C" ? d : n-1-d;
    if(zarr_format == 3 && isV3Json(zarray)){
        const V3Codecs c = parseV3Codecs(chunkCodecsV3(zarray), n, false);
        if(!c.axisOrder.empty()) axes = c.axisOrder;
    }
    return axes;
}

const std::string zarr::get_compressor() const{
    if(cname == "none") return "none";
    // (zlib streams decode the same way as gzip)
    if(cname == "gzip" || id == "gzip" || id == "zlib") return "gzip";
    if(id == "zstd") return "zstd";
    if(id.empty() || id == "blosc") return "blosc";
    return id;
}

bool zarr::get_chunkChecksum() const{
    if(zarr_format != 3 || !isV3Json(zarray)) return false;
    return parseV3Codecs(chunkCodecsV3(zarray), shape.size(), false).checksum;
}

const uint64_t zarr::get_ShardPosition(const std::vector<uint64_t> &cAV) const{
    uint64_t pos = 0;
    for(uint64_t d = 0; d < cAV.size(); d++) pos = pos*shards[d] + cAV[d];
    return pos;
}

// Position of a chunk inside its shard, last axis fastest
const uint64_t zarr::get_chunkShardPosition(const std::vector<uint64_t> &cAV) const{
    uint64_t pos = 0;
    for(uint64_t d = 0; d < cAV.size(); d++) pos = pos*chunksPerShard[d] + (cAV[d]%chunksPerShard[d]);
    return pos;
}

// Parse a chunk name ("x.y.z", "a/b/c/d", ...) into its per-axis indices
const std::vector<uint64_t> zarr::get_chunkAxisVals(const std::string &fileName) const{
    std::vector<uint64_t> cAV;
    // A 0-dimensional array's single chunk ("0") has no indices
    if(shape.empty()) return cAV;
    cAV.reserve(shape.size());
    const char* p = fileName.c_str();
    while(true){
        char* end;
        cAV.push_back(strtoull(p, &end, 10));
        if(*end == '\0' || end == p) break;
        p = end+1;
    }
    return cAV;
}

void zarr::set_chunkInfo(const std::vector<uint64_t> &startCoords,
                         const std::vector<uint64_t> &endCoords)
{
    normalizeDims();
    const uint64_t n = shape.size();
    // A 0-dimensional array (scalar) is a single chunk named "0"
    if(n == 0){
        numChunks = 1;
        chunkNames = std::vector<std::string>(1, "0");
        return;
    }

    // Chunks (or shards when sharded) that the region touches along each axis. A
    // shard spans its shape rounded up to a multiple of the inner chunk shape
    if(shard) set_shardData();
    std::vector<uint64_t> first(n), count(n);
    uint64_t numOuter = 1;
    for(uint64_t d = 0; d < n; d++){
        const uint64_t extent = shard ? chunksPerShard[d]*chunk_shape[d] : chunks[d];
        first[d] = startCoords[d]/extent;
        const uint64_t last = endCoords[d]/extent + (endCoords[d]%extent ? 1 : 0);
        count[d] = last > first[d] ? last-first[d] : 0;
        numOuter *= count[d];
    }

    // Default behavior for when chunks are not sharded: every chunk, last axis fastest
    if(!shard){
        numChunks = numOuter;
        chunkNames = std::vector<std::string>(numChunks);
        #pragma omp parallel for
        for(uint64_t i = 0; i < numChunks; i++){
            chunkNames[i] = gridIndexName(i,first,count,dimension_separator);
        }
    }
    // Sharding: every inner chunk of every shard, shard by shard (both last axis fastest)
    else{
        numChunks = numOuter*numChunksPerShard;
        chunkNames = std::vector<std::string>(numChunks);
        #pragma omp parallel for
        for(uint64_t i = 0; i < numChunks; i++){
            uint64_t o = i/numChunksPerShard, c = i%numChunksPerShard;
            uint64_t local[64];
            std::vector<uint64_t> heap;
            uint64_t* idx = local;
            if(n > 64){ heap.resize(n); idx = heap.data(); }
            for(int64_t d = (int64_t)n-1; d >= 0; d--){
                const uint64_t s = first[d] + o%count[d];
                o /= count[d];
                idx[d] = s*chunksPerShard[d] + c%chunksPerShard[d];
                c /= chunksPerShard[d];
            }
            std::string name;
            for(uint64_t d = 0; d < n; d++){
                if(d) name += dimension_separator;
                name += std::to_string(idx[d]);
            }
            chunkNames[i] = name;
        }
    }
}

const std::string &zarr::get_chunkNames(const uint64_t &index) const{
    return chunkNames[index];
}

const uint64_t &zarr::get_numChunks() const{
    return numChunks;
}

const std::string &zarr::get_errString() const{
    return errString;
}

void zarr::set_errString(const std::string &errString){
    this->errString = errString;
}

uint64_t zarr::get_ndims() const{
    return shape.size();
}

void zarr::normalizeDims(){
    const uint64_t n = shape.size();
    if(chunks.size() > n) chunks.resize(n);
    while(chunks.size() < n) chunks.push_back(chunks.size() < 3 ? 256 : 1);
    if(subfolders.size() > n) subfolders.resize(n);
    while(subfolders.size() < n) subfolders.push_back(0);
    if(chunk_shape.size() > n) chunk_shape.resize(n);
    while(chunk_shape.size() < n) chunk_shape.push_back(1);
}
