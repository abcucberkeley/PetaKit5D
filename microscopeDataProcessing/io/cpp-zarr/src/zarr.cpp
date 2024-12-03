#include <cstdint>
#include <omp.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <stdarg.h>
#include <sys/time.h>
#else
#include <uuid/uuid.h>
#endif
#include <fstream>
#include "zarr.h"
#include "helperfunctions.h"

// Create a blank zarr object with default values
zarr::zarr() :
fileName(""), chunks({256,256,256}), blocksize(0),
clevel(5), cname("lz4"), id("blosc"), shuffle(1), dtype("<u2"),
dimension_separator("."), fill_value("0"), filters({}), order("F"), 
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
    if(!fileExists(fileName+"/.zarray")){
        throw std::string("metadataFileMissing:"+fileName);
        //mexErrMsgIdAndTxt("zarr:zarrayError","Metadata file in \"%s\" is missing. Does the file exist?",fileName.c_str());
    }
    std::ifstream f(fileName+"/.zarray");
    zarray = json::parse(f);

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
        try{
            // Try blosc compression types
            cname = zarray.at("compressor").at("cname");
            clevel = zarray.at("compressor").at("clevel");
            blocksize = zarray.at("compressor").at("blocksize");
            id = zarray.at("compressor").at("id");
            shuffle = zarray.at("compressor").at("shuffle");
        }
        catch(...){
            // Try gzip
            clevel = zarray.at("compressor").at("level");
            cname = zarray.at("compressor").at("id");
            blocksize = 0;
            id = "";
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
        if(zarray.at("fill_value").empty()) fill_value = "0";
        else{
            if(zarray.at("fill_value").type() == json::value_t::number_integer || 
               zarray.at("fill_value").type() == json::value_t::number_unsigned || 
               zarray.at("fill_value").type() == json::value_t::number_float)
            {
                fill_value = std::to_string((int64_t)zarray.at("fill_value"));
            }
            else fill_value = zarray.at("fill_value");
            // TODO: Make NaN actually NaN here and in other functions
            if(fill_value == "null" || fill_value == "NaN") fill_value = "0";
            else if(fill_value == "Infinity") fill_value = std::to_string(std::numeric_limits<int64_t>::max());
            else if(fill_value == "-Infinity") fill_value = std::to_string(std::numeric_limits<int64_t>::min());
        }
        //filters = "";

        order = zarray.at("order");
        shape = zarray.at("shape").get<std::vector<uint64_t>>();
        zarr_format = zarray.at("zarr_format");
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

// Write the current Metadata to the .zarray and create subfolders if needed
void zarr::write_zarray(){
    createSubfolders();
    set_jsonValues();
    write_jsonValues();
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
    this->dtype = dtype;
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
    if(dtype.size() != 3) return 0;
    else if(dtype[2] == '1') return 1;
    else if(dtype[2] == '2') return 2;
    else if(dtype[2] == '4') return 4;
    else if(dtype[2] == '8') return 8;
    else return 0;
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
    else throw std::string("unsupportedCompressor"); 
    
    // dimension_separator only if dimension_separator is "/"
    if(dimension_separator == "/") zarray["dimension_separator"] = dimension_separator;

    zarray["dtype"] = dtype;
    if(fill_value == "NaN") zarray["fill_value"] = fill_value;
    else if(fill_value == "null") zarray["fill_value"] = nullptr;
    else if(fill_value == "Infinity") zarray["fill_value"] = std::numeric_limits<int64_t>::max();
    else if(fill_value == "-Infinity") zarray["fill_value"] = std::numeric_limits<int64_t>::min();
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

    rename(fnFull.c_str(),fileNameFinal.c_str());
}

const std::string zarr::get_subfoldersString(const std::vector<uint64_t> &cAV) const{
    if(subfolders[0] == 0 && subfolders[1] == 0 && subfolders[2] == 0) return "";

    std::vector<uint64_t> currVals = {0,0,0};
    if(subfolders[0] > 0) currVals[0] = cAV[0]/subfolders[0];
    if(subfolders[1] > 0) currVals[1] = cAV[1]/subfolders[1];
    if(subfolders[2] > 0) currVals[2] = cAV[2]/subfolders[2];

    return std::string(std::to_string(currVals[0])+"_"+
                       std::to_string(currVals[1])+"_"+
                       std::to_string(currVals[2]));
}

void zarr::set_subfolders(const std::vector<uint64_t> &subfolders){
    this->subfolders = subfolders;
    zarray["subfolders"] = subfolders;
}

void zarr::set_shardData(){
    uint64_t prod = 1;
    shards = std::vector<uint64_t>(3);
    chunksPerShard = std::vector<uint64_t>(3);
    for(uint64_t i = 0; i < 3; i++){
        chunksPerShard[i] = fastCeilDiv(chunks[i],chunk_shape[i]);
        shards[i] = ceil((double)shape[i]/(double)chunks[i]);
    }
    prod = 1;
    for (const auto& i : shards){
        prod *= i;
    }
    numShards = prod;
    numChunksPerShard = fastCeilDiv(chunks[0],chunk_shape[0])*fastCeilDiv(chunks[1],chunk_shape[1])*fastCeilDiv(chunks[2],chunk_shape[2]);
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
    set_shardData();
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

// Create subfolder "chunks"
void zarr::createSubfolders(){
    // dimension_separator subfolders
    if(dimension_separator == "/"){
        set_chunkInfo({0,0,0},shape);
        #pragma omp parallel for
        for(uint64_t i = 0; i < chunkNames.size(); i++){
            makeDimensionFolders(fileName+"/"+chunkNames[i]);
        }
    }

    // If all elements are zero then we don't make subfolders
    if(std::all_of(subfolders.begin(),
                   subfolders.end(),
                   [](int i){return !i;}))
    {
        return;
    }

    std::vector<uint64_t> nChunks;
    if(!shard){
        nChunks = {fastCeilDiv(shape[0],chunks[0]),
                   fastCeilDiv(shape[1],chunks[1]),
                   fastCeilDiv(shape[2],chunks[2])};
    }
    else{
        set_shardData();
        // Use shards instead
        nChunks = {shards[0],
                   shards[1],
                   shards[2]};
    }
    std::vector<uint64_t> nSubfolders = {1,1,1};
    for(uint64_t i = 0; i < nSubfolders.size(); i++){
        if(subfolders[i] > 0){
            nSubfolders[i] = fastCeilDiv(nChunks[i],subfolders[i]);
        }
    }

    // Create subfolders
    #pragma omp parallel for collapse(3)
    for(uint64_t x = 0; x < nSubfolders[0]; x++){
        for(uint64_t y = 0; y < nSubfolders[1]; y++){
            for(uint64_t z = 0; z < nSubfolders[2]; z++){
                std::string currName(fileName+"/"+std::to_string(x)+"_"+
                                     std::to_string(y)+"_"+std::to_string(z));
                mkdirRecursive(currName.c_str());
            }
        }
    }
}

const std::string zarr::chunkNameToShardName(const std::string &chunkName) const{
    std::vector<uint64_t> cAV = get_chunkAxisVals(chunkName);
    /*
    return std::to_string(cAV[0]/(uint64_t)ceil((double)chunks[0]/(double)chunk_shape[0]))+
                        dimension_separator+std::to_string(cAV[1]/(uint64_t)ceil((double)chunks[1]/(double)chunk_shape[1]))+
                        dimension_separator+std::to_string(cAV[2]/(uint64_t)ceil((double)chunks[2]/(double)chunk_shape[2]));
    */
    return std::to_string(cAV[0]/chunksPerShard[0])+
                          dimension_separator+std::to_string(cAV[1]/chunksPerShard[1])+
                          dimension_separator+std::to_string(cAV[2]/chunksPerShard[2]);
}

const std::vector<uint64_t> zarr::chunkToShard(const std::vector<uint64_t> &cAV) const{
    uint64_t x = ceil((double)cAV[0]/ceil((double)chunks[0]/(double)chunk_shape[0]));
    if(x) x--;
    uint64_t y = ceil((double)cAV[1]/ceil((double)chunks[1]/(double)chunk_shape[1]));
    if(y) y--;
    uint64_t z = ceil((double)cAV[2]/ceil((double)chunks[2]/(double)chunk_shape[2]));
    if(z) z--;

    return {x,y,z};
}

const uint64_t zarr::get_ShardPosition(const std::vector<uint64_t> &cAV) const{
    return (cAV[0]*(shards[1]*shards[2]))+(cAV[1]*shards[2])+cAV[2];
}

const uint64_t zarr::get_chunkShardPosition(const std::vector<uint64_t> &cAV) const{
    return (cAV[2]%chunksPerShard[2]) +
        ((cAV[1]%chunksPerShard[1])*chunksPerShard[2]) +
        ((cAV[0]%chunksPerShard[0])*chunksPerShard[2]*chunksPerShard[1]);
}

const std::vector<uint64_t> zarr::get_chunkAxisVals(const std::string &fileName) const{
    std::vector<uint64_t> cAV(3);
    char* ptr;
    cAV[0] = strtol(fileName.c_str(), &ptr, 10);
    ptr++;
    cAV[1] = strtol(ptr, &ptr, 10);
    ptr++;
    cAV[2] = strtol(ptr, &ptr, 10);
    return cAV;
}

void zarr::set_chunkInfo(const std::vector<uint64_t> &startCoords,
                         const std::vector<uint64_t> &endCoords)
{
    //std::cout << xChunks << yChunks << zChunks << " numChunks: " << numChunks << std::endl;
    
    // Defualt behavior for when chunks are not sharded
    if(!shard){
        uint64_t xStartAligned = startCoords[0]-(startCoords[0]%chunks[0]);
        uint64_t yStartAligned = startCoords[1]-(startCoords[1]%chunks[1]);
        uint64_t zStartAligned = startCoords[2]-(startCoords[2]%chunks[2]);
        uint64_t xStartChunk = (xStartAligned/chunks[0]);
        uint64_t yStartChunk = (yStartAligned/chunks[1]);
        uint64_t zStartChunk = (zStartAligned/chunks[2]);
    
        uint64_t xEndAligned = endCoords[0];
        uint64_t yEndAligned = endCoords[1];
        uint64_t zEndAligned = endCoords[2];
    
        if(xEndAligned%chunks[0]) xEndAligned = endCoords[0]-(endCoords[0]%chunks[0])+chunks[0];
        if(yEndAligned%chunks[1]) yEndAligned = endCoords[1]-(endCoords[1]%chunks[1])+chunks[1];
        if(zEndAligned%chunks[2]) zEndAligned = endCoords[2]-(endCoords[2]%chunks[2])+chunks[2];
        uint64_t xEndChunk = (xEndAligned/chunks[0]);
        uint64_t yEndChunk = (yEndAligned/chunks[1]);
        uint64_t zEndChunk = (zEndAligned/chunks[2]);
    
        uint64_t xChunks = (xEndChunk-xStartChunk);
        uint64_t yChunks = (yEndChunk-yStartChunk);
        uint64_t zChunks = (zEndChunk-zStartChunk);
        numChunks = xChunks*yChunks*zChunks;

        chunkNames = std::vector<std::string>(numChunks);
        #pragma omp parallel for collapse(3)
        for(uint64_t x = xStartChunk; x < xEndChunk; x++){
            for(uint64_t y = yStartChunk; y < yEndChunk; y++){
                for(uint64_t z = zStartChunk; z < zEndChunk; z++){
                    uint64_t currFile = (z-zStartChunk)+((y-yStartChunk)*zChunks)+((x-xStartChunk)*yChunks*zChunks);
                    chunkNames[currFile] = std::to_string(x)+dimension_separator+std::to_string(y)+dimension_separator+std::to_string(z);
                }
            }
        }
    }
    // Sharding
    else{
        

        uint64_t xStartAligned = startCoords[0]-(startCoords[0]%chunks[0]);
        uint64_t yStartAligned = startCoords[1]-(startCoords[1]%chunks[1]);
        uint64_t zStartAligned = startCoords[2]-(startCoords[2]%chunks[2]);
        uint64_t xStartChunk = (xStartAligned/chunks[0]);
        uint64_t yStartChunk = (yStartAligned/chunks[1]);
        uint64_t zStartChunk = (zStartAligned/chunks[2]);
    
        uint64_t xEndAligned = endCoords[0];
        uint64_t yEndAligned = endCoords[1];
        uint64_t zEndAligned = endCoords[2];
    
        if(xEndAligned%chunks[0]) xEndAligned = endCoords[0]-(endCoords[0]%chunks[0])+chunks[0];
        if(yEndAligned%chunks[1]) yEndAligned = endCoords[1]-(endCoords[1]%chunks[1])+chunks[1];
        if(zEndAligned%chunks[2]) zEndAligned = endCoords[2]-(endCoords[2]%chunks[2])+chunks[2];
        uint64_t xEndChunk = (xEndAligned/chunks[0]);
        uint64_t yEndChunk = (yEndAligned/chunks[1]);
        uint64_t zEndChunk = (zEndAligned/chunks[2]);
    
        uint64_t xChunks = (xEndChunk-xStartChunk);
        uint64_t yChunks = (yEndChunk-yStartChunk);
        uint64_t zChunks = (zEndChunk-zStartChunk);
        numChunks = xChunks*yChunks*zChunks;

        set_shardData();

        std::vector<uint64_t> startShard = {xStartChunk,yStartChunk,zStartChunk};
        std::vector<uint64_t> endShard = {xEndChunk,yEndChunk,zEndChunk};
        uint64_t numShardsI = (endShard[0]-startShard[0])*(endShard[1]-startShard[1])*(endShard[2]-startShard[2])*numChunksPerShard;
        numChunks = numShardsI;
        


        chunkNames = std::vector<std::string>(numShardsI);
        // Add parallel later
        uint64_t currFile = 0;
        
        for(uint64_t x = startShard[0]; x < endShard[0]; x++){
            for(uint64_t y = startShard[1]; y < endShard[1]; y++){
                for(uint64_t z = startShard[2]; z < endShard[2]; z++){
                    uint64_t xStart = x * fastCeilDiv(chunks[0],chunk_shape[0]);
                    uint64_t xEnd = (x + 1) * fastCeilDiv(chunks[0],chunk_shape[0]);
                
                    uint64_t yStart = y * fastCeilDiv(chunks[1],chunk_shape[1]);
                    uint64_t yEnd = (y + 1) * fastCeilDiv(chunks[1],chunk_shape[1]);
                
                    uint64_t zStart = z * fastCeilDiv(chunks[2],chunk_shape[2]);
                    uint64_t zEnd = (z + 1) * fastCeilDiv(chunks[2],chunk_shape[2]);
        
                    for (uint64_t xI = xStart; xI < xEnd; xI++) {
                        for (uint64_t yI = yStart; yI < yEnd; yI++) {
                            for (uint64_t zI = zStart; zI < zEnd; zI++) {
                                //uint64_t currFile = (z-zStartChunk)+((y-yStartChunk)*zChunks)+((x-xStartChunk)*yChunks*zChunks);
                                chunkNames[currFile] = std::to_string(xI)+dimension_separator+std::to_string(yI)+dimension_separator+std::to_string(zI);
                                //std::cout << currFile << " " << chunkNames[currFile] << " " << chunkNameToShardName(chunkNames[currFile]) << " " << get_chunkShardPosition(get_chunkAxisVals(chunkNames[currFile])) << std::endl;
                                currFile++;
                            }
                        }
                    }
                }
            }
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
