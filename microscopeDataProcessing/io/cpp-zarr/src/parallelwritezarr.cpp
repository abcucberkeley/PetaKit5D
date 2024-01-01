#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <omp.h>
#include <stddef.h>
#ifdef _WIN32
#include <sys/time.h>
#else
#include <uuid/uuid.h>
#endif
#include <sys/stat.h>
#include <fstream>
#include <algorithm>
#include "blosc.h"
#include "parallelreadzarr.h"
#include "parallelwritezarr.h"
#include "helperfunctions.h"
#include "zarr.h"
#include "zlib.h"

uint32_t crc32c(const uint8_t* data, size_t length) {
    uint32_t crc = 0xFFFFFFFF;
    // CRC32C
    const uint32_t polynomial = 0x82F63B78;

    for (size_t i = 0; i < length; ++i) {
        crc ^= data[i];
        for (size_t j = 0; j < 8; ++j) {
            crc = (crc >> 1) ^ (-(crc & 1) & polynomial);
        }
    }

    return ~crc;
}

uint8_t parallelWriteZarr(zarr &Zarr, void* zarrArr,
                          const std::vector<uint64_t> &startCoords,
                          const std::vector<uint64_t> &endCoords,
                          const std::vector<uint64_t> &writeShape,
                          const uint64_t bits, const bool useUuid,
                          const bool crop, const bool sparse){
    const uint64_t bytes = (bits/8);

    int32_t numWorkers = omp_get_max_threads();

    int32_t nBloscThreads = 1;
    if(numWorkers>Zarr.get_numChunks()){
        nBloscThreads = std::ceil(((double)numWorkers)/((double)Zarr.get_numChunks()));
        numWorkers = Zarr.get_numChunks();
    }

    int32_t batchSize = (Zarr.get_numChunks()-1)/numWorkers+1;
    if(Zarr.get_shard()){
        // batchSize has to align to a shard
        if(batchSize <= Zarr.get_numChunksPerShard()) batchSize = Zarr.get_numChunksPerShard();
        else batchSize += (Zarr.get_numChunksPerShard()-(batchSize%Zarr.get_numChunksPerShard()));

        // The chunk size is actually the inner chunk size now
        Zarr.set_chunks({Zarr.get_chunk_shape(0),Zarr.get_chunk_shape(1),Zarr.get_chunk_shape(2)});
    }

    const uint64_t s = Zarr.get_chunks(0)*Zarr.get_chunks(1)*Zarr.get_chunks(2);
    const uint64_t sB = s*bytes;

    const std::string uuid(generateUUID());

    void* zeroChunkUnc = NULL;
    if(sparse){
        zeroChunkUnc = calloc(s,bytes);
    }

    int err = 0;
    std::string errString;
    #pragma omp parallel for
    for(int32_t w = 0; w < numWorkers; w++){
        void* chunkUnC = malloc(sB);
        void* chunkC = malloc(sB+BLOSC_MAX_OVERHEAD);
        void* cRegion = nullptr;
        int64_t currChunk = -1;
        uint64_t shardFooterSize;
        uint64_t* shardFooter = nullptr;
        std::vector<uint64_t> cAV;
        if(Zarr.get_shard()){
            shardFooterSize = Zarr.get_numChunksPerShard()*2;
            shardFooter = (uint64_t*)malloc(shardFooterSize*sizeof(uint64_t));
        }
        uint64_t lastF = 0;
        bool unWritten = true;
        for(int64_t f = w*batchSize; f < (w+1)*batchSize; f++){
            if(f>=Zarr.get_numChunks()  || err) break;
            lastF = f;

            if(Zarr.get_shard()){
                unWritten = true;
                currChunk++;
                std::vector<uint64_t> pAV = Zarr.get_chunkAxisVals(Zarr.get_chunkNames(f));
                bool pad = pAV[0] > endCoords[0]/Zarr.get_chunk_shape(0) ||
                    pAV[1] > endCoords[1]/Zarr.get_chunk_shape(1) ||
                    pAV[2] > endCoords[2]/Zarr.get_chunk_shape(2);
                
                if(currChunk == Zarr.get_numChunksPerShard() || pad){
                    if(pad){
                        shardFooter[currChunk*2] = std::numeric_limits<uint64_t>::max();
                        shardFooter[(currChunk*2)+1] = std::numeric_limits<uint64_t>::max();

                        // Edge case for when we can't write because padded chunks are in the middle
                        if(currChunk != Zarr.get_numChunksPerShard()){
                            continue;
                        }
                    }
                    unWritten = false;
                    
                    // calculate CRC32C
                    uint32_t shardFooterCRC32C = crc32c(reinterpret_cast<uint8_t*>(shardFooter), shardFooterSize*sizeof(uint64_t));
                    // Not sure how sharding interacts with the subfolders at the moment (cAV needs to be converted)
                    const std::string subfolderName = Zarr.get_subfoldersString(cAV);
                    std::string fileName(Zarr.get_fileName()+"/"+subfolderName+"/"+Zarr.chunkNameToShardName(Zarr.get_chunkNames(f-1)));
                    std::string fileNameFinal;
                    if(useUuid){
                        fileNameFinal = std::string(fileName);
                        fileName.append(uuid);
                    }
                    std::ofstream file(fileName, std::ios::binary | std::ios::app);

                    if(!file.is_open()){
                            #pragma omp critical
                            {
                                err = 1;
                                errString = "Check permissions or filepath. Cannot write to path: "+
                                    fileName+"\n";
                            }
                            break;
                    }
                    file.write(reinterpret_cast<char*>(shardFooter),shardFooterSize*sizeof(uint64_t));
                    file.write(reinterpret_cast<char*>(&shardFooterCRC32C),sizeof(uint32_t));
                    file.close();
                    if(useUuid){
                        rename(fileName.c_str(),fileNameFinal.c_str());
                    }
                    /*
                    if(pad){
                        f += (Zarr.get_numChunksPerShard()-currChunk-1);
                        currChunk = -1;
                        continue;
                    }
                    */
                    currChunk = 0;
                }
            }
            
            cAV = Zarr.get_chunkAxisVals(Zarr.get_chunkNames(f));
            cRegion = nullptr;

            if(crop && ((((cAV[0])*Zarr.get_chunks(0)) < startCoords[0] || ((cAV[0]+1)*Zarr.get_chunks(0) > endCoords[0] && endCoords[0] < Zarr.get_shape(0)))
                        || (((cAV[1])*Zarr.get_chunks(1)) < startCoords[1] || ((cAV[1]+1)*Zarr.get_chunks(1) > endCoords[1] && endCoords[1] < Zarr.get_shape(1)))
                        || (((cAV[2])*Zarr.get_chunks(2)) < startCoords[2] || ((cAV[2]+1)*Zarr.get_chunks(2) > endCoords[2] && endCoords[2] < Zarr.get_shape(2))))){
                cRegion = parallelReadZarrWriteWrapper(Zarr, crop,
                                                  {((cAV[0])*Zarr.get_chunks(0)),
                                                   ((cAV[1])*Zarr.get_chunks(1)),
                                                   ((cAV[2])*Zarr.get_chunks(2))},
                                                  {(cAV[0]+1)*Zarr.get_chunks(0),
                                                   (cAV[1]+1)*Zarr.get_chunks(1),
                                                   (cAV[2]+1)*Zarr.get_chunks(2)});
                if(!cRegion){
                    err = 1;
                    errString = "Error in Writer Read. Chunk: "+Zarr.get_chunkNames(f)+"\n";
                    break;
                }
            }
            if(Zarr.get_order() == "F"){
                for(int64_t z = cAV[2]*Zarr.get_chunks(2); z < (cAV[2]+1)*Zarr.get_chunks(2); z++){
                    if(z>=endCoords[2]){
                        if(crop){
                            if((cAV[2]+1)*Zarr.get_chunks(2) > Zarr.get_shape(2)){
                                memcpy((uint8_t*)chunkUnC+((((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)cRegion+((((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),((Zarr.get_shape(2)-z)*Zarr.get_chunks(0)*Zarr.get_chunks(1))*bytes);
                                uint64_t zRest = ((cAV[2]+1)*Zarr.get_chunks(2))-Zarr.get_shape(2);
                                memset((uint8_t*)chunkUnC+(((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1))*bytes),stoi(Zarr.get_fill_value()),(zRest*(Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes);
                            }
                            else{
                                memcpy((uint8_t*)chunkUnC+((((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)cRegion+((((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),((((cAV[2]+1)*Zarr.get_chunks(2))-z)*Zarr.get_chunks(0)*Zarr.get_chunks(1))*bytes);
                            }
                        }
                        else{
                            uint64_t zRest = ((cAV[2]+1)*Zarr.get_chunks(2))-z;
                            memset((uint8_t*)chunkUnC+(((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1))*bytes),stoi(Zarr.get_fill_value()),(zRest*(Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes);
                        }
                        break;
                    }
                    else if(z<startCoords[2]){
                        if(crop){
                            memcpy((uint8_t*)chunkUnC+(((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1))*bytes),(uint8_t*)cRegion+(((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1))*bytes),((startCoords[2]-z)*Zarr.get_chunks(0)*Zarr.get_chunks(1))*bytes);
                        }
                        else{
                            memset((uint8_t*)chunkUnC+(((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1))*bytes),stoi(Zarr.get_fill_value()),((startCoords[2]-z)*(Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes);
                        }
                        z = startCoords[2]-1;
                        continue;
                    }
                    for(int64_t y = cAV[1]*Zarr.get_chunks(1); y < (cAV[1]+1)*Zarr.get_chunks(1); y++){
                        if(y>=endCoords[1]){
                            if(crop){
                                if((cAV[1]+1)*Zarr.get_chunks(1) > Zarr.get_shape(1)){
                                    memcpy((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)cRegion+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),((Zarr.get_shape(1)-y)*Zarr.get_chunks(0))*bytes);
                                    uint64_t yRest = ((cAV[1]+1)*Zarr.get_chunks(1))-Zarr.get_shape(1);
                                    memset((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),stoi(Zarr.get_fill_value()),(yRest*(Zarr.get_chunks(0)))*bytes);
                                }
                                else{
                                    memcpy((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)cRegion+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),((((cAV[1]+1)*Zarr.get_chunks(1))-y)*Zarr.get_chunks(0))*bytes);
                                }
                            }
                            else{
                                uint64_t yRest = ((cAV[1]+1)*Zarr.get_chunks(1))-y;
                                memset((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),stoi(Zarr.get_fill_value()),(yRest*Zarr.get_chunks(0))*bytes);
                            }
                            break;
                        }
                        else if(y<startCoords[1]){
                            if(crop){
                                memcpy((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)cRegion+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),((startCoords[1]-y)*Zarr.get_chunks(0))*bytes);
                            }
                            else{
                                memset((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),stoi(Zarr.get_fill_value()),(startCoords[1]-y)*bytes);
                            }
                            y = startCoords[1]-1;
                            continue;
                        }

                        if(((cAV[0]*Zarr.get_chunks(0)) < startCoords[0] && ((cAV[0]+1)*Zarr.get_chunks(0)) > startCoords[0]) || (cAV[0]+1)*Zarr.get_chunks(0)>endCoords[0]){
                            if(((cAV[0]*Zarr.get_chunks(0)) < startCoords[0] && ((cAV[0]+1)*Zarr.get_chunks(0)) > startCoords[0]) && (cAV[0]+1)*Zarr.get_chunks(0)>endCoords[0]){
                                if(crop){
                                    memcpy((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)cRegion+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(startCoords[0]%Zarr.get_chunks(0))*bytes);
                                    memcpy((uint8_t*)chunkUnC+(((startCoords[0]%Zarr.get_chunks(0))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)zarrArr+((((cAV[0]*Zarr.get_chunks(0))-startCoords[0]+(startCoords[0]%Zarr.get_chunks(0)))+((y-startCoords[1])*writeShape[0])+((z-startCoords[2])*writeShape[0]*writeShape[1]))*bytes),((endCoords[0]%Zarr.get_chunks(0))-(startCoords[0]%Zarr.get_chunks(0)))*bytes);
                                    memcpy((uint8_t*)chunkUnC+(((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))+(endCoords[0]%Zarr.get_chunks(0)))*bytes),(uint8_t*)cRegion+(((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))+(endCoords[0]%Zarr.get_chunks(0)))*bytes),(Zarr.get_chunks(0)-(endCoords[0]%Zarr.get_chunks(0)))*bytes);
                                }
                                else{
                                    memset((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),stoi(Zarr.get_fill_value()),(startCoords[0]%Zarr.get_chunks(0))*bytes);
                                    memcpy((uint8_t*)chunkUnC+(((startCoords[0]%Zarr.get_chunks(0))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)zarrArr+((((cAV[0]*Zarr.get_chunks(0))-startCoords[0]+(startCoords[0]%Zarr.get_chunks(0)))+((y-startCoords[1])*writeShape[0])+((z-startCoords[2])*writeShape[0]*writeShape[1]))*bytes),((endCoords[0]%Zarr.get_chunks(0))-(startCoords[0]%Zarr.get_chunks(0)))*bytes);
                                    memset((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))+(endCoords[0]%Zarr.get_chunks(0))*bytes),stoi(Zarr.get_fill_value()),(Zarr.get_chunks(0)-(endCoords[0]%Zarr.get_chunks(0)))*bytes);
                                }
                            }
                            else if((cAV[0]+1)*Zarr.get_chunks(0)>endCoords[0]){
                                if(crop){
                                    memcpy((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)zarrArr+((((cAV[0]*Zarr.get_chunks(0))-startCoords[0])+((y-startCoords[1])*writeShape[0])+((z-startCoords[2])*writeShape[0]*writeShape[1]))*bytes),(endCoords[0]-(cAV[0]*Zarr.get_chunks(0)))*bytes);

                                    if((cAV[0]+1)*Zarr.get_chunks(0) > Zarr.get_shape(0)){
                                        memcpy((uint8_t*)chunkUnC+((((endCoords[0]-(cAV[0]*Zarr.get_chunks(0))))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)cRegion+((((endCoords[0]-(cAV[0]*Zarr.get_chunks(0))))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(Zarr.get_shape(0)-endCoords[0])*bytes);
                                        uint64_t xRest = ((cAV[0]+1)*Zarr.get_chunks(0))-Zarr.get_shape(0);
                                        memset((uint8_t*)chunkUnC+(((Zarr.get_shape(0)-(cAV[0]*Zarr.get_chunks(0)))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),stoi(Zarr.get_fill_value()),(xRest)*bytes);
                                    }
                                    else{
                                        memcpy((uint8_t*)chunkUnC+((((endCoords[0]-(cAV[0]*Zarr.get_chunks(0))))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)cRegion+((((endCoords[0]-(cAV[0]*Zarr.get_chunks(0))))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(((cAV[0]+1)*Zarr.get_chunks(0))-endCoords[0])*bytes);
                                    }
                                }
                                else{
                                    memcpy((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)zarrArr+((((cAV[0]*Zarr.get_chunks(0))-startCoords[0])+((y-startCoords[1])*writeShape[0])+((z-startCoords[2])*writeShape[0]*writeShape[1]))*bytes),(endCoords[0]%Zarr.get_chunks(0))*bytes);
                                    memset((uint8_t*)chunkUnC+(((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))+(endCoords[0]%Zarr.get_chunks(0)))*bytes),stoi(Zarr.get_fill_value()),(Zarr.get_chunks(0)-(endCoords[0]%Zarr.get_chunks(0)))*bytes);
                                }
                            }
                            else if((cAV[0]*Zarr.get_chunks(0)) < startCoords[0] && ((cAV[0]+1)*Zarr.get_chunks(0)) > startCoords[0]){
                                if(crop){
                                    memcpy((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)cRegion+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(startCoords[0]%Zarr.get_chunks(0))*bytes);
                                    memcpy((uint8_t*)chunkUnC+(((startCoords[0]%Zarr.get_chunks(0))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)zarrArr+((((cAV[0]*Zarr.get_chunks(0))-startCoords[0]+(startCoords[0]%Zarr.get_chunks(0)))+((y-startCoords[1])*writeShape[0])+((z-startCoords[2])*writeShape[0]*writeShape[1]))*bytes),(Zarr.get_chunks(0)-(startCoords[0]%Zarr.get_chunks(0)))*bytes);
                                }
                                else{
                                    memset((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),stoi(Zarr.get_fill_value()),(startCoords[0]%Zarr.get_chunks(0))*bytes);
                                    memcpy((uint8_t*)chunkUnC+(((startCoords[0]%Zarr.get_chunks(0))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)zarrArr+((((cAV[0]*Zarr.get_chunks(0))-startCoords[0]+(startCoords[0]%Zarr.get_chunks(0)))+((y-startCoords[1])*writeShape[0])+((z-startCoords[2])*writeShape[0]*writeShape[1]))*bytes),(Zarr.get_chunks(0)-(startCoords[0]%Zarr.get_chunks(0)))*bytes);
                                }
                            }
                        }
                        else{
                            memcpy((uint8_t*)chunkUnC+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(uint8_t*)zarrArr+((((cAV[0]*Zarr.get_chunks(0))-startCoords[0])+((y-startCoords[1])*writeShape[0])+((z-startCoords[2])*writeShape[0]*writeShape[1]))*bytes),Zarr.get_chunks(0)*bytes);
                        }
                    }
                }
            }
            else if (Zarr.get_order() == "C"){
                for(int64_t x = cAV[0]*Zarr.get_chunks(2); x < (cAV[0]+1)*Zarr.get_chunks(0); x++){
                    for(int64_t y = cAV[1]*Zarr.get_chunks(1); y < (cAV[1]+1)*Zarr.get_chunks(1); y++){
                        for(int64_t z = cAV[2]*Zarr.get_chunks(2); z < (cAV[2]+1)*Zarr.get_chunks(2); z++){
                            switch(bytes){
                                case 1:
                                    if(x>=endCoords[0] || x<startCoords[0] || y>= endCoords[1] || y<startCoords[1] || z>=endCoords[2] || z<startCoords[2]){
                                        ((uint8_t*)chunkUnC)[(((z%Zarr.get_chunks(2))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(2))+((x%Zarr.get_chunks(0))*Zarr.get_chunks(2)*Zarr.get_chunks(1)))*bytes)] = stoi(Zarr.get_fill_value());
                                        continue;
                                    }
                                        ((uint8_t*)chunkUnC)[(((z%Zarr.get_chunks(2))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(2))+((x%Zarr.get_chunks(0))*Zarr.get_chunks(2)*Zarr.get_chunks(1)))*bytes)] = ((uint8_t*)zarrArr)[((x+(y*writeShape[0])+(z*writeShape[0]*writeShape[1]))*bytes)];
                                        break;
                                case 2:
                                    if(x>=endCoords[0] || x<startCoords[0] || y>= endCoords[1] || y<startCoords[1] || z>=endCoords[2] || z<startCoords[2]){
                                        ((uint16_t*)chunkUnC)[(((z%Zarr.get_chunks(2))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(2))+((x%Zarr.get_chunks(0))*Zarr.get_chunks(2)*Zarr.get_chunks(1)))*bytes)] = stoi(Zarr.get_fill_value());
                                        continue;
                                    }
                                        ((uint16_t*)chunkUnC)[(((z%Zarr.get_chunks(2))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(2))+((x%Zarr.get_chunks(0))*Zarr.get_chunks(2)*Zarr.get_chunks(1)))*bytes)] = ((uint16_t*)zarrArr)[((x+(y*writeShape[0])+(z*writeShape[0]*writeShape[1]))*bytes)];
                                        break;
                                case 4:
                                    if(x>=endCoords[0] || x<startCoords[0] || y>= endCoords[1] || y<startCoords[1] || z>=endCoords[2] || z<startCoords[2]){
                                        ((float*)chunkUnC)[(((z%Zarr.get_chunks(2))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(2))+((x%Zarr.get_chunks(0))*Zarr.get_chunks(2)*Zarr.get_chunks(1)))*bytes)] = stoi(Zarr.get_fill_value());
                                        continue;
                                    }
                                        ((float*)chunkUnC)[(((z%Zarr.get_chunks(2))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(2))+((x%Zarr.get_chunks(0))*Zarr.get_chunks(2)*Zarr.get_chunks(1)))*bytes)] = ((float*)zarrArr)[((x+(y*writeShape[0])+(z*writeShape[0]*writeShape[1]))*bytes)];
                                        break;
                                case 8:
                                    if(x>=endCoords[0] || x<startCoords[0] || y>= endCoords[1] || y<startCoords[1] || z>=endCoords[2] || z<startCoords[2]){
                                        ((double*)chunkUnC)[(((z%Zarr.get_chunks(2))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(2))+((x%Zarr.get_chunks(0))*Zarr.get_chunks(2)*Zarr.get_chunks(1)))*bytes)] = stoi(Zarr.get_fill_value());
                                        continue;
                                    }
                                        ((double*)chunkUnC)[(((z%Zarr.get_chunks(2))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(2))+((x%Zarr.get_chunks(0))*Zarr.get_chunks(2)*Zarr.get_chunks(1)))*bytes)] = ((double*)zarrArr)[((x+(y*writeShape[0])+(z*writeShape[0]*writeShape[1]))*bytes)];
                                        break;
                            }

                        }
                    }
                }
            }

            if(sparse){
                const bool allZeros = memcmp(zeroChunkUnc,chunkUnC,sB);
                if(!allZeros){
                    if(Zarr.get_shard()){
                        shardFooter[(currChunk*2)] = std::numeric_limits<uint64_t>::max();
                        shardFooter[(currChunk*2)+1] = std::numeric_limits<uint64_t>::max();
                    } 
                    free(cRegion);
                    cRegion = nullptr;
                    continue;
                }

            }

            // Use the same blosc compress as Zarr
            const std::string subfolderName = Zarr.get_subfoldersString(cAV);
            int64_t csize = 0;

            if(Zarr.get_cname() != "gzip"){
                /*
                if(numWorkers<=Zarr.get_numChunks()){
                    csize = blosc_compress_ctx(Zarr.get_clevel(), BLOSC_SHUFFLE, bytes, sB, chunkUnC, chunkC, sB+BLOSC_MAX_OVERHEAD,Zarr.get_cname().c_str(),0,1);
                }
                else{
                    csize = blosc_compress_ctx(Zarr.get_clevel(), BLOSC_SHUFFLE, bytes, sB, chunkUnC, chunkC, sB+BLOSC_MAX_OVERHEAD,Zarr.get_cname().c_str(),0,numWorkers);
                }
                */
                csize = blosc_compress_ctx(Zarr.get_clevel(), BLOSC_SHUFFLE, bytes, sB, chunkUnC, chunkC, sB+BLOSC_MAX_OVERHEAD,Zarr.get_cname().c_str(),0,nBloscThreads);
            }
            else{
                csize = sB+BLOSC_MAX_OVERHEAD;
                z_stream stream;
                stream.zalloc = Z_NULL;
                stream.zfree = Z_NULL;
                stream.opaque = Z_NULL;

                stream.next_in = (uint8_t*)chunkUnC;
                stream.next_out = (uint8_t*)chunkC;

                stream.avail_in = sB;
                stream.avail_out = csize;
                int cErr = deflateInit2(&stream, Zarr.get_clevel(), Z_DEFLATED, MAX_WBITS + 16, MAX_MEM_LEVEL, Z_DEFAULT_STRATEGY);
                if(cErr){
                    #pragma omp critical
                    {
                        err = 1;
                        errString = "Compression error. Error code: "+
                            std::to_string(cErr)+" ChunkName: "+
                            Zarr.get_fileName()+"/"+subfolderName+"/"+
                            Zarr.get_chunkNames(f)+"\n";
                    }
                    break;
                }

                cErr = deflate(&stream, Z_FINISH);

                if(cErr != Z_STREAM_END){
                    #pragma omp critical
                    {
                        err = 1;
                        errString = "Compression error. Error code: "+
                            std::to_string(cErr)+" ChunkName: "+
                            Zarr.get_fileName()+"/"+subfolderName+"/"+
                            Zarr.get_chunkNames(f)+"\n";                }
                    break;
                }

                if(deflateEnd(&stream)){
                    #pragma omp critical
                    {
                        err = 1;
                        errString = "Compression error. Error code: "+
                            std::to_string(cErr)+" ChunkName: "+
                            Zarr.get_fileName()+"/"+subfolderName+"/"+
                            Zarr.get_chunkNames(f)+"\n";
                    }
                    break;
                }
                csize = csize - stream.avail_out;
            }
            
            // Default write
            if(!Zarr.get_shard()){
                std::string fileName(Zarr.get_fileName()+"/"+subfolderName+"/"+Zarr.get_chunkNames(f));
                std::string fileNameFinal;
                if(useUuid){
                    fileNameFinal = std::string(fileName);
                    fileName.append(uuid);
                }
                    std::ofstream file(fileName, std::ios::binary | std::ios::trunc);
    
                    if(!file.is_open()){
                        #pragma omp critical
                        {
                            err = 1;
                            errString = "Check permissions or filepath. Cannot write to path: "+
                                fileName+"\n";
                        }
                        break;
                    }
                    file.write(reinterpret_cast<char*>(chunkC),csize);
                    file.close();
                if(useUuid){
                    rename(fileName.c_str(),fileNameFinal.c_str());
                }
            }
            // Sharding
            else{
                std::string fileName(Zarr.get_fileName()+"/"+subfolderName+"/"+Zarr.chunkNameToShardName(Zarr.get_chunkNames(f)));
                if(useUuid){    
                    fileName.append(uuid);
                }
                std::ofstream file;

                if(currChunk > 0) {
                    int64_t currInd = currChunk;
                    while(currInd >= 0 && (shardFooter[((currInd-1)*2)] == std::numeric_limits<uint64_t>::max() &&
                       shardFooter[((currInd-1)*2)+1] == std::numeric_limits<uint64_t>::max())){
                       currInd--;
                    }
                    uint64_t shardOffset = 0; 
                    if(currInd){
                        shardOffset = shardFooter[((currInd-1)*2)]+
                                      shardFooter[((currInd-1)*2)+1];
                    }
                    shardFooter[(currChunk*2)] = shardOffset;
                    shardFooter[(currChunk*2)+1] = csize;
                    file = std::ofstream(fileName, std::ios::binary | std::ios::app);

                }
                else{
                    shardFooter[0] = 0;
                    shardFooter[1] = csize;
                    file = std::ofstream(fileName, std::ios::binary | std::ios::trunc);
                }
                if(!file.is_open()){
                    #pragma omp critical
                    {
                        err = 1;
                        errString = "Check permissions or filepath. Cannot write to path: "+
                            fileName+"\n";
                    }
                    break;
                }
                file.write(reinterpret_cast<char*>(chunkC),csize);
                file.close();
            
            }
            free(cRegion);
            cRegion = nullptr;
        }

        if(Zarr.get_shard() && unWritten && (((uint64_t)((w)*batchSize)))<(Zarr.get_numChunks())){
            // calculate CRC32C
            uint32_t shardFooterCRC32C = crc32c(reinterpret_cast<uint8_t*>(shardFooter), shardFooterSize*sizeof(uint64_t));
            uint64_t f = lastF;
            std::vector<uint64_t> pAV = Zarr.get_chunkAxisVals(Zarr.get_chunkNames(f));

            bool pad = pAV[0] > endCoords[0]/Zarr.get_chunk_shape(0) ||
                    pAV[1] > endCoords[1]/Zarr.get_chunk_shape(1) ||
                    pAV[2] > endCoords[2]/Zarr.get_chunk_shape(2);

            if(pad){
                for(uint64_t i = currChunk; i < Zarr.get_numChunksPerShard(); i++){
                    shardFooter[i*2] = std::numeric_limits<uint64_t>::max();
                    shardFooter[(i*2)+1] = std::numeric_limits<uint64_t>::max();
                }
            }
            const std::string subfolderName = Zarr.get_subfoldersString(cAV);
            std::string fileName(Zarr.get_fileName()+"/"+subfolderName+"/"+Zarr.chunkNameToShardName(Zarr.get_chunkNames(f)));
            std::string fileNameFinal;
            if(useUuid){
                fileNameFinal = std::string(fileName);
                fileName.append(uuid);
            }
            
            std::ofstream file(fileName, std::ios::binary | std::ios::app);

            if(!file.is_open()){
                    #pragma omp critical
                    {
                        err = 1;
                        errString = "Check permissions or filepath. Cannot write to path: "+
                            fileName+"\n";
                    }
                    continue;
            }
            file.write(reinterpret_cast<char*>(shardFooter),shardFooterSize*sizeof(uint64_t));
            file.write(reinterpret_cast<char*>(&shardFooterCRC32C),sizeof(uint32_t));
            file.close();
            if(useUuid){
                rename(fileName.c_str(),fileNameFinal.c_str());
            }
        }
        free(shardFooter);
        free(chunkUnC);
        free(chunkC);

    }
    free(zeroChunkUnc);

    if(err) {
        Zarr.set_errString(errString);
        return 1;
    }
    else return 0;
}
