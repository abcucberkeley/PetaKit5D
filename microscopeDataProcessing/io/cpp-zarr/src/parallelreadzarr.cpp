#include <fstream>
#include <stdint.h>
#include <omp.h>
#include "parallelreadzarr.h"
#include "blosc2.h"
//#include "blosc.h"
#include "zarr.h"
#include "helperfunctions.h"
#include "zlib.h"

//#include <iostream>

// zarrArr should be initialized to all zeros if you have empty chunks
uint8_t parallelReadZarr(zarr &Zarr, void* zarrArr,
                         const std::vector<uint64_t> &startCoords, 
                         const std::vector<uint64_t> &endCoords,
                         const std::vector<uint64_t> &readShape,
                         const uint64_t bits,
                         const bool useCtx,
                         const bool sparse)
{
    //std::cout << "fileName: " << Zarr.get_fileName() << std::endl;
    //std::cout << "startCoords: " << startCoords[0] << " " << startCoords[1] << " " << startCoords[2] << std::endl;
    //std::cout << "endCoords: " << endCoords[0] << " " << endCoords[1] << " " << endCoords[2] << std::endl;
    //std::cout << "readShape: " << readShape[0] << " " << readShape[1] << " " << readShape[2] << std::endl;
    void* zarrArrC = nullptr;
    const uint64_t bytes = (bits/8);
    
    int32_t numWorkers = omp_get_max_threads();

    // nBloscThreads used when using blosc_ctx
    uint32_t nBloscThreads = 1;
    if(!useCtx){
        blosc2_init();
        blosc2_set_nthreads(numWorkers);
    }

    // no ctx
    /*
    if(numWorkers>Zarr.get_numChunks()){
        blosc_set_nthreads(std::ceil(((double)numWorkers)/((double)Zarr.get_numChunks())));
        numWorkers = Zarr.get_numChunks();
    }
    else {
        blosc_set_nthreads(numWorkers);
    }
    */
    
    // ctx
    /*
    if(numWorkers>Zarr.get_numChunks()){
        nBloscThreads = std::ceil(((double)numWorkers)/((double)Zarr.get_numChunks()));
        numWorkers = Zarr.get_numChunks();
    }
    */

    // The chunk size is actually the inner chunk size if the zarr file is sharded
    if(Zarr.get_shard()){
        Zarr.set_chunks({Zarr.get_chunk_shape(0),Zarr.get_chunk_shape(1),Zarr.get_chunk_shape(2)});
    }
    
    const int32_t batchSize = (Zarr.get_numChunks()-1)/numWorkers+1;
    const uint64_t s = Zarr.get_chunks(0)*Zarr.get_chunks(1)*Zarr.get_chunks(2);
    const uint64_t sB = s*bytes;

    // If C->F order then we need a temp C order array
    if(Zarr.get_order() == "C") zarrArrC = calloc(readShape[0]*readShape[1]*readShape[2],bytes);
    
    void* zeroChunkUnc = NULL;
    if(sparse){
        zeroChunkUnc = calloc(s,bytes);
    }

    int err = 0;
    std::string errString;
    //return -1;
    #pragma omp parallel for
    for(int32_t w = 0; w < numWorkers; w++){
        void* bufferDest = operator new(sB);
        void* buffer = NULL;
        std::streamsize lastFileLen = 0;
        int64_t dsize = -1;
        int uncErr = 0;
        for(int64_t f = w*batchSize; f < (w+1)*batchSize; f++){
            if(f>=Zarr.get_numChunks() || err) break;
            const std::vector<uint64_t> cAV = Zarr.get_chunkAxisVals(Zarr.get_chunkNames(f));
            const std::string subfolderName = Zarr.get_subfoldersString(cAV);

            std::string fileName;
            
            if(!Zarr.get_shard()){
                fileName = (Zarr.get_fileName()+"/"+subfolderName+"/"+Zarr.get_chunkNames(f));
            }
            else{
                // Can change this to the check for zeros maybe
                /*
                bool pad = (cAV[0] >= ceil((double)readShape[0]/(double)Zarr.get_chunk_shape(0)) ||
                    cAV[1] >= ceil((double)readShape[1]/(double)Zarr.get_chunk_shape(1)) ||
                    cAV[2] >= ceil((double)readShape[2]/(double)Zarr.get_chunk_shape(2)));
                */
                /*
                bool pad = (cAV[0] >= ceil((double)endCoords[0]/(double)Zarr.get_chunk_shape(0)) ||
                    cAV[1] >= ceil((double)endCoords[1]/(double)Zarr.get_chunk_shape(1)) ||
                    cAV[2] >= ceil((double)endCoords[2]/(double)Zarr.get_chunk_shape(2)) ||
                    cAV[0] < ceil((double)startCoords[0]/(double)Zarr.get_chunk_shape(0)) ||
                    cAV[1] < ceil((double)startCoords[1]/(double)Zarr.get_chunk_shape(1)) ||
                    cAV[2] < ceil((double)startCoords[2]/(double)Zarr.get_chunk_shape(2)));
                */
                bool pad = cAV[0] > endCoords[0]/Zarr.get_chunk_shape(0) ||
                    cAV[1] > endCoords[1]/Zarr.get_chunk_shape(1) ||
                    cAV[2] > endCoords[2]/Zarr.get_chunk_shape(2) ||
                    cAV[0] < startCoords[0]/Zarr.get_chunk_shape(0) ||
                    cAV[1] < startCoords[1]/Zarr.get_chunk_shape(1) ||
                    cAV[2] < startCoords[2]/Zarr.get_chunk_shape(2);
                /*
                    std::cout <<  "cAV[0]: " << cAV[0] << " endC0: " << endCoords[0]/Zarr.get_chunk_shape(0) <<
                        " cAV[1]: " << cAV[1] << " endC1: " << endCoords[1]/Zarr.get_chunk_shape(1) <<
                        " cAv[2]: " << cAV[2] << " rnfC2: " << endCoords[2]/Zarr.get_chunk_shape(2) <<
                        " startC0: " << startCoords[0]/Zarr.get_chunk_shape(0) << 
                        " startC1: " << startCoords[1]/Zarr.get_chunk_shape(1) <<
                        " startC2: " << startCoords[2]/Zarr.get_chunk_shape(2) << std::endl;
                */
                if(pad) {
                    //std::cout << "PADDING" << std::endl;
                    /*
                    std::cout << "PADDING: "<< "cAV[0]: " << cAV[0] << " Zarr.get_chunks(0): " << Zarr.get_chunks(0) <<
                        " startCoords[0]: " << startCoords[0] << " startCoords[1]: " << startCoords[1] <<
                        " readShape[0]: " << readShape[0] << " startCoords[2]: " << startCoords[2] <<
                        " readShape[1]: " << readShape[1] << " bytes: " << bytes << std::endl;
                    */
                    continue;
                }
                fileName = Zarr.get_fileName()+"/"+subfolderName+"/"+Zarr.chunkNameToShardName(Zarr.get_chunkNames(f));
            }
            // If we cannot open the file then set to all zeros
            // Can make this better by checking the errno
            std::ifstream file(fileName, std::ios::binary);
            if(!file.is_open()){
                continue;
                //memset(bufferDest,0,sB);
            }
            else{
                std::streamsize fileLen; 
                if(!Zarr.get_shard()){
                    file.seekg(0, std::ios::end);
                    fileLen = file.tellg();
                    if(lastFileLen < fileLen){
                        operator delete(buffer);
                        buffer = operator new(fileLen);
                        lastFileLen = fileLen;
                    }
                    file.seekg(0, std::ios::beg);
                    file.read(reinterpret_cast<char*>(buffer), fileLen);
                    file.close();
                }
                // Sharding
                else{
                    uint64_t currChunkShardPosition = Zarr.get_chunkShardPosition(cAV);
                    uint64_t offsetNBytes[2];
                    
                    file.seekg(-(int64_t)(((Zarr.get_numChunksPerShard()*2*sizeof(uint64_t))+4)-(currChunkShardPosition*2*sizeof(uint64_t))), std::ios::end);
                    file.read(reinterpret_cast<char*>(offsetNBytes), sizeof(offsetNBytes));
                    
                    // All zeros or skippable chunk
                    if(offsetNBytes[0]== std::numeric_limits<uint64_t>::max() &&
                       offsetNBytes[1] == std::numeric_limits<uint64_t>::max()){
                        file.close();
                        continue;
                    }
                    fileLen = offsetNBytes[1];
                    //std::cout << "offset: " << offsetNBytes[0] << " lastFileLen: " << lastFileLen << " fileLen: " << fileLen << std::endl;
                    if(lastFileLen < fileLen){
                        operator delete(buffer);
                        buffer = operator new(fileLen);
                        lastFileLen = fileLen;
                    }
                    file.seekg(offsetNBytes[0], std::ios::beg);
                    file.read(reinterpret_cast<char*>(buffer), fileLen);

                    file.close();
                }
                
                // Decompress
                if(Zarr.get_cname() != "gzip"){
                    if(!useCtx){
                        dsize = blosc2_decompress(buffer, fileLen, bufferDest, sB);
                    }
                    else{
                        blosc2_context *dctx;
                        blosc2_dparams dparams = {(int16_t)nBloscThreads,NULL,NULL,NULL};
                        dctx = blosc2_create_dctx(dparams);
                        dsize = blosc2_decompress_ctx(dctx, buffer, fileLen, bufferDest, sB);
                        blosc2_free_ctx(dctx);
                    }
                }
                else{
                    dsize = sB;
                    z_stream stream;
                    stream.zalloc = Z_NULL;
                    stream.zfree = Z_NULL;
                    stream.opaque = Z_NULL;
                    stream.avail_in = (uInt)fileLen;
                    stream.avail_out = (uInt)dsize;
                    while(stream.avail_in > 0){
    
                        dsize = sB;
    
                        stream.next_in = (uint8_t*)buffer+(fileLen-stream.avail_in);
                        stream.next_out = (uint8_t*)bufferDest+(sB-stream.avail_out);
    
                        uncErr = inflateInit2(&stream, 32);
                        if(uncErr){
                        #pragma omp critical
                        {
                        err = 1;
                        errString = "Decompression error. Error code: "+
                                     std::to_string(uncErr)+" ChunkName: "+
                                     Zarr.get_fileName()+"/"+subfolderName+"/"+
                                     Zarr.get_chunkNames(f)+"\n";
                        }
                        break;
                        }
        
                        uncErr = inflate(&stream, Z_NO_FLUSH);
        
                        if(uncErr != Z_STREAM_END){
                        #pragma omp critical
                        {
                        err = 1;
                        errString = "Decompression error. Error code: "+
                                     std::to_string(uncErr)+" ChunkName: "+
                                     Zarr.get_fileName()+"/"+subfolderName+"/"+
                                     Zarr.get_chunkNames(f)+"\n";
                        }
                        break;
                        }
                    }
                    if(inflateEnd(&stream)){
                        #pragma omp critical
                        {
                        err = 1;
                        errString = "Decompression error. Error code: "+
                                     std::to_string(uncErr)+" ChunkName: "+
                                     Zarr.get_fileName()+"/"+subfolderName+"/"+
                                     Zarr.get_chunkNames(f)+"\n";
                        }
                        break;
                    }
                }
                
                
                if(dsize < 0){
                    #pragma omp critical
                    {
                    err = 1;
                    errString = "Decompression error. Error code: "+
                                     std::to_string(uncErr)+" ChunkName: "+
                                     Zarr.get_fileName()+"/"+subfolderName+"/"+
                                     Zarr.get_chunkNames(f)+"\n";
                    }
                    break;
                }
            }
            if(sparse){
                // If the chunk is all zeros (memcmp == 0) then we skip it
                const bool allZeros = memcmp(zeroChunkUnc,bufferDest,sB);
                if(!allZeros) continue;
            }
            
            // F->F
            /*
             std::cout << "cAV[0]: " << cAV[0] << " Zarr.get_chunks(0): " << Zarr.get_chunks(0) <<
                 " startCoords[0]: " << startCoords[0] << " startCoords[1]: " << startCoords[1] <<
                 " readShape[0]: " << readShape[0] << " startCoords[2]: " << startCoords[2] <<
                 " readShape[1]: " << readShape[1] << " bytes: " << bytes << std::endl;
        */
            
            if(Zarr.get_order() == "F"){  
                for(int64_t y = cAV[1]*Zarr.get_chunks(1); y < (cAV[1]+1)*Zarr.get_chunks(1); y++){
                    if(y>=endCoords[1]) break;
                    else if(y<startCoords[1]) continue;
                    for(int64_t z = cAV[2]*Zarr.get_chunks(2); z < (cAV[2]+1)*Zarr.get_chunks(2); z++){
                        if(z>=endCoords[2]) break;
                        else if(z<startCoords[2]) continue;
                        if(((cAV[0]*Zarr.get_chunks(0)) < startCoords[0] && ((cAV[0]+1)*Zarr.get_chunks(0)) > startCoords[0]) || (cAV[0]+1)*Zarr.get_chunks(0)>endCoords[0]){
                            if(((cAV[0]*Zarr.get_chunks(0)) < startCoords[0] && ((cAV[0]+1)*Zarr.get_chunks(0)) > startCoords[0]) && (cAV[0]+1)*Zarr.get_chunks(0)>endCoords[0]){
                                memcpy((uint8_t*)zarrArr+((((cAV[0]*Zarr.get_chunks(0))-startCoords[0]+(startCoords[0]%Zarr.get_chunks(0)))+((y-startCoords[1])*readShape[0])+((z-startCoords[2])*readShape[0]*readShape[1]))*bytes),(uint8_t*)bufferDest+(((startCoords[0]%Zarr.get_chunks(0))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),((endCoords[0]%Zarr.get_chunks(0))-(startCoords[0]%Zarr.get_chunks(0)))*bytes);
                            }
                            else if((cAV[0]+1)*Zarr.get_chunks(0)>endCoords[0]){
                                memcpy((uint8_t*)zarrArr+((((cAV[0]*Zarr.get_chunks(0))-startCoords[0])+((y-startCoords[1])*readShape[0])+((z-startCoords[2])*readShape[0]*readShape[1]))*bytes),(uint8_t*)bufferDest+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(endCoords[0]%Zarr.get_chunks(0))*bytes);
                            }
                            else if((cAV[0]*Zarr.get_chunks(0)) < startCoords[0] && ((cAV[0]+1)*Zarr.get_chunks(0)) > startCoords[0]){
                                memcpy((uint8_t*)zarrArr+((((cAV[0]*Zarr.get_chunks(0)-startCoords[0]+(startCoords[0]%Zarr.get_chunks(0))))+((y-startCoords[1])*readShape[0])+((z-startCoords[2])*readShape[0]*readShape[1]))*bytes),(uint8_t*)bufferDest+(((startCoords[0]%Zarr.get_chunks(0))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),(Zarr.get_chunks(0)-(startCoords[0]%Zarr.get_chunks(0)))*bytes);
                            }
                        }
                        else{
                            memcpy((uint8_t*)zarrArr+((((cAV[0]*Zarr.get_chunks(0))-startCoords[0])+((y-startCoords[1])*readShape[0])+((z-startCoords[2])*readShape[0]*readShape[1]))*bytes),(uint8_t*)bufferDest+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(0))+((z%Zarr.get_chunks(2))*Zarr.get_chunks(0)*Zarr.get_chunks(1)))*bytes),Zarr.get_chunks(0)*bytes);
                        }
                    }
                }
                
            }
            // C->C (x and z are flipped) then we flip to F below
            else if (Zarr.get_order() == "C"){
                for(int64_t y = cAV[1]*Zarr.get_chunks(1); y < (cAV[1]+1)*Zarr.get_chunks(1); y++){
                    if(y>=endCoords[1]) break;
                    else if(y<startCoords[1]) continue;
                    for(int64_t z = cAV[0]*Zarr.get_chunks(0); z < (cAV[0]+1)*Zarr.get_chunks(0); z++){
                        if(z>=endCoords[0]) break;
                        else if(z<startCoords[0]) continue;
                        if(((cAV[2]*Zarr.get_chunks(2)) < startCoords[2] && ((cAV[2]+1)*Zarr.get_chunks(2)) > startCoords[2]) || (cAV[2]+1)*Zarr.get_chunks(2)>endCoords[2]){
                            if(((cAV[2]*Zarr.get_chunks(2)) < startCoords[2] && ((cAV[2]+1)*Zarr.get_chunks(2)) > startCoords[2]) && (cAV[2]+1)*Zarr.get_chunks(2)>endCoords[2]){
                                memcpy((uint8_t*)zarrArrC+((((cAV[2]*Zarr.get_chunks(2))-startCoords[2]+(startCoords[2]%Zarr.get_chunks(2)))+((y-startCoords[1])*readShape[2])+((z-startCoords[0])*readShape[2]*readShape[1]))*bytes),(uint8_t*)bufferDest+(((startCoords[2]%Zarr.get_chunks(2))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(2))+((z%Zarr.get_chunks(0))*Zarr.get_chunks(2)*Zarr.get_chunks(1)))*bytes),((endCoords[2]%Zarr.get_chunks(2))-(startCoords[2]%Zarr.get_chunks(2)))*bytes);
                            }
                            else if((cAV[2]+1)*Zarr.get_chunks(2)>endCoords[2]){
                                memcpy((uint8_t*)zarrArrC+((((cAV[2]*Zarr.get_chunks(2))-startCoords[2])+((y-startCoords[1])*readShape[2])+((z-startCoords[0])*readShape[2]*readShape[1]))*bytes),(uint8_t*)bufferDest+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(2))+((z%Zarr.get_chunks(0))*Zarr.get_chunks(2)*Zarr.get_chunks(1)))*bytes),(endCoords[2]%Zarr.get_chunks(2))*bytes);
                            }
                            else if((cAV[2]*Zarr.get_chunks(2)) < startCoords[2] && ((cAV[2]+1)*Zarr.get_chunks(2)) > startCoords[2]){
                                memcpy((uint8_t*)zarrArrC+((((cAV[2]*Zarr.get_chunks(2)-startCoords[2]+(startCoords[2]%Zarr.get_chunks(2))))+((y-startCoords[1])*readShape[2])+((z-startCoords[0])*readShape[2]*readShape[1]))*bytes),(uint8_t*)bufferDest+(((startCoords[2]%Zarr.get_chunks(2))+((y%Zarr.get_chunks(1))*Zarr.get_chunks(2))+((z%Zarr.get_chunks(0))*Zarr.get_chunks(2)*Zarr.get_chunks(1)))*bytes),(Zarr.get_chunks(2)-(startCoords[2]%Zarr.get_chunks(2)))*bytes);
                            }
                        }
                        else{
                            memcpy((uint8_t*)zarrArrC+((((cAV[2]*Zarr.get_chunks(2))-startCoords[2])+((y-startCoords[1])*readShape[2])+((z-startCoords[0])*readShape[2]*readShape[1]))*bytes),(uint8_t*)bufferDest+((((y%Zarr.get_chunks(1))*Zarr.get_chunks(2))+((z%Zarr.get_chunks(0))*Zarr.get_chunks(2)*Zarr.get_chunks(1)))*bytes),Zarr.get_chunks(2)*bytes);
                        }
                    }
                }  
                
            }
            
        }
        operator delete(bufferDest);
        operator delete(buffer);
    }
    if(!useCtx){
        blosc2_destroy();
    }
    free(zeroChunkUnc);

    if(err){
        Zarr.set_errString(errString);
        free(zarrArrC);
        return 1;
    }
    else if (Zarr.get_order() == "C"){
        // This transpose can potentially be optimized more        
        #pragma omp parallel for collapse(3)
        for(size_t j = 0; j < readShape[1]; j++) {
            for(size_t i = 0; i < readShape[0]; i++) {
                for(size_t k = 0; k < readShape[2]; k++) {
                    switch(bytes){
                        case 1:
                            *(((uint8_t*)zarrArr)+(k*readShape[0]*readShape[1])+(j*readShape[0])+i) = *((uint8_t*)zarrArrC+(i*readShape[1]*readShape[2])+(j*readShape[2])+k);
                            break;
                        case 2:
                            *(((uint16_t*)zarrArr)+(k*readShape[0]*readShape[1])+(j*readShape[0])+i) = *((uint16_t*)zarrArrC+(i*readShape[1]*readShape[2])+(j*readShape[2])+k);
                            break;
                        case 4:
                            *(((float*)zarrArr)+(k*readShape[0]*readShape[1])+(j*readShape[0])+i) = *((float*)zarrArrC+(i*readShape[1]*readShape[2])+(j*readShape[2])+k);
                            break;
                        case 8:
                            *(((double*)zarrArr)+(k*readShape[0]*readShape[1])+(j*readShape[0])+i) = *((double*)zarrArrC+(i*readShape[1]*readShape[2])+(j*readShape[2])+k);
                            break;
                    }
                }
            }
        }
        free(zarrArrC);
    }
    return 0;
}

// TODO: FIX MEMORY LEAKS
// Wrapper used by parallelWriteZarr
void* parallelReadZarrWriteWrapper(zarr Zarr, const bool &crop,
                              std::vector<uint64_t> startCoords, 
                              std::vector<uint64_t> endCoords){
   
    if(!crop){
        startCoords[0] = 0;
        startCoords[1] = 0;
        startCoords[2] = 0;
        endCoords[0] = Zarr.get_shape(0);
        endCoords[1] = Zarr.get_shape(1);
        endCoords[2] = Zarr.get_shape(2);
    }

    
    std::vector<uint64_t> readShape = {endCoords[0]-startCoords[0],
                                       endCoords[1]-startCoords[1],
                                       endCoords[2]-startCoords[2]};

    Zarr.set_chunkInfo(startCoords, endCoords);
    uint8_t err = 0;
    uint64_t readSize = readShape[0]*readShape[1]*readShape[2];
    if(Zarr.get_dtype() == "<u1"){
        uint64_t bits = 8;
        uint8_t* zarrArr = nullptr;
        if(stoi(Zarr.get_fill_value())){
            zarrArr = (uint8_t*)malloc(readSize*sizeof(uint8_t));
            memset(zarrArr,stoi(Zarr.get_fill_value()),readSize*sizeof(uint8_t));
        }
        else zarrArr = (uint8_t*)calloc(readSize,sizeof(uint8_t));
        err = parallelReadZarr(Zarr, (void*)zarrArr,startCoords,endCoords,readShape,bits,true);
        if(err){
            free(zarrArr);
            return NULL;
        }
        else return (void*)zarrArr;
    }
    else if(Zarr.get_dtype() == "<u2"){
        uint64_t bits = 16;
        uint16_t* zarrArr = nullptr;
        if(stoi(Zarr.get_fill_value())){
            zarrArr = (uint16_t*)malloc(readSize*(uint64_t)(sizeof(uint16_t)));
            memset(zarrArr,stoi(Zarr.get_fill_value()),readSize*sizeof(uint16_t));
        }
        else zarrArr = (uint16_t*)calloc(readSize,(uint64_t)(sizeof(uint16_t)));
        err = parallelReadZarr(Zarr, (void*)zarrArr,startCoords,endCoords,readShape,bits,true);
        if(err){
            free(zarrArr);
            return NULL;
        }
        else return (void*)zarrArr;
    }
    else if(Zarr.get_dtype() == "<f4"){
        uint64_t bits = 32;
        float* zarrArr = nullptr;
        if(stoi(Zarr.get_fill_value())){
            zarrArr = (float*)malloc(readSize*(sizeof(float)));
            memset(zarrArr,stoi(Zarr.get_fill_value()),readSize*sizeof(float));
        }
        else zarrArr = (float*)calloc(readSize,(sizeof(float)));
        err = parallelReadZarr(Zarr, (void*)zarrArr,startCoords,endCoords,readShape,bits,true);
        if(err){
            free(zarrArr);
            return NULL;
        }
        else return (void*)zarrArr;
    }
    else if(Zarr.get_dtype() == "<f8"){
        uint64_t bits = 64;
        double* zarrArr = nullptr;
        if(stoi(Zarr.get_fill_value())){
            zarrArr = (double*)malloc(readSize*(sizeof(double)));
            memset(zarrArr,stoi(Zarr.get_fill_value()),readSize*sizeof(double));
        }
        else zarrArr = (double*)calloc(readSize,(sizeof(double)));
        err = parallelReadZarr(Zarr, (void*)zarrArr,startCoords,endCoords,readShape,bits,true);
        if(err){
            free(zarrArr);
            return NULL;
        }
        else return (void*)zarrArr;
    }
    else{
        return NULL;
    }
}

void* readZarrParallelHelper(const char* folderName, uint64_t startX, uint64_t startY, uint64_t startZ, uint64_t endX, uint64_t endY, uint64_t endZ, uint8_t imageJIm){
    zarr Zarr(folderName);
    void* zarrArr = parallelReadZarrWriteWrapper(Zarr, true,
                              {startX, startY, startZ},
                              {endX, endY, endZ});
    // May need to add a check for if the data is f order or c order for ImageJ
    // For the c order images I have tested, we also have to do this flip for now
    if(imageJIm /*&& (order == 'F' || order == 'f')*/){
        void* zarrArrC = malloc(Zarr.get_shape(0)*Zarr.get_shape(1)*Zarr.get_shape(2)*Zarr.dtypeBytes());
		#pragma omp parallel for
        for(uint64_t k = 0; k < Zarr.get_shape(2); k++){
            for(uint64_t j = 0; j < Zarr.get_shape(1); j++){
                for(uint64_t i = 0; i < Zarr.get_shape(0); i++){
                    switch(Zarr.dtypeBytes()){
                        case 1:
                            ((uint8_t*)zarrArrC)[j+(i*Zarr.get_shape(1))+(k*Zarr.get_shape(1)*Zarr.get_shape(0))] = ((uint8_t*)zarrArr)[i+(j*Zarr.get_shape(0))+(k*Zarr.get_shape(1)*Zarr.get_shape(0))];
                            break;
                        case 2:
                            ((uint16_t*)zarrArrC)[j+(i*Zarr.get_shape(1))+(k*Zarr.get_shape(1)*Zarr.get_shape(0))] = ((uint16_t*)zarrArr)[i+(j*Zarr.get_shape(0))+(k*Zarr.get_shape(1)*Zarr.get_shape(0))];
                            break;
                        case 4:
                            ((float*)zarrArrC)[j+(i*Zarr.get_shape(1))+(k*Zarr.get_shape(1)*Zarr.get_shape(0))] = ((float*)zarrArr)[i+(j*Zarr.get_shape(0))+(k*Zarr.get_shape(1)*Zarr.get_shape(0))];
                            break;
                        case 8:
                            ((double*)zarrArrC)[j+(i*Zarr.get_shape(1))+(k*Zarr.get_shape(1)*Zarr.get_shape(0))] = ((double*)zarrArr)[i+(j*Zarr.get_shape(0))+(k*Zarr.get_shape(1)*Zarr.get_shape(0))];
                            break;
                    }
                }
            }
        }
		free(zarrArr);
        return zarrArrC;
    }
    return zarrArr;
}

