#include "parallelReadZarr.h"
#include <stdio.h>
#include <dirent.h>
#include <blosc2.h>
#include <cjson/cJSON.h>
#include <omp.h>
#include <stdlib.h>
#include "helperFunctions.h"

#include "zlib.h"

//mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/include/' '-I/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/lib64' -lblosc2 '-L/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/lib64' -lcjson zarrMex.c

void parallelReadZarr(void* zarr, char* folderName,uint64_t startX, uint64_t startY, uint64_t startZ, uint64_t endX, uint64_t endY,uint64_t endZ,uint64_t chunkXSize,uint64_t chunkYSize,uint64_t chunkZSize,uint64_t shapeX,uint64_t shapeY,uint64_t shapeZ, uint64_t bits, char order, char* cname){
    char fileSepS[2];
    const char fileSep =
#ifdef _WIN32
            '\\';
#else
    '/';
#endif
    
    uint64_t bytes = (bits/8);
    fileSepS[0] = fileSep;
    fileSepS[1] = '\0';
    
    /* Initialize the Blosc compressor */
    int32_t numWorkers = omp_get_max_threads();
    
    struct chunkInfo cI = getChunkInfo(folderName, startX, startY, startZ, endX, endY, endZ,chunkXSize,chunkYSize,chunkZSize);
    //if(!cI.chunkNames) mexErrMsgIdAndTxt("zarr:inputError","File \"%s\" cannot be opened",folderName);
    
    int32_t batchSize = (cI.numChunks-1)/numWorkers+1;
    uint64_t s = chunkXSize*chunkYSize*chunkZSize;
    uint64_t sB = s*bytes;
    int32_t w;
    int err = 0;
    char errString[10000];
    
#pragma omp parallel for //if(numWorkers<=cI.numChunks)
    for(w = 0; w < numWorkers; w++){
        void* bufferDest = mallocDynamic(s,bits);
        uint64_t lastFileLen = 0;
        char *buffer = NULL;
        for(int64_t f = w*batchSize; f < (w+1)*batchSize; f++){
            if(f>=cI.numChunks || err) break;
            struct chunkAxisVals cAV = getChunkAxisVals(cI.chunkNames[f]);
            
            //malloc +2 for null term and filesep
            char *fileName = malloc(strlen(folderName)+strlen(cI.chunkNames[f])+2);
            fileName[0] = '\0';
            strcat(fileName,folderName);
            strcat(fileName,fileSepS);
            strcat(fileName,cI.chunkNames[f]);
            
            FILE *fileptr = fopen(fileName, "rb");
            if(!fileptr){
#pragma omp critical
                {
                err = 1;
                sprintf(errString,"Could not open file: %s\n",fileName);
                }
                break;
            }
            free(fileName);
            
            fseek(fileptr, 0, SEEK_END);
            long filelen = ftell(fileptr);
            rewind(fileptr);
            if(lastFileLen < filelen){
                free(buffer);
                buffer = (char*) malloc(filelen);
                lastFileLen = filelen;
            }
            fread(buffer, filelen, 1, fileptr);
            fclose(fileptr);
            
            // Decompress
            int dsize = -1;
            int uncErr = 0;
            if(strcmp(cname,"gzip")){
                //dsize = blosc2_decompress(buffer, filelen, bufferDest, sB);
                blosc2_context *dctx;
                blosc2_dparams dparams = BLOSC2_DPARAMS_DEFAULTS;
                dctx = blosc2_create_dctx(dparams);
                
                dsize = blosc2_decompress_ctx(dctx, buffer, filelen,bufferDest, sB);
                blosc2_free_ctx(dctx);
            }
            else{
                dsize = sB;
                z_stream stream;
                stream.zalloc = Z_NULL;
                stream.zfree = Z_NULL;
                stream.opaque = Z_NULL;
                stream.avail_in = (uInt)filelen;
                stream.avail_out = (uInt)sB;
                while(stream.avail_in > 0){

                    dsize = sB;

                    stream.next_in = (uint8_t*)buffer+(filelen-stream.avail_in);
                    stream.next_out = (uint8_t*)bufferDest+(sB-stream.avail_out);

                    uncErr = inflateInit2(&stream, 32);
                    if(uncErr){
                    #pragma omp critical
                    {
                    err = 1;
                    sprintf(errString,"Decompression error. Error code: %d ChunkName: %s/%s\n",uncErr,folderName,cI.chunkNames[f]);
                    }
                    break;
                    }
    
                    uncErr = inflate(&stream, Z_NO_FLUSH);
    
                    if(uncErr != Z_STREAM_END){
                    #pragma omp critical
                    {
                    err = 1;
                    sprintf(errString,"Decompression error. Error code: %d ChunkName: %s/%s\n",uncErr,folderName,cI.chunkNames[f]);
                    }
                    break;
                    }
                }
                if(inflateEnd(&stream)){
                    #pragma omp critical
                    {
                    err = 1;
                    sprintf(errString,"Decompression error. Error code: %d ChunkName: %s/%s\n",uncErr,folderName,cI.chunkNames[f]);
                    }
                    break;
                }
            }
            
            
            if(dsize < 0){
#pragma omp critical
                {
                err = 1;
                sprintf(errString,"Decompression error. Error code: %d ChunkName: %s/%s\n",dsize,folderName,cI.chunkNames[f]);
                }
                break;
            }
            
            //printf("ChunkName: %s\n",cI.chunkNames[f]);
            //printf("w: %d b: %d\n",w,f);
            if(order == 'F'){
                for(int64_t z = cAV.z*chunkZSize; z < (cAV.z+1)*chunkZSize; z++){
                    if(z>=endZ) break;
                    else if(z<startZ) continue;
                    for(int64_t y = cAV.y*chunkYSize; y < (cAV.y+1)*chunkYSize; y++){
                        if(y>=endY) break;
                        else if(y<startY) continue;
                        if(((cAV.x*chunkXSize) < startX && ((cAV.x+1)*chunkXSize) > startX) || (cAV.x+1)*chunkXSize>endX){
                            if(((cAV.x*chunkXSize) < startX && ((cAV.x+1)*chunkXSize) > startX) && (cAV.x+1)*chunkXSize>endX){
                                memcpy((uint8_t*)zarr+((((cAV.x*chunkXSize)-startX+(startX%chunkXSize))+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))*bytes),(uint8_t*)bufferDest+(((startX%chunkXSize)+((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),((endX%chunkXSize)-(startX%chunkXSize))*bytes);
                            }
                            else if((cAV.x+1)*chunkXSize>endX){
                                memcpy((uint8_t*)zarr+((((cAV.x*chunkXSize)-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))*bytes),(uint8_t*)bufferDest+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(endX%chunkXSize)*bytes);
                            }
                            else if((cAV.x*chunkXSize) < startX && ((cAV.x+1)*chunkXSize) > startX){
                                memcpy((uint8_t*)zarr+((((cAV.x*chunkXSize-startX+(startX%chunkXSize)))+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))*bytes),(uint8_t*)bufferDest+(((startX%chunkXSize)+((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(chunkXSize-(startX%chunkXSize))*bytes);
                            }
                        }
                        else{
                            memcpy((uint8_t*)zarr+((((cAV.x*chunkXSize)-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))*bytes),(uint8_t*)bufferDest+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),chunkXSize*bytes);
                        }
                    }
                }
            }
            else if (order == 'C'){
                for(int64_t x = cAV.x*chunkXSize; x < (cAV.x+1)*chunkXSize; x++){
                    if(x>=endX) break;
                    else if(x<startX) continue;
                    for(int64_t y = cAV.y*chunkYSize; y < (cAV.y+1)*chunkYSize; y++){
                        if(y>=endY) break;
                        else if(y<startY) continue;
                        for(int64_t z = cAV.z*chunkZSize; z < (cAV.z+1)*chunkZSize; z++){
                            if(z>=endZ) break;
                            else if(z<startZ) continue;
                            switch(bytes){
                                case 1:
                                    ((uint8_t*)zarr)[((x-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))] = ((uint8_t*)bufferDest)[((z%chunkZSize)+((y%chunkYSize)*chunkZSize)+((x%chunkXSize)*chunkZSize*chunkYSize))];
                                    break;
                                case 2:
                                    ((uint16_t*)zarr)[((x-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))] = ((uint16_t*)bufferDest)[((z%chunkZSize)+((y%chunkYSize)*chunkZSize)+((x%chunkXSize)*chunkZSize*chunkYSize))];
                                    break;
                                case 4:
                                    ((float*)zarr)[((x-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))] = ((float*)bufferDest)[((z%chunkZSize)+((y%chunkYSize)*chunkZSize)+((x%chunkXSize)*chunkZSize*chunkYSize))];
                                    break;
                                case 8:
                                    ((double*)zarr)[((x-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))] = ((double*)bufferDest)[((z%chunkZSize)+((y%chunkYSize)*chunkZSize)+((x%chunkXSize)*chunkZSize*chunkYSize))];
                                    break;
                            }
                        }
                    }
                }
                
            }
            
        }
        free(bufferDest);
        free(buffer);
    }
#pragma omp parallel for
    for(int i = 0; i < cI.numChunks; i++){
        free(cI.chunkNames[i]);
    }
    free(cI.chunkNames);
    
    if(err){
        printf("zarr:threadError: %s\n",errString);
    }
}
uint64_t dTypeToBits(char* dtype){
    
    if(dtype[1] == 'u' && dtype[2] == '1'){
        return 8;
    }
    else if(dtype[1] == 'u' && dtype[2] == '2'){
        return 16;
    }
    else if(dtype[1] == 'f' && dtype[2] == '4'){
        return 32;
    }
    else if(dtype[1] == 'f' && dtype[2] == '8'){
        return 64;
    }
    else{
        return 0;
    }
    
}

void* parallelReadZarrWrapper(char* folderName,uint8_t crop, uint64_t startX, uint64_t startY, uint64_t startZ, uint64_t endX, uint64_t endY,uint64_t endZ){
    
    uint64_t shapeX = 0;
    uint64_t shapeY = 0;
    uint64_t shapeZ = 0;
    uint64_t chunkXSize = 0;
    uint64_t chunkYSize = 0;
    uint64_t chunkZSize = 0;
    char dtype[4];
    char order;
    char* cname = NULL;
    uint64_t clevel = 5;
    setValuesFromJSON(folderName,&chunkXSize,&chunkYSize,&chunkZSize,dtype,&order,&shapeX,&shapeY,&shapeZ,&cname,&clevel);
    
    if(!crop){
        startX = 0;
        startY = 0;
        startZ = 0;
        endX = shapeX;
        endY = shapeY;
        endZ = shapeZ;
    }
    else{
        startX--;
        startY--;
        startZ--;
    }
    /*
     * if(endX > shapeX || endY > shapeY || endZ > shapeZ){
     * printf("Upper bound is invalid\n");
     * return NULL;
     * }*/
    uint64_t dim[3];
    shapeX = endX-startX;
    shapeY = endY-startY;
    shapeZ = endZ-startZ;
    dim[0] = shapeX;
    dim[1] = shapeY;
    dim[2] = shapeZ;
    
    if(dtype[1] == 'u' && dtype[2] == '1'){
        uint64_t bits = 8;
        uint8_t* zarr = (uint8_t*)malloc(sizeof(uint8_t)*shapeX*shapeY*shapeZ);
        parallelReadZarr((void*)zarr,folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ,bits,order,cname);
        return (void*)zarr;
    }
    else if(dtype[1] == 'u' && dtype[2] == '2'){
        uint64_t bits = 16;
        uint16_t* zarr = (uint16_t*)malloc((uint64_t)(sizeof(uint16_t)*shapeX*shapeY*shapeZ));
        parallelReadZarr((void*)zarr,folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ,bits,order,cname);
        return (void*)zarr;
    }
    else if(dtype[1] == 'f' && dtype[2] == '4'){
        uint64_t bits = 32;
        float* zarr = (float*)malloc(sizeof(float)*shapeX*shapeY*shapeZ);
        parallelReadZarr((void*)zarr,folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ,bits,order,cname);
        return (void*)zarr;
    }
    else if(dtype[1] == 'f' && dtype[2] == '8'){
        uint64_t bits = 64;
        double* zarr = (double*)malloc(sizeof(double)*shapeX*shapeY*shapeZ);
        parallelReadZarr((void*)zarr,folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ,bits,order,cname);
        return (void*)zarr;
    }
    else{
        return NULL;
    }
}
