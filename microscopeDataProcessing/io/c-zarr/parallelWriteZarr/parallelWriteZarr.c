#include "parallelWriteZarr.h"
#include "parallelReadZarr.h"
#include "helperFunctions.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <blosc.h>
#include <cjson/cJSON.h>
#include <omp.h>
#include <uuid/uuid.h>
#include <sys/stat.h>
#include "mex.h"

//mex -v COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/lib64' -lblosc2 zarrMex.c

void parallelWriteZarrMex(void* zarr, char* folderName,uint64_t startX, uint64_t startY, uint64_t startZ, uint64_t endX, uint64_t endY,uint64_t endZ,uint64_t chunkXSize,uint64_t chunkYSize,uint64_t chunkZSize,uint64_t shapeX,uint64_t shapeY,uint64_t shapeZ,uint64_t origShapeX,uint64_t origShapeY,uint64_t origShapeZ, uint64_t bits, char order, uint8_t useUuid, uint8_t crop){
    char fileSepS[2];
    const char fileSep =
#ifdef _WIN32
            '\\';
#else
    '/';
#endif
    fileSepS[0] = fileSep;
    fileSepS[1] = '\0';
    //printf("%s startxyz: %d %d %d endxyz: %d %d %d chunkxyz: %d %d %d shapexyz: %d %d %d bits: %d\n",folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ,bits);
    
    
    uint64_t bytes = (bits/8);
    
    /* Initialize the Blosc compressor */
    int32_t numWorkers = omp_get_max_threads();
    
    struct chunkInfo cI = getChunkInfo(folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize);
    //if(!cI.chunkNames) mexErrMsgIdAndTxt("zarr:inputError","File \"%s\" cannot be opened",folderName);
    char** chunkNamesUuid = malloc(cI.numChunks*sizeof(char*));
    
    int32_t batchSize = (cI.numChunks-1)/numWorkers+1;
    uint64_t s = chunkXSize*chunkYSize*chunkZSize;
    uint64_t sB = s*bytes;
    int32_t w;
    
    uint64_t xRest = 0;
    uint64_t yRest = 0;
    uint64_t zRest = 0;
    
    
    uuid_t binuuid;
    uuid_generate_random(binuuid);
    char *uuid = malloc(37);
    uuid_unparse(binuuid, uuid);
    #pragma omp parallel for //if(numWorkers<=cI.numChunks)
    for(w = 0; w < numWorkers; w++){
        void* chunkUnC = mallocDynamic(s,bits);
        void* chunkC = malloc(sB+BLOSC_MAX_OVERHEAD);
        for(int64_t f = w*batchSize; f < (w+1)*batchSize; f++){
            if(f>=cI.numChunks) break;
            struct chunkAxisVals cAV = getChunkAxisVals(cI.chunkNames[f]);
            void* cRegion = NULL;
            if(crop && ((((cAV.x)*chunkXSize) < startX || ((cAV.x+1)*chunkXSize > endX && endX < origShapeX))
            || (((cAV.y)*chunkYSize) < startY || ((cAV.y+1)*chunkYSize > endY && endY < origShapeY))
            || (((cAV.z)*chunkZSize) < startZ || ((cAV.z+1)*chunkZSize > endZ && endZ < origShapeZ)))){
                cRegion = parallelReadZarrWrapper(folderName, crop, ((cAV.x)*chunkXSize)+1, ((cAV.y)*chunkYSize)+1, ((cAV.z)*chunkZSize)+1, (cAV.x+1)*chunkXSize, (cAV.y+1)*chunkYSize, (cAV.z+1)*chunkZSize);
            }
            if(order == 'F'){
                for(int64_t z = cAV.z*chunkZSize; z < (cAV.z+1)*chunkZSize; z++){
                    if(z>=endZ){
                        if(crop){
                            if((cAV.z+1)*chunkZSize > origShapeZ){
                                memcpy((uint8_t*)chunkUnC+((((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)cRegion+((((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),((origShapeZ-z)*chunkXSize*chunkYSize)*bytes);
                                uint64_t zRest = ((cAV.z+1)*chunkZSize)-origShapeZ;
                                memset((uint8_t*)chunkUnC+(((z%chunkZSize)*chunkXSize*chunkYSize)*bytes),0,(zRest*(chunkXSize*chunkYSize))*bytes);
                            }
                            else{
                                memcpy((uint8_t*)chunkUnC+((((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)cRegion+((((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),((((cAV.z+1)*chunkZSize)-z)*chunkXSize*chunkYSize)*bytes);
                            }
                        }
                        else{
                            uint64_t zRest = ((cAV.z+1)*chunkZSize)-z;
                            memset((uint8_t*)chunkUnC+(((z%chunkZSize)*chunkXSize*chunkYSize)*bytes),0,(zRest*(chunkXSize*chunkYSize))*bytes);
                        }
                        break;
                    }
                    else if(z<startZ){
                        if(crop){
                            memcpy((uint8_t*)chunkUnC+(((z%chunkZSize)*chunkXSize*chunkYSize)*bytes),(uint8_t*)cRegion+(((z%chunkZSize)*chunkXSize*chunkYSize)*bytes),((startZ-z)*chunkXSize*chunkYSize)*bytes);
                        }
                        else{
                            memset((uint8_t*)chunkUnC+(((z%chunkZSize)*chunkXSize*chunkYSize)*bytes),0,((startZ-z)*(chunkXSize*chunkYSize))*bytes);
                        }
                        z = startZ-1;
                        continue;
                    }
                    for(int64_t y = cAV.y*chunkYSize; y < (cAV.y+1)*chunkYSize; y++){
                        if(y>=endY){
                            if(crop){
                                if((cAV.y+1)*chunkYSize > origShapeY){
                                    memcpy((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)cRegion+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),((origShapeY-y)*chunkXSize)*bytes);
                                    uint64_t yRest = ((cAV.y+1)*chunkYSize)-origShapeY;
                                    memset((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),0,(yRest*(chunkXSize))*bytes);
                                }
                                else{
                                    memcpy((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)cRegion+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),((((cAV.y+1)*chunkYSize)-y)*chunkXSize)*bytes);
                                }
                            }
                            else{
                                uint64_t yRest = ((cAV.y+1)*chunkYSize)-y;
                                memset((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),0,(yRest*chunkXSize)*bytes);
                            }
                            break;
                        }
                        else if(y<startY){
                            if(crop){
                                memcpy((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)cRegion+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),((startY-y)*chunkXSize)*bytes);
                            }
                            else{
                                memset((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),0,(startY-y)*bytes);
                            }
                            y = startY-1;
                            continue;
                        }
                        
                        if(((cAV.x*chunkXSize) < startX && ((cAV.x+1)*chunkXSize) > startX) || (cAV.x+1)*chunkXSize>endX){
                            if(((cAV.x*chunkXSize) < startX && ((cAV.x+1)*chunkXSize) > startX) && (cAV.x+1)*chunkXSize>endX){
                                if(crop){
                                    memcpy((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)cRegion+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(startX%chunkXSize)*bytes);
                                    memcpy((uint8_t*)chunkUnC+(((startX%chunkXSize)+((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)zarr+((((cAV.x*chunkXSize)-startX+(startX%chunkXSize))+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))*bytes),((endX%chunkXSize)-(startX%chunkXSize))*bytes);
                                    memcpy((uint8_t*)chunkUnC+(((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))+(endX%chunkXSize))*bytes),(uint8_t*)cRegion+(((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))+(endX%chunkXSize))*bytes),(chunkXSize-(endX%chunkXSize))*bytes);
                                }
                                else{
                                    memset((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),0,(startX%chunkXSize)*bytes);
                                    memcpy((uint8_t*)chunkUnC+(((startX%chunkXSize)+((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)zarr+((((cAV.x*chunkXSize)-startX+(startX%chunkXSize))+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))*bytes),((endX%chunkXSize)-(startX%chunkXSize))*bytes);
                                    memset((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))+(endX%chunkXSize)*bytes),0,(chunkXSize-(endX%chunkXSize))*bytes);
                                }
                            }
                            else if((cAV.x+1)*chunkXSize>endX){
                                if(crop){
                                    memcpy((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)zarr+((((cAV.x*chunkXSize)-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))*bytes),(endX-(cAV.x*chunkXSize))*bytes);
                                    
                                    if((cAV.x+1)*chunkXSize > origShapeX){
                                        memcpy((uint8_t*)chunkUnC+((((endX-(cAV.x*chunkXSize)))+((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)cRegion+((((endX-(cAV.x*chunkXSize)))+((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(origShapeX-endX)*bytes);
                                        uint64_t xRest = ((cAV.x+1)*chunkXSize)-origShapeX;
                                        memset((uint8_t*)chunkUnC+(((origShapeX-(cAV.x*chunkXSize))+((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),0,(xRest)*bytes);
                                    }
                                    else{
                                        memcpy((uint8_t*)chunkUnC+((((endX-(cAV.x*chunkXSize)))+((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)cRegion+((((endX-(cAV.x*chunkXSize)))+((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(((cAV.x+1)*chunkXSize)-endX)*bytes);
                                    }
                                }
                                else{
                                    memcpy((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)zarr+((((cAV.x*chunkXSize)-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))*bytes),(endX%chunkXSize)*bytes);
                                    memset((uint8_t*)chunkUnC+(((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))+(endX%chunkXSize))*bytes),0,(chunkXSize-(endX%chunkXSize))*bytes);
                                }
                            }
                            else if((cAV.x*chunkXSize) < startX && ((cAV.x+1)*chunkXSize) > startX){
                                if(crop){
                                    memcpy((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)cRegion+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(startX%chunkXSize)*bytes);
                                    memcpy((uint8_t*)chunkUnC+(((startX%chunkXSize)+((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)zarr+((((cAV.x*chunkXSize)-startX+(startX%chunkXSize))+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))*bytes),(chunkXSize-(startX%chunkXSize))*bytes);
                                }
                                else{
                                    memset((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),0,(startX%chunkXSize)*bytes);
                                    memcpy((uint8_t*)chunkUnC+(((startX%chunkXSize)+((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)zarr+((((cAV.x*chunkXSize)-startX+(startX%chunkXSize))+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))*bytes),(chunkXSize-(startX%chunkXSize))*bytes);
                                }
                            }
                        }
                        else{
                            memcpy((uint8_t*)chunkUnC+((((y%chunkYSize)*chunkXSize)+((z%chunkZSize)*chunkXSize*chunkYSize))*bytes),(uint8_t*)zarr+((((cAV.x*chunkXSize)-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY))*bytes),chunkXSize*bytes);
                        }
                    }
                }
            }
            else if (order == 'C'){
                for(int64_t x = cAV.x*chunkZSize; x < (cAV.x+1)*chunkXSize; x++){
                    for(int64_t y = cAV.y*chunkYSize; y < (cAV.y+1)*chunkYSize; y++){
                        for(int64_t z = cAV.z*chunkZSize; z < (cAV.z+1)*chunkZSize; z++){
                            if(x>=endX || x<startX || y>= endY || y<startY || z>=endZ || z<startZ){
                                ((uint8_t*)chunkUnC)[(((z%chunkZSize)+((y%chunkYSize)*chunkZSize)+((x%chunkXSize)*chunkZSize*chunkYSize))*bytes)] = 0;
                                continue;
                            }
                            ((uint8_t*)chunkUnC)[(((z%chunkZSize)+((y%chunkYSize)*chunkZSize)+((x%chunkXSize)*chunkZSize*chunkYSize))*bytes)] = ((uint8_t*)zarr)[((x+(y*shapeX)+(z*shapeX*shapeY))*bytes)];
                        }
                    }
                }
            }
            /* Compress with clevel=5 and shuffle active  */
            //char* compressor = blosc_get_compressor();
            //printf("Thread: %d Compressor: %s\n",w,compressor);
            
            // Use the same blosc compress as Zarr
            //int64_t csize = blosc_compress(5, BLOSC_SHUFFLE, bytes, sB, chunkUnC, chunkC, sB+BLOSC_MAX_OVERHEAD);
            int64_t csize = blosc_compress_ctx(5, BLOSC_SHUFFLE, bytes, sB, chunkUnC, chunkC, sB+BLOSC_MAX_OVERHEAD,"lz4",0,1);
            
            //malloc +2 for null term and filesep
            char *fileName = malloc(strlen(folderName)+1+strlen(cI.chunkNames[f])+36+1);
            fileName[0] = '\0';
            strcat(fileName,folderName);
            strcat(fileName,fileSepS);
            strcat(fileName,cI.chunkNames[f]);
            
            //FILE *fileptr = fopen(fileName, "r+b");
            
            if(useUuid){
                strcat(fileName,uuid);
                char* fileNameFinal = strndup(fileName,strlen(folderName)+1+strlen(cI.chunkNames[f]));
                FILE *fileptr = fopen(fileName, "w+b");
                fwrite(chunkC,csize,1,fileptr);
                fclose(fileptr);
                rename(fileName,fileNameFinal);
                free(fileNameFinal);
            }
            else{
                
                FILE *fileptr = fopen(fileName, "w+b");
                fwrite(chunkC,csize,1,fileptr);
                fclose(fileptr);
            }
            free(fileName);
            free(cRegion);
        }
        free(chunkUnC);
        free(chunkC);
        
    }
    
    free(uuid);
#pragma omp parallel for
    for(int i = 0; i < cI.numChunks; i++){
        free(cI.chunkNames[i]);
    }
    free(cI.chunkNames);
    
    /* After using it, destroy the Blosc environment */
    //blosc_destroy();
}

// TODO: FIX MEMORY LEAKS
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    uint64_t startX = 0;
    uint64_t startY = 0;
    uint64_t startZ = 0;
    uint64_t endX = 0;
    uint64_t endY = 0;
    uint64_t endZ = 0;
    uint8_t crop = 0;
    if(nrhs < 3) mexErrMsgIdAndTxt("zarr:inputError","This functions requires at least 3 arguments");
    else if(nrhs == 4){
        if(mxGetN(prhs[3]) == 6){
            
            crop = 1;
            startX = (uint64_t)*(mxGetPr(prhs[3]))-1;
            startY = (uint64_t)*((mxGetPr(prhs[3])+1))-1;
            startZ = (uint64_t)*((mxGetPr(prhs[3])+2))-1;
            endX = (uint64_t)*((mxGetPr(prhs[3])+3));
            endY = (uint64_t)*((mxGetPr(prhs[3])+4));
            endZ = (uint64_t)*((mxGetPr(prhs[3])+5));
            
            if(startX+1 < 1 || startY+1 < 1 || startZ+1 < 1) mexErrMsgIdAndTxt("zarr:inputError","Lower bounds must be at least 1");
        }
        else if(mxGetN(prhs[3]) != 3) mexErrMsgIdAndTxt("zarr:inputError","Input range is not 6 or 3");
    }
    else if (nrhs > 4) mexErrMsgIdAndTxt("zarr:inputError","Number of input arguments must be 4 or less");
    if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("zarr:inputError","The first argument must be a string");
    char* folderName = mxArrayToString(prhs[0]);
    uint8_t useUuid = (uint8_t)*(mxGetPr(prhs[2]));
    uint64_t shapeX = 0;
    uint64_t shapeY = 0;
    uint64_t shapeZ = 0;
    uint64_t chunkXSize = 0;
    uint64_t chunkYSize = 0;
    uint64_t chunkZSize = 0;
    char dtype[4];
    char order;
    void* zarrC = NULL;
    
    
    mxClassID mDType = mxGetClassID(prhs[1]);
    dtype[0] = '<';
    if(mDType == mxUINT8_CLASS){
        dtype[1] = 'u';
        dtype[2] = '1';
        
    }
    else if(mDType == mxUINT16_CLASS){
        dtype[1] = 'u';
        dtype[2] = '2';
        
    }
    else if(mDType == mxSINGLE_CLASS){
        dtype[1] = 'f';
        dtype[2] = '4';
        
    }
    else if(mDType == mxDOUBLE_CLASS){
        dtype[1] = 'f';
        dtype[2] = '8';
        
    }
    dtype[3] = '\0';
    chunkXSize = 256;
    chunkYSize = 256;
    chunkZSize = 256;
    order = 'F';
    
    char* zArray = ".zarray";
    char* fnFull = (char*)malloc(strlen(folderName)+9);
    fnFull[0] = '\0';
    char fileSepS[2];
    fileSepS[0] = '/';
    fileSepS[1] = '\0';
    
    strcat(fnFull,folderName);
    strcat(fnFull,fileSepS);
    strcat(fnFull,zArray);
    
    if(!crop){
        uint64_t* dims = (uint64_t*)mxGetDimensions(prhs[1]);
        shapeX = dims[0];
        shapeY = dims[1];
        shapeZ = dims[2];
        chunkXSize = (uint64_t)*(mxGetPr(prhs[3]));
        chunkYSize = (uint64_t)*((mxGetPr(prhs[3])+1));
        chunkZSize = (uint64_t)*((mxGetPr(prhs[3])+2));

        FILE* f = fopen(fnFull,"r");
        if(f) fclose(f);
        else{
            mkdir(folderName, 0775);
            chmod(folderName, 0775);
        }
        
        
        setJSONValues(folderName,&chunkXSize,&chunkYSize,&chunkZSize,dtype,&order,&shapeX,&shapeY, &shapeZ);
        setValuesFromJSON(folderName,&chunkXSize,&chunkYSize,&chunkZSize,dtype,&order,&shapeX,&shapeY,&shapeZ);
        //}
        
    }
    else{
        shapeX = endX;
        shapeY = endY;
        shapeZ = endZ;
        
        FILE* f = fopen(fnFull,"r");
        if(f) fclose(f);
        else {
            mkdir(folderName, 0775);
            chmod(folderName, 0775);
            setJSONValues(folderName,&chunkXSize,&chunkYSize,&chunkZSize,dtype,&order,&shapeX,&shapeY, &shapeZ);
        }
        
        char dtypeT[4];
        for(int i = 0; i < 4; i++) dtypeT[i] = dtype[i];
        
        setValuesFromJSON(folderName,&chunkXSize,&chunkYSize,&chunkZSize,dtype,&order,&shapeX,&shapeY,&shapeZ);
        
        if(dtypeT[2] != dtype[2]){
            uint64_t size = (endX-startX)*(endY-startY)*(endZ-startZ);
            
            uint64_t bitsT = 0;
            if(dtypeT[1] == 'u' && dtypeT[2] == '1') bitsT = 8;
            else if(dtypeT[1] == 'u' && dtypeT[2] == '2') bitsT = 16;
            else if(dtypeT[1] == 'f' && dtypeT[2] == '4') bitsT = 32;
            else if(dtypeT[1] == 'f' && dtypeT[2] == '8') bitsT = 64;
            else mexErrMsgIdAndTxt("tiff:dataTypeError","Cannont convert to passed in data type. Data type not suppported");
            
            
            if(dtype[1] == 'u' && dtype[2] == '1'){
                zarrC = malloc(size*sizeof(uint8_t));
                if(bitsT == 16){
                    uint16_t* zarrT = (uint16_t*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((uint8_t*)zarrC)[i] = (uint8_t)zarrT[i];
                    }
                }
                else if(bitsT == 32){
                    float* zarrT = (float*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((uint8_t*)zarrC)[i] = (uint8_t)zarrT[i];
                    }
                }
                else if(bitsT == 64){
                    double* zarrT = (double*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((uint8_t*)zarrC)[i] = (uint8_t)zarrT[i];
                    }
                }
            }
            else if(dtype[1] == 'u' && dtype[2] == '2'){
                zarrC = malloc(size*sizeof(uint16_t));
                if(bitsT == 8){
                    uint8_t* zarrT = (uint8_t*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((uint16_t*)zarrC)[i] = (uint16_t)zarrT[i];
                    }
                }
                else if (bitsT == 32){
                    float* zarrT = (float*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((uint16_t*)zarrC)[i] = (uint16_t)zarrT[i];
                    }
                }
                else if (bitsT == 64){
                    double* zarrT = (double*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((uint16_t*)zarrC)[i] = (uint16_t)zarrT[i];
                    }
                }
            }
            else if(dtype[1] == 'f' && dtype[2] == '4'){
                zarrC = malloc(size*sizeof(float));
            	if(bitsT == 8){
                    uint8_t* zarrT = (uint8_t*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((float*)zarrC)[i] = (float)zarrT[i];
                    }
                }
                else if(bitsT == 16){
                    uint16_t* zarrT = (uint16_t*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((float*)zarrC)[i] = (float)zarrT[i];
                    }
                }
                else if(bitsT == 64){
                    double* zarrT = (double*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((float*)zarrC)[i] = (float)zarrT[i];
                    }
                }
            }
            else if(dtype[1] == 'f' && dtype[2] == '8'){
                zarrC = malloc(size*sizeof(double));
                if(bitsT == 8){
                    uint8_t* zarrT = (uint8_t*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((double*)zarrC)[i] = (double)zarrT[i];
                    }
                }
                else if(bitsT == 16){
                    uint16_t* zarrT = (uint16_t*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((double*)zarrC)[i] = (double)zarrT[i];
                    }
                }
                else if(bitsT == 32){
                    float* zarrT = (float*)mxGetPr(prhs[1]);
                    #pragma omp parallel for
                    for(uint64_t i = 0; i < size; i++){
                        ((double*)zarrC)[i] = (double)zarrT[i];
                    }
                }
            }
            else{
                mexErrMsgIdAndTxt("tiff:dataTypeError","Cannont convert to passed in data type. Data type not suppported");
            }
        }
    }
    
    free(fnFull);
    
    uint64_t origShapeX = shapeX;
    uint64_t origShapeY = shapeY;
    uint64_t origShapeZ = shapeZ;
    if(endX > shapeX || endY > shapeY || endZ > shapeZ) mexErrMsgIdAndTxt("zarr:inputError","Upper bound is invalid");
    if(!crop){
        endX = shapeX;
        endY = shapeY;
        endZ = shapeZ;
        startX = 0;
        startY = 0;
        startZ = 0;
    }
    uint64_t dim[3];
    shapeX = endX-startX;
    shapeY = endY-startY;
    shapeZ = endZ-startZ;
    dim[0] = shapeX;
    dim[1] = shapeY;
    dim[2] = shapeZ;
    
    if(dtype[1] == 'u' && dtype[2] == '1'){
        uint64_t bits = 8;
        uint8_t* zarr;
        if(zarrC) zarr = (uint8_t*)zarrC;
        else zarr =  (uint8_t*)mxGetPr(prhs[1]);
        parallelWriteZarrMex((void*)zarr,folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ, origShapeX, origShapeY,origShapeZ, bits,order,useUuid,crop);
    }
    else if(dtype[1] == 'u' && dtype[2] == '2'){
        uint64_t bits = 16;
        uint16_t* zarr;
        if(zarrC) zarr = (uint16_t*)zarrC;
        else zarr = (uint16_t*)mxGetPr(prhs[1]);
        parallelWriteZarrMex((void*)zarr,folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ, origShapeX, origShapeY,origShapeZ, bits,order,useUuid,crop);
    }
    else if(dtype[1] == 'f' && dtype[2] == '4'){
        uint64_t bits = 32;
        float* zarr = (float*)mxGetPr(prhs[1]);
        parallelWriteZarrMex((void*)zarr,folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ, origShapeX, origShapeY,origShapeZ, bits,order,useUuid,crop);
    }
    else if(dtype[1] == 'f' && dtype[2] == '8'){
        uint64_t bits = 64;
        double* zarr = (double*)mxGetPr(prhs[1]);
        parallelWriteZarrMex((void*)zarr,folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ, origShapeX, origShapeY,origShapeZ, bits,order,useUuid,crop);
    }
    else{
        mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
    }
    
    
}

