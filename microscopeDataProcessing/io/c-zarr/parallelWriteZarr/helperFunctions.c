#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <cjson/cJSON.h>
#include <omp.h>
#include <uuid/uuid.h>
#include "mex.h"
#include "helperFunctions.h"

const char fileSep =
#ifdef _WIN32
    '\\';
#else
    '/';
#endif

void* mallocDynamic(uint64_t x, uint64_t bits){
    switch(bits){
        case 8:
            return malloc(x*sizeof(uint8_t));
        case 16:
            return malloc(x*sizeof(uint16_t));
        case 32:
            return malloc(x*sizeof(float));
        case 64:
            return malloc(x*sizeof(double));
        default:
            printf("Image is not 8/16 bit, single, or double. Using single.");
            return malloc(x*sizeof(float));
    }
}

struct chunkAxisVals getChunkAxisVals(char* fileName){
    struct chunkAxisVals cAV;
    char* ptr;
    cAV.x = strtol(fileName, &ptr, 10);
    ptr++;
    cAV.y = strtol(ptr, &ptr, 10);
    ptr++;
    cAV.z = strtol(ptr, &ptr, 10);
    return cAV;
}

struct chunkInfo getChunkInfo(char* folderName, uint64_t startX, uint64_t startY, uint64_t startZ, uint64_t endX, uint64_t endY,uint64_t endZ,uint64_t chunkXSize,uint64_t chunkYSize,uint64_t chunkZSize){
    struct chunkInfo cI;
    cI.numChunks = 0;
    cI.chunkNames = NULL;

    uint64_t xStartAligned = startX-(startX%chunkXSize);
    uint64_t yStartAligned = startY-(startY%chunkYSize);
    uint64_t zStartAligned = startZ-(startZ%chunkZSize);
    uint64_t xStartChunk = (xStartAligned/chunkXSize);
    uint64_t yStartChunk = (yStartAligned/chunkYSize);
    uint64_t zStartChunk = (zStartAligned/chunkZSize);

    uint64_t xEndAligned = endX;
    uint64_t yEndAligned = endY;
    uint64_t zEndAligned = endZ;

    if(xEndAligned%chunkXSize) xEndAligned = endX-(endX%chunkXSize)+chunkXSize;
    if(yEndAligned%chunkYSize) yEndAligned = endY-(endY%chunkYSize)+chunkYSize;
    if(zEndAligned%chunkZSize) zEndAligned = endZ-(endZ%chunkZSize)+chunkZSize;
    uint64_t xEndChunk = (xEndAligned/chunkXSize);
    uint64_t yEndChunk = (yEndAligned/chunkYSize);
    uint64_t zEndChunk = (zEndAligned/chunkZSize);

    uint64_t xChunks = (xEndChunk-xStartChunk);
    uint64_t yChunks = (yEndChunk-yStartChunk);
    uint64_t zChunks = (zEndChunk-zStartChunk);

    uint64_t file_count = xChunks*yChunks*zChunks;

    char** chunkNames = malloc(file_count*sizeof(char*));
    #pragma omp parallel for collapse(3)
    for(uint64_t x = xStartChunk; x < xEndChunk; x++){
        for(uint64_t y = yStartChunk; y < yEndChunk; y++){
            for(uint64_t z = zStartChunk; z < zEndChunk; z++){
                uint64_t currFile = (z-zStartChunk)+((y-yStartChunk)*zChunks)+((x-xStartChunk)*yChunks*zChunks);
                asprintf(&chunkNames[currFile],"%llu.%llu.%llu",x,y,z);
            }
        }
    }
    cI.chunkNames = chunkNames;
    cI.numChunks = file_count;
    return cI;
}

void setChunkShapeFromJSON(cJSON *json, uint64_t *x, uint64_t *y, uint64_t *z){
    *x = json->child->valueint;
    *y = json->child->next->valueint;
    *z = json->child->next->next->valueint;
}

void setDTypeFromJSON(cJSON *json, char* dtype){
    dtype[0] = json->valuestring[0];
    dtype[1] = json->valuestring[1];
    dtype[2] = json->valuestring[2];
    dtype[3] = json->valuestring[3];
}

void setOrderFromJSON(cJSON *json, char* order){
    *order = json->valuestring[0];
}

void setShapeFromJSON(cJSON *json, uint64_t *x, uint64_t *y, uint64_t *z){
    *x = json->child->valueint;
    *y = json->child->next->valueint;
    *z = json->child->next->next->valueint;
}

void setValuesFromJSON(char* fileName,uint64_t *chunkXSize,uint64_t *chunkYSize,uint64_t *chunkZSize,char* dtype,char* order,uint64_t *shapeX,uint64_t *shapeY,uint64_t *shapeZ){

    char* zArray = ".zarray";
    char* fnFull = (char*)malloc(strlen(fileName)+9);
    fnFull[0] = '\0';
    char fileSepS[2];
    fileSepS[0] = fileSep;
    fileSepS[1] = '\0';

    strcat(fnFull,fileName);
    strcat(fnFull,fileSepS);
    strcat(fnFull,zArray);

    FILE *fileptr = fopen(fnFull, "rb");
    if(!fileptr) mexErrMsgIdAndTxt("zarr:inputError","Failed to open JSON File: %s\n",fnFull);
    free(fnFull);

    fseek(fileptr, 0, SEEK_END);
    long filelen = ftell(fileptr);
    rewind(fileptr);
    char* buffer = (char *)malloc((filelen));
    fread(buffer, filelen, 1, fileptr);
    fclose(fileptr);
    cJSON *json = cJSON_ParseWithLength(buffer,filelen);
    uint8_t flags[4] = {0,0,0,0};

    while(!(flags[0] && flags[1] && flags[2] && flags[3])){
        if(!json->string){
            json = json->child;
            continue;
        }
        else if(!strcmp(json->string,"chunks")){
            setChunkShapeFromJSON(json, chunkXSize,chunkYSize,chunkZSize);
            flags[0] = 1;
        }
        else if(!strcmp(json->string,"dtype")){
            setDTypeFromJSON(json, dtype);
            flags[1] = 1;
        }
        else if(!strcmp(json->string,"order")){
            setOrderFromJSON(json, order);
            flags[2] = 1;
        }
        else if(!strcmp(json->string,"shape")){
            setShapeFromJSON(json, shapeX,shapeY,shapeZ);
            flags[3] = 1;
        }
        json = json->next;
    }
    cJSON_Delete(json);
}

void setJSONValues(char* fileName,uint64_t *chunkXSize,uint64_t *chunkYSize,uint64_t *chunkZSize,char* dtype,char* order,uint64_t *shapeX,uint64_t *shapeY,uint64_t *shapeZ){
    
    // Overflows for ints greater than 32 bit for chunkSizes and shape
    cJSON* zArray = cJSON_CreateObject();
    const int chunkSizes[3] = {*chunkXSize,*chunkYSize,*chunkZSize};
    cJSON* chunks = cJSON_CreateIntArray(chunkSizes, 3);
    cJSON_AddItemToObject(zArray, "chunks", chunks);

    cJSON* compressor = cJSON_CreateObject();
    cJSON_AddItemToObject(zArray, "compressor", compressor);

    cJSON_AddNumberToObject(compressor, "blocksize", 0);
    cJSON_AddNumberToObject(compressor, "clevel", 5);
    cJSON_AddStringToObject(compressor, "cname", "lz4");
    cJSON_AddStringToObject(compressor, "id", "blosc");
    cJSON_AddNumberToObject(compressor, "shuffle", 1);

    cJSON_AddStringToObject(zArray, "dtype", dtype);
    cJSON_AddNumberToObject(zArray, "fill_value", 0);
    cJSON_AddNullToObject(zArray, "filters");
    cJSON_AddStringToObject(zArray, "order", "F");

    const int shapeSizes[3] = {*shapeX,*shapeY,*shapeZ};
    cJSON* shape = cJSON_CreateIntArray(shapeSizes, 3);
    cJSON_AddItemToObject(zArray, "shape", shape);
    cJSON_AddNumberToObject(zArray, "zarr_format", 2);
    
    uuid_t binuuid;
    uuid_generate_random(binuuid);
    char *uuid = malloc(37);
    uuid_unparse(binuuid, uuid);

    char* zArrayS = ".zarray";
    char* fnFull = (char*)malloc(strlen(fileName)+8+36+1);
    fnFull[0] = '\0';
    char fileSepS[2];
    fileSepS[0] = fileSep;
    fileSepS[1] = '\0';

    strcat(fnFull,fileName);
    strcat(fnFull,fileSepS);
    strcat(fnFull,zArrayS);
    char* fileNameFinal = strndup(fnFull,strlen(fileName)+1+8);
    strcat(fnFull,uuid);
    
    char* string = cJSON_Print(zArray);

    FILE *fileptr = fopen(fnFull, "w+");
    if(!fileptr) mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s\n",fnFull);
    fprintf(fileptr,"%s",string);
    fclose(fileptr);
    
    
    rename(fnFull,fileNameFinal);
    //file = fopen(fileNameFinal, "r");
    //if(!file) rename(fileName,fileNameFinal);
    
    cJSON_Delete(zArray);
    free(fnFull);
    free(uuid);
    free(fileNameFinal);
}
