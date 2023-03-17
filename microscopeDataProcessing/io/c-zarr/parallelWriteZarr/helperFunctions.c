#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <cjson/cJSON.h>
#include <omp.h>
#ifdef _WIN32
#include <stdarg.h> 
#include <sys/time.h>
#else
#include <uuid/uuid.h>
#endif
#include <sys/stat.h>
#include "mex.h"
#include "helperFunctions.h"

#ifndef _WIN32
#include <wordexp.h>
char* expandTilde(char* path) {
    wordexp_t expPath;
    wordexp(path, &expPath, 0);
    return expPath.we_wordv[0];
}
#endif

const char fileSep =
#ifdef _WIN32
    '\\';
#else
    '/';
#endif

#ifdef _WIN32
char* strndup (const char *s, size_t n)
{
  size_t len = strnlen (s, n);
  char *new = (char *) malloc (len + 1);
  if (new == NULL)
    return NULL;
  new[len] = '\0';
  return (char *) memcpy (new, s, len);
}

int _vscprintf_so(const char * format, va_list pargs) {
    int retval;
    va_list argcopy;
    va_copy(argcopy, pargs);
    retval = vsnprintf(NULL, 0, format, argcopy);
    va_end(argcopy);
    return retval;
}

int vasprintf(char **strp, const char *fmt, va_list ap) {
    int len = _vscprintf_so(fmt, ap);
    if (len == -1) return -1;
    char *str = malloc((size_t) len + 1);
    if (!str) return -1;
    int r = vsnprintf(str, len + 1, fmt, ap); /* "secure" version of vsprintf */
    if (r == -1) return free(str), -1;
    *strp = str;
    return r;
}

int asprintf(char *strp[], const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    int r = vasprintf(strp, fmt, ap);
    va_end(ap);
    return r;
}
#endif

static void mkdirRecursive(const char *dir) {
    char tmp[8192];
    char *p = NULL;
    size_t len;
    #ifdef _WIN32
    char fileSep = '\\';
    #else
    char fileSep = '/';
    #endif
    int status;
    snprintf(tmp, sizeof(tmp),"%s",dir);
    len = strlen(tmp);
    if (tmp[len - 1] == fileSep)
        tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++){
        if (*p == fileSep) {
            *p = 0;

            #ifdef _WIN32
            mkdir(tmp);
            #else
            mkdir(tmp, 0775);
            #endif

            chmod(tmp, 0775);
            *p = fileSep;
        }
    }
    #ifdef _WIN32
    mkdir(tmp);
    #else
    mkdir(tmp, 0775);
    #endif
    chmod(tmp, 0775);
}

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

char* getSubfolderString(struct chunkAxisVals *cAV, uint64_t subfolderSizeX, uint64_t subfolderSizeY, uint64_t subfolderSizeZ){
    if(subfolderSizeX == 0 && subfolderSizeY == 0 && subfolderSizeZ == 0) return NULL;
    
    uint64_t currX = 0;
    uint64_t currY = 0;
    uint64_t currZ = 0;
    if(subfolderSizeX > 0) currX = cAV->x/subfolderSizeX;
    if(subfolderSizeY > 0) currY = cAV->y/subfolderSizeY;
    if(subfolderSizeZ > 0) currZ = cAV->z/subfolderSizeZ;

    char* currName = NULL;
    asprintf(&currName,"%llu_%llu_%llu",currX,currY,currZ);
    return currName;

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

void setCnameFromJSON(cJSON *json, char** cname){
    cJSON *jsonItem = json->child;
    

    while(jsonItem){
        if(!strcmp(jsonItem->string,"cname")){
            *cname = strdup(jsonItem->valuestring);
            return;
        }
        // For gzip
        if(!strcmp(jsonItem->string,"id") && !strcmp(jsonItem->valuestring,"gzip")){     
            *cname = strdup(jsonItem->valuestring);
            //mexErrMsgIdAndTxt("zarr:zarrayError","Compressor: \"%s\" is not currently supportedEXTRAIN\n",cname);
            return;
        }
        jsonItem = jsonItem->next;
    } 
    mexErrMsgIdAndTxt("zarr:zarrayError","Compressor: \"%s\" is not currently supported\n",*cname);
    
}

void setClevelFromJSON(cJSON *json, uint64_t* clevel){
    cJSON *jsonItem = json->child;
    
    while(jsonItem){
        if(!strcmp(jsonItem->string,"clevel")){
            *clevel = jsonItem->valueint;
            return;
        }
        // For gzip
        if(!strcmp(jsonItem->string,"level")){     
            *clevel = jsonItem->valueint;
            return;
        }
        jsonItem = jsonItem->next;
    } 
    mexErrMsgIdAndTxt("zarr:zarrayError","Compression level not found in .zarray file\n");
    
}

void setSubfoldersFromJSON(cJSON *json, uint64_t *subfolderSizeX, uint64_t *subfolderSizeY, uint64_t *subfolderSizeZ){
    *subfolderSizeX = json->child->valueint;
    *subfolderSizeY = json->child->next->valueint;
    *subfolderSizeZ = json->child->next->next->valueint;

}

void setValuesFromJSON(char* fileName,uint64_t *chunkXSize,uint64_t *chunkYSize,uint64_t *chunkZSize,char* dtype,char* order,uint64_t *shapeX,uint64_t *shapeY,uint64_t *shapeZ,char** cname, uint64_t* clevel, uint64_t *subfolderSizeX, uint64_t *subfolderSizeY, uint64_t *subfolderSizeZ){

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

    while(json){
        if(!json->string){
            json = json->child;
            continue;
        }
        else if(!strcmp(json->string,"chunks")){
            setChunkShapeFromJSON(json, chunkXSize,chunkYSize,chunkZSize);
        }
        else if(!strcmp(json->string,"dtype")){
            setDTypeFromJSON(json, dtype);
        }
        else if(!strcmp(json->string,"order")){
            setOrderFromJSON(json, order);
        }
        else if(!strcmp(json->string,"shape")){
            setShapeFromJSON(json, shapeX,shapeY,shapeZ);
        }
        else if(!strcmp(json->string,"compressor")){
            setCnameFromJSON(json, cname);
            setClevelFromJSON(json, clevel);
        }
        else if(!strcmp(json->string,"subfolders")){
            setSubfoldersFromJSON(json, subfolderSizeX, subfolderSizeY, subfolderSizeZ);
        }
        json = json->next;
    }
    cJSON_Delete(json);
}

void setJSONValues(char* fileName,uint64_t *chunkXSize,uint64_t *chunkYSize,uint64_t *chunkZSize,char* dtype,char* order,uint64_t *shapeX,uint64_t *shapeY,uint64_t *shapeZ, char* cname, uint64_t* clevel, uint64_t *subfolderSizeX, uint64_t *subfolderSizeY, uint64_t *subfolderSizeZ){
    
    // Overflows for ints greater than 32 bit for chunkSizes and shape
    cJSON* zArray = cJSON_CreateObject();
    const int chunkSizes[3] = {*chunkXSize,*chunkYSize,*chunkZSize};
    cJSON* chunks = cJSON_CreateIntArray(chunkSizes, 3);
    cJSON_AddItemToObject(zArray, "chunks", chunks);

    cJSON* compressor = cJSON_CreateObject();
    cJSON_AddItemToObject(zArray, "compressor", compressor);
    
    if(!strcmp(cname,"lz4") || !strcmp(cname,"blosclz") || !strcmp(cname,"lz4hc") || !strcmp(cname,"zlib") || !strcmp(cname,"zstd")){
        cJSON_AddNumberToObject(compressor, "blocksize", 0);
        cJSON_AddNumberToObject(compressor, "clevel", *clevel);
        cJSON_AddStringToObject(compressor, "cname", cname);
        cJSON_AddStringToObject(compressor, "id", "blosc");
        cJSON_AddNumberToObject(compressor, "shuffle", 1);
    }
    else if(!strcmp(cname,"gzip")){
        cJSON_AddStringToObject(compressor, "id", cname);
        cJSON_AddNumberToObject(compressor, "level", *clevel);
    }
    else mexErrMsgIdAndTxt("zarr:zarrayError","Compressor: \"%s\" is not currently supported\n",cname);

    cJSON_AddStringToObject(zArray, "dtype", dtype);
    cJSON_AddNumberToObject(zArray, "fill_value", 0);
    cJSON_AddNullToObject(zArray, "filters");
    cJSON_AddStringToObject(zArray, "order", "F");

    const int shapeSizes[3] = {*shapeX,*shapeY,*shapeZ};
    cJSON* shape = cJSON_CreateIntArray(shapeSizes, 3);
    cJSON_AddItemToObject(zArray, "shape", shape);
    cJSON_AddNumberToObject(zArray, "zarr_format", 2);

    const int subfolderSizeSizes[3] = {*subfolderSizeX,*subfolderSizeY,*subfolderSizeZ};
    cJSON* subfolderSize = cJSON_CreateIntArray(subfolderSizeSizes, 3);
    cJSON_AddItemToObject(zArray, "subfolders", subfolderSize);
    
    uint64_t uuidLen;
    #ifdef _WIN32
    uuidLen = 5;
    char *uuid = malloc(uuidLen+1);
    char *seedArr = malloc(1000);
    struct timeval cSeed;
    gettimeofday(&cSeed,NULL);
    int nChars = sprintf(seedArr,"%d%d",cSeed.tv_sec,cSeed.tv_usec);
    int aSeed = 0;
    char* ptr;
    if(nChars > 9)
        aSeed = strtol(seedArr+nChars-10, &ptr, 9);
    else aSeed = strtol(seedArr, &ptr, 9);
    srand(aSeed);
    sprintf(uuid,"%.5d",rand() % 99999);
    free(seedArr);
    #else
    uuidLen = 36;
    uuid_t binuuid;
    uuid_generate_random(binuuid);
    char *uuid = malloc(uuidLen+1);
    uuid_unparse(binuuid, uuid);
    #endif

    char* zArrayS = ".zarray";
    char* fnFull = (char*)malloc(strlen(fileName)+8+uuidLen+1);
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

uint64_t fastCeilDiv(uint64_t num, uint64_t denom){
    return 1 + ((num - 1) / denom);
}

void createSubfolders(char* folderName, uint64_t shapeX, uint64_t shapeY, uint64_t shapeZ, uint64_t chunkXSize, uint64_t chunkYSize, uint64_t chunkZSize, uint64_t subfolderSizeX, uint64_t subfolderSizeY, uint64_t subfolderSizeZ){
    if(subfolderSizeX == 0 && subfolderSizeY == 0 && subfolderSizeZ == 0) return;
    uint64_t nChunksX = fastCeilDiv(shapeX,chunkXSize);
    uint64_t nChunksY = fastCeilDiv(shapeY,chunkYSize);
    uint64_t nChunksZ = fastCeilDiv(shapeZ,chunkZSize);

    

    
    uint64_t nSubfoldersX = 1;
    uint64_t nSubfoldersY = 1;
    uint64_t nSubfoldersZ = 1;

    if(subfolderSizeX > 0) nSubfoldersX = fastCeilDiv(nChunksX,subfolderSizeX);
    if(subfolderSizeY > 0) nSubfoldersY = fastCeilDiv(nChunksY,subfolderSizeY);
    if(subfolderSizeZ > 0) nSubfoldersZ = fastCeilDiv(nChunksZ,subfolderSizeZ);

    //printf("%d,%d,%d\n",nSubfoldersX,nSubfoldersY,nSubfoldersZ);
    //printf("%d,%d,%d\n",nChunksX,nChunksY,nChunksZ);
    //printf("%d,%d,%d\n",subfolderSizeX,subfolderSizeY,subfolderSizeZ);
    //mexErrMsgIdAndTxt("zarr:zarrayError","bork");
    #pragma omp parallel for collapse(3)
    for(uint64_t x = 0; x < nSubfoldersX; x++){
        for(uint64_t y = 0; y < nSubfoldersY; y++){
            for(uint64_t z = 0; z < nSubfoldersZ; z++){
                char* currName = NULL;
                asprintf(&currName,"%s/%llu_%llu_%llu",folderName,x,y,z);
                mkdirRecursive(currName);
                free(currName);
            }
        }
    }
}
