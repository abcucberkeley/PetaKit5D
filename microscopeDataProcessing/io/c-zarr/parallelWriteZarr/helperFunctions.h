#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H
#include <stdint.h>
#include <cjson/cJSON.h>

#ifdef _WIN32
char* strndup (const char *s, size_t n);

#endif  

void* mallocDynamic(uint64_t x, uint64_t bits);

struct chunkInfo{
    char** chunkNames;
    int64_t numChunks;
};

struct chunkAxisVals{
    uint64_t x;
    uint64_t y;
    uint64_t z;
};

struct chunkAxisVals getChunkAxisVals(char* fileName);

struct chunkInfo getChunkInfo(char* folderName, uint64_t startX, uint64_t startY, uint64_t startZ, uint64_t endX, uint64_t endY,uint64_t endZ,uint64_t chunkXSize,uint64_t chunkYSize,uint64_t chunkZSize);

void setChunkShapeFromJSON(cJSON *json, uint64_t *x, uint64_t *y, uint64_t *z);

void setDTypeFromJSON(cJSON *json, char* dtype);

void setOrderFromJSON(cJSON *json, char* order);

void setShapeFromJSON(cJSON *json, uint64_t *x, uint64_t *y, uint64_t *z);

void setValuesFromJSON(char* fileName,uint64_t *chunkXSize,uint64_t *chunkYSize,uint64_t *chunkZSize,char* dtype,char* order,uint64_t *shapeX,uint64_t *shapeY,uint64_t *shapeZ,char** cname);

void setJSONValues(char* fileName,uint64_t *chunkXSize,uint64_t *chunkYSize,uint64_t *chunkZSize,char* dtype,char* order,uint64_t *shapeX,uint64_t *shapeY,uint64_t *shapeZ, char* cname);

#endif // HELPERFUNCTIONS_H
