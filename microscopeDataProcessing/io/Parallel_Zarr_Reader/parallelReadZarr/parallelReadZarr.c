#include <stdio.h>
#include <stdint.h>
#include <dirent.h>
#include "cBlosc2/include/blosc2.h"
#include "cJSON/include/cjson/cJSON.h"
#include <omp.h>
#include <stdlib.h>
#include "mex.h"

//mex -v COPTIMFLAGS="-O3 -fwrapv -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/include/' '-I/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/include/' '-L/global/home/groups/software/sl-7.x86_64/modules/cBlosc/2.0.4/lib64' -lblosc2 -L/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/lib64' -lcjson zarrMex.c


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

struct chunkInfo{
    char** chunkNames;
    int64_t numChunks;
};

struct chunkAxisVals{
    uint64_t x;
    uint64_t y;
    uint64_t z;
};

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
    int file_count = 0;
    DIR * dirp;
    struct dirent * entry;
    struct chunkInfo cI;
    cI.numChunks = 0;
    cI.chunkNames = NULL;

    dirp = opendir(folderName);
    if(!dirp){
        printf("Failed to open dir\n");
        return cI;
    }

    while ((entry = readdir(dirp)) != NULL) {
        if (entry->d_name[0] != '.') { /* If the entry is a regular file */
            struct chunkAxisVals cAV = getChunkAxisVals(entry->d_name);
            if((cAV.x+1)*chunkXSize < startX || (cAV.y+1)*chunkYSize < startY || (cAV.z+1)*chunkZSize < startZ) continue;
            if((cAV.x)*chunkXSize >= endX || (cAV.y)*chunkYSize >= endY || (cAV.z)*chunkZSize >= endZ) continue;
            file_count++;
        }
    }
    rewinddir(dirp);
    char** chunkNames = malloc(file_count*sizeof(char*));
    int currDir = 0;
    while ((entry = readdir(dirp)) != NULL) {
        if (entry->d_name[0] != '.') { /* If the entry is a regular file */
            struct chunkAxisVals cAV = getChunkAxisVals(entry->d_name);
            if((cAV.x+1)*chunkXSize < startX || (cAV.y+1)*chunkYSize < startY || (cAV.z+1)*chunkZSize < startZ) continue;
            if((cAV.x)*chunkXSize >= endX || (cAV.y)*chunkYSize >= endY || (cAV.z)*chunkZSize >= endZ) continue;
            chunkNames[currDir] = malloc(strlen(entry->d_name)+1);
            strcpy(chunkNames[currDir],entry->d_name);
            currDir++;
        }
    }
    
    closedir(dirp);
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

void setShapeFromJSON(cJSON *json, uint64_t *x, uint64_t *y, uint64_t *z){
    *x = json->child->valueint;
    *y = json->child->next->valueint;
    *z = json->child->next->next->valueint;
}

void setValuesFromJSON(char* fileName,uint64_t *chunkXSize,uint64_t *chunkYSize,uint64_t *chunkZSize,char* dtype,uint64_t *shapeX,uint64_t *shapeY,uint64_t *shapeZ){

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
    uint8_t flags[3] = {0,0,0};

    while(!(flags[0] && flags[1] && flags[2])){
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
        else if(!strcmp(json->string,"shape")){
            setShapeFromJSON(json, shapeX,shapeY,shapeZ);
            flags[2] = 1;
        }
        json = json->next;
    }
    cJSON_Delete(json);
}

void parallelReadZarrMex(void* zarr, char* folderName,uint64_t startX, uint64_t startY, uint64_t startZ, uint64_t endX, uint64_t endY,uint64_t endZ,uint64_t chunkXSize,uint64_t chunkYSize,uint64_t chunkZSize,uint64_t shapeX,uint64_t shapeY,uint64_t shapeZ, uint64_t bits){
    char fileSepS[2];
    fileSepS[0] = fileSep;
    fileSepS[1] = '\0';

    /* Initialize the Blosc compressor */
    int32_t numWorkers = omp_get_max_threads();
    blosc_init();
    blosc_set_nthreads(numWorkers);
    
    struct chunkInfo cI = getChunkInfo(folderName, startX, startY, startZ, endX, endY, endZ,chunkXSize,chunkYSize,chunkZSize);
    if(!cI.chunkNames) mexErrMsgIdAndTxt("zarr:inputError","File \"%s\" cannot be opened",folderName);
    //printf("nWorkers: %d\n",omp_get_max_threads());
    //mexErrMsgIdAndTxt("zarr:inputError","File \"%s\" cannot be opened",folderName);
    //printf("numChunks: %d, ChunkSizeX: %d ChunkSizeY: %d ChunkSizeZ: %d \n",cI.numChunks,chunkXSize,chunkYSize,chunkZSize);
    //mexErrMsgIdAndTxt("zarr:inputError","TESTING");
    
    int32_t batchSize = (cI.numChunks-1)/numWorkers+1;
    uint64_t s = chunkXSize*chunkYSize*chunkZSize;
    int32_t w;
    int err = 0;
    char errString[10000];
    #pragma omp parallel for
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
                //char *buffer = (char *)malloc((filelen));
                free(buffer);
                buffer = (char*) malloc(filelen);
                lastFileLen = filelen;
            }
            fread(buffer, filelen, 1, fileptr);
            fclose(fileptr);

            // Decompress
            int dsize = -1;
            switch(bits){
                case 8:
                    dsize = blosc2_decompress(buffer, filelen,bufferDest, s*sizeof(uint8_t));
                    break;
                case 16:
                    dsize = blosc2_decompress(buffer, filelen,bufferDest, s*sizeof(uint16_t));
                    break;
                case 32:
                    dsize = blosc2_decompress(buffer, filelen,bufferDest, s*sizeof(float));
                    break;
                case 64:
                    dsize = blosc2_decompress(buffer, filelen,bufferDest, s*sizeof(double));
                    break;
            }
            
            if(dsize < 0){
                #pragma omp critical
                {
                    err = 1;
                    sprintf(errString,"Decompression error. Error code: %d\n",dsize);
                }
                break;
            }
            
            //printf("ChunkName: %s\n",cI.chunkNames[f]);
            //printf("w: %d b: %d\n",w,f);
            for(int64_t x = cAV.x*chunkXSize; x < (cAV.x+1)*chunkXSize; x++){
                if(x>=endX) break;
                else if(x<startX) continue;
                for(int64_t y = cAV.y*chunkYSize; y < (cAV.y+1)*chunkYSize; y++){
                    if(y>=endY) break;
                    else if(y<startY) continue;
                    for(int64_t z = cAV.z*chunkZSize; z < (cAV.z+1)*chunkZSize; z++){
                        if(z>=endZ) break;
                        else if(z<startZ) continue;
                        switch(bits){
                            case 8:
                                ((uint8_t*)zarr)[(x-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY)] = ((uint8_t*)bufferDest)[(z%chunkZSize)+((y%chunkYSize)*chunkZSize)+((x%chunkXSize)*chunkZSize*chunkYSize)];
                                break;
                            case 16:
                                ((uint16_t*)zarr)[(x-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY)] = ((uint16_t*)bufferDest)[(z%chunkZSize)+((y%chunkYSize)*chunkZSize)+((x%chunkXSize)*chunkZSize*chunkYSize)];
                                break;
                            case 32:
                                ((float*)zarr)[(x-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY)] = ((float*)bufferDest)[(z%chunkZSize)+((y%chunkYSize)*chunkZSize)+((x%chunkXSize)*chunkZSize*chunkYSize)];
                                break;
                            case 64:
                                ((double*)zarr)[(x-startX)+((y-startY)*shapeX)+((z-startZ)*shapeX*shapeY)] = ((double*)bufferDest)[(z%chunkZSize)+((y%chunkYSize)*chunkZSize)+((x%chunkXSize)*chunkZSize*chunkYSize)];
                                break;
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

    /* After using it, destroy the Blosc environment */
    blosc_destroy();
    
    if(err) mexErrMsgIdAndTxt("zarr:threadError",errString);
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
    if(!nrhs) mexErrMsgIdAndTxt("zarr:inputError","This functions requires at least one argument");
    else if(nrhs == 2){
        if(mxGetN(prhs[1]) != 6) mexErrMsgIdAndTxt("zarr:inputError","Input range is not 6");
        startX = (uint64_t)*(mxGetPr(prhs[1]))-1;
        startY = (uint64_t)*((mxGetPr(prhs[1])+1))-1;
        startZ = (uint64_t)*((mxGetPr(prhs[1])+2))-1;
        endX = (uint64_t)*((mxGetPr(prhs[1])+3));
        endY = (uint64_t)*((mxGetPr(prhs[1])+4));
        endZ = (uint64_t)*((mxGetPr(prhs[1])+5));

        if(startX+1 < 1 || startY+1 < 1 || startZ+1 < 1) mexErrMsgIdAndTxt("zarr:inputError","Lower bounds must be at least 1");
    }
    else if (nrhs > 2) mexErrMsgIdAndTxt("zarr:inputError","Number of input arguments must be 1 or 2");
    if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("zarr:inputError","The first argument must be a string");
    char* folderName = mxArrayToString(prhs[0]);
    uint64_t shapeX = 0;
    uint64_t shapeY = 0;
    uint64_t shapeZ = 0;
    uint64_t chunkXSize = 0;
    uint64_t chunkYSize = 0;
    uint64_t chunkZSize = 0;
    char dtype[4];
    setValuesFromJSON(folderName,&chunkXSize,&chunkYSize,&chunkZSize,dtype,&shapeX,&shapeY,&shapeZ);
    if(endX > shapeX || endY > shapeY || endZ > shapeZ) mexErrMsgIdAndTxt("zarr:inputError","Upper bound is invalid");
    if(nrhs == 1){
        endX = shapeX;
        endY = shapeY;
        endZ = shapeZ;
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
        plhs[0] = mxCreateNumericArray(3,dim,mxUINT8_CLASS, mxREAL);
        uint8_t* zarr = (uint8_t*)mxGetPr(plhs[0]);
        parallelReadZarrMex((void*)zarr,folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ,bits);
    }
    else if(dtype[1] == 'u' && dtype[2] == '2'){
        uint64_t bits = 16;
        plhs[0] = mxCreateNumericArray(3,dim,mxUINT16_CLASS, mxREAL);
        uint16_t* zarr = (uint16_t*)mxGetPr(plhs[0]);
        parallelReadZarrMex((void*)zarr,folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ,bits);
    }
    else if(dtype[1] == 'f' && dtype[2] == '4'){
        uint64_t bits = 32;
        plhs[0] = mxCreateNumericArray(3,dim,mxSINGLE_CLASS, mxREAL);
        float* zarr = (float*)mxGetPr(plhs[0]);
        parallelReadZarrMex((void*)zarr,folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ,bits);
    }
    else if(dtype[1] == 'f' && dtype[2] == '8'){
        uint64_t bits = 64;
        plhs[0] = mxCreateNumericArray(3,dim,mxDOUBLE_CLASS, mxREAL);
        double* zarr = (double*)mxGetPr(plhs[0]);
        parallelReadZarrMex((void*)zarr,folderName,startX,startY,startZ,endX,endY,endZ,chunkXSize,chunkYSize,chunkZSize,shapeX,shapeY,shapeZ,bits);
    }
    else{
        mexErrMsgIdAndTxt("tiff:dataTypeError","Data type not suppported");
    }

    //plhs[0] = mxCreateNumericArray(3,dim,mxUINT16_CLASS, mxREAL);
    //uint16_t* out = (uint16_t*)mxGetPr(plhs[0]);


}
