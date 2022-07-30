#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <errno.h>
//#include <blosc.h>
#include <cjson/cJSON.h>
#include <omp.h>
#ifdef __linux__
#include <uuid/uuid.h>
#endif
#ifdef _WIN32
#include <sys/time.h>
#endif
#include <sys/stat.h>
#include "mex.h"

//compile
//mex -v COPTIMFLAGS="-O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/include/' -L'/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/lib64' -lcjson -luuid createZarrFile.c

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
#endif

static void mkdirRecursive(const char *dir) {
    char tmp[8192];
    char *p = NULL;
    size_t len;
    #ifdef __linux__
    char fileSep = '/';
    #endif
    #ifdef _WIN32
    char fileSep = '\\';
    #endif
    int status;
    snprintf(tmp, sizeof(tmp),"%s",dir);
    len = strlen(tmp);
    if (tmp[len - 1] == fileSep)
        tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++){
        if (*p == fileSep) {
            *p = 0;

            #ifdef __linux__
            mkdir(tmp, 0775);
            #endif

            #ifdef _WIN32
            mkdir(tmp);
            #endif

            chmod(tmp, 0775);
            *p = fileSep;
        }
    }
    #ifdef __linux__
    mkdir(tmp, 0775);
    #endif
    #ifdef _WIN32
    mkdir(tmp);
    #endif
    chmod(tmp, 0775);
}

void setJSONValues(char* fileName,uint64_t *chunkXSize,uint64_t *chunkYSize,uint64_t *chunkZSize,char* dtype,char* order,uint64_t *shapeX,uint64_t *shapeY,uint64_t *shapeZ, char* cname, uint64_t level){
    
    // Overflows for ints greater than 32 bit for chunkSizes and shape
    cJSON* zArray = cJSON_CreateObject();
    const int chunkSizes[3] = {*chunkXSize,*chunkYSize,*chunkZSize};
    cJSON* chunks = cJSON_CreateIntArray(chunkSizes, 3);
    cJSON_AddItemToObject(zArray, "chunks", chunks);

    cJSON* compressor = cJSON_CreateObject();
    cJSON_AddItemToObject(zArray, "compressor", compressor);
    if(!strcmp(cname,"lz4") || !strcmp(cname,"blosclz") || !strcmp(cname,"lz4hc") || !strcmp(cname,"zlib") || !strcmp(cname,"zstd")){
        cJSON_AddNumberToObject(compressor, "blocksize", 0);
        cJSON_AddNumberToObject(compressor, "clevel", level);
        cJSON_AddStringToObject(compressor, "cname", cname);
        cJSON_AddStringToObject(compressor, "id", "blosc");
        cJSON_AddNumberToObject(compressor, "shuffle", 1);
    }
    else if(!strcmp(cname,"gzip")){
        cJSON_AddStringToObject(compressor, "id", cname);
        cJSON_AddNumberToObject(compressor, "level", level);
    }
    else mexErrMsgIdAndTxt("zarr:zarrayError","Compressor: \"%s\" is not currently supported\n",cname);

    cJSON_AddStringToObject(zArray, "dtype", dtype);
    cJSON_AddNumberToObject(zArray, "fill_value", 0);
    cJSON_AddNullToObject(zArray, "filters");
    cJSON_AddStringToObject(zArray, "order", order);

    const int shapeSizes[3] = {*shapeX,*shapeY,*shapeZ};
    cJSON* shape = cJSON_CreateIntArray(shapeSizes, 3);
    cJSON_AddItemToObject(zArray, "shape", shape);
    cJSON_AddNumberToObject(zArray, "zarr_format", 2);
    
    uint64_t uuidLen;
    #ifdef __linux__
    uuidLen = 36;
    uuid_t binuuid;
    uuid_generate_random(binuuid);
    char *uuid = malloc(uuidLen+1);
    uuid_unparse(binuuid, uuid);
    #endif
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

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    // Needed variables for .zarray file
    uint64_t chunkXSize = 256;
    uint64_t chunkYSize = 256;
    uint64_t chunkZSize = 256;
    char* cname = "lz4";
    char* dtype = "<u2";
    char* order = "F";
    uint64_t shapeX = 0;
    uint64_t shapeY = 0;
    uint64_t shapeZ = 0;
    uint64_t level = 5;
    if(nrhs < 1) mexErrMsgIdAndTxt("zarr:inputError","This functions requires at least 1 argument\n");
    if(nrhs == 2) mexErrMsgIdAndTxt("zarr:inputError","This functions does not accept only 2 arguments\n");
    if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("zarr:inputError","The first argument must be a string\n");
    char* folderName = mxArrayToString(prhs[0]);

    for(int i = 1; i < nrhs; i+=2){
        if(i+1 == nrhs) mexErrMsgIdAndTxt("zarr:inputError","Mismatched argument pair for input number %d\n");
        if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("zarr:inputError","The argument in input location %d is not a string\n",i+1);
        char* currInput = mxArrayToString(prhs[i]);

        if(!strcmp(currInput,"chunks")){
            if(mxGetN(prhs[i+1]) != 3) mexErrMsgIdAndTxt("zarr:inputError","chunks must be an array of 3 numbers\n");
            chunkXSize = (uint64_t)*(mxGetPr(prhs[i+1]));
            chunkYSize = (uint64_t)*((mxGetPr(prhs[i+1])+1));
            chunkZSize = (uint64_t)*((mxGetPr(prhs[i+1])+2));
        }
        else if(!strcmp(currInput,"cname")){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","cname must be a string\n");
            cname = mxArrayToString(prhs[i+1]);
            if(!strcmp(currInput,cname) || !strcmp(currInput,cname)) mexErrMsgIdAndTxt("zarr:inputError","cname supported compressors are lz4 and zlib\n");
        }
        else if(!strcmp(currInput,"dtype")){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","dtype must be a string\n");
            dtype = mxArrayToString(prhs[i+1]);
            int dtypeLen = strlen(dtype);

            if(dtypeLen == 2){
                char* dtypeCopy = strdup(dtype);
                dtype = (char*)malloc(4);
                dtype[0] = '<';
                dtype[1] = dtypeCopy[0];
                dtype[2] = dtypeCopy[1];
                dtype[3] = '\0';
                free(dtypeCopy);
            }
            else if(dtypeLen != 3)mexErrMsgIdAndTxt("zarr:inputError","dtype must be of length 2 or 3\n");
            
            if(dtype[0] != '<' && dtype[0] != '>'){
                mexErrMsgIdAndTxt("zarr:inputError","The first character of dtype must be \"<\" or \">\"\n");
            }

            if(dtype[1] == 'u' || dtype[1] == 'f'){
                if(dtype[1] == 'u'){
                    if(dtype[2] != '1' && dtype[2] != '2'){
                        mexErrMsgIdAndTxt("zarr:inputError","For dtype u, the third character of dtype must be \"1\" or \"2\"\n");
                    }
                }
                else{
                    if(dtype[2] != '4' && dtype[2] != '8'){
                        mexErrMsgIdAndTxt("zarr:inputError","For dtype f, the third character of dtype must be \"4\" or \"8\"\n");
                    }
                }
            }
            else mexErrMsgIdAndTxt("zarr:inputError","The second character of dtype must be \"u\" or \"f\"\n");
        }
        else if(!strcmp(currInput,"order")){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","order must be a string\n");
            order = mxArrayToString(prhs[i+1]);
            if(strlen(order) > 1) mexErrMsgIdAndTxt("zarr:inputError","order must be of length 1\n");
            if(order[0] != 'F' || order[0] != 'C' || order[0] != 'f' || order[0] != 'c'){
                if(order[0] == 'f') order = "F";
                else if(order [0] == 'c') order = "C";
            }
            else mexErrMsgIdAndTxt("zarr:inputError","order must be \"F\" or \"C\"\n");
        }
        else if(!strcmp(currInput,"shape")){
            if(mxGetN(prhs[i+1]) != 3) mexErrMsgIdAndTxt("zarr:inputError","shape must be an array of 3 numbers\n");
            shapeX = (uint64_t)*(mxGetPr(prhs[i+1]));
            shapeY = (uint64_t)*((mxGetPr(prhs[i+1])+1));
            shapeZ = (uint64_t)*((mxGetPr(prhs[i+1])+2));
        }
        else if(!strcmp(currInput,"level")){
            if(!mxIsNumeric(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","level must be a numerical value\n");
            level = (uint64_t)*(mxGetPr(prhs[i+1]));
        }
        else{
            mexErrMsgIdAndTxt("zarr:inputError","The argument \"%s\" does not match the name of any supported input name.\n \
                              Currently Supported Names: chunks, cname, dtype, order, shape\n",currInput);
        }
    }
    



    char* zArray = ".zarray";
    char* fnFull = (char*)malloc(strlen(folderName)+9);
    fnFull[0] = '\0';
    char fileSepS[2];
    fileSepS[0] = '/';
    fileSepS[1] = '\0';
    
    strcat(fnFull,folderName);
    strcat(fnFull,fileSepS);
    strcat(fnFull,zArray);
    
    FILE* f = fopen(fnFull,"r");
    if(f) fclose(f);
    else{
        #ifdef __linux__
        mkdirRecursive(folderName);
        #endif
        #ifdef _WIN32
        mkdirRecursive(folderName);
        #endif
        chmod(folderName, 0775);
    }

    setJSONValues(folderName,&chunkXSize,&chunkYSize,&chunkZSize,dtype,order,&shapeX,&shapeY,&shapeZ,cname,level);


    /////////
}